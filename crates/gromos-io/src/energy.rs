//! GROMOS energy trajectory (.tre) — gromosXX's `TIMESTEP` + `ENERGY03` + `VOLUMEPRESSURE03` layout.
//!
//! The writer emits what the native binary emits (`out_configuration.cc::_print_energyred` via
//! `_print_energyred_helper`, and `_print_volumepressurered`), so gromos++ `ene_ana` reads the file
//! with the current `ene_ana.md++.lib` (`totene`, `totkin`, `totpot`, `totlj`, `totcrf`, …):
//!
//! * after `ENEVERSION`, the profile's self-description — `# energy-schema` plus one
//!   `# energy-layout` line per declaration, and `# energy-provenance` when the caller gives it
//!   ([`with_provenance`](EnergyWriter::with_provenance)) — all comment lines to gromosXX and
//!   gromos++ (`docs/src/reference/energy-library.md`).
//! * `ENERGY03`: the body [`write_energyred_body`] prints for one [`Energy`] record, shared with
//!   the `.trg` writer: the 52 `# totals` in gromosXX order, `# baths` (kinetic, centre-of-mass
//!   and internal/rotational kinetic energy per bath), `# bonded` and `# special` per energy
//!   group, `# nonbonded` per group pair (gromosXX's folded cells), and the EDS / GaMD /
//!   precalculated-λ / AB-dihedral sections.
//! * `VOLUMEPRESSURE03`: total mass, temperature per bath (total; COM / internal split and the
//!   scaling factor are not tracked: 0, 0, 1), volume, the box vectors, the pressure (the two
//!   scalars gromosXX prints after it — virial and kinetic contributions — are not tracked: 0),
//!   the pressure tensor (not tracked: zeros), the virial tensor and the molecular kinetic-energy
//!   tensor.
//!
//! Numbers use Rust's `e` exponent (`1.5e-2`, not `1.500000000e-02`); every GROMOS reader parses it.
//! [`read_energy_trajectory_native`] reads this layout (and gromosXX's own files); [`EnergyReader`]
//! still reads the one-block `ENERTRJ` layout this crate wrote before 0.0.34.

use crate::energy_traj::{builtin_schema, self_description_lines, Provenance};
use crate::IoError;
use gromos_core::configuration::Energy;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

/// Energy component identifiers (GROMOS convention)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyBlock {
    Time = 0,
    Kinetic = 1,
    Potential = 2,
    Total = 3,
    Temperature = 4,
    Volume = 5,
    Pressure = 6,

    // Bonded
    Bond = 7,
    Angle = 8,
    ImproperDihedral = 9,
    ProperDihedral = 10,

    // Nonbonded
    LennardJones = 11,
    ElectrostaticReal = 12,
    ElectrostaticReciprocal = 13,
    ElectrostaticSelf = 14,

    // Constraints
    Shake = 15,
    DistanceRestraint = 16,

    // Kinetic components
    KineticTrans = 17,
    KineticRot = 18,
    KineticInt = 19,
}

/// Energy data for a single timestep
#[derive(Debug, Clone)]
pub struct EnergyFrame {
    pub time: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub temperature: f64,
    pub volume: f64,
    pub pressure: f64,

    // Bonded energies
    pub bond: f64,
    pub angle: f64,
    pub improper: f64,
    pub dihedral: f64,

    // Nonbonded energies
    pub lj: f64,
    pub coul_real: f64,
    pub coul_recip: f64,
    pub coul_self: f64,

    // Constraints
    pub shake: f64,
    pub restraint: f64,

    // Additional components (optional)
    pub extra: Vec<f64>,
    /// MD step number (`TIMESTEP` block).
    pub step: usize,
    /// Sum of every special-interaction energy, restraints included — gromosXX's
    /// `special_total`, `ENER[21]`; the total is kinetic + potential + this.
    pub special: f64,
    /// SASA solvation energy (`ENER[22]`).
    pub sasa: f64,
    /// Total mass of the system (`VOLUMEPRESSURE03`).
    pub mass: f64,
    /// Box vectors K, L, M as rows.
    pub box_vectors: [[f64; 3]; 3],
    /// Virial tensor, `[a][b]`.
    pub virial_tensor: [[f64; 3]; 3],
    /// Molecular kinetic-energy tensor, `[a][b]`.
    pub kinetic_tensor: [[f64; 3]; 3],
}

impl Default for EnergyFrame {
    fn default() -> Self {
        Self {
            time: 0.0,
            kinetic: 0.0,
            potential: 0.0,
            total: 0.0,
            temperature: 0.0,
            volume: 0.0,
            pressure: 0.0,
            bond: 0.0,
            angle: 0.0,
            improper: 0.0,
            dihedral: 0.0,
            lj: 0.0,
            coul_real: 0.0,
            coul_recip: 0.0,
            coul_self: 0.0,
            shake: 0.0,
            restraint: 0.0,
            extra: Vec::new(),
            step: 0,
            special: 0.0,
            sasa: 0.0,
            mass: 0.0,
            box_vectors: [[0.0; 3]; 3],
            virial_tensor: [[0.0; 3]; 3],
            kinetic_tensor: [[0.0; 3]; 3],
        }
    }
}

impl EnergyFrame {
    /// Create energy frame from components
    pub fn new(time: f64, kinetic: f64, potential: f64, temperature: f64) -> Self {
        Self {
            time,
            kinetic,
            potential,
            total: kinetic + potential,
            temperature,
            ..Default::default()
        }
    }

    /// Update total energy (sum of kinetic and potential)
    pub fn update_total(&mut self) {
        self.total = self.kinetic + self.potential;
    }

    /// Update potential energy (sum of all components)
    pub fn update_potential(&mut self) {
        self.potential = self.bond
            + self.angle
            + self.improper
            + self.dihedral
            + self.lj
            + self.coul_real
            + self.coul_recip
            + self.coul_self
            + self.shake
            + self.restraint;
    }

    /// The scalar summary of the runtime [`Energy`] accumulator: which total lands in which
    /// field (`pyo3-gromos`'s `Simulation.run()` columns). The frame does not keep the record:
    /// a `.tre` frame of the engine is written from the record itself
    /// ([`EnergyWriter::write_energy`]), so per-bath and per-group values have one source.
    pub fn from_energy(
        energy: &Energy,
        time: f64,
        temperature: f64,
        volume: f64,
        pressure: f64,
    ) -> Self {
        Self {
            time,
            kinetic: energy.kinetic_total,
            potential: energy.potential_total,
            total: energy.total(),
            temperature,
            volume,
            pressure,
            bond: energy.bond_total,
            angle: energy.angle_total,
            improper: energy.improper_total,
            dihedral: energy.dihedral_total,
            lj: energy.lj_total,
            coul_real: energy.crf_total,
            coul_recip: energy.ls_kspace_total,
            coul_self: 0.0,
            shake: energy.constraint_total,
            restraint: energy.distanceres_total,
            special: energy.special_total + energy.distanceres_total,
            sasa: energy.sasa_total,
            ..Default::default()
        }
    }

    /// What the `VOLUMEPRESSURE03` block of this frame prints.
    pub fn volume_pressure(&self) -> VolumePressure {
        VolumePressure {
            mass: self.mass,
            temperature: self.temperature,
            volume: self.volume,
            pressure: self.pressure,
            box_vectors: self.box_vectors,
            virial_tensor: self.virial_tensor,
            kinetic_tensor: self.kinetic_tensor,
        }
    }

    /// The energy record derived from the scalar fields, for `n_baths` baths and `n_groups`
    /// energy groups (group 0 carries the bonded terms, the per-bath and per-pair sections are
    /// zero; the total is recomputed as kinetic + potential + special, where `special` includes
    /// `restraint`).
    pub fn to_energy(&self, n_baths: usize, n_groups: usize) -> Energy {
        let mut e = Energy::new(n_baths.max(1), n_groups.max(1));
        e.kinetic_total = self.kinetic;
        e.potential_total = self.potential;
        e.bond_total = self.bond;
        e.angle_total = self.angle;
        e.improper_total = self.improper;
        e.dihedral_total = self.dihedral;
        e.lj_total = self.lj;
        e.crf_total = self.coul_real;
        e.ls_kspace_total = self.coul_recip;
        e.constraint_total = self.shake;
        e.distanceres_total = self.restraint;
        e.special_total = self.special - self.restraint;
        e.sasa_total = self.sasa;
        e.bond_energy[0] = self.bond;
        e.angle_energy[0] = self.angle;
        e.improper_energy[0] = self.improper;
        e.dihedral_energy[0] = self.dihedral;
        e
    }
}

/// What a `VOLUMEPRESSURE03` block prints: the system's mass, its temperature (one row per
/// bath, the same value), the volume and box, the pressure, and the virial and molecular
/// kinetic-energy tensors (`[a][b]`).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct VolumePressure {
    pub mass: f64,
    pub temperature: f64,
    pub volume: f64,
    pub pressure: f64,
    /// Box vectors K, L, M as rows.
    pub box_vectors: [[f64; 3]; 3],
    pub virial_tensor: [[f64; 3]; 3],
    pub kinetic_tensor: [[f64; 3]; 3],
}

/// gromosXX's `_print_energyred_helper` (`out_configuration.cc:3068-3267`): the body of one
/// `ENERGY03` (energies) or `FREEENERDERIVS03` (λ derivatives) block, from `# totals` to the
/// `# ABdih` row. The section sizes are the record's: its baths (`kinetic_energy`) and energy
/// groups (`bond_energy`). `n_accel_groups` is the GaMD acceleration-group count of the record:
/// gromosXX sizes the energies for its default IGROUPS = 1 even without GaMD (one zero row in
/// every `.tre`) and never sizes the derivatives (none in a `.trg`).
pub fn write_energyred_body<W: Write>(
    w: &mut W,
    e: &Energy,
    n_accel_groups: usize,
) -> std::io::Result<()> {
    fn row<W: Write>(w: &mut W, vals: &[f64]) -> std::io::Result<()> {
        for v in vals {
            write!(w, "{v:>18.9e}")?;
        }
        writeln!(w)
    }
    let n = e.bond_energy.len();
    let at = |v: &Vec<f64>, i: usize| v.get(i).copied().unwrap_or(0.0);
    let cell = |m: &Vec<Vec<f64>>, i: usize, j: usize| {
        m.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0)
    };
    // gromosXX `calculate_totals` folds every pair into [max][min]: [j][i] + [i][j] (energy.cc:433-470)
    let folded = |m: &Vec<Vec<f64>>, i: usize, j: usize| {
        if i == j {
            cell(m, i, i)
        } else {
            cell(m, j, i) + cell(m, i, j)
        }
    };
    let folded_total = |m: &Vec<Vec<f64>>| {
        let mut t = 0.0;
        for i in 0..n {
            for j in i..n {
                t += folded(m, i, j);
            }
        }
        t
    };

    let total = e.total();
    let mut totals = [0.0f64; N_TOTALS];
    totals[0] = total;
    totals[1] = e.kinetic_total;
    totals[2] = e.potential_total;
    totals[3] =
        e.bond_total + e.angle_total + e.dihedral_total + e.improper_total + e.cross_dihedral_total;
    totals[4] = e.bond_total;
    totals[5] = e.angle_total;
    totals[6] = e.improper_total;
    totals[7] = e.dihedral_total;
    totals[8] = e.cross_dihedral_total;
    totals[9] = e.lj_total + e.crf_total + e.ls_realspace_total + e.ls_kspace_total;
    totals[10] = e.lj_total;
    totals[11] = e.crf_total;
    totals[12] = e.ls_total;
    totals[13] = e.ls_realspace_total + e.ls_kspace_total;
    totals[14] = e.ls_realspace_total;
    totals[15] = e.ls_kspace_total;
    totals[20] = e.special_total + e.distanceres_total;
    totals[21] = e.sasa_total;
    totals[23] = e.constraint_total;
    totals[24] = e.distanceres_total;
    totals[27] = e.posrest_total;
    totals[39] = e.entropy_term;
    totals[46] = total + folded_total(&e.shift_extra_orig);
    totals[47] = total + folded_total(&e.shift_extra_phys);
    writeln!(w, "# totals")?;
    for v in totals {
        row(w, &[v])?;
    }

    let nb = e.kinetic_energy.len();
    writeln!(w, "# baths")?;
    writeln!(w, "{nb}")?;
    for b in 0..nb {
        row(
            w,
            &[
                e.kinetic_energy[b],
                at(&e.com_kinetic_energy, b),
                at(&e.ir_kinetic_energy, b),
            ],
        )?;
    }

    writeln!(w, "# bonded")?;
    writeln!(w, "{n}")?;
    for g in 0..n {
        row(
            w,
            &[
                e.bond_energy[g],
                at(&e.angle_energy, g),
                at(&e.improper_energy, g),
                at(&e.dihedral_energy, g),
                at(&e.crossdihedral_energy, g),
            ],
        )?;
    }

    writeln!(w, "# nonbonded")?;
    for i in 0..n {
        for j in i..n {
            row(
                w,
                &[
                    folded(&e.lj_energy, i, j),
                    folded(&e.crf_energy, i, j),
                    e.ls_total,
                    folded(&e.ls_k_energy, i, j),
                    folded(&e.shift_extra_orig, i, j),
                    folded(&e.shift_extra_phys, i, j),
                ],
            )?;
        }
    }

    writeln!(w, "# special")?;
    for g in 0..n {
        row(
            w,
            &[
                at(&e.constraints_energy, g),
                at(&e.posrest_energy, g),
                at(&e.distanceres_energy, g),
                at(&e.disfieldres_energy, g),
                at(&e.angrest_energy, g),
                at(&e.dihrest_energy, g),
                at(&e.sasa_energy, g),
                at(&e.sasa_volume_energy, g),
                0.0,
                at(&e.rdc_energy, g),
                0.0,
                0.0,
                0.0,
            ],
        )?;
    }

    writeln!(w, "# eds\n# numstates\n0")?;
    writeln!(
        w,
        "           # total         nonbonded           special            offset        total_orig       total_phys"
    )?;
    writeln!(w, "# gamd\n# numaccelgroups\n{n_accel_groups}")?;
    writeln!(
        w,
        "      # E_dihedral       E_potential        K_dihedral       K_potential         dV_group"
    )?;
    for _ in 0..n_accel_groups {
        row(w, &[0.0; 5])?;
    }
    let pc = &e.precalc;
    let n_lambdas = pc.nr_lambdas();
    writeln!(w, "# precalclam\n# nr_lambdas\n{n_lambdas}")?;
    writeln!(
        w,
        "          # A_e_lj            B_e_lj           A_e_crf           B_e_crf        AB_kinetic           AB_bond          AB_angle       AB_improper         AB_disres         AB_angres         AB_dihres        AB_disfld"
    )?;
    for k in 0..n_lambdas {
        row(w, &pc.row(k))?;
    }
    writeln!(w, "# ABdih\n      # A_dihedral       B_dihedral")?;
    row(w, &[pc.a_dihedral, pc.b_dihedral])
}

/// Writer for gromosXX-layout energy trajectories (`.tre`).
pub struct EnergyWriter {
    writer: BufWriter<File>,
    frame_count: usize,
    n_baths: usize,
    n_groups: usize,
}

/// The gromosXX layout this writer emits; the `ENEVERSION` block and the built-in schema
/// (`energy_traj::builtin_schema`) are both keyed by it.
pub const ENE_VERSION: &str = "2023-04-15";

/// Number of entries in an `ENERGY03` `# totals` section (gromosXX `Energy`, 2023-04-15 layout).
const N_TOTALS: usize = 52;

impl EnergyWriter {
    /// Create the file with its `TITLE` and `ENEVERSION` blocks. Per-bath / per-group sections
    /// default to one bath and one energy group; set the run's counts with
    /// [`with_layout`](Self::with_layout).
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "TITLE")?;
        writeln!(writer, "\t{}", title)?;
        writeln!(writer, "\n\tenergy trajectory")?;
        writeln!(writer, "END")?;
        writeln!(writer, "ENEVERSION")?;
        writeln!(writer, "\t{ENE_VERSION}")?;
        writeln!(writer, "END")?;
        let schema = builtin_schema(ENE_VERSION, "ENERTRJ").expect("the layout this writer emits");
        for line in self_description_lines(&schema) {
            writeln!(writer, "{line}")?;
        }
        Ok(Self {
            writer,
            frame_count: 0,
            n_baths: 1,
            n_groups: 1,
        })
    }

    /// Record what the run was (`# energy-provenance` lines). Must precede the first frame.
    pub fn with_provenance(mut self, provenance: &Provenance) -> Result<Self, IoError> {
        debug_assert_eq!(
            self.frame_count, 0,
            "provenance goes before the first frame"
        );
        for line in provenance.lines() {
            writeln!(self.writer, "{line}")?;
        }
        Ok(self)
    }

    /// The number of temperature baths and energy groups of the run.
    pub fn with_layout(mut self, n_baths: usize, n_groups: usize) -> Self {
        self.n_baths = n_baths.max(1);
        self.n_groups = n_groups.max(1);
        self
    }

    /// Append one frame: `TIMESTEP`, `ENERGY03` and `VOLUMEPRESSURE03` blocks.
    pub fn write_frame(&mut self, frame: &EnergyFrame) -> Result<(), IoError> {
        let energy = frame.to_energy(self.n_baths, self.n_groups);
        self.write_energy(frame.step, frame.time, &energy, &frame.volume_pressure())
    }

    /// Append one frame of the engine: a `TIMESTEP` block, the `ENERGY03` block printed from
    /// `energy` (the single source of its body), and the `VOLUMEPRESSURE03` block from `vp`.
    pub fn write_energy(
        &mut self,
        step: usize,
        time: f64,
        energy: &Energy,
        vp: &VolumePressure,
    ) -> Result<(), IoError> {
        let nb = self.n_baths;
        let w = &mut self.writer;
        let num = |v: f64| format!("{:>18.9e}", v);
        let row = |vals: &[f64]| vals.iter().map(|&v| num(v)).collect::<Vec<_>>().join("");
        writeln!(w, "TIMESTEP")?;
        writeln!(w, "{:>15}{:>15.9}", step, time)?;
        writeln!(w, "END")?;

        writeln!(w, "ENERGY03")?;
        write_energyred_body(w, energy, 1)?;
        writeln!(w, "END")?;

        writeln!(w, "VOLUMEPRESSURE03")?;
        writeln!(w, "# mass")?;
        writeln!(w, "{}", num(vp.mass))?;
        writeln!(w, "# temperature")?;
        writeln!(w, "{nb}")?;
        for _ in 0..nb {
            writeln!(w, "{}", row(&[vp.temperature, 0.0, 0.0, 1.0]))?;
        }
        writeln!(w, "# volume")?;
        writeln!(w, "{}", num(vp.volume))?;
        for r in &vp.box_vectors {
            writeln!(w, "{}", row(r))?;
        }
        writeln!(w, "# pressure")?;
        writeln!(w, "{}", num(vp.pressure))?;
        writeln!(w, "{}", num(0.0))?;
        writeln!(w, "{}", num(0.0))?;
        for _ in 0..3 {
            writeln!(w, "{}", row(&[0.0; 3]))?;
        }
        for r in &vp.virial_tensor {
            writeln!(w, "{}", row(r))?;
        }
        for r in &vp.kinetic_tensor {
            writeln!(w, "{}", row(r))?;
        }
        writeln!(w, "END")?;
        self.frame_count += 1;
        Ok(())
    }

    /// Kept for callers of the pre-0.0.34 API; the native layout has no trailing block.
    pub fn finalize(&mut self) -> Result<(), IoError> {
        self.flush()
    }

    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush().map_err(IoError::Io)
    }

    pub fn frame_count(&self) -> usize {
        self.frame_count
    }
}

/// Read a native-layout energy trajectory (this crate's writer since 0.0.34, or gromosXX's own
/// files): the `# totals` of every `ENERGY03` block, the time from the preceding `TIMESTEP`,
/// volume / pressure / temperature from `VOLUMEPRESSURE03`. Comment lines are skipped.
pub fn read_energy_trajectory_native<P: AsRef<Path>>(path: P) -> Result<Vec<EnergyFrame>, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|e| IoError::FileNotFound(format!("{}: {e}", path.as_ref().display())))?;
    #[derive(PartialEq)]
    enum Block {
        None,
        Timestep,
        Energy,
        VolPres,
    }
    let mut frames: Vec<EnergyFrame> = Vec::new();
    let mut block = Block::None;
    let mut section = String::new();
    let mut step = 0usize;
    let mut time = 0.0;
    let mut totals: Vec<f64> = Vec::new();
    let mut vp: Vec<(String, Vec<f64>)> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(IoError::Io)?;
        let t = line.trim();
        match t {
            "TIMESTEP" => {
                block = Block::Timestep;
                continue;
            },
            "ENERGY03" => {
                block = Block::Energy;
                section.clear();
                totals.clear();
                continue;
            },
            "VOLUMEPRESSURE03" => {
                block = Block::VolPres;
                section.clear();
                vp.clear();
                continue;
            },
            "END" => {
                match block {
                    Block::Energy => {
                        let at = |i: usize| totals.get(i).copied().unwrap_or(0.0);
                        frames.push(EnergyFrame {
                            step,
                            time,
                            total: at(0),
                            kinetic: at(1),
                            potential: at(2),
                            bond: at(4),
                            angle: at(5),
                            improper: at(6),
                            dihedral: at(7),
                            lj: at(10),
                            coul_real: at(11),
                            coul_recip: at(13),
                            special: at(20),
                            sasa: at(21),
                            shake: at(23),
                            restraint: at(24),
                            ..Default::default()
                        });
                    },
                    Block::VolPres => {
                        if let Some(f) = frames.last_mut() {
                            for (name, vals) in &vp {
                                match (name.as_str(), vals.first()) {
                                    ("volume", Some(v)) => f.volume = *v,
                                    ("pressure", Some(v)) => f.pressure = *v,
                                    ("mass", Some(v)) => f.mass = *v,
                                    ("temperature", _)
                                        // count, then one line of 4 per bath: the first bath's total
                                        if vals.len() >= 2 => {
                                            f.temperature = vals[1];
                                        },
                                    _ => {},
                                }
                            }
                        }
                    },
                    _ => {},
                }
                block = Block::None;
                continue;
            },
            _ => {},
        }
        if t.is_empty() {
            continue;
        }
        match block {
            Block::Timestep => {
                let vals: Vec<f64> = t
                    .split_whitespace()
                    .filter_map(|x| x.parse().ok())
                    .collect();
                if vals.len() >= 2 {
                    step = vals[0] as usize;
                    time = vals[1];
                }
            },
            Block::Energy => {
                if let Some(name) = t.strip_prefix('#') {
                    section = name.trim().to_string();
                    continue;
                }
                if section.is_empty() || section == "totals" {
                    totals.extend(t.split_whitespace().filter_map(|x| x.parse::<f64>().ok()));
                }
            },
            Block::VolPres => {
                if let Some(name) = t.strip_prefix('#') {
                    vp.push((name.trim().to_string(), Vec::new()));
                    continue;
                }
                if let Some((_, vals)) = vp.last_mut() {
                    vals.extend(t.split_whitespace().filter_map(|x| x.parse::<f64>().ok()));
                }
            },
            Block::None => {},
        }
    }
    Ok(frames)
}

/// Write a single energy frame to a file (convenience function)
pub fn write_energy_frame<P: AsRef<Path>>(
    path: P,
    frame: &EnergyFrame,
    title: &str,
) -> Result<(), IoError> {
    let mut writer = EnergyWriter::new(path, title)?;
    writer.write_frame(frame)?;
    writer.finalize()?;
    Ok(())
}

/// GROMOS energy file reader
pub struct EnergyReader {
    reader: Box<dyn BufRead>,
    title: String,
    frames_read: usize,
}

impl EnergyReader {
    /// Open an energy trajectory file for reading
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, IoError> {
        let mut reader = crate::gz::open_text(path)?;

        // Read TITLE block
        let title = Self::read_title_block(&mut reader)?;

        Ok(Self {
            reader,
            title,
            frames_read: 0,
        })
    }

    /// Get the title from the energy file
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Get the number of frames read so far
    pub fn frames_read(&self) -> usize {
        self.frames_read
    }

    /// Read a single energy frame
    ///
    /// Returns Ok(Some(frame)) if a frame was read, Ok(None) if end of file reached
    pub fn read_frame(&mut self) -> Result<Option<EnergyFrame>, IoError> {
        let mut line = String::new();

        loop {
            line.clear();
            let bytes_read = self.reader.read_line(&mut line)?;

            if bytes_read == 0 {
                // End of file
                return Ok(None);
            }

            let trimmed = line.trim();

            // Skip empty lines and comments
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // Check for END block
            if trimmed == "END" {
                return Ok(None);
            }

            // Check for ENERTRJ block
            if trimmed == "ENERTRJ" {
                continue;
            }

            // Parse energy values
            let values: Result<Vec<f64>, _> = trimmed
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect();

            match values {
                Ok(vals) if vals.len() >= 17 => {
                    // Successfully parsed energy frame
                    let frame = EnergyFrame {
                        time: vals[0],
                        kinetic: vals[1],
                        potential: vals[2],
                        total: vals[3],
                        temperature: vals[4],
                        volume: vals[5],
                        pressure: vals[6],
                        bond: vals[7],
                        angle: vals[8],
                        improper: vals[9],
                        dihedral: vals[10],
                        lj: vals[11],
                        coul_real: vals[12],
                        coul_recip: vals[13],
                        coul_self: vals[14],
                        shake: vals[15],
                        restraint: vals[16],
                        extra: if vals.len() > 17 {
                            vals[17..].to_vec()
                        } else {
                            Vec::new()
                        },
                        ..Default::default()
                    };

                    self.frames_read += 1;
                    return Ok(Some(frame));
                },
                Ok(vals) if vals.is_empty() => {
                    // Empty line after splitting, continue
                    continue;
                },
                Ok(vals) => {
                    return Err(IoError::ParseError(format!(
                        "Incomplete energy frame: expected at least 17 values, got {}",
                        vals.len()
                    )));
                },
                Err(e) => {
                    return Err(IoError::ParseError(format!(
                        "Failed to parse energy values: {}",
                        e
                    )));
                },
            }
        }
    }

    /// Read all energy frames from the file
    pub fn read_all_frames(&mut self) -> Result<Vec<EnergyFrame>, IoError> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_frame()? {
            frames.push(frame);
        }
        Ok(frames)
    }

    /// Read TITLE block
    fn read_title_block(reader: &mut dyn BufRead) -> Result<String, IoError> {
        let mut line = String::new();
        let mut title_lines = Vec::new();
        let mut in_title = false;

        loop {
            line.clear();
            let bytes_read = reader.read_line(&mut line)?;

            if bytes_read == 0 {
                return Err(IoError::ParseError(
                    "Unexpected EOF while reading TITLE".to_string(),
                ));
            }

            let trimmed = line.trim();

            if trimmed == "TITLE" {
                in_title = true;
                continue;
            }

            if trimmed == "END" {
                if in_title {
                    return Ok(title_lines.join(" ").trim().to_string());
                }
                continue;
            }

            if in_title {
                title_lines.push(trimmed.to_string());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The shared body prints gromosXX's layout: folded group pairs, the global lattice-sum
    /// total in column 3, `total + shift_extra` in slots 47/48, the GaMD count it is given.
    #[test]
    fn the_energyred_body_is_gromosxx_layout() {
        let mut e = Energy::new(2, 2);
        e.kinetic_total = 10.0;
        e.potential_total = -4.0;
        e.lj_energy[0][1] = 1.0;
        e.lj_energy[1][0] = 2.0;
        e.crf_energy[1][1] = -3.0;
        e.kinetic_energy = vec![6.0, 4.0];
        e.com_kinetic_energy = vec![1.0, 0.5];
        e.ir_kinetic_energy = vec![5.0, 3.5];
        e.posrest_energy[1] = 0.25;
        let mut out = Vec::new();
        write_energyred_body(&mut out, &e, 1).unwrap();
        let text = String::from_utf8(out).unwrap();
        let rows = |section: &str| -> Vec<Vec<f64>> {
            text.lines()
                .skip_while(|l| l.trim() != section)
                .skip(1)
                .take_while(|l| !l.trim_start().starts_with('#'))
                .map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect())
                .collect()
        };
        let totals = rows("# totals");
        assert_eq!(totals.len(), 52);
        assert_eq!(
            (totals[0][0], totals[46][0], totals[47][0]),
            (6.0, 6.0, 6.0)
        );
        assert_eq!(
            rows("# baths"),
            vec![vec![2.0], vec![6.0, 1.0, 5.0], vec![4.0, 0.5, 3.5]]
        );
        let nb = rows("# nonbonded");
        assert_eq!(nb.len(), 3); // (0,0) (0,1) (1,1)
        assert_eq!((nb[1][0], nb[2][1]), (3.0, -3.0));
        assert_eq!(rows("# special")[1][1], 0.25);
        assert!(text.contains("# numaccelgroups\n1\n"));
        let mut out = Vec::new();
        write_energyred_body(&mut out, &e, 0).unwrap();
        assert!(String::from_utf8(out)
            .unwrap()
            .contains("# numaccelgroups\n0\n"));
    }

    #[test]
    fn test_energy_frame_creation() {
        let frame = EnergyFrame::new(0.0, 100.0, -200.0, 300.0);
        assert_eq!(frame.time, 0.0);
        assert_eq!(frame.kinetic, 100.0);
        assert_eq!(frame.potential, -200.0);
        assert_eq!(frame.total, -100.0);
        assert_eq!(frame.temperature, 300.0);
    }

    #[test]
    fn test_energy_writer_creation() {
        let temp_file = "/tmp/test_energy.tre";
        let writer = EnergyWriter::new(temp_file, "Test energy");
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_energy_frame() {
        let temp_file = "/tmp/test_energy2.tre";
        let mut writer = EnergyWriter::new(temp_file, "Test").unwrap();

        let mut frame = EnergyFrame::new(0.0, 150.0, -300.0, 298.15);
        frame.bond = -50.0;
        frame.angle = -30.0;
        frame.lj = -150.0;
        frame.coul_real = -70.0;
        frame.update_potential();

        let result = writer.write_frame(&frame);
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);

        let result = writer.finalize();
        assert!(result.is_ok());
    }

    /// What this writer emits is read back under tier a — the file's own layout, verified —
    /// with the official library, and satisfies the first-frame identity with a non-zero
    /// special energy (the slot it lands in is the one `ene_ana`'s `totspecial` reads).
    #[test]
    fn written_file_describes_itself_and_adds_up() {
        use crate::energy_traj::{read_preamble, EnergyTraj, Tier, OFFICIAL_LIBRARY};
        let path =
            std::env::temp_dir().join(format!("gromos_tre_roundtrip_{}.tre", std::process::id()));
        let mut writer = EnergyWriter::new(&path, "roundtrip")
            .unwrap()
            .with_layout(2, 2)
            .with_provenance(&Provenance {
                writer: "gromos-rs test".into(),
                topology: Some(("sha256:00".into(), "x.top".into())),
                energy_groups: vec![(1, 10), (11, 30)],
            })
            .unwrap();
        let mut frame = EnergyFrame::new(0.002, 150.0, -300.0, 298.15);
        frame.step = 1;
        frame.special = 7.5;
        frame.restraint = 2.5;
        frame.shake = 0.25;
        frame.total = frame.kinetic + frame.potential + frame.special;
        frame.volume = 27.0;
        writer.write_frame(&frame).unwrap();
        writer.finalize().unwrap();

        let mut reader = crate::gz::TextReader::open(&path).unwrap();
        let pre = read_preamble(&mut reader).unwrap();
        assert_eq!(pre.version.as_deref(), Some(ENE_VERSION));
        assert_eq!(pre.schema, builtin_schema(ENE_VERSION, "ENERTRJ"));
        assert_eq!(pre.writer(), Some("gromos-rs test"));
        assert_eq!(
            pre.provenance[2],
            ("energy-groups".to_string(), "1-10 11-30".to_string())
        );

        let mut lib = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
        assert_eq!(
            lib.bind("ENERTRJ", &pre, "official", "roundtrip.tre")
                .unwrap(),
            Tier::SelfDescribed
        );
        assert!(lib.read_frame(&mut reader, "ENERTRJ").unwrap());
        assert_eq!(lib.drain_warnings(), Vec::<String>::new());
        assert_eq!(lib.value("totene").unwrap(), -142.5);
        assert_eq!(lib.value("totspecial").unwrap(), 7.5);
        assert_eq!(lib.value("totdisres").unwrap(), 2.5);
        assert_eq!(lib.value("totconstraint").unwrap(), 0.25);
        assert_eq!(lib.value("boxvol").unwrap(), 27.0);

        let back = read_energy_trajectory_native(&path).unwrap();
        assert_eq!(
            (back[0].special, back[0].restraint, back[0].shake),
            (7.5, 2.5, 0.25)
        );
        std::fs::remove_file(&path).ok();
    }

    /// The PRECALCLAM rows and the ABdih row of a `.tre` and a `.trg` written by gromos-rs read
    /// back through the library-bound reader (`ext_ti_ana`'s input) with the values written.
    #[test]
    fn precalclam_rows_read_back_through_the_library() {
        use crate::energy_traj::{read_preamble, EnergyTraj, OFFICIAL_LIBRARY};
        use crate::free_energy::FreeEnergyWriter;
        let mut e = EnergyFrame::new(0.0, 1.0, 2.0, 300.0).to_energy(1, 2);
        e.precalc.resize(3, 2);
        for k in 0..3 {
            let pc = &mut e.precalc;
            pc.a_lj_total[k] = 10.0 + k as f64;
            pc.b_lj_total[k] = 20.0 + k as f64;
            pc.a_crf_total[k] = 30.0 + k as f64;
            pc.b_crf_total[k] = 40.0 + k as f64;
            pc.ab_kinetic[k] = 50.0 + k as f64;
            pc.ab_bond[k] = 60.0 + k as f64;
            pc.ab_angle[k] = 70.0 + k as f64;
            pc.ab_improper[k] = 80.0 + k as f64;
            pc.ab_disres[k] = 90.0 + k as f64;
        }
        e.precalc.a_dihedral = 1.5;
        e.precalc.b_dihedral = -2.5;
        let expect = |rows: &[Vec<f64>]| {
            assert_eq!(rows.len(), 3);
            for (k, row) in rows.iter().enumerate() {
                assert_eq!(row, &e.precalc.row(k).to_vec(), "row {k}");
            }
        };
        let dir = std::env::temp_dir();
        let (tre, trg) = (
            dir.join(format!("gromos_precalc_{}.tre", std::process::id())),
            dir.join(format!("gromos_precalc_{}.trg", std::process::id())),
        );
        let mut w = EnergyWriter::new(&tre, "precalc")
            .unwrap()
            .with_layout(1, 2);
        w.write_energy(0, 0.0, &e, &VolumePressure::default())
            .unwrap();
        w.finalize().unwrap();
        let mut fw = FreeEnergyWriter::new(&trg, "precalc")
            .unwrap()
            .with_layout(1, 2);
        fw.write_derivatives(0, 0.0, 0.25, &e).unwrap();
        fw.flush().unwrap();
        drop(fw);

        for (path, ft, block, abdih) in [
            (&tre, "ENERTRJ", "PRECALCLAM", "ABDIH"),
            (&trg, "FRENERTRJ", "FREEPRECALCLAM", "FREEABDIH"),
        ] {
            let mut reader = crate::gz::TextReader::open(path).unwrap();
            let pre = read_preamble(&mut reader).unwrap();
            let mut lib = EnergyTraj::from_library(OFFICIAL_LIBRARY).unwrap();
            lib.bind(ft, &pre, "official", "precalc").unwrap();
            assert!(lib.read_frame(&mut reader, ft).unwrap(), "{ft}");
            assert_eq!(lib.size("NUM_LAMBDAS"), Some(3), "{ft}");
            expect(lib.table(block).unwrap());
            assert_eq!(lib.table(abdih).unwrap(), &[vec![1.5, -2.5]], "{ft}");
            std::fs::remove_file(path).ok();
        }
    }
}
