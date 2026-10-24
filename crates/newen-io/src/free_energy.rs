//! Free energy trajectory (.trg) — gromosXX's `TIMESTEP` + `FREEENERDERIVS03` layout.
//!
//! The writer emits what the native binary emits (`out_configuration.cc::_print_free_energyred`
//! via `_print_energyred_helper`): a `TIMESTEP` block, then `FREEENERDERIVS03` with `# lambda`
//! and the body the `.tre` writer prints too ([`write_energyred_body`]), from the record of λ
//! derivatives — so gromos++ `ene_ana`/`ext_ti_ana` read the file, and the reference suite
//! compares dH/dλ value by value. Numbers are printed with Rust's `e` exponent (`1.5e-2`, not
//! `1.500000000e-02`), which every GROMOS reader parses.
//!
//! After `ENEVERSION` the writer adds the profile's self-description (`# energy-schema`,
//! `# energy-layout`, optionally `# energy-provenance`): comment lines to every other reader
//! (`docs/src/reference/energy-library.md`).
//!
//! The reader goes through the energy library, so a file that does not match the layout it
//! claims is refused instead of yielding dH/dλ from the wrong slot; the one-line `FREEENERGY03`
//! block older files carry is still read, on its own path.
//!
//! dH/dλ (kJ/mol) is the thermodynamic-integration integrand at the current λ. Integrate
//! ⟨dH/dλ⟩_λ over λ with `ext_ti_ana` to get ΔG.

use crate::energy::{write_energyred_body, ENE_VERSION};
use crate::energy_traj::{
    builtin_schema, read_preamble, self_description_lines, EnergyTraj, Provenance, OFFICIAL_LIBRARY,
};
use crate::gz::TextReader;
use crate::IoError;
use newen_core::configuration::Energy;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// One frame of free-energy derivative data.
#[derive(Debug, Clone, Default)]
pub struct FreeEnergyFrame {
    /// MD step number (the `TIMESTEP` block).
    pub step: usize,
    pub time: f64,
    pub lambda: f64,
    // Per-component dH/dλ (kJ/mol); zero if not computed separately.
    pub dhdl_bond: f64,
    pub dhdl_angle: f64,
    pub dhdl_improper: f64,
    pub dhdl_dihedral: f64,
    pub dhdl_lj: f64,
    pub dhdl_crf: f64,
    pub dhdl_special: f64,
    /// Total dH/dλ — the quantity consumed by ext_ti_ana / TI integration.
    pub dhdl_total: f64,
}

impl FreeEnergyFrame {
    pub fn new(time: f64, lambda: f64, dhdl_total: f64) -> Self {
        Self {
            time,
            lambda,
            dhdl_total,
            ..Default::default()
        }
    }

    /// The record of λ derivatives derived from the scalar fields, for `n_baths` baths and
    /// `n_groups` energy groups (potential = total − special, no kinetic derivative).
    pub fn to_energy(&self, n_baths: usize, n_groups: usize) -> Energy {
        let mut d = Energy::new(n_baths.max(1), n_groups.max(1));
        d.bond_total = self.dhdl_bond;
        d.angle_total = self.dhdl_angle;
        d.improper_total = self.dhdl_improper;
        d.dihedral_total = self.dhdl_dihedral;
        d.lj_total = self.dhdl_lj;
        d.crf_total = self.dhdl_crf;
        d.special_total = self.dhdl_special;
        d.potential_total = self.dhdl_total - self.dhdl_special;
        d
    }
}

/// A `.trg` file as read: its frames and whatever the layout check had to say about them.
#[derive(Debug, Clone, Default)]
pub struct FreeEnergyTrajectory {
    pub frames: Vec<FreeEnergyFrame>,
    /// Layout warnings, in the words the profile gives them
    /// (`docs/src/reference/energy-library.md`) — an unverified layout, a version mismatch.
    /// Empty for a file whose layout is established.
    pub warnings: Vec<String>,
}

/// Read every frame of a `.trg` file, in either layout:
///
/// * gromosXX's `FREEENERDERIVS03`, read through the energy library like `ene_ana` reads it:
///   the file's own self-description if it has one, else the built-in layout for its
///   `ENEVERSION`, else the official library on trust with a warning. A file that does not
///   match the layout it claims is an error here rather than dH/dλ from the wrong slot
///   (`docs/src/reference/energy-library.md`).
/// * the one-line `FREEENERGY03` block this crate wrote before 0.0.33 (time, λ, per-term
///   dH/dλ, total), so older files keep reading. It predates the profile and gets none of
///   its checks; nothing writes it any more.
///
/// Returns frames in order.
pub fn read_free_energy_trajectory<P: AsRef<Path>>(
    path: P,
) -> Result<FreeEnergyTrajectory, IoError> {
    let path = path.as_ref();
    let mut reader = TextReader::open(path)
        .map_err(|e| IoError::FileNotFound(format!("{}: {e}", path.display())))?;
    let fail = |e: String| IoError::FormatError(format!("{}: {e}", path.display()));

    let preamble = read_preamble(&mut reader).map_err(fail)?;
    if peek_block(&mut reader).map_err(IoError::Io)?.as_deref() == Some("FREEENERGY03") {
        return Ok(FreeEnergyTrajectory {
            frames: read_legacy_free_energy(&mut reader).map_err(IoError::Io)?,
            warnings: Vec::new(),
        });
    }

    let mut etrj = EnergyTraj::from_library(OFFICIAL_LIBRARY)
        .map_err(|e| IoError::FormatError(format!("built-in energy library: {e}")))?;
    etrj.bind(
        "FRENERTRJ",
        &preamble,
        "<built-in library>",
        &path.display().to_string(),
    )
    .map_err(fail)?;

    let mut frames = Vec::new();
    let mut warnings = etrj.drain_warnings();
    while etrj.read_frame(&mut reader, "FRENERTRJ").map_err(fail)? {
        let at = |i: usize| etrj.value(&format!("FREEENER[{i}]")).unwrap_or(0.0);
        frames.push(FreeEnergyFrame {
            step: etrj.value("TIME[1]").unwrap_or(0.0) as usize,
            time: etrj.value("TIME[2]").unwrap_or(0.0),
            lambda: etrj.value("RLAM[1]").unwrap_or(0.0),
            dhdl_total: at(1),
            dhdl_bond: at(5),
            dhdl_angle: at(6),
            dhdl_improper: at(7),
            dhdl_dihedral: at(8),
            dhdl_lj: at(11),
            dhdl_crf: at(12),
            dhdl_special: at(21),
        });
        warnings.extend(etrj.drain_warnings());
    }
    Ok(FreeEnergyTrajectory { frames, warnings })
}

/// The name of the block the reader is positioned on, put back for the next reader.
fn peek_block(reader: &mut TextReader) -> std::io::Result<Option<String>> {
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let t = line.trim();
        if !t.is_empty() {
            reader.unread_line(&line);
            return Ok(Some(t.to_string()));
        }
    }
    Ok(None)
}

/// The one-line `FREEENERGY03` block: `time λ bond angle improper dihedral lj crf total`, or
/// the three-value `time λ total` form. Predates the profile; no layout check applies.
fn read_legacy_free_energy(reader: &mut TextReader) -> std::io::Result<Vec<FreeEnergyFrame>> {
    let mut frames = Vec::new();
    let mut inside = false;
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let t = line.trim();
        match t {
            "FREEENERGY03" => inside = true,
            "END" => inside = false,
            _ if !inside || t.is_empty() || t.starts_with('#') => {},
            _ => {
                let v: Vec<f64> = t
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if v.len() >= 9 {
                    frames.push(FreeEnergyFrame {
                        step: 0,
                        time: v[0],
                        lambda: v[1],
                        dhdl_bond: v[2],
                        dhdl_angle: v[3],
                        dhdl_improper: v[4],
                        dhdl_dihedral: v[5],
                        dhdl_lj: v[6],
                        dhdl_crf: v[7],
                        dhdl_special: 0.0,
                        dhdl_total: v[8],
                    });
                } else if v.len() == 3 {
                    frames.push(FreeEnergyFrame::new(v[0], v[1], v[2]));
                }
            },
        }
    }
    Ok(frames)
}

/// Writer for gromosXX-layout free-energy trajectory files (.trg).
pub struct FreeEnergyWriter {
    writer: BufWriter<File>,
    n_baths: usize,
    n_groups: usize,
}

impl FreeEnergyWriter {
    /// Create the file with its `TITLE` and `ENEVERSION` blocks. The per-bath / per-group
    /// sections default to one bath and one energy group; set the run's counts with
    /// [`with_layout`](Self::with_layout).
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "TITLE")?;
        writeln!(writer, "\t{}", title)?;
        writeln!(writer, "\n\tfree energy trajectory")?;
        writeln!(writer, "END")?;
        writeln!(writer, "ENEVERSION")?;
        writeln!(writer, "\t{ENE_VERSION}")?;
        writeln!(writer, "END")?;
        let schema =
            builtin_schema(ENE_VERSION, "FRENERTRJ").expect("the layout this writer emits");
        for line in self_description_lines(&schema) {
            writeln!(writer, "{line}")?;
        }
        Ok(Self {
            writer,
            n_baths: 1,
            n_groups: 1,
        })
    }

    /// Record what the run was (`# energy-provenance` lines). Must precede the first frame.
    pub fn with_provenance(mut self, provenance: &Provenance) -> Result<Self, IoError> {
        for line in provenance.lines() {
            writeln!(self.writer, "{line}")?;
        }
        Ok(self)
    }

    /// The number of temperature baths and energy groups of the run (sizes of the per-bath
    /// and per-group sections, as in the native file).
    pub fn with_layout(mut self, n_baths: usize, n_groups: usize) -> Self {
        self.n_baths = n_baths.max(1);
        self.n_groups = n_groups.max(1);
        self
    }

    /// Append one frame: a `TIMESTEP` block and a `FREEENERDERIVS03` block.
    pub fn write_frame(&mut self, frame: &FreeEnergyFrame) -> Result<(), IoError> {
        let derivatives = frame.to_energy(self.n_baths, self.n_groups);
        self.write_derivatives(frame.step, frame.time, frame.lambda, &derivatives)
    }

    /// Append one frame of the engine: a `TIMESTEP` block and a `FREEENERDERIVS03` block printed
    /// from `derivatives` (the engine's `perturbed_energy_derivatives`, the single source of its
    /// body).
    pub fn write_derivatives(
        &mut self,
        step: usize,
        time: f64,
        lambda: f64,
        derivatives: &Energy,
    ) -> Result<(), IoError> {
        let w = &mut self.writer;
        writeln!(w, "TIMESTEP")?;
        writeln!(w, "{:>15}{:>15.9}", step, time)?;
        writeln!(w, "END")?;
        writeln!(w, "FREEENERDERIVS03")?;
        writeln!(w, "# lambda")?;
        writeln!(w, "{:>18.9e}", lambda)?;
        write_energyred_body(w, derivatives, &[])?;
        writeln!(w, "END")?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush().map_err(IoError::Io)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames_survive_write_then_read() {
        let p = std::env::temp_dir().join(format!("newen_trg_test_{}.trg", std::process::id()));
        {
            let mut w = FreeEnergyWriter::new(&p, "round trip")
                .unwrap()
                .with_layout(2, 3);
            let mut f = FreeEnergyFrame::new(0.002, 0.25, -12.5);
            f.dhdl_lj = -10.0;
            f.dhdl_crf = -2.5;
            w.write_frame(&f).unwrap();
            w.write_frame(&FreeEnergyFrame::new(0.004, 0.25, -11.0))
                .unwrap();
            w.flush().unwrap();
        }
        let traj = read_free_energy_trajectory(&p).unwrap();
        let frames = traj.frames;
        // written with the profile's self-description, so the layout is established, not assumed
        assert!(traj.warnings.is_empty(), "{:?}", traj.warnings);
        assert_eq!(frames.len(), 2);
        assert!((frames[0].time - 0.002).abs() < 1e-12);
        assert!((frames[0].lambda - 0.25).abs() < 1e-12);
        assert!((frames[0].dhdl_total + 12.5).abs() < 1e-9);
        assert!((frames[0].dhdl_lj + 10.0).abs() < 1e-9);
        assert!((frames[1].dhdl_total + 11.0).abs() < 1e-9);
        std::fs::remove_file(p).ok();
    }

    /// The one-line block this crate wrote before 0.0.33: no `ENEVERSION`, no schema, its own
    /// path through the reader.
    #[test]
    fn the_pre_0_0_33_block_still_reads() {
        let p = std::env::temp_dir().join(format!("newen_trg_legacy_{}.trg", std::process::id()));
        std::fs::write(
            &p,
            "TITLE\n\tan old file\nEND\nFREEENERGY03\n# time lambda ...\n\
             0.0 0.5 1.0 2.0 3.0 4.0 5.0 6.0 21.0\n\
             0.002 0.5 1.0 2.0 3.0 4.0 5.0 6.0 23.0\nEND\n",
        )
        .unwrap();
        let traj = read_free_energy_trajectory(&p).unwrap();
        assert_eq!(traj.frames.len(), 2);
        assert!((traj.frames[1].dhdl_total - 23.0).abs() < 1e-9);
        assert!((traj.frames[0].dhdl_lj - 5.0).abs() < 1e-9);
        std::fs::remove_file(p).ok();
    }

    /// A `.trg` that claims a layout it does not have is refused, where the fixed-slot reader
    /// used to hand back dH/dλ from whatever landed in the slot.
    #[test]
    fn a_shifted_table_is_refused_not_misread() {
        let p = std::env::temp_dir().join(format!("newen_trg_short_{}.trg", std::process::id()));
        {
            let mut w = FreeEnergyWriter::new(&p, "one total short").unwrap();
            w.write_frame(&FreeEnergyFrame::new(0.002, 0.25, -12.5))
                .unwrap();
            w.flush().unwrap();
        }
        let text = std::fs::read_to_string(&p).unwrap();
        let mut values = 0;
        let short: String = text
            .lines()
            .filter(|l| {
                let keep = l.trim().starts_with('#') || l.trim().parse::<f64>().is_err() || {
                    values += 1;
                    values != 53 // one value of FREEENER dropped (RLAM is value 1)
                };
                keep
            })
            .fold(String::new(), |mut acc, l| {
                acc.push_str(l);
                acc.push('\n');
                acc
            });
        std::fs::write(&p, short).unwrap();

        let err = read_free_energy_trajectory(&p).unwrap_err().to_string();
        assert!(err.contains("falls inside subblock FREEENER"), "{err}");
        std::fs::remove_file(p).ok();
    }
}
