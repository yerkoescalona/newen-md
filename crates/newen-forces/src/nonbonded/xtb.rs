//! `XtbInteraction` — a real `xtb` (GFN-xTB) semiempirical QM engine as a [`PotentialProvider`].
//!
//! Runs `xtb` as a subprocess per [`PotentialProvider::contribute`] call (file-in/file-out),
//! writing an `.xyz` geometry and parsing the `.engrad` output for energy and gradient. `xtb`'s
//! own xyz parser accepts a bare atomic number in the element column (verified against a real
//! run), so atom identity is passed the same way [`crate::nonbonded::SchNetInteraction`] takes
//! `elements` — a caller-supplied table, no symbol lookup needed here.
//!
//! Two embedding modes:
//! - [`Embedding::None`] (default): the region's own isolated-cluster energy only.
//! - [`Embedding::Electrostatic`]: **real** electrostatic embedding — MM point charges within
//!   [`Self::with_cutoff`]'s radius are written to xtb's `pcharge` file (format and *units*
//!   verified empirically against a real run: coordinates in **Bohr**, matching every other raw
//!   xtb I/O file — the man page documents the columns but not the unit, and guessing wrong here
//!   would have been a silent, hard-to-detect error). xtb's own SCF is polarized by them, and
//!   returns MM-atom forces directly via `pcgrad` (one row per point charge, same order,
//!   Eh/bohr) — this is Poliak et al. 2025's coupling path (a): *the QM program hands back MM
//!   forces itself*, which is what `xtb_worker.cc` in real gromosXX does
//!   (`xtb_setExternalCharges`/`xtb_getPCGradient`). No charge-derivative ("dq/dR") correction
//!   term is needed because xtb differentiates the whole coupled energy itself — confirmed by
//!   reading gromosXX's actual C++ source, which implements no such term anywhere.
//!
//! **Double-counting hazard.** This is an *alternative* to
//! [`crate::nonbonded::ElectrostaticEmbedding`] (Poliak path (c): the host computes Coulomb from
//! caller-supplied QM charges), not something to combine with it — registering both on the same
//! QM/MM boundary computes the coupling twice. Path (a) here is strictly more correct when the
//! engine supports it (real SCF response, not a fixed-charge approximation); path (c) remains
//! necessary for engines that don't return MM gradients (e.g. MOPAC).
//!
//! [`Embedding::Mechanical`] is still rejected — nothing in this crate implements it.

use newen_core::configuration::Configuration;
use newen_core::math::{Mat3, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::SpatialIndex;
use newen_core::topology::Topology;
use newen_core::units::{BOHR_TO_NM, HARTREE_TO_KJMOL, NM_TO_ANGSTROM};
// Only the test module spawns the real `xtb` binary; the library goes through
// `qm_subprocess`, so an unconditional import is reported unused.
use std::io::Write;
use std::path::PathBuf;
#[cfg(test)]
use std::process::Command;

use super::qm_subprocess;
use crate::provider::{Contribution, Embedding, PotentialProvider, ProviderError, ProviderExtra};

/// Runs a real `xtb` binary as a subprocess to get a QM region's isolated energy and forces.
pub struct XtbInteraction {
    binary: String,
    gfn: u8,
    charge: i32,
    multiplicity: i32,
    /// Atomic number per global atom index — the caller's job to map from `Topology`, same
    /// convention as `SchNetInteraction::elements`.
    atomic_numbers: Vec<i64>,
    work_dir: PathBuf,
    embedding: Embedding,
    /// GROMOS `RCUTQM` (nm) — cutoff for gathering MM point charges when `embedding` is
    /// `Electrostatic`. Unused otherwise. Must be set via `with_cutoff` before an `Electrostatic`
    /// call; `contribute()` errors clearly rather than silently using a zero cutoff.
    cutoff: f64,
    /// Wall-clock bound per `xtb` call (`with_timeout`); `None` waits forever.
    timeout: Option<std::time::Duration>,
}

impl XtbInteraction {
    /// `work_dir` is created if missing and reused (overwritten) on every `contribute()` call —
    /// no temp-directory dependency, same convention as the older `qmmm.rs::XTBEngine`.
    pub fn new(
        work_dir: impl Into<PathBuf>,
        gfn: u8,
        charge: i32,
        multiplicity: i32,
        atomic_numbers: Vec<i64>,
    ) -> Result<Self, ProviderError> {
        let work_dir = qm_subprocess::ensure_work_dir(work_dir)?;
        Ok(Self {
            binary: "xtb".to_string(),
            gfn,
            charge,
            multiplicity,
            atomic_numbers,
            work_dir,
            embedding: Embedding::None,
            cutoff: 0.0,
            timeout: None,
        })
    }

    /// Bound every `xtb` call: past `timeout` the subprocess is killed and `contribute()` fails
    /// with a clear error instead of stalling the MD step.
    pub fn with_timeout(mut self, timeout: std::time::Duration) -> Self {
        self.timeout = Some(timeout);
        self
    }

    /// Override the `xtb` executable (default: `"xtb"`, resolved via `PATH`).
    pub fn with_binary(mut self, binary: impl Into<String>) -> Self {
        self.binary = binary.into();
        self
    }

    /// Declare a non-default embedding scheme. `Electrostatic` is real (see module docs);
    /// `Mechanical` is accepted here but rejected at `contribute()` time with a clear error,
    /// same pattern as `SchNetInteraction::with_embedding`.
    pub fn with_embedding(mut self, embedding: Embedding) -> Self {
        self.embedding = embedding;
        self
    }

    /// Set the MM point-charge gathering cutoff (GROMOS `RCUTQM`, nm) for `Embedding::Electrostatic`.
    pub fn with_cutoff(mut self, cutoff: f64) -> Self {
        self.cutoff = cutoff;
        self
    }

    fn write_xyz(&self, path: &std::path::Path, atom_indices: &[usize], conf: &Configuration) {
        let mut file = std::fs::File::create(path).expect("xyz path is in a writable work_dir");
        writeln!(file, "{}", atom_indices.len()).unwrap();
        writeln!(file, "newen XtbInteraction").unwrap();
        for &g in atom_indices {
            let p = conf.current().pos[g] * NM_TO_ANGSTROM;
            writeln!(
                file,
                "{} {:.10} {:.10} {:.10}",
                self.atomic_numbers[g], p.x, p.y, p.z
            )
            .unwrap();
        }
    }

    /// Write xtb's `pcharge` file: count line, then `charge x y z` per MM atom, coordinates in
    /// **Bohr** (verified empirically — see module docs).
    fn write_pcharge(
        &self,
        path: &std::path::Path,
        mm_atoms: &[usize],
        topo: &Topology,
        conf: &Configuration,
    ) {
        let mut file = std::fs::File::create(path).expect("pcharge path is in a writable work_dir");
        writeln!(file, "{}", mm_atoms.len()).unwrap();
        for &i in mm_atoms {
            let p = conf.current().pos[i] / BOHR_TO_NM;
            writeln!(
                file,
                "{:.10} {:.10} {:.10} {:.10}",
                topo.charge[i], p.x, p.y, p.z
            )
            .unwrap();
        }
    }

    /// Gather MM atoms within `self.cutoff` of `region`, deduplicated — mirrors
    /// `ElectrostaticEmbedding::contribute`'s exact cross-boundary-pair gathering pattern.
    fn gather_mm_atoms(
        &self,
        region: &AtomSelection,
        atom_indices: &[usize],
        topo: &Topology,
        neigh: &dyn SpatialIndex,
    ) -> Vec<usize> {
        let n_atoms = topo.num_atoms();
        let mut in_region = vec![false; n_atoms];
        for &i in atom_indices {
            in_region[i] = true;
        }
        let mut touched = vec![false; n_atoms];
        let mut mm_atoms = Vec::new();
        for (lo, hi, _shift) in neigh.neighbor_pairs(region, self.cutoff) {
            let mm = match (in_region[lo], in_region[hi]) {
                (true, false) => hi,
                (false, true) => lo,
                _ => continue, // QM-QM (xtb's own business) — never both-outside per SpatialIndex's contract.
            };
            if !touched[mm] {
                touched[mm] = true;
                mm_atoms.push(mm);
            }
        }
        mm_atoms
    }
}

/// Parse xtb's `pcgrad`: one `gx gy gz` row (Eh/bohr) per point charge, same order as written to
/// `pcharge`.
fn parse_pcgrad(path: &std::path::Path, expected_n: usize) -> Result<Vec<Vec3>, ProviderError> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| ProviderError::ComputationFailed(format!("failed to read pcgrad: {e}")))?;
    let mut gradients = Vec::new();
    for line in content.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let parts: Vec<&str> = t.split_whitespace().collect();
        if parts.len() != 3 {
            return Err(ProviderError::ComputationFailed(format!(
                "pcgrad line has {} fields, expected 3: {t}",
                parts.len()
            )));
        }
        let parse = |s: &str| {
            s.parse::<f64>().map_err(|e| {
                ProviderError::ComputationFailed(format!("failed to parse pcgrad component: {e}"))
            })
        };
        gradients.push(Vec3::new(
            parse(parts[0])?,
            parse(parts[1])?,
            parse(parts[2])?,
        ));
    }
    if gradients.len() != expected_n {
        return Err(ProviderError::ComputationFailed(format!(
            "pcgrad has {} rows, expected {expected_n} (one per MM point charge)",
            gradients.len()
        )));
    }
    Ok(gradients)
}

/// Parse xtb's `.engrad` (ORCA-style) output: atom count, energy (Eh), then the flattened
/// gradient (Eh/bohr). Comment lines start with `#`; everything else is one value per line.
fn parse_engrad(
    path: &std::path::Path,
    expected_n: usize,
) -> Result<(f64, Vec<Vec3>), ProviderError> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| ProviderError::ComputationFailed(format!("failed to read .engrad: {e}")))?;
    let mut values = content
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));

    let parse_next = |values: &mut dyn Iterator<Item = &str>,
                      what: &str|
     -> Result<f64, ProviderError> {
        values
            .next()
            .ok_or_else(|| ProviderError::ComputationFailed(format!(".engrad missing {what}")))?
            .parse::<f64>()
            .map_err(|e| ProviderError::ComputationFailed(format!("failed to parse {what}: {e}")))
    };

    let n = parse_next(&mut values, "atom count")? as usize;
    if n != expected_n {
        return Err(ProviderError::ComputationFailed(format!(
            ".engrad atom count {n} does not match region size {expected_n}"
        )));
    }
    let energy = parse_next(&mut values, "energy")?;

    let mut gradient = Vec::with_capacity(n);
    for _ in 0..n {
        let gx = parse_next(&mut values, "gradient component")?;
        let gy = parse_next(&mut values, "gradient component")?;
        let gz = parse_next(&mut values, "gradient component")?;
        gradient.push(Vec3::new(gx, gy, gz));
    }

    Ok((energy, gradient))
}

impl PotentialProvider for XtbInteraction {
    fn contribute(
        &mut self,
        region: &AtomSelection,
        topo: &Topology,
        conf: &Configuration,
        neigh: &dyn SpatialIndex,
    ) -> Result<Contribution, ProviderError> {
        match self.embedding {
            Embedding::None | Embedding::Electrostatic => {},
            Embedding::Mechanical => {
                return Err(ProviderError::ComputationFailed(
                    "Mechanical embedding is not implemented for XtbInteraction — nothing in \
                     this crate implements it yet"
                        .to_string(),
                ));
            },
        }
        if self.embedding == Embedding::Electrostatic && self.cutoff <= 0.0 {
            return Err(ProviderError::ComputationFailed(
                "Electrostatic embedding requires a positive cutoff — call with_cutoff() first"
                    .to_string(),
            ));
        }

        let atom_indices: Vec<usize> = region.indices().to_vec();
        let n = atom_indices.len();
        if n == 0 {
            return Err(ProviderError::InvalidRegion("empty region".to_string()));
        }
        if atom_indices.iter().any(|&g| g >= self.atomic_numbers.len()) {
            return Err(ProviderError::InvalidRegion(
                "atomic_numbers table shorter than region's atom indices".to_string(),
            ));
        }

        let xyz_path = self.work_dir.join("region.xyz");
        let engrad_path = self.work_dir.join("region.engrad");
        let pcharge_path = self.work_dir.join("pcharge");
        let pcgrad_path = self.work_dir.join("pcgrad");
        // Remove any stale output before running, so a crashed xtb can't leave a prior result
        // silently readable as this call's answer. `pcharge` is removed unconditionally too —
        // xtb reads it automatically if present, so a leftover file from a prior Electrostatic
        // call must never silently leak into a later None/Mechanical (or empty-MM-list) call.
        qm_subprocess::remove_stale(&engrad_path);
        qm_subprocess::remove_stale(&pcgrad_path);
        qm_subprocess::remove_stale(&pcharge_path);
        self.write_xyz(&xyz_path, &atom_indices, conf);

        let mm_atoms = if self.embedding == Embedding::Electrostatic {
            let mm_atoms = self.gather_mm_atoms(region, &atom_indices, topo, neigh);
            // Matches gromosXX's own guard (`xtb_worker.cc`): an empty external-charge list is
            // skipped rather than written, since xtb errors on a zero-charge pcharge file.
            if !mm_atoms.is_empty() {
                self.write_pcharge(&pcharge_path, &mm_atoms, topo, conf);
            }
            mm_atoms
        } else {
            Vec::new()
        };

        let gfn_arg = self.gfn.to_string();
        let chrg_arg = self.charge.to_string();
        let uhf_arg = (self.multiplicity - 1).to_string();
        qm_subprocess::run_subprocess_with_timeout(
            &self.binary,
            &self.work_dir,
            &[
                "region.xyz",
                "--gfn",
                &gfn_arg,
                "--grad",
                "--chrg",
                &chrg_arg,
                "--uhf",
                &uhf_arg,
            ],
            self.timeout,
        )?;

        let (energy_eh, gradient_eh_bohr) = parse_engrad(&engrad_path, n)?;

        let energy = energy_eh * HARTREE_TO_KJMOL;
        let force_conversion = HARTREE_TO_KJMOL / BOHR_TO_NM;
        // Force = -gradient (xtb, like every QM engine, reports dE/dx, not the force).
        let mut forces: Vec<(usize, Vec3)> = atom_indices
            .iter()
            .zip(gradient_eh_bohr.iter())
            .map(|(&global, &g)| (global, -g * force_conversion))
            .collect();

        if !mm_atoms.is_empty() {
            let mm_gradients = parse_pcgrad(&pcgrad_path, mm_atoms.len())?;
            forces.extend(
                mm_atoms
                    .iter()
                    .zip(mm_gradients.iter())
                    .map(|(&global, &g)| (global, -g * force_conversion)),
            );
        }

        Ok(Contribution {
            energy,
            forces,
            virial: Mat3::ZERO, // not computed yet — no real-system caller needs it yet.
            extra: ProviderExtra::default(),
        })
    }

    fn name(&self) -> &str {
        "xtb (GFN-xTB, subprocess)"
    }

    fn embedding(&self) -> Embedding {
        self.embedding
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::{Periodicity, Vacuum};
    use newen_core::spatial_index::ConfigurationSpatialIndex;

    fn xtb_available() -> bool {
        Command::new("xtb")
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }

    fn water_positions_nm() -> Vec<Vec3> {
        // Same geometry (in Å, converted to nm) verified against a real xtb 6.7.1 GFN2 run on
        // this machine: total energy -5.018180941704 Eh.
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0758602, 0.0, 0.0504284),
            Vec3::new(0.0758602, 0.0, -0.0504284),
        ]
    }

    fn water_topo_and_conf() -> (Topology, Configuration) {
        let topo = Topology::new();
        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().pos = water_positions_nm();
        (topo, conf)
    }

    #[test]
    fn water_monomer_energy_matches_pinned_xtb_oracle() {
        if !xtb_available() {
            eprintln!("skipping: xtb not found on PATH");
            return;
        }
        let work_dir = std::env::temp_dir().join("newen_xtb_test_water_energy");
        let mut interaction =
            XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1]).expect("work_dir creatable");

        let (topo, conf) = water_topo_and_conf();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let c = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("xtb calculation should succeed");

        let expected_energy = -5.018180941704 * HARTREE_TO_KJMOL;
        assert!(
            (c.energy - expected_energy).abs() < 1e-6,
            "energy {} vs pinned oracle {}",
            c.energy,
            expected_energy
        );
        assert_eq!(c.forces.len(), 3);
    }

    /// `Mechanical` must fail loudly, not be silently ignored — same guard `schnet.rs` has for
    /// the same reason. (`Electrostatic` is real now — see the tests below.)
    #[test]
    fn mechanical_embedding_is_rejected_not_silently_ignored() {
        let work_dir = std::env::temp_dir().join("newen_xtb_test_embedding_rejection");
        let mut interaction = XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1])
            .expect("work_dir creatable")
            .with_embedding(Embedding::Mechanical);
        assert_eq!(
            PotentialProvider::embedding(&interaction),
            Embedding::Mechanical
        );

        let (topo, conf) = water_topo_and_conf();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let err = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect_err("mechanical embedding is not implemented and must not silently pass");
        assert!(
            format!("{err}").contains("not implemented"),
            "error should say the scheme is unimplemented, got: {err}"
        );
    }

    /// `Electrostatic` without a cutoff must fail loudly rather than silently gather zero MM
    /// atoms (which would look identical to `Embedding::None` — a much harder bug to notice).
    #[test]
    fn electrostatic_embedding_without_a_cutoff_is_rejected() {
        let work_dir = std::env::temp_dir().join("newen_xtb_test_embedding_no_cutoff");
        let mut interaction = XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1])
            .expect("work_dir creatable")
            .with_embedding(Embedding::Electrostatic);

        let (topo, conf) = water_topo_and_conf();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let err = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect_err("missing cutoff must be a loud error, not a silent zero-MM-atoms run");
        assert!(format!("{err}").contains("cutoff"), "got: {err}");
    }

    /// The tier that actually validates *this provider's* sign/unit-conversion pipeline
    /// (parsing, Hartree->kJ/mol, Bohr->nm, force = -gradient), independent of xtb's own physics
    /// — same rationale as `schnet.rs::model_forces_match_finite_differences`.
    #[test]
    fn forces_match_finite_differences() {
        if !xtb_available() {
            eprintln!("skipping: xtb not found on PATH");
            return;
        }
        let work_dir = std::env::temp_dir().join("newen_xtb_test_fd");
        let mut interaction =
            XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1]).expect("work_dir creatable");

        let positions = water_positions_nm();
        let (topo, _) = water_topo_and_conf();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(Vacuum);

        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().pos = positions.clone();
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let contribution = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("xtb calculation should succeed");

        let mut energy_at = |positions: &[Vec3]| -> f64 {
            let mut conf = Configuration::new(3, 1, 1);
            conf.current_mut().pos = positions.to_vec();
            let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
            interaction
                .contribute(&region, &topo, &conf, &index)
                .unwrap()
                .energy
        };

        let h = 1e-4;
        for atom in 0..3 {
            for axis in 0..3 {
                let mut plus = positions.clone();
                let mut minus = positions.clone();
                let delta = match axis {
                    0 => Vec3::new(h, 0.0, 0.0),
                    1 => Vec3::new(0.0, h, 0.0),
                    _ => Vec3::new(0.0, 0.0, h),
                };
                plus[atom] += delta;
                minus[atom] -= delta;

                let finite_diff_force = -(energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
                let (_, model_force) = contribution.forces[atom];
                let model_component = match axis {
                    0 => model_force.x,
                    1 => model_force.y,
                    _ => model_force.z,
                };

                assert!(
                    (finite_diff_force - model_component).abs() < 5.0,
                    "atom {atom} axis {axis}: finite-diff {finite_diff_force} vs xtb {model_component}"
                );
            }
        }
    }

    /// A water monomer (QM zone, atoms 0-2) plus one external MM point charge (atom 3) at
    /// `mm_pos`, carrying `mm_charge`. Real topology fields populated the way `contribute()`
    /// (via `gather_mm_atoms`) actually reads them: `topo.charge`/`topo.num_atoms()`.
    fn water_plus_mm_atom_topo_and_conf(mm_charge: f64, mm_pos: Vec3) -> (Topology, Configuration) {
        use newen_core::topology::MolTypeAtom;
        let masses = [15.9994, 1.008, 1.008, 35.45];
        let mut topo = Topology::new();
        topo.charge = vec![0.0, 0.0, 0.0, mm_charge];
        topo.iac = vec![0; 4];
        topo.mass = masses.to_vec();
        topo.inverse_mass = masses.iter().map(|m| 1.0 / m).collect();
        topo.exclusions = vec![Vec::new(); 4];
        topo.moltypes[0].atoms = (0..4)
            .map(|i| MolTypeAtom {
                name: format!("A{i}"),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 0,
                mass: masses[i],
                charge: topo.charge[i],
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            })
            .collect();

        let mut positions = water_positions_nm();
        positions.push(mm_pos);
        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos = positions;
        (topo, conf)
    }

    /// **Real electronic polarization, not a formality.** Path (c) (`ElectrostaticEmbedding`)
    /// can only ever add a fixed classical Coulomb term to a fixed QM energy — by construction,
    /// it cannot change what the QM SCF itself converges to. This is the one check that path (a)
    /// actually does something structurally different: the MM charge is fed into xtb's own
    /// Hamiltonian, so the QM energy itself must shift, and the MM atom must receive a real force
    /// back (`pcgrad`) — not just accumulate a bystander Coulomb term computed by our own code.
    #[test]
    fn electrostatic_embedding_changes_qm_energy_and_puts_force_on_the_mm_atom() {
        if !xtb_available() {
            eprintln!("skipping: xtb not found on PATH");
            return;
        }
        let mm_pos = Vec3::new(0.3, 0.0, 0.0);
        let (topo, conf) = water_plus_mm_atom_topo_and_conf(-0.8, mm_pos);
        let region = AtomSelection::from_indices(vec![0, 1, 2], 4).unwrap();
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let none_work_dir = std::env::temp_dir().join("newen_xtb_test_embedding_none_cmp");
        let mut none_interaction =
            XtbInteraction::new(none_work_dir, 2, 0, 1, vec![8, 1, 1]).expect("work_dir creatable");
        let baseline = none_interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("baseline xtb calculation should succeed");

        let elst_work_dir = std::env::temp_dir().join("newen_xtb_test_embedding_elst_cmp");
        let mut elst_interaction = XtbInteraction::new(elst_work_dir, 2, 0, 1, vec![8, 1, 1])
            .expect("work_dir creatable")
            .with_embedding(Embedding::Electrostatic)
            .with_cutoff(1.0);
        let embedded = elst_interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("embedded xtb calculation should succeed");

        assert!(
            (embedded.energy - baseline.energy).abs() > 1.0,
            "an external -0.8e charge 0.3nm away should measurably change the QM energy via real \
             SCF polarization: baseline={} embedded={}",
            baseline.energy,
            embedded.energy
        );

        let mm_force = embedded.forces.iter().find(|&&(i, _)| i == 3);
        assert!(
            mm_force.is_some(),
            "the MM atom (index 3, outside XtbInteraction's own region) must receive a real \
             force via pcgrad — this is the actual QM->MM coupling"
        );
    }

    /// The tier that validates *this embedding path's* sign/unit-conversion pipeline
    /// (`pcharge` in Bohr, `pcgrad` sign/units), independent of xtb's own physics — same
    /// rationale as `forces_match_finite_differences`, extended to cover the MM atom's force
    /// too, which only exists under `Electrostatic`.
    #[test]
    fn electrostatic_embedding_mm_and_qm_forces_match_finite_differences() {
        if !xtb_available() {
            eprintln!("skipping: xtb not found on PATH");
            return;
        }
        let mm_pos = Vec3::new(0.3, 0.0, 0.0);
        let (topo, _) = water_plus_mm_atom_topo_and_conf(-0.8, mm_pos);
        let region = AtomSelection::from_indices(vec![0, 1, 2], 4).unwrap();
        let periodicity = Periodicity::Vacuum(Vacuum);

        let work_dir = std::env::temp_dir().join("newen_xtb_test_embedding_fd");
        let mut interaction = XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1])
            .expect("work_dir creatable")
            .with_embedding(Embedding::Electrostatic)
            .with_cutoff(1.0);

        let mut positions = water_positions_nm();
        positions.push(mm_pos);

        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos = positions.clone();
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let contribution = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("xtb calculation should succeed");

        let mut energy_at = |positions: &[Vec3]| -> f64 {
            let mut conf = Configuration::new(4, 1, 1);
            conf.current_mut().pos = positions.to_vec();
            let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
            interaction
                .contribute(&region, &topo, &conf, &index)
                .unwrap()
                .energy
        };

        let h = 1e-4;
        for atom in 0..4 {
            // includes the MM atom (index 3)
            for axis in 0..3 {
                let mut plus = positions.clone();
                let mut minus = positions.clone();
                let delta = match axis {
                    0 => Vec3::new(h, 0.0, 0.0),
                    1 => Vec3::new(0.0, h, 0.0),
                    _ => Vec3::new(0.0, 0.0, h),
                };
                plus[atom] += delta;
                minus[atom] -= delta;

                let finite_diff_force = -(energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
                let model_force = contribution
                    .forces
                    .iter()
                    .find(|&&(i, _)| i == atom)
                    .map(|&(_, f)| f)
                    .unwrap_or(Vec3::ZERO);
                let model_component = match axis {
                    0 => model_force.x,
                    1 => model_force.y,
                    _ => model_force.z,
                };

                assert!(
                    (finite_diff_force - model_component).abs() < 5.0,
                    "atom {atom} axis {axis}: finite-diff {finite_diff_force} vs xtb {model_component}"
                );
            }
        }
    }
}
