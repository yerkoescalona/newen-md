//! `MopacInteraction` — a real MOPAC semiempirical QM engine as a [`PotentialProvider`].
//!
//! Same subprocess file-in/file-out shape as [`crate::nonbonded::XtbInteraction`] (shared
//! mechanics in `qm_subprocess.rs`), but a different engine and file format. `method` is a plain
//! field, so PM7, AM1, PM3 and PM6 all come from this one implementation — MOPAC keyword names,
//! not separate ports.
//!
//! **Oracle note, found empirically, not assumed.** The BuRNN tutorial ships 1722 real MOPAC
//! `.aux` files (`MOPAC_VERSION=MOPAC2016.22.067L`, method `PM7`) — but re-running one archived
//! geometry through this machine's installed MOPAC 23.1.2 showed a single-point energy 165
//! kcal/mol (11%) away from the archived value, even though each version's *own* independently
//! optimized minimum agrees to ~1.2% (confirmed by three real runs: single-point at the archived
//! geometry, full re-optimization from it, and single-point at that re-optimized geometry — the
//! last two agree to 5 decimal places, proving MOPAC 23 is internally self-consistent and the
//! mismatch is a real cross-version PM7 surface difference, not a bug here). So those files
//! aren't usable as tight-tolerance single-point oracles. This provider's own tests instead pin
//! real values from *this machine's* installed MOPAC — the same pattern `XtbInteraction` uses
//! (`nonbonded/xtb.rs`) — using the archived geometries as realistic test inputs, not as
//! output-value oracles.
//!
//! **Embedding: [`Embedding::None`] only.** gromosXX's real `mopac_worker.cc` doesn't implement
//! `parse_mm_gradients` (unlike `xtb_worker.cc`) — MOPAC never returns MM-atom forces, so real
//! electrostatic embedding (path (a), see `nonbonded/xtb.rs`) isn't available for this engine.
//! Pairing `MopacInteraction` with QM/MM coupling means path (c)
//! ([`crate::nonbonded::ElectrostaticEmbedding`]).

use newen_core::configuration::Configuration;
use newen_core::math::{Mat3, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::SpatialIndex;
use newen_core::topology::Topology;
use newen_core::units::{KCAL_TO_KJ, NM_TO_ANGSTROM};
use std::io::Write;
use std::path::PathBuf;

use super::qm_subprocess;
use crate::provider::{Contribution, Embedding, PotentialProvider, ProviderError, ProviderExtra};

/// MOPAC semiempirical method — a keyword difference only, not a separate implementation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MopacMethod {
    Pm7,
    Am1,
    Pm3,
    Pm6,
}

impl MopacMethod {
    fn keyword(self) -> &'static str {
        match self {
            MopacMethod::Pm7 => "PM7",
            MopacMethod::Am1 => "AM1",
            MopacMethod::Pm3 => "PM3",
            MopacMethod::Pm6 => "PM6",
        }
    }
}

/// Element symbol for atomic numbers 1-54 (H through Xe). MOPAC's `.mop` format requires a real
/// element symbol — unlike xtb's xyz parser, it does not also accept a bare atomic number.
fn element_symbol(atomic_number: i64) -> Result<&'static str, ProviderError> {
    const SYMBOLS: [&str; 54] = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    ];
    usize::try_from(atomic_number)
        .ok()
        .and_then(|n| n.checked_sub(1))
        .and_then(|i| SYMBOLS.get(i))
        .copied()
        .ok_or_else(|| {
            ProviderError::InvalidRegion(format!(
                "unsupported atomic number {atomic_number} (element_symbol table covers 1-54)"
            ))
        })
}

/// MOPAC spin-state keyword for a multiplicity (1 = singlet, MOPAC's default — no keyword needed).
fn spin_keyword(multiplicity: i32) -> Result<Option<&'static str>, ProviderError> {
    match multiplicity {
        1 => Ok(None),
        2 => Ok(Some("DOUBLET")),
        3 => Ok(Some("TRIPLET")),
        4 => Ok(Some("QUARTET")),
        5 => Ok(Some("QUINTET")),
        6 => Ok(Some("SEXTET")),
        7 => Ok(Some("SEPTET")),
        _ => Err(ProviderError::InvalidRegion(format!(
            "unsupported multiplicity {multiplicity} (MOPAC spin keywords cover singlet-septet)"
        ))),
    }
}

/// Runs a real `mopac` binary as a subprocess to get a QM region's isolated energy and forces.
pub struct MopacInteraction {
    binary: String,
    method: MopacMethod,
    charge: i32,
    multiplicity: i32,
    /// Atomic number per global atom index — the caller's job to map from `Topology`, same
    /// convention as `XtbInteraction::atomic_numbers`.
    atomic_numbers: Vec<i64>,
    work_dir: PathBuf,
    embedding: Embedding,
}

impl MopacInteraction {
    /// `work_dir` is created if missing and reused (overwritten) on every `contribute()` call.
    pub fn new(
        work_dir: impl Into<PathBuf>,
        method: MopacMethod,
        charge: i32,
        multiplicity: i32,
        atomic_numbers: Vec<i64>,
    ) -> Result<Self, ProviderError> {
        let work_dir = qm_subprocess::ensure_work_dir(work_dir)?;
        Ok(Self {
            binary: "mopac".to_string(),
            method,
            charge,
            multiplicity,
            atomic_numbers,
            work_dir,
            embedding: Embedding::None,
        })
    }

    /// Override the `mopac` executable (default: `"mopac"`, resolved via `PATH`).
    pub fn with_binary(mut self, binary: impl Into<String>) -> Self {
        self.binary = binary.into();
        self
    }

    /// Declare a non-default embedding scheme. Accepted here but rejected at `contribute()` time
    /// with a clear error — see module docs for why (MOPAC never returns MM gradients).
    pub fn with_embedding(mut self, embedding: Embedding) -> Self {
        self.embedding = embedding;
        self
    }

    /// Write the `.mop` input: keyword line, title, blank, then per-atom `El x 1 y 1 z 1`. The
    /// `1` flags mark every coordinate as an optimization variable — required for MOPAC's `GRAD`
    /// keyword to report anything at all ("Keyword GRADIENTS used, but geometry has no
    /// variables" otherwise, confirmed against a real run); `1SCF` still prevents any actual
    /// optimization from moving the geometry, matching gromosXX's own tutorial MOPAC input
    /// convention (`train_dataset_tutorial/mopac.py`).
    fn write_mop(
        &self,
        path: &std::path::Path,
        atom_indices: &[usize],
        conf: &Configuration,
    ) -> Result<(), ProviderError> {
        let mut file = std::fs::File::create(path).map_err(|e| {
            ProviderError::ComputationFailed(format!("failed to create .mop input: {e}"))
        })?;
        let io = |e: std::io::Error| {
            ProviderError::ComputationFailed(format!("failed to write .mop input: {e}"))
        };
        let mut keywords = format!(
            "{} 1SCF GRAD AUX(PRECISION=9) CHARGE={}",
            self.method.keyword(),
            self.charge
        );
        if let Some(spin) = spin_keyword(self.multiplicity)? {
            keywords.push(' ');
            keywords.push_str(spin);
        }
        writeln!(file, "{keywords}").map_err(io)?;
        writeln!(file, "newen MopacInteraction").map_err(io)?;
        writeln!(file).map_err(io)?;
        for &g in atom_indices {
            let symbol = element_symbol(self.atomic_numbers[g])?;
            let p = conf.current().pos[g] * NM_TO_ANGSTROM;
            writeln!(
                file,
                "{:>2}{:14.9} 1{:14.9} 1{:14.9} 1",
                symbol, p.x, p.y, p.z
            )
            .map_err(io)?;
        }
        Ok(())
    }
}

/// Parse a single scalar `KEY=value` line (e.g. `HEAT_OF_FORMATION:KCAL/MOL=...`). MOPAC writes
/// Fortran `D` exponents (`-0.161088149248928D+04`), which `f64::from_str` rejects — normalized
/// to `E` first.
fn parse_scalar(content: &str, key: &str) -> Result<f64, ProviderError> {
    for line in content.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix(key) {
            if let Some(val) = rest.strip_prefix('=') {
                let normalized = val.trim().replace(['D', 'd'], "E");
                return normalized.parse::<f64>().map_err(|e| {
                    ProviderError::ComputationFailed(format!("failed to parse {key}: {e}"))
                });
            }
        }
    }
    Err(ProviderError::ComputationFailed(format!(
        "{key} not found in .aux output"
    )))
}

/// Parse a `KEY[count]=` block: the header line gives the value count, then that many
/// whitespace-separated (possibly line-wrapped) floats follow.
fn parse_array(content: &str, key_prefix: &str) -> Result<Vec<f64>, ProviderError> {
    let mut lines = content.lines();
    let mut count = None;
    for line in lines.by_ref() {
        let t = line.trim();
        if t.starts_with(key_prefix) {
            if let (Some(start), Some(end)) = (t.find('['), t.find(']')) {
                count = t[start + 1..end].trim().parse::<usize>().ok();
            }
            break;
        }
    }
    let count = count.ok_or_else(|| {
        ProviderError::ComputationFailed(format!("{key_prefix} block not found in .aux output"))
    })?;

    let mut values = Vec::with_capacity(count);
    for line in lines {
        if values.len() >= count {
            break;
        }
        for tok in line.split_whitespace() {
            if values.len() >= count {
                break;
            }
            let normalized = tok.replace(['D', 'd'], "E");
            let v: f64 = normalized.parse().map_err(|e| {
                ProviderError::ComputationFailed(format!(
                    "failed to parse {key_prefix} value '{tok}': {e}"
                ))
            })?;
            values.push(v);
        }
    }
    if values.len() != count {
        return Err(ProviderError::ComputationFailed(format!(
            "{key_prefix}: expected {count} values, got {}",
            values.len()
        )));
    }
    Ok(values)
}

impl PotentialProvider for MopacInteraction {
    fn contribute(
        &mut self,
        region: &AtomSelection,
        _topo: &Topology,
        conf: &Configuration,
        _neigh: &dyn SpatialIndex,
    ) -> Result<Contribution, ProviderError> {
        match self.embedding {
            Embedding::None => {},
            Embedding::Mechanical | Embedding::Electrostatic => {
                return Err(ProviderError::ComputationFailed(format!(
                    "{:?} embedding is not implemented for MopacInteraction — MOPAC never \
                     returns MM-atom forces (confirmed against gromosXX's real \
                     mopac_worker.cc), so real embedding isn't available for this engine; pair \
                     it with ElectrostaticEmbedding instead",
                    self.embedding
                )));
            },
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

        let mop_path = self.work_dir.join("region.mop");
        let aux_path = self.work_dir.join("region.aux");
        qm_subprocess::remove_stale(&aux_path);
        self.write_mop(&mop_path, &atom_indices, conf)?;

        qm_subprocess::run_subprocess(&self.binary, &self.work_dir, &["region.mop"])?;

        let content = std::fs::read_to_string(&aux_path).map_err(|e| {
            ProviderError::ComputationFailed(format!("failed to read .aux output: {e}"))
        })?;

        let energy_kcal = parse_scalar(&content, "HEAT_OF_FORMATION:KCAL/MOL")?;
        let gradients_kcal_ang = parse_array(&content, "GRADIENTS:KCAL/MOL/ANGSTROM")?;
        if gradients_kcal_ang.len() != 3 * n {
            return Err(ProviderError::ComputationFailed(format!(
                "expected {} gradient components, got {}",
                3 * n,
                gradients_kcal_ang.len()
            )));
        }

        let energy = energy_kcal * KCAL_TO_KJ;
        // kcal/mol/Å -> kJ/mol/nm: KCAL_TO_KJ for the energy unit, NM_TO_ANGSTROM (10.0) because
        // a "per Å" quantity is 10x larger expressed "per nm" (same derivation as XtbInteraction's
        // Hartree/bohr -> kJ/mol/nm conversion).
        let force_conversion = KCAL_TO_KJ * NM_TO_ANGSTROM;
        let forces = atom_indices
            .iter()
            .enumerate()
            .map(|(i, &global)| {
                let g = Vec3::new(
                    gradients_kcal_ang[3 * i],
                    gradients_kcal_ang[3 * i + 1],
                    gradients_kcal_ang[3 * i + 2],
                );
                // Force = -gradient, same convention as XtbInteraction / every QM engine here.
                (global, -g * force_conversion)
            })
            .collect();

        Ok(Contribution {
            energy,
            forces,
            virial: Mat3::ZERO, // not computed yet — same honest gap as XtbInteraction.
            extra: ProviderExtra::default(),
        })
    }

    fn name(&self) -> &str {
        "MOPAC (subprocess)"
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

    fn mopac_available() -> bool {
        std::process::Command::new("mopac")
            .output()
            .map(|o| o.status.success() || o.status.code() == Some(0))
            .unwrap_or(false)
    }

    fn water_positions_nm() -> Vec<Vec3> {
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
    fn water_monomer_energy_matches_pinned_mopac_oracle() {
        if !mopac_available() {
            eprintln!("skipping: mopac not found on PATH");
            return;
        }
        let work_dir = std::env::temp_dir().join("newen_mopac_test_water_energy");
        let mut interaction =
            MopacInteraction::new(work_dir, MopacMethod::Pm7, 0, 1, vec![8, 1, 1])
                .expect("work_dir creatable");

        let (topo, conf) = water_topo_and_conf();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let c = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("mopac calculation should succeed");

        // Pinned against a real run of this machine's own installed MOPAC 23.1.2 (not the
        // archived MOPAC2016 tutorial data — see module docs on why those aren't a valid oracle).
        let expected_energy = -27.6185382947224 * KCAL_TO_KJ;
        assert!(
            (c.energy - expected_energy).abs() < 1e-4,
            "energy {} vs pinned oracle {}",
            c.energy,
            expected_energy
        );
        assert_eq!(c.forces.len(), 3);
    }

    /// `Mechanical`/`Electrostatic` must fail loudly, not be silently ignored — same guard
    /// `xtb.rs`/`schnet.rs` have for the same reason.
    #[test]
    fn unsupported_embedding_is_rejected_not_silently_ignored() {
        let work_dir = std::env::temp_dir().join("newen_mopac_test_embedding_rejection");
        let mut interaction =
            MopacInteraction::new(work_dir, MopacMethod::Pm7, 0, 1, vec![8, 1, 1])
                .expect("work_dir creatable")
                .with_embedding(Embedding::Electrostatic);
        assert_eq!(
            PotentialProvider::embedding(&interaction),
            Embedding::Electrostatic
        );

        let (topo, conf) = water_topo_and_conf();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let err = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect_err("electrostatic embedding is not implemented and must not silently pass");
        assert!(
            format!("{err}").contains("not implemented"),
            "error should say the scheme is unimplemented, got: {err}"
        );
    }

    /// The tier that actually validates *this provider's* sign/unit-conversion pipeline
    /// (D-exponent normalization, kcal/mol/Å -> kJ/mol/nm, force = -gradient), independent of
    /// MOPAC's own physics — same rationale as `xtb.rs::forces_match_finite_differences`.
    #[test]
    fn forces_match_finite_differences() {
        if !mopac_available() {
            eprintln!("skipping: mopac not found on PATH");
            return;
        }
        let work_dir = std::env::temp_dir().join("newen_mopac_test_fd");
        let mut interaction =
            MopacInteraction::new(work_dir, MopacMethod::Pm7, 0, 1, vec![8, 1, 1])
                .expect("work_dir creatable");

        let positions = water_positions_nm();
        let (topo, _) = water_topo_and_conf();
        let region = AtomSelection::all(3);
        let periodicity = Periodicity::Vacuum(Vacuum);

        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().pos = positions.clone();
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let contribution = interaction
            .contribute(&region, &topo, &conf, &index)
            .expect("mopac calculation should succeed");

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
                    "atom {atom} axis {axis}: finite-diff {finite_diff_force} vs mopac {model_component}"
                );
            }
        }
    }
}
