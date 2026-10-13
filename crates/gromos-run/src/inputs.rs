//! Auxiliary run inputs and execution options.
//!
//! An `.imd` file describes a run only together with its auxiliary files (3.9 A5);
//! [`RunInputs`] carries those. [`RunOptions`] carries what is *not* physics and therefore has
//! no IMD field (A15, `execution` group of the future recipe): the parallel policy.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

/// Paths of the optional GROMOS input files — the `md` binary's `@pttopo`, `@posresspec`,
/// `@refpos` and `@distrest` flags.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct RunInputs {
    /// Perturbation topology (`.ptp`); only read when `NTG != 0`.
    pub pttopo: Option<PathBuf>,
    /// Position-restraint atom specification (`POSRESSPEC`); required when `NTPOR > 0`.
    pub posresspec: Option<PathBuf>,
    /// Reference positions for position restraints; required when `NTPORB >= 1`.
    pub refpos: Option<PathBuf>,
    /// Distance-restraint specification; required when `NTDIR != 0`.
    pub distrest: Option<PathBuf>,
}

/// Whether the nonbonded kernels run on the rayon pool.
///
/// A platform concern, not physics: the parallel kernels are reference-verified at the suite
/// tolerance (1e-8) and bit-identical run-to-run at a fixed thread count, but differ from the
/// serial kernels in the last digit (BENCHMARKING.md, "Kernel determinism").
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ParallelPolicy {
    /// What the `md` binary has always done: parallel above 100 atoms.
    #[default]
    Auto,
    /// Never parallel — the setting for exact (`==`) cross-front-end comparisons.
    Serial,
    /// Always parallel.
    Parallel,
}

impl ParallelPolicy {
    /// Resolve the policy for a system of `n_atoms` atoms.
    pub fn resolve(self, n_atoms: usize) -> bool {
        match self {
            ParallelPolicy::Auto => n_atoms > 100,
            ParallelPolicy::Serial => false,
            ParallelPolicy::Parallel => true,
        }
    }
}

/// Execution options that are not part of the physics description.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct RunOptions {
    pub parallel: ParallelPolicy,
    /// MPI pair decomposition `(rank, size)` for the force field (`md` under `mpirun`); the
    /// caller installs the reduction on the built `Forcefield`. `None`: single process.
    pub pair_partition: Option<(usize, usize)>,
    /// Unmodelled `.imd` blocks this caller accepts as passthrough (3.9 A17). The
    /// binary allows `GAMD`/`EDS` because it applies them itself; the Python binding allows none.
    pub passthrough: crate::recipe::PassthroughPolicy,
    /// Recipe fields no `.imd` block carries, set by dotted path after `from_imd`
    /// ([`RunRecipe::set`](crate::RunRecipe::set); `md @recipe_set key=value`).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub recipe_overrides: Vec<(String, String)>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parallel_policy_resolves_as_documented() {
        assert!(!ParallelPolicy::Auto.resolve(100));
        assert!(ParallelPolicy::Auto.resolve(101));
        assert!(!ParallelPolicy::Serial.resolve(1_000_000));
        assert!(ParallelPolicy::Parallel.resolve(1));
    }

    #[test]
    fn run_options_serialise_with_snake_case_policies() {
        let o = RunOptions {
            parallel: ParallelPolicy::Serial,
            ..Default::default()
        };
        let text = toml::to_string(&o).unwrap();
        assert!(text.contains("parallel = \"serial\""), "{text}");
        let back: RunOptions = toml::from_str(&text).unwrap();
        assert_eq!(back, o);
    }
}
