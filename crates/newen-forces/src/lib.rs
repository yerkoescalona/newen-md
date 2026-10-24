//! newen Forces: All force field calculations
//!
//! This crate provides force calculation routines:
//! - Bonded forces (bonds, angles, dihedrals, impropers)
//! - Nonbonded forces (Lennard-Jones, Coulomb)
//! - Electrostatics (PME, reaction field)
//! - Restraints (position, distance, dihedral, NMR)
//! - Enhanced sampling (local elevation, QM/MM, polarization)
//!
//! All calculations use f64 for molecular dynamics precision.
//!
//! The nonbonded kernels take their inputs as separate slices (positions, charges, types,
//! parameters, pairlist, storage, …) — the flat signature is the point of a kernel, and it is
//! what gromosXX's inner loops take too. Clippy's argument-count and index-loop lints are
//! therefore off for this crate (`Cargo.toml` `[lints.clippy]` for the loop lint).
#![allow(clippy::too_many_arguments, clippy::chunks_exact_to_as_chunks)]

pub mod bonded;
pub mod energy;
pub mod nonbonded;
pub mod orchestrator;
pub mod orchestrator_algorithm;
pub mod provider;
pub mod restraints;
pub mod zones;

// Advanced force modules
pub mod local_elevation;
pub mod polarization;
pub mod qmmm;

// GPU acceleration (feature-gated)
#[cfg(feature = "gpu")]
pub mod gpu;

// Re-export main types
pub use bonded::{
    calculate_angle_forces, calculate_bonded_forces, calculate_dihedral_forces,
    calculate_improper_dihedral_forces, ForceEnergy, ForceEnergyLambda,
};
pub use nonbonded::{lj_crf_interaction, CRFParameters, ForceStorage, LJParamMatrix, LJParameters};
pub use orchestrator::ProviderOrchestrator;
pub use provider::{Contribution, Embedding, PotentialProvider, ProviderError, ProviderExtra};
pub use zones::{PairOwner, Zone, ZonePartition};
