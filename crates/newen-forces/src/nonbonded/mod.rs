//! Nonbonded interactions: Lennard-Jones and Coulomb (CRF)

pub mod eds;
pub mod embedding;
pub mod ewald;
pub mod excluded;
pub mod innerloops;
pub mod interaction;
pub mod lattice_sum;
pub mod ls_real_space;
pub mod mopac;
pub mod p3m;
pub mod params;
pub mod perturbed;
mod qm_subprocess;
#[cfg(feature = "ml")]
pub mod schnet;
pub mod spme;
pub mod xtb;

pub use eds::{
    eds_one_four_interactions, eds_pair_interactions, eds_rf_excluded_interactions,
    eds_ti_one_four_interactions, eds_ti_pair_interactions, eds_ti_rf_excluded_interactions,
    split_eds_pairs, EdsAtoms, EdsStorage,
};
pub use embedding::ElectrostaticEmbedding;
pub use excluded::*;
pub use innerloops::*;
pub use interaction::LjCrfInteraction;
pub use mopac::{MopacInteraction, MopacMethod};
pub use params::*;
pub use perturbed::*;
#[cfg(feature = "ml")]
pub use schnet::SchNetInteraction;
pub use xtb::XtbInteraction;

/// MPI-parallel nonbonded force calculation
///
/// # Status
/// Stub — pending a dedicated `newen-mpi` crate that holds `MpiControl` /
/// `MpiNonbondedMaster`. `newen-forces` cannot import from `newen-integrators`
/// (would be circular). Wire up when the MPI support crate is extracted.
///
/// Ported from `md++/src/interaction/nonbonded/interaction/mpi_nonbonded_master.cc`
#[cfg(feature = "use-mpi")]
pub mod mpi_nonbonded {
    use super::*;

    /// Calculate nonbonded interactions with MPI parallelisation (stub).
    ///
    /// Signature mirrors `MPI_Nonbonded_Master::calculate_interactions`.
    /// Returns a zeroed `ForceStorage` until the MPI layer is wired.
    pub fn calculate_mpi<BC: newen_core::math::BoundaryCondition>(
        _positions: &[newen_core::math::Vec3],
        _charges: &[f64],
        _iac: &[u32],
        _pairlist: &[(u32, u32)],
        _lj_params: &LJParamMatrix,
        _crf: &CRFParameters,
        _periodicity: &BC,
    ) -> ForceStorage {
        // TODO: rewire — needs newen-mpi crate with MpiControl + MpiNonbondedMaster
        ForceStorage::new(_positions.len())
    }
}
