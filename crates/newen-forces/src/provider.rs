//! `PotentialProvider` — the additive-provider shape for force/energy computation.
//!
//! Classical LJ+CRF, bonded terms, and (eventually) QM/MM and ML potentials all implement
//! this one trait. If classical terms are *already* providers, QM/MM and ML potentials are
//! just more providers later — additive, touching no existing physics (`architecture.md`
//! "the provider pattern"; the P1–P9 constraints in this crate's CONTEXT.md).
//!
//! This module intentionally does **not** hand a provider `&mut Configuration`: that would
//! give any provider — including a future external QM/ML one — unchecked mutable access to
//! every atom's force/position/velocity, defeating the "contributes to an arbitrary but
//! *scoped* atom subset" invariant (P5). Instead a provider reads state and
//! returns a scattered [`Contribution`]; the caller (which owns `Configuration`/`Energy`)
//! accumulates it and can validate which indices were touched before doing so.

use newen_core::configuration::Configuration;
use newen_core::math::{Mat3, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::SpatialIndex;
use newen_core::topology::Topology;

/// A force/energy term, expressed as an additive contribution over shared state.
///
/// Classical terms and future QM/MM/ML engines are all implementors. The shape — data in,
/// scattered contribution out, stateful, fallible — is what must be right early; see
/// This crate's CONTEXT.md §"The provider constraints P1–P9" for what this is designed
/// against. Deliberately not
/// solved yet (see the plan this landed from): async/cancellable evaluation for slow or
/// external providers, and a shared typed-units boundary (each provider currently owns its
/// own unit conversion).
pub trait PotentialProvider {
    /// Compute this term's contribution for `region`.
    ///
    /// `region` is an arbitrary atom subset — the whole system for classical terms, a QM
    /// zone plus its embedding MM atoms for a QM/MM provider. `&mut self` allows stateful
    /// providers (SCF iteration, wavefunction reuse between steps). Fallible because
    /// external engines can fail (crash, timeout, non-convergence).
    fn contribute(
        &mut self,
        region: &AtomSelection,
        topo: &Topology,
        conf: &Configuration,
        neigh: &dyn SpatialIndex,
    ) -> Result<Contribution, ProviderError>;

    /// Human-readable name for logging/diagnostics (e.g. `"LJ+CRF"`, `"MACE-OFF"`).
    fn name(&self) -> &str;

    /// How this provider treats atoms **outside** `region` (P2.7 Step 2).
    ///
    /// [`SpatialIndex::neighbor_pairs`] deliberately returns pairs with only one endpoint in
    /// `region`, so every provider has to decide what to do with the environment. Declaring it
    /// makes "this provider ignores its surroundings" a stated property rather than an
    /// accident — the P2.6 bug was exactly an unstated assumption here. Defaults to
    /// [`Embedding::None`], which is correct for any provider whose region is the whole system.
    fn embedding(&self) -> Embedding {
        Embedding::None
    }
}

/// How a [`PotentialProvider`] couples to atoms outside its region.
///
/// Grounded in the schemes GROMOS actually implements (Poliak et al. 2025, *J. Comput. Chem.*
/// 46:e70053), not invented: mechanical embedding sees only environment *geometry*, while
/// electrostatic embedding hands the environment's point charges to the engine so they polarize
/// it. See P2.7 for which of these are implemented.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Embedding {
    /// Environment is ignored entirely; cross-boundary neighbor pairs are dropped.
    ///
    /// Correct when the region *is* the whole system, and the honest description of an
    /// isolated-region ML potential that was never given environment inputs.
    #[default]
    None,
    /// Environment geometry participates, but its charges do not polarize the region.
    Mechanical,
    /// Environment point charges polarize the region, and forces flow back onto them
    /// (P5). Note GROMOS's own EE polarizes the QM zone by the MM atoms but
    /// **not** vice versa — mutual polarization is a separate, harder scheme.
    Electrostatic,
}

/// Scattered output of a [`PotentialProvider`]: energy plus sparse per-atom forces.
///
/// `forces` lists only the atoms actually touched — which may be wider than `region` (an
/// electrostatically embedded QM provider also puts forces on nearby MM atoms), but the
/// caller can still check every index against the set it expects before accumulating.
#[derive(Debug, Clone)]
pub struct Contribution {
    /// Energy contributed by this provider (kJ/mol).
    pub energy: f64,
    /// Sparse per-atom force contributions: `(global atom index, force in kJ/mol/nm)`.
    pub forces: Vec<(usize, Vec3)>,
    /// Virial tensor contribution (kJ/mol), for pressure coupling.
    pub virial: Mat3,
    /// Extensible return channel (uncertainty, etc.) — see [`ProviderExtra`].
    pub extra: ProviderExtra,
}

impl Contribution {
    /// A zero contribution — useful as a starting accumulator or for providers with nothing
    /// to add this step.
    pub fn zero() -> Self {
        Self {
            energy: 0.0,
            forces: Vec::new(),
            virial: Mat3::ZERO,
            extra: ProviderExtra::default(),
        }
    }
}

/// Extensible extra return data alongside energy/forces.
///
/// `uncertainty` is a placeholder for the ensemble/active-learning hook (P7:
/// disagreement across a committee of models triggers a DFT recompute or flags the step for
/// active learning). A single `f64` cannot represent a committee's per-model spread — this
/// field is expected to widen (e.g. to `Vec<f64>`) when the first ML provider needing it
/// actually lands; it is not a locked-in shape.
#[derive(Debug, Clone, Default)]
pub struct ProviderExtra {
    /// Optional scalar uncertainty estimate, if this provider can produce one.
    pub uncertainty: Option<f64>,
}

/// Why a [`PotentialProvider::contribute`] call failed.
#[derive(Debug, Clone)]
pub enum ProviderError {
    /// The underlying computation failed (e.g. external engine crashed, non-convergence).
    ComputationFailed(String),
    /// The requested region was invalid for this provider (e.g. empty, or outside its domain).
    InvalidRegion(String),
}

impl std::fmt::Display for ProviderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProviderError::ComputationFailed(msg) => {
                write!(f, "provider computation failed: {msg}")
            },
            ProviderError::InvalidRegion(msg) => write!(f, "invalid region for provider: {msg}"),
        }
    }
}

impl std::error::Error for ProviderError {}
