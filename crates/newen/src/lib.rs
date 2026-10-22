//! newen: High-performance Rust molecular dynamics engine
//!
//! This is the facade crate that re-exports everything from the sub-crates.
//! For most users, this is the only crate you need to depend on.
//!
//! ## Features
//!
//! - **Pure Rust MD engine**: Complete translation, not just FFI wrappers
//! - **SIMD vectorization**: Automatic vectorization with glam + wide
//! - **Fearless concurrency**: Parallel algorithms with Rayon
//! - **Memory safety**: Zero-cost abstractions without runtime overhead
//! - **Modular design**: Easy to extend and customize
//!
//! ## Architecture
//!
//! ```text
//! newen
//! ├── core          - Math types, topology, configuration
//! ├── forces        - Force calculations (bonded + nonbonded)
//! ├── integrators   - Time integration algorithms
//! ├── io            - File I/O
//! └── analysis      - Trajectory analysis
//! ```

// Re-export sub-crates as modules (Polars style)
pub use newen_core as core;
pub use newen_forces as forces;
pub use newen_integrators as integrators;
pub use newen_io as io;
// NOTE: newen-analysis is NOT re-exported here to avoid circular dependencies.
// Use `newen_analysis` directly when needed.

// Also re-export with common names for compatibility with existing code
pub mod math {
    pub use newen_core::math::*;
    pub use newen_core::{Mat3, Vec3};
}

pub mod topology {
    pub use newen_core::topology::*;
    pub use newen_core::Topology;
}

pub mod configuration {
    pub use newen_core::configuration::*;
    pub use newen_core::{Configuration, Energy, State};
}

pub mod selection {
    pub use newen_core::selection::*;
}

pub mod random {
    pub use newen_core::random::*;
}

pub mod pairlist {
    pub use newen_core::pairlist::*;
}

// Compatibility aliases for old module names
pub mod algorithm {
    pub use newen_core::algorithm::{Algorithm, AlgorithmSequence, SimulationState};
    pub use newen_integrators::algorithms;
    pub use newen_integrators::algorithms::*;
    pub use newen_integrators::barostats;
    pub use newen_integrators::barostats::*;
    pub use newen_integrators::constraints;
    pub use newen_integrators::constraints::*;
    pub use newen_integrators::thermostats;
    pub use newen_integrators::thermostats::*;
    pub use newen_integrators::virtual_atoms;
    pub use newen_integrators::virtual_atoms::*;
}

pub mod interaction {
    // Both bonded::* and nonbonded::* expose a `perturbed` submodule; suppress
    // the re-export ambiguity warning — callers use the namespaced paths.
    pub use newen_forces::bonded;
    #[allow(ambiguous_glob_reexports)]
    pub use newen_forces::bonded::*;
    pub use newen_forces::nonbonded;
    #[allow(ambiguous_glob_reexports)]
    pub use newen_forces::nonbonded::*;
    pub use newen_forces::restraints;
}

pub mod integrator {
    pub use newen_integrators::integrator::*;
    pub use newen_integrators::{
        Integrator, LeapFrog, SteepestDescent, StochasticDynamics, VelocityVerlet,
    };
}

// Enhanced sampling re-exports
pub use newen_integrators::eds;
pub use newen_integrators::fep;
pub use newen_integrators::gamd;
pub use newen_integrators::remd;
pub use newen_integrators::replica;

// Validation
pub use newen_core::validation;

// Re-export main types at crate root for convenience
pub use newen_core::{
    BoundaryCondition, Configuration, Energy, LJParameters, Mat3, Periodicity, State, Topology,
    Vec3,
};

pub use newen_integrators::{
    Integrator, LeapFrog, SteepestDescent, StochasticDynamics, VelocityVerlet,
};

#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
