#![warn(missing_docs)]
//! GROMOS-RS Core: Math types, topology, and fundamental data structures
//!
//! This crate provides the foundational types used throughout GROMOS-RS:
//! - Floating-point precision trait (`Float`)
//! - Generic vector and matrix types (`Vector3<F>`, `Matrix3<F>`)
//! - Boundary conditions and periodicity
//! - Atom and topology definitions
//! - System configuration and state
//! - Pairlist generation
//! - Selection language
//!
//! # Precision
//!
//! All types are generic over the floating-point precision. Use `f64` for
//! CPU-based MD simulations (better energy conservation) or `f32` for
//! GPU acceleration and reduced memory bandwidth.
//!
//! ```rust
//! use gromos_core::{Vector3, Vec3f, Vec3d};
//!
//! // Type aliases for convenience
//! let pos_f32: Vec3f = Vector3::new(1.0, 2.0, 3.0);
//! let pos_f64: Vec3d = Vector3::new(1.0, 2.0, 3.0);
//!
//! // Convert between precisions
//! let converted = pos_f32.to_f64();
//! ```

pub mod algorithm;
pub mod configuration;
pub mod float;
pub mod gather;
pub mod lambdas;
pub mod math;
pub mod matrix;
pub mod pairlist;
pub mod random;
pub mod selection;
pub mod spatial_index;
pub mod stat;
pub mod topology;
pub mod units;
pub mod validation;
pub mod vector;

// FFI for C bindings (feature-gated)
#[cfg(feature = "ffi")]
pub mod ffi;

// Re-export Float trait and precision types
pub use float::{Float, F32, F64};

// Re-export generic vector and matrix types
pub use matrix::{Mat3d, Mat3f, Matrix3};
pub use vector::{Vec3d, Vec3f, Vector3, VEC3D_ZERO, VEC3F_ZERO};

// Legacy compatibility - keep the old Vec3/Mat3 names pointing to the glam-based types
// These will be phased out in favor of the generic types
pub use math::{BoundaryCondition, Mat3, Periodicity, Rectangular, Triclinic, Vacuum, Vec3};

// Re-export configuration types
pub use configuration::{
    AedsSearch, Box, BoxType, Configuration, EdsEndStates, Energy, PrecalcEnergies, State,
    StochasticVariables,
};

// Re-export the individual-λ table
pub use lambdas::{IndividualLambdas, InteractionLambda, PrecalcLambdas};

// Re-export topology types
pub use topology::{
    EdsPerturbedAtom, EdsSolute, LJParameters, NetCharge, PairInteraction, PerturbedAngle,
    PerturbedAtom, PerturbedAtomPair, PerturbedBond, PerturbedDihedral, PerturbedSolute, SoftAngle,
    SoftBond, SoftImproper, Topology,
};

// Re-export pairlist types
pub use pairlist::{
    CellListPairlistAlgorithm, Pairlist, PairlistAlgorithm, PairlistContainer,
    StandardPairlistAlgorithm,
};

// Re-export selection
pub use selection::AtomSelection;
