//! GROMOS-RS Analysis: Trajectory analysis tools
//!
//! This crate provides analysis functions:
//! - RDF (radial distribution function)
//! - RMSD, RMSF
//! - Hydrogen bonds
//! - Gyration radius
//! - Diffusion coefficients
//! - Clustering
//! - And more...

// Re-export gromos-core types for convenience within this crate
pub use gromos_core::{configuration, math, topology, Configuration, Mat3, State, Topology, Vec3};

/// gromos++-style `Arguments` (lives in gromos-io so the tools share it).
pub use gromos_io::args;
/// `@pbc` parsing and gathering (lives in gromos-io so the tools share it).
pub use gromos_io::pbc;
pub mod diffusion;
pub mod distribution;
pub mod eds;
pub mod fit;
pub mod gyration;
pub mod hbond;
pub mod lnexp;
pub mod pb;
pub mod property;
pub mod rdf;
pub mod rmsd;
pub mod time;

pub mod programs;
pub use programs::{find, PROGRAMS};

// Re-export main functions
pub use diffusion::{calculate_diffusion_coefficient, calculate_msd};
pub use gyration::calculate_gyration_radius;
pub use rdf::calculate_rdf;
pub use rmsd::calculate_rmsd;

/// The Boltzmann constant gromos++ uses when no topology supplies one (`gmath/Physics.cc`):
/// the free-energy programs (`dg_ener`, `dfmult`) reproduce gromos++ with this value.
pub const GROMOSPP_BOLTZMANN: f64 = 0.00831451;
