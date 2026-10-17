//! newen Analysis: Trajectory analysis tools
//!
//! This crate provides analysis functions:
//! - RDF (radial distribution function, gromos++ `rdf`)
//! - RMSD, RMSF
//! - Hydrogen bonds
//! - Gyration radius
//! - Diffusion coefficients
//! - Clustering
//! - And more...

// Re-export newen-core types for convenience within this crate
pub use newen_core::{configuration, math, topology, Configuration, Mat3, State, Topology, Vec3};

/// gromos++-style `Arguments` (lives in newen-io so the tools share it).
pub use newen_io::args;
/// `@pbc` parsing and gathering (lives in newen-io so the tools share it).
pub use newen_io::pbc;
pub mod diffusion;
pub mod distribution;
pub mod eds;
pub mod energy_sequence;
pub mod fit;
pub mod free_energy;
pub mod gyration;
pub mod hbond;
pub mod lnexp;
pub mod pb;
pub mod property;
pub mod rdf;
pub mod rmsd;
pub mod ti;
pub mod time;

pub mod programs;
pub use programs::{find, PROGRAMS};

// Re-export main functions
pub use diffusion::{calculate_diffusion_coefficient, calculate_msd};
pub use gyration::calculate_gyration_radius;
pub use rdf::Rdf;
pub use rmsd::calculate_rmsd;

/// The Boltzmann constant gromos++ uses when no topology supplies one (`gmath/Physics.cc`):
/// the free-energy programs (`dg_ener`, `dfmult`) reproduce gromos++ with this value.
pub const GROMOSPP_BOLTZMANN: f64 = 0.00831451;

/// gromos++'s Boltzmann constant in double precision, evaluated as `gmath/Physics.cc:51-52` does
/// (`1.380658e-23 · 6.0221367e23 / 1000`, 0.008314511211948601): `ext_ti_ana`, `bar` and `mbar`
/// reproduce gromos++ to the last bit with it. [`GROMOSPP_BOLTZMANN`] is its 6-digit print.
pub const GROMOSPP_KB: f64 = 1.380658e-23 * 6.0221367e23 / 1000.0;

/// [`GROMOSPP_KB`] with gromos++'s once-per-process warning that no topology supplied it
/// (`PhysConst::get_boltzmann`, `Physics.cc:137-142, 244-250`).
pub fn gromospp_kb_warn() -> f64 {
    static WARNED: std::sync::Once = std::sync::Once::new();
    WARNED.call_once(|| {
        eprintln!("WARNING: No topology given or no value for boltzmann in topology file.");
        eprintln!("         Setting boltzmann to hardcoded value from gmath/Physics.cc:");
        eprintln!(
            "         boltzmann = {}",
            newen_io::table::cpp_g(GROMOSPP_KB, 6)
        );
    });
    GROMOSPP_KB
}
