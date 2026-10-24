//! newen Integrators: Time integration algorithms
//!
//! This crate provides integration algorithms:
//! - Velocity Verlet
//! - Leap-frog
//! - Stochastic dynamics (Langevin)
//! - Steepest descent minimization
//! - Constraints (SHAKE, LINCS)
//! - Thermostats and barostats
//! - Virtual atoms
//! - Enhanced sampling: REMD, EDS, GAMD, FEP

pub mod algorithms;
pub mod barostats;
pub mod constraints;
pub mod integrator;
pub mod thermostats;
pub mod virtual_atoms;

// Enhanced sampling modules
pub mod eds;
pub mod fep;
pub mod gamd;
pub mod remd;
pub mod replica;

// MPI support (optional)

// Re-export main types
pub use barostats::{berendsen_barostat, calculate_virial, parrinello_rahman_barostat};
pub use barostats::{BerendsenBarostatParameters, ParrinelloRahmanBarostatParameters};
pub use constraints::{
    lincs, settle, shake, shake_buffered, shake_positions, shake_velocities, ConstraintResult,
    LincsParameters, NtcMode, ShakeBuffers, ShakeParameters,
};
pub use integrator::{Integrator, LeapFrog, SteepestDescent, StochasticDynamics, VelocityVerlet};
pub use thermostats::{andersen_thermostat, berendsen_thermostat, nose_hoover_thermostat};
pub use thermostats::{
    AndersenThermostatParameters, BerendsenThermostatParameters, NoseHooverThermostatParameters,
};
pub use virtual_atoms::{VirtualAtom, VirtualAtomsManager};
