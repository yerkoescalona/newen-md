//! Molecular dynamics algorithms
//!
//! This module contains various MD algorithms:
//! - Constraints (SHAKE, SETTLE, M-SHAKE)
//! - Thermostats (Berendsen, Nosé-Hoover, Andersen)
//! - Barostats (Berendsen, Parrinello-Rahman)
//! - Virtual Atoms (for coarse-grained and united-atom models)

pub mod barostats;
pub mod constraints;
pub mod thermostats;
pub mod virtual_atoms;

// Re-export commonly used items
pub use constraints::{m_shake, settle, shake, shake_buffered, ConstraintResult, ShakeBuffers, ShakeParameters};

pub use thermostats::{
    andersen_thermostat, berendsen_thermostat, nose_hoover_thermostat,
    AndersenThermostatParameters, BerendsenThermostatParameters, NoseHooverThermostatParameters,
};

pub use barostats::{
    berendsen_barostat, calculate_virial, parrinello_rahman_barostat, BerendsenBarostatParameters,
    ParrinelloRahmanBarostatParameters,
};
