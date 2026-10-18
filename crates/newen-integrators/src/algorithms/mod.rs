//! Algorithm implementations for the MD sequence.
//!
//! Each struct implements `newen_core::algorithm::Algorithm` and can be
//! added to an `AlgorithmSequence` to build the MD step.
//!
//! The standard GROMOS leap-frog MD sequence is:
//! ```text
//!   1. Forcefield (nonbonded + bonded forces)
//!   2. LeapFrogVelocity (exchange_state + velocity update)
//!   3. Thermostat (optional: Berendsen, Nosé-Hoover)
//!   4. LeapFrogPosition (position update)
//!   5. Constraints (optional: SHAKE)
//!   6. TemperatureCalculation (kinetic energy)
//!   7. PressureCalculation (optional)
//!   8. Barostat (optional: Berendsen)
//!   9. EnergyCalculation (total energy)
//! ```

mod berendsen_barostat;
mod berendsen_thermostat;
mod energy;
mod forcefield;
mod lattice_shift;
mod leap_frog;
mod lincs_algorithm;
mod multibath;
mod nosehoover_thermostat;
mod pressure_calculation;
mod remove_com_motion;
mod rottrans;
mod settle_algorithm;
mod shake_algorithm;
mod steepest_descent;
mod temperature;

pub use berendsen_barostat::{BerendsenBarostat, BerendsenBarostatParams};
pub use berendsen_thermostat::{BerendsenThermostat, BerendsenThermostatParams};
pub use energy::EnergyCalculation;
pub use forcefield::Forcefield;
pub use lattice_shift::LatticeShift;
pub use leap_frog::{LeapFrogPosition, LeapFrogVelocity};
pub use lincs_algorithm::LincsAlgorithm;
pub use multibath::{in_bath, BathRange};
pub use nosehoover_thermostat::{restored_zeta, NoseHooverThermostat, NoseHooverThermostatParams};
pub use pressure_calculation::{PressureCalculation, VirialType};
pub use remove_com_motion::RemoveCOMMotion;
pub use rottrans::RottransConstraints;
pub use settle_algorithm::SettleAlgorithm;
pub use shake_algorithm::ShakeAlgorithm;
pub use steepest_descent::SteepestDescentAlgorithm;
pub use temperature::TemperatureCalculation;
