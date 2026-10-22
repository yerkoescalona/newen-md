//! Temperature and the thermostat: how many kinetic degrees of freedom a system has, and how they
//! are shared between temperature baths.
//!
//! Theory: `spec/NewenTheory/Dynamics/DegreesOfFreedom.lean`; models and bridges:
//! `spec/NewenSpec/Thermo/Dof.lean`, `spec/NewenSpec/Thermo/Bath.lean`.

pub mod bath;
pub mod dof;
