//! Free-energy corrections that are formulas of the simulation's set-up, not analyses of its
//! frames.

pub mod charge_correction;

pub use charge_correction::{
    charge_correction, ChargeCorrection, ChargeCorrectionInput, Scheme, SolventModel,
};
