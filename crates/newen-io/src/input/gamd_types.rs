//! GAMD type definitions for I/O
//!
//! These types are used for parsing GAMD input blocks.

/// GAMD boost type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GamdBoostType {
    /// Boost total potential energy
    TotalPotential,
    /// Boost dihedral potential energy
    DihedralPotential,
    /// Dual boost (both total and dihedral)
    Dual,
}

impl Default for GamdBoostType {
    fn default() -> Self {
        GamdBoostType::TotalPotential
    }
}

/// GAMD boost form (how boost is applied)
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BoostForm {
    #[default]
    Constant,
    Linear,
    Harmonic,
}

/// GAMD search mode for parameter estimation
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum SearchMode {
    #[default]
    Off,
    RunningAverage,
    AllTime,
}

/// Threshold type for boost calculation
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ThresholdType {
    #[default]
    Lower,
    Upper,
    Dual,
}

/// GAMD parameters from input file  
#[derive(Debug, Clone, Default)]
pub struct GamdParameters {
    pub enabled: bool,
    pub boost_type: GamdBoostType,
    pub sigma0: f64,
    pub equilibration_steps: u64,
    pub update_interval: u64,
}
