//! Replica exchange type definitions for I/O
//!
//! These types are used for parsing replica exchange input blocks.

/// Replica exchange type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReplicaExchangeType {
    /// Temperature REMD
    Temperature,
    /// Hamiltonian REMD (λ-REMD)
    Hamiltonian,
    /// REST2 (replica exchange with solute tempering)
    REST2,
}

impl Default for ReplicaExchangeType {
    fn default() -> Self {
        ReplicaExchangeType::Temperature
    }
}

/// Replica exchange parameters from input file
#[derive(Debug, Clone, Default)]
pub struct ReplicaParameters {
    pub enabled: bool,
    pub exchange_type: ReplicaExchangeType,
    pub num_replicas: usize,
    pub exchange_interval: u64,
    pub temperatures: Vec<f64>,
    pub lambda_values: Vec<f64>,
}
