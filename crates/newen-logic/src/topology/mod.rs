//! Topology bookkeeping: which kind of nonbonded pair two atoms form, and which group an atom
//! belongs to.
//!
//! Theory: `spec/NewenTheory/Topology/Pairs.lean`, `spec/NewenTheory/Topology/Groups.lean`;
//! models and bridges: `spec/NewenSpec/Topology/Pairs.lean`, `spec/NewenSpec/Topology/Groups.lean`.

pub mod groups;
pub mod pairs;
