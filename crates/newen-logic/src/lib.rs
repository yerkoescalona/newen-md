//! Integer decision logic of newen with a Lean specification.
//!
//! Every public function here has a model in `spec/NewenSpec/` (named by its `/// spec:` line),
//! a Lean theorem that the model means what the doc comment says, and, when Aeneas can translate
//! it, a Lean proof that this Rust function computes the model. See `spec/README.md`.
#![forbid(unsafe_code)]

pub mod schedule;
pub mod thermo;
pub mod topology;
