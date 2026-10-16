//! gromos-run — the one place that turns "what to run" into a running `AlgorithmSequence`.
//!
//! L3 orchestration (see `.claude/architecture.md`). Before this crate existed, the `md` binary
//! and the Python binding each carried their own copy of the IMD → algorithm-sequence assembly and
//! had drifted (3.8). Both now call the functions here, so the gromosXX reference suite
//! (which drives the binary) and the Python suite (which drives the binding) exercise the same
//! code from two sides — 3.9, complete 2026-08-29.
//!
//! The description of a run is plain data, [`RunRecipe`]: the content of a GROMOS `.imd` grouped
//! by concern, plus additive force [`TermSpec`]s and auxiliary [`RunInputs`]. A `.imd` file and a
//! Python `Recipe` are two front-ends to it ([`RunRecipe::from_imd_with`] / [`RunRecipe::to_imd`]
//! are lossless in both directions).
//!
//! Pipeline:
//!
//! ```text
//! .imd ──read_imd──▶ ImdParameters ──RunRecipe::from_imd_with──▶ RunRecipe  ◀── TOML/JSON, Python
//!                                                                    │
//!   Topology + Coordinates + RunInputs ──prepare_system──▶ Prepared    │
//!   (pttopo merge, truncated-octahedron transform, NSM from            │
//!    the coordinate file, initial velocities, validation)              │
//!                                                                      ▼
//!                                   build_plan (recipe, topology) ──▶ Vec<AlgorithmSpec>   stage 1
//!                                   validate_plan                      (GROMOS step order)  stage 2
//!                                   instantiate (plan only, G8)   ──▶ AlgorithmSequence     stage 3
//!                                   start (init + step 0)              then the caller loops
//! ```
//!
//! [`build_sequence_from_recipe`] / [`build_sequence_from_plan`] / [`build_sequence_from_imd`]
//! run stages 1–3 in one call; [`load_bundle`] / [`write_bundle`] read and write a run directory
//! (`input.toml`, `run.recipe.toml`, `run.imd`).
//!
//! Rules (crate contract, enforced by `just lint`'s G6 greps):
//! - No `println!`, no `process::exit`: everything the binary prints comes back as data
//!   ([`Prepared::notes`], [`BuildSummary`], [`Diagnostics`]) or as a [`RunError`].
//! - No second builder, no second IMD reader: a caller that needs a different sequence edits the
//!   plan; it never pushes algorithms itself, and it opens parameter files only through
//!   [`read_imd`].
//! - Defaults come from `ImdParameters::default()` and mean "block absent" — what gromosXX does
//!   without the block (G7/A18).
//! - Adding a force term = a [`TermSpec`] variant, its registry arms in `plan.rs`, and one arm in
//!   `build.rs::instantiate_orchestrator` (measured with `xtb`: three files).

pub mod build;
pub mod bundle;
pub mod constraints;
pub mod dof;
pub mod error;
pub mod inputs;
#[cfg(test)]
mod legacy_builder;
pub mod plan;
pub mod prepare;
pub mod recipe;

pub use build::{
    build_sequence_from_imd, build_sequence_from_plan, build_sequence_from_recipe, instantiate,
    periodicity_of, start, BarostatSummary, BuildSummary, Built, Instantiated, ThermostatSummary,
};
pub use bundle::{load_bundle, read_imd, write_bundle, RunBundle};
pub use constraints::ConstraintSelection;
pub use dof::{total_dof, BathDof};
pub use error::RunError;
pub use inputs::{ParallelPolicy, RunInputs, RunOptions};
pub use plan::{
    build_plan, plan_from_json, plan_to_json, validate_plan, AedsPlan, AlgorithmSpec, EdsPlan,
};
pub use prepare::{prepare_system, Coordinates, PrepareNote, Prepared};
pub use recipe::{Aeds, Diagnostics, Eds, PassthroughPolicy, RunRecipe, TermSpec, RECIPE_VERSION};
