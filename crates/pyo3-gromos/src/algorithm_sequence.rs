//! `AlgorithmSequence` — a deprecation shim (one release) over `recipe.plan(system)`.
//!
//! 3.9 step 4 deleted the descriptor path: the twelve `Forcefield`/`LeapFrog…`
//! pyclasses, the `AlgorithmDescriptor` enum and `resolve_algorithm_sequence` — the second
//! IMD→sequence builder. The MD step is `gromos.Plan`, built by `gromos-run` from a recipe;
//! this file builds nothing. What survives, warning, are the preset names, each returning that
//! `Plan`:
//!
//! ```python
//! plan = AlgorithmSequence.from_parameters(topo, params)    # -> Plan (warns)
//! sim = Simulation.from_sequence(topo, conf, params, plan)  # (warns)
//! # today:
//! plan = recipe.plan(system)
//! sim = Simulation(system, recipe, plan=plan)
//! ```

use pyo3::prelude::*;

use gromos_run::recipe::{PassthroughPolicy, RunRecipe};
use gromos_run::{build_plan, validate_plan};

use super::parameters::{deprecated, PyInputParameters};
use super::recipe::{run_err, PyPlan};
use super::topology::PyTopology;

/// Deprecated: use `recipe.plan(system)`. Cannot be instantiated; its presets return a `Plan`.
#[pyclass(name = "AlgorithmSequence", module = "gromos")]
pub struct PyAlgorithmSequence;

/// The plan of `params` on `topo` — what `Recipe.from_imd(...).plan(system)` builds, minus the
/// coordinates: without them NSM comes from the parameters (`prepare_system` takes it from the
/// coordinate file).
fn plan_for(topo: &PyTopology, params: &PyInputParameters) -> PyResult<PyPlan> {
    let (recipe, _) =
        RunRecipe::from_imd_with(&params.inner, &PassthroughPolicy::default()).map_err(run_err)?;
    let mut topology = topo.inner.clone();
    let unsolvated = topology.num_atoms() == topology.num_solute_atoms();
    if unsolvated && params.inner.nsm > 0 && !topology.solvent_atom_template.is_empty() {
        topology.solvate(params.inner.nsm);
    }
    let specs = build_plan(&recipe, &topology, &topo.physical_constants).map_err(run_err)?;
    validate_plan(&specs).map_err(run_err)?;
    Ok(PyPlan { specs })
}

fn preset(
    py: Python<'_>,
    what: &str,
    topo: &PyTopology,
    params: &PyInputParameters,
) -> PyResult<PyPlan> {
    deprecated(
        py,
        &format!("AlgorithmSequence.{what}(topo, params)"),
        "Recipe.from_imd(path) (or Recipe.nve/nvt/npt/minimize) and recipe.plan(system)",
    )?;
    plan_for(topo, params)
}

#[pymethods]
impl PyAlgorithmSequence {
    /// Deprecated: the plan of `params` (its ensemble comes from `params`, not from the name).
    #[staticmethod]
    fn nve(py: Python<'_>, topo: &PyTopology, params: &PyInputParameters) -> PyResult<PyPlan> {
        preset(py, "nve", topo, params)
    }

    /// Deprecated: the plan of `params`.
    #[staticmethod]
    fn nvt(py: Python<'_>, topo: &PyTopology, params: &PyInputParameters) -> PyResult<PyPlan> {
        preset(py, "nvt", topo, params)
    }

    /// Deprecated: the plan of `params`.
    #[staticmethod]
    fn npt(py: Python<'_>, topo: &PyTopology, params: &PyInputParameters) -> PyResult<PyPlan> {
        preset(py, "npt", topo, params)
    }

    /// Deprecated: the plan of `params`.
    #[staticmethod]
    fn minimize(py: Python<'_>, topo: &PyTopology, params: &PyInputParameters) -> PyResult<PyPlan> {
        preset(py, "minimize", topo, params)
    }

    /// Deprecated: the plan of `params` — use `recipe.plan(system)`.
    #[staticmethod]
    fn from_parameters(
        py: Python<'_>,
        topo: &PyTopology,
        params: &PyInputParameters,
    ) -> PyResult<PyPlan> {
        preset(py, "from_parameters", topo, params)
    }
}

pub fn register_algorithm_sequence(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAlgorithmSequence>()?;
    Ok(())
}
