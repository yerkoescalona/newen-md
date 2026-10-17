//! The recipe model from Python — `Recipe`, `Term`, `Algorithm`, `Plan` (3.9 step 3).
//!
//! These classes hold the **Rust data types** (`RunRecipe`, `TermSpec`, `AlgorithmSpec`); the
//! bridge is serde through `pythonize`, so a new field or variant on the Rust side is accepted
//! from Python with no binding code, and a typo'd key is an error naming it
//! (`deny_unknown_fields`). Nothing here holds physics or defaults — 3.9 A12/G7.
//!
//! ```python
//! recipe = Recipe.nvt(dt=0.002, steps=1000, temperature=300.0, constraints="hbonds")
//! recipe = recipe.update(outputs={"energy_every": 10}).with_term(Term("schnet", ...))
//! plan = recipe.plan(system)            # the MD step as data, editable
//! sim = Simulation(system, recipe)      # or Simulation(system, recipe, plan=plan)
//! recipe.to_imd(n_atoms=system.n_atoms) # what gromosXX would run
//! ```

use pyo3::create_exception;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyType};
use pythonize::{depythonize_bound, pythonize};
use serde::de::DeserializeOwned;
use serde::Serialize;

use newen_run::plan::{build_plan, validate_plan, AlgorithmSpec};
use newen_run::recipe::{PassthroughPolicy, RunRecipe, TermSpec, RECIPE_VERSION};
use newen_run::{prepare_system, Coordinates};

use super::system::PySystem;

create_exception!(
    newen,
    RecipeError,
    PyValueError,
    "A run recipe, term or algorithm is invalid, or the run cannot be represented as one."
);
create_exception!(
    newen,
    PlanError,
    PyValueError,
    "An algorithm plan violates the GROMOS step-order invariants."
);
create_exception!(
    newen,
    MissingFeatureError,
    PyRuntimeError,
    "A term or algorithm needs a cargo feature this build does not include."
);
create_exception!(
    newen,
    RunError,
    PyRuntimeError,
    "The engine failed to initialise or to take a step."
);

/// Map every `newen_run::RunError` onto the `newen.exceptions` hierarchy. The classes subclass
/// the builtin exceptions the binding raised before (ValueError / RuntimeError / IOError), so
/// existing `except` clauses keep working.
pub(crate) fn run_err(e: newen_run::RunError) -> PyErr {
    use newen_run::RunError as E;
    match e {
        E::Io { .. } => PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()),
        E::MissingInput { .. }
        | E::AtomCountMismatch { .. }
        | E::SolventCount { .. }
        | E::Inconsistent(_)
        | E::Validation { .. }
        | E::UnknownBlock { .. }
        | E::Recipe(_)
        | E::Serde(_) => RecipeError::new_err(e.to_string()),
        E::InvalidPlan(_) => PlanError::new_err(e.to_string()),
        E::MissingFeature { .. } => MissingFeatureError::new_err(e.to_string()),
        E::Init(_) | E::Step { .. } => RunError::new_err(e.to_string()),
    }
}

fn to_py<'py, T: Serialize + ?Sized>(py: Python<'py>, value: &T) -> PyResult<Bound<'py, PyAny>> {
    pythonize(py, value)
        .map(|o| o.into_bound(py))
        .map_err(|e| RecipeError::new_err(e.to_string()))
}

fn from_py<T: DeserializeOwned>(obj: &Bound<'_, PyAny>, what: &str) -> PyResult<T> {
    depythonize_bound(obj.clone()).map_err(|e| RecipeError::new_err(format!("{what}: {e}")))
}

/// `{"kind": kind, **params}` — the tagged form serde expects.
fn tagged<'py>(
    py: Python<'py>,
    kind: &str,
    params: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("kind", kind)?;
    if let Some(p) = params {
        for (k, v) in p.iter() {
            d.set_item(k, v)?;
        }
    }
    Ok(d)
}

fn enabled_features() -> Vec<&'static str> {
    let mut v = Vec::new();
    if cfg!(feature = "ml") {
        v.push("ml");
    }
    v
}

fn feature_available(feature: Option<&str>) -> bool {
    feature.is_none_or(|f| enabled_features().contains(&f))
}

// ---------------------------------------------------------------------------------------------
// Term
// ---------------------------------------------------------------------------------------------

/// An additive force provider (`newen.terms()` lists the kinds and their parameters).
#[pyclass(name = "Term", module = "newen")]
#[derive(Clone)]
pub struct PyTerm {
    pub(crate) inner: TermSpec,
}

#[pymethods]
impl PyTerm {
    /// `Term(kind, **params)` — validated against the Rust `TermSpec`: an unknown kind or
    /// parameter raises `RecipeError` naming it.
    #[new]
    #[pyo3(signature = (kind, **params))]
    fn new(py: Python<'_>, kind: &str, params: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        if !TermSpec::KINDS.contains(&kind) {
            return Err(RecipeError::new_err(format!(
                "unknown term kind {kind:?}; known kinds: {}",
                TermSpec::KINDS.join(", ")
            )));
        }
        let d = tagged(py, kind, params)?;
        Ok(Self {
            inner: from_py(&d, &format!("Term({kind:?})"))?,
        })
    }

    #[getter]
    fn kind(&self) -> &'static str {
        self.inner.name()
    }

    /// The cargo feature this term needs, if any (`None` = always available).
    #[getter]
    fn feature(&self) -> Option<&'static str> {
        self.inner.feature()
    }

    /// Whether this build can instantiate the term.
    #[getter]
    fn available(&self) -> bool {
        feature_available(self.inner.feature())
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &self.inner)
    }

    #[staticmethod]
    fn from_dict(d: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: from_py(d, "Term.from_dict")?,
        })
    }

    fn __eq__(&self, other: &PyTerm) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!("Term({})", self.to_dict(py)?.repr()?))
    }
}

// ---------------------------------------------------------------------------------------------
// Algorithm
// ---------------------------------------------------------------------------------------------

/// One algorithm of the MD step, fully resolved (`newen.algorithms()` lists the kinds).
#[pyclass(name = "Algorithm", module = "newen")]
#[derive(Clone)]
pub struct PyAlgorithm {
    pub(crate) inner: AlgorithmSpec,
}

#[pymethods]
impl PyAlgorithm {
    /// `Algorithm(kind, **params)` — validated against the Rust `AlgorithmSpec`.
    #[new]
    #[pyo3(signature = (kind, **params))]
    fn new(py: Python<'_>, kind: &str, params: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        if !AlgorithmSpec::KINDS.contains(&kind) {
            return Err(RecipeError::new_err(format!(
                "unknown algorithm kind {kind:?}; known kinds: {}",
                AlgorithmSpec::KINDS.join(", ")
            )));
        }
        let d = tagged(py, kind, params)?;
        Ok(Self {
            inner: from_py(&d, &format!("Algorithm({kind:?})"))?,
        })
    }

    #[getter]
    fn kind(&self) -> &'static str {
        self.inner.name()
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &self.inner)
    }

    #[staticmethod]
    fn from_dict(d: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: from_py(d, "Algorithm.from_dict")?,
        })
    }

    fn __eq__(&self, other: &PyAlgorithm) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!("Algorithm({})", self.to_dict(py)?.repr()?))
    }
}

// ---------------------------------------------------------------------------------------------
// Plan
// ---------------------------------------------------------------------------------------------

/// The MD step as an ordered, editable list of `Algorithm`s (stage 1 of the builder). Every
/// value is resolved; `validate()` (also run by `Simulation`) enforces the GROMOS ordering
/// invariants. Targets for the edit verbs are an index or a kind name (a kind that occurs more
/// than once must be addressed by index).
#[pyclass(name = "Plan", module = "newen")]
#[derive(Clone)]
pub struct PyPlan {
    pub(crate) specs: Vec<AlgorithmSpec>,
}

impl PyPlan {
    fn resolve(&self, target: &Bound<'_, PyAny>) -> PyResult<usize> {
        if let Ok(i) = target.extract::<isize>() {
            let n = self.specs.len() as isize;
            let idx = if i < 0 { n + i } else { i };
            if idx < 0 || idx >= n {
                return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(format!(
                    "index {i} out of range for a plan of {n} algorithms"
                )));
            }
            return Ok(idx as usize);
        }
        let kind: String = target.extract().map_err(|_| {
            PlanError::new_err("target must be an index (int) or an algorithm kind (str)")
        })?;
        let hits: Vec<usize> = self
            .specs
            .iter()
            .enumerate()
            .filter(|(_, s)| s.name() == kind)
            .map(|(i, _)| i)
            .collect();
        match hits.as_slice() {
            [i] => Ok(*i),
            [] => Err(PlanError::new_err(format!(
                "no {kind:?} in the plan (kinds: {})",
                self.kinds().join(", ")
            ))),
            many => Err(PlanError::new_err(format!(
                "{kind:?} occurs {} times (indices {many:?}); address it by index",
                many.len()
            ))),
        }
    }
}

#[pymethods]
impl PyPlan {
    #[new]
    #[pyo3(signature = (algorithms=None))]
    fn new(algorithms: Option<Vec<PyAlgorithm>>) -> Self {
        Self {
            specs: algorithms
                .unwrap_or_default()
                .into_iter()
                .map(|a| a.inner)
                .collect(),
        }
    }

    /// Kind names in step order.
    #[getter]
    fn kinds(&self) -> Vec<&'static str> {
        self.specs.iter().map(|s| s.name()).collect()
    }

    fn __len__(&self) -> usize {
        self.specs.len()
    }

    fn __getitem__(&self, target: &Bound<'_, PyAny>) -> PyResult<PyAlgorithm> {
        let i = self.resolve(target)?;
        Ok(PyAlgorithm {
            inner: self.specs[i].clone(),
        })
    }

    fn __contains__(&self, kind: &str) -> bool {
        self.specs.iter().any(|s| s.name() == kind)
    }

    /// The algorithms as a list (copies; edit through the verbs below).
    fn to_list(&self) -> Vec<PyAlgorithm> {
        self.specs
            .iter()
            .cloned()
            .map(|inner| PyAlgorithm { inner })
            .collect()
    }

    fn insert(&mut self, index: usize, algorithm: &PyAlgorithm) -> PyResult<()> {
        if index > self.specs.len() {
            return Err(PlanError::new_err(format!(
                "index {index} out of range for a plan of {} algorithms",
                self.specs.len()
            )));
        }
        self.specs.insert(index, algorithm.inner.clone());
        Ok(())
    }

    fn insert_after(&mut self, target: &Bound<'_, PyAny>, algorithm: &PyAlgorithm) -> PyResult<()> {
        let i = self.resolve(target)?;
        self.specs.insert(i + 1, algorithm.inner.clone());
        Ok(())
    }

    fn insert_before(
        &mut self,
        target: &Bound<'_, PyAny>,
        algorithm: &PyAlgorithm,
    ) -> PyResult<()> {
        let i = self.resolve(target)?;
        self.specs.insert(i, algorithm.inner.clone());
        Ok(())
    }

    fn remove(&mut self, target: &Bound<'_, PyAny>) -> PyResult<PyAlgorithm> {
        let i = self.resolve(target)?;
        Ok(PyAlgorithm {
            inner: self.specs.remove(i),
        })
    }

    fn replace(&mut self, target: &Bound<'_, PyAny>, algorithm: &PyAlgorithm) -> PyResult<()> {
        let i = self.resolve(target)?;
        self.specs[i] = algorithm.inner.clone();
        Ok(())
    }

    /// Check the GROMOS ordering invariants; raises `PlanError` naming the violation.
    fn validate(&self) -> PyResult<()> {
        validate_plan(&self.specs).map_err(run_err)
    }

    fn to_json(&self) -> PyResult<String> {
        newen_run::plan_to_json(&self.specs).map_err(run_err)
    }

    #[staticmethod]
    fn from_json(text: &str) -> PyResult<Self> {
        Ok(Self {
            specs: newen_run::plan_from_json(text).map_err(run_err)?,
        })
    }

    /// The plan as a list of dicts (one per algorithm, `kind` + parameters).
    fn to_dicts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &self.specs)
    }

    fn __repr__(&self) -> String {
        format!("Plan([{}])", self.kinds().join(" → "))
    }
}

// ---------------------------------------------------------------------------------------------
// Recipe
// ---------------------------------------------------------------------------------------------

/// A run description — the same data an `.imd` file carries, grouped by concern, plus the
/// additive `terms` and the auxiliary `inputs`. Immutable: every edit returns a new `Recipe`.
#[pyclass(name = "Recipe", module = "newen")]
#[derive(Clone)]
pub struct PyRecipe {
    pub(crate) inner: RunRecipe,
    /// What `from_imd` defaulted or passed through (empty for recipes built in memory).
    pub(crate) diagnostics: Vec<String>,
}

impl PyRecipe {
    fn wrap(inner: RunRecipe) -> Self {
        Self {
            inner,
            diagnostics: Vec::new(),
        }
    }

    fn merged(&self, py: Python<'_>, groups: Option<&Bound<'_, PyDict>>) -> PyResult<RunRecipe> {
        let Some(groups) = groups else {
            return Ok(self.inner.clone());
        };
        let mut base =
            serde_json::to_value(&self.inner).map_err(|e| RecipeError::new_err(e.to_string()))?;
        let patch: serde_json::Value = from_py(groups.as_any(), "Recipe.update")?;
        let _ = py;
        merge(&mut base, patch);
        serde_json::from_value(base)
            .map_err(|e| RecipeError::new_err(format!("Recipe.update: {e}")))
    }

    fn group<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        let v =
            serde_json::to_value(&self.inner).map_err(|e| RecipeError::new_err(e.to_string()))?;
        to_py(py, &v[name])
    }
}

/// Deep-merge `patch` into `base`: objects merge key by key, everything else is replaced.
fn merge(base: &mut serde_json::Value, patch: serde_json::Value) {
    match (base, patch) {
        (serde_json::Value::Object(b), serde_json::Value::Object(p)) => {
            for (k, v) in p {
                match b.get_mut(&k) {
                    Some(slot) => merge(slot, v),
                    None => {
                        b.insert(k, v);
                    },
                }
            }
        },
        (b, p) => *b = p,
    }
}

#[pymethods]
impl PyRecipe {
    /// `Recipe(**groups)` — gromosXX's defaults (every optional block absent) patched with the
    /// given groups, e.g. `Recipe(control={"steps": 100, "dt": 0.002})`.
    #[new]
    #[pyo3(signature = (**groups))]
    fn new(py: Python<'_>, groups: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let base = Self::wrap(RunRecipe::default());
        Ok(Self::wrap(base.merged(py, groups)?))
    }

    /// The recipe of a GROMOS `.imd` file. Unmodelled blocks are an error unless listed in
    /// `allow_passthrough`. Absent optional blocks are reported in `.diagnostics`.
    #[staticmethod]
    #[pyo3(signature = (path, allow_passthrough=None))]
    pub(crate) fn from_imd(path: &str, allow_passthrough: Option<Vec<String>>) -> PyResult<Self> {
        let imd = newen_run::read_imd(path).map_err(run_err)?;
        let policy = PassthroughPolicy::allow(allow_passthrough.unwrap_or_default());
        let (inner, diag) = RunRecipe::from_imd_with(&imd, &policy).map_err(run_err)?;
        Ok(Self {
            inner,
            diagnostics: diag.notes,
        })
    }

    #[staticmethod]
    fn from_dict(d: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self::wrap(from_py(d, "Recipe.from_dict")?))
    }

    #[staticmethod]
    fn from_toml(text: &str) -> PyResult<Self> {
        Ok(Self::wrap(RunRecipe::from_toml(text).map_err(run_err)?))
    }

    #[staticmethod]
    fn from_json(text: &str) -> PyResult<Self> {
        Ok(Self::wrap(RunRecipe::from_json(text).map_err(run_err)?))
    }

    /// NVE: no thermostat, no barostat. `constraints`: "none" | "hbonds" | "allbonds".
    #[staticmethod]
    #[pyo3(signature = (dt, steps, constraints="none"))]
    fn nve(dt: f64, steps: usize, constraints: &str) -> PyResult<Self> {
        Ok(Self::wrap(
            RunRecipe::nve(dt, steps, constraints).map_err(run_err)?,
        ))
    }

    /// NVT: Berendsen thermostat at `temperature` K (tau 0.1 ps).
    #[staticmethod]
    #[pyo3(signature = (dt, steps, temperature, constraints="none"))]
    fn nvt(dt: f64, steps: usize, temperature: f64, constraints: &str) -> PyResult<Self> {
        Ok(Self::wrap(
            RunRecipe::nvt(dt, steps, temperature, constraints).map_err(run_err)?,
        ))
    }

    /// NPT: Berendsen thermostat + barostat (`pressure` in bar).
    #[staticmethod]
    #[pyo3(signature = (dt, steps, temperature, pressure, constraints="none"))]
    fn npt(
        dt: f64,
        steps: usize,
        temperature: f64,
        pressure: f64,
        constraints: &str,
    ) -> PyResult<Self> {
        Ok(Self::wrap(
            RunRecipe::npt(dt, steps, temperature, pressure, constraints).map_err(run_err)?,
        ))
    }

    /// Steepest-descent energy minimisation.
    #[staticmethod]
    fn minimize(steps: usize) -> Self {
        Self::wrap(RunRecipe::minimize(steps))
    }

    // --- views -------------------------------------------------------------------------------

    #[getter]
    fn version(&self) -> u32 {
        self.inner.version
    }

    #[getter]
    fn title(&self) -> String {
        self.inner.title.clone()
    }

    #[getter]
    fn diagnostics(&self) -> Vec<String> {
        self.diagnostics.clone()
    }

    #[getter]
    fn terms(&self) -> Vec<PyTerm> {
        self.inner
            .forcefield
            .terms
            .iter()
            .cloned()
            .map(|inner| PyTerm { inner })
            .collect()
    }

    #[getter]
    fn system<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "system")
    }
    #[getter]
    fn control<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "control")
    }
    #[getter]
    fn boundary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "boundary")
    }
    #[getter]
    fn forcefield<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "forcefield")
    }
    #[getter]
    fn constraints<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "constraints")
    }
    #[getter]
    fn ensemble<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "ensemble")
    }
    #[getter]
    fn minimisation<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "minimisation")
    }
    #[getter]
    fn perturbation<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "perturbation")
    }
    /// The `EDS` block as a dict (`enabled`, `form`, `num_states`, `s`, `pairs`, `eir`), or
    /// `None` without one.
    #[getter]
    fn eds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "eds")
    }
    /// The `AEDS` block as a dict (`enabled`, `form`, `num_states`, `emax`, `emin`, `eir`,
    /// `ntiaedss`, `restremin`, `bmaxtype`, `bmax`, `asteps`, `bsteps`), or `None` without one.
    /// `energy_columns`, `prevalence` and `target_state` (no `.imd` block; listed only when not
    /// `"intended"`) reproduce gromosXX's A-EDS defects (`"gromosxx-post-update"`,
    /// `"gromosxx-plain-exp"`, `"gromosxx-state1-default"`).
    #[getter]
    fn aeds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "aeds")
    }
    #[getter]
    fn outputs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "outputs")
    }
    #[getter]
    fn execution<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "execution")
    }
    #[getter]
    fn inputs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "inputs")
    }
    #[getter]
    fn passthrough<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.group(py, "passthrough")
    }

    fn is_minimization(&self) -> bool {
        self.inner.is_minimization()
    }

    // --- edits (each returns a new Recipe) ---------------------------------------------------

    /// Patch groups: `recipe.update(control={"steps": 50}, ensemble={"thermostat": None})`.
    /// Dicts merge key by key; everything else is replaced. Unknown keys raise `RecipeError`.
    #[pyo3(signature = (**groups))]
    fn update(&self, py: Python<'_>, groups: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        Ok(Self {
            inner: self.merged(py, groups)?,
            diagnostics: self.diagnostics.clone(),
        })
    }

    /// Append an additive term (`newen.Term`).
    fn with_term(&self, term: &PyTerm) -> Self {
        let mut r = self.clone();
        r.inner.forcefield.terms.push(term.inner.clone());
        r
    }

    fn without_terms(&self) -> Self {
        let mut r = self.clone();
        r.inner.forcefield.terms.clear();
        r
    }

    /// Set the auxiliary input files (the `md` binary's `@pttopo`/`@posresspec`/`@refpos`/
    /// `@distrest`). `None` leaves a path unchanged.
    #[pyo3(signature = (pttopo=None, posresspec=None, refpos=None, distrest=None))]
    fn with_inputs(
        &self,
        pttopo: Option<String>,
        posresspec: Option<String>,
        refpos: Option<String>,
        distrest: Option<String>,
    ) -> Self {
        let mut r = self.clone();
        let i = &mut r.inner.inputs;
        if let Some(p) = pttopo {
            i.pttopo = Some(p.into());
        }
        if let Some(p) = posresspec {
            i.posresspec = Some(p.into());
        }
        if let Some(p) = refpos {
            i.refpos = Some(p.into());
        }
        if let Some(p) = distrest {
            i.distrest = Some(p.into());
        }
        r
    }

    /// The run with the given execution policy (`parallel`: "auto" | "serial" | "parallel").
    fn with_execution(&self, py: Python<'_>, parallel: &str) -> PyResult<Self> {
        let d = PyDict::new_bound(py);
        let e = PyDict::new_bound(py);
        e.set_item("parallel", parallel)?;
        d.set_item("execution", e)?;
        self.update(py, Some(&d))
    }

    // --- serialisation -----------------------------------------------------------------------

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        to_py(py, &self.inner)
    }

    fn to_toml(&self) -> PyResult<String> {
        self.inner.to_toml().map_err(run_err)
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json().map_err(run_err)
    }

    /// The GROMOS `.imd` text for this run (gromosXX reads it). `n_atoms` fills the atom counts
    /// a recipe built in memory does not know (MULTIBATH DOFSET, FORCE NRE).
    #[pyo3(signature = (n_atoms=None))]
    fn to_imd(&self, n_atoms: Option<usize>) -> String {
        self.inner.to_imd_string(n_atoms)
    }

    /// Write the `.imd` text to `path`.
    #[pyo3(signature = (path, n_atoms=None))]
    fn save_imd(&self, path: &str, n_atoms: Option<usize>) -> PyResult<()> {
        std::fs::write(path, self.inner.to_imd_string(n_atoms)).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("cannot write '{path}': {e}"))
        })
    }

    /// Write a run bundle (`input.toml`, `run.recipe.toml`, `run.imd`) into `directory`.
    fn to_bundle(
        &self,
        directory: &str,
        system: &PySystem,
        topology_path: &str,
        configuration_path: &str,
    ) -> PyResult<String> {
        let path = newen_run::write_bundle(
            std::path::Path::new(directory),
            &self.inner,
            std::path::Path::new(topology_path),
            std::path::Path::new(configuration_path),
            Some(system.topology.inner.num_atoms()),
        )
        .map_err(run_err)?;
        Ok(path.to_string_lossy().into_owned())
    }

    /// The recipe named by a bundle (`input.toml`): its `recipe` file if present, else its
    /// `parameters` file; the bundle's auxiliary paths become `inputs`.
    #[staticmethod]
    #[pyo3(signature = (path, allow_passthrough=None))]
    fn from_bundle(path: &str, allow_passthrough: Option<Vec<String>>) -> PyResult<Self> {
        let policy = PassthroughPolicy::allow(allow_passthrough.unwrap_or_default());
        let (_, inner, diag) = newen_run::load_bundle(path, &policy).map_err(run_err)?;
        Ok(Self {
            inner,
            diagnostics: diag.notes,
        })
    }

    // --- the plan ----------------------------------------------------------------------------

    /// Stage 1 of the builder for this recipe on `system`: the MD step as an editable `Plan`,
    /// already validated. Pass it back with `Simulation(system, recipe, plan=plan)`.
    fn plan(&self, system: &PySystem) -> PyResult<PyPlan> {
        let imd = self.inner.to_imd();
        let coords = Coordinates {
            positions: system.configuration.pos_data.clone(),
            velocities: system.configuration.vel_data.clone(),
            box_dims: system.configuration.box_dims,
            pertdata: system.configuration.pertdata,
            aedssearch: system.configuration.aedssearch.clone(),
        };
        let prepared = prepare_system(
            &imd,
            system.topology.inner.clone(),
            system.topology.physical_constants,
            coords,
            &self.inner.inputs,
        )
        .map_err(run_err)?;
        let specs = build_plan(
            &self.inner,
            &prepared.topology,
            &prepared.physical_constants,
        )
        .map_err(run_err)?;
        validate_plan(&specs).map_err(run_err)?;
        Ok(PyPlan { specs })
    }

    fn __eq__(&self, other: &PyRecipe) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> String {
        let r = &self.inner;
        // A bath with TAU < 0 is gromosXX's "no coupling": it is kept in the recipe (the
        // `.imd` round-trips) but it does not make the run NVT.
        let coupled = r
            .ensemble
            .thermostat
            .as_ref()
            .is_some_and(|t| t.baths.iter().any(|b| b.tau > 0.0));
        let ensemble = match (coupled, &r.ensemble.barostat) {
            (false, None) => "NVE",
            (true, None) => "NVT",
            (_, Some(_)) => "NPT",
        };
        format!(
            "Recipe({}, steps={}, dt={}, terms={}, {})",
            if r.is_minimization() { "EM" } else { ensemble },
            r.control.steps,
            r.control.dt,
            r.forcefield.terms.len(),
            r.title.lines().next().unwrap_or_default()
        )
    }

    /// Pickle support: the recipe round-trips through its TOML text.
    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, (String,))> {
        let cls = py.get_type_bound::<PyRecipe>();
        Ok((cls.getattr("from_toml")?, (self.to_toml()?,)))
    }
}

// ---------------------------------------------------------------------------------------------
// Registries and module
// ---------------------------------------------------------------------------------------------

fn strip_kind<'py>(py: Python<'py>, example: &impl Serialize) -> PyResult<Bound<'py, PyAny>> {
    let d = to_py(py, example)?;
    let d = d.downcast_into::<PyDict>()?;
    d.del_item("kind")?;
    Ok(d.into_any())
}

/// The term kinds this build knows: `[{"kind", "params" (example values), "feature",
/// "available"}, …]`.
#[pyfunction]
pub fn terms<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
    let out = PyList::empty_bound(py);
    for ex in TermSpec::examples() {
        let d = PyDict::new_bound(py);
        d.set_item("kind", ex.name())?;
        d.set_item("params", strip_kind(py, &ex)?)?;
        d.set_item("feature", ex.feature())?;
        d.set_item("available", feature_available(ex.feature()))?;
        out.append(d)?;
    }
    Ok(out)
}

/// The algorithm kinds and their ordering rules: `[{"kind", "params", "rules"}, …]`.
#[pyfunction]
pub fn algorithms<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
    let out = PyList::empty_bound(py);
    for ex in AlgorithmSpec::examples() {
        let d = PyDict::new_bound(py);
        d.set_item("kind", ex.name())?;
        d.set_item("params", strip_kind(py, &ex)?)?;
        let r = AlgorithmSpec::rules(ex.name());
        let rules = PyDict::new_bound(py);
        rules.set_item("unique", r.unique)?;
        rules.set_item("required", r.required)?;
        rules.set_item("first", r.first)?;
        rules.set_item("last", r.last)?;
        rules.set_item("must_follow", r.must_follow)?;
        rules.set_item("after", r.after)?;
        rules.set_item("before", r.before)?;
        rules.set_item("excludes", r.excludes)?;
        rules.set_item("requires", r.requires)?;
        d.set_item("rules", rules)?;
        out.append(d)?;
    }
    Ok(out)
}

/// What this extension was built with: version, recipe schema version, cargo features.
#[pyfunction]
pub fn build_info<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("version", env!("CARGO_PKG_VERSION"))?;
    d.set_item("recipe_version", RECIPE_VERSION)?;
    d.set_item("features", enabled_features())?;
    Ok(d)
}

pub fn register_recipe(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add_class::<PyRecipe>()?;
    m.add_class::<PyTerm>()?;
    m.add_class::<PyAlgorithm>()?;
    m.add_class::<PyPlan>()?;
    m.add_function(wrap_pyfunction!(terms, m)?)?;
    m.add_function(wrap_pyfunction!(algorithms, m)?)?;
    m.add_function(wrap_pyfunction!(build_info, m)?)?;
    m.add("RecipeError", py.get_type_bound::<RecipeError>())?;
    m.add("PlanError", py.get_type_bound::<PlanError>())?;
    m.add(
        "MissingFeatureError",
        py.get_type_bound::<MissingFeatureError>(),
    )?;
    m.add("RunError", py.get_type_bound::<RunError>())?;
    let _: &Bound<'_, PyType> = &py.get_type_bound::<PyRecipe>();
    Ok(())
}
