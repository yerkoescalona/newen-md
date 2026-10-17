//! `InputParameters` — the deprecated, pre-recipe front-end (3.9 step 3; removed next
//! release). Every constructor warns and names its `Recipe` replacement; `Simulation` translates
//! an `InputParameters` into a recipe (`RunRecipe::from_imd_with`), so it is not a second path.

use pyo3::prelude::*;

use newen_io::imd::ImdParameters;

use super::recipe::{run_err, RecipeError};

/// Simulation parameters loaded from a GROMOS input file (.imd/.in).
///
/// Contains all MD control parameters: timestep, number of steps,
/// cutoffs, thermostat/barostat settings, constraints, etc.
///
/// # Example (Python)
///
/// ```python
/// params = InputParameters("run.imd")
/// print(f"dt={params.dt} ps, steps={params.nstlim}")
/// print(f"T={params.temperature} K, cutoff={params.cutoff} nm")
/// print(f"solvent molecules: {params.nsm}")
/// ```
#[pyclass(name = "InputParameters")]
pub struct PyInputParameters {
    pub(crate) inner: ImdParameters,
}

/// `InputParameters` is the pre-recipe front-end (3.9 step 3): every constructor
/// warns once per call and names its `Recipe` replacement. One deprecation cycle (A9).
pub(crate) fn deprecated(py: Python<'_>, what: &str, instead: &str) -> PyResult<()> {
    PyErr::warn_bound(
        py,
        &py.get_type_bound::<pyo3::exceptions::PyDeprecationWarning>(),
        &format!("{what} is deprecated and will be removed; use {instead}"),
        1,
    )
}

#[pymethods]
impl PyInputParameters {
    /// Load parameters from a GROMOS input file.
    ///
    /// Deprecated: use `Recipe.from_imd(path)`.
    ///
    /// Args:
    ///     input_file: Path to input file (.imd, .in)
    #[new]
    fn new(py: Python<'_>, input_file: &str) -> PyResult<Self> {
        deprecated(py, "InputParameters(path)", "Recipe.from_imd(path)")?;
        Self::read(input_file)
    }

    /// Load parameters from a GROMOS input file (alias of constructor).
    ///
    /// Deprecated: use `Recipe.from_imd(path)`.
    #[staticmethod]
    fn from_file(py: Python<'_>, input_file: &str) -> PyResult<Self> {
        deprecated(
            py,
            "InputParameters.from_file(path)",
            "Recipe.from_imd(path)",
        )?;
        Self::read(input_file)
    }

    /// NVE (microcanonical) parameters: no thermostat or barostat.
    ///
    /// Deprecated: use `Recipe.nve(dt, steps, constraints)`.
    ///
    /// `constraints`: `"none"` (default) | `"hbonds"` | `"allbonds"` — SHAKE-constrain
    /// solute bonds. A constrained system (e.g. one with solute H-bonds) run with
    /// `constraints="none"` will silently diverge; matching a constrained reference
    /// `.in` file needs `"hbonds"`/`"allbonds"` here.
    #[staticmethod]
    #[pyo3(signature = (dt, steps, constraints="none"))]
    fn nve(py: Python<'_>, dt: f64, steps: usize, constraints: &str) -> PyResult<Self> {
        deprecated(
            py,
            "InputParameters.nve(...)",
            "Recipe.nve(dt, steps, constraints)",
        )?;
        let inner = ImdParameters::nve(dt, steps, constraints).map_err(RecipeError::new_err)?;
        Ok(Self { inner })
    }

    /// NVT (canonical) parameters: Berendsen thermostat at `temperature` K.
    ///
    /// Deprecated: use `Recipe.nvt(dt, steps, temperature, constraints)`.
    #[staticmethod]
    #[pyo3(signature = (dt, steps, temperature, constraints="none"))]
    fn nvt(
        py: Python<'_>,
        dt: f64,
        steps: usize,
        temperature: f64,
        constraints: &str,
    ) -> PyResult<Self> {
        deprecated(
            py,
            "InputParameters.nvt(...)",
            "Recipe.nvt(dt, steps, temperature, constraints)",
        )?;
        let inner = ImdParameters::nvt(dt, steps, temperature, constraints)
            .map_err(RecipeError::new_err)?;
        Ok(Self { inner })
    }

    /// NPT (isothermal-isobaric) parameters: Berendsen thermostat + barostat.
    ///
    /// Deprecated: use `Recipe.npt(dt, steps, temperature, pressure, constraints)`.
    #[staticmethod]
    #[pyo3(signature = (dt, steps, temperature, pressure, constraints="none"))]
    fn npt(
        py: Python<'_>,
        dt: f64,
        steps: usize,
        temperature: f64,
        pressure: f64,
        constraints: &str,
    ) -> PyResult<Self> {
        deprecated(
            py,
            "InputParameters.npt(...)",
            "Recipe.npt(dt, steps, temperature, pressure, constraints)",
        )?;
        let inner = ImdParameters::npt(dt, steps, temperature, pressure, constraints)
            .map_err(RecipeError::new_err)?;
        Ok(Self { inner })
    }

    /// Steepest-descent energy minimization parameters.
    ///
    /// Deprecated: use `Recipe.minimize(steps)`.
    #[staticmethod]
    fn steepest_descent(py: Python<'_>, steps: usize) -> PyResult<Self> {
        deprecated(
            py,
            "InputParameters.steepest_descent(steps)",
            "Recipe.minimize(steps)",
        )?;
        Ok(Self {
            inner: ImdParameters::steepest_descent(steps),
        })
    }

    /// Time step in picoseconds.
    #[getter]
    fn dt(&self) -> f64 {
        self.inner.dt
    }

    /// Number of simulation steps.
    #[getter]
    fn nstlim(&self) -> usize {
        self.inner.nstlim
    }

    /// Number of solvent molecules.
    #[getter]
    fn nsm(&self) -> usize {
        self.inner.nsm
    }

    /// Long-range nonbonded cutoff (nm).
    #[getter]
    fn cutoff(&self) -> f64 {
        self.inner.rcutl
    }

    /// Short-range pairlist cutoff (nm).
    #[getter]
    fn rcutp(&self) -> f64 {
        self.inner.rcutp
    }

    /// Target temperature (K) of the first temperature bath.
    #[getter]
    fn temperature(&self) -> f64 {
        if !self.inner.temp_bath.is_empty() && !self.inner.temp_bath[0].temp0.is_empty() {
            self.inner.temp_bath[0].temp0[0]
        } else {
            0.0
        }
    }

    /// SHAKE constraint mode (1=none, 2=H-bonds, 3=all bonds).
    #[getter]
    fn ntc(&self) -> i32 {
        self.inner.ntc
    }

    /// SHAKE constraint mode as the `nve`/`nvt`/`npt` convenience string
    /// (`"none"`/`"hbonds"`/`"allbonds"`) — the inverse of the `constraints=`
    /// factory argument.
    #[getter]
    fn constraints(&self) -> &'static str {
        match self.inner.ntc {
            3 => "allbonds",
            2 => "hbonds",
            _ => "none",
        }
    }

    /// Boundary condition type (0=vacuum, 1=rectangular).
    #[getter]
    fn ntb(&self) -> i32 {
        self.inner.ntb
    }

    /// Pairlist update frequency (steps).
    #[getter]
    fn nsnb(&self) -> usize {
        self.inner.nsnb
    }

    /// Trajectory write frequency (steps).
    #[getter]
    fn ntwx(&self) -> usize {
        self.inner.ntwx
    }

    /// Energy write frequency (steps).
    #[getter]
    fn ntwe(&self) -> usize {
        self.inner.ntwe
    }

    fn __repr__(&self) -> String {
        format!(
            "InputParameters(dt={}, nstlim={}, nsm={}, cutoff={:.3}, T={:.1})",
            self.inner.dt,
            self.inner.nstlim,
            self.inner.nsm,
            self.inner.rcutl,
            self.temperature(),
        )
    }
}

impl PyInputParameters {
    fn read(input_file: &str) -> PyResult<Self> {
        let imd = newen_run::read_imd(input_file).map_err(run_err)?;
        Ok(Self { inner: imd })
    }
}

pub fn register_parameters(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyInputParameters>()?;
    Ok(())
}
