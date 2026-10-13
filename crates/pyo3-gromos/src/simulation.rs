//! Interactive Simulation API for Python
//!
//! A `Simulation` is built from a `System` (topology + coordinates) and a `Recipe` (what to run)
//! through the same `gromos-run` functions the `md` binary calls (`prepare_system` →
//! `build_plan` → `validate_plan` → `instantiate` → `start`), so the gromosXX reference suite
//! (which drives the binary) and the Python suite guard one code path from two sides
//! (3.9). This file contains **no run assembly of its own**.
//!
//! # Example (Python)
//!
//! ```python
//! from gromos import Recipe, Simulation, System
//!
//! system = System.from_files("system.topo", "initial.cnf")
//! recipe = Recipe.from_imd("run.imd") # or Recipe.nvt(0.002, 1000, 300.0)
//! sim = Simulation(system, recipe)
//! sim.step(100)
//! print(sim.energies, sim.positions)
//!
//! plan = recipe.plan(system)                              # the MD step as data
//! plan.remove("remove_com")
//! sim = Simulation(system, recipe, plan=plan)
//! ```
//!
//! The pre-recipe forms — `Simulation(topo, conf, InputParameters(...))`, the `distrest=`/
//! `posresspec=`/`refpos=`/`ml_*=` keyword arguments and `Simulation.from_sequence` — still
//! work for one deprecation cycle, warn, and are translated into a recipe.

use numpy::PyArray2;
use pyo3::prelude::*;

use gromos_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::{BoxType, Configuration},
    math::{truncoct_triclinic_rotmat, Vec3},
    units::PhysicalConstants,
    Topology,
};
use gromos_io::{
    coordinate::read_coordinates,
    topology::{build_topology, read_topology_file},
};
use gromos_run::plan::AlgorithmSpec;
use gromos_run::recipe::{PassthroughPolicy, RunRecipe};
use gromos_run::{
    build_sequence_from_plan, build_sequence_from_recipe, load_bundle, prepare_system, start,
    Coordinates,
};

use super::parameters::{deprecated, PyInputParameters};
use super::py_conf::PyConfiguration;
use super::recipe::{run_err, MissingFeatureError, PyPlan, PyRecipe};
use super::system::PySystem;
use super::topology::PyTopology;
use super::PyEnergy;

/// The raw pieces every constructor reduces to.
struct Parts {
    topo: Topology,
    physical_constants: PhysicalConstants,
    positions: Vec<Vec3>,
    velocities: Vec<Vec3>,
    box_dims: Vec3,
    /// λ of the configuration file's PERTDATA block (NRDGL=1).
    pertdata: Option<f64>,
}

impl Parts {
    fn from_objects(topo: &PyTopology, conf: &PyConfiguration) -> Self {
        Self {
            topo: topo.inner.clone(),
            physical_constants: topo.physical_constants,
            positions: conf.pos_data.clone(),
            velocities: conf.vel_data.clone(),
            box_dims: conf.box_dims,
            pertdata: conf.pertdata,
        }
    }

    fn from_system(system: &PySystem) -> Self {
        Self::from_objects(&system.topology, &system.configuration)
    }

    fn from_files(topo_file: &str, conf_file: &str) -> PyResult<Self> {
        let parsed = read_topology_file(topo_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read topology '{}': {}",
                topo_file, e
            ))
        })?;
        let physical_constants = parsed.physical_constants;
        let topo = build_topology(parsed);
        let coord_data = read_coordinates(conf_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read coordinates '{}': {}",
                conf_file, e
            ))
        })?;
        Ok(Self {
            topo,
            physical_constants,
            positions: coord_data.positions,
            velocities: coord_data.velocities,
            box_dims: coord_data.box_dims,
            pertdata: coord_data.pertdata,
        })
    }

    fn into_pieces(self) -> (Topology, PhysicalConstants, Coordinates) {
        let coords = Coordinates {
            positions: self.positions,
            velocities: self.velocities,
            box_dims: self.box_dims,
            pertdata: self.pertdata,
        };
        (self.topo, self.physical_constants, coords)
    }
}

/// A recipe as a constructor argument: a `Recipe`, or (deprecated) an `InputParameters`.
fn recipe_arg(arg: &Bound<'_, PyAny>) -> PyResult<(RunRecipe, Vec<String>)> {
    if let Ok(r) = arg.extract::<PyRef<PyRecipe>>() {
        return Ok((r.inner.clone(), r.diagnostics.clone()));
    }
    if let Ok(p) = arg.extract::<PyRef<PyInputParameters>>() {
        deprecated(
            arg.py(),
            "Simulation(..., InputParameters)",
            "Simulation(system, Recipe.from_imd(path))",
        )?;
        let (recipe, diag) =
            RunRecipe::from_imd_with(&p.inner, &PassthroughPolicy::default()).map_err(run_err)?;
        return Ok((recipe, diag.notes));
    }
    Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
        "expected a Recipe (or a deprecated InputParameters)",
    ))
}

/// Deprecated keyword arguments, translated into the recipe (one deprecation cycle, A9).
#[allow(clippy::too_many_arguments)]
fn apply_legacy_kwargs(
    py: Python<'_>,
    recipe: &mut RunRecipe,
    distrest: Option<String>,
    posresspec: Option<String>,
    refpos: Option<String>,
    ml_potential: Option<&Bound<'_, PyAny>>,
    ml_region: Option<String>,
    ml_buffer: Option<String>,
) -> PyResult<()> {
    if distrest.is_some() || posresspec.is_some() || refpos.is_some() {
        deprecated(
            py,
            "Simulation(..., distrest=/posresspec=/refpos=)",
            "recipe.with_inputs(distrest=..., posresspec=..., refpos=...)",
        )?;
        if let Some(p) = distrest {
            recipe.inputs.distrest = Some(p.into());
        }
        if let Some(p) = posresspec {
            recipe.inputs.posresspec = Some(p.into());
        }
        if let Some(p) = refpos {
            recipe.inputs.refpos = Some(p.into());
        }
    }
    match (ml_potential, ml_region) {
        (None, None) => {},
        (Some(_), None) | (None, Some(_)) => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "ml_potential and ml_region must be given together",
            ))
        },
        (Some(potential), Some(region)) => {
            deprecated(
                py,
                "Simulation(..., ml_potential=, ml_region=, ml_buffer=)",
                "recipe.with_term(Term(\"schnet\", model=..., cutoff=..., elements=..., \
                 region=..., buffer=...))",
            )?;
            #[cfg(feature = "ml")]
            {
                use gromos_run::recipe::{Coupling, TermSpec};
                let potential = potential
                    .extract::<PyRef<crate::ml_potential::PySchNetPotential>>()
                    .map_err(|_| {
                        PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                            "ml_potential must be a SchNetPotential",
                        )
                    })?;
                recipe.forcefield.terms.push(TermSpec::Schnet {
                    model: potential.model_path.clone(),
                    cutoff: potential.cutoff,
                    elements: potential.elements.clone(),
                    region,
                    buffer: ml_buffer,
                    coupling: Coupling::Delta,
                });
            }
            #[cfg(not(feature = "ml"))]
            {
                let _ = (potential, region, ml_buffer);
                return Err(MissingFeatureError::new_err(
                    "ML support was not compiled into this build (rebuild pyo3-gromos/py-gromos \
                     with --features ml)",
                ));
            }
        },
    }
    Ok(())
}

/// The one way a `Simulation` comes into being: prepare the system, build (or take) the
/// plan, instantiate it, run `init` + step 0.
fn build_from_recipe(
    parts: Parts,
    recipe: RunRecipe,
    plan: Option<Vec<AlgorithmSpec>>,
    diagnostics: Vec<String>,
) -> PyResult<PySimulation> {
    let imd = recipe.to_imd();
    let (topo, physical_constants, coords) = parts.into_pieces();
    let prepared =
        prepare_system(&imd, topo, physical_constants, coords, &recipe.inputs).map_err(run_err)?;
    let built = match plan {
        Some(plan) => build_sequence_from_plan(&recipe, plan, &prepared),
        None => build_sequence_from_recipe(&recipe, &prepared),
    }
    .map_err(run_err)?;
    let gromos_run::Built {
        sequence: mut md_sequence,
        recipe,
        plan,
        summary,
        ..
    } = built;
    let gromos_run::Prepared {
        topology: topo,
        configuration: mut conf,
        ..
    } = prepared;
    let n_atoms = topo.num_atoms();
    let dt = recipe.control.dt;

    let mut sim_state = SimulationState::new(dt, recipe.control.steps);
    start(&mut md_sequence, &topo, &mut conf, &sim_state).map_err(run_err)?;
    sim_state.advance();

    Ok(PySimulation {
        topology: topo,
        configuration: conf,
        md_sequence,
        sim_state,
        dt,
        n_atoms,
        total_dof: summary.total_dof,
        recipe,
        plan,
        diagnostics,
    })
}

// ============================================================================
// PySimulation
// ============================================================================

/// A GROMOS-RS molecular dynamics simulation.
///
/// Interactive simulation object inspired by OpenMM's Simulation class,
/// using GROMOS naming conventions and file formats.
///
/// Create from a `System` and a `Recipe` (or from file paths). Call `step(n)` to advance and
/// access properties like `positions`, `forces`, `energies`, `temperature`.
#[pyclass(name = "Simulation", unsendable)]
pub struct PySimulation {
    topology: Topology,
    configuration: Configuration,
    md_sequence: AlgorithmSequence,
    sim_state: SimulationState,
    dt: f64,
    n_atoms: usize,
    /// Kinetic degrees of freedom (constraint- and NDFMIN-aware), computed once
    /// at build time by `gromos_run::total_dof` and reused for `temperature` —
    /// the same value the thermostat (if any) couples to.
    total_dof: f64,
    /// The run recipe this simulation was built from.
    recipe: RunRecipe,
    /// The algorithm plan that was instantiated (a frozen snapshot).
    plan: Vec<AlgorithmSpec>,
    /// What `from_imd` defaulted or passed through.
    diagnostics: Vec<String>,
}

#[pymethods]
impl PySimulation {
    /// Create a new simulation.
    ///
    /// Forms:
    ///   - `Simulation(system, recipe)` — `System` + `Recipe`
    ///   - `Simulation(topology, configuration, recipe)` — pre-loaded objects
    ///   - `Simulation("system.topo", "initial.cnf", "run.imd")` — file paths
    ///   - `plan=` — an edited `Plan` from `recipe.plan(system)` (validated, then run)
    ///
    /// Deprecated (one release, with a warning): an `InputParameters` in place of the
    /// recipe; the `distrest`/`posresspec`/`refpos` restraint-file keywords (use
    /// `recipe.with_inputs(...)`); the `ml_potential`/`ml_region`/`ml_buffer` keywords (use
    /// `recipe.with_term(Term("schnet", ...))`).
    #[new]
    #[pyo3(signature = (arg1, arg2, arg3=None, *, plan=None, distrest=None, posresspec=None, refpos=None, ml_potential=None, ml_region=None, ml_buffer=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        arg1: &Bound<'_, PyAny>,
        arg2: &Bound<'_, PyAny>,
        arg3: Option<&Bound<'_, PyAny>>,
        plan: Option<PyRef<'_, PyPlan>>,
        distrest: Option<String>,
        posresspec: Option<String>,
        refpos: Option<String>,
        ml_potential: Option<&Bound<'_, PyAny>>,
        ml_region: Option<String>,
        ml_buffer: Option<String>,
    ) -> PyResult<Self> {
        let (parts, mut recipe, diagnostics) = match arg3 {
            None => {
                let system = arg1.extract::<PyRef<PySystem>>().map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "With two arguments, first must be a System object",
                    )
                })?;
                let (recipe, diagnostics) = recipe_arg(arg2).map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                        "With two arguments, second must be a Recipe (or InputParameters)",
                    )
                })?;
                (Parts::from_system(&system), recipe, diagnostics)
            },
            Some(arg3) => {
                if let (Ok(topo_file), Ok(conf_file), Ok(input_file)) = (
                    arg1.extract::<String>(),
                    arg2.extract::<String>(),
                    arg3.extract::<String>(),
                ) {
                    let parts = Parts::from_files(&topo_file, &conf_file)?;
                    let recipe = PyRecipe::from_imd(&input_file, None)?;
                    (parts, recipe.inner, recipe.diagnostics)
                } else {
                    let topo = arg1.extract::<PyRef<PyTopology>>().map_err(|_| {
                        PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                            "First argument must be a file path (str) or Topology object",
                        )
                    })?;
                    let conf = arg2.extract::<PyRef<PyConfiguration>>().map_err(|_| {
                        PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                            "Second argument must be a file path (str) or Configuration object",
                        )
                    })?;
                    let (recipe, diagnostics) = recipe_arg(arg3).map_err(|_| {
                        PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                            "Third argument must be a file path (str), a Recipe, or \
                             InputParameters",
                        )
                    })?;
                    (Parts::from_objects(&topo, &conf), recipe, diagnostics)
                }
            },
        };
        apply_legacy_kwargs(
            py,
            &mut recipe,
            distrest,
            posresspec,
            refpos,
            ml_potential,
            ml_region,
            ml_buffer,
        )?;
        let plan = plan.map(|p| p.specs.clone());
        build_from_recipe(parts, recipe, plan, diagnostics)
    }

    /// Create a simulation from file paths (alternative to constructor).
    ///
    /// Example:
    ///     sim = Simulation.from_files("system.topo", "initial.cnf", "run.imd")
    #[staticmethod]
    fn from_files(topo_file: &str, conf_file: &str, input_file: &str) -> PyResult<Self> {
        let parts = Parts::from_files(topo_file, conf_file)?;
        let recipe = PyRecipe::from_imd(input_file, None)?;
        build_from_recipe(parts, recipe.inner, None, recipe.diagnostics)
    }

    /// Create a simulation from a run bundle (`input.toml` naming topology, coordinates and
    /// a recipe or parameter file — the reference systems' layout).
    #[staticmethod]
    #[pyo3(signature = (path, allow_passthrough=None))]
    fn from_bundle(path: &str, allow_passthrough: Option<Vec<String>>) -> PyResult<Self> {
        let policy = PassthroughPolicy::allow(allow_passthrough.unwrap_or_default());
        let (bundle, recipe, diag) = load_bundle(path, &policy).map_err(run_err)?;
        let parts = Parts::from_files(
            &bundle.topology.to_string_lossy(),
            &bundle.configuration.to_string_lossy(),
        )?;
        build_from_recipe(parts, recipe, None, diag.notes)
    }

    /// Run the simulation for `n_steps` MD steps.
    ///
    /// Example:
    ///     sim.step(1000)  # advance 1000 steps
    fn step(&mut self, n_steps: usize) -> PyResult<()> {
        for _ in 0..n_steps {
            self.md_sequence
                .run_step(&self.topology, &mut self.configuration, &self.sim_state)
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Error at step {}: {}",
                        self.sim_state.step, e
                    ))
                })?;
            self.sim_state.advance();
        }
        Ok(())
    }

    /// Run `n_steps` MD steps, sampling energies every `ene_freq` steps.
    ///
    /// Returns an `(n_frames, 13)` numpy array with columns
    /// `[time, kinetic, potential, total, volume, pressure, bond, angle, improper, dihedral, lj, coulomb]`
    /// — the same component order as the `.tre` file written by the `md` binary
    /// (see `gromos_io::energy::EnergyFrame`), so a `run()` array and a `.tre`
    /// dump of the same trajectory line up column-for-column.
    /// Frame 0 is the state before any of these steps ran; subsequent frames
    /// are sampled after every `ene_freq`-th step.
    ///
    /// Example:
    ///     energies = sim.run(1000, ene_freq=100)  # (11, 13) array
    #[pyo3(signature = (n_steps, ene_freq=100))]
    fn run<'py>(
        &mut self,
        py: Python<'py>,
        n_steps: usize,
        ene_freq: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let mut rows: Vec<Vec<f64>> = vec![self.energy_row()];
        for i in 0..n_steps {
            self.md_sequence
                .run_step(&self.topology, &mut self.configuration, &self.sim_state)
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Error at step {}: {}",
                        self.sim_state.step, e
                    ))
                })?;
            self.sim_state.advance();
            if (i + 1) % ene_freq == 0 {
                rows.push(self.energy_row());
            }
        }
        Ok(PyArray2::from_vec2_bound(py, &rows)?)
    }

    // -- State getters -------------------------------------------------------

    /// Current simulation time in picoseconds.
    #[getter]
    fn time(&self) -> f64 {
        self.sim_state.time
    }

    /// Current step number.
    #[getter]
    fn current_step(&self) -> usize {
        self.sim_state.step
    }

    /// Time step size in picoseconds.
    #[getter]
    fn dt(&self) -> f64 {
        self.dt
    }

    /// Number of atoms in the system.
    #[getter]
    fn n_atoms(&self) -> usize {
        self.n_atoms
    }

    /// Names of algorithms in the MD sequence.
    ///
    /// Example:
    ///     print(sim.algorithm_names)
    ///     # ['Forcefield', 'LeapFrogVelocity', 'LeapFrogPosition', ...]
    #[getter]
    fn algorithm_names(&self) -> Vec<String> {
        self.md_sequence
            .algorithm_names()
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    /// The recipe this simulation was built from.
    #[getter]
    fn recipe(&self) -> PyRecipe {
        PyRecipe {
            inner: self.recipe.clone(),
            diagnostics: self.diagnostics.clone(),
        }
    }

    /// The plan that was instantiated — a frozen snapshot, editing it does not change this
    /// simulation.
    #[getter]
    fn plan(&self) -> PyPlan {
        PyPlan {
            specs: self.plan.clone(),
        }
    }

    /// The effective run recipe as TOML — every value the engine actually used, the same
    /// text the `md` binary writes next to its `.tre`.
    #[getter]
    fn recipe_toml(&self) -> PyResult<String> {
        self.recipe.to_toml().map_err(run_err)
    }

    /// The algorithm plan as JSON (one entry per algorithm, fully resolved).
    #[getter]
    fn plan_json(&self) -> PyResult<String> {
        gromos_run::plan_to_json(&self.plan).map_err(run_err)
    }

    /// Notes from reading the parameters: which optional blocks were absent (and what that
    /// means), which blocks passed through. Empty when nothing was defaulted.
    #[getter]
    fn diagnostics(&self) -> Vec<String> {
        self.diagnostics.clone()
    }

    // -- Phase-space getters -------------------------------------------------

    /// Current positions as an Nx3 numpy array (nm).
    ///
    /// Returns positions from the "old" state, which holds the most
    /// recently completed step's data (GROMOS convention).
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = state.pos.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Current velocities as an Nx3 numpy array (nm/ps).
    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = state.vel.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Current forces as an Nx3 numpy array (kJ/mol/nm).
    ///
    /// For a truncated-octahedron box, GROMOS rotates FREEFORCERED output
    /// back into the original (cube) frame via
    /// `truncoct_triclinic_rotmat(false)` — positions/velocities stay in the
    /// triclinic frame, only forces are rotated back — mirrored here so this
    /// getter agrees with the `md` binary's `.trf` output (no ROTTRANS block
    /// support here, so `rmat(phi,theta,psi)` is identity, same as `md.rs`).
    #[getter]
    fn forces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let state = self.configuration.old();
        let data: Vec<Vec<f64>> = if state.box_config.box_type == BoxType::TruncatedOctahedral {
            let rot = truncoct_triclinic_rotmat(false);
            state
                .force
                .iter()
                .map(|f| {
                    let r = rot * *f;
                    vec![r.x, r.y, r.z]
                })
                .collect()
        } else {
            state.force.iter().map(|v| vec![v.x, v.y, v.z]).collect()
        };
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    // -- Energy getters ------------------------------------------------------

    /// Current energies (total, kinetic, potential, and components).
    #[getter]
    fn energies(&self) -> PyEnergy {
        let e = &self.configuration.old().energies;
        PyEnergy {
            total: e.total(),
            kinetic: e.kinetic_total,
            potential: e.potential_total,
            bond: e.bond_total,
            angle: e.angle_total,
            dihedral: e.dihedral_total,
            improper: e.improper_total,
            lj: e.lj_total,
            coulomb: e.crf_total,
        }
    }

    /// Total energy (kJ/mol).
    #[getter]
    fn total_energy(&self) -> f64 {
        self.configuration.old().energies.total()
    }

    /// Potential energy (kJ/mol).
    #[getter]
    fn potential_energy(&self) -> f64 {
        self.configuration.old().energies.potential_total
    }

    /// The energy of every additive term of the recipe (`Term("xtb", …)`, …) at the last
    /// step, keyed by the term's registry name (`"xtb"`, or `"xtb:1"` when a kind occurs more
    /// than once), in plan order. Their sum is what the terms add to `total_energy`;
    /// `potential_energy` is the classical force field alone. Empty without terms.
    #[getter]
    fn term_energies<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        let d = pyo3::types::PyDict::new_bound(py);
        for (name, energy) in &self.configuration.old().energies.term_energies {
            d.set_item(name, energy)?;
        }
        Ok(d)
    }

    /// dH/dλ (kJ/mol) at the last step — the thermodynamic-integration integrand: the whole
    /// derivative (kinetic, potential, constraints and restraints; gromosXX `energy.cc:597`),
    /// what the first value of a `.trg` frame holds. 0 unless the recipe's `perturbation` is
    /// enabled (a `pttopo` in `inputs`, `NTG=1`).
    #[getter]
    fn dhdl(&self) -> f64 {
        self.configuration
            .old()
            .perturbed_energy_derivatives
            .total()
    }

    /// dH/dλ split by term at the last step (kJ/mol): `lj`, `crf`, `bonded`, `kinetic` (the
    /// λ-dependent masses), `constraints` (perturbed SHAKE) and `distanceres` (perturbed
    /// distance restraints).
    #[getter]
    fn dhdl_terms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        let e = &self.configuration.old().perturbed_energy_derivatives;
        let d = pyo3::types::PyDict::new_bound(py);
        d.set_item("lj", e.lj_total)?;
        d.set_item("crf", e.crf_total)?;
        d.set_item(
            "bonded",
            e.bond_total + e.angle_total + e.improper_total + e.dihedral_total,
        )?;
        d.set_item("kinetic", e.kinetic_total)?;
        d.set_item("constraints", e.constraint_total)?;
        d.set_item("distanceres", e.distanceres_total)?;
        Ok(d)
    }

    /// EDS reference-state energy V_R (kJ/mol) at the last step (`.tre` slot 35, part of
    /// `total_energy` through the special energy, not of `potential_energy`); 0 without EDS.
    #[getter]
    fn eds_vr(&self) -> f64 {
        self.configuration.old().energies.eds_vr
    }

    /// EDS end-state energies V_i (kJ/mol) at the last step, one per state (the `.tre`'s
    /// `# eds` rows); empty without EDS.
    #[getter]
    fn eds_vi(&self) -> Vec<f64> {
        self.configuration.old().energies.eds_vi.clone()
    }

    /// Kinetic energy (kJ/mol).
    #[getter]
    fn kinetic_energy(&self) -> f64 {
        self.configuration.old().energies.kinetic_total
    }

    /// Current temperature (K), computed from kinetic energy.
    ///
    /// Uses the same constraint- and NDFMIN-aware degrees-of-freedom count the
    /// thermostat couples to (`total_dof`), not a bare `3*n_atoms` — a
    /// constrained system's true DOF is lower, and using the unconstrained
    /// count here used to silently disagree with what the thermostat targets.
    #[getter]
    fn temperature(&self) -> f64 {
        let n_dof = self.total_dof.round() as usize;
        self.configuration.old().temperature(n_dof)
    }

    /// Current box volume (nm³).
    #[getter]
    fn volume(&self) -> f64 {
        self.configuration.old().box_config.volume()
    }

    /// Current instantaneous pressure (kJ mol⁻¹ nm⁻³): trace / 3 of the pressure tensor
    /// `2/V·(E_kin + ½·Σ r⊗f)`, the number the barostat couples to and gromosXX writes.
    ///
    /// The tensor is only computed by `PressureCalculation` in the algorithm sequence —
    /// automatic for NPT (`Recipe.npt(...)`), absent for NVE/NVT, where this is 0.
    #[getter]
    fn pressure(&self) -> f64 {
        self.configuration.old().pressure()
    }

    // -- Topology getters ----------------------------------------------------

    /// Number of solute atoms.
    #[getter]
    fn n_solute_atoms(&self) -> usize {
        self.topology.num_solute_atoms()
    }

    /// Number of solvent atoms.
    #[getter]
    fn n_solvent_atoms(&self) -> usize {
        self.n_atoms - self.topology.num_solute_atoms()
    }

    /// Create a simulation with a custom *descriptor* sequence.
    ///
    /// Deprecated: use `plan = recipe.plan(system)`, edit it, and
    /// `Simulation(system, recipe, plan=plan)`.
    #[staticmethod]
    fn from_sequence(
        py: Python<'_>,
        topo: &PyTopology,
        conf: &PyConfiguration,
        params: &PyInputParameters,
        sequence: &PyPlan,
    ) -> PyResult<Self> {
        deprecated(
            py,
            "Simulation.from_sequence(topo, conf, params, sequence)",
            "Simulation(system, recipe, plan=recipe.plan(system))",
        )?;
        let (recipe, diagnostics) =
            RunRecipe::from_imd_with(&params.inner, &PassthroughPolicy::default())
                .map_err(run_err)?;
        build_from_recipe(
            Parts::from_objects(topo, conf),
            recipe,
            Some(sequence.specs.clone()),
            diagnostics.notes,
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "Simulation(n_atoms={}, step={}, time={:.3} ps, E_tot={:.6e} kJ/mol)",
            self.n_atoms,
            self.sim_state.step,
            self.sim_state.time,
            self.configuration.old().energies.total(),
        )
    }
}

// Private helpers
impl PySimulation {
    /// `[time, kinetic, potential, total, volume, pressure, bond, angle, improper, dihedral, lj,
    /// coulomb, dhdl]` for the current state — the `EnergyFrame` layout plus dH/dλ, used by
    /// `run()` to
    /// build the energy timeseries. Temperature is left at 0.0: unlike volume
    /// and pressure it needs the degrees-of-freedom count, which this row does
    /// not carry (use the `temperature` getter).
    ///
    /// The energies are those of the step just taken, evaluated at the positions that step
    /// started from, so the row is stamped with that step's time: `sim_state` has already
    /// advanced past it. gromosXX's `.tre` stamps step 0 with t = 0 the same way.
    fn energy_row(&self) -> Vec<f64> {
        let state = self.configuration.old();
        let volume = state.box_config.volume();
        let pressure = state.pressure();
        let frame = gromos_io::energy::EnergyFrame::from_energy(
            &state.energies,
            (self.sim_state.step - 1) as f64 * self.dt,
            0.0,
            volume,
            pressure,
        );
        vec![
            frame.time,
            frame.kinetic,
            frame.potential,
            frame.total,
            frame.volume,
            frame.pressure,
            frame.bond,
            frame.angle,
            frame.improper,
            frame.dihedral,
            frame.lj,
            frame.coul_real,
            state.perturbed_energy_derivatives.total(),
        ]
    }
}

/// Register simulation bindings in the Python module.
pub fn register_simulation(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySimulation>()?;
    Ok(())
}
