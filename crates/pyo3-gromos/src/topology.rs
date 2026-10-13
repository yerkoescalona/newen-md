//! Python wrapper for GROMOS molecular topology.

use numpy::PyArray1;
use pyo3::prelude::*;

use gromos_core::{units::PhysicalConstants, Topology};
use gromos_io::topology::{
    build_topology, parse_topology, read_topology_file, write_parsed_topology, ParsedTopology,
};

/// Molecular topology loaded from a GROMOS topology file.
///
/// Contains atom types, masses, charges, bonded parameters, LJ parameters,
/// exclusions, and molecule definitions.
///
/// # Example (Python)
///
/// ```python
/// topo = Topology("system.topo")
/// print(topo.n_atoms)          # solute atoms (before solvation)
/// print(topo.n_solute_atoms)
/// print(topo.masses[:5])       # first 5 masses
/// topo.solvate(216)            # add 216 solvent molecules
/// print(topo.n_atoms)          # now includes solvent
/// ```
#[pyclass(name = "Topology")]
#[derive(Debug)]
pub struct PyTopology {
    pub(crate) inner: Topology,
    pub(crate) physical_constants: PhysicalConstants,
    /// The topology as its file blocks, kept beside the built form: `gromos.tools` works at block
    /// level (`com_top`, `red_top`, `ion` all take a `ParsedTopology`) and `build_topology`
    /// consumes it with no way back.
    pub(crate) parsed: ParsedTopology,
}

#[pymethods]
impl PyTopology {
    /// Load a topology from a GROMOS topology file.
    ///
    /// Args:
    ///     topo_file: Path to topology file (.topo, .top)
    #[new]
    fn new(topo_file: &str) -> PyResult<Self> {
        Self::from_file(topo_file)
    }

    /// Load a topology from a GROMOS topology file.
    #[staticmethod]
    pub fn from_file(topo_file: &str) -> PyResult<Self> {
        let parsed = read_topology_file(topo_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read topology '{}': {}",
                topo_file, e
            ))
        })?;
        Ok(Self::from_parsed(parsed))
    }

    /// Read a topology from the text of a GROMOS topology file.
    ///
    /// What `str(topology)` writes, this reads — so a topology a program printed needs no file.
    #[staticmethod]
    pub fn from_string(text: &str) -> PyResult<Self> {
        let parsed = parse_topology(text.as_bytes()).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Failed to parse topology: {e}"
            ))
        })?;
        Ok(Self::from_parsed(parsed))
    }

    /// The topology as a GROMOS topology file (gromos++ `OutTopology` layout).
    fn __str__(&self) -> String {
        write_parsed_topology(&self.parsed, "gromos-rs")
    }

    /// Write the topology to a GROMOS topology file.
    fn write(&self, path: &str, title: Option<&str>) -> PyResult<()> {
        std::fs::write(
            path,
            write_parsed_topology(&self.parsed, title.unwrap_or("gromos-rs")),
        )
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to write '{path}': {e}"))
        })
    }

    /// Total number of atoms (solute + solvent).
    #[getter]
    fn n_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    /// Number of solute atoms.
    #[getter]
    fn n_solute_atoms(&self) -> usize {
        self.inner.num_solute_atoms()
    }

    /// Number of solvent atoms.
    #[getter]
    fn n_solvent_atoms(&self) -> usize {
        self.inner.num_atoms() - self.inner.num_solute_atoms()
    }

    /// Atomic masses as a numpy array.
    #[getter]
    fn masses<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec_bound(py, self.inner.mass.clone())
    }

    /// Atomic charges as a numpy array.
    #[getter]
    fn charges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec_bound(py, self.inner.charge.clone())
    }

    /// Total integer charge of the system (e).
    #[getter]
    pub fn charge(&self) -> i32 {
        self.inner.charge.iter().sum::<f64>().round() as i32
    }

    /// Add solvent molecules to the topology.
    ///
    /// Args:
    ///     nsm: Number of solvent molecules to add
    fn solvate(&mut self, nsm: usize) {
        self.inner.solvate(nsm);
    }

    /// Merge a GROMOS perturbation topology (`.ptp`) for free-energy runs (FEP/TI) and EDS.
    ///
    /// Exactly what the `md` binary does with `@pttopo`: the perturbed solute terms are stored
    /// and the bonds/angles/dihedrals they perturb are removed from the regular solute, and an
    /// `MPERTATOM` block's EDS end states are merged. Call before `solvate()` and before
    /// building a `Simulation`; the run's `InputParameters` must have `NTG != 0` for the
    /// perturbed terms and an `EDS` block with `EDS = 1` for the end states (a topology
    /// carrying either without them is rejected at construction).
    ///
    /// Args:
    ///     pttopo_file: Path to the perturbation topology (.ptp)
    fn apply_perturbation(&mut self, pttopo_file: &str) -> PyResult<()> {
        let pt = gromos_io::ptp::read_pttopo(pttopo_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read perturbation topology '{}': {}",
                pttopo_file, e
            ))
        })?;
        self.inner
            .apply_perturbation(pt.perturbed_solute, pt.is_perturbed);
        if let Some(m) = pt.eds {
            self.inner
                .apply_eds_perturbation(m.num_states, m.state_names, m.atoms)
                .map_err(PyErr::new::<pyo3::exceptions::PyValueError, _>)?;
        }
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!(
            "Topology(n_atoms={}, n_solute={}, n_solvent={})",
            self.inner.num_atoms(),
            self.inner.num_solute_atoms(),
            self.inner.num_atoms() - self.inner.num_solute_atoms(),
        )
    }

    fn __len__(&self) -> usize {
        self.inner.num_atoms()
    }
}

impl PyTopology {
    /// Wrap a parsed topology, keeping both forms (`gromos.tools` produces these).
    pub(crate) fn from_parsed(parsed: ParsedTopology) -> Self {
        let physical_constants = parsed.physical_constants;
        let inner = build_topology(parsed.clone());
        Self {
            inner,
            physical_constants,
            parsed,
        }
    }
}

pub fn register_topology(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTopology>()?;
    Ok(())
}
