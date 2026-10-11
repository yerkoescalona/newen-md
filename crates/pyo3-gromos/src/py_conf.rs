//! Python wrapper for GROMOS configuration (system state).

use numpy::PyArray2;
use pyo3::prelude::*;

use gromos_core::math::Vec3;
use gromos_io::coordinate::{format_g96, parse_g96_labeled, read_coordinates, read_g96_labeled};
use gromos_io::G96Atom;

/// System configuration: positions, velocities, and simulation box.
///
/// Loaded from a GROMOS coordinate file (.cnf/.g96). Represents the
/// initial state of the system before simulation.
///
/// # Example (Python)
///
/// ```python
/// conf = Configuration("initial.cnf")
/// print(conf.n_atoms)
/// print(conf.positions.shape)      # (N, 3) numpy array
/// print(conf.box_dimensions)       # (Lx, Ly, Lz) in nm
/// ```
#[pyclass(name = "Configuration")]
#[derive(Debug)]
pub struct PyConfiguration {
    pub(crate) pos_data: Vec<Vec3>,
    pub(crate) vel_data: Vec<Vec3>,
    pub(crate) box_dims: Vec3,
    /// Residue and atom names, when the source had them. Every coordinate tool needs them
    /// (`sim_box` labels its solvent `SOLV`, `ion` names the ion it places), so a `Configuration`
    /// that came from a `.cnf` or from `gromos.tools` carries them along.
    pub(crate) labels: Option<Vec<G96Atom>>,
    /// λ of the file's PERTDATA block (a run with NRDGL=1 starts from it).
    pub(crate) pertdata: Option<f64>,
}

#[pymethods]
impl PyConfiguration {
    /// Load a configuration from a GROMOS coordinate file.
    ///
    /// Args:
    ///     conf_file: Path to coordinate file (.cnf, .g96)
    #[new]
    fn new(conf_file: &str) -> PyResult<Self> {
        Self::from_file(conf_file)
    }

    /// Load a configuration from a GROMOS coordinate file.
    #[staticmethod]
    pub fn from_file(conf_file: &str) -> PyResult<Self> {
        let coord_data = read_coordinates(conf_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read coordinates '{}': {}",
                conf_file, e
            ))
        })?;
        // A second pass for the labels: `read_coordinates` gives the velocities, `read_g96_labeled`
        // the residue/atom names. Both are needed and neither yields the other.
        let labels = read_g96_labeled(conf_file).ok().map(|d| d.atoms);
        Ok(Self {
            pos_data: coord_data.positions,
            vel_data: coord_data.velocities,
            box_dims: coord_data.box_dims,
            labels,
            pertdata: coord_data.pertdata,
        })
    }

    /// Read a configuration from the text of a GROMOS coordinate file.
    ///
    /// What `str(configuration)` writes, this reads — so coordinates a program printed need no
    /// file. Velocities are not part of the G96 POSITION block and come back empty.
    #[staticmethod]
    pub fn from_string(text: &str) -> PyResult<Self> {
        let data = parse_g96_labeled(text.as_bytes()).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Failed to parse coordinates: {e}"
            ))
        })?;
        Ok(Self::from_atoms(data.atoms, data.box_dims))
    }

    /// The configuration as a GROMOS96 coordinate file.
    fn __str__(&self) -> String {
        self.to_g96(None)
    }

    /// The configuration as GROMOS96 text, with an optional TITLE.
    #[pyo3(signature = (title = None))]
    pub fn to_g96(&self, title: Option<&str>) -> String {
        let box_dims = (self.box_dims != Vec3::ZERO).then_some(self.box_dims);
        format_g96(
            title.unwrap_or("gromos-rs"),
            &self.atoms_or_placeholders(),
            box_dims,
        )
    }

    /// Number of atoms with coordinates.
    #[getter]
    fn n_atoms(&self) -> usize {
        self.pos_data.len()
    }

    /// Positions as an Nx3 numpy array (nm).
    #[getter]
    pub fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let data: Vec<Vec<f64>> = self.pos_data.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Velocities as an Nx3 numpy array (nm/ps).
    #[getter]
    pub fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let data: Vec<Vec<f64>> = self.vel_data.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Simulation box dimensions (Lx, Ly, Lz) in nm.
    #[getter]
    pub fn box_dimensions(&self) -> (f64, f64, f64) {
        (self.box_dims.x, self.box_dims.y, self.box_dims.z)
    }

    fn __repr__(&self) -> String {
        format!(
            "Configuration(n_atoms={}, box=({:.3}, {:.3}, {:.3}))",
            self.pos_data.len(),
            self.box_dims.x,
            self.box_dims.y,
            self.box_dims.z,
        )
    }
}

impl PyConfiguration {
    /// Wrap labelled atoms (what every coordinate tool returns).
    pub(crate) fn from_atoms(atoms: Vec<G96Atom>, box_dims: Vec3) -> Self {
        Self {
            pos_data: atoms.iter().map(|a| a.pos).collect(),
            vel_data: Vec::new(),
            box_dims,
            labels: Some(atoms),
            pertdata: None,
        }
    }

    /// The atoms with their labels; positions alone get placeholder names, so a configuration
    /// built from arrays still serialises rather than failing.
    pub(crate) fn atoms_or_placeholders(&self) -> Vec<G96Atom> {
        match &self.labels {
            Some(atoms) if atoms.len() == self.pos_data.len() => atoms
                .iter()
                .zip(&self.pos_data)
                .map(|(a, &pos)| G96Atom { pos, ..a.clone() })
                .collect(),
            _ => self
                .pos_data
                .iter()
                .enumerate()
                .map(|(i, &pos)| G96Atom {
                    res_num: i + 1,
                    res_name: "UNK".to_string(),
                    atom_name: "X".to_string(),
                    atom_num: i + 1,
                    pos,
                })
                .collect(),
        }
    }
}

pub fn register_configuration(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyConfiguration>()?;
    Ok(())
}
