//! Python binding for `XtbInteraction` — a real QM engine, directly callable from Python.
//!
//! Grew out of a direct follow-up to P3.7 (ML potential binding): `SchNetPotential` could be
//! attached to a `Simulation`, but there was no way to get a real QM reference value from
//! Python to compare it against — the QM-vs-ML comparison `qm_vs_ml_comparison.rs` does
//! (real `XtbInteraction` vs a trained `SchNetInteraction`, RMSE check) only existed in Rust.
//! `XtbPotential.evaluate()` closes that gap: a standalone, direct call (not wired into a
//! `Simulation`'s step loop — this is for comparison/reference use, not production QM/MM), same
//! `Embedding::None` isolated-cluster scope `qm_vs_ml_comparison.rs` uses.
//!
//! Not feature-gated — unlike `SchNetPotential`, `XtbInteraction` is a subprocess wrapper around
//! the real `xtb` binary, no `libtorch`/`--features ml` needed.

use numpy::{PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::Topology;
use newen_forces::nonbonded::XtbInteraction;
use newen_forces::provider::PotentialProvider;

/// A real `xtb` (GFN-xTB) QM engine, callable directly on a set of positions for comparison
/// against a trained `SchNetPotential` — see module docs.
#[pyclass(name = "XtbPotential")]
#[derive(Debug, Clone)]
pub struct PyXtbPotential {
    work_dir: String,
    gfn: u8,
    charge: i32,
    multiplicity: i32,
    elements: Vec<i64>,
}

#[pymethods]
impl PyXtbPotential {
    /// Args:
    ///     work_dir: scratch directory for xtb's input/output files (created if missing,
    ///         reused/overwritten on every `evaluate()` call).
    ///     elements: atomic number per atom.
    ///     gfn: 1 or 2 (GFN1-xTB / GFN2-xTB), default 2.
    ///     charge: net molecular charge, default 0.
    ///     multiplicity: spin multiplicity (1 = singlet), default 1.
    #[new]
    #[pyo3(signature = (work_dir, elements, gfn=2, charge=0, multiplicity=1))]
    fn new(work_dir: String, elements: Vec<i64>, gfn: u8, charge: i32, multiplicity: i32) -> Self {
        Self {
            work_dir,
            gfn,
            charge,
            multiplicity,
            elements,
        }
    }

    /// Real xtb single-point energy (kJ/mol) and forces (kJ/mol/nm) for `positions` (nm, an
    /// Nx3 array). The whole array is treated as one isolated QM cluster in vacuum
    /// (`Embedding::None`) — no environment/embedding support here, matching
    /// `qm_vs_ml_comparison.rs`'s own scope.
    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        positions: PyReadonlyArray2<'py, f64>,
    ) -> PyResult<(f64, Bound<'py, PyArray2<f64>>)> {
        let pos_arr = positions.as_array();
        if pos_arr.shape().len() != 2 || pos_arr.shape()[1] != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "positions must be an Nx3 array",
            ));
        }
        let n = pos_arr.shape()[0];
        if n != self.elements.len() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "positions has {n} atoms but elements has {}",
                self.elements.len()
            )));
        }

        let mut interaction = XtbInteraction::new(
            self.work_dir.clone(),
            self.gfn,
            self.charge,
            self.multiplicity,
            self.elements.clone(),
        )
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        let mut conf = Configuration::new(n, 1, 1);
        conf.current_mut().pos = (0..n)
            .map(|i| Vec3::new(pos_arr[[i, 0]], pos_arr[[i, 1]], pos_arr[[i, 2]]))
            .collect();

        let topo = Topology::new();
        let region = AtomSelection::all(n);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let contribution = interaction
            .contribute(&region, &topo, &conf, &index)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        let mut forces = vec![vec![0.0; 3]; n];
        for (i, f) in contribution.forces {
            forces[i] = vec![f.x, f.y, f.z];
        }

        Ok((contribution.energy, PyArray2::from_vec2_bound(py, &forces)?))
    }

    fn __repr__(&self) -> String {
        format!(
            "XtbPotential(gfn={}, charge={}, multiplicity={}, n_elements={})",
            self.gfn,
            self.charge,
            self.multiplicity,
            self.elements.len()
        )
    }
}

pub fn register_qm_potential(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyXtbPotential>()?;
    Ok(())
}
