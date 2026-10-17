//! Python binding for ML potentials (P3.7) — connecting a trained model to a
//! `Simulation`'s zone definition by name, not by hand-counted atom index.
//!
//! `SchNetPotential` is deliberately a *recipe*, not an eagerly-loaded model: the real
//! `SchNetInteraction::load()` call needs `elements` validated against a region size that's only
//! known once a `Topology` is available (inside `Simulation`'s own construction, via the
//! `ml_potential=`/`ml_region=`/`ml_buffer=` kwargs — see `simulation.rs`), so this just holds
//! the path/cutoff/elements until then.
//!
//! `elements` is caller-supplied, not derived from `Topology` — `Topology` has no atomic-number
//! field (only `iac`, a force-field type index, and `mass`), and no such derivation exists
//! anywhere in this codebase yet; every real QM/ML provider this session (`XtbInteraction`,
//! `MopacInteraction`, the real `t_06` SchNet test) gets atomic numbers from an external source.
//! Building a `Topology` -> atomic-number inference is real, separate, speculative work.

use numpy::{PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::Topology;
use newen_forces::nonbonded::SchNetInteraction;
use newen_forces::orchestrator::ProviderOrchestrator;
use newen_forces::orchestrator_algorithm::ProviderOrchestratorAlgorithm;
use newen_forces::provider::PotentialProvider;
use newen_forces::zones::ZonePartition;

use crate::topology::PyTopology;

/// A trained SchNetPack 2 TorchScript model, to be attached to a `Simulation`'s QM/ML zone via
/// `ml_potential=`/`ml_region=`/`ml_buffer=` kwargs.
#[pyclass(name = "SchNetPotential")]
#[derive(Debug, Clone)]
pub struct PySchNetPotential {
    pub(crate) model_path: String,
    pub(crate) cutoff: f64,
    pub(crate) elements: Vec<i64>,
}

#[pymethods]
impl PySchNetPotential {
    /// Args:
    ///     model_path: path to a TorchScript-exported SchNetPack 2 model (see
    ///         `scripts/export_toy_schnet.py`/`scripts/train_qmmm_schnet.py`).
    ///     cutoff: radial cutoff (nm) for the model's own neighbor graph.
    ///     elements: atomic number per atom in the region this potential will be attached to
    ///         (caller-supplied — see module docs for why this isn't derived automatically).
    #[new]
    fn new(model_path: String, cutoff: f64, elements: Vec<i64>) -> Self {
        Self {
            model_path,
            cutoff,
            elements,
        }
    }

    /// Real (loaded) SchNetPack model energy (kJ/mol) and forces (kJ/mol/nm) for `positions`
    /// (nm, an Nx3 array) — a direct, standalone call (not wired into a `Simulation`), for
    /// comparison against `XtbPotential.evaluate()` on the same positions. Same isolated-cluster
    /// scope as `XtbPotential`/`qm_vs_ml_comparison.rs` — `Embedding::None`, vacuum, no PBC.
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

        let mut interaction = self.load()?;

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
            "SchNetPotential(model_path={:?}, cutoff={}, n_elements={})",
            self.model_path,
            self.cutoff,
            self.elements.len()
        )
    }
}

impl PySchNetPotential {
    /// Build the real, loaded `SchNetInteraction` this recipe describes.
    pub(crate) fn load(&self) -> PyResult<SchNetInteraction> {
        SchNetInteraction::load(&self.model_path, self.cutoff, self.elements.clone())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    }
}

/// Resolve an inner/buffer zone split by selector string against a real topology, for testing
/// against a known-good hand-built `ZonePartition` (e.g. `zone_partition_reference.rs`'s `t_06`
/// fixture) — the same resolution `Simulation`'s `ml_region=`/`ml_buffer=` kwargs use internally.
///
/// Returns `(inner_indices, buffer_indices, outer_indices)`.
#[pyfunction]
#[pyo3(signature = (topology, inner, buffer=None))]
pub fn resolve_zone_partition(
    topology: &PyTopology,
    inner: &str,
    buffer: Option<&str>,
) -> PyResult<(Vec<usize>, Vec<usize>, Vec<usize>)> {
    let partition = ZonePartition::from_selections(&topology.inner, inner, buffer)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
    Ok((
        partition.atoms_in(newen_forces::zones::Zone::Inner),
        partition.atoms_in(newen_forces::zones::Zone::Buffer),
        partition.atoms_in(newen_forces::zones::Zone::Outer),
    ))
}

pub fn register_ml_potential(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySchNetPotential>()?;
    m.add_function(wrap_pyfunction!(resolve_zone_partition, m)?)?;
    Ok(())
}
