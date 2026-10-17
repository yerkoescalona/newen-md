//! PySystem — paired Topology + Configuration with atom-count validation.

use numpy::PyArray2;
use pyo3::prelude::*;

use newen_io::g96::write_g96;

use crate::py_conf::PyConfiguration;
use crate::recipe::RecipeError;
use crate::topology::PyTopology;

/// A molecular system: topology paired with a matching configuration.
///
/// Validates at construction that both objects describe the same number of
/// atoms, mirroring the gromosXX invariant that `.top` and `.cnf` must agree. A solute-only
/// topology with a configuration of solute plus whole solvent molecules is solvated here.
///
/// # Example (Python)
///
/// ```python
/// system = System.from_files("system.topo", "initial.cnf")
/// print(system.n_atoms)          # 648
/// print(system.charge)           # 0
/// print(system.positions.shape)  # (648, 3)
/// ```
#[pyclass(name = "System", module = "newen")]
#[derive(Debug)]
pub struct PySystem {
    pub(crate) topology: PyTopology,
    pub(crate) configuration: PyConfiguration,
    /// The (topology, coordinates) files of a `System.from_files` system, absolute: how it
    /// pickles (a worker process reads the same files). `None` for a system built in memory.
    source: Option<(String, String)>,
}

#[pymethods]
impl PySystem {
    /// Create a System from pre-loaded Topology and Configuration objects.
    ///
    /// Solvates a solute-only topology to the configuration's size; raises `RecipeError` (a
    /// `ValueError`) if the atom counts cannot be made to match.
    #[new]
    fn new(
        topology: &Bound<'_, PyTopology>,
        configuration: &Bound<'_, PyConfiguration>,
    ) -> PyResult<Self> {
        let topo_ref = topology.borrow();
        let conf_ref = configuration.borrow();
        Self::paired(
            PyTopology {
                inner: topo_ref.inner.clone(),
                physical_constants: topo_ref.physical_constants,
                parsed: topo_ref.parsed.clone(),
            },
            PyConfiguration {
                labels: conf_ref.labels.clone(),
                pos_data: conf_ref.pos_data.clone(),
                vel_data: conf_ref.vel_data.clone(),
                box_dims: conf_ref.box_dims,
                pertdata: conf_ref.pertdata,
                aedssearch: conf_ref.aedssearch.clone(),
            },
            None,
        )
    }

    /// Load a System directly from topology and coordinate file paths.
    #[staticmethod]
    fn from_files(topo_file: &str, conf_file: &str) -> PyResult<Self> {
        let topo = PyTopology::from_file(topo_file)?;
        let conf = PyConfiguration::from_file(conf_file)?;
        let absolute = |p: &str| {
            std::path::absolute(p).map_or_else(|_| p.to_string(), |a| a.display().to_string())
        };
        Self::paired(topo, conf, Some((absolute(topo_file), absolute(conf_file))))
    }

    /// Pickling: a `System.from_files` system is `System.from_files` on the same (absolute)
    /// paths, which is how `newen.free_energy.run_windows(parallel=N)` hands it to its worker
    /// processes. A system built in memory has no such form and does not pickle.
    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, (String, String))> {
        let (topo, conf) = self.source.clone().ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "only a System made by System.from_files pickles (as its two files); this one \
                 was built in memory",
            )
        })?;
        let from_files = py.get_type_bound::<Self>().getattr("from_files")?;
        Ok((from_files, (topo, conf)))
    }

    // ── Properties ────────────────────────────────────────────────────────────

    /// Total number of atoms.
    #[getter]
    fn n_atoms(&self) -> usize {
        self.topology.inner.num_atoms()
    }

    /// Total integer charge of the system (e).
    #[getter]
    fn charge(&self) -> i32 {
        self.topology.charge()
    }

    /// Positions as an Nx3 numpy array (nm).
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.configuration.positions(py)
    }

    /// Velocities as an Nx3 numpy array (nm/ps).
    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.configuration.velocities(py)
    }

    /// Simulation box dimensions (Lx, Ly, Lz) in nm.
    #[getter]
    fn r#box(&self) -> (f64, f64, f64) {
        self.configuration.box_dimensions()
    }

    /// The underlying Topology object.
    #[getter]
    fn topology(&self) -> PyTopology {
        PyTopology {
            inner: self.topology.inner.clone(),
            physical_constants: self.topology.physical_constants,
            parsed: self.topology.parsed.clone(),
        }
    }

    /// The underlying Configuration object.
    #[getter]
    fn configuration(&self) -> PyConfiguration {
        PyConfiguration {
            labels: self.configuration.labels.clone(),
            pos_data: self.configuration.pos_data.clone(),
            vel_data: self.configuration.vel_data.clone(),
            box_dims: self.configuration.box_dims,
            pertdata: self.configuration.pertdata,
            aedssearch: self.configuration.aedssearch.clone(),
        }
    }

    /// Write the configuration to a GROMOS coordinate file (.cnf).
    fn write(&self, path: &str) -> PyResult<()> {
        let pos = &self.configuration.pos_data;
        let vel = &self.configuration.vel_data;
        let box_v = self.configuration.box_dims;

        let vels = if vel.iter().any(|v| *v != newen_core::math::Vec3::ZERO) {
            Some(vel.as_slice())
        } else {
            None
        };
        let box_opt = if box_v.x > 0.0 { Some(box_v) } else { None };

        write_g96(
            path,
            "newen",
            pos,
            vels,
            box_opt,
            Some(&self.topology.inner),
        )
        .map_err(PyErr::new::<pyo3::exceptions::PyIOError, _>)
    }

    // ── Preparation verbs ───────────────────────────────────────────
    //
    // Sugar over `newen.tools`: each calls the same `newen_tools::programs` verb the module
    // function calls, and returns a *new* System — an input system is never mutated, so a
    // preparation chain is re-runnable.

    /// Solvate this system in a box of pre-equilibrated solvent (`newen.tools.sim_box`).
    ///
    /// The topology follows the coordinates: it is extended by the solvent molecules that were
    /// kept, so the result is a System, not a loose configuration.
    #[pyo3(signature = (solvent, *, minwall = None, thresh = 0.23, boxsize = false))]
    fn solvate(
        &self,
        solvent: &PyConfiguration,
        minwall: Option<Vec<f64>>,
        thresh: f64,
        boxsize: bool,
    ) -> PyResult<Self> {
        use newen_tools::programs::sim_box;

        let opts = sim_box::SimBoxOptions {
            minwall: minwall.unwrap_or_default(),
            thresh,
            use_boxsize: boxsize,
        };
        let solute_box = (self.configuration.box_dims != newen_core::math::Vec3::ZERO)
            .then_some(self.configuration.box_dims);
        let result = sim_box::solvate(
            &self.configuration.atoms_or_placeholders(),
            &self.topology.parsed.masses,
            solute_box,
            &solvent.atoms_or_placeholders(),
            solvent.box_dims,
            &opts,
        )
        .map_err(crate::programs::ProgramError::new_err)?;

        let mut topology = PyTopology::from_parsed(self.topology.parsed.clone());
        topology.inner.solvate(result.n_solvent_molecules);
        Ok(Self {
            topology,
            configuration: PyConfiguration::from_atoms(result.atoms, result.box_dims),
            source: None,
        })
    }

    /// Replace solvent molecules by counter-ions (`newen.tools.ion`), topology included.
    ///
    /// gromos++'s `ion` changes only the coordinates and leaves the topology to the caller; this
    /// does both halves — it builds each ion's topology from `forcefield` (an `(mtb, ifp)` pair,
    /// the same files `make_top` takes) and combines it with the solute's, so the result is a
    /// System whose two halves agree.
    #[pyo3(signature = (*, forcefield, positive = None, negative = None, potential = None, mindist = 0.0))]
    fn add_ions(
        &self,
        forcefield: (String, String),
        positive: Option<(usize, String)>,
        negative: Option<(usize, String)>,
        potential: Option<f64>,
        mindist: f64,
    ) -> PyResult<Self> {
        use newen_tools::programs::{com_top, ion, make_top};
        let err = crate::programs::ProgramError::new_err;

        let opts = ion::IonOptions {
            positive: positive.clone(),
            negative: negative.clone(),
            potential,
            mindist,
        };
        let data = newen_io::LabeledCoordinateData {
            atoms: self.configuration.atoms_or_placeholders(),
            box_dims: self.configuration.box_dims,
            box_type: 1,
        };
        let result = ion::add_ions(&self.topology.parsed, &data, &opts).map_err(err)?;

        // Each ion needs its own topology entry, in the order `ion` placed them: positive first.
        let (mtb_path, ifp_path) = forcefield;
        let bb = newen_io::mtb::read_mtb_file(&mtb_path)
            .map_err(|e| err(format!("Error reading MTB file '{mtb_path}': {e:?}")))?;
        let ffp = newen_io::ifp::read_ifp_file(&ifp_path)
            .map_err(|e| err(format!("Error reading IFP file '{ifp_path}': {e:?}")))?;

        let mut ion_tops = Vec::new();
        for (count, name) in [positive, negative].into_iter().flatten() {
            if count == 0 {
                continue;
            }
            let built = make_top::make_top(&bb, &ffp, std::slice::from_ref(&name), None)
                .map_err(|e| err(format!("building the topology of '{name}': {e}")))?;
            ion_tops.push((count, built.topology));
        }

        let mut parts: Vec<(usize, &newen_io::topology::ParsedTopology)> =
            vec![(1, &self.topology.parsed)];
        parts.extend(ion_tops.iter().map(|(n, t)| (*n, t)));
        let combined = com_top::com_top(&parts, 1, 1).map_err(err)?;

        let mut topology = PyTopology::from_parsed(combined);
        let n_solvent_atoms = result.atoms.len() - topology.inner.num_solute_atoms();
        let per = topology.parsed.solvent_atoms.len().max(1);
        topology.inner.solvate(n_solvent_atoms / per);
        Ok(Self {
            topology,
            configuration: PyConfiguration::from_atoms(result.atoms, result.box_dims),
            source: None,
        })
    }

    /// Run every consistency check over this system's topology (`newen.tools.check_top`).
    ///
    /// Returns `(errors, warnings, report)`; `newen.tools.check_top` wraps it as a
    /// `CheckReport`.
    fn check(&self) -> (Vec<String>, Vec<String>, String) {
        // check_top reads a topology file, solute only; the System's own may be solvated.
        let as_file = newen_io::topology::build_topology(self.topology.parsed.clone());
        let report = newen_tools::programs::check_top::check_top(&as_file);
        (report.errors, report.warnings, report.text)
    }

    fn __repr__(&self) -> String {
        let (lx, ly, lz) = self.configuration.box_dimensions();
        format!(
            "System(n_atoms={}, charge={}, box=({:.3}, {:.3}, {:.3}))",
            self.n_atoms(),
            self.charge(),
            lx,
            ly,
            lz,
        )
    }
}

impl PySystem {
    /// Pair a topology with its configuration. A solute-only topology with a configuration of
    /// solute plus whole solvent molecules is solvated here, by the routine `prepare_system`
    /// uses (`Topology::solvate`), so the System's atom counts agree from the start.
    fn paired(
        mut topology: PyTopology,
        configuration: PyConfiguration,
        source: Option<(String, String)>,
    ) -> PyResult<Self> {
        let n_conf = configuration.pos_data.len();
        validate_atom_count_match(&topology.inner, n_conf).map_err(RecipeError::new_err)?;
        if topology.inner.num_atoms() != n_conf {
            newen_io::topology::solvate_to_atoms(&mut topology.inner, n_conf)
                .map_err(|e| RecipeError::new_err(e.to_string()))?;
        }
        Ok(Self {
            topology,
            configuration,
            source,
        })
    }
}

pub fn register_system(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySystem>()?;
    Ok(())
}

/// Rust-only validation logic extracted for unit testing without pyo3 linker symbols.
/// The configuration must hold the topology's atoms — or, for a topology that is not yet
/// solvated, the solute plus a whole number of solvent molecules: `prepare_system` solvates
/// from the coordinate count (the `md` binary's rule; the recipe's NSM is only a hint), so a
/// `System` of a solute topology and a solvated coordinate file is a valid one.
pub(crate) fn validate_atom_count_match(
    topo: &newen_core::Topology,
    n_conf: usize,
) -> Result<(), String> {
    atom_count_ok(
        topo.num_atoms(),
        topo.num_solute_atoms(),
        topo.solvent_atom_template.len(),
        n_conf,
    )
}

fn atom_count_ok(
    n_topo: usize,
    n_solute: usize,
    per_solvent: usize,
    n_conf: usize,
) -> Result<(), String> {
    if n_topo == n_conf {
        return Ok(());
    }
    let unsolvated = n_topo == n_solute && per_solvent > 0;
    if unsolvated && n_conf >= n_solute && (n_conf - n_solute).is_multiple_of(per_solvent) {
        return Ok(());
    }
    Err(if unsolvated {
        format!(
            "Topology has {n_solute} solute atoms and {per_solvent}-atom solvent molecules, but \
             the configuration has {n_conf} atoms — not solute plus whole solvent molecules"
        )
    } else {
        format!("Topology has {n_topo} atoms but configuration has {n_conf} — they must match")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_count_match_ok() {
        assert!(atom_count_ok(648, 0, 3, 648).is_ok());
    }

    #[test]
    fn atom_count_mismatch_error() {
        let err = atom_count_ok(648, 0, 3, 3).unwrap_err();
        assert!(err.contains("must match"), "message: {err}");
    }

    #[test]
    fn unsolvated_topology_accepts_solute_plus_whole_solvent_molecules() {
        // nacl_1water_box: 2 solute atoms, 3-atom SPC solvent, 5 atoms in the coordinates.
        assert!(atom_count_ok(2, 2, 3, 5).is_ok());
        assert!(atom_count_ok(2, 2, 3, 2).is_ok());
        let err = atom_count_ok(2, 2, 3, 6).unwrap_err();
        assert!(err.contains("whole solvent molecules"), "message: {err}");
        // A solvated topology (all 648 atoms present) is still an exact match only.
        assert!(atom_count_ok(648, 0, 3, 651).is_err());
    }
}
