//! Python wrapper for the library-driven energy-trajectory reader.
//!
//! `Simulation.run()` hands back the scalar totals of a run in progress (`EnergyTimeseries`).
//! This reads a finished `.tre`/`.trg` the way `ene_ana` does — through the energy library, so
//! the file's layout is established before a value is believed (`EnergyTraj::bind`, profile
//! tiers in `docs/src/reference/energy-library.md`) — and keeps the *shape* the library
//! declares: per-bath and per-energy-group tables, not just the totals.
//!
//! A subblock comes back as `(n_frames, rows, cols)`; `NONBONDED` has one row per
//! energy-group pair, in gromosXX's order, which `group_pairs` names.

use newen_io::energy_traj::{
    group_pair, read_preamble, EnergyTraj, OFFICIAL_LIBRARY, PROFILE_VERSION,
};
use newen_io::gz::TextReader;
use numpy::{PyArray1, PyArray2, PyArray3, PyArrayMethods};
use pyo3::exceptions::{PyIOError, PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::BTreeMap;

/// One subblock over every frame: `values[frame][row][col]`.
struct Table {
    rows: usize,
    cols: usize,
    values: Vec<f64>,
}

/// An energy trajectory read into memory, as the library declares it.
///
/// # Example (Python)
///
/// ```python
/// tre = EnergyTrajectory("md.tre")
/// tre.times                     # (n_frames,)
/// tre["totene"]                 # any library variable, or "ENER[1]"
/// tre.table("NONBONDED")        # (n_frames, n_pairs, 6)
/// tre.group_pairs               # [(0, 0), (0, 1), ...] — the rows of that table
/// tre.to_long("NONBONDED")      # tidy columns, ready for a group_by
/// ```
#[pyclass(name = "EnergyTrajectory")]
pub struct PyEnergyTrajectory {
    path: String,
    file_type: String,
    version: Option<String>,
    warnings: Vec<String>,
    times: Vec<f64>,
    steps: Vec<f64>,
    tables: BTreeMap<String, Table>,
    /// Library variables evaluated per frame, on request.
    properties: BTreeMap<String, Vec<f64>>,
    n_energy_groups: usize,
    n_baths: usize,
}

fn io_err(e: impl std::fmt::Display) -> PyErr {
    PyErr::new::<PyIOError, _>(e.to_string())
}

#[pymethods]
impl PyEnergyTrajectory {
    /// Read a `.tre` (or a `.trg` with `free_energy=True`).
    ///
    /// Args:
    ///     path: the trajectory, plain or gzipped.
    ///     library: an `ene_ana` library file; the built-in one when omitted.
    ///     properties: library variables to evaluate per frame (`"totene"`, `"ENER[1]"`, ...).
    ///     free_energy: read the `FRENERTRJ` layout of a `.trg`.
    #[new]
    #[pyo3(signature = (path, library=None, properties=None, free_energy=false))]
    fn new(
        path: &str,
        library: Option<&str>,
        properties: Option<Vec<String>>,
        free_energy: bool,
    ) -> PyResult<Self> {
        let lib_text = match library {
            Some(p) => std::fs::read_to_string(p).map_err(|e| io_err(format!("{p}: {e}")))?,
            None => OFFICIAL_LIBRARY.to_string(),
        };
        let lib_name = library.unwrap_or("<built-in library>");
        let mut etrj = EnergyTraj::from_library(&lib_text)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("{lib_name}: {e}")))?;
        // MASS/NUMMOL come from a topology in ene_ana; without one the variables that use
        // them evaluate to 0 rather than failing the whole read.
        etrj.add_constant("MASS", 0.0);
        etrj.add_constant("NUMMOL", 0.0);

        let file_type = if free_energy { "FRENERTRJ" } else { "ENERTRJ" };
        let mut reader = TextReader::open(path).map_err(|e| io_err(format!("{path}: {e}")))?;
        let preamble = read_preamble(&mut reader)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("{path}: {e}")))?;
        let version = preamble.version.clone();
        etrj.bind(file_type, &preamble, lib_name, path)
            .map_err(PyErr::new::<PyValueError, _>)?;

        let props = properties.unwrap_or_default();
        let mut out = PyEnergyTrajectory {
            path: path.to_string(),
            file_type: file_type.to_string(),
            version,
            warnings: etrj.drain_warnings(),
            times: Vec::new(),
            steps: Vec::new(),
            tables: BTreeMap::new(),
            properties: props.iter().map(|p| (p.clone(), Vec::new())).collect(),
            n_energy_groups: 0,
            n_baths: 0,
        };

        while etrj
            .read_frame(&mut reader, file_type)
            .map_err(|e| PyErr::new::<PyValueError, _>(format!("{path}: {e}")))?
        {
            out.warnings.extend(etrj.drain_warnings());
            out.times.push(etrj.value("TIME[2]").unwrap_or(f64::NAN));
            out.steps.push(etrj.value("TIME[1]").unwrap_or(f64::NAN));
            for (name, series) in out.properties.iter_mut() {
                let v = etrj
                    .value(name)
                    .map_err(|e| PyErr::new::<PyKeyError, _>(format!("{path}: {e}")))?;
                series.push(v);
            }
            for name in etrj.subblock_names() {
                // a library declares both file types; the other one's subblocks stay empty,
                // as does one whose `size` is 0 in this run
                let Some(rows) = etrj.table(name).filter(|r| !r.is_empty()) else {
                    continue;
                };
                let (nr, nc) = (rows.len(), rows.first().map_or(0, |r| r.len()));
                let t = out.tables.entry(name.to_string()).or_insert(Table {
                    rows: nr,
                    cols: nc,
                    values: Vec::new(),
                });
                // a size that changed between frames is already an error in the reader
                for row in rows {
                    t.values.extend_from_slice(row);
                }
            }
            out.n_energy_groups = etrj.size("NUM_ENERGY_GROUPS").unwrap_or(0);
            out.n_baths = etrj.size("NUM_BATHS").unwrap_or(0);
        }
        Ok(out)
    }

    /// Number of frames read.
    fn __len__(&self) -> usize {
        self.times.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "EnergyTrajectory({:?}, {} frames, {}, ENEVERSION {}, {} energy groups)",
            self.path,
            self.times.len(),
            self.file_type,
            self.version.as_deref().unwrap_or("NONE"),
            self.n_energy_groups,
        )
    }

    /// A library variable or a direct `SUBBLOCK[i][j]` reference, over every frame.
    fn __getitem__<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let series = self.properties.get(name).ok_or_else(|| {
            PyErr::new::<PyKeyError, _>(format!(
                "{name:?} was not requested; pass it in `properties=` when opening {:?}",
                self.path
            ))
        })?;
        Ok(PyArray1::from_slice_bound(py, series))
    }

    /// The simulated time of each frame, ps.
    #[getter]
    fn times<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice_bound(py, &self.times)
    }

    /// The step number of each frame.
    #[getter]
    fn steps<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice_bound(py, &self.steps)
    }

    /// What the layout check had to say — empty when the file established its own layout.
    #[getter]
    fn warnings(&self) -> Vec<String> {
        self.warnings.clone()
    }

    /// The `ENEVERSION` the file carries, if any.
    #[getter]
    fn version(&self) -> Option<String> {
        self.version.clone()
    }

    /// The profile version this reader implements.
    #[getter]
    fn profile_version(&self) -> u32 {
        PROFILE_VERSION
    }

    #[getter]
    fn n_energy_groups(&self) -> usize {
        self.n_energy_groups
    }

    #[getter]
    fn n_baths(&self) -> usize {
        self.n_baths
    }

    /// The subblock names the library declares.
    #[getter]
    fn subblocks(&self) -> Vec<String> {
        self.tables.keys().cloned().collect()
    }

    /// The energy-group pair each row of a `matrix_`-sized subblock (`NONBONDED`) belongs to,
    /// 0-based, in gromosXX's write order.
    #[getter]
    fn group_pairs(&self) -> Vec<(usize, usize)> {
        (0..)
            .map_while(|r| group_pair(r, self.n_energy_groups))
            .collect()
    }

    /// One subblock over every frame: `(n_frames, rows, cols)`.
    fn table<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let t = self.tables.get(name).ok_or_else(|| {
            PyErr::new::<PyKeyError, _>(format!(
                "{name:?} is not declared by the library; have {:?}",
                self.tables.keys().collect::<Vec<_>>()
            ))
        })?;
        let n = self.times.len();
        let flat = PyArray1::from_slice_bound(py, &t.values);
        flat.reshape([n, t.rows, t.cols])
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
    }

    /// A tidy table of one subblock — `frame`, `time`, `row`, `col`, `group_i`, `group_j`,
    /// `value` — as a dict of equal-length arrays, ready for `polars.DataFrame(...)` and a
    /// `group_by`. `group_i`/`group_j` are `-1` for a subblock that is not per group pair.
    fn to_long<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyDict>> {
        let t = self.tables.get(name).ok_or_else(|| {
            PyErr::new::<PyKeyError, _>(format!("{name:?} is not declared by the library"))
        })?;
        let per_frame = t.rows * t.cols;
        let n = per_frame * self.times.len();
        let pairs = t.rows == self.group_pairs().len() && self.n_energy_groups > 0;

        let mut frame = Vec::with_capacity(n);
        let mut time = Vec::with_capacity(n);
        let mut row = Vec::with_capacity(n);
        let mut col = Vec::with_capacity(n);
        let mut gi = Vec::with_capacity(n);
        let mut gj = Vec::with_capacity(n);
        for (f, t0) in self.times.iter().enumerate() {
            for r in 0..t.rows {
                let (i, j) = if pairs {
                    group_pair(r, self.n_energy_groups)
                        .map_or((-1, -1), |(a, b)| (a as i64, b as i64))
                } else {
                    (-1, -1)
                };
                for c in 0..t.cols {
                    frame.push(f as i64);
                    time.push(*t0);
                    row.push(r as i64);
                    col.push(c as i64);
                    gi.push(i);
                    gj.push(j);
                }
            }
        }

        let out = PyDict::new_bound(py);
        out.set_item("frame", PyArray1::from_vec_bound(py, frame))?;
        out.set_item("time", PyArray1::from_vec_bound(py, time))?;
        out.set_item("row", PyArray1::from_vec_bound(py, row))?;
        out.set_item("col", PyArray1::from_vec_bound(py, col))?;
        out.set_item("group_i", PyArray1::from_vec_bound(py, gi))?;
        out.set_item("group_j", PyArray1::from_vec_bound(py, gj))?;
        out.set_item("value", PyArray1::from_slice_bound(py, &t.values))?;
        Ok(out)
    }

    /// Every requested property as a dict of arrays.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let out = PyDict::new_bound(py);
        out.set_item("time", PyArray1::from_slice_bound(py, &self.times))?;
        for (name, series) in &self.properties {
            out.set_item(name, PyArray1::from_slice_bound(py, series))?;
        }
        Ok(out)
    }

    /// The last frame of a subblock as `(rows, cols)` — the matrix a `group_by` would rebuild.
    fn last<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let t = self.tables.get(name).ok_or_else(|| {
            PyErr::new::<PyKeyError, _>(format!("{name:?} is not declared by the library"))
        })?;
        let per_frame = t.rows * t.cols;
        let start = t.values.len().saturating_sub(per_frame);
        let flat = PyArray1::from_slice_bound(py, &t.values[start..]);
        flat.reshape([t.rows, t.cols])
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
    }
}

pub fn register_energy_traj(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyEnergyTrajectory>()?;
    Ok(())
}
