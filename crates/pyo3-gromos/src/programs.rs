//! The gromos++-style programs, in process.
//!
//! `gromos_tools::PROGRAMS` is the registry; each entry's `run(argv, out)` is the same code the
//! binary runs, so `gromos.tools.run("sim_box", ...)` and `sim_box @topo ...` on the command line
//! cannot produce different output. Nothing here parses arguments or formats results — argv comes
//! from Python already in gromos++ `@key` form and goes through `gromos_io::normalize`, exactly as
//! `run_main` does for the command line.

use pyo3::create_exception;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

create_exception!(
    gromos,
    ProgramError,
    PyRuntimeError,
    "A GROMOS program failed. `.stdout` carries whatever it had written before failing — for \
     `check_top` that is the report itself."
);

/// Raise `ProgramError` with the message and whatever the program managed to write.
fn program_error(py: Python<'_>, message: String, stdout: String) -> PyErr {
    let err = ProgramError::new_err(message);
    if let Ok(value) = err.value_bound(py).getattr("args") {
        let _ = value; // keep `args` as the message alone; stdout is a named attribute
    }
    let _ = err.value_bound(py).setattr("stdout", stdout);
    err
}

/// Run a program by name and return what it wrote to stdout.
///
/// `argv` is the argument list *without* the program name, in gromos++ form
/// (`["@topo", "a.top", "@pbc", "r"]`); `@f` argument files are expanded here, as on the command
/// line. Diagnostics the program writes to stderr stay on the process stderr.
#[pyfunction]
#[pyo3(signature = (name, argv))]
fn _run_program(py: Python<'_>, name: &str, argv: Vec<String>) -> PyResult<String> {
    let program = gromos_tools::find(name)
        .ok_or_else(|| ProgramError::new_err(format!("no GROMOS program named '{name}'")))?;

    let mut full = vec![name.to_string()];
    full.extend(argv);
    let full = gromos_io::normalize(full).map_err(|e| program_error(py, e, String::new()))?;

    // The program is pure Rust and touches no Python object: let other threads run.
    let (result, buf) = py.allow_threads(|| {
        let mut buf: Vec<u8> = Vec::new();
        let result = (program.run)(&full, &mut buf);
        (result, buf)
    });

    let text = String::from_utf8_lossy(&buf).into_owned();
    match result {
        Ok(()) => Ok(text),
        Err(e) => Err(program_error(py, e, text)),
    }
}

/// Every program this build provides: name, accepted `@key`s, usage text and status.
#[pyfunction]
fn _programs(py: Python<'_>) -> PyResult<Bound<'_, PyList>> {
    let rows = PyList::empty_bound(py);
    for p in gromos_tools::PROGRAMS {
        let d = PyDict::new_bound(py);
        d.set_item("name", p.name)?;
        d.set_item("known", p.known.to_vec())?;
        d.set_item("usage", p.usage)?;
        d.set_item("status", p.status.as_str())?;
        rows.append(d)?;
    }
    Ok(rows)
}

// ── The preparation chain, typed ────────────────────────────────────────────────────────────
//
// One binding per verb in `gromos_tools::programs`. These take and return objects, so a workflow
// needs no files between steps; `crates/gromos-tools/tests/typed_verbs.rs` pins each verb's result
// to what the program's text says, so the two forms cannot drift.

use crate::py_conf::PyConfiguration;
use crate::topology::PyTopology;
use gromos_tools::programs::{
    build_box, check_top, com_top, ion, make_top, pdb2g96, red_top, sim_box,
};

fn verb_err(e: String) -> PyErr {
    ProgramError::new_err(e)
}

/// `sim_box`: solvate a solute in a box of pre-equilibrated solvent.
#[pyfunction]
#[pyo3(signature = (topology, solute, solvent, *, minwall = None, thresh = 0.23, boxsize = false))]
fn _sim_box(
    topology: &PyTopology,
    solute: &PyConfiguration,
    solvent: &PyConfiguration,
    minwall: Option<Vec<f64>>,
    thresh: f64,
    boxsize: bool,
) -> PyResult<PyConfiguration> {
    let opts = sim_box::SimBoxOptions {
        minwall: minwall.unwrap_or_default(),
        thresh,
        use_boxsize: boxsize,
    };
    let solute_box = (solute.box_dims != gromos_core::math::Vec3::ZERO).then_some(solute.box_dims);
    let result = sim_box::solvate(
        &solute.atoms_or_placeholders(),
        &topology.parsed.masses,
        solute_box,
        &solvent.atoms_or_placeholders(),
        solvent.box_dims,
        &opts,
    )
    .map_err(verb_err)?;
    Ok(PyConfiguration::from_atoms(result.atoms, result.box_dims))
}

/// `ion`: replace solvent molecules by counter-ions.
#[pyfunction]
#[pyo3(signature = (topology, system, *, positive = None, negative = None, potential = None, mindist = 0.0))]
fn _ion(
    topology: &PyTopology,
    system: &PyConfiguration,
    positive: Option<(usize, String)>,
    negative: Option<(usize, String)>,
    potential: Option<f64>,
    mindist: f64,
) -> PyResult<PyConfiguration> {
    let opts = ion::IonOptions {
        positive,
        negative,
        potential,
        mindist,
    };
    let data = gromos_io::LabeledCoordinateData {
        atoms: system.atoms_or_placeholders(),
        box_dims: system.box_dims,
        box_type: 1,
    };
    let result = ion::add_ions(&topology.parsed, &data, &opts).map_err(verb_err)?;
    Ok(PyConfiguration::from_atoms(result.atoms, result.box_dims))
}

/// `build_box`: replicate one molecule on a grid at a target density.
#[pyfunction]
#[pyo3(signature = (topology, molecule, *, nsm, dens))]
fn _build_box(
    topology: &PyTopology,
    molecule: &PyConfiguration,
    nsm: usize,
    dens: f64,
) -> PyResult<PyConfiguration> {
    let result = build_box::build_box(
        &molecule.atoms_or_placeholders(),
        &topology.parsed.masses,
        nsm,
        dens,
    )
    .map_err(verb_err)?;
    let side = result.box_side;
    Ok(PyConfiguration::from_atoms(
        result.atoms,
        gromos_core::math::Vec3::new(side, side, side),
    ))
}

/// `com_top`: combine topologies. `parts` is `[(repeat, topology), …]`.
#[pyfunction]
#[pyo3(signature = (parts, *, param = 1, solv = 1))]
fn _com_top(
    parts: Vec<(usize, PyRef<'_, PyTopology>)>,
    param: usize,
    solv: usize,
) -> PyResult<PyTopology> {
    let refs: Vec<(usize, &gromos_io::topology::ParsedTopology)> =
        parts.iter().map(|(n, t)| (*n, &t.parsed)).collect();
    let combined = com_top::com_top(&refs, param, solv).map_err(verb_err)?;
    Ok(PyTopology::from_parsed(combined))
}

/// `red_top`: keep only atoms `first..=last` (1-based, inclusive).
#[pyfunction]
#[pyo3(signature = (topology, *, first, last))]
fn _red_top(topology: &PyTopology, first: usize, last: usize) -> PyResult<PyTopology> {
    let reduced = red_top::reduce(&topology.parsed, first, last).map_err(verb_err)?;
    Ok(PyTopology::from_parsed(reduced))
}

/// `pdb2g96`: put a PDB's coordinates into the topology's atom order.
#[pyfunction]
#[pyo3(signature = (topology, pdb, *, lib = None, gch = false))]
fn _pdb2g96(
    topology: &PyTopology,
    pdb: &str,
    lib: Option<&str>,
    gch: bool,
) -> PyResult<(PyConfiguration, usize)> {
    let structure = gromos_io::pdb::PDBStructure::read_pdb(pdb)
        .map_err(|e| verb_err(format!("Error reading PDB '{pdb}': {e}")))?;
    let library = match lib {
        Some(path) => pdb2g96::Library::read(path).map_err(verb_err)?,
        None => pdb2g96::Library::new(),
    };
    let result = pdb2g96::pdb2g96(&topology.parsed, &structure, &library, gch).map_err(verb_err)?;
    let missing = result.missing;
    Ok((
        PyConfiguration::from_atoms(result.atoms, gromos_core::math::Vec3::ZERO),
        missing,
    ))
}

/// `make_top`: build a topology from building blocks.
#[pyfunction]
#[pyo3(signature = (*, build, param, seq, solv = None))]
fn _make_top(
    build: Vec<String>,
    param: &str,
    seq: Vec<String>,
    solv: Option<&str>,
) -> PyResult<PyTopology> {
    let mut bb: Option<gromos_io::mtb::BuildingBlocks> = None;
    for path in &build {
        let b = gromos_io::mtb::read_mtb_file(path)
            .map_err(|e| verb_err(format!("Error reading MTB file '{path}': {e:?}")))?;
        match bb.as_mut() {
            None => bb = Some(b),
            Some(first) => {
                first.solute_blocks.extend(b.solute_blocks);
                first.end_blocks.extend(b.end_blocks);
                first.solvent_blocks.extend(b.solvent_blocks);
            },
        }
    }
    let bb = bb.ok_or_else(|| verb_err("at least one building-block file is required".into()))?;
    let ffp = gromos_io::ifp::read_ifp_file(param)
        .map_err(|e| verb_err(format!("Error reading IFP file '{param}': {e:?}")))?;

    let result = make_top::make_top(&bb, &ffp, &seq, solv).map_err(verb_err)?;
    Ok(PyTopology::from_parsed(result.topology))
}

/// `check_top`: run every consistency check. Returns `(errors, warnings, report)`.
#[pyfunction]
fn _check_top(topology: &PyTopology) -> (Vec<String>, Vec<String>, String) {
    let report = check_top::check_top(&topology.inner);
    (report.errors, report.warnings, report.text)
}

// ── Analysis programs whose results are numbers ─────────────────────────────────────────────

/// `dGslv_pbsolv` in process (`gromos_analysis::programs::dgslv_pbsolv::evaluate`, the code the
/// binary runs): the per-atom potentials and free energies of every solve as numbers, with the
/// table and the `@verbose` log as text. `argv` is in `@key` form, without the program name.
#[pyfunction]
#[pyo3(signature = (argv))]
fn _dgslv_pbsolv(py: Python<'_>, argv: Vec<String>) -> PyResult<Bound<'_, PyDict>> {
    use gromos_analysis::programs::dgslv_pbsolv;
    use numpy::PyArray1;

    let mut full = vec![dgslv_pbsolv::PROGRAM.name.to_string()];
    full.extend(argv);
    let full = gromos_io::normalize(full).map_err(|e| program_error(py, e, String::new()))?;
    let (result, buf) = py.allow_threads(|| {
        let mut buf: Vec<u8> = Vec::new();
        let result = dgslv_pbsolv::evaluate(&full, &mut buf);
        (result, buf)
    });
    let stdout = String::from_utf8_lossy(&buf).into_owned();
    let (r, log) = result.map_err(|e| program_error(py, e, stdout.clone()))?;

    let d = PyDict::new_bound(py);
    let array = |v: &[f64]| PyArray1::from_vec_bound(py, v.to_vec());
    d.set_item(
        "atom",
        r.atoms.iter().map(|a| a.atom + 1).collect::<Vec<_>>(),
    )?;
    d.set_item(
        "residue",
        r.atoms.iter().map(|a| a.resnum + 1).collect::<Vec<_>>(),
    )?;
    d.set_item(
        "name",
        r.atoms.iter().map(|a| a.name.clone()).collect::<Vec<_>>(),
    )?;
    d.set_item(
        "charge",
        array(&r.atoms.iter().map(|a| a.charge).collect::<Vec<_>>()),
    )?;
    let iterations = PyDict::new_bound(py);
    for (key, s) in [
        ("npbc_slv", &r.fd.npbc_slv),
        ("npbc_vac", &r.fd.npbc_vac),
        ("pbc_slv", &r.fd.pbc_slv),
        ("pbc_vac", &r.fd.pbc_vac),
    ] {
        d.set_item(key, array(&s.potentials))?;
        d.set_item(format!("dg_{key}"), s.dg)?;
        iterations.set_item(key, s.iterations)?;
    }
    d.set_item("dg_npbc", r.fd.dg_npbc)?;
    d.set_item("dg_pbc", r.fd.dg_pbc)?;
    match &r.fft {
        Some(fft) => {
            for (key, s) in [("fft_ls_pbc", &fft.ls), ("fft_rf_pbc", &fft.rf)] {
                d.set_item(key, array(&s.potentials))?;
                d.set_item(format!("dg_{key}"), s.dg)?;
                iterations.set_item(key, s.steps)?;
            }
        },
        None => {
            for key in ["fft_ls_pbc", "fft_rf_pbc"] {
                d.set_item(key, py.None())?;
                d.set_item(format!("dg_{key}"), py.None())?;
            }
        },
    }
    d.set_item("iterations", iterations)?;
    d.set_item("documented_correction", r.documented_correction())?;
    d.set_item("stdout", stdout)?;
    d.set_item("log", log)?;
    Ok(d)
}

pub fn register_programs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_dgslv_pbsolv, m)?)?;
    m.add_function(wrap_pyfunction!(_run_program, m)?)?;
    m.add_function(wrap_pyfunction!(_programs, m)?)?;
    m.add_function(wrap_pyfunction!(_sim_box, m)?)?;
    m.add_function(wrap_pyfunction!(_ion, m)?)?;
    m.add_function(wrap_pyfunction!(_build_box, m)?)?;
    m.add_function(wrap_pyfunction!(_com_top, m)?)?;
    m.add_function(wrap_pyfunction!(_red_top, m)?)?;
    m.add_function(wrap_pyfunction!(_pdb2g96, m)?)?;
    m.add_function(wrap_pyfunction!(_make_top, m)?)?;
    m.add_function(wrap_pyfunction!(_check_top, m)?)?;
    m.add("ProgramError", m.py().get_type_bound::<ProgramError>())?;
    Ok(())
}
