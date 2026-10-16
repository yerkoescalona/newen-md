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

// ── Perturbation topologies (`gromos.free_energy.Perturbation`) ─────────────────────────────

/// 0-based atom indices from a gromos++ atom specifier (`"1:1"`, `"1:CM"`) or a sequence of
/// 0-based indices.
fn atom_indices(topology: &PyTopology, atoms: &Bound<'_, PyAny>) -> PyResult<Vec<usize>> {
    if let Ok(spec) = atoms.extract::<String>() {
        let sel = gromos_core::selection::AtomSelection::from_string(&spec, &topology.inner)
            .map_err(|e| verb_err(format!("atoms {spec:?}: {e}")))?;
        return Ok(sel.indices().to_vec());
    }
    atoms.extract::<Vec<usize>>().map_err(|_| {
        PyErr::new::<pyo3::exceptions::PyTypeError, _>(
            "atoms: an atom specifier (str) or a sequence of 0-based atom indices",
        )
    })
}

fn molecule_charges(charges: Vec<(usize, gromos_core::NetCharge)>) -> Vec<(usize, f64, f64)> {
    charges.into_iter().map(|(m, q)| (m, q.a, q.b)).collect()
}

/// `make_pt_top` from `topology` (state A) to its `make_pt_top::mutate`d copy (state B):
/// returns `(ptp text, [(molecule, Q_A, Q_B)])`, molecules 0-based. `iac_b` is the topology
/// file's (1-based) IAC.
#[pyfunction]
#[pyo3(signature = (topology, atoms, *, charge_b = None, iac_b = None, mass_b = None, alpha_lj = 1.51, alpha_crf = 0.5, title = None))]
#[allow(clippy::too_many_arguments)]
fn _make_pt_top_mutate(
    topology: &PyTopology,
    atoms: &Bound<'_, PyAny>,
    charge_b: Option<Vec<f64>>,
    iac_b: Option<Vec<usize>>,
    mass_b: Option<Vec<f64>>,
    alpha_lj: f64,
    alpha_crf: f64,
    title: Option<&str>,
) -> PyResult<(String, Vec<(usize, f64, f64)>)> {
    use gromos_tools::programs::make_pt_top::{self, Mutation};
    let atoms = atom_indices(topology, atoms)?;
    let iac = iac_b
        .map(|v| {
            v.into_iter()
                .map(|t| {
                    t.checked_sub(1)
                        .ok_or_else(|| verb_err("iac_b: IAC numbers start at 1".into()))
                })
                .collect::<PyResult<Vec<usize>>>()
        })
        .transpose()?;
    let mutation = Mutation {
        charge: charge_b.as_deref(),
        iac: iac.as_deref(),
        mass: mass_b.as_deref(),
    };
    let topo_b = make_pt_top::mutate(&topology.inner, &atoms, &mutation).map_err(verb_err)?;
    let (text, charges) = make_pt_top::make_pt_top(
        &topology.inner,
        &topo_b,
        alpha_lj,
        alpha_crf,
        title.unwrap_or("gromos.free_energy.Perturbation.mutate"),
    )
    .map_err(verb_err)?;
    Ok((text, molecule_charges(charges)))
}

fn parse_ptp(text: &str) -> PyResult<gromos_io::ptp::PerturbedTopology> {
    gromos_io::ptp::parse_pttopo(text).map_err(|e| verb_err(format!("perturbation topology: {e}")))
}

/// The solute's net charge in state A and B with the perturbation `text` merged into a copy of
/// `topology`: what `Topology.apply_perturbation` returns.
#[pyfunction]
fn _ptp_net_charge(topology: &PyTopology, text: &str) -> PyResult<(f64, f64)> {
    let pt = parse_ptp(text)?;
    let mut topo = topology.inner.clone();
    let q = topo.apply_perturbation(pt.perturbed_solute, pt.is_perturbed);
    Ok((q.a, q.b))
}

/// `pt_top`: state B of `topology` with the perturbation `text` applied from atom
/// `first_atom` (0-based), and `[(molecule, Q_A, Q_B)]`.
#[pyfunction]
#[pyo3(signature = (topology, text, *, first_atom = 0))]
fn _pt_top(
    topology: &PyTopology,
    text: &str,
    first_atom: usize,
) -> PyResult<(PyTopology, Vec<(usize, f64, f64)>)> {
    let pt = parse_ptp(text)?;
    let (b, charges) =
        gromos_tools::programs::pt_top::pt_top(topology.parsed.clone(), &pt, first_atom)
            .map_err(verb_err)?;
    Ok((PyTopology::from_parsed(b), molecule_charges(charges)))
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

// ── Free-energy estimators (item 1.8) ───────────────────────────────────────────────────────

/// `ext_ti_ana` in process (`gromos_analysis::programs::ext_ti_ana`): the predictions of every
/// softness combination as arrays, and the files when `write` is set. `argv` is in `@key` form
/// without the program name; without `@library` the built-in library is used. `seed` fixes the
/// bootstrap stream (gromos++ takes the wall clock).
#[pyfunction]
#[pyo3(signature = (argv, *, seed = None, write = true))]
fn _ext_ti_ana(
    py: Python<'_>,
    argv: Vec<String>,
    seed: Option<u32>,
    write: bool,
) -> PyResult<Bound<'_, PyDict>> {
    use gromos_analysis::programs::ext_ti_ana;
    use numpy::{PyArray1, PyArray2};

    let mut full = vec![ext_ti_ana::PROGRAM.name.to_string()];
    full.extend(argv);
    let full = gromos_io::normalize(full).map_err(|e| program_error(py, e, String::new()))?;
    let result = py.allow_threads(|| -> Result<_, String> {
        let args =
            gromos_io::args::Arguments::from_vec(full, ext_ti_ana::KNOWN, ext_ti_ana::USAGE)?;
        let mut opts = ext_ti_ana::Options::from_args(&args, false)?;
        opts.seed = seed;
        let analysis = ext_ti_ana::analyse(&opts)?;
        let files = if write {
            ext_ti_ana::write_files(&opts, &analysis)?
        } else {
            Vec::new()
        };
        Ok((analysis, files))
    });
    let (a, files) = result.map_err(|e| program_error(py, e, String::new()))?;

    let rows = |cols: &[Vec<f64>]| -> PyResult<Bound<'_, PyArray2<f64>>> {
        PyArray2::from_vec2_bound(py, cols)
            .map_err(|e| PyRuntimeError::new_err(format!("ragged table: {e}")))
    };
    let predictions = PyList::empty_bound(py);
    for p in &a.predictions {
        let d = PyDict::new_bound(py);
        d.set_item("labels", p.labels)?;
        d.set_item("slj", p.slj.to_vec())?;
        d.set_item("scrf", p.scrf.to_vec())?;
        d.set_item("plam", PyArray1::from_slice_bound(py, &p.plam))?;
        d.set_item("dhdl", PyArray1::from_slice_bound(py, &p.dhdl))?;
        match &p.error {
            Some(e) => d.set_item("error", PyArray1::from_slice_bound(py, e))?,
            None => d.set_item("error", py.None())?,
        }
        d.set_item("contributing_frames", p.contributing_frames.clone())?;
        match &p.bar {
            Some(b) => {
                let bd = PyDict::new_bound(py);
                bd.set_item("plam", PyArray1::from_slice_bound(py, &b.plam))?;
                bd.set_item("e_s", PyArray1::from_slice_bound(py, &b.e_s))?;
                bd.set_item("e_p", rows(&b.e_p)?)?;
                d.set_item("bar", bd)?;
            },
            None => d.set_item("bar", py.None())?,
        }
        match &p.dhdl_data {
            Some(h) => {
                let hd = PyDict::new_bound(py);
                hd.set_item("plam", PyArray1::from_slice_bound(py, &h.plam))?;
                hd.set_item("de", rows(&h.de)?)?;
                d.set_item("dhdl_data", hd)?;
            },
            None => d.set_item("dhdl_data", py.None())?,
        }
        predictions.append(d)?;
    }
    let out = PyDict::new_bound(py);
    out.set_item("n_frames", a.n_frames)?;
    out.set_item("predictions", predictions)?;
    out.set_item("files", files)?;
    Ok(out)
}

/// `bar` in process (`gromos_analysis::programs::bar::bar`) on `ext_ti_ana @bar_data` files:
/// per interval and in total, kJ/mol, with the table the program prints. `parity` is `@parity`
/// (`"gromospp"` reproduces gromos++ 9cd9f6b3's error; `None` the corrected one).
#[pyfunction]
#[pyo3(signature = (files, temperature, *, max_iterations = 500, convergence = 1e-5, bootstrap = 0, seed = None, printdist = None, parity = None))]
#[allow(clippy::too_many_arguments)]
fn _bar(
    py: Python<'_>,
    files: Vec<String>,
    temperature: f64,
    max_iterations: usize,
    convergence: f64,
    bootstrap: usize,
    seed: Option<u32>,
    printdist: Option<String>,
    parity: Option<String>,
) -> PyResult<Bound<'_, PyDict>> {
    use gromos_analysis::programs::bar;
    use gromos_analysis::programs::free_energy::Parity;

    let parity = parity
        .as_deref()
        .map_or(Ok(Parity::default()), Parity::parse)
        .map_err(|e| program_error(py, e, String::new()))?;
    let result = py.allow_threads(|| -> Result<_, (String, String)> {
        let data = files
            .iter()
            .map(|p| bar::read_bar_file(p))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| (e, String::new()))?;
        let opts = bar::BarOptions {
            temp: temperature,
            max_iterations,
            convergence,
            bootstrap,
            seed,
            printdist: printdist.map(std::path::PathBuf::from),
            parity,
        };
        let mut buf: Vec<u8> = Vec::new();
        let r = bar::bar(&data, &opts, &mut buf);
        let text = String::from_utf8_lossy(&buf).into_owned();
        r.map(|r| (r, text.clone())).map_err(|e| (e, text))
    });
    let (r, stdout) = result.map_err(|(e, text)| program_error(py, e, text))?;
    bar_result_dict(py, &r, stdout)
}

/// The dict `_bar` returns for a `bar::BarResult` (kJ/mol) and the table it printed.
fn bar_result_dict<'py>(
    py: Python<'py>,
    r: &gromos_analysis::programs::bar::BarResult,
    stdout: String,
) -> PyResult<Bound<'py, PyDict>> {
    use gromos_analysis::programs::bar;

    let d = PyDict::new_bound(py);
    let col = |f: fn(&bar::Interval) -> f64| r.intervals.iter().map(f).collect::<Vec<f64>>();
    d.set_item("lam_i", col(|i| i.lam_i))?;
    d.set_item("lam_j", col(|i| i.lam_j))?;
    d.set_item("dg", col(|i| i.dg))?;
    d.set_item("error", col(|i| i.error))?;
    d.set_item("error_corrected", col(|i| i.error_corrected))?;
    d.set_item(
        "bootstrap",
        r.intervals.iter().map(|i| i.bootstrap).collect::<Vec<_>>(),
    )?;
    d.set_item("overlap", col(|i| i.overlap))?;
    d.set_item(
        "iterations",
        r.intervals.iter().map(|i| i.iterations).collect::<Vec<_>>(),
    )?;
    d.set_item("dg_total", r.dg_total)?;
    d.set_item("error_total", r.error_total)?;
    d.set_item("error_total_corrected", r.error_total_corrected)?;
    d.set_item("bootstrap_total", r.bootstrap_total)?;
    d.set_item("stdout", stdout)?;
    Ok(d)
}

/// `ext_ti_merge` in process (`gromos_analysis::programs::ext_ti_merge::evaluate`): the merged
/// curve, its integral and the text the program prints. `argv` in `@key` form.
#[pyfunction]
#[pyo3(signature = (argv))]
fn _ext_ti_merge(py: Python<'_>, argv: Vec<String>) -> PyResult<Bound<'_, PyDict>> {
    use gromos_analysis::programs::ext_ti_merge;
    use numpy::PyArray1;

    let mut full = vec![ext_ti_merge::PROGRAM.name.to_string()];
    full.extend(argv);
    let full = gromos_io::normalize(full).map_err(|e| program_error(py, e, String::new()))?;
    let (result, buf) = py.allow_threads(|| {
        let mut buf: Vec<u8> = Vec::new();
        let result = ext_ti_merge::evaluate(&full, &mut buf);
        (result, buf)
    });
    let stdout = String::from_utf8_lossy(&buf).into_owned();
    let m = result.map_err(|e| program_error(py, e, stdout.clone()))?;
    let d = PyDict::new_bound(py);
    d.set_item("lambda", PyArray1::from_slice_bound(py, &m.lambda))?;
    d.set_item("dhdl", PyArray1::from_slice_bound(py, &m.dhdl))?;
    match &m.error {
        Some(e) => d.set_item("error", PyArray1::from_slice_bound(py, e))?,
        None => d.set_item("error", py.None())?,
    }
    d.set_item("integral", m.integral)?;
    d.set_item("integral_error", m.integral_error)?;
    d.set_item("integral_error_corrected", m.integral_error_corrected)?;
    d.set_item("stdout", stdout)?;
    Ok(d)
}

/// The analytic charging-free-energy correction terms
/// (`gromos_analysis::free_energy::charge_correction`): a dict of the terms in kJ/mol.
#[pyfunction]
#[pyo3(signature = (
    scheme, charge_a, charge_b, box_edge, cutoff, eps_rf, solvent, *, eps_solvent, eps_model,
    ion_radius, density, quadrupole_trace, four_pi_eps_i
))]
#[allow(clippy::too_many_arguments)]
fn _charge_correction<'py>(
    py: Python<'py>,
    scheme: &str,
    charge_a: f64,
    charge_b: f64,
    box_edge: f64,
    cutoff: f64,
    eps_rf: f64,
    solvent: &str,
    eps_solvent: f64,
    eps_model: Option<f64>,
    ion_radius: f64,
    density: f64,
    quadrupole_trace: Option<f64>,
    four_pi_eps_i: f64,
) -> PyResult<Bound<'py, PyDict>> {
    use gromos_analysis::free_energy::charge_correction as cc;
    use pyo3::exceptions::PyValueError;

    let solvent = *cc::SolventModel::find(solvent).ok_or_else(|| {
        let known: Vec<_> = cc::SOLVENTS.iter().map(|m| m.name).collect();
        PyValueError::new_err(format!(
            "unknown solvent '{solvent}' ({})",
            known.join(", ")
        ))
    })?;
    let input = cc::ChargeCorrectionInput {
        scheme: scheme.parse().map_err(PyValueError::new_err)?,
        charge_a,
        charge_b,
        box_edge,
        cutoff,
        eps_rf,
        eps_solvent,
        eps_model,
        ion_radius,
        density,
        solvent,
        quadrupole_trace,
        four_pi_eps_i,
    };
    let c = py
        .allow_threads(|| cc::charge_correction(&input))
        .map_err(PyValueError::new_err)?;
    let d = PyDict::new_bound(py);
    d.set_item("a1", c.a1)?;
    d.set_item("a2", c.a2)?;
    d.set_item("b", c.b)?;
    d.set_item("c1", c.c1)?;
    d.set_item("c2", c.c2)?;
    d.set_item("total", c.total())?;
    d.set_item("beyond_pb", c.beyond_pb())?;
    d.set_item("delta_phi_odl", c.delta_phi_odl)?;
    d.set_item("quadrupole_trace", c.quadrupole_trace)?;
    Ok(d)
}

/// MBAR (`gromos_analysis::programs::mbar::mbar`) on potentials `u_kn` in kJ/mol (K × N, the
/// samples of state 0 first) at `temperature`: free energies and their errors in kJ/mol.
#[pyfunction]
#[pyo3(signature = (u_kn, n_k, temperature, *, max_iterations = 10000, tolerance = 1e-10))]
fn _mbar<'py>(
    py: Python<'py>,
    u_kn: numpy::PyReadonlyArray2<'py, f64>,
    n_k: Vec<usize>,
    temperature: f64,
    max_iterations: usize,
    tolerance: f64,
) -> PyResult<Bound<'py, PyDict>> {
    use gromos_analysis::programs::mbar;

    let kbt = gromos_analysis::GROMOSPP_KB * temperature;
    let view = u_kn.as_array();
    let rows: Vec<Vec<f64>> = view
        .rows()
        .into_iter()
        .map(|r| r.iter().map(|u| u / kbt).collect())
        .collect();
    let r = py
        .allow_threads(|| mbar::mbar(&rows, &n_k, max_iterations, tolerance))
        .map_err(|e| program_error(py, e, String::new()))?;
    mbar_result_dict(py, &r, kbt)
}

/// The dict `_mbar` returns for an `mbar::MbarResult` (kT units) at `kbt`: kJ/mol.
fn mbar_result_dict<'py>(
    py: Python<'py>,
    r: &gromos_analysis::programs::mbar::MbarResult,
    kbt: f64,
) -> PyResult<Bound<'py, PyDict>> {
    use numpy::{PyArray1, PyArray2};

    let scale = |m: &[Vec<f64>]| -> Vec<Vec<f64>> {
        m.iter()
            .map(|r| r.iter().map(|x| x * kbt).collect())
            .collect()
    };
    let d = PyDict::new_bound(py);
    let f: Vec<f64> = r.f.iter().map(|x| x * kbt).collect();
    d.set_item("f", PyArray1::from_vec_bound(py, f))?;
    let to2 = |m: Vec<Vec<f64>>| {
        PyArray2::from_vec2_bound(py, &m).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    };
    d.set_item("dg", to2(scale(&r.delta_f))?)?;
    d.set_item("ddg", to2(scale(&r.d_delta_f))?)?;
    d.set_item("kT", kbt)?;
    d.set_item("self_consistent_iterations", r.self_consistent_iterations)?;
    d.set_item("newton_iterations", r.newton_iterations)?;
    Ok(d)
}

/// One simulated λ's extended-TI frames, as `gromos.free_energy` passes them: a dict with
/// `lam`, `energies` and `derivatives` (frames × grid λ × 12) and `abdih`, `free_abdih`
/// (frames × 2) — `Simulation.precalclam` of every sampled step, stacked.
#[derive(FromPyObject)]
struct PrecalcWindow<'py> {
    #[pyo3(item)]
    lam: f64,
    #[pyo3(item)]
    energies: numpy::PyReadonlyArray3<'py, f64>,
    #[pyo3(item)]
    derivatives: numpy::PyReadonlyArray3<'py, f64>,
    #[pyo3(item)]
    abdih: numpy::PyReadonlyArray2<'py, f64>,
    #[pyo3(item)]
    free_abdih: numpy::PyReadonlyArray2<'py, f64>,
}

/// The frames of `w` after the first `skip`.
fn precalc_frames(
    w: &PrecalcWindow<'_>,
    skip: usize,
) -> PyResult<Vec<gromos_analysis::programs::ext_ti_ana::PrecalcFrame>> {
    use gromos_analysis::programs::ext_ti_ana::PrecalcFrame;

    let (e, d) = (w.energies.as_array(), w.derivatives.as_array());
    let (a, fa) = (w.abdih.as_array(), w.free_abdih.as_array());
    let n = e.shape()[0];
    if d.shape()[0] != n || a.shape()[0] != n || fa.shape()[0] != n {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "lambda {}: energies, derivatives, abdih and free_abdih need one row per frame",
            w.lam
        )));
    }
    if a.shape()[1] != 2 || fa.shape()[1] != 2 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "lambda {}: abdih and free_abdih need two columns (A, B)",
            w.lam
        )));
    }
    let table = |t: numpy::ndarray::ArrayView2<'_, f64>| -> Vec<Vec<f64>> {
        t.outer_iter().map(|r| r.to_vec()).collect()
    };
    Ok((skip..n)
        .map(|f| PrecalcFrame {
            energies: table(e.index_axis(numpy::ndarray::Axis(0), f)),
            derivatives: table(d.index_axis(numpy::ndarray::Axis(0), f)),
            abdih: [a[[f, 0]], a[[f, 1]]],
            free_abdih: [fa[[f, 0]], fa[[f, 1]]],
        })
        .collect())
}

/// `ext_ti_ana @bar_data` of every window at every window's λ, from frames held in memory
/// (`ext_ti_ana::bar_file`): what `bar` and `mbar` read. `grid` is PRECALCLAM's
/// `(NRLAM, MINLAM, MAXLAM)`, `nlam` the simulations' NLAM.
fn precalc_bar_files(
    py: Python<'_>,
    windows: &[PrecalcWindow<'_>],
    grid: (usize, f64, f64),
    nlam: f64,
    temperature: f64,
    skip: usize,
) -> PyResult<Vec<gromos_analysis::programs::bar::BarFile>> {
    use gromos_analysis::programs::ext_ti_ana;

    let lams: Vec<f64> = windows.iter().map(|w| w.lam).collect();
    let frames = windows
        .iter()
        .map(|w| precalc_frames(w, skip))
        .collect::<PyResult<Vec<_>>>()?;
    py.allow_threads(|| {
        lams.iter()
            .zip(&frames)
            .map(|(&lam, f)| {
                let o = ext_ti_ana::Options::for_frames(
                    temperature,
                    grid.0,
                    grid.1,
                    grid.2,
                    lam,
                    nlam,
                    lams.clone(),
                )?;
                ext_ti_ana::bar_file(&o, f).map_err(|e| format!("window lambda = {lam}: {e}"))
            })
            .collect::<Result<Vec<_>, String>>()
    })
    .map_err(|e| program_error(py, e, String::new()))
}

/// The per-frame energies BAR and MBAR use, from in-memory PRECALCLAM frames: per window
/// `slam`, `plam` (the window λ), `e_s` (at `slam`) and `e_p` (at each `plam`), kJ/mol.
#[pyfunction]
#[pyo3(signature = (windows, grid, nlam, temperature, skip = 0))]
fn _precalc_bar_data<'py>(
    py: Python<'py>,
    windows: Vec<PrecalcWindow<'py>>,
    grid: (usize, f64, f64),
    nlam: f64,
    temperature: f64,
    skip: usize,
) -> PyResult<Bound<'py, PyList>> {
    use numpy::{PyArray1, PyArray2};

    let files = precalc_bar_files(py, &windows, grid, nlam, temperature, skip)?;
    let out = PyList::empty_bound(py);
    for f in &files {
        let d = PyDict::new_bound(py);
        d.set_item("slam", f.slam)?;
        d.set_item("plam", PyArray1::from_slice_bound(py, &f.plam))?;
        d.set_item("e_s", PyArray1::from_slice_bound(py, &f.e_s))?;
        d.set_item(
            "e_p",
            PyArray2::from_vec2_bound(py, &f.e_p)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))?,
        )?;
        out.append(d)?;
    }
    Ok(out)
}

/// `bar` (`gromos_analysis::programs::bar::bar`, corrected error) over in-memory PRECALCLAM
/// frames: the dict `_bar` returns.
#[pyfunction]
#[pyo3(signature = (windows, grid, nlam, temperature, skip = 0, *, max_iterations = 500, convergence = 1e-5))]
#[allow(clippy::too_many_arguments)]
fn _precalc_bar<'py>(
    py: Python<'py>,
    windows: Vec<PrecalcWindow<'py>>,
    grid: (usize, f64, f64),
    nlam: f64,
    temperature: f64,
    skip: usize,
    max_iterations: usize,
    convergence: f64,
) -> PyResult<Bound<'py, PyDict>> {
    use gromos_analysis::programs::bar;

    let files = precalc_bar_files(py, &windows, grid, nlam, temperature, skip)?;
    let opts = bar::BarOptions {
        temp: temperature,
        max_iterations,
        convergence,
        bootstrap: 0,
        seed: None,
        printdist: None,
        parity: Default::default(),
    };
    let result = py.allow_threads(|| {
        let mut buf: Vec<u8> = Vec::new();
        let r = bar::bar(&files, &opts, &mut buf);
        (r, String::from_utf8_lossy(&buf).into_owned())
    });
    let r = result
        .0
        .map_err(|e| program_error(py, e, result.1.clone()))?;
    bar_result_dict(py, &r, result.1)
}

/// MBAR (`gromos_analysis::programs::mbar`) over in-memory PRECALCLAM frames, every window's
/// energies at every window's λ: the dict `_mbar` returns, plus `lambdas` (ascending).
#[pyfunction]
#[pyo3(signature = (windows, grid, nlam, temperature, skip = 0, *, max_iterations = 10000, tolerance = 1e-10))]
#[allow(clippy::too_many_arguments)]
fn _precalc_mbar<'py>(
    py: Python<'py>,
    windows: Vec<PrecalcWindow<'py>>,
    grid: (usize, f64, f64),
    nlam: f64,
    temperature: f64,
    skip: usize,
    max_iterations: usize,
    tolerance: f64,
) -> PyResult<Bound<'py, PyDict>> {
    use gromos_analysis::programs::mbar;

    let files = precalc_bar_files(py, &windows, grid, nlam, temperature, skip)?;
    let kbt = gromos_analysis::GROMOSPP_KB * temperature;
    let (lam, r) = py
        .allow_threads(|| -> Result<_, String> {
            let (lam, u_kn, n_k) = mbar::reduced_potentials(&files, kbt)?;
            Ok((lam, mbar::mbar(&u_kn, &n_k, max_iterations, tolerance)?))
        })
        .map_err(|e| program_error(py, e, String::new()))?;
    let d = mbar_result_dict(py, &r, kbt)?;
    d.set_item("lambdas", lam)?;
    Ok(d)
}

/// The plain multi-window TI of `gromos_analysis::ti::trapezoid_from_trg`: ⟨∂H/∂λ⟩ per `.trg`
/// and the trapezoid ΔG (kJ/mol).
#[pyfunction]
#[pyo3(signature = (files, lambdas = None, skip = 0))]
fn _ti_trapezoid(
    py: Python<'_>,
    files: Vec<String>,
    lambdas: Option<Vec<f64>>,
    skip: usize,
) -> PyResult<Bound<'_, PyDict>> {
    let r = py
        .allow_threads(|| gromos_analysis::ti::trapezoid_from_trg(&files, lambdas.as_deref(), skip))
        .map_err(|e| program_error(py, e, String::new()))?;
    ti_result_dict(py, &r)
}

/// Plain multi-window TI over in-memory dH/dλ series, one per λ, each after dropping its
/// first `skip` samples: the estimator of `_ti_trapezoid` (`gromos_analysis::ti`).
#[pyfunction]
#[pyo3(signature = (lambdas, series, skip=0))]
fn _ti_trapezoid_series(
    py: Python<'_>,
    lambdas: Vec<f64>,
    series: Vec<Vec<f64>>,
    skip: usize,
) -> PyResult<Bound<'_, PyDict>> {
    if lambdas.len() != series.len() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "{} lambda values for {} series: give one per series",
            lambdas.len(),
            series.len()
        )));
    }
    let windows = lambdas
        .iter()
        .zip(&series)
        .filter_map(|(&l, s)| {
            gromos_analysis::ti::TiWindow::from_samples(l, s.iter().skip(skip).copied())
        })
        .collect();
    let r = gromos_analysis::ti::trapezoid(windows);
    ti_result_dict(py, &r)
}

fn ti_result_dict<'py>(
    py: Python<'py>,
    r: &gromos_analysis::ti::TiResult,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    let w = &r.windows;
    d.set_item("lambda", w.iter().map(|x| x.lambda).collect::<Vec<_>>())?;
    d.set_item("mean", w.iter().map(|x| x.mean).collect::<Vec<_>>())?;
    d.set_item("ee", w.iter().map(|x| x.ee).collect::<Vec<_>>())?;
    d.set_item("n", w.iter().map(|x| x.n).collect::<Vec<_>>())?;
    d.set_item("dg", r.dg)?;
    d.set_item("dg_error", r.dg_error)?;
    Ok(d)
}

pub fn register_programs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_ext_ti_ana, m)?)?;
    m.add_function(wrap_pyfunction!(_bar, m)?)?;
    m.add_function(wrap_pyfunction!(_ext_ti_merge, m)?)?;
    m.add_function(wrap_pyfunction!(_mbar, m)?)?;
    m.add_function(wrap_pyfunction!(_ti_trapezoid, m)?)?;
    m.add_function(wrap_pyfunction!(_ti_trapezoid_series, m)?)?;
    m.add_function(wrap_pyfunction!(_precalc_bar_data, m)?)?;
    m.add_function(wrap_pyfunction!(_precalc_bar, m)?)?;
    m.add_function(wrap_pyfunction!(_precalc_mbar, m)?)?;
    m.add_function(wrap_pyfunction!(_dgslv_pbsolv, m)?)?;
    m.add_function(wrap_pyfunction!(_charge_correction, m)?)?;
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
    m.add_function(wrap_pyfunction!(_make_pt_top_mutate, m)?)?;
    m.add_function(wrap_pyfunction!(_ptp_net_charge, m)?)?;
    m.add_function(wrap_pyfunction!(_pt_top, m)?)?;
    m.add("ProgramError", m.py().get_type_bound::<ProgramError>())?;
    Ok(())
}
