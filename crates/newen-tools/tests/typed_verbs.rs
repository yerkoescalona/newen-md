//! The typed verb and the program's text are one computation.
//!
//! Each preparation-chain program is `run` = read files + verb + serialise. These tests call the
//! verb directly, serialise its result with the program's own writer, and require that to equal
//! what `run` writes — so `newen.tools.sim_box(...) -> Configuration` and
//! `newen.tools.run("sim_box", ...) -> str` can never describe different systems.

use newen_io::coordinate::format_g96;
use newen_io::topology::{read_topology_file, write_parsed_topology};
use newen_io::{read_g96_labeled, G96Atom};
use newen_tools::programs::{build_box, com_top, ion, make_top, pdb2g96, red_top, sim_box};
use std::path::{Path, PathBuf};

fn data(tool: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(tool)
}

fn shared() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../newen-md/tests/gromosXX_references/shared")
}

fn s(p: &Path) -> String {
    p.to_str().unwrap().to_string()
}

/// What `run` writes for `argv`.
fn text_of(program: &str, args: &[String]) -> String {
    let p = newen_tools::find(program).expect("registered");
    let mut argv = vec![program.to_string()];
    argv.extend(args.iter().cloned());
    let argv = newen_io::normalize(argv).expect("normalize");
    let mut buf = Vec::new();
    (p.run)(&argv, &mut buf).unwrap_or_else(|e| panic!("{program}: {e}"));
    String::from_utf8(buf).expect("utf-8")
}

fn argv(items: &[&str]) -> Vec<String> {
    items.iter().map(|s| s.to_string()).collect()
}

#[test]
fn sim_box_solvate_equals_its_text() {
    let topo_path = shared().join("meoh_spc.top");
    let pos_path = data("sim_box").join("meoh.cnf");
    let solvent_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../newen-md/tests/gromosXX_references/water_1000_spc_gridcell/water_1000_spc.cnf");

    let topo = read_topology_file(&topo_path).unwrap();
    let solute = read_g96_labeled(&pos_path).unwrap();
    let solvent = read_g96_labeled(&solvent_path).unwrap();
    let opts = sim_box::SimBoxOptions {
        minwall: vec![1.0],
        ..Default::default()
    };

    let result = sim_box::solvate(
        &solute.atoms,
        &topo.masses,
        None,
        &solvent.atoms,
        solvent.box_dims,
        &opts,
    )
    .unwrap();

    // 352 SPC molecules in a 2.2 nm box — the reference figure (tests/data/sim_box/README.md).
    assert_eq!(result.n_solvent_molecules, 352);
    assert_eq!(result.atoms.len(), solute.atoms.len() + 352 * 3);

    let title = sim_box::title(&s(&pos_path), &s(&solvent_path), result.n_solvent_molecules);
    let typed = format_g96(&title, &result.atoms, Some(result.box_dims));

    let text = text_of(
        "sim_box",
        &argv(&[
            "@topo",
            &s(&topo_path),
            "@pbc",
            "r",
            "@pos",
            &s(&pos_path),
            "@solvent",
            &s(&solvent_path),
            "@minwall",
            "1.0",
        ]),
    );
    assert_eq!(typed, text);
}

#[test]
fn ion_add_ions_equals_its_text() {
    let topo_path = shared().join("aladip.topo");
    let pos_path = shared().join("aladip.conf");

    let topo = read_topology_file(&topo_path).unwrap();
    let system = read_g96_labeled(&pos_path).unwrap();
    let opts = ion::IonOptions {
        positive: Some((1, "NA+".to_string())),
        potential: Some(1.4),
        ..Default::default()
    };

    let result = ion::add_ions(&topo, &system, &opts).unwrap();
    assert_eq!(result.n_positive, 1);
    // one water (3 atoms) became one ion (1 atom)
    assert_eq!(result.atoms.len(), system.atoms.len() - 2);

    let typed = ion::format(&result, &s(&pos_path), &opts);
    let text = text_of(
        "ion",
        &argv(&[
            "@topo",
            &s(&topo_path),
            "@pbc",
            "r",
            "@positive",
            "1",
            "NA+",
            "@potential",
            "1.4",
            "@pos",
            &s(&pos_path),
        ]),
    );
    assert_eq!(typed, text);
}

#[test]
fn com_top_equals_its_text() {
    let a_path = shared().join("aladip.topo");
    let b_path = shared().join("ch4_spc.top");
    let a = read_topology_file(&a_path).unwrap();
    let b = read_topology_file(&b_path).unwrap();

    let combined = com_top::com_top(&[(1, &a), (1, &b)], 1, 1).unwrap();
    assert_eq!(combined.n_atoms, a.n_atoms + b.n_atoms);

    let title = format!(
        "COM_TOP: Combined topology using:\n  1 : {}\n  2 : {}\n",
        s(&a_path),
        s(&b_path)
    );
    let typed = write_parsed_topology(&combined, &title);
    let text = text_of(
        "com_top",
        &argv(&[
            "@topo",
            &s(&a_path),
            &s(&b_path),
            "@param",
            "1",
            "@solv",
            "1",
        ]),
    );
    assert_eq!(typed, text);
}

#[test]
fn red_top_reduce_equals_its_text() {
    let topo_path = shared().join("aladip.topo");
    let topo = read_topology_file(&topo_path).unwrap();

    let reduced = red_top::reduce(&topo, 1, 6).unwrap();
    assert_eq!(reduced.n_atoms, 6);

    let typed = write_parsed_topology(
        &reduced,
        &format!("red_top: atoms 1-6 of {}", s(&topo_path)),
    );
    let text = text_of(
        "red_top",
        &argv(&["@topo", &s(&topo_path), "@atoms", "1:1-6"]),
    );
    assert_eq!(typed, text);
}

#[test]
fn red_top_rejects_a_range_outside_the_topology() {
    let topo = read_topology_file(shared().join("aladip.topo")).unwrap();
    assert!(red_top::reduce(&topo, 1, topo.n_atoms + 1).is_err());
    assert!(red_top::reduce(&topo, 0, 3).is_err());
    assert!(red_top::reduce(&topo, 5, 4).is_err());
}

#[test]
fn build_box_equals_its_text() {
    let topo_path = shared().join("meoh_spc.top");
    let pos_path = data("sim_box").join("meoh.cnf");
    let topo = read_topology_file(&topo_path).unwrap();
    let molecule: Vec<G96Atom> = read_g96_labeled(&pos_path).unwrap().atoms;

    let result = build_box::build_box(&molecule, &topo.masses, 2, 790.0).unwrap();
    assert_eq!(result.n_molecules, 8);
    assert_eq!(result.atoms.len(), molecule.len() * 8);

    let typed = build_box::format(&result, &s(&pos_path), 790.0);
    let text = text_of(
        "build_box",
        &argv(&[
            "@topo",
            &s(&topo_path),
            "@pos",
            &s(&pos_path),
            "@nsm",
            "2",
            "@dens",
            "790",
        ]),
    );
    assert_eq!(typed, text);
}

#[test]
fn pdb2g96_equals_its_text() {
    let topo_path = shared().join("aladip.topo");
    let pdb_path = data("pdb2g96").join("aladip.pdb");

    let topo = read_topology_file(&topo_path).unwrap();
    let pdb = newen_io::pdb::PDBStructure::read_pdb(&pdb_path).unwrap();
    let result = pdb2g96::pdb2g96(&topo, &pdb, &pdb2g96::Library::new(), false).unwrap();
    assert_eq!(result.atoms.len(), topo.n_atoms);

    let typed = pdb2g96::format(&result, &s(&pdb_path));
    let text = text_of(
        "pdb2g96",
        &argv(&["@topo", &s(&topo_path), "@pdb", &s(&pdb_path)]),
    );
    assert_eq!(typed, text);
}

#[test]
fn make_top_reads_back_the_topology_it_writes() {
    let mtb = newen_io::mtb::read_mtb_file(data("make_top").join("54a7.mtb")).unwrap();
    let ifp = newen_io::ifp::read_ifp_file(data("make_top").join("54a7.ifp")).unwrap();
    let seq: Vec<String> = ["NH3+", "ALA", "GLY", "COO-"]
        .iter()
        .map(|s| s.to_string())
        .collect();

    let result = make_top::make_top(&mtb, &ifp, &seq, Some("H2O")).unwrap();

    // The typed topology is the text, parsed — same atoms, same bonded terms, same solvent.
    assert!(result.topology.n_atoms > 0);
    assert!(!result.topology.solvent_atoms.is_empty());
    assert_eq!(result.topology.masses.len(), result.topology.n_atoms);
    assert_eq!(result.topology.charges.len(), result.topology.n_atoms);

    let text = text_of(
        "make_top",
        &argv(&[
            "@build",
            &s(&data("make_top").join("54a7.mtb")),
            "@param",
            &s(&data("make_top").join("54a7.ifp")),
            "@seq",
            "NH3+",
            "ALA",
            "GLY",
            "COO-",
            "@solv",
            "H2O",
        ]),
    );
    assert_eq!(result.text, text);
}

/// The missing-solvent check must fire before any output exists — gromos++ (and our own binary
/// until 3.10) exited half-way through writing, leaving a truncated topology on stdout.
#[test]
fn make_top_rejects_a_missing_solvent_before_writing() {
    let mtb = newen_io::mtb::read_mtb_file(data("make_top").join("54a7.mtb")).unwrap();
    let ifp = newen_io::ifp::read_ifp_file(data("make_top").join("54a7.ifp")).unwrap();
    let seq = vec!["ALA".to_string()];

    let err = make_top::make_top(&mtb, &ifp, &seq, Some("NO_SUCH_SOLVENT")).unwrap_err();
    assert!(err.contains("NO_SUCH_SOLVENT"), "{err}");

    let p = newen_tools::find("make_top").unwrap();
    let mut buf = Vec::new();
    let argv = argv(&[
        "make_top",
        "--build",
        &s(&data("make_top").join("54a7.mtb")),
        "--param",
        &s(&data("make_top").join("54a7.ifp")),
        "--seq",
        "ALA",
        "--solv",
        "NO_SUCH_SOLVENT",
    ]);
    assert!((p.run)(&argv, &mut buf).is_err());
    assert!(buf.is_empty(), "a failed make_top still wrote a topology");
}

/// `check_top` is the one program whose output survives its own failure: the report *is* the
/// answer, and `run` writes it before returning the error.
#[test]
fn check_top_writes_its_report_even_when_it_fails() {
    let topo_path = shared().join("aladip.topo");
    let parsed = read_topology_file(&topo_path).unwrap();
    let built = newen_io::topology::build_topology(parsed);
    let report = newen_tools::programs::check_top::check_top(&built);

    assert!(!report.text.is_empty());
    assert_eq!(report.has_errors(), !report.errors.is_empty());

    let p = newen_tools::find("check_top").unwrap();
    let argv = argv(&["check_top", "--topo", &s(&topo_path)]);
    let mut buf = Vec::new();
    let result = (p.run)(&argv, &mut buf);
    let text = String::from_utf8(buf).unwrap();

    assert!(text.contains(&report.text), "run did not write the report");
    assert_eq!(result.is_err(), report.has_errors());
}
