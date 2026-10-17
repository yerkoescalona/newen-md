//! AtomSelection on topologies as read from file — the parsed counterpart of the unit tests in
//! `newen-core/src/selection.rs`, which build the same layouts by hand.

use std::path::PathBuf;

use newen_core::selection::{AtomSelection, SelectionError};
use newen_io::topology::{build_topology, read_topology_file};

fn reference_topology(system: &str) -> newen_core::Topology {
    reference_topology_at(&format!("{system}/{system}.topo"))
}

/// A topology under `newen-md/tests/gromosXX_references/`.
fn reference_topology_at(rel: &str) -> newen_core::Topology {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../newen-md/tests/gromosXX_references")
        .join(rel);
    build_topology(read_topology_file(&path).expect("reference topology reads"))
}

/// water_216_box.topo: 216 three-atom SOL molecules, all solute (SOLVENTATOM is empty).
#[test]
fn water_216_solute_molecules() {
    let t = reference_topology("water_216_box");
    let count = |spec: &str| AtomSelection::from_string(spec, &t).map(|s| s.len());

    assert_eq!(count("a:OW").unwrap(), 216);
    assert_eq!(count("1-216:1").unwrap(), 216);
    assert_eq!(count("2:1").unwrap(), 1);
    assert!(matches!(
        AtomSelection::from_string("s:OW", &t),
        Err(SelectionError::NoSolvent)
    ));
    assert_eq!(t.atom_name(3), Some("OW"));
    assert_eq!(t.residue_nr(647), Some(216));
}

// ── `a:` against gromos++ `atominfo` on aladip.topo (`tests/data/atominfo`) ──────────────────

fn atominfo_reference(case: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/atominfo")
        .join(format!("{case}.gromospp.out"));
    std::fs::read_to_string(path).expect("gromos++ atominfo reference reads")
}

/// The data rows of an `atominfo` ATOMS block, whitespace-split. gromos++ adds two `#` header
/// lines and writes a compressed specifier in TITLE; the rows themselves are the comparison.
fn atom_rows(out: &str) -> Vec<Vec<String>> {
    out.lines()
        .skip_while(|l| l.trim() != "ATOMS")
        .skip(1)
        .take_while(|l| l.trim() != "END")
        .filter(|l| !l.trim_start().starts_with('#'))
        .map(|l| l.split_whitespace().map(String::from).collect())
        .collect()
}

const ATOMINFO_CASES: &[(&str, &str)] = &[
    ("a:OW", "a_OW"),
    ("a:N", "a_N"),
    ("a:a", "a_a"),
    ("1:a", "1_a"),
];

/// Our `atominfo` binary prints the same rows as gromos++ wherever gromos++ selects atoms.
#[test]
fn atominfo_rows_match_gromospp() {
    let topo = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../newen-md/tests/gromosXX_references/shared/aladip.topo");
    for &(spec, case) in ATOMINFO_CASES {
        let theirs = atom_rows(&atominfo_reference(case));
        let out = std::process::Command::new(env!("CARGO_BIN_EXE_atominfo"))
            .arg("@topo")
            .arg(&topo)
            .args(["@atomspec", spec])
            .output()
            .expect("run atominfo");
        if theirs.is_empty() {
            // gromos++ prints an empty ATOMS block; ours refuses an empty selection.
            assert!(!out.status.success(), "{spec}: expected an empty selection");
            continue;
        }
        assert!(
            out.status.success(),
            "{spec}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert_eq!(
            atom_rows(&String::from_utf8_lossy(&out.stdout)),
            theirs,
            "{spec}"
        );
    }
}

/// gromos++ resolves `a:` over solute molecules whether or not solvent is present. With aladip
/// solvated to the 72 atoms of `aladip_solvated`, ours selects the same GROMOS atom numbers.
#[test]
fn a_prefix_on_solvated_aladip_matches_gromospp() {
    let mut t = reference_topology_at("shared/aladip.topo");
    t.solvate(20);
    assert_eq!(t.num_atoms(), 72);
    for &(spec, case) in ATOMINFO_CASES {
        let theirs: Vec<usize> = atom_rows(&atominfo_reference(case))
            .iter()
            .map(|row| row[1].parse().expect("GROMOS number"))
            .collect();
        let ours: Vec<usize> = match AtomSelection::from_string(spec, &t) {
            Ok(s) => s.iter().map(|i| i + 1).collect(),
            Err(SelectionError::EmptySelection) => Vec::new(),
            Err(e) => panic!("{spec}: {e}"),
        };
        assert_eq!(ours, theirs, "{spec}");
    }
    assert_eq!(AtomSelection::from_string("s:OW", &t).unwrap().len(), 20);
}
