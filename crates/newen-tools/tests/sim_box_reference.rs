//! `sim_box` against gromos++: same box, same solvent count, same positions
//! (`tests/data/sim_box/README.md`).

use std::path::{Path, PathBuf};
use std::process::Command;

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sim_box")
}

fn references() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../newen-md/tests/gromosXX_references")
}

struct Cnf {
    atoms: Vec<(String, String, [f64; 3])>, // residue name, atom name, position
    box_dims: [f64; 3],
}

fn parse(text: &str) -> Cnf {
    let mut atoms = Vec::new();
    let mut box_dims = [0.0; 3];
    let mut block = "";
    let mut genbox_line = 0;
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if block.is_empty() {
            block = t;
            genbox_line = 0;
            continue;
        }
        if t == "END" {
            block = "";
            continue;
        }
        match block {
            "POSITION" => {
                let res = line[5..11].trim().to_string();
                let name = line[11..17].trim().to_string();
                let v: Vec<f64> = line[24..]
                    .split_whitespace()
                    .map(|x| x.parse().unwrap())
                    .collect();
                atoms.push((res, name, [v[0], v[1], v[2]]));
            },
            "GENBOX" => {
                genbox_line += 1;
                if genbox_line == 2 {
                    let v: Vec<f64> = t.split_whitespace().map(|x| x.parse().unwrap()).collect();
                    box_dims = [v[0], v[1], v[2]];
                }
            },
            _ => {},
        }
    }
    Cnf { atoms, box_dims }
}

fn sorted_positions(atoms: &[(String, String, [f64; 3])], res: &str) -> Vec<[f64; 3]> {
    let mut v: Vec<[f64; 3]> = atoms
        .iter()
        .filter(|(r, _, _)| r == res)
        .map(|(_, _, p)| *p)
        .collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v
}

#[test]
fn methanol_in_spc_matches_gromospp() {
    let out = Command::new(env!("CARGO_BIN_EXE_sim_box"))
        .arg("@topo")
        .arg(references().join("shared/meoh_spc.top"))
        .arg("@pbc")
        .arg("r")
        .arg("@pos")
        .arg(data().join("meoh.cnf"))
        .arg("@solvent")
        .arg(references().join("water_1000_spc_gridcell/water_1000_spc.cnf"))
        .arg("@minwall")
        .arg("1.0")
        .output()
        .expect("run sim_box");
    assert!(
        out.status.success(),
        "sim_box failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let ours = parse(&String::from_utf8_lossy(&out.stdout));
    let theirs =
        parse(&std::fs::read_to_string(data().join("meoh_spc_minwall1.gromospp.cnf")).unwrap());

    for d in 0..3 {
        assert!(
            (ours.box_dims[d] - theirs.box_dims[d]).abs() < 1e-8,
            "box {:?} vs gromos++ {:?}",
            ours.box_dims,
            theirs.box_dims
        );
    }
    assert_eq!(ours.atoms.len(), theirs.atoms.len(), "atom count");
    let waters = theirs.atoms.iter().filter(|(r, _, _)| r == "SOLV").count() / 3;
    assert_eq!(waters, 352, "the reference carries 352 waters");

    // solute, in file order
    for (a, b) in ours.atoms.iter().zip(&theirs.atoms).take(3) {
        assert_eq!(a.1, b.1);
        for d in 0..3 {
            assert!((a.2[d] - b.2[d]).abs() < 1e-8, "solute {a:?} vs {b:?}");
        }
    }
    // solvent, as a set (gromos++ and we may emit the molecules in a different order)
    let (mine, ref_) = (
        sorted_positions(&ours.atoms, "SOLV"),
        sorted_positions(&theirs.atoms, "SOLV"),
    );
    assert_eq!(mine.len(), ref_.len());
    for (a, b) in mine.iter().zip(&ref_) {
        for d in 0..3 {
            assert!(
                (a[d] - b[d]).abs() < 1e-8,
                "solvent position {a:?} vs gromos++ {b:?}"
            );
        }
    }
}
