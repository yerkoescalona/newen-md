//! The analysis programs against gromos++ output (`tests/data/README.md`): same arguments, every
//! number compared to 1e-6 relative (`-nan` literally), every label exactly.

use std::path::{Path, PathBuf};
use std::process::Command;

fn data(tool: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(tool)
}

fn shared() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
}

fn run(bin: &str, args: &[String], cwd: Option<&Path>) -> String {
    let mut cmd = Command::new(bin);
    cmd.args(args);
    if let Some(d) = cwd {
        cmd.current_dir(d);
    }
    let out = cmd.output().expect("run");
    assert!(
        out.status.success(),
        "{bin} failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// Token-wise comparison: numbers within `rel`, everything else exact.
fn assert_same_text(ours: &str, theirs: &str, rel: f64) {
    let ol: Vec<&str> = ours.lines().collect();
    let tl: Vec<&str> = theirs.lines().collect();
    assert_eq!(
        ol.len(),
        tl.len(),
        "line count\n--- ours ---\n{ours}\n--- gromos++ ---\n{theirs}"
    );
    for (n, (a, b)) in ol.iter().zip(&tl).enumerate() {
        let at: Vec<&str> = a.split_whitespace().collect();
        let bt: Vec<&str> = b.split_whitespace().collect();
        assert_eq!(
            at.len(),
            bt.len(),
            "line {}: {a:?} vs gromos++ {b:?}",
            n + 1
        );
        for (x, y) in at.iter().zip(&bt) {
            match (x.parse::<f64>(), y.parse::<f64>()) {
                (Ok(fx), Ok(fy)) if fx.is_finite() && fy.is_finite() => {
                    assert!(
                        (fx - fy).abs() <= rel * fx.abs().max(fy.abs()).max(1e-3),
                        "line {}: {x} vs gromos++ {y}\n{a}\n{b}",
                        n + 1
                    );
                },
                _ => assert_eq!(x, y, "line {}: {a:?} vs gromos++ {b:?}", n + 1),
            }
        }
    }
}

#[test]
fn tser_matches_gromospp() {
    let out = run(
        env!("CARGO_BIN_EXE_tser"),
        &[
            "@topo",
            &shared().join("shared/aladip.topo").to_string_lossy(),
            "@pbc",
            "r",
            "@prop",
            "d%1:1,2 a%1:1,2,3 t%1:1,2,3,4 tp%1:2,3,4,5%60",
            "@traj",
            &shared()
                .join("aladip_solvated/expected/trajectory.trc")
                .to_string_lossy(),
            "@dist",
            "5",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        None,
    );
    let theirs =
        std::fs::read_to_string(data("tser").join("aladip_solvated.gromospp.out")).unwrap();
    assert_same_text(&out, &theirs, 1e-5);
}

#[test]
fn mdf_matches_gromospp() {
    let dir = std::env::temp_dir().join(format!(
        "gromos_mdf_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    run(
        env!("CARGO_BIN_EXE_mdf"),
        &[
            "@topo",
            &shared().join("shared/aladip.topo").to_string_lossy(),
            "@pbc",
            "r",
            "@centre",
            "1:1",
            "1:6",
            "@with",
            "1:3-12",
            "@traj",
            &shared()
                .join("aladip_solvated/expected/trajectory.trc")
                .to_string_lossy(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        Some(&dir),
    );
    // The `:` is gromos++'s atom specifier, in its output filenames and so in ours.
    for name in ["MIN_1:1", "MIN_1:6"] {
        let ours = std::fs::read_to_string(dir.join(format!("{name}.dat"))).unwrap();
        let theirs =
            std::fs::read_to_string(data("mdf").join(format!("{name}.gromospp.dat"))).unwrap();
        assert_same_text(&ours, &theirs, 1e-8);
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn dg_ener_matches_gromospp() {
    let d = data("dg_ener");
    for (extra, reference) in [
        (vec![], "dg_ener.gromospp.out"),
        (vec!["@col", "3", "3"], "dg_ener_col3.gromospp.out"),
    ] {
        let mut args: Vec<String> = [
            "@temp",
            "300",
            "@stateA",
            &d.join("eA.dat").to_string_lossy(),
            "@stateB",
            &d.join("eB.dat").to_string_lossy(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        args.extend(extra.iter().map(|s| s.to_string()));
        let out = run(env!("CARGO_BIN_EXE_dg_ener"), &args, None);
        assert_same_text(
            &out,
            &std::fs::read_to_string(d.join(reference)).unwrap(),
            2e-5,
        );
    }
}

#[test]
fn dfmult_matches_gromospp() {
    let d = data("dfmult");
    let out = run(
        env!("CARGO_BIN_EXE_dfmult"),
        &[
            "@temp",
            "300",
            "@stateR",
            &d.join("eR.dat").to_string_lossy(),
            "@endstates",
            &d.join("eX1.dat").to_string_lossy(),
            &d.join("eX2.dat").to_string_lossy(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("dfmult.gromospp.out")).unwrap(),
        1e-6,
    );
}

#[test]
fn matrix_overlap_matches_gromospp() {
    let d = data("matrix_overlap");
    let out = run(
        env!("CARGO_BIN_EXE_matrix_overlap"),
        &[
            "@m1",
            &d.join("m1.dat").to_string_lossy(),
            "@m2",
            &d.join("m2.dat").to_string_lossy(),
            "@dim",
            "3",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("matrix_overlap.gromospp.out")).unwrap(),
        1e-5,
    );
}

fn topo() -> String {
    shared()
        .join("shared/aladip.topo")
        .to_string_lossy()
        .into_owned()
}

fn trc() -> String {
    shared()
        .join("aladip_solvated/expected/trajectory.trc")
        .to_string_lossy()
        .into_owned()
}

fn strs(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn bilayer_dist_matches_gromospp() {
    let d = data("bilayer_dist");
    let out = run(
        env!("CARGO_BIN_EXE_bilayer_dist"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:a",
            "@selection",
            "s:OW",
            "@grid",
            "10",
            "@traj",
            &trc(),
        ]),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("aladip_grid10.gromospp.out")).unwrap(),
        1e-5,
    );
    let out = run(
        env!("CARGO_BIN_EXE_bilayer_dist"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:a",
            "@selection",
            "s:OW",
            "@grid",
            "10",
            "@density",
            "@traj",
            &trc(),
        ]),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("aladip_grid10_density.gromospp.out")).unwrap(),
        1e-5,
    );
}

#[test]
fn bilayer_oparam_matches_gromospp() {
    let out = run(
        env!("CARGO_BIN_EXE_bilayer_oparam"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:2-11",
            "@refvec",
            "0",
            "0",
            "1",
            "@traj",
            &trc(),
        ]),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(data("bilayer_oparam").join("aladip_2-11.gromospp.out")).unwrap(),
        1e-4,
    );
}

#[test]
fn jval_matches_gromospp() {
    let d = data("jval");
    let jv = d.join("aladip.jval").to_string_lossy().into_owned();
    for (extra, reference) in [
        (vec!["@timeseries"], "timeseries.gromospp.out"),
        (vec!["@rmsd"], "rmsd.gromospp.out"),
        (vec![], "averages.gromospp.out"),
    ] {
        let mut args = strs(&["@topo", &topo(), "@pbc", "r", "@jval", &jv]);
        args.extend(strs(&extra));
        args.extend(strs(&["@traj", &trc()]));
        let out = run(env!("CARGO_BIN_EXE_jval"), &args, None);
        assert_same_text(
            &out,
            &std::fs::read_to_string(d.join(reference)).unwrap(),
            1e-5,
        );
    }
}

#[test]
fn edyn_matches_gromospp_and_projects_on_eigenvectors() {
    let dir = std::env::temp_dir().join(format!(
        "gromos_edyn_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    run(
        env!("CARGO_BIN_EXE_edyn"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:a",
            "@eigenvalues",
            "1",
            "2",
            "@traj",
            &trc(),
        ]),
        Some(&dir),
    );
    let d = data("edyn");
    for (ours, theirs, rel) in [
        ("COVAR.out", "COVAR.gromospp.out", 1e-4),
        ("COVATOM.out", "COVATOM.gromospp.out", 1e-4),
        ("EIFLUC.out", "EIFLUC.gromospp.out", 1e-4),
    ] {
        assert_same_text(
            &std::fs::read_to_string(dir.join(ours)).unwrap(),
            &std::fs::read_to_string(d.join(theirs)).unwrap(),
            rel,
        );
    }
    // eigenvalues: the nine non-zero ones (10 frames → rank ≤ 9) to 1e-4, the rest are numerical zeros
    let ours = std::fs::read_to_string(dir.join("EIVAL.out")).unwrap();
    let theirs = std::fs::read_to_string(d.join("EIVAL.gromospp.out")).unwrap();
    for (a, b) in ours.lines().zip(theirs.lines()).skip(1).take(9) {
        let (x, y): (f64, f64) = (
            a.split_whitespace().nth(1).unwrap().parse().unwrap(),
            b.split_whitespace().nth(1).unwrap().parse().unwrap(),
        );
        assert!((x - y).abs() <= 1e-4 * y.abs(), "eigenvalue {a} vs {b}");
    }
    // projections on eigenvector k fluctuate with variance = eigenvalue k
    let essdyn = std::fs::read_to_string(dir.join("ESSDYN.out")).unwrap();
    let eig: Vec<f64> = ours
        .lines()
        .skip(1)
        .map(|l| l.split_whitespace().nth(1).unwrap().parse().unwrap())
        .collect();
    for (k, line) in essdyn.lines().filter(|l| !l.starts_with('#')).enumerate() {
        let fluct: f64 = line.split_whitespace().nth(2).unwrap().parse().unwrap();
        assert!(
            (fluct * fluct - eig[k]).abs() <= 1e-4 * eig[k].max(1e-12) + 1e-12,
            "EV {}: fluctuation² {} vs eigenvalue {}",
            k + 1,
            fluct * fluct,
            eig[k]
        );
    }
    std::fs::remove_dir_all(dir).ok();
}

fn positions_of(text: &str) -> Vec<Vec<[f64; 3]>> {
    let mut frames = Vec::new();
    let mut cur: Option<Vec<[f64; 3]>> = None;
    for line in text.lines() {
        if line.starts_with("POSITION") {
            cur = Some(Vec::new());
            continue;
        }
        if let Some(c) = cur.as_mut() {
            if line.starts_with("END") {
                frames.push(cur.take().unwrap());
            } else if !line.starts_with('#') && line.len() > 24 {
                let v: Vec<f64> = line[24..]
                    .split_whitespace()
                    .map(|x| x.parse().unwrap())
                    .collect();
                c.push([v[0], v[1], v[2]]);
            }
        }
    }
    frames
}

#[test]
fn gca_matches_gromospp() {
    let d = data("gca");
    let dihedral = std::fs::read_to_string(d.join("dihedral.txt"))
        .unwrap()
        .trim()
        .to_string();
    let conf = shared()
        .join("shared/aladip.conf")
        .to_string_lossy()
        .into_owned();
    for (props, mobile, reference) in [
        (
            format!("t%1:{dihedral}%30%-60%60"),
            None,
            "torsion_scan.gromospp.cnf",
        ),
        (
            "d%1:1,2%0.2 a%1:1,2,3%100".to_string(),
            Some("first"),
            "dist_angle_first.gromospp.cnf",
        ),
        ("a%1:1,2,3%100".to_string(), None, "angle_last.gromospp.cnf"),
    ] {
        let mut args = strs(&["@topo", &topo(), "@pbc", "r", "@prop", &props]);
        if let Some(m) = mobile {
            args.extend(strs(&["@mobile", m]));
        }
        args.extend(strs(&["@traj", &conf]));
        let out = run(env!("CARGO_BIN_EXE_gca"), &args, None);
        let ours = positions_of(&out);
        let theirs = positions_of(&std::fs::read_to_string(d.join(reference)).unwrap());
        assert_eq!(ours.len(), theirs.len(), "{reference}: frame count");
        for (f, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            assert_eq!(a.len(), b.len());
            for (k, (x, y)) in a.iter().zip(b).enumerate() {
                for dim in 0..3 {
                    assert!(
                        (x[dim] - y[dim]).abs() < 1e-7,
                        "{reference} frame {f} atom {}: {x:?} vs gromos++ {y:?}",
                        k + 1
                    );
                }
            }
        }
    }
}

fn tmp_dir(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "gromos_{tag}_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn read(path: &std::path::Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

#[test]
fn eds_update_1_matches_gromospp() {
    let d = data("eds_update");
    let dir = tmp_dir("eds1");
    let f = |n: &str| d.join(n).to_string_lossy().into_owned();
    let (er, ea, eb, ec, tree) = (
        f("eR.dat"),
        f("eA.dat"),
        f("eB.dat"),
        f("eC.dat"),
        f("tree.dat"),
    );
    let cases: Vec<(Vec<&str>, &str, Option<&str>)> = vec![
        (
            vec!["@form", "1", "@s", "0.03"],
            "u1_form1.gromospp.out",
            None,
        ),
        (
            vec!["@form", "2", "@s", "0.03", "0.04", "0.05"],
            "u1_form2.gromospp.out",
            None,
        ),
        (
            vec![
                "@form",
                "3",
                "@s",
                "0.03",
                "0.04",
                "@update_tree",
                "0",
                "@tree",
                &tree,
            ],
            "u1_form3t0.gromospp.out",
            Some("u1_form3t0.tree_new.gromospp.dat"),
        ),
        (
            vec![
                "@form",
                "3",
                "@s",
                "0.03",
                "0.04",
                "@update_tree",
                "1",
                "@tree",
                &tree,
            ],
            "u1_form3t1.gromospp.out",
            Some("u1_form3t1.tree_new.gromospp.dat"),
        ),
    ];
    for (extra, reference, tree_ref) in cases {
        let mut args = strs(&[
            "@temp", "300", "@numstat", "3", "@vr", &er, "@vy", &ea, &eb, &ec, "@EiR", "0", "-5.4",
            "4.3",
        ]);
        args.extend(strs(&extra));
        let out = run(env!("CARGO_BIN_EXE_eds_update_1"), &args, Some(&dir));
        assert_same_text(&out, &read(&d.join(reference)), 1e-5);
        if let Some(t) = tree_ref {
            assert_same_text(&read(&dir.join("tree_new.dat")), &read(&d.join(t)), 1e-5);
        }
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn eds_update_2_matches_gromospp() {
    let d = data("eds_update");
    let f = |n: &str| d.join(n).to_string_lossy().into_owned();
    let (er, ea, eb) = (f("eR.dat"), f("eA.dat"), f("eB.dat"));
    for (extra, reference) in [
        (
            vec!["@update", "1", "@eunder", "-100"],
            "u2_upd1.gromospp.out",
        ),
        (
            vec!["@update", "1", "@eunder", "0"],
            "u2_upd1b.gromospp.out",
        ),
        (
            vec![
                "@update", "2", "@eunder", "5", "@etrans", "2", "@scale", "1.0",
            ],
            "u2_upd2.gromospp.out",
        ),
    ] {
        let mut args = strs(&[
            "@temp", "300", "@vr", &er, "@vy", &ea, &eb, "@s", "0.03", "@s_old", "0.06", "@EiR",
            "0", "-5.4",
        ]);
        args.extend(strs(&extra));
        let out = run(env!("CARGO_BIN_EXE_eds_update_2"), &args, None);
        assert_same_text(&out, &read(&d.join(reference)), 1e-5);
    }
}

#[test]
fn jepot_matches_gromospp() {
    let d = data("jepot");
    let f = |n: &str| d.join(n).to_string_lossy().into_owned();
    let jv = data("jval")
        .join("aladip.jval")
        .to_string_lossy()
        .into_owned();
    let conf = shared()
        .join("shared/aladip.conf")
        .to_string_lossy()
        .into_owned();
    let (fin, res, one_jval, one_trs) = (f("fin.cnf"), f("res.trs"), f("one.jval"), f("one.trs"));
    let topo_path = topo();
    let cases: Vec<(Vec<&str>, &str)> = vec![
        (
            vec!["@jval", &jv, "@K", "0.5", "@ngrid", "8", "@fin", &fin],
            "fin.gromospp.out",
        ),
        (
            vec!["@jval", &jv, "@K", "0.5", "@ngrid", "8", "@restraj", &res],
            "ts_all.gromospp.out",
        ),
        (
            vec![
                "@jval",
                &jv,
                "@K",
                "0.5",
                "@ngrid",
                "8",
                "@restraj",
                &res,
                "@timespec",
                "SPEC",
                "@timepts",
                "2",
                "@time",
                "5",
                "0.5",
            ],
            "ts_spec.gromospp.out",
        ),
        (
            vec![
                "@jval", &one_jval, "@K", "0.5", "@ngrid", "8", "@angles", "CURR", "@topo",
                &topo_path, "@pbc", "r", "@postraj", &conf, "@restraj", &one_trs,
            ],
            "curr.gromospp.out",
        ),
    ];
    for (args, reference) in cases {
        let out = run(env!("CARGO_BIN_EXE_jepot"), &strs(&args), None);
        assert_same_text(&out, &read(&d.join(reference)), 1e-6);
    }
}

#[test]
fn pocket_matches_gromospp() {
    let d = data("pocket");
    let trc_path = trc();
    let center = read(&d.join("center.txt")).trim().to_string();
    let conf = shared()
        .join("shared/aladip.conf")
        .to_string_lossy()
        .into_owned();
    let cases: Vec<(Vec<&str>, &str, Vec<&str>)> = vec![
        (
            vec![
                "@protein",
                "1:a",
                "@radius",
                "0.8",
                "@vec_number_factor",
                "2",
                "@volume_and_area",
                "@final_vector_coords",
                "@traj",
                &conf,
            ],
            "p1",
            vec!["charges", "lengths", "volume", "area", "vector_coords"],
        ),
        (
            vec![
                "@protein",
                "1:a",
                "@reject",
                "1:1-3",
                "@radius",
                "0.8",
                "@vec_number_factor",
                "3",
                "@radH",
                "0.05",
                "@hemisphere",
                "@volume_and_area",
                "@traj",
                &trc_path,
            ],
            "p2",
            vec!["charges", "lengths", "volume", "area"],
        ),
    ];
    for (extra, tag, files) in cases {
        let dir = tmp_dir("pocket");
        let mut args = strs(&["@topo", &topo(), "@pbc", "r", "@center", &center]);
        args.extend(strs(&extra));
        let out = run(env!("CARGO_BIN_EXE_pocket"), &args, Some(&dir));
        assert_same_text(&out, &read(&d.join(format!("{tag}.gromospp.out"))), 1e-6);
        for f in files {
            assert_same_text(
                &read(&dir.join(format!("{f}.txt"))),
                &read(&d.join(format!("{tag}_{f}.gromospp.txt"))),
                1e-6,
            );
        }
        std::fs::remove_dir_all(dir).ok();
    }
}

#[test]
fn dfgrid_matches_gromospp() {
    let d = data("dfgrid");
    let dir = tmp_dir("dfgrid");
    let args = strs(&[
        "@topo",
        &topo(),
        "@pbc",
        "r",
        "@atom",
        "1:1",
        "@distatoms",
        "1:12",
        "1:6",
        "@gridspacing",
        "0.3",
        "@proteinoffset",
        "15",
        "@proteincutoff",
        "0.25",
        "@proteinatoms",
        "1:a",
        "@max",
        "3",
        "@smooth",
        "1",
        "@protect",
        "0.3",
        "@frames",
        "0",
        "2",
        "@traj",
        &trc(),
    ]);
    let out = run(env!("CARGO_BIN_EXE_dfgrid"), &args, Some(&dir));
    assert_same_text(&out, &read(&d.join("distatoms.gromospp.out")), 1e-5);
    for f in ["grid00000", "grid00002"] {
        // the title carries the trajectory path: compare from TIMESTEP on
        let ours = read(&dir.join(format!("{f}.cnf")));
        let theirs = read(&d.join(format!("{f}.gromospp.cnf")));
        let from = |s: &str| s[s.find("TIMESTEP").unwrap()..].to_string();
        assert_same_text(&from(&ours), &from(&theirs), 1e-7);
    }
    std::fs::remove_dir_all(dir).ok();
}

/// `rmsd` and `rmsf` take two independent atom sets — the frame is superimposed using
/// `@atomsfit` and the deviation/fluctuation measured over `@atomsrmsd`/`@atomsrmsf` — and a
/// reference that is either `@ref` or the trajectory's first frame. Each combination is
/// compared against gromos++ (`tests/data/{rmsd,rmsf}/README` rows).
fn structural_case(bin: &str, tool: &str, expected: &str, extra: &[&str]) {
    let topo = shared().join("shared/aladip.topo");
    let traj = shared().join("aladip_solvated/expected/trajectory.trc");
    let mut args: Vec<String> = vec![
        "@topo".into(),
        topo.to_string_lossy().into_owned(),
        "@pbc".into(),
        "r".into(),
    ];
    args.extend(extra.iter().map(|s| {
        if *s == "%REF%" {
            shared()
                .join("aladip_solvated/expected/final.conf")
                .to_string_lossy()
                .into_owned()
        } else {
            (*s).to_string()
        }
    }));
    args.push("@traj".into());
    args.push(traj.to_string_lossy().into_owned());
    let out = run(bin, &args, None);
    let theirs = std::fs::read_to_string(data(tool).join(expected)).unwrap();
    assert_same_text(&out, &theirs, 1e-5);
}

#[test]
fn rmsd_matches_gromospp() {
    let bin = env!("CARGO_BIN_EXE_rmsd");
    structural_case(bin, "rmsd", "fit_all.gromospp.out", &["@atomsrmsd", "1:a"]);
    structural_case(
        bin,
        "rmsd",
        "fit_separate.gromospp.out",
        &["@atomsrmsd", "1:3-8", "@atomsfit", "1:a"],
    );
    structural_case(
        bin,
        "rmsd",
        "ref_conf.gromospp.out",
        &["@atomsrmsd", "1:a", "@ref", "%REF%", "@time", "5", "0.5"],
    );
}

#[test]
fn rmsf_matches_gromospp() {
    let bin = env!("CARGO_BIN_EXE_rmsf");
    structural_case(bin, "rmsf", "fit_all.gromospp.out", &["@atomsrmsf", "1:a"]);
    structural_case(
        bin,
        "rmsf",
        "fit_separate.gromospp.out",
        &["@atomsrmsf", "1:3-8", "@atomsfit", "1:a", "@ref", "%REF%"],
    );
}

/// `ene_ana` reads a `.tre` according to `@library`'s declaration of it, so the reference cases
/// cover a plain run, `@time`, a `@topo`-derived constant (`densit` needs MASS), direct
/// `SUBBLOCK[i]` references and a two-file sequence. stdout *and* every `<prop>.dat` are
/// compared (`tests/data/ene_ana` rows in the README).
fn ene_ana_case(name: &str, extra: &[&str]) {
    let dir = std::env::temp_dir().join(format!(
        "gromos_ene_ana_{name}_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();

    // no `@library`: the built-in one, which `gromos-io` derives from the same 2023-04-15 layout
    // md++ ships (`official_library_is_generated_from_the_layout` keeps the two equal)
    let tre = shared().join("aladip_solvated/expected/energies.tre");
    let mut args: Vec<String> = vec!["@en_files".into(), tre.to_string_lossy().into_owned()];
    args.extend(extra.iter().map(|s| {
        match *s {
            "%TRE%" => tre.to_string_lossy().into_owned(),
            "%TOPO%" => shared()
                .join("shared/aladip.topo")
                .to_string_lossy()
                .into_owned(),
            other => other.to_string(),
        }
    }));
    let out = run(env!("CARGO_BIN_EXE_ene_ana"), &args, Some(&dir));

    let expected = data("ene_ana").join(name);
    let theirs = std::fs::read_to_string(expected.join("stdout.gromospp.out")).unwrap();
    // gromos++'s NaN warning goes to stderr and is part of the stored output
    let warning = "# WARNING: One of the values is a NaN,\n#   the data provided are not enough to \n#   give a sensible error estimate\n";
    let ours = if theirs.starts_with(warning) {
        format!("{warning}{out}")
    } else {
        out
    };
    assert_same_text(&ours, &theirs, 1e-5);

    for entry in std::fs::read_dir(&expected).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_some_and(|e| e == "dat") {
            let file = path.file_name().unwrap();
            let ours = std::fs::read_to_string(dir.join(file))
                .unwrap_or_else(|e| panic!("{name}: our {file:?} missing: {e}"));
            let theirs = std::fs::read_to_string(&path).unwrap();
            assert_same_text(&ours, &theirs, 1e-5);
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn ene_ana_matches_gromospp() {
    ene_ana_case("plain", &["@prop", "totene", "totpot", "totkin"]);
    ene_ana_case("usertime", &["@prop", "totene", "@time", "5", "0.5"]);
    ene_ana_case(
        "topo_densit",
        &["@prop", "densit", "totene", "@topo", "%TOPO%"],
    );
    ene_ana_case("direct_refs", &["@prop", "ENER[1]", "ENER[11]"]);
    ene_ana_case("two_files", &["@en_files", "%TRE%", "@prop", "totene"]);
}

/// One `dGslv_pbsolv` case (`tests/data/dgslv_pbsolv/README.md`): stdout and the `@verbose` log
/// byte-identical to gromos++'s. The finite-difference solver stops on an L1 residual test, so a
/// change in floating-point order shows up as a different iteration count and different digits.
fn dgslv_pbsolv_case(case: &str, args: &[&str]) {
    let dir = tmp_dir(&format!("dgslv_pbsolv_{case}"));
    let log_name = format!("{case}.log");
    let mut argv: Vec<String> = args
        .iter()
        .map(|a| {
            a.replace("%SHARED%", &shared().join("shared").to_string_lossy())
                .replace("%DATA%", &data("dgslv_pbsolv").to_string_lossy())
        })
        .collect();
    argv.extend(["@verbose".to_string(), log_name.clone()]);
    let out = run(env!("CARGO_BIN_EXE_dGslv_pbsolv"), &argv, Some(&dir));
    let reference = data("dgslv_pbsolv");
    assert_eq!(
        out,
        read(&reference.join(format!("{case}.stdout"))),
        "{case}: stdout"
    );

    let ours = read(&dir.join(&log_name));
    let theirs = read(&reference.join(&log_name));
    // The lines that carry results first, so a failure names what moved.
    let results = |log: &str| -> Vec<String> {
        log.lines()
            .filter(|l| {
                [
                    "# CONVERGED AFTER",
                    "# NO CONVERGENCE AFTER",
                    "# DG RETURNED",
                    "# DGRESULT",
                    "# gridstart",
                ]
                .iter()
                .any(|p| l.starts_with(p))
            })
            .map(str::to_string)
            .collect()
    };
    assert_eq!(results(&ours), results(&theirs), "{case}: log results");
    for (n, (a, b)) in ours.lines().zip(theirs.lines()).enumerate() {
        assert_eq!(a, b, "{case}: log line {}", n + 1);
    }
    assert_eq!(ours, theirs, "{case}: log");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn dgslv_pbsolv_matches_gromospp() {
    let coord = [
        "@topo",
        "%SHARED%/aladip.topo",
        "@coord",
        "%SHARED%/aladip.conf",
    ];
    let with = |extra: &[&'static str]| -> Vec<&'static str> {
        coord.iter().copied().chain(extra.iter().copied()).collect()
    };
    dgslv_pbsolv_case(
        "c1",
        &with(&[
            "@pbc",
            "r",
            "@atoms",
            "1:a",
            "@atomsTOcharge",
            "1:a",
            "@schemeELEC",
            "LS",
            "@epsSOLV",
            "66.6",
            "@gridspacing",
            "0.1",
            "@probeIAC",
            "4",
        ]),
    );
    dgslv_pbsolv_case(
        "sigma",
        &with(&[
            "@atoms",
            "1:a",
            "@atomsTOcharge",
            "1:4,5",
            "@schemeELEC",
            "LS",
            "@epsSOLV",
            "66.6",
            "@gridspacing",
            "0.1",
            "@NPBCsize",
            "5.0",
            "@increasegrid",
            "2",
            "2",
            "2",
            "@rminORsigma",
            "1",
            "@radscal",
            "1.1",
            "@probeIAC",
            "4",
        ]),
    );
    dgslv_pbsolv_case(
        "order",
        &with(&[
            "@pbc",
            "r",
            "@atoms",
            "1:5,3,1",
            "1:7",
            "@atomsTOcharge",
            "1:9,2",
            "@schemeELEC",
            "ls",
            "@epsSOLV",
            "66.6",
            "@epsNPBC",
            "66.6",
            "@gridspacing",
            "0.13",
            "@probeIAC",
            "5",
            "@probeRAD",
            "0.12",
            "@maxiter",
            "40",
        ]),
    );
    dgslv_pbsolv_case(
        "pqr",
        &[
            "@pqr",
            "%DATA%/aladip.pqr",
            "@coordinates",
            "3.767055681",
            "3.767055681",
            "3.767055681",
            "@nogridpoints",
            "24",
            "@gridspacing",
            "0.15",
            "@atoms",
            "1-12",
            "@atomsTOcharge",
            "2,3,5",
            "@schemeELEC",
            "LS",
            "@epsSOLV",
            "66.6",
            "@radscal",
            "1.1",
        ],
    );
}
