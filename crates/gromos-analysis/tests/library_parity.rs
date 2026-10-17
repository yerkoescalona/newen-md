//! The library and the binary are one program.
//!
//! Every registered program is run twice — once through `gromos_analysis::programs::<name>::run`
//! into a buffer, once as the compiled binary — and the two stdouts must be byte-identical. The
//! reference tests say our output matches gromos++'s; these say the in-process path a Python
//! caller takes is the *same* output, not a second implementation of it.
//!
//! A program joins the registry with a case here. Programs that write their results to files in
//! the working directory (`ene_ana`, `mdf`, `dfgrid`, `edyn`, `pocket`, `eds_update_1`, `ext_ti_ana`,
//! `frameout`) are given a missing required argument instead, so the parity check exercises the
//! failure path without touching the working tree — the same choice `gromos-tools` makes.

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

fn s(p: PathBuf) -> String {
    p.to_string_lossy().into_owned()
}

/// One case: the program, the binary's path (`CARGO_BIN_EXE_*`), the arguments, and whether it is
/// a real run or a deliberate failure. Without `succeeds`, a case whose fixture went missing would
/// fail on both sides and still pass the parity check.
struct Case {
    name: &'static str,
    bin: &'static str,
    args: Vec<String>,
    succeeds: bool,
}

fn case(name: &'static str, bin: &'static str, args: &[&str], succeeds: bool) -> Case {
    Case {
        name,
        bin,
        args: args.iter().map(|s| s.to_string()).collect(),
        succeeds,
    }
}

/// `run(argv, out)` from the library must write exactly what the binary writes to stdout, and
/// fail exactly when the binary exits non-zero.
fn assert_same_as_binary(c: &Case) {
    let program =
        gromos_analysis::find(c.name).unwrap_or_else(|| panic!("{} not registered", c.name));

    // The binary: argv[0] is the program name, as a shell would pass it.
    let out = Command::new(c.bin)
        .args(&c.args)
        .output()
        .unwrap_or_else(|e| panic!("{} failed to start: {e}", c.name));

    // The library: the same argv, normalised the way `run_main` normalises it.
    let mut argv = vec![c.name.to_string()];
    argv.extend(c.args.iter().cloned());
    let argv = gromos_io::normalize(argv).expect("normalize");
    let mut buf: Vec<u8> = Vec::new();
    let result = (program.run)(&argv, &mut buf);

    assert_eq!(
        result.is_ok(),
        c.succeeds,
        "{}: expected the run to {}, got {:?}",
        c.name,
        if c.succeeds { "succeed" } else { "fail" },
        result
    );
    assert_eq!(
        result.is_ok(),
        out.status.success(),
        "{}: library {:?} but the binary exited {}",
        c.name,
        result,
        out.status
    );
    assert_eq!(
        String::from_utf8_lossy(&buf),
        String::from_utf8_lossy(&out.stdout),
        "{}: library stdout differs from the binary's",
        c.name
    );
}

fn cases() -> Vec<Case> {
    let topo = s(shared().join("shared/aladip.topo"));
    let conf = s(shared().join("shared/aladip.conf"));
    let trc = s(shared().join("aladip_solvated/expected/trajectory.trc"));
    let jval_file = s(data("jval").join("aladip.jval"));

    let eds = |n: &str| s(data("eds_update").join(n));
    let (e_r, e_a, e_b) = (eds("eR.dat"), eds("eA.dat"), eds("eB.dat"));
    let chain = |f: &str| s(data("ext_ti_ana").join("chain").join(f));
    let (bar0, bar25, bar50) = (
        chain("bar_data_l0.00.dat"),
        chain("bar_data_l0.25.dat"),
        chain("bar_data_l0.50.dat"),
    );
    let (pred0, pred25) = (chain("predicted_l0.00.dat"), chain("predicted_l0.25.dat"));
    let (m1, m2) = (
        s(data("matrix_overlap").join("m1.dat")),
        s(data("matrix_overlap").join("m2.dat")),
    );
    let (mult_r, mult_x1, mult_x2) = (
        s(data("dfmult").join("eR.dat")),
        s(data("dfmult").join("eX1.dat")),
        s(data("dfmult").join("eX2.dat")),
    );
    let (ener_a, ener_b) = (
        s(data("dg_ener").join("eA.dat")),
        s(data("dg_ener").join("eB.dat")),
    );
    let fin = s(data("jepot").join("fin.cnf"));

    let real: &[(&str, &str, Vec<&str>)] = &[
        // ── Structural ───────────────────────────────────────────────
        (
            "atominfo",
            env!("CARGO_BIN_EXE_atominfo"),
            vec!["@topo", &topo, "@atomspec", "1:a"],
        ),
        (
            "rmsd",
            env!("CARGO_BIN_EXE_rmsd"),
            vec![
                "@topo",
                &topo,
                "@pbc",
                "r",
                "@atomsrmsd",
                "1:a",
                "@traj",
                &trc,
            ],
        ),
        (
            "rmsf",
            env!("CARGO_BIN_EXE_rmsf"),
            vec![
                "@topo",
                &topo,
                "@pbc",
                "r",
                "@atomsrmsf",
                "1:a",
                "@traj",
                &trc,
            ],
        ),
        // ── Energy ───────────────────────────────────────────────────
        (
            "ener",
            env!("CARGO_BIN_EXE_ener"),
            vec!["@topo", &topo, "@traj", &conf],
        ),
        // ── Distribution ─────────────────────────────────────────────
        (
            "bilayer_dist",
            env!("CARGO_BIN_EXE_bilayer_dist"),
            vec![
                "@topo",
                &topo,
                "@pbc",
                "r",
                "@atoms",
                "1:a",
                "@selection",
                "s:OW",
                "@grid",
                "10",
                "@traj",
                &trc,
            ],
        ),
        (
            "bilayer_oparam",
            env!("CARGO_BIN_EXE_bilayer_oparam"),
            vec![
                "@topo", &topo, "@pbc", "r", "@atoms", "1:2-11", "@refvec", "0", "0", "1", "@traj",
                &trc,
            ],
        ),
        (
            "rdf",
            env!("CARGO_BIN_EXE_rdf"),
            vec![
                "@topo", &topo, "@pbc", "r", "@centre", "1:1", "@with", "s:OW", "@cut", "1.5",
                "@grid", "50", "@traj", &trc,
            ],
        ),
        // ── Free energy ──────────────────────────────────────────────
        (
            "dfmult",
            env!("CARGO_BIN_EXE_dfmult"),
            vec![
                "@temp",
                "300",
                "@stateR",
                &mult_r,
                "@endstates",
                &mult_x1,
                &mult_x2,
            ],
        ),
        (
            "dg_ener",
            env!("CARGO_BIN_EXE_dg_ener"),
            vec!["@temp", "300", "@stateA", &ener_a, "@stateB", &ener_b],
        ),
        (
            "eds_update_2",
            env!("CARGO_BIN_EXE_eds_update_2"),
            vec![
                "@temp", "300", "@vr", &e_r, "@vy", &e_a, &e_b, "@s", "0.03", "@s_old", "0.06",
                "@EiR", "0", "-5.4", "@update", "1", "@eunder", "0",
            ],
        ),
        (
            "bar",
            env!("CARGO_BIN_EXE_bar"),
            vec!["@files", &bar25, &bar0, &bar50, "@temp", "300"],
        ),
        (
            "ext_ti_merge",
            env!("CARGO_BIN_EXE_ext_ti_merge"),
            vec!["@files", &pred25, &pred0],
        ),
        (
            "mbar",
            env!("CARGO_BIN_EXE_mbar"),
            vec!["@files", &bar0, &bar25, &bar50, "@temp", "300"],
        ),
        // ── Trajectory manipulation ──────────────────────────────────
        (
            "tser",
            env!("CARGO_BIN_EXE_tser"),
            vec![
                "@topo",
                &topo,
                "@pbc",
                "r",
                "@prop",
                "d%1:1,2 a%1:1,2,3 t%1:1,2,3,4 tp%1:2,3,4,5%60",
                "@traj",
                &trc,
                "@dist",
                "5",
            ],
        ),
        // ── NOE / NMR ────────────────────────────────────────────────
        (
            "jepot",
            env!("CARGO_BIN_EXE_jepot"),
            vec![
                "@jval", &jval_file, "@K", "0.5", "@ngrid", "8", "@fin", &fin,
            ],
        ),
        (
            "jval",
            env!("CARGO_BIN_EXE_jval"),
            vec![
                "@topo", &topo, "@pbc", "r", "@jval", &jval_file, "@traj", &trc,
            ],
        ),
        (
            "nhoparam",
            env!("CARGO_BIN_EXE_nhoparam"),
            vec!["@topo", &topo, "@traj", &conf, "@atoms", "1:a"],
        ),
        // ── Special ──────────────────────────────────────────────────
        (
            "gca",
            env!("CARGO_BIN_EXE_gca"),
            vec![
                "@topo",
                &topo,
                "@pbc",
                "r",
                "@prop",
                "a%1:1,2,3%100",
                "@traj",
                &conf,
            ],
        ),
        (
            "matrix_overlap",
            env!("CARGO_BIN_EXE_matrix_overlap"),
            vec!["@m1", &m1, "@m2", &m2, "@dim", "3"],
        ),
        // A coarse grid keeps the four Poisson–Boltzmann solves fast; no @verbose, no file.
        (
            "dGslv_pbsolv",
            env!("CARGO_BIN_EXE_dGslv_pbsolv"),
            vec![
                "@topo",
                &topo,
                "@coord",
                &conf,
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
                "0.3",
                "@probeIAC",
                "4",
            ],
        ),
    ];

    // These write their results into the working directory. A missing required argument makes
    // them fail before they open anything, which is what the parity check needs from them.
    let failing: &[(&str, &str, Vec<&str>)] = &[
        (
            "ene_ana",
            env!("CARGO_BIN_EXE_ene_ana"),
            vec!["@prop", "totene"],
        ),
        ("mdf", env!("CARGO_BIN_EXE_mdf"), vec!["@topo", &topo]),
        ("dfgrid", env!("CARGO_BIN_EXE_dfgrid"), vec!["@topo", &topo]),
        ("edyn", env!("CARGO_BIN_EXE_edyn"), vec!["@topo", &topo]),
        ("pocket", env!("CARGO_BIN_EXE_pocket"), vec!["@topo", &topo]),
        (
            "eds_update_1",
            env!("CARGO_BIN_EXE_eds_update_1"),
            vec!["@temp", "300"],
        ),
        (
            "frameout",
            env!("CARGO_BIN_EXE_frameout"),
            vec!["@topo", &topo],
        ),
        (
            "ext_ti_ana",
            env!("CARGO_BIN_EXE_ext_ti_ana"),
            vec!["@en_files", &bar0, "@library", &bar0],
        ),
    ];

    let mut v: Vec<Case> = real
        .iter()
        .map(|(name, bin, args)| case(name, bin, args, true))
        .collect();
    v.extend(
        failing
            .iter()
            .map(|(name, bin, args)| case(name, bin, args, false)),
    );
    v
}

#[test]
fn library_run_matches_the_binary() {
    for c in cases() {
        assert_same_as_binary(&c);
    }
}

/// Every registered program has at least one case above, so adding one to the registry without a
/// parity case fails here rather than going unchecked.
#[test]
fn every_registered_program_has_a_case() {
    let covered: Vec<&str> = cases().iter().map(|c| c.name).collect();
    for p in gromos_analysis::PROGRAMS {
        assert!(
            covered.contains(&p.name),
            "{} is registered but has no library_parity case",
            p.name
        );
    }
}

/// `--help` is part of what a caller sees, so the library and the binary must answer it the same.
#[test]
fn help_matches_the_binary() {
    for c in cases() {
        assert_same_as_binary(&Case {
            name: c.name,
            bin: c.bin,
            args: vec!["--help".to_string()],
            succeeds: false, // `Arguments` answers --help with `Err(usage)`
        });
    }
}

/// Every program here parses with `Arguments`, so `--help` is its registered usage text: Python's
/// `programs()` and the command line cannot describe a program differently.
#[test]
fn programs_answer_help_with_their_usage() {
    for p in gromos_analysis::PROGRAMS {
        let argv = vec![p.name.to_string(), "--help".to_string()];
        let mut buf = Vec::new();
        let err = (p.run)(&argv, &mut buf).expect_err("--help is an Err(usage)");
        assert_eq!(err, p.usage, "{}: --help is not the usage text", p.name);
        assert!(buf.is_empty(), "{}: --help wrote to stdout", p.name);
    }
}
