//! The library and the binary are one program.
//!
//! Every registered program is run twice — once through `newen_tools::programs::<name>::run`
//! into a buffer, once as the compiled binary — and the two stdouts must be byte-identical. The
//! reference tests say our output matches gromos++'s; these say the in-process path a Python
//! caller takes is the *same* output, not a second implementation of it.
//!
//! A program joins the registry with a case here. `--help` is checked the same way, because the
//! usage text is part of what a caller sees.

use std::path::{Path, PathBuf};
use std::process::Command;

fn data(tool: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(tool)
}

fn shared() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../newen-md/tests/gromosXX_references/shared")
}

/// One case: the program, the binary's path (`CARGO_BIN_EXE_*`) and the arguments.
struct Case {
    name: &'static str,
    bin: &'static str,
    args: Vec<String>,
}

fn case(name: &'static str, bin: &'static str, args: &[&str]) -> Case {
    Case {
        name,
        bin,
        args: args.iter().map(|s| s.to_string()).collect(),
    }
}

/// `run(argv, out)` from the library must write exactly what the binary writes to stdout, and
/// fail exactly when the binary exits non-zero.
fn assert_same_as_binary(c: &Case) {
    let program = newen_tools::find(c.name).unwrap_or_else(|| panic!("{} not registered", c.name));

    // The binary: argv[0] is the program name, as a shell would pass it.
    let out = Command::new(c.bin)
        .args(&c.args)
        .output()
        .unwrap_or_else(|e| panic!("{} failed to start: {e}", c.name));

    // The library: the same argv, normalised the way `run_main` normalises it.
    let mut argv = vec![c.name.to_string()];
    argv.extend(c.args.iter().cloned());
    let argv = newen_io::normalize(argv).expect("normalize");
    let mut buf: Vec<u8> = Vec::new();
    let result = (program.run)(&argv, &mut buf);

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
    let meoh2_top = data("explode").join("meoh2.top");
    let meoh2_cnf = data("explode").join("meoh2.cnf");
    let dup_cnf = data("duplicate").join("meoh2_dup.cnf");
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();

    let mut v = vec![
        case(
            "duplicate",
            env!("CARGO_BIN_EXE_duplicate"),
            &["@topo", &s(&meoh2_top), "@pos", &s(&dup_cnf), "@pbc", "r"],
        ),
        case(
            "duplicate",
            env!("CARGO_BIN_EXE_duplicate"),
            &[
                "@topo",
                &s(&meoh2_top),
                "@pos",
                &s(&dup_cnf),
                "@pbc",
                "r",
                "@write",
            ],
        ),
        case(
            "explode",
            env!("CARGO_BIN_EXE_explode"),
            &[
                "@topo",
                &s(&meoh2_top),
                "@pos",
                &s(&meoh2_cnf),
                "@nsm",
                "2",
                "@dist",
                "1.5",
            ],
        ),
        case(
            "pt_top",
            env!("CARGO_BIN_EXE_pt_top"),
            &[
                "@topo",
                &s(&shared().join("ch4_spc.top")),
                "@pttopo",
                &s(&data("make_pt_top").join("ch4_to_B.gromospp.ptp")),
                "@type",
                "TOPO",
            ],
        ),
        case(
            "make_sasa_top",
            env!("CARGO_BIN_EXE_make_sasa_top"),
            &[
                "@topo",
                &s(&shared().join("aladip.topo")),
                "@sasaspec",
                &s(&data("make_sasa_top").join("aladip.sasaspec")),
            ],
        ),
        case(
            "make_sasa_top",
            env!("CARGO_BIN_EXE_make_sasa_top"),
            &[
                "@topo",
                &s(&shared().join("aladip.topo")),
                "@sasaspec",
                &s(&data("make_sasa_top").join("aladip.sasaspec")),
                "@noH",
            ],
        ),
        case(
            "sim_box",
            env!("CARGO_BIN_EXE_sim_box"),
            &[
                "@topo",
                &s(&shared().join("meoh_spc.top")),
                "@pbc",
                "r",
                "@pos",
                &s(&data("sim_box").join("meoh.cnf")),
                "@solvent",
                &s(&Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("../newen-md/tests/gromosXX_references/water_1000_spc_gridcell/water_1000_spc.cnf")
                    .to_path_buf()),
                "@minwall",
                "1.0",
            ],
        ),
        case(
            "red_top",
            env!("CARGO_BIN_EXE_red_top"),
            &[
                "@topo",
                &s(&shared().join("aladip.topo")),
                "@atoms",
                "1:1-6",
            ],
        ),
        case(
            "com_top",
            env!("CARGO_BIN_EXE_com_top"),
            &[
                "@topo",
                &s(&shared().join("aladip.topo")),
                &s(&shared().join("ch4_spc.top")),
                "@param",
                "1",
                "@solv",
                "1",
            ],
        ),
        case(
            "ion",
            env!("CARGO_BIN_EXE_ion"),
            &[
                "@topo",
                &s(&shared().join("aladip.topo")),
                "@pbc",
                "r",
                "@positive",
                "1",
                "NA+",
                "@potential",
                "1.4",
                "@pos",
                &s(&shared().join("aladip.conf")),
            ],
        ),
        case(
            "build_box",
            env!("CARGO_BIN_EXE_build_box"),
            &[
                "@topo",
                &s(&data("sim_box").join("../../../newen-md/tests/gromosXX_references/shared/meoh_spc.top")),
                "@pos",
                &s(&data("sim_box").join("meoh.cnf")),
                "@nsm",
                "2",
                "@dens",
                "790",
            ],
        ),
        case(
            "make_top",
            env!("CARGO_BIN_EXE_make_top"),
            &[
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
            ],
        ),
        case(
            "check_top",
            env!("CARGO_BIN_EXE_check_top"),
            &["@topo", &s(&shared().join("aladip.topo"))],
        ),
        case(
            "pdb2g96",
            env!("CARGO_BIN_EXE_pdb2g96"),
            &[
                "@topo",
                &s(&shared().join("aladip.topo")),
                "@pdb",
                &s(&data("pdb2g96").join("aladip.pdb")),
            ],
        ),
        // A missing required argument: the failure path must agree too.
        case("explode", env!("CARGO_BIN_EXE_explode"), &["@nsm", "2"]),
        // A solvent block that does not exist: the check that used to fire mid-output.
        case(
            "make_top",
            env!("CARGO_BIN_EXE_make_top"),
            &[
                "@build",
                &s(&data("make_top").join("54a7.mtb")),
                "@param",
                &s(&data("make_top").join("54a7.ifp")),
                "@seq",
                "ALA",
                "@solv",
                "NO_SUCH_SOLVENT",
            ],
        ),
        case(
            "ion",
            env!("CARGO_BIN_EXE_ion"),
            &["@topo", &s(&shared().join("aladip.topo")), "@pbc", "r"],
        ),
        case("red_top", env!("CARGO_BIN_EXE_red_top"), &["@natoms", "6"]),
        // Out of range: the typed verb's own validation, surfaced through run.
        case(
            "red_top",
            env!("CARGO_BIN_EXE_red_top"),
            &["@topo", &s(&shared().join("aladip.topo")), "@natoms", "999"],
        ),
    ];

    // The rest of the tools, on the shared aladip system. A `.cnf` reads as a one-frame
    // trajectory, so the `@traj` programs run on it too.
    let topo = s(&shared().join("aladip.topo"));
    let conf = s(&shared().join("aladip.conf"));
    let ch4 = s(&shared().join("ch4_spc.top"));

    let rest: &[(&str, &str, Vec<&str>)] = &[
        (
            "check_box",
            env!("CARGO_BIN_EXE_check_box"),
            vec!["@traj", &conf],
        ),
        (
            "bin_box",
            env!("CARGO_BIN_EXE_bin_box"),
            vec!["@topo", &topo, "@traj", &conf, "@bins", "2", "2", "2"],
        ),
        (
            "unify_box",
            env!("CARGO_BIN_EXE_unify_box"),
            vec!["@topo", &topo, "@traj", &conf, "@box", "3", "3", "3"],
        ),
        ("gch", env!("CARGO_BIN_EXE_gch"), vec!["@traj", &conf]),
        (
            "traj2pdb",
            env!("CARGO_BIN_EXE_traj2pdb"),
            vec!["@traj", &conf],
        ),
        (
            "copy_box",
            env!("CARGO_BIN_EXE_copy_box"),
            vec!["@topo", &topo, "@pos", &conf, "@dir", "x:1"],
        ),
        (
            "inbox",
            env!("CARGO_BIN_EXE_inbox"),
            vec!["@topo", &topo, "@pos", &conf, "@pbc", "r"],
        ),
        (
            "ran_box",
            env!("CARGO_BIN_EXE_ran_box"),
            vec!["@topo", &topo, "@box", "3", "3", "3", "@seed", "7"],
        ),
        (
            "ran_solvation",
            env!("CARGO_BIN_EXE_ran_solvation"),
            vec![
                "@solute", &conf, "@solvent", &conf, "@nsolv", "4", "@box", "3", "3", "3",
            ],
        ),
        (
            "amber2gromos",
            env!("CARGO_BIN_EXE_amber2gromos"),
            vec!["@amber", "prmtop", "@natoms", "5"],
        ),
        (
            "top2build",
            env!("CARGO_BIN_EXE_top2build"),
            vec!["@topo", &topo],
        ),
        (
            "con_top",
            env!("CARGO_BIN_EXE_con_top"),
            vec!["@topo", &topo, &ch4],
        ),
        (
            "link_top",
            env!("CARGO_BIN_EXE_link_top"),
            vec!["@topo1", &topo, "@topo2", &ch4],
        ),
        (
            "pert_top",
            env!("CARGO_BIN_EXE_pert_top"),
            vec!["@topo1", &topo, "@topo2", &topo],
        ),
        (
            "addvirt_top",
            env!("CARGO_BIN_EXE_addvirt_top"),
            vec!["@topo", &topo, "@nvirt", "3"],
        ),
        (
            "prep_eds",
            env!("CARGO_BIN_EXE_prep_eds"),
            vec!["@topo", &topo, "@states", "3"],
        ),
        (
            "prep_noe",
            env!("CARGO_BIN_EXE_prep_noe"),
            vec!["@pairs", "no_such_file"],
        ),
        (
            "prep_xray",
            env!("CARGO_BIN_EXE_prep_xray"),
            vec!["@topo", &topo],
        ),
        // These three write files; run them with a missing required argument so the parity
        // check exercises the failure path without touching the working tree.
        (
            "make_pt_top",
            env!("CARGO_BIN_EXE_make_pt_top"),
            vec!["@topo", &topo],
        ),
        (
            "prep_posres",
            env!("CARGO_BIN_EXE_prep_posres"),
            vec!["@topo", &topo],
        ),
        (
            "mk_script",
            env!("CARGO_BIN_EXE_mk_script"),
            vec!["@sys", "x"],
        ),
    ];
    for (name, bin, args) in rest {
        v.push(case(name, bin, args));
    }
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
    for p in newen_tools::PROGRAMS {
        assert!(
            covered.contains(&p.name),
            "{} is registered but has no library_parity case",
            p.name
        );
    }
}

/// `--help` is part of what a caller sees, and the two parser kinds answer it differently: an
/// `Arguments` program returns `Err(usage)` (stderr, exit 1), a clap program prints its help to
/// stdout and exits 0. Either is fine — the library and the binary must just agree.
#[test]
fn help_matches_the_binary() {
    for c in cases() {
        assert_same_as_binary(&Case {
            name: c.name,
            bin: c.bin,
            args: vec!["--help".to_string()],
        });
    }
}

/// An `Arguments` program's `--help` is its registered usage text, so `programs()` in Python and
/// the command line cannot describe a program differently.
#[test]
fn arguments_programs_answer_help_with_their_usage() {
    for p in newen_tools::PROGRAMS {
        let argv = vec![p.name.to_string(), "--help".to_string()];
        let mut buf = Vec::new();
        match (p.run)(&argv, &mut buf) {
            Err(e) => {
                assert_eq!(e, p.usage, "{}: --help is not the usage text", p.name);
                assert!(buf.is_empty(), "{}: --help wrote to stdout", p.name);
            },
            // a clap program: help went to stdout and the run is over
            Ok(()) => assert!(
                !buf.is_empty(),
                "{}: --help produced neither an error nor output",
                p.name
            ),
        }
    }
}
