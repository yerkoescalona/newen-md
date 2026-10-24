//! `make_top` and `com_top` against gromos++: the topology we write must read back identical to
//! the one gromos++ writes for the same arguments (`tests/data/{make_top,com_top}/README.md`).
//! Both files go through `newen_io::read_topology_file`, so formatting is irrelevant and every
//! block that matters for a simulation is compared.

use std::path::{Path, PathBuf};
use std::process::Command;

use newen_io::topology::{read_topology_file, ParsedTopology};

fn data(tool: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(tool)
}

fn references() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../newen-md/tests/gromosXX_references")
}

fn run(bin: &str, args: &[String]) -> ParsedTopology {
    let out = Command::new(bin).args(args).output().expect("run tool");
    assert!(
        out.status.success(),
        "{bin} failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    static COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let tmp = std::env::temp_dir().join(format!(
        "newen_tools_{}_{}_{}.top",
        Path::new(bin).file_name().unwrap().to_string_lossy(),
        std::process::id(),
        COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    std::fs::write(&tmp, &out.stdout).unwrap();
    let parsed = read_topology_file(&tmp).unwrap_or_else(|e| {
        panic!(
            "our output is not readable: {e}\n{}",
            String::from_utf8_lossy(&out.stdout)
        )
    });
    std::fs::remove_file(tmp).ok();
    parsed
}

/// gromos++ writes its parameters with six significant digits (`%e`); ours carry more.
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-5 * a.abs().max(b.abs()).max(1.0)
}

fn assert_same(ours: &ParsedTopology, theirs: &ParsedTopology) {
    eprintln!(
        "ours: {} atoms, {} bonds, {} angles, {} impropers, {} dihedrals, {} LJ, {} types, {} res | theirs: {} atoms, {} bonds, {} angles, {} impropers, {} dihedrals, {} LJ, {} types, {} res",
        ours.n_atoms, ours.bonds.len(), ours.angles.len(), ours.improper_dihedrals.len(), ours.proper_dihedrals.len(), ours.lj_parameters.len(), ours.atom_type_names.len(), ours.residue_names.len(),
        theirs.n_atoms, theirs.bonds.len(), theirs.angles.len(), theirs.improper_dihedrals.len(), theirs.proper_dihedrals.len(), theirs.lj_parameters.len(), theirs.atom_type_names.len(), theirs.residue_names.len()
    );
    assert_eq!(ours.n_atoms, theirs.n_atoms, "atom count");
    assert_eq!(ours.atom_names, theirs.atom_names, "atom names");
    assert_eq!(
        ours.residue_numbers, theirs.residue_numbers,
        "residue numbers"
    );
    assert_eq!(ours.residue_names, theirs.residue_names, "residue names");
    assert_eq!(
        ours.atom_type_names, theirs.atom_type_names,
        "atom type names"
    );
    assert_eq!(ours.iac, theirs.iac, "IAC");
    assert_eq!(
        ours.chargegroup_codes, theirs.chargegroup_codes,
        "charge groups"
    );
    for (i, (a, b)) in ours.masses.iter().zip(&theirs.masses).enumerate() {
        assert!(close(*a, *b), "mass of atom {}: {a} vs {b}", i + 1);
    }
    for (i, (a, b)) in ours.charges.iter().zip(&theirs.charges).enumerate() {
        assert!(close(*a, *b), "charge of atom {}: {a} vs {b}", i + 1);
    }
    let sorted = |v: &Vec<Vec<usize>>| -> Vec<Vec<usize>> {
        v.iter()
            .map(|l| {
                let mut l = l.clone();
                l.sort();
                l
            })
            .collect()
    };
    assert_eq!(
        sorted(&ours.exclusions),
        sorted(&theirs.exclusions),
        "exclusions"
    );
    assert_eq!(
        sorted(&ours.one_four_pairs),
        sorted(&theirs.one_four_pairs),
        "1-4 pairs"
    );
    let sorted3 = |v: &Vec<(usize, usize, usize)>| {
        let mut v = v.clone();
        v.sort();
        v
    };
    assert_eq!(sorted3(&ours.bonds), sorted3(&theirs.bonds), "bonds");
    let mut a = ours.angles.clone();
    let mut b = theirs.angles.clone();
    a.sort();
    b.sort();
    assert_eq!(a, b, "angles");
    let mut a = ours.improper_dihedrals.clone();
    let mut b = theirs.improper_dihedrals.clone();
    a.sort();
    b.sort();
    assert_eq!(a, b, "improper dihedrals");
    let mut a = ours.proper_dihedrals.clone();
    let mut b = theirs.proper_dihedrals.clone();
    a.sort();
    b.sort();
    assert_eq!(a, b, "dihedrals");
    assert_eq!(
        ours.bond_parameters.len(),
        theirs.bond_parameters.len(),
        "bond types"
    );
    for (i, (a, b)) in ours
        .bond_parameters
        .iter()
        .zip(&theirs.bond_parameters)
        .enumerate()
    {
        assert!(
            close(a.k_quartic, b.k_quartic)
                && close(a.k_harmonic, b.k_harmonic)
                && close(a.r0, b.r0),
            "bond type {}: {a:?} vs {b:?}",
            i + 1
        );
    }
    assert_eq!(
        ours.angle_parameters.len(),
        theirs.angle_parameters.len(),
        "angle types"
    );
    for (i, (a, b)) in ours
        .angle_parameters
        .iter()
        .zip(&theirs.angle_parameters)
        .enumerate()
    {
        assert!(
            close(a.k_cosine, b.k_cosine)
                && close(a.k_harmonic, b.k_harmonic)
                && close(a.theta0, b.theta0),
            "angle type {}: {a:?} vs {b:?}",
            i + 1
        );
    }
    assert_eq!(
        ours.dihedral_parameters.len(),
        theirs.dihedral_parameters.len(),
        "dihedral types"
    );
    for (i, (a, b)) in ours
        .dihedral_parameters
        .iter()
        .zip(&theirs.dihedral_parameters)
        .enumerate()
    {
        assert!(
            close(a.k, b.k) && close(a.pd.cos(), b.pd.cos()) && a.m == b.m,
            "dihedral type {}: {a:?} vs {b:?}",
            i + 1
        );
    }
    assert_eq!(
        ours.improper_dihedral_parameters.len(),
        theirs.improper_dihedral_parameters.len(),
        "improper types"
    );
    for (i, (a, b)) in ours
        .improper_dihedral_parameters
        .iter()
        .zip(&theirs.improper_dihedral_parameters)
        .enumerate()
    {
        assert!(
            close(a.k, b.k) && (a.q0 - b.q0).abs() < 1e-5,
            "improper type {}: {a:?} vs {b:?}",
            i + 1
        );
    }
    assert_eq!(
        ours.lj_parameters.len(),
        theirs.lj_parameters.len(),
        "LJ pairs"
    );
    for (k, b) in &theirs.lj_parameters {
        let a = ours
            .lj_parameters
            .get(k)
            .unwrap_or_else(|| panic!("LJ pair {k:?} missing"));
        assert!(
            close(a.c6, b.c6)
                && close(a.c12, b.c12)
                && close(a.cs6, b.cs6)
                && close(a.cs12, b.cs12),
            "LJ pair {k:?}: {a:?} vs {b:?}"
        );
    }
    assert_eq!(
        ours.solute_molecules, theirs.solute_molecules,
        "solute molecules"
    );
    assert_eq!(
        ours.temperature_groups, theirs.temperature_groups,
        "temperature groups"
    );
    assert_eq!(
        ours.pressure_groups, theirs.pressure_groups,
        "pressure groups"
    );
    assert_eq!(
        ours.solvent_atoms.len(),
        theirs.solvent_atoms.len(),
        "solvent atoms"
    );
    for (a, b) in ours.solvent_atoms.iter().zip(&theirs.solvent_atoms) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.iac, b.iac, "solvent IAC");
        assert!(close(a.mass, b.mass) && close(a.charge, b.charge));
    }
    assert_eq!(
        ours.solvent_constraints.len(),
        theirs.solvent_constraints.len(),
        "solvent constraints"
    );
}

#[test]
fn make_top_peptide_matches_gromospp() {
    let d = data("make_top");
    let s = |p: &str| d.join(p).to_string_lossy().into_owned();
    let args: Vec<String> = [
        "@build",
        &s("54a7.mtb"),
        "@param",
        &s("54a7.ifp"),
        "@seq",
        "NH3+",
        "ALA",
        "GLY",
        "COO-",
        "@solv",
        "H2O",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let ours = run(env!("CARGO_BIN_EXE_make_top"), &args);
    let theirs = read_topology_file(d.join("nh3_ala_gly_coo.gromospp.top")).unwrap();
    assert_eq!(theirs.n_atoms, 14);
    assert_same(&ours, &theirs);
}

#[test]
fn make_top_two_building_block_files_match_gromospp() {
    let d = data("make_top");
    let s = |p: &str| d.join(p).to_string_lossy().into_owned();
    let args: Vec<String> = [
        "@build",
        &s("54a7.mtb"),
        &s("MeOH_exp_H.mtb"),
        "@param",
        &s("54a7.ifp"),
        "@seq",
        "MeOH",
        "@solv",
        "H2O",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let ours = run(env!("CARGO_BIN_EXE_make_top"), &args);
    let theirs = read_topology_file(d.join("meoh_exp_h.gromospp.top")).unwrap();
    assert_eq!(theirs.n_atoms, 6);
    assert_same(&ours, &theirs);
}

#[test]
fn com_top_matches_gromospp() {
    let r = references();
    let args: Vec<String> = [
        "--topo",
        &r.join("shared/aladip.topo").to_string_lossy(),
        &r.join("shared/ch4_spc.top").to_string_lossy(),
        "--param",
        "1",
        "--solv",
        "1",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let ours = run(env!("CARGO_BIN_EXE_com_top"), &args);
    let theirs = read_topology_file(data("com_top").join("aladip_ch4.gromospp.top")).unwrap();
    assert_eq!(theirs.n_atoms, 13);
    assert_same(&ours, &theirs);
}

#[test]
fn red_top_matches_gromospp() {
    let r = references();
    let args: Vec<String> = [
        "@topo",
        &r.join("shared/aladip.topo").to_string_lossy(),
        "@atoms",
        "1:1-6",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let ours = run(env!("CARGO_BIN_EXE_red_top"), &args);
    let theirs = read_topology_file(data("red_top").join("aladip_1_6.gromospp.top")).unwrap();
    assert_eq!(theirs.n_atoms, 6);
    assert_same(&ours, &theirs);
}

#[test]
fn make_pt_top_matches_gromospp() {
    let r = references();
    let d = data("make_pt_top");
    let out = std::env::temp_dir().join(format!("newen_make_pt_top_{}.ptp", std::process::id()));
    let args: Vec<String> = [
        "@topo",
        &r.join("shared/ch4_spc.top").to_string_lossy(),
        &d.join("ch4_B.top").to_string_lossy(),
        "@softpar",
        "1.51",
        "0.5",
        "@out",
        &out.to_string_lossy(),
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let status = Command::new(env!("CARGO_BIN_EXE_make_pt_top"))
        .args(&args)
        .output()
        .expect("run");
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let ours = newen_io::ptp::read_pttopo(&out).expect("our .ptp reads back");
    let theirs =
        newen_io::ptp::read_pttopo(d.join("ch4_to_B.gromospp.ptp")).expect("gromos++ .ptp");
    std::fs::remove_file(out).ok();
    let (ours, theirs) = (&ours.perturbed_solute, &theirs.perturbed_solute);
    assert_eq!(ours.atoms.len(), 1);
    assert_eq!(ours.atoms.len(), theirs.atoms.len());
    for (a, b) in ours.atoms.iter().zip(&theirs.atoms) {
        assert_eq!((a.seq, a.a_iac, a.b_iac), (b.seq, b.a_iac, b.b_iac));
        assert!(close(a.a_charge, b.a_charge) && close(a.b_charge, b.b_charge));
        assert!(close(a.a_mass, b.a_mass) && close(a.b_mass, b.b_mass));
        assert!(close(a.lj_soft, b.lj_soft) && close(a.crf_soft, b.crf_soft));
    }
}

#[test]
fn pt_top_matches_gromospp() {
    let r = references();
    let args: Vec<String> = [
        "@topo",
        &r.join("shared/ch4_spc.top").to_string_lossy(),
        "@pttopo",
        &data("make_pt_top")
            .join("ch4_to_B.gromospp.ptp")
            .to_string_lossy(),
        "@type",
        "TOPO",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let ours = run(env!("CARGO_BIN_EXE_pt_top"), &args);
    let theirs = read_topology_file(data("pt_top").join("ch4_stateB.gromospp.top")).unwrap();
    assert_eq!(theirs.iac[0], 17, "state B IAC 18 applied to atom 1");
    assert_same(&ours, &theirs);
}

/// A +1 → 0 perturbation: `pt_top` puts each perturbed molecule's net charge in A and B into
/// its title, `make_pt_top` prints it (1.8).
#[test]
fn pt_top_and_make_pt_top_print_the_net_charges() {
    use newen_tools::programs::topology::{make_pt_top, pt_top};
    let dir = std::env::temp_dir().join(format!("newen_net_charge_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let top = references().join("shared/ch4_spc.top");
    let ptp = |name: &str, qa: f64, qb: f64| {
        let path = dir.join(name);
        std::fs::write(
            &path,
            format!(
                "TITLE\ncharge\nEND\nPERTATOMPARAM\n1\n    1 1 CM 17 16.043 {qa} 17 16.043 {qb} \
                 0.0 0.0\nEND\n"
            ),
        )
        .unwrap();
        path.to_string_lossy().into_owned()
    };
    let call = |run: fn(&[String], &mut dyn std::io::Write) -> Result<(), String>,
                args: &[&str]| {
        let mut argv = vec!["prog".to_string()];
        argv.extend(args.iter().map(|s| s.to_string()));
        let argv = newen_io::normalize(argv).expect("normalize");
        let mut out = Vec::new();
        run(&argv, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    };
    let top_s = top.to_string_lossy();

    let minus = ptp("minus.ptp", 1.0, 0.0);
    let text = call(
        pt_top::run,
        &["@topo", &top_s, "@pttopo", &minus, "@type", "TOPO"],
    );
    assert!(
        text.contains("net charge of molecule 1: A 1.000 B 0.000 (ΔQ = -1.000)"),
        "{text}"
    );

    // the state-B topology of a 0 → +1 perturbation is the +1 topology
    let plus = ptp("plus.ptp", 0.0, 1.0);
    let plus_top = dir.join("plus.top");
    let text = call(
        pt_top::run,
        &["@topo", &top_s, "@pttopo", &plus, "@type", "TOPO"],
    );
    assert!(text.contains("A 0.000 B 1.000 (ΔQ = 1.000)"), "{text}");
    std::fs::write(&plus_top, text).unwrap();
    let out_ptp = dir.join("made.ptp");
    let text = call(
        make_pt_top::run,
        &[
            "@topo",
            &plus_top.to_string_lossy(),
            &top_s,
            "@out",
            &out_ptp.to_string_lossy(),
        ],
    );
    assert!(
        text.contains("net charge of molecule 1: A 1.000 B 0.000 (ΔQ = -1.000)"),
        "{text}"
    );
    std::fs::remove_dir_all(dir).ok();
}

/// The typed verbs `newen.free_energy.Perturbation` calls: `mutate` + `make_pt_top` write the
/// `.ptp` `make_pt_top @topo A B` writes from the same two states, and `pt_top` applies it back.
#[test]
fn mutate_make_pt_top_and_pt_top_verbs_round_trip() {
    use newen_io::ptp::parse_pttopo;
    use newen_io::topology::build_topology;
    use newen_tools::programs::topology::make_pt_top::{self, Mutation};
    use newen_tools::programs::topology::pt_top;

    let parsed = read_topology_file(references().join("shared/ch4_spc.top")).unwrap();
    let a = build_topology(parsed.clone());
    let plus = Mutation {
        charge: Some(&[1.0]),
        ..Default::default()
    };
    let b = make_pt_top::mutate(&a, &[0], &plus).unwrap();
    assert_eq!((a.charge[0], b.charge[0]), (0.0, 1.0));
    assert_eq!((a.iac[0], a.mass[0]), (b.iac[0], b.mass[0]));
    let (text, charges) = make_pt_top::make_pt_top(&a, &b, 0.0, 0.0, "plus").unwrap();
    assert_eq!(charges.len(), 1);
    assert_eq!(
        (charges[0].0, charges[0].1.a, charges[0].1.b),
        (0, 0.0, 1.0)
    );
    // the file writer the `make_pt_top` binary used to call writes the same text
    let file = std::env::temp_dir().join(format!("newen_verbs_{}.ptp", std::process::id()));
    newen_io::PtpWriter::new(&file)
        .unwrap()
        .write(&a, &b, 0.0, 0.0, "plus")
        .unwrap();
    assert_eq!(text, std::fs::read_to_string(&file).unwrap());
    std::fs::remove_file(&file).ok();

    let pt = parse_pttopo(&text).unwrap();
    assert_eq!(pt.perturbed_solute.atoms.len(), 1);
    let (state_b, charges) = pt_top::pt_top(parsed, &pt, 0).unwrap();
    assert_eq!(state_b.charges[0], 1.0);
    assert_eq!((charges[0].1.a, charges[0].1.b), (0.0, 1.0));

    let err = |atoms: &[usize], m: Mutation| make_pt_top::mutate(&a, atoms, &m).unwrap_err();
    assert!(err(&[], plus).contains("no atoms"));
    assert!(err(&[1], plus).contains("not a solute atom"));
    assert!(err(&[0, 0], plus).contains("twice"));
    assert!(err(&[0], Mutation::default()).contains("nothing to change"));
    let too_many = Mutation {
        charge: Some(&[1.0, 2.0]),
        ..Default::default()
    };
    assert!(err(&[0], too_many).contains("one per atom"));
    let bad_type = Mutation {
        iac: Some(&[999]),
        ..Default::default()
    };
    assert!(err(&[0], bad_type).contains("atom types"));
}

/// gromos++'s output text, token by token (numbers within `rel`); the TITLE block is skipped
/// because ours keeps the input title's trailing blanks.
fn assert_same_text(ours: &str, theirs: &str, rel: f64) {
    let toks = |s: &str| -> Vec<String> {
        s.lines()
            .skip_while(|l| l.trim() != "END")
            .flat_map(|l| l.split_whitespace().map(str::to_string).collect::<Vec<_>>())
            .collect()
    };
    let (a, b) = (toks(ours), toks(theirs));
    assert_eq!(a.len(), b.len(), "token count {} vs {}", a.len(), b.len());
    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        if x == y {
            continue;
        }
        match (x.parse::<f64>(), y.parse::<f64>()) {
            (Ok(fx), Ok(fy)) => assert!(
                (fx - fy).abs() <= rel * fx.abs().max(fy.abs()) + 1e-12,
                "token {i}: {x} vs gromos++ {y}"
            ),
            _ => panic!("token {i}: '{x}' vs gromos++ '{y}'"),
        }
    }
}

#[test]
fn make_sasa_top_matches_gromospp() {
    let d = data("make_sasa_top");
    let topo = references()
        .join("shared/aladip.topo")
        .to_string_lossy()
        .into_owned();
    let spec = d.join("aladip.sasaspec").to_string_lossy().into_owned();
    for (extra, reference) in [
        (vec![], "aladip_sasa.gromospp.top"),
        (vec!["@noH"], "aladip_sasa_noH.gromospp.top"),
    ] {
        let mut args: Vec<String> = ["@topo", &topo, "@sasaspec", &spec]
            .iter()
            .map(|s| s.to_string())
            .collect();
        args.extend(extra.iter().map(|s| s.to_string()));
        let out = Command::new(env!("CARGO_BIN_EXE_make_sasa_top"))
            .args(&args)
            .output()
            .expect("run make_sasa_top");
        assert!(
            out.status.success(),
            "make_sasa_top failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let ours = String::from_utf8_lossy(&out.stdout).into_owned();
        let theirs = std::fs::read_to_string(d.join(reference)).unwrap();
        assert_same_text(&ours, &theirs, 1e-5);
        // and the topology part reads back
        let tmp = std::env::temp_dir().join(format!(
            "newen_sasa_{}_{}.top",
            std::process::id(),
            extra.len()
        ));
        std::fs::write(&tmp, &ours).unwrap();
        let parsed = read_topology_file(&tmp).unwrap();
        std::fs::remove_file(tmp).ok();
        assert_eq!(parsed.n_atoms, 12);
    }
}
