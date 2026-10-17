//! A *null* perturbation must be a no-op — the invariance that pinned every FEP defect of
//! 2026-08-29/30 (1.7). gromosXX has it by construction; here it is explicit.
//!
//! A `.ptp` whose atoms have identical A and B states and zero soft-core parameters describes
//! the unperturbed system: energies and forces must be bit-identical to a run without any
//! perturbation topology, for a perturbed atom anywhere in the topology (the pairlist routing
//! bug only showed for atoms that were not the first). With a non-zero α the same `.ptp` must
//! change the energy: soft-core is not a no-op even at equal end states. A `LAMBDAS` block that
//! is on but gives every interaction kind the identity polynomial (`L = λ`, `L′ = 1`) must leave
//! the null perturbation bit-identical too.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use newen_io::read_energy_trajectory_native;

fn refs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

fn md_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_md"))
}

/// PERTATOMPARAM for the alanine dipeptide: (atom, name, iac, mass, charge) as in the topology.
const ATOMS: &[(usize, &str, usize, f64, f64)] = &[
    (1, "CB", 14, 15.0350, 0.000),
    (4, "N", 5, 14.0067, -0.280),
    (11, "H", 18, 1.0080, 0.280),
];

fn ptp(atoms: &[(usize, &str, usize, f64, f64)], alpha_lj: f64, alpha_crf: f64) -> String {
    ptp_with_pairs(atoms, alpha_lj, alpha_crf, &[])
}

/// [`ptp`] with `PERTATOMPAIR` rows `(i, j, A type, B type)`, 1-based, types as in the file.
fn ptp_with_pairs(
    atoms: &[(usize, &str, usize, f64, f64)],
    alpha_lj: f64,
    alpha_crf: f64,
    pairs: &[(usize, usize, u8, u8)],
) -> String {
    let mut s = String::from("TITLE\nnull perturbation\nEND\nPERTATOMPARAM\n");
    s += &format!("{}\n", atoms.len());
    for (nr, name, iac, mass, q) in atoms {
        s += &format!(
            "{nr:5} 1 {name:<4} {iac:3} {mass:9.4} {q:8.3} {iac:3} {mass:9.4} {q:8.3} {alpha_lj:5.2} {alpha_crf:5.2}\n"
        );
    }
    s += "END\n";
    s += &format!("PERTATOMPAIR\n{}\n", pairs.len());
    for (i, j, a, b) in pairs {
        s += &format!("{i:5} {j:5} {a} {b}\n");
    }
    s += "END\n";
    for block in [
        "PERTBONDSTRETCH",
        "PERTBONDANGLE",
        "PERTIMPROPERDIH",
        "PERTPROPERDIH",
    ] {
        s += &format!("{block}\n   0\nEND\n");
    }
    s
}

fn run(pttopo: Option<&Path>, tag: &str) -> (Vec<f64>, String) {
    run_input(
        pttopo,
        &refs().join("aladip_vacuum_fep/aladip_vacuum_fep.in"),
        tag,
    )
}

fn run_input(pttopo: Option<&Path>, input: &Path, tag: &str) -> (Vec<f64>, String) {
    run_with(pttopo, input, tag, &[])
}

fn run_with(pttopo: Option<&Path>, input: &Path, tag: &str, extra: &[&str]) -> (Vec<f64>, String) {
    let (e, forces, _) = run_logged(pttopo, input, tag, extra);
    (e, forces)
}

/// [`run_with`], plus what `md` printed (stdout, then stderr).
fn run_logged(
    pttopo: Option<&Path>,
    input: &Path,
    tag: &str,
    extra: &[&str],
) -> (Vec<f64>, String, String) {
    let r = refs();
    let out = std::env::temp_dir().join(format!("newen_fep_null_{tag}_{}", std::process::id()));
    fs::create_dir_all(&out).unwrap();
    let mut cmd = Command::new(md_bin());
    cmd.arg("@topo")
        .arg(r.join("shared/aladip.topo"))
        .arg("@conf")
        .arg(r.join("aladip_vacuum/aladip_vacuum.conf"))
        .arg("@input")
        .arg(input)
        .arg("@fin")
        .arg(out.join("fin.cnf"))
        .arg("@tre")
        .arg(out.join("e.tre"))
        .arg("@trc")
        .arg(out.join("t.trc"))
        .arg("@trf")
        .arg(out.join("f.trf"))
        .arg("@trg")
        .arg(out.join("g.trg"));
    if let Some(p) = pttopo {
        cmd.arg("@pttopo").arg(p);
    }
    cmd.args(extra);
    let status = cmd.output().expect("run md");
    assert!(
        status.status.success(),
        "md failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    let frames = read_energy_trajectory_native(out.join("e.tre")).unwrap();
    let f = &frames[0];
    let forces = fs::read_to_string(out.join("f.trf")).unwrap();
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );
    (
        vec![
            f.total,
            f.kinetic,
            f.potential,
            f.bond,
            f.angle,
            f.improper,
            f.dihedral,
            f.lj,
            f.coul_real,
        ],
        forces,
        log,
    )
}

#[test]
fn null_perturbation_is_a_no_op_and_soft_core_is_not() {
    let dir = std::env::temp_dir().join(format!("newen_fep_null_ptp_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let (reference, ref_forces) = run(None, "none");

    // Identical end states, α = 0: nothing may change — energies and the force trajectory.
    let null = dir.join("null.ptp");
    fs::write(&null, ptp(ATOMS, 0.0, 0.0)).unwrap();
    let (e, forces) = run(Some(&null), "null");
    assert_eq!(
        e, reference,
        "a null perturbation changed the energies: {e:?} vs {reference:?}"
    );
    assert_eq!(forces, ref_forces, "a null perturbation changed the forces");

    // Identical end states, α ≠ 0: the soft-core cutoff/distance terms act — energy must move.
    let soft = dir.join("soft.ptp");
    fs::write(&soft, ptp(ATOMS, 0.3, 0.6)).unwrap();
    let (e, _) = run(Some(&soft), "soft");
    assert_ne!(
        e[8], reference[8],
        "soft-core CRF with α_crf = 0.6 left the energy unchanged"
    );
    assert_ne!(
        e[7], reference[7],
        "soft-core LJ with α_lj = 0.3 left the energy unchanged"
    );
}

/// A null perturbation keeps the net charge: `md` reports ΔQ = 0 and does not warn about a
/// finite-size correction.
#[test]
fn null_perturbation_reports_no_net_charge_change() {
    let dir = std::env::temp_dir().join(format!("newen_fep_null_dq_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let null = dir.join("null.ptp");
    fs::write(&null, ptp(ATOMS, 0.0, 0.0)).unwrap();
    let input = refs().join("aladip_vacuum_fep/aladip_vacuum_fep.in");
    let (_, _, log) = run_logged(Some(&null), &input, "dq", &[]);
    assert!(
        log.contains("Perturbation net charge: A 0.000 B 0.000 (ΔQ = 0.000)"),
        "{log}"
    );
    assert!(!log.contains("finite-size correction"), "{log}");
}

/// `aladip_vacuum_fep.in` with a `LAMBDAS` block that is on and sets every kind gromosXX lets a
/// row set to the identity polynomial `0 0 0 1 0`.
fn identity_lambdas_input(dir: &Path) -> PathBuf {
    let text = fs::read_to_string(refs().join("aladip_vacuum_fep/aladip_vacuum_fep.in")).unwrap();
    let mut block = String::from("LAMBDAS\n# NTIL\n  on\n# NTLI NILG1 NILG2 ALI BLI CLI DLI ELI\n");
    for kind in [
        "bond",
        "angle",
        "dihedral",
        "improper",
        "vdw",
        "vdw_soft",
        "crf",
        "crf_soft",
        "distanceres",
        "distancefield",
        "dihedralres",
        "mass",
    ] {
        block += &format!("  {kind} 1 1 0 0 0 1 0\n");
    }
    block += "END\n";
    let path = dir.join("identity_lambdas.in");
    fs::write(&path, text + &block).unwrap();
    path
}

#[test]
fn identity_individual_lambdas_keep_the_null_perturbation_a_no_op() {
    let dir = std::env::temp_dir().join(format!("newen_fep_null_lam_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let (reference, ref_forces) = run(None, "lam_none");
    let null = dir.join("null.ptp");
    fs::write(&null, ptp(ATOMS, 0.0, 0.0)).unwrap();
    let input = identity_lambdas_input(&dir);
    let (e, forces) = run_input(Some(&null), &input, "lam_null");
    assert_eq!(
        e, reference,
        "identity LAMBDAS on a null perturbation changed the energies: {e:?} vs {reference:?}"
    );
    assert_eq!(
        forces, ref_forces,
        "identity LAMBDAS on a null perturbation changed the forces"
    );
}

/// The numbers of a force trajectory, in order.
fn numbers(text: &str) -> Vec<f64> {
    text.split_whitespace()
        .filter_map(|t| t.parse::<f64>().ok())
        .collect()
}

/// `PERTATOMPAIR` rows on a perturbed atom whose A and B types are the pair's topology type
/// (4-10 a 1-4 pair, 4-11 a normal one, both charged) take the pair out of the regular loops
/// into the atom-pair term: with equal end states and α = 0 that is still the unperturbed
/// system, to rounding, under the intended kernel. gromosXX's shifted call is not.
#[test]
fn pertatompair_rows_with_equal_states_keep_the_null_perturbation_a_no_op() {
    let dir = std::env::temp_dir().join(format!("newen_fep_null_pair_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let input = refs().join("aladip_vacuum_fep/aladip_vacuum_fep.in");
    let (reference, ref_forces) = run(None, "pair_none");
    let null = dir.join("null_pairs.ptp");
    fs::write(
        &null,
        ptp_with_pairs(ATOMS, 0.0, 0.0, &[(4, 10, 2, 2), (4, 11, 1, 1)]),
    )
    .unwrap();

    let (e, forces) = run_with(Some(&null), &input, "pair_null", &[]);
    for (k, (a, b)) in e.iter().zip(&reference).enumerate() {
        assert!(
            (a - b).abs() <= 1e-10 * b.abs().max(1.0),
            "energy {k}: {a} vs {b} with null PERTATOMPAIR rows"
        );
    }
    let (f, f_ref) = (numbers(&forces), numbers(&ref_forces));
    assert_eq!(f.len(), f_ref.len());
    for (a, b) in f.iter().zip(&f_ref) {
        assert!(
            (a - b).abs() <= 1e-8,
            "force {a} vs {b} with null PERTATOMPAIR rows"
        );
    }

    let (shifted, _) = run_with(
        Some(&null),
        &input,
        "pair_null_shifted",
        &[
            "@recipe_set",
            "perturbation.atom_pair_kernel=gromosxx-shifted",
        ],
    );
    assert!(
        (shifted[8] - reference[8]).abs() > 1e-3,
        "the shifted kernel left the Coulomb energy unchanged"
    );
}
