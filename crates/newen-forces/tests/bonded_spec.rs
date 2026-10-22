//! The bonded kernels against their Lean model, bit for bit, and against their own energy by
//! finite differences.
//!
//! `spec/vectors/bonded.json` holds `NewenSpec.Bonded.quarticBondFloat`, `harmonicBondFloat`,
//! `cosineAngleFloat`, `harmonicAngleFloat`, `dihedralChebyshevFloat`, `dihedralNewFloat` and
//! `improperFloat` evaluated at Lean's `Float` (IEEE 754 binary64) on random and edge inputs,
//! written by `lake exe spec-vectors` (checked in, so this test needs no Lean). Every number is
//! the hex string of its 64 bits. The kernels use `+ − × ÷` and `sqrt`, which IEEE 754 rounds
//! correctly, and `acos`, `cos`, `sin`, which it does not: Lean's `Float` and Rust's `f64` call the
//! same C library for those, so on one machine the Rust must return the same bits as Lean. A NaN
//! matches any NaN (Lean has a single NaN).
//!
//! The finite-difference tests check the physics whatever the bits: on 200 random rows per
//! kernel, the central difference of the kernel's energy along every coordinate of every atom
//! (`h = 1e-6` nm) is minus the force the kernel returns, to 1e-5 of the largest force component
//! of the row plus 1e-6.

use newen_core::math::Vec3;
use newen_forces::bonded::{
    cosine_angle, dihedral_chebyshev, dihedral_new, harmonic_angle, harmonic_bond,
    improper_dihedral, quartic_bond, NEAR_LINEAR, ZERO_DISTANCE, ZERO_DISTANCE2, ZERO_SINE,
};

/// `tinyFloat` (`1e-10`) and `piFloat` in `spec/NewenSpec/Bonded.lean`.
const TINY_BITS: u64 = 0x3ddb_7cdf_d9d7_bdbb;
const PI_BITS: u64 = 0x4009_21fb_5444_2d18;

/// The three sections of `bonded.json`, each a list of rows of floats.
struct Rows {
    bonds: Vec<Vec<f64>>,
    angles: Vec<Vec<f64>>,
    dihedrals: Vec<Vec<f64>>,
}

fn lean_rows() -> Rows {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/bonded.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let mut rows = Rows {
        bonds: Vec::new(),
        angles: Vec::new(),
        dihedrals: Vec::new(),
    };
    let mut section = "";
    for line in text.lines().map(str::trim) {
        for name in ["bonds", "angles", "dihedrals"] {
            if line.starts_with(&format!("\"{name}\"")) {
                section = name;
            }
        }
        let line = line.trim_end_matches(',');
        if !(line.starts_with('[') && line.ends_with(']')) {
            continue;
        }
        let fields: Vec<f64> = line[1..line.len() - 1]
            .split(',')
            .map(|f| float(f.trim().trim_matches('"')))
            .collect();
        match section {
            "bonds" => rows.bonds.push(fields),
            "angles" => rows.angles.push(fields),
            "dihedrals" => rows.dihedrals.push(fields),
            _ => panic!("row outside a section: {line}"),
        }
    }
    rows
}

fn float(field: &str) -> f64 {
    let hex = field.strip_prefix("0x").expect("hex bits in bonded.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in bonded.json"))
}

/// The vector at `row[i..i + 3]`.
fn vec3(row: &[f64], i: usize) -> Vec3 {
    Vec3::new(row[i], row[i + 1], row[i + 2])
}

/// An energy and its force vectors as one list of floats, in the layout of the rows.
fn flat(energy: f64, forces: &[Vec3]) -> Vec<f64> {
    let mut out = vec![energy];
    for f in forces {
        out.extend([f.x, f.y, f.z]);
    }
    out
}

/// The same bits, or both NaN (see the module comment).
fn same(rust: f64, lean: f64) -> bool {
    rust.to_bits() == lean.to_bits() || (rust.is_nan() && lean.is_nan())
}

/// Distance in units in the last place, through the sign-magnitude to two's-complement map.
fn ulps(a: f64, b: f64) -> u128 {
    let ordered = |x: f64| {
        let b = x.to_bits() as i64;
        if b < 0 {
            i64::MIN as i128 - b as i128
        } else {
            b as i128
        }
    };
    (ordered(a) - ordered(b)).unsigned_abs()
}

/// Compare one kernel's Rust outputs with its Lean outputs; on a difference, describe the row.
fn check(inputs: &[f64], rust: &[f64], lean: &[f64], mismatches: &mut Vec<String>) {
    assert_eq!(rust.len(), lean.len());
    if rust.iter().zip(lean).all(|(&r, &l)| same(r, l)) {
        return;
    }
    let worst = rust
        .iter()
        .zip(lean)
        .map(|(&r, &l)| ulps(r, l))
        .max()
        .unwrap_or(0);
    mismatches.push(format!(
        "  inputs={inputs:?}: rust={rust:?} lean={lean:?} ({worst} ulp)"
    ));
}

fn report(what: &str, rows: usize, mismatches: &[String]) {
    assert!(
        mismatches.is_empty(),
        "{what}: {} of {rows} rows differ from the Lean model; the first ones:\n{}",
        mismatches.len(),
        mismatches
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// Row layouts (see `layout_*` in `bonded.json`).
const BOND_INPUTS: usize = 6;
const ANGLE_INPUTS: usize = 10;
const DIHEDRAL_INPUTS: usize = 15;

fn bond_rows() -> Vec<Vec<f64>> {
    let rows = lean_rows().bonds;
    assert!(rows.len() >= 1_000, "{} bonds rows", rows.len());
    assert!(rows.iter().all(|r| r.len() == BOND_INPUTS + 8));
    rows
}

fn angle_rows() -> Vec<Vec<f64>> {
    let rows = lean_rows().angles;
    assert!(rows.len() >= 1_000, "{} angles rows", rows.len());
    assert!(rows.iter().all(|r| r.len() == ANGLE_INPUTS + 20));
    rows
}

fn dihedral_rows() -> Vec<Vec<f64>> {
    let rows = lean_rows().dihedrals;
    assert!(rows.len() >= 700, "{} dihedrals rows", rows.len());
    assert!(rows.iter().all(|r| r.len() == DIHEDRAL_INPUTS + 39));
    rows
}

/// `m` of a dihedrals row, written as a float.
fn multiplicity(row: &[f64]) -> u32 {
    let m = row[12];
    assert!(m >= 0.0 && m.fract() == 0.0, "m = {m}");
    m as u32
}

#[test]
fn constants_are_the_bits_the_lean_model_uses() {
    for (name, c) in [
        ("ZERO_DISTANCE", ZERO_DISTANCE),
        ("ZERO_DISTANCE2", ZERO_DISTANCE2),
        ("ZERO_SINE", ZERO_SINE),
        ("NEAR_LINEAR", NEAR_LINEAR),
    ] {
        assert_eq!(c.to_bits(), TINY_BITS, "{name}");
    }
    assert_eq!(std::f64::consts::PI.to_bits(), PI_BITS);
}

#[test]
fn quartic_bond_equals_lean_model_bit_for_bit() {
    let rows = bond_rows();
    let mut mismatches = Vec::new();
    for row in &rows {
        let (e, f) = quartic_bond(vec3(row, 0), row[3], row[4]);
        check(
            &row[..BOND_INPUTS],
            &flat(e, &[f]),
            &row[6..10],
            &mut mismatches,
        );
    }
    report("quartic_bond", rows.len(), &mismatches);
}

#[test]
fn harmonic_bond_equals_lean_model_bit_for_bit() {
    let rows = bond_rows();
    let mut mismatches = Vec::new();
    for row in &rows {
        let (e, f) = harmonic_bond(vec3(row, 0), row[5], row[4]);
        check(
            &row[..BOND_INPUTS],
            &flat(e, &[f]),
            &row[10..14],
            &mut mismatches,
        );
    }
    report("harmonic_bond", rows.len(), &mismatches);
}

#[test]
fn cosine_angle_equals_lean_model_bit_for_bit() {
    let rows = angle_rows();
    let mut mismatches = Vec::new();
    for row in &rows {
        let (e, fi, fj, fk) = cosine_angle(vec3(row, 0), vec3(row, 3), row[6], row[7]);
        check(
            &row[..ANGLE_INPUTS],
            &flat(e, &[fi, fj, fk]),
            &row[10..20],
            &mut mismatches,
        );
    }
    report("cosine_angle", rows.len(), &mismatches);
}

#[test]
fn harmonic_angle_equals_lean_model_bit_for_bit() {
    let rows = angle_rows();
    let mut mismatches = Vec::new();
    for row in &rows {
        let (e, fi, fj, fk) = harmonic_angle(vec3(row, 0), vec3(row, 3), row[8], row[9]);
        check(
            &row[..ANGLE_INPUTS],
            &flat(e, &[fi, fj, fk]),
            &row[20..30],
            &mut mismatches,
        );
    }
    report("harmonic_angle", rows.len(), &mismatches);
}

#[test]
fn dihedral_chebyshev_equals_lean_model_bit_for_bit() {
    let rows = dihedral_rows();
    let mut mismatches = Vec::new();
    for row in &rows {
        let (e, fi, fj, fk, fl) = dihedral_chebyshev(
            vec3(row, 0),
            vec3(row, 3),
            vec3(row, 6),
            row[9],
            row[10],
            multiplicity(row),
        );
        let rust = flat(e, &[fi, fj, fk, fl]);
        check(
            &row[..DIHEDRAL_INPUTS],
            &rust,
            &row[15..28],
            &mut mismatches,
        );
    }
    report("dihedral_chebyshev", rows.len(), &mismatches);
}

#[test]
fn dihedral_new_equals_lean_model_bit_for_bit() {
    let rows = dihedral_rows();
    let mut mismatches = Vec::new();
    for row in &rows {
        let (e, fi, fj, fk, fl) = dihedral_new(
            vec3(row, 0),
            vec3(row, 3),
            vec3(row, 6),
            row[9],
            row[11],
            row[12],
        );
        let rust = flat(e, &[fi, fj, fk, fl]);
        check(
            &row[..DIHEDRAL_INPUTS],
            &rust,
            &row[28..41],
            &mut mismatches,
        );
    }
    report("dihedral_new", rows.len(), &mismatches);
}

#[test]
fn improper_dihedral_equals_lean_model_bit_for_bit() {
    let rows = dihedral_rows();
    let mut mismatches = Vec::new();
    for row in &rows {
        let (e, fi, fj, fk, fl) =
            improper_dihedral(vec3(row, 0), vec3(row, 3), vec3(row, 6), row[13], row[14]);
        let rust = flat(e, &[fi, fj, fk, fl]);
        check(
            &row[..DIHEDRAL_INPUTS],
            &rust,
            &row[41..54],
            &mut mismatches,
        );
    }
    report("improper_dihedral", rows.len(), &mismatches);
}

/// The generated rows reach every guard and branch of the kernels.
#[test]
fn rows_cover_the_edges() {
    let bonds = bond_rows();
    assert!(bonds.iter().any(|r| vec3(r, 0).length() == 0.0));
    assert!(bonds.iter().any(|r| vec3(r, 0).length() == ZERO_DISTANCE));
    assert!(bonds
        .iter()
        .any(|r| vec3(r, 0).length() < ZERO_DISTANCE && vec3(r, 0).length() > 0.0));

    let angles = angle_rows();
    let angle_cos = |r: &Vec<f64>| {
        let (a, b) = (vec3(r, 0), vec3(r, 3));
        a.dot(b) / (a.length() * b.length())
    };
    assert!(angles.iter().any(|r| vec3(r, 3).length() == 0.0));
    assert!(angles.iter().any(|r| angle_cos(r) >= 1.0 - 1e-15));
    assert!(angles.iter().any(|r| angle_cos(r) <= -1.0 + 1e-15));
    // |sin θ| < ZERO_SINE with θ0 = π (the linear branch) and with θ0 ≠ π (no term)
    let linear = |r: &Vec<f64>| angle_cos(r).clamp(-1.0, 1.0).acos().sin().abs() < ZERO_SINE;
    assert!(angles
        .iter()
        .any(|r| linear(r) && r[9] == std::f64::consts::PI));
    assert!(angles
        .iter()
        .any(|r| linear(r) && r[9] != std::f64::consts::PI));
    // the linear branch at θ ≈ 0 with θ0 = π (energy k(1 + cos θ), as gromosXX)
    assert!(angles
        .iter()
        .any(|r| linear(r) && angle_cos(r) > 0.0 && r[9] == std::f64::consts::PI && r[20] != 0.0));

    let dihedrals = dihedral_rows();
    let normals = |r: &Vec<f64>| {
        let (rij, rkj, rkl) = (vec3(r, 0), vec3(r, 3), vec3(r, 6));
        let (mj, nk) = (rij.cross(rkj), rkj.cross(rkl));
        (
            mj.length_squared(),
            nk.length_squared(),
            rkj.length_squared(),
        )
    };
    assert!(dihedrals.iter().any(|r| vec3(r, 3).length() == 0.0));
    assert!(dihedrals.iter().any(|r| {
        let (_, _, dkj2) = normals(r);
        dkj2 > 0.0 && dkj2 < ZERO_DISTANCE2
    }));
    // a linear bond angle: the NTBDN = 0 kernel drops f_i, the Chebyshev kernel does not
    assert!(dihedrals.iter().any(|r| {
        let (dmj2, _, dkj2) = normals(r);
        dmj2 > 0.0 && dmj2 < NEAR_LINEAR * dkj2 && r[29..32].iter().all(|&f| f == 0.0)
    }));
    // j, k, l nearly collinear with i off the axis: both kernels drop f_l and keep f_i
    let zero = |f: &[f64]| f.iter().all(|&x| x == 0.0);
    let l_side = |r: &Vec<f64>| {
        let (dmj2, dnk2, dkj2) = normals(r);
        dnk2.sqrt() >= ZERO_DISTANCE2 && dnk2 < NEAR_LINEAR * dkj2 && dmj2 >= NEAR_LINEAR * dkj2
    };
    assert!(dihedrals
        .iter()
        .any(|r| l_side(r) && zero(&r[38..41]) && !zero(&r[29..32])));
    assert!(dihedrals
        .iter()
        .any(|r| l_side(r) && zero(&r[51..54]) && !zero(&r[42..45])));
    // ζ0 outside [−π, π], and the wrap of ζ into [ζ0 − π, ζ0 + π] turning twice
    let pi = std::f64::consts::PI;
    let raw_zeta = |r: &Vec<f64>| {
        let (rij, rkj, rkl) = (vec3(r, 0), vec3(r, 3), vec3(r, 6));
        let (mj, nk) = (rij.cross(rkj), rkj.cross(rkl));
        let z = (mj.dot(nk) / (mj.length() * nk.length()))
            .clamp(-1.0, 1.0)
            .acos();
        if rij.dot(nk) < 0.0 {
            -z
        } else {
            z
        }
    };
    assert!(dihedrals.iter().any(|r| r[14].abs() > pi && r[41] != 0.0));
    assert!(dihedrals
        .iter()
        .any(|r| (raw_zeta(r) - r[14]).abs() > 3.0 * pi && r[41] != 0.0));
    assert!(dihedrals.iter().any(|r| multiplicity(r) == 7));
    for m in 0..=6 {
        assert!(dihedrals.iter().any(|r| multiplicity(r) == m), "m = {m}");
    }
}

// ─── Finite differences ───────────────────────────────────────────────────────

/// Central-difference step (nm).
const H: f64 = 1e-6;

/// Rows checked per kernel.
const FD_ROWS: usize = 200;

/// For `N` atoms at `x`, every coordinate's central difference of `energy` against `−force`.
/// Returns a description of the worst coordinate when it misses the tolerance.
fn fd_mismatch<const N: usize>(
    x: [Vec3; N],
    energy: impl Fn(&[Vec3; N]) -> f64,
    force: [Vec3; N],
) -> Option<String> {
    let scale = force
        .iter()
        .flat_map(|f| [f.x.abs(), f.y.abs(), f.z.abs()])
        .fold(0.0, f64::max);
    let tol = 1e-5 * scale + 1e-6;
    for a in 0..N {
        for c in 0..3 {
            let shifted = |h: f64| {
                let mut y = x;
                y[a][c] += h;
                energy(&y)
            };
            let fd = -(shifted(H) - shifted(-H)) / (2.0 * H);
            let err = (fd - force[a][c]).abs();
            // a NaN error fails too
            if err.is_nan() || err > tol {
                return Some(format!(
                    "atom {a} coordinate {c}: -dE/dx = {fd:e}, force {:e} (error {err:e}, tolerance {tol:e})",
                    force[a][c]
                ));
            }
        }
    }
    None
}

/// Run the finite-difference check on the first `FD_ROWS` rows `select` accepts.
fn fd_check<const N: usize>(
    what: &str,
    rows: &[Vec<f64>],
    select: impl Fn(&[f64]) -> bool,
    case: impl Fn(&[f64]) -> Option<String>,
) {
    let chosen: Vec<&Vec<f64>> = rows.iter().filter(|r| select(r)).take(FD_ROWS).collect();
    assert_eq!(
        chosen.len(),
        FD_ROWS,
        "{what}: too few well-conditioned rows"
    );
    let failures: Vec<String> = chosen
        .iter()
        .filter_map(|r| case(r).map(|m| format!("  row={r:?}: {m}")))
        .collect();
    assert!(
        failures.is_empty(),
        "{what}: {} of {FD_ROWS} rows fail the finite-difference check:\n{}",
        failures.len(),
        failures
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// `|sin|` of the angle between two vectors.
fn sine(a: Vec3, b: Vec3) -> f64 {
    a.cross(b).length() / (a.length() * b.length())
}

/// Atoms `i, j` with `r = x_i − x_j`, `x_j = 0`.
fn bond_atoms(row: &[f64]) -> [Vec3; 2] {
    [vec3(row, 0), Vec3::ZERO]
}

/// Atoms `i, j, k` with `x_j = 0`.
fn angle_atoms(row: &[f64]) -> [Vec3; 3] {
    [vec3(row, 0), Vec3::ZERO, vec3(row, 3)]
}

/// Atoms `i, j, k, l` with `x_j = 0`, `x_l = x_k − r_kl`.
fn dihedral_atoms(row: &[f64]) -> [Vec3; 4] {
    [
        vec3(row, 0),
        Vec3::ZERO,
        vec3(row, 3),
        vec3(row, 3) - vec3(row, 6),
    ]
}

/// Bond angles of the dihedral at least ~6° away from 0 and π.
fn well_conditioned_dihedral(row: &[f64]) -> bool {
    let (rij, rkj, rkl) = (vec3(row, 0), vec3(row, 3), vec3(row, 6));
    sine(rij, rkj) > 0.1 && sine(rkj, rkl) > 0.1
}

#[test]
fn quartic_bond_force_is_minus_the_energy_gradient() {
    fd_check::<2>(
        "quartic_bond",
        &bond_rows(),
        |_| true,
        |row| {
            let (k, b0) = (row[3], row[4]);
            let x = bond_atoms(row);
            let (_, f) = quartic_bond(x[0] - x[1], k, b0);
            fd_mismatch(x, |y| quartic_bond(y[0] - y[1], k, b0).0, [f, -f])
        },
    );
}

#[test]
fn harmonic_bond_force_is_minus_the_energy_gradient() {
    fd_check::<2>(
        "harmonic_bond",
        &bond_rows(),
        |_| true,
        |row| {
            let (k, b0) = (row[5], row[4]);
            let x = bond_atoms(row);
            let (_, f) = harmonic_bond(x[0] - x[1], k, b0);
            fd_mismatch(x, |y| harmonic_bond(y[0] - y[1], k, b0).0, [f, -f])
        },
    );
}

#[test]
fn cosine_angle_force_is_minus_the_energy_gradient() {
    let select = |row: &[f64]| sine(vec3(row, 0), vec3(row, 3)) > 0.05;
    fd_check::<3>("cosine_angle", &angle_rows(), select, |row| {
        let (k, cos0) = (row[6], row[7]);
        let x = angle_atoms(row);
        let e = |y: &[Vec3; 3]| cosine_angle(y[0] - y[1], y[2] - y[1], k, cos0);
        let (_, fi, fj, fk) = e(&x);
        fd_mismatch(x, |y| e(y).0, [fi, fj, fk])
    });
}

#[test]
fn harmonic_angle_force_is_minus_the_energy_gradient() {
    let select = |row: &[f64]| sine(vec3(row, 0), vec3(row, 3)) > 0.05;
    fd_check::<3>("harmonic_angle", &angle_rows(), select, |row| {
        let (k, theta0) = (row[8], row[9]);
        let x = angle_atoms(row);
        let e = |y: &[Vec3; 3]| harmonic_angle(y[0] - y[1], y[2] - y[1], k, theta0);
        let (_, fi, fj, fk) = e(&x);
        fd_mismatch(x, |y| e(y).0, [fi, fj, fk])
    });
}

#[test]
fn dihedral_chebyshev_force_is_minus_the_energy_gradient() {
    fd_check::<4>(
        "dihedral_chebyshev",
        &dihedral_rows(),
        well_conditioned_dihedral,
        |row| {
            let (k, cos_delta, m) = (row[9], row[10], multiplicity(row));
            let x = dihedral_atoms(row);
            let e = |y: &[Vec3; 4]| {
                dihedral_chebyshev(y[0] - y[1], y[2] - y[1], y[2] - y[3], k, cos_delta, m)
            };
            let (_, fi, fj, fk, fl) = e(&x);
            fd_mismatch(x, |y| e(y).0, [fi, fj, fk, fl])
        },
    );
}

#[test]
fn dihedral_new_force_is_minus_the_energy_gradient() {
    fd_check::<4>(
        "dihedral_new",
        &dihedral_rows(),
        well_conditioned_dihedral,
        |row| {
            let (k, delta, m) = (row[9], row[11], row[12]);
            let x = dihedral_atoms(row);
            let e =
                |y: &[Vec3; 4]| dihedral_new(y[0] - y[1], y[2] - y[1], y[2] - y[3], k, delta, m);
            let (_, fi, fj, fk, fl) = e(&x);
            fd_mismatch(x, |y| e(y).0, [fi, fj, fk, fl])
        },
    );
}

#[test]
fn improper_dihedral_force_is_minus_the_energy_gradient() {
    // away from ζ = ζ0 ± π, where the wrap into [ζ0 − π, ζ0 + π] makes the energy jump
    let select = |row: &[f64]| {
        let e = improper_dihedral(vec3(row, 0), vec3(row, 3), vec3(row, 6), row[13], row[14]).0;
        let d_zeta = (2.0 * e / row[13]).sqrt();
        well_conditioned_dihedral(row) && d_zeta < std::f64::consts::PI - 0.05
    };
    fd_check::<4>("improper_dihedral", &dihedral_rows(), select, |row| {
        let (k, zeta0) = (row[13], row[14]);
        let x = dihedral_atoms(row);
        let e = |y: &[Vec3; 4]| improper_dihedral(y[0] - y[1], y[2] - y[1], y[2] - y[3], k, zeta0);
        let (_, fi, fj, fk, fl) = e(&x);
        fd_mismatch(x, |y| e(y).0, [fi, fj, fk, fl])
    });
}
