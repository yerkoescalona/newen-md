//! The perturbed bonded kernels against their Lean model, bit for bit, and their λ derivatives
//! and forces against finite differences of their own energies.
//!
//! `spec/vectors/perturbed_bonded.json` holds the `NewenSpec.PerturbedBonded` models evaluated
//! at Lean's `Float` on random inputs, written by `lake exe spec-vectors`. Every float is the hex
//! string of its 64 bits; multiplicities are plain integers. The kernels use only `+ − × ÷`
//! (`powi` is a fixed sequence of multiplications), so the Rust must return the same bits as
//! Lean. A NaN matches any NaN.

use newen_forces::bonded::perturbed::{
    chebyshev, perturbed_angle_term, perturbed_dihedral_term, perturbed_improper_term,
    perturbed_quartic_bond_term, soft_perturbed_term,
};

const SECTIONS: [&str; 6] = [
    "chebyshev",
    "quartic",
    "angle",
    "improper",
    "dihedral",
    "soft",
];

fn section(name: &str) -> Vec<Vec<f64>> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/perturbed_bonded.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let mut current = "";
    let mut rows = Vec::new();
    for line in text.lines().map(str::trim) {
        if let Some(s) = SECTIONS
            .iter()
            .find(|s| line.starts_with(&format!("\"{s}\": [")))
        {
            current = s;
        }
        let line = line.trim_end_matches(',');
        if current != name || !(line.starts_with('[') && line.ends_with(']')) {
            continue;
        }
        rows.push(
            line[1..line.len() - 1]
                .split(',')
                .map(|f| field(f.trim()))
                .collect(),
        );
    }
    assert!(rows.len() >= 200, "{} rows in section {name}", rows.len());
    rows
}

fn field(f: &str) -> f64 {
    match f.trim_matches('"').strip_prefix("0x") {
        Some(hex) => f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex")),
        None => f.parse::<i32>().expect("integer field") as f64,
    }
}

fn same(rust: f64, lean: f64) -> bool {
    rust.to_bits() == lean.to_bits() || (rust.is_nan() && lean.is_nan())
}

fn assert_rows(name: &str, kernel: impl Fn(&[f64]) -> Vec<f64>, n_out: usize) {
    let rows = section(name);
    let mismatches: Vec<String> = rows
        .iter()
        .filter_map(|row| {
            let (inputs, lean) = row.split_at(row.len() - n_out);
            let rust = kernel(inputs);
            (!rust.iter().zip(lean).all(|(&r, &l)| same(r, l)))
                .then(|| format!("  inputs={inputs:?}: rust={rust:?} lean={lean:?}"))
        })
        .collect();
    assert!(
        mismatches.is_empty(),
        "{name}: {} of {} rows differ from the Lean model; the first ones:\n{}",
        mismatches.len(),
        rows.len(),
        mismatches[..mismatches.len().min(20)].join("\n")
    );
}

#[test]
fn chebyshev_equals_lean_model_bit_for_bit() {
    assert_rows(
        "chebyshev",
        |r| {
            let (c, d) = chebyshev(r[0] as i32, r[1]);
            vec![c, d]
        },
        2,
    );
}

#[test]
fn quartic_bond_equals_lean_model_bit_for_bit() {
    assert_rows(
        "quartic",
        |r| {
            let (e, f, de) = perturbed_quartic_bond_term(r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
            vec![e, f, de]
        },
        3,
    );
}

#[test]
fn angle_equals_lean_model_bit_for_bit() {
    assert_rows(
        "angle",
        |r| {
            let (e, df, de) = perturbed_angle_term(r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
            vec![e, df, de]
        },
        3,
    );
}

#[test]
fn improper_equals_lean_model_bit_for_bit() {
    assert_rows(
        "improper",
        |r| {
            let (e, mkd, de, z) = perturbed_improper_term(r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
            vec![e, mkd, de, z]
        },
        4,
    );
}

fn dihedral(r: &[f64]) -> newen_forces::bonded::perturbed::PerturbedDihedralTerm {
    perturbed_dihedral_term(
        r[0],
        (r[1], r[2], r[3] as i32),
        (r[4], r[5], r[6] as i32),
        r[7],
        r[8],
    )
}

#[test]
fn dihedral_equals_lean_model_bit_for_bit() {
    assert_rows(
        "dihedral",
        |r| {
            let t = dihedral(r);
            vec![t.e_a, t.e_b, t.g_a, t.g_b, t.energy, t.de_dl]
        },
        6,
    );
}

#[test]
fn soft_term_equals_lean_model_bit_for_bit() {
    assert_rows(
        "soft",
        |r| {
            let (e, g, de) = soft_perturbed_term(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
            vec![e, g, de]
        },
        3,
    );
}

/// The generated rows cover every multiplicity, both end states, absent states, several turns.
#[test]
fn rows_cover_the_edges() {
    let cheb = section("chebyshev");
    for m in 0..=7 {
        assert!(cheb.iter().any(|r| r[0] == m as f64));
    }
    for name in ["quartic", "angle", "improper", "soft"] {
        let rows = section(name);
        let lam = if name == "soft" { 6 } else { 5 };
        assert!(rows.iter().any(|r| r[lam] == 0.0), "{name}: λ = 0");
        assert!(rows.iter().any(|r| r[lam] == 1.0), "{name}: λ = 1");
        assert!(rows.iter().any(|r| r[lam + 1] != 1.0), "{name}: λ′ ≠ 1");
        assert!(
            rows.iter().any(|r| r[1] == 0.0 || r[2] == 0.0),
            "{name}: K = 0"
        );
    }
    let imp = section("improper");
    assert!(imp.iter().any(|r| r[0].abs() > 2.0 * std::f64::consts::PI));
}

const H: f64 = 1e-6;

/// `|a − b| ≤ 1e-6·|a| + 1e-9·scale`.
fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-6 * a.abs() + 1e-9 * scale
}

/// Soft rows whose softening `α·Δ²` is at most 10: beyond it the energy varies on a λ scale
/// so short that the truncation error of a central difference with `h = 1e-6` exceeds 1e-6.
fn soft_well_conditioned(r: &[f64]) -> bool {
    let d = r[0] - ((1.0 - r[6]) * r[3] + r[6] * r[4]);
    r[5] * d * d <= 10.0
}

/// A kernel's `(energy, dE/dλ)` at a row whose λ is moved by `δ·λ′`.
fn check_dvdl(name: &str, lam: usize, energy_de: impl Fn(&[f64]) -> (f64, f64)) {
    let rows: Vec<_> = section(name)
        .into_iter()
        .filter(|r| name != "soft" || soft_well_conditioned(r))
        .collect();
    assert!(rows.len() >= 150, "{name}: {} rows", rows.len());
    for row in rows {
        let (e, de) = energy_de(&row);
        let mut r = row.clone();
        r[lam] = row[lam] + H * row[lam + 1];
        let ep = energy_de(&r).0;
        r[lam] = row[lam] - H * row[lam + 1];
        let em = energy_de(&r).0;
        let fd = (ep - em) / (2.0 * H);
        assert!(
            close(de, fd, 1.0 + e.abs()),
            "{name} row={row:?}: de={de:e} fd={fd:e}"
        );
    }
}

/// dE/dλ of every bonded kernel equals the central difference of its energy, with λ moving
/// along λ′ (the chain rule of individual λ).
#[test]
fn bonded_dvdl_is_the_derivative_of_the_energy() {
    check_dvdl("quartic", 5, |r| {
        let (e, _, de) = perturbed_quartic_bond_term(r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
        (e, de)
    });
    check_dvdl("angle", 5, |r| {
        let (e, _, de) = perturbed_angle_term(r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
        (e, de)
    });
    check_dvdl("dihedral", 7, |r| {
        let t = dihedral(r);
        (t.energy, t.de_dl)
    });
    check_dvdl("soft", 6, |r| {
        let (e, _, de) = soft_perturbed_term(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
        (e, de)
    });
    // the wrapped ζ is a step function of λ; check only rows where ζ stays inside its window
    let rows: Vec<_> = section("improper")
        .into_iter()
        .filter(|r| {
            let z = |l: f64| perturbed_improper_term(r[0], r[1], r[2], r[3], r[4], l, r[6]).3;
            z(r[5] + H * r[6]) == z(r[5] - H * r[6])
        })
        .collect();
    assert!(rows.len() >= 250);
    for row in rows {
        let at =
            |l: f64| perturbed_improper_term(row[0], row[1], row[2], row[3], row[4], l, row[6]);
        let (e, _, de, _) = at(row[5]);
        let fd = (at(row[5] + H * row[6]).0 - at(row[5] - H * row[6]).0) / (2.0 * H);
        assert!(
            close(de, fd, 1.0 + e.abs()),
            "improper row={row:?}: de={de:e} fd={fd:e}"
        );
    }
}

/// The force factors are minus the derivatives of the energies in the coordinate.
#[test]
fn bonded_forces_are_minus_the_gradient() {
    for row in section("quartic") {
        let b = row[0].sqrt();
        let e = |x: f64| {
            perturbed_quartic_bond_term(x * x, row[1], row[2], row[3], row[4], row[5], row[6])
        };
        let (en, f, _) = e(b);
        let fd = (e(b + H).0 - e(b - H).0) / (2.0 * H);
        assert!(
            close(-f * b, fd, 1.0 + en.abs()),
            "quartic row={row:?}: -f·b={:e} fd={fd:e}",
            -f * b
        );
    }
    for row in section("angle") {
        let e = |x: f64| perturbed_angle_term(x, row[1], row[2], row[3], row[4], row[5], row[6]);
        let (en, df, _) = e(row[0]);
        let fd = (e(row[0] + H).0 - e(row[0] - H).0) / (2.0 * H);
        assert!(
            close(-df, fd, 1.0 + en.abs()),
            "angle row={row:?}: -df={:e} fd={fd:e}",
            -df
        );
    }
    for row in section("soft")
        .into_iter()
        .filter(|r| soft_well_conditioned(r))
    {
        let e =
            |x: f64| soft_perturbed_term(x, row[1], row[2], row[3], row[4], row[5], row[6], row[7]);
        let (en, g, _) = e(row[0]);
        let fd = (e(row[0] + H).0 - e(row[0] - H).0) / (2.0 * H);
        assert!(
            close(g, fd, 1.0 + en.abs()),
            "soft row={row:?}: g={g:e} fd={fd:e}"
        );
    }
}

/// `chebyshev(m, cos φ) = (cos(mφ), −sin(mφ)·m / (−sin φ))` for `m = 1 … 6`.
#[test]
fn chebyshev_is_cos_m_phi() {
    for i in 0..200 {
        let phi = -3.1 + 6.2 * (i as f64 + 0.5) / 200.0;
        for m in 1..=6 {
            let (c, d) = chebyshev(m, phi.cos());
            let mf = m as f64;
            assert!((c - (mf * phi).cos()).abs() < 1e-12, "m={m} φ={phi}");
            if phi.sin().abs() > 0.1 {
                let dd = mf * (mf * phi).sin() / phi.sin();
                assert!((d - dd).abs() < 1e-10 * (1.0 + dd.abs()), "m={m} φ={phi}");
            }
        }
    }
}
