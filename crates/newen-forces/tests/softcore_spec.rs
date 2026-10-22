//! The perturbed (soft-core) nonbonded kernels against their Lean model, bit for bit, and their
//! λ derivatives and forces against finite differences of their own energies.
//!
//! `spec/vectors/softcore.json` holds the `NewenSpec.SoftCore` models evaluated at Lean's `Float`
//! (IEEE 754 binary64) on random inputs, written by `lake exe spec-vectors` (checked in, so this
//! test needs no Lean). Every float is the hex string of its 64 bits; NLAM is a plain integer.
//! The kernels use `+ − × ÷` and `sqrt`, which IEEE 754 rounds correctly in Rust and in Lean (C)
//! alike, and the λ table calls the C library's `pow` on both sides (Rust's `powf` and Lean's
//! `Float.pow` both link the platform libm), so the Rust must return the same bits as Lean. A NaN
//! matches any NaN (Lean has a single NaN).
//!
//! Each row carries the four inputs of `CRFParameters::new` and the nine inputs of
//! `PerturbedLambdaParams::from_lambda`; the Rust rebuilds both from them.

use newen_forces::nonbonded::{
    lj_crf_soft_interaction_ext, perturbed_lj_crf_pair, perturbed_rf_excluded_pair,
    perturbed_rf_self_term, rf_soft_interaction, rf_soft_interaction_ext, CRFParameters,
    PerturbedLambdaParams,
};

/// `guardFloat` in `spec/NewenSpec/SoftCore.lean`: the `1e-10` of the pair kernel's guards.
const GUARD_BITS: u64 = 0x3ddb_7cdf_d9d7_bdbb;

const SECTIONS: [&str; 7] = [
    "lambda_params",
    "perturbed_lj_crf_pair",
    "rf_soft_interaction",
    "perturbed_rf_excluded_pair",
    "perturbed_rf_self_term",
    "lj_crf_soft_interaction_ext",
    "rf_soft_interaction_ext",
];

/// The rows of one section of `softcore.json`; NLAM is stored as an `f64`.
fn section(name: &str) -> Vec<Vec<f64>> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/softcore.json"
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
    assert!(!rows.is_empty(), "no rows in section {name}");
    rows
}

fn field(f: &str) -> f64 {
    match f.trim_matches('"').strip_prefix("0x") {
        Some(hex) => f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex")),
        None => f.parse::<i32>().expect("integer field") as f64,
    }
}

fn crf_of(inputs: &[f64]) -> CRFParameters {
    CRFParameters::new(inputs[0], inputs[1], inputs[2], inputs[3])
}

/// `from_lambda` of the nine table fields `lj, ljs, crf, crfs, dlj, dljs, dcrf, dcrfs, nlam`.
fn table_of(t: &[f64]) -> PerturbedLambdaParams {
    PerturbedLambdaParams::from_lambda(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8] as i32)
}

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

fn check(row: &[f64], rust: &[f64], lean: &[f64], mismatches: &mut Vec<String>) {
    if rust.iter().zip(lean).all(|(&r, &l)| same(r, l)) {
        return;
    }
    let worst = rust
        .iter()
        .zip(lean)
        .map(|(&r, &l)| ulps(r, l))
        .max()
        .unwrap_or(0);
    let inputs = &row[..row.len() - lean.len()];
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

fn params_fields(p: &PerturbedLambdaParams) -> [f64; 16] {
    [
        p.a_lj_lambda_n,
        p.b_lj_lambda_n,
        p.a_crf_lambda_n,
        p.b_crf_lambda_n,
        p.a_lj_lambda_n_1,
        p.b_lj_lambda_n_1,
        p.a_crf_lambda_n_1,
        p.b_crf_lambda_n_1,
        p.a_ljs_lambda,
        p.b_ljs_lambda,
        p.a_ljs_lambda2,
        p.b_ljs_lambda2,
        p.a_crfs_lambda,
        p.b_crfs_lambda,
        p.a_crfs_lambda2,
        p.b_crfs_lambda2,
    ]
}

/// One row of `perturbed_lj_crf_pair` with the table replaced by `t`.
fn pair_with(row: &[f64], t: &PerturbedLambdaParams, r2: f64) -> (f64, f64, f64, f64, f64) {
    perturbed_lj_crf_pair(
        r2,
        row[1],
        row[2],
        row[3],
        row[4],
        row[5],
        row[6],
        row[7],
        row[8],
        t,
        &crf_of(&row[9..13]),
        row[13],
        row[14],
    )
}

#[test]
fn guard_is_the_bits_the_lean_model_uses() {
    assert_eq!(1e-10_f64.to_bits(), GUARD_BITS);
}

#[test]
fn lambda_params_equal_lean_model_bit_for_bit() {
    let rows = section("lambda_params");
    assert!(rows.len() >= 200, "{} lambda_params rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 25, "lambda_params row {row:?}");
        let p = table_of(&row[..9]);
        assert_eq!(p.lambda_exp, row[8] as i32);
        check(row, &params_fields(&p), &row[9..], &mut mismatches);
    }
    report(
        "PerturbedLambdaParams::from_lambda",
        rows.len(),
        &mismatches,
    );
}

#[test]
fn perturbed_lj_crf_pair_equals_lean_model_bit_for_bit() {
    let rows = section("perturbed_lj_crf_pair");
    assert!(rows.len() >= 500, "{} pair rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 29, "pair row {row:?}");
        let (f, e_lj, e_crf, de_lj, de_crf) = pair_with(row, &table_of(&row[15..24]), row[0]);
        check(
            row,
            &[f, e_lj, e_crf, de_lj, de_crf],
            &row[24..],
            &mut mismatches,
        );
    }
    report("perturbed_lj_crf_pair", rows.len(), &mismatches);
}

type RfKernel =
    fn(f64, f64, f64, f64, &PerturbedLambdaParams, &CRFParameters, f64) -> (f64, f64, f64);

fn rf_rows_equal_lean(name: &str, kernel: RfKernel) {
    let rows = section(name);
    assert!(rows.len() >= 200, "{} {name} rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 21, "{name} row {row:?}");
        let (f, e, de) = kernel(
            row[0],
            row[1],
            row[2],
            row[3],
            &table_of(&row[9..18]),
            &crf_of(&row[4..8]),
            row[8],
        );
        check(row, &[f, e, de], &row[18..], &mut mismatches);
    }
    report(name, rows.len(), &mismatches);
}

#[test]
fn rf_soft_interaction_equals_lean_model_bit_for_bit() {
    rf_rows_equal_lean("rf_soft_interaction", rf_soft_interaction);
}

#[test]
fn perturbed_rf_excluded_pair_equals_lean_model_bit_for_bit() {
    rf_rows_equal_lean("perturbed_rf_excluded_pair", perturbed_rf_excluded_pair);
}

#[test]
fn perturbed_rf_self_term_equals_lean_model_bit_for_bit() {
    let rows = section("perturbed_rf_self_term");
    assert!(rows.len() >= 100, "{} self rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 18, "self row {row:?}");
        let (e, de) = perturbed_rf_self_term(
            row[0],
            row[1],
            &table_of(&row[7..16]),
            &crf_of(&row[2..6]),
            row[6],
        );
        check(row, &[e, de], &row[16..], &mut mismatches);
    }
    report("perturbed_rf_self_term", rows.len(), &mismatches);
}

#[test]
fn lj_crf_soft_interaction_ext_equals_lean_model_bit_for_bit() {
    let rows = section("lj_crf_soft_interaction_ext");
    assert!(rows.len() >= 200, "{} ext rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 23, "ext row {row:?}");
        let (e, de) = lj_crf_soft_interaction_ext(
            row[0],
            row[1],
            row[2],
            row[3],
            row[4],
            row[5],
            row[6],
            row[7],
            row[8],
            &crf_of(&row[9..13]),
            row[13],
            row[14],
        );
        let rust: Vec<f64> = e.iter().chain(&de).copied().collect();
        check(row, &rust, &row[15..], &mut mismatches);
    }
    report("lj_crf_soft_interaction_ext", rows.len(), &mismatches);
}

#[test]
fn rf_soft_interaction_ext_equals_lean_model_bit_for_bit() {
    let rows = section("rf_soft_interaction_ext");
    assert!(rows.len() >= 100, "{} rf ext rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 14, "rf ext row {row:?}");
        let (e, de) = rf_soft_interaction_ext(
            row[0],
            row[1],
            row[2],
            row[3],
            &crf_of(&row[4..8]),
            row[8],
            row[9],
        );
        check(
            row,
            &[e[0], e[1], de[0], de[1]],
            &row[10..],
            &mut mismatches,
        );
    }
    report("rf_soft_interaction_ext", rows.len(), &mismatches);
}

/// The generated rows cover both NLAM, the end states, every α, dummy states, the guards.
#[test]
fn rows_cover_the_edges() {
    let pair = section("perturbed_lj_crf_pair");
    let guard = f64::from_bits(GUARD_BITS);
    for n in [1.0, 2.0] {
        assert!(pair.iter().any(|r| r[23] == n));
    }
    for lam in [0.0, 0.25, 0.5, 0.75, 1.0] {
        assert!(pair.iter().any(|r| r[15..19].iter().all(|&l| l == lam)));
    }
    for alpha in [0.0, 0.5, 1.0, 1.51] {
        assert!(pair.iter().any(|r| r[7] == alpha));
        assert!(pair.iter().any(|r| r[8] == alpha));
    }
    // individual λ: different λ per kind, derivatives other than 1
    assert!(pair.iter().any(|r| r[15] != r[16] && r[19] != 1.0));
    // a dummy state A, an uncharged state B
    assert!(pair
        .iter()
        .any(|r| r[1] == 0.0 && r[2] == 0.0 && r[5] == 0.0));
    assert!(pair.iter().any(|r| r[4] > 0.0 && r[6] == 0.0));
    // below the distance guard, at the C6 guard
    assert!(pair.iter().any(|r| r[0] < guard));
    assert!(pair.iter().any(|r| r[1] == guard || r[3] == guard));
    assert!(pair
        .iter()
        .any(|r| (r[1] > 0.0 && r[1] < guard) || (r[3] > 0.0 && r[3] < guard)));
}

/// `h` of the central differences.
const H: f64 = 1e-6;

/// A row's table at `λ + δ`, each individual λ moved along its own slope: `Λ_k + δ·Λ_k′`.
fn table_shifted(t: &[f64], delta: f64) -> PerturbedLambdaParams {
    let mut s = t.to_vec();
    for k in 0..4 {
        s[k] = t[k] + delta * t[k + 4];
    }
    table_of(&s)
}

/// `|a − b| ≤ 1e-6·|a| + 1e-9·scale`: relative to the derivative, plus the rounding floor of a
/// central difference of energies of size `scale` with `h = 1e-6`.
fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-6 * a.abs() + 1e-9 * scale
}

/// Rows away from the soft-core singularity and the guards: `r² ≥ 0.01` nm², a reaction-field
/// cutoff, and each `C6` zero or at least `1e-6` (a smaller one makes `C12/C6` so large that the
/// energy's second λ derivative swamps a central difference).
fn well_conditioned(row: &[f64]) -> bool {
    let c6_ok = |c6: f64| c6 == 0.0 || c6 >= 1e-6;
    row[0] >= 0.01 && row[9] > 0.0 && c6_ok(row[1]) && c6_ok(row[3])
}

/// dV/dλ of the pair kernel equals the central difference of its energy in λ, with each
/// individual λ moving along its own slope (the chain rule the kernel applies).
#[test]
fn perturbed_lj_crf_pair_dvdl_is_the_derivative_of_its_energy() {
    let rows = section("perturbed_lj_crf_pair");
    let mut checked = 0;
    let mut failures = Vec::new();
    for row in rows.iter().filter(|r| well_conditioned(r)) {
        let t = &row[15..24];
        let (_, e_lj, e_crf, de_lj, de_crf) = pair_with(row, &table_of(t), row[0]);
        let (_, lj_p, crf_p, _, _) = pair_with(row, &table_shifted(t, H), row[0]);
        let (_, lj_m, crf_m, _, _) = pair_with(row, &table_shifted(t, -H), row[0]);
        let fd_lj = (lj_p - lj_m) / (2.0 * H);
        let fd_crf = (crf_p - crf_m) / (2.0 * H);
        checked += 1;
        if !close(de_lj, fd_lj, 1.0 + e_lj.abs()) || !close(de_crf, fd_crf, 1.0 + e_crf.abs()) {
            failures.push(format!(
                "  row={row:?}: de_lj={de_lj:e} fd={fd_lj:e}; de_crf={de_crf:e} fd={fd_crf:e}"
            ));
        }
    }
    assert!(checked >= 250, "{checked} well-conditioned rows");
    assert!(
        failures.is_empty(),
        "{} of {checked} rows:\n{}",
        failures.len(),
        failures[..failures.len().min(10)].join("\n")
    );
}

/// The scalar force of the pair kernel is `−dV/dr / r`, against a central difference in `r`.
#[test]
fn perturbed_lj_crf_pair_force_is_minus_the_gradient() {
    let rows = section("perturbed_lj_crf_pair");
    let mut checked = 0;
    let mut failures = Vec::new();
    for row in rows.iter().filter(|r| well_conditioned(r)) {
        let t = table_of(&row[15..24]);
        let r = row[0].sqrt();
        let (f, e_lj, e_crf, _, _) = pair_with(row, &t, row[0]);
        let v = |x: f64| {
            let (_, lj, crf, _, _) = pair_with(row, &t, x * x);
            lj + crf
        };
        let fd = (v(r + H) - v(r - H)) / (2.0 * H);
        checked += 1;
        if !close(-f * r, fd, 1.0 + e_lj.abs() + e_crf.abs()) {
            failures.push(format!("  row={row:?}: -f·r={:e} fd={fd:e}", -f * r));
        }
    }
    assert!(checked >= 250, "{checked} well-conditioned rows");
    assert!(
        failures.is_empty(),
        "{} of {checked} rows:\n{}",
        failures.len(),
        failures[..failures.len().min(10)].join("\n")
    );
}

/// dV/dλ of both excluded-pair kernels and of the self term against central differences.
#[test]
fn rf_kernels_dvdl_is_the_derivative_of_their_energy() {
    for (name, kernel) in [
        ("rf_soft_interaction", rf_soft_interaction as RfKernel),
        ("perturbed_rf_excluded_pair", perturbed_rf_excluded_pair),
    ] {
        let mut checked = 0;
        for row in section(name).iter().filter(|r| r[4] > 0.0) {
            let t = &row[9..18];
            let crf = crf_of(&row[4..8]);
            let at = |tab: &PerturbedLambdaParams| {
                kernel(row[0], row[1], row[2], row[3], tab, &crf, row[8])
            };
            let (_, e, de) = at(&table_of(t));
            let fd = (at(&table_shifted(t, H)).1 - at(&table_shifted(t, -H)).1) / (2.0 * H);
            assert!(
                close(de, fd, 1.0 + e.abs()),
                "{name} row={row:?}: de={de:e} fd={fd:e}"
            );
            checked += 1;
        }
        assert!(checked >= 100, "{name}: {checked} rows");
    }
    for row in section("perturbed_rf_self_term") {
        let t = &row[7..16];
        let crf = crf_of(&row[2..6]);
        let at =
            |tab: &PerturbedLambdaParams| perturbed_rf_self_term(row[0], row[1], tab, &crf, row[6]);
        let (e, de) = at(&table_of(t));
        let fd = (at(&table_shifted(t, H)).0 - at(&table_shifted(t, -H)).0) / (2.0 * H);
        assert!(
            close(de, fd, 1.0 + e.abs()),
            "self row={row:?}: de={de:e} fd={fd:e}"
        );
    }
}

/// The ext kernels: `A_de = ∂E_A/∂λ`, `B_de = −∂E_B/∂λ`, against central differences.
#[test]
fn ext_kernels_derivatives_are_the_derivatives_of_their_energies() {
    let mut checked = 0;
    for row in section("lj_crf_soft_interaction_ext")
        .iter()
        .filter(|r| well_conditioned(r))
    {
        let at = |lam: f64| {
            lj_crf_soft_interaction_ext(
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                row[5],
                row[6],
                row[7],
                row[8],
                &crf_of(&row[9..13]),
                row[13],
                lam,
            )
        };
        let (e, de) = at(row[14]);
        let (ep, _) = at(row[14] + H);
        let (em, _) = at(row[14] - H);
        for k in 0..4 {
            let fd = (ep[k] - em[k]) / (2.0 * H);
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            assert!(
                close(sign * de[k], fd, 1.0 + e[k].abs()),
                "ext column {k} row={row:?}: de={:e} fd={fd:e}",
                sign * de[k]
            );
        }
        checked += 1;
    }
    assert!(checked >= 150, "{checked} ext rows");
    for row in section("rf_soft_interaction_ext")
        .iter()
        .filter(|r| r[4] > 0.0)
    {
        let at = |lam: f64| {
            rf_soft_interaction_ext(
                row[0],
                row[1],
                row[2],
                row[3],
                &crf_of(&row[4..8]),
                row[8],
                lam,
            )
        };
        let (e, de) = at(row[9]);
        let (ep, _) = at(row[9] + H);
        let (em, _) = at(row[9] - H);
        for k in 0..2 {
            let fd = (ep[k] - em[k]) / (2.0 * H);
            let sign = if k == 0 { 1.0 } else { -1.0 };
            assert!(
                close(sign * de[k], fd, 1.0 + e[k].abs()),
                "rf ext column {k} row={row:?}: de={:e} fd={fd:e}",
                sign * de[k]
            );
        }
    }
}
