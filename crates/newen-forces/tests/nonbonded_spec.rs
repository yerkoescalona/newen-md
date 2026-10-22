//! The Lennard-Jones and reaction-field pair kernels against their Lean model, bit for bit.
//!
//! `spec/vectors/nonbonded.json` holds `NewenSpec.Nonbonded.crfParametersFloat`,
//! `ljCrfPairFloat`, `rfExcludedPairFloat` and `rfSelfTermFloat` evaluated at Lean's `Float`
//! (IEEE 754 binary64) on random inputs, written by `lake exe spec-vectors` (checked in, so this
//! test needs no Lean). Every number is the hex string of its 64 bits. The kernels use only
//! `+ − × ÷` and `sqrt`, which IEEE 754 rounds correctly in Rust and in Lean (C) alike, so the
//! Rust must return the same bits as Lean. A NaN matches any NaN (Lean has a single NaN).
//!
//! Each row carries the four inputs of `CRFParameters::new`; the pair rows rebuild the
//! parameters from them on the Rust side, so a pair row also exercises `CRFParameters::new`.

use newen_forces::nonbonded::{
    lj_crf_pair, rf_excluded_pair, rf_self_term, CRFParameters, ZERO_DISTANCE2,
};

/// `zeroDistance2Float` in `spec/NewenSpec/Nonbonded.lean`.
const ZERO_DISTANCE2_BITS: u64 = 0x3ddb_7cdf_d9d7_bdbb;

/// The four sections of `nonbonded.json`, each a list of rows of floats.
struct Rows {
    crf_parameters: Vec<Vec<f64>>,
    lj_crf_pair: Vec<Vec<f64>>,
    rf_excluded_pair: Vec<Vec<f64>>,
    rf_self_term: Vec<Vec<f64>>,
}

fn lean_rows() -> Rows {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/nonbonded.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let mut rows = Rows {
        crf_parameters: Vec::new(),
        lj_crf_pair: Vec::new(),
        rf_excluded_pair: Vec::new(),
        rf_self_term: Vec::new(),
    };
    let mut section = "";
    for line in text.lines().map(str::trim) {
        for name in [
            "crf_parameters",
            "lj_crf_pair",
            "rf_excluded_pair",
            "rf_self_term",
        ] {
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
            "crf_parameters" => rows.crf_parameters.push(fields),
            "lj_crf_pair" => rows.lj_crf_pair.push(fields),
            "rf_excluded_pair" => rows.rf_excluded_pair.push(fields),
            "rf_self_term" => rows.rf_self_term.push(fields),
            _ => panic!("row outside a section: {line}"),
        }
    }
    rows
}

fn float(field: &str) -> f64 {
    let hex = field
        .strip_prefix("0x")
        .expect("hex bits in nonbonded.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in nonbonded.json"))
}

/// `CRFParameters::new` of a row's `cutoff, epsilon, rf_epsilon, rf_kappa`.
fn crf_of(inputs: &[f64]) -> CRFParameters {
    CRFParameters::new(inputs[0], inputs[1], inputs[2], inputs[3])
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

/// Compare a row's Rust outputs with its Lean outputs; on a difference, describe the row.
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

#[test]
fn constants_are_the_bits_the_lean_model_uses() {
    assert_eq!(ZERO_DISTANCE2.to_bits(), ZERO_DISTANCE2_BITS);
}

#[test]
fn crf_parameters_equal_lean_model_bit_for_bit() {
    let rows = lean_rows().crf_parameters;
    assert!(rows.len() >= 500, "{} crf_parameters rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 8, "crf_parameters row {row:?}");
        let p = crf_of(&row[..4]);
        let rust = [p.crf_cut, p.crf_2cut3i, p.crf_cut3i, p.cutoff_sq];
        check(row, &rust, &row[4..], &mut mismatches);
    }
    report("CRFParameters::new", rows.len(), &mismatches);
}

#[test]
fn lj_crf_pair_equals_lean_model_bit_for_bit() {
    let rows = lean_rows().lj_crf_pair;
    assert!(rows.len() >= 5_000, "{} lj_crf_pair rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 11, "lj_crf_pair row {row:?}");
        let (r2, c6, c12, q_prod) = (row[0], row[1], row[2], row[3]);
        let (f_mag, e_lj, e_crf) = lj_crf_pair(r2, c6, c12, q_prod, &crf_of(&row[4..8]));
        check(row, &[f_mag, e_lj, e_crf], &row[8..], &mut mismatches);
    }
    report("lj_crf_pair", rows.len(), &mismatches);
}

#[test]
fn rf_excluded_pair_equals_lean_model_bit_for_bit() {
    let rows = lean_rows().rf_excluded_pair;
    assert!(rows.len() >= 1_000, "{} rf_excluded_pair rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 8, "rf_excluded_pair row {row:?}");
        let (f_mag, e) = rf_excluded_pair(row[0], row[1], &crf_of(&row[2..6]));
        check(row, &[f_mag, e], &row[6..], &mut mismatches);
    }
    report("rf_excluded_pair", rows.len(), &mismatches);
}

#[test]
fn rf_self_term_equals_lean_model_bit_for_bit() {
    let rows = lean_rows().rf_self_term;
    assert!(rows.len() >= 200, "{} rf_self_term rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 7, "rf_self_term row {row:?}");
        let e = rf_self_term(row[0], row[1], &crf_of(&row[2..6]));
        check(row, &[e], &row[6..], &mut mismatches);
    }
    report("rf_self_term", rows.len(), &mismatches);
}

/// The generated rows cover the guard, the cutoff and both reaction-field sentinels.
#[test]
fn rows_cover_the_edges() {
    let rows = lean_rows();
    let lj = &rows.lj_crf_pair;
    assert!(lj.iter().any(|r| r[0] < ZERO_DISTANCE2));
    assert!(lj.iter().any(|r| r[0] == ZERO_DISTANCE2));
    assert!(lj.iter().any(|r| r[4] > 0.0 && r[0] == r[4] * r[4]));
    let crf = &rows.crf_parameters;
    assert!(crf.iter().any(|r| r[0] == 0.0));
    assert!(crf.iter().any(|r| r[0] > 0.0 && r[2] == 0.0));
    assert!(crf
        .iter()
        .any(|r| r[0] > 0.0 && r[2] == r[1] && r[3] == 0.0));
}
