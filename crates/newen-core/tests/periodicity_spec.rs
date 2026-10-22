//! The periodic-boundary kernels against their Lean model, bit for bit.
//!
//! `spec/vectors/periodicity.json` holds the `*Float` definitions of `NewenSpec.Periodicity`
//! evaluated at Lean's `Float` (IEEE 754 binary64) on random inputs, written by
//! `lake exe spec-vectors` (checked in, so this test needs no Lean). Every number is the hex string
//! of its 64 bits. The kernels use `+ − × ÷`, comparisons, `floor` and round-to-nearest-even
//! (`f64::round_ties_even`, C `rint`): all exact or correctly rounded in Rust and in Lean (C)
//! alike, so the Rust must return the same bits as Lean. A NaN matches any NaN (Lean has a single
//! NaN).
//!
//! Each row carries the box, not its precomputed reciprocals or inverse: the Rust builds
//! `Rectangular::new` / `Triclinic::new` from it, so every row also exercises the constructor, and
//! calls the kernels through `Periodicity`, the way the force loops do.

use std::collections::HashMap;

use newen_core::math::{
    rectangular_nearest_image, rectangular_put_into_box, triclinic_nearest_image,
    triclinic_put_into_box, truncoct_triclinic_box, vacuum_nearest_image, BoundaryCondition, Mat3,
    Periodicity, Rectangular, Triclinic, Vacuum, Vec3,
};

/// `pointFiveFloat` in `spec/NewenSpec/Periodicity.lean`.
const POINT_FIVE_BITS: u64 = 0x3fe0_0000_0000_0000;

/// The sections of `periodicity.json`, each a list of rows of floats.
fn lean_rows() -> HashMap<String, Vec<Vec<f64>>> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/periodicity.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let mut sections: HashMap<String, Vec<Vec<f64>>> = HashMap::new();
    let mut section = String::new();
    for line in text.lines().map(str::trim) {
        if let Some(name) = line.strip_prefix('"').and_then(|l| l.strip_suffix("\": [")) {
            section = name.to_string();
            continue;
        }
        let line = line.trim_end_matches(',');
        if !(line.starts_with('[') && line.ends_with(']')) {
            continue;
        }
        let fields: Vec<f64> = line[1..line.len() - 1]
            .split(',')
            .map(|f| float(f.trim().trim_matches('"')))
            .collect();
        assert!(!section.is_empty(), "row outside a section: {line}");
        sections.entry(section.clone()).or_default().push(fields);
    }
    sections
}

fn section(rows: &HashMap<String, Vec<Vec<f64>>>, name: &str, len: usize) -> Vec<Vec<f64>> {
    let rows = rows
        .get(name)
        .unwrap_or_else(|| panic!("no section {name} in periodicity.json"))
        .clone();
    for row in &rows {
        assert_eq!(row.len(), len, "{name} row {row:?}");
    }
    rows
}

fn float(field: &str) -> f64 {
    let hex = field
        .strip_prefix("0x")
        .expect("hex bits in periodicity.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in periodicity.json"))
}

fn vec3(f: &[f64]) -> Vec3 {
    Vec3::new(f[0], f[1], f[2])
}

/// A matrix from its columns `a, b, c`.
fn mat3(f: &[f64]) -> Mat3 {
    Mat3::from_cols(vec3(&f[0..3]), vec3(&f[3..6]), vec3(&f[6..9]))
}

fn parts(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn mat_parts(m: Mat3) -> Vec<f64> {
    [m.x_axis, m.y_axis, m.z_axis]
        .into_iter()
        .flat_map(parts)
        .collect()
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
    assert_eq!(0.5_f64.to_bits(), POINT_FIVE_BITS);
}

#[test]
fn vacuum_equals_lean_model_bit_for_bit() {
    let rows = section(&lean_rows(), "vacuum_nearest_image", 9);
    assert!(rows.len() >= 200, "{} vacuum rows", rows.len());
    let pbc = Periodicity::Vacuum(Vacuum);
    let mut mismatches = Vec::new();
    for row in &rows {
        let (ri, rj) = (vec3(&row[0..3]), vec3(&row[3..6]));
        let lean = &row[6..9];
        check(
            row,
            &parts(pbc.nearest_image(ri, rj)),
            lean,
            &mut mismatches,
        );
        check(
            row,
            &parts(vacuum_nearest_image(ri, rj)),
            lean,
            &mut mismatches,
        );
        // put_into_box is the identity
        check(
            row,
            &parts(pbc.put_into_box(ri)),
            &row[0..3],
            &mut mismatches,
        );
    }
    report("vacuum", rows.len(), &mismatches);
}

#[test]
fn rectangular_new_equals_lean_model_bit_for_bit() {
    let rows = section(&lean_rows(), "rectangular_new", 9);
    assert!(rows.len() >= 100, "{} rectangular_new rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        let rect = Rectangular::new(vec3(&row[0..3]));
        let rust = [parts(rect.half_box), parts(rect.inv_box)].concat();
        check(row, &rust, &row[3..9], &mut mismatches);
    }
    report("Rectangular::new", rows.len(), &mismatches);
}

#[test]
fn rectangular_nearest_image_equals_lean_model_bit_for_bit() {
    let rows = section(&lean_rows(), "rectangular_nearest_image", 12);
    assert!(rows.len() >= 1_000, "{} rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        let rect = Rectangular::new(vec3(&row[0..3]));
        let (ri, rj) = (vec3(&row[3..6]), vec3(&row[6..9]));
        let lean = &row[9..12];
        let pbc = Periodicity::Rectangular(rect);
        check(
            row,
            &parts(pbc.nearest_image(ri, rj)),
            lean,
            &mut mismatches,
        );
        let pure = rectangular_nearest_image(rect.box_size, rect.inv_box, ri, rj);
        check(row, &parts(pure), lean, &mut mismatches);
    }
    report("rectangular nearest_image", rows.len(), &mismatches);
}

#[test]
fn rectangular_put_into_box_equals_lean_model_bit_for_bit() {
    let rows = section(&lean_rows(), "rectangular_put_into_box", 9);
    assert!(rows.len() >= 500, "{} rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        let rect = Rectangular::new(vec3(&row[0..3]));
        let pos = vec3(&row[3..6]);
        let lean = &row[6..9];
        let pbc = Periodicity::Rectangular(rect);
        check(row, &parts(pbc.put_into_box(pos)), lean, &mut mismatches);
        let pure = rectangular_put_into_box(rect.box_size, rect.inv_box, pos);
        check(row, &parts(pure), lean, &mut mismatches);
    }
    report("rectangular put_into_box", rows.len(), &mismatches);
}

#[test]
fn truncoct_boxes_are_the_ones_newen_builds() {
    let rows = section(&lean_rows(), "truncoct_box", 10);
    let mut mismatches = Vec::new();
    for row in &rows {
        let m = truncoct_triclinic_box(Mat3::from_diagonal(Vec3::splat(row[0])), true);
        check(row, &mat_parts(m), &row[1..10], &mut mismatches);
    }
    report("truncoct_triclinic_box", rows.len(), &mismatches);
}

#[test]
fn triclinic_new_equals_lean_model_bit_for_bit() {
    let rows = section(&lean_rows(), "triclinic_new", 18);
    assert!(rows.len() >= 150, "{} triclinic_new rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        let tri = Triclinic::new(mat3(&row[0..9]));
        check(
            row,
            &mat_parts(tri.inv_box_matrix),
            &row[9..18],
            &mut mismatches,
        );
    }
    report("Triclinic::new", rows.len(), &mismatches);
}

#[test]
fn triclinic_nearest_image_equals_lean_model_bit_for_bit() {
    let rows = section(&lean_rows(), "triclinic_nearest_image", 18);
    assert!(rows.len() >= 1_000, "{} rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        let tri = Triclinic::new(mat3(&row[0..9]));
        let (ri, rj) = (vec3(&row[9..12]), vec3(&row[12..15]));
        let lean = &row[15..18];
        let pbc = Periodicity::Triclinic(tri);
        check(
            row,
            &parts(pbc.nearest_image(ri, rj)),
            lean,
            &mut mismatches,
        );
        let pure = triclinic_nearest_image(tri.box_matrix, ri, rj);
        check(row, &parts(pure), lean, &mut mismatches);
    }
    report("triclinic nearest_image", rows.len(), &mismatches);
}

#[test]
fn triclinic_put_into_box_equals_lean_model_bit_for_bit() {
    let rows = section(&lean_rows(), "triclinic_put_into_box", 15);
    assert!(rows.len() >= 500, "{} rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        let tri = Triclinic::new(mat3(&row[0..9]));
        let pos = vec3(&row[9..12]);
        let lean = &row[12..15];
        let pbc = Periodicity::Triclinic(tri);
        check(row, &parts(pbc.put_into_box(pos)), lean, &mut mismatches);
        // the trait method, as the generic callers use it
        check(row, &parts(tri.put_into_box(pos)), lean, &mut mismatches);
        let pure = triclinic_put_into_box(tri.box_matrix, tri.inv_box_matrix, pos);
        check(row, &parts(pure), lean, &mut mismatches);
    }
    report("triclinic put_into_box", rows.len(), &mismatches);
}

/// The generated rows cover exact half-box ties, far separations, negative coordinates, NaN and
/// infinity, faces of the triclinic cell, and every box shape.
#[test]
fn rows_cover_the_edges() {
    let rows = lean_rows();
    let rect = section(&rows, "rectangular_nearest_image", 12);
    let sep = |r: &Vec<f64>, d: usize| r[3 + d] - r[6 + d];
    // a separation of exactly half an edge, of either sign, and an odd multiple of it
    assert!(rect.iter().any(|r| sep(r, 0) == 0.5 * r[0]));
    assert!(rect.iter().any(|r| sep(r, 0) == -0.5 * r[0]));
    assert!(rect.iter().any(|r| sep(r, 1).abs() == 1.5 * r[1]));
    // far outside the box, negative, NaN, infinite
    assert!(rect.iter().any(|r| sep(r, 2).abs() > 20.0 * r[2]));
    assert!(rect.iter().any(|r| r[3] < -r[0]));
    assert!(rect.iter().any(|r| r[3..9].iter().any(|x| x.is_nan())));
    assert!(rect.iter().any(|r| r[3..9].iter().any(|x| x.is_infinite())));
    // anisotropic boxes
    assert!(rect.iter().any(|r| r[0] != r[1]));

    let put = section(&rows, "rectangular_put_into_box", 9);
    assert!(put.iter().any(|r| r[3] == -2.0 * r[0]));
    assert!(put.iter().any(|r| r[3] == r[0]));

    let tri = section(&rows, "triclinic_nearest_image", 18);
    let truncoct = section(&rows, "truncoct_box", 10);
    // truncated octahedra, rectangular, skewed (|b_x| > a_x/2) and monoclinic boxes
    assert!(tri
        .iter()
        .any(|r| truncoct.iter().any(|t| t[1..10] == r[0..9])));
    assert!(tri
        .iter()
        .any(|r| r[3] == 0.0 && r[6] == 0.0 && r[7] == 0.0));
    assert!(tri.iter().any(|r| r[3].abs() > 0.5 * r[0]));
    assert!(tri
        .iter()
        .any(|r| r[3] == 0.0 && r[6] != 0.0 && r[7] == 0.0));
    // separations exactly on the faces z = ±c_z/2 of the cell, and far ones
    let tsep = |r: &Vec<f64>, d: usize| r[9 + d] - r[12 + d];
    assert!(tri.iter().any(|r| tsep(r, 2) == 0.5 * r[8]));
    assert!(tri.iter().any(|r| tsep(r, 2) == -0.5 * r[8]));
    assert!(tri.iter().any(|r| tsep(r, 2).abs() > 20.0 * r[8]));
    assert!(tri.iter().any(|r| r[9..15].iter().any(|x| x.is_nan())));
}
