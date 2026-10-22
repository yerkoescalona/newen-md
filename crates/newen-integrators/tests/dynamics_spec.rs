//! The dynamics kernels against their Lean model, bit for bit.
//!
//! `spec/vectors/dynamics.json` holds the models of `NewenSpec.Dynamics` evaluated at Lean's
//! `Float` (IEEE 754 binary64) on random inputs, written by `lake exe spec-vectors` (checked in,
//! so this test needs no Lean). Every number is the hex string of its 64 bits. The leap-frog,
//! kinetic-energy, temperature and pressure kernels use only `+ − × ÷`, which IEEE 754 rounds
//! correctly in Rust and in Lean (C) alike, so the Rust must return the same bits as Lean. The
//! barostat calls `powf`, which is not required to be correctly rounded: Rust's `f64::powf` and
//! Lean's `Float.pow` both call the C library's `pow`, so its rows agree when both run against
//! the same `libm` (the platform's).
//!
//! One exception: Lean's `Float` has a single NaN, so a NaN matches any NaN; every other value
//! must match bit for bit.

use newen_core::math::{Mat3, Vec3};
use newen_integrators::algorithms::{
    berendsen_barostat_scale, group_kinetic_energy, kinetic_energy, kinetic_temperature,
    leap_frog_position, leap_frog_velocity, pressure_tensor, scalar_pressure,
};

/// `halfFloat` in `spec/NewenSpec/Dynamics.lean`: the `0.5` of the kinetic energies and the
/// pressure tensor.
const HALF_BITS: u64 = 0x3fe0_0000_0000_0000;
/// `thirdFloat` in `spec/NewenSpec/Dynamics.lean`: the `1.0 / 3.0` of the barostat.
const THIRD_BITS: u64 = 0x3fd5_5555_5555_5555;

/// The rows of one section of `dynamics.json`, each a list of fields (hex bits or a count).
fn lean_rows(section: &str) -> Vec<Vec<String>> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/dynamics.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let header = format!("\"{section}\": [");
    let mut rows = Vec::new();
    let mut inside = false;
    for line in text.lines().map(str::trim) {
        if line.starts_with('"') {
            inside = line == header;
            continue;
        }
        let line = line.trim_end_matches(',');
        if !(inside && line.starts_with('[') && line.ends_with(']')) {
            continue;
        }
        rows.push(
            line[1..line.len() - 1]
                .split(',')
                .map(|f| f.trim().trim_matches('"').to_string())
                .collect(),
        );
    }
    rows
}

fn float(field: &str) -> f64 {
    let hex = field.strip_prefix("0x").expect("hex bits in dynamics.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in dynamics.json"))
}

fn floats(fields: &[String]) -> Vec<f64> {
    fields.iter().map(|f| float(f)).collect()
}

fn vec3(v: &[f64]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

/// A `Mat3` from nine numbers, column by column.
fn mat3(v: &[f64]) -> Mat3 {
    Mat3::from_cols(vec3(&v[0..3]), vec3(&v[3..6]), vec3(&v[6..9]))
}

fn mat3_values(m: Mat3) -> [f64; 9] {
    m.to_cols_array()
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

/// Compares the outputs of one row; on a mismatch, records the inputs, both outputs and the
/// largest ulp distance.
fn check(mismatches: &mut Vec<String>, inputs: &str, rust: &[f64], lean: &[f64]) {
    assert_eq!(rust.len(), lean.len());
    if rust.iter().zip(lean).all(|(&r, &l)| same(r, l)) {
        return;
    }
    let worst = rust
        .iter()
        .zip(lean)
        .filter(|(&r, &l)| !same(r, l))
        .map(|(&r, &l)| ulps(r, l))
        .max()
        .unwrap_or(0);
    mismatches.push(format!(
        "  {inputs}: rust={rust:?} lean={lean:?} ({worst} ulp)"
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
    assert_eq!(0.5_f64.to_bits(), HALF_BITS);
    assert_eq!((1.0_f64 / 3.0).to_bits(), THIRD_BITS);
}

#[test]
fn leap_frog_equals_lean_model_bit_for_bit() {
    let rows = lean_rows("leap_frog");
    assert!(rows.len() >= 1_000, "{} leap-frog rows", rows.len());
    let (mut vel_bad, mut pos_bad) = (Vec::new(), Vec::new());
    for row in &rows {
        assert_eq!(row.len(), 17, "leap-frog row {row:?}");
        let v = floats(row);
        let (vel, force, im, dt, pos) =
            (vec3(&v[0..3]), vec3(&v[3..6]), v[6], v[7], vec3(&v[8..11]));
        let (vel_lean, pos_lean) = (vec3(&v[11..14]), vec3(&v[14..17]));
        let vel_rust = leap_frog_velocity(vel, force, im, dt);
        check(
            &mut vel_bad,
            &format!("v={vel:?} f={force:?} im={im:e} dt={dt}"),
            &vel_rust.to_array(),
            &vel_lean.to_array(),
        );
        // the position step from the model's new velocity, so each kernel is tested on its own
        let pos_rust = leap_frog_position(pos, vel_lean, dt);
        check(
            &mut pos_bad,
            &format!("x={pos:?} v={vel_lean:?} dt={dt}"),
            &pos_rust.to_array(),
            &pos_lean.to_array(),
        );
    }
    report("leap_frog_velocity", rows.len(), &vel_bad);
    report("leap_frog_position", rows.len(), &pos_bad);
}

#[test]
fn kinetic_energies_equal_lean_model_bit_for_bit() {
    let rows = lean_rows("kinetic");
    assert!(rows.len() >= 400, "{} kinetic rows", rows.len());
    let (mut total_bad, mut group_bad) = (Vec::new(), Vec::new());
    for row in &rows {
        let n: usize = row[0].parse().expect("atom count in dynamics.json");
        assert_eq!(row.len(), 1 + 7 * n + 6, "kinetic row {row:?}");
        let v = floats(&row[1..]);
        let mass: Vec<f64> = (0..n).map(|a| v[7 * a]).collect();
        let vel_new: Vec<Vec3> = (0..n).map(|a| vec3(&v[7 * a + 1..7 * a + 4])).collect();
        let vel_old: Vec<Vec3> = (0..n).map(|a| vec3(&v[7 * a + 4..7 * a + 7])).collect();
        let lean = &v[7 * n..];
        let inputs = format!("m={mass:?} v_new={vel_new:?} v_old={vel_old:?}");

        let (e, e_new) = kinetic_energy(&mass, &vel_new, &vel_old);
        check(&mut total_bad, &inputs, &[e, e_new], &lean[0..2]);
        let g = group_kinetic_energy(&mass, &vel_new, &vel_old);
        check(
            &mut group_bad,
            &inputs,
            &[g.com, g.total, g.new_com, g.new_total],
            &lean[2..6],
        );
    }
    report("kinetic_energy", rows.len(), &total_bad);
    report("group_kinetic_energy", rows.len(), &group_bad);
}

#[test]
fn temperature_equals_lean_model_bit_for_bit() {
    let rows = lean_rows("temperature");
    assert!(rows.len() >= 500, "{} temperature rows", rows.len());
    let mut bad = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 4, "temperature row {row:?}");
        let v = floats(row);
        check(
            &mut bad,
            &format!("ekin={:e} dof={} kB={:e}", v[0], v[1], v[2]),
            &[kinetic_temperature(v[0], v[1], v[2])],
            &v[3..4],
        );
    }
    report("kinetic_temperature", rows.len(), &bad);
    // a negative and a zero dof, both at temperature 0
    let v: Vec<Vec<f64>> = rows.iter().map(|r| floats(r)).collect();
    assert!(v.iter().any(|r| r[1] < 0.0 && r[0] > 0.0 && r[3] == 0.0));
    assert!(v.iter().any(|r| r[1] == 0.0 && r[0] > 0.0 && r[3] == 0.0));
}

#[test]
fn pressure_equals_lean_model_bit_for_bit() {
    let rows = lean_rows("pressure");
    assert!(rows.len() >= 500, "{} pressure rows", rows.len());
    let (mut tensor_bad, mut scalar_bad) = (Vec::new(), Vec::new());
    for row in &rows {
        assert_eq!(row.len(), 29, "pressure row {row:?}");
        let v = floats(row);
        let (kinetic, virial, volume) = (mat3(&v[0..9]), mat3(&v[9..18]), v[18]);
        let (p_lean, scalar_lean) = (mat3(&v[19..28]), v[28]);
        let p = pressure_tensor(kinetic, virial, volume);
        check(
            &mut tensor_bad,
            &format!("K={kinetic:?} W={virial:?} V={volume:e}"),
            &mat3_values(p),
            &mat3_values(p_lean),
        );
        check(
            &mut scalar_bad,
            &format!("P={p_lean:?}"),
            &[scalar_pressure(p_lean)],
            &[scalar_lean],
        );
    }
    report("pressure_tensor", rows.len(), &tensor_bad);
    report("scalar_pressure", rows.len(), &scalar_bad);
}

#[test]
fn berendsen_barostat_equals_lean_model_bit_for_bit() {
    let rows = lean_rows("barostat");
    assert!(rows.len() >= 1_000, "{} barostat rows", rows.len());
    let mut bad = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 6, "barostat row {row:?}");
        let v = floats(row);
        let (p, p0, beta, tau, dt) = (v[0], v[1], v[2], v[3], v[4]);
        check(
            &mut bad,
            &format!("P={p:e} P0={p0:e} beta={beta:e} tau={tau} dt={dt}"),
            &[berendsen_barostat_scale(p, p0, beta, tau, dt)],
            &v[5..6],
        );
    }
    report("berendsen_barostat_scale", rows.len(), &bad);
}

/// At the reference pressure the box is left alone, exactly (`barostat_fixed` over ℝ).
#[test]
fn barostat_at_the_reference_pressure_is_exactly_one() {
    for p in [0.0, 0.06102, 1.0, -3.5] {
        assert_eq!(berendsen_barostat_scale(p, p, 4.575e-4, 0.5, 0.002), 1.0);
    }
}

/// Constant force: leap-frog is exact (`leapFrog_iterate_exact` over ℝ), here to rounding.
#[test]
fn leap_frog_is_exact_for_a_constant_force() {
    let (im, dt, a) = (0.5, 0.01, Vec3::new(3.0, -2.0, 1.0));
    let (x0, u) = (Vec3::new(1.0, 2.0, 3.0), Vec3::new(-0.5, 0.25, 1.0));
    let (mut x, mut v) = (x0, u - a * (dt / 2.0));
    let n = 100;
    for _ in 0..n {
        v = leap_frog_velocity(v, a / im, im, dt);
        x = leap_frog_position(x, v, dt);
    }
    let t = n as f64 * dt;
    let exact = x0 + u * t + a * (t * t / 2.0);
    assert!((x - exact).length() < 1e-12, "x = {x:?}, exact {exact:?}");
}
