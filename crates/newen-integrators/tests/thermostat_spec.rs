//! The thermostat kernels against their Lean model, bit for bit.
//!
//! `spec/vectors/thermostat.json` holds `NewenSpec.Thermostat.berendsenModel` and
//! `noseHooverModel` evaluated at Lean's `Float` (IEEE 754 binary64) on random inputs, written by
//! `lake exe spec-vectors` (checked in, so this test needs no Lean). Every number is the hex
//! string of its 64 bits. The kernels use only `+ − × ÷` and `sqrt`, which IEEE 754 rounds
//! correctly in Rust and in Lean (C) alike, so the Rust must return the same bits as Lean.
//!
//! One exception: Lean's `Float` has a single NaN (`Float.toBits` writes every NaN as
//! `0x7ff8000000000000`), while x86 produces `0xfff8000000000000` for `∞ − ∞`. A NaN therefore
//! matches any NaN; every other value must match bit for bit.

use newen_core::units::kB;
use newen_integrators::algorithms::{berendsen_scale, nose_hoover_scale, ZERO_TEMPERATURE};

/// `kBFloat` in `spec/NewenSpec/Thermostat.lean`.
const KB_BITS: u64 = 0x3f81_0725_3848_c6c3;
/// `zeroTemperatureFloat` in `spec/NewenSpec/Thermostat.lean`.
const ZERO_TEMPERATURE_BITS: u64 = 0x3d71_9799_812d_ea11;

/// The Berendsen rows and the Nosé–Hoover rows, each a list of fields (hex bits or `nhc`).
fn lean_rows() -> (Vec<Vec<String>>, Vec<Vec<String>>) {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/thermostat.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let (mut berendsen, mut nose_hoover) = (Vec::new(), Vec::new());
    let mut section = "";
    for line in text.lines().map(str::trim) {
        if line.starts_with("\"berendsen\"") {
            section = "berendsen";
        } else if line.starts_with("\"nose_hoover\"") {
            section = "nose_hoover";
        }
        let line = line.trim_end_matches(',');
        if !(line.starts_with('[') && line.ends_with(']')) {
            continue;
        }
        let fields = line[1..line.len() - 1]
            .split(',')
            .map(|f| f.trim().trim_matches('"').to_string())
            .collect();
        match section {
            "berendsen" => berendsen.push(fields),
            "nose_hoover" => nose_hoover.push(fields),
            _ => panic!("row outside a section: {line}"),
        }
    }
    (berendsen, nose_hoover)
}

fn float(field: &str) -> f64 {
    let hex = field
        .strip_prefix("0x")
        .expect("hex bits in thermostat.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in thermostat.json"))
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
    assert_eq!(kB.to_bits(), KB_BITS);
    assert_eq!(ZERO_TEMPERATURE.to_bits(), ZERO_TEMPERATURE_BITS);
}

#[test]
fn berendsen_equals_lean_model_bit_for_bit() {
    let (rows, _) = lean_rows();
    assert!(rows.len() >= 5_000, "{} Berendsen rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 6, "Berendsen row {row:?}");
        let v: Vec<f64> = row.iter().map(|f| float(f)).collect();
        let (ekin, dof, temperature, tau, dt, lean) = (v[0], v[1], v[2], v[3], v[4], v[5]);
        let rust = berendsen_scale(ekin, dof, temperature, tau, dt);
        if !same(rust, lean) {
            mismatches.push(format!(
                "  ekin={ekin:e} dof={dof} T0={temperature} tau={tau} dt={dt}: \
                 rust={rust:e} lean={lean:e} ({} ulp)",
                ulps(rust, lean)
            ));
        }
    }
    report("berendsen_scale", rows.len(), &mismatches);
}

#[test]
fn nose_hoover_equals_lean_model_bit_for_bit() {
    let (_, rows) = lean_rows();
    assert!(rows.len() >= 2_500, "{} Nosé–Hoover rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        let nhc: usize = row[5].parse().expect("nhc in thermostat.json");
        assert_eq!(row.len(), 7 + 2 * nhc, "Nosé–Hoover row {row:?}");
        let v: Vec<f64> = row[..5].iter().map(|f| float(f)).collect();
        let (ekin, dof, temperature, tau, dt) = (v[0], v[1], v[2], v[3], v[4]);
        let zeta_in: Vec<f64> = row[6..6 + nhc].iter().map(|f| float(f)).collect();
        let zeta_lean: Vec<f64> = row[6 + nhc..6 + 2 * nhc].iter().map(|f| float(f)).collect();
        let scale_lean = float(&row[6 + 2 * nhc]);

        let mut zeta = zeta_in.clone();
        let scale = nose_hoover_scale(&mut zeta, ekin, dof, temperature, tau, dt);
        let zeta_ok = zeta.iter().zip(&zeta_lean).all(|(&r, &l)| same(r, l));
        if !(zeta_ok && same(scale, scale_lean)) {
            let worst = zeta
                .iter()
                .zip(&zeta_lean)
                .map(|(&r, &l)| ulps(r, l))
                .chain([ulps(scale, scale_lean)])
                .max()
                .unwrap_or(0);
            mismatches.push(format!(
                "  ekin={ekin:e} dof={dof} T0={temperature} tau={tau} dt={dt} zeta_in={zeta_in:?}: \
                 rust=({zeta:?}, {scale:e}) lean=({zeta_lean:?}, {scale_lean:e}) ({worst} ulp)"
            ));
        }
    }
    report("nose_hoover_scale", rows.len(), &mismatches);
}

/// τ = Δt: one Berendsen step reaches the target temperature (as `berendsen_relax` says over ℝ).
#[test]
fn berendsen_tau_equal_to_dt_reaches_the_target() {
    let (dof, t0, dt) = (3.0, 300.0, 0.002);
    let ekin = 0.5 * dof * kB * 450.0;
    let scale = berendsen_scale(ekin, dof, t0, dt, dt);
    let t_new = scale * scale * 2.0 * ekin / (dof * kB);
    assert!((t_new - t0).abs() < 1e-9, "T after one step = {t_new}");
}

#[test]
fn negative_tau_leaves_zeta_and_velocities_alone() {
    assert_eq!(berendsen_scale(10.0, 3.0, 300.0, -1.0, 0.002), 1.0);
    let mut zeta = vec![0.25, -1.5, 3.0];
    let scale = nose_hoover_scale(&mut zeta, 10.0, 3.0, 300.0, -1.0, 0.002);
    assert_eq!(scale, 1.0);
    assert_eq!(zeta, [0.25, -1.5, 3.0]);
}

/// No kinetic energy: the bath falls back to its reference temperature, so Berendsen does not
/// scale and Nosé–Hoover does not change ζ.
#[test]
fn a_cold_bath_falls_back_to_the_reference_temperature() {
    assert_eq!(berendsen_scale(0.0, 3.0, 300.0, 0.1, 0.002), 1.0);
    let mut zeta = vec![0.0];
    assert_eq!(
        nose_hoover_scale(&mut zeta, 0.0, 3.0, 300.0, 0.1, 0.002),
        1.0
    );
    assert_eq!(zeta, [0.0]);
}
