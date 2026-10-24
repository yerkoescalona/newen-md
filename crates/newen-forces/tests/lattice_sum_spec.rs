//! The lattice-sum kernels against their Lean model, bit for bit, and against their physics.
//!
//! `spec/vectors/lattice_sum.json` holds `NewenSpec.LatticeSum.shapeSwitchFloat`,
//! `shapeFourierFloat`, `ljLsPairFloat`, `lsExcludedPairFloat` and `selfConstantsFloat` evaluated
//! at Lean's `Float` (IEEE 754 binary64), written by `lake exe spec-vectors` (checked in, so this
//! test needs no Lean). Every number is the hex string of its 64 bits. Besides `+ − × ÷` and
//! `sqrt`, which IEEE 754 rounds correctly, the kernels call `erfc`, `exp`, `sin` and `cos`, which
//! it does not; Rust and Lean both call the platform's C library for them, so the bits agree on
//! one machine. A NaN matches any NaN (Lean has a single NaN).
//!
//! The other tests compare the Rust with itself and with GROMOS Vol. 2's tables: `η'` is the
//! derivative of `η`, `γ̂'` of `γ̂`, the pair force of the pair energy (finite differences);
//! `η(0) = 1` for every shape; and `γ̂`, `A1`, `A3` follow from the tabulated shape `γ` and
//! switch function `η` (quadrature).

use newen_core::units::{four_pi_eps_i, EPS0_INV};
use newen_forces::nonbonded::lattice_sum::{
    lj_ls_pair, ls_excluded_pair, self_constants, shape_fourier, shape_switch, ChargeShape,
};

/// `piFloat` in `spec/NewenSpec/LatticeSum.lean`.
const PI_BITS: u64 = 0x4009_21fb_5444_2d18;
/// gromosXX `math::four_pi_eps_i = 138.9354` (`math/math.cc:43`).
const FOUR_PI_EPS_I_BITS: u64 = 0x4061_5dee_cbfb_15b5;
/// gromosXX `math::eps0_i = 1745.9137` (`math/math.cc:45`).
const EPS0_INV_BITS: u64 = 0x409b_47a7_a0f9_096c;

/// Every NSHAPE, `−1..=10`.
fn all_shapes() -> Vec<ChargeShape> {
    (-1..=10)
        .map(|n| ChargeShape::from_nshape(n).unwrap())
        .collect()
}

/// The five sections of `lattice_sum.json`, each a list of rows of floats.
#[derive(Default)]
struct Rows {
    shape_switch: Vec<Vec<f64>>,
    shape_fourier: Vec<Vec<f64>>,
    lj_ls_pair: Vec<Vec<f64>>,
    ls_excluded_pair: Vec<Vec<f64>>,
    self_constants: Vec<Vec<f64>>,
}

const SECTIONS: [&str; 5] = [
    "shape_switch",
    "shape_fourier",
    "lj_ls_pair",
    "ls_excluded_pair",
    "self_constants",
];

fn lean_rows() -> Rows {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/lattice_sum.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let mut rows = Rows::default();
    let mut section = "";
    for line in text.lines().map(str::trim) {
        for name in SECTIONS {
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
            "shape_switch" => rows.shape_switch.push(fields),
            "shape_fourier" => rows.shape_fourier.push(fields),
            "lj_ls_pair" => rows.lj_ls_pair.push(fields),
            "ls_excluded_pair" => rows.ls_excluded_pair.push(fields),
            "self_constants" => rows.self_constants.push(fields),
            _ => panic!("row outside a section: {line}"),
        }
    }
    rows
}

fn float(field: &str) -> f64 {
    let hex = field
        .strip_prefix("0x")
        .expect("hex bits in lattice_sum.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in lattice_sum.json"))
}

/// The shape of a row's `nshape` field.
fn shape_of(nshape: f64) -> ChargeShape {
    assert_eq!(nshape, nshape.trunc(), "nshape {nshape}");
    ChargeShape::from_nshape(nshape as i32).unwrap_or_else(|| panic!("nshape {nshape}"))
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
fn constants_are_the_bits_the_lean_model_and_gromosxx_use() {
    assert_eq!(std::f64::consts::PI.to_bits(), PI_BITS);
    assert_eq!(four_pi_eps_i.to_bits(), FOUR_PI_EPS_I_BITS);
    assert_eq!(EPS0_INV.to_bits(), EPS0_INV_BITS);
}

#[test]
fn nshape_maps_to_the_twelve_shapes() {
    assert_eq!(ChargeShape::from_nshape(-1), Some(ChargeShape::Gaussian));
    for n in 0..=10 {
        assert_eq!(
            ChargeShape::from_nshape(n),
            Some(ChargeShape::Polynomial(n as u8))
        );
    }
    for n in [-2, 11, 100, i32::MIN] {
        assert_eq!(ChargeShape::from_nshape(n), None);
    }
}

#[test]
fn shape_switch_equals_lean_model_bit_for_bit() {
    let rows = lean_rows().shape_switch;
    assert!(rows.len() >= 2_000, "{} shape_switch rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 4, "shape_switch row {row:?}");
        let (eta, d_eta) = shape_switch(shape_of(row[0]), row[1]);
        check(row, &[eta, d_eta], &row[2..], &mut mismatches);
    }
    report("shape_switch", rows.len(), &mismatches);
}

#[test]
fn shape_fourier_equals_lean_model_bit_for_bit() {
    let rows = lean_rows().shape_fourier;
    assert!(rows.len() >= 2_000, "{} shape_fourier rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 4, "shape_fourier row {row:?}");
        let (g, dg) = shape_fourier(shape_of(row[0]), row[1]);
        check(row, &[g, dg], &row[2..], &mut mismatches);
    }
    report("shape_fourier", rows.len(), &mismatches);
}

#[test]
fn lj_ls_pair_equals_lean_model_bit_for_bit() {
    let rows = lean_rows().lj_ls_pair;
    assert!(rows.len() >= 2_000, "{} lj_ls_pair rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 10, "lj_ls_pair row {row:?}");
        let (r2, c6, c12, q_eps, shape) = (row[0], row[1], row[2], row[3], shape_of(row[4]));
        let (width, inv_width) = (row[5], row[6]);
        assert_eq!(inv_width.to_bits(), (1.0 / width).to_bits());
        let (f_mag, e_lj, e_ls) = lj_ls_pair(r2, c6, c12, q_eps, shape, inv_width);
        check(row, &[f_mag, e_lj, e_ls], &row[7..], &mut mismatches);
    }
    report("lj_ls_pair", rows.len(), &mismatches);
}

#[test]
fn ls_excluded_pair_equals_lean_model_bit_for_bit() {
    let rows = lean_rows().ls_excluded_pair;
    assert!(rows.len() >= 1_000, "{} ls_excluded_pair rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 7, "ls_excluded_pair row {row:?}");
        let (r2, q_eps, shape, width, inv_width) =
            (row[0], row[1], shape_of(row[2]), row[3], row[4]);
        assert_eq!(inv_width.to_bits(), (1.0 / width).to_bits());
        let (f_mag, e_ls) = ls_excluded_pair(r2, q_eps, shape, inv_width);
        check(row, &[f_mag, e_ls], &row[5..], &mut mismatches);
    }
    report("ls_excluded_pair", rows.len(), &mismatches);
}

#[test]
fn self_constants_equal_lean_model_bit_for_bit() {
    let rows = lean_rows().self_constants;
    assert!(rows.len() >= 150, "{} self_constants rows", rows.len());
    let mut mismatches = Vec::new();
    for row in &rows {
        assert_eq!(row.len(), 5, "self_constants row {row:?}");
        let (a1, a3) = self_constants(shape_of(row[0]), row[1], row[2]);
        check(row, &[a1, a3], &row[3..], &mut mismatches);
    }
    report("self_constants", rows.len(), &mismatches);
}

/// The generated rows cover every shape, the support edge on both sides, `κ = 0` and pairs at
/// `r² = a²`.
#[test]
fn rows_cover_the_edges() {
    let rows = lean_rows();
    for nshape in -1..=10 {
        let n = nshape as f64;
        let sw: Vec<_> = rows.shape_switch.iter().filter(|r| r[0] == n).collect();
        for xi in [0.0, 1.0, 1.0 - f64::EPSILON / 2.0, 1.0 + f64::EPSILON] {
            assert!(sw.iter().any(|r| r[1] == xi), "nshape {nshape}, xi {xi}");
        }
        assert!(rows.shape_fourier.iter().any(|r| r[0] == n && r[1] == 0.0));
        assert!(rows.self_constants.iter().any(|r| r[0] == n));
    }
    assert!(rows.lj_ls_pair.iter().any(|r| r[0] == r[5] * r[5]));
    assert!(rows.ls_excluded_pair.iter().any(|r| r[0] == r[3] * r[3]));
}

/// `η'` is the derivative of `η`: central differences at `h = 10⁻⁶` inside the support (the
/// Gaussian on `[0.01, 3]`), and a vanishing slope at the edge `ξ = 1`.
#[test]
fn switch_derivative_is_derivative_of_switch() {
    let h = 1e-6;
    for shape in all_shapes() {
        let top = match shape {
            ChargeShape::Gaussian => 3.0,
            ChargeShape::Polynomial(_) => 0.99,
        };
        let mut xi = 0.01;
        while xi <= top {
            let fd = (shape_switch(shape, xi + h).0 - shape_switch(shape, xi - h).0) / (2.0 * h);
            let d_eta = shape_switch(shape, xi).1;
            assert!(
                (fd - d_eta).abs() <= 1e-7 * d_eta.abs().max(1.0),
                "{shape:?} at xi = {xi}: eta' = {d_eta}, finite difference {fd}"
            );
            xi += 0.01;
        }
        if let ChargeShape::Polynomial(_) = shape {
            let fd = (shape_switch(shape, 1.0 + h).0 - shape_switch(shape, 1.0 - h).0) / (2.0 * h);
            assert_eq!(
                shape_switch(shape, 1.0),
                (0.0, 0.0),
                "{shape:?} at the edge"
            );
            assert!(fd.abs() <= 1e-5, "{shape:?}: slope {fd} across the edge");
            assert_eq!(
                shape_switch(shape, 1.5),
                (0.0, 0.0),
                "{shape:?} beyond the edge"
            );
        }
    }
}

/// `γ̂'` is the derivative of `γ̂`. The polynomial shapes cancel catastrophically for small `κ`
/// (`κ^(N+3)` in the denominator), so they are checked for `κ ≥ 5`, the Gaussian from `0.01`.
#[test]
fn fourier_derivative_is_derivative_of_fourier() {
    for shape in all_shapes() {
        let bottom = match shape {
            ChargeShape::Gaussian => 0.01,
            ChargeShape::Polynomial(_) => 5.0,
        };
        let mut kappa = bottom;
        while kappa <= 40.0 {
            let h = 1e-5 * kappa;
            let fd =
                (shape_fourier(shape, kappa + h).0 - shape_fourier(shape, kappa - h).0) / (2.0 * h);
            let d_g = shape_fourier(shape, kappa).1;
            assert!(
                (fd - d_g).abs() <= 1e-7,
                "{shape:?} at kappa = {kappa}: gamma_hat' = {d_g}, finite difference {fd}"
            );
            kappa += 0.37;
        }
    }
}

/// Whether a central difference at `r` would straddle the support edge `r = a`, where `η''` of a
/// polynomial shape jumps and the difference is off by `O(h)`.
fn near_edge(r: f64, inv_width: f64) -> bool {
    (r * inv_width - 1.0).abs() < 1e-4
}

/// The pair force is `−(dE/dr)/r` of the pair energy, LJ and lattice sum together, and so is the
/// excluded-pair force: central differences in `r` on the generated pairs with `r² ≥ 0.01` nm²,
/// away from the support edge (`switch_derivative_is_derivative_of_switch` covers the edge).
#[test]
fn pair_forces_are_minus_gradients() {
    let rows = lean_rows();
    let mut checked = 0;
    for row in &rows.lj_ls_pair {
        let (r2, c6, c12, q_eps, shape, inv_width) =
            (row[0], row[1], row[2], row[3], shape_of(row[4]), row[6]);
        if r2 < 0.01 {
            continue;
        }
        let r = r2.sqrt();
        if near_edge(r, inv_width) {
            continue;
        }
        let h = 1e-7 * r;
        let energy = |r: f64| {
            let (_, e_lj, e_ls) = lj_ls_pair(r * r, c6, c12, q_eps, shape, inv_width);
            e_lj + e_ls
        };
        let fd = -(energy(r + h) - energy(r - h)) / (2.0 * h) / r;
        let (f_mag, e_lj, e_ls) = lj_ls_pair(r2, c6, c12, q_eps, shape, inv_width);
        let scale = f_mag.abs() + (e_lj.abs() + e_ls.abs()) / r2;
        assert!(
            (fd - f_mag).abs() <= 1e-6 * scale,
            "lj_ls_pair {row:?}: f_mag = {f_mag}, finite difference {fd}"
        );
        checked += 1;
    }
    for row in &rows.ls_excluded_pair {
        let (r2, q_eps, shape, inv_width) = (row[0], row[1], shape_of(row[2]), row[4]);
        if r2 < 0.01 {
            continue;
        }
        let r = r2.sqrt();
        if near_edge(r, inv_width) {
            continue;
        }
        let h = 1e-7 * r;
        let energy = |r: f64| ls_excluded_pair(r * r, q_eps, shape, inv_width).1;
        let fd = -(energy(r + h) - energy(r - h)) / (2.0 * h) / r;
        let (f_mag, e_ls) = ls_excluded_pair(r2, q_eps, shape, inv_width);
        let scale = f_mag.abs() + e_ls.abs() / r2 + q_eps.abs() / r2;
        assert!(
            (fd - f_mag).abs() <= 1e-6 * scale,
            "ls_excluded_pair {row:?}: f_mag = {f_mag}, finite difference {fd}"
        );
        checked += 1;
    }
    assert!(checked >= 1_500, "{checked} rows checked");
}

/// `η(0) = 1` for every shape (the full Coulomb term at contact), and the order-0 polynomial
/// against its expanded form `1 − 3ξ/2 + ξ³/2`.
#[test]
fn switch_values_from_vol2() {
    for shape in all_shapes() {
        assert_eq!(shape_switch(shape, 0.0).0, 1.0, "{shape:?}: eta(0)");
    }
    for xi in [0.1, 0.25, 0.5, 0.75, 0.9] {
        let (eta, d_eta) = shape_switch(ChargeShape::Polynomial(0), xi);
        assert!((eta - (1.0 - 1.5 * xi + 0.5 * xi * xi * xi)).abs() <= 1e-15);
        assert!((d_eta - (1.5 * xi * xi - 1.5)).abs() <= 1e-15);
    }
}

/// `π γ(ξ)` of Vol. 2, table `ls_charge_shaping`.
fn pi_gamma(shape: ChargeShape, xi: f64) -> f64 {
    let m = 1.0 - xi;
    match shape {
        ChargeShape::Gaussian => (-xi * xi).exp() / std::f64::consts::PI.sqrt(),
        _ if m < 0.0 => 0.0,
        ChargeShape::Polynomial(order) => match order {
            0 => 3.0 / 4.0,
            1 => 3.0 * m,
            2 => 15.0 / 2.0 * m.powi(2),
            3 => 15.0 / 4.0 * m.powi(2) * (1.0 + 2.0 * xi),
            4 => 105.0 / 16.0 * m.powi(3) * (3.0 * xi + 1.0),
            5 => 21.0 / 2.0 * m.powi(4) * (4.0 * xi + 1.0),
            6 => 63.0 / 8.0 * m.powi(4) * (5.0 * xi * xi + 4.0 * xi + 1.0),
            7 => 45.0 / 4.0 * m.powi(5) * (8.0 * xi * xi + 5.0 * xi + 1.0),
            8 => 165.0 / 32.0 * m.powi(6) * (35.0 * xi * xi + 18.0 * xi + 3.0),
            9 => 165.0 / 64.0 * m.powi(6) * (64.0 * xi.powi(3) + 69.0 * xi * xi + 30.0 * xi + 5.0),
            10 => {
                2145.0 / 128.0 * m.powi(7) * (21.0 * xi.powi(3) + 19.0 * xi * xi + 7.0 * xi + 1.0)
            },
            _ => unreachable!(),
        },
    }
}

/// Simpson's rule on `[0, b]` with `n` (even) intervals.
fn simpson(f: impl Fn(f64) -> f64, b: f64, n: usize) -> f64 {
    let h = b / n as f64;
    let inner: f64 = (1..n)
        .map(|i| f(i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 })
        .sum();
    (f(0.0) + inner + f(b)) * h / 3.0
}

/// The integration range: the support, or far enough into the Gaussian's tail.
fn support(shape: ChargeShape) -> f64 {
    match shape {
        ChargeShape::Gaussian => 8.0,
        ChargeShape::Polynomial(_) => 1.0,
    }
}

/// The tabulated `γ̂` is the Fourier transform of the tabulated shape (Vol. 2,
/// `latsum_gamma_ft` at `a = 1`): `γ̂(κ) = 4/κ ∫ ξ sin(κξ) πγ(ξ) dξ`, and the shape carries unit
/// charge, `4 ∫ ξ² πγ(ξ) dξ = 1`. Checked where the closed form does not cancel (`κ ≥ 2`).
#[test]
fn fourier_coefficient_is_transform_of_shape() {
    for shape in all_shapes() {
        let b = support(shape);
        let charge = 4.0 * simpson(|xi| xi * xi * pi_gamma(shape, xi), b, 4000);
        assert!((charge - 1.0).abs() <= 1e-10, "{shape:?}: charge {charge}");
        for kappa in [2.0, 3.5, 5.0, 10.0, 20.0, 35.0] {
            let quad =
                4.0 / kappa * simpson(|xi| xi * (kappa * xi).sin() * pi_gamma(shape, xi), b, 8000);
            let (g, _) = shape_fourier(shape, kappa);
            assert!(
                (g - quad).abs() <= 1e-6,
                "{shape:?} at kappa = {kappa}: gamma_hat = {g}, quadrature {quad}"
            );
        }
    }
}

/// The self-term constants follow from `η` (Vol. 2, `latsum_a1_def`, `latsum_a3_def`):
/// `−V A1/(π a²) = 4 ∫ ξ η(ξ) dξ` and `−a A3 = −η'(0)`, checked at `a = 1`, `V = 1`.
#[test]
fn self_constants_follow_from_switch() {
    for shape in all_shapes() {
        let (a1, a3) = self_constants(shape, 1.0, 1.0);
        let a1_tilde = -a1 / std::f64::consts::PI;
        let quad = 4.0 * simpson(|xi| xi * shape_switch(shape, xi).0, support(shape), 4000);
        assert!(
            (a1_tilde - quad).abs() <= 1e-10,
            "{shape:?}: table A1 {a1_tilde}, quadrature {quad}"
        );
        let d_eta0 = shape_switch(shape, 0.0).1;
        assert!(
            (a3 - d_eta0).abs() <= 1e-15,
            "{shape:?}: A3 = {a3}, eta'(0) = {d_eta0}"
        );
    }
}
