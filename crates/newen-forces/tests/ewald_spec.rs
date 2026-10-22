//! The Ewald-sum kernels against their Lean model, bit for bit, and against their physics.
//!
//! `spec/vectors/ewald.json` holds `NewenSpec.LatticeSum.reciprocalBoxFloat`, `kVectorsFloat`,
//! `ewaldKspaceFloat`, `a2CubicFloat`, `a2ShellsFloat`, `chargeSumsFloat`, `selfATermsFloat`,
//! `selfAVirialsFloat` and `surfaceFloat` evaluated at Lean's `Float` on random small systems,
//! written by `lake exe spec-vectors` (checked in, so this test needs no Lean). Every number is
//! the hex string of its 64 bits; counts, flags and NSHAPE are floats. `sin`, `cos`, `exp` and
//! `erfc` come from the platform's C library on both sides.
//!
//! The other tests compare the Rust with itself: the k-space and surface forces are minus the
//! gradients of their energies, their virials minus half the derivative of the energy under a
//! deformation of box and positions, the k-space energy does not change under a common
//! translation, nor (in a triclinic box) when one atom moves by a lattice vector, and the wave
//! vectors are dual to the box edges.

use newen_core::math::Vec3;
use newen_core::units::EPS0_INV;
use newen_forces::nonbonded::ewald::{ewald_kspace, k_vectors, KVector};
use newen_forces::nonbonded::lattice_sum::{
    a2_cubic, a2_shells, charge_sums, lattice_vector, reciprocal_box, self_and_a_terms,
    self_and_a_virials, self_constants, surface_term, ChargeShape, Sym3, WIGNER_CUBE,
};

const SECTIONS: [&str; 9] = [
    "reciprocal_box",
    "k_vectors",
    "kspace",
    "a2_cubic",
    "a2_shells",
    "charge_sums",
    "self_and_a_terms",
    "self_and_a_virials",
    "surface_term",
];

/// The sections of `ewald.json`, each a list of rows of floats.
fn lean_rows() -> Vec<(String, Vec<Vec<f64>>)> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../spec/vectors/ewald.json");
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let mut out: Vec<(String, Vec<Vec<f64>>)> = Vec::new();
    for line in text.lines().map(str::trim) {
        for name in SECTIONS {
            if line.starts_with(&format!("\"{name}\"")) {
                out.push((name.to_string(), Vec::new()));
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
        out.last_mut()
            .expect("row outside a section")
            .1
            .push(fields);
    }
    out
}

fn section(name: &str) -> Vec<Vec<f64>> {
    lean_rows()
        .into_iter()
        .find(|(n, _)| n == name)
        .unwrap_or_else(|| panic!("no section {name}"))
        .1
}

fn float(field: &str) -> f64 {
    let hex = field.strip_prefix("0x").expect("hex bits in ewald.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in ewald.json"))
}

fn shape_of(nshape: f64) -> ChargeShape {
    ChargeShape::from_nshape(nshape as i32).unwrap_or_else(|| panic!("nshape {nshape}"))
}

/// The same bits, or both NaN (Lean has a single NaN).
fn same(rust: f64, lean: f64) -> bool {
    rust.to_bits() == lean.to_bits() || (rust.is_nan() && lean.is_nan())
}

fn compare(what: &str, i: usize, rust: &[f64], lean: &[f64], mismatches: &mut Vec<String>) {
    assert_eq!(rust.len(), lean.len(), "{what} row {i}: lengths");
    if let Some(j) = (0..rust.len()).find(|&j| !same(rust[j], lean[j])) {
        mismatches.push(format!(
            "  {what} row {i}, output {j}: rust {:e} lean {:e}",
            rust[j], lean[j]
        ));
    }
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

fn v3(row: &[f64]) -> [f64; 3] {
    [row[0], row[1], row[2]]
}

fn kvec_fields(kv: &KVector) -> [f64; 10] {
    [
        kv.k[0],
        kv.k[1],
        kv.k[2],
        kv.k2,
        kv.k2i,
        kv.ak,
        kv.gamma_hat,
        kv.k2i_gamma_hat,
        kv.ak_gamma_hat_prime,
        0.0,
    ]
}

/// The k-space case at the start of a `k_vectors` or `kspace` row: `(nk, kcut, shape, width,
/// recip)`, 15 fields.
fn k_case(row: &[f64]) -> ([i64; 3], f64, ChargeShape, f64, [[f64; 3]; 3]) {
    (
        [row[0] as i64, row[1] as i64, row[2] as i64],
        row[3],
        shape_of(row[4]),
        row[5],
        [v3(&row[6..]), v3(&row[9..]), v3(&row[12..])],
    )
}

/// `n` atoms `(r, q)` from `4n` fields.
fn atoms(fields: &[f64], n: usize) -> (Vec<Vec3>, Vec<f64>) {
    (0..n)
        .map(|i| {
            let a = &fields[4 * i..];
            (Vec3::new(a[0], a[1], a[2]), a[3])
        })
        .unzip()
}

#[test]
fn reciprocal_box_equals_lean_model_bit_for_bit() {
    let rows = section("reciprocal_box");
    assert!(rows.len() >= 100);
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row.len(), 20);
        let edges = [v3(&row[1..]), v3(&row[4..]), v3(&row[7..])];
        let (volume, recip) = reciprocal_box(&edges, row[0] == 1.0);
        let mut rust = vec![volume];
        rust.extend(recip.iter().flatten());
        compare("reciprocal_box", i, &rust, &row[10..], &mut mismatches);
    }
    report("reciprocal_box", rows.len(), &mismatches);
}

#[test]
fn k_vectors_equal_lean_model_bit_for_bit() {
    let rows = section("k_vectors");
    assert!(rows.len() >= 20);
    let mut mismatches = Vec::new();
    let mut total = 0;
    for (i, row) in rows.iter().enumerate() {
        let (nk, kcut, shape, width, recip) = k_case(row);
        let ks = k_vectors(nk, kcut, &recip, shape, width);
        let count = row[15] as usize;
        assert_eq!(ks.len(), count, "k_vectors row {i}: count");
        assert_eq!(row.len(), 16 + 9 * count);
        let rust: Vec<f64> = ks
            .iter()
            .flat_map(|kv| kvec_fields(kv)[..9].to_vec())
            .collect();
        compare("k_vectors", i, &rust, &row[16..], &mut mismatches);
        total += count;
    }
    assert!(total > 2_000, "{total} wave vectors");
    report("k_vectors", rows.len(), &mismatches);
}

#[test]
fn ewald_kspace_equals_lean_model_bit_for_bit() {
    let rows = section("kspace");
    assert!(rows.len() >= 150);
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let (nk, kcut, shape, width, recip) = k_case(row);
        let (volume, eps0_inv, n) = (row[15], row[16], row[17] as usize);
        let (pos, q) = atoms(&row[18..], n);
        let out = &row[18 + 4 * n..];
        assert_eq!(out.len(), 14 + 3 * n);
        let ks = k_vectors(nk, kcut, &recip, shape, width);
        let r = ewald_kspace(&pos, &q, &ks, volume, eps0_inv);
        let mut rust = vec![r.energy, r.a2_tilde];
        rust.extend(r.virial);
        rust.extend(r.a2_tilde_derivative);
        rust.extend(r.forces.iter().flat_map(|f| [f.x, f.y, f.z]));
        compare("kspace", i, &rust, out, &mut mismatches);
    }
    report("kspace", rows.len(), &mismatches);
}

#[test]
fn a2_cubic_equals_lean_model_bit_for_bit() {
    let rows = section("a2_cubic");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row[0].to_bits(), WIGNER_CUBE.to_bits());
        let (a2, d) = a2_cubic(row[0], row[1], row[2], row[3]);
        compare("a2_cubic", i, &[a2, d], &row[4..], &mut mismatches);
    }
    report("a2_cubic", rows.len(), &mismatches);
}

#[test]
fn a2_shells_equals_lean_model_bit_for_bit() {
    let rows = section("a2_shells");
    assert!(rows.len() >= 40);
    let (mut converged, mut stopped) = (0, 0);
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let recip = [v3(&row[2..]), v3(&row[5..]), v3(&row[8..])];
        let r = a2_shells(
            shape_of(row[0]),
            row[1],
            &recip,
            row[11],
            row[12],
            row[13] as u32,
        );
        let mut rust = vec![r.a2, r.shells as f64, if r.converged { 1.0 } else { 0.0 }];
        rust.extend(r.derivative);
        compare("a2_shells", i, &rust, &row[14..], &mut mismatches);
        if r.converged {
            converged += 1;
        } else {
            stopped += 1;
        }
    }
    assert!(
        converged > 5 && stopped > 5,
        "{converged} converged, {stopped} stopped"
    );
    report("a2_shells", rows.len(), &mismatches);
}

#[test]
fn charge_sums_equal_lean_model_bit_for_bit() {
    let rows = section("charge_sums");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let n = row[0] as usize;
        let (s, s2) = charge_sums(&row[1..1 + n]);
        compare("charge_sums", i, &[s, s2], &row[1 + n..], &mut mismatches);
    }
    report("charge_sums", rows.len(), &mismatches);
}

#[test]
fn self_and_a_terms_equal_lean_model_bit_for_bit() {
    let rows = section("self_and_a_terms");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let (e_self, e_a) =
            self_and_a_terms(row[0], row[1], row[2], row[3], row[4], row[5], row[6]);
        compare(
            "self_and_a_terms",
            i,
            &[e_self, e_a],
            &row[7..],
            &mut mismatches,
        );
    }
    report("self_and_a_terms", rows.len(), &mismatches);
}

#[test]
fn self_and_a_virials_equal_lean_model_bit_for_bit() {
    let rows = section("self_and_a_virials");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let sym = |f: &[f64]| -> Sym3 { [f[0], f[1], f[2], f[3], f[4], f[5]] };
        let (sv, cv) = self_and_a_virials(
            row[0],
            row[1],
            row[2],
            &sym(&row[3..]),
            &sym(&row[9..]),
            row[15],
        );
        let rust: Vec<f64> = sv.iter().chain(&cv).copied().collect();
        compare("self_and_a_virials", i, &rust, &row[16..], &mut mismatches);
    }
    report("self_and_a_virials", rows.len(), &mismatches);
}

#[test]
fn surface_term_equals_lean_model_bit_for_bit() {
    let rows = section("surface_term");
    let mut mismatches = Vec::new();
    let mut with_dielectric = 0;
    for (i, row) in rows.iter().enumerate() {
        let centre = Vec3::new(row[0], row[1], row[2]);
        let (eps_ls, volume, eps0_inv, n) = (row[3], row[4], row[5], row[6] as usize);
        let (pos, q) = atoms(&row[7..], n);
        let out = &row[7 + 4 * n..];
        let r = surface_term(&pos, &q, centre, eps_ls, volume, eps0_inv);
        let mut rust = vec![r.energy];
        rust.extend(r.virial);
        rust.extend(r.forces.iter().flat_map(|f| [f.x, f.y, f.z]));
        compare("surface_term", i, &rust, out, &mut mismatches);
        if eps_ls != 0.0 {
            with_dielectric += 1;
        }
    }
    assert!(with_dielectric > 50);
    report("surface_term", rows.len(), &mismatches);
}

// ─── Physics, without Lean ──────────────────────────────────────────────────────────────────

/// A small deterministic generator (SplitMix64).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    }

    fn uniform(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next()
    }
}

/// A triclinic box in GROMOS's lower-triangular form, `n` atoms in it with charges summing to
/// zero.
fn triclinic_system(rng: &mut Rng, n: usize) -> ([[f64; 3]; 3], Vec<Vec3>, Vec<f64>) {
    let (la, lb, lc) = (
        rng.uniform(1.8, 2.6),
        rng.uniform(1.8, 2.6),
        rng.uniform(1.8, 2.6),
    );
    let deg = std::f64::consts::PI / 180.0;
    let (al, be, ga) = (
        rng.uniform(70.0, 110.0) * deg,
        rng.uniform(70.0, 110.0) * deg,
        rng.uniform(70.0, 110.0) * deg,
    );
    let cx = lc * be.cos();
    let cy = lc * (al.cos() - be.cos() * ga.cos()) / ga.sin();
    let cz = (lc * lc - cx * cx - cy * cy).sqrt();
    let edges = [
        [la, 0.0, 0.0],
        [lb * ga.cos(), lb * ga.sin(), 0.0],
        [cx, cy, cz],
    ];
    let mut q: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let total: f64 = q[..n - 1].iter().sum();
    q[n - 1] = -total;
    let pos = (0..n)
        .map(|_| {
            let (fa, fb, fc) = (rng.next(), rng.next(), rng.next());
            Vec3::new(
                fa * edges[0][0] + fb * edges[1][0] + fc * edges[2][0],
                fa * edges[0][1] + fb * edges[1][1] + fc * edges[2][1],
                fa * edges[0][2] + fb * edges[1][2] + fc * edges[2][2],
            )
        })
        .collect();
    (edges, pos, q)
}

fn kspace_energy(edges: &[[f64; 3]; 3], pos: &[Vec3], q: &[f64], shape: ChargeShape) -> f64 {
    let (volume, recip) = reciprocal_box(edges, false);
    let ks = k_vectors([4, 4, 4], 1e9, &recip, shape, 0.5);
    ewald_kspace(pos, q, &ks, volume, EPS0_INV).energy
}

/// The k-space and surface forces are minus the gradients of their energies (central
/// differences, `h = 10⁻⁶` nm, to `10⁻⁵` relative or absolute: the energies are some hundred
/// kJ/mol), for a Gaussian and a polynomial shape.
#[test]
fn kspace_and_surface_forces_are_minus_gradients() {
    let mut rng = Rng(7);
    for shape in [ChargeShape::Gaussian, ChargeShape::Polynomial(3)] {
        let (edges, pos, q) = triclinic_system(&mut rng, 6);
        let (volume, recip) = reciprocal_box(&edges, false);
        let ks = k_vectors([4, 4, 4], 1e9, &recip, shape, 0.5);
        let forces = ewald_kspace(&pos, &q, &ks, volume, EPS0_INV).forces;
        let centre = Vec3::new(1.0, 1.1, 0.9);
        let surface = |p: &[Vec3]| surface_term(p, &q, centre, 61.0, volume, EPS0_INV);
        let surface_forces = surface(&pos).forces;
        let h = 1e-6;
        for i in 0..pos.len() {
            for c in 0..3 {
                let shifted = |d: f64| {
                    let mut p = pos.clone();
                    p[i][c] += d;
                    p
                };
                let (plus, minus) = (shifted(h), shifted(-h));
                let fd = -(ewald_kspace(&plus, &q, &ks, volume, EPS0_INV).energy
                    - ewald_kspace(&minus, &q, &ks, volume, EPS0_INV).energy)
                    / (2.0 * h);
                let f = forces[i][c];
                assert!(
                    (fd - f).abs() <= 1e-5 * f.abs().max(1.0),
                    "{shape:?}: k-space force on atom {i}, {c}: {f:e} vs FD {fd:e}"
                );
                let fd = -(surface(&plus).energy - surface(&minus).energy) / (2.0 * h);
                let f = surface_forces[i][c];
                assert!(
                    (fd - f).abs() <= 1e-5 * f.abs().max(1.0),
                    "surface force on atom {i}, {c}: {f:e} vs FD {fd:e}"
                );
            }
        }
    }
}

/// `Ξ = ½ ∂E/∂ε` for the deformation `r → (1 + ε) r` of positions and box edges together: the
/// k-space virial (with the same integer triples) and the surface virial, off-diagonal elements
/// included. (gromosXX adds the diagonal of the surface virial without the factor −2 of every
/// other term; newen stores `−2 Ξ` for all.)
#[test]
fn kspace_and_surface_virials_are_half_the_deformation_derivative() {
    let mut rng = Rng(11);
    let (edges, pos, q) = triclinic_system(&mut rng, 5);
    let shape = ChargeShape::Gaussian;
    let deform = |e: [[f64; 3]; 3], a: usize, b: usize, h: f64| {
        let m = |v: [f64; 3]| {
            let mut w = v;
            w[a] += h * v[b];
            if a != b {
                w[b] += h * v[a];
            }
            w
        };
        let pv = |p: &Vec3| {
            let v = m([p.x, p.y, p.z]);
            Vec3::new(v[0], v[1], v[2])
        };
        (
            [m(e[0]), m(e[1]), m(e[2])],
            pos.iter().map(pv).collect::<Vec<_>>(),
            pv(&Vec3::new(1.0, 1.1, 0.9)),
        )
    };
    let energies = |e: [[f64; 3]; 3], p: &[Vec3], c: Vec3| {
        let (volume, recip) = reciprocal_box(&e, false);
        let ks = k_vectors([4, 4, 4], 1e9, &recip, shape, 0.5);
        (
            ewald_kspace(p, &q, &ks, volume, EPS0_INV).energy,
            surface_term(p, &q, c, 61.0, volume, EPS0_INV).energy,
        )
    };
    let (volume, recip) = reciprocal_box(&edges, false);
    let ks = k_vectors([4, 4, 4], 1e9, &recip, shape, 0.5);
    let xi_k = ewald_kspace(&pos, &q, &ks, volume, EPS0_INV).virial;
    let xi_s = surface_term(&pos, &q, Vec3::new(1.0, 1.1, 0.9), 61.0, volume, EPS0_INV).virial;
    let h = 1e-6;
    for (idx, (a, b)) in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
        .into_iter()
        .enumerate()
    {
        let (ep, pp, cp) = deform(edges, a, b, h);
        let (em, pm, cm) = deform(edges, a, b, -h);
        let (kp, sp) = energies(ep, &pp, cp);
        let (km, sm) = energies(em, &pm, cm);
        // a symmetric off-diagonal deformation moves both ε_ab and ε_ba
        let factor = if a == b { 0.5 } else { 0.25 };
        for (what, xi, fd) in [
            ("k space", xi_k[idx], factor * (kp - km) / (2.0 * h)),
            ("surface", xi_s[idx], factor * (sp - sm) / (2.0 * h)),
        ] {
            assert!(
                (fd - xi).abs() <= 1e-5 * xi.abs().max(1.0),
                "{what} virial ({a}, {b}): {xi:e} vs FD {fd:e}"
            );
        }
    }
}

/// The k-space energy does not change when every atom moves by the same vector (within
/// `10⁻¹⁰` relative: the rounding differs).
#[test]
fn kspace_energy_is_translation_invariant() {
    let mut rng = Rng(3);
    for shape in [ChargeShape::Gaussian, ChargeShape::Polynomial(0)] {
        let (edges, pos, q) = triclinic_system(&mut rng, 8);
        let e0 = kspace_energy(&edges, &pos, &q, shape);
        let d = Vec3::new(0.37, -0.21, 0.55);
        let moved: Vec<Vec3> = pos.iter().map(|p| *p + d).collect();
        let e1 = kspace_energy(&edges, &moved, &q, shape);
        assert!((e1 - e0).abs() <= 1e-10 * e0.abs(), "{e0:e} vs {e1:e}");
    }
}

/// In a triclinic box the wave vectors are dual to the edges, `k · a = 2π l_x` and so on, so one
/// atom moved by a lattice vector leaves the k-space energy unchanged (gromosXX's wave vectors,
/// the transposed product, fail both).
#[test]
fn triclinic_wave_vectors_are_dual_to_the_box_and_the_energy_is_periodic() {
    let mut rng = Rng(5);
    let (edges, pos, q) = triclinic_system(&mut rng, 6);
    let (_, recip) = reciprocal_box(&edges, false);
    let two_pi = 2.0 * std::f64::consts::PI;
    for l in [[1, 0, 0], [0, 1, 0], [0, 0, 1], [2, -1, 3]] {
        let k = lattice_vector(&recip, l);
        for (j, e) in edges.iter().enumerate() {
            let dot = k[0] * e[0] + k[1] * e[1] + k[2] * e[2];
            assert!(
                (dot - two_pi * l[j] as f64).abs() < 1e-12,
                "l {l:?}, edge {j}: {dot}"
            );
        }
    }
    let e0 = kspace_energy(&edges, &pos, &q, ChargeShape::Gaussian);
    let mut moved = pos.clone();
    for c in 0..3 {
        moved[2][c] += edges[0][c] - edges[1][c] + 2.0 * edges[2][c];
    }
    let e1 = kspace_energy(&edges, &moved, &q, ChargeShape::Gaussian);
    assert!((e1 - e0).abs() <= 1e-10 * e0.abs(), "{e0:e} vs {e1:e}");
}

/// The cube's `A2` from the Wigner constant agrees with the shell sum to the shell sum's
/// tolerance, for the Gaussian.
#[test]
fn a2_cubic_agrees_with_the_shell_sum() {
    for (edge, width) in [(2.0, 0.4), (1.8652, 0.4), (3.0, 0.6)] {
        let volume = edge * edge * edge;
        let (a1, a3) = self_constants(ChargeShape::Gaussian, width, volume);
        let (cubic, _) = a2_cubic(WIGNER_CUBE, edge, a1, a3);
        let edges = [[edge, 0.0, 0.0], [0.0, edge, 0.0], [0.0, 0.0, edge]];
        let (v, recip) = reciprocal_box(&edges, true);
        let shells = a2_shells(ChargeShape::Gaussian, width, &recip, v, 1e-12, 100);
        assert!(shells.converged);
        assert!(
            (shells.a2 - cubic).abs() <= 1e-7 * cubic.abs(),
            "edge {edge}: shells {} vs cubic {cubic}",
            shells.a2
        );
    }
}
