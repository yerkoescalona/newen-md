//! The SPME kernels against their Lean model, and SPME against its physics.
//!
//! `spec/vectors/spme.json` holds the `NewenSpec.SPME` functions evaluated at Lean's `Float`
//! (written by `lake exe spec-vectors`, checked in, so this test needs no Lean). Every number is
//! the hex string of its 64 bits. The per-point kernels must give the same bits.
//!
//! gromosXX has no SPME, so there is no gromosXX reference. The physics oracles: the forces are
//! the exact gradient of the mesh energy (central differences agree to rounding), SPME agrees
//! with the direct Ewald sum (bit for bit the Lean Ewald model, `ewald_spec.rs`) within its error
//! estimate, and with P3M on the same grid. The measured numbers are printed
//! (`cargo test -- --nocapture`).

use newen_core::math::Vec3;
use newen_forces::nonbonded::ewald::{ewald_kspace, k_vectors};
use newen_forces::nonbonded::lattice_sum::{reciprocal_box, ChargeShape};
use newen_forces::nonbonded::p3m::{
    a2_tilde_p3m, assign_charges, assign_squared_charges, field_ik, influence_function,
    interpolate_forces, kspace_energy, potential_hat, Fft3, InfluenceInputs, Mesh,
};
use newen_forces::nonbonded::spme::{
    bspline, bspline_and_deriv, euler_factors, euler_spline_factor, spline_weights, spme_assign,
    spme_assign_squared_charges, spme_forces, spme_influence_at, spme_influence_function,
};

const SECTIONS: [&str; 9] = [
    "bspline",
    "spline_weights",
    "euler_spline_factor",
    "influence_at",
    "influence_function",
    "assign",
    "assign_squared_charges",
    "forces",
    "spme_energy",
];

/// The sections of `spme.json`, each a list of rows of floats.
fn lean_rows() -> Vec<(String, Vec<Vec<f64>>)> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../spec/vectors/spme.json");
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
    let hex = field.strip_prefix("0x").expect("hex bits in spme.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in spme.json"))
}

/// The same bits, or both NaN (Lean has a single NaN).
fn same(rust: f64, lean: f64) -> bool {
    rust.to_bits() == lean.to_bits() || (rust.is_nan() && lean.is_nan())
}

/// Compare a whole row's outputs; record the first difference.
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
    assert!(rows > 0, "{what}: no rows");
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

/// Reads a row field by field.
struct Cursor<'a> {
    row: &'a [f64],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn new(row: &'a [f64]) -> Self {
        Self { row, at: 0 }
    }
    fn f(&mut self) -> f64 {
        self.at += 1;
        self.row[self.at - 1]
    }
    fn n(&mut self) -> usize {
        self.f() as usize
    }
    fn i(&mut self) -> i64 {
        self.f() as i64
    }
    fn v3(&mut self) -> [f64; 3] {
        [self.f(), self.f(), self.f()]
    }
    fn take(&mut self, n: usize) -> Vec<f64> {
        self.at += n;
        self.row[self.at - n..self.at].to_vec()
    }
    fn rest(&self) -> &'a [f64] {
        &self.row[self.at..]
    }
    fn grid(&mut self) -> [usize; 3] {
        [self.n(), self.n(), self.n()]
    }
    /// `(rectangular, edges)`.
    fn box_(&mut self) -> (bool, [[f64; 3]; 3]) {
        (self.f() == 1.0, [self.v3(), self.v3(), self.v3()])
    }
    fn atoms(&mut self) -> (Vec<Vec3>, Vec<f64>) {
        let n = self.n();
        let mut r = Vec::new();
        let mut q = Vec::new();
        for _ in 0..n {
            let v = self.v3();
            r.push(Vec3::new(v[0], v[1], v[2]));
            q.push(self.f());
        }
        (r, q)
    }
    /// An influence case: grid, order, shape, width, box.
    fn case(&mut self) -> ([usize; 3], usize, ChargeShape, f64, bool, [[f64; 3]; 3]) {
        let grid = self.grid();
        let n = self.n();
        let shape = ChargeShape::from_nshape(self.f() as i32).expect("nshape");
        let width = self.f();
        let (rect, edges) = self.box_();
        (grid, n, shape, width, rect, edges)
    }
}

#[test]
fn bspline_and_its_derivative_match_lean_bit_for_bit() {
    let rows = section("bspline");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let (n, u) = (row[0] as usize, row[1]);
        let (m, d) = bspline_and_deriv(n, u);
        compare(
            "bspline",
            i,
            &[bspline(n, u), m, d],
            &row[2..],
            &mut mismatches,
        );
    }
    report("bspline", rows.len(), &mismatches);
}

#[test]
fn spline_weights_and_euler_factors_match_lean_bit_for_bit() {
    let rows = section("spline_weights");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let (n, u) = (row[0] as usize, row[1]);
        let (t, d) = spline_weights(n, u);
        let mut rust = t.clone();
        rust.extend(&d);
        compare("spline_weights", i, &rust, &row[2..], &mut mismatches);
        // a partition of unity, and its derivative sums to zero (to rounding)
        let sum: f64 = t.iter().sum();
        let dsum: f64 = d.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14, "n {n} u {u}: Σθ = {sum}");
        assert!(dsum.abs() < 1e-14, "n {n} u {u}: Σdθ = {dsum}");
        // the table is the recursion of `bspline` (to rounding: the arguments are rounded
        // differently)
        let w = u - u.floor();
        for j in 0..n {
            let (m, dm) = bspline_and_deriv(n, w + j as f64);
            assert!(
                (t[j] - m).abs() < 1e-14 && (d[j] - dm).abs() < 1e-13,
                "n {n} u {u} j {j}"
            );
        }
    }
    report("spline_weights", rows.len(), &mismatches);

    let rows = section("euler_spline_factor");
    for (i, row) in rows.iter().enumerate() {
        let rust = euler_spline_factor(row[0] as usize, row[1] as usize, row[2] as i64);
        compare(
            "euler_spline_factor",
            i,
            &[rust],
            &row[3..],
            &mut mismatches,
        );
    }
    report("euler_spline_factor", rows.len(), &mismatches);
}

#[test]
fn influence_function_matches_lean_bit_for_bit() {
    let rows = section("influence_at");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let (grid, n, shape, width, rect, edges) = c.case();
        let (_, recip) = reciprocal_box(&edges, rect);
        let l = [c.i(), c.i(), c.i()];
        let (g, gam) = spme_influence_at(l, grid, n, shape, width, &recip);
        let mut rust = vec![g];
        rust.extend(gam);
        compare("influence_at", i, &rust, c.rest(), &mut mismatches);
    }
    report("influence_at", rows.len(), &mismatches);

    let rows = section("influence_function");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let (grid, n, shape, width, rect, edges) = c.case();
        let (_, recip) = reciprocal_box(&edges, rect);
        let (g, gam) = spme_influence_function(grid, shape, width, &recip, &euler_factors(grid, n));
        let mut rust = g.clone();
        rust.extend(gam.iter().flatten());
        compare("influence_function", i, &rust, c.rest(), &mut mismatches);
    }
    report("influence_function", rows.len(), &mismatches);
}

#[test]
fn assignment_and_forces_match_lean_bit_for_bit() {
    for (name, squared) in [("assign", false), ("assign_squared_charges", true)] {
        let rows = section(name);
        let mut mismatches = Vec::new();
        for (i, row) in rows.iter().enumerate() {
            let mut c = Cursor::new(row);
            let grid = c.grid();
            let n = c.n();
            let (_, edges) = c.box_();
            let (r, q) = c.atoms();
            let m = Mesh::new(&edges, grid);
            let rust = if squared {
                spme_assign_squared_charges(&r, &q, &m, n)
            } else {
                spme_assign(&r, &q, &m, n)
            };
            compare(name, i, &rust, c.rest(), &mut mismatches);
            if !squared {
                // the grid carries the charges: Σρ V_cell = Σq (to rounding)
                let total: f64 = rust.iter().sum::<f64>() * m.cell_volume;
                let sum_q: f64 = q.iter().sum();
                assert!((total - sum_q).abs() < 1e-12, "row {i}: {total} vs {sum_q}");
            }
        }
        report(name, rows.len(), &mismatches);
    }

    let rows = section("forces");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let grid = c.grid();
        let n = c.n();
        let (_, edges) = c.box_();
        let (r, q) = c.atoms();
        let e0 = c.f();
        let m = Mesh::new(&edges, grid);
        let conv = c.take(m.len());
        let f = spme_forces(&r, &q, &m, n, &conv, e0);
        let rust: Vec<f64> = f.iter().flat_map(|v| [v.x, v.y, v.z]).collect();
        compare("forces", i, &rust, c.rest(), &mut mismatches);
    }
    report("forces", rows.len(), &mismatches);
}

#[test]
fn spme_kspace_energy_matches_the_lean_pipeline_within_1e_11() {
    let rows = section("spme_energy");
    let mut worst = 0.0f64;
    for row in &rows {
        let mut c = Cursor::new(row);
        let (grid, n, shape, width, rect, edges) = c.case();
        let e0 = c.f();
        let (r, q) = c.atoms();
        let lean = c.f();
        let (_, recip) = reciprocal_box(&edges, rect);
        let (ghat, _) =
            spme_influence_function(grid, shape, width, &recip, &euler_factors(grid, n));
        let mesh = Mesh::new(&edges, grid);
        let mut re = spme_assign(&r, &q, &mesh, n);
        let mut im = vec![0.0; re.len()];
        Fft3::new(grid).forward(&mut re, &mut im);
        let (e, _) = kspace_energy(&ghat, None, &re, &im, mesh.cell_volume, e0);
        let err = ((e - lean) / lean).abs();
        worst = worst.max(err);
        assert!(err < 1e-11, "SPME energy {e} vs Lean {lean}: {err:e}");
    }
    assert!(!rows.is_empty());
    println!("spme k-space energy: largest relative difference from Lean {worst:.2e}");
}

/// A small neutral system in a box: positions inside, charges.
fn system(edges: &[[f64; 3]; 3], n: usize, seed: u64) -> (Vec<Vec3>, Vec<f64>) {
    let mut s = seed;
    let mut next = move || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut r = Vec::new();
    let mut q = Vec::new();
    for i in 0..n {
        let f = [next(), next(), next()];
        let p: [f64; 3] =
            std::array::from_fn(|a| f[0] * edges[0][a] + f[1] * edges[1][a] + f[2] * edges[2][a]);
        r.push(Vec3::new(p[0], p[1], p[2]));
        q.push(if i % 2 == 0 { 0.8 } else { -0.8 } * (0.5 + next()));
    }
    let total: f64 = q.iter().sum();
    *q.last_mut().unwrap() -= total;
    (r, q)
}

const FPE: f64 = 138.9354;
const WIDTH: f64 = 0.4;

fn eps0_inv() -> f64 {
    FPE * 4.0 * std::f64::consts::PI
}

/// The SPME pipeline of one box: `Ĝ` and the grid.
struct SpmeCase {
    mesh: Mesh,
    n: usize,
    ghat: Vec<f64>,
}

impl SpmeCase {
    fn new(edges: &[[f64; 3]; 3], rect: bool, grid: usize, n: usize, shape: ChargeShape) -> Self {
        let (_, recip) = reciprocal_box(edges, rect);
        let g = [grid; 3];
        let (ghat, _) = spme_influence_function(g, shape, WIDTH, &recip, &euler_factors(g, n));
        Self {
            mesh: Mesh::new(edges, g),
            n,
            ghat,
        }
    }

    /// The mesh energy `½ V_cell ε₀⁻¹ Σ Ĝ |ρ̂|²/N`.
    fn energy(&self, r: &[Vec3], q: &[f64]) -> f64 {
        let mut re = spme_assign(r, q, &self.mesh, self.n);
        let mut im = vec![0.0; re.len()];
        Fft3::new(self.mesh.grid).forward(&mut re, &mut im);
        kspace_energy(
            &self.ghat,
            None,
            &re,
            &im,
            self.mesh.cell_volume,
            eps0_inv(),
        )
        .0
    }

    /// `(E, F, E_Ã2)`: the mesh energy, the forces, and the `Ã2` part of the A term.
    fn evaluate(&self, r: &[Vec3], q: &[f64]) -> (f64, Vec<Vec3>, f64) {
        let mut fft = Fft3::new(self.mesh.grid);
        let mut re = spme_assign(r, q, &self.mesh, self.n);
        let mut im = vec![0.0; re.len()];
        fft.forward(&mut re, &mut im);
        let (e, _) = kspace_energy(
            &self.ghat,
            None,
            &re,
            &im,
            self.mesh.cell_volume,
            eps0_inv(),
        );
        for i in 0..re.len() {
            re[i] *= self.ghat[i];
            im[i] *= self.ghat[i];
        }
        fft.backward(&mut re, &mut im);
        let f = spme_forces(r, q, &self.mesh, self.n, &re, eps0_inv());
        let sum_q2: f64 = q.iter().map(|x| x * x).sum();
        let mut rq = spme_assign_squared_charges(r, q, &self.mesh, self.n);
        let mut rq_im = vec![0.0; rq.len()];
        fft.forward(&mut rq, &mut rq_im);
        let (a2t, _) = a2_tilde_p3m(&self.ghat, None, &rq, sum_q2);
        let e_self = -a2t * sum_q2 * eps0_inv() / (8.0 * std::f64::consts::PI);
        (e, f, e_self)
    }
}

fn rms(f: &[Vec3]) -> f64 {
    (f.iter().map(|v| v.length_squared()).sum::<f64>() / f.len() as f64).sqrt()
}

#[test]
fn spme_forces_are_the_exact_gradient_of_the_mesh_energy() {
    let edges = [[2.0, 0.0, 0.0], [0.4, 2.1, 0.0], [-0.3, 0.2, 1.9]];
    let (r, q) = system(&edges, 12, 3);
    for (n, grid) in [(4, 16), (6, 16), (8, 16), (12, 16), (4, 32)] {
        let case = SpmeCase::new(&edges, false, grid, n, ChargeShape::Gaussian);
        let (_, f, _) = case.evaluate(&r, &q);
        let h = 1e-5;
        let mut worst = 0.0f64;
        for i in 0..r.len() {
            for a in 0..3 {
                let mut rp = r.clone();
                let mut rm = r.clone();
                rp[i][a] += h;
                rm[i][a] -= h;
                let g = -(case.energy(&rp, &q) - case.energy(&rm, &q)) / (2.0 * h);
                worst = worst.max((f[i][a] - g).abs());
            }
        }
        let rel = worst / rms(&f);
        println!(
            "SPME n {n} grid {grid}³: max |F − (−∇E)| {worst:.2e} (central differences, \
             h = 1e-5 nm), relative to rms F {rel:.2e}"
        );
        assert!(rel < 1e-6, "n {n}: {rel:e}");
    }
}

/// `Ξ = ½ ∂E/∂ε` for the deformation `r → (1 + ε) r` of positions and box edges together (the
/// grid coordinates of every charge stay put): the k-space virial of `kspace_energy` with the
/// SPME `Γ̂`, off-diagonal elements included, as for the Ewald sum.
#[test]
fn spme_virial_is_half_the_deformation_derivative() {
    let edges = [[2.0, 0.0, 0.0], [0.5, 1.9, 0.0], [-0.4, 0.3, 2.2]];
    let (r, q) = system(&edges, 10, 5);
    let n = 6;
    let grid = [16; 3];
    let eval = |e: &[[f64; 3]; 3], p: &[Vec3], with_gamma: bool| {
        let (_, recip) = reciprocal_box(e, false);
        let (ghat, gam) = spme_influence_function(
            grid,
            ChargeShape::Gaussian,
            WIDTH,
            &recip,
            &euler_factors(grid, n),
        );
        let mesh = Mesh::new(e, grid);
        let mut re = spme_assign(p, &q, &mesh, n);
        let mut im = vec![0.0; re.len()];
        Fft3::new(grid).forward(&mut re, &mut im);
        kspace_energy(
            &ghat,
            with_gamma.then_some(&gam[..]),
            &re,
            &im,
            mesh.cell_volume,
            eps0_inv(),
        )
    };
    let deform = |a: usize, b: usize, h: f64| {
        let m = |v: [f64; 3]| {
            let mut w = v;
            w[a] += h * v[b];
            if a != b {
                w[b] += h * v[a];
            }
            w
        };
        let e = [m(edges[0]), m(edges[1]), m(edges[2])];
        let p: Vec<Vec3> = r
            .iter()
            .map(|x| {
                let v = m([x.x, x.y, x.z]);
                Vec3::new(v[0], v[1], v[2])
            })
            .collect();
        (e, p)
    };
    let (_, xi) = eval(&edges, &r, true);
    let h = 1e-6;
    let mut worst = 0.0f64;
    for (idx, (a, b)) in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
        .into_iter()
        .enumerate()
    {
        let (ep, pp) = deform(a, b, h);
        let (em, pm) = deform(a, b, -h);
        let factor = if a == b { 0.5 } else { 0.25 };
        let fd = factor * (eval(&ep, &pp, false).0 - eval(&em, &pm, false).0) / (2.0 * h);
        let rel = (fd - xi[idx]).abs() / xi[idx].abs().max(1.0);
        worst = worst.max(rel);
        assert!(rel < 1e-6, "virial ({a}, {b}): {:e} vs FD {fd:e}", xi[idx]);
    }
    println!("SPME k-space virial vs ½ ∂E/∂ε: largest relative difference {worst:.1e}");
}

#[test]
fn a_single_charge_meets_itself_through_the_grid_as_a2_tilde_says() {
    let edges = [[2.0, 0.0, 0.0], [0.5, 1.9, 0.0], [-0.4, 0.3, 2.2]];
    for n in [4, 6, 8] {
        let case = SpmeCase::new(&edges, false, 16, n, ChargeShape::Gaussian);
        let r = [Vec3::new(0.37, 1.21, 0.93)];
        let q = [0.83];
        let (e, _, e_self) = case.evaluate(&r, &q);
        // E_γ of one charge is its self-interaction through the grid, −(Ã2 part of the A term)
        assert!(
            (e + e_self).abs() < 1e-12 * e.abs(),
            "n {n}: {e} vs {e_self}"
        );
    }
}

fn p3m_inputs(edges: &[[f64; 3]; 3], rect: bool, grid: usize, order: u8) -> (InfluenceInputs, f64) {
    let (volume, recip) = reciprocal_box(edges, rect);
    (
        InfluenceInputs {
            order,
            nfdord: 0,
            nalias: 2,
            shape: ChargeShape::Gaussian,
            width: WIDTH,
            mesh: Mesh::new(edges, [grid; 3]),
            recip,
        },
        volume,
    )
}

/// The P3M k-space energy, forces and the `Ã2` part of the A term (ik differentiation).
fn p3m(inp: &InfluenceInputs, volume: f64, r: &[Vec3], q: &[f64]) -> (f64, Vec<Vec3>, f64) {
    let inf = influence_function(inp, volume, false);
    let mesh = inp.mesh;
    let mut fft = Fft3::new(mesh.grid);
    let mut re = assign_charges(&mesh, inp.order, r, q);
    let mut im = vec![0.0; re.len()];
    fft.forward(&mut re, &mut im);
    let (e, _) = kspace_energy(&inf.ghat, None, &re, &im, mesh.cell_volume, eps0_inv());
    let sum_q2: f64 = q.iter().map(|x| x * x).sum();
    let mut rq = assign_squared_charges(&mesh, inp.order, r, q);
    let mut rq_im = vec![0.0; rq.len()];
    fft.forward(&mut rq, &mut rq_im);
    let (a2t, _) = a2_tilde_p3m(&inf.ghat, None, &rq, sum_q2);
    let e_self = -a2t * sum_q2 * eps0_inv() / (8.0 * std::f64::consts::PI);
    potential_hat(&inf.ghat, &mut re, &mut im, eps0_inv());
    let mut fld = field_ik(&mesh, &inp.recip, &re, &im);
    for (a, b) in fld.iter_mut() {
        fft.backward(a, b);
    }
    let f = interpolate_forces(&mesh, inp.order, r, q, [&fld[0].0, &fld[1].0, &fld[2].0]);
    (e, f, e_self)
}

/// The boxes of the comparisons: cubic, rectangular, triclinic.
fn boxes() -> [([[f64; 3]; 3], bool, &'static str); 3] {
    [
        (
            [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
            true,
            "cubic",
        ),
        (
            [[1.9, 0.0, 0.0], [0.0, 2.3, 0.0], [0.0, 0.0, 2.1]],
            true,
            "rectangular",
        ),
        (
            [[2.0, 0.0, 0.0], [0.5, 1.9, 0.0], [-0.4, 0.3, 2.2]],
            false,
            "triclinic",
        ),
    ]
}

#[test]
fn spme_agrees_with_the_ewald_sum_and_with_p3m() {
    let mut worst: std::collections::BTreeMap<(usize, usize), (f64, f64, f64)> = Default::default();
    for (edges, rect, name) in boxes() {
        let (volume, recip) = reciprocal_box(&edges, rect);
        let kv = k_vectors([14, 14, 14], 1e9, &recip, ChargeShape::Gaussian, WIDTH);
        for (atoms, seed) in [(8, 11), (16, 7)] {
            let (r, q) = system(&edges, atoms, seed);
            let ew = ewald_kspace(&r, &q, &kv, volume, eps0_inv());
            let sum_q2: f64 = q.iter().map(|x| x * x).sum();
            let ew_total =
                ew.energy - ew.a2_tilde * sum_q2 * eps0_inv() / (8.0 * std::f64::consts::PI);
            for grid in [16, 32] {
                for n in [4, 6, 8] {
                    let case = SpmeCase::new(&edges, rect, grid, n, ChargeShape::Gaussian);
                    let (e, f, e_self) = case.evaluate(&r, &q);
                    let df: Vec<Vec3> = f.iter().zip(&ew.forces).map(|(a, b)| *a - *b).collect();
                    let rms_df = rms(&df) / rms(&ew.forces);
                    let de = ((e + e_self - ew_total) / ew.energy).abs();
                    // the bound: B-spline interpolation of order n differentiated analytically
                    // is accurate to O((h/a)^(n−1)), h the longest grid spacing
                    let h = edges
                        .iter()
                        .map(|e| (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt())
                        .fold(0.0f64, f64::max)
                        / grid as f64;
                    let est = 0.3 * (h / WIDTH).powi(n as i32 - 1);
                    println!(
                        "{name} {atoms} charges grid {grid}³ n {n}: relative energy error \
                         {de:.1e}, relative rms force error {rms_df:.1e} (bound {est:.1e})"
                    );
                    assert!(rms_df < est, "rms {rms_df:e} vs bound {est:e}");
                    assert!(de < 0.1 * est, "energy {de:e} vs bound {est:e}");
                    let w = worst.entry((grid, n)).or_default();
                    *w = (w.0.max(de), w.1.max(rms_df), w.2.max(est));
                }
                // P3M of the same order (NASORD 4, ik, NALIAS 2) on the same grid
                let p = 4u8;
                let case = SpmeCase::new(&edges, rect, grid, p as usize, ChargeShape::Gaussian);
                let (e, f, e_self) = case.evaluate(&r, &q);
                let (inp, vol) = p3m_inputs(&edges, rect, grid, p);
                let (ep, fp, ep_self) = p3m(&inp, vol, &r, &q);
                let scale = rms(&ew.forces);
                let (ds, dp) = (
                    rms(&df_of(&f, &ew.forces)) / scale,
                    rms(&df_of(&fp, &ew.forces)) / scale,
                );
                let dsp = rms(&df_of(&f, &fp)) / scale;
                let (es, ep_err) = (
                    ((e + e_self - ew_total) / ew.energy).abs(),
                    ((ep + ep_self - ew_total) / ew.energy).abs(),
                );
                let esp = ((e + e_self - (ep + ep_self)) / ew.energy).abs();
                println!(
                    "{name} {atoms} charges grid {grid}³ order {p}: P3M vs Ewald energy \
                     {ep_err:.1e}, forces {dp:.1e}; SPME vs P3M energy {esp:.1e}, forces {dsp:.1e}"
                );
                // the two meshes differ by no more than their errors against Ewald
                assert!(dsp <= 1.5 * (ds + dp));
                assert!(esp <= 1.5 * (es + ep_err));
            }
        }
    }
    for ((grid, n), (de, df, est)) in worst {
        println!(
            "worst over boxes and systems, grid {grid}³ n {n}: relative energy error {de:.1e}, \
             relative rms force error {df:.1e}, bound {est:.1e}"
        );
    }
}

fn df_of(a: &[Vec3], b: &[Vec3]) -> Vec<Vec3> {
    a.iter().zip(b).map(|(x, y)| *x - *y).collect()
}

#[test]
#[should_panic(expected = "NSPORD must be even")]
fn an_odd_order_is_refused() {
    use newen_forces::nonbonded::ewald::LatticeSumParams;
    use newen_forces::nonbonded::spme::{Spme, SpmeParams};
    let params = LatticeSumParams {
        shape: ChargeShape::Gaussian,
        width: WIDTH,
        nk: [1, 1, 1],
        kcut: 1.0,
        na2clc: 3,
        tola2: 1e-10,
        eps_ls: 0.0,
        eps0_inv: eps0_inv(),
        four_pi_eps_i: FPE,
    };
    let _ = Spme::new(
        params,
        SpmeParams {
            grid: [16; 3],
            order: 5,
        },
    );
}
