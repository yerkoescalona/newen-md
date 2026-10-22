//! The P3M kernels against their Lean model, and P3M against its physics.
//!
//! `spec/vectors/p3m.json` holds the `NewenSpec.P3M` functions evaluated at Lean's `Float`
//! (written by `lake exe spec-vectors`, checked in, so this test needs no Lean). Every number is
//! the hex string of its 64 bits. The per-point kernels must give the same bits. The transforms
//! do not: PhastFT's butterflies are not the Lean model's direct sum, so `Fft3` and the whole
//! k-space energy are compared with `dftFloat` within a tolerance, and the measured maximum is
//! printed (`cargo test -- --nocapture`).
//!
//! Without Lean: P3M against the direct Ewald sum of phase 2 on small systems, within the rms
//! force-error estimate, and the P3M forces against central differences of the P3M energy.

use newen_core::math::Vec3;
use newen_forces::nonbonded::ewald::{ewald_kspace, k_vectors};
use newen_forces::nonbonded::lattice_sum::{reciprocal_box, ChargeShape, Sym3};
use newen_forces::nonbonded::p3m::{
    a2_tilde_p3m, assign_charges, assign_squared_charges, assignment_1d, assignment_weights,
    fd_coefficients, field_fd, field_ik, force_error_estimate, influence_at,
    influence_box_correction, influence_function, interpolate_forces, kspace_energy,
    nearest_grid_point, potential_hat, quality_at, selfterm_fp, sinc_power, Fft3, InfluenceInputs,
    Mesh, M3,
};

const SECTIONS: [&str; 20] = [
    "assignment_1d",
    "weights",
    "selfterm_fp",
    "fd_coefficients",
    "sinc_power",
    "mesh",
    "influence_at",
    "quality_at",
    "box_correction",
    "force_error",
    "influence_function",
    "assign_charges",
    "assign_squared_charges",
    "kspace_energy",
    "a2_tilde",
    "field_ik",
    "field_fd",
    "interpolate_forces",
    "dft",
    "p3m_energy",
];

/// The sections of `p3m.json`, each a list of rows of floats.
fn lean_rows() -> Vec<(String, Vec<Vec<f64>>)> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../spec/vectors/p3m.json");
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
    let hex = field.strip_prefix("0x").expect("hex bits in p3m.json");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("64-bit hex in p3m.json"))
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
    fn sym(&mut self) -> Sym3 {
        [self.f(), self.f(), self.f(), self.f(), self.f(), self.f()]
    }
    fn m3(&mut self) -> M3 {
        [self.v3(), self.v3(), self.v3()]
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
    /// An influence case: the inputs and the box volume.
    fn case(&mut self) -> (InfluenceInputs, f64) {
        let grid = self.grid();
        let order = self.n() as u8;
        let nfdord = self.n() as u8;
        let nalias = self.i();
        let shape = ChargeShape::from_nshape(self.f() as i32).expect("nshape");
        let width = self.f();
        let (rect, edges) = self.box_();
        let (volume, recip) = reciprocal_box(&edges, rect);
        (
            InfluenceInputs {
                order,
                nfdord,
                nalias,
                shape,
                width,
                mesh: Mesh::new(&edges, grid),
                recip,
            },
            volume,
        )
    }
}

fn sym_vec(s: &Sym3) -> Vec<f64> {
    s.to_vec()
}

#[test]
fn assignment_1d_matches_lean_bit_for_bit() {
    let rows = section("assignment_1d");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let rust = assignment_1d(row[0] as u8, row[1]);
        compare("assignment_1d", i, &[rust], &row[2..], &mut mismatches);
    }
    report("assignment_1d", rows.len(), &mismatches);
}

#[test]
fn weights_and_nearest_point_match_lean_bit_for_bit() {
    let rows = section("weights");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let (p, u) = (row[0] as u8, row[1]);
        let mut rust = vec![nearest_grid_point(p, u) as f64];
        rust.extend(assignment_weights(p, u));
        // the nearest point as an integer (Lean's truncation of -0.0 is 0)
        let mut lean = row[2..].to_vec();
        lean[0] = lean[0] as i64 as f64;
        compare("weights", i, &rust, &lean, &mut mismatches);
        // the weights of a stencil add up to one away from the half-integers of NGP
        let sum: f64 = assignment_weights(p, u).iter().sum();
        let half = (u - u.floor() - 0.5).abs() == 0.0;
        // (an even order truncates toward zero, so a negative u loses weight: the kernels see
        // in-box positions only)
        if u >= 0.0 && !(p == 1 && half) {
            assert!((sum - 1.0).abs() < 1e-14, "p {p} u {u}: Σw = {sum}");
        }
    }
    report("weights", rows.len(), &mismatches);
}

#[test]
fn selfterm_fp_and_fd_coefficients_and_sinc_power_match_lean_bit_for_bit() {
    let mut mismatches = Vec::new();
    let rows = section("selfterm_fp");
    for (i, row) in rows.iter().enumerate() {
        // gromosXX's f_p is even in n: check both signs
        for sign in [1, -1] {
            let rust = selfterm_fp(row[0] as u8, sign * row[1] as i64, row[2]);
            compare("selfterm_fp", i, &[rust], &row[3..], &mut mismatches);
        }
    }
    report("selfterm_fp", rows.len(), &mismatches);
    let rows = section("fd_coefficients");
    for (i, row) in rows.iter().enumerate() {
        compare(
            "fd_coefficients",
            i,
            &fd_coefficients(row[0] as u8),
            &row[1..],
            &mut mismatches,
        );
    }
    report("fd_coefficients", rows.len(), &mismatches);
    let rows = section("sinc_power");
    for (i, row) in rows.iter().enumerate() {
        let rust = sinc_power(row[1], row[0] as u8);
        compare("sinc_power", i, &[rust], &row[2..], &mut mismatches);
    }
    report("sinc_power", rows.len(), &mismatches);
}

#[test]
fn mesh_matches_lean_bit_for_bit() {
    let rows = section("mesh");
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let (_, edges) = c.box_();
        let grid = c.grid();
        let m = Mesh::new(&edges, grid);
        let mut rust: Vec<f64> = m.h.iter().flatten().copied().collect();
        rust.extend(m.h_inv.iter().flatten());
        rust.push(m.cell_volume);
        compare("mesh", i, &rust, c.rest(), &mut mismatches);
    }
    report("mesh", rows.len(), &mismatches);
}

#[test]
fn influence_function_matches_lean_bit_for_bit() {
    let mut mismatches = Vec::new();
    let rows = section("influence_at");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let (inp, _) = c.case();
        let l = [c.i(), c.i(), c.i()];
        let wg = c.f() == 1.0;
        let pt = influence_at(&inp, l, wg);
        let mut rust = vec![pt.ghat];
        rust.extend(pt.gammahat);
        rust.push(pt.q);
        compare("influence_at", i, &rust, c.rest(), &mut mismatches);
    }
    report("influence_at", rows.len(), &mismatches);

    let rows = section("quality_at");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let (inp, _) = c.case();
        let l = [c.i(), c.i(), c.i()];
        let g = c.f();
        let (a, b) = quality_at(&inp, l, g);
        compare("quality_at", i, &[a, b], c.rest(), &mut mismatches);
    }
    report("quality_at", rows.len(), &mismatches);

    let rows = section("influence_function");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let (inp, volume) = c.case();
        let n = c.n();
        let inf = influence_function(&inp, volume, true);
        assert_eq!(inf.ghat.len(), n);
        let mut rust = inf.ghat.clone();
        rust.extend(inf.gammahat.iter().flat_map(sym_vec));
        rust.push(inf.q);
        compare("influence_function", i, &rust, c.rest(), &mut mismatches);
    }
    report("influence_function", rows.len(), &mismatches);

    let rows = section("box_correction");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let g = c.f();
        let gam = c.sym();
        let t0 = c.m3();
        let dl = c.m3();
        let rust = influence_box_correction(g, &gam, &t0, &dl);
        compare("box_correction", i, &[rust], c.rest(), &mut mismatches);
    }
    report("box_correction", rows.len(), &mismatches);

    let rows = section("force_error");
    for (i, row) in rows.iter().enumerate() {
        let rust = force_error_estimate(row[0], row[1], row[2], row[3] as usize, row[4]);
        compare("force_error", i, &[rust], &row[5..], &mut mismatches);
    }
    report("force_error", rows.len(), &mismatches);
}

#[test]
fn assignment_matches_lean_bit_for_bit() {
    for (name, squared) in [("assign_charges", false), ("assign_squared_charges", true)] {
        let rows = section(name);
        let mut mismatches = Vec::new();
        for (i, row) in rows.iter().enumerate() {
            let mut c = Cursor::new(row);
            let grid = c.grid();
            let p = c.n() as u8;
            let (_, edges) = c.box_();
            let (r, q) = c.atoms();
            let m = Mesh::new(&edges, grid);
            let rust = if squared {
                assign_squared_charges(&m, p, &r, &q)
            } else {
                assign_charges(&m, p, &r, &q)
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
}

#[test]
fn grid_kernels_match_lean_bit_for_bit() {
    let mut mismatches = Vec::new();
    let rows = section("kspace_energy");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let n = c.n();
        let wg = c.f() == 1.0;
        let vc = c.f();
        let e0 = c.f();
        let (mut g, mut gam, mut re, mut im) = (vec![], vec![], vec![], vec![]);
        for _ in 0..n {
            g.push(c.f());
            gam.push(c.sym());
            re.push(c.f());
            im.push(c.f());
        }
        let (e, vir) = kspace_energy(&g, wg.then_some(&gam[..]), &re, &im, vc, e0);
        potential_hat(&g, &mut re, &mut im, e0);
        let mut rust = vec![e];
        rust.extend(vir);
        for k in 0..n {
            rust.extend([re[k], im[k]]);
        }
        compare("kspace_energy", i, &rust, c.rest(), &mut mismatches);
    }
    report("kspace_energy", rows.len(), &mismatches);

    let rows = section("a2_tilde");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let n = c.n();
        let wg = c.f() == 1.0;
        let s2 = c.f();
        let (mut g, mut gam, mut rq) = (vec![], vec![], vec![]);
        for _ in 0..n {
            g.push(c.f());
            gam.push(c.sym());
            rq.push(c.f());
        }
        let (a2t, d) = a2_tilde_p3m(&g, wg.then_some(&gam[..]), &rq, s2);
        let mut rust = vec![a2t];
        rust.extend(d);
        compare("a2_tilde", i, &rust, c.rest(), &mut mismatches);
    }
    report("a2_tilde", rows.len(), &mismatches);

    let rows = section("field_ik");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let grid = c.grid();
        let (rect, edges) = c.box_();
        let (_, recip) = reciprocal_box(&edges, rect);
        let m = Mesh::new(&edges, grid);
        let (mut re, mut im) = (vec![], vec![]);
        for _ in 0..m.len() {
            re.push(c.f());
            im.push(c.f());
        }
        let e = field_ik(&m, &recip, &re, &im);
        let mut rust = vec![];
        for k in 0..m.len() {
            rust.extend([
                e[0].0[k], e[1].0[k], e[2].0[k], e[0].1[k], e[1].1[k], e[2].1[k],
            ]);
        }
        compare("field_ik", i, &rust, c.rest(), &mut mismatches);
    }
    report("field_ik", rows.len(), &mismatches);

    let rows = section("field_fd");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let grid = c.grid();
        let q = c.n() as u8;
        let (_, edges) = c.box_();
        let m = Mesh::new(&edges, grid);
        let phi = c.take(m.len());
        let e = field_fd(&m, q, &phi);
        let mut rust = vec![];
        for k in 0..m.len() {
            rust.extend([e[0][k], e[1][k], e[2][k]]);
        }
        compare("field_fd", i, &rust, c.rest(), &mut mismatches);
    }
    report("field_fd", rows.len(), &mismatches);

    let rows = section("interpolate_forces");
    for (i, row) in rows.iter().enumerate() {
        let mut c = Cursor::new(row);
        let grid = c.grid();
        let p = c.n() as u8;
        let (_, edges) = c.box_();
        let (r, q) = c.atoms();
        let m = Mesh::new(&edges, grid);
        let ex = c.take(m.len());
        let ey = c.take(m.len());
        let ez = c.take(m.len());
        let f = interpolate_forces(&m, p, &r, &q, [&ex, &ey, &ez]);
        let rust: Vec<f64> = f.iter().flat_map(|v| [v.x, v.y, v.z]).collect();
        compare("interpolate_forces", i, &rust, c.rest(), &mut mismatches);
    }
    report("interpolate_forces", rows.len(), &mismatches);
}

/// `max |a − b| / max |b|`.
fn rel_error(a: &[f64], b: &[f64]) -> f64 {
    let scale = b.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let diff = a
        .iter()
        .zip(b)
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    diff / scale
}

#[test]
fn fft3_matches_the_direct_transform_within_1e_12() {
    let rows = section("dft");
    let mut worst = 0.0f64;
    for row in &rows {
        let mut c = Cursor::new(row);
        let grid = c.grid();
        let seed = c.f();
        let sign = c.f();
        let n = grid[0] * grid[1] * grid[2];
        // the input of `dftInput`: C sin and cos, as Lean's `Float.sin`, `Float.cos`
        let mut re: Vec<f64> = (0..n).map(|i| (0.37 * i as f64 + seed).sin()).collect();
        let mut im: Vec<f64> = (0..n)
            .map(|i| (0.91 * i as f64 + 2.0 * seed).cos())
            .collect();
        let mut fft = Fft3::new(grid);
        if sign < 0.0 {
            fft.forward(&mut re, &mut im);
        } else {
            fft.backward(&mut re, &mut im);
        }
        let lean = c.rest();
        let mut rust = re.clone();
        rust.extend(&im);
        let err = rel_error(&rust, lean);
        println!("fft3 {grid:?} sign {sign}: max relative difference {err:.2e}");
        assert!(err < 1e-12, "fft3 {grid:?}: {err:e}");
        worst = worst.max(err);
    }
    println!("fft3: largest relative difference from the direct transform {worst:.2e}");
}

#[test]
fn fft3_backward_of_forward_is_n_times_the_grid() {
    let grid = [8, 4, 16];
    let n = 8 * 4 * 16;
    let re0: Vec<f64> = (0..n)
        .map(|i| ((i * 7919) % 101) as f64 / 50.0 - 1.0)
        .collect();
    let im0: Vec<f64> = (0..n)
        .map(|i| ((i * 104729) % 97) as f64 / 48.0 - 1.0)
        .collect();
    let (mut re, mut im) = (re0.clone(), im0.clone());
    let mut fft = Fft3::new(grid);
    fft.forward(&mut re, &mut im);
    fft.backward(&mut re, &mut im);
    for i in 0..n {
        assert!((re[i] / n as f64 - re0[i]).abs() < 1e-14);
        assert!((im[i] / n as f64 - im0[i]).abs() < 1e-14);
    }
}

#[test]
fn p3m_kspace_energy_matches_the_lean_pipeline_within_1e_11() {
    let rows = section("p3m_energy");
    let mut worst = 0.0f64;
    for row in &rows {
        let mut c = Cursor::new(row);
        let (inp, volume) = c.case();
        let e0 = c.f();
        let (r, q) = c.atoms();
        let lean = c.f();
        let inf = influence_function(&inp, volume, false);
        let mut re = assign_charges(&inp.mesh, inp.order, &r, &q);
        let mut im = vec![0.0; re.len()];
        Fft3::new(inp.mesh.grid).forward(&mut re, &mut im);
        let (e, _) = kspace_energy(&inf.ghat, None, &re, &im, inp.mesh.cell_volume, e0);
        let err = ((e - lean) / lean).abs();
        worst = worst.max(err);
        assert!(err < 1e-11, "P3M energy {e} vs Lean {lean}: {err:e}");
    }
    println!("p3m k-space energy: largest relative difference from Lean {worst:.2e}");
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

/// The P3M k-space energy, forces, rms force-error estimate, and the `Ã2` part of the A term
/// (`P3m::evaluate` without the other constant terms).
fn p3m_kspace(
    inp: &InfluenceInputs,
    inf: &newen_forces::nonbonded::p3m::Influence,
    volume: f64,
    r: &[Vec3],
    q: &[f64],
) -> (f64, Vec<Vec3>, f64, f64) {
    let eps0_inv = FPE * 4.0 * std::f64::consts::PI;
    let sum_q2: f64 = q.iter().map(|x| x * x).sum();
    let error = force_error_estimate(inf.q, sum_q2, volume, r.len(), FPE);
    let mesh = inp.mesh;
    let mut fft = Fft3::new(mesh.grid);
    let mut re = assign_charges(&mesh, inp.order, r, q);
    let mut im = vec![0.0; re.len()];
    fft.forward(&mut re, &mut im);
    let (e, _) = kspace_energy(&inf.ghat, None, &re, &im, mesh.cell_volume, eps0_inv);
    // the A term's Ã2 part, `−Ã2 S̃² / (8π ε₀)`, with the P3M Ã2 of this configuration
    let mut rq = assign_squared_charges(&mesh, inp.order, r, q);
    let mut rq_im = vec![0.0; rq.len()];
    fft.forward(&mut rq, &mut rq_im);
    let (a2t, _) = a2_tilde_p3m(&inf.ghat, None, &rq, sum_q2);
    let e_self = -a2t * sum_q2 * eps0_inv / (8.0 * std::f64::consts::PI);
    potential_hat(&inf.ghat, &mut re, &mut im, eps0_inv);
    let forces = if inp.nfdord == 0 {
        let mut f = field_ik(&mesh, &inp.recip, &re, &im);
        for (a, b) in f.iter_mut() {
            fft.backward(a, b);
        }
        interpolate_forces(&mesh, inp.order, r, q, [&f[0].0, &f[1].0, &f[2].0])
    } else {
        fft.backward(&mut re, &mut im);
        let f = field_fd(&mesh, inp.nfdord, &re);
        interpolate_forces(&mesh, inp.order, r, q, [&f[0], &f[1], &f[2]])
    };
    (e, forces, error, e_self)
}

fn inputs(
    edges: &[[f64; 3]; 3],
    rect: bool,
    grid: usize,
    order: u8,
    nfdord: u8,
) -> (InfluenceInputs, f64) {
    let (volume, recip) = reciprocal_box(edges, rect);
    (
        InfluenceInputs {
            order,
            nfdord,
            nalias: 2,
            shape: ChargeShape::Gaussian,
            width: 0.4,
            mesh: Mesh::new(edges, [grid; 3]),
            recip,
        },
        volume,
    )
}

fn rms(f: &[Vec3]) -> f64 {
    (f.iter().map(|v| v.length_squared()).sum::<f64>() / f.len() as f64).sqrt()
}

#[test]
fn p3m_agrees_with_the_direct_ewald_sum_within_its_error_estimate() {
    let eps0_inv = FPE * 4.0 * std::f64::consts::PI;
    let cubic = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
    let rect = [[1.9, 0.0, 0.0], [0.0, 2.3, 0.0], [0.0, 0.0, 2.1]];
    let tric = [[2.0, 0.0, 0.0], [0.5, 1.9, 0.0], [-0.4, 0.3, 2.2]];
    for (edges, is_rect) in [(cubic, true), (rect, true), (tric, false)] {
        let (r, q) = system(&edges, 16, 7);
        let (volume, recip) = reciprocal_box(&edges, is_rect);
        let kv = k_vectors([14, 14, 14], 1e9, &recip, ChargeShape::Gaussian, 0.4);
        let ew = ewald_kspace(&r, &q, &kv, volume, eps0_inv);
        let sum_q2: f64 = q.iter().map(|x| x * x).sum();
        for (grid, order, nfdord) in [(32, 3, 0), (32, 5, 2), (16, 4, 0), (32, 2, 4)] {
            let (inp, vol) = inputs(&edges, is_rect, grid, order, nfdord);
            let inf = influence_function(&inp, vol, false);
            let (e, f, est, e_self) = p3m_kspace(&inp, &inf, vol, &r, &q);
            let df: Vec<Vec3> = f.iter().zip(&ew.forces).map(|(a, b)| *a - *b).collect();
            let rms_df = rms(&df);
            // each method's k-space energy with its own Ã2 part of the A term
            let ew_total =
                ew.energy - ew.a2_tilde * sum_q2 * eps0_inv / (8.0 * std::f64::consts::PI);
            // relative to the k-space energy itself, which the A term nearly cancels here
            let de = ((e + e_self - ew_total) / ew.energy).abs();
            println!(
                "box {:?} grid {grid}³ p {order} q {nfdord}: E_γ {e:.6} vs Ewald {:.6}; \
                 E_γ + E_Ã2 {:.8} vs Ewald {ew_total:.8} (diff / E_γ {de:.1e}); rms force difference \
                 {rms_df:.2e}, estimate {est:.2e}, rms Ewald force {:.2e}",
                edges[1],
                ew.energy,
                e + e_self,
                rms(&ew.forces)
            );
            assert!(rms_df < 3.0 * est, "rms {rms_df:e} vs estimate {est:e}");
            // cloud-in-cell (p = 2) is the least accurate assignment
            let bound = if order <= 2 { 2e-3 } else { 1e-4 };
            assert!(de < bound, "energy {} vs Ewald {ew_total}", e + e_self);
        }
    }
}

#[test]
fn p3m_forces_are_minus_the_gradient_of_its_energy_to_grid_accuracy() {
    let edges = [[2.0, 0.0, 0.0], [0.4, 2.1, 0.0], [0.0, 0.0, 1.9]];
    let (r, q) = system(&edges, 12, 3);
    for (grid, order, nfdord) in [(32, 3, 0), (32, 5, 2), (32, 4, 3)] {
        let rect = false;
        let (inp, vol) = inputs(&edges, rect, grid, order, nfdord);
        let inf = influence_function(&inp, vol, false);
        let (_, f, est, _) = p3m_kspace(&inp, &inf, vol, &r, &q);
        let h = 1e-5;
        let mut diff = Vec::new();
        for i in 0..r.len() {
            let mut g = [0.0; 3];
            for (a, ga) in g.iter_mut().enumerate() {
                let mut rp = r.clone();
                let mut rm = r.clone();
                rp[i][a] += h;
                rm[i][a] -= h;
                // the energy with the Ã2 part of the A term, which removes each charge's
                // interaction with itself through the grid
                let (ep, _, _, sp) = p3m_kspace(&inp, &inf, vol, &rp, &q);
                let (em, _, _, sm) = p3m_kspace(&inp, &inf, vol, &rm, &q);
                let (ep, em) = (ep + sp, em + sm);
                *ga = -(ep - em) / (2.0 * h);
            }
            diff.push(f[i] - Vec3::new(g[0], g[1], g[2]));
        }
        let rms_d = rms(&diff);
        println!(
            "grid {grid}³ p {order} q {nfdord}: rms (F − (−∇E)) {rms_d:.2e}, estimate {est:.2e}, \
             rms F {:.2e}",
            rms(&f)
        );
        // The P3M energy is not translation invariant at the grid scale: its error δE varies
        // over a cell h, so its gradient misses the forces by about δE/h, more than the rms
        // force error of the forces themselves (which the Ewald test bounds).
        assert!(
            rms_d < 5e-3 * rms(&f),
            "rms {rms_d:e} vs rms F {:e}",
            rms(&f)
        );
    }
}
