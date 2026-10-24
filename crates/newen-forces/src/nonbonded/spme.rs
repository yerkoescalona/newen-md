//! Smooth particle-mesh Ewald (NONBONDED NLRELE = 5, newen only, experimental): the
//! reciprocal-space part of the lattice sum on a grid with cardinal B-splines.
//!
//! Essmann, Perera, Berkowitz, Darden, Lee and Pedersen, J. Chem. Phys. 103, 8577 (1995): the
//! charges are spread onto the grid with the cardinal B-spline `M_n` of even order `n` = NSPORD
//! (eq. 4.1), the grid is Fourier transformed, multiplied by `Ĝ(k) = γ̂(a|k|)/k² · B(m)` with the
//! Euler exponential-spline factors `B = Π_α |b_α(m_α)|²` (eqs. 4.4, 4.7, 4.8), and transformed
//! back once; the forces are the derivatives of the splines against that convolution (eq. 4.9),
//! so they are the exact gradient of the mesh energy. With the Gaussian shape, `γ̂(κ) = e^{−κ²/4}`
//! with `κ = a|k|` is Essmann's `exp(−π²m²/β²)` for `β = 1/a`; with the eleven polynomial shapes
//! of GROMOS it generalises SPME to GROMOS's charge shaping.
//!
//! gromosXX has no SPME (its reader stops NLRELE at 4, the shifted reaction field, and reads
//! NSPORD into a field nothing uses), so there is no gromosXX oracle: the oracles are the Lean
//! Ewald model and P3M (`tests/spme_spec.rs`). The grid, the transforms (PhastFT, powers of two),
//! the real-space loops, the energy formula, the `Ã2` formula and the self, A and surface terms
//! are P3M's (`p3m.rs`, `ewald.rs`).

use super::ewald::{
    add_constant_terms, self_constants_of_box, LatticeSumParams, LsEnergies, SelfConstants,
};
use super::lattice_sum::{
    charge_sums, lattice_vector, reciprocal_box, shape_fourier, sym_tensor, ChargeShape, Sym3,
};
use super::p3m::{a2_tilde_p3m, kspace_energy, product, wave_index, Fft3, Mesh};
use newen_core::math::Vec3;

/// The smallest and largest spline order (NSPORD) newen accepts; the order must be even.
pub const MIN_ORDER: usize = 4;
pub const MAX_ORDER: usize = 12;

/// The cardinal B-spline `M_n(u)` by its recursion (Essmann eq. 4.1 and the line after it):
///
/// ```text
/// M_1(u) = 1 for 0 ≤ u < 1, else 0
/// M_k(u) = (u · M_{k−1}(u) + (k − u) · M_{k−1}(u − 1)) / (k − 1)
/// ```
///
/// `2^(n−1)` evaluations of `M_1`: for tests and the spec; the grid kernels use
/// [`spline_weights`].
///
/// spec: spme.bspline
pub fn bspline(n: usize, u: f64) -> f64 {
    assert!(n >= 1, "B-spline order 0");
    if n == 1 {
        return if (0.0..1.0).contains(&u) { 1.0 } else { 0.0 };
    }
    let k = n as f64;
    let a = bspline(n - 1, u);
    let b = bspline(n - 1, u - 1.0);
    (u * a + (k - u) * b) / (k - 1.0)
}

/// `(M_n(u), M_n'(u))` with `M_n'(u) = M_{n−1}(u) − M_{n−1}(u − 1)` (Essmann after eq. 4.1),
/// `n ≥ 2`.
///
/// spec: spme.bspline_and_deriv
pub fn bspline_and_deriv(n: usize, u: f64) -> (f64, f64) {
    assert!(n >= 2, "B-spline derivative of order {n}");
    (bspline(n, u), bspline(n - 1, u) - bspline(n - 1, u - 1.0))
}

/// The weights `θ_j = M_n(w + j)` and derivatives `dθ_j = M_{n−1}(w + j) − M_{n−1}(w + j − 1)`,
/// `j = 0 … n − 1`, a charge at grid coordinate `u` gives the grid points `⌊u⌋ − j`, with
/// `w = u − ⌊u⌋`: the recursion of [`bspline`] on the whole table at once (`n ≥ 2`),
///
/// ```text
/// t = (1, 0, …, 0)                                         (M_1(w + j))
/// per k = 2 … n:  (at k = n first dθ_j = t_j − t_{j−1}, t_{−1} = 0)
///   per j = k − 1 … 0:  x = w + j,  t_j = (x · t_j + (k − x) · t_{j−1}) / (k − 1)
/// ```
///
/// spec: spme.spline_weights
pub fn spline_weights(n: usize, u: f64) -> (Vec<f64>, Vec<f64>) {
    assert!(n >= 2, "B-spline order {n}");
    let w = u - u.floor();
    let mut t = vec![0.0; n];
    t[0] = 1.0;
    let mut d = vec![0.0; n];
    for k in 2..=n {
        if k == n {
            d[0] = t[0];
            for j in 1..n {
                d[j] = t[j] - t[j - 1];
            }
        }
        let kf = k as f64;
        let div = kf - 1.0;
        for j in (0..k).rev() {
            let x = w + j as f64;
            let below = if j > 0 { t[j - 1] } else { 0.0 };
            t[j] = (x * t[j] + (kf - x) * below) / div;
        }
    }
    (t, d)
}

/// The Euler exponential-spline factor `|b(m)|²` of an axis of `K` grid points for the wave
/// index `m` (Essmann eqs. 4.4, 4.7): `|e^{2πi(n−1)m/K}| = 1`, so
///
/// ```text
/// per j = 0 … n − 2:  θ = M_n(j + 1) (bspline),  x = 2π · m · j / K,
///   re += θ cos x,   im += θ sin x
/// |b(m)|² = 1 / (re · re + im · im)
/// ```
///
/// The denominator never vanishes for even `n` (it does at `2m = K` for odd `n`).
///
/// spec: spme.euler_spline_factor
pub fn euler_spline_factor(n: usize, k_alpha: usize, m: i64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut re = 0.0;
    let mut im = 0.0;
    for j in 0..n - 1 {
        let theta = bspline(n, (j + 1) as f64);
        let x = two_pi * m as f64 * j as f64 / k_alpha as f64;
        re += theta * x.cos();
        im += theta * x.sin();
    }
    1.0 / (re * re + im * im)
}

/// `Ĝ` and `Γ̂` at the wave vector `k` for the Euler factor `b = B(m)`:
///
/// ```text
/// k² = k_x k_x + k_y k_y + k_z k_z,  k2i = 1/k²,  ak = a · sqrt(k²),  (γ̂, γ̂') = shape_fourier(ak)
/// Ĝ = k2i · γ̂ · b
/// Γ̂ = s k ⊗ k with s = (ak γ̂' − 2γ̂) · (k2i · k2i) · b
/// ```
///
/// `Γ̂` is the derivative of `Ĝ` the k-space virial uses (P3M's, here exact: `B` depends on the
/// integer `m` only); `Ĝ + Γ̂` is the Ewald sum's `γ̂/k² 1 + (akγ̂' − 2γ̂)/k⁴ k⊗k` times `B`.
fn influence_from(k: [f64; 3], b: f64, shape: ChargeShape, width: f64) -> (f64, Sym3) {
    let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
    let k2i = 1.0 / k2;
    let ak = width * k2.sqrt();
    let (fc, fc_deriv) = shape_fourier(shape, ak);
    let ghat = k2i * fc * b;
    let s = (ak * fc_deriv - 2.0 * fc) * (k2i * k2i) * b;
    (ghat, sym_tensor([s * k[0], s * k[1], s * k[2]], k))
}

/// The SPME influence function `Ĝ(k_l) = γ̂(a|k_l|)/k_l² · B(l)` and its derivative tensor `Γ̂` at
/// the integer triple `l ≠ 0` of a grid of `grid` points (Essmann eqs. 4.7, 4.8 with GROMOS's
/// shape `γ̂`): `k_l = l_x b₁ + l_y b₂ + l_z b₃`, `B = |b_x(l_x)|² · |b_y(l_y)|² · |b_z(l_z)|²`
/// ([`euler_spline_factor`]). No alias sums, no differentiation operator.
///
/// spec: spme.influence_at
pub fn spme_influence_at(
    l: [i64; 3],
    grid: [usize; 3],
    n: usize,
    shape: ChargeShape,
    width: f64,
    reciprocal: &[[f64; 3]; 3],
) -> (f64, Sym3) {
    let k = lattice_vector(reciprocal, l);
    let b = euler_spline_factor(n, grid[0], l[0])
        * euler_spline_factor(n, grid[1], l[1])
        * euler_spline_factor(n, grid[2], l[2]);
    influence_from(k, b, shape, width)
}

/// The Euler factors of every storage position of each axis, `|b_α(wave_index(i))|²`.
pub fn euler_factors(grid: [usize; 3], n: usize) -> [Vec<f64>; 3] {
    grid.map(|k| {
        (0..k)
            .map(|i| euler_spline_factor(n, k, wave_index(i, k)))
            .collect()
    })
}

/// `Ĝ` and `Γ̂` on the whole grid (storage order of [`Mesh::index`], zero at `l = 0`), each point
/// as [`spme_influence_at`] with the Euler factors `euler` of [`euler_factors`].
pub fn spme_influence_function(
    grid: [usize; 3],
    shape: ChargeShape,
    width: f64,
    reciprocal: &[[f64; 3]; 3],
    euler: &[Vec<f64>; 3],
) -> (Vec<f64>, Vec<Sym3>) {
    let [nx, ny, nz] = grid;
    let mut ghat = vec![0.0; nx * ny * nz];
    let mut gammahat = vec![[0.0; 6]; nx * ny * nz];
    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                let l = [wave_index(x, nx), wave_index(y, ny), wave_index(z, nz)];
                if l == [0, 0, 0] {
                    continue;
                }
                let k = lattice_vector(reciprocal, l);
                let b = euler[0][x] * euler[1][y] * euler[2][z];
                let i = z + nz * (y + ny * x);
                (ghat[i], gammahat[i]) = influence_from(k, b, shape, width);
            }
        }
    }
    (ghat, gammahat)
}

/// The spline weights and derivatives of a position along the three grid axes, and the grid
/// point `⌊u⌋` the weight `j` is counted back from.
fn atom_splines(mesh: &Mesh, n: usize, r: Vec3) -> ([i64; 3], [(Vec<f64>, Vec<f64>); 3]) {
    let u = mesh.grid_coordinates(r);
    (u.map(|c| c.floor() as i64), u.map(|c| spline_weights(n, c)))
}

/// The charge density on the grid (Essmann eq. 4.6, divided by the cell volume as P3M's
/// [`assign_charges`](super::p3m::assign_charges)):
///
/// ```text
/// per atom i (in order):  u = H⁻ᵀ r_i,  θ_α = spline_weights(n, u_α)
///   per j_x, j_y, j_z = 0 … n − 1 (x slowest):
///     ρ(⌊u⌋ − j) += (θ_x[j_x] · θ_y[j_y] · θ_z[j_z] / V_cell) · q_i
/// ```
///
/// spec: spme.assign
pub fn spme_assign(positions: &[Vec3], charges: &[f64], mesh: &Mesh, n: usize) -> Vec<f64> {
    let mut rho = vec![0.0; mesh.len()];
    for (r, &q) in positions.iter().zip(charges) {
        let (base, [(wx, _), (wy, _), (wz, _)]) = atom_splines(mesh, n, *r);
        for jx in 0..n {
            let px = base[0] - jx as i64;
            for jy in 0..n {
                let py = base[1] - jy as i64;
                for jz in 0..n {
                    let pz = base[2] - jz as i64;
                    let w = wx[jx] * wy[jy] * wz[jz] / mesh.cell_volume;
                    rho[mesh.index([px, py, pz])] += w * q;
                }
            }
        }
    }
    rho
}

/// The overlap `f(d) = Σ_j θ_j θ_{j−d}` of a charge's weights along one axis with themselves
/// shifted by `d = −(n − 1) … n − 1` grid points (`j` ascending, the terms with `j − d` outside
/// `0 … n − 1` left out).
fn overlaps(theta: &[f64]) -> Vec<f64> {
    let n = theta.len() as i64;
    (-(n - 1)..n)
        .map(|d| {
            let mut s = 0.0;
            for j in 0..n {
                let i = j - d;
                if (0..n).contains(&i) {
                    s += theta[j as usize] * theta[i as usize];
                }
            }
            s
        })
        .collect()
}

/// The squared charges on the grid for the SPME `Ã2`, as P3M's
/// [`assign_squared_charges`](super::p3m::assign_squared_charges): the autocorrelation of each
/// charge's weights, around the grid origin,
///
/// ```text
/// per atom i (in order):  θ_α = spline_weights(n, u_α),  f_α = overlaps(θ_α)
///   per d_x, d_y, d_z in −(n − 1) … n − 1 (x slowest):
///     R(d) += (f_x(d_x) · f_y(d_y) · f_z(d_z) · (1/V_cell)) · q_i²
/// ```
///
/// so that `Re R̂ = Σ_i q_i² |Ŵ_i|²/V_cell` and [`a2_tilde_p3m`] gives each charge's interaction
/// with itself through the grid.
///
/// spec: spme.assign_squared_charges
pub fn spme_assign_squared_charges(
    positions: &[Vec3],
    charges: &[f64],
    mesh: &Mesh,
    n: usize,
) -> Vec<f64> {
    let mut grid = vec![0.0; mesh.len()];
    let cell_volume_i = 1.0 / mesh.cell_volume;
    let bound = n as i64 - 1;
    for (r, &q) in positions.iter().zip(charges) {
        let q2 = q * q;
        let (_, [(wx, _), (wy, _), (wz, _)]) = atom_splines(mesh, n, *r);
        let (fx, fy, fz) = (overlaps(&wx), overlaps(&wy), overlaps(&wz));
        for dx in -bound..=bound {
            let w_x = fx[(dx + bound) as usize];
            for dy in -bound..=bound {
                let w_y = fy[(dy + bound) as usize];
                for dz in -bound..=bound {
                    let w_z = fz[(dz + bound) as usize];
                    let w = w_x * w_y * w_z * cell_volume_i;
                    grid[mesh.index([dx, dy, dz])] += w * q2;
                }
            }
        }
    }
    grid
}

/// The forces from the splines' derivatives (Essmann eq. 4.9), for the convolution
/// `c = Re backward(Ĝ ρ̂)` (unnormalised, as [`Fft3::backward`]) on the grid:
///
/// ```text
/// per atom i:  u = H⁻ᵀ r_i,  (θ_α, dθ_α) = spline_weights(n, u_α),  g = 0
///   per j_x, j_y, j_z (x slowest):  c = c(⌊u⌋ − j)
///     g_x += dθ_x θ_y θ_z c;   g_y += θ_x dθ_y θ_z c;   g_z += θ_x θ_y dθ_z c
///   F_i = (−(q_i · eps0_inv) / N) · H⁻¹ g
/// ```
///
/// `g` is the gradient of `Σ_grid W_i c` in grid coordinates and `∂u_β/∂r_α = (H⁻¹)_αβ`. With
/// `E = ½ V_cell eps0_inv Σ_l Ĝ_l |ρ̂_l|²/N` ([`kspace_energy`]) and `ρ = Σ_i q_i W_i / V_cell`,
/// `F_i = −∂E/∂r_i` exactly, for the symmetric part of `Ĝ`, which the real part keeps.
///
/// spec: spme.forces
pub fn spme_forces(
    positions: &[Vec3],
    charges: &[f64],
    mesh: &Mesh,
    n: usize,
    convolution: &[f64],
    eps0_inv: f64,
) -> Vec<Vec3> {
    let n_points = mesh.len() as f64;
    positions
        .iter()
        .zip(charges)
        .map(|(r, &q)| {
            let (base, [(wx, dx), (wy, dy), (wz, dz)]) = atom_splines(mesh, n, *r);
            let mut g = [0.0; 3];
            for jx in 0..n {
                let px = base[0] - jx as i64;
                for jy in 0..n {
                    let py = base[1] - jy as i64;
                    for jz in 0..n {
                        let pz = base[2] - jz as i64;
                        let c = convolution[mesh.index([px, py, pz])];
                        g[0] += dx[jx] * wy[jy] * wz[jz] * c;
                        g[1] += wx[jx] * dy[jy] * wz[jz] * c;
                        g[2] += wx[jx] * wy[jy] * dz[jz] * c;
                    }
                }
            }
            let grad = product(&mesh.h_inv, g);
            let s = -(q * eps0_inv) / n_points;
            Vec3::new(s * grad[0], s * grad[1], s * grad[2])
        })
        .collect()
}

/// `Ĝ ρ̂` in place of `ρ̂`, the transform the backward FFT turns into the convolution of
/// [`spme_forces`].
fn multiply(ghat: &[f64], re: &mut [f64], im: &mut [f64]) {
    for i in 0..ghat.len() {
        re[i] *= ghat[i];
        im[i] *= ghat[i];
    }
}

/// The SPME settings of a run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpmeParams {
    /// NGX, NGY, NGZ (powers of two)
    pub grid: [usize; 3],
    /// NSPORD: even, 4 ..= 12
    pub order: usize,
}

/// The influence function of a box.
#[derive(Debug, Clone)]
struct InfluenceState {
    edges: [[f64; 3]; 3],
    ghat: Vec<f64>,
    gammahat: Vec<Sym3>,
}

/// The reciprocal-space, self, A and surface terms of SPME for one run: P3M's pipeline
/// ([`P3m::evaluate`](super::p3m::P3m::evaluate)) with B-spline assignment, the SPME influence
/// function, recomputed exactly whenever the box changes (no quality check, no first-order
/// correction), and the forces from the splines' derivatives after one backward transform.
#[derive(Debug, Clone)]
pub struct Spme {
    pub params: LatticeSumParams,
    pub mesh_params: SpmeParams,
    fft: Fft3,
    euler: [Vec<f64>; 3],
    influence: Option<InfluenceState>,
    constants: Option<SelfConstants>,
}

impl Spme {
    /// Panics unless every grid dimension is a power of two and the order is even in
    /// `4 ..= 12` (the run's checks refuse others).
    pub fn new(params: LatticeSumParams, mesh_params: SpmeParams) -> Self {
        let n = mesh_params.order;
        assert!(
            n.is_multiple_of(2) && (MIN_ORDER..=MAX_ORDER).contains(&n),
            "SPME order {n}: NSPORD must be even, 4 ..= 12"
        );
        Self {
            params,
            fft: Fft3::new(mesh_params.grid),
            euler: euler_factors(mesh_params.grid, n),
            mesh_params,
            influence: None,
            constants: None,
        }
    }

    /// Add the reciprocal-space, self, A and surface terms to `storage` (the virial as `−2 Ξ`)
    /// and return their energies, as [`P3m::evaluate`](super::p3m::P3m::evaluate); `Ã2`
    /// (NA2CLC 3, and 4 read as 3) from the squared-charge grid every step.
    ///
    /// spec: spme.kspace
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &mut self,
        positions: &[Vec3],
        charges: &[f64],
        edges: [[f64; 3]; 3],
        rectangular: bool,
        periodicity: &newen_core::math::Periodicity,
        storage: &mut super::ForceStorage,
        virial: bool,
    ) -> Result<LsEnergies, String> {
        let p = self.params;
        let n = self.mesh_params.order;
        if self.constants.as_ref().is_none_or(|c| c.edges != edges) {
            self.constants = Some(self_constants_of_box(&p, edges, rectangular)?);
        }
        let constants = self.constants.clone().expect("set above");
        let mesh = Mesh::new(&edges, self.mesh_params.grid);
        if self.influence.as_ref().is_none_or(|s| s.edges != edges) {
            let (_, recip) = reciprocal_box(&edges, rectangular);
            let (ghat, gammahat) =
                spme_influence_function(mesh.grid, p.shape, p.width, &recip, &self.euler);
            self.influence = Some(InfluenceState {
                edges,
                ghat,
                gammahat,
            });
        }
        let inf = self.influence.as_ref().expect("set above");
        let gammahat = virial.then_some(inf.gammahat.as_slice());
        let (_, sum_q2) = charge_sums(charges);

        // the positions put into the box, as P3M
        let [a, b, c] = edges;
        let o = Vec3::new(
            (a[0] + b[0] + c[0]) / 2.0,
            (a[1] + b[1] + c[1]) / 2.0,
            (a[2] + b[2] + c[2]) / 2.0,
        );
        let in_box: Vec<Vec3> = positions
            .iter()
            .map(|&r| periodicity.nearest_image(r, o) + o)
            .collect();
        let len = mesh.len();

        // Ã2 from the squared-charge grid (NA2CLC 3, and 4 read as 3)
        let mut a2_tilde = (0.0, [0.0; 6]);
        if matches!(p.na2clc, 3 | 4) {
            let mut rq_re = spme_assign_squared_charges(&in_box, charges, &mesh, n);
            let mut rq_im = vec![0.0; len];
            self.fft.forward(&mut rq_re, &mut rq_im);
            a2_tilde = a2_tilde_p3m(&inf.ghat, gammahat, &rq_re, sum_q2);
        }

        // the charge density, its transform, the energy, one backward transform, the forces
        let mut re = spme_assign(&in_box, charges, &mesh, n);
        let mut im = vec![0.0; len];
        self.fft.forward(&mut re, &mut im);
        let (energy, kspace_virial) =
            kspace_energy(&inf.ghat, gammahat, &re, &im, mesh.cell_volume, p.eps0_inv);
        multiply(&inf.ghat, &mut re, &mut im);
        self.fft.backward(&mut re, &mut im);
        let forces = spme_forces(&in_box, charges, &mesh, n, &re, p.eps0_inv);
        for (f, df) in storage.forces.iter_mut().zip(&forces) {
            *f += *df;
        }
        if virial {
            super::ewald::add_minus_two(storage, &kspace_virial);
        }
        let mut e = add_constant_terms(
            &p, &constants, positions, charges, edges, a2_tilde, storage, virial,
        );
        e.kspace = energy;
        Ok(e)
    }
}
