//! P3M (NONBONDED NLRELE = 3): the reciprocal-space part of the lattice sum on a grid.
//!
//! GROMOS Vol. 2 (`assign_p3m`, `optim_influ`, `force_from_field_p3m`): the charges are spread
//! onto a grid of `N_x × N_y × N_z` points with an assignment function `w_p` of order `p`
//! (NASORD), the grid charge density is Fourier transformed, multiplied by the optimal influence
//! function `Ĝ` and transformed back to a potential, the field on the grid comes from `ik` or a
//! finite-difference operator of order NFDORD, and is interpolated back to the atoms with the same
//! `w_p`. The grid cell is `H` with rows `a/N_x`, `b/N_y`, `c/N_z`; a position `r` has grid
//! coordinates `u = H⁻ᵀ r` (the dual vectors of the rows of `H`).
//!
//! The per-point kernels follow gromosXX (`interaction/nonbonded/interaction/latticesum.{h,cc}`,
//! `configuration/influence_function.cc`) in their loop and operation order and are checked bit
//! for bit against their Lean model (`spec/NewenSpec/P3M.lean`). The FFT is PhastFT, gromosXX's
//! is FFTW: the two agree to rounding, not to the bit, so everything after the first transform is
//! checked within a tolerance (`tests/p3m_spec.rs`). PhastFT transforms powers of two only, so
//! every grid dimension must be one (gromosXX accepts any even size).

use super::ewald::{
    add_constant_terms, self_constants_of_box, LatticeSumParams, LsEnergies, SelfConstants,
};
use super::lattice_sum::{
    charge_sums, lattice_vector, reciprocal_box, shape_fourier, sym_tensor, ChargeShape, Sym3,
};
use newen_core::math::Vec3;
use phastft::planner::{Direction, PlannerDit64};

/// C `pow` from the platform's libm, the function gromosXX's `fourier_charge_assignment_1d` calls
/// (and Lean's `Float.pow`).
mod libm {
    extern "C" {
        fn pow(x: f64, y: f64) -> f64;
    }

    /// `pow(x, y)`.
    #[inline]
    pub fn c_pow(x: f64, y: f64) -> f64 {
        // SAFETY: a pure C math function of two doubles.
        unsafe { pow(x, y) }
    }
}

/// A 3 × 3 matrix, `m[row][column]`.
pub type M3 = [[f64; 3]; 3];

/// gromosXX `math::epsilon`: below this `|D̂|²` the influence function is set to zero.
pub const EPSILON: f64 = 0.000000000001;

/// The largest assignment order (NASORD) and finite-difference order (NFDORD) gromosXX has.
pub const MAX_ORDER: u8 = 5;

/// The one-dimensional charge-assignment function `w_p(ξ)` of order `p` = NASORD, gromosXX
/// `Lattice_Sum::charge_assignment_1d` (`latticesum.h:508-565`) in its operation order, with
/// `|ξ|` as `a`:
///
/// ```text
/// p = 1:  1 for a < 1/2                                  (nearest grid point)
/// p = 2:  1 − a for a < 1                                (cloud in cell)
/// p = 3:  0.75 − ξ·ξ for a < 1/2;  0.125 t² with t = 2a − 3 for a < 3/2
/// p = 4:  1/6 (3a³ − 6a² + 4) for a < 1;  −1/6 t³ with t = a − 2 for a < 2
/// p = 5:  1/192 (48ξ⁴ − 120ξ² + 115);  −1/96 (16ξ⁴ − 80a³ + 120ξ² − 20a − 55);
///         1/384 t⁴ with t = 2a − 5           (a < 1/2, 3/2, 5/2)
/// zero outside
/// ```
///
/// The B-spline of order `p` (Vol. 2 `assign_p3m`): `Σ_n w_p(x − n) = 1`. Every branch reads `|ξ|`
/// or `ξ²`, so `w_p(−ξ)` has the bits of `w_p(ξ)`.
///
/// spec: p3m.assignment_1d
pub fn assignment_1d(order: u8, xi: f64) -> f64 {
    let absf = xi.abs();
    match order {
        1 => {
            if absf < 0.5 {
                1.0
            } else {
                0.0
            }
        },
        2 => {
            if absf < 1.0 {
                1.0 - absf
            } else {
                0.0
            }
        },
        3 => {
            if absf < 0.5 {
                0.75 - xi * xi
            } else if absf < 1.5 {
                let term = 2.0 * absf - 3.0;
                0.125 * term * term
            } else {
                0.0
            }
        },
        4 => {
            if absf < 1.0 {
                let absf2 = absf * absf;
                1.0 / 6.0 * (3.0 * absf2 * absf - 6.0 * absf2 + 4.0)
            } else if absf < 2.0 {
                let term = absf - 2.0;
                -1.0 / 6.0 * term * term * term
            } else {
                0.0
            }
        },
        5 => {
            if absf < 0.5 {
                let xi2 = xi * xi;
                let xi4 = xi2 * xi2;
                1.0 / 192.0 * (48.0 * xi4 - 120.0 * xi2 + 115.0)
            } else if absf < 1.5 {
                let xi2 = xi * xi;
                let xi4 = xi2 * xi2;
                let absf3 = absf * absf * absf;
                -1.0 / 96.0 * (16.0 * xi4 - 80.0 * absf3 + 120.0 * xi2 - 20.0 * absf - 55.0)
            } else if absf < 2.5 {
                let tmp = 2.0 * absf - 5.0;
                let tmp2 = tmp * tmp;
                1.0 / 384.0 * tmp2 * tmp2
            } else {
                0.0
            }
        },
        _ => panic!("assignment order {order}: only NASORD 1..=5 exist"),
    }
}

/// The stencil of an assignment of order `p`: grid offsets `−(p − 1)/2 ..= p/2` (integer
/// division) around the nearest grid point.
#[inline]
pub fn stencil(order: u8) -> (i64, i64) {
    let p = i64::from(order);
    (-(p - 1) / 2, p / 2)
}

/// The grid point a stencil of order `p` is centred on, for the grid coordinate `u`: the nearest
/// integer (`rint`, halves to even) for odd `p`, `u` truncated toward zero (C `int(u)`) for even
/// `p` (`latticesum.cc:96-106`).
///
/// spec: p3m.nearest_grid_point
#[inline]
pub fn nearest_grid_point(order: u8, u: f64) -> i64 {
    if order % 2 == 1 {
        u.round_ties_even() as i64
    } else {
        u.trunc() as i64
    }
}

/// `w_p(u − (n₀ + d))` for the stencil offsets `d = −(p − 1)/2 ..= p/2` in order, `n₀` the
/// [`nearest_grid_point`]; the unused entries are zero. gromosXX computes `u − point` when it
/// assigns and `point − u` when it interpolates (`latticesum.cc:113, 747`); the weights have the
/// same bits (see [`assignment_1d`]).
///
/// spec: p3m.assignment_weights
pub fn assignment_weights(order: u8, u: f64) -> [f64; 5] {
    let n0 = nearest_grid_point(order, u);
    let (lower, upper) = stencil(order);
    let mut w = [0.0; 5];
    for d in lower..=upper {
        w[(d - lower) as usize] = assignment_1d(order, u - (n0 + d) as f64);
    }
    w
}

/// `det m` as gromosXX's `math::det` (`gmath.h:1197-1205`) sums it.
pub fn det3(m: &M3) -> f64 {
    m[0][1] * m[1][2] * m[2][0] - m[0][2] * m[1][1] * m[2][0] + m[0][2] * m[1][0] * m[2][1]
        - m[0][0] * m[1][2] * m[2][1]
        - m[0][1] * m[1][0] * m[2][2]
        + m[0][0] * m[1][1] * m[2][2]
}

/// `m⁻¹` as gromosXX's `math::inverse` (`gmath.h:1211-1223`): the adjugate times `1/det m`.
pub fn inverse3(m: &M3) -> M3 {
    let d = 1.0 / det3(m);
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d,
        ],
    ]
}

/// `mᵀ`.
pub fn transpose3(m: &M3) -> M3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// `Σ_i m[i][j] v_i` for each `j`: gromosXX's `math::product(Matrix, Vec)` (`gmath.h:925-932`),
/// which multiplies by the transpose.
#[inline]
pub fn product_t(m: &M3, v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
        m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
        m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2],
    ]
}

/// `m v`: `math::product` of `mᵀ`, as gromosXX writes `product(H_trans, k)` and
/// `product(H_inv_trans, D)`.
#[inline]
pub fn product(m: &M3, v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// The grid of a P3M run in a box: `H` (rows `a/N_x`, `b/N_y`, `c/N_z`, gromosXX
/// `math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz)`), `H⁻¹`, and the cell volume
/// `det Hᵀ` (`latticesum.cc:57-64`).
///
/// spec: p3m.mesh
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mesh {
    pub grid: [usize; 3],
    pub h: M3,
    pub h_inv: M3,
    pub cell_volume: f64,
}

impl Mesh {
    pub fn new(edges: &[[f64; 3]; 3], grid: [usize; 3]) -> Self {
        let mut h = [[0.0; 3]; 3];
        for i in 0..3 {
            let n = grid[i] as f64;
            for j in 0..3 {
                h[i][j] = edges[i][j] / n;
            }
        }
        Self {
            grid,
            h,
            h_inv: inverse3(&h),
            cell_volume: det3(&transpose3(&h)),
        }
    }

    /// The number of grid points.
    pub fn len(&self) -> usize {
        self.grid[0] * self.grid[1] * self.grid[2]
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The grid coordinates `H⁻ᵀ r` of a position, `math::product(H_inv, r)`.
    #[inline]
    pub fn grid_coordinates(&self, r: Vec3) -> [f64; 3] {
        product_t(&self.h_inv, [r.x, r.y, r.z])
    }

    /// The storage index of the grid point `p`, periodic (gromosXX `Mesh::operator()`,
    /// `mesh.h:151-158`): `z + N_z (y + N_y x)`, `x` slowest.
    #[inline]
    pub fn index(&self, p: [i64; 3]) -> usize {
        let [nx, ny, nz] = self.grid.map(|n| n as i64);
        let x = p[0].rem_euclid(nx);
        let y = p[1].rem_euclid(ny);
        let z = p[2].rem_euclid(nz);
        (z + nz * (y + ny * x)) as usize
    }
}

/// The wave-vector index `l` of the storage position `i` along an axis of `n` points: `i` for
/// `i ≤ n/2`, `i − n` above (`influence_function.cc:134-139`), so that `l ∈ (−n/2, n/2]`.
#[inline]
pub fn wave_index(i: usize, n: usize) -> i64 {
    if i > n / 2 {
        i as i64 - n as i64
    } else {
        i as i64
    }
}

/// The charge density on the grid, gromosXX `Lattice_Sum::calculate_charge_density`
/// (`latticesum.cc:46-139`) in its loop and operation order, for positions in the box:
///
/// ```text
/// per atom i (in order):  u = H⁻ᵀ r_i,  n₀ = nearest_grid_point(p, u)
///   per d_x, d_y, d_z in the stencil (x slowest):
///     ρ(n₀ + d) += (w_p(u_x − n_x) · w_p(u_y − n_y) · w_p(u_z − n_z) / V_cell) · q_i
/// ```
///
/// GROMOS Vol. 2: `assign_p3m`.
///
/// spec: p3m.assign_charges
pub fn assign_charges(mesh: &Mesh, order: u8, positions: &[Vec3], charges: &[f64]) -> Vec<f64> {
    let mut rho = vec![0.0; mesh.len()];
    let (lower, upper) = stencil(order);
    for (r, &q) in positions.iter().zip(charges) {
        let u = mesh.grid_coordinates(*r);
        let n0 = u.map(|c| nearest_grid_point(order, c));
        for dx in lower..=upper {
            let px = n0[0] + dx;
            let w_x = assignment_1d(order, u[0] - px as f64);
            for dy in lower..=upper {
                let py = n0[1] + dy;
                let w_y = assignment_1d(order, u[1] - py as f64);
                for dz in lower..=upper {
                    let pz = n0[2] + dz;
                    let w_z = assignment_1d(order, u[2] - pz as f64);
                    let w = w_x * w_y * w_z / mesh.cell_volume;
                    rho[mesh.index([px, py, pz])] += w * q;
                }
            }
        }
    }
    rho
}

/// gromosXX's `f_p(n; s)` for the P3M `Ã2` (`Lattice_Sum::p3m_selfterm_fp`,
/// `latticesum.h:672-769`): the overlap `Σ_m w_p(s − m) w_p(s − m − n)` of an assignment of order
/// `p` with itself shifted by `n` grid points, for a charge at offset `s` from its stencil centre,
/// a polynomial in `s²`; zero is never returned for `|n| ≤ p − 1`.
///
/// spec: p3m.selfterm_fp
pub fn selfterm_fp(order: u8, n: i64, s: f64) -> f64 {
    let s2 = s * s;
    match (order, n.abs()) {
        (1, _) => 1.0,
        (2, 0) => 1.0 / 2.0 * (1.0 + 4.0 * s2),
        (2, 1) => 1.0 / 4.0 * (1.0 - 4.0 * s2),
        (3, 0) => 1.0 / 32.0 * (19.0 + s2 * (-24.0 + s2 * 48.0)),
        (3, 1) => 1.0 / 16.0 * (3.0 + s2 * (8.0 - s2 * 16.0)),
        (3, 2) => {
            let term = 1.0 - 4.0 * s2;
            1.0 / 64.0 * term * term
        },
        (4, 0) => 1.0 / 576.0 * (265.0 + s2 * (204.0 + s2 * (-528.0 + s2 * 320.0))),
        (4, 1) => 1.0 / 2304.0 * (575.0 + s2 * (-564.0 + s2 * (1488.0 + s2 * -960.0))),
        (4, 2) => 1.0 / 1152.0 * (23.0 + s2 * (84.0 + s2 * (-240.0 + s2 * 192.0))),
        (4, 3) => {
            let term = 1.0 - 4.0 * s2;
            1.0 / 2304.0 * term * term * term
        },
        (5, 0) => {
            1.0 / 73728.0
                * (32227.0 + s2 * (-9520.0 + s2 * (28960.0 + s2 * (-29440.0 + s2 * 8960.0))))
        },
        (5, 1) => {
            1.0 / 18432.0 * (4389.0 + s2 * (1792.0 + s2 * (-5472.0 + s2 * (5632.0 + s2 * -1792.0))))
        },
        (5, 2) => {
            1.0 / 36864.0 * (1559.0 + s2 * (-1456.0 + s2 * (4512.0 + s2 * (-4864.0 + s2 * 1792.0))))
        },
        (5, 3) => {
            1.0 / 18432.0 * (19.0 + s2 * (128.0 + s2 * (-416.0 + s2 * (512.0 + s2 * -256.0))))
        },
        (5, 4) => {
            let term = 1.0 - 4.0 * s2;
            let term2 = term * term;
            1.0 / 147456.0 * term2 * term2
        },
        _ => panic!("selfterm_fp: order {order}, |n| = {}", n.abs()),
    }
}

/// The squared charges on the grid for the P3M `Ã2`, gromosXX
/// `Lattice_Sum::calculate_squared_charge_grid` (`latticesum.cc:170-258`) in its loop and
/// operation order: around the grid origin, not around the atom,
///
/// ```text
/// per atom i (in order):  u = H⁻ᵀ r_i,  n₀ = nearest_grid_point(p, u)
///   s = u − n₀ (odd p),  s = u − n₀ − 1/2 (even p)
///   per n_x, n_y, n_z in −(p − 1) ..= p − 1 (x slowest):
///     R(n) += (f_p(n_x; s_x) · f_p(n_y; s_y) · f_p(n_z; s_z) · (1/V_cell)) · q_i²
/// ```
///
/// GROMOS Vol. 2: `u_ga_slf_def_p3m`, `f_g_func`.
///
/// spec: p3m.assign_squared_charges
pub fn assign_squared_charges(
    mesh: &Mesh,
    order: u8,
    positions: &[Vec3],
    charges: &[f64],
) -> Vec<f64> {
    let mut grid = vec![0.0; mesh.len()];
    let cell_volume_i = 1.0 / mesh.cell_volume;
    let bound = i64::from(order) - 1;
    let odd = order % 2 == 1;
    for (r, &q) in positions.iter().zip(charges) {
        let q2 = q * q;
        let u = mesh.grid_coordinates(*r);
        let n0 = u.map(|c| nearest_grid_point(order, c));
        let s: [f64; 3] = std::array::from_fn(|k| {
            if odd {
                u[k] - n0[k] as f64
            } else {
                u[k] - n0[k] as f64 - 0.5
            }
        });
        for dx in -bound..=bound {
            let w_x = selfterm_fp(order, dx, s[0]);
            for dy in -bound..=bound {
                let w_y = selfterm_fp(order, dy, s[1]);
                for dz in -bound..=bound {
                    let w_z = selfterm_fp(order, dz, s[2]);
                    let w = w_x * w_y * w_z * cell_volume_i;
                    grid[mesh.index([dx, dy, dz])] += w * q2;
                }
            }
        }
    }
    grid
}

/// The coefficients `c_1 … c_q` of the finite-difference operator of order `q` = NFDORD,
/// gromosXX `Lattice_Sum::finite_difference_constants` (`latticesum.h:598-633`), unused entries
/// zero; order 0 (`ik` differentiation) has none. GROMOS Vol. 2 table `tab:fd_coef`:
/// `Σ_j c_j (f(x + j) − f(x − j)) / (2j)` is `f'(x)` for polynomials of degree `≤ 2q`.
///
/// spec: p3m.fd_coefficients
pub fn fd_coefficients(nfdord: u8) -> [f64; 5] {
    match nfdord {
        0 => [0.0; 5],
        1 => [1.0, 0.0, 0.0, 0.0, 0.0],
        2 => [4.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0],
        3 => [3.0 / 2.0, -3.0 / 5.0, 0.1, 0.0, 0.0],
        4 => [8.0 / 5.0, -4.0 / 5.0, 8.0 / 35.0, -1.0 / 35.0, 0.0],
        5 => [
            5.0 / 3.0,
            -20.0 / 21.0,
            5.0 / 14.0,
            -5.0 / 63.0,
            1.0 / 126.0,
        ],
        _ => panic!("finite-difference order {nfdord}: only NFDORD 0..=5 exist"),
    }
}

/// `P̂_p(x) = (sin(x/2) / (x/2))^p`, the Fourier transform of `w_p` at `x`, gromosXX
/// `Lattice_Sum::fourier_charge_assignment_1d` (`latticesum.h:583-589`):
/// `1` at `x = 0`, else `pow(2/x · sin(x/2), p)` with C `pow`.
///
/// spec: p3m.sinc_power
pub fn sinc_power(x: f64, p: u8) -> f64 {
    if x == 0.0 {
        1.0
    } else {
        libm::c_pow(2.0 / x * (x / 2.0).sin(), f64::from(p))
    }
}

/// What the influence function of a grid needs: the method's orders and the box, through the
/// grid matrix and the reciprocal vectors of [`reciprocal_box`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InfluenceInputs {
    /// NASORD
    pub order: u8,
    /// NFDORD (0: `ik` differentiation)
    pub nfdord: u8,
    /// NALIAS
    pub nalias: i64,
    pub shape: ChargeShape,
    /// ASHAPE
    pub width: f64,
    pub mesh: Mesh,
    /// `2π/V (b × c, c × a, a × b)`
    pub recip: M3,
}

/// The sums over the alias vectors `m` of one wave vector (`influence_function.cc:173-233`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AliasSums {
    /// `Σ_m k_m γ̂(a k_m)/k_m² P̂²(k_m)`
    pub numerator: [f64; 3],
    /// `Σ_m P̂²(k_m)` (`1` for p = 1)
    pub denominator: f64,
    /// `Σ_m γ̂(a k_m)²/k_m²`
    pub k2i_fourier2: f64,
    /// The numerator of `Γ̂` (zero unless asked for).
    pub gamma_numerator: Sym3,
}

/// `D̂(k)`, the Fourier transform of the differentiation operator, and `|D̂|²`
/// (`influence_function.cc:156-171`): `k` itself for `ik` (NFDORD 0), else
/// `H⁻¹ (Σ_j c_j sin(j x_α) / j)_α` with `x = H k`. GROMOS Vol. 2: `finite_diff`, `exact_diff`.
pub fn d_hat(inp: &InfluenceInputs, k_l: [f64; 3]) -> ([f64; 3], f64) {
    let d = if inp.nfdord == 0 {
        k_l
    } else {
        let c = fd_coefficients(inp.nfdord);
        let th = product(&inp.mesh.h, k_l);
        let mut d = [0.0; 3];
        for j in 1..=inp.nfdord {
            let jf = f64::from(j);
            let f = c[usize::from(j) - 1] / jf;
            for a in 0..3 {
                d[a] += (jf * th[a]).sin() * f;
            }
        }
        product(&inp.mesh.h_inv, d)
    };
    (d, d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
}

/// The sums over the alias vectors `k_m = k_l + 2π H⁻¹ m`, `m ∈ [−NALIAS, NALIAS]³` (`m_x`
/// slowest), gromosXX `Influence_Function::calculate` (`influence_function.cc:173-224`) in its
/// operation order: per `m`
///
/// ```text
/// k_m = k_l + 2π · H⁻¹ m,  k² = |k_m|²,  k2i = 1/k²,  ak = a · sqrt(k²)
/// (γ̂, γ̂') = shape_fourier(ak),   x = H k_m,   P̂ = P̂_p(x_x) · P̂_p(x_y) · P̂_p(x_z)
/// num += k_m · ((k2i · γ̂) · P̂²);  den += P̂²;  Σγ̂²/k² += γ̂ · (k2i · γ̂)
/// with_gamma:  Γnum += (P̂² · k2i) · ((k⊗D + D⊗k − (2 k·D) (k2i k⊗k + |D|⁻² D⊗D)) γ̂
///                                   + (ak γ̂' k2i k·D) k⊗k)
/// ```
///
/// and `den = 1` for p = 1. GROMOS Vol. 2: `optim_influ`, `alias_def`, `big_gammahat_def`.
pub fn alias_sums(
    inp: &InfluenceInputs,
    k_l: [f64; 3],
    d: [f64; 3],
    abs2_d: f64,
    with_gamma: bool,
) -> AliasSums {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut num = [0.0; 3];
    let mut den = 0.0;
    let mut k2if2 = 0.0;
    let mut gnum: Sym3 = [0.0; 6];
    let na = inp.nalias;
    for m1 in -na..=na {
        for m2 in -na..=na {
            for m3 in -na..=na {
                let hm = product(&inp.mesh.h_inv, [m1 as f64, m2 as f64, m3 as f64]);
                let k = [
                    k_l[0] + two_pi * hm[0],
                    k_l[1] + two_pi * hm[1],
                    k_l[2] + two_pi * hm[2],
                ];
                let k_2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                let k_2i = 1.0 / k_2;
                let ak = inp.width * k_2.sqrt();
                let (fc, fc_deriv) = shape_fourier(inp.shape, ak);
                let th = product(&inp.mesh.h, k);
                let p_hat = sinc_power(th[0], inp.order)
                    * sinc_power(th[1], inp.order)
                    * sinc_power(th[2], inp.order);
                let p_hat_2 = p_hat * p_hat;
                let k2i_fourier = k_2i * fc;
                let s = k2i_fourier * p_hat_2;
                num[0] += k[0] * s;
                num[1] += k[1] * s;
                num[2] += k[2] * s;
                den += p_hat_2;
                k2if2 += fc * k2i_fourier;
                if with_gamma {
                    let k_dot_d = k[0] * d[0] + k[1] * d[1] + k[2] * d[2];
                    let k_x_d = sym_tensor(k, d);
                    let d_x_k = sym_tensor(d, k);
                    let k_x_k = sym_tensor(k, k);
                    let d_x_d = sym_tensor(d, d);
                    let d_2i = 1.0 / abs2_d;
                    let two_kd = 2.0 * k_dot_d;
                    let pre = p_hat_2 * k_2i;
                    let deriv = ak * fc_deriv * k_2i * k_dot_d;
                    for j in 0..6 {
                        let bracket = k_2i * k_x_k[j] + d_2i * d_x_d[j];
                        let t = (k_x_d[j] + d_x_k[j] - two_kd * bracket) * fc + deriv * k_x_k[j];
                        gnum[j] += pre * t;
                    }
                }
            }
        }
    }
    if inp.order == 1 {
        den = 1.0;
    }
    AliasSums {
        numerator: num,
        denominator: den,
        k2i_fourier2: k2if2,
        gamma_numerator: gnum,
    }
}

/// The influence function at one wave vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InfluencePoint {
    /// `Ĝ(k_l)`
    pub ghat: f64,
    /// `Γ̂(k_l)`, its derivative tensor (zero unless asked for)
    pub gammahat: Sym3,
    /// The contribution of `k_l` to `Q` (before `16π²/V`).
    pub q: f64,
}

/// The optimal influence function at the wave vector `k_l` of the integer triple `l ≠ 0`,
/// gromosXX `Influence_Function::calculate` (`influence_function.cc:140-287`) for one `l`:
///
/// ```text
/// k_l = l_x b₁ + l_y b₂ + l_z b₃;   (D̂, |D̂|²) = d_hat(k_l);   sums = alias_sums(k_l, D̂)
/// |D̂|² < 1e-12:  Ĝ = 0, Γ̂ = 0, q = Σγ̂²/k²
/// else:          num = D̂ · Σ k_m γ̂/k_m² P̂²,   den = |D̂|² (Σ P̂²)²
///                Ĝ = num / den,   q = Σγ̂²/k² − num · Ĝ,   Γ̂ = Γnum / den
/// ```
///
/// gromosXX builds `k_l` with the transpose of the reciprocal basis, as for the Ewald sum (the
/// same vector in a rectangular box). GROMOS Vol. 2: `optim_influ`, `big_gammahat_def`,
/// `q_calc_2`.
///
/// spec: p3m.influence_at
pub fn influence_at(inp: &InfluenceInputs, l: [i64; 3], with_gamma: bool) -> InfluencePoint {
    let k_l = lattice_vector(&inp.recip, l);
    let (d, abs2_d) = d_hat(inp, k_l);
    let sums = alias_sums(inp, k_l, d, abs2_d, with_gamma);
    if abs2_d < EPSILON {
        return InfluencePoint {
            ghat: 0.0,
            gammahat: [0.0; 6],
            q: sums.k2i_fourier2,
        };
    }
    let n = sums.numerator;
    let numerator = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
    let denominator = abs2_d * sums.denominator * sums.denominator;
    let ghat = numerator / denominator;
    InfluencePoint {
        ghat,
        gammahat: if with_gamma {
            sums.gamma_numerator.map(|x| x / denominator)
        } else {
            [0.0; 6]
        },
        q: sums.k2i_fourier2 - numerator * ghat,
    }
}

/// The contribution of `l ≠ 0` to `Q` for the influence-function value `ghat` in the current
/// box, gromosXX `Influence_Function::evaluate_quality` (`influence_function.cc:410-503`) for
/// one `l`: `Σγ̂²/k²` if finite differences give `|D̂|² < 1e-12`, else
/// `(Ĝ² |D̂|² (ΣP̂²)²) + (Σγ̂²/k² − 2Ĝ D̂ · Σ k_m γ̂/k_m² P̂²)`, the two added in turn.
///
/// spec: p3m.quality_at
pub fn quality_at(inp: &InfluenceInputs, l: [i64; 3], ghat: f64) -> (f64, f64) {
    let k_l = lattice_vector(&inp.recip, l);
    let (d, abs2_d) = d_hat(inp, k_l);
    let sums = alias_sums(inp, k_l, d, abs2_d, false);
    if inp.nfdord != 0 && abs2_d < EPSILON {
        return (sums.k2i_fourier2, 0.0);
    }
    let den = sums.denominator;
    let n = sums.numerator;
    let first = ghat * ghat * abs2_d * den * den;
    let second = sums.k2i_fourier2 - 2.0 * ghat * (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]);
    (first, second)
}

/// The influence function of a grid.
#[derive(Debug, Clone, PartialEq)]
pub struct Influence {
    /// `Ĝ` per grid point (storage order of [`Mesh::index`]; `Ĝ(0) = 0`)
    pub ghat: Vec<f64>,
    /// `Γ̂` per grid point (empty unless asked for)
    pub gammahat: Vec<Sym3>,
    /// `Q = 16π²/V Σ_l q_l`
    pub q: f64,
}

/// The wave indices of the grid in gromosXX's order (`influence_function.cc:132-147`), one list
/// per storage `x` plane: `l_y` from `1 − N_y/2` to `N_y/2`, `l_z` likewise inside it.
fn wave_planes(grid: [usize; 3]) -> impl Iterator<Item = Vec<[i64; 3]>> {
    let [nx, ny, nz] = grid;
    let (ny2, nz2) = ((ny / 2) as i64, (nz / 2) as i64);
    (0..nx).map(move |x| {
        let l0 = wave_index(x, nx);
        let mut plane = Vec::with_capacity(ny * nz);
        for l1 in 1 - ny2..=ny2 {
            for l2 in 1 - nz2..=nz2 {
                plane.push([l0, l1, l2]);
            }
        }
        plane
    })
}

/// The optimal influence function on the whole grid and `Q`, gromosXX
/// `Influence_Function::calculate` (`influence_function.cc:91-344`): [`influence_at`] for every
/// `l ≠ 0` in gromosXX's order, `q` summed per `x` plane and the planes in turn, then
/// `Q = q · (16π · π / V)`.
pub fn influence_function(inp: &InfluenceInputs, volume: f64, with_gamma: bool) -> Influence {
    let mesh = inp.mesh;
    let n = mesh.len();
    let mut ghat = vec![0.0; n];
    let mut gammahat = if with_gamma {
        vec![[0.0; 6]; n]
    } else {
        Vec::new()
    };
    let mut q = 0.0;
    for plane in wave_planes(mesh.grid) {
        let mut q_plane = 0.0;
        for l in plane {
            if l == [0, 0, 0] {
                continue;
            }
            let pt = influence_at(inp, l, with_gamma);
            let i = mesh.index(l);
            ghat[i] = pt.ghat;
            if with_gamma {
                gammahat[i] = pt.gammahat;
            }
            q_plane += pt.q;
        }
        q += q_plane;
    }
    q *= 16.0 * std::f64::consts::PI * std::f64::consts::PI / volume;
    Influence { ghat, gammahat, q }
}

/// `Q` of the influence function `ghat` in the box of `inp`, gromosXX
/// `Influence_Function::evaluate_quality` (`influence_function.cc:355-526`): [`quality_at`] for
/// every `l ≠ 0`, summed as [`influence_function`] sums.
pub fn influence_quality(inp: &InfluenceInputs, volume: f64, ghat: &[f64]) -> f64 {
    let mesh = inp.mesh;
    let mut q = 0.0;
    for plane in wave_planes(mesh.grid) {
        let mut q_plane = 0.0;
        for l in plane {
            if l == [0, 0, 0] {
                continue;
            }
            let (a, b) = quality_at(inp, l, ghat[mesh.index(l)]);
            q_plane += a;
            q_plane += b;
        }
        q += q_plane;
    }
    q * (16.0 * std::f64::consts::PI * std::f64::consts::PI / volume)
}

/// The estimated rms force error of the influence function, gromosXX
/// (`influence_function.cc:311`): `fpe · S̃² · sqrt(Q / (V · N_atoms))` with `fpe = 1/(4πε₀)`.
/// GROMOS Vol. 2: `deserno`.
///
/// spec: p3m.force_error
pub fn force_error_estimate(q: f64, sum_q2: f64, volume: f64, n_atoms: usize, fpe: f64) -> f64 {
    fpe * sum_q2 * (q / (volume * n_atoms as f64)).sqrt()
}

/// `Ĝ − tr(Γ̂ · ᵗL₀⁻¹ · (ᵗL − ᵗL₀))`, the influence function corrected to first order for a box
/// that differs from the one it was computed in, gromosXX `Influence_Function::operator()`
/// (`influence_function.h:140-148`), with `ᵗL` the matrix whose rows are the box edges and the
/// matrix products summed from zero in index order. GROMOS Vol. 2: `gammahat_upd`.
///
/// spec: p3m.influence_box_correction
pub fn influence_box_correction(ghat: f64, gammahat: &Sym3, tl0_inv: &M3, delta_tl: &M3) -> f64 {
    let g = super::lattice_sum::sym_to_matrix(gammahat);
    let mut p = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                p[i][j] += g[i][k] * tl0_inv[k][j];
            }
        }
    }
    let mut t = [0.0; 3];
    for (i, ti) in t.iter_mut().enumerate() {
        for k in 0..3 {
            *ti += p[i][k] * delta_tl[k][i];
        }
    }
    ghat - (t[0] + t[1] + t[2])
}

/// The k-space energy and virial of P3M, gromosXX `Lattice_Sum::calculate_potential_and_energy`
/// (`latticesum.cc:384-422`) in its loop and operation order, over the transformed charge density
/// `ρ̂` (unnormalised forward transform) in storage order:
///
/// ```text
/// |ρ̂|²/N = (Re ρ̂ · Re ρ̂ + Im ρ̂ · Im ρ̂) / N;   E += Ĝ · |ρ̂|²/N;   Ξ += (Γ̂ + Ĝ 1) · |ρ̂|²/N
/// E ·= ½ · V_cell · eps0_inv  (as ((½ E) V_cell) eps0_inv);   Ξ ·= (−¼ eps0_inv) V_cell
/// ```
///
/// with `N` the number of grid points; the `1/V` of the energy comes from `V_cell/N` and the
/// unnormalised transform. Without `gammahat` the virial is zero. GROMOS Vol. 2:
/// `e_gamma_p3m`, `reci_space_virial`.
///
/// spec: p3m.kspace_energy
pub fn kspace_energy(
    ghat: &[f64],
    gammahat: Option<&[Sym3]>,
    rho_re: &[f64],
    rho_im: &[f64],
    cell_volume: f64,
    eps0_inv: f64,
) -> (f64, Sym3) {
    let n = ghat.len() as f64;
    let mut energy = 0.0;
    let mut virial: Sym3 = [0.0; 6];
    for i in 0..ghat.len() {
        let abs2 = (rho_re[i] * rho_re[i] + rho_im[i] * rho_im[i]) / n;
        energy += ghat[i] * abs2;
        if let Some(g) = gammahat {
            let mut t = g[i];
            t[0] += ghat[i];
            t[3] += ghat[i];
            t[5] += ghat[i];
            for j in 0..6 {
                virial[j] += t[j] * abs2;
            }
        }
    }
    let s = -0.25 * eps0_inv * cell_volume;
    (0.5 * energy * cell_volume * eps0_inv, virial.map(|x| x * s))
}

/// `φ̂ = ((eps0_inv · Ĝ) · ρ̂) / √N`, the transformed potential (`latticesum.cc:392`), in place
/// of `ρ̂`.
///
/// spec: p3m.potential_hat
pub fn potential_hat(ghat: &[f64], re: &mut [f64], im: &mut [f64], eps0_inv: f64) {
    let sqrt_n = (ghat.len() as f64).sqrt();
    for i in 0..ghat.len() {
        let f = eps0_inv * ghat[i];
        re[i] = f * re[i] / sqrt_n;
        im[i] = f * im[i] / sqrt_n;
    }
}

/// The transformed field by `ik` differentiation, gromosXX `Lattice_Sum::calculate_electric_field`
/// with NFDORD = 0 (`latticesum.cc:544-580`): per grid point with wave vector
/// `k_l = l_x b₁ + l_y b₂ + l_z b₃` (`l` from [`wave_index`], the Nyquist planes `l = N/2`
/// included), `Ê_α = k_α · (−i φ̂ / √N)` with `−i φ̂ = (Im φ̂, −Re φ̂)`. Returns
/// `[(Re Ê_α, Im Ê_α); 3]`; three backward transforms give the field, of which the real part is
/// used. Taking the real part makes the Nyquist planes harmless: it keeps the Hermitian part of
/// `Ê`, which is `−i k' φ̂` with the Nyquist components of `l` set to zero.
///
/// spec: p3m.field_ik
pub fn field_ik(
    mesh: &Mesh,
    recip: &M3,
    phi_re: &[f64],
    phi_im: &[f64],
) -> [(Vec<f64>, Vec<f64>); 3] {
    let n = mesh.len();
    let sqrt_n = (n as f64).sqrt();
    let mut out: [(Vec<f64>, Vec<f64>); 3] = std::array::from_fn(|_| (vec![0.0; n], vec![0.0; n]));
    let [nx, ny, nz] = mesh.grid;
    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                let l = [wave_index(x, nx), wave_index(y, ny), wave_index(z, nz)];
                let k = lattice_vector(recip, l);
                let i = z + nz * (y + ny * x);
                let mre = phi_im[i] / sqrt_n;
                let mim = -phi_re[i] / sqrt_n;
                for a in 0..3 {
                    out[a].0[i] = k[a] * mre;
                    out[a].1[i] = k[a] * mim;
                }
            }
        }
    }
    out
}

/// The field on the grid by finite differences of order `q` = NFDORD from the potential `φ`
/// (real part of the backward transform of `φ̂`), after gromosXX
/// `Lattice_Sum::calculate_electric_field` (`latticesum.cc:581-640`): per grid point
///
/// ```text
/// g_α = Σ_{j=1..q} (c_j / √N) · (φ(n − j e_α) − φ(n + j e_α)) / (2 j)      (j in order)
/// E = H⁻¹ g
/// ```
///
/// gromosXX divides the difference along grid axis `α` by `2 j |h_α|` (the length of the
/// cell edge) and stores it as the Cartesian component `α`, which is the gradient only when the
/// cell edges are orthogonal; newen maps the per-axis derivatives through `H⁻¹` (`Σ_α g_α h*_α`
/// with `h*_α` the dual vectors), the operator whose transform is the `D̂` of the influence
/// function. In a rectangular box the two differ only in rounding. GROMOS Vol. 2: `finite_diff`.
///
/// spec: p3m.field_fd
pub fn field_fd(mesh: &Mesh, nfdord: u8, phi: &[f64]) -> [Vec<f64>; 3] {
    let n = mesh.len();
    let sqrt_n = (n as f64).sqrt();
    let c = fd_coefficients(nfdord);
    let mut out: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; n]);
    let [nx, ny, nz] = mesh.grid;
    for x in 0..nx as i64 {
        for y in 0..ny as i64 {
            for z in 0..nz as i64 {
                let mut g = [0.0; 3];
                for j in 1..=i64::from(nfdord) {
                    let two_j = 2.0 * j as f64;
                    let d = [
                        (phi[mesh.index([x - j, y, z])] - phi[mesh.index([x + j, y, z])]) / two_j,
                        (phi[mesh.index([x, y - j, z])] - phi[mesh.index([x, y + j, z])]) / two_j,
                        (phi[mesh.index([x, y, z - j])] - phi[mesh.index([x, y, z + j])]) / two_j,
                    ];
                    let f = c[(j - 1) as usize] / sqrt_n;
                    for a in 0..3 {
                        g[a] += f * d[a];
                    }
                }
                let e = product(&mesh.h_inv, g);
                let i = mesh.index([x, y, z]);
                for a in 0..3 {
                    out[a][i] = e[a];
                }
            }
        }
    }
    out
}

/// The force on every atom from the field on the grid, gromosXX `Lattice_Sum::calculate_force`
/// (`latticesum.cc:690-778`) in its loop and operation order, for positions in the box:
///
/// ```text
/// per atom i:  u = H⁻ᵀ r_i,  n₀ = nearest_grid_point(p, u),  E_i = 0
///   per d_x, d_y, d_z (x slowest):  E_i += (w_p(n_x − u_x) w_p(n_y − u_y) w_p(n_z − u_z)) E(n₀ + d)
///   F_i = q_i E_i
/// ```
///
/// GROMOS Vol. 2: `force_from_field_p3m`, `force_interpol`.
///
/// spec: p3m.interpolate_forces
pub fn interpolate_forces(
    mesh: &Mesh,
    order: u8,
    positions: &[Vec3],
    charges: &[f64],
    field: [&[f64]; 3],
) -> Vec<Vec3> {
    let (lower, upper) = stencil(order);
    positions
        .iter()
        .zip(charges)
        .map(|(r, &q)| {
            let u = mesh.grid_coordinates(*r);
            let n0 = u.map(|c| nearest_grid_point(order, c));
            let mut e = [0.0; 3];
            for dx in lower..=upper {
                let px = n0[0] + dx;
                let w_x = assignment_1d(order, px as f64 - u[0]);
                for dy in lower..=upper {
                    let py = n0[1] + dy;
                    let w_y = assignment_1d(order, py as f64 - u[1]);
                    for dz in lower..=upper {
                        let pz = n0[2] + dz;
                        let w_z = assignment_1d(order, pz as f64 - u[2]);
                        let w = w_x * w_y * w_z;
                        let i = mesh.index([px, py, pz]);
                        e[0] += w * field[0][i];
                        e[1] += w * field[1][i];
                        e[2] += w * field[2][i];
                    }
                }
            }
            Vec3::new(q * e[0], q * e[1], q * e[2])
        })
        .collect()
}

/// `Ã2` of P3M and its derivative, gromosXX `Lattice_Sum::calculate_p3m_selfterm`
/// (`latticesum.cc:440-489`) in its loop and operation order, over the real part of the
/// transformed squared-charge grid `R̂`:
///
/// ```text
/// S += Ĝ · Re R̂;   D += (Γ̂ + Ĝ 1) · Re R̂            (grid points in storage order)
/// Ã2 = 4π · S / N / S̃²;   dÃ2 = ((4π / N) / S̃²) · D
/// ```
///
/// GROMOS Vol. 2: `u_ga_slf_def_p3m`.
///
/// spec: p3m.a2_tilde
pub fn a2_tilde_p3m(
    ghat: &[f64],
    gammahat: Option<&[Sym3]>,
    rq_re: &[f64],
    sum_q2: f64,
) -> (f64, Sym3) {
    let n = ghat.len() as f64;
    let mut sum = 0.0;
    let mut deriv: Sym3 = [0.0; 6];
    for i in 0..ghat.len() {
        sum += ghat[i] * rq_re[i];
        if let Some(g) = gammahat {
            let mut t = g[i];
            t[0] += ghat[i];
            t[3] += ghat[i];
            t[5] += ghat[i];
            for j in 0..6 {
                deriv[j] += t[j] * rq_re[i];
            }
        }
    }
    let four_pi = 4.0 * std::f64::consts::PI;
    let s = four_pi / n / sum_q2;
    (four_pi * sum / n / sum_q2, deriv.map(|x| s * x))
}

/// The three-dimensional discrete Fourier transform of a grid stored as `z + N_z (y + N_y x)`,
/// by one-dimensional PhastFT transforms along `z`, `y` and `x`, unnormalised in both directions
/// as FFTW's (`configuration/mesh.cc:62-68, 142-147`): forward `Σ_n f_n e^{−2πi k·n/N}`, backward
/// `Σ_k f̂_k e^{+2πi k·n/N}`. Every dimension must be a power of two.
pub struct Fft3 {
    grid: [usize; 3],
    planners: [PlannerDit64; 3],
    line_re: Vec<f64>,
    line_im: Vec<f64>,
}

impl std::fmt::Debug for Fft3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Fft3").field("grid", &self.grid).finish()
    }
}

impl Clone for Fft3 {
    fn clone(&self) -> Self {
        Self::new(self.grid)
    }
}

impl Fft3 {
    /// Panics unless every dimension is a power of two.
    pub fn new(grid: [usize; 3]) -> Self {
        assert!(
            grid.iter().all(|n| n.is_power_of_two()),
            "P3M grid {grid:?}: PhastFT transforms powers of two only"
        );
        let longest = grid.iter().copied().max().unwrap_or(1);
        Self {
            grid,
            planners: grid.map(PlannerDit64::new),
            line_re: vec![0.0; longest],
            line_im: vec![0.0; longest],
        }
    }

    /// Forward transform in place.
    ///
    /// spec: p3m.fft3
    pub fn forward(&mut self, re: &mut [f64], im: &mut [f64]) {
        self.transform(re, im, Direction::Forward);
    }

    /// Backward (inverse, unnormalised) transform in place.
    pub fn backward(&mut self, re: &mut [f64], im: &mut [f64]) {
        self.transform(re, im, Direction::Inverse);
    }

    fn transform(&mut self, re: &mut [f64], im: &mut [f64], direction: Direction) {
        let [nx, ny, nz] = self.grid;
        assert_eq!(re.len(), nx * ny * nz);
        assert_eq!(im.len(), nx * ny * nz);
        // stride between consecutive points along each axis
        let strides = [ny * nz, nz, 1];
        for axis in [2, 1, 0] {
            let n = self.grid[axis];
            let stride = strides[axis];
            // PhastFT's inverse divides by n; FFTW's backward transform does not. n is a power
            // of two, so multiplying by it again is exact.
            let rescale = if direction == Direction::Inverse {
                n as f64
            } else {
                1.0
            };
            let (lr, li) = (&mut self.line_re[..n], &mut self.line_im[..n]);
            for start in line_starts(self.grid, axis) {
                for t in 0..n {
                    lr[t] = re[start + t * stride];
                    li[t] = im[start + t * stride];
                }
                phastft::fft_f64_dit_with_planner(lr, li, direction, &self.planners[axis]);
                for t in 0..n {
                    re[start + t * stride] = lr[t] * rescale;
                    im[start + t * stride] = li[t] * rescale;
                }
            }
        }
    }
}

/// The storage index of the first point of every grid line along `axis`.
fn line_starts(grid: [usize; 3], axis: usize) -> Vec<usize> {
    let [nx, ny, nz] = grid;
    let mut out = Vec::with_capacity(nx * ny * nz / grid[axis]);
    match axis {
        2 => {
            for x in 0..nx {
                for y in 0..ny {
                    out.push(nz * (y + ny * x));
                }
            }
        },
        1 => {
            for x in 0..nx {
                for z in 0..nz {
                    out.push(z + nz * ny * x);
                }
            }
        },
        _ => {
            for y in 0..ny {
                for z in 0..nz {
                    out.push(z + nz * y);
                }
            }
        },
    }
    out
}

/// The P3M settings of a run (NONBONDED block, GROMOS Vol. 4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct P3mParams {
    /// NGX, NGY, NGZ (powers of two)
    pub grid: [usize; 3],
    /// NASORD, 1..=5
    pub order: u8,
    /// NFDORD, 0..=5 (0: `ik` differentiation)
    pub nfdord: u8,
    /// NALIAS
    pub nalias: i64,
    /// NQEVAL: evaluate the force error every NQEVAL steps (0: never after step 0)
    pub nqeval: usize,
    /// FACCUR: the largest rms force error accepted (kJ mol⁻¹ nm⁻¹)
    pub faccur: f64,
}

/// The influence function of a run and the box it was computed in.
#[derive(Debug, Clone)]
struct InfluenceState {
    ghat: Vec<f64>,
    gammahat: Vec<Sym3>,
    edges0: [[f64; 3]; 3],
    /// `(ᵗL₀)⁻¹`
    tl0_inv: M3,
}

/// The reciprocal-space, self, A and surface terms of P3M for one run. The influence function is
/// computed at the first evaluation (gromosXX: step 0) and refused if its rms force error
/// exceeds FACCUR; every NQEVAL steps its error is evaluated in the current box and, if too
/// large, it is computed again (`nonbonded_outerloop.cc:2109-2167`). In a box that differs from
/// the one it was computed in (pressure scaling) it is corrected to first order
/// ([`influence_box_correction`]).
#[derive(Debug, Clone)]
pub struct P3m {
    pub params: LatticeSumParams,
    pub mesh_params: P3mParams,
    fft: Fft3,
    influence: Option<InfluenceState>,
    constants: Option<SelfConstants>,
    /// The rms force error of the influence function, as last estimated.
    pub force_error: f64,
}

impl P3m {
    /// Panics unless every grid dimension is a power of two (the run's checks refuse others).
    pub fn new(params: LatticeSumParams, mesh_params: P3mParams) -> Self {
        Self {
            params,
            fft: Fft3::new(mesh_params.grid),
            mesh_params,
            influence: None,
            constants: None,
            force_error: 0.0,
        }
    }

    fn inputs(&self, edges: &[[f64; 3]; 3], recip: M3) -> InfluenceInputs {
        let m = self.mesh_params;
        InfluenceInputs {
            order: m.order,
            nfdord: m.nfdord,
            nalias: m.nalias,
            shape: self.params.shape,
            width: self.params.width,
            mesh: Mesh::new(edges, m.grid),
            recip,
        }
    }

    fn too_large(&self) -> String {
        format!(
            "P3M: RMS force error is still too big after reevaluation of the influence \
             function ({:e} > FACCUR = {}). Increase the number of grid points, the rms force \
             error threshold, or the charge width parameter.",
            self.force_error, self.mesh_params.faccur
        )
    }

    /// Compute the influence function in the box `edges` (and the force error), or refuse it.
    fn compute_influence(
        &mut self,
        inp: &InfluenceInputs,
        edges: [[f64; 3]; 3],
        volume: f64,
        sum_q2: f64,
        n_atoms: usize,
    ) -> Result<(), String> {
        let inf = influence_function(inp, volume, true);
        self.force_error =
            force_error_estimate(inf.q, sum_q2, volume, n_atoms, self.params.four_pi_eps_i);
        self.influence = Some(InfluenceState {
            ghat: inf.ghat,
            gammahat: inf.gammahat,
            edges0: edges,
            tl0_inv: inverse3(&edges),
        });
        if self.force_error > self.mesh_params.faccur {
            return Err(self.too_large());
        }
        Ok(())
    }

    /// Add the reciprocal-space, self, A and surface terms to `storage` (the virial as `−2 Ξ`)
    /// and return their energies; as [`EwaldSum::evaluate`](super::ewald::EwaldSum::evaluate)
    /// for the Ewald sum. `step` is the MD step (0 first). gromosXX
    /// `_ls_p3m_kspace_outerloop` (`nonbonded_outerloop.cc:2076-2267`); `Ã2` is computed every
    /// step (gromosXX: at step 0, under pressure scaling and on energy-output steps).
    ///
    /// spec: p3m.kspace
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &mut self,
        step: usize,
        positions: &[Vec3],
        charges: &[f64],
        edges: [[f64; 3]; 3],
        rectangular: bool,
        periodicity: &newen_core::math::Periodicity,
        storage: &mut super::ForceStorage,
        virial: bool,
    ) -> Result<LsEnergies, String> {
        let p = self.params;
        let m = self.mesh_params;
        let stale = self.constants.as_ref().is_none_or(|c| c.edges != edges);
        if stale {
            self.constants = Some(self_constants_of_box(&p, edges, rectangular)?);
        }
        let constants = self.constants.clone().expect("set above");
        let (volume, recip) = reciprocal_box(&edges, rectangular);
        let inp = self.inputs(&edges, recip);
        let mesh = inp.mesh;
        let (_, sum_q2) = charge_sums(charges);

        // the influence function: at the first step, then every NQEVAL steps if its error grew
        if self.influence.is_none() {
            self.compute_influence(&inp, edges, volume, sum_q2, positions.len())?;
        } else if m.nqeval > 0 && step > 0 && step.is_multiple_of(m.nqeval) {
            let ghat = self.corrected_ghat(&edges);
            let q = influence_quality(&inp, volume, &ghat);
            self.force_error =
                force_error_estimate(q, sum_q2, volume, positions.len(), p.four_pi_eps_i);
            if self.force_error > m.faccur {
                self.compute_influence(&inp, edges, volume, sum_q2, positions.len())?;
            }
        }
        let ghat = self.corrected_ghat(&edges);
        let gammahat = virial.then(|| self.influence.as_ref().expect("computed").gammahat.clone());

        // the positions put into the box, as gromosXX (`nonbonded_outerloop.cc:2088-2094`)
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
        let n = mesh.len();

        // Ã2 by P3M (NA2CLC 3, and 4 which gromosXX reads as 3)
        let mut a2_tilde = (0.0, [0.0; 6]);
        if matches!(p.na2clc, 3 | 4) {
            let mut rq_re = assign_squared_charges(&mesh, m.order, &in_box, charges);
            let mut rq_im = vec![0.0; n];
            self.fft.forward(&mut rq_re, &mut rq_im);
            a2_tilde = a2_tilde_p3m(&ghat, gammahat.as_deref(), &rq_re, sum_q2);
        }

        // the charge density, its transform, the energy and the potential
        let mut re = assign_charges(&mesh, m.order, &in_box, charges);
        let mut im = vec![0.0; n];
        self.fft.forward(&mut re, &mut im);
        let (energy, kspace_virial) = kspace_energy(
            &ghat,
            gammahat.as_deref(),
            &re,
            &im,
            mesh.cell_volume,
            p.eps0_inv,
        );
        potential_hat(&ghat, &mut re, &mut im, p.eps0_inv);

        // the field on the grid, then the forces
        let forces = if m.nfdord == 0 {
            let mut e = field_ik(&mesh, &recip, &re, &im);
            for (er, ei) in e.iter_mut() {
                self.fft.backward(er, ei);
            }
            interpolate_forces(
                &mesh,
                m.order,
                &in_box,
                charges,
                [&e[0].0, &e[1].0, &e[2].0],
            )
        } else {
            self.fft.backward(&mut re, &mut im);
            let e = field_fd(&mesh, m.nfdord, &re);
            interpolate_forces(&mesh, m.order, &in_box, charges, [&e[0], &e[1], &e[2]])
        };
        for (f, df) in storage.forces.iter_mut().zip(&forces) {
            *f += *df;
        }
        if virial {
            let mtx = super::lattice_sum::sym_to_matrix(&kspace_virial);
            for (row, m_row) in storage.virial.iter_mut().zip(mtx) {
                for (v, x) in row.iter_mut().zip(m_row) {
                    *v += -2.0 * x;
                }
            }
        }
        let mut e = add_constant_terms(
            &p, &constants, positions, charges, edges, a2_tilde, storage, virial,
        );
        e.kspace = energy;
        Ok(e)
    }

    /// `Ĝ` in the box `edges`: as computed, or corrected to first order if the box changed since.
    fn corrected_ghat(&self, edges: &[[f64; 3]; 3]) -> Vec<f64> {
        let inf = self.influence.as_ref().expect("computed");
        if inf.edges0 == *edges {
            return inf.ghat.clone();
        }
        let mut delta = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                delta[i][j] = edges[i][j] - inf.edges0[i][j];
            }
        }
        inf.ghat
            .iter()
            .zip(&inf.gammahat)
            .map(|(&g, gam)| influence_box_correction(g, gam, &inf.tl0_inv, &delta))
            .collect()
    }
}
