//! Lattice-sum electrostatics (NLRELE 2 Ewald, 3 P3M): the scalar kernels both methods share.
//!
//! GROMOS splits the Coulomb energy of a periodic system into a real-space sum over pairs, in
//! which each charge is screened by a charge-shaping function of width `a`, a reciprocal-space
//! sum, and constant self and surface terms (GROMOS Vol. 2, `chg_fe_all`). The functions here are
//! the pieces that do not depend on the reciprocal-space method: the switch function `η(ξ)`,
//! `ξ = r/a`, of each charge shape and its derivative, the Fourier coefficients `γ̂(κ)` of the
//! shape, the real-space pair terms, the self-term constants `A1`, `A2`, `A3`, the reciprocal
//! lattice of the box, the self and `A` terms and the surface term. They follow gromosXX
//! operation by operation and are checked bit for bit against their Lean model
//! (`spec/NewenSpec/LatticeSum.lean`). The Ewald sum over wave vectors is in `ewald.rs`.

use newen_core::math::Vec3;

/// C `erfc` from the platform's libm (glibc on Linux), the function gromosXX calls; Rust's
/// standard library has none.
mod libm {
    extern "C" {
        fn erfc(x: f64) -> f64;
    }

    /// `erfc(x)`.
    #[inline]
    pub fn c_erfc(x: f64) -> f64 {
        // SAFETY: a pure C math function of one double.
        unsafe { erfc(x) }
    }
}

/// The charge-shaping function, LONGRANGE `NSHAPE` (GROMOS Vol. 2, table `ls_charge_shaping`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChargeShape {
    /// `NSHAPE = −1`: Gaussian, `η(ξ) = erfc(ξ)`.
    Gaussian,
    /// `NSHAPE = 0..=10`: the truncated polynomial of that order, zero beyond `ξ = 1`.
    Polynomial(u8),
}

impl ChargeShape {
    /// The shape of `NSHAPE`, or `None` outside `−1..=10`.
    pub fn from_nshape(nshape: i32) -> Option<Self> {
        match nshape {
            -1 => Some(ChargeShape::Gaussian),
            0..=10 => Some(ChargeShape::Polynomial(nshape as u8)),
            _ => None,
        }
    }
}

fn unknown_order(order: u8) -> ! {
    panic!("charge shape of order {order}: only NSHAPE -1..=10 exist (ChargeShape::from_nshape)")
}

/// `(η(ξ), η'(ξ))` of the switch function at `ξ = r/a`, gromosXX `charge_shape_switch`
/// (`latticesum.h:86-252`) in its operation order. With `m = 1 − ξ`:
///
/// ```text
/// Gaussian:  η = erfc(ξ),   η' = exp(−ξ·ξ) · −2 / sqrt(π)
/// m < 0:     (0, 0)                                  (the Heaviside factor H(1 − ξ))
/// order 0:   η = 1/2 · m² · (ξ + 2),          η' = 3/2 · (ξ·ξ − 1)
/// order 1:   η = m³ · (ξ + 1),                η' = −2 · m² · (1 + 2ξ)
/// …          (orders 2 to 10 as in Vol. 2, tables ls_switch_function and ls_switch_prime,
///             the polynomials in Horner form)
/// ```
///
/// GROMOS Vol. 2: `eta_function_def`. The real-space interaction of a pair is `η(r/a)/r`.
///
/// spec: lattice_sum.shape_switch
pub fn shape_switch(shape: ChargeShape, xi: f64) -> (f64, f64) {
    let order = match shape {
        ChargeShape::Gaussian => {
            let eta = libm::c_erfc(xi);
            let d_eta = (-xi * xi).exp() * -2.0 / std::f64::consts::PI.sqrt();
            return (eta, d_eta);
        },
        ChargeShape::Polynomial(order) => order,
    };
    let m = 1.0 - xi;
    if m < 0.0 {
        return (0.0, 0.0);
    }
    match order {
        0 => {
            let m2 = m * m;
            (1.0 / 2.0 * m2 * (xi + 2.0), 3.0 / 2.0 * (xi * xi - 1.0))
        },
        1 => {
            let m2 = m * m;
            let m3 = m2 * m;
            (m3 * (xi + 1.0), -2.0 * m2 * (1.0 + 2.0 * xi))
        },
        2 => {
            let m2 = m * m;
            let m3 = m2 * m;
            let m4 = m2 * m2;
            (
                1.0 / 2.0 * m4 * (3.0 * xi + 2.0),
                -5.0 / 2.0 * m3 * (1.0 + 3.0 * xi),
            )
        },
        3 => {
            let m2 = m * m;
            let m3 = m2 * m;
            let m4 = m2 * m2;
            (
                1.0 / 4.0 * m4 * ((4.0 * xi + 7.0) * xi + 4.0),
                -3.0 / 4.0 * m3 * ((8.0 * xi + 9.0) * xi + 3.0),
            )
        },
        4 => {
            let m2 = m * m;
            let m4 = m2 * m2;
            let m5 = m4 * m;
            (
                1.0 / 8.0 * m5 * ((15.0 * xi + 19.0) * xi + 8.0),
                -21.0 / 8.0 * m4 * (1.0 + xi * (4.0 + 5.0 * xi)),
            )
        },
        5 => {
            let m2 = m * m;
            let m4 = m2 * m2;
            let m5 = m4 * m;
            let m6 = m4 * m2;
            (
                m6 * ((3.0 * xi + 3.0) * xi + 1.0),
                -3.0 * m5 * (1.0 + xi * (5.0 + 8.0 * xi)),
            )
        },
        6 => {
            let m2 = m * m;
            let m5 = m2 * m2 * m;
            let m6 = m5 * m;
            (
                1.0 / 16.0 * m6 * (((35.0 * xi + 66.0) * xi + 51.0) * xi + 16.0),
                -9.0 / 16.0 * m5 * (5.0 + xi * (25.0 + xi * (47.0 + 35.0 * xi))),
            )
        },
        7 => {
            let m2 = m * m;
            let m6 = m2 * m2 * m2;
            let m7 = m6 * m;
            (
                1.0 / 8.0 * m7 * (((32.0 * xi + 49.0) * xi + 31.0) * xi + 8.0),
                -5.0 / 8.0 * m6 * (5.0 + xi * (30.0 + xi * (69.0 + 64.0 * xi))),
            )
        },
        8 => {
            let m3 = m * m * m;
            let m7 = m3 * m3 * m;
            let m8 = m7 * m;
            (
                1.0 / 16.0 * m8 * (((105.0 * xi + 136.0) * xi + 73.0) * xi + 16.0),
                -55.0 / 16.0 * m7 * (1.0 + 3.0 * xi) * (1.0 + xi * (4.0 + 7.0 * xi)),
            )
        },
        9 => {
            let m3 = m * m * m;
            let m7 = m3 * m3 * m;
            let m8 = m7 * m;
            (
                1.0 / 32.0 * m8 * ((((160.0 * xi + 335.0) * xi + 312.0) * xi + 151.0) * xi + 32.0),
                -15.0 / 32.0 * m7 * (7.0 + xi * (49.0 + xi * (141.0 + xi * (203.0 + 128.0 * xi)))),
            )
        },
        10 => {
            let m2 = m * m;
            let m4 = m2 * m2;
            let m8 = m4 * m4;
            let m9 = m8 * m;
            (
                1.0 / 128.0
                    * m9
                    * ((((1155.0 * xi + 2075.0) * xi + 1665.0) * xi + 697.0) * xi + 128.0),
                -65.0 / 128.0
                    * m8
                    * (7.0 + xi * (56.0 + 3.0 * xi * (62.0 + xi * (104.0 + 77.0 * xi)))),
            )
        },
        _ => unknown_order(order),
    }
}

/// `(γ̂(κ), γ̂'(κ))`, the Fourier coefficient of the charge shape at `κ = a·k` and its
/// derivative, gromosXX `charge_shape_fourier` (`latticesum.h:309-467`) in its operation order,
/// with `C = cos κ`, `S = sin κ`, `kn = κⁿ`:
///
/// ```text
/// Gaussian:  q = κ · κ · −0.25,   γ̂ = exp(q),   γ̂' = −0.5 · κ · exp(q)
/// order 0:   γ̂  = 3 · (S − κC) / k3,            γ̂' = 3 · (3κC + (k2 − 3) S) / k4
/// order 1:   γ̂  = 12 · (2 − 2C − κS) / k4,      γ̂' = 12 · (−8 − (k2 − 8) C + 5κS) / k5
/// …          (orders 2 to 10 as in Vol. 2, tables ls_fcoef and ls_fcoef_prime)
/// ```
///
/// The polynomial shapes divide by `κ^(N+3)`: `κ = 0` gives NaN (Vol. 2 sets `γ̂(0) = 1`; the
/// `k = 0` term is never summed). GROMOS Vol. 2: `latsum_gamma_ft`.
///
/// spec: lattice_sum.shape_fourier
pub fn shape_fourier(shape: ChargeShape, kappa: f64) -> (f64, f64) {
    let order = match shape {
        ChargeShape::Gaussian => {
            let k2_quarter = kappa * kappa * -0.25;
            return (k2_quarter.exp(), -0.5 * kappa * k2_quarter.exp());
        },
        ChargeShape::Polynomial(order) => order,
    };
    match order {
        0 => {
            let k2 = kappa * kappa;
            let k3 = k2 * kappa;
            let s = kappa.sin();
            let c = kappa.cos();
            let kappa_c = kappa * c;
            let k4 = k2 * k2;
            (
                3.0 * (s - kappa_c) / k3,
                3.0 * (3.0 * kappa_c + (k2 - 3.0) * s) / k4,
            )
        },
        1 => {
            let k2 = kappa * kappa;
            let k4 = k2 * k2;
            let c = kappa.cos();
            let kappa_s = kappa * kappa.sin();
            let k5 = k4 * kappa;
            (
                12.0 * (2.0 - 2.0 * c - kappa_s) / k4,
                12.0 * (-8.0 - (k2 - 8.0) * c + 5.0 * kappa_s) / k5,
            )
        },
        2 => {
            let k2 = kappa * kappa;
            let k5 = k2 * k2 * kappa;
            let c = kappa.cos();
            let s = kappa.sin();
            let k6 = k5 * kappa;
            (
                60.0 * (2.0 * kappa + kappa * c - 3.0 * s) / k5,
                60.0 * (-8.0 * kappa - 7.0 * kappa * c - (k2 - 15.0) * s) / k6,
            )
        },
        3 => {
            let k2 = kappa * kappa;
            let k6 = k2 * k2 * k2;
            let c = kappa.cos();
            let kappa_s = kappa * kappa.sin();
            let k7 = k6 * kappa;
            (
                90.0 * (8.0 + (k2 - 8.0) * c - 5.0 * kappa_s) / k6,
                90.0 * (-48.0 - (9.0 * k2 - 48.0) * c - (k2 - 33.0) * kappa_s) / k7,
            )
        },
        4 => {
            let k2 = kappa * kappa;
            let k7 = k2 * k2 * k2 * kappa;
            let kappa_c = kappa * kappa.cos();
            let s = kappa.sin();
            let k8 = k7 * kappa;
            (
                630.0 * (8.0 * kappa + 7.0 * kappa_c + (k2 - 15.0) * s) / k7,
                630.0 * (-48.0 * kappa + (k2 - 57.0) * kappa_c - 3.0 * (4.0 * k2 - 35.0) * s) / k8,
            )
        },
        5 => {
            let k2 = kappa * kappa;
            let k4 = k2 * k2;
            let k8 = k4 * k4;
            let c = kappa.cos();
            let kappa_s = kappa * kappa.sin();
            let k9 = k8 * kappa;
            (
                5040.0 * (4.0 * (k2 - 6.0) - (k2 - 24.0) * c + 9.0 * kappa_s) / k8,
                5040.0 * (-24.0 * (k2 - 8.0) + 3.0 * (5.0 * k2 - 64.0) * c + (k2 - 87.0) * kappa_s)
                    / k9,
            )
        },
        6 => {
            let k2 = kappa * kappa;
            let k4 = k2 * k2;
            let k9 = k4 * k4 * kappa;
            let kappa_c = kappa * kappa.cos();
            let s = kappa.sin();
            let k10 = k9 * kappa;
            (
                7560.0 * (48.0 * kappa - (k2 - 57.0) * kappa_c + 3.0 * (4.0 * k2 - 35.0) * s) / k9,
                7560.0
                    * (-384.0 * kappa
                        + 3.0 * (6.0 * k2 - 187.0) * kappa_c
                        + (k4 - 141.0 * k2 + 945.0) * s)
                    / k10,
            )
        },
        7 => {
            let k2 = kappa * kappa;
            let k4 = k2 * k2;
            let k10 = k4 * k4 * k2;
            let c = kappa.cos();
            let kappa_s = kappa * kappa.sin();
            let k11 = k10 * kappa;
            (
                75600.0 * (24.0 * (k2 - 8.0) - 3.0 * (5.0 * k2 - 64.0) * c - (k2 - 87.0) * kappa_s)
                    / k10,
                75600.0
                    * (-192.0 * (k2 - 10.0) - (k4 - 207.0 * k2 + 1920.0) * c
                        + (22.0 * k2 - 975.0) * kappa_s)
                    / k11,
            )
        },
        8 => {
            let k2 = kappa * kappa;
            let k4 = k2 * k2;
            let k11 = k4 * k4 * k2 * kappa;
            let kappa_c = kappa * kappa.cos();
            let s = kappa.sin();
            let k12 = k11 * kappa;
            (
                831600.0
                    * (8.0 * kappa * (k2 - 24.0) + (k2 - 123.0) * kappa_c
                        - (18.0 * k2 - 315.0) * s)
                    / k11,
                831600.0
                    * (-64.0 * kappa * (k2 - 30.0)
                        - (26.0 * k2 - 1545.0) * kappa_c
                        - (k4 - 285.0 * k2 + 3465.0) * s)
                    / k12,
            )
        },
        9 => {
            let k2 = kappa * kappa;
            let k4 = k2 * k2;
            let k12 = k4 * k4 * k4;
            let c = kappa.cos();
            let kappa_s = kappa * kappa.sin();
            let k13 = k12 * kappa;
            (
                1247400.0
                    * (192.0 * (k2 - 10.0) + (k4 - 207.0 * k2 + 1920.0) * c
                        - (22.0 * k2 - 975.0) * kappa_s)
                    / k12,
                1247400.0
                    * (-1920.0 * (k2 - 12.0)
                        - 15.0 * (2.0 * k4 - 203.0 * k2 + 1536.0) * c
                        - (k4 - 405.0 * k2 + 12645.0) * kappa_s)
                    / k13,
            )
        },
        10 => {
            let k2 = kappa * kappa;
            let k4 = k2 * k2;
            let k13 = k4 * k4 * k4 * kappa;
            let kappa_c = kappa * kappa.cos();
            let s = kappa.sin();
            let k14 = k13 * kappa;
            (
                16216200.0
                    * (64.0 * kappa * (k2 - 30.0)
                        + (26.0 * k2 - 1545.0) * kappa_c
                        + (k4 - 285.0 * k2 + 3465.0) * s)
                    / k13,
                16216200.0
                    * (-640.0 * kappa * (k2 - 36.0) + (k4 - 545.0 * k2 + 22005.0) * kappa_c
                        - 5.0 * (7.0 * k4 - 936.0 * k2 + 9009.0) * s)
                    / k14,
            )
        },
        _ => unknown_order(order),
    }
}

/// One LJ + real-space lattice-sum pair at squared distance `r2`: `(f_mag, e_lj, e_ls)`,
/// gromosXX `lj_ls_interaction` (`nonbonded_term.cc:375-387`) in its operation order, with
/// `q_eps = q_i q_j / (4πε₀)` and `inv_width = 1/a`:
///
/// ```text
/// dist = sqrt(r2),  dist2i = 1 / r2,  dist6i = dist2i · dist2i · dist2i,  disti = 1 / dist
/// c12_dist6i = c12 · dist6i
/// e_lj  = (c12_dist6i − c6) · dist6i                             C12/r¹² − C6/r⁶
/// (η, η') = shape_switch(shape, inv_width · dist)
/// e_ls  = q_eps · disti · η                                      q η(r/a) / r
/// f_ls  = dist2i · (e_ls − q_eps · η' · inv_width)               −(d e_ls/dr) / r
/// f_mag = (c12_dist6i + c12_dist6i − c6) · 6 · dist6i · dist2i + f_ls
/// ```
///
/// The force on atom i is `r_vec · f_mag`. There is no zero-distance guard (gromosXX asserts
/// `r² ≠ 0`). GROMOS Vol. 2: `e_eta_approx` (first sum), `e_eta_approx_force`.
///
/// spec: lattice_sum.lj_ls_pair
#[inline]
pub fn lj_ls_pair(
    r2: f64,
    c6: f64,
    c12: f64,
    q_eps: f64,
    shape: ChargeShape,
    inv_width: f64,
) -> (f64, f64, f64) {
    let dist = r2.sqrt();
    let dist2i = 1.0 / r2;
    let dist6i = dist2i * dist2i * dist2i;
    let disti = 1.0 / dist;
    let c12_dist6i = c12 * dist6i;
    let ai_dist = inv_width * dist;

    let e_lj = (c12_dist6i - c6) * dist6i;
    let (eta, d_eta) = shape_switch(shape, ai_dist);
    let e_ls = q_eps * disti * eta;
    let f_ls = dist2i * (e_ls - q_eps * d_eta * inv_width);
    let f_mag = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i + f_ls;
    (f_mag, e_lj, e_ls)
}

/// The real-space lattice-sum term of an excluded pair at squared distance `r2`:
/// `(f_mag, e_ls)`, gromosXX `ls_excluded_interaction` (`nonbonded_term.cc:393-411`) in its
/// operation order:
///
/// ```text
/// dist = sqrt(r2),  dist2i = 1 / r2,  disti = 1 / dist
/// (η, η') = shape_switch(shape, inv_width · dist)
/// e_ls  = q_eps · disti · (η − 1)                                q (η(r/a) − 1) / r
/// f_mag = dist2i · (e_ls − q_eps · η' · inv_width)
/// ```
///
/// An excluded pair has no Coulomb interaction, but the reciprocal-space sum includes it; this
/// term removes it again. GROMOS Vol. 2: `e_eta_approx` (second sum).
///
/// spec: lattice_sum.ls_excluded_pair
#[inline]
pub fn ls_excluded_pair(r2: f64, q_eps: f64, shape: ChargeShape, inv_width: f64) -> (f64, f64) {
    let dist = r2.sqrt();
    let dist2i = 1.0 / r2;
    let disti = 1.0 / dist;
    let ai_dist = inv_width * dist;

    let (eta, d_eta) = shape_switch(shape, ai_dist);
    let e_ls = q_eps * disti * (eta - 1.0);
    let f_mag = dist2i * (e_ls - q_eps * d_eta * inv_width);
    (f_mag, e_ls)
}

/// The self-term constants `(A1, A3)` of a charge shape of width `a` (`width`) in a box of
/// volume `V`, gromosXX `Nonbonded_Outerloop::ls_self_outerloop`
/// (`nonbonded_outerloop.cc:2291-2353`) in its operation order:
///
/// ```text
/// (ã1, ã3) from the table (Gaussian: ã1 = 1, ã3 = 2 · 1 / sqrt(π); order 0: 2/5, 3/2; …)
/// A1 = ã1 · (−π · a · a / V)
/// A3 = ã3 · (−1 / a)
/// ```
///
/// GROMOS Vol. 2: `latsum_a1_def`, `latsum_a3_def`, table `ls_acst` (`−V π⁻¹ a⁻² A1`, `−a A3`).
///
/// spec: lattice_sum.self_constants
pub fn self_constants(shape: ChargeShape, width: f64, volume: f64) -> (f64, f64) {
    let (a1, a3) = match shape {
        ChargeShape::Gaussian => (1.0, 2.0 * 1.0 / std::f64::consts::PI.sqrt()),
        ChargeShape::Polynomial(order) => match order {
            0 => (2.0 / 5.0, 3.0 / 2.0),
            1 => (4.0 / 15.0, 2.0),
            2 => (4.0 / 21.0, 5.0 / 2.0),
            3 => (3.0 / 14.0, 9.0 / 4.0),
            4 => (1.0 / 6.0, 21.0 / 8.0),
            5 => (2.0 / 15.0, 3.0),
            6 => (8.0 / 55.0, 45.0 / 16.0),
            7 => (4.0 / 33.0, 25.0 / 8.0),
            8 => (4.0 / 39.0, 55.0 / 16.0),
            9 => (10.0 / 91.0, 105.0 / 32.0),
            10 => (2.0 / 21.0, 455.0 / 128.0),
            _ => unknown_order(order),
        },
    };
    (
        a1 * (-std::f64::consts::PI * width * width / volume),
        a3 * (-1.0 / width),
    )
}

/// A symmetric 3×3 tensor as gromosXX's `math::SymmetricMatrix` stores it: `[xx, xy, xz, yy, yz,
/// zz]`.
pub type Sym3 = [f64; 6];

/// `v1 ⊗ v2` for the upper triangle, gromosXX `math::symmetric_tensor_product`.
#[inline]
pub fn sym_tensor(v1: [f64; 3], v2: [f64; 3]) -> Sym3 {
    [
        v1[0] * v2[0],
        v1[0] * v2[1],
        v1[0] * v2[2],
        v1[1] * v2[1],
        v1[1] * v2[2],
        v1[2] * v2[2],
    ]
}

#[inline]
fn sym_add(m: &mut Sym3, t: &Sym3) {
    for (a, b) in m.iter_mut().zip(t) {
        *a += b;
    }
}

#[inline]
fn add_to_diagonal(m: &mut Sym3, d: f64) {
    m[0] += d;
    m[3] += d;
    m[5] += d;
}

/// The full matrix of a symmetric tensor, `[a][b]`.
pub fn sym_to_matrix(m: &Sym3) -> [[f64; 3]; 3] {
    [[m[0], m[1], m[2]], [m[1], m[3], m[4]], [m[2], m[4], m[5]]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The volume `V` of the box with edge vectors `a, b, c` and the vectors `2π/V · (b×c, c×a,
/// a×b)` of its reciprocal lattice, gromosXX `math::volume` (`math/volume.cc:33-52`) and
/// `KSpace_Utils::l_to_k_matrix` (`kspace.cc:128-136`) in their operation order:
///
/// ```text
/// rectangular:  V = |a| · |b| · |c|          (|v| = sqrt(v_x² + v_y² + v_z²))
/// triclinic:    V = (a × b) · c
/// s = 2π / V,   reciprocal = ((b × c)·s, (c × a)·s, (a × b)·s)
/// ```
///
/// The wave vector of the integer triple `l` is `l_x b₁ + l_y b₂ + l_z b₃`
/// ([`lattice_vector`]), so that `k · a = 2π l_x`. GROMOS Vol. 2 after `e_gamma_ewald`
/// (`k = 2π l` in the reciprocal basis).
///
/// spec: lattice_sum.reciprocal_box
pub fn reciprocal_box(edges: &[[f64; 3]; 3], rectangular: bool) -> (f64, [[f64; 3]; 3]) {
    let [a, b, c] = *edges;
    let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let volume = if rectangular {
        norm(a) * norm(b) * norm(c)
    } else {
        let ab = cross(a, b);
        ab[0] * c[0] + ab[1] * c[1] + ab[2] * c[2]
    };
    let s = 2.0 * std::f64::consts::PI / volume;
    let scale = |v: [f64; 3]| [v[0] * s, v[1] * s, v[2] * s];
    (
        volume,
        [scale(cross(b, c)), scale(cross(c, a)), scale(cross(a, b))],
    )
}

/// `k = l_x b₁ + l_y b₂ + l_z b₃` for the reciprocal vectors `b₁, b₂, b₃` of [`reciprocal_box`],
/// component by component as `(b₁ l_x + b₂ l_y) + b₃ l_z`.
///
/// gromosXX (`kspace.cc:88`, `nonbonded_outerloop.cc:2449`) multiplies with the transpose
/// (`k_x = b₁·l`), the same vector for a rectangular box, where the reciprocal matrix is
/// diagonal, and not a reciprocal-lattice vector for a triclinic one.
#[inline]
pub fn lattice_vector(recip: &[[f64; 3]; 3], l: [i64; 3]) -> [f64; 3] {
    let [u, v, w] = recip;
    let (x, y, z) = (l[0] as f64, l[1] as f64, l[2] as f64);
    [
        u[0] * x + v[0] * y + w[0] * z,
        u[1] * x + v[1] * y + w[1] * z,
        u[2] * x + v[2] * y + w[2] * z,
    ]
}

/// The Wigner constant of the simple cubic lattice, gromosXX `nonbonded_outerloop.cc:2379`
/// (GROMOS Vol. 2 prints `−2.83792748`, two digits transposed).
pub const WIGNER_CUBE: f64 = -2.83729748;

/// `A2` in a cube of edge `L`, and its derivative (the isotropic diagonal of `∂A2/∂ε`), gromosXX
/// `ls_self_outerloop` (`nonbonded_outerloop.cc:2376-2385`):
///
/// ```text
/// A2 = wigner / L − A1 − A3,     dA2 = wigner / (3 L) − A1
/// ```
///
/// GROMOS Vol. 2: `A1 + A2 + A3 = ξ_EW / L` for a cube (after `e_slf`).
///
/// spec: lattice_sum.a2_cubic
pub fn a2_cubic(wigner: f64, edge: f64, a1: f64, a3: f64) -> (f64, f64) {
    (wigner / edge - a1 - a3, wigner / (3.0 * edge) - a1)
}

/// The result of [`a2_shells`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct A2Shells {
    /// `A2 = 4π/V · Σ_{l ≠ 0} γ̂(a|k|)/k²`.
    pub a2: f64,
    /// The number of shells summed.
    pub shells: u32,
    /// Whether the last shell changed the sum by at most `tola2` relative.
    pub converged: bool,
    /// `4π/V · Σ (2 (a|k| γ̂' − 2γ̂)/k⁴ k⊗k + γ̂/k² 1)`, the derivative gromosXX's virial uses.
    pub derivative: Sym3,
}

/// `A2` summed over cubic shells of growing `l_max` in the integer lattice, gromosXX
/// `ls_self_outerloop` (`nonbonded_outerloop.cc:2387-2537`) in its loop order. Shell `n` holds
/// the triples with `max |l_i| = n`, enumerated as three faces: `l_x = n` with
/// `|l_y|, |l_z| ≤ n`; `l_y = n` with `|l_x| ≤ n − 1`, `|l_z| ≤ n`; `l_z = n` with
/// `|l_x|, |l_y| ≤ n − 1`. Those are half the shell; `−l` gives the same term, so the face sum is
/// doubled:
///
/// ```text
/// per l:      k = lattice_vector(l),  k2 = |k|²,  ak = sqrt(k2) · a
///             term += γ̂(ak) / k2
///             S += (2 (ak γ̂'(ak) − 2 γ̂(ak)) / (k2 k2)) k ⊗ k
/// per shell:  term ·= 2;  A2 += term;  S += term 1;  stop when |term / A2| ≤ tola2
/// A2 ·= 4π / V,   derivative = S · 4π / V
/// ```
///
/// gromosXX loops until the tolerance is met; here at most `max_shells` shells are summed and
/// `converged` says whether it was met. GROMOS Vol. 2: `latsum_a2_calc`.
///
/// spec: lattice_sum.a2_shells
pub fn a2_shells(
    shape: ChargeShape,
    width: f64,
    recip: &[[f64; 3]; 3],
    volume: f64,
    tola2: f64,
    max_shells: u32,
) -> A2Shells {
    let mut sum_gammahat: Sym3 = [0.0; 6];
    let mut l_max: i64 = 0;
    let mut a2 = 0.0;
    let converged = loop {
        l_max += 1;
        let mut term = 0.0;
        for coord in 0..3 {
            let (coord_a, coord_b, bound_a, bound_b) = match coord {
                0 => (1, 2, l_max, l_max),
                1 => (0, 2, l_max - 1, l_max),
                _ => (0, 1, l_max - 1, l_max - 1),
            };
            for l_a in -bound_a..=bound_a {
                let mut l = [0i64; 3];
                l[coord] = l_max;
                l[coord_a] = l_a;
                for l_b in -bound_b..=bound_b {
                    l[coord_b] = l_b;
                    let k = lattice_vector(recip, l);
                    let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                    let ak = k2.sqrt() * width;
                    let (gamma_hat, gamma_hat_prime) = shape_fourier(shape, ak);
                    term += gamma_hat / k2;
                    let iso = 2.0 * (ak * gamma_hat_prime - 2.0 * gamma_hat) / (k2 * k2);
                    sym_add(
                        &mut sum_gammahat,
                        &sym_tensor([iso * k[0], iso * k[1], iso * k[2]], k),
                    );
                }
            }
        }
        term *= 2.0;
        a2 += term;
        add_to_diagonal(&mut sum_gammahat, term);
        let tolerance = (term / a2).abs();
        // gromosXX `while (tolerance > required_precision)`: a NaN tolerance stops too
        let another_shell = tolerance > tola2;
        if !another_shell {
            break true;
        }
        if l_max >= i64::from(max_shells) {
            break false;
        }
    };
    let s = 4.0 * std::f64::consts::PI / volume;
    A2Shells {
        a2: a2 * s,
        shells: l_max as u32,
        converged,
        derivative: sum_gammahat.map(|x| x * s),
    }
}

/// `(Σ q_i, Σ q_i²)` in atom order, gromosXX `Topology::squared_sum_charges` (before its final
/// square) and `sum_squared_charges` (`topology.cc:1572-1586`).
pub fn charge_sums(charges: &[f64]) -> (f64, f64) {
    let mut sum = 0.0;
    let mut sum_sq = 0.0;
    for &q in charges {
        sum += q;
        sum_sq += q * q;
    }
    (sum, sum_sq)
}

/// The self term and the `A` term of the lattice sum, gromosXX `ls_self_outerloop`
/// (`nonbonded_outerloop.cc:2588, 2605`) in its operation order, with `S = Σ q_i` (`sum_q`),
/// `S̃² = Σ q_i²` (`sum_q2`) and `eps0_inv = 1/ε₀`:
///
/// ```text
/// ls_self   = (A1 + A2 + A3) · S̃² · eps0_inv / (8π)
/// ls_a_term = (A1 · S·S − (A1 + Ã2) · S̃²) · eps0_inv / (8π)
/// ```
///
/// GROMOS Vol. 2: `e_slf`, `e_a_def`.
///
/// spec: lattice_sum.self_and_a_terms
pub fn self_and_a_terms(
    a1: f64,
    a2: f64,
    a3: f64,
    a2_tilde: f64,
    sum_q: f64,
    sum_q2: f64,
    eps0_inv: f64,
) -> (f64, f64) {
    let s2 = sum_q * sum_q;
    let eight_pi = 8.0 * std::f64::consts::PI;
    (
        (a1 + a2 + a3) * sum_q2 * eps0_inv / eight_pi,
        (a1 * s2 - (a1 + a2_tilde) * sum_q2) * eps0_inv / eight_pi,
    )
}

/// The virials (gromosXX's `Ξ`, before the storage factor −2) of the self term and of the `A`
/// term, gromosXX `ls_self_outerloop` (`nonbonded_outerloop.cc:2357-2366, 2384, 2519, 2555-2607`):
///
/// ```text
/// self     = −¼ fpe A1 S̃² 1 + (−¼ S̃² fpe) dA2
/// constant = −¼ fpe A1 (S² − S̃²) 1 − (−¼ S̃² fpe) dÃ2
/// ```
///
/// with `fpe = 1/(4πε₀)` and `dA2`, `dÃ2` the derivatives of `A2` and `Ã2` (the NA2CLC choice is
/// the caller's).
pub fn self_and_a_virials(
    a1: f64,
    sum_q: f64,
    sum_q2: f64,
    d_a2: &Sym3,
    d_a2_tilde: &Sym3,
    four_pi_eps_i: f64,
) -> (Sym3, Sym3) {
    let s2 = sum_q * sum_q;
    let mut self_virial: Sym3 = [0.0; 6];
    add_to_diagonal(&mut self_virial, -0.25 * four_pi_eps_i * a1 * sum_q2);
    let mut constant_virial: Sym3 = [0.0; 6];
    add_to_diagonal(
        &mut constant_virial,
        -0.25 * four_pi_eps_i * a1 * (s2 - sum_q2),
    );
    let f = -0.25 * sum_q2 * four_pi_eps_i;
    for i in 0..6 {
        self_virial[i] += f * d_a2[i];
        constant_virial[i] -= f * d_a2_tilde[i];
    }
    (self_virial, constant_virial)
}

/// The result of [`surface_term`].
#[derive(Debug, Clone, PartialEq)]
pub struct Surface {
    pub energy: f64,
    pub forces: Vec<Vec3>,
    /// `Ξ = ½ p M ⊗ M − ¼ p |M|² 1`, the virial of the term (GROMOS sign, before the storage
    /// factor −2).
    pub virial: Sym3,
}

/// The surface term of the box dipole, gromosXX `_ls_surface_outerloop`
/// (`nonbonded_outerloop.cc:2652-2688`) in its operation order, with the unwrapped positions,
/// the box centre `c` and `eps0_inv = 1/ε₀`:
///
/// ```text
/// EPSLS = 0 (conducting boundary):  nothing
/// M   = Σ_i q_i (r_i − c)
/// p   = eps0_inv / ((EPSLS · 2 + 1) · V)
/// E   = ½ |M|² p,     F_i = −p q_i M
/// Ξ_ab = ½ p M_a M_b + δ_ab (−¼ p |M|²)
/// ```
///
/// gromosXX computes the diagonal of `Ξ` only and adds it to its virial storage as it is, where
/// every other lattice-sum term adds `−2 Ξ`; this returns the whole tensor and the caller stores
/// `−2 Ξ`. GROMOS Vol. 2: `e_srf`, `e_srf_forces`.
///
/// spec: lattice_sum.surface_term
pub fn surface_term(
    positions: &[Vec3],
    charges: &[f64],
    centre: Vec3,
    eps_ls: f64,
    volume: f64,
    eps0_inv: f64,
) -> Surface {
    let n = positions.len();
    if eps_ls == 0.0 {
        return Surface {
            energy: 0.0,
            forces: vec![Vec3::ZERO; n],
            virial: [0.0; 6],
        };
    }
    let mut m = [0.0f64; 3];
    for (r, &q) in positions.iter().zip(charges) {
        m[0] += q * (r.x - centre.x);
        m[1] += q * (r.y - centre.y);
        m[2] += q * (r.z - centre.z);
    }
    let prefactor = eps0_inv / ((eps_ls * 2.0 + 1.0) * volume);
    let abs2 = m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
    let energy = 0.5 * abs2 * prefactor;
    let forces = charges
        .iter()
        .map(|&q| {
            let s = -prefactor * q;
            Vec3::new(s * m[0], s * m[1], s * m[2])
        })
        .collect();
    let isotropic = -0.25 * prefactor * abs2;
    let half_p = 0.5 * prefactor;
    let mut virial = sym_tensor([half_p * m[0], half_p * m[1], half_p * m[2]], m);
    // gromosXX's diagonal: isotropic + ½ p M_i M_i
    virial[0] += isotropic;
    virial[3] += isotropic;
    virial[5] += isotropic;
    Surface {
        energy,
        forces,
        virial,
    }
}
