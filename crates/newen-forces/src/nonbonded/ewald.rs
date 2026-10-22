//! The direct Ewald sum (NONBONDED NLRELE = 2): the reciprocal-space part of the lattice sum as
//! an explicit sum over wave vectors.
//!
//! GROMOS Vol. 2 (`e_gamma_ewald`, `ewald_reci_force`): with `C(k) = Σ_i q_i cos(k·r_i)` and
//! `S(k) = Σ_i q_i sin(k·r_i)`, the energy is `(2ε₀V)⁻¹ Σ_k γ̂(a|k|)/k² (C(k)² + S(k)²)` over the
//! wave vectors `k ≠ 0` of the box's reciprocal lattice inside the cutoff. [`k_vectors`] builds
//! that set once per box, [`ewald_kspace`] evaluates the sum. Both follow gromosXX
//! (`configuration/kspace.cc`, `Nonbonded_Outerloop::_ls_ewald_kspace_outerloop`) in its loop
//! and operation order and are checked bit for bit against their Lean model
//! (`spec/NewenSpec/LatticeSum.lean`).

use super::lattice_sum::{lattice_vector, shape_fourier, sym_tensor, ChargeShape, Sym3};
use newen_core::math::Vec3;

/// One wave vector of the Ewald sum with the quantities gromosXX caches for it
/// (`configuration::KSpace_Element`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KVector {
    pub k: [f64; 3],
    /// `|k|²`
    pub k2: f64,
    /// `1/k²`
    pub k2i: f64,
    /// `a|k|`
    pub ak: f64,
    /// `γ̂(a|k|)`
    pub gamma_hat: f64,
    /// `γ̂(a|k|)/k²`, as `k2i · γ̂`
    pub k2i_gamma_hat: f64,
    /// `a|k| γ̂'(a|k|)`
    pub ak_gamma_hat_prime: f64,
}

/// The wave vectors of the Ewald sum, gromosXX `configuration::calculate_k_space`
/// (`kspace.cc:48-126`) in its loop order: every integer triple `l` with
/// `|l_x| ≤ NKX, |l_y| ≤ NKY, |l_z| ≤ NKZ` (`l_x` slowest), `k = lattice_vector(l)`, kept when
///
/// ```text
/// k2 = k_x·k_x + k_y·k_y + k_z·k_z,    k2 ≤ KCUT  and  k2 ≠ 0
/// ```
///
/// then `|k| = sqrt(k2)`, `ak = |k| · a`, `k2i = 1/k2`, `(γ̂, γ̂') = shape_fourier(ak)`,
/// `k2i_gamma_hat = k2i · γ̂`, `ak_gamma_hat_prime = ak · γ̂'`.
///
/// gromosXX compares `|k|²` with KCUT itself (`kspace.cc:72, 90`, where the variable is named
/// `k_cut2` and commented "k cutoff squared"); GROMOS Vol. 2 cuts at `|k| ≤ k_c` and Vol. 4 calls
/// KCUT the "k-space cutoff". This follows the code. An empty set is the caller's error
/// (gromosXX `kspace.cc:122-125`).
///
/// spec: ewald.k_vectors
pub fn k_vectors(
    nk: [i64; 3],
    kcut: f64,
    recip: &[[f64; 3]; 3],
    shape: ChargeShape,
    width: f64,
) -> Vec<KVector> {
    let mut out = Vec::new();
    for lx in -nk[0]..=nk[0] {
        for ly in -nk[1]..=nk[1] {
            for lz in -nk[2]..=nk[2] {
                let k = lattice_vector(recip, [lx, ly, lz]);
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                if k2 <= kcut && k2 != 0.0 {
                    let abs_k = k2.sqrt();
                    let ak = abs_k * width;
                    let k2i = 1.0 / k2;
                    let (gamma_hat, gamma_hat_prime) = shape_fourier(shape, ak);
                    out.push(KVector {
                        k,
                        k2,
                        k2i,
                        ak,
                        gamma_hat,
                        k2i_gamma_hat: k2i * gamma_hat,
                        ak_gamma_hat_prime: ak * gamma_hat_prime,
                    });
                }
            }
        }
    }
    out
}

/// The result of [`ewald_kspace`].
#[derive(Debug, Clone, PartialEq)]
pub struct KSpace {
    /// `E_γ`
    pub energy: f64,
    /// The force on every atom.
    pub forces: Vec<Vec3>,
    /// `Ã2 = 4π/V Σ_k γ̂/k²`, the A2 of this set of wave vectors (NA2CLC = 1).
    pub a2_tilde: f64,
    /// The virial `Ξ` of the term (GROMOS sign, before the storage factor −2).
    pub virial: Sym3,
    /// `4π/V · (Σ_k (a|k|γ̂' − 2γ̂)/k⁴ k⊗k + Σ_k γ̂/k² 1)`, the derivative of `Ã2` its virial uses.
    pub a2_tilde_derivative: Sym3,
}

/// The reciprocal-space energy, forces, `Ã2` and virial of the direct Ewald sum, gromosXX
/// `_ls_ewald_kspace_outerloop` (`nonbonded_outerloop.cc:1946-2059`) in its loop and operation
/// order, for positions already put into the box and `eps0_inv = 1/ε₀`:
///
/// ```text
/// per k (in order):
///   C = Σ_i q_i cos(k·r_i),  S = Σ_i q_i sin(k·r_i)          (atoms in order, k·r = k_x r_x + …)
///   f_i += k · (k2i_gamma_hat · (C sin(k·r_i) − S cos(k·r_i)))
///   E   += k2i_gamma_hat · (C·C + S·S);   Ã2 += k2i_gamma_hat
///   T    = ((ak_gamma_hat_prime − 2 γ̂) · (k2i · k2i)) k ⊗ k
///   D += T;   Ξ += (C·C + S·S) T;   Ξ += (C·C + S·S) k2i_gamma_hat 1
/// f_i ·= q_i · (eps0_inv / V);   E ·= (eps0_inv / V) · ½
/// D += Ã2 1;   Ã2 ·= 4π / V;   Ξ ·= −¼ (eps0_inv / V);   dÃ2 = (4π / V) D
/// ```
///
/// GROMOS Vol. 2: `e_gamma_ewald`, `ewald_reci_force`, `vir_ew`, `a_2_tilde_ew`.
///
/// spec: ewald.kspace
pub fn ewald_kspace(
    positions: &[Vec3],
    charges: &[f64],
    kvecs: &[KVector],
    volume: f64,
    eps0_inv: f64,
) -> KSpace {
    let n = positions.len();
    let eps_volume_i = eps0_inv / volume;
    let mut energy = 0.0;
    let mut forces = vec![[0.0f64; 3]; n];
    let mut sin_kr = vec![0.0f64; n];
    let mut cos_kr = vec![0.0f64; n];
    let mut a2_tilde = 0.0;
    let mut virial: Sym3 = [0.0; 6];
    let mut sum_gammahat: Sym3 = [0.0; 6];

    for kv in kvecs {
        let k = kv.k;
        let mut c_k = 0.0;
        let mut s_k = 0.0;
        for (i, r) in positions.iter().enumerate() {
            let r_dot_k = k[0] * r.x + k[1] * r.y + k[2] * r.z;
            sin_kr[i] = r_dot_k.sin();
            cos_kr[i] = r_dot_k.cos();
            c_k += charges[i] * cos_kr[i];
            s_k += charges[i] * sin_kr[i];
        }
        for (i, f) in forces.iter_mut().enumerate() {
            let s = kv.k2i_gamma_hat * (c_k * sin_kr[i] - s_k * cos_kr[i]);
            f[0] += k[0] * s;
            f[1] += k[1] * s;
            f[2] += k[2] * s;
        }
        let ewald_factor = c_k * c_k + s_k * s_k;
        energy += kv.k2i_gamma_hat * ewald_factor;
        a2_tilde += kv.k2i_gamma_hat;
        let iso = (kv.ak_gamma_hat_prime - 2.0 * kv.gamma_hat) * (kv.k2i * kv.k2i);
        let term = sym_tensor([iso * k[0], iso * k[1], iso * k[2]], k);
        let diagonal = ewald_factor * kv.k2i_gamma_hat;
        for j in 0..6 {
            sum_gammahat[j] += term[j];
            virial[j] += ewald_factor * term[j];
        }
        virial[0] += diagonal;
        virial[3] += diagonal;
        virial[5] += diagonal;
    }

    let forces = forces
        .iter()
        .zip(charges)
        .map(|(f, &q)| {
            let s = q * eps_volume_i;
            Vec3::new(f[0] * s, f[1] * s, f[2] * s)
        })
        .collect();
    energy *= eps_volume_i * 0.5;
    sum_gammahat[0] += a2_tilde;
    sum_gammahat[3] += a2_tilde;
    sum_gammahat[5] += a2_tilde;
    let four_pi_v = 4.0 * std::f64::consts::PI / volume;
    a2_tilde *= four_pi_v;
    let scale = -0.25 * eps_volume_i;
    KSpace {
        energy,
        forces,
        a2_tilde,
        virial: virial.map(|x| x * scale),
        a2_tilde_derivative: sum_gammahat.map(|x| four_pi_v * x),
    }
}

/// The lattice-sum settings of a run (NONBONDED block, GROMOS Vol. 4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatticeSumParams {
    /// NSHAPE
    pub shape: ChargeShape,
    /// ASHAPE (nm)
    pub width: f64,
    /// NKX, NKY, NKZ
    pub nk: [i64; 3],
    /// KCUT, compared with `|k|²` (see [`k_vectors`])
    pub kcut: f64,
    /// NA2CLC: 0 no A2, 1 `Ã2` from the wave vectors, 2 A2 numerically, 3 both
    pub na2clc: i32,
    /// TOLA2
    pub tola2: f64,
    /// EPSLS (0: conducting boundary)
    pub eps_ls: f64,
    /// `1/ε₀` (gromosXX `math::eps0_i`)
    pub eps0_inv: f64,
    /// `1/(4πε₀)` (gromosXX `math::four_pi_eps_i`)
    pub four_pi_eps_i: f64,
}

/// The shells [`a2_shells`](super::lattice_sum::a2_shells) may sum before the run is refused
/// (gromosXX has no bound).
pub const A2_MAX_SHELLS: u32 = 300;

/// The energies of the lattice-sum terms of one step.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct LsEnergies {
    pub realspace: f64,
    pub kspace: f64,
    pub a_term: f64,
    pub self_term: f64,
    pub surface: f64,
}

impl LsEnergies {
    /// Real space + k space + A term (gromosXX `ls_pair_total`).
    pub fn pair(&self) -> f64 {
        self.realspace + self.kspace + self.a_term
    }

    /// Pair + self + surface (gromosXX `ls_total`).
    pub fn total(&self) -> f64 {
        self.pair() + self.self_term + self.surface
    }
}

/// What depends on the box only: the wave vectors and `A1`, `A2`, `A3`.
#[derive(Debug, Clone)]
struct BoxTerms {
    edges: [[f64; 3]; 3],
    volume: f64,
    kvecs: Vec<KVector>,
    a1: f64,
    a3: f64,
    /// The numerical A2 and its derivative (NA2CLC 2, 3), zero otherwise.
    a2: f64,
    d_a2: Sym3,
}

/// The reciprocal-space, self, A and surface terms of the direct Ewald sum for one run: the
/// wave vectors and the box constants are rebuilt whenever the box changes (gromosXX: at step 0
/// and every step under pressure scaling, `nonbonded_set.cc:168-171, 207-210`).
#[derive(Debug, Clone)]
pub struct EwaldSum {
    pub params: LatticeSumParams,
    box_terms: Option<BoxTerms>,
}

impl EwaldSum {
    pub fn new(params: LatticeSumParams) -> Self {
        Self {
            params,
            box_terms: None,
        }
    }

    /// The number of wave vectors of the current box (0 before the first evaluation).
    pub fn num_k_vectors(&self) -> usize {
        self.box_terms.as_ref().map_or(0, |b| b.kvecs.len())
    }

    fn box_terms(&mut self, edges: [[f64; 3]; 3], rectangular: bool) -> Result<&BoxTerms, String> {
        use super::lattice_sum::{
            a2_cubic, a2_shells, reciprocal_box, self_constants, WIGNER_CUBE,
        };
        let stale = self.box_terms.as_ref().is_none_or(|b| b.edges != edges);
        if stale {
            let p = self.params;
            let (volume, recip) = reciprocal_box(&edges, rectangular);
            let kvecs = k_vectors(p.nk, p.kcut, &recip, p.shape, p.width);
            if kvecs.is_empty() {
                return Err(
                    "calculate k space: kspace is empty. Please increase cutoffs \
                            (NONBONDED NKX, NKY, NKZ, KCUT)."
                        .to_string(),
                );
            }
            let (a1, a3) = self_constants(p.shape, p.width, volume);
            let (mut a2, mut d_a2) = (0.0, [0.0; 6]);
            if matches!(p.na2clc, 2 | 3) {
                // gromosXX `ls_self_outerloop:2370-2373`: a rectangular box with equal edges
                let cube = rectangular && edges[0][0] == edges[1][1] && edges[0][0] == edges[2][2];
                if cube {
                    let (value, d) = a2_cubic(WIGNER_CUBE, edges[0][0], a1, a3);
                    a2 = value;
                    d_a2 = [d, 0.0, 0.0, d, 0.0, d];
                } else {
                    let s = a2_shells(p.shape, p.width, &recip, volume, p.tola2, A2_MAX_SHELLS);
                    if !s.converged {
                        return Err(format!(
                            "lattice sum: A2 did not reach the relative tolerance TOLA2 = {} in \
                             {} shells of wave vectors; increase TOLA2 or use NA2CLC = 1",
                            p.tola2, s.shells
                        ));
                    }
                    a2 = s.a2;
                    d_a2 = s.derivative;
                }
            }
            self.box_terms = Some(BoxTerms {
                edges,
                volume,
                kvecs,
                a1,
                a3,
                a2,
                d_a2,
            });
        }
        Ok(self.box_terms.as_ref().expect("set above"))
    }

    /// Add the reciprocal-space, self, A and surface terms to `storage` (forces, and the virial
    /// as `W = −2 Ξ` when `virial`) and return their energies (`realspace` is the caller's).
    ///
    /// `positions` are the configuration's (unwrapped) positions; the k-space term puts them
    /// into the box first, as gromosXX does (`periodicity.put_into_positive_box`,
    /// `nonbonded_outerloop.cc:1975-1979`, `periodicity.cc:45-54`).
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
        use super::lattice_sum::{charge_sums, self_and_a_terms, self_and_a_virials, surface_term};
        let p = self.params;
        let terms = self.box_terms(edges, rectangular)?.clone();

        // the box centre, `(a + b + c) / 2`
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
        let ks = ewald_kspace(&in_box, charges, &terms.kvecs, terms.volume, p.eps0_inv);

        // NA2CLC (gromosXX `ls_self_outerloop:2558-2585`)
        let (a2, a2_tilde, d_a2, d_a2_tilde) = match p.na2clc {
            0 => (0.0, 0.0, [0.0; 6], [0.0; 6]),
            1 => (
                ks.a2_tilde,
                ks.a2_tilde,
                ks.a2_tilde_derivative,
                ks.a2_tilde_derivative,
            ),
            2 => (terms.a2, terms.a2, terms.d_a2, terms.d_a2),
            _ => (terms.a2, ks.a2_tilde, terms.d_a2, ks.a2_tilde_derivative),
        };
        let (sum_q, sum_q2) = charge_sums(charges);
        let (self_term, a_term) =
            self_and_a_terms(terms.a1, a2, terms.a3, a2_tilde, sum_q, sum_q2, p.eps0_inv);

        // the surface term: unwrapped positions, centre `a/2 + b/2 + c/2`
        let centre = Vec3::new(
            a[0] / 2.0 + b[0] / 2.0 + c[0] / 2.0,
            a[1] / 2.0 + b[1] / 2.0 + c[1] / 2.0,
            a[2] / 2.0 + b[2] / 2.0 + c[2] / 2.0,
        );
        let srf = surface_term(
            positions,
            charges,
            centre,
            p.eps_ls,
            terms.volume,
            p.eps0_inv,
        );

        for (i, f) in storage.forces.iter_mut().enumerate() {
            *f += ks.forces[i];
        }
        for (i, f) in storage.forces.iter_mut().enumerate() {
            *f += srf.forces[i];
        }
        if virial {
            let (self_virial, constant_virial) =
                self_and_a_virials(terms.a1, sum_q, sum_q2, &d_a2, &d_a2_tilde, p.four_pi_eps_i);
            for xi in [&ks.virial, &self_virial, &constant_virial, &srf.virial] {
                let m = super::lattice_sum::sym_to_matrix(xi);
                for (row, m_row) in storage.virial.iter_mut().zip(m) {
                    for (v, x) in row.iter_mut().zip(m_row) {
                        *v += -2.0 * x;
                    }
                }
            }
        }
        Ok(LsEnergies {
            realspace: 0.0,
            kspace: ks.energy,
            a_term,
            self_term,
            surface: srf.energy,
        })
    }
}
