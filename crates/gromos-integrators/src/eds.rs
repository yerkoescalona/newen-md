//! Enveloping distribution sampling: the reference-state Hamiltonian (gromosXX `algorithm::EDS`,
//! `algorithm/integration/eds.cc`, non-accelerated forms).
//!
//! The force field leaves the energy V_i, forces and virial of every end state in the
//! configuration (`Energy::eds_vi`, `Configuration::eds`); this algorithm, placed right after the
//! force field, computes the reference energy V_R and the state weights p_i = ∂V_R/∂V_i, adds
//! Σ p_i F_i to the forces and Σ p_i W_i to the virial, and books V_R as `eds_vr` (part of the
//! special energy, not of the potential energy). The integrator then sees the mixed forces.
//!
//! The three forms are the documented ones (`simulation/parameter.h:624-658`):
//!
//! - FORM 1 (single s): V_R = −1/(βs) · ln Σ_i exp(−βs(V_i − E_i^R));
//! - FORM 2 (multi s): V_R = −1/β · ln[ 1/(N−1) · Σ_{i<j} (e^{−βs_ij(V_i−E_i)} + e^{−βs_ij(V_j−E_j)})^{1/s_ij} ];
//! - FORM 3 (pairs): the same over the N−1 listed pairs, normalised by N/(2(N−1)).
//!
//! FORM 1 reproduces gromosXX 247c95c4 bit for bit (its running pairwise fold, in state order).
//! For FORMS 2 and 3 gromosXX's released weights are not the derivative of its V_R (zero-seeded
//! log-sum-exp, a missing `− max`, multi-s skipping the pairs (0, j ≥ 2), item 1.8); the weights
//! here are the documented ones, what gromosXX computes with `D2_eds_weights.patch`.

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::math::{Mat3, Vec3};
use gromos_core::topology::Topology;

/// The functional form of the reference Hamiltonian, with its smoothness parameters.
#[derive(Debug, Clone, PartialEq)]
pub enum EdsForm {
    /// FORM 1: one s for every state.
    SingleS { s: f64 },
    /// FORM 2: one s per state pair, upper-triangle row-major: (0,1), (0,2) … (1,2) …
    MultiS { s: Vec<f64> },
    /// FORM 3: N−1 pairs of states (0-based) with their s.
    PairS { pairs: Vec<(usize, usize, f64)> },
}

/// gromosXX's pairwise log-sum-exp step `max(a, b) + log(1 + exp(min(a, b) − max(a, b)))`
/// (`eds.cc:409-424`), started from −∞ (the patched multi-s/pair-s seed).
#[inline]
fn fold(acc: f64, x: f64) -> f64 {
    if acc == f64::NEG_INFINITY {
        // max(−∞, x) + log(1 + exp(−∞)) = x exactly
        return x;
    }
    let hi = if acc < x { x } else { acc };
    let lo = if x < acc { x } else { acc };
    hi + (1.0 + (lo - hi).exp()).ln()
}

/// Two states i, j of a pair term: `(LSE / s, A + LSE(1/s − 1), B + LSE(1/s − 1))` with
/// A = −βs(V_i − E_i), B = −βs(V_j − E_j) and LSE = ln(e^A + e^B) — the pair's contribution
/// to ln Σ and the two states' log-weights (`eds.cc:536-551`).
#[inline]
fn pair_term(beta: f64, s: f64, vi: &[f64], eir: &[f64], i: usize, j: usize) -> (f64, f64, f64) {
    let part_a = -beta * s * (vi[i] - eir[i]);
    let part_b = -beta * s * (vi[j] - eir[j]);
    let lse = fold(part_a, part_b);
    let elem = lse / s;
    let elem2 = lse * (1.0 / s - 1.0);
    (elem, part_a + elem2, part_b + elem2)
}

/// The reference energy V_R of end-state energies `vi` with offsets `eir` at β = 1/(k_B T),
/// and the weight p_i = ∂V_R/∂V_i of every state in `weights` (Σ p_i = 1).
pub fn reference_hamiltonian(
    form: &EdsForm,
    beta: f64,
    eir: &[f64],
    vi: &[f64],
    weights: &mut [f64],
) -> f64 {
    let n = vi.len();
    let (vr, sum) = match form {
        EdsForm::SingleS { s } => {
            // the running fold in state order, from the pair (0, 1) on (`eds.cc:409-424`)
            let mut sum = f64::NEG_INFINITY;
            for st in 0..n {
                let part = -beta * s * (vi[st] - eir[st]);
                weights[st] = part;
                sum = fold(sum, part);
            }
            (-1.0 / (beta * s) * sum, sum)
        },
        EdsForm::MultiS { s } => {
            weights.fill(f64::NEG_INFINITY);
            let mut sum = f64::NEG_INFINITY;
            for i in 0..n.saturating_sub(1) {
                for j in i + 1..n {
                    let k = (i * (2 * n - i - 1)) / 2 + j - i - 1;
                    let (elem, wi, wj) = pair_term(beta, s[k], vi, eir, i, j);
                    sum = fold(sum, elem);
                    weights[i] = fold(weights[i], wi);
                    weights[j] = fold(weights[j], wj);
                }
            }
            (-1.0 / beta * (sum + (1.0 / (n - 1) as f64).ln()), sum)
        },
        EdsForm::PairS { pairs } => {
            weights.fill(f64::NEG_INFINITY);
            let mut sum = f64::NEG_INFINITY;
            for &(i, j, s) in pairs {
                let (elem, wi, wj) = pair_term(beta, s, vi, eir, i, j);
                sum = fold(sum, elem);
                weights[i] = fold(weights[i], wi);
                weights[j] = fold(weights[j], wj);
            }
            let norm = (n as f64 / ((n - 1) * 2) as f64).ln();
            (-1.0 / beta * (sum + norm), sum)
        },
    };
    for w in weights.iter_mut() {
        *w = (*w - sum).exp();
    }
    vr
}

/// The EDS step: mixes the end states into the forces, virial and energies of the current
/// state (gromosXX `algorithm::EDS::apply`, `eds.cc:380-716`).
#[derive(Debug, Clone)]
pub struct Eds {
    pub form: EdsForm,
    /// E_i^R, one per state.
    pub eir: Vec<f64>,
    /// 1/(k_B T), T of the first bath (gromosXX `multibath.bath(0).temperature`), k_B from the
    /// topology's PHYSICALCONSTANTS.
    pub beta: f64,
    weights: Vec<f64>,
}

impl Eds {
    pub fn new(form: EdsForm, eir: Vec<f64>, beta: f64) -> Self {
        let n = eir.len();
        Self {
            form,
            eir,
            beta,
            weights: vec![0.0; n],
        }
    }

    pub fn num_states(&self) -> usize {
        self.eir.len()
    }
}

impl Algorithm for Eds {
    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        let n = self.num_states();
        if conf.eds.num_states() != n || conf.current().energies.eds_vi.len() != n {
            return Err(format!(
                "EDS with {n} states, but the force field left {} end states",
                conf.eds.num_states()
            ));
        }
        let vr = reference_hamiltonian(
            &self.form,
            self.beta,
            &self.eir,
            &conf.current().energies.eds_vi,
            &mut self.weights,
        );
        let (state, eds) = conf.current_and_eds_mut();
        let mut d_vr = 0.0;
        for (s, &p) in self.weights.iter().enumerate() {
            for (f, fs) in state.force.iter_mut().zip(&eds.forces[s]) {
                *f += *fs * p;
            }
            d_vr += p * state.perturbed_energy_derivatives.eds_vi[s];
            let v = &eds.virials[s];
            let w = &mut state.virial_tensor;
            *w = Mat3::from_cols(
                w.x_axis + Vec3::new(v[0][0], v[1][0], v[2][0]) * p,
                w.y_axis + Vec3::new(v[0][1], v[1][1], v[2][1]) * p,
                w.z_axis + Vec3::new(v[0][2], v[1][2], v[2][2]) * p,
            );
        }
        state.energies.eds_vr = vr;
        state.energies.special_total += vr;
        // both records' special totals include eds_vr (gromosXX `energy.cc:592-596`)
        state.perturbed_energy_derivatives.eds_vr = d_vr;
        state.perturbed_energy_derivatives.special_total += d_vr;
        Ok(())
    }

    fn name(&self) -> &str {
        "EDS"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BETA: f64 = 1.0 / (0.00831441 * 300.0);

    /// −1/β · ln Σ_k exp(−β V_eff,k) by brute force over the documented formula's terms.
    fn lse(terms: &[f64]) -> f64 {
        let m = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        m + terms.iter().map(|t| (t - m).exp()).sum::<f64>().ln()
    }

    /// The single-s fold is gromosXX's: pairwise, in state order, never a global max shift.
    #[test]
    fn single_s_is_gromosxx_running_fold_in_state_order() {
        let vi = [15.39, -16.08, 3.5];
        let eir = [0.0, 10.0, -5.0];
        let s = 0.3;
        let parts: Vec<f64> = (0..3).map(|k| -BETA * s * (vi[k] - eir[k])).collect();
        let mut sum = parts[0].max(parts[1])
            + (1.0 + (parts[0].min(parts[1]) - parts[0].max(parts[1])).exp()).ln();
        sum = sum.max(parts[2]) + (1.0 + (sum.min(parts[2]) - sum.max(parts[2])).exp()).ln();
        let mut w = [0.0; 3];
        let vr = reference_hamiltonian(&EdsForm::SingleS { s }, BETA, &eir, &vi, &mut w);
        assert_eq!(vr.to_bits(), (-1.0 / (BETA * s) * sum).to_bits());
        for k in 0..3 {
            assert_eq!(w[k].to_bits(), (parts[k] - sum).exp().to_bits());
        }
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    }

    /// FORMS 2 and 3 give the documented V_R, and their weights are ∂V_R/∂V_i (summing to 1).
    #[test]
    fn multi_s_and_pair_s_follow_the_documented_formulas() {
        let vi = [12.0, -20.0, 4.0, -1.0];
        let eir = [0.0, 10.0, -5.0, 2.0];
        let n = vi.len();
        let s_pairs = [0.3, 0.2, 0.1, 0.25, 0.4, 0.15];
        let multi = EdsForm::MultiS {
            s: s_pairs.to_vec(),
        };
        let pairs = EdsForm::PairS {
            pairs: vec![(0, 1, 0.3), (1, 2, 0.1), (3, 1, 0.2)],
        };
        let documented = |form: &EdsForm, vi: &[f64]| -> f64 {
            let term = |i: usize, j: usize, s: f64| {
                lse(&[-BETA * s * (vi[i] - eir[i]), -BETA * s * (vi[j] - eir[j])]) / s
            };
            match form {
                EdsForm::MultiS { s } => {
                    let mut terms = Vec::new();
                    let mut k = 0;
                    for i in 0..n {
                        for j in i + 1..n {
                            terms.push(term(i, j, s[k]));
                            k += 1;
                        }
                    }
                    -1.0 / BETA * (lse(&terms) + (1.0 / (n - 1) as f64).ln())
                },
                EdsForm::PairS { pairs } => {
                    let terms: Vec<f64> = pairs.iter().map(|&(i, j, s)| term(i, j, s)).collect();
                    -1.0 / BETA * (lse(&terms) + (n as f64 / (2 * (n - 1)) as f64).ln())
                },
                EdsForm::SingleS { .. } => unreachable!(),
            }
        };
        for form in [multi, pairs] {
            let mut w = [0.0; 4];
            let vr = reference_hamiltonian(&form, BETA, &eir, &vi, &mut w);
            assert!((vr - documented(&form, &vi)).abs() < 1e-10, "{form:?}");
            assert!(
                (w.iter().sum::<f64>() - 1.0).abs() < 1e-14,
                "{form:?}: {w:?}"
            );
            for k in 0..n {
                let h = 1e-5;
                let (mut up, mut dn) = (vi, vi);
                up[k] += h;
                dn[k] -= h;
                let dvr = (documented(&form, &up) - documented(&form, &dn)) / (2.0 * h);
                assert!(
                    (w[k] - dvr).abs() < 1e-7,
                    "{form:?} state {k}: {} vs {dvr}",
                    w[k]
                );
            }
        }
    }

    /// A state no listed pair names has no weight.
    #[test]
    fn pair_s_state_outside_every_pair_has_zero_weight() {
        let mut w = [0.0; 3];
        let form = EdsForm::PairS {
            pairs: vec![(0, 1, 0.5), (1, 0, 0.5)],
        };
        reference_hamiltonian(&form, BETA, &[0.0; 3], &[1.0, 2.0, 3.0], &mut w);
        assert_eq!(w[2], 0.0);
    }
}
