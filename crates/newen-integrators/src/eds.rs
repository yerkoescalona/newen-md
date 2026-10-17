//! Enveloping distribution sampling: the reference-state Hamiltonian (gromosXX `algorithm::EDS`,
//! `algorithm/integration/eds.cc`), plain and accelerated.
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
//!
//! Accelerated EDS (the `AEDS` block, `eds.cc:57-378`) boosts the single-s (s = 1) mixed energy
//! V_mix between E_min and E_max: V_R* = V_mix − (V_mix − E_min)²/(2(E_max − E_min)) inside the
//! window, V_mix − (E_max − E_min)/2 above it, V_mix below it; the EDS part of the forces, virial
//! and dV/dλ is scaled by 1 − (V_mix − E_min)/(E_max − E_min) inside the window. FORMS 2-6 then
//! search E_max, E_min and the offsets E_i^R every step ([`search_step`]); the search state lives
//! in the configuration ([`AedsSearch`], the `AEDSSEARCH` block) so a run continues it. Three
//! gromosXX defects are corrected by default and reproduced on request (item 1.8): the energy
//! columns written after the update ([`AedsColumns`], D16), the FORM 5/6 prevalences as plain
//! exponentials ([`AedsPrevalence`], D25) and a never-visited state 1 as the E_min target
//! ([`AedsTarget`], D32).

use newen_core::algorithm::{Algorithm, SimulationState};
pub use newen_core::configuration::AedsSearch;
use newen_core::configuration::Configuration;
use newen_core::math::{Mat3, Vec3};
use newen_core::topology::Topology;

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

/// The forms of accelerated EDS (the `AEDS` block's FORM, gromosXX `simulation/parameter.h:659-683`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AedsForm {
    /// FORM 1: E_max, E_min and E^R fixed.
    Fixed,
    /// FORM 2: search the offsets E^R (memory-averaged free energies).
    SearchEir,
    /// FORM 3: search E_max and E_min.
    SearchEmaxEmin,
    /// FORM 4: search E_max, E_min and E^R.
    SearchAll,
    /// FORM 5: search E_max and E_min; offsets by the prevalence-driven ("advanced") update.
    AdvancedAll,
    /// FORM 6: the advanced offset update with E_max and E_min fixed.
    AdvancedEir,
}

impl AedsForm {
    /// The form of the block's FORM value (1-6).
    pub fn from_gromos(form: i32) -> Option<Self> {
        Some(match form {
            1 => Self::Fixed,
            2 => Self::SearchEir,
            3 => Self::SearchEmaxEmin,
            4 => Self::SearchAll,
            5 => Self::AdvancedAll,
            6 => Self::AdvancedEir,
            _ => return None,
        })
    }

    /// Any parameter search (FORMS 2-6): these keep and write the `AEDSSEARCH` state.
    pub fn searches(self) -> bool {
        self != Self::Fixed
    }

    /// FORMS 2, 4, 5, 6: the memory-averaged end-state free energies are updated.
    pub fn searches_offsets(self) -> bool {
        !matches!(self, Self::Fixed | Self::SearchEmaxEmin)
    }

    /// FORMS 2 and 4: the offsets are set to those free energies.
    pub fn conventional_offsets(self) -> bool {
        matches!(self, Self::SearchEir | Self::SearchAll)
    }

    /// FORMS 3, 4, 5: E_max and E_min are searched.
    pub fn searches_emax_emin(self) -> bool {
        matches!(
            self,
            Self::SearchEmaxEmin | Self::SearchAll | Self::AdvancedAll
        )
    }

    /// FORMS 5 and 6: the prevalence-driven offset update.
    pub fn advanced(self) -> bool {
        matches!(self, Self::AdvancedAll | Self::AdvancedEir)
    }
}

/// The E_max, E_min and offsets an A-EDS energy frame carries (item 1.8 D16).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(rename_all = "kebab-case")
)]
pub enum AedsColumns {
    /// The ones this step's V_mix and V_R* were computed with, so a frame's V_R* follows from
    /// its own e_i, e_ir, E_max and E_min.
    #[default]
    Intended,
    /// gromosXX 247c95c4: the ones after this step's search update (`eds.cc:368-372`), which
    /// the next step uses.
    GromosxxPostUpdate,
}

/// How FORMS 5 and 6 compute the prevalences p_i = e^{x_i}/Σ_j e^{x_j}, x_i = −β(V_i − E_i^R)
/// (item 1.8 D25).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(rename_all = "kebab-case")
)]
pub enum AedsPrevalence {
    /// Every exponent shifted by the largest one (log-sum-exp): finite for any energies.
    #[default]
    Intended,
    /// gromosXX 247c95c4: plain exponentials (`eds.cc:339-343`); once an x_i exceeds ≈ 709.78
    /// (or all are below ≈ −745) the prevalences are NaN and no offset moves.
    GromosxxPlainExp,
}

/// Which state's mean and σ set E_min (G and σ_G, item 1.8 D32).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(rename_all = "kebab-case")
)]
pub enum AedsTarget {
    /// The lowest mean of V_i − E_i^R among the states visited more than once, state 1 among
    /// them; state 1 when none is.
    #[default]
    Intended,
    /// gromosXX 247c95c4 (`eds.cc:252-261`): state 1 is the starting candidate whatever its
    /// visit count, so a state 1 never visited (mean 0, σ 0) stays the target until a sampled
    /// state has a mean ≤ 0.
    GromosxxState1Default,
}

/// The fixed settings of an A-EDS search (the `AEDS` block after FORM).
#[derive(Debug, Clone, PartialEq)]
pub struct AedsParams {
    pub form: AedsForm,
    /// BMAXTYPE: 1 BMAX in kJ/mol, 2 BMAX in standard deviations of the target state's V_i.
    pub bmax_type: u32,
    /// BMAX
    pub bmax: f64,
    /// ASTEPS: memory time (steps) at the start of the run; in FORMS 5 and 6 also the offset
    /// step, ASTEPS·Δt.
    pub asteps: u32,
    /// BSTEPS: memory time (steps) at the end of the run.
    pub bsteps: u32,
    /// The parameters the energy frames carry.
    pub columns: AedsColumns,
    /// The FORM 5/6 prevalences.
    pub prevalence: AedsPrevalence,
    /// The state that sets G and σ_G.
    pub target: AedsTarget,
}

/// The accelerated reference energy V_R* of the mixed energy `v_mix` and the factor that scales
/// the EDS part of the forces (gromosXX `eds.cc:122-134`, V2 eq. 23301-23307). The factor is 1
/// outside the window (`x · 1.0` is exact, so the forces there are plain EDS's).
pub fn boosted_reference(v_mix: f64, emax: f64, emin: f64) -> (f64, f64) {
    if v_mix <= emin {
        (v_mix, 1.0)
    } else if v_mix >= emax {
        (v_mix - 0.5 * (emax - emin), 1.0)
    } else {
        let demix = v_mix - emin;
        let kfac = 1.0 / (emax - emin);
        (v_mix - 0.5 * kfac * demix * demix, 1.0 - kfac * demix)
    }
}

/// The memory time τ (steps) at step `step` of a run of `nstlim` steps: ASTEPS ramped linearly
/// to BSTEPS over the current run (`eds.cc:186`), so a search chained over several jobs restarts
/// the ramp in each. The difference is unsigned, as in gromosXX: with
/// BSTEPS < ASTEPS it wraps to about 4·10⁹ (the recipe refuses that input, item 1.8 D13).
pub fn memory_time(asteps: u32, bsteps: u32, step: usize, nstlim: usize) -> f64 {
    f64::from(asteps) + f64::from(bsteps.wrapping_sub(asteps)) * step as f64 / nstlim as f64
}

/// The first step of a fresh search (`eds.cc:107-120`): FORMS 3, 4, 5 start E_max, E_min and the
/// E_max candidate at V_mix; FORMS 2, 4, 5, 6 seed the averages with the input offsets.
pub fn initialise_search(form: AedsForm, s: &mut AedsSearch, v_mix: f64, beta: f64) {
    if form.searches_emax_emin() {
        s.emax = v_mix;
        s.emin = v_mix;
        s.searchemax = v_mix;
    }
    if form.searches_offsets() {
        for i in 0..s.num_states() {
            s.lnexpde[i] = -(s.eir[i] - s.eir[0]) * beta;
        }
    }
}

/// The boosted energy of one end state against its own window [E_min + E_i^R, E_max + E_i^R]
/// (`eds.cc:189-201`).
fn boosted_state(v: f64, emax: f64, emin: f64) -> f64 {
    if v <= emin {
        v
    } else if v >= emax {
        v - 0.5 * (emax - emin)
    } else {
        let demix = v - emin;
        let kfac = 1.0 / (emax - emin);
        v - 0.5 * kfac * demix * demix
    }
}

/// The memory-averaged free energies of the end states and, in FORMS 2 and 4, the offsets
/// (`eds.cc:184-221`): ⟨e^{−β(V_i* − V_R*)}⟩ relaxes with memory time `tau` (steps), folded in log
/// space; E_i^R = F_i − F_0 in FORMS 2 and 4 (FORMS 5 and 6 only keep F_i).
pub fn offset_update(
    form: AedsForm,
    s: &mut AedsSearch,
    vi: &[f64],
    v_r: f64,
    beta: f64,
    tau: f64,
) {
    let decay = (-1.0 / tau).exp();
    for i in 0..s.num_states() {
        let estar = boosted_state(vi[i], s.emax + s.eir[i], s.emin + s.eir[i]);
        let expde = -beta * (estar - v_r);
        let mut l = s.lnexpde[i] + (decay / (1.0 - decay)).ln();
        let (hi, lo) = if l < expde { (expde, l) } else { (l, expde) };
        l = hi + (1.0 + (lo - hi).exp()).ln();
        l += (1.0 - decay).ln();
        s.lnexpde[i] = l;
        s.statefren[i] = -1.0 / beta * l;
        if form.conventional_offsets() {
            s.eir[i] = s.statefren[i] - s.statefren[0];
        }
    }
}

/// What [`track_state`] found this step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateTrack {
    /// The current state: the lowest V_i − E_i^R.
    pub state: usize,
    /// All states have been visited since the last round trip (and the flags were reset).
    pub round_trip: bool,
    /// The state with the lowest mean of V_i − E_i^R among those visited more than once (with
    /// [`AedsTarget::GromosxxState1Default`] state 0 always a candidate); state 0 when none
    /// qualifies.
    pub target: usize,
    /// The target's mean of V_i − E_i^R (E_min on the first step of a fresh search).
    pub globminavg: f64,
    /// The target's standard deviation of V_i.
    pub globminfluc: f64,
}

/// State tracking, round trips and the running averages (`eds.cc:223-266`, `check_round_trip`
/// :731-750). The offsets are the ones just updated, as gromosXX reads them. The mean of
/// V_i − E_i^R averages over the offsets of every step it saw, as in gromosXX (`:246`).
pub fn track_state(s: &mut AedsSearch, vi: &[f64], rule: AedsTarget) -> StateTrack {
    let n = s.num_states();
    let mut state = 0;
    let mut min = vi[0] - s.eir[0];
    for i in 1..n {
        if vi[i] - s.eir[i] < min {
            min = vi[i] - s.eir[i];
            state = i;
        }
    }
    s.visited[state] = true;
    s.visitcounts[state] = s.visitcounts[state].wrapping_add(1);
    let round_trip = s.visited.iter().all(|&v| v);
    if round_trip {
        s.emaxcounts = s.emaxcounts.wrapping_add(1);
        s.visited.fill(false);
    }
    let count = f64::from(s.visitcounts[state]);
    let v = vi[state];
    let mean = s.avgenergy[state] + (v - s.avgenergy[state]) / count;
    s.eiravgenergy[state] += (v - s.eiravgenergy[state] - s.eir[state]) / count;
    if s.visitcounts[state] > 1 {
        s.bigs[state] += (v - mean) * (v - s.avgenergy[state]);
        s.stdevenergy[state] = (s.bigs[state] / f64::from(s.visitcounts[state] - 1)).sqrt();
    }
    s.avgenergy[state] = mean;
    let mut target = 0;
    let mut lowest = s.eiravgenergy[0];
    match rule {
        AedsTarget::Intended => {
            let mut sampled = false;
            for i in 0..n {
                if s.visitcounts[i] > 1 && (!sampled || lowest >= s.eiravgenergy[i]) {
                    target = i;
                    lowest = s.eiravgenergy[i];
                    sampled = true;
                }
            }
        },
        AedsTarget::GromosxxState1Default => {
            for i in 1..n {
                if lowest >= s.eiravgenergy[i] && s.visitcounts[i] > 1 {
                    target = i;
                    lowest = s.eiravgenergy[i];
                }
            }
        },
    }
    let globminavg = if s.initialise {
        // gromosXX: "prevent rounding error in the first simulation step" (`eds.cc:263-266`)
        s.emin
    } else {
        s.eiravgenergy[target]
    };
    StateTrack {
        state,
        round_trip,
        target,
        globminavg,
        globminfluc: s.stdevenergy[target],
    }
}

/// The E_max/E_min search of FORMS 3, 4, 5 (`eds.cc:270-330`): E_max is the highest V_mix at a
/// state change until the first round trip, then the mean of the per-round-trip maxima; E_min
/// puts the boosted transition energy at G + b_max (G = the target's mean), quadratic when that
/// falls below G, held at G until the first round trip with RESTREMIN 1.
///
/// gromosXX divides by b_max at `eds.cc:314`; with BMAXTYPE 2 and a target without variance
/// (b_max = 0) E_min becomes −∞ and V_R* −∞, then NaN (item 1.8 D15). That case is an error here
/// whenever the RESTREMIN clamp does not catch it.
pub fn emax_emin_update(
    p: &AedsParams,
    s: &mut AedsSearch,
    v_mix: f64,
    t: &StateTrack,
) -> Result<(), String> {
    if t.state as u32 != s.oldstate && s.searchemax < v_mix {
        s.searchemax = v_mix;
    }
    if s.emaxcounts == 0 {
        s.emax = s.searchemax;
    }
    if t.round_trip {
        if s.emaxcounts == 1 {
            s.emax = 0.0;
            s.fullemin = true;
        }
        s.emax += (s.searchemax - s.emax) / f64::from(s.emaxcounts);
        s.searchemax = t.globminavg;
    }
    let g = t.globminavg;
    let bmax = match p.bmax_type {
        1 => p.bmax,
        2 => s.stdevenergy[t.target] * p.bmax,
        _ => 0.0,
    };
    if s.emax - g <= bmax {
        s.emin = s.emax;
    } else {
        s.emin = 2.0 * (g + bmax) - s.emax;
        if s.emin < g {
            s.emin =
                (-s.emax * s.emax + 2.0 * s.emax * bmax + 2.0 * s.emax * g - g * g) / (2.0 * bmax);
        }
        if !s.fullemin && s.emin < g {
            s.emin = g;
            if s.emin > s.emax {
                s.emin = s.emax;
            }
        }
    }
    if !s.emin.is_finite() {
        return Err(format!(
            "AEDS: E_min = {} (b_max = {bmax}: the target state {} has no variance yet, and \
             gromosXX divides by b_max, eds.cc:314, then runs into NaN energies); use BMAXTYPE 1 \
             or RESTREMIN 1",
            s.emin,
            t.target + 1
        ));
    }
    s.oldstate = t.state as u32;
    Ok(())
}

/// The prevalence-driven offset update of FORMS 5 and 6 (`eds.cc:332-363`): every state holding
/// at least 1/N of the weight this step has its offset lowered by `step_size`·p_i/(n_rt + 1)
/// with n_rt the round trips so far, then the offsets are centred on zero. `step_size` is
/// ASTEPS·Δt, as gromosXX has it (ASTEPS read as a rate per ps; item 1.8, 0.0.60's D17). The
/// prevalences are a log-sum-exp ([`AedsPrevalence::Intended`]) or gromosXX's plain
/// exponentials (`eds.cc:339-343`), which overflow once some V_i − E_i^R is below about −709/β.
pub fn advanced_offset_update(
    s: &mut AedsSearch,
    vi: &[f64],
    v_r: f64,
    beta: f64,
    step_size: f64,
    prevalence_rule: AedsPrevalence,
    weights: &mut [f64],
) {
    let n = s.num_states();
    let shift = match prevalence_rule {
        AedsPrevalence::Intended => (0..n)
            .map(|i| -beta * (vi[i] - s.eir[i]))
            .fold(f64::NEG_INFINITY, f64::max),
        AedsPrevalence::GromosxxPlainExp => 0.0,
    };
    let mut total = 0.0;
    for i in 0..n {
        weights[i] = match prevalence_rule {
            AedsPrevalence::Intended => (-beta * (vi[i] - s.eir[i]) - shift).exp(),
            AedsPrevalence::GromosxxPlainExp => (-beta * (vi[i] - s.eir[i])).exp(),
        };
        total += weights[i];
    }
    let damping = f64::from(s.emaxcounts.wrapping_add(1));
    for i in 0..n {
        let prevalence = weights[i] / total;
        if prevalence >= 1.0 / n as f64 {
            if vi[i] - v_r <= s.statefren[i] + 1.0 / beta {
                s.framecounts[i] = s.framecounts[i].wrapping_add(1);
            }
            s.eir[i] -= step_size * (prevalence / damping);
        }
    }
    let mean = s.eir.iter().fold(0.0, |a, &e| a + e) / n as f64;
    for e in s.eir.iter_mut() {
        *e -= mean;
    }
}

/// One step of the A-EDS parameter search on this step's end-state energies `vi`, mixed energy
/// `v_mix` and accelerated reference `v_r` (`eds.cc:178-366`), at step `step` of a run of
/// `nstlim` steps of `dt`. Returns `(G, σ_G)` in the forms that search E_max/E_min (the energy
/// trajectory's `eds_globmin`, `eds_globminfluc`).
#[allow(clippy::too_many_arguments)]
pub fn search_step(
    p: &AedsParams,
    s: &mut AedsSearch,
    vi: &[f64],
    v_mix: f64,
    v_r: f64,
    beta: f64,
    step: usize,
    nstlim: usize,
    dt: f64,
    scratch: &mut [f64],
) -> Result<Option<(f64, f64)>, String> {
    if p.form.searches_offsets() {
        let tau = memory_time(p.asteps, p.bsteps, step, nstlim);
        offset_update(p.form, s, vi, v_r, beta, tau);
    }
    let t = track_state(s, vi, p.target);
    let mut globmin = None;
    if p.form.searches_emax_emin() {
        emax_emin_update(p, s, v_mix, &t)?;
        globmin = Some((t.globminavg, t.globminfluc));
    }
    if p.form.advanced() {
        let step_size = f64::from(p.asteps) * dt;
        advanced_offset_update(s, vi, v_r, beta, step_size, p.prevalence, scratch);
    }
    Ok(globmin)
}

/// The accelerated part of an [`Eds`] step./// The accelerated part of an [`Eds`] step.
#[derive(Debug, Clone)]
struct Accelerated {
    params: AedsParams,
    /// The state a fresh search (or FORM 1) starts from.
    fresh: AedsSearch,
    /// NTIAEDSS = 0 in a search form: continue the configuration's state instead.
    restart: bool,
}

/// The EDS step: mixes the end states into the forces, virial and energies of the current
/// state (gromosXX `algorithm::EDS::apply`, `eds.cc:380-716`).
#[derive(Debug, Clone)]
pub struct Eds {
    pub form: EdsForm,
    /// E_i^R, one per state (A-EDS: the starting offsets; the run's are in the configuration).
    pub eir: Vec<f64>,
    /// 1/(k_B T), T of the first bath (gromosXX `multibath.bath(0).temperature`), k_B from the
    /// topology's PHYSICALCONSTANTS.
    pub beta: f64,
    weights: Vec<f64>,
    accelerated: Option<Accelerated>,
    vi: Vec<f64>,
}

impl Eds {
    pub fn new(form: EdsForm, eir: Vec<f64>, beta: f64) -> Self {
        let n = eir.len();
        Self {
            form,
            eir,
            beta,
            weights: vec![0.0; n],
            accelerated: None,
            vi: vec![0.0; n],
        }
    }

    /// Accelerated EDS (the `AEDS` block): single s = 1, starting from `emax`, `emin` and `eir`
    /// with `fullemin` (RESTREMIN 0), or — `restart`, NTIAEDSS = 0 in a search form — from the
    /// configuration's `AEDSSEARCH` state.
    pub fn accelerated(
        params: AedsParams,
        emax: f64,
        emin: f64,
        eir: Vec<f64>,
        fullemin: bool,
        restart: bool,
        beta: f64,
    ) -> Self {
        let restart = restart && params.form.searches();
        let fresh = AedsSearch::new(emax, emin, eir.clone(), fullemin);
        let mut eds = Self::new(EdsForm::SingleS { s: 1.0 }, eir, beta);
        eds.accelerated = Some(Accelerated {
            params,
            fresh,
            restart,
        });
        eds
    }

    pub fn num_states(&self) -> usize {
        self.eir.len()
    }

    /// Plain EDS (`eds.cc:380-716`).
    fn apply_plain(&mut self, conf: &mut Configuration) {
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
    }

    /// Accelerated EDS (`eds.cc:57-378`) with the configuration's search state `s`.
    fn apply_accelerated(
        &mut self,
        conf: &mut Configuration,
        s: &mut AedsSearch,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let params = &self.accelerated.as_ref().expect("accelerated").params;
        self.vi.clone_from(&conf.current().energies.eds_vi);
        let v_mix =
            reference_hamiltonian(&self.form, self.beta, &s.eir, &self.vi, &mut self.weights);
        if s.initialise {
            initialise_search(params.form, s, v_mix, self.beta);
        }
        let (v_r, fkfac) = boosted_reference(v_mix, s.emax, s.emin);
        let (state, eds) = conf.current_and_eds_mut();
        let mut d_vr = 0.0;
        for (k, &p) in self.weights.iter().enumerate() {
            for (f, fs) in state.force.iter_mut().zip(&eds.forces[k]) {
                *f += *fs * p * fkfac;
            }
            d_vr += p * state.perturbed_energy_derivatives.eds_vi[k] * fkfac;
            let v = &eds.virials[k];
            let w = &mut state.virial_tensor;
            *w = Mat3::from_cols(
                w.x_axis + Vec3::new(v[0][0], v[1][0], v[2][0]) * p * fkfac,
                w.y_axis + Vec3::new(v[0][1], v[1][1], v[2][1]) * p * fkfac,
                w.z_axis + Vec3::new(v[0][2], v[1][2], v[2][2]) * p * fkfac,
            );
        }
        let e = &mut state.energies;
        e.eds_vmix = v_mix;
        e.eds_vr = v_r;
        e.special_total += v_r;
        if params.columns == AedsColumns::Intended {
            e.eds_emax = s.emax;
            e.eds_emin = s.emin;
            e.eds_eir.clone_from(&s.eir);
        }
        state.perturbed_energy_derivatives.eds_vr = d_vr;
        state.perturbed_energy_derivatives.special_total += d_vr;

        if params.form.searches() {
            let globmin = search_step(
                params,
                s,
                &self.vi,
                v_mix,
                v_r,
                self.beta,
                sim.step,
                sim.n_steps,
                sim.dt,
                &mut self.weights,
            )?;
            if let Some((g, fluc)) = globmin {
                let e = &mut conf.current_mut().energies;
                e.eds_globmin = g;
                e.eds_globminfluc = fluc;
            }
        }
        if params.columns == AedsColumns::GromosxxPostUpdate {
            // gromosXX writes them after the update (`eds.cc:368-372`): the ones the next step
            // uses (item 1.8 D16)
            let e = &mut conf.current_mut().energies;
            e.eds_emax = s.emax;
            e.eds_emin = s.emin;
            e.eds_eir.clone_from(&s.eir);
        }
        s.initialise = false;
        Ok(())
    }
}

impl Algorithm for Eds {
    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    fn init(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        let Some(a) = &self.accelerated else {
            return Ok(());
        };
        if !a.restart {
            conf.aeds = Some(a.fresh.clone());
            return Ok(());
        }
        // NTIAEDSS = 0: the configuration's AEDSSEARCH block overrides the input's E_max, E_min
        // and E^R (`in_configuration.cc:1399-1455`)
        match conf.aeds.as_mut() {
            None => Err("AEDS NTIAEDSS = 0: no AEDSSEARCH block in configuration".to_string()),
            Some(s) if s.num_states() != self.num_states() => Err(format!(
                "AEDS with {} states, but the AEDSSEARCH block has {}",
                self.num_states(),
                s.num_states()
            )),
            Some(s) => {
                s.initialise = false;
                Ok(())
            },
        }
    }

    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let n = self.num_states();
        if conf.eds.num_states() != n || conf.current().energies.eds_vi.len() != n {
            return Err(format!(
                "EDS with {n} states, but the force field left {} end states",
                conf.eds.num_states()
            ));
        }
        if self.accelerated.is_none() {
            self.apply_plain(conf);
            return Ok(());
        }
        let mut s = conf
            .aeds
            .take()
            .ok_or_else(|| "A-EDS without its search state (init not run)".to_string())?;
        let result = self.apply_accelerated(conf, &mut s, sim);
        conf.aeds = Some(s);
        result
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

    /// The P1 probe of gromosXX (aladip, e_1 = 15.39360093, e_2 = −16.07623145, EIR 0 10): V_mix
    /// below, inside and above the window and the printed V_R*.
    #[test]
    fn boosted_reference_reproduces_the_three_placements() {
        let beta = 1.0 / (0.00831441 * 300.0);
        let mut w = [0.0; 2];
        let v_mix = reference_hamiltonian(
            &EdsForm::SingleS { s: 1.0 },
            beta,
            &[0.0, 10.0],
            &[15.39360093, -16.07623145],
            &mut w,
        );
        assert!((v_mix - -26.0762316).abs() < 1e-7, "{v_mix}");
        let below = boosted_reference(v_mix, 50.0, -10.0);
        assert_eq!(below, (v_mix, 1.0));
        let above = boosted_reference(v_mix, -40.0, -80.0);
        assert_eq!(above, (v_mix - 20.0, 1.0));
        let (vr, k) = boosted_reference(v_mix, 40.0, -60.0);
        assert!((vr - -31.83034191).abs() < 1e-7, "{vr}");
        assert!((k - 0.660762316).abs() < 1e-8, "{k}");
        // C¹ at both borders, and the factor is dV_R*/dV_mix
        for (emax, emin) in [(40.0, -60.0), (-10.0, -30.0)] {
            for x in [emin, emax, 0.5 * (emin + emax)] {
                let h = 1e-6;
                let d = (boosted_reference(x + h, emax, emin).0
                    - boosted_reference(x - h, emax, emin).0)
                    / (2.0 * h);
                let mid = boosted_reference(x, emax, emin).1;
                let expect = if x <= emin || x >= emax { 1.0 } else { mid };
                assert!(
                    (d - expect).abs() < 1e-6 || x == emin || x == emax,
                    "{x}: {d}"
                );
            }
        }
    }

    /// τ ramps from ASTEPS to BSTEPS over the run; BSTEPS < ASTEPS wraps as gromosXX's unsigned
    /// difference does (refused by the recipe).
    #[test]
    fn memory_time_ramps_and_wraps_like_gromosxx() {
        assert_eq!(memory_time(20, 100, 0, 200), 20.0);
        assert_eq!(memory_time(20, 100, 100, 200), 60.0);
        assert_eq!(memory_time(50, 50, 7, 10), 50.0);
        let wrapped = memory_time(1000, 100, 1, 200);
        assert!((wrapped - (1000.0 + (u32::MAX as f64 - 899.0) / 200.0)).abs() < 1e-6);
    }

    fn search(n: usize, fullemin: bool) -> AedsSearch {
        let mut s = AedsSearch::new(0.0, 0.0, vec![0.0; n], fullemin);
        s.initialise = false;
        s
    }

    /// Round trips reset the visited flags; the target is the visited state with the lowest mean
    /// V_i − E_i^R, state 0 by default. Both rules agree here: state 0 is sampled.
    #[test]
    fn track_state_counts_round_trips_and_picks_the_target() {
        for rule in [AedsTarget::Intended, AedsTarget::GromosxxState1Default] {
            let mut s = search(2, true);
            let t = track_state(&mut s, &[5.0, 6.0], rule);
            assert_eq!((t.state, t.round_trip, t.target), (0, false, 0));
            let t = track_state(&mut s, &[7.0, 2.5], rule);
            assert_eq!((t.state, t.round_trip), (1, true));
            assert_eq!(s.emaxcounts, 1);
            assert_eq!(s.visited, vec![false, false]);
            // state 1 visited once only: not a target yet
            assert_eq!(t.target, 0);
            let t = track_state(&mut s, &[7.0, 1.5], rule);
            assert_eq!((t.state, t.target), (1, 1));
            assert_eq!(t.globminavg, 2.0);
            assert!((t.globminfluc - (0.5f64).sqrt()).abs() < 1e-15);
        }
    }

    /// D32: a state 0 never visited (mean 0, σ 0) is gromosXX's target while the sampled
    /// state's mean is above 0; the intended rule takes the sampled state, mean and σ by hand.
    /// With no state sampled twice both fall back to state 0.
    #[test]
    fn target_is_a_sampled_state_not_an_unvisited_state_1() {
        let energies = [[9.0, 6.0], [9.0, 8.0], [9.0, 4.0]];
        let mut ours = search(2, true);
        let mut gxx = search(2, true);
        let mut last = None;
        for (k, vi) in energies.iter().enumerate() {
            let a = track_state(&mut ours, vi, AedsTarget::Intended);
            let b = track_state(&mut gxx, vi, AedsTarget::GromosxxState1Default);
            assert_eq!((a.state, b.state), (1, 1));
            if k == 0 {
                // one visit of state 1: nothing qualifies, both take state 0
                assert_eq!((a.target, a.globminavg), (0, 0.0));
                assert_eq!((b.target, b.globminavg), (0, 0.0));
            }
            last = Some((a, b));
        }
        let (a, b) = last.unwrap();
        assert_eq!(ours.visitcounts, vec![0, 3]);
        // gromosXX: the unvisited state 0, G = 0, σ_G = 0
        assert_eq!((b.target, b.globminavg, b.globminfluc), (0, 0.0, 0.0));
        // intended: state 1, G = mean(6, 8, 4) = 6, σ_G = sample σ of (6, 8, 4) = 2
        assert_eq!(a.target, 1);
        assert!((a.globminavg - 6.0).abs() < 1e-15, "{}", a.globminavg);
        assert!((a.globminfluc - 2.0).abs() < 1e-15, "{}", a.globminfluc);
        // a sampled state with a mean ≤ 0 is the target under both rules
        let mut ours = search(2, true);
        let mut gxx = search(2, true);
        for vi in [[9.0, -1.0], [9.0, -3.0]] {
            track_state(&mut ours, &vi, AedsTarget::Intended);
            track_state(&mut gxx, &vi, AedsTarget::GromosxxState1Default);
        }
        let a = track_state(&mut ours, &[9.0, -2.0], AedsTarget::Intended);
        let b = track_state(&mut gxx, &[9.0, -2.0], AedsTarget::GromosxxState1Default);
        assert_eq!((a.target, b.target), (1, 1));
        assert_eq!(a.globminavg, b.globminavg);
    }

    /// The four E_min branches (`eds.cc:296-325`) and the D15 refusal.
    #[test]
    fn emin_branches_and_the_zero_variance_refusal() {
        let p = AedsParams {
            form: AedsForm::SearchEmaxEmin,
            bmax_type: 1,
            bmax: 1.0,
            asteps: 10,
            bsteps: 10,
            columns: AedsColumns::Intended,
            prevalence: AedsPrevalence::Intended,
            target: AedsTarget::Intended,
        };
        let t = |g: f64| StateTrack {
            state: 0,
            round_trip: false,
            target: 0,
            globminavg: g,
            globminfluc: 0.0,
        };
        let run = |emax: f64, g: f64, fullemin: bool, p: &AedsParams| {
            let mut s = search(2, fullemin);
            s.emaxcounts = 1;
            s.emax = emax;
            emax_emin_update(p, &mut s, -100.0, &t(g)).map(|_| s.emin)
        };
        // emax − G ≤ b: no acceleration
        assert_eq!(run(0.5, 0.0, true, &p).unwrap(), 0.5);
        // 2(G + b) − emax ≥ G
        assert_eq!(run(1.5, 0.0, true, &p).unwrap(), 0.5);
        // below G: the quadratic emax − (emax − G)²/(2b)
        assert_eq!(run(4.0, 0.0, true, &p).unwrap(), -4.0);
        // RESTREMIN 1 before the first round trip: held at G
        assert_eq!(run(4.0, 0.0, false, &p).unwrap(), 0.0);
        // BMAXTYPE 2 with σ = 0: b = 0, gromosXX's −∞ is refused …
        let p2 = AedsParams { bmax_type: 2, ..p };
        let e = run(4.0, 0.0, true, &p2).unwrap_err();
        assert!(e.contains("BMAXTYPE 1 or RESTREMIN 1"), "{e}");
        // … unless the clamp catches it, as in gromosXX
        assert_eq!(run(4.0, 0.0, false, &p2).unwrap(), 0.0);
    }

    /// FORMS 5/6: the dominant state's offset goes down by ASTEPS·Δt·p/(round trips + 1), then
    /// the offsets are centred on zero.
    #[test]
    fn advanced_offsets_lower_the_dominant_state_and_centre() {
        let beta = 1.0 / (0.00831441 * 300.0);
        let mut s = search(3, true);
        s.eir = vec![0.0, 1.0, -1.0];
        let mut scratch = [0.0; 3];
        for rule in [AedsPrevalence::Intended, AedsPrevalence::GromosxxPlainExp] {
            let mut s = s.clone();
            advanced_offset_update(
                &mut s,
                &[0.0, 50.0, 50.0],
                0.0,
                beta,
                0.04,
                rule,
                &mut scratch,
            );
            assert!(s.eir.iter().sum::<f64>().abs() < 1e-15);
            let shift = 0.04 * 1.0; // prevalence of state 0 ≈ 1
            assert!(
                (s.eir[0] - (-shift - (-shift) / 3.0)).abs() < 1e-9,
                "{:?}",
                s.eir
            );
            assert!(s.eir[1] - s.eir[2] - 2.0 < 1e-12);
        }
    }

    /// D25: the log-sum-exp prevalences equal gromosXX's plain ones (to rounding) wherever those
    /// are finite, and stay finite where e^{−β(V_i − E_i^R)} overflows (or all underflow): the
    /// dominant state's offset then moves by the hand-computed step·(1 − 1/N), the others by
    /// step/N, where gromosXX's NaN prevalences move nothing.
    #[test]
    fn log_sum_exp_prevalences_match_where_finite_and_survive_overflow() {
        let beta = 1.0 / (0.00831441 * 300.0);
        let step = 0.04;
        let run = |vi: &[f64], eir: &[f64], rule| {
            let mut s = search(vi.len(), true);
            s.eir = eir.to_vec();
            let mut scratch = vec![0.0; vi.len()];
            advanced_offset_update(&mut s, vi, 0.0, beta, step, rule, &mut scratch);
            s
        };
        // moderate energies: both finite, equal to rounding
        for (vi, eir) in [
            (vec![3.0, -2.0, 1.5], vec![0.0, 0.5, -1.0]),
            (vec![-700.0, -690.0, -705.0], vec![0.0, 0.0, 0.0]),
            (vec![15.39, -16.08, 3.5, 0.2], vec![0.0, 10.0, -5.0, 2.0]),
        ] {
            let a = run(&vi, &eir, AedsPrevalence::Intended);
            let b = run(&vi, &eir, AedsPrevalence::GromosxxPlainExp);
            for (x, y) in a.eir.iter().zip(&b.eir) {
                assert!((x - y).abs() <= 1e-15 * y.abs().max(1.0), "{x} vs {y}");
            }
            assert_eq!(a.framecounts, b.framecounts);
        }
        // V_1 − E_1^R = −2000 kJ/mol: −β·(−2000) = 801.8 > 709.78, e^x = ∞
        let vi = [-2000.0, 0.0, 0.0];
        assert_eq!((-beta * vi[0]).exp(), f64::INFINITY);
        let plain = run(&vi, &[0.0; 3], AedsPrevalence::GromosxxPlainExp);
        assert_eq!(plain.eir, vec![0.0; 3], "gromosXX: frozen");
        let ours = run(&vi, &[0.0; 3], AedsPrevalence::Intended);
        let (lead, rest) = (-step * (1.0 - 1.0 / 3.0), step / 3.0);
        assert!((ours.eir[0] - lead).abs() < 1e-15, "{:?}", ours.eir);
        assert!((ours.eir[1] - rest).abs() < 1e-15 && (ours.eir[2] - rest).abs() < 1e-15);
        // all above +745/β: every e^x underflows to 0, 0/0 = NaN in gromosXX
        let vi = [2000.0, 2010.0, 2030.0];
        let plain = run(&vi, &[0.0; 3], AedsPrevalence::GromosxxPlainExp);
        assert_eq!(plain.eir, vec![0.0; 3], "gromosXX: frozen");
        let ours = run(&vi, &[0.0; 3], AedsPrevalence::Intended);
        // p_1 = 1/(1 + e^{−10β} + e^{−30β}) = 1 − 1.8e-2
        let p1 = 1.0 / (1.0 + (-10.0 * beta).exp() + (-30.0 * beta).exp());
        assert!(
            (ours.eir[0] - -step * p1 * (2.0 / 3.0)).abs() < 1e-15,
            "{:?}",
            ours.eir
        );
        assert!(ours.eir.iter().all(|e| e.is_finite()));
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
