//! Perturbed (free-energy) LJ+CRF interactions with soft-core potentials

use super::innerloops::with_groups;
use super::lj_crf_interaction;
use super::{
    CRFParameters, EnergyGroups, GroupEnergies, LJParamMatrix, LJParameters, PrecalcNonbonded,
};
use newen_core::lambdas::{IndividualLambdas, InteractionLambda};
use newen_core::math::{BoundaryCondition, Vec3};
use newen_core::topology::{PairInteraction, PerturbedAtomPair};

/// Lambda-dependent parameters for perturbed nonbonded interactions
///
/// Direct translation from gromosXX perturbed_nonbonded_term.h/.cc
#[derive(Debug, Clone)]
pub struct PerturbedLambdaParams {
    /// Lambda for LJ interaction, state A: (1-λ)^n
    pub a_lj_lambda_n: f64,
    /// Lambda for LJ interaction, state B: λ^n
    pub b_lj_lambda_n: f64,
    /// Lambda for CRF interaction, state A: (1-λ)^n
    pub a_crf_lambda_n: f64,
    /// Lambda for CRF interaction, state B: λ^n
    pub b_crf_lambda_n: f64,

    /// Lambda derivative for LJ, state A: (1-λ)^(n-1) · dΛ/dλ (the individual-λ slope)
    pub a_lj_lambda_n_1: f64,
    /// Lambda derivative for LJ, state B: λ^(n-1) · dΛ/dλ
    pub b_lj_lambda_n_1: f64,
    /// Lambda derivative for CRF, state A: (1-λ)^(n-1) · dΛ/dλ
    pub a_crf_lambda_n_1: f64,
    /// Lambda derivative for CRF, state B: λ^(n-1) · dΛ/dλ
    pub b_crf_lambda_n_1: f64,

    /// (1-λ_ljs) · dΛ_ljs/dλ: the slope of state B's LJ softness
    pub a_ljs_lambda: f64,
    /// λ_ljs · dΛ_ljs/dλ: the slope of state A's LJ softness
    pub b_ljs_lambda: f64,
    /// (1-λ_ljs)²: the LJ softness of state B
    pub a_ljs_lambda2: f64,
    /// λ_ljs²: the LJ softness of state A
    pub b_ljs_lambda2: f64,

    /// (1-λ_crfs) · dΛ_crfs/dλ: the slope of state B's CRF softness
    pub a_crfs_lambda: f64,
    /// λ_crfs · dΛ_crfs/dλ: the slope of state A's CRF softness
    pub b_crfs_lambda: f64,
    /// (1-λ_crfs)²: the CRF softness of state B
    pub a_crfs_lambda2: f64,
    /// λ_crfs²: the CRF softness of state A
    pub b_crfs_lambda2: f64,

    /// Lambda exponent (n)
    pub lambda_exp: i32,
}

impl PerturbedLambdaParams {
    /// gromosXX `set_lambda` (`perturbed_nonbonded_term.cc:322-347`), whose `pow(x, n)` is the
    /// double-precision `pow`: `powf`, not `powi` (they differ in the last bits for n ≥ 2).
    ///
    /// spec: nonbonded.perturbed_lambda_params
    pub fn from_lambda(
        lj_lambda: f64,
        ljs_lambda: f64,
        crf_lambda: f64,
        crfs_lambda: f64,
        lj_lambda_derivative: f64,
        ljs_lambda_derivative: f64,
        crf_lambda_derivative: f64,
        crfs_lambda_derivative: f64,
        n: i32,
    ) -> Self {
        Self {
            // State A: (1-λ)^n
            a_lj_lambda_n: (1.0 - lj_lambda).powf(n as f64),
            a_crf_lambda_n: (1.0 - crf_lambda).powf(n as f64),

            // State B: λ^n
            b_lj_lambda_n: lj_lambda.powf(n as f64),
            b_crf_lambda_n: crf_lambda.powf(n as f64),

            // Derivatives: (1-λ)^(n-1) · dΛ/dλ
            a_lj_lambda_n_1: (1.0 - lj_lambda).powf((n - 1) as f64) * lj_lambda_derivative,
            a_crf_lambda_n_1: (1.0 - crf_lambda).powf((n - 1) as f64) * crf_lambda_derivative,

            // Derivatives: λ^(n-1) · dΛ/dλ
            b_lj_lambda_n_1: lj_lambda.powf((n - 1) as f64) * lj_lambda_derivative,
            b_crf_lambda_n_1: crf_lambda.powf((n - 1) as f64) * crf_lambda_derivative,

            // Softness slopes: (1-λ_s) · dΛ_s/dλ
            a_ljs_lambda: (1.0 - ljs_lambda) * ljs_lambda_derivative,
            a_crfs_lambda: (1.0 - crfs_lambda) * crfs_lambda_derivative,

            // Softness slopes: λ_s · dΛ_s/dλ
            b_ljs_lambda: ljs_lambda * ljs_lambda_derivative,
            b_crfs_lambda: crfs_lambda * crfs_lambda_derivative,

            // Squared soft-core lambdas
            a_ljs_lambda2: (1.0 - ljs_lambda) * (1.0 - ljs_lambda),
            a_crfs_lambda2: (1.0 - crfs_lambda) * (1.0 - crfs_lambda),
            b_ljs_lambda2: ljs_lambda * ljs_lambda,
            b_crfs_lambda2: crfs_lambda * crfs_lambda,

            lambda_exp: n,
        }
    }

    /// One parameter set per energy-group pair, row-major `[g1][g2]` over the groups of
    /// `lambdas`, from the individual λ of the four nonbonded kinds and NLAM `n` — what gromosXX's
    /// innerloops pass to `set_lambda` for a pair in groups `(g1, g2)`
    /// (`perturbed_nonbonded_innerloop.cc:191-200`). `set_lambda` is a pure function of its
    /// inputs, so one evaluation per pair of groups gives the per-pair values exactly.
    pub fn fill_table(lambdas: &IndividualLambdas, n: i32, table: &mut Vec<Self>) {
        use InteractionLambda::{Crf, CrfSoftness, Lj, LjSoftness};
        let groups = lambdas.n_groups();
        table.clear();
        for g1 in 0..groups {
            for g2 in 0..groups {
                let (lj, dlj) = lambdas.get(Lj, g1, g2);
                let (ljs, dljs) = lambdas.get(LjSoftness, g1, g2);
                let (crf, dcrf) = lambdas.get(Crf, g1, g2);
                let (crfs, dcrfs) = lambdas.get(CrfSoftness, g1, g2);
                table.push(Self::from_lambda(
                    lj, ljs, crf, crfs, dlj, dljs, dcrf, dcrfs, n,
                ));
            }
        }
    }
}

/// The parameter set of the pair `(i, j)`: its energy-group cell when there are several groups
/// (the table is then `n_groups²` long), the only entry otherwise.
#[inline(always)]
fn pair_params<'a, const GROUPS: bool>(
    lps: &'a [PerturbedLambdaParams],
    egroups: &EnergyGroups<'_>,
    i: usize,
    j: usize,
) -> (usize, &'a PerturbedLambdaParams) {
    let k = if GROUPS { egroups.cell(i, j) } else { 0 };
    (k, &lps[k])
}

/// Perturbed LJ+CRF interaction with soft-core potentials
///
/// Direct translation from eds_pert_lj_crf_interaction() in perturbed_nonbonded_term.cc
///
/// # Arguments
/// * `r` - Distance vector
/// * `a_c6, a_c12` - State A LJ parameters
/// * `b_c6, b_c12` - State B LJ parameters
/// * `a_q, b_q` - State A and B charge products
/// * `alpha_lj, alpha_crf` - Soft-core parameters
/// * `lambda_params` - Lambda values and derivatives
/// * `crf` - CRF parameters
///
/// # Returns
/// * `force_magnitude` - Scalar force magnitude
/// * `e_lj` - Lennard-Jones energy
/// * `e_crf` - Coulomb reaction field energy
/// * `de_lj` - Lambda derivative of LJ energy (dH/dλ)
/// * `de_crf` - Lambda derivative of CRF energy (dH/dλ)
#[inline]
pub fn perturbed_lj_crf_interaction(
    r: Vec3,
    a_c6: f64,
    a_c12: f64,
    b_c6: f64,
    b_c12: f64,
    a_q: f64,
    b_q: f64,
    alpha_lj: f64,
    alpha_crf: f64,
    lambda_params: &PerturbedLambdaParams,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
) -> (f64, f64, f64, f64, f64) {
    perturbed_lj_crf_interaction_scaled(
        r,
        a_c6,
        a_c12,
        b_c6,
        b_c12,
        a_q,
        b_q,
        alpha_lj,
        alpha_crf,
        lambda_params,
        crf,
        four_pi_eps_i,
        1.0,
    )
}

/// [`perturbed_lj_crf_interaction`] with gromosXX's `coulomb_scaling` on the `1/r` Coulomb
/// part (`perturbed_nonbonded_term.cc:387-510`); a scaling of 1 gives the same bits.
#[inline]
pub fn perturbed_lj_crf_interaction_scaled(
    r: Vec3,
    a_c6: f64,
    a_c12: f64,
    b_c6: f64,
    b_c12: f64,
    a_q: f64,
    b_q: f64,
    alpha_lj: f64,
    alpha_crf: f64,
    lambda_params: &PerturbedLambdaParams,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
    coulomb_scaling: f64,
) -> (f64, f64, f64, f64, f64) {
    perturbed_lj_crf_pair(
        r.length_squared(),
        a_c6,
        a_c12,
        b_c6,
        b_c12,
        a_q,
        b_q,
        alpha_lj,
        alpha_crf,
        lambda_params,
        crf,
        four_pi_eps_i,
        coulomb_scaling,
    )
}

/// One perturbed LJ + CRF pair with soft-core at squared distance `r2`: the kernel of
/// [`perturbed_lj_crf_interaction_scaled`]. Returns `(force, e_lj, e_crf, de_lj, de_crf)`, the
/// force a scalar along the pair vector.
///
/// spec: nonbonded.perturbed_lj_crf_pair
#[inline]
pub fn perturbed_lj_crf_pair(
    r2: f64,
    a_c6: f64,
    a_c12: f64,
    b_c6: f64,
    b_c12: f64,
    a_q: f64,
    b_q: f64,
    alpha_lj: f64,
    alpha_crf: f64,
    lambda_params: &PerturbedLambdaParams,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
    coulomb_scaling: f64,
) -> (f64, f64, f64, f64, f64) {
    let r6 = r2 * r2 * r2;
    let r4 = r2 * r2;

    // C12/C6 of the soft-core length; 0 for C6 = 0, as gromosXX
    let a_c126 = if a_c6 != 0.0 { a_c12 / a_c6 } else { 0.0 };
    let b_c126 = if b_c6 != 0.0 { b_c12 / b_c6 } else { 0.0 };

    // ===== STATE A SOFT-CORE DISTANCES =====
    // State A is softened by λ_soft² (it vanishes as λ → 1)
    // Electrostatics: r_soft² = r² + α_crf · λ_crfs²
    let a_dist2soft = r2 + alpha_crf * lambda_params.b_crfs_lambda2;
    let a_distisoft = 1.0 / a_dist2soft.sqrt();
    let a_dist3isoft = a_distisoft / a_dist2soft;

    // LJ: r⁶_soft = r⁶ + α_lj · λ_ljs² · C12/C6
    let a_dist6soft = r6 + alpha_lj * lambda_params.b_ljs_lambda2 * a_c126;
    let a_dist6isoft = 1.0 / a_dist6soft;

    // ===== STATE B SOFT-CORE DISTANCES =====
    // State B is softened by (1 − λ_soft)² (it vanishes as λ → 0)
    let b_dist2soft = r2 + alpha_crf * lambda_params.a_crfs_lambda2;
    let b_distisoft = 1.0 / b_dist2soft.sqrt();
    let b_dist3isoft = b_distisoft / b_dist2soft;

    let b_dist6soft = r6 + alpha_lj * lambda_params.a_ljs_lambda2 * b_c126;
    let b_dist6isoft = 1.0 / b_dist6soft;

    // ===== SOFT-CORE CUTOFF MODIFICATIONS =====
    let cut2 = crf.cutoff_sq; // RCRF², the soft-core reaction-field radius
    let a_cut2soft = cut2 + alpha_crf * lambda_params.b_crfs_lambda2;
    let b_cut2soft = cut2 + alpha_crf * lambda_params.a_crfs_lambda2;

    let a_cut2soft3 = a_cut2soft * a_cut2soft * a_cut2soft;
    let b_cut2soft3 = b_cut2soft * b_cut2soft * b_cut2soft;

    // CRF constants with soft-core cutoff
    // GROMOS: A_crf_2cut3i = (crf/2) / sqrt(A_cut2soft³). `crf.crf_2cut3i` is crf/(2·cutoff³),
    // so the softened form is crf_2cut3i · cutoff³ / sqrt(cut2soft³) — until 0.0.34 this divided
    // crf_2cut3i by sqrt(cut2soft³) once more (an extra 1/cutoff³ on the r² reaction-field
    // term of every perturbed pair; invisible with zero charges, 0.16 kJ/mol on methanol).
    let cut3 = cut2 * cut2.sqrt();
    let a_crf_2cut3i = crf.crf_2cut3i * cut3 / a_cut2soft3.sqrt();
    let b_crf_2cut3i = crf.crf_2cut3i * cut3 / b_cut2soft3.sqrt();

    let a_crf_cut3i = 2.0 * a_crf_2cut3i;
    let b_crf_cut3i = 2.0 * b_crf_2cut3i;

    // Derivative terms for soft-core
    let a_crf_pert = 3.0 * a_crf_2cut3i / a_cut2soft;
    let b_crf_pert = 3.0 * b_crf_2cut3i / b_cut2soft;

    // ===== FORCE CALCULATION =====
    // CRF force
    let mut force =
        (lambda_params.a_crf_lambda_n * a_q * (coulomb_scaling * a_dist3isoft + a_crf_cut3i)
            + lambda_params.b_crf_lambda_n * b_q * (coulomb_scaling * b_dist3isoft + b_crf_cut3i))
            * four_pi_eps_i;

    // LJ attractive force: -6 * C6 / r⁸
    force += -6.0
        * (lambda_params.a_lj_lambda_n * a_c6 * a_dist6isoft * a_dist6isoft
            + lambda_params.b_lj_lambda_n * b_c6 * b_dist6isoft * b_dist6isoft)
        * r4;

    // LJ repulsive force: 12 * C12 / r¹⁴
    force += 12.0
        * (lambda_params.a_lj_lambda_n * a_c12 * a_dist6isoft * a_dist6isoft * a_dist6isoft
            + lambda_params.b_lj_lambda_n * b_c12 * b_dist6isoft * b_dist6isoft * b_dist6isoft)
        * r4;

    // ===== ENERGY CALCULATION =====
    let a_e_lj = (a_c12 * a_dist6isoft - a_c6) * a_dist6isoft;
    let b_e_lj = (b_c12 * b_dist6isoft - b_c6) * b_dist6isoft;

    let crf_cut = crf.crf_cut; // (1 - crf/2) / cutoff  [was wrongly crf_cut3i before]
    let a_e_crf = a_q * (coulomb_scaling * a_distisoft - a_crf_2cut3i * r2 - crf_cut);
    let b_e_crf = b_q * (coulomb_scaling * b_distisoft - b_crf_2cut3i * r2 - crf_cut);

    // Total energy: E = (1-λ)^n * E_A + λ^n * E_B
    let e_lj = lambda_params.a_lj_lambda_n * a_e_lj + lambda_params.b_lj_lambda_n * b_e_lj;
    let e_crf = (lambda_params.a_crf_lambda_n * a_e_crf + lambda_params.b_crf_lambda_n * b_e_crf)
        * four_pi_eps_i;

    // ===== LAMBDA DERIVATIVES (dH/dλ) =====
    // LJ derivative: soft-core contribution + direct lambda contribution
    let de_lj = -2.0
        * alpha_lj
        * (lambda_params.a_lj_lambda_n
            * lambda_params.b_ljs_lambda
            * a_c126
            * a_dist6isoft
            * a_dist6isoft
            * (2.0 * a_c12 * a_dist6isoft - a_c6)
            - lambda_params.b_lj_lambda_n
                * lambda_params.a_ljs_lambda
                * b_c126
                * b_dist6isoft
                * b_dist6isoft
                * (2.0 * b_c12 * b_dist6isoft - b_c6))
        + (lambda_params.lambda_exp as f64)
            * (lambda_params.b_lj_lambda_n_1 * b_e_lj - lambda_params.a_lj_lambda_n_1 * a_e_lj);

    // CRF derivative: soft-core contribution + direct lambda contribution
    let de_crf = -(lambda_params.a_crf_lambda_n
        * a_q
        * lambda_params.b_crfs_lambda
        * (coulomb_scaling * a_dist3isoft - a_crf_pert * r2)
        - lambda_params.b_crf_lambda_n
            * b_q
            * lambda_params.a_crfs_lambda
            * (coulomb_scaling * b_dist3isoft - b_crf_pert * r2))
        * four_pi_eps_i
        * alpha_crf
        + (lambda_params.lambda_exp as f64)
            * (lambda_params.b_crf_lambda_n_1 * b_e_crf - lambda_params.a_crf_lambda_n_1 * a_e_crf)
            * four_pi_eps_i;

    (force, e_lj, e_crf, de_lj, de_crf)
}

// ─── Extended TI (PRECALCLAM) ──────────────────────────────────────────────────

/// gromosXX `lj_crf_soft_interaction_ext` (`perturbed_nonbonded_term.cc:614-716`): the A- and
/// B-state LJ and RF energies of a perturbed pair at λ = `lam`, unweighted, with softness
/// `λ_soft = lam` (no NLAM, no individual λ) and their soft-core derivatives with
/// `A_de = +∂A/∂λ`, `B_de = −∂B/∂λ` (extended TI, JCTC 2016, 12, 4476, Eqs. 9/12; item 1.8).
/// `r2` is the squared distance; the Coulomb scaling is 1, as at every gromosXX call site (the
/// 1-4 one included).
///
/// Returns `([A_e_lj, B_e_lj, A_e_crf, B_e_crf], [A_de_lj, B_de_lj, A_de_crf, B_de_crf])`.
///
/// spec: nonbonded.lj_crf_soft_interaction_ext
#[inline]
pub fn lj_crf_soft_interaction_ext(
    r2: f64,
    a_c6: f64,
    a_c12: f64,
    b_c6: f64,
    b_c12: f64,
    a_q: f64,
    b_q: f64,
    alpha_lj: f64,
    alpha_crf: f64,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
    lam: f64,
) -> ([f64; 4], [f64; 4]) {
    let a_ljs_lambda2 = (1.0 - lam) * (1.0 - lam);
    let b_ljs_lambda2 = lam * lam;
    let (a_crfs_lambda2, b_crfs_lambda2) = (a_ljs_lambda2, b_ljs_lambda2);
    let a_c126 = if a_c6 != 0.0 { a_c12 / a_c6 } else { 0.0 };
    let b_c126 = if b_c6 != 0.0 { b_c12 / b_c6 } else { 0.0 };

    let a_dist2soft = r2 + alpha_crf * b_crfs_lambda2;
    let b_dist2soft = r2 + alpha_crf * a_crfs_lambda2;
    let a_distisoft = 1.0 / a_dist2soft.sqrt();
    let b_distisoft = 1.0 / b_dist2soft.sqrt();
    let a_dist3isoft = a_distisoft / a_dist2soft;
    let b_dist3isoft = b_distisoft / b_dist2soft;

    let dist6 = r2 * r2 * r2;
    let a_dist6isoft = 1.0 / (dist6 + alpha_lj * b_ljs_lambda2 * a_c126);
    let b_dist6isoft = 1.0 / (dist6 + alpha_lj * a_ljs_lambda2 * b_c126);

    let cut2 = crf.cutoff_sq;
    let a_cut2soft = cut2 + alpha_crf * b_crfs_lambda2;
    let b_cut2soft = cut2 + alpha_crf * a_crfs_lambda2;
    // crf/2 / sqrt(cut2soft³), with crf/2 = crf_2cut3i·cut³ as in `perturbed_lj_crf_interaction`
    let cut3 = cut2 * cut2.sqrt();
    let a_crf_2cut3i = crf.crf_2cut3i * cut3 / (a_cut2soft * a_cut2soft * a_cut2soft).sqrt();
    let b_crf_2cut3i = crf.crf_2cut3i * cut3 / (b_cut2soft * b_cut2soft * b_cut2soft).sqrt();
    let a_crf_pert = 3.0 * a_crf_2cut3i / a_cut2soft;
    let b_crf_pert = 3.0 * b_crf_2cut3i / b_cut2soft;

    let e = [
        (a_c12 * a_dist6isoft - a_c6) * a_dist6isoft,
        (b_c12 * b_dist6isoft - b_c6) * b_dist6isoft,
        a_q * (a_distisoft - a_crf_2cut3i * r2 - crf.crf_cut) * four_pi_eps_i,
        b_q * (b_distisoft - b_crf_2cut3i * r2 - crf.crf_cut) * four_pi_eps_i,
    ];
    let de = [
        -2.0 * alpha_lj
            * lam
            * a_c126
            * a_dist6isoft
            * a_dist6isoft
            * (2.0 * a_c12 * a_dist6isoft - a_c6),
        -2.0 * alpha_lj
            * (1.0 - lam)
            * b_c126
            * b_dist6isoft
            * b_dist6isoft
            * (2.0 * b_c12 * b_dist6isoft - b_c6),
        -a_q * lam * (a_dist3isoft - a_crf_pert * r2) * four_pi_eps_i * alpha_crf,
        -b_q * (1.0 - lam) * (b_dist3isoft - b_crf_pert * r2) * four_pi_eps_i * alpha_crf,
    ];
    (e, de)
}

/// gromosXX `rf_soft_interaction_ext` (`perturbed_nonbonded_term.cc:853-881`): the reaction-field
/// term of an excluded pair (or, at `r2 = 0`, the self term) in states A and B at λ = `lam`.
/// Unlike `rf_soft_interaction` there is no self-term correction.
///
/// Returns `([A_e_rf, B_e_rf], [A_de_rf, B_de_rf])`, the derivatives in the sign convention of
/// [`lj_crf_soft_interaction_ext`].
///
/// spec: nonbonded.rf_soft_interaction_ext
#[inline]
pub fn rf_soft_interaction_ext(
    r2: f64,
    a_q: f64,
    b_q: f64,
    alpha_crf: f64,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
    lam: f64,
) -> ([f64; 2], [f64; 2]) {
    let cut2 = crf.cutoff_sq;
    let a_cut2soft = cut2 + alpha_crf * lam * lam;
    let b_cut2soft = cut2 + alpha_crf * (1.0 - lam) * (1.0 - lam);
    let cut3 = cut2 * cut2.sqrt();
    let a_crf_2cut3i = crf.crf_2cut3i * cut3 / (a_cut2soft * a_cut2soft * a_cut2soft).sqrt();
    let b_crf_2cut3i = crf.crf_2cut3i * cut3 / (b_cut2soft * b_cut2soft * b_cut2soft).sqrt();
    let a_crf_pert = 3.0 * a_crf_2cut3i / a_cut2soft;
    let b_crf_pert = 3.0 * b_crf_2cut3i / b_cut2soft;
    (
        [
            a_q * (-a_crf_2cut3i * r2 - crf.crf_cut) * four_pi_eps_i,
            b_q * (-b_crf_2cut3i * r2 - crf.crf_cut) * four_pi_eps_i,
        ],
        [
            a_q * lam * a_crf_pert * r2 * four_pi_eps_i * alpha_crf,
            b_q * (1.0 - lam) * b_crf_pert * r2 * four_pi_eps_i * alpha_crf,
        ],
    )
}

impl PrecalcNonbonded {
    /// Book [`lj_crf_soft_interaction_ext`] of one pair at every λ_k into `cell`.
    #[inline]
    pub fn add_lj_crf_ext(
        &mut self,
        cell: usize,
        r2: f64,
        lj: (f64, f64, f64, f64),
        a_q: f64,
        b_q: f64,
        alpha_lj: f64,
        alpha_crf: f64,
        crf: &CRFParameters,
        four_pi_eps_i: f64,
    ) {
        let (a_c6, a_c12, b_c6, b_c12) = lj;
        for k in 0..self.lambdas.len() {
            let (e, de) = lj_crf_soft_interaction_ext(
                r2,
                a_c6,
                a_c12,
                b_c6,
                b_c12,
                a_q,
                b_q,
                alpha_lj,
                alpha_crf,
                crf,
                four_pi_eps_i,
                self.lambdas[k],
            );
            self.book(k, cell, e, de);
        }
    }

    /// Book `factor` × [`rf_soft_interaction_ext`] of one pair (or self term) at every λ_k into
    /// the CRF columns of `cell`.
    #[inline]
    pub fn add_rf_ext(
        &mut self,
        cell: usize,
        r2: f64,
        a_q: f64,
        b_q: f64,
        alpha_crf: f64,
        factor: f64,
        crf: &CRFParameters,
        four_pi_eps_i: f64,
    ) {
        for k in 0..self.lambdas.len() {
            let ([a_e, b_e], [a_de, b_de]) = rf_soft_interaction_ext(
                r2,
                a_q,
                b_q,
                alpha_crf,
                crf,
                four_pi_eps_i,
                self.lambdas[k],
            );
            self.book(
                k,
                cell,
                [0.0, 0.0, factor * a_e, factor * b_e],
                [0.0, 0.0, factor * a_de, factor * b_de],
            );
        }
    }

    /// Book λ-independent energies (an unperturbed pair of a `PERTATOMPAIR`) into every λ_k,
    /// without a derivative (gromosXX `perturbed_nonbonded_pair.cc:410-420, 600-612`).
    #[inline]
    pub fn add_constant(&mut self, cell: usize, e: [f64; 4]) {
        for k in 0..self.lambdas.len() {
            self.book(k, cell, e, [0.0; 4]);
        }
    }
}

// ─── Step 3: perturbed nonbonded correction functions ─────────────────────────

/// Per-atom perturbation data for fast lookup during nonbonded force calculation.
#[derive(Debug, Clone)]
pub struct PertAtomInfo {
    pub a_iac: usize,
    pub b_iac: usize,
    pub a_charge: f64,
    pub b_charge: f64,
    pub alpha_lj: f64,
    pub alpha_crf: f64,
}

/// Delta result applied on top of the state-A result already in `nonbonded_storage` — or, from
/// [`perturbed_atom_pair_interactions`], the pairs' whole contribution.
pub struct PertNBCorrection {
    pub forces: Vec<Vec3>,
    pub delta_e_lj: f64,
    pub delta_e_crf: f64,
    /// dH/dλ of the corrected pairs (LJ + CRF).
    pub dhdl: f64,
    /// The LJ part of `dhdl` — kept apart so the `.trg` file carries the same per-term
    /// derivatives gromosXX writes (FREEENERDERIVS03 totals: LJ, CRF).
    pub dhdl_lj: f64,
    /// The CRF part of `dhdl`.
    pub dhdl_crf: f64,
    /// Virial `r ⊗ f` of the pairs at the current λ (atom pairs, 1-4, excluded pairs); from
    /// [`perturbed_pairlist_correction`] the virial of the force change `f(λ) − f_A`.
    pub virial: [[f64; 3]; 3],
    /// Virial `r ⊗ f_A` of the state-A forces the regular 1-4 and excluded-pair loops already
    /// put into the nonbonded storage for these pairs (1-4 and excluded corrections only).
    pub virial_a: [[f64; 3]; 3],
}

/// `m[a][b] += r(a)·f(b)`.
#[inline]
fn add_outer(m: &mut [[f64; 3]; 3], r: Vec3, f: Vec3) {
    let (rv, fv) = ([r.x, r.y, r.z], [f.x, f.y, f.z]);
    for (x, row) in m.iter_mut().enumerate() {
        for (y, v) in row.iter_mut().enumerate() {
            *v += rv[x] * fv[y];
        }
    }
}

impl PertNBCorrection {
    pub fn new(n_atoms: usize) -> Self {
        Self {
            forces: vec![Vec3::ZERO; n_atoms],
            delta_e_lj: 0.0,
            delta_e_crf: 0.0,
            dhdl: 0.0,
            dhdl_lj: 0.0,
            dhdl_crf: 0.0,
            virial: [[0.0; 3]; 3],
            virial_a: [[0.0; 3]; 3],
        }
    }
}

/// Correct pairlist energies/forces for pairs involving ≥1 perturbed atom.
///
/// The regular innerloop ran with state-A charges/IAC.  For each perturbed pair
/// this function subtracts the state-A contribution and adds the full dual-topology
/// perturbed result, returning the net delta.
///
/// Alpha convention (faithful to GROMOS):
///   single perturbed atom  →  that atom's α
///   both perturbed         →  average of the two atoms' α values
///
/// With more than one energy group the energy change and the dH/dλ of each pair also go to
/// cell `[eg(i)][eg(j)]` of `sink` (gromosXX `perturbed_nonbonded_innerloop.cc:391-403`).
/// With `precalc` (PRECALCLAM on and this step one it is evaluated on) each pair's A/B energies
/// at every λ_k go there too (`:315-370`, `:547-600`).
pub fn perturbed_pairlist_correction<BC: BoundaryCondition>(
    pos: &[Vec3],
    pairlist: &[(u32, u32)],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    iac: &[u32],
    lj: &LJParamMatrix,
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    four_pi_eps_i: f64,
    egroups: EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: Option<&mut PrecalcNonbonded>,
) -> PertNBCorrection {
    match precalc {
        Some(pc) => with_groups!(egroups, G => pairlist_correction_impl::<BC, G, true>(
            pos,
            pairlist,
            pert,
            charges,
            iac,
            lj,
            crf,
            lps,
            bc,
            four_pi_eps_i,
            &egroups,
            sink,
            pc,
        )),
        None => with_groups!(egroups, G => pairlist_correction_impl::<BC, G, false>(
            pos,
            pairlist,
            pert,
            charges,
            iac,
            lj,
            crf,
            lps,
            bc,
            four_pi_eps_i,
            &egroups,
            sink,
            &mut PrecalcNonbonded::default(),
        )),
    }
}

fn pairlist_correction_impl<BC: BoundaryCondition, const GROUPS: bool, const PRECALC: bool>(
    pos: &[Vec3],
    pairlist: &[(u32, u32)],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    iac: &[u32],
    lj: &LJParamMatrix,
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    four_pi_eps_i: f64,
    egroups: &EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: &mut PrecalcNonbonded,
) -> PertNBCorrection {
    let n = pos.len();
    let mut out = PertNBCorrection::new(n);

    for &(ii, jj) in pairlist {
        let i = ii as usize;
        let j = jj as usize;
        // gromosXX's insert_pair routes a pair to the perturbed pairlist when *either* atom is
        // perturbed (`if is_perturbed(i) … else if is_perturbed(j) → (j, i)`). Until 0.0.34 this
        // skipped pairs whose lower-index atom was unperturbed, so a perturbed atom that is not
        // at the start of the topology lost the soft-core treatment of most of its pairs
        // (invisible for CH4 — atom 1 — and methanol — atoms 1–3).
        if pert[i].is_none() && pert[j].is_none() {
            continue;
        }
        let pi = pert[i].as_ref();
        let pj = pert[j].as_ref();
        let (k, lp) = pair_params::<GROUPS>(lps, egroups, i, j);

        let a_iac_i = pi.map_or(iac[i] as usize, |p| p.a_iac);
        let a_iac_j = pj.map_or(iac[j] as usize, |p| p.a_iac);
        let b_iac_i = pi.map_or(iac[i] as usize, |p| p.b_iac);
        let b_iac_j = pj.map_or(iac[j] as usize, |p| p.b_iac);
        let a_q_i = pi.map_or(charges[i], |p| p.a_charge);
        let a_q_j = pj.map_or(charges[j], |p| p.a_charge);
        let b_q_i = pi.map_or(charges[i], |p| p.b_charge);
        let b_q_j = pj.map_or(charges[j], |p| p.b_charge);

        let lj_a = lj.get(a_iac_i, a_iac_j);
        let lj_b = lj.get(b_iac_i, b_iac_j);
        let a_q = a_q_i * a_q_j;
        let b_q = b_q_i * b_q_j;

        let (alpha_lj, alpha_crf) = match (pi, pj) {
            (Some(pi), Some(pj)) => (
                (pi.alpha_lj + pj.alpha_lj) * 0.5,
                (pi.alpha_crf + pj.alpha_crf) * 0.5,
            ),
            (Some(p), None) | (None, Some(p)) => (p.alpha_lj, p.alpha_crf),
            _ => unreachable!(),
        };

        let r = bc.nearest_image(pos[i], pos[j]);
        // lj_crf_interaction expects q_prod already scaled by four_pi_eps_i
        let (f_a, e_lj_a, e_crf_a) =
            lj_crf_interaction(r, lj_a.c6, lj_a.c12, a_q * four_pi_eps_i, crf);
        let (f_p, e_lj_p, e_crf_p, de_lj, de_crf) = perturbed_lj_crf_interaction(
            r,
            lj_a.c6,
            lj_a.c12,
            lj_b.c6,
            lj_b.c12,
            a_q,
            b_q,
            alpha_lj,
            alpha_crf,
            lp,
            crf,
            four_pi_eps_i,
        );

        let df = r * (f_p - f_a);
        out.forces[i] += df;
        out.forces[j] -= df;
        add_outer(&mut out.virial, r, df);
        out.delta_e_lj += e_lj_p - e_lj_a;
        out.delta_e_crf += e_crf_p - e_crf_a;
        out.dhdl += de_lj + de_crf;
        out.dhdl_lj += de_lj;
        out.dhdl_crf += de_crf;
        if GROUPS {
            sink.book(k, e_lj_p - e_lj_a, e_crf_p - e_crf_a, de_lj, de_crf);
        }
        if PRECALC {
            precalc.add_lj_crf_ext(
                k,
                r.length_squared(),
                (lj_a.c6, lj_a.c12, lj_b.c6, lj_b.c12),
                a_q,
                b_q,
                alpha_lj,
                alpha_crf,
                crf,
                four_pi_eps_i,
            );
        }
    }
    out
}

/// Correct RF self-energy for perturbed atoms.
///
/// The regular `rf_excluded_interactions` computed `-0.5*q_A²*FPEPSI*crf_cut`.
/// GROMOS adds `0.5*e_rf` from `rf_soft_interaction(r=0, selfterm_correction=true)`.
/// The net extra contribution is:
///   ΔE_self = 0.5*[(1−(1-λ)^n)*q_A² − λ^n*q_B²]*FPEPSI*crf_cut
///   Δ(dH/dλ)  = 0.5*n*[(1-λ)^(n-1)*q_A² − λ^(n-1)*q_B²]*FPEPSI*crf_cut
///
/// With more than one energy group each atom's terms also go to cell `[eg(i)][eg(i)]` of `sink`
/// (gromosXX `perturbed_nonbonded_innerloop.cc:3452-3453`); with `precalc`, half its
/// `rf_soft_interaction_ext` at r = 0 at every λ_k (`:3390-3418`).
pub fn perturbed_self_energy_correction(
    pert: &[Option<PertAtomInfo>],
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    four_pi_eps_i: f64,
    egroups: EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: Option<&mut PrecalcNonbonded>,
) -> (f64, f64) {
    let mut none = PrecalcNonbonded::default();
    let pc = precalc.unwrap_or(&mut none);
    with_groups!(egroups, G => self_energy_correction_impl::<G>(
        pert,
        crf,
        lps,
        four_pi_eps_i,
        &egroups,
        sink,
        pc,
    ))
}

fn self_energy_correction_impl<const GROUPS: bool>(
    pert: &[Option<PertAtomInfo>],
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    four_pi_eps_i: f64,
    egroups: &EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: &mut PrecalcNonbonded,
) -> (f64, f64) {
    let mut delta_e = 0.0;
    let mut dhdl = 0.0;

    for (i, pi) in pert
        .iter()
        .enumerate()
        .filter_map(|(i, x)| x.as_ref().map(|p| (i, p)))
    {
        // the self term takes the individual λ of the atom's own group (innerloop.cc:3367-3377)
        let (k, lp) = pair_params::<GROUPS>(lps, egroups, i, i);
        let qa2 = pi.a_charge * pi.a_charge;
        let qb2 = pi.b_charge * pi.b_charge;

        let (e_rf, de_rf) = perturbed_rf_self_term(qa2, qb2, lp, crf, four_pi_eps_i);
        delta_e += 0.5 * e_rf;
        dhdl += 0.5 * de_rf;
        if GROUPS {
            sink.book(k, 0.0, 0.5 * e_rf, 0.0, 0.5 * de_rf);
        }
        if precalc.is_active() {
            precalc.add_rf_ext(k, 0.0, qa2, qb2, pi.alpha_crf, 0.5, crf, four_pi_eps_i);
        }
    }
    (delta_e, dhdl)
}

/// The reaction-field self term of a perturbed atom with squared charges `qa2`, `qb2`, as the
/// change from the state-A self term the unperturbed loop booked: gromosXX `rf_soft_interaction`
/// at r = 0 with `selfterm_correction` (`perturbed_nonbonded_term.cc:718-765`), where the soft
/// terms vanish. Returns `(e_rf, de_rf)`; the caller books half of each.
///
/// spec: nonbonded.perturbed_rf_self_term
#[inline]
pub fn perturbed_rf_self_term(
    qa2: f64,
    qb2: f64,
    lp: &PerturbedLambdaParams,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
) -> (f64, f64) {
    let n = lp.lambda_exp as f64;
    // e_rf = [-(1-λ)^n*qa2 - λ^n*qb2 + qa2] * crf_cut * FPEPSI
    let e_rf =
        (-(lp.a_crf_lambda_n * qa2) - lp.b_crf_lambda_n * qb2 + qa2) * crf.crf_cut * four_pi_eps_i;
    // de_rf = n*((1-λ)^(n-1)*qa2 - λ^(n-1)*qb2) * crf_cut * FPEPSI
    let de_rf =
        n * (lp.a_crf_lambda_n_1 * qa2 - lp.b_crf_lambda_n_1 * qb2) * crf.crf_cut * four_pi_eps_i;
    (e_rf, de_rf)
}

/// The soft reaction-field term of an excluded perturbed pair at squared distance `r2`, as the
/// excluded-pair correction computes it: gromosXX `rf_soft_interaction` without the self-term
/// correction, with the softened `crf/2 / sqrt(cut2soft³)` formed as
/// `crf_2cut3i · sqrt((cut2/cut2soft)³)`. Returns `(force, e_crf, de_crf)`, the force a scalar
/// along the pair vector.
///
/// spec: nonbonded.perturbed_rf_excluded_pair
#[inline]
pub fn perturbed_rf_excluded_pair(
    r2: f64,
    a_q: f64,
    b_q: f64,
    alpha_crf: f64,
    lp: &PerturbedLambdaParams,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
) -> (f64, f64, f64) {
    // only the quadratic crf_2cut3i term is softened, not the crf_cut constant
    let cut2 = crf.cutoff_sq;
    let a_cut2soft = cut2 + alpha_crf * lp.b_crfs_lambda2;
    let b_cut2soft = cut2 + alpha_crf * lp.a_crfs_lambda2;
    let a_crf_2cut3i = crf.crf_2cut3i * (cut2 / a_cut2soft).powi(3).sqrt();
    let b_crf_2cut3i = crf.crf_2cut3i * (cut2 / b_cut2soft).powi(3).sqrt();
    let a_crf_cut3i = 2.0 * a_crf_2cut3i;
    let b_crf_cut3i = 2.0 * b_crf_2cut3i;
    let a_crf_pert = 3.0 * a_crf_2cut3i / a_cut2soft;
    let b_crf_pert = 3.0 * b_crf_2cut3i / b_cut2soft;

    // no 1/r term for excluded pairs
    let a_e_crf = a_q * (-a_crf_2cut3i * r2 - crf.crf_cut);
    let b_e_crf = b_q * (-b_crf_2cut3i * r2 - crf.crf_cut);

    let e_crf = (lp.a_crf_lambda_n * a_e_crf + lp.b_crf_lambda_n * b_e_crf) * four_pi_eps_i;
    let force = (lp.a_crf_lambda_n * a_q * a_crf_cut3i + lp.b_crf_lambda_n * b_q * b_crf_cut3i)
        * four_pi_eps_i;

    let n_exp = lp.lambda_exp as f64;
    let de_crf = ((lp.a_crf_lambda_n * a_q * lp.b_crfs_lambda * a_crf_pert
        - lp.b_crf_lambda_n * b_q * lp.a_crfs_lambda * b_crf_pert)
        * r2
        * alpha_crf
        + n_exp * (lp.b_crf_lambda_n_1 * b_e_crf - lp.a_crf_lambda_n_1 * a_e_crf))
        * four_pi_eps_i;
    (force, e_crf, de_crf)
}

/// Correct excluded-pair RF for pairs where ≥1 atom is perturbed.
///
/// Regular code computed state-A CRF for ALL exclusions.  This adds
/// the delta (perturbed_RF − stateA_RF) for each affected pair.
///
/// With more than one energy group each pair also goes to cell `[eg(i)][eg(j)]` of `sink`
/// (gromosXX `perturbed_nonbonded_innerloop.cc:3652-3653`); with `precalc`, its
/// `rf_soft_interaction_ext` at every λ_k (`:3659-3687`, `:3781-3805`).
pub fn perturbed_excluded_correction<BC: BoundaryCondition>(
    exclusions: &[Vec<usize>],
    pos: &[Vec3],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    n_solute: usize,
    four_pi_eps_i: f64,
    out: &mut PertNBCorrection,
    egroups: EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: Option<&mut PrecalcNonbonded>,
) {
    let mut none = PrecalcNonbonded::default();
    let pc = precalc.unwrap_or(&mut none);
    with_groups!(egroups, G => excluded_correction_impl::<BC, G>(
        exclusions,
        pos,
        pert,
        charges,
        crf,
        lps,
        bc,
        n_solute,
        four_pi_eps_i,
        out,
        &egroups,
        sink,
        pc,
    ))
}

fn excluded_correction_impl<BC: BoundaryCondition, const GROUPS: bool>(
    exclusions: &[Vec<usize>],
    pos: &[Vec3],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    n_solute: usize,
    four_pi_eps_i: f64,
    out: &mut PertNBCorrection,
    egroups: &EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: &mut PrecalcNonbonded,
) {
    for i in 0..n_solute {
        let pi = pert[i].as_ref();
        for &j in &exclusions[i] {
            if j <= i {
                continue;
            }
            let pj = pert[j].as_ref();
            if pi.is_none() && pj.is_none() {
                continue;
            }
            let (k, lp) = pair_params::<GROUPS>(lps, egroups, i, j);

            let a_q_i = pi.map_or(charges[i], |p| p.a_charge);
            let b_q_i = pi.map_or(charges[i], |p| p.b_charge);
            let a_q_j = pj.map_or(charges[j], |p| p.a_charge);
            let b_q_j = pj.map_or(charges[j], |p| p.b_charge);
            let alpha_crf = match (pi, pj) {
                (Some(pi), Some(pj)) => (pi.alpha_crf + pj.alpha_crf) * 0.5,
                (Some(p), None) | (None, Some(p)) => p.alpha_crf,
                _ => unreachable!(),
            };

            let r = bc.nearest_image(pos[i], pos[j]);
            let r2 = r.length_squared();
            let a_q = a_q_i * a_q_j;
            let b_q = b_q_i * b_q_j;

            let (f_p, e_crf_p, de_crf) =
                perturbed_rf_excluded_pair(r2, a_q, b_q, alpha_crf, lp, crf, four_pi_eps_i);

            // State-A: what the regular excluded loop computed (same formula, no 1/r)
            let e_crf_a = a_q * four_pi_eps_i * (-crf.crf_2cut3i * r2 - crf.crf_cut);
            let f_a = a_q * four_pi_eps_i * crf.crf_cut3i;

            let df = r * (f_p - f_a);
            out.forces[i] += df;
            out.forces[j] -= df;
            add_outer(&mut out.virial, r, r * f_p);
            add_outer(&mut out.virial_a, r, r * f_a);
            out.delta_e_crf += e_crf_p - e_crf_a;
            out.dhdl += de_crf;
            out.dhdl_crf += de_crf;
            if GROUPS {
                sink.book(k, 0.0, e_crf_p - e_crf_a, 0.0, de_crf);
            }
            if precalc.is_active() {
                precalc.add_rf_ext(k, r2, a_q, b_q, alpha_crf, 1.0, crf, four_pi_eps_i);
            }
        }
    }
}

/// Correct 1-4 nonbonded contributions for pairs involving ≥1 perturbed atom.
///
/// Regular loop used state-A cs6/cs12 and charges.  Applies
/// delta = (1-λ)^n*stateA + λ^n*stateB − stateA for each perturbed pair.
///
/// With more than one energy group each pair also goes to cell `[eg(i)][eg(j)]` of `sink`
/// (gromosXX `perturbed_nonbonded_innerloop.cc:2257-2271`); with `precalc`, its
/// `lj_crf_soft_interaction_ext` with cs6/cs12 at every λ_k (`:2183-2236`, `:2407-2460`; the
/// ext call takes no Coulomb scaling).
pub fn perturbed_one_four_correction<BC: BoundaryCondition>(
    one_four_pairs: &[Vec<usize>],
    pos: &[Vec3],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    iac: &[u32],
    lj: &LJParamMatrix,
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    four_pi_eps_i: f64,
    out: &mut PertNBCorrection,
    egroups: EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: Option<&mut PrecalcNonbonded>,
) {
    let mut none = PrecalcNonbonded::default();
    let pc = precalc.unwrap_or(&mut none);
    with_groups!(egroups, G => one_four_correction_impl::<BC, G>(
        one_four_pairs,
        pos,
        pert,
        charges,
        iac,
        lj,
        crf,
        lps,
        bc,
        four_pi_eps_i,
        out,
        &egroups,
        sink,
        pc,
    ))
}

fn one_four_correction_impl<BC: BoundaryCondition, const GROUPS: bool>(
    one_four_pairs: &[Vec<usize>],
    pos: &[Vec3],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    iac: &[u32],
    lj: &LJParamMatrix,
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    four_pi_eps_i: f64,
    out: &mut PertNBCorrection,
    egroups: &EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: &mut PrecalcNonbonded,
) {
    for i in 0..one_four_pairs.len() {
        let pi = pert[i].as_ref();
        for &j in &one_four_pairs[i] {
            if j <= i {
                continue;
            }
            let pj = pert[j].as_ref();
            if pi.is_none() && pj.is_none() {
                continue;
            }
            let (k, lp) = pair_params::<GROUPS>(lps, egroups, i, j);

            let a_iac_i = pi.map_or(iac[i] as usize, |p| p.a_iac);
            let b_iac_i = pi.map_or(iac[i] as usize, |p| p.b_iac);
            let a_iac_j = pj.map_or(iac[j] as usize, |p| p.a_iac);
            let b_iac_j = pj.map_or(iac[j] as usize, |p| p.b_iac);
            let a_q_i = pi.map_or(charges[i], |p| p.a_charge);
            let b_q_i = pi.map_or(charges[i], |p| p.b_charge);
            let a_q_j = pj.map_or(charges[j], |p| p.a_charge);
            let b_q_j = pj.map_or(charges[j], |p| p.b_charge);

            let lj_a = lj.get(a_iac_i, a_iac_j);
            let lj_b = lj.get(b_iac_i, b_iac_j);
            let a_q = a_q_i * a_q_j;
            let b_q = b_q_i * b_q_j;

            let (alpha_lj, alpha_crf) = match (pi, pj) {
                (Some(pi), Some(pj)) => (
                    (pi.alpha_lj + pj.alpha_lj) * 0.5,
                    (pi.alpha_crf + pj.alpha_crf) * 0.5,
                ),
                (Some(p), None) | (None, Some(p)) => (p.alpha_lj, p.alpha_crf),
                _ => unreachable!(),
            };

            let r = bc.nearest_image(pos[i], pos[j]);
            let r2 = r.length_squared();
            if r2 < 1e-10 {
                continue;
            }
            let inv_r2 = 1.0 / r2;
            let inv_r6 = inv_r2 * inv_r2 * inv_r2;
            let inv_r = inv_r2.sqrt();

            // State-A 1-4 what the regular one_four_interaction_loop computed (no soft-core)
            let e_lj_a = (lj_a.cs12 * inv_r6 - lj_a.cs6) * inv_r6;
            let f_lj_a = (12.0 * lj_a.cs12 * inv_r6 - 6.0 * lj_a.cs6) * inv_r6 * inv_r2;
            let e_crf_a = a_q * four_pi_eps_i * (inv_r - crf.crf_2cut3i * r2 - crf.crf_cut);
            let f_crf_a = a_q * four_pi_eps_i * (inv_r * inv_r2 + crf.crf_cut3i);

            // Perturbed 1-4: call full dual-topology kernel with cs6/cs12 (faithful to
            // GROMOS eds_pert_lj_crf_interaction which uses cs6/cs12 for 1-4 pairs)
            let (f_p, e_lj_p, e_crf_p, de_lj, de_crf) = perturbed_lj_crf_interaction(
                r,
                lj_a.cs6,
                lj_a.cs12,
                lj_b.cs6,
                lj_b.cs12,
                a_q,
                b_q,
                alpha_lj,
                alpha_crf,
                lp,
                crf,
                four_pi_eps_i,
            );

            let df = r * (f_p - (f_lj_a + f_crf_a));
            out.forces[i] += df;
            out.forces[j] -= df;
            add_outer(&mut out.virial, r, r * f_p);
            add_outer(&mut out.virial_a, r, r * (f_lj_a + f_crf_a));
            out.delta_e_lj += e_lj_p - e_lj_a;
            out.delta_e_crf += e_crf_p - e_crf_a;
            out.dhdl += de_lj + de_crf;
            out.dhdl_lj += de_lj;
            out.dhdl_crf += de_crf;
            if GROUPS {
                sink.book(k, e_lj_p - e_lj_a, e_crf_p - e_crf_a, de_lj, de_crf);
            }
            if precalc.is_active() {
                precalc.add_lj_crf_ext(
                    k,
                    r2,
                    (lj_a.cs6, lj_a.cs12, lj_b.cs6, lj_b.cs12),
                    a_q,
                    b_q,
                    alpha_lj,
                    alpha_crf,
                    crf,
                    four_pi_eps_i,
                );
            }
        }
    }
}

/// gromosXX `rf_soft_interaction` (`perturbed_nonbonded_term.cc:718-765`) without the self-term
/// correction: the soft reaction-field term of an excluded perturbed pair at the current λ.
/// `r2` is the squared distance.
///
/// Returns `(force, e_rf, de_rf)`, the force a scalar along the pair vector.
///
/// spec: nonbonded.rf_soft_interaction
#[inline]
pub fn rf_soft_interaction(
    r2: f64,
    a_q: f64,
    b_q: f64,
    alpha_crf: f64,
    lp: &PerturbedLambdaParams,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
) -> (f64, f64, f64) {
    let cut2 = crf.cutoff_sq;
    let a_cut2soft = cut2 + alpha_crf * lp.b_crfs_lambda2;
    let b_cut2soft = cut2 + alpha_crf * lp.a_crfs_lambda2;
    // crf/2 / sqrt(cut2soft³), with crf/2 = crf_2cut3i·cut³ as in `perturbed_lj_crf_interaction`
    let cut3 = cut2 * cut2.sqrt();
    let a_crf_2cut3i = crf.crf_2cut3i * cut3 / (a_cut2soft * a_cut2soft * a_cut2soft).sqrt();
    let b_crf_2cut3i = crf.crf_2cut3i * cut3 / (b_cut2soft * b_cut2soft * b_cut2soft).sqrt();
    let a_crf_pert = 3.0 * a_crf_2cut3i / a_cut2soft;
    let b_crf_pert = 3.0 * b_crf_2cut3i / b_cut2soft;

    let force = (lp.a_crf_lambda_n * a_q * 2.0 * a_crf_2cut3i
        + lp.b_crf_lambda_n * b_q * 2.0 * b_crf_2cut3i)
        * four_pi_eps_i;
    let a_e_rf = a_q * (-a_crf_2cut3i * r2 - crf.crf_cut);
    let b_e_rf = b_q * (-b_crf_2cut3i * r2 - crf.crf_cut);
    let e_rf = (lp.a_crf_lambda_n * a_e_rf + lp.b_crf_lambda_n * b_e_rf) * four_pi_eps_i;
    let de_rf = ((lp.a_crf_lambda_n * a_q * lp.b_crfs_lambda * a_crf_pert
        - lp.b_crf_lambda_n * b_q * lp.a_crfs_lambda * b_crf_pert)
        * r2
        * alpha_crf
        + (lp.b_crf_lambda_n_1 * b_e_rf - lp.a_crf_lambda_n_1 * a_e_rf) * lp.lambda_exp as f64)
        * four_pi_eps_i;
    (force, e_rf, de_rf)
}

/// The `PERTATOMPAIR` pairs (gromosXX `Perturbed_Nonbonded_Pair`,
/// `perturbed_nonbonded_pair.cc:118-1612`). [`Topology::apply_perturbation`] took them out of
/// the exclusion and 1-4 lists and off the pairlist, so this is the only place they are
/// computed: each state by its own interaction type (excluded, normal, 1-4), with the A/B
/// parameters of perturbed atoms and soft-core, or λ-scaled plain terms when neither atom is
/// perturbed, then summed `A + B`.
///
/// The force, the energies and the dH/dλ go to `out`; with more than one energy group also to
/// cell `[eg(i)][eg(j)]` of `sink` (`:1583-1611`); with `precalc` the pair's ext energies at
/// every λ_k (`:343-365, 410-420, 479-525, 600-612, 685-731, 803-815` and the B mirrors).
/// `coulomb_scaling` is the 1-4 Coulomb scaling of a 1-4 state with a perturbed atom.
///
/// [`Topology::apply_perturbation`]: newen_core::Topology::apply_perturbation
pub fn perturbed_atom_pair_interactions<BC: BoundaryCondition>(
    atom_pairs: &[PerturbedAtomPair],
    pos: &[Vec3],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    iac: &[u32],
    lj: &LJParamMatrix,
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    four_pi_eps_i: f64,
    rf_excluded: bool,
    coulomb_scaling: f64,
    out: &mut PertNBCorrection,
    egroups: EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: Option<&mut PrecalcNonbonded>,
) {
    let mut none = PrecalcNonbonded::default();
    let pc = precalc.unwrap_or(&mut none);
    with_groups!(egroups, G => atom_pair_impl::<BC, G>(
        atom_pairs,
        pos,
        pert,
        charges,
        iac,
        lj,
        crf,
        lps,
        bc,
        four_pi_eps_i,
        rf_excluded,
        coulomb_scaling,
        out,
        &egroups,
        sink,
        pc,
    ))
}

/// The parameters of one end state of a perturbed atom pair.
struct PairState {
    kind: PairInteraction,
    /// state B (else A)
    b: bool,
    lj: LJParameters,
    q: f64,
}

/// One end state of a `PERTATOMPAIR` pair — the body of gromosXX's `switch (A_type)` /
/// `switch (B_type)` for `lj_crf_func` (`perturbed_nonbonded_pair.cc:328-868`, `:875-1414`).
/// Returns `[force, e_lj, e_crf, de_lj, de_crf]`, the force a scalar along `r`, already
/// λ-weighted.
#[inline]
fn atom_pair_state(
    s: &PairState,
    perturbed: bool,
    r: Vec3,
    alpha_lj: f64,
    alpha_crf: f64,
    lp: &PerturbedLambdaParams,
    crf: &CRFParameters,
    four_pi_eps_i: f64,
    rf_excluded: bool,
    coulomb_scaling: f64,
    cell: usize,
    precalc: &mut PrecalcNonbonded,
) -> [f64; 5] {
    let r2 = r.length_squared();
    let n = lp.lambda_exp as f64;
    // λ weight of the state, its derivative with gromosXX's sign (−n(1−λ)^(n−1), +nλ^(n−1)),
    // and the state's column of the ext record
    let (lj_n, crf_n, lj_d, crf_d) = if s.b {
        (
            lp.b_lj_lambda_n,
            lp.b_crf_lambda_n,
            n * lp.b_lj_lambda_n_1,
            n * lp.b_crf_lambda_n_1,
        )
    } else {
        (
            lp.a_lj_lambda_n,
            lp.a_crf_lambda_n,
            -n * lp.a_lj_lambda_n_1,
            -n * lp.a_crf_lambda_n_1,
        )
    };
    let ext = |e_lj: f64, e_crf: f64| {
        if s.b {
            [0.0, e_lj, 0.0, e_crf]
        } else {
            [e_lj, 0.0, e_crf, 0.0]
        }
    };
    let (a_q, b_q) = if s.b { (0.0, s.q) } else { (s.q, 0.0) };

    match (s.kind, perturbed) {
        (PairInteraction::Excluded, _) if !rf_excluded => [0.0; 5],
        (PairInteraction::Excluded, true) => {
            // pair.cc:338-365, 882-909
            let (f, e_crf, de_crf) =
                rf_soft_interaction(r2, a_q, b_q, alpha_crf, lp, crf, four_pi_eps_i);
            if precalc.is_active() {
                precalc.add_rf_ext(cell, r2, a_q, b_q, alpha_crf, 1.0, crf, four_pi_eps_i);
            }
            [f, 0.0, e_crf, 0.0, de_crf]
        },
        (PairInteraction::Excluded, false) => {
            // pair.cc:402-420, 931-948: gromosXX `rf_interaction`, scaled here
            let qe = s.q * four_pi_eps_i;
            let f = qe * crf.crf_cut3i;
            let e_crf = qe * (-crf.crf_2cut3i * r2 - crf.crf_cut);
            if precalc.is_active() {
                precalc.add_constant(cell, ext(0.0, e_crf));
            }
            [crf_n * f, 0.0, crf_n * e_crf, 0.0, crf_d * e_crf]
        },
        (kind, true) => {
            // pair.cc:463-525, 670-731, 986-1047, 1194-1253
            let (c6, c12) = match kind {
                PairInteraction::OneFour => (s.lj.cs6, s.lj.cs12),
                _ => (s.lj.c6, s.lj.c12),
            };
            let (a_lj, b_lj) = if s.b {
                ((0.0, 0.0), (c6, c12))
            } else {
                ((c6, c12), (0.0, 0.0))
            };
            if precalc.is_active() {
                precalc.add_lj_crf_ext(
                    cell,
                    r2,
                    (a_lj.0, a_lj.1, b_lj.0, b_lj.1),
                    a_q,
                    b_q,
                    alpha_lj,
                    alpha_crf,
                    crf,
                    four_pi_eps_i,
                );
            }
            // the outputs in their own slots (D1), Coulomb scaling 1 for a normal state and the
            // 1-4 one for a 1-4 state (perturbed_nonbonded_innerloop.cc:2102)
            let scaling = match kind {
                PairInteraction::OneFour => coulomb_scaling,
                _ => 1.0,
            };
            let (f, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction_scaled(
                r,
                a_lj.0,
                a_lj.1,
                b_lj.0,
                b_lj.1,
                a_q,
                b_q,
                alpha_lj,
                alpha_crf,
                lp,
                crf,
                four_pi_eps_i,
                scaling,
            );
            [f, e_lj, e_crf, de_lj, de_crf]
        },
        (kind, false) => {
            // pair.cc:582-612, 788-815, 1105-1131, 1319-1345
            let (c6, c12) = match kind {
                PairInteraction::OneFour => (s.lj.cs6, s.lj.cs12),
                _ => (s.lj.c6, s.lj.c12),
            };
            let (f, e_lj, e_crf) = lj_crf_interaction(r, c6, c12, s.q * four_pi_eps_i, crf);
            if precalc.is_active() {
                precalc.add_constant(cell, ext(e_lj, e_crf));
            }
            // both parts of the force take the LJ λ (pair.cc:592, 798, 1115, 1329)
            [
                lj_n * f,
                lj_n * e_lj,
                crf_n * e_crf,
                lj_d * e_lj,
                crf_d * e_crf,
            ]
        },
    }
}

fn atom_pair_impl<BC: BoundaryCondition, const GROUPS: bool>(
    atom_pairs: &[PerturbedAtomPair],
    pos: &[Vec3],
    pert: &[Option<PertAtomInfo>],
    charges: &[f64],
    iac: &[u32],
    lj: &LJParamMatrix,
    crf: &CRFParameters,
    lps: &[PerturbedLambdaParams],
    bc: &BC,
    four_pi_eps_i: f64,
    rf_excluded: bool,
    coulomb_scaling: f64,
    out: &mut PertNBCorrection,
    egroups: &EnergyGroups<'_>,
    sink: &mut GroupEnergies,
    precalc: &mut PrecalcNonbonded,
) {
    for ap in atom_pairs {
        let (i, j) = (ap.i, ap.j);
        let (k, lp) = pair_params::<GROUPS>(lps, egroups, i, j);
        let pi = pert.get(i).and_then(Option::as_ref);
        let pj = pert.get(j).and_then(Option::as_ref);

        // pair.cc:179-298
        let a_iac_i = pi.map_or(iac[i] as usize, |p| p.a_iac);
        let a_iac_j = pj.map_or(iac[j] as usize, |p| p.a_iac);
        let b_iac_i = pi.map_or(iac[i] as usize, |p| p.b_iac);
        let b_iac_j = pj.map_or(iac[j] as usize, |p| p.b_iac);
        let a_q = pi.map_or(charges[i], |p| p.a_charge) * pj.map_or(charges[j], |p| p.a_charge);
        let b_q = pi.map_or(charges[i], |p| p.b_charge) * pj.map_or(charges[j], |p| p.b_charge);
        let (alpha_lj, alpha_crf) = match (pi, pj) {
            (Some(pi), Some(pj)) => (
                (pi.alpha_lj + pj.alpha_lj) / 2.0,
                (pi.alpha_crf + pj.alpha_crf) / 2.0,
            ),
            (Some(p), None) | (None, Some(p)) => (p.alpha_lj, p.alpha_crf),
            (None, None) => (0.0, 0.0),
        };
        let perturbed = pi.is_some() || pj.is_some();

        let r = bc.nearest_image(pos[i], pos[j]);
        let state = |s: PairState, precalc: &mut PrecalcNonbonded| {
            atom_pair_state(
                &s,
                perturbed,
                r,
                alpha_lj,
                alpha_crf,
                lp,
                crf,
                four_pi_eps_i,
                rf_excluded,
                coulomb_scaling,
                k,
                precalc,
            )
        };
        let a = state(
            PairState {
                kind: ap.a_type,
                b: false,
                lj: lj.get(a_iac_i, a_iac_j),
                q: a_q,
            },
            precalc,
        );
        let b = state(
            PairState {
                kind: ap.b_type,
                b: true,
                lj: lj.get(b_iac_i, b_iac_j),
                q: b_q,
            },
            precalc,
        );
        // pair.cc:1532-1539
        let [f, e_lj, e_crf, de_lj, de_crf] = std::array::from_fn(|c| b[c] + a[c]);

        let fv = r * f;
        out.forces[i] += fv;
        out.forces[j] -= fv;
        let (rv, fa) = ([r.x, r.y, r.z], [fv.x, fv.y, fv.z]);
        for (x, row) in out.virial.iter_mut().enumerate() {
            for (y, v) in row.iter_mut().enumerate() {
                *v += rv[x] * fa[y];
            }
        }
        out.delta_e_lj += e_lj;
        out.delta_e_crf += e_crf;
        out.dhdl += de_lj + de_crf;
        out.dhdl_lj += de_lj;
        out.dhdl_crf += de_crf;
        if GROUPS {
            sink.book(k, e_lj, e_crf, de_lj, de_crf);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use newen_core::math::Vec3;

    #[test]
    fn test_perturbed_lambda_params() {
        // Test lambda parameter calculation
        let lambda = 0.5;
        let lambda_soft = 0.5;
        let derivative = 1.0;
        let n = 2;

        let params = PerturbedLambdaParams::from_lambda(
            lambda,
            lambda_soft, // lj_lambda, ljs_lambda
            lambda,
            lambda_soft, // crf_lambda, crfs_lambda
            derivative,
            derivative, // lj_lambda_derivative, ljs_lambda_derivative
            derivative,
            derivative, // crf_lambda_derivative, crfs_lambda_derivative
            n,
        );

        // Check λ^n calculations
        assert_relative_eq!(params.b_lj_lambda_n, 0.25, epsilon = 1e-10); // 0.5²
        assert_relative_eq!(params.a_lj_lambda_n, 0.25, epsilon = 1e-10); // (1-0.5)²

        // Check λ^(n-1) calculations
        assert_relative_eq!(params.b_lj_lambda_n_1, 0.5, epsilon = 1e-10); // 0.5¹ * 1.0
        assert_relative_eq!(params.a_lj_lambda_n_1, 0.5, epsilon = 1e-10); // (1-0.5)¹ * 1.0

        // Check soft-core lambda squared
        assert_relative_eq!(params.b_ljs_lambda2, 0.25, epsilon = 1e-10); // 0.5²
        assert_relative_eq!(params.a_ljs_lambda2, 0.25, epsilon = 1e-10); // (1-0.5)²
    }

    #[test]
    fn test_perturbed_interaction_at_lambda_0() {
        // At λ=0, should get pure state A
        let r = Vec3::new(1.0, 0.0, 0.0);

        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.002; // Different from state A
        let b_c12 = 0.0002;

        let a_q = 0.25; // q_i * q_j for state A
        let b_q = 0.0; // No charge in state B

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let lambda_params = PerturbedLambdaParams::from_lambda(
            0.0, 0.0, // lj_lambda=0, ljs_lambda=0
            0.0, 0.0, // crf_lambda=0, crfs_lambda=0
            1.0, 1.0, // derivatives
            1.0, 1.0, 1, // n=1
        );

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let (_force, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction(
            r,
            a_c6,
            a_c12,
            b_c6,
            b_c12,
            a_q,
            b_q,
            alpha_lj,
            alpha_crf,
            &lambda_params,
            &crf,
            newen_core::units::four_pi_eps_i,
        );

        // At λ=0: should be dominated by state A
        // Energy should be finite
        assert!(e_lj.is_finite(), "LJ energy should be finite");
        assert!(e_crf.is_finite(), "CRF energy should be finite");

        // Derivatives should be finite
        assert!(de_lj.is_finite(), "dH/dλ for LJ should be finite");
        assert!(de_crf.is_finite(), "dH/dλ for CRF should be finite");

        println!(
            "λ=0: E_lj={:.6}, E_crf={:.6}, dE_lj/dλ={:.6}, dE_crf/dλ={:.6}",
            e_lj, e_crf, de_lj, de_crf
        );
    }

    /// The ext kernels' sign convention (`A_de = +∂A/∂λ`, `B_de = −∂B/∂λ`) and their relation
    /// to the λ-weighted kernel at λ_k = λ, NLAM = 1, `L = λ`:
    /// `E = (1−λ)·A + λ·B`, `dE/dλ = (1−λ)·A_de − λ·B_de + B − A`.
    #[test]
    fn ext_kernels_are_the_unweighted_end_states_of_the_soft_core_kernel() {
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let fpe = newen_core::units::four_pi_eps_i;
        let (r2, a_q, b_q, alj, acrf) = (0.09, -0.3, 0.2, 0.5, 0.6);
        let (a_c6, a_c12, b_c6, b_c12) = (2.6e-3, 2.7e-6, 1.1e-3, 0.9e-6);
        let ext = |lam: f64| {
            lj_crf_soft_interaction_ext(
                r2, a_c6, a_c12, b_c6, b_c12, a_q, b_q, alj, acrf, &crf, fpe, lam,
            )
        };
        let lam = 0.3;
        let (e, de) = ext(lam);
        let eps = 1e-6;
        let (ep, _) = ext(lam + eps);
        let (em, _) = ext(lam - eps);
        for (c, sign) in [(0, 1.0), (1, -1.0), (2, 1.0), (3, -1.0)] {
            let fd = (ep[c] - em[c]) / (2.0 * eps);
            assert!(
                (sign * fd - de[c]).abs() <= 1e-6 * fd.abs().max(1e-3),
                "column {c}: fd {fd} vs {}",
                de[c]
            );
        }
        let lp = PerturbedLambdaParams::from_lambda(lam, lam, lam, lam, 1.0, 1.0, 1.0, 1.0, 1);
        let r = Vec3::new(r2.sqrt(), 0.0, 0.0);
        let (_, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction(
            r, a_c6, a_c12, b_c6, b_c12, a_q, b_q, alj, acrf, &lp, &crf, fpe,
        );
        let close = |x: f64, y: f64| (x - y).abs() <= 1e-12 * x.abs().max(1.0);
        assert!(close(e_lj, (1.0 - lam) * e[0] + lam * e[1]));
        assert!(close(e_crf, (1.0 - lam) * e[2] + lam * e[3]));
        assert!(close(
            de_lj,
            (1.0 - lam) * de[0] - lam * de[1] + e[1] - e[0]
        ));
        assert!(close(
            de_crf,
            (1.0 - lam) * de[2] - lam * de[3] + e[3] - e[2]
        ));

        // the excluded-pair RF term, likewise
        let rf = |lam: f64| rf_soft_interaction_ext(r2, a_q, b_q, acrf, &crf, fpe, lam);
        let (_, d) = rf(lam);
        let ((p, _), (m, _)) = (rf(lam + eps), rf(lam - eps));
        assert!(((p[0] - m[0]) / (2.0 * eps) - d[0]).abs() <= 1e-6);
        assert!((-(p[1] - m[1]) / (2.0 * eps) - d[1]).abs() <= 1e-6);
    }

    /// What one `PERTATOMPAIR` pair (atom 0 perturbed, atom 1 not) books through
    /// [`perturbed_atom_pair_interactions`].
    struct PairOut {
        force: Vec3,
        e_lj: f64,
        e_crf: f64,
        dhdl: f64,
    }

    const PAIR_LJ: [[LJParameters; 2]; 2] = [
        [
            LJParameters {
                c6: 2.6e-3,
                c12: 2.7e-6,
                cs6: 2.1e-3,
                cs12: 1.9e-6,
            },
            LJParameters {
                c6: 1.8e-3,
                c12: 1.6e-6,
                cs6: 1.2e-3,
                cs12: 1.1e-6,
            },
        ],
        [
            LJParameters {
                c6: 1.8e-3,
                c12: 1.6e-6,
                cs6: 1.2e-3,
                cs12: 1.1e-6,
            },
            LJParameters {
                c6: 1.1e-3,
                c12: 0.9e-6,
                cs6: 0.8e-3,
                cs12: 0.7e-6,
            },
        ],
    ];
    const PAIR_Q: (f64, f64, f64) = (0.4, -0.3, -0.5); // atom 0 in A, in B; atom 1

    fn pair_crf() -> CRFParameters {
        CRFParameters::new(1.4, 1.0, 61.0, 0.0)
    }

    fn one_pair(
        types: (PairInteraction, PairInteraction),
        pos0: Vec3,
        lam: f64,
        n: i32,
        coulomb_scaling: f64,
    ) -> PairOut {
        let pos = [pos0, Vec3::ZERO];
        let pert = [
            Some(PertAtomInfo {
                a_iac: 0,
                b_iac: 1,
                a_charge: PAIR_Q.0,
                b_charge: PAIR_Q.1,
                alpha_lj: 0.5,
                alpha_crf: 0.6,
            }),
            None,
        ];
        let lj = LJParamMatrix::from_nested(&PAIR_LJ.map(Vec::from));
        let lps = [PerturbedLambdaParams::from_lambda(
            lam, lam, lam, lam, 1.0, 1.0, 1.0, 1.0, n,
        )];
        let pairs = [PerturbedAtomPair {
            i: 0,
            j: 1,
            a_type: types.0,
            b_type: types.1,
        }];
        let mut out = PertNBCorrection::new(2);
        perturbed_atom_pair_interactions(
            &pairs,
            &pos,
            &pert,
            &[PAIR_Q.0, PAIR_Q.2],
            &[0, 0],
            &lj,
            &pair_crf(),
            &lps,
            &newen_core::math::Vacuum,
            newen_core::units::four_pi_eps_i,
            true,
            coulomb_scaling,
            &mut out,
            EnergyGroups::single(),
            &mut GroupEnergies::new(1),
            None,
        );
        PairOut {
            force: out.forces[0],
            e_lj: out.delta_e_lj,
            e_crf: out.delta_e_crf,
            dhdl: out.dhdl,
        }
    }

    const PAIR_POS: Vec3 = Vec3::new(0.21, -0.17, 0.23);

    /// At λ = 0 only the A state acts and the soft core vanishes, so the intended kernel of a
    /// normal (or 1-4) A state is the plain LJ+CRF kernel; likewise the B state at λ = 1.
    #[test]
    fn atom_pair_kernel_matches_plain_kernel_at_endpoints() {
        use PairInteraction::{Excluded, Normal, OneFour};
        let fpe = newen_core::units::four_pi_eps_i;
        let crf = pair_crf();
        let close = |x: f64, y: f64| (x - y).abs() <= 1e-12 * x.abs().max(1e-6);
        let cases = [
            (
                (Normal, OneFour),
                0.0,
                PAIR_LJ[0][0].c6,
                PAIR_LJ[0][0].c12,
                PAIR_Q.0,
            ),
            (
                (Normal, Excluded),
                0.0,
                PAIR_LJ[0][0].c6,
                PAIR_LJ[0][0].c12,
                PAIR_Q.0,
            ),
            (
                (OneFour, Normal),
                0.0,
                PAIR_LJ[0][0].cs6,
                PAIR_LJ[0][0].cs12,
                PAIR_Q.0,
            ),
            (
                (OneFour, Normal),
                1.0,
                PAIR_LJ[1][0].c6,
                PAIR_LJ[1][0].c12,
                PAIR_Q.1,
            ),
            (
                (Excluded, Normal),
                1.0,
                PAIR_LJ[1][0].c6,
                PAIR_LJ[1][0].c12,
                PAIR_Q.1,
            ),
            (
                (Normal, OneFour),
                1.0,
                PAIR_LJ[1][0].cs6,
                PAIR_LJ[1][0].cs12,
                PAIR_Q.1,
            ),
        ];
        for (types, lam, c6, c12, q0) in cases {
            let got = one_pair(types, PAIR_POS, lam, 1, 1.0);
            let (f, e_lj, e_crf) = lj_crf_interaction(PAIR_POS, c6, c12, q0 * PAIR_Q.2 * fpe, &crf);
            let want = PAIR_POS * f;
            let tag = format!("{types:?} at λ = {lam}");
            assert!(close(got.e_lj, e_lj), "{tag}: e_lj {} vs {e_lj}", got.e_lj);
            assert!(
                close(got.e_crf, e_crf),
                "{tag}: e_crf {} vs {e_crf}",
                got.e_crf
            );
            for (a, b) in [got.force.x, got.force.y, got.force.z]
                .into_iter()
                .zip([want.x, want.y, want.z])
            {
                assert!(close(a, b), "{tag}: force {a} vs {b}");
            }
        }
    }

    /// The intended pair term is consistent: its force is −∂V/∂r and its dH/dλ is ∂V/∂λ
    /// (central differences), a 1-4 state with a non-unit Coulomb scaling included.
    #[test]
    fn atom_pair_kernel_force_and_dhdl_are_derivatives_of_its_energy() {
        use PairInteraction::{Excluded, Normal, OneFour};
        let v = |o: &PairOut| o.e_lj + o.e_crf;
        let check = |what: &str, analytic: f64, fd: f64| {
            assert!(
                (analytic - fd).abs() <= 1e-6 * analytic.abs().max(1.0),
                "{what}: analytic {analytic} vs finite difference {fd}"
            );
        };
        for types in [
            (Normal, OneFour),
            (OneFour, Normal),
            (OneFour, OneFour),
            (Excluded, Normal),
        ] {
            for scaling in [1.0, 0.5] {
                for n in [1, 2] {
                    for lam in [0.25, 0.5, 0.75] {
                        let run = |p: Vec3, l: f64| one_pair(types, p, l, n, scaling);
                        let tag = format!("{types:?} scaling {scaling} n {n} λ {lam}");
                        let at = run(PAIR_POS, lam);
                        let h = 1e-6;
                        for (d, (axis, f)) in [
                            (Vec3::X, at.force.x),
                            (Vec3::Y, at.force.y),
                            (Vec3::Z, at.force.z),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            let fd = -(v(&run(PAIR_POS + axis * h, lam))
                                - v(&run(PAIR_POS - axis * h, lam)))
                                / (2.0 * h);
                            check(&format!("{tag}: force {d}"), f, fd);
                        }
                        let hl = 1e-6;
                        let fd = (v(&run(PAIR_POS, lam + hl)) - v(&run(PAIR_POS, lam - hl)))
                            / (2.0 * hl);
                        check(&format!("{tag}: dH/dλ"), at.dhdl, fd);
                    }
                }
            }
        }
    }

    #[test]
    fn test_perturbed_interaction_at_lambda_1() {
        // At λ=1, should get pure state B
        let r = Vec3::new(1.0, 0.0, 0.0);

        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.002; // Different from state A
        let b_c12 = 0.0002;

        let a_q = 0.0; // No charge in state A
        let b_q = 0.25; // q_i * q_j for state B

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let lambda_params = PerturbedLambdaParams::from_lambda(
            1.0, 1.0, // lj_lambda=1, ljs_lambda=1
            1.0, 1.0, // crf_lambda=1, crfs_lambda=1
            1.0, 1.0, // derivatives
            1.0, 1.0, 1, // n=1
        );

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let (_force, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction(
            r,
            a_c6,
            a_c12,
            b_c6,
            b_c12,
            a_q,
            b_q,
            alpha_lj,
            alpha_crf,
            &lambda_params,
            &crf,
            newen_core::units::four_pi_eps_i,
        );

        // At λ=1: should be dominated by state B
        assert!(e_lj.is_finite(), "LJ energy should be finite");
        assert!(e_crf.is_finite(), "CRF energy should be finite");
        assert!(de_lj.is_finite(), "dH/dλ for LJ should be finite");
        assert!(de_crf.is_finite(), "dH/dλ for CRF should be finite");

        println!(
            "λ=1: E_lj={:.6}, E_crf={:.6}, dE_lj/dλ={:.6}, dE_crf/dλ={:.6}",
            e_lj, e_crf, de_lj, de_crf
        );
    }

    #[test]
    fn test_perturbed_interaction_lambda_sweep() {
        // Test lambda sweep from 0 to 1
        let r = Vec3::new(0.5, 0.0, 0.0);

        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.002;
        let b_c12 = 0.0003;

        let a_q = 0.25;
        let b_q = -0.25;

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let lambdas = [0.0, 0.25, 0.5, 0.75, 1.0];
        let mut energies = Vec::new();
        let mut derivatives = Vec::new();

        for &lambda in &lambdas {
            let lambda_params = PerturbedLambdaParams::from_lambda(
                lambda, lambda, lambda, lambda, 1.0, 1.0, 1.0, 1.0, 1,
            );

            let (_, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_interaction(
                r,
                a_c6,
                a_c12,
                b_c6,
                b_c12,
                a_q,
                b_q,
                alpha_lj,
                alpha_crf,
                &lambda_params,
                &crf,
                newen_core::units::four_pi_eps_i,
            );

            energies.push(e_lj + e_crf);
            derivatives.push(de_lj + de_crf);

            // All should be finite
            assert!(e_lj.is_finite(), "Energy should be finite at λ={}", lambda);
            assert!(
                de_lj.is_finite(),
                "Derivative should be finite at λ={}",
                lambda
            );
        }

        println!("Lambda sweep:");
        for (i, lambda) in lambdas.iter().enumerate() {
            println!(
                "  λ={:.2}: E={:.6}, dE/dλ={:.6}",
                lambda, energies[i], derivatives[i]
            );
        }

        // Check smoothness: no huge jumps
        for i in 0..energies.len() - 1 {
            let energy_change = (energies[i + 1] - energies[i]).abs();
            assert!(
                energy_change < 100.0,
                "Energy jump too large between λ={} and λ={}",
                lambdas[i],
                lambdas[i + 1]
            );
        }
    }

    #[test]
    fn test_soft_core_prevents_singularity() {
        // Test that soft-core prevents singularities at close distances
        let r = Vec3::new(0.1, 0.0, 0.0); // Very close distance

        // Particle appearing (state A=dummy, state B=real)
        let a_c6 = 0.0; // Dummy particle
        let a_c12 = 0.0;
        let b_c6 = 0.001; // Real particle
        let b_c12 = 0.0001;

        let a_q = 0.0;
        let b_q = 0.25;

        let alpha_lj = 1.0; // Strong soft-core
        let alpha_crf = 1.0;

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
            cutoff_sq: 1.4_f64.powi(2),
        };

        // At λ=0 (particle not yet present), soft-core should prevent singularity
        let lambda_params_0 =
            PerturbedLambdaParams::from_lambda(0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1);

        let (_, e_lj_0, e_crf_0, _, _) = perturbed_lj_crf_interaction(
            r,
            a_c6,
            a_c12,
            b_c6,
            b_c12,
            a_q,
            b_q,
            alpha_lj,
            alpha_crf,
            &lambda_params_0,
            &crf,
            newen_core::units::four_pi_eps_i,
        );

        // At λ=1 (particle fully present), will have high energy but should be finite
        let lambda_params_1 =
            PerturbedLambdaParams::from_lambda(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1);

        let (_, e_lj_1, e_crf_1, _, _) = perturbed_lj_crf_interaction(
            r,
            a_c6,
            a_c12,
            b_c6,
            b_c12,
            a_q,
            b_q,
            alpha_lj,
            alpha_crf,
            &lambda_params_1,
            &crf,
            newen_core::units::four_pi_eps_i,
        );

        // Energies should be finite (no NaN or Inf)
        assert!(
            e_lj_0.is_finite(),
            "Soft-core should prevent LJ singularity at λ=0"
        );
        assert!(
            e_crf_0.is_finite(),
            "Soft-core should prevent CRF singularity at λ=0"
        );
        assert!(e_lj_1.is_finite(), "Energy should be finite at λ=1");
        assert!(e_crf_1.is_finite(), "Energy should be finite at λ=1");

        // At λ=0 with soft-core, energy should be much lower than at λ=1
        println!("Soft-core test at r=0.1 nm:");
        println!("  λ=0: E_lj={:.6}, E_crf={:.6}", e_lj_0, e_crf_0);
        println!("  λ=1: E_lj={:.6}, E_crf={:.6}", e_lj_1, e_crf_1);
    }

    #[test]
    fn test_lambda_exponent() {
        // Test different lambda exponents
        let r = Vec3::new(1.0, 0.0, 0.0);
        let lambda = 0.5;

        // Use different A and B state parameters to see perturbation effect
        let a_c6 = 0.001;
        let a_c12 = 0.0001;
        let b_c6 = 0.002; // Different from A state
        let b_c12 = 0.0002; // Different from A state
        let a_q = 0.1;
        let b_q = -0.1; // Different charge

        let alpha_lj = 0.5;
        let alpha_crf = 0.5;

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431 / 2.0,
            crf_cut3i: 0.364431 / 2.0 / 1.4,
            cutoff_sq: 1.4_f64.powi(2),
        };

        // Test n=1 (linear coupling)
        let params_n1 = PerturbedLambdaParams::from_lambda(
            lambda, lambda, lambda, lambda, 1.0, 1.0, 1.0, 1.0, 1,
        );

        let (_, e_lj_n1, _e_crf_n1, de_lj_n1, _de_crf_n1) = perturbed_lj_crf_interaction(
            r,
            a_c6,
            a_c12,
            b_c6,
            b_c12,
            a_q,
            b_q,
            alpha_lj,
            alpha_crf,
            &params_n1,
            &crf,
            newen_core::units::four_pi_eps_i,
        );

        // Test n=2 (quadratic coupling)
        let params_n2 = PerturbedLambdaParams::from_lambda(
            lambda, lambda, lambda, lambda, 1.0, 1.0, 1.0, 1.0, 2,
        );

        let (_, e_lj_n2, _e_crf_n2, de_lj_n2, _de_crf_n2) = perturbed_lj_crf_interaction(
            r,
            a_c6,
            a_c12,
            b_c6,
            b_c12,
            a_q,
            b_q,
            alpha_lj,
            alpha_crf,
            &params_n2,
            &crf,
            newen_core::units::four_pi_eps_i,
        );

        // Energies should differ due to different λ^n weighting
        assert_ne!(
            e_lj_n1, e_lj_n2,
            "Energy should differ for different exponents"
        );
        assert_ne!(
            de_lj_n1, de_lj_n2,
            "Derivative should differ for different exponents"
        );

        println!("Lambda exponent comparison at λ=0.5:");
        println!("  n=1: E_lj={:.6}, dE_lj/dλ={:.6}", e_lj_n1, de_lj_n1);
        println!("  n=2: E_lj={:.6}, dE_lj/dλ={:.6}", e_lj_n2, de_lj_n2);
    }

    /// Two perturbed atoms on top of each other: with α > 0 and C6 ≠ 0 in both states the
    /// soft-core pair has a finite energy and dV/dλ at r = 0 (a state with C6 = 0 has no
    /// soft-core length, so there it diverges, as in gromosXX), and at NLAM = 1 it is the λ-mix
    /// of the ext kernel's end states.
    #[test]
    fn overlapping_soft_core_pair_is_finite() {
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let fpe = newen_core::units::four_pi_eps_i;
        let (a_c6, a_c12, b_c6, b_c12, a_q, b_q) = (2.6e-3, 2.7e-6, 1.1e-3, 0.9e-6, -0.3, 0.2);
        let (alpha_lj, alpha_crf, lam) = (0.5, 0.5, 0.4);
        let lp = PerturbedLambdaParams::from_lambda(lam, lam, lam, lam, 1.0, 1.0, 1.0, 1.0, 1);
        for r2 in [0.0, 1e-12, 9e-11] {
            let (force, e_lj, e_crf, de_lj, de_crf) = perturbed_lj_crf_pair(
                r2, a_c6, a_c12, b_c6, b_c12, a_q, b_q, alpha_lj, alpha_crf, &lp, &crf, fpe, 1.0,
            );
            for v in [force, e_lj, e_crf, de_lj, de_crf] {
                assert!(v.is_finite(), "r2 = {r2}: {v}");
            }
            assert!(
                e_lj != 0.0 && e_crf != 0.0,
                "r2 = {r2}: the pair is not dropped"
            );
            let (e, _) = lj_crf_soft_interaction_ext(
                r2, a_c6, a_c12, b_c6, b_c12, a_q, b_q, alpha_lj, alpha_crf, &crf, fpe, lam,
            );
            assert_relative_eq!(e_lj, (1.0 - lam) * e[0] + lam * e[1], max_relative = 1e-12);
            assert_relative_eq!(e_crf, (1.0 - lam) * e[2] + lam * e[3], max_relative = 1e-12);
        }
    }

    /// Any nonzero C6 enters the soft-core length through C12/C6, as in gromosXX: a C6 of
    /// 1e-12 is softened, not treated as zero.
    #[test]
    fn a_tiny_c6_is_softened() {
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let (c6, c12, alpha_lj, lam, r2) = (1e-12, 1e-6, 0.5, 0.5, 0.09);
        let lp = PerturbedLambdaParams::from_lambda(lam, lam, lam, lam, 1.0, 1.0, 1.0, 1.0, 1);
        let (_, e_lj, _, _, _) = perturbed_lj_crf_pair(
            r2, c6, c12, 0.0, 0.0, 0.0, 0.0, alpha_lj, 0.0, &lp, &crf, 1.0, 1.0,
        );
        let dist6 = r2 * r2 * r2 + alpha_lj * lam * lam * (c12 / c6);
        let softened = (1.0 - lam) * (c12 / dist6 - c6) / dist6;
        assert_relative_eq!(e_lj, softened, max_relative = 1e-12);
        let r6 = r2 * r2 * r2;
        let plain = (1.0 - lam) * (c12 / r6 - c6) / r6;
        assert!(e_lj < 1e-6 * plain, "{e_lj} vs unsoftened {plain}");
    }
}
