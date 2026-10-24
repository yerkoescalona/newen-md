//! RF excluded-pair corrections and 1-4 interaction loop

use super::innerloops::with_groups;
use super::{CRFParameters, EnergyGroups, ForceStorage, LJParamMatrix};
use newen_core::math::{BoundaryCondition, Vec3};

/// Reaction-field correction of one excluded solute pair at squared distance `r2`:
/// `(f_mag, e)` with `f_mag = q · crf_cut3i` (`q C_RF/R³`) and
/// `e = q · (−crf_2cut3i · r2 − crf_cut)` (`q (−C_RF r²/(2R³) − (1 − C_RF/2)/R)`), the CRF pair
/// term without its Coulomb `1/r` part. The force on atom i is `r_vec · f_mag`.
///
/// spec: nonbonded.rf_excluded_pair
#[inline]
pub fn rf_excluded_pair(r2: f64, q_prod: f64, crf: &CRFParameters) -> (f64, f64) {
    let e = q_prod * (-crf.crf_2cut3i * r2 - crf.crf_cut);
    let f_mag = q_prod * crf.crf_cut3i;
    (f_mag, e)
}

/// Reaction-field self term of one solute atom with charge `q`:
/// `−0.5 · q · q · four_pi_eps_i · crf_cut`, half the excluded-pair energy of the atom with
/// itself at `r = 0`.
///
/// spec: nonbonded.rf_self_term
#[inline]
pub fn rf_self_term(q: f64, four_pi_eps_i: f64, crf: &CRFParameters) -> f64 {
    -0.5 * q * q * four_pi_eps_i * crf.crf_cut
}

/// Calculate reaction-field excluded-pair corrections (energy + forces).
///
/// In GROMOS, excluded atom pairs (bonded neighbours) are removed from the
/// pairlist.  With `NSLFEXCL=1` (default), they receive only the RF
/// **correction** — NOT the full Coulomb 1/r term.
///
/// **Solute atoms** (GROMOS `RF_excluded_interaction_innerloop`):
///   Self-term: `E_self = -0.5 * qi² * FPEPSI * crf_cut`
///   Excluded pairs: full `rf_interaction` with force + energy:
///     Force  = q_prod * FPEPSI * crf_cut3i * r_vec
///     Energy = q_prod * FPEPSI * (-crf_2cut3i * r² - crf_cut)
///
/// **Solvent atoms** (GROMOS `RF_solvent_interaction_innerloop`):
///   NO self-term (distance-independent parts cancel for neutral charge groups)
///   NO forces (rigid molecules)
///   Only energy: `E = -qi*qj * FPEPSI * crf_2cut3i * r²` (distance-dependent part only)
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn rf_excluded_interactions<BC: BoundaryCondition>(
    charges: &[f64],
    exclusions: &[Vec<usize>],
    positions: &[Vec3],
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    num_solute_atoms: usize,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => rf_excluded_impl::<BC, G>(
        charges,
        exclusions,
        positions,
        crf,
        boundary,
        four_pi_eps_i,
        storage,
        num_solute_atoms,
        &egroups,
    ));
}

/// [`rf_excluded_interactions`]; with `GROUPS` each term also goes to its group pair: the self
/// term to `[eg(i)][eg(i)]`, an excluded pair to `[eg(i)][eg(j)]` (gromosXX
/// `nonbonded_innerloop.cc:1706-1707, 1758-1759, 2129-2157`).
fn rf_excluded_impl<BC: BoundaryCondition, const GROUPS: bool>(
    charges: &[f64],
    exclusions: &[Vec<usize>],
    positions: &[Vec3],
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    num_solute_atoms: usize,
    egroups: &EnergyGroups<'_>,
) {
    // === Solute atoms ===
    let mut e_self = 0.0;
    let mut e_excl_solute = 0.0;
    let mut n_solute_excl = 0usize;

    for i in 0..num_solute_atoms {
        let term = rf_self_term(charges[i], four_pi_eps_i, crf);
        e_self += term;
        storage.e_crf += term;
        if GROUPS {
            storage.groups.crf[egroups.cell(i, i)] += term;
        }

        // Excluded pair RF correction (energy + forces)
        for &j in exclusions[i].iter() {
            if j > i {
                let r = boundary.nearest_image(positions[i], positions[j]);
                let r2 = r.length_squared();
                let q_prod = charges[i] * charges[j] * four_pi_eps_i;

                let (f_mag, term) = rf_excluded_pair(r2, q_prod, crf);
                e_excl_solute += term;
                storage.e_crf += term;
                if GROUPS {
                    storage.groups.crf[egroups.cell(i, j)] += term;
                }
                n_solute_excl += 1;

                let force = r * f_mag;
                storage.forces[i] += force;
                storage.forces[j] -= force;

                // Virial: GROMOS convention virial_tensor(b, a) += r(b) * force(a)
                let rv = [r.x, r.y, r.z];
                let fv = [force.x, force.y, force.z];
                for a in 0..3 {
                    for b in 0..3 {
                        storage.virial[a][b] += rv[a] * fv[b];
                    }
                }
            }
        }
    }

    // === Solvent atoms ===
    // GROMOS: no self-term, no forces, only distance-dependent energy
    let mut e_excl_solvent = 0.0;
    let mut n_solvent_excl = 0usize;

    for i in num_solute_atoms..charges.len() {
        for &j in exclusions[i].iter() {
            if j > i {
                let r = boundary.nearest_image(positions[i], positions[j]);
                let r2 = r.length_squared();

                // Only: -qi*qj * FPEPSI * crf_2cut3i * r²
                let term = -charges[i] * charges[j] * four_pi_eps_i * crf.crf_2cut3i * r2;
                e_excl_solvent += term;
                storage.e_crf += term;
                if GROUPS {
                    storage.groups.crf[egroups.cell(i, j)] += term;
                }
                n_solvent_excl += 1;

                // NO force for solvent excluded pairs (rigid molecules)
            }
        }
    }
    log::debug!(
        "  RF solute: self={:.10e}, excl={:.10e} ({} pairs)",
        e_self,
        e_excl_solute,
        n_solute_excl
    );
    log::debug!(
        "  RF solvent: excl={:.10e} ({} pairs)",
        e_excl_solvent,
        n_solvent_excl
    );
    log::debug!(
        "  RF total_corr={:.10e}",
        e_self + e_excl_solute + e_excl_solvent
    );
}

/// Calculate 1-4 nonbonded interactions (LJ with cs6/cs12 + CRF with coulomb scaling).
///
/// In GROMOS, 1-4 pairs are excluded from the pairlist and computed separately
/// using scaled LJ parameters (cs6, cs12) and optionally scaled Coulomb.
/// No cutoff is applied — all 1-4 pairs are always computed.
///
/// # Arguments
/// * `one_four_pairs` - Per-atom lists of 1-4 partner indices (only j > i processed)
/// * `positions` - Current atomic positions
/// * `charges` - Per-atom partial charges
/// * `iac` - Integer atom codes (atom types)
/// * `lj_params` - LJ parameter matrix (indexed by atom type)
/// * `crf` - CRF parameters
/// * `boundary` - Boundary condition for nearest-image calculation
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
/// * `storage` - Force/energy accumulator
/// * `coulomb_scaling` - Scaling factor for Coulomb 1-4 interactions (1.0 for GROMOS, 0.5 for AMBER)
pub fn one_four_interaction_loop<BC: BoundaryCondition>(
    one_four_pairs: &[Vec<usize>],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    coulomb_scaling: f64,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => one_four_impl::<BC, G>(
        one_four_pairs,
        positions,
        charges,
        iac,
        lj_params,
        crf,
        boundary,
        four_pi_eps_i,
        storage,
        coulomb_scaling,
        &egroups,
    ));
}

/// [`one_four_interaction_loop`]; with `GROUPS` each pair also goes to `[eg(i)][eg(j)]`
/// (gromosXX `nonbonded_innerloop.cc:1383-1390`).
fn one_four_impl<BC: BoundaryCondition, const GROUPS: bool>(
    one_four_pairs: &[Vec<usize>],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    boundary: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    coulomb_scaling: f64,
    egroups: &EnergyGroups<'_>,
) {
    let mut e_lj_14 = 0.0;
    let mut e_crf_14 = 0.0;
    let mut n_pairs = 0usize;

    for i in 0..one_four_pairs.len() {
        for &j in &one_four_pairs[i] {
            if j <= i {
                continue;
            }

            let r = boundary.nearest_image(positions[i], positions[j]);
            let r2 = r.length_squared();

            if r2 < 1e-10 {
                continue;
            }

            let type_i = iac[i] as usize;
            let type_j = iac[j] as usize;
            let lj = lj_params.get(type_i, type_j);

            let inv_r2 = 1.0 / r2;
            let inv_r6 = inv_r2 * inv_r2 * inv_r2;

            // LJ with 1-4 parameters (cs6, cs12)
            let e_lj = (lj.cs12 * inv_r6 - lj.cs6) * inv_r6;
            let f_lj = (12.0 * lj.cs12 * inv_r6 - 6.0 * lj.cs6) * inv_r6 * inv_r2;

            // CRF with coulomb scaling
            let q_prod = charges[i] * charges[j] * four_pi_eps_i;
            let inv_r = inv_r2.sqrt();
            let e_crf = q_prod * (inv_r * coulomb_scaling - crf.crf_2cut3i * r2 - crf.crf_cut);
            let f_crf = q_prod * (inv_r * coulomb_scaling * inv_r2 + crf.crf_cut3i);

            let f_total = f_lj + f_crf;
            let force = r * f_total;

            storage.forces[i] += force;
            storage.forces[j] -= force;

            e_lj_14 += e_lj;
            e_crf_14 += e_crf;
            if GROUPS {
                let k = egroups.cell(i, j);
                storage.groups.lj[k] += e_lj;
                storage.groups.crf[k] += e_crf;
            }

            let rv = [r.x, r.y, r.z];
            let fv = [force.x, force.y, force.z];
            for a in 0..3 {
                for b in 0..3 {
                    storage.virial[a][b] += rv[a] * fv[b];
                }
            }

            n_pairs += 1;
        }
    }

    storage.e_lj += e_lj_14;
    storage.e_crf += e_crf_14;
    log::debug!(
        "  1-4 interactions: e_lj={:.10e}, e_crf={:.10e} ({} pairs)",
        e_lj_14,
        e_crf_14,
        n_pairs
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::Vacuum;

    /// Two excluded solute atoms: the self terms plus the excluded-pair reaction-field term,
    /// forces equal and opposite along the pair vector (gromosXX `RF_excluded_interaction`).
    #[test]
    fn excluded_pair_energy_and_forces_follow_the_reaction_field_formula() {
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let charges = [0.5, -0.3];
        let exclusions = vec![vec![1], vec![]];
        let r = Vec3::new(0.1, 0.2, 0.0);
        let positions = [Vec3::ZERO, r];
        let fpepsi = 138.9354;
        let mut storage = ForceStorage::new(2);
        rf_excluded_interactions(
            &charges,
            &exclusions,
            &positions,
            &crf,
            &Vacuum,
            fpepsi,
            &mut storage,
            2,
            EnergyGroups::single(),
        );
        let q_prod = charges[0] * charges[1] * fpepsi;
        let expected = -0.5 * (charges[0].powi(2) + charges[1].powi(2)) * fpepsi * crf.crf_cut
            + q_prod * (-crf.crf_2cut3i * r.length_squared() - crf.crf_cut);
        assert!(
            (storage.e_crf - expected).abs() < 1e-12,
            "{} vs {expected}",
            storage.e_crf
        );
        let f = (positions[0] - positions[1]) * (q_prod * crf.crf_cut3i);
        assert!((storage.forces[0] - f).length() < 1e-12);
        assert!((storage.forces[1] + f).length() < 1e-12);
    }

    #[test]
    fn only_solute_atoms_get_self_terms() {
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let charges = [1.0, 1.0];
        let exclusions = vec![vec![], vec![]];
        let positions = [Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)];
        let mut storage = ForceStorage::new(2);
        rf_excluded_interactions(
            &charges,
            &exclusions,
            &positions,
            &crf,
            &Vacuum,
            1.0,
            &mut storage,
            1,
            EnergyGroups::single(),
        );
        assert!((storage.e_crf - (-0.5 * crf.crf_cut)).abs() < 1e-15);
    }
}
