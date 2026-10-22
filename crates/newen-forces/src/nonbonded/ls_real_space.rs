//! The real-space loops of the lattice sum (NLRELE 2, 3): the pairlist pairs, the excluded pairs
//! and the 1-4 pairs, each through the kernels of `lattice_sum.rs`, as gromosXX's
//! `_ls_real_outerloop` and `one_four_interaction_innerloop` (`lj_ls_func` case) run them.

use super::lattice_sum::{lj_ls_pair, ls_excluded_pair, ChargeShape};
use super::{ForceStorage, LJParamMatrix};
use newen_core::math::{BoundaryCondition, Vec3};

/// The charge shape of a lattice-sum run and `1/a`.
#[derive(Debug, Clone, Copy)]
pub struct ShapeWidth {
    pub shape: ChargeShape,
    pub inv_width: f64,
}

/// Add one pair's force `r · f_mag` to atoms `i` and `j` and its virial `r ⊗ f`.
#[inline]
fn add_pair(storage: &mut ForceStorage, i: usize, j: usize, r: Vec3, f_mag: f64, virial: bool) {
    let force = Vec3::new(f_mag * r.x, f_mag * r.y, f_mag * r.z);
    storage.forces[i] += force;
    storage.forces[j] -= force;
    if virial {
        let rv = [r.x, r.y, r.z];
        let fv = [force.x, force.y, force.z];
        for a in 0..3 {
            for b in 0..3 {
                storage.virial[a][b] += rv[a] * fv[b];
            }
        }
    }
}

#[inline]
fn r2_of(r: Vec3) -> f64 {
    r.x * r.x + r.y * r.y + r.z * r.z
}

/// LJ + real-space lattice-sum term of every pair of a pairlist, gromosXX `lj_ls_real_innerloop`
/// (`nonbonded_innerloop.cc:2515-2600`). The LJ energy goes to `storage.e_lj`, the real-space
/// lattice-sum energy to `storage.e_crf` (the electrostatic slot of the storage).
#[allow(clippy::too_many_arguments)]
pub fn lj_ls_pairs<BC: BoundaryCondition>(
    pairs: &[(u32, u32)],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    sw: ShapeWidth,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    virial: bool,
) {
    for &(i, j) in pairs {
        let (i, j) = (i as usize, j as usize);
        let r = periodicity.nearest_image(positions[i], positions[j]);
        let lj = lj_params.get(iac[i] as usize, iac[j] as usize);
        let q_eps = charges[i] * charges[j] * four_pi_eps_i;
        let (f_mag, e_lj, e) = lj_ls_pair(r2_of(r), lj.c6, lj.c12, q_eps, sw.shape, sw.inv_width);
        add_pair(storage, i, j, r, f_mag, virial);
        storage.e_lj += e_lj;
        storage.e_crf += e;
    }
}

/// Lennard-Jones only, for the long-range shell of a twin-range lattice-sum run: gromosXX runs
/// `lj_crf_innerloop` there with `lj_ls_func`, which computes the LJ term alone
/// (`nonbonded_innerloop.cc:429-437`, `lj_interaction`). Through [`lj_ls_pair`] with no charge:
/// its LJ part is `lj_interaction` operation by operation, and a zero charge adds `+0.0` to the
/// force.
#[allow(clippy::too_many_arguments)]
pub fn lj_only_pairs<BC: BoundaryCondition>(
    pairs: &[(u32, u32)],
    positions: &[Vec3],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    sw: ShapeWidth,
    periodicity: &BC,
    storage: &mut ForceStorage,
    virial: bool,
) {
    for &(i, j) in pairs {
        let (i, j) = (i as usize, j as usize);
        let r = periodicity.nearest_image(positions[i], positions[j]);
        let lj = lj_params.get(iac[i] as usize, iac[j] as usize);
        let (f_mag, e_lj, _) = lj_ls_pair(r2_of(r), lj.c6, lj.c12, 0.0, sw.shape, sw.inv_width);
        add_pair(storage, i, j, r, f_mag, virial);
        storage.e_lj += e_lj;
    }
}

/// The real-space term of every excluded pair, gromosXX `ls_real_excluded_innerloop`
/// (`nonbonded_innerloop.cc:2603-2680`): the solute exclusions and the pairs inside a solvent
/// molecule, here both from `exclusions` (each pair once, `j > i`). The energy goes to
/// `storage.e_crf`.
#[allow(clippy::too_many_arguments)]
pub fn ls_excluded_pairs<BC: BoundaryCondition>(
    exclusions: &[Vec<usize>],
    positions: &[Vec3],
    charges: &[f64],
    sw: ShapeWidth,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    virial: bool,
) {
    for (i, excl) in exclusions.iter().enumerate() {
        for &j in excl {
            if j <= i {
                continue;
            }
            let r = periodicity.nearest_image(positions[i], positions[j]);
            let q_eps = charges[i] * charges[j] * four_pi_eps_i;
            let (f_mag, e) = ls_excluded_pair(r2_of(r), q_eps, sw.shape, sw.inv_width);
            add_pair(storage, i, j, r, f_mag, virial);
            storage.e_crf += e;
        }
    }
}

/// The 1-4 pairs of a lattice-sum run, gromosXX `one_four_interaction_innerloop`, case
/// `lj_ls_func` (`nonbonded_innerloop.cc:1309-1335`): [`lj_ls_pair`] with the 1-4 parameters
/// CS6, CS12 and the full charge product. The LJ energy goes to `storage.e_lj`, the real-space
/// lattice-sum energy to `storage.e_crf`.
#[allow(clippy::too_many_arguments)]
pub fn one_four_ls_pairs<BC: BoundaryCondition>(
    one_four_pairs: &[Vec<usize>],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    sw: ShapeWidth,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    virial: bool,
) {
    for (i, partners) in one_four_pairs.iter().enumerate() {
        for &j in partners {
            if j <= i {
                continue;
            }
            let r = periodicity.nearest_image(positions[i], positions[j]);
            let lj = lj_params.get(iac[i] as usize, iac[j] as usize);
            let q_eps = charges[i] * charges[j] * four_pi_eps_i;
            let (f_mag, e_lj, e) =
                lj_ls_pair(r2_of(r), lj.cs6, lj.cs12, q_eps, sw.shape, sw.inv_width);
            add_pair(storage, i, j, r, f_mag, virial);
            storage.e_lj += e_lj;
            storage.e_crf += e;
        }
    }
}
