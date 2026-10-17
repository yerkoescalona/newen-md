//! LJ+CRF inner loop kernels: serial, parallel, CG-grouped, and solvent variants

use super::{CGPairGroup, CRFParameters, EnergyGroups, ForceStorage, GroupEnergies, LJParamMatrix};
use gromos_core::math::{BoundaryCondition, Vec3};

use rayon::prelude::*;
use wide::{f64x4, CmpLt};

/// Core LJ + CRF interaction calculation (hot path!)
///
/// # Arguments
/// * `r` - Distance vector from i to j
/// * `c6` - LJ C6 coefficient
/// * `c12` - LJ C12 coefficient
/// * `q_prod` - Charge product qi * qj * four_pi_eps_i
/// * `crf` - CRF parameters
///
/// # Returns
/// * `force_magnitude` - Scalar force magnitude
/// * `e_lj` - Lennard-Jones energy
/// * `e_crf` - Coulomb reaction field energy
#[inline(always)]
pub fn lj_crf_interaction(
    r: Vec3,
    c6: f64,
    c12: f64,
    q_prod: f64,
    crf: &CRFParameters,
) -> (f64, f64, f64) {
    let r2 = r.length_squared();

    // Early exit for zero distance (should not happen, but safety)
    if r2 < 1e-10 {
        return (0.0, 0.0, 0.0);
    }

    let inv_r2 = 1.0 / r2;
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;

    // Lennard-Jones: E_lj = C12/r^12 - C6/r^6
    let e_lj = (c12 * inv_r6 - c6) * inv_r6;

    // Lennard-Jones force: F = 12*C12/r^14 - 6*C6/r^8
    let f_lj = (12.0 * c12 * inv_r6 - 6.0 * c6) * inv_r6 * inv_r2;

    // Coulomb Reaction Field
    // E_crf = q * (1/r - crf_2cut3i * r² - crf_cut)
    // F_crf = q * (1/r³ + crf_cut3i)
    let inv_r = inv_r2.sqrt();
    let e_crf = q_prod * (inv_r - crf.crf_2cut3i * r2 - crf.crf_cut);
    let f_crf = q_prod * (inv_r * inv_r2 + crf.crf_cut3i);

    let force_magnitude = f_lj + f_crf;

    (force_magnitude, e_lj, e_crf)
}

/// Pairs evaluated per SIMD batch — one AVX2 `f64x4` register.
const BATCH: usize = 4;

/// [`lj_crf_interaction`] for four pairs at once.
///
/// Takes `r²` rather than `r` so the caller decides how the separation is
/// formed (minimum image vs. charge-group shift). Every expression is written
/// in the same association order as the scalar kernel, and `wide` lowers `/`
/// and `sqrt` to the correctly-rounded `vdivpd`/`vsqrtpd`, so each lane is
/// bit-identical to the scalar result for the same pair. That is what lets the
/// batched and remainder paths coexist without changing any reference output.
///
/// The divide and the square root are the expensive part of the pair
/// arithmetic; doing four of each per instruction is the point.
#[inline(always)]
fn lj_crf_interaction_x4(
    r2: [f64; BATCH],
    c6: [f64; BATCH],
    c12: [f64; BATCH],
    q_prod: [f64; BATCH],
    crf: &CRFParameters,
) -> ([f64; BATCH], [f64; BATCH], [f64; BATCH]) {
    let r2 = f64x4::from(r2);
    let c6 = f64x4::from(c6);
    let c12 = f64x4::from(c12);
    let q_prod = f64x4::from(q_prod);

    let inv_r2 = f64x4::splat(1.0) / r2;
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;

    let e_lj = (c12 * inv_r6 - c6) * inv_r6;
    let f_lj = (f64x4::splat(12.0) * c12 * inv_r6 - f64x4::splat(6.0) * c6) * inv_r6 * inv_r2;

    let inv_r = inv_r2.sqrt();
    let e_crf = q_prod * (inv_r - f64x4::splat(crf.crf_2cut3i) * r2 - f64x4::splat(crf.crf_cut));
    let f_crf = q_prod * (inv_r * inv_r2 + f64x4::splat(crf.crf_cut3i));

    let force = f_lj + f_crf;

    // Same zero-distance guard as the scalar kernel, as a lane mask.
    let near_zero = r2.cmp_lt(f64x4::splat(1e-10));
    let zero = f64x4::splat(0.0);
    (
        near_zero.blend(zero, force).to_array(),
        near_zero.blend(zero, e_lj).to_array(),
        near_zero.blend(zero, e_crf).to_array(),
    )
}

/// Run `$body` with the const `$g` set to whether `$groups` books per energy group, so each
/// kernel picks its instantiation once, outside the pair loop.
macro_rules! with_groups {
    ($groups:expr, $g:ident => $body:expr) => {
        if $groups.is_grouped() {
            const $g: bool = true;
            $body
        } else {
            const $g: bool = false;
            $body
        }
    };
}
pub(crate) use with_groups;

/// Energies of consecutive pairs in the same energy-group cell, summed in registers and added
/// to the cell when the pairs move on to another one (or the kernel ends), so the group
/// booking is not a read-modify-write of memory per pair.
struct GroupRun {
    cell: usize,
    lj: f64,
    crf: f64,
}

impl GroupRun {
    const NONE: usize = usize::MAX;

    #[inline(always)]
    fn new() -> Self {
        Self {
            cell: Self::NONE,
            lj: 0.0,
            crf: 0.0,
        }
    }

    #[inline(always)]
    fn add(&mut self, cell: usize, e_lj: f64, e_crf: f64, groups: &mut GroupEnergies) {
        if cell == self.cell {
            self.lj += e_lj;
            self.crf += e_crf;
        } else {
            self.flush(groups);
            *self = Self {
                cell,
                lj: e_lj,
                crf: e_crf,
            };
        }
    }

    #[inline(always)]
    fn flush(&mut self, groups: &mut GroupEnergies) {
        if self.cell != Self::NONE {
            groups.lj[self.cell] += self.lj;
            groups.crf[self.cell] += self.crf;
            self.cell = Self::NONE;
        }
    }
}

/// Scatter one evaluated pair into the accumulator. Shared by the batched and
/// remainder paths so the accumulation order is identical to the scalar loop.
/// With `GROUPS` the pair's energies also go to cell `[eg(i)][eg(j)]` (gromosXX
/// `nonbonded_innerloop.cc:521-524`), through `run`.
#[inline(always)]
fn accumulate_pair<const VIRIAL: bool, const GROUPS: bool>(
    storage: &mut ForceStorage,
    egroups: &EnergyGroups<'_>,
    run: &mut GroupRun,
    i: usize,
    j: usize,
    r: Vec3,
    f_mag: f64,
    e_lj: f64,
    e_crf: f64,
) {
    let force = r * f_mag;
    storage.forces[i] += force;
    storage.forces[j] -= force;

    storage.e_lj += e_lj;
    storage.e_crf += e_crf;
    if GROUPS {
        run.add(egroups.cell(i, j), e_lj, e_crf, &mut storage.groups);
    }

    if VIRIAL {
        for a in 0..3 {
            for b in 0..3 {
                storage.virial[a][b] += r[a] * force[b];
            }
        }
    }
}

/// The nonbonded innerloop over a slice of atom pairs.
///
/// This is the single source of truth for the LJ + CRF pair physics. The flat,
/// charge-group-grouped, and solvent kernels differ only in how the separation
/// vector is formed and where the pair's `(c6, c12, q_prod)` come from, which
/// they pass in as `separation` and `params`. Pairs are evaluated four at a
/// time through [`lj_crf_interaction_x4`], with the tail through the scalar
/// [`lj_crf_interaction`]; both are scattered by [`accumulate_pair`] in pair
/// order, so results match the original one-pair-at-a-time loops exactly.
///
/// Only the divide/sqrt arithmetic is vectorised; the separations stay as
/// per-pair `Vec3`s. A structure-of-arrays variant that also vectorised the
/// minimum image and the force products was measured slower (the AoS↔SoA
/// shuffles cost more than the dozen scalar ops they replaced), so this is
/// the deliberate shape, not an oversight.
#[inline(always)]
fn process_pair_slice<const VIRIAL: bool, const GROUPS: bool>(
    pairs: &[(u32, u32)],
    separation: impl Fn(usize, usize) -> Vec3,
    params: impl Fn(usize, usize) -> (f64, f64, f64),
    crf: &CRFParameters,
    egroups: &EnergyGroups<'_>,
    run: &mut GroupRun,
    storage: &mut ForceStorage,
) {
    let mut batches = pairs.chunks_exact(BATCH);

    for batch in &mut batches {
        let mut r = [Vec3::ZERO; BATCH];
        let mut r2 = [0.0; BATCH];
        let mut c6 = [0.0; BATCH];
        let mut c12 = [0.0; BATCH];
        let mut q_prod = [0.0; BATCH];

        for (lane, &(i, j)) in batch.iter().enumerate() {
            let i = i as usize;
            let j = j as usize;
            r[lane] = separation(i, j);
            r2[lane] = r[lane].length_squared();
            let (pair_c6, pair_c12, pair_q) = params(i, j);
            c6[lane] = pair_c6;
            c12[lane] = pair_c12;
            q_prod[lane] = pair_q;
        }

        let (f_mag, e_lj, e_crf) = lj_crf_interaction_x4(r2, c6, c12, q_prod, crf);

        for (lane, &(i, j)) in batch.iter().enumerate() {
            accumulate_pair::<VIRIAL, GROUPS>(
                storage,
                egroups,
                run,
                i as usize,
                j as usize,
                r[lane],
                f_mag[lane],
                e_lj[lane],
                e_crf[lane],
            );
        }
    }

    for &(i, j) in batches.remainder() {
        let i = i as usize;
        let j = j as usize;
        let r = separation(i, j);
        let (c6, c12, q_prod) = params(i, j);
        let (f_mag, e_lj, e_crf) = lj_crf_interaction(r, c6, c12, q_prod, crf);
        accumulate_pair::<VIRIAL, GROUPS>(storage, egroups, run, i, j, r, f_mag, e_lj, e_crf);
    }
}

/// CG-grouped kernel: compute nearest_image once per charge-group pair,
/// then process all atom pairs in that group using simple subtraction + shift.
///
/// This eliminates per-atom-pair nearest_image calls (expensive floor/branch ops)
/// and replaces them with cheap additions. For 3-atom CGs (water), this is ~3x fewer
/// nearest_image calls.
#[inline]
fn process_pairs_cg_grouped<BC: BoundaryCondition, const VIRIAL: bool, const GROUPS: bool>(
    pairs: &[(u32, u32)],
    groups: &[CGPairGroup],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    let mut run = GroupRun::new();
    for group in groups {
        let ri = positions[group.ref_atom_i as usize];
        let rj = positions[group.ref_atom_j as usize];

        let r_ref = periodicity.nearest_image(ri, rj);
        let tx = r_ref.x - ri.x + rj.x;
        let ty = r_ref.y - ri.y + rj.y;
        let tz = r_ref.z - ri.z + rj.z;

        let block = &pairs[group.start as usize..group.end as usize];
        process_pair_slice::<VIRIAL, GROUPS>(
            block,
            |i, j| {
                Vec3::new(
                    positions[i].x + tx - positions[j].x,
                    positions[i].y + ty - positions[j].y,
                    positions[i].z + tz - positions[j].z,
                )
            },
            |i, j| {
                let lj = lj_params.get(iac[i] as usize, iac[j] as usize);
                (lj.c6, lj.c12, charges[i] * charges[j] * four_pi_eps_i)
            },
            crf,
            &egroups,
            &mut run,
            storage,
        );
    }
    run.flush(&mut storage.groups);
}

/// Shared kernel: process a slice of atom pairs into a ForceStorage.
///
/// This is the single source of truth for the nonbonded innerloop physics.
/// Both serial and parallel paths call this — no code duplication.
#[inline]
fn process_pairs<BC: BoundaryCondition, const VIRIAL: bool, const GROUPS: bool>(
    pairs: &[(u32, u32)],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    let mut run = GroupRun::new();
    process_pair_slice::<VIRIAL, GROUPS>(
        pairs,
        |i, j| periodicity.nearest_image(positions[i], positions[j]),
        |i, j| {
            let lj = lj_params.get(iac[i] as usize, iac[j] as usize);
            (lj.c6, lj.c12, charges[i] * charges[j] * four_pi_eps_i)
        },
        crf,
        &egroups,
        &mut run,
        storage,
    );
    run.flush(&mut storage.groups);
}

/// Shared kernel: process a slice of solvent molecule pairs into a ForceStorage.
///
/// Single source of truth for solvent-solvent interactions with shared PBC shift.
/// Precomputes LJ params and charge products for the N×N atom type combinations
/// within a solvent molecule (GROMOS optimization: avoids repeated lookups).
#[inline]
fn process_solvent_pairs<BC: BoundaryCondition, const VIRIAL: bool, const GROUPS: bool>(
    pairs: &[(u32, u32)],
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    if pairs.is_empty() {
        return;
    }
    let mut run = GroupRun::new();

    let first_mol = pairs[0].0 as usize;
    let n = atoms_per_solvent;
    // Stack-allocated arrays for up to 4-site water models (covers SPC, SPC/E, TIP4P)
    assert!(
        n <= 4,
        "solvent with >4 atoms not supported in optimized path"
    );
    let mut lj_c6 = [[0.0f64; 4]; 4];
    let mut lj_c12 = [[0.0f64; 4]; 4];
    let mut q_prod = [[0.0f64; 4]; 4];

    for ai in 0..n {
        let type_i = iac[first_mol + ai] as usize;
        let qi = charges[first_mol + ai];
        for aj in 0..n {
            let type_j = iac[first_mol + aj] as usize;
            let lj = lj_params.get(type_i, type_j);
            lj_c6[ai][aj] = lj.c6;
            lj_c12[ai][aj] = lj.c12;
            q_prod[ai][aj] = qi * charges[first_mol + aj] * four_pi_eps_i;
        }
    }

    // One sentinel pair expands to the n×n atom pairs of the two molecules,
    // enumerated in (atom_i, atom_j) order so accumulation matches the
    // original nested loop. n ≤ 4 is asserted above, hence the fixed 16.
    let mut block = [(0u32, 0u32); 16];
    for &(i_first, j_first) in pairs {
        let i_first = i_first as usize;
        let j_first = j_first as usize;

        let pos_i0 = positions[i_first];
        let pos_j0 = positions[j_first];

        let r_first = periodicity.nearest_image(pos_i0, pos_j0);
        let tx = r_first.x - pos_i0.x + pos_j0.x;
        let ty = r_first.y - pos_i0.y + pos_j0.y;
        let tz = r_first.z - pos_i0.z + pos_j0.z;

        let mut count = 0;
        for atom_i in 0..n {
            for atom_j in 0..n {
                block[count] = ((i_first + atom_i) as u32, (j_first + atom_j) as u32);
                count += 1;
            }
        }

        process_pair_slice::<VIRIAL, GROUPS>(
            &block[..count],
            |i, j| {
                Vec3::new(
                    positions[i].x + tx - positions[j].x,
                    positions[i].y + ty - positions[j].y,
                    positions[i].z + tz - positions[j].z,
                )
            },
            |i, j| {
                let (ai, aj) = (i - i_first, j - j_first);
                (lj_c6[ai][aj], lj_c12[ai][aj], q_prod[ai][aj])
            },
            crf,
            &egroups,
            &mut run,
            storage,
        );
    }
    run.flush(&mut storage.groups);
}

/// Minimum pairlist size to justify parallel overhead.
const PARALLEL_THRESHOLD: usize = 2048;

/// Inner loop for nonbonded interactions (solute-solute and solute-solvent).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => process_pairs::<BC, true, G>(
        pairlist,
        positions,
        charges,
        iac,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        storage,
        egroups,
    ));
}

/// Inner loop without virial computation (for NVE/NVT without pressure coupling).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => process_pairs::<BC, false, G>(
        pairlist,
        positions,
        charges,
        iac,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        storage,
        egroups,
    ));
}

/// CG-grouped innerloop with virial (for NPT).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop_cg_grouped<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => process_pairs_cg_grouped::<BC, true, G>(
        pairlist,
        groups,
        positions,
        charges,
        iac,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        storage,
        egroups,
    ));
}

/// CG-grouped innerloop without virial (for NVE/NVT).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop_cg_grouped_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => process_pairs_cg_grouped::<BC, false, G>(
        pairlist,
        groups,
        positions,
        charges,
        iac,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        storage,
        egroups,
    ));
}

/// Parallel CG-grouped innerloop with virial.
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop_cg_grouped_parallel<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    with_groups!(egroups, G => lj_crf_innerloop_cg_grouped_parallel_virial::<BC, true, G>(
        positions,
        charges,
        iac,
        pairlist,
        groups,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        n_atoms,
        egroups,
    ))
}

/// Parallel CG-grouped innerloop without virial.
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop_cg_grouped_parallel_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    with_groups!(egroups, G => lj_crf_innerloop_cg_grouped_parallel_virial::<BC, false, G>(
        positions,
        charges,
        iac,
        pairlist,
        groups,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        n_atoms,
        egroups,
    ))
}

fn lj_crf_innerloop_cg_grouped_parallel_virial<
    BC: BoundaryCondition,
    const VIRIAL: bool,
    const GROUPS: bool,
>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    groups: &[CGPairGroup],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    let mut result = ForceStorage::with_groups(n_atoms, egroups.n_groups());
    if groups.len() < 64 {
        process_pairs_cg_grouped::<BC, VIRIAL, GROUPS>(
            pairlist,
            groups,
            positions,
            charges,
            iac,
            lj_params,
            crf,
            periodicity,
            four_pi_eps_i,
            &mut result,
            egroups,
        );
        return result;
    }

    let pc_fpi = four_pi_eps_i;
    result = groups
        .par_chunks(
            groups
                .len()
                .div_ceil(rayon::current_num_threads().max(1))
                .max(1),
        )
        .fold(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut local_storage, chunk| {
                process_pairs_cg_grouped::<BC, VIRIAL, GROUPS>(
                    pairlist,
                    chunk,
                    positions,
                    charges,
                    iac,
                    lj_params,
                    crf,
                    periodicity,
                    pc_fpi,
                    &mut local_storage,
                    egroups,
                );
                local_storage
            },
        )
        .reduce(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut a, b| {
                a.merge(&b);
                a
            },
        );

    result
}

/// Parallel nonbonded innerloop — same physics as `lj_crf_innerloop`.
///
/// Uses Rayon fold/reduce with thread-local ForceStorage buffers (GROMOS pattern:
/// per-thread Nonbonded_Set with private force arrays, sequential reduction).
///
/// Automatically falls back to serial for small pairlists.
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop_parallel<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    with_groups!(egroups, G => lj_crf_innerloop_parallel_virial::<BC, true, G>(
        positions,
        charges,
        iac,
        pairlist,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        n_atoms,
        egroups,
    ))
}

/// Parallel nonbonded innerloop without virial (for NVE/NVT).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn lj_crf_innerloop_parallel_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    with_groups!(egroups, G => lj_crf_innerloop_parallel_virial::<BC, false, G>(
        positions,
        charges,
        iac,
        pairlist,
        lj_params,
        crf,
        periodicity,
        four_pi_eps_i,
        n_atoms,
        egroups,
    ))
}

fn lj_crf_innerloop_parallel_virial<
    BC: BoundaryCondition,
    const VIRIAL: bool,
    const GROUPS: bool,
>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    let mut result = ForceStorage::with_groups(n_atoms, egroups.n_groups());
    if pairlist.len() < PARALLEL_THRESHOLD {
        process_pairs::<BC, VIRIAL, GROUPS>(
            pairlist,
            positions,
            charges,
            iac,
            lj_params,
            crf,
            periodicity,
            four_pi_eps_i,
            &mut result,
            egroups,
        );
        return result;
    }

    // One chunk per thread, not a fixed 1024. Each fold accumulator is a full
    // n_atoms force buffer that has to be allocated, zeroed, and merged back, so
    // chunking finer than the thread count multiplies that cost for no added
    // parallelism — at 62k pairs the fixed size produced 61 buffers per step and
    // dominated the kernel. Matches the CG-grouped parallel path's sizing.
    let pc_fpi = four_pi_eps_i;
    let chunk = pairlist
        .len()
        .div_ceil(rayon::current_num_threads().max(1))
        .max(1);
    result = pairlist
        .par_chunks(chunk)
        .fold(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut acc, chunk| {
                process_pairs::<BC, VIRIAL, GROUPS>(
                    chunk,
                    positions,
                    charges,
                    iac,
                    lj_params,
                    crf,
                    periodicity,
                    pc_fpi,
                    &mut acc,
                    egroups,
                );
                acc
            },
        )
        .reduce(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut a, b| {
                a.merge(&b);
                a
            },
        );
    result
}

/// Solvent-solvent innerloop with shared PBC shift (GROMOS convention).
///
/// For each solvent-solvent CG pair, the PBC shift is computed once from
/// the first atoms (typically O-O for water) and reused for all atom pairs.
///
/// The pairlist stores pairs of first-atom indices (one per solvent molecule).
/// `atoms_per_solvent` defines the molecule size (e.g. 3 for water).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn solvent_innerloop<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => process_solvent_pairs::<BC, true, G>(
        pairlist,
        positions,
        charges,
        iac,
        lj_params,
        crf,
        periodicity,
        atoms_per_solvent,
        four_pi_eps_i,
        storage,
        egroups,
    ));
}

/// Solvent-solvent innerloop without virial computation (for NVE/NVT).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn solvent_innerloop_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    four_pi_eps_i: f64,
    storage: &mut ForceStorage,
    egroups: EnergyGroups<'_>,
) {
    storage.check_groups(&egroups);
    with_groups!(egroups, G => process_solvent_pairs::<BC, false, G>(
        pairlist,
        positions,
        charges,
        iac,
        lj_params,
        crf,
        periodicity,
        atoms_per_solvent,
        four_pi_eps_i,
        storage,
        egroups,
    ));
}

/// Parallel solvent-solvent innerloop — same physics as `solvent_innerloop`.
///
/// Uses Rayon fold/reduce with thread-local buffers. Falls back to serial
/// for small pairlists.
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn solvent_innerloop_parallel<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    let mut result = ForceStorage::with_groups(n_atoms, egroups.n_groups());
    if pairlist.len() < PARALLEL_THRESHOLD {
        with_groups!(egroups, G => process_solvent_pairs::<BC, true, G>(
            pairlist,
            positions,
            charges,
            iac,
            lj_params,
            crf,
            periodicity,
            atoms_per_solvent,
            four_pi_eps_i,
            &mut result,
            egroups,
        ));
        return result;
    }

    // One chunk per thread; see the note in lj_crf_innerloop_parallel_virial.
    let pc_fpi = four_pi_eps_i;
    let chunk = pairlist
        .len()
        .div_ceil(rayon::current_num_threads().max(1))
        .max(1);
    result = pairlist
        .par_chunks(chunk)
        .fold(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut acc, chunk| {
                with_groups!(egroups, G => process_solvent_pairs::<BC, true, G>(
                    chunk,
                    positions,
                    charges,
                    iac,
                    lj_params,
                    crf,
                    periodicity,
                    atoms_per_solvent,
                    pc_fpi,
                    &mut acc,
                    egroups,
                ));
                acc
            },
        )
        .reduce(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut a, b| {
                a.merge(&b);
                a
            },
        );
    result
}

/// Parallel solvent-solvent innerloop without virial (for NVE/NVT).
///
/// # Arguments
/// * `constants` - Physical constants (controls `four_pi_eps_i` value)
pub fn solvent_innerloop_parallel_novirial<BC: BoundaryCondition>(
    positions: &[Vec3],
    charges: &[f64],
    iac: &[u32],
    pairlist: &[(u32, u32)],
    lj_params: &LJParamMatrix,
    crf: &CRFParameters,
    periodicity: &BC,
    atoms_per_solvent: usize,
    four_pi_eps_i: f64,
    n_atoms: usize,
    egroups: EnergyGroups<'_>,
) -> ForceStorage {
    let mut result = ForceStorage::with_groups(n_atoms, egroups.n_groups());
    if pairlist.len() < PARALLEL_THRESHOLD {
        with_groups!(egroups, G => process_solvent_pairs::<BC, false, G>(
            pairlist,
            positions,
            charges,
            iac,
            lj_params,
            crf,
            periodicity,
            atoms_per_solvent,
            four_pi_eps_i,
            &mut result,
            egroups,
        ));
        return result;
    }

    // One chunk per thread; see the note in lj_crf_innerloop_parallel_virial.
    let pc_fpi = four_pi_eps_i;
    let chunk = pairlist
        .len()
        .div_ceil(rayon::current_num_threads().max(1))
        .max(1);
    result = pairlist
        .par_chunks(chunk)
        .fold(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut acc, chunk| {
                with_groups!(egroups, G => process_solvent_pairs::<BC, false, G>(
                    chunk,
                    positions,
                    charges,
                    iac,
                    lj_params,
                    crf,
                    periodicity,
                    atoms_per_solvent,
                    pc_fpi,
                    &mut acc,
                    egroups,
                ));
                acc
            },
        )
        .reduce(
            || ForceStorage::with_groups(n_atoms, egroups.n_groups()),
            |mut a, b| {
                a.merge(&b);
                a
            },
        );
    result
}

#[cfg(test)]
mod tests {
    use super::super::params::LJParameters;
    use super::*;
    use approx::assert_relative_eq;
    use gromos_core::math::{Rectangular, Vacuum};

    #[test]
    fn test_lj_interaction() {
        let r = Vec3::new(1.0, 0.0, 0.0);
        let c6 = 0.001;
        let c12 = 0.0001;
        let q_prod = 0.0;

        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.0,
            crf_cut3i: 0.0,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let (f, e_lj, _e_crf) = lj_crf_interaction(r, c6, c12, q_prod, &crf);

        let expected_e_lj = c12 - c6;
        assert_relative_eq!(e_lj, expected_e_lj, epsilon = 1e-9);
        assert!(f < 0.0, "Force should be attractive");
    }

    #[test]
    fn test_innerloop_simple() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let charges = vec![0.5, -0.5];
        let iac = vec![0, 0];
        let pairlist = vec![(0, 1)];

        let lj_params = LJParamMatrix::from_nested(&[vec![LJParameters {
            c6: 0.001,
            c12: 0.0001,
            cs6: 0.001,
            cs12: 0.0001,
        }]]);
        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.364431,
            crf_cut3i: 0.364431 / 2.0,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let periodicity = Vacuum;
        let mut storage = ForceStorage::new(2);

        lj_crf_innerloop(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            gromos_core::units::four_pi_eps_i,
            &mut storage,
            EnergyGroups::single(),
        );

        assert_relative_eq!(storage.forces[0].x, -storage.forces[1].x, epsilon = 1e-6);
        assert_relative_eq!(storage.forces[0].y, -storage.forces[1].y, epsilon = 1e-6);
        assert_relative_eq!(storage.forces[0].z, -storage.forces[1].z, epsilon = 1e-6);
        assert!(storage.e_lj != 0.0);
        assert!(storage.e_crf != 0.0);
    }

    #[test]
    fn test_periodic_boundary() {
        let box_size = Vec3::splat(10.0);
        let periodicity = Rectangular::new(box_size);

        let positions = vec![Vec3::new(9.5, 0.0, 0.0), Vec3::new(0.5, 0.0, 0.0)];
        let charges = vec![0.0, 0.0];
        let iac = vec![0, 0];
        let pairlist = vec![(0, 1)];

        let lj_params = LJParamMatrix::from_nested(&[vec![LJParameters {
            c6: 0.001,
            c12: 0.0001,
            cs6: 0.001,
            cs12: 0.0001,
        }]]);
        let crf = CRFParameters {
            crf_cut: 1.4,
            crf_2cut3i: 0.0,
            crf_cut3i: 0.0,
            cutoff_sq: 1.4_f64.powi(2),
        };

        let mut storage = ForceStorage::new(2);

        lj_crf_innerloop(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            gromos_core::units::four_pi_eps_i,
            &mut storage,
            EnergyGroups::single(),
        );

        assert!(
            storage.forces[0].length() < 10.0,
            "Force too large - PBC not working"
        );
    }
}
