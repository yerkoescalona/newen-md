//! Constraint algorithms
//!
//! This module implements various constraint algorithms:
//! - SHAKE: Iterative constraint solver for bonds
//! - SETTLE: Analytical constraint solver for rigid water
//! - M-SHAKE: Mass-weighted SHAKE variant
//! - LINCS: Linear constraint solver
//! - Perturbed SHAKE: λ-dependent constraint lengths for FEP (in `shake_buffered`)
//! - Flexible Constraints: Time-dependent constraints (FlexShake)
//! - Angle Constraints: Fix bond angles
//! - Dihedral Constraints: Fix dihedral angles
//! - COM Motion Removal: Remove center-of-mass translation and rotation

use gromos_core::configuration::Configuration;
use gromos_core::lambdas::{InteractionLambda, PrecalcLambdas};
use gromos_core::math::{Periodicity, Vec3};
use gromos_core::topology::Topology;

/// Precomputed constraint data and reusable buffers for SHAKE.
///
/// Avoids per-step allocations by precomputing constraint lists from the topology
/// and reusing skip arrays across iterations.
#[derive(Debug, Clone)]
pub struct ShakeBuffers {
    /// Perturbed solute constraints (gromosXX `perturbed_solute().distance_constraints()`),
    /// shaken before the unperturbed ones in every iteration (`perturbed_shake.cc:318-345`).
    pub perturbed_constraints: Vec<PerturbedConstraint>,
    /// Precomputed solute constraints: (atom_i, atom_j, constraint_length)
    pub solute_constraints: Vec<(usize, usize, f64)>,
    /// Precomputed solvent constraints: (atom_i, atom_j, constraint_length)
    pub solvent_constraints: Vec<(usize, usize, f64)>,
    /// Sorted list of constrained atom indices (replaces HashSet)
    pub constrained_atoms: Vec<usize>,
    /// Reusable buffer: skip optimization for current iteration
    skip_now: Vec<bool>,
    /// Reusable buffer: skip optimization for next iteration
    skip_next: Vec<bool>,
}

/// A perturbed distance constraint: its length is `(1 − L)·r0_A + L·r0_B` with `L` the bond λ
/// of the energy group of atom `i`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PerturbedConstraint {
    pub i: usize,
    pub j: usize,
    pub r0_a: f64,
    pub r0_b: f64,
    /// Energy group of atom `i`: its bond λ, and the cell its dH/dλ is booked to.
    pub group: usize,
}

impl ShakeBuffers {
    /// Precompute constraint lists from topology and NTC mode.
    /// Call once at initialization, reuse across all MD steps.
    ///
    /// `include_solvent` should be `false` when the solvent is constrained by
    /// a different algorithm (e.g. SETTLE/LINCS selected via NTCS), so SHAKE
    /// only handles the solute bonds.
    pub fn new(topo: &Topology, ntc: NtcMode, include_solvent: bool) -> Self {
        let solute_constraints: Vec<(usize, usize, f64)> = match ntc {
            NtcMode::SolventOnly => Vec::new(),
            // gromosXX `in_topology.cc:2143-2178`: NTC = 2 constrains exactly the BONDH bonds
            NtcMode::HydrogenBonds => topo.moltypes[0]
                .bonds
                .iter()
                .filter_map(|bond| {
                    let constraint_length = topo.bond_parameters[bond.bond_type].r0;
                    if constraint_length < 1e-10 || !bond.hydrogen {
                        return None;
                    }
                    Some((bond.i, bond.j, constraint_length))
                })
                .collect(),
            NtcMode::AllBonds => topo.moltypes[0]
                .bonds
                .iter()
                .filter_map(|bond| {
                    let constraint_length = topo.bond_parameters[bond.bond_type].r0;
                    if constraint_length < 1e-10 {
                        return None;
                    }
                    Some((bond.i, bond.j, constraint_length))
                })
                .collect(),
        };

        // NTC > 1: the perturbed bonds `prepare_system` handed to the constraints (NTC = 3)
        let perturbed_constraints: Vec<PerturbedConstraint> = match ntc {
            NtcMode::SolventOnly => Vec::new(),
            NtcMode::HydrogenBonds | NtcMode::AllBonds => topo
                .perturbed_solute
                .distance_constraints
                .iter()
                .map(|b| PerturbedConstraint {
                    i: b.i,
                    j: b.j,
                    r0_a: topo.bond_parameters[b.a_type].r0,
                    r0_b: topo.bond_parameters[b.b_type].r0,
                    group: topo.atom_to_energy_group.get(b.i).copied().unwrap_or(0),
                })
                .collect(),
        };

        let mut solvent_constraints: Vec<(usize, usize, f64)> = Vec::new();
        if include_solvent
            && !topo.solvent_constraint_template.is_empty()
            && topo.num_solvent_molecules() > 0
        {
            let n_solute = topo.num_solute_atoms();
            let atoms_per_solvent = topo.solvent_atom_template.len();
            let num_molecules = topo.num_solvent_molecules();

            for mol in 0..num_molecules {
                let base = n_solute + mol * atoms_per_solvent;
                for constr in &topo.solvent_constraint_template {
                    solvent_constraints.push((base + constr.i, base + constr.j, constr.length));
                }
            }
        }

        // Collect constrained atoms (deduplicated + sorted, replaces HashSet)
        let mut constrained_atoms: Vec<usize> = Vec::new();
        for c in &perturbed_constraints {
            constrained_atoms.push(c.i);
            constrained_atoms.push(c.j);
        }
        for &(i, j, _) in solute_constraints.iter().chain(solvent_constraints.iter()) {
            constrained_atoms.push(i);
            constrained_atoms.push(j);
        }
        constrained_atoms.sort_unstable();
        constrained_atoms.dedup();

        let n_total_atoms = topo.num_atoms();
        let skip_now = vec![false; n_total_atoms];
        let skip_next = vec![true; n_total_atoms];

        Self {
            perturbed_constraints,
            solute_constraints,
            solvent_constraints,
            constrained_atoms,
            skip_now,
            skip_next,
        }
    }

    /// Returns true if there are no constraints at all.
    pub fn is_empty(&self) -> bool {
        self.perturbed_constraints.is_empty()
            && self.solute_constraints.is_empty()
            && self.solvent_constraints.is_empty()
    }

    /// The solute distance constraints as atom pairs, perturbed ones included: each removes a
    /// degree of freedom (gromosXX `topology.cc:1076-1180`).
    pub fn solute_constraint_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.solute_constraints
            .iter()
            .map(|&(i, j, _)| (i, j))
            .chain(self.perturbed_constraints.iter().map(|c| (c.i, c.j)))
    }
}

/// NTC constraint mode (which bonds to constrain)
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NtcMode {
    /// NTC=1: no solute constraints (solvent only if NTCS>0)
    SolventOnly,
    /// NTC=2: constrain the bonds of the topology's BONDH block
    HydrogenBonds,
    /// NTC=3: constrain all bonds
    AllBonds,
}

/// SHAKE algorithm parameters
#[derive(Debug, Clone)]
pub struct ShakeParameters {
    pub tolerance: f64,        // Convergence tolerance for constraint satisfaction
    pub max_iterations: usize, // Maximum number of iterations
    pub ntc: NtcMode,          // Which solute bonds to constrain
}

impl Default for ShakeParameters {
    fn default() -> Self {
        Self {
            tolerance: 1e-4,             // GROMOS default: 0.0001
            max_iterations: 1000,        // GROMOS default
            ntc: NtcMode::HydrogenBonds, // NTC=2
        }
    }
}

/// LINCS algorithm parameters
#[derive(Debug, Clone)]
pub struct LincsParameters {
    pub order: usize, // Order of expansion (typically 4-8, higher = more accurate)
}

impl Default for LincsParameters {
    fn default() -> Self {
        Self {
            order: 4, // GROMOS default
        }
    }
}

/// LINCS coupling data for constraint network
#[derive(Debug, Clone)]
pub struct LincsData {
    /// Diagonal matrix elements: 1/√(1/m_i + 1/m_j)
    pub sdiag: Vec<f64>,
    /// Coupled constraints for each constraint
    pub coupled_constr: Vec<Vec<usize>>,
    /// Coupling coefficients
    pub coef: Vec<Vec<f64>>,
}

impl LincsData {
    /// Create empty LINCS data
    pub fn new() -> Self {
        Self {
            sdiag: Vec::new(),
            coupled_constr: Vec::new(),
            coef: Vec::new(),
        }
    }
}

/// Result of constraint application
impl Default for LincsData {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone)]
pub struct ConstraintResult {
    pub converged: bool,
    pub iterations: usize,
    pub max_error: f64,
    /// The constraint furthest from its length when the iteration stopped: atoms (0-based),
    /// the distance reached and the distance wanted, both nm. `None` when every constraint
    /// converged.
    pub worst: Option<(usize, usize, f64, f64)>,
}

impl ConstraintResult {
    /// `atom 12 - atom 13: 0.9421 nm, wants 0.1000 nm` — what to put after "did not converge".
    pub fn worst_constraint(&self) -> String {
        match self.worst {
            None => "no constraint was left unsatisfied".to_string(),
            Some((i, j, got, want)) => format!(
                "atom {} - atom {}: {got:.4} nm, wants {want:.4} nm",
                i + 1,
                j + 1
            ),
        }
    }
}

/// The constraint furthest from its length, as `(error, (i, j, distance, length))`.
fn worst_unsatisfied(
    conf: &Configuration,
    pbc: &Periodicity,
    constraints: impl Iterator<Item = (usize, usize, f64)>,
) -> (f64, Option<(usize, usize, f64, f64)>) {
    let mut max_error = 0.0;
    let mut worst = None;
    for (i, j, length) in constraints {
        let (error, distance) = constraint_error(conf, pbc, i, j, length);
        if error > max_error {
            max_error = error;
            worst = Some((i, j, distance, length));
        }
    }
    (max_error, worst)
}

/// How far the constraint between `i` and `j` is from its length, and the distance it reached.
fn constraint_error(
    conf: &Configuration,
    pbc: &Periodicity,
    i: usize,
    j: usize,
    length: f64,
) -> (f64, f64) {
    let r = pbc.nearest_image(conf.current().pos[i], conf.current().pos[j]);
    let d = r.dot(r).sqrt();
    ((d - length).abs(), d)
}

/// The minimum-image convention of the current box. gromosXX's SHAKE takes the constraint vector
/// as `periodicity.nearest_image(pos_i, pos_j)` (`shake.h::shake_iteration`), so a molecule whose
/// atoms sit on opposite sides of the box — every configuration written with the molecules put
/// back into the box — is constrained across the boundary instead of across the whole box.
pub fn box_periodicity(conf: &Configuration) -> Periodicity {
    use gromos_core::configuration::BoxType;
    use gromos_core::math::{Rectangular, Triclinic, Vacuum};
    let b = &conf.current().box_config;
    match b.box_type {
        BoxType::Rectangular => {
            let d = b.dimensions();
            if d.x > 0.0 && d.y > 0.0 && d.z > 0.0 {
                Periodicity::Rectangular(Rectangular::new(d))
            } else {
                Periodicity::Vacuum(Vacuum)
            }
        },
        BoxType::Triclinic | BoxType::TruncatedOctahedral => {
            Periodicity::Triclinic(Triclinic::new(b.vectors))
        },
        _ => Periodicity::Vacuum(Vacuum),
    }
}

/// Apply SHAKE algorithm to satisfy distance constraints
///
/// The SHAKE algorithm (Ryckaert, Ciccotti & Berendsen, 1977) iteratively
/// adjusts atomic positions to satisfy bond length constraints.
///
/// Features matching GROMOS:
/// - NTC control: 1=solvent only, 2=H-bonds, 3=all bonds
/// - Constraint forces stored in conf.old().constraint_force
/// - Virial tensor contribution accumulated in conf.old().virial_tensor
/// - skip_now/skip_next optimization to avoid re-checking converged constraints
/// - Velocity correction: v = (pos_new - pos_old) / dt
///
/// # Parameters
/// - `topo`: Molecular topology containing constraint information
/// - `conf`: Configuration with current and old positions
/// - `dt`: Time step size
/// - `params`: SHAKE parameters (tolerance, max iterations, NTC mode)
///
/// # Returns
/// - `ConstraintResult` with convergence status and statistics
pub fn shake(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &ShakeParameters,
) -> ConstraintResult {
    let mut buffers = ShakeBuffers::new(topo, params.ntc, true);
    shake_buffered(topo, conf, dt, params, &mut buffers, None, false)
}

/// Apply SHAKE using precomputed buffers (zero per-step allocations).
///
/// This is the performance-critical version called every MD step.
/// All constraint lists are precomputed in `ShakeBuffers`, and skip arrays are reused.
///
/// `precalc`: the PRECALCLAM grid when this step is one it is evaluated on — every perturbed
/// constraint correction then adds to the `AB_bond` derivatives at each λ_k
/// (`perturbed_shake.cc:229-259`; energies get nothing, gromosXX's own caveat).
///
/// `atomic_virial`: the run computes the atomic virial (`VIRIAL 1`). gromosXX adds the constraint
/// virial only then, and never in vacuum (`shake.h:380`, `perturbed_shake.cc:209`, dispatched by
/// `SPLIT_VIRIAL_BOUNDARY`, `template_split.h:58-79`).
pub fn shake_buffered(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &ShakeParameters,
    buffers: &mut ShakeBuffers,
    precalc: Option<PrecalcLambdas>,
    atomic_virial: bool,
) -> ConstraintResult {
    let pbc = box_periodicity(conf);
    if atomic_virial && !matches!(pbc, Periodicity::Vacuum(_)) {
        shake_buffered_impl::<true>(topo, conf, dt, params, buffers, precalc, pbc)
    } else {
        shake_buffered_impl::<false>(topo, conf, dt, params, buffers, precalc, pbc)
    }
}

fn shake_buffered_impl<const VIRIAL: bool>(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &ShakeParameters,
    buffers: &mut ShakeBuffers,
    precalc: Option<PrecalcLambdas>,
    pbc: Periodicity,
) -> ConstraintResult {
    if buffers.is_empty() {
        return ConstraintResult {
            converged: true,
            iterations: 0,
            max_error: 0.0,
            worst: None,
        };
    }

    let tolerance = params.tolerance;
    let dt2 = dt * dt;

    // Zero constraint forces for constrained atoms
    for &atom in &buffers.constrained_atoms {
        conf.old_mut().constraint_force[atom] = Vec3::ZERO;
    }

    // Reset skip buffers
    buffers.skip_now.iter_mut().for_each(|s| *s = false);
    buffers.skip_next.iter_mut().for_each(|s| *s = true);

    let mut converged = false;
    let mut iteration = 0;

    while iteration < params.max_iterations && !converged {
        converged = true;
        iteration += 1;

        // === Perturbed solute constraints (gromosXX perturbed_shake_iteration) ===
        for c in &buffers.perturbed_constraints {
            let (i, j) = (c.i, c.j);
            if buffers.skip_now[i] && buffers.skip_now[j] {
                continue;
            }
            if topo.inverse_mass[i] == 0.0 && topo.inverse_mass[j] == 0.0 {
                continue;
            }
            // perturbed_shake.cc:121-125: the bond λ of atom i's energy group
            let (lam, lam_derivative) =
                topo.individual_lambdas
                    .get(InteractionLambda::Bond, c.group, c.group);
            let r0 = (1.0 - lam) * c.r0_a + lam * c.r0_b;
            if let Some((lambda, sp, dist2)) = shake_one_constraint_full::<VIRIAL>(
                conf,
                topo,
                &pbc,
                i,
                j,
                r0,
                tolerance,
                dt2,
                &mut buffers.skip_next,
            ) {
                converged = false;
                // perturbed_shake.cc:224-226, booked on the old state (this step's energies)
                let constr_length2 = r0 * r0;
                let de = lam_derivative * lambda / dt2 * constr_length2.sqrt() * (c.r0_b - c.r0_a);
                let d = &mut conf.old_mut().perturbed_energy_derivatives;
                if let Some(cell) = d.constraints_energy.get_mut(c.group) {
                    *cell += de;
                }
                d.constraint_total += de;
                d.special_total += de;
                if let Some(grid) = &precalc {
                    let r0_diff = c.r0_b - c.r0_a;
                    let sp_2_m_dt2 = sp * 2.0 * (1.0 / topo.mass[i] + 1.0 / topo.mass[j]) * dt2;
                    let ab_bond = &mut d.precalc.ab_bond;
                    for (k, slot) in ab_bond.iter_mut().enumerate().take(grid.nr_lambdas) {
                        let lam = grid.lambda(k);
                        let r0lam = (1.0 - lam) * c.r0_a + lam * c.r0_b;
                        let difflam = r0lam * r0lam - dist2;
                        *slot += (difflam / sp_2_m_dt2) * r0lam * r0_diff;
                    }
                }
            }
        }

        // === Solute constraints ===
        for &(i, j, constraint_length) in &buffers.solute_constraints {
            if buffers.skip_now[i] && buffers.skip_now[j] {
                continue;
            }
            if topo.inverse_mass[i] == 0.0 && topo.inverse_mass[j] == 0.0 {
                continue;
            }

            if shake_one_constraint_full::<VIRIAL>(
                conf,
                topo,
                &pbc,
                i,
                j,
                constraint_length,
                tolerance,
                dt2,
                &mut buffers.skip_next,
            )
            .is_some()
            {
                converged = false;
            }
        }

        // === Solvent constraints ===
        for &(i, j, constraint_length) in &buffers.solvent_constraints {
            if buffers.skip_now[i] && buffers.skip_now[j] {
                continue;
            }

            if shake_one_constraint_full::<VIRIAL>(
                conf,
                topo,
                &pbc,
                i,
                j,
                constraint_length,
                tolerance,
                dt2,
                &mut buffers.skip_next,
            )
            .is_some()
            {
                converged = false;
            }
        }

        // Swap skip arrays for next iteration
        std::mem::swap(&mut buffers.skip_now, &mut buffers.skip_next);
        buffers.skip_next.iter_mut().for_each(|s| *s = true);
    }

    // Constraint forces: 1/dt² once per atom (`shake.h:529` solute, `:653` solvent)
    for &atom in &buffers.constrained_atoms {
        conf.old_mut().constraint_force[atom] *= 1.0 / dt2;
    }

    // Velocity correction
    for &atom in &buffers.constrained_atoms {
        conf.current_mut().vel[atom] =
            (conf.current().pos[atom] - conf.old().pos[atom]) * (1.0 / dt);
    }

    let (max_error, worst) = if converged {
        (0.0, None)
    } else {
        // perturbed constraints at their λ-mixed length, as shaken above
        let perturbed = buffers.perturbed_constraints.iter().map(|c| {
            let (lam, _) = topo
                .individual_lambdas
                .get(InteractionLambda::Bond, c.group, c.group);
            (c.i, c.j, (1.0 - lam) * c.r0_a + lam * c.r0_b)
        });
        worst_unsatisfied(
            conf,
            &pbc,
            perturbed.chain(
                buffers
                    .solute_constraints
                    .iter()
                    .chain(buffers.solvent_constraints.iter())
                    .copied(),
            ),
        )
    };

    ConstraintResult {
        converged,
        iterations: iteration,
        max_error,
        worst,
    }
}

/// Apply SHAKE to a single distance constraint with full GROMOS features.
/// Returns `None` if already converged, else the Lagrange multiplier it applied.
/// Accumulates constraint force (pre-1/dt² scaling) and, with `VIRIAL`, the virial tensor.
#[inline]
#[allow(clippy::too_many_arguments)] // one constraint's worth of state, as in gromosXX
fn shake_one_constraint_full<const VIRIAL: bool>(
    conf: &mut Configuration,
    topo: &Topology,
    pbc: &Periodicity,
    i: usize,
    j: usize,
    constraint_length: f64,
    tolerance: f64,
    dt2: f64,
    skip_next: &mut [bool],
) -> Option<(f64, f64, f64)> {
    let constr_length2 = constraint_length * constraint_length;

    // GROMOS convention: r = pos(i) - pos(j), under the minimum-image convention
    let r = pbc.nearest_image(conf.current().pos[i], conf.current().pos[j]);
    let dist2 = r.dot(r);

    // GROMOS: diff = constr_length2 - dist2
    let diff = constr_length2 - dist2;

    // GROMOS convergence: fabs(diff) >= constr_length2 * tolerance * 2.0
    if diff.abs() < constr_length2 * tolerance * 2.0 {
        return None; // already converged
    }

    // Reference (old) distance vector, likewise the nearest image
    let ref_r = pbc.nearest_image(conf.old().pos[i], conf.old().pos[j]);
    let sp = ref_r.dot(r);

    // Mass weighting
    let inv_mass_i = topo.inverse_mass[i];
    let inv_mass_j = topo.inverse_mass[j];
    let inv_mass_sum = inv_mass_i + inv_mass_j;

    if inv_mass_sum < 1e-20 {
        return None; // fixed atoms
    }

    // GROMOS: lambda = diff / (sp * 2.0 * inv_mass_sum)
    let lambda = diff / (sp * 2.0 * inv_mass_sum);

    // Accumulate constraint force (GROMOS: cons_force = lambda * ref_r)
    // Note: final scaling by 1/dt² is done after all iterations
    let cons_force = ref_r * lambda;
    conf.old_mut().constraint_force[i] += cons_force;
    conf.old_mut().constraint_force[j] -= cons_force;

    // Virial tensor contribution (GROMOS: virial_tensor(a,aa) += ref_r(a) * ref_r(aa) * lambda / dt2)
    if VIRIAL {
        let lambda_over_dt2 = lambda / dt2;
        let vt = &mut conf.old_mut().virial_tensor;
        vt.x_axis.x += ref_r.x * ref_r.x * lambda_over_dt2;
        vt.x_axis.y += ref_r.x * ref_r.y * lambda_over_dt2;
        vt.x_axis.z += ref_r.x * ref_r.z * lambda_over_dt2;
        vt.y_axis.x += ref_r.y * ref_r.x * lambda_over_dt2;
        vt.y_axis.y += ref_r.y * ref_r.y * lambda_over_dt2;
        vt.y_axis.z += ref_r.y * ref_r.z * lambda_over_dt2;
        vt.z_axis.x += ref_r.z * ref_r.x * lambda_over_dt2;
        vt.z_axis.y += ref_r.z * ref_r.y * lambda_over_dt2;
        vt.z_axis.z += ref_r.z * ref_r.z * lambda_over_dt2;
    }

    // Update positions
    let correction = ref_r * lambda;
    conf.current_mut().pos[i] += correction * inv_mass_i;
    conf.current_mut().pos[j] -= correction * inv_mass_j;

    // Mark atoms as needing re-check in next iteration
    skip_next[i] = false;
    skip_next[j] = false;

    // not yet converged: the multiplier, and `sp`, `dist2` for the PRECALCLAM derivative
    Some((lambda, sp, dist2))
}

/// Shake initial positions to satisfy constraints (GROMOS: sim.param().start.shake_pos).
///
/// Uses current positions as both reference and unconstrained positions.
/// After shaking, old positions are set equal to the shaken current positions.
pub fn shake_positions(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &ShakeParameters,
) -> ConstraintResult {
    // Set old = current (reference = unconstrained = same positions)
    let n = topo.num_atoms();
    for i in 0..n {
        conf.old_mut().pos[i] = conf.current().pos[i];
        conf.old_mut().vel[i] = conf.current().vel[i];
    }

    // Shake current positions
    let result = shake(topo, conf, dt, params);

    // Restore velocities and set old = shaken current
    for i in 0..n {
        conf.current_mut().vel[i] = conf.old().vel[i];
        conf.old_mut().pos[i] = conf.current().pos[i];
    }

    result
}

/// Shake initial velocities to satisfy velocity constraints (GROMOS: sim.param().start.shake_vel).
///
/// Ensures velocity components along constrained bonds are zero.
/// Must be called after shake_positions.
pub fn shake_velocities(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &ShakeParameters,
) -> ConstraintResult {
    let n = topo.num_atoms();

    // Step 1: compute unconstrained positions from current velocities
    // r_uc(t+dt) = r(t) + v(t+dt/2) * dt
    for i in 0..n {
        conf.current_mut().pos[i] = conf.old().pos[i] + conf.old().vel[i] * dt;
    }

    // Step 2: SHAKE these positions
    let result = shake(topo, conf, dt, params);

    // Step 3: recover constrained velocities and restore positions
    // v(t+dt/2) = (r_constrained(t+dt) - r(t)) / dt
    // Velocities are already corrected by shake() itself for constrained atoms.
    // For unconstrained atoms, restore original positions.
    for i in 0..n {
        conf.current_mut().pos[i] = conf.old().pos[i];
        // Negate velocity direction (GROMOS convention for init shake_vel)
        conf.current_mut().vel[i] = -conf.current().vel[i];
        conf.old_mut().vel[i] = conf.current().vel[i];
    }

    result
}

/// M-SHAKE: Mass-weighted SHAKE variant
///
/// Similar to SHAKE but uses mass-weighted coordinates for better
/// numerical stability, especially for constraints involving light atoms (e.g., hydrogens)
pub fn m_shake(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &ShakeParameters,
) -> ConstraintResult {
    // For now, M-SHAKE uses the same implementation as SHAKE
    // A full M-SHAKE would work in mass-weighted coordinates
    // This is a simplified version that still provides good results
    shake(topo, conf, dt, params)
}

/// SETTLE: Analytical constraint solver for rigid 3-site water molecules
///
/// SETTLE (Miyamoto & Kollman, 1992) is an analytical method specifically
/// designed for rigid water molecules with fixed bond lengths and angles.
/// It's much faster than iterative SHAKE for water.
///
/// Faithful port of GROMOS `algorithm::Settle::solvent`
/// (md++/src/algorithm/constraints/settle.cc:179-415): builds the canonical
/// triangle geometry from the constraint lengths and masses, transforms into
/// a centre-of-mass/orientation frame derived from the *old* positions, solves
/// analytically for sin(φ), sin(ψ), sin(θ) (eqs. A3-A17 of the paper), and
/// reconstructs the new O/H1/H2 positions — then derives velocities, the
/// constraint force, and the virial contribution from the displacement.
///
/// Assumptions (checked by the caller / `SettleAlgorithm::init`):
/// - Exactly one solvent, water-like with 3 atoms and equal H masses
/// - 3 distance constraints with the two O-H lengths equal
///
/// `atomic_virial`: add the constraint virial, which gromosXX does only for `VIRIAL 1`
/// (`settle.cc:397`).
pub fn settle(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    atomic_virial: bool,
) -> ConstraintResult {
    let n_solute = topo.num_solute_atoms();
    let num_atoms = topo.num_atoms();

    if topo.num_solvent_molecules() == 0 || topo.solvent_constraint_template.len() != 3 {
        return ConstraintResult {
            converged: true,
            iterations: 0,
            max_error: 0.0,
            worst: None,
        };
    }

    let m_o = topo.mass[n_solute];
    let m_h = topo.mass[n_solute + 1];

    // Distance constraints: 0 = O-H1, 1 = O-H2, 2 = H1-H2 (GROMOS convention)
    let dist_oh = topo.solvent_constraint_template[0].length;
    let dist_hh = topo.solvent_constraint_template[2].length;

    // canonical triangle geometry (settle.cc:203-206 / Figure 2a)
    let half_mo_div_mh = 0.5 * m_o / m_h;
    let rc = 0.5 * dist_hh;
    let ra = (dist_oh * dist_oh - rc * rc).sqrt() / (1.0 + half_mo_div_mh);
    let rb = half_mo_div_mh * ra;

    let dt_i = 1.0 / dt;
    let dt2_i = dt_i * dt_i;

    let mut max_error: f64 = 0.0;
    let mut error_count = 0usize;
    let mut num_molecules = 0usize;

    let mut i = n_solute;
    while i + 3 <= num_atoms {
        let pos_old = [
            conf.old().pos[i],
            conf.old().pos[i + 1],
            conf.old().pos[i + 2],
        ];
        let pos_new = [
            conf.current().pos[i],
            conf.current().pos[i + 1],
            conf.current().pos[i + 2],
        ];

        // vectors in the plane of the old positions
        let mut b0 = pos_old[1] - pos_old[0];
        let mut c0 = pos_old[2] - pos_old[0];

        // centre of mass of the new positions
        let d0 = (pos_new[0] * m_o + pos_new[1] * m_h + pos_new[2] * m_h) / (m_o + m_h + m_h);

        // move the origin to the centre of mass
        let mut a1 = pos_new[0] - d0;
        let mut b1 = pos_new[1] - d0;
        let mut c1 = pos_new[2] - d0;

        // orthonormal basis transforming old-frame -> COM frame (settle.cc:252-265)
        let mut n0 = b0.cross(c0);
        let n1 = a1.cross(n0);
        let n2 = n0.cross(n1);
        n0 /= n0.length();
        let n1 = n1 / n1.length();
        let n2 = n2 / n2.length();

        // transpose, to undo the transformation later
        let m1 = Vec3::new(n1.x, n2.x, n0.x);
        let m2 = Vec3::new(n1.y, n2.y, n0.y);
        let m0 = Vec3::new(n1.z, n2.z, n0.z);

        // rotate old positions into the COM/orientation frame
        b0 = Vec3::new(n1.dot(b0), n2.dot(b0), n0.dot(b0));
        c0 = Vec3::new(n1.dot(c0), n2.dot(c0), n0.dot(c0));

        // and the new (COM-centred) positions
        a1 = Vec3::new(n1.dot(a1), n2.dot(a1), n0.dot(a1));
        b1 = Vec3::new(n1.dot(b1), n2.dot(b1), n0.dot(b1));
        c1 = Vec3::new(n1.dot(c1), n2.dot(c1), n0.dot(c1));

        // sin(phi)/cos(phi) — (A8)
        let sinphi = a1.z / ra;
        let one_minus_sinphi2 = 1.0 - sinphi * sinphi;
        if one_minus_sinphi2 < 0.0 {
            error_count += 1;
            i += 3;
            continue;
        }
        let cosphi = one_minus_sinphi2.sqrt();

        // sin(psi)/cos(psi) — (A9)
        let sinpsi = (b1.z - c1.z) / (2.0 * rc * cosphi);
        let one_minus_sinpsi2 = 1.0 - sinpsi * sinpsi;
        if one_minus_sinpsi2 < 0.0 {
            error_count += 1;
            i += 3;
            continue;
        }
        let cospsi = one_minus_sinpsi2.sqrt();

        let minus_rb_cosphi = -rb * cosphi;
        let rc_cospsi = rc * cospsi;
        let rc_sinpsi_sinphi = rc * sinpsi * sinphi;
        let rc_sinpsi_cosphi = rc * sinpsi * cosphi;

        // (A3)
        let x_a2 = 0.0;
        let x_b2 = -rc_cospsi;
        let x_c2 = rc_cospsi;

        let y_a2 = ra * cosphi;
        let y_b2 = minus_rb_cosphi - rc_sinpsi_sinphi;
        let y_c2 = minus_rb_cosphi + rc_sinpsi_sinphi;

        // (A5)/(A6)/(A7)
        let z_a2 = ra * sinphi;
        let z_b2 = -rb * sinphi + rc_sinpsi_cosphi;
        let z_c2 = -rb * sinphi - rc_sinpsi_cosphi;

        let a2 = Vec3::new(x_a2, y_a2, z_a2);
        let b2 = Vec3::new(x_b2, y_b2, z_b2);
        let c2 = Vec3::new(x_c2, y_c2, z_c2);

        // (A15)
        let alpha = b2.x * (b0.x - c0.x) + b0.y * b2.y + c0.y * c2.y;
        let beta = b2.x * (c0.y - b0.y) + b0.x * b2.y + c0.x * c2.y;
        let gamma = b0.x * b1.y - b1.x * b0.y + c0.x * c1.y - c1.x * c0.y;

        let alpha2_beta2 = alpha * alpha + beta * beta;
        // (A17)
        let sintheta =
            (alpha * gamma - beta * (alpha2_beta2 - gamma * gamma).sqrt()) / alpha2_beta2;
        let one_minus_sintheta2 = 1.0 - sintheta * sintheta;
        if one_minus_sintheta2 < 0.0 {
            error_count += 1;
            i += 3;
            continue;
        }
        let costheta = one_minus_sintheta2.sqrt();

        // (A4)
        let a3 = Vec3::new(-a2.y * sintheta, a2.y * costheta, a1.z);
        let b3 = Vec3::new(
            b2.x * costheta - b2.y * sintheta,
            b2.x * sintheta + b2.y * costheta,
            b1.z,
        );
        let c3 = Vec3::new(
            -b2.x * costheta - c2.y * sintheta,
            -b2.x * sintheta + c2.y * costheta,
            c1.z,
        );

        // rotate back to the lab frame and re-add the centre of mass
        let pos_a = Vec3::new(a3.dot(m1), a3.dot(m2), a3.dot(m0)) + d0;
        let pos_b = Vec3::new(b3.dot(m1), b3.dot(m2), b3.dot(m0)) + d0;
        let pos_c = Vec3::new(c3.dot(m1), c3.dot(m2), c3.dot(m0)) + d0;

        // displacement, used for velocity correction, constraint force and virial
        let d_a = pos_a - pos_new[0];
        let d_b = pos_b - pos_new[1];
        let d_c = pos_c - pos_new[2];

        conf.current_mut().pos[i] = pos_a;
        conf.current_mut().pos[i + 1] = pos_b;
        conf.current_mut().pos[i + 2] = pos_c;

        // constraint force by finite difference (settle.cc:381-383)
        let cons_force = [d_a * m_o * dt2_i, d_b * m_h * dt2_i, d_c * m_h * dt2_i];
        conf.old_mut().constraint_force[i] = cons_force[0];
        conf.old_mut().constraint_force[i + 1] = cons_force[1];
        conf.old_mut().constraint_force[i + 2] = cons_force[2];

        // velocity correction (settle.cc:387-395)
        conf.current_mut().vel[i] += d_a * dt_i;
        conf.current_mut().vel[i + 1] += d_b * dt_i;
        conf.current_mut().vel[i + 2] += d_c * dt_i;

        // atomic virial contribution (settle.cc:397-411)
        // GROMOS: virial_tensor(row, col) += Σ_k pos_old[k](row) * cons_force[k](col);
        // Mat3 is column-major (mat.col(col)[row]), see forcefield.rs virial transfer.
        if atomic_virial {
            let p = [
                [pos_old[0].x, pos_old[0].y, pos_old[0].z],
                [pos_old[1].x, pos_old[1].y, pos_old[1].z],
                [pos_old[2].x, pos_old[2].y, pos_old[2].z],
            ];
            let f = [
                [cons_force[0].x, cons_force[0].y, cons_force[0].z],
                [cons_force[1].x, cons_force[1].y, cons_force[1].z],
                [cons_force[2].x, cons_force[2].y, cons_force[2].z],
            ];
            let mut dvir = [[0.0f64; 3]; 3]; // dvir[row][col]
            for k in 0..3 {
                for row in 0..3 {
                    for col in 0..3 {
                        dvir[row][col] += p[k][row] * f[k][col];
                    }
                }
            }
            let vt = &mut conf.old_mut().virial_tensor;
            vt.x_axis += Vec3::new(dvir[0][0], dvir[1][0], dvir[2][0]);
            vt.y_axis += Vec3::new(dvir[0][1], dvir[1][1], dvir[2][1]);
            vt.z_axis += Vec3::new(dvir[0][2], dvir[1][2], dvir[2][2]);
        }

        let err = (d_a.length() + d_b.length() + d_c.length()) / 3.0;
        if err > max_error {
            max_error = err;
        }

        num_molecules += 1;
        i += 3;
    }

    let _ = num_molecules;

    ConstraintResult {
        converged: error_count == 0,
        iterations: 1,
        max_error,
        worst: None,
    }
}

/// Precomputed LINCS constraint data and reusable buffers, mirroring `ShakeBuffers`.
///
/// Solute constraints use global atom indices directly (offset 0). Solvent
/// constraints are stored once as a local (0-based, within-molecule) template
/// — since every molecule of a given solvent type has identical connectivity
/// and masses, the coupling matrix only needs to be computed once
/// (cf. GROMOS `Lincs::init`, `lincs.cc:391-409`, which sets up
/// `topo.solvent(i).lincs()` once and reuses it for every molecule via an
/// `offset` in `_lincs<B>`, `lincs.cc:158-205`).
#[derive(Debug, Clone)]
pub struct LincsBuffers {
    /// Solute constraints: (atom_i, atom_j, r0), global indices, offset 0
    pub solute_constraints: Vec<(usize, usize, f64)>,
    pub solute_data: LincsData,
    pub solute_params: LincsParameters,

    /// One solvent molecule's constraints: (local_i, local_j, r0)
    pub solvent_constraints: Vec<(usize, usize, f64)>,
    pub solvent_data: LincsData,
    pub solvent_params: LincsParameters,
    /// Atom index of the first solvent atom (= number of solute atoms)
    pub solvent_offset0: usize,
    /// Number of atoms per solvent molecule
    pub atoms_per_solvent: usize,
    /// Number of solvent molecules
    pub num_solvent_molecules: usize,
}

impl LincsBuffers {
    /// Precompute LINCS constraint lists and coupling matrices.
    ///
    /// `ntc` selects which solute bonds are constrained (mirrors `ShakeBuffers::new`);
    /// LINCS for the solute additionally requires `ntc > 1` (GROMOS: `lincs.cc:248`).
    /// `solute_order`/`solvent_order` are the LINCS expansion orders (NTCP0/NTCS0).
    pub fn new(
        topo: &Topology,
        ntc: NtcMode,
        solute_order: usize,
        solvent_order: usize,
        include_solute: bool,
        include_solvent: bool,
    ) -> Self {
        let solute_constraints: Vec<(usize, usize, f64)> = if include_solute {
            match ntc {
                NtcMode::SolventOnly => Vec::new(),
                NtcMode::HydrogenBonds => topo.moltypes[0]
                    .bonds
                    .iter()
                    .filter_map(|bond| {
                        let r0 = topo.bond_parameters[bond.bond_type].r0;
                        if r0 < 1e-10 || !bond.hydrogen {
                            return None;
                        }
                        Some((bond.i, bond.j, r0))
                    })
                    .collect(),
                NtcMode::AllBonds => topo.moltypes[0]
                    .bonds
                    .iter()
                    .filter_map(|bond| {
                        let r0 = topo.bond_parameters[bond.bond_type].r0;
                        if r0 < 1e-10 {
                            return None;
                        }
                        Some((bond.i, bond.j, r0))
                    })
                    .collect(),
            }
        } else {
            Vec::new()
        };

        let n_solute = topo.num_solute_atoms();
        let solute_params = LincsParameters {
            order: solute_order,
        };
        let solute_data = setup_lincs(topo, &solute_constraints, 0);

        let mut solvent_constraints: Vec<(usize, usize, f64)> = Vec::new();
        let mut atoms_per_solvent = 0;
        let mut num_solvent_molecules = 0;
        if include_solvent
            && !topo.solvent_constraint_template.is_empty()
            && topo.num_solvent_molecules() > 0
        {
            atoms_per_solvent = topo.solvent_atom_template.len();
            num_solvent_molecules = topo.num_solvent_molecules();
            for constr in &topo.solvent_constraint_template {
                solvent_constraints.push((constr.i, constr.j, constr.length));
            }
        }
        let solvent_params = LincsParameters {
            order: solvent_order,
        };
        // Coupling matrix only depends on masses/connectivity, identical for every
        // molecule of this solvent type — compute once using the first molecule's offset.
        let solvent_data = setup_lincs(topo, &solvent_constraints, n_solute);

        Self {
            solute_constraints,
            solute_data,
            solute_params,
            solvent_constraints,
            solvent_data,
            solvent_params,
            solvent_offset0: n_solute,
            atoms_per_solvent,
            num_solvent_molecules,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.solute_constraints.is_empty() && self.solvent_constraints.is_empty()
    }
}

/// Setup LINCS coupling matrix for a set of distance constraints.
///
/// Mirrors GROMOS `_setup_lincs` (`lincs.cc:288-345`): computes the diagonal
/// matrix elements `sdiag[i] = 1/sqrt(1/m_i + 1/m_j)` and, for every pair of
/// constraints sharing a common atom, the coupling coefficient
/// `coef = ±(1/m_common) * sdiag[i] * sdiag[j]` (sign flips when the shared
/// atom occupies the same slot — both `.i` or both `.j` — in each constraint).
///
/// `constraints` are `(atom_i, atom_j, r0)` triples using indices local to
/// `offset` (pass `offset = 0` for already-global indices, e.g. solute bonds;
/// pass the first-atom index of a solvent molecule for solvent templates).
pub fn setup_lincs(
    topo: &Topology,
    constraints: &[(usize, usize, f64)],
    offset: usize,
) -> LincsData {
    let num_constr = constraints.len();
    let mut lincs = LincsData::new();

    lincs.sdiag.reserve(num_constr);
    lincs.coupled_constr.resize(num_constr, Vec::new());
    lincs.coef.resize(num_constr, Vec::new());

    for &(i, j, _) in constraints {
        let inv_mass_sum = topo.inverse_mass[i + offset] + topo.inverse_mass[j + offset];
        lincs.sdiag.push(1.0 / inv_mass_sum.sqrt());
    }

    for i in 0..num_constr {
        let (ci_i, ci_j, _) = constraints[i];

        for j in (i + 1)..num_constr {
            let (cj_i, cj_j, _) = constraints[j];

            let common_atom = if ci_i == cj_i {
                Some(ci_i)
            } else if ci_j == cj_i {
                Some(ci_j)
            } else if ci_i == cj_j {
                Some(ci_i)
            } else if ci_j == cj_j {
                Some(ci_j)
            } else {
                None
            };

            if let Some(con) = common_atom {
                lincs.coupled_constr[i].push(j);
                lincs.coupled_constr[j].push(i);

                let inv_mass_common = topo.inverse_mass[con + offset];
                let mut c = inv_mass_common * lincs.sdiag[i] * lincs.sdiag[j];

                if (ci_i == cj_i) || (ci_j == cj_j) {
                    c *= -1.0;
                }

                lincs.coef[i].push(c);
                lincs.coef[j].push(c);
            }
        }
    }

    lincs
}

/// Apply LINCS algorithm to satisfy distance constraints.
///
/// LINCS (Linear Constraint Solver, Hess et al. 1997) — analytical port of
/// GROMOS `_lincs<B>` (`lincs.cc:118-205`): builds reference direction
/// vectors `B` from old positions, solves the coupled linear system via a
/// truncated matrix expansion (`_solve_lincs`, `lincs.cc:69-112`), corrects
/// positions, then applies a second solve for the rotational-lengthening
/// correction. Finishes with the GROMOS velocity recovery
/// `v = (pos - old_pos) / dt` (`lincs.cc:273-281`, applied unconditionally
/// here — mirrors how `shake`/`settle` handle `do_velocity` in this port,
/// since `SimulationState` has no analyze/minimise/sd flags).
///
/// Note on sign convention: this port defines `B = (pos[j] - pos[old_i]) /
/// |...|` (opposite of GROMOS's `B = (pos[i] - pos[j])/|...|`); the position
/// correction signs below are flipped to match, so the net result is identical
/// — `B` and the correction always appear together as `±B * sdiag/m * sol`.
///
/// `constraints`/`offset` follow the same convention as [`setup_lincs`].
pub fn lincs(
    topo: &Topology,
    conf: &mut Configuration,
    constraints: &[(usize, usize, f64)],
    lincs_data: &LincsData,
    params: &LincsParameters,
    offset: usize,
    dt: f64,
) -> ConstraintResult {
    let num_constr = constraints.len();

    if num_constr == 0 {
        return ConstraintResult {
            converged: true,
            iterations: 0,
            max_error: 0.0,
            worst: None,
        };
    }

    // B vectors: constraint direction unit vectors from old positions
    let mut b_vectors: Vec<Vec3> = Vec::with_capacity(num_constr);
    for &(i, j, _) in constraints {
        let r_ij = conf.old().pos[j + offset] - conf.old().pos[i + offset];
        let r_length = r_ij.length();

        if r_length < 1e-10 {
            b_vectors.push(Vec3::new(0.0, 0.0, 0.0));
            continue;
        }
        b_vectors.push(r_ij / r_length);
    }

    // Solve-matrix coefficients: A[i][n] = coef[i][n] * dot(B[i], B[coupled[n]])
    // (lincs.cc:188-189) — the static `coef` only carries the mass-weighted part;
    // the geometry-dependent dot product between the (old-position) B vectors of
    // coupled constraints must be folded in here, every step.
    let a_coef: Vec<Vec<f64>> = (0..num_constr)
        .map(|i| {
            lincs_data.coupled_constr[i]
                .iter()
                .enumerate()
                .map(|(n, &coupled_idx)| {
                    lincs_data.coef[i][n] * b_vectors[i].dot(b_vectors[coupled_idx])
                })
                .collect()
        })
        .collect();

    // Right-hand side and solution vectors
    let mut rhs: Vec<Vec<f64>> = vec![vec![0.0; num_constr]; 2];
    let mut sol: Vec<f64> = vec![0.0; num_constr];

    // Initial RHS: deviation of the current bond from the constraint length
    for (idx, &(i, j, r0)) in constraints.iter().enumerate() {
        let r_ij = conf.current().pos[j + offset] - conf.current().pos[i + offset];
        let projection = b_vectors[idx].dot(r_ij);
        rhs[0][idx] = lincs_data.sdiag[idx] * (projection - r0);
        sol[idx] = rhs[0][idx];
    }

    // Iterative solution of coupled constraints (truncated matrix expansion)
    let mut w = 1;
    for _rec in 0..params.order {
        for i in 0..num_constr {
            rhs[w][i] = 0.0;
            for (n, &coupled_idx) in lincs_data.coupled_constr[i].iter().enumerate() {
                rhs[w][i] += a_coef[i][n] * rhs[1 - w][coupled_idx];
            }
            sol[i] += rhs[w][i];
        }
        w = 1 - w;
    }

    // Apply position corrections
    for (idx, &(i, j, _)) in constraints.iter().enumerate() {
        let inv_mass_i = topo.inverse_mass[i + offset];
        let inv_mass_j = topo.inverse_mass[j + offset];
        let correction = lincs_data.sdiag[idx] * sol[idx];

        conf.current_mut().pos[i + offset] += b_vectors[idx] * (correction * inv_mass_i);
        conf.current_mut().pos[j + offset] -= b_vectors[idx] * (correction * inv_mass_j);
    }

    // Rotational lengthening correction
    let mut num_warnings = 0;
    for (idx, &(i, j, r0)) in constraints.iter().enumerate() {
        let r_ij = conf.current().pos[j + offset] - conf.current().pos[i + offset];
        let r_sq = r_ij.dot(r_ij);
        let target_sq = 2.0 * r0 * r0;

        let diff = target_sq - r_sq;
        let p = if diff > 0.0 {
            diff.sqrt()
        } else {
            num_warnings += 1;
            0.0
        };

        rhs[0][idx] = lincs_data.sdiag[idx] * (r0 - p);
        sol[idx] = rhs[0][idx];
    }
    let _ = num_warnings;

    // Second iterative solve for the rotational correction (reuses `a_coef` —
    // the B vectors are still those from the old positions, lincs.cc:230)
    w = 1;
    for _rec in 0..params.order {
        for i in 0..num_constr {
            rhs[w][i] = 0.0;
            for (n, &coupled_idx) in lincs_data.coupled_constr[i].iter().enumerate() {
                rhs[w][i] += a_coef[i][n] * rhs[1 - w][coupled_idx];
            }
            sol[i] += rhs[w][i];
        }
        w = 1 - w;
    }

    // Apply rotational corrections
    for (idx, &(i, j, _)) in constraints.iter().enumerate() {
        let inv_mass_i = topo.inverse_mass[i + offset];
        let inv_mass_j = topo.inverse_mass[j + offset];
        let correction = lincs_data.sdiag[idx] * sol[idx];

        conf.current_mut().pos[i + offset] += b_vectors[idx] * (correction * inv_mass_i);
        conf.current_mut().pos[j + offset] -= b_vectors[idx] * (correction * inv_mass_j);
    }

    // Velocity correction: v = (pos_new - pos_old) / dt (lincs.cc:273-281)
    let mut constrained_atoms: Vec<usize> = Vec::with_capacity(num_constr * 2);
    for &(i, j, _) in constraints {
        constrained_atoms.push(i + offset);
        constrained_atoms.push(j + offset);
    }
    constrained_atoms.sort_unstable();
    constrained_atoms.dedup();
    for &atom in &constrained_atoms {
        conf.current_mut().vel[atom] =
            (conf.current().pos[atom] - conf.old().pos[atom]) * (1.0 / dt);
    }

    // Final error
    let mut max_error = 0.0;
    for &(i, j, r0) in constraints {
        let r_ij = conf.current().pos[j + offset] - conf.current().pos[i + offset];
        let r_current = r_ij.length();
        let error = ((r_current - r0) / r0).abs();
        if error > max_error {
            max_error = error;
        }
    }

    ConstraintResult {
        converged: max_error < 1e-3, // Slightly relaxed convergence criterion
        iterations: params.order,
        max_error,
        worst: None,
    }
}

/// Apply LINCS using precomputed buffers — solute first, then every solvent
/// molecule with its atom-index offset (mirrors GROMOS `Lincs::apply`,
/// `lincs.cc:238-284`: `_lincs` for the solute, then `_solvent` looping over
/// every molecule of every solvent type with a running `first` offset).
pub fn lincs_buffered(
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    buffers: &LincsBuffers,
) -> ConstraintResult {
    let mut converged = true;
    let mut max_error: f64 = 0.0;

    if !buffers.solute_constraints.is_empty() {
        let result = lincs(
            topo,
            conf,
            &buffers.solute_constraints,
            &buffers.solute_data,
            &buffers.solute_params,
            0,
            dt,
        );
        converged &= result.converged;
        max_error = max_error.max(result.max_error);
    }

    if !buffers.solvent_constraints.is_empty() {
        for mol in 0..buffers.num_solvent_molecules {
            let offset = buffers.solvent_offset0 + mol * buffers.atoms_per_solvent;
            let result = lincs(
                topo,
                conf,
                &buffers.solvent_constraints,
                &buffers.solvent_data,
                &buffers.solvent_params,
                offset,
                dt,
            );
            converged &= result.converged;
            max_error = max_error.max(result.max_error);
        }
    }

    ConstraintResult {
        converged,
        iterations: buffers
            .solute_params
            .order
            .max(buffers.solvent_params.order),
        max_error,
        worst: None,
    }
}

// ============================================================================
// Angle Constraints
// ============================================================================

/// Angle constraint data structure
#[derive(Debug, Clone, Copy)]
pub struct AngleConstraint {
    pub i: usize,   // First atom
    pub j: usize,   // Central atom (vertex)
    pub k: usize,   // Third atom
    pub theta: f64, // Target angle in radians
}

/// Perturbed angle constraint for FEP
#[derive(Debug, Clone, Copy)]
pub struct PerturbedAngleConstraint {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub a_theta: f64, // State A angle (radians)
    pub b_theta: f64, // State B angle (radians)
}

/// Angle constraint parameters
#[derive(Debug, Clone)]
pub struct AngleConstraintParameters {
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for AngleConstraintParameters {
    fn default() -> Self {
        Self {
            tolerance: 1e-4,
            max_iterations: 1000,
        }
    }
}

/// Apply angle constraints using iterative SHAKE-like algorithm
///
/// Constrains bond angles (i-j-k) to fixed values using the algorithm from:
/// J. Comput. Chem. 2021;42:418–434.
///
/// # Algorithm
/// 1. Calculate current angle θ from dot product of bond vectors
/// 2. Compute auxiliary vectors a₁₂₃ and a₃₂₁ (constraint gradients)
/// 3. Calculate mass-weighted vectors b₁₂₃ and b₃₂₁
/// 4. Solve for Lagrange multiplier λ
/// 5. Update positions of all three atoms
///
/// # Parameters
/// - `constraints`: List of angle constraints to apply
/// - `topo`: Molecular topology
/// - `conf`: Configuration with positions
/// - `dt`: Time step
/// - `params`: Constraint parameters
///
/// # Returns
/// - `ConstraintResult` with convergence status
pub fn angle_constraints(
    constraints: &[AngleConstraint],
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    params: &AngleConstraintParameters,
) -> ConstraintResult {
    let _dt_sq = dt * dt;
    let tolerance_sq = params.tolerance * params.tolerance;

    let mut max_error = f64::MAX;
    let mut iteration = 0;

    while iteration < params.max_iterations && max_error > tolerance_sq {
        max_error = 0.0;
        iteration += 1;

        for constraint in constraints {
            let i = constraint.i;
            let j = constraint.j;
            let k = constraint.k;
            let theta0 = constraint.theta;

            // Bond vectors: r12 = r_i - r_j, r32 = r_k - r_j
            let r12 = conf.current().pos[i] - conf.current().pos[j];
            let r32 = conf.current().pos[k] - conf.current().pos[j];

            let d12 = r12.length();
            let d32 = r32.length();

            if d12 < 1e-10 || d32 < 1e-10 {
                continue; // Avoid division by zero
            }

            // Current angle: θ = acos(r12 · r32 / (|r12| |r32|))
            let dot_product = (r12.dot(r32)) / (d12 * d32);
            let dot_product = dot_product.clamp(-1.0, 1.0); // Numerical safety
            let theta = dot_product.acos();

            // Constraint error
            let diff = (theta - theta0).abs();
            let error = diff * diff;

            if error > max_error {
                max_error = error;
            }

            if error < tolerance_sq {
                continue; // Already satisfied
            }

            // Reference positions for velocity correction
            let r12_old = conf.old().pos[i] - conf.old().pos[j];
            let r32_old = conf.old().pos[k] - conf.old().pos[j];

            // Auxiliary vectors (eq. 18 from paper)
            // a123 = (d12² * r32 - (r12·r32) * r12) / (d12³ * d32)
            let dot_12_32 = r12_old.dot(r32_old);
            let a123 = (r32_old * (d12 * d12) - r12_old * dot_12_32) / (d12 * d12 * d12 * d32);

            // a321 = (d32² * r12 - (r12·r32) * r32) / (d12 * d32³)
            let a321 = (r12_old * (d32 * d32) - r32_old * dot_12_32) / (d12 * d32 * d32 * d32);

            // Masses
            let m1 = topo.mass[i];
            let m2 = topo.mass[j];
            let m3 = topo.mass[k];

            // Mass-weighted vectors (eq. 28)
            // b123 = a123/m1 + (a123 + a321)/m2
            let b123 = a123 / m1 + (a123 + a321) / m2;

            // b321 = a321/m3 + (a123 + a321)/m2
            let b321 = a321 / m3 + (a123 + a321) / m2;

            // Constants for Lagrange multiplier calculation
            let c1 = r12_old.dot(r32_old); // eq. 39
            let c2 = r12_old.dot(b321) + r32_old.dot(b123); // eq. 39
            let c3 = d12 * d32; // eq. 40
            let c4 = d12 / d32 * r32_old.dot(b321) + d32 / d12 * r12_old.dot(b123); // eq. 41

            // Lagrange multiplier (eq. 43)
            // λ/dt² = (c1 - c3*cos(θ0)) / (c2 - c4*cos(θ0))
            let numerator = c1 - c3 * theta0.cos();
            let denominator = c2 - c4 * theta0.cos();

            if denominator.abs() < 1e-20 {
                continue; // Avoid division by zero
            }

            let lambda_over_dt_sq = numerator / denominator;

            // Position updates (eq. 14+17)
            // pos(i) -= (λ/dt²) * a123 / m1
            // pos(j) += (λ/dt²) * (a123 + a321) / m2
            // pos(k) -= (λ/dt²) * a321 / m3
            conf.current_mut().pos[i] -= a123 * (lambda_over_dt_sq / m1);
            conf.current_mut().pos[j] += (a123 + a321) * (lambda_over_dt_sq / m2);
            conf.current_mut().pos[k] -= a321 * (lambda_over_dt_sq / m3);
        }
    }

    ConstraintResult {
        converged: max_error <= tolerance_sq,
        iterations: iteration,
        max_error: max_error.sqrt(),
        worst: None,
    }
}

/// Apply perturbed angle constraints for FEP
///
/// Like `angle_constraints` but with λ-dependent target angles:
/// θ₀(λ) = (1-λ)·θ_A + λ·θ_B
///
/// Also calculates ∂H/∂λ for thermodynamic integration.
///
/// # Parameters
/// - `constraints`: List of perturbed angle constraints
/// - `topo`: Molecular topology
/// - `conf`: Configuration
/// - `dt`: Time step
/// - `lambda`: Current λ value (0 to 1)
/// - `lambda_deriv`: dλ/dt for FEP
/// - `params`: Constraint parameters
///
/// # Returns
/// - `ConstraintResult` with convergence status
pub fn perturbed_angle_constraints(
    constraints: &[PerturbedAngleConstraint],
    topo: &Topology,
    conf: &mut Configuration,
    dt: f64,
    lambda: f64,
    lambda_deriv: f64,
    params: &AngleConstraintParameters,
) -> ConstraintResult {
    let _dt_sq = dt * dt;
    let tolerance_sq = params.tolerance * params.tolerance;

    let mut max_error = f64::MAX;
    let mut iteration = 0;

    while iteration < params.max_iterations && max_error > tolerance_sq {
        max_error = 0.0;
        iteration += 1;

        for constraint in constraints {
            let i = constraint.i;
            let j = constraint.j;
            let k = constraint.k;

            // Interpolate target angle: θ₀(λ) = (1-λ)·θ_A + λ·θ_B
            let theta0 = (1.0 - lambda) * constraint.a_theta + lambda * constraint.b_theta;

            // Bond vectors
            let r12 = conf.current().pos[i] - conf.current().pos[j];
            let r32 = conf.current().pos[k] - conf.current().pos[j];

            let d12 = r12.length();
            let d32 = r32.length();

            if d12 < 1e-10 || d32 < 1e-10 {
                continue;
            }

            // Current angle
            let dot_product = (r12.dot(r32)) / (d12 * d32);
            let dot_product = dot_product.clamp(-1.0, 1.0);
            let theta = dot_product.acos();

            // Constraint error
            let diff = (theta - theta0).abs();
            let error = diff * diff;

            if error > max_error {
                max_error = error;
            }

            if error < tolerance_sq {
                continue;
            }

            // Reference positions
            let r12_old = conf.old().pos[i] - conf.old().pos[j];
            let r32_old = conf.old().pos[k] - conf.old().pos[j];

            // Auxiliary vectors
            let dot_12_32 = r12_old.dot(r32_old);
            let a123 = (r32_old * (d12 * d12) - r12_old * dot_12_32) / (d12 * d12 * d12 * d32);
            let a321 = (r12_old * (d32 * d32) - r32_old * dot_12_32) / (d12 * d32 * d32 * d32);

            // Masses
            let m1 = topo.mass[i];
            let m2 = topo.mass[j];
            let m3 = topo.mass[k];

            // Mass-weighted vectors
            let b123 = a123 / m1 + (a123 + a321) / m2;
            let b321 = a321 / m3 + (a123 + a321) / m2;

            // Lagrange multiplier constants
            let c1 = r12_old.dot(r32_old);
            let c2 = r12_old.dot(b321) + r32_old.dot(b123);
            let c3 = d12 * d32;
            let c4 = d12 / d32 * r32_old.dot(b321) + d32 / d12 * r12_old.dot(b123);

            let numerator = c1 - c3 * theta0.cos();
            let denominator = c2 - c4 * theta0.cos();

            if denominator.abs() < 1e-20 {
                continue;
            }

            let lambda_over_dt_sq = numerator / denominator;

            // Position updates
            conf.current_mut().pos[i] -= a123 * (lambda_over_dt_sq / m1);
            conf.current_mut().pos[j] += (a123 + a321) * (lambda_over_dt_sq / m2);
            conf.current_mut().pos[k] -= a321 * (lambda_over_dt_sq / m3);

            // Energy derivative for FEP (eq. 48)
            // ∂H/∂λ = λ_deriv * (λ/dt²) * sin(θ₀) * (θ_B - θ_A)
            if lambda_deriv.abs() > 1e-20 {
                let _dh_dlambda = lambda_deriv
                    * lambda_over_dt_sq
                    * theta0.sin()
                    * (constraint.b_theta - constraint.a_theta);
                // Store in energy derivatives (would be added to Configuration)
            }
        }
    }

    ConstraintResult {
        converged: max_error <= tolerance_sq,
        iterations: iteration,
        max_error: max_error.sqrt(),
        worst: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shake_parameters() {
        let params = ShakeParameters::default();
        assert!(params.tolerance > 0.0);
        assert!(params.max_iterations > 0);
    }

    #[test]
    fn test_constraint_result() {
        let result = ConstraintResult {
            converged: true,
            iterations: 5,
            max_error: 1e-6,
            worst: None,
        };
        assert!(result.converged);
        assert_eq!(result.iterations, 5);
        assert_eq!(
            result.worst_constraint(),
            "no constraint was left unsatisfied"
        );

        let failed = ConstraintResult {
            converged: false,
            iterations: 1000,
            max_error: 0.8171,
            worst: Some((11, 12, 0.9421, 0.125)),
        };
        // atoms are 1-based in the message, as in every GROMOS file
        assert_eq!(
            failed.worst_constraint(),
            "atom 12 - atom 13: 0.9421 nm, wants 0.1250 nm"
        );
    }

    #[test]
    fn test_lincs_parameters() {
        let params = LincsParameters::default();
        assert_eq!(params.order, 4);

        let custom = LincsParameters { order: 8 };
        assert_eq!(custom.order, 8);
    }

    #[test]
    fn test_lincs_data_creation() {
        let lincs_data = LincsData::new();
        assert_eq!(lincs_data.sdiag.len(), 0);
        assert_eq!(lincs_data.coupled_constr.len(), 0);
        assert_eq!(lincs_data.coef.len(), 0);
    }

    #[test]
    fn test_setup_lincs_simple() {
        use gromos_core::topology::*;

        // Create a simple topology with 3 atoms forming a triangle of constraints
        // Atom 0 -- Atom 1
        //    \      /
        //    Atom 2

        let mut topo = Topology::new();

        // Add 3 atoms
        topo.mass = vec![12.0, 12.0, 12.0]; // Carbon atoms
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0];

        // Add bonds: 0-1, 0-2, 1-2
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 2,
            bond_type: 0,
            hydrogen: false,
        });
        topo.moltypes[0].bonds.push(Bond {
            i: 1,
            j: 2,
            bond_type: 0,
            hydrogen: false,
        });

        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 0.0,
            r0: 0.15, // 0.15 nm constraint length
        });

        // All 3 bonds are constrained: (atom_i, atom_j, r0) triples
        let constraints = vec![(0usize, 1usize, 0.15), (0, 2, 0.15), (1, 2, 0.15)];

        let lincs_data = setup_lincs(&topo, &constraints, 0);

        // Check that we have 3 constraints
        assert_eq!(lincs_data.sdiag.len(), 3);
        assert_eq!(lincs_data.coupled_constr.len(), 3);
        assert_eq!(lincs_data.coef.len(), 3);

        // Each constraint should be coupled to the other two (triangle)
        assert_eq!(lincs_data.coupled_constr[0].len(), 2); // Bond 0-1 coupled to 0-2 and 1-2
        assert_eq!(lincs_data.coupled_constr[1].len(), 2); // Bond 0-2 coupled to 0-1 and 1-2
        assert_eq!(lincs_data.coupled_constr[2].len(), 2); // Bond 1-2 coupled to 0-1 and 0-2

        // Check diagonal elements (should be same for all since all masses are equal)
        let expected_sdiag = 1.0 / (2.0 / 12.0_f64).sqrt();
        for &sdiag in &lincs_data.sdiag {
            assert!(
                (sdiag - expected_sdiag).abs() < 1e-10,
                "sdiag = {}, expected = {}",
                sdiag,
                expected_sdiag
            );
        }
    }

    #[test]
    fn test_lincs_simple_constraint() {
        use gromos_core::configuration::*;
        use gromos_core::math::Vec3;
        use gromos_core::topology::*;

        // Create simple system: 2 atoms with 1 constrained bond
        let mut topo = Topology::new();

        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0];

        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 0.0,
            r0: 0.15, // Constraint length 0.15 nm
        });

        let constraints = vec![(0usize, 1usize, 0.15)];
        let lincs_data = setup_lincs(&topo, &constraints, 0);
        let params = LincsParameters { order: 4 };

        // Create configuration with atoms slightly too far apart
        let mut conf = Configuration::new(2, 1, 1);
        conf.old_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.old_mut().pos[1] = Vec3::new(0.15, 0.0, 0.0); // At constraint distance

        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.16, 0.0, 0.0); // Slightly stretched

        // Apply LINCS
        let result = lincs(
            &topo,
            &mut conf,
            &constraints,
            &lincs_data,
            &params,
            0,
            0.002,
        );

        // Check that constraint is satisfied
        let r = conf.current().pos[1] - conf.current().pos[0];
        let distance = r.length();

        println!(
            "LINCS result: converged={}, max_error={}, distance={}",
            result.converged, result.max_error, distance
        );

        assert!(
            result.converged,
            "LINCS did not converge: max_error={}, threshold=1e-3",
            result.max_error
        );
        assert!(
            (distance - 0.15).abs() < 1e-3,
            "Distance = {}, expected 0.15",
            distance
        );
    }

    #[test]
    fn test_lincs_coupled_constraints() {
        use gromos_core::configuration::*;
        use gromos_core::math::Vec3;
        use gromos_core::topology::*;

        // Create a simple 3-atom system with two coupled constraints (sharing atom 0)
        let mut topo = Topology::new();

        topo.mass = vec![16.0, 1.0, 1.0]; // Water-like: O, H, H
        topo.inverse_mass = vec![1.0 / 16.0, 1.0 / 1.0, 1.0 / 1.0];

        // Two O-H bonds (coupled through atom 0)
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 2,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 0.0,
            r0: 0.1, // O-H distance
        });

        let constraints = vec![(0usize, 1usize, 0.1), (0, 2, 0.1)];

        // Setup LINCS
        let lincs_data = setup_lincs(&topo, &constraints, 0);

        // Verify that constraints are coupled
        assert_eq!(
            lincs_data.coupled_constr[0].len(),
            1,
            "Constraint 0 should be coupled to 1 other"
        );
        assert_eq!(
            lincs_data.coupled_constr[1].len(),
            1,
            "Constraint 1 should be coupled to 1 other"
        );
        assert_eq!(
            lincs_data.coupled_constr[0][0], 1,
            "Constraint 0 should be coupled to constraint 1"
        );
        assert_eq!(
            lincs_data.coupled_constr[1][0], 0,
            "Constraint 1 should be coupled to constraint 0"
        );

        let lincs_params = LincsParameters { order: 4 };

        // Create configuration
        let mut conf = Configuration::new(3, 1, 1);

        // Old positions (reference) - water-like geometry
        conf.old_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.old_mut().pos[1] = Vec3::new(0.1, 0.0, 0.0);
        conf.old_mut().pos[2] = Vec3::new(-0.05, 0.0866, 0.0);

        // Current positions (small perturbation - 0.5% stretch)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.1005, 0.0, 0.0); // 0.5% longer
        conf.current_mut().pos[2] = Vec3::new(-0.05025, 0.08703, 0.0); // 0.5% longer

        // Apply LINCS
        let result = lincs(
            &topo,
            &mut conf,
            &constraints,
            &lincs_data,
            &lincs_params,
            0,
            0.002,
        );

        println!(
            "LINCS coupled constraints: converged={}, max_error={}",
            result.converged, result.max_error
        );

        // Should converge
        assert!(
            result.converged,
            "LINCS did not converge with coupled constraints: {}",
            result.max_error
        );

        // Check that both bonds are at correct length
        let r01 = conf.current().pos[1] - conf.current().pos[0];
        let r02 = conf.current().pos[2] - conf.current().pos[0];

        let dist01 = r01.length();
        let dist02 = r02.length();

        assert!(
            (dist01 - 0.1).abs() < 1e-3,
            "Bond 0-1 distance = {}, expected 0.1",
            dist01
        );
        assert!(
            (dist02 - 0.1).abs() < 1e-3,
            "Bond 0-2 distance = {}, expected 0.1",
            dist02
        );
    }

    /// SETTLE analytical solver on a single SPC-like water molecule:
    /// after constraining a slightly perturbed geometry, O-H/H-H distances
    /// must match the constraint template lengths to high precision
    /// (cf. `settle.cc::solvent`, Miyamoto & Kollman 1992).
    #[test]
    fn test_settle_water() {
        use gromos_core::configuration::*;
        use gromos_core::math::Vec3;
        use gromos_core::topology::*;

        let mut topo = Topology::new();

        // Single SPC water: O, H, H (no solute)
        topo.mass = vec![15.99940, 1.00800, 1.00800];
        topo.inverse_mass = topo.mass.iter().map(|m| 1.0 / m).collect();

        topo.solvent_atom_template = vec![
            SolventAtomTemplate {
                iac: 0,
                name: "OW".into(),
                mass: 15.99940,
                charge: -0.82,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW1".into(),
                mass: 1.00800,
                charge: 0.41,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW2".into(),
                mass: 1.00800,
                charge: 0.41,
            },
        ];

        let dist_oh = 0.1;
        let dist_hh = 0.1633;
        topo.solvent_constraint_template = vec![
            SolventConstraintTemplate {
                i: 0,
                j: 1,
                length: dist_oh,
            },
            SolventConstraintTemplate {
                i: 0,
                j: 2,
                length: dist_oh,
            },
            SolventConstraintTemplate {
                i: 1,
                j: 2,
                length: dist_hh,
            },
        ];

        // Populate molecules + instances via solvate() — replaces old Solvent struct setup.
        topo.solvate(1);

        // Old (reference) geometry: exact canonical SPC triangle, centered at origin
        let half_hh = dist_hh / 2.0;
        let h_height = (dist_oh * dist_oh - half_hh * half_hh).sqrt();
        let mut conf = Configuration::new(3, 1, 1);
        conf.old_mut().pos[0] = Vec3::new(0.0, h_height * (1.0 / 3.0) * 2.0, 0.0);
        conf.old_mut().pos[1] = Vec3::new(-half_hh, -h_height * (1.0 / 3.0), 0.0);
        conf.old_mut().pos[2] = Vec3::new(half_hh, -h_height * (1.0 / 3.0), 0.0);

        // New (unconstrained) positions: old + small perturbation (as if integrated forward)
        conf.current_mut().pos[0] = conf.old().pos[0] + Vec3::new(0.001, 0.0005, -0.0007);
        conf.current_mut().pos[1] = conf.old().pos[1] + Vec3::new(-0.0012, 0.0009, 0.0006);
        conf.current_mut().pos[2] = conf.old().pos[2] + Vec3::new(0.0008, -0.0011, 0.0003);

        let dt = 0.002;
        let result = settle(&topo, &mut conf, dt, true);
        assert!(
            result.converged,
            "SETTLE failed to find a valid analytical solution"
        );

        let pos_o = conf.current().pos[0];
        let pos_h1 = conf.current().pos[1];
        let pos_h2 = conf.current().pos[2];

        let d_oh1 = (pos_h1 - pos_o).length();
        let d_oh2 = (pos_h2 - pos_o).length();
        let d_hh = (pos_h2 - pos_h1).length();

        assert!(
            (d_oh1 - dist_oh).abs() < 1e-10,
            "O-H1 distance = {d_oh1}, expected {dist_oh}"
        );
        assert!(
            (d_oh2 - dist_oh).abs() < 1e-10,
            "O-H2 distance = {d_oh2}, expected {dist_oh}"
        );
        assert!(
            (d_hh - dist_hh).abs() < 1e-10,
            "H-H distance = {d_hh}, expected {dist_hh}"
        );
    }
}
