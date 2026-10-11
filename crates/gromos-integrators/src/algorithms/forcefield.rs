//! Forcefield algorithm — wraps the complete force calculation.
//!
//! Equivalent to GROMOS `algorithm::Forcefield`:
//! 1. Update pairlist (if needed)
//! 2. Calculate bonded forces (bonds, angles, dihedrals)
//! 3. Calculate nonbonded forces (LJ + CRF)
//! 4. Assemble total forces and energies onto Configuration

use gromos_core::algorithm::{Algorithm, SimulationState};
use gromos_core::configuration::Configuration;
use gromos_core::math::{Mat3, Periodicity, Vec3};
use gromos_core::pairlist::{PairlistAlgorithm, PairlistContainer, CG_RADIUS_MARGIN};
use gromos_core::topology::Topology;

use gromos_core::lambdas::InteractionLambda;
use gromos_forces::bonded::{
    calculate_bonded_forces_ntf, calculate_perturbed_bonded_forces, BondedGroupEnergies, TermLambda,
};
use gromos_forces::nonbonded::{
    lj_crf_innerloop, lj_crf_innerloop_cg_grouped, lj_crf_innerloop_cg_grouped_novirial,
    lj_crf_innerloop_cg_grouped_parallel, lj_crf_innerloop_cg_grouped_parallel_novirial,
    lj_crf_innerloop_novirial, lj_crf_innerloop_parallel, lj_crf_innerloop_parallel_novirial,
    one_four_interaction_loop, perturbed_atom_pair_correction, perturbed_excluded_correction,
    perturbed_one_four_correction, perturbed_pairlist_correction, perturbed_self_energy_correction,
    rf_excluded_interactions, solvent_innerloop, solvent_innerloop_novirial,
    solvent_innerloop_parallel, solvent_innerloop_parallel_novirial, CGPairGroup, CRFParameters,
    EnergyGroups, ForceStorage, LJParamMatrix, LJParameters, PertAtomInfo, PertNBCorrection,
    PerturbedLambdaParams,
};
use gromos_forces::restraints::{
    DistanceRestraints, PerturbedDistanceRestraints, PositionRestraints,
};

use crate::algorithms::pressure_calculation::VirialType;
use rayon::prelude::*;

/// Accumulated wall time per force-calculation sub-phase.
///
/// The gromos-rs counterpart of the gromosXX `TIMING` block's `NonBonded`
/// sub-entries, so the two engines' profiles can be compared line by line.
/// Populated only when [`Forcefield::enable_timing`] is set.
#[derive(Debug, Clone, Copy, Default)]
pub struct ForcefieldTiming {
    /// Pairlist construction (amortised over NSNB steps).
    pub pairlist: std::time::Duration,
    /// Bonded terms: bonds, angles, impropers, dihedrals.
    pub bonded: std::time::Duration,
    /// Short-range solute-solute nonbonded innerloop.
    pub solute: std::time::Duration,
    /// Short-range solvent-solvent nonbonded innerloop.
    pub solvent: std::time::Duration,
    /// Rebuilding charge-group pair metadata for the grouped kernels.
    /// Has no gromosXX counterpart.
    pub cg_groups: std::time::Duration,
    /// Twin-range long-range interactions.
    pub longrange: std::time::Duration,
    /// Reaction-field self-energy and excluded-pair Coulomb.
    pub rf_excluded: std::time::Duration,
}

/// Forcefield algorithm that performs the complete force calculation each step.
///
/// Owns all the nonbonded interaction parameters and pairlist state.
/// After `apply()`, `conf.current()` has updated forces and potential energies.
pub struct Forcefield {
    /// Flat Lennard-Jones parameter matrix (cache-friendly contiguous layout)
    pub lj_params: LJParamMatrix,
    /// Coulomb reaction field parameters
    pub crf_params: CRFParameters,
    /// Coulomb prefactor 1/(4πε₀) in GROMOS units (kJ mol⁻¹ nm e⁻²).
    pub four_pi_eps_i: f64,
    /// MPI domain decomposition by pair: `(rank, size)` keeps only the pairs whose first atom
    /// index ≡ rank (mod size) after every pairlist update, so each rank evaluates its share
    /// and [`Self::set_nonbonded_reducer`] sums the shares. `None`: every pair, no reduction.
    pub pair_partition: Option<(usize, usize)>,
    /// Whether the pairlist currently held has been reduced to this rank's share (the initial
    /// pairlist is generated before the first `apply`).
    partition_applied: bool,
    /// Contributions every rank computes in full (excluded-pair and self reaction field, 1-4,
    /// the perturbed corrections that are not pairlist-based) — kept out of the reduction.
    local_storage: ForceStorage,
    /// Sums `nonbonded_storage` across ranks (MPI); called once per step after the pair terms.
    nonbonded_reducer: Option<Box<dyn FnMut(&mut ForceStorage) + Send>>,
    /// NSLFEXCL: add the reaction-field term of the excluded pairs and the self term
    /// (gromosXX `sim.param().nonbonded.rf_excluded`, which gates `RF_excluded_outerloop`,
    /// the perturbed variant, and the excluded-state term of perturbed atom pairs).
    pub rf_excluded: bool,
    /// Boundary conditions
    pub periodicity: Periodicity,
    /// Pairlist container
    pub pairlist: PairlistContainer,
    /// Pairlist construction algorithm
    pub pairlist_algorithm: PairlistAlgorithm,
    /// Whether to use quartic bond potentials (vs harmonic)
    pub use_quartic_bonds: bool,
    /// COVALENTFORM: which functional form each covalent term takes (gromosXX's NTBBH/NTBAH/NTBDN).
    pub covalent_form: gromos_forces::bonded::CovalentForm,
    /// Whether to run nonbonded in parallel
    pub parallel_nonbonded: bool,
    /// NTF flags: which bonded terms to include (GROMOS FORCE block)
    pub ntf_bond: bool,
    pub ntf_angle: bool,
    pub ntf_improper: bool,
    pub ntf_dihedral: bool,
    /// Number of atoms per solvent molecule (e.g. 3 for water)
    pub atoms_per_solvent: usize,
    /// Virial type: None, Atomic, or Molecular (controls KE tensor + virial correction)
    pub virial_type: VirialType,
    /// Reusable nonbonded force storage (avoids allocation per step)
    nonbonded_storage: ForceStorage,
    /// CG-pair group metadata for the solute shortrange pairlist
    solute_cg_groups: Vec<CGPairGroup>,
    /// CG-pair group metadata for the long-range solute pairlist.
    solute_long_cg_groups: Vec<CGPairGroup>,
    /// Reusable long-range force accumulator, so a twin-range update does not
    /// allocate a fresh `ForceStorage` every NSNB steps.
    lr_storage: ForceStorage,
    /// Whether `solute_cg_groups`/`solute_long_cg_groups` reflect the current pairlist.
    cg_groups_built: bool,
    /// Cached charge vector
    charges: Vec<f64>,
    /// Cached IAC vector
    iac_u32: Vec<u32>,
    /// Cached long-range forces (twin-range: reused between pairlist updates)
    longrange_forces: Vec<Vec3>,
    /// Cached long-range LJ energy
    longrange_e_lj: f64,
    /// Cached long-range CRF energy
    longrange_e_crf: f64,
    /// Cached long-range virial tensor (twin-range: reused between pairlist updates)
    /// Must be added to nonbonded_storage.virial each step, just like forces/energies.
    longrange_virial: [[f64; 3]; 3],
    /// Number of energy groups (`Topology::energy_groups`); above one, every term is also
    /// booked per group (pair) as gromosXX does.
    n_energy_groups: usize,
    /// Bonded energies per energy group, unperturbed and perturbed terms together.
    bonded_groups: BondedGroupEnergies,
    /// Perturbed bonded dH/dλ per energy group.
    bonded_dhdl_groups: BondedGroupEnergies,
    /// Scratch for a restraint energy per energy group.
    restraint_groups: Vec<f64>,
    /// Scratch for a restraint dH/dλ per energy group.
    restraint_dhdl_groups: Vec<f64>,
    /// Whether twin-range is active (RCUTP < RCUTL)
    twin_range_active: bool,
    /// Whether long-range forces have been computed at least once
    longrange_computed: bool,
    /// Whether CG-grouped solute innerloop is safe (box large enough relative to cutoff + CG size)
    use_cg_grouped_solute: bool,
    /// Accumulated sub-phase wall time; only written when `timing_enabled`.
    timing: ForcefieldTiming,
    /// Whether to accumulate `timing`. Off by default — each phase costs two
    /// `Instant::now()` reads per step, which is pure overhead in production.
    timing_enabled: bool,
    /// Position restraints (empty = none)
    pub position_restraints: PositionRestraints,
    /// Distance restraints (empty = none)
    pub distance_restraints: DistanceRestraints,
    /// Perturbed distance restraints (empty = none)
    pub perturbed_distance_restraints: PerturbedDistanceRestraints,
    /// The perturbed nonbonded λ parameters of every energy-group pair, rebuilt each step from
    /// the topology's individual λ (`Topology::individual_lambdas`, `lambda_exp`).
    lambda_table: Vec<PerturbedLambdaParams>,
    /// Global soft-core alpha for LJ (ALPHLJ from imd PERTURBATION block).
    /// Per-atom lj_soft from .ptp is multiplied by this to get effective alpha_lj,
    /// matching GROMOS in_perturbation.cc:1308 `lj_soft *= param.perturbation.soft_vdw`.
    pub global_alphlj: f64,
    /// Global soft-core alpha for CRF (ALPHC from imd PERTURBATION block).
    pub global_alphc: f64,
    /// Per-atom perturbation data; None for non-perturbed atoms.
    /// Populated when a pttopo is loaded (topo.perturbed_solute.atoms is non-empty).
    pert_info: Vec<Option<PertAtomInfo>>,
}

impl Forcefield {
    pub fn new(
        lj_params: Vec<Vec<LJParameters>>,
        crf_params: CRFParameters,
        periodicity: Periodicity,
        pairlist: PairlistContainer,
        pairlist_algorithm: PairlistAlgorithm,
    ) -> Self {
        let twin_range_active = pairlist.short_range_cutoff < pairlist.long_range_cutoff - 1e-10;
        Self {
            lj_params: LJParamMatrix::from_nested(&lj_params),
            crf_params,
            four_pi_eps_i: gromos_core::units::four_pi_eps_i,
            pair_partition: None,
            partition_applied: false,
            local_storage: ForceStorage::new(0),
            nonbonded_reducer: None,
            rf_excluded: true,
            periodicity,
            pairlist,
            pairlist_algorithm,
            use_quartic_bonds: true,
            covalent_form: gromos_forces::bonded::CovalentForm::default(),
            parallel_nonbonded: false,
            ntf_bond: true,
            ntf_angle: true,
            ntf_improper: true,
            ntf_dihedral: true,
            atoms_per_solvent: 3,
            virial_type: VirialType::None,
            nonbonded_storage: ForceStorage::new(0),
            solute_cg_groups: Vec::new(),
            solute_long_cg_groups: Vec::new(),
            lr_storage: ForceStorage::new(0),
            cg_groups_built: false,
            charges: Vec::new(),
            iac_u32: Vec::new(),
            longrange_forces: Vec::new(),
            longrange_e_lj: 0.0,
            longrange_e_crf: 0.0,
            longrange_virial: [[0.0; 3]; 3],
            n_energy_groups: 1,
            bonded_groups: BondedGroupEnergies::default(),
            bonded_dhdl_groups: BondedGroupEnergies::default(),
            restraint_groups: Vec::new(),
            restraint_dhdl_groups: Vec::new(),
            twin_range_active,
            longrange_computed: false,
            use_cg_grouped_solute: false, // computed at init when box size is known
            timing: ForcefieldTiming::default(),
            timing_enabled: false,
            position_restraints: PositionRestraints::new(),
            distance_restraints: DistanceRestraints::new(),
            perturbed_distance_restraints: PerturbedDistanceRestraints::new(),
            lambda_table: Vec::new(),
            global_alphlj: 1.0,
            global_alphc: 1.0,
            pert_info: Vec::new(),
        }
    }

    /// Enable accumulation of the per-sub-phase breakdown in [`Self::timing`].
    ///
    /// Off by default: timing costs two `Instant::now()` reads per phase per step.
    pub fn enable_timing(&mut self, enabled: bool) {
        self.timing_enabled = enabled;
    }

    /// Accumulated wall time per force sub-phase.
    ///
    /// All-zero unless [`Self::enable_timing`] was called.
    pub fn timing(&self) -> ForcefieldTiming {
        self.timing
    }

    /// Build per-atom perturbation info from the topology's perturbed_solute.
    ///
    /// GROMOS multiplies the per-atom lj_soft/crf_soft from the .ptp file by the
    /// global ALPHLJ/ALPHC from the .imd PERTURBATION block before storing them
    /// (in_perturbation.cc:1308: `lj_soft *= param.perturbation.soft_vdw`).
    /// We replicate that here: effective alpha_lj = pa.lj_soft * global_alphlj.
    fn build_pert_info(
        topo: &Topology,
        global_alphlj: f64,
        global_alphc: f64,
    ) -> Vec<Option<PertAtomInfo>> {
        let n = topo.inverse_mass.len();
        let mut info: Vec<Option<PertAtomInfo>> = vec![None; n];
        for pa in &topo.perturbed_solute.atoms {
            if pa.seq < n {
                // topo.iac and PerturbedAtom.a_iac / b_iac are both 0-indexed:
                //   topo.iac   → topology.rs reader subtracts 1 from file's 1-indexed value
                //   pa.a_iac   → ptp.rs reader subtracts 1 from pttopo's 1-indexed value
                // No conversion needed here.
                info[pa.seq] = Some(PertAtomInfo {
                    a_iac: pa.a_iac,
                    b_iac: pa.b_iac,
                    a_charge: pa.a_charge,
                    b_charge: pa.b_charge,
                    alpha_lj: pa.lj_soft * global_alphlj,
                    alpha_crf: pa.crf_soft * global_alphc,
                });
            }
        }
        info
    }

    /// Convert topology LJ parameters to the nonbonded-crate format.
    pub fn convert_lj_parameters(topo: &Topology) -> Vec<Vec<LJParameters>> {
        topo.lj_parameters
            .iter()
            .map(|row| row.iter().map(Into::into).collect())
            .collect()
    }
}

impl Forcefield {
    /// Allocate the per-topology buffers (atoms, energy groups); called once, from `init`.
    fn size_for(&mut self, topo: &Topology) {
        let n = topo.inverse_mass.len();
        let n_groups = topo.energy_groups.len().max(1);
        if n_groups > 1 {
            assert_eq!(
                topo.atom_to_energy_group.len(),
                n,
                "every atom needs an energy group"
            );
        }
        self.n_energy_groups = n_groups;
        self.nonbonded_storage = ForceStorage::with_groups(n, n_groups);
        self.local_storage = ForceStorage::with_groups(n, n_groups);
        self.lr_storage = ForceStorage::with_groups(0, n_groups);
        self.bonded_groups = BondedGroupEnergies::new(n_groups);
        self.bonded_dhdl_groups = BondedGroupEnergies::new(n_groups);
        self.restraint_groups = vec![0.0; n_groups];
        self.restraint_dhdl_groups = vec![0.0; n_groups];
        self.charges = topo.charge.clone();
        self.iac_u32 = topo.iac.iter().map(|&i| i as u32).collect();
        self.longrange_forces = vec![Vec3::ZERO; n];
        self.pert_info = Self::build_pert_info(topo, self.global_alphlj, self.global_alphc);
    }

    /// Install the cross-rank sum of the pair terms (MPI). See [`Self::pair_partition`].
    pub fn set_nonbonded_reducer(&mut self, f: Box<dyn FnMut(&mut ForceStorage) + Send>) {
        self.nonbonded_reducer = Some(f);
    }
}

impl Algorithm for Forcefield {
    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    fn init(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        self.size_for(topo);
        // gromosXX `Distance_Restraint_Interaction::init` and its perturbed twin: the r⁻³
        // averages start from the start-up configuration
        self.distance_restraints
            .init_averages(conf, &self.periodicity);
        self.perturbed_distance_restraints
            .init_averages(conf, &self.periodicity);
        Ok(())
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        let n_atoms = topo.inverse_mass.len();

        // Lazy init if init() was not called
        if self.nonbonded_storage.forces.len() != n_atoms {
            self.size_for(topo);
        }
        let egroups = if self.n_energy_groups > 1 {
            EnergyGroups::new(&topo.atom_to_energy_group, self.n_energy_groups)
        } else {
            EnergyGroups::single()
        };

        // IMPORTANT: Update periodicity from current box each step.
        // Under NPT, the barostat scales the box at the end of the previous step,
        // so self.periodicity (which caches box_size, half_box, inv_box) must be
        // refreshed before any force/pairlist/virial calculation uses it.
        // Without this, all PBC nearest_image calls use stale box dimensions.
        if !matches!(self.periodicity, Periodicity::Vacuum(_)) {
            use gromos_core::configuration::BoxType;
            use gromos_core::math::Rectangular;
            let box_cfg = &conf.current().box_config;
            match box_cfg.box_type {
                BoxType::Rectangular => {
                    let dims = box_cfg.dimensions();
                    if dims.x > 0.0 && dims.y > 0.0 && dims.z > 0.0 {
                        self.periodicity = Periodicity::Rectangular(Rectangular::new(dims));
                    }
                },
                BoxType::Triclinic | BoxType::TruncatedOctahedral => {
                    use gromos_core::math::Triclinic;
                    self.periodicity = Periodicity::Triclinic(Triclinic::new(box_cfg.vectors));
                },
                _ => {},
            }
        }

        // Determine if CG-grouped solute innerloop is safe for this box.
        // Safety condition: half_box_min > cutoff_long + max_cg_diameter
        // When this fails, atom pairs near the half-box boundary can get
        // different periodic image rounding than their CG reference atoms.
        let cg_data_available =
            topo.atom_to_chargegroup.len() == topo.num_atoms() && !topo.chargegroups.is_empty();
        self.use_cg_grouped_solute = cg_data_available
            && match &self.periodicity {
                Periodicity::Vacuum(_) => true, // no PBC → shift is always 0, safe
                Periodicity::Rectangular(rect) => {
                    let half_min = rect.half_box.x.min(rect.half_box.y).min(rect.half_box.z);
                    let cutoff = self.pairlist.long_range_cutoff + self.pairlist.skin;
                    half_min > cutoff + CG_RADIUS_MARGIN
                },
                Periodicity::Triclinic(_) => false, // be conservative for triclinic
            };

        // --- prepare_virial: compute kinetic energy tensor (GROMOS: before force calc) ---
        // Uses current velocities and positions; stored in current().kinetic_energy_tensor
        if self.virial_type != VirialType::None {
            let ke_tensor =
                if self.virial_type == VirialType::Molecular && !topo.pressure_groups.is_empty() {
                    compute_molecular_ke_tensor(topo, conf, &self.periodicity)
                } else {
                    compute_atomic_ke_tensor(topo, conf, n_atoms)
                };
            conf.current_mut().kinetic_energy_tensor = ke_tensor;
        }

        // --- 1. Update pairlist if needed ---
        // The initial pairlist was generated before the first `apply`: reduce it to this
        // rank's share once (later updates are reduced right after the update below).
        if let (Some((rank, size)), false) = (self.pair_partition, self.partition_applied) {
            for list in [
                &mut self.pairlist.solute_short,
                &mut self.pairlist.solute_long,
                &mut self.pairlist.solvent_short,
                &mut self.pairlist.solvent_long,
            ] {
                gromos_core::pairlist::keep_partition(list, rank, size);
            }
            self.partition_applied = true;
            self.cg_groups_built = false;
        }
        let pairlist_updated = self.pairlist.needs_update();
        if pairlist_updated {
            let t_pl = self.timing_enabled.then(std::time::Instant::now);
            self.pairlist_algorithm
                .update(topo, conf, &mut self.pairlist, &self.periodicity);
            if let Some((rank, size)) = self.pair_partition {
                for list in [
                    &mut self.pairlist.solute_short,
                    &mut self.pairlist.solute_long,
                    &mut self.pairlist.solvent_short,
                    &mut self.pairlist.solvent_long,
                ] {
                    gromos_core::pairlist::keep_partition(list, rank, size);
                }
            }
            if let Some(t) = t_pl {
                self.timing.pairlist += t.elapsed();
            }
        }
        self.pairlist.step();

        // Charge-group pair metadata for the grouped kernels, rebuilt whenever the
        // pairlist changes (and once on the first call, since the initial pairlist
        // is generated before `apply`). The kernels read the pairlist's own `u32`
        // lists directly; there is no per-step copy.
        let t_groups = self.timing_enabled.then(std::time::Instant::now);
        if pairlist_updated || !self.cg_groups_built {
            self.solute_cg_groups.clear();
            self.solute_long_cg_groups.clear();
            if self.use_cg_grouped_solute {
                build_cg_groups(
                    topo,
                    &self.pairlist.solute_short,
                    &mut self.solute_cg_groups,
                );
                // Same grouping for the long-range list; the safety condition is
                // keyed on the long-range cutoff, so it covers these pairs too.
                build_cg_groups(
                    topo,
                    &self.pairlist.solute_long,
                    &mut self.solute_long_cg_groups,
                );
            }
            self.cg_groups_built = true;
        }
        if let Some(t) = t_groups {
            self.timing.cg_groups += t.elapsed();
        }

        // --- 2. Calculate bonded forces (conditional on NTF flags) ---
        let any_bonded = self.ntf_bond || self.ntf_angle || self.ntf_improper || self.ntf_dihedral;
        let t_bonded = self.timing_enabled.then(std::time::Instant::now);
        self.bonded_groups.clear();
        self.bonded_dhdl_groups.clear();
        let bonded_result = if any_bonded {
            calculate_bonded_forces_ntf(
                topo,
                conf,
                &self.periodicity,
                self.use_quartic_bonds,
                self.covalent_form,
                self.ntf_bond,
                self.ntf_angle,
                self.ntf_dihedral,
                self.ntf_improper,
                &mut self.bonded_groups,
            )
        } else {
            gromos_forces::bonded::ForceEnergy::new(n_atoms)
        };
        if let Some(t) = t_bonded {
            self.timing.bonded += t.elapsed();
        }

        // --- 2b. Perturbed bonded forces (FEP/TI, NTG != 0) ---
        // Forces are accumulated in the returned ForceEnergyLambda and applied in step 4.
        let perturbed_bonded = if !topo.perturbed_solute.is_empty() {
            let pr = calculate_perturbed_bonded_forces(
                topo,
                conf,
                &self.periodicity,
                &mut self.bonded_groups,
                &mut self.bonded_dhdl_groups,
            );
            log::debug!(
                "  Perturbed bonded: E={:.10e} dH/dλ={:.10e}",
                pr.energy,
                pr.lambda_derivative
            );
            Some(pr)
        } else {
            None
        };

        // --- 3. Calculate nonbonded forces ---
        self.nonbonded_storage.clear();
        self.local_storage.clear();
        let need_virial = self.virial_type != VirialType::None;

        // Short-range solute interactions
        let t_solute = self.timing_enabled.then(std::time::Instant::now);
        if self.use_cg_grouped_solute {
            // CG-grouped: compute nearest_image once per CG pair (fast path for large boxes)
            if self.parallel_nonbonded {
                let result = if need_virial {
                    lj_crf_innerloop_cg_grouped_parallel(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist.solute_short,
                        &self.solute_cg_groups,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.four_pi_eps_i,
                        n_atoms,
                        egroups,
                    )
                } else {
                    lj_crf_innerloop_cg_grouped_parallel_novirial(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist.solute_short,
                        &self.solute_cg_groups,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.four_pi_eps_i,
                        n_atoms,
                        egroups,
                    )
                };
                self.nonbonded_storage.merge(&result);
            } else if need_virial {
                lj_crf_innerloop_cg_grouped(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solute_short,
                    &self.solute_cg_groups,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            } else {
                lj_crf_innerloop_cg_grouped_novirial(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solute_short,
                    &self.solute_cg_groups,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            }
        } else {
            // Per-atom nearest_image fallback (exact, for small boxes)
            if self.parallel_nonbonded {
                let result = if need_virial {
                    lj_crf_innerloop_parallel(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist.solute_short,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.four_pi_eps_i,
                        n_atoms,
                        egroups,
                    )
                } else {
                    lj_crf_innerloop_parallel_novirial(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist.solute_short,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.four_pi_eps_i,
                        n_atoms,
                        egroups,
                    )
                };
                self.nonbonded_storage.merge(&result);
            } else if need_virial {
                lj_crf_innerloop(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solute_short,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            } else {
                lj_crf_innerloop_novirial(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solute_short,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            }
        }
        if let Some(t) = t_solute {
            self.timing.solute += t.elapsed();
        }
        let e_lj_after_solute = self.nonbonded_storage.e_lj;
        let e_crf_after_solute = self.nonbonded_storage.e_crf;
        log::debug!(
            "  After solute innerloop: e_lj={:.10e}, e_crf={:.10e}",
            e_lj_after_solute,
            e_crf_after_solute
        );

        // Short-range solvent-solvent interactions (shared PBC shift)
        let t_solvent = self.timing_enabled.then(std::time::Instant::now);
        if !self.pairlist.solvent_short.is_empty() {
            if self.parallel_nonbonded {
                let result = if need_virial {
                    solvent_innerloop_parallel(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist.solvent_short,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.atoms_per_solvent,
                        self.four_pi_eps_i,
                        n_atoms,
                        egroups,
                    )
                } else {
                    solvent_innerloop_parallel_novirial(
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.pairlist.solvent_short,
                        &self.lj_params,
                        &self.crf_params,
                        &self.periodicity,
                        self.atoms_per_solvent,
                        self.four_pi_eps_i,
                        n_atoms,
                        egroups,
                    )
                };
                self.nonbonded_storage.merge(&result);
            } else if need_virial {
                solvent_innerloop(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solvent_short,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.atoms_per_solvent,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            } else {
                solvent_innerloop_novirial(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solvent_short,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.atoms_per_solvent,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            }
        }
        if let Some(t) = t_solvent {
            self.timing.solvent += t.elapsed();
        }
        let e_lj_after_solvent = self.nonbonded_storage.e_lj;
        let e_crf_after_solvent = self.nonbonded_storage.e_crf;
        log::debug!("  After solvent innerloop: e_lj={:.10e}, e_crf={:.10e} (delta_lj={:.10e}, delta_crf={:.10e})",
            e_lj_after_solvent, e_crf_after_solvent,
            e_lj_after_solvent - e_lj_after_solute, e_crf_after_solvent - e_crf_after_solute);

        // Twin-range long-range forces: recalculate on pairlist update, reuse otherwise
        let t_lr = self.timing_enabled.then(std::time::Instant::now);
        if self.twin_range_active {
            // Recalculate long-range if pairlist was updated OR first step
            if pairlist_updated || !self.longrange_computed {
                // Recalculate long-range forces and cache them. Virial variants only when
                // pressure coupling will read the tensor — the same rule the
                // short-range path applies; the virial loop is nine extra
                // multiply-adds per pair..
                // Buffer is reused across updates rather than reallocated.
                if self.lr_storage.forces.len() != n_atoms {
                    self.lr_storage = ForceStorage::with_groups(n_atoms, self.n_energy_groups);
                } else {
                    self.lr_storage.clear();
                }
                // Detached so the innerloop calls below can borrow `self`'s
                // pairlists and parameters while writing into it.
                let mut lr_storage = std::mem::replace(&mut self.lr_storage, ForceStorage::new(0));

                // Long-range solute interactions: the same kernel selection as the
                // short-range block (grouped vs flat, parallel vs serial, virial or
                // not). The long-range block ran serial kernels until 2026-08-29,
                // which left this phase flat when threads were added.
                if !self.pairlist.solute_long.is_empty() {
                    let pos = &conf.current().pos;
                    if self.use_cg_grouped_solute {
                        if self.parallel_nonbonded {
                            let partial = if need_virial {
                                lj_crf_innerloop_cg_grouped_parallel(
                                    pos,
                                    &self.charges,
                                    &self.iac_u32,
                                    &self.pairlist.solute_long,
                                    &self.solute_long_cg_groups,
                                    &self.lj_params,
                                    &self.crf_params,
                                    &self.periodicity,
                                    self.four_pi_eps_i,
                                    n_atoms,
                                    egroups,
                                )
                            } else {
                                lj_crf_innerloop_cg_grouped_parallel_novirial(
                                    pos,
                                    &self.charges,
                                    &self.iac_u32,
                                    &self.pairlist.solute_long,
                                    &self.solute_long_cg_groups,
                                    &self.lj_params,
                                    &self.crf_params,
                                    &self.periodicity,
                                    self.four_pi_eps_i,
                                    n_atoms,
                                    egroups,
                                )
                            };
                            lr_storage.merge(&partial);
                        } else if need_virial {
                            lj_crf_innerloop_cg_grouped(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solute_long,
                                &self.solute_long_cg_groups,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.four_pi_eps_i,
                                &mut lr_storage,
                                egroups,
                            );
                        } else {
                            lj_crf_innerloop_cg_grouped_novirial(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solute_long,
                                &self.solute_long_cg_groups,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.four_pi_eps_i,
                                &mut lr_storage,
                                egroups,
                            );
                        }
                    } else if self.parallel_nonbonded {
                        let partial = if need_virial {
                            lj_crf_innerloop_parallel(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solute_long,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.four_pi_eps_i,
                                n_atoms,
                                egroups,
                            )
                        } else {
                            lj_crf_innerloop_parallel_novirial(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solute_long,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.four_pi_eps_i,
                                n_atoms,
                                egroups,
                            )
                        };
                        lr_storage.merge(&partial);
                    } else if need_virial {
                        lj_crf_innerloop(
                            pos,
                            &self.charges,
                            &self.iac_u32,
                            &self.pairlist.solute_long,
                            &self.lj_params,
                            &self.crf_params,
                            &self.periodicity,
                            self.four_pi_eps_i,
                            &mut lr_storage,
                            egroups,
                        );
                    } else {
                        lj_crf_innerloop_novirial(
                            pos,
                            &self.charges,
                            &self.iac_u32,
                            &self.pairlist.solute_long,
                            &self.lj_params,
                            &self.crf_params,
                            &self.periodicity,
                            self.four_pi_eps_i,
                            &mut lr_storage,
                            egroups,
                        );
                    }
                }
                let lr_solute_e_lj = lr_storage.e_lj;
                let lr_solute_e_crf = lr_storage.e_crf;

                // Long-range solvent-solvent interactions. The pairlist says which
                // representation it stored: sentinel molecule pairs (one shared
                // image, solvent kernel with parameter tables) when that is provably
                // exact for this box, else expanded atom pairs with a per-atom
                // minimum image. Parallel and virial selection as for short-range.
                if !self.pairlist.solvent_long.is_empty() {
                    let pos = &conf.current().pos;
                    if self.pairlist.solvent_long_is_sentinel {
                        if self.parallel_nonbonded {
                            let partial = if need_virial {
                                solvent_innerloop_parallel(
                                    pos,
                                    &self.charges,
                                    &self.iac_u32,
                                    &self.pairlist.solvent_long,
                                    &self.lj_params,
                                    &self.crf_params,
                                    &self.periodicity,
                                    self.atoms_per_solvent,
                                    self.four_pi_eps_i,
                                    n_atoms,
                                    egroups,
                                )
                            } else {
                                solvent_innerloop_parallel_novirial(
                                    pos,
                                    &self.charges,
                                    &self.iac_u32,
                                    &self.pairlist.solvent_long,
                                    &self.lj_params,
                                    &self.crf_params,
                                    &self.periodicity,
                                    self.atoms_per_solvent,
                                    self.four_pi_eps_i,
                                    n_atoms,
                                    egroups,
                                )
                            };
                            lr_storage.merge(&partial);
                        } else if need_virial {
                            solvent_innerloop(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solvent_long,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.atoms_per_solvent,
                                self.four_pi_eps_i,
                                &mut lr_storage,
                                egroups,
                            );
                        } else {
                            solvent_innerloop_novirial(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solvent_long,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.atoms_per_solvent,
                                self.four_pi_eps_i,
                                &mut lr_storage,
                                egroups,
                            );
                        }
                    } else if self.parallel_nonbonded {
                        let partial = if need_virial {
                            lj_crf_innerloop_parallel(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solvent_long,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.four_pi_eps_i,
                                n_atoms,
                                egroups,
                            )
                        } else {
                            lj_crf_innerloop_parallel_novirial(
                                pos,
                                &self.charges,
                                &self.iac_u32,
                                &self.pairlist.solvent_long,
                                &self.lj_params,
                                &self.crf_params,
                                &self.periodicity,
                                self.four_pi_eps_i,
                                n_atoms,
                                egroups,
                            )
                        };
                        lr_storage.merge(&partial);
                    } else if need_virial {
                        lj_crf_innerloop(
                            pos,
                            &self.charges,
                            &self.iac_u32,
                            &self.pairlist.solvent_long,
                            &self.lj_params,
                            &self.crf_params,
                            &self.periodicity,
                            self.four_pi_eps_i,
                            &mut lr_storage,
                            egroups,
                        );
                    } else {
                        lj_crf_innerloop_novirial(
                            pos,
                            &self.charges,
                            &self.iac_u32,
                            &self.pairlist.solvent_long,
                            &self.lj_params,
                            &self.crf_params,
                            &self.periodicity,
                            self.four_pi_eps_i,
                            &mut lr_storage,
                            egroups,
                        );
                    }
                }
                log::debug!(
                    "  LR solute: e_lj={:.10e}, e_crf={:.10e} ({} pairs)",
                    lr_solute_e_lj,
                    lr_solute_e_crf,
                    self.pairlist.solute_long.len()
                );
                log::debug!(
                    "  LR solvent: e_lj={:.10e}, e_crf={:.10e} ({} pairs)",
                    lr_storage.e_lj - lr_solute_e_lj,
                    lr_storage.e_crf - lr_solute_e_crf,
                    self.pairlist.solvent_long.len()
                );

                // Cache the long-range results
                self.longrange_forces.clear();
                self.longrange_forces.extend_from_slice(&lr_storage.forces);
                self.longrange_e_lj = lr_storage.e_lj;
                self.longrange_e_crf = lr_storage.e_crf;
                self.longrange_virial = lr_storage.virial;
                self.longrange_computed = true;
                self.lr_storage = lr_storage;
                log::debug!("  Long-range (recomputed): lr_e_lj={:.10e}, lr_e_crf={:.10e}, solute_long={}, solvent_long={}",
                    self.longrange_e_lj, self.longrange_e_crf,
                    self.pairlist.solute_long.len(), self.pairlist.solvent_long.len());
            }

            // Add cached long-range forces, energies, and virial to current step.
            // IMPORTANT: The virial from long-range pairs MUST be included here.
            // Without it, the virial trace is ~1000 kJ/mol too small, causing
            // wildly wrong pressures and barostat box-collapse under NPT.
            for i in 0..n_atoms {
                self.nonbonded_storage.forces[i] += self.longrange_forces[i];
            }
            self.nonbonded_storage.e_lj += self.longrange_e_lj;
            self.nonbonded_storage.e_crf += self.longrange_e_crf;
            // per-group long-range energies: `lr_storage` holds them until the next update
            self.nonbonded_storage.groups.merge(&self.lr_storage.groups);
            for a in 0..3 {
                for b in 0..3 {
                    self.nonbonded_storage.virial[a][b] += self.longrange_virial[a][b];
                }
            }
        } else {
            // No twin-range: evaluate long-range pairs as short-range (same as before)
            if !self.pairlist.solute_long.is_empty() {
                lj_crf_innerloop(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solute_long,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            }
            if !self.pairlist.solvent_long.is_empty() {
                solvent_innerloop(
                    &conf.current().pos,
                    &self.charges,
                    &self.iac_u32,
                    &self.pairlist.solvent_long,
                    &self.lj_params,
                    &self.crf_params,
                    &self.periodicity,
                    self.atoms_per_solvent,
                    self.four_pi_eps_i,
                    &mut self.nonbonded_storage,
                    egroups,
                );
            }
        }

        if let Some(t) = t_lr {
            self.timing.longrange += t.elapsed();
        }

        // RF self-energy and excluded-pair Coulomb (energy + forces) — NSLFEXCL
        if self.rf_excluded {
            let t_rf = self.timing_enabled.then(std::time::Instant::now);
            rf_excluded_interactions(
                &self.charges,
                &topo.exclusions,
                &conf.current().pos,
                &self.crf_params,
                &self.periodicity,
                self.four_pi_eps_i,
                &mut self.local_storage,
                topo.num_solute_atoms(),
                egroups,
            );
            if let Some(t) = t_rf {
                self.timing.rf_excluded += t.elapsed();
            }
        }

        // 1-4 interactions: LJ with cs6/cs12 + CRF with coulomb scaling
        // GROMOS: computed separately from pairlist, no cutoff
        if !topo.one_four_pairs.is_empty() {
            one_four_interaction_loop(
                &topo.one_four_pairs,
                &conf.current().pos,
                &self.charges,
                &self.iac_u32,
                &self.lj_params,
                &self.crf_params,
                &self.periodicity,
                self.four_pi_eps_i,
                &mut self.local_storage,
                1.0, // coulomb_scaling: 1.0 for standard GROMOS
                egroups,
            );
        }

        // --- 3b. Perturbed nonbonded corrections (FEP/TI) ---
        // The four steps below mirror GROMOS's perturbed_nonbonded_outerloop:
        //   perturbed_pairlist_correction   → replaces state-A contribution for perturbed pairs
        //   perturbed_self_energy_correction → corrects RF self for perturbed atoms
        //   perturbed_excluded_correction    → corrects RF excluded pairs for perturbed atoms
        //   perturbed_one_four_correction    → corrects 1-4 for perturbed atoms
        // A `.ptp` may carry only PERTATOMPAIR entries (no perturbed atom): the atom-pair
        // correction must still run (`pert_info` is then all `None`, which every other
        // correction treats as "nothing to do").
        let has_perturbed_atoms =
            !self.pert_info.is_empty() && self.pert_info.iter().any(|x| x.is_some());
        let has_atom_pairs = !topo.perturbed_solute.atom_pairs.is_empty();
        let (pert_nb_dhdl, pert_nb_dhdl_lj, pert_nb_dhdl_crf) =
            if has_perturbed_atoms || has_atom_pairs {
                if topo.individual_lambdas.n_groups() != self.n_energy_groups {
                    return Err(format!(
                        "the individual-λ table covers {} energy groups, the topology has {}: \
                         call Topology::update_for_lambda after assigning the energy groups",
                        topo.individual_lambdas.n_groups(),
                        self.n_energy_groups
                    ));
                }
                PerturbedLambdaParams::fill_table(
                    &topo.individual_lambdas,
                    topo.lambda_exp,
                    &mut self.lambda_table,
                );
                let lp = &self.lambda_table;
                // 3b-1: pairlist correction (delta over what state-A innerloop already computed)
                // Applied to BOTH short-range and long-range pairlists:
                // the regular innerloops computed state-A for all pairs; here we replace the
                // contribution for perturbed pairs with the full soft-core dual-topology result.
                let corr_pl = perturbed_pairlist_correction(
                    &conf.current().pos,
                    &self.pairlist.solute_short,
                    &self.pert_info,
                    &self.charges,
                    &self.iac_u32,
                    &self.lj_params,
                    &self.crf_params,
                    lp,
                    &self.periodicity,
                    self.four_pi_eps_i,
                    egroups,
                    &mut self.nonbonded_storage.groups,
                );
                self.nonbonded_storage.e_lj += corr_pl.delta_e_lj;
                self.nonbonded_storage.e_crf += corr_pl.delta_e_crf;
                // pairlist-based dH/dλ travels with the (reduced) pair storage; the locals
                // below collect the terms every rank computes in full
                self.nonbonded_storage.dhdl_lj += corr_pl.dhdl_lj;
                self.nonbonded_storage.dhdl_crf += corr_pl.dhdl_crf;
                let mut dhdl_nb = 0.0;
                let mut dhdl_nb_lj = 0.0;
                let mut dhdl_nb_crf = 0.0;
                for i in 0..n_atoms {
                    self.nonbonded_storage.forces[i] += corr_pl.forces[i];
                }

                // 3b-1b: same correction for ALL remaining pairlists that may contain
                // perturbed pairs. Solute-solvent long-range pairs are in solute_long.
                let pls: &[(&str, &[(u32, u32)])] = &[
                    ("solute_long", &self.pairlist.solute_long),
                    ("solvent_short", &self.pairlist.solvent_short),
                    ("solvent_long", &self.pairlist.solvent_long),
                ];
                for (_name, pl) in pls {
                    if pl.is_empty() {
                        continue;
                    }
                    let corr = perturbed_pairlist_correction(
                        &conf.current().pos,
                        pl,
                        &self.pert_info,
                        &self.charges,
                        &self.iac_u32,
                        &self.lj_params,
                        &self.crf_params,
                        lp,
                        &self.periodicity,
                        self.four_pi_eps_i,
                        egroups,
                        &mut self.nonbonded_storage.groups,
                    );
                    self.nonbonded_storage.e_lj += corr.delta_e_lj;
                    self.nonbonded_storage.e_crf += corr.delta_e_crf;
                    self.nonbonded_storage.dhdl_lj += corr.dhdl_lj;
                    self.nonbonded_storage.dhdl_crf += corr.dhdl_crf;
                    for i in 0..n_atoms {
                        self.nonbonded_storage.forces[i] += corr.forces[i];
                    }
                }
                log::debug!(
                    "  Pert corr solute_short: ({} pairs, Δe_lj={:.6e}, Δe_crf={:.6e})",
                    self.pairlist.solute_short.len(),
                    corr_pl.delta_e_lj,
                    corr_pl.delta_e_crf
                );
                log::debug!(
                    "  Perturbed NB: Δe_lj_pl={:.6e} Δe_crf_pl={:.6e} Δe_self=TBD dH/dλ_nb={:.6e}",
                    corr_pl.delta_e_lj,
                    corr_pl.delta_e_crf,
                    dhdl_nb
                );

                // 3b-2: RF self-energy correction for perturbed atoms (NSLFEXCL, like the
                // unperturbed term: gromosXX `perturbed_nonbonded_set.cc` gates the perturbed
                // RF-excluded outerloop, self term included, on `rf_excluded`)
                let mut de_self = 0.0;
                if self.rf_excluded {
                    let (e, dhdl_self) = perturbed_self_energy_correction(
                        &self.pert_info,
                        &self.crf_params,
                        lp,
                        self.four_pi_eps_i,
                        egroups,
                        &mut self.local_storage.groups,
                    );
                    de_self = e;
                    self.local_storage.e_crf += de_self;
                    dhdl_nb += dhdl_self;
                    dhdl_nb_crf += dhdl_self;
                }

                // 3b-3: excluded-pair CRF correction (separate accumulator to avoid double-add)
                if self.rf_excluded {
                    use gromos_forces::nonbonded::PertNBCorrection;
                    let mut corr_ex = PertNBCorrection::new(n_atoms);
                    perturbed_excluded_correction(
                        &topo.exclusions,
                        &conf.current().pos,
                        &self.pert_info,
                        &self.charges,
                        &self.crf_params,
                        lp,
                        &self.periodicity,
                        topo.num_solute_atoms(),
                        self.four_pi_eps_i,
                        &mut corr_ex,
                        egroups,
                        &mut self.local_storage.groups,
                    );
                    self.local_storage.e_crf += corr_ex.delta_e_crf;
                    dhdl_nb += corr_ex.dhdl;
                    dhdl_nb_crf += corr_ex.dhdl_crf;
                    for i in 0..n_atoms {
                        self.local_storage.forces[i] += corr_ex.forces[i];
                    }
                }

                // 3b-4: 1-4 correction (separate accumulator)
                if !topo.one_four_pairs.is_empty() {
                    use gromos_forces::nonbonded::PertNBCorrection;
                    let mut corr14 = PertNBCorrection::new(n_atoms);
                    perturbed_one_four_correction(
                        &topo.one_four_pairs,
                        &conf.current().pos,
                        &self.pert_info,
                        &self.charges,
                        &self.iac_u32,
                        &self.lj_params,
                        &self.crf_params,
                        lp,
                        &self.periodicity,
                        self.four_pi_eps_i,
                        &mut corr14,
                        egroups,
                        &mut self.local_storage.groups,
                    );
                    log::debug!(
                        "  1-4 corr: Δe_lj={:.6e} Δe_crf={:.6e}",
                        corr14.delta_e_lj,
                        corr14.delta_e_crf
                    );
                    self.local_storage.e_lj += corr14.delta_e_lj;
                    self.local_storage.e_crf += corr14.delta_e_crf;
                    dhdl_nb += corr14.dhdl;
                    dhdl_nb_lj += corr14.dhdl_lj;
                    dhdl_nb_crf += corr14.dhdl_crf;
                    for i in 0..n_atoms {
                        self.local_storage.forces[i] += corr14.forces[i];
                    }
                }

                // 3b-5: PERTATOMPAIR correction (special pair types replacing 1-4 LJ)
                if !topo.perturbed_solute.atom_pairs.is_empty() {
                    let mut corr_ap = PertNBCorrection::new(n_atoms);
                    perturbed_atom_pair_correction(
                        &topo.perturbed_solute.atom_pairs,
                        &topo.one_four_pairs,
                        &conf.current().pos,
                        &self.charges,
                        &self.iac_u32,
                        &self.lj_params,
                        &self.crf_params,
                        lp,
                        &self.periodicity,
                        self.four_pi_eps_i,
                        self.rf_excluded,
                        &mut corr_ap,
                        egroups,
                        &mut self.local_storage.groups,
                    );
                    self.local_storage.e_lj += corr_ap.delta_e_lj;
                    self.local_storage.e_crf += corr_ap.delta_e_crf;
                    dhdl_nb += corr_ap.dhdl;
                    dhdl_nb_lj += corr_ap.dhdl_lj;
                    dhdl_nb_crf += corr_ap.dhdl_crf;
                    for i in 0..n_atoms {
                        self.local_storage.forces[i] += corr_ap.forces[i];
                    }
                    log::debug!(
                        "  AtomPair corr: Δe_lj={:.6e} Δe_crf={:.6e}",
                        corr_ap.delta_e_lj,
                        corr_ap.delta_e_crf
                    );
                }

                log::debug!(
                "  Perturbed NB: Δe_lj_pl={:.6e} Δe_crf_pl={:.6e} Δe_self={:.6e} dH/dλ_nb={:.6e}",
                corr_pl.delta_e_lj,
                corr_pl.delta_e_crf,
                de_self,
                dhdl_nb
            );
                log::debug!(
                    "  nb_storage after all: e_lj={:.6e} e_crf={:.6e}",
                    self.nonbonded_storage.e_lj,
                    self.nonbonded_storage.e_crf
                );

                (dhdl_nb, dhdl_nb_lj, dhdl_nb_crf)
            } else {
                (0.0, 0.0, 0.0)
            };

        // --- 4. Assemble forces and energies ---
        // Sum the pair terms across MPI ranks (identity without a reducer), then add the terms
        // every rank computed in full.
        if let Some(reduce) = self.nonbonded_reducer.as_mut() {
            reduce(&mut self.nonbonded_storage);
        }
        let pair_dhdl_lj = self.nonbonded_storage.dhdl_lj;
        let pair_dhdl_crf = self.nonbonded_storage.dhdl_crf;
        let local = std::mem::replace(&mut self.local_storage, ForceStorage::new(0));
        self.nonbonded_storage.merge(&local);
        self.local_storage = local;
        let pert_nb_dhdl = pert_nb_dhdl + pair_dhdl_lj + pair_dhdl_crf;
        let pert_nb_dhdl_lj = pert_nb_dhdl_lj + pair_dhdl_lj;
        let pert_nb_dhdl_crf = pert_nb_dhdl_crf + pair_dhdl_crf;
        {
            let state = conf.current_mut();
            let has_bonded = !bonded_result.forces.is_empty();
            for i in 0..n_atoms {
                let bonded_f = if has_bonded {
                    bonded_result.forces[i]
                } else {
                    gromos_core::math::Vec3::ZERO
                };
                let pert_f = perturbed_bonded
                    .as_ref()
                    .and_then(|p| p.forces.get(i).copied())
                    .unwrap_or(gromos_core::math::Vec3::ZERO);
                state.force[i] = bonded_f + pert_f + self.nonbonded_storage.forces[i];
            }

            // Perturbed (FEP) bonded terms are booked with their unperturbed counterparts, as
            // gromosXX does (before 0.0.34 they were all folded into bond_total — the totals
            // were right, the columns were not).
            let pb = perturbed_bonded.as_ref();
            state.energies.bond_total =
                bonded_result.bond_energy + pb.map_or(0.0, |p| p.bond_energy);
            state.energies.angle_total =
                bonded_result.angle_energy + pb.map_or(0.0, |p| p.angle_energy);
            state.energies.dihedral_total =
                bonded_result.dihedral_energy + pb.map_or(0.0, |p| p.dihedral_energy);
            state.energies.improper_total =
                bonded_result.improper_energy + pb.map_or(0.0, |p| p.improper_energy);
            let dhdl_bonded = perturbed_bonded
                .as_ref()
                .map_or(0.0, |p| p.lambda_derivative);
            state.energies.lj_total = self.nonbonded_storage.e_lj;
            state.energies.crf_total = self.nonbonded_storage.e_crf;
            state.energies.update_potential_total();
            // dH/dλ per term, as gromosXX books it (soft terms count as their kind)
            let d = &mut state.perturbed_energy_derivatives;
            d.bond_total = pb.map_or(0.0, |p| p.bond_dhdl);
            d.angle_total = pb.map_or(0.0, |p| p.angle_dhdl);
            d.improper_total = pb.map_or(0.0, |p| p.improper_dhdl);
            d.dihedral_total = pb.map_or(0.0, |p| p.dihedral_dhdl);
            d.lj_total = pert_nb_dhdl_lj;
            d.crf_total = pert_nb_dhdl_crf;
            d.potential_total = dhdl_bonded + pert_nb_dhdl;
            // perturbed SHAKE adds its dH/dλ here after the step's position update
            // (gromosXX `perturbed_shake.cc:224`, on the state the force field just filled)
            d.constraint_total = 0.0;
            d.special_total = 0.0;
            d.constraints_energy.fill(0.0);

            // Per energy group (more than one group; with one, `EnergyCalculation` books the
            // totals into group 0). Both records are double-buffered and never cleared, so
            // every cell is assigned here.
            if egroups.is_grouped() {
                let n = self.n_energy_groups;
                let g = &self.nonbonded_storage.groups;
                let (e, d) = (&mut state.energies, &mut state.perturbed_energy_derivatives);
                if e.lj_energy.len() != n || d.lj_energy.len() != n {
                    let baths = e.kinetic_energy.len();
                    e.resize(n, baths);
                    d.resize(n, baths);
                }
                for i in 0..n {
                    for j in 0..n {
                        let k = i * n + j;
                        e.lj_energy[i][j] = g.lj[k];
                        e.crf_energy[i][j] = g.crf[k];
                        d.lj_energy[i][j] = g.dhdl_lj[k];
                        d.crf_energy[i][j] = g.dhdl_crf[k];
                    }
                }
                let (be, bd) = (&self.bonded_groups, &self.bonded_dhdl_groups);
                e.bond_energy.copy_from_slice(&be.bond);
                e.angle_energy.copy_from_slice(&be.angle);
                e.improper_energy.copy_from_slice(&be.improper);
                e.dihedral_energy.copy_from_slice(&be.dihedral);
                d.bond_energy.copy_from_slice(&bd.bond);
                d.angle_energy.copy_from_slice(&bd.angle);
                d.improper_energy.copy_from_slice(&bd.improper);
                d.dihedral_energy.copy_from_slice(&bd.dihedral);
            }

            // Transfer virial tensor to configuration (nonbonded + bonded)
            // GROMOS convention: virial[a][b] = Σ r[a] * force[b]
            let vir = &self.nonbonded_storage.virial;
            let bvir = &bonded_result.virial;
            state.virial_tensor = Mat3::from_cols(
                Vec3::new(
                    vir[0][0] + bvir[0][0],
                    vir[1][0] + bvir[1][0],
                    vir[2][0] + bvir[2][0],
                ),
                Vec3::new(
                    vir[0][1] + bvir[0][1],
                    vir[1][1] + bvir[1][1],
                    vir[2][1] + bvir[2][1],
                ),
                Vec3::new(
                    vir[0][2] + bvir[0][2],
                    vir[1][2] + bvir[1][2],
                    vir[2][2] + bvir[2][2],
                ),
            );

            // Debug: show force magnitudes
            let f_max = state
                .force
                .iter()
                .map(|f| f.length())
                .fold(0.0_f64, f64::max);
            log::debug!(
                "  Bond: {:.10e}  LJ: {:.10e}  CRF: {:.10e}",
                bonded_result.energy,
                self.nonbonded_storage.e_lj,
                self.nonbonded_storage.e_crf
            );
            log::debug!(
                "  Max |force|: {:.10e}, solute_pairs: {}, solvent_pairs: {}",
                f_max,
                self.pairlist.solute_short.len(),
                self.pairlist.solvent_short.len()
            );
        }

        // --- 4b. Position restraints (special interaction) ---
        if !self.position_restraints.restraints.is_empty() {
            let e_posres = if egroups.is_grouped() {
                self.restraint_groups.fill(0.0);
                let e = self.position_restraints.calculate_all_grouped(
                    conf,
                    &self.periodicity,
                    &topo.atom_to_energy_group,
                    &mut self.restraint_groups,
                );
                conf.current_mut()
                    .energies
                    .posrest_energy
                    .copy_from_slice(&self.restraint_groups);
                e
            } else {
                self.position_restraints
                    .calculate_all(conf, &self.periodicity)
            };
            let e = &mut conf.current_mut().energies;
            e.special_total = e_posres;
            e.posrest_total = e_posres;
            log::debug!(
                "  Position restraints: {:.10e} kJ/mol ({} atoms)",
                e_posres,
                self.position_restraints.restraints.len()
            );
        }

        // --- 4c. Distance restraints ---
        // Energies and dH/dλ are assigned (not accumulated) each step: both records are
        // double-buffered and never cleared.
        let any_distres = !self.distance_restraints.restraints.is_empty()
            || !self.perturbed_distance_restraints.restraints.is_empty();
        if any_distres {
            let grouped = egroups.is_grouped();
            if grouped {
                self.restraint_groups.fill(0.0);
                self.restraint_dhdl_groups.fill(0.0);
            }
            let mut e_dr = 0.0;
            if !self.distance_restraints.restraints.is_empty() {
                e_dr = if grouped {
                    self.distance_restraints.calculate_all_grouped(
                        conf,
                        &self.periodicity,
                        &topo.atom_to_energy_group,
                        &mut self.restraint_groups,
                    )
                } else {
                    self.distance_restraints
                        .calculate_all(conf, &self.periodicity)
                };
                log::debug!("  Distance restraints: {:.10e} kJ/mol", e_dr);
            }
            let (mut e_pdr, mut d_pdr) = (0.0, 0.0);
            if !self.perturbed_distance_restraints.restraints.is_empty() {
                let lam = TermLambda::table(topo, InteractionLambda::Disres);
                (e_pdr, d_pdr) = if grouped {
                    self.perturbed_distance_restraints.calculate_all_grouped(
                        conf,
                        &self.periodicity,
                        lam,
                        &topo.atom_to_energy_group,
                        &mut self.restraint_groups,
                        &mut self.restraint_dhdl_groups,
                    )
                } else {
                    self.perturbed_distance_restraints
                        .calculate_all(conf, &self.periodicity, lam)
                };
                log::debug!(
                    "  Perturbed distance restraints: {:.10e} kJ/mol, dH/dλ {:.10e}",
                    e_pdr,
                    d_pdr
                );
            }
            let state = conf.current_mut();
            state.energies.distanceres_total = e_dr + e_pdr;
            state.perturbed_energy_derivatives.distanceres_total = d_pdr;
            if grouped {
                state
                    .energies
                    .distanceres_energy
                    .copy_from_slice(&self.restraint_groups);
                state
                    .perturbed_energy_derivatives
                    .distanceres_energy
                    .copy_from_slice(&self.restraint_dhdl_groups);
            }
        }

        // --- 5. atomic_to_molecular_virial: correct virial from atomic to molecular ---
        // GROMOS: last interaction in forcefield, subtracts intramolecular virial contributions
        if self.virial_type == VirialType::Molecular && !topo.pressure_groups.is_empty() {
            apply_molecular_virial_correction(topo, conf, &self.periodicity);
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "Forcefield"
    }

    fn enable_timing(&mut self, enabled: bool) {
        self.timing_enabled = enabled;
    }

    fn sub_timings(&self) -> Vec<(&'static str, std::time::Duration)> {
        vec![
            ("pairlist update", self.timing.pairlist),
            ("CG pair groups", self.timing.cg_groups),
            ("bonded", self.timing.bonded),
            ("shortrange solute", self.timing.solute),
            ("shortrange solvent-solvent", self.timing.solvent),
            ("longrange", self.timing.longrange),
            ("RF excluded interaction", self.timing.rf_excluded),
        ]
    }
}

/// Group a charge-group-ordered atom pairlist into runs of pairs sharing the
/// same charge-group pair.
///
/// The CG-grouped innerloop uses one `nearest_image` per group instead of one
/// per atom pair — for SPC water that is one per nine pairs. Relies on the
/// pairlist being emitted in charge-group order, which every
/// `PairlistAlgorithm` guarantees.
fn build_cg_groups(topo: &Topology, pairs: &[(u32, u32)], groups: &mut Vec<CGPairGroup>) {
    groups.clear();
    if pairs.is_empty() {
        return;
    }

    // Runs are found per contiguous chunk of the pairlist in parallel; a run
    // that straddles a chunk boundary comes out as two groups with the same
    // charge-group pair, and is stitched back together below. Group boundaries
    // are therefore exactly the serial ones, so the grouped kernel sees the
    // same blocks and produces the same numbers.
    let n_threads = rayon::current_num_threads().max(1);
    let chunk = pairs.len().div_ceil(4 * n_threads).max(1);
    let partials: Vec<Vec<CGPairGroup>> = pairs
        .par_chunks(chunk)
        .enumerate()
        .map(|(c, slice)| build_cg_groups_serial(topo, slice, (c * chunk) as u32))
        .collect();

    let atom_to_cg = &topo.atom_to_chargegroup;
    let cg_pair_of = |g: &CGPairGroup| {
        let (i, j) = pairs[g.start as usize];
        (atom_to_cg[i as usize], atom_to_cg[j as usize])
    };
    for partial in partials {
        let mut rest = partial.as_slice();
        if let (Some(last), Some(first)) = (groups.last_mut(), rest.first()) {
            if cg_pair_of(last) == cg_pair_of(first) {
                last.end = first.end;
                rest = &rest[1..];
            }
        }
        groups.extend_from_slice(rest);
    }
}

/// Group a slice of pairs into runs of identical charge-group pairs;
/// `offset` is the slice's start index in the full pairlist.
fn build_cg_groups_serial(topo: &Topology, pairs: &[(u32, u32)], offset: u32) -> Vec<CGPairGroup> {
    let mut groups = Vec::new();
    if pairs.is_empty() {
        return groups;
    }

    let atom_to_cg = &topo.atom_to_chargegroup;
    let mut start = offset;
    let (first_i, first_j) = pairs[0];
    let mut cur_cg_i = atom_to_cg[first_i as usize];
    let mut cur_cg_j = atom_to_cg[first_j as usize];
    let mut ref_atom_i = topo.chargegroups[cur_cg_i].atoms[0] as u32;
    let mut ref_atom_j = topo.chargegroups[cur_cg_j].atoms[0] as u32;

    for (idx, &(ii, jj)) in pairs.iter().enumerate().skip(1) {
        let cg_i = atom_to_cg[ii as usize];
        let cg_j = atom_to_cg[jj as usize];
        if cg_i != cur_cg_i || cg_j != cur_cg_j {
            groups.push(CGPairGroup {
                ref_atom_i,
                ref_atom_j,
                start,
                end: offset + idx as u32,
            });
            start = offset + idx as u32;
            cur_cg_i = cg_i;
            cur_cg_j = cg_j;
            ref_atom_i = topo.chargegroups[cg_i].atoms[0] as u32;
            ref_atom_j = topo.chargegroups[cg_j].atoms[0] as u32;
        }
    }

    groups.push(CGPairGroup {
        ref_atom_i,
        ref_atom_j,
        start,
        end: offset + pairs.len() as u32,
    });
    groups
}

// === Virial helper functions (GROMOS: prepare_virial.cc) ===

/// Compute atomic kinetic energy tensor: KE_ij = 0.5 * Σ_k m_k * v_k_i * v_k_j
/// Uses conf.current().vel (GROMOS convention: computed before force calc)
fn compute_atomic_ke_tensor(topo: &Topology, conf: &Configuration, n_atoms: usize) -> Mat3 {
    let vel = &conf.current().vel;
    let mut ke = [[0.0f64; 3]; 3];
    for k in 0..n_atoms {
        let m = topo.mass[k];
        let v = vel[k];
        let vv = [v.x, v.y, v.z];
        for a in 0..3 {
            for b in 0..3 {
                ke[a][b] += 0.5 * m * vv[a] * vv[b];
            }
        }
    }
    Mat3::from_cols(
        Vec3::new(ke[0][0], ke[1][0], ke[2][0]),
        Vec3::new(ke[0][1], ke[1][1], ke[2][1]),
        Vec3::new(ke[0][2], ke[1][2], ke[2][2]),
    )
}

/// Compute molecular kinetic energy tensor using COM velocities per pressure group.
/// GROMOS: prepare_virial.cc _centre_of_mass with chain-gathering.
/// Uses conf.current().vel and conf.current().pos.
fn compute_molecular_ke_tensor(
    topo: &Topology,
    conf: &Configuration,
    _periodicity: &Periodicity,
) -> Mat3 {
    let vel = &conf.current().vel;
    let mut ke = [[0.0f64; 3]; 3];

    for pg in &topo.pressure_groups {
        let mut mv = [0.0f64; 3]; // mass-weighted velocity sum
        let mut m_total = 0.0f64;

        for i in pg.clone() {
            let m = topo.mass[i];
            m_total += m;
            mv[0] += m * vel[i].x;
            mv[1] += m * vel[i].y;
            mv[2] += m * vel[i].z;
        }
        if m_total > 0.0 {
            for a in 0..3 {
                for b in 0..3 {
                    ke[a][b] += 0.5 * mv[a] * mv[b] / m_total;
                }
            }
        }
    }
    Mat3::from_cols(
        Vec3::new(ke[0][0], ke[1][0], ke[2][0]),
        Vec3::new(ke[0][1], ke[1][1], ke[2][1]),
        Vec3::new(ke[0][2], ke[1][2], ke[2][2]),
    )
}

/// Apply molecular virial correction: transform atomic virial to molecular virial.
/// GROMOS: atomic_to_molecular_virial in prepare_virial.cc
///
/// For each pressure group: compute COM position (chain-gathered), then for each atom:
///   r = nearest_image(pos_atom, com_pos)
///   corrP(b, a) += force(a) * r(b)
/// Finally: virial -= corrP
///
/// Operates on conf.current() (GROMOS convention).
fn apply_molecular_virial_correction(
    topo: &Topology,
    conf: &mut Configuration,
    periodicity: &Periodicity,
) {
    let pos = &conf.current().pos;
    let force = &conf.current().force;

    let mut corr = [[0.0f64; 3]; 3];

    for pg in &topo.pressure_groups {
        // Chain-gather COM position (GROMOS: _centre_of_mass)
        let first = pg.start;
        let mut com = [0.0f64; 3];
        let mut m_total = 0.0f64;
        let mut prev = pos[first];

        for i in pg.clone() {
            let m = topo.mass[i];
            m_total += m;
            // nearest_image returns (pos[i] - prev) with PBC
            let p = periodicity.nearest_image(pos[i], prev);
            // gathered position = p + prev
            let gathered = p + prev;
            com[0] += m * gathered.x;
            com[1] += m * gathered.y;
            com[2] += m * gathered.z;
            prev = gathered; // chain: advance reference to gathered position
        }
        if m_total > 0.0 {
            com[0] /= m_total;
            com[1] /= m_total;
            com[2] /= m_total;
        }

        let com_vec = Vec3::new(com[0], com[1], com[2]);

        // Accumulate correction: corrP(b,a) += f(a) * (pos - com)(b)
        for i in pg.clone() {
            let r = periodicity.nearest_image(pos[i], com_vec);
            let f = force[i];
            let rv = [r.x, r.y, r.z];
            let fv = [f.x, f.y, f.z];
            for a in 0..3 {
                for b in 0..3 {
                    corr[b][a] += fv[a] * rv[b];
                }
            }
        }
    }

    // virial -= corrP
    let corr_mat = Mat3::from_cols(
        Vec3::new(corr[0][0], corr[1][0], corr[2][0]),
        Vec3::new(corr[0][1], corr[1][1], corr[2][1]),
        Vec3::new(corr[0][2], corr[1][2], corr[2][2]),
    );
    let state = conf.current_mut();
    state.virial_tensor = Mat3::from_cols(
        state.virial_tensor.x_axis - corr_mat.x_axis,
        state.virial_tensor.y_axis - corr_mat.y_axis,
        state.virial_tensor.z_axis - corr_mat.z_axis,
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::math::{Vacuum, Vec3};
    use gromos_core::pairlist::StandardPairlistAlgorithm;
    use gromos_forces::nonbonded::CRFParameters;

    /// Test that the forcefield algorithm computes correct LJ forces for the pair_lj system.
    #[test]
    fn test_forcefield_pair_lj() {
        // Build topology matching pair_lj reference:
        // 2 argon atoms, sigma=0.3405nm, epsilon=0.9961 kJ/mol
        // C6=6.2647e-3, C12=9.847e-6
        let mut topo = Topology::new();
        topo.mass = vec![39.948, 39.948];
        topo.inverse_mass = vec![1.0 / 39.948, 1.0 / 39.948];
        topo.charge = vec![0.0, 0.0];
        topo.iac = vec![0, 0];

        // Exclusions: neither atom excludes the other
        topo.exclusions = vec![Vec::new(), Vec::new()];

        // LJ parameter matrix (1x1 for single atom type)
        let lj = LJParameters {
            c6: 6.2647e-3,
            c12: 9.847e-6,
            cs6: 6.2647e-3,
            cs12: 9.847e-6,
        };
        let lj_params = vec![vec![lj]];

        // CRF params (epsilon=1, rf_epsilon=1 => pure Coulomb cutoff, but charges=0)
        let crf_params = CRFParameters::new(1.4, 1.0, 1.0, 0.0);

        let periodicity = Periodicity::Vacuum(Vacuum);

        // Pairlist: manually set to contain the single pair
        let mut pairlist = PairlistContainer::new(1.4, 1.4, 0.0);
        pairlist.solute_short = vec![(0, 1)];
        pairlist.update_frequency = 1000; // Don't auto-update in test

        let pairlist_algorithm = PairlistAlgorithm::Standard(StandardPairlistAlgorithm::new(false));

        let mut ff = Forcefield::new(
            lj_params,
            crf_params,
            periodicity,
            pairlist,
            pairlist_algorithm,
        );

        // Configuration: 2 atoms at r=0.35 nm
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.35, 0.0, 0.0)];
        conf.current_mut().vel = vec![Vec3::ZERO, Vec3::ZERO];
        conf.current_mut().force = vec![Vec3::ZERO, Vec3::ZERO];

        let sim = SimulationState::new(0.002, 10);
        ff.init(&topo, &mut conf, &sim).unwrap();
        ff.apply(&topo, &mut conf, &sim).unwrap();

        // Check forces: should be attractive/repulsive along x
        let f0 = conf.current().force[0];
        let f1 = conf.current().force[1];

        // Forces should be equal and opposite (Newton's 3rd law)
        assert!(
            (f0.x + f1.x).abs() < 1e-10,
            "Forces not equal and opposite: f0.x={}, f1.x={}",
            f0.x,
            f1.x
        );
        assert!(f0.y.abs() < 1e-10);
        assert!(f0.z.abs() < 1e-10);

        // At r=0.35nm for Ar-Ar: net force should be repulsive (r < sigma=0.3405*2^(1/6)=0.382)
        // Actually at r=0.35 which is just below the LJ minimum, let's check:
        // F = -(12*C12/r^13 - 6*C6/r^7)
        // For atom 0: F_x should be negative (pushed towards -x)
        let r = 0.35_f64;
        let r6 = r.powi(6);
        let r7 = r6 * r;
        let r12 = r6 * r6;
        let r13 = r12 * r;
        let c6 = 6.2647e-3;
        let c12 = 9.847e-6;
        let expected_fx = -(12.0 * c12 / r13 - 6.0 * c6 / r7);

        assert!(
            (f0.x - expected_fx).abs() < 1e-6,
            "f0.x = {:.10}, expected = {:.10}",
            f0.x,
            expected_fx
        );

        // Check potential energy
        let expected_e_lj = c12 / r12 - c6 / r6;
        let e_lj = conf.current().energies.lj_total;
        assert!(
            (e_lj - expected_e_lj).abs() < 1e-10,
            "E_LJ = {:.10e}, expected = {:.10e}",
            e_lj,
            expected_e_lj
        );
    }
}
