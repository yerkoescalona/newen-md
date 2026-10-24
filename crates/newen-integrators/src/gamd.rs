use crate::constraints::box_periodicity;
use crate::constraints::{shake, ShakeParameters};
use crate::integrator::Integrator;
use crate::thermostats::{berendsen_thermostat, BerendsenThermostatParameters};
/// Gaussian Accelerated Molecular Dynamics (GaMD)
///
/// GaMD enhances sampling by adding a boost potential that smoothens
/// the energy landscape, accelerating transitions over barriers without
/// requiring prior knowledge of transition pathways.
///
/// ## Theory
///
/// The boost potential is applied when the system potential energy
/// is below a threshold E:
///
/// ```text
/// ΔV(r) = ½K(E - V(r))²  when V(r) < E
/// ΔV(r) = 0              when V(r) ≥ E
/// ```
///
/// Where:
/// - K = force constant (computed from statistics)
/// - E = energy threshold (Vmax for lower bound, derived for upper bound)
/// - V(r) = system potential energy
///
/// The boosted force becomes:
/// ```text
/// F_boost = F_orig × [1 + K(E - V)]
/// ```
///
/// ## Boost Modes
///
/// - **Dihedral boost**: Accelerate only dihedral energies
/// - **Total boost**: Accelerate entire potential (including dihedrals)
/// - **Dual boost**: Separate acceleration for dihedrals and non-dihedrals
///
/// ## Parameter Computation
///
/// ### Force Constant K:
/// ```text
/// K = k0 / (Vmax - Vmin)
/// ```
///
/// ### k0 (threshold factor, 0 < k0 ≤ 1):
/// **Lower bound:** k0 = min(1.0, (σ₀/σᵥ) × (Vmax - Vmin) / (Vmax - Vmean))
/// **Upper bound:** k0 = (1 - σ₀/σᵥ) × (Vmax - Vmin) / (Vmean - Vmin)
///
/// ### Energy Threshold E:
/// **Lower bound:** E = Vmax
/// **Upper bound:** E = Vmin + (Vmax - Vmin) / k0
///
/// ## Workflow
///
/// 1. **CMD Search Phase**: Classical MD, collect statistics (no boost)
/// 2. **GaMD Search Phase**: Accelerated MD with adaptive parameters
/// 3. **Production Phase**: Fixed parameters for production runs
///
/// ## References
/// - Miao et al., J. Chem. Theory Comput. 2015, 11, 3584-3595
/// - Pang et al., J. Chem. Theory Comput. 2017, 13, 9-19
use newen_core::configuration::Configuration;
use newen_core::math::Rectangular;
use newen_core::math::Vec3;
use newen_core::pairlist::{PairlistAlgorithm, PairlistContainer, StandardPairlistAlgorithm};
use newen_core::topology::Topology;
use newen_forces::bonded::calculate_bonded_forces;
use newen_forces::nonbonded::{
    lj_crf_innerloop, lj_crf_innerloop_parallel, CRFParameters, EnergyGroups, ForceStorage,
    LJParamMatrix, LJParameters as NBLJParams,
};

/// GaMD search mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMode {
    /// Production: use fixed parameters, no statistics update
    NoSearch,
    /// CMD search: classical MD, collect statistics without boost
    CmdSearch,
    /// GaMD search: accelerated MD with adaptive parameters
    GamdSearch,
}

/// GaMD boost form
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoostForm {
    /// Boost only dihedral potential
    DihedralBoost,
    /// Boost total potential (including dihedrals)
    TotalBoost,
    /// Separate boost for dihedral and non-dihedral
    DualBoost,
}

/// Energy threshold type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThresholdType {
    /// Lower bound: E = Vmax
    LowerBound,
    /// Upper bound: E = Vmin + (Vmax - Vmin) / k0
    UpperBound,
}

/// GaMD statistics for a single energy component
#[derive(Debug, Clone)]
pub struct GamdStatistics {
    /// Maximum energy observed
    pub v_max: f64,
    /// Minimum energy observed
    pub v_min: f64,
    /// Mean energy
    pub v_mean: f64,
    /// Standard deviation of energy
    pub sigma_v: f64,
    /// Welford M2 accumulator for variance calculation
    pub m2: f64,
    /// Number of steps collected
    pub steps: usize,
}

impl GamdStatistics {
    /// Create new statistics tracker
    pub fn new() -> Self {
        GamdStatistics {
            v_max: f64::NEG_INFINITY,
            v_min: f64::INFINITY,
            v_mean: 0.0,
            sigma_v: 0.0,
            m2: 0.0,
            steps: 0,
        }
    }

    /// Update statistics with new energy value using Welford's algorithm
    ///
    /// This is a numerically stable online algorithm for computing
    /// mean and variance without storing all samples.
    pub fn update(&mut self, energy: f64) {
        self.steps += 1;

        // Update max/min
        if energy > self.v_max {
            self.v_max = energy;
        }
        if energy < self.v_min {
            self.v_min = energy;
        }

        // Welford's online algorithm for mean and variance
        let delta = energy - self.v_mean;
        self.v_mean += delta / self.steps as f64;
        let delta2 = energy - self.v_mean;
        self.m2 += delta * delta2;

        // Update standard deviation
        if self.steps > 1 {
            self.sigma_v = (self.m2 / self.steps as f64).sqrt();
        }
    }

    /// Reset statistics
    pub fn reset(&mut self) {
        self.v_max = f64::NEG_INFINITY;
        self.v_min = f64::INFINITY;
        self.v_mean = 0.0;
        self.sigma_v = 0.0;
        self.m2 = 0.0;
        self.steps = 0;
    }

    /// Check if statistics are ready for parameter calculation
    pub fn is_ready(&self) -> bool {
        self.steps > 10 && self.v_max > self.v_min && self.sigma_v > 0.0
    }
}

impl Default for GamdStatistics {
    fn default() -> Self {
        Self::new()
    }
}

/// GaMD acceleration parameters
#[derive(Debug, Clone)]
pub struct GamdParameters {
    /// Search mode
    pub search_mode: SearchMode,

    /// Boost form
    pub boost_form: BoostForm,

    /// Threshold type
    pub threshold_type: ThresholdType,

    /// Dihedral statistics
    pub dih_stats: GamdStatistics,

    /// Total potential statistics
    pub tot_stats: GamdStatistics,

    /// Dihedral: maximum standard deviation (sigma0)
    pub sigma0_dih: f64,

    /// Total potential: maximum standard deviation (sigma0)
    pub sigma0_tot: f64,

    /// Dihedral: k0 parameter (0 < k0 <= 1)
    pub k0_dih: f64,

    /// Total: k0 parameter (0 < k0 <= 1)
    pub k0_tot: f64,

    /// Dihedral: force constant K
    pub k_dih: f64,

    /// Total: force constant K
    pub k_tot: f64,

    /// Dihedral: energy threshold E
    pub e_dih: f64,

    /// Total: energy threshold E
    pub e_tot: f64,

    /// Current boost potential (for output)
    pub boost_potential: f64,

    /// Equilibration steps (statistics reset after equilibration)
    pub equilibration_steps: usize,

    /// Statistics window size (0 = use all data, >0 = sliding window)
    pub window_size: usize,

    /// Current step counter
    pub current_step: usize,
}

impl GamdParameters {
    /// Create new GaMD parameters
    pub fn new(
        search_mode: SearchMode,
        boost_form: BoostForm,
        threshold_type: ThresholdType,
        sigma0_dih: f64,
        sigma0_tot: f64,
    ) -> Self {
        GamdParameters {
            search_mode,
            boost_form,
            threshold_type,
            dih_stats: GamdStatistics::new(),
            tot_stats: GamdStatistics::new(),
            sigma0_dih,
            sigma0_tot,
            k0_dih: 0.0,
            k0_tot: 0.0,
            k_dih: 0.0,
            k_tot: 0.0,
            e_dih: 0.0,
            e_tot: 0.0,
            boost_potential: 0.0,
            equilibration_steps: 0,
            window_size: 0,
            current_step: 0,
        }
    }

    /// Calculate GaMD parameters from statistics
    ///
    /// Returns (k0, k, E) for the given statistics and sigma0
    fn calculate_parameters(&self, stats: &GamdStatistics, sigma0: f64) -> Option<(f64, f64, f64)> {
        if !stats.is_ready() {
            return None;
        }

        let v_max = stats.v_max;
        let v_min = stats.v_min;
        let v_mean = stats.v_mean;
        let sigma_v = stats.sigma_v;

        // Avoid division by zero
        if (v_max - v_min).abs() < 1e-10 || sigma_v < 1e-10 {
            return None;
        }

        let (k0, e) = match self.threshold_type {
            ThresholdType::LowerBound => {
                // Lower bound: E = Vmax
                // k0 = (σ₀/σᵥ) × (Vmax - Vmin) / (Vmax - Vmean)
                let k0 = (sigma0 / sigma_v) * (v_max - v_min) / (v_max - v_mean);
                let k0 = k0.min(1.0); // Clamp to 1.0
                (k0, v_max)
            },
            ThresholdType::UpperBound => {
                // Upper bound: E = Vmin + (Vmax - Vmin) / k0
                // k0 = (1 - σ₀/σᵥ) × (Vmax - Vmin) / (Vmean - Vmin)
                let mut k0 = (1.0 - sigma0 / sigma_v) * (v_max - v_min) / (v_mean - v_min);

                // Check if k0 is valid (0 < k0 <= 1)
                if k0 <= 0.0 || k0 > 1.0 {
                    // Fall back to lower bound
                    k0 = (sigma0 / sigma_v) * (v_max - v_min) / (v_max - v_mean);
                    k0 = k0.min(1.0);
                    (k0, v_max)
                } else {
                    let e = v_min + (v_max - v_min) / k0;
                    (k0, e)
                }
            },
        };

        // Calculate force constant K = k0 / (Vmax - Vmin)
        let k = k0 / (v_max - v_min);

        Some((k0, k, e))
    }

    /// Update dihedral parameters
    pub fn update_dihedral_parameters(&mut self) {
        if let Some((k0, k, e)) = self.calculate_parameters(&self.dih_stats, self.sigma0_dih) {
            self.k0_dih = k0;
            self.k_dih = k;
            self.e_dih = e;
        }
    }

    /// Update total potential parameters
    pub fn update_total_parameters(&mut self) {
        if let Some((k0, k, e)) = self.calculate_parameters(&self.tot_stats, self.sigma0_tot) {
            self.k0_tot = k0;
            self.k_tot = k;
            self.e_tot = e;
        }
    }

    /// Update all parameters
    pub fn update_all_parameters(&mut self) {
        self.update_dihedral_parameters();
        self.update_total_parameters();
    }

    /// Update statistics with current energies
    pub fn update_statistics(&mut self, dihedral_energy: f64, total_energy: f64) {
        self.current_step += 1;

        // Check for equilibration reset
        if self.equilibration_steps > 0 && self.current_step == self.equilibration_steps {
            self.dih_stats.reset();
            self.tot_stats.reset();
        }

        // Check for window reset
        if self.window_size > 0 && self.dih_stats.steps.is_multiple_of(self.window_size) {
            self.dih_stats.reset();
            self.tot_stats.reset();
        }

        // Update statistics
        if self.boost_form == BoostForm::DihedralBoost || self.boost_form == BoostForm::DualBoost {
            self.dih_stats.update(dihedral_energy);
        }

        if self.boost_form == BoostForm::TotalBoost || self.boost_form == BoostForm::DualBoost {
            let non_dihedral_energy = total_energy - dihedral_energy;
            self.tot_stats.update(non_dihedral_energy);
        }

        // For TotalBoost, use the full potential
        if self.boost_form == BoostForm::TotalBoost {
            self.tot_stats.update(total_energy);
        }
    }

    /// Apply GaMD boost to forces
    ///
    /// Returns the total boost potential added
    pub fn apply_boost(
        &mut self,
        configuration: &mut Configuration,
        dihedral_energy: f64,
        dihedral_forces: &[Vec3],
        total_energy: f64,
        total_forces: &[Vec3],
    ) -> f64 {
        self.boost_potential = 0.0;

        // Only apply boost in GaMD search or production modes
        if self.search_mode == SearchMode::CmdSearch {
            return 0.0;
        }

        let n_atoms = configuration.current().pos.len();

        match self.boost_form {
            BoostForm::DihedralBoost => {
                // Boost only dihedral term
                let v_e = self.e_dih - dihedral_energy;
                if v_e > 0.0 {
                    let prefactor = self.k_dih * v_e;
                    self.boost_potential = 0.5 * prefactor * v_e;

                    // Apply boost: F_boost = F_orig * k(E - V)
                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] += dihedral_forces[i] * prefactor;
                    }
                }
            },
            BoostForm::TotalBoost => {
                // Boost entire potential (including dihedrals)
                let v_e = self.e_tot - total_energy;
                if v_e > 0.0 {
                    let prefactor = self.k_tot * v_e;
                    self.boost_potential = 0.5 * prefactor * v_e;

                    // Apply boost to all forces
                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] +=
                            (dihedral_forces[i] + total_forces[i]) * prefactor;
                    }
                }
            },
            BoostForm::DualBoost => {
                // Boost dihedral and non-dihedral independently
                let non_dihedral_energy = total_energy - dihedral_energy;

                // Non-dihedral boost
                let v_e_tot = self.e_tot - non_dihedral_energy;
                if v_e_tot > 0.0 {
                    let prefactor_tot = self.k_tot * v_e_tot;
                    self.boost_potential += 0.5 * prefactor_tot * v_e_tot;

                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] += total_forces[i] * prefactor_tot;
                    }
                }

                // Dihedral boost
                let v_e_dih = self.e_dih - dihedral_energy;
                if v_e_dih > 0.0 {
                    let prefactor_dih = self.k_dih * v_e_dih;
                    self.boost_potential += 0.5 * prefactor_dih * v_e_dih;

                    for i in 0..n_atoms {
                        configuration.current_mut().force[i] += dihedral_forces[i] * prefactor_dih;
                    }
                }
            },
        }

        self.boost_potential
    }
}

/// GaMD MD Runner
///
/// Manages GaMD molecular dynamics simulation with separate
/// tracking of dihedral and non-dihedral forces
pub struct GamdRunner {
    /// GaMD parameters
    pub gamd: GamdParameters,

    /// Pairlist for nonbonded interactions
    pairlist: PairlistContainer,

    /// Pairlist algorithm
    pairlist_algorithm: PairlistAlgorithm,

    /// Use parallel force calculation
    parallel_forces: bool,

    /// CRF parameters for nonbonded interactions
    crf_params: CRFParameters,
    four_pi_eps_i: f64,

    /// Current MD step
    current_step: usize,

    /// Stored dihedral forces
    dihedral_forces: Vec<Vec3>,

    /// Stored non-dihedral forces
    other_forces: Vec<Vec3>,

    /// Dihedral energy
    dihedral_energy: f64,

    /// Non-dihedral energy
    other_energy: f64,
}

impl GamdRunner {
    /// Create a new GaMD runner
    pub fn new(gamd: GamdParameters, cutoff: f64, parallel_forces: bool) -> Self {
        let pairlist = PairlistContainer::new(cutoff, cutoff, 0.2);
        let pairlist_algorithm = PairlistAlgorithm::Standard(StandardPairlistAlgorithm::new(true));

        // Calculate CRF parameters
        let eps = 1.0;
        let eps_rf = 61.0;
        let crf = (eps_rf - eps) * (2.0 * eps + eps_rf)
            / (eps_rf + 2.0 * eps)
            / (cutoff * cutoff * cutoff);

        let crf_params = CRFParameters {
            crf_cut: cutoff,
            crf_2cut3i: crf / (2.0 * cutoff * cutoff * cutoff),
            crf_cut3i: (1.0 - crf / 2.0) / cutoff,
            cutoff_sq: cutoff * cutoff,
        };

        GamdRunner {
            gamd,
            pairlist,
            pairlist_algorithm,
            parallel_forces,
            crf_params,
            four_pi_eps_i: newen_core::units::four_pi_eps_i,
            current_step: 0,
            dihedral_forces: Vec::new(),
            other_forces: Vec::new(),
            dihedral_energy: 0.0,
            other_energy: 0.0,
        }
    }

    /// Update pairlist if needed
    fn update_pairlist_if_needed(&mut self, topology: &Topology, configuration: &Configuration) {
        if self.pairlist.needs_update() {
            let box_dims = configuration.current().box_config.dimensions();
            let periodicity = Rectangular::new(box_dims);

            self.pairlist_algorithm.update(
                topology,
                configuration,
                &mut self.pairlist,
                &periodicity,
            );
            self.pairlist.reset_counter();
        } else {
            self.pairlist.step();
        }
    }

    /// Calculate all forces and energies
    fn calculate_forces(&mut self, topology: &Topology, configuration: &mut Configuration) {
        let n_atoms = topology.num_atoms();

        // Initialize force storage
        self.dihedral_forces = vec![Vec3::ZERO; n_atoms];
        self.other_forces = vec![Vec3::ZERO; n_atoms];
        self.dihedral_energy = 0.0;
        self.other_energy = 0.0;

        // Update pairlist
        self.update_pairlist_if_needed(topology, configuration);

        // Clear configuration forces
        for force in configuration.current_mut().force.iter_mut() {
            *force = Vec3::ZERO;
        }

        // Calculate bonded forces (bonds, angles, dihedrals)
        let bonded_result = calculate_bonded_forces(
            topology,
            configuration,
            &box_periodicity(configuration),
            true,
        );

        // For simplicity, assume all bonded energy is non-dihedral
        // (In a full implementation, you'd separate dihedral contributions)
        // TODO: Separate dihedral forces from bonded calculation
        for (i, &force) in bonded_result.forces.iter().enumerate() {
            self.other_forces[i] += force;
        }
        self.other_energy += bonded_result.energy;

        // Calculate nonbonded forces
        let box_dims = configuration.current().box_config.dimensions();
        let periodicity = Rectangular::new(box_dims);

        // Convert LJ parameters
        let lj_params_nested: Vec<Vec<NBLJParams>> = topology
            .lj_parameters
            .iter()
            .map(|row| {
                row.iter()
                    .map(|params| NBLJParams {
                        c6: params.c6,
                        c12: params.c12,
                        cs6: params.cs6,
                        cs12: params.cs12,
                    })
                    .collect()
            })
            .collect();
        let lj_params = LJParamMatrix::from_nested(&lj_params_nested);

        // Convert pairlist to (u32, u32) format
        let pairlist_short: Vec<(u32, u32)> = self.pairlist.solute_short.to_vec();

        // Convert charge and iac

        let iac_u32: Vec<u32> = topology.iac.iter().map(|&i| i as u32).collect();

        let nonbonded_storage = if self.parallel_forces {
            lj_crf_innerloop_parallel(
                &configuration.current().pos,
                &topology.charge,
                &iac_u32,
                &pairlist_short,
                &lj_params,
                &self.crf_params,
                &periodicity,
                self.four_pi_eps_i,
                n_atoms,
                EnergyGroups::single(),
            )
        } else {
            let mut storage = ForceStorage::new(n_atoms);
            lj_crf_innerloop(
                &configuration.current().pos,
                &topology.charge,
                &iac_u32,
                &pairlist_short,
                &lj_params,
                &self.crf_params,
                &periodicity,
                self.four_pi_eps_i,
                &mut storage,
                EnergyGroups::single(),
            );
            storage
        };

        // Add nonbonded forces to other_forces
        for (i, &force) in nonbonded_storage.forces.iter().enumerate() {
            self.other_forces[i] += force;
        }
        self.other_energy += nonbonded_storage.e_lj + nonbonded_storage.e_crf;

        // Copy forces to configuration (will be modified by GaMD boost)
        for i in 0..n_atoms {
            configuration.current_mut().force[i] = self.dihedral_forces[i] + self.other_forces[i];
        }

        // Update configuration energies
        configuration.current_mut().energies.bond_total = bonded_result.bond_energy;
        configuration.current_mut().energies.angle_total = bonded_result.angle_energy;
        configuration.current_mut().energies.dihedral_total = bonded_result.dihedral_energy;
        configuration.current_mut().energies.improper_total = bonded_result.improper_energy;
        configuration.current_mut().energies.lj_total = nonbonded_storage.e_lj;
        configuration.current_mut().energies.crf_total = nonbonded_storage.e_crf;
        configuration.current_mut().energies.potential_total =
            self.dihedral_energy + self.other_energy;
    }

    /// Perform one GaMD MD step
    pub fn md_step<I: Integrator>(
        &mut self,
        topology: &Topology,
        configuration: &mut Configuration,
        integrator: &mut I,
        dt: f64,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        // Calculate forces
        self.calculate_forces(topology, configuration);

        let total_energy = self.dihedral_energy + self.other_energy;

        // Update statistics (in all search modes)
        self.gamd
            .update_statistics(self.dihedral_energy, total_energy);

        // Update parameters (in GaMD search mode)
        if self.gamd.search_mode == SearchMode::GamdSearch {
            self.gamd.update_all_parameters();
        }

        // Apply GaMD boost
        let boost = self.gamd.apply_boost(
            configuration,
            self.dihedral_energy,
            &self.dihedral_forces,
            total_energy,
            &self.other_forces,
        );

        // Add boost to potential energy
        configuration.current_mut().energies.potential_total += boost;

        // Integrate
        integrator.step(dt, topology, configuration);

        // Apply constraints (SHAKE)
        if let Some(shake_p) = shake_params {
            shake(topology, configuration, dt, shake_p);
        }

        // Calculate kinetic energy
        configuration
            .current_mut()
            .calculate_kinetic_energy(&topology.mass);

        // Apply thermostat
        if let Some(params) = thermostat_params {
            berendsen_thermostat(topology, configuration, dt, params);
        }

        self.current_step += 1;
    }

    /// Run multiple GaMD MD steps
    pub fn run_md_steps<I: Integrator>(
        &mut self,
        topology: &Topology,
        configuration: &mut Configuration,
        integrator: &mut I,
        n_steps: usize,
        dt: f64,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        for _ in 0..n_steps {
            self.md_step(
                topology,
                configuration,
                integrator,
                dt,
                thermostat_params,
                shake_params,
            );
        }
    }

    /// Get current boost potential
    pub fn boost_potential(&self) -> f64 {
        self.gamd.boost_potential
    }

    /// Get dihedral statistics
    pub fn dih_stats(&self) -> &GamdStatistics {
        &self.gamd.dih_stats
    }

    /// Get total potential statistics
    pub fn tot_stats(&self) -> &GamdStatistics {
        &self.gamd.tot_stats
    }
}

// ─── gromosXX's GAMD algorithm (`algorithm/integration/gamd.cc`) ───────────────────────────
//
// The types above are newen's own GaMD runner (the `gamd` binary). What follows is gromosXX's
// algorithm, run by the force field after every force evaluation: per interaction group (index 0
// is the non-accelerated group) the statistics of the dihedral energy V_D (proper + improper
// dihedrals) and of the potential energy V_T (bonds, angles and the nonbonded pair terms), and
// the boost ΔV = k(V − E)²/2 for V < E with the force scaling F += k(V − E)·F_term.

/// One interaction group's GaMD statistics (gromosXX `gamd_struct`, `parameter.h:3364-3440`),
/// in the column order of the `GAMDSTAT` configuration block (`out_configuration.cc:4021-4031`).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct GamdGroupStats {
    /// Welford sum of squared deviations of V_D.
    pub m2_d: f64,
    /// Welford sum of squared deviations of V_T.
    pub m2_t: f64,
    pub vmax_d: f64,
    pub vmax_t: f64,
    pub vmean_d: f64,
    pub vmean_t: f64,
    pub vmin_d: f64,
    pub vmin_t: f64,
    /// Standard deviation of V_D.
    pub sigma_d: f64,
    /// Standard deviation of V_T.
    pub sigma_t: f64,
}

impl GamdGroupStats {
    /// The `GAMDSTAT` row: M2D M2T VmaxD VmaxT VmeanD VmeanT VminD VminT sigmaVD sigmaVT.
    pub fn to_row(&self) -> [f64; 10] {
        [
            self.m2_d,
            self.m2_t,
            self.vmax_d,
            self.vmax_t,
            self.vmean_d,
            self.vmean_t,
            self.vmin_d,
            self.vmin_t,
            self.sigma_d,
            self.sigma_t,
        ]
    }

    /// A group from its `GAMDSTAT` row.
    pub fn from_row(r: [f64; 10]) -> Self {
        Self {
            m2_d: r[0],
            m2_t: r[1],
            vmax_d: r[2],
            vmax_t: r[3],
            vmean_d: r[4],
            vmean_t: r[5],
            vmin_d: r[6],
            vmin_t: r[7],
            sigma_d: r[8],
            sigma_t: r[9],
        }
    }
}

/// The GAMD settings of a run (gromosXX `read_GAMD`, `in_parameter.cc:4148-4320`). The
/// per-group vectors have one entry per interaction group, index 0 being the non-accelerated
/// group (0 from the input, as gromosXX sizes them IGROUPS + 1).
#[derive(Debug, Clone, PartialEq)]
pub struct GamdSettings {
    pub search: SearchMode,
    pub form: BoostForm,
    pub thresh: ThresholdType,
    /// NTIGAMDS = 1: start the search statistics from scratch instead of `GAMDSTAT`.
    pub ntisearch: bool,
    /// DIHSTD: σ₀ of the dihedral boost.
    pub dihstd: f64,
    /// TOTSTD: σ₀ of the potential-energy boost.
    pub totstd: f64,
    /// EQSTEPS: steps before the statistics are collected.
    pub equilibration: i64,
    /// GWINDOW: steps between statistics resets (0 = never).
    pub window: i64,
    pub ed: Vec<f64>,
    pub et: Vec<f64>,
    pub kd: Vec<f64>,
    pub kt: Vec<f64>,
}

/// What one interaction group collects from a force evaluation (gromosXX books every term to
/// the group of its atoms' acceleration groups). Empty force vectors stand for zero forces.
#[derive(Debug, Clone, Default)]
pub struct GamdGroupTerms {
    /// V_D: proper and improper dihedral energy (gromosXX `gamd_dihedral_total`).
    pub v_dihedral: f64,
    /// V_T: bond, angle and nonbonded energy (gromosXX `gamd_potential_total`).
    pub v_potential: f64,
    /// The forces of the dihedral terms (gromosXX `special().gamd.dihe_force`).
    pub f_dihedral: Vec<Vec3>,
    /// The forces of the potential terms (gromosXX `special().gamd.total_force`).
    pub f_potential: Vec<Vec3>,
    /// Their virials, in the layout of `State::virial_tensor`.
    pub vir_dihedral: newen_core::math::Mat3,
    pub vir_potential: newen_core::math::Mat3,
}

/// gromosXX's GAMD algorithm and its state: settings, the k0 and statistics of every
/// interaction group, `stepsdone`, and this step's energy record.
///
/// The force field books its terms to two interaction groups: [`Self::group`], the one every
/// pair of the run's acceleration groups maps to (the caller refuses an assignment that splits
/// them), and [`Self::excluded_group`], where gromosXX books the reaction-field excluded solute
/// pairs. The other groups see V = 0 and no forces, as in gromosXX.
#[derive(Debug, Clone)]
pub struct Gamd {
    pub settings: GamdSettings,
    /// The interaction group every term is booked to (0: none is accelerated).
    pub group: usize,
    /// The interaction group of the reaction-field excluded solute pairs: gromosXX looks it up
    /// with the key (0, group of j) — its `gamdi` is shadowed (`nonbonded_innerloop.cc:1694,
    /// 1711-1723`) — so they usually leave the accelerated group.
    pub excluded_group: usize,
    pub k0d: Vec<f64>,
    pub k0t: Vec<f64>,
    pub stats: Vec<GamdGroupStats>,
    /// gromosXX `stepsdone`: steps since the statistics started (written to `GAMDSTAT`).
    pub stepsdone: i64,
    /// This step's `# gamd` rows of the energy record, per interaction group: ED, ET, KD, KT and
    /// ΔV (gromosXX `gamd_ED`, `gamd_ET`, `gamd_KD`, `gamd_KT`, `gamd_DV`).
    pub rows: Vec<[f64; 5]>,
    /// gromosXX's warnings (`io::messages`), in the order they were raised.
    pub warnings: Vec<String>,
}

const K0_WARNING: &str = "gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ";

impl Gamd {
    /// The algorithm for `settings` (vectors of `n` = IGROUPS + 1 entries), with the statistics
    /// a run continues from (`GAMDSTAT`; zero for a fresh start). Runs gromosXX's `GAMD::init`
    /// (`gamd.cc:44-97`): without a cMD search and with NTIGAMDS = 0, a group whose KD and KT are
    /// both 0 takes E and k from the statistics.
    pub fn new(
        settings: GamdSettings,
        (group, excluded_group): (usize, usize),
        stats: Vec<GamdGroupStats>,
        stepsdone: i64,
    ) -> Self {
        let n = settings.ed.len();
        let mut g = Self {
            settings,
            group,
            excluded_group,
            k0d: vec![0.0; n],
            k0t: vec![0.0; n],
            stats,
            stepsdone,
            rows: vec![[0.0; 5]; n],
            warnings: Vec::new(),
        };
        g.stats.resize(n, GamdGroupStats::default());
        let (search, ntisearch) = (g.settings.search, g.settings.ntisearch);
        if search != SearchMode::CmdSearch && !ntisearch {
            for gg in 1..n {
                if g.settings.kd[gg] == 0.0 && g.settings.kt[gg] == 0.0 {
                    match g.settings.form {
                        BoostForm::DihedralBoost => g.update_dihedral(gg),
                        BoostForm::TotalBoost => g.update_potential(gg),
                        BoostForm::DualBoost => {
                            g.update_dihedral(gg);
                            g.update_potential(gg);
                        },
                    }
                }
            }
        } else if search == SearchMode::GamdSearch && ntisearch {
            g.warnings.push(
                "GAMD search should be started from an cMD search with ntisearch = 0. With \
                 ntisearch = 1 the statistics will be estimated only from the actual GAMD search \
                 run"
                .to_string(),
            );
        }
        g
    }

    /// Number of interaction groups, the non-accelerated one included.
    pub fn n_groups(&self) -> usize {
        self.rows.len()
    }

    /// ΔV of this step, summed over the groups (gromosXX `gamd_DV_total`, ENER[51]).
    pub fn dv_total(&self) -> f64 {
        self.rows.iter().map(|r| r[4]).sum()
    }

    /// E and k of the dihedral boost of group `gg` from its statistics.
    fn update_dihedral(&mut self, gg: usize) {
        let st = self.stats[gg];
        let s = &mut self.settings;
        if calc_gamd_e_k(
            s.thresh,
            s.dihstd,
            st.vmax_d,
            st.vmin_d,
            st.vmean_d,
            st.sigma_d,
            &mut self.k0d[gg],
            &mut s.kd[gg],
            &mut s.ed[gg],
        ) {
            self.warnings.push(K0_WARNING.to_string());
            s.thresh = ThresholdType::LowerBound;
        }
    }

    /// E and k of the potential-energy boost of group `gg` from its statistics.
    fn update_potential(&mut self, gg: usize) {
        let st = self.stats[gg];
        let s = &mut self.settings;
        if calc_gamd_e_k(
            s.thresh,
            s.totstd,
            st.vmax_t,
            st.vmin_t,
            st.vmean_t,
            st.sigma_t,
            &mut self.k0t[gg],
            &mut s.kt[gg],
            &mut s.et[gg],
        ) {
            self.warnings.push(K0_WARNING.to_string());
            s.thresh = ThresholdType::LowerBound;
        }
    }

    /// One step of gromosXX's `GAMD::apply` (`gamd.cc:104-411`) after the force evaluation of
    /// step `step`, `terms[gg]` being what interaction group `gg` collected (missing groups
    /// collected nothing): the statistics (from EQSTEPS on), then — except in a cMD search — the
    /// boost on `force` and `virial`. Returns ΔV (to be booked in the special energy).
    pub fn apply(
        &mut self,
        step: usize,
        terms: &[GamdGroupTerms],
        force: &mut [Vec3],
        virial: &mut newen_core::math::Mat3,
    ) -> f64 {
        let n = self.n_groups();
        let none = GamdGroupTerms::default();
        let group_terms = |gg: usize| terms.get(gg).unwrap_or(&none);
        let energies = |gg: usize| {
            let t = group_terms(gg);
            (t.v_dihedral, t.v_potential)
        };
        let eq = self.settings.equilibration;
        let window = self.settings.window;
        if step as i64 >= eq {
            self.stepsdone += 1;
            let steps = self.stepsdone;
            for gg in 1..n {
                // reset the statistics after the equilibration (`gamd.cc:122-127`)
                let st = &mut self.stats[gg];
                if eq == step as i64 && eq > 0 {
                    st.vmax_t = st.vmean_t;
                    st.vmin_t = st.vmean_t;
                    st.vmax_d = st.vmean_d;
                    st.vmin_d = st.vmean_d;
                    st.m2_t = 0.0;
                    st.m2_d = 0.0;
                    st.sigma_t = 0.0;
                    st.sigma_d = 0.0;
                }
                let (vd, vt) = energies(gg);
                let reset = window > 0 && steps % window == 0;
                let collect = |v: f64, st: &mut GamdGroupStats, dihedral: bool| {
                    let (vmax, vmin, vmean, m2, sigma) = if dihedral {
                        (
                            &mut st.vmax_d,
                            &mut st.vmin_d,
                            &mut st.vmean_d,
                            &mut st.m2_d,
                            &mut st.sigma_d,
                        )
                    } else {
                        (
                            &mut st.vmax_t,
                            &mut st.vmin_t,
                            &mut st.vmean_t,
                            &mut st.m2_t,
                            &mut st.sigma_t,
                        )
                    };
                    if reset {
                        // a window reset keeps Vmax/Vmin and restarts mean and M2 at V
                        if v > *vmax {
                            *vmax = v;
                        }
                        if v < *vmin {
                            *vmin = v;
                        }
                        *m2 = 0.0;
                        *sigma = 0.0;
                        *vmean = v;
                    } else {
                        calc_gamd_std_mean(v, steps, vmax, vmin, vmean, m2, sigma);
                    }
                };
                match self.settings.search {
                    SearchMode::NoSearch => {},
                    SearchMode::CmdSearch | SearchMode::GamdSearch => {
                        let st = &mut self.stats[gg];
                        match self.settings.form {
                            BoostForm::DihedralBoost => collect(vd, st, true),
                            BoostForm::TotalBoost => collect(vt + vd, st, false),
                            BoostForm::DualBoost => {
                                collect(vd, st, true);
                                collect(vt, st, false);
                            },
                        }
                        if self.settings.search == SearchMode::GamdSearch {
                            match self.settings.form {
                                BoostForm::DihedralBoost => self.update_dihedral(gg),
                                BoostForm::TotalBoost => self.update_potential(gg),
                                BoostForm::DualBoost => {
                                    self.update_dihedral(gg);
                                    self.update_potential(gg);
                                },
                            }
                        }
                    },
                }
            }
        }

        // the energy record is zeroed every step (`energy.cc:164-173`)
        for r in self.rows.iter_mut() {
            *r = [0.0; 5];
        }
        if self.settings.search != SearchMode::CmdSearch {
            for gg in 1..n {
                let (vd, vt) = energies(gg);
                let split = group_terms(gg);
                let s = &self.settings;
                let mut dv = 0.0;
                let add = |f: &[Vec3], p: f64, force: &mut [Vec3]| {
                    for (fi, ti) in force.iter_mut().zip(f) {
                        *fi += *ti * p;
                    }
                };
                match s.form {
                    BoostForm::DihedralBoost => {
                        let ve = vd - s.ed[gg];
                        if ve < 0.0 {
                            let p = s.kd[gg] * ve;
                            dv = p * ve / 2.0;
                            add(&split.f_dihedral, p, force);
                            *virial += split.vir_dihedral * p;
                        }
                    },
                    BoostForm::TotalBoost => {
                        let ve = vt + vd - s.et[gg];
                        if ve < 0.0 {
                            let p = s.kt[gg] * ve;
                            dv = p * ve / 2.0;
                            let zero = Vec3::ZERO;
                            for (i, fi) in force.iter_mut().enumerate() {
                                *fi += *split.f_dihedral.get(i).unwrap_or(&zero) * p;
                                *fi += *split.f_potential.get(i).unwrap_or(&zero) * p;
                            }
                            *virial += split.vir_dihedral * p;
                            *virial += split.vir_potential * p;
                        }
                    },
                    BoostForm::DualBoost => {
                        let vet = vt - s.et[gg];
                        if vet < 0.0 {
                            let p = s.kt[gg] * vet;
                            dv = p * vet / 2.0;
                            add(&split.f_potential, p, force);
                            *virial += split.vir_potential * p;
                        }
                        let ved = vd - s.ed[gg];
                        if ved < 0.0 {
                            let p = s.kd[gg] * ved;
                            // gromosXX assigns, so a dihedral boost replaces the potential
                            // one in ΔV while both act on the forces (`gamd.cc:364, 381`)
                            dv = p * ved / 2.0;
                            add(&split.f_dihedral, p, force);
                            *virial += split.vir_dihedral * p;
                        }
                    },
                }
                let row = [s.ed[gg], s.et[gg], s.kd[gg], s.kt[gg], dv];
                self.rows[gg] = row;
            }
        }
        self.dv_total()
    }
}

/// gromosXX `GAMD::calc_gamd_std_mean` (`gamd.cc:415-433`): Vmax/Vmin (both start at the first
/// V when still 0) and Welford's running mean and M2 over `step` samples.
pub fn calc_gamd_std_mean(
    v: f64,
    step: i64,
    vmax: &mut f64,
    vmin: &mut f64,
    vmean: &mut f64,
    m2: &mut f64,
    sigma: &mut f64,
) {
    if *vmax == 0.0 && *vmin == 0.0 {
        *vmin = v;
        *vmax = v;
    }
    if v > *vmax {
        *vmax = v;
    } else if v < *vmin {
        *vmin = v;
    }
    let delta = v - *vmean;
    *vmean += delta / step as f64;
    *m2 += delta * (v - *vmean);
    *sigma = (*m2 / step as f64).sqrt();
}

/// gromosXX `GAMD::calc_gamd_E_K` (`gamd.cc:435-468`): k0, k and E from the statistics, for the
/// lower bound E = Vmax or the upper bound E = Vmin + (Vmax − Vmin)/k0. An upper-bound k0
/// outside (0, 1] falls back to the lower bound and returns `true` (the caller then switches
/// THRESH to the lower bound). Nothing changes while Vmax, Vmin and Vmean are all 0.
#[allow(clippy::too_many_arguments)]
pub fn calc_gamd_e_k(
    thresh: ThresholdType,
    sigma0: f64,
    vmax: f64,
    vmin: f64,
    vmean: f64,
    sigma: f64,
    k0: &mut f64,
    k: &mut f64,
    e: &mut f64,
) -> bool {
    if vmax == 0.0 && vmin == 0.0 && vmean == 0.0 {
        return false;
    }
    let lower = |k0: &mut f64, k: &mut f64, e: &mut f64| {
        *k0 = (sigma0 / sigma) * (vmax - vmin) / (vmax - vmean);
        if *k0 > 1.0 {
            *k0 = 1.0;
        }
        *e = vmax;
        *k = *k0 / (vmax - vmin);
    };
    match thresh {
        ThresholdType::UpperBound => {
            *k0 = (1.0 - sigma0 / sigma) * (vmax - vmin) / (vmean - vmin);
            if *k0 > 1.0 || *k0 <= 0.0 {
                lower(k0, k, e);
                true
            } else {
                *e = vmin + (vmax - vmin) / *k0;
                *k = *k0 / (vmax - vmin);
                false
            }
        },
        ThresholdType::LowerBound => {
            lower(k0, k, e);
            false
        },
    }
}

#[cfg(test)]
mod gromos_gamd_tests {
    use super::*;
    use newen_core::math::Mat3;

    fn settings(search: SearchMode, form: BoostForm) -> GamdSettings {
        GamdSettings {
            search,
            form,
            thresh: ThresholdType::LowerBound,
            ntisearch: true,
            dihstd: 24.79,
            totstd: 24.79,
            equilibration: 0,
            window: 0,
            ed: vec![0.0, 30.0],
            et: vec![0.0, -100.0],
            kd: vec![0.0, 0.005],
            kt: vec![0.0, 0.0005],
        }
    }

    fn split(vd: f64, vt: f64) -> Vec<GamdGroupTerms> {
        let own = GamdGroupTerms {
            v_dihedral: vd,
            v_potential: vt,
            f_dihedral: vec![Vec3::new(1.0, 0.0, 0.0)],
            f_potential: vec![Vec3::new(0.0, 2.0, 0.0)],
            vir_dihedral: Mat3::IDENTITY,
            vir_potential: Mat3::IDENTITY * 2.0,
        };
        vec![GamdGroupTerms::default(), own]
    }

    /// FORM 3 (dihedral boost), fixed parameters: ΔV = kD(V_D − ED)²/2 and F += kD(V_D − ED)·F_D.
    #[test]
    fn dihedral_boost_scales_the_dihedral_forces() {
        let mut g = Gamd::new(
            settings(SearchMode::NoSearch, BoostForm::DihedralBoost),
            (1, 0),
            vec![],
            0,
        );
        let mut f = vec![Vec3::ZERO];
        let mut vir = Mat3::ZERO;
        let dv = g.apply(0, &split(10.0, -50.0), &mut f, &mut vir);
        let p = 0.005 * (10.0 - 30.0);
        assert_eq!(dv, p * (10.0 - 30.0) / 2.0);
        assert_eq!(f[0], Vec3::new(p, 0.0, 0.0));
        assert_eq!(vir, Mat3::IDENTITY * p);
        assert_eq!(g.rows[0], [0.0; 5]);
        assert_eq!(g.rows[1], [30.0, -100.0, 0.005, 0.0005, dv]);
    }

    /// FORM 1 (dual boost): both boosts act on the forces, but ΔV keeps only the dihedral one
    /// (gromosXX assigns `gamd_DV` twice).
    #[test]
    fn dual_boost_reports_only_the_dihedral_dv() {
        let mut g = Gamd::new(
            settings(SearchMode::NoSearch, BoostForm::DualBoost),
            (1, 0),
            vec![],
            0,
        );
        let mut f = vec![Vec3::ZERO];
        let mut vir = Mat3::ZERO;
        let dv = g.apply(0, &split(10.0, -150.0), &mut f, &mut vir);
        let pt = 0.0005 * (-150.0 + 100.0);
        let ph = 0.005 * (10.0 - 30.0);
        assert_eq!(dv, ph * (10.0 - 30.0) / 2.0);
        assert_eq!(f[0], Vec3::new(ph, 2.0 * pt, 0.0));
    }

    /// A cMD search collects statistics and applies no boost; GAMDSTAT's numbers follow
    /// Welford's update with Vmax/Vmin starting at the first sample.
    #[test]
    fn cmd_search_collects_statistics_only() {
        let mut g = Gamd::new(
            settings(SearchMode::CmdSearch, BoostForm::DihedralBoost),
            (1, 0),
            vec![],
            0,
        );
        let mut f = vec![Vec3::ZERO];
        let mut vir = Mat3::ZERO;
        for (step, v) in [4.0, 8.0, 6.0].into_iter().enumerate() {
            assert_eq!(g.apply(step, &split(v, 0.0), &mut f, &mut vir), 0.0);
        }
        assert_eq!(f[0], Vec3::ZERO);
        let st = g.stats[1];
        assert_eq!(g.stepsdone, 3);
        assert_eq!((st.vmax_d, st.vmin_d, st.vmean_d), (8.0, 4.0, 6.0));
        assert!((st.m2_d - 8.0).abs() < 1e-12);
        assert!((st.sigma_d - (8.0f64 / 3.0).sqrt()).abs() < 1e-12);
        assert_eq!(g.rows, vec![[0.0; 5]; 2]);
    }

    /// The lower bound: E = Vmax, k0 = min(1, σ₀/σ·(Vmax − Vmin)/(Vmax − Vmean)); an invalid
    /// upper-bound k0 falls back to it and reports the switch.
    #[test]
    fn energy_threshold_and_force_constant() {
        let (mut k0, mut k, mut e) = (0.0, 0.0, 0.0);
        assert!(!calc_gamd_e_k(
            ThresholdType::LowerBound,
            1.0,
            10.0,
            0.0,
            5.0,
            4.0,
            &mut k0,
            &mut k,
            &mut e
        ));
        assert_eq!((k0, e), (0.5, 10.0));
        assert_eq!(k, 0.05);
        assert!(calc_gamd_e_k(
            ThresholdType::UpperBound,
            8.0,
            10.0,
            0.0,
            5.0,
            4.0,
            &mut k0,
            &mut k,
            &mut e
        ));
        assert_eq!((k0, e), (1.0, 10.0));
        k0 = 7.0;
        assert!(!calc_gamd_e_k(
            ThresholdType::UpperBound,
            1.0,
            0.0,
            0.0,
            0.0,
            4.0,
            &mut k0,
            &mut k,
            &mut e
        ));
        assert_eq!(k0, 7.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_statistics_update() {
        let mut stats = GamdStatistics::new();

        // Add some samples
        stats.update(100.0);
        stats.update(110.0);
        stats.update(90.0);
        stats.update(105.0);

        assert_eq!(stats.steps, 4);
        assert_eq!(stats.v_max, 110.0);
        assert_eq!(stats.v_min, 90.0);
        assert!((stats.v_mean - 101.25).abs() < 1e-10);
        assert!(stats.sigma_v > 0.0);
    }

    #[test]
    fn test_parameter_calculation_lower_bound() {
        let params = GamdParameters::new(
            SearchMode::GamdSearch,
            BoostForm::TotalBoost,
            ThresholdType::LowerBound,
            1.0,
            1.0,
        );

        let mut stats = GamdStatistics::new();
        for i in 0..100 {
            stats.update(-100.0 + i as f64);
        }

        if let Some((k0, k, e)) = params.calculate_parameters(&stats, 1.0) {
            assert!(k0 > 0.0 && k0 <= 1.0);
            assert!(k > 0.0);
            assert!((e - stats.v_max).abs() < 1e-10); // Lower bound: E = Vmax
        } else {
            panic!("Parameter calculation failed");
        }
    }

    #[test]
    fn test_parameter_calculation_upper_bound() {
        let params = GamdParameters::new(
            SearchMode::GamdSearch,
            BoostForm::TotalBoost,
            ThresholdType::UpperBound,
            1.0,
            1.0,
        );

        let mut stats = GamdStatistics::new();
        for i in 0..100 {
            stats.update(-100.0 + i as f64);
        }

        if let Some((k0, k, e)) = params.calculate_parameters(&stats, 1.0) {
            assert!(k0 > 0.0 && k0 <= 1.0);
            assert!(k > 0.0);
            assert!(e >= stats.v_max); // Upper bound: E >= Vmax (or E = Vmax if fallback)
        } else {
            panic!("Parameter calculation failed");
        }
    }

    #[test]
    fn test_boost_potential() {
        let mut params = GamdParameters::new(
            SearchMode::GamdSearch,
            BoostForm::TotalBoost,
            ThresholdType::LowerBound,
            1.0,
            1.0,
        );

        // Set up parameters manually
        params.k_tot = 0.01;
        params.e_tot = 100.0;

        // Create dummy configuration
        let mut config = Configuration::new(10, 3, 1);
        for i in 0..10 {
            config.current_mut().force[i] = Vec3::new(1.0, 0.0, 0.0);
        }

        let dih_forces = vec![Vec3::ZERO; 10];
        let tot_forces = vec![Vec3::new(1.0, 0.0, 0.0); 10];

        // Apply boost with V = 90 < E = 100
        let boost = params.apply_boost(&mut config, 0.0, &dih_forces, 90.0, &tot_forces);

        // Check boost potential: ΔV = 0.5 * k * (E - V)^2
        let expected = 0.5 * 0.01 * (100.0 - 90.0) * (100.0 - 90.0);
        assert!((boost - expected).abs() < 1e-10);

        // Check forces were modified
        let prefactor = 0.01 * (100.0 - 90.0);
        let expected_force = 1.0 + 1.0 * prefactor;
        assert!((config.current().force[0].x - expected_force).abs() < 1e-5);
    }

    #[test]
    fn test_cmd_search_no_boost() {
        let mut params = GamdParameters::new(
            SearchMode::CmdSearch,
            BoostForm::TotalBoost,
            ThresholdType::LowerBound,
            1.0,
            1.0,
        );

        params.k_tot = 0.01;
        params.e_tot = 100.0;

        let mut config = Configuration::new(10, 3, 1);
        let dih_forces = vec![Vec3::ZERO; 10];
        let tot_forces = vec![Vec3::new(1.0, 0.0, 0.0); 10];

        // CMD search should not apply boost
        let boost = params.apply_boost(&mut config, 0.0, &dih_forces, 90.0, &tot_forces);
        assert_eq!(boost, 0.0);
    }
}
