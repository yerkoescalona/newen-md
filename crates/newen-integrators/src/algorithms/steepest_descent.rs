//! Steepest Descent energy minimization algorithm.
//!
//! Wraps the SteepestDescent integrator as an Algorithm for the AlgorithmSequence.
//! This replaces both LeapFrogVelocity and LeapFrogPosition in the sequence.
//!
//! In GROMOS, steepest descent is an integrator that replaces the leap-frog
//! steps. The rest of the sequence (Forcefield, EnergyCalculation) stays the same.

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::math::Vec3;
use newen_core::topology::Topology;

/// Steepest Descent Algorithm for AlgorithmSequence.
///
/// Performs energy minimization instead of MD integration.
/// Replaces LeapFrogVelocity + LeapFrogPosition in the sequence.
#[derive(Debug, Clone)]
pub struct SteepestDescentAlgorithm {
    /// Current step size (nm)
    pub step_size: f64,
    /// Initial step size (dx0)
    pub initial_step_size: f64,
    /// Maximum step size (dxm)
    pub max_step_size: f64,
    /// Energy convergence tolerance (dele) in kJ/mol
    pub energy_tolerance: f64,
    /// Minimum steps before convergence check (nmin)
    pub min_steps: usize,
    /// Maximum force magnitude limit (flim, 0.0 = unlimited)
    pub force_limit: f64,
    /// Whether minimization has converged
    converged: bool,
}

impl SteepestDescentAlgorithm {
    /// Create with GROMOS default parameters
    pub fn new() -> Self {
        Self {
            step_size: 0.01,
            initial_step_size: 0.01,
            max_step_size: 0.05,
            energy_tolerance: 0.1,
            min_steps: 1,
            force_limit: 0.0,
            converged: false,
        }
    }

    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.energy_tolerance = tolerance;
        self
    }

    pub fn with_step_sizes(mut self, initial: f64, max: f64) -> Self {
        self.initial_step_size = initial;
        self.step_size = initial;
        self.max_step_size = max;
        self
    }

    pub fn with_force_limit(mut self, limit: f64) -> Self {
        self.force_limit = limit;
        self
    }

    pub fn with_min_steps(mut self, steps: usize) -> Self {
        self.min_steps = steps;
        self
    }

    /// Check if minimization has converged
    pub fn is_converged(&self) -> bool {
        self.converged
    }

    /// Apply force limiting if configured
    fn apply_force_limit(&self, forces: &mut [Vec3]) {
        if self.force_limit > 0.0 {
            for force in forces.iter_mut() {
                let magnitude = force.length();
                if magnitude > self.force_limit {
                    *force *= self.force_limit / magnitude;
                }
            }
        }
    }

    /// Calculate 1/|F| for normalized step direction
    fn calculate_force_norm(&self, forces: &[Vec3]) -> f64 {
        let f_squared: f64 = forces
            .iter()
            .map(|f| f.x * f.x + f.y * f.y + f.z * f.z)
            .sum();

        if f_squared < 1e-15 {
            1.0
        } else {
            1.0 / f_squared.sqrt()
        }
    }
}

impl Default for SteepestDescentAlgorithm {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for SteepestDescentAlgorithm {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        if self.converged {
            return Ok(());
        }

        let n_atoms = topo.num_atoms();

        // GROMOS: only check convergence if past nmin steps
        // Uses sim.steps() (the global step counter), not an internal counter
        if sim.step > self.min_steps {
            // GROMOS uses potential_total + special_total for convergence
            let ecur =
                conf.current().energies.potential_total + conf.current().energies.special_total;
            let eold = conf.old().energies.potential_total + conf.old().energies.special_total;

            let energy_change = (ecur - eold).abs();

            log::debug!(
                "SD step {}: E_pot={:.6e}, E_old={:.6e}, dE={:.6e}, dele={:.6e}",
                sim.step,
                ecur,
                eold,
                energy_change,
                self.energy_tolerance
            );

            if energy_change < self.energy_tolerance {
                self.converged = true;
                log::info!(
                    "Steepest descent converged at step {}: dE = {:.6e} kJ/mol",
                    sim.step,
                    energy_change
                );
                return Ok(());
            }

            // Adaptive step sizing (GROMOS: *1.2 if decrease, *0.5 if increase)
            if ecur < eold {
                self.step_size = (self.step_size * 1.2).min(self.max_step_size);
            } else {
                self.step_size *= 0.5;
            }
        } else {
            // GROMOS: reset step size to dx0 for initial steps
            self.step_size = self.initial_step_size;
        }

        // Apply force limiting (GROMOS: per-atom |f| > flim → scale down)
        self.apply_force_limit(&mut conf.current_mut().force);

        // Calculate <f|f>^-0.5 (normalized force magnitude)
        let f_norm = self.calculate_force_norm(&conf.current().force);

        // NaN check (GROMOS: gisnan)
        if f_norm.is_nan() {
            return Err("Force is NaN during steepest descent".to_string());
        }

        // exchange_state: swap old <-> current
        conf.exchange_state();
        // Copy box to current (GROMOS convention)
        conf.current_mut().box_config = conf.old().box_config.clone();

        // Update positions: r_new = r_old + step_size * f_norm * force_old
        let step_scale = self.step_size * f_norm;

        for i in 0..n_atoms {
            conf.current_mut().pos[i] = conf.old().pos[i] + conf.old().force[i] * step_scale;
        }

        // Zero velocities in both states (GROMOS: conf.old().vel = 0; conf.current().vel = 0)
        conf.old_mut().vel.fill(Vec3::ZERO);
        conf.current_mut().vel.fill(Vec3::ZERO);
        conf.current_mut().clear_forces();
        conf.current_mut().vel.fill(Vec3::ZERO);
        conf.current_mut().clear_forces();

        Ok(())
    }

    fn name(&self) -> &str {
        "Steepest-Descent"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::topology::MolTypeAtom;

    /// A topology of free atoms with the given masses (one solute molecule, no bonds;
    /// `Topology::new()` already registers the solute instance).
    fn bare(masses: &[f64]) -> Topology {
        let mut topo = Topology::new();
        topo.mass = masses.to_vec();
        topo.inverse_mass = masses.iter().map(|m| 1.0 / m).collect();
        topo.moltypes[0].atoms = masses
            .iter()
            .map(|m| MolTypeAtom {
                name: "X".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: *m,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            })
            .collect();
        topo
    }

    #[test]
    fn first_step_moves_along_the_normalised_force_by_the_initial_step_size() {
        let topo = bare(&[1.0, 1.0]);
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().force[0] = Vec3::new(3.0, 0.0, 0.0);
        conf.current_mut().force[1] = Vec3::new(0.0, 4.0, 0.0); // |F| = 5
        let sim = SimulationState::new(0.002, 10); // step 0 ≤ min_steps
        let mut sd = SteepestDescentAlgorithm::new().with_step_sizes(0.01, 0.05);
        sd.apply(&topo, &mut conf, &sim).unwrap();
        // displacement = step_size · F / |F|
        assert!((conf.current().pos[0].x - 0.01 * 3.0 / 5.0).abs() < 1e-15);
        assert!((conf.current().pos[1].y - 0.01 * 4.0 / 5.0).abs() < 1e-15);
        assert_eq!(conf.current().vel[0], Vec3::ZERO);
        assert!(!sd.is_converged());
    }

    #[test]
    fn converges_when_the_energy_change_drops_below_the_tolerance() {
        let topo = bare(&[1.0]);
        let mut conf = Configuration::new(1, 1, 1);
        conf.current_mut().energies.potential_total = -10.0;
        conf.old_mut().energies.potential_total = -10.00001;
        let mut sim = SimulationState::new(0.002, 10);
        sim.step = 5;
        let mut sd = SteepestDescentAlgorithm::new().with_tolerance(0.001);
        sd.apply(&topo, &mut conf, &sim).unwrap();
        assert!(sd.is_converged());
    }

    #[test]
    fn a_nan_force_is_reported_as_an_error() {
        let topo = bare(&[1.0]);
        let mut conf = Configuration::new(1, 1, 1);
        conf.current_mut().force[0] = Vec3::new(f64::NAN, 0.0, 0.0);
        let sim = SimulationState::new(0.002, 10);
        assert!(SteepestDescentAlgorithm::new()
            .apply(&topo, &mut conf, &sim)
            .is_err());
    }
}
