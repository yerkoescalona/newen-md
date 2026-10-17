//! Leap-Frog velocity and position algorithms.
//!
//! Direct translation of:
//! - md++/src/algorithm/integration/leap_frog.cc (Leap_Frog_Velocity, Leap_Frog_Position)

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::topology::Topology;
use rayon::prelude::*;

/// Leap-Frog velocity update algorithm.
///
/// Equivalent to GROMOS `Leap_Frog_Velocity::apply()`:
/// 1. `exchange_state()` - swap current/old pointers
/// 2. `v_current = v_old + f_old * dt / m`
///
/// After apply(), `conf.current().vel` contains v(t+dt/2),
/// and `conf.old()` contains the previous state (positions at t, forces at t).
#[derive(Debug, Clone)]
pub struct LeapFrogVelocity {
    pub parallel: bool,
}

impl LeapFrogVelocity {
    pub fn new() -> Self {
        Self { parallel: false }
    }
}

impl Default for LeapFrogVelocity {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for LeapFrogVelocity {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        // GROMOS: exchange_state() first, then velocity update
        conf.exchange_state();
        // Copy the box to current state (GROMOS does this)
        conf.current_mut().box_config = conf.old().box_config.clone();

        let n_atoms = topo.inverse_mass.len();
        let dt = sim.dt;

        if self.parallel {
            let old_vel = conf.old().vel.clone();
            let old_force = conf.old().force.clone();

            conf.current_mut()
                .vel
                .par_iter_mut()
                .enumerate()
                .for_each(|(i, vel_new)| {
                    let accel = old_force[i] * topo.inverse_mass[i];
                    *vel_new = old_vel[i] + accel * dt;
                });
        } else {
            for i in 0..n_atoms {
                let accel = conf.old().force[i] * topo.inverse_mass[i];
                conf.current_mut().vel[i] = conf.old().vel[i] + accel * dt;
            }
        }

        let max_f = conf
            .old()
            .force
            .iter()
            .map(|f| f.length())
            .fold(0.0_f64, f64::max);
        let max_v = conf
            .current()
            .vel
            .iter()
            .map(|v| v.length())
            .fold(0.0_f64, f64::max);
        log::debug!("  max|f_old|={:.6e}  max|v_new|={:.6e}", max_f, max_v);

        Ok(())
    }

    fn name(&self) -> &str {
        "Leap_Frog_Velocity"
    }
}

/// Leap-Frog position update algorithm.
///
/// Equivalent to GROMOS `Leap_Frog_Position::apply()`:
/// `r_current = r_old + v_current * dt`
///
/// Must be called AFTER `LeapFrogVelocity` (which sets up current velocities).
#[derive(Debug, Clone)]
pub struct LeapFrogPosition {
    pub parallel: bool,
}

impl LeapFrogPosition {
    pub fn new() -> Self {
        Self { parallel: false }
    }
}

impl Default for LeapFrogPosition {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for LeapFrogPosition {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let n_atoms = topo.inverse_mass.len();
        let dt = sim.dt;

        if self.parallel {
            let old_pos = conf.old().pos.clone();
            let cur_vel = conf.current().vel.clone();

            conf.current_mut()
                .pos
                .par_iter_mut()
                .enumerate()
                .for_each(|(i, pos_new)| {
                    *pos_new = old_pos[i] + cur_vel[i] * dt;
                });
        } else {
            for i in 0..n_atoms {
                conf.current_mut().pos[i] = conf.old().pos[i] + conf.current().vel[i] * dt;
            }
        }

        let max_dr = (0..n_atoms)
            .map(|i| (conf.current().pos[i] - conf.old().pos[i]).length())
            .fold(0.0_f64, f64::max);
        log::debug!("  max|dr|={:.6e}", max_dr);

        Ok(())
    }

    fn name(&self) -> &str {
        "Leap_Frog_Position"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::Vec3;

    /// Helper to create a simple 2-atom topology (like pair_lj)
    fn make_pair_topology(mass: f64) -> Topology {
        let mut topo = Topology::new();
        // Add 2 atoms with given mass
        topo.mass = vec![mass, mass];
        topo.inverse_mass = vec![1.0 / mass, 1.0 / mass];
        topo.charge = vec![0.0, 0.0];
        topo.iac = vec![0, 0];
        topo
    }

    /// Helper to create a configuration for 2 atoms
    fn make_pair_configuration(pos: [Vec3; 2], vel: [Vec3; 2], force: [Vec3; 2]) -> Configuration {
        let mut conf = Configuration::new(2, 1, 1);
        // Set initial state in current
        conf.current_mut().pos = pos.to_vec();
        conf.current_mut().vel = vel.to_vec();
        conf.current_mut().force = force.to_vec();
        // Copy to old so exchange_state finds the right data
        conf.copy_current_to_old();
        conf
    }

    #[test]
    fn test_leapfrog_velocity_from_rest() {
        // Two argon atoms at rest with a force pushing them apart
        // F = [-41.486, 0, 0] on atom 0, [41.486, 0, 0] on atom 1
        // mass = 39.948 amu
        // dt = 0.002 ps
        // Expected: v = F*dt/m = [-41.486 * 0.002 / 39.948, 0, 0]
        //         = [-0.002077, 0, 0] for atom 0

        let mass = 39.948;
        let dt = 0.002;
        let force_x = -41.486240407;

        let topo = make_pair_topology(mass);
        let mut conf = make_pair_configuration(
            [Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.35, 0.0, 0.0)],
            [Vec3::ZERO, Vec3::ZERO],
            [Vec3::new(force_x, 0.0, 0.0), Vec3::new(-force_x, 0.0, 0.0)],
        );
        let sim = SimulationState::new(dt, 10);

        let mut vel_step = LeapFrogVelocity::new();
        vel_step.apply(&topo, &mut conf, &sim).unwrap();

        // After velocity step, current has new velocities
        let expected_vx = force_x * dt / mass; // = -0.002077...
        let v0 = conf.current().vel[0];
        let v1 = conf.current().vel[1];

        assert!(
            (v0.x - expected_vx).abs() < 1e-10,
            "v0.x = {}, expected {}",
            v0.x,
            expected_vx
        );
        assert!(
            (v1.x - (-expected_vx)).abs() < 1e-10,
            "v1.x = {}, expected {}",
            v1.x,
            -expected_vx
        );
        assert_eq!(v0.y, 0.0);
        assert_eq!(v0.z, 0.0);
    }

    #[test]
    fn test_leapfrog_position_update() {
        // After velocity step: v = [-0.002077, 0, 0] for atom 0
        // Position step: r_new = r_old + v*dt
        // r0_new = 0 + (-0.002077) * 0.002 = -4.154e-6

        let mass = 39.948;
        let dt = 0.002;
        let force_x = -41.486240407;
        let vx = force_x * dt / mass;

        let topo = make_pair_topology(mass);

        // Set up state as it would be AFTER velocity_step:
        // - current has new velocities
        // - old has old positions
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().vel = vec![Vec3::new(vx, 0.0, 0.0), Vec3::new(-vx, 0.0, 0.0)];
        conf.current_mut().pos = vec![Vec3::ZERO, Vec3::new(0.35, 0.0, 0.0)]; // will be overwritten

        // old state has the original positions
        // We need to manually set this up since position_step reads from old().pos
        // After exchange_state in velocity_step, old() has the original positions
        // Let's simulate that by setting state properly:
        // In the dual-state model, after velocity_step's exchange:
        //   current_idx was swapped, so old() = original state
        // We'll set the other state:
        conf.copy_current_to_old(); // wrong - we need old to have original pos
                                    // Actually let's just use the full sequence
        let mut conf2 = make_pair_configuration(
            [Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.35, 0.0, 0.0)],
            [Vec3::ZERO, Vec3::ZERO],
            [Vec3::new(force_x, 0.0, 0.0), Vec3::new(-force_x, 0.0, 0.0)],
        );
        let sim = SimulationState::new(dt, 10);

        // Run velocity step first
        let mut vel_step = LeapFrogVelocity::new();
        vel_step.apply(&topo, &mut conf2, &sim).unwrap();

        // Now run position step
        let mut pos_step = LeapFrogPosition::new();
        pos_step.apply(&topo, &mut conf2, &sim).unwrap();

        // Check new positions
        let r0 = conf2.current().pos[0];
        let r1 = conf2.current().pos[1];

        let expected_r0x = 0.0 + vx * dt;
        let expected_r1x = 0.35 + (-vx) * dt;

        assert!(
            (r0.x - expected_r0x).abs() < 1e-12,
            "r0.x = {}, expected {}",
            r0.x,
            expected_r0x
        );
        assert!(
            (r1.x - expected_r1x).abs() < 1e-12,
            "r1.x = {}, expected {}",
            r1.x,
            expected_r1x
        );
    }

    #[test]
    fn test_leapfrog_kinetic_energy_after_velocity_step() {
        // GROMOS convention for E_kin:
        // E_kin = 0.5 * sum_i(m_i * (|v_new_i|^2 + |v_old_i|^2) / 2)
        // This averages the kinetic energy between old and new velocities.
        //
        // At step 0: v_old = 0, v_new = F*dt/m
        // E_kin = 0.5 * sum_i(m_i * |v_new_i|^2 / 2) = quarter of naive sum
        //
        // Reference (pair_lj step 0): E_kin = 8.616742481e-05 kJ/mol

        let mass = 39.948;
        let dt = 0.002;
        let force_x = -41.486240407;

        let topo = make_pair_topology(mass);
        let mut conf = make_pair_configuration(
            [Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.35, 0.0, 0.0)],
            [Vec3::ZERO, Vec3::ZERO],
            [Vec3::new(force_x, 0.0, 0.0), Vec3::new(-force_x, 0.0, 0.0)],
        );
        let sim = SimulationState::new(dt, 10);

        let mut vel_step = LeapFrogVelocity::new();
        vel_step.apply(&topo, &mut conf, &sim).unwrap();

        // GROMOS formula: E_kin = 0.5 * sum(m_i * (|v_new|^2 + |v_old|^2) / 2)
        // After exchange_state in velocity_step, old() has original v=0
        let mut e_kin = 0.0;
        for i in 0..2 {
            let v_new = conf.current().vel[i];
            let v_old = conf.old().vel[i]; // should be zero at step 0
            e_kin += 0.5 * mass * (v_new.length_squared() + v_old.length_squared()) / 2.0;
        }

        let expected_ekin = 8.616742481e-05;
        let rel_error = (e_kin - expected_ekin).abs() / expected_ekin;
        assert!(
            rel_error < 1e-6,
            "E_kin = {:.10e}, expected {:.10e}, rel_error = {:.2e}",
            e_kin,
            expected_ekin,
            rel_error
        );
    }
}
