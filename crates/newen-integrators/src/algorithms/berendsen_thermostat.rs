//! Berendsen weak-coupling thermostat algorithm.
//!
//! Equivalent to GROMOS `algorithm::Berendsen_Thermostat`.
//! Rescales velocities to weakly couple the system to a heat bath.
//!
//! Placed between LeapFrogVelocity and LeapFrogPosition in the sequence.
//! Uses the "new" kinetic energy from the previous step's TemperatureCalculation
//! (stored in energies.kinetic_energy_new, equivalent to GROMOS multibath.bath.ekin).
//!
//! Source: md++/src/algorithm/temperature/berendsen_thermostat.cc
//!         md++/src/algorithm/temperature/thermostat.cc (scale method)

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::topology::Topology;
use newen_core::units::kB;

use super::multibath::{BathRange, ZERO_TEMPERATURE};

/// The Berendsen scaling factor of one bath (gromosXX `Berendsen_Thermostat::calc_scaling`):
///
/// ```text
///   T_free = 2·E_kin / (N_df·k_B)        (0 when N_df ≤ 0; T₀ when below ZERO_TEMPERATURE)
///   λ = sqrt(1 + dt/τ · (T₀/T_free − 1))   τ > 0
///   λ = sqrt(T₀/T_free)                    τ = 0, instantaneous
///   λ = 1                                  τ < 0 (no coupling), or T_free = T₀ = 0
/// ```
///
/// spec: thermostat.berendsen_scale
pub fn berendsen_scale(ekin: f64, dof: f64, temperature: f64, tau: f64, dt: f64) -> f64 {
    if tau < 0.0 {
        return 1.0;
    }
    let mut free_temp = if dof > 0.0 {
        2.0 * ekin / (dof * kB)
    } else {
        0.0
    };
    if free_temp < ZERO_TEMPERATURE {
        free_temp = temperature;
    }
    if free_temp < ZERO_TEMPERATURE {
        1.0
    } else if tau == 0.0 {
        (temperature / free_temp).sqrt()
    } else {
        (1.0 + dt / tau * (temperature / free_temp - 1.0)).sqrt()
    }
}

/// Berendsen thermostat parameters for a single temperature bath.
#[derive(Debug, Clone)]
pub struct BerendsenThermostatParams {
    /// Reference temperature T0 (K)
    pub temperature: f64,
    /// Coupling time constant τ (ps). τ < 0 → no coupling, τ = 0 → instantaneous
    pub tau: f64,
    /// Degrees of freedom for this bath
    pub dof: f64,
}

/// Berendsen weak-coupling thermostat.
///
/// GROMOS sequence position: after Leap_Frog_Velocity, before Leap_Frog_Position.
///
/// The scaling factor is:
///   λ = sqrt(1 + dt/τ · (T₀/T_free - 1))
/// where T_free = 2·E_kin / (DOF · k_B) using E_kin from the previous step.
#[derive(Debug, Clone)]
pub struct BerendsenThermostat {
    pub params: Vec<BerendsenThermostatParams>,
    /// The DOFSET ranges that map atoms to `params` (gromosXX `bath_index`), in file order.
    pub ranges: Vec<BathRange>,
    /// Scaling factor per bath of the current step.
    scale: Vec<f64>,
}

impl BerendsenThermostat {
    /// Create a Berendsen thermostat from resolved baths: reference temperature, coupling time and
    /// degrees of freedom per bath, and the DOFSET ranges that assign the atoms to them.
    pub fn new(params: Vec<BerendsenThermostatParams>, ranges: Vec<BathRange>) -> Self {
        Self {
            scale: vec![1.0; params.len()],
            params,
            ranges,
        }
    }

    /// The scaling factor of every bath at the last step (gromosXX `bath_struct::scale`).
    pub fn scales(&self) -> &[f64] {
        &self.scale
    }

    /// Create a new Berendsen thermostat with a single bath covering all atoms.
    pub fn new_single_bath(temperature: f64, tau: f64, dof: f64, n_atoms: usize) -> Self {
        Self {
            params: vec![BerendsenThermostatParams {
                temperature,
                tau,
                dof,
            }],
            ranges: vec![BathRange {
                last_atom: n_atoms - 1,
                com_bath: 0,
                ir_bath: 0,
            }],
            scale: vec![1.0],
        }
    }
}

impl Algorithm for BerendsenThermostat {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        // After exchange_state in LeapFrogVelocity:
        //   current().energies.kinetic_energy_new = E_kin_new from previous step
        // (set by TemperatureCalculation into old().energies, then swapped by exchange_state)

        // gromosXX `Berendsen_Thermostat::calc_scaling`: one factor per bath
        self.scale.resize(self.params.len(), 1.0);
        for (bath_idx, params) in self.params.iter().enumerate() {
            // Read E_kin from the previous step (stored by TemperatureCalculation). A multi-bath
            // run scales each bath on its own kinetic energy, as gromosXX does
            // (`bath_struct::ekin`); with one bath the per-bath vector is empty and the total is
            // the same number.
            let ekin = conf
                .current()
                .energies
                .kinetic_energy_new_bath
                .get(bath_idx)
                .copied()
                .unwrap_or(conf.current().energies.kinetic_energy_new);

            let scale = berendsen_scale(ekin, params.dof, params.temperature, params.tau, sim.dt);

            log::debug!(
                "  Berendsen bath {}: E_kin_new={:.6e}, T0={:.2}, scale={:.10}",
                bath_idx,
                ekin,
                params.temperature,
                scale
            );

            self.scale[bath_idx] = scale;
        }

        // gromosXX `Thermostat::scale` (`thermostat.cc:45-136`): every DOFSET range is scaled
        // by its bath. A COM bath different from the IR bath is refused at plan validation.
        // gromosXX's `continue` on an unscaled range also skips its `last_atom = r_it->last_atom`,
        // so the next range would scale this range's atoms too; that is not reproduced.
        let mut first = 0;
        for r in &self.ranges {
            let scale = self.scale[r.com_bath];
            if (scale - 1.0).abs() >= 1e-15 {
                for v in &mut conf.current_mut().vel[first..=r.last_atom] {
                    *v *= scale;
                }
            }
            first = r.last_atom + 1;
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "Berendsen_Thermostat"
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::Vec3;

    /// One atom, three degrees of freedom, τ = Δt: one step scales the velocities so that the
    /// new kinetic energy is exactly the target temperature's.
    #[test]
    fn tau_equal_to_dt_reaches_the_target_temperature_in_one_step() {
        let topo = Topology::new();
        let mut conf = Configuration::new(1, 1, 1);
        conf.current_mut().vel[0] = Vec3::new(1.0, 2.0, 2.0); // |v|² = 9
        let mass = 4.0;
        let ekin = 0.5 * mass * 9.0;
        conf.current_mut().energies.kinetic_energy_new = ekin;
        let dof = 3.0;
        let t_free = 2.0 * ekin / (dof * kB);
        let t0 = 300.0;
        let dt = 0.002;
        let sim = SimulationState::new(dt, 1);
        let mut bath = BerendsenThermostat::new_single_bath(t0, dt, dof, 1);
        bath.apply(&topo, &mut conf, &sim).unwrap();
        let scale = (t0 / t_free).sqrt();
        let v = conf.current().vel[0];
        assert!((v.x - scale).abs() < 1e-12 && (v.y - 2.0 * scale).abs() < 1e-12);
        let ekin_after = 0.5 * mass * v.length_squared();
        assert!((2.0 * ekin_after / (dof * kB) - t0).abs() < 1e-9);
    }

    #[test]
    fn negative_tau_disables_the_bath() {
        let topo = Topology::new();
        let mut conf = Configuration::new(1, 1, 1);
        conf.current_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().energies.kinetic_energy_new = 0.5;
        let sim = SimulationState::new(0.002, 1);
        let mut bath = BerendsenThermostat::new_single_bath(300.0, -1.0, 3.0, 1);
        bath.apply(&topo, &mut conf, &sim).unwrap();
        assert_eq!(conf.current().vel[0], Vec3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn each_dofset_range_is_scaled_by_its_own_bath() {
        // DOFSET `1 2 2  2 1 1`: atom 0 couples to bath 2, atom 1 to bath 1 (only bath 2 coupled)
        let topo = Topology::new();
        let mut conf = Configuration::new(2, 2, 1);
        conf.current_mut().vel = vec![Vec3::new(1.0, 0.0, 0.0); 2];
        conf.current_mut().energies.kinetic_energy_new_bath = vec![0.5, 0.5];
        let sim = SimulationState::new(0.002, 1);
        let bath = |tau| BerendsenThermostatParams {
            temperature: 300.0,
            tau,
            dof: 3.0,
        };
        let ranges = vec![
            BathRange {
                last_atom: 0,
                com_bath: 1,
                ir_bath: 1,
            },
            BathRange {
                last_atom: 1,
                com_bath: 0,
                ir_bath: 0,
            },
        ];
        let mut t = BerendsenThermostat::new(vec![bath(-1.0), bath(0.002)], ranges);
        t.apply(&topo, &mut conf, &sim).unwrap();
        let v = &conf.current().vel;
        assert_ne!(v[0].x, 1.0);
        assert_eq!(v[1].x, 1.0);
    }

    #[test]
    fn a_system_already_at_the_target_is_left_alone() {
        let topo = Topology::new();
        let mut conf = Configuration::new(1, 1, 1);
        let dof = 3.0;
        let t0 = 300.0;
        conf.current_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().energies.kinetic_energy_new = 0.5 * dof * kB * t0;
        let sim = SimulationState::new(0.002, 1);
        let mut bath = BerendsenThermostat::new_single_bath(t0, 0.1, dof, 1);
        bath.apply(&topo, &mut conf, &sim).unwrap();
        assert!((conf.current().vel[0].x - 1.0).abs() < 1e-12);
    }
}
