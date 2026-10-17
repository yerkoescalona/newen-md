//! Pressure tensor calculation algorithm.
//!
//! Equivalent to GROMOS `algorithm::Pressure_Calculation`.
//! Computes the pressure tensor from the kinetic energy tensor and virial tensor:
//!   P_tensor = (KE_tensor + 0.5 * virial_tensor) * (2.0 / V)
//!
//! The virial is stored as the raw outer product r⊗f (without -0.5 prefactor).
//! The 0.5 factor is applied here.
//!
//! KE tensor and virial correction are now computed inside the Forcefield
//! (prepare_virial + atomic_to_molecular_virial), matching GROMOS.
//! After exchange_state, they are in conf.old().
//!
//! Source: md++/src/algorithm/pressure/pressure_calculation.cc

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::math::Mat3;
use newen_core::topology::Topology;

/// Virial type: how the virial and KE tensor are computed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VirialType {
    None,
    Atomic,
    Molecular,
}

/// Calculates the pressure tensor from KE tensor and virial tensor.
///
/// GROMOS sequence position: after TemperatureCalculation, before BerendsenBarostat.
/// Reads from conf.old() (after exchange_state).
///
/// The KE tensor and molecular virial correction have already been computed
/// by the Forcefield algorithm (on conf.current(), which is now conf.old()).
///
/// Formula: P_tensor = (KE_tensor + 0.5 * virial_tensor) * (2.0 / V)
#[derive(Debug, Clone)]
pub struct PressureCalculation {
    pub virial_type: VirialType,
}

impl PressureCalculation {
    pub fn new(virial_type: VirialType) -> Self {
        Self { virial_type }
    }
}

impl Algorithm for PressureCalculation {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        if self.virial_type == VirialType::None {
            return Ok(());
        }

        let volume = conf.old().box_config.volume();

        if volume < 1e-30 {
            return Ok(());
        }

        // KE tensor and virial tensor are already in conf.old()
        // (computed by Forcefield on current(), then swapped by exchange_state)
        let ke_tensor = conf.old().kinetic_energy_tensor;
        let virial_tensor = conf.old().virial_tensor;

        // P_tensor = (KE_tensor + 0.5 * virial_tensor) * (2.0 / V)
        // tot_cg_factor = 1.0 for non-coarse-grained systems
        let factor = 2.0 / volume;
        let half_virial = Mat3::from_cols(
            virial_tensor.x_axis * 0.5,
            virial_tensor.y_axis * 0.5,
            virial_tensor.z_axis * 0.5,
        );
        let sum = Mat3::from_cols(
            ke_tensor.x_axis + half_virial.x_axis,
            ke_tensor.y_axis + half_virial.y_axis,
            ke_tensor.z_axis + half_virial.z_axis,
        );
        let p_tensor = Mat3::from_cols(
            sum.x_axis * factor,
            sum.y_axis * factor,
            sum.z_axis * factor,
        );

        // Store pressure tensor in old() (GROMOS convention)
        conf.old_mut().pressure_tensor = p_tensor;

        Ok(())
    }

    fn name(&self) -> &str {
        "Pressure_Calculation"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::configuration::Box as SimBox;
    use newen_core::math::Vec3;

    #[test]
    fn pressure_is_two_over_volume_times_kinetic_plus_half_virial() {
        let topo = Topology::new();
        let mut conf = Configuration::new(1, 1, 1);
        conf.old_mut().box_config = SimBox::rectangular(2.0, 2.0, 2.0); // V = 8
        conf.old_mut().kinetic_energy_tensor = Mat3::from_cols(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.0, 0.0, 3.0),
        );
        conf.old_mut().virial_tensor = Mat3::from_cols(
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, -2.0),
        );
        let sim = SimulationState::new(0.002, 1);
        PressureCalculation::new(VirialType::Atomic)
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        let p = conf.old().pressure_tensor;
        // P = 2/V (E_kin + ½ W): xx = 2/8·(1+2) = 0.75, yy = 2/8·2 = 0.5, zz = 2/8·(3−1) = 0.5
        assert!((p.x_axis.x - 0.75).abs() < 1e-14);
        assert!((p.y_axis.y - 0.5).abs() < 1e-14);
        assert!((p.z_axis.z - 0.5).abs() < 1e-14);
    }

    #[test]
    fn no_virial_means_no_pressure_update() {
        let topo = Topology::new();
        let mut conf = Configuration::new(1, 1, 1);
        conf.old_mut().box_config = SimBox::rectangular(2.0, 2.0, 2.0);
        conf.old_mut().kinetic_energy_tensor = Mat3::IDENTITY;
        let before = conf.old().pressure_tensor;
        let sim = SimulationState::new(0.002, 1);
        PressureCalculation::new(VirialType::None)
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        assert_eq!(conf.old().pressure_tensor, before);
    }
}
