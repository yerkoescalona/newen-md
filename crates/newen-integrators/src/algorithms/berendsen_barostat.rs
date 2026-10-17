//! Berendsen barostat algorithm (weak pressure coupling).
//!
//! Equivalent to GROMOS `algorithm::Berendsen_Barostat`.
//! Isotropic scaling: mu = (1 - comp * dt / tau * (P0 - P))^(1/3)
//! Then scale box and positions by mu.
//!
//! GROMOS sequence position: after PressureCalculation, before EnergyCalculation.
//! Reads pressure_tensor from conf.old(), scales conf.current() positions and box.
//!
//! Source: md++/src/algorithm/pressure/berendsen_barostat.cc

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::topology::Topology;

/// Berendsen barostat parameters.
#[derive(Debug, Clone)]
pub struct BerendsenBarostatParams {
    /// Reference pressure P0 (kJ/(mol·nm³)) — GROMOS uses internal units
    pub pressure0: f64,
    /// Isothermal compressibility κ (kJ/(mol·nm³))⁻¹
    pub compressibility: f64,
    /// Coupling time constant τ_P (ps)
    pub tau: f64,
}

/// Berendsen weak-coupling barostat (isotropic).
///
/// Scaling factor: μ = (1 - κ · dt/τ · (P₀ - P))^(1/3)
/// Scale box: box *= μ
/// Scale positions: pos[i] *= μ
#[derive(Debug, Clone)]
pub struct BerendsenBarostat {
    pub params: BerendsenBarostatParams,
}

impl BerendsenBarostat {
    pub fn new(params: BerendsenBarostatParams) -> Self {
        Self { params }
    }
}

impl Algorithm for BerendsenBarostat {
    fn apply(
        &mut self,
        _topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        // Read pressure tensor from PressureCalculation (stored in old())
        let pressure = conf.old().pressure_tensor;

        // Isotropic pressure: P = (P_xx + P_yy + P_zz) / 3
        let total_pressure = (pressure.x_axis.x + pressure.y_axis.y + pressure.z_axis.z) / 3.0;

        // Scaling factor: mu = (1 - comp * dt / tau * (P0 - P))^(1/3)
        let mu = (1.0
            - self.params.compressibility * sim.dt / self.params.tau
                * (self.params.pressure0 - total_pressure))
            .powf(1.0 / 3.0);

        // Scale the box
        {
            let box_vecs = conf.current().box_config.vectors;
            let scaled = newen_core::math::Mat3::from_cols(
                box_vecs.x_axis * mu,
                box_vecs.y_axis * mu,
                box_vecs.z_axis * mu,
            );
            conf.current_mut().box_config.vectors = scaled;
            conf.current_mut().box_config.inv_vectors = scaled.inverse();
        }

        // Scale the positions
        let n = conf.current().pos.len();
        for i in 0..n {
            conf.current_mut().pos[i] *= mu;
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "Berendsen_Barostat"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::configuration::Box as SimBox;
    use newen_core::math::{Mat3, Vec3};

    fn system(pressure: f64) -> (Topology, Configuration) {
        let topo = Topology::new();
        let mut conf = Configuration::new(1, 1, 1);
        conf.current_mut().box_config = SimBox::rectangular(2.0, 2.0, 2.0);
        conf.current_mut().pos[0] = Vec3::new(1.0, 0.5, 0.25);
        conf.old_mut().pressure_tensor = Mat3::IDENTITY * pressure;
        (topo, conf)
    }

    #[test]
    fn at_the_reference_pressure_nothing_moves() {
        let (topo, mut conf) = system(0.06102);
        let sim = SimulationState::new(0.002, 1);
        let mut alg = BerendsenBarostat::new(BerendsenBarostatParams {
            pressure0: 0.06102,
            compressibility: 4.575e-4,
            tau: 0.5,
        });
        alg.apply(&topo, &mut conf, &sim).unwrap();
        assert!((conf.current().box_config.vectors.x_axis.x - 2.0).abs() < 1e-15);
        assert_eq!(conf.current().pos[0], Vec3::new(1.0, 0.5, 0.25));
    }

    #[test]
    fn above_the_reference_pressure_the_box_and_positions_grow_together() {
        let (topo, mut conf) = system(1.0);
        let sim = SimulationState::new(0.002, 1);
        let (p0, comp, tau) = (0.0, 4.575e-4, 0.5);
        let mut alg = BerendsenBarostat::new(BerendsenBarostatParams {
            pressure0: p0,
            compressibility: comp,
            tau,
        });
        alg.apply(&topo, &mut conf, &sim).unwrap();
        let mu = (1.0 - comp * 0.002 / tau * (p0 - 1.0)).powf(1.0 / 3.0);
        assert!(mu > 1.0);
        assert!((conf.current().box_config.vectors.x_axis.x - 2.0 * mu).abs() < 1e-15);
        assert!((conf.current().pos[0].x - mu).abs() < 1e-15);
        assert!((conf.current().pos[0].z - 0.25 * mu).abs() < 1e-15);
        // the inverse box is kept consistent
        let prod = conf.current().box_config.vectors * conf.current().box_config.inv_vectors;
        assert!((prod.x_axis.x - 1.0).abs() < 1e-12 && prod.x_axis.y.abs() < 1e-12);
    }
}
