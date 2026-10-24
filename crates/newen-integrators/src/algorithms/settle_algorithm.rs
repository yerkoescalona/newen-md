//! SETTLE constraint algorithm wrapper.
//!
//! Equivalent to GROMOS `algorithm::Settle` — analytical 3-site rigid water
//! constraint solver (Miyamoto & Kollman, 1992). Solvent-only.
//!
//! Source: md++/src/algorithm/constraints/settle.cc

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::topology::Topology;

use crate::constraints::settle;

/// SETTLE constraint algorithm for the MD sequence.
///
/// Applied after the position update to analytically enforce the rigid
/// 3-site water geometry. Equivalent to GROMOS's `Settle::apply`.
pub struct SettleAlgorithm {
    /// The run computes the atomic virial (`VIRIAL 1`): only then does SETTLE add its virial.
    pub atomic_virial: bool,
}

impl SettleAlgorithm {
    pub fn new() -> Self {
        Self {
            atomic_virial: false,
        }
    }
}

impl Default for SettleAlgorithm {
    fn default() -> Self {
        Self::new()
    }
}

impl Algorithm for SettleAlgorithm {
    fn init(
        &mut self,
        topo: &Topology,
        _conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        // GROMOS Settle::init validation (settle.cc:53-118)
        if topo.num_solvent_molecules() == 0 {
            return Err("SETTLE does only work if 1 solvent.".to_string());
        }
        if topo.atoms_per_solvent() != 3 {
            return Err("SETTLE does only work with water like molecules (3 atoms).".to_string());
        }
        let n_solute = topo.num_solute_atoms();
        if topo.mass[n_solute + 1] != topo.mass[n_solute + 2] {
            return Err(
                "SETTLE does only work with water like molecules (wrong masses).".to_string(),
            );
        }
        if topo.solvent_constraint_template.len() != 3 {
            return Err(
                "SETTLE does only work with water like molecules (3 distance constraints)."
                    .to_string(),
            );
        }
        if (topo.solvent_constraint_template[0].length - topo.solvent_constraint_template[1].length)
            .abs()
            > 1e-12
        {
            return Err(
                "SETTLE does only work with water like molecules (distance constraints wrong)."
                    .to_string(),
            );
        }
        Ok(())
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let result = settle(topo, conf, sim.dt, self.atomic_virial);
        if result.converged {
            Ok(())
        } else {
            Err("SETTLE failed: sin(angle) > 1.0 in analytical solve".to_string())
        }
    }

    fn name(&self) -> &str {
        "SettleConstraints"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::Vec3;
    use newen_core::topology::{SolventAtomTemplate, SolventConstraintTemplate};

    const D_OH: f64 = 0.1;
    const D_HH: f64 = 0.1633;

    fn one_spc_water() -> Topology {
        let mut topo = Topology::new();
        topo.mass = vec![15.9994, 1.008, 1.008];
        topo.inverse_mass = topo.mass.iter().map(|m| 1.0 / m).collect();
        topo.solvent_atom_template = vec![
            SolventAtomTemplate {
                iac: 0,
                name: "OW".into(),
                mass: 15.9994,
                charge: -0.82,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW1".into(),
                mass: 1.008,
                charge: 0.41,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW2".into(),
                mass: 1.008,
                charge: 0.41,
            },
        ];
        topo.solvent_constraint_template = vec![
            SolventConstraintTemplate {
                i: 0,
                j: 1,
                length: D_OH,
            },
            SolventConstraintTemplate {
                i: 0,
                j: 2,
                length: D_OH,
            },
            SolventConstraintTemplate {
                i: 1,
                j: 2,
                length: D_HH,
            },
        ];
        topo.solvate(1);
        topo
    }

    /// Canonical SPC triangle as the old geometry, a small displacement as the new one.
    fn perturbed_water() -> Configuration {
        let mut conf = Configuration::new(3, 1, 1);
        let half = D_HH / 2.0;
        let height = (D_OH * D_OH - half * half).sqrt();
        conf.old_mut().pos[0] = Vec3::new(0.0, height * (2.0 / 3.0), 0.0);
        conf.old_mut().pos[1] = Vec3::new(-half, -height / 3.0, 0.0);
        conf.old_mut().pos[2] = Vec3::new(half, -height / 3.0, 0.0);
        conf.current_mut().pos[0] = conf.old().pos[0] + Vec3::new(0.001, 0.0005, -0.0007);
        conf.current_mut().pos[1] = conf.old().pos[1] + Vec3::new(-0.0012, 0.0009, 0.0006);
        conf.current_mut().pos[2] = conf.old().pos[2] + Vec3::new(0.0008, -0.0011, 0.0003);
        conf
    }

    #[test]
    fn apply_restores_the_water_geometry() {
        let topo = one_spc_water();
        let mut conf = perturbed_water();
        let sim = SimulationState::new(0.002, 10);
        let mut alg = SettleAlgorithm::new();
        alg.init(&topo, &mut conf, &sim).unwrap();
        alg.apply(&topo, &mut conf, &sim).unwrap();
        let p = &conf.current().pos;
        assert!(((p[1] - p[0]).length() - D_OH).abs() < 1e-10);
        assert!(((p[2] - p[0]).length() - D_OH).abs() < 1e-10);
        assert!(((p[2] - p[1]).length() - D_HH).abs() < 1e-10);
        assert_eq!(alg.name(), "SettleConstraints");
    }

    #[test]
    fn init_rejects_a_system_without_water() {
        let mut topo = Topology::new();
        topo.mass = vec![12.0];
        let mut conf = Configuration::new(1, 1, 1);
        let sim = SimulationState::new(0.002, 10);
        let err = SettleAlgorithm::new()
            .init(&topo, &mut conf, &sim)
            .unwrap_err();
        assert!(err.contains("1 solvent"), "{err}");
    }

    #[test]
    fn init_rejects_unequal_hydrogen_masses() {
        let mut topo = one_spc_water();
        topo.mass[2] = 2.014;
        let mut conf = perturbed_water();
        let sim = SimulationState::new(0.002, 10);
        let err = SettleAlgorithm::new()
            .init(&topo, &mut conf, &sim)
            .unwrap_err();
        assert!(err.contains("masses"), "{err}");
    }
}
