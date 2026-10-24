//! LINCS constraint algorithm wrapper.
//!
//! Equivalent to GROMOS `algorithm::Lincs` — linear constraint solver
//! (Hess et al. 1997), selectable for the solute (NTCP) and/or solvent (NTCS).
//!
//! Source: md++/src/algorithm/constraints/lincs.cc

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::topology::Topology;

use crate::constraints::{lincs_buffered, LincsBuffers, NtcMode};

/// LINCS constraint algorithm for the MD sequence.
///
/// `init()` precomputes the coupling matrices for the solute and/or solvent
/// (whichever is selected via `include_solute`/`include_solvent`), `apply()`
/// runs LINCS for both groups every step (GROMOS `Lincs::apply`,
/// `lincs.cc:238-284`).
pub struct LincsAlgorithm {
    ntc: NtcMode,
    solute_order: usize,
    solvent_order: usize,
    include_solute: bool,
    include_solvent: bool,
    buffers: Option<LincsBuffers>,
}

impl LincsAlgorithm {
    pub fn new(
        ntc: NtcMode,
        solute_order: usize,
        solvent_order: usize,
        include_solute: bool,
        include_solvent: bool,
    ) -> Self {
        Self {
            ntc,
            solute_order,
            solvent_order,
            include_solute,
            include_solvent,
            buffers: None,
        }
    }
}

impl Algorithm for LincsAlgorithm {
    fn init(
        &mut self,
        topo: &Topology,
        _conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        self.buffers = Some(LincsBuffers::new(
            topo,
            self.ntc,
            self.solute_order,
            self.solvent_order,
            self.include_solute,
            self.include_solvent,
        ));
        Ok(())
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        if let Some(ref buffers) = self.buffers {
            // GROMOS `Lincs::apply` (lincs.cc:238-288) always returns success —
            // LINCS is a fixed-order analytical expansion, not an iterative solver
            // that can fail to converge, so we don't gate on `max_error` here
            // (unlike SHAKE, whose iteration can genuinely fail).
            let _ = lincs_buffered(topo, conf, sim.dt, buffers);
            Ok(())
        } else {
            Err("LincsAlgorithm::apply called before init".to_string())
        }
    }

    fn name(&self) -> &str {
        "LincsConstraints"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::math::Vec3;
    use newen_core::topology::{Bond, BondParameters, MolTypeAtom};

    fn two_atom_bond_topology() -> Topology {
        let mut topo = Topology::new();
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0];
        topo.moltypes[0].atoms = vec![
            MolTypeAtom {
                name: "C1".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
            MolTypeAtom {
                name: "C2".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
        ];
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 0.0,
            r0: 0.15,
        });
        topo
    }

    /// Old positions at the constraint length, new ones stretched by 0.01 nm.
    fn stretched_pair() -> Configuration {
        let mut conf = Configuration::new(2, 1, 1);
        conf.old_mut().pos[1] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.16, 0.0, 0.0);
        conf
    }

    #[test]
    fn apply_restores_the_bond_length() {
        let topo = two_atom_bond_topology();
        let mut conf = stretched_pair();
        let sim = SimulationState::new(0.002, 10);
        let mut alg = LincsAlgorithm::new(NtcMode::AllBonds, 4, 4, true, false);
        alg.init(&topo, &mut conf, &sim).unwrap();
        alg.apply(&topo, &mut conf, &sim).unwrap();
        let d = (conf.current().pos[1] - conf.current().pos[0]).length();
        assert!((d - 0.15).abs() < 1e-6, "distance after LINCS {d}");
        assert_eq!(alg.name(), "LincsConstraints");
    }

    #[test]
    fn apply_before_init_is_an_error() {
        let topo = two_atom_bond_topology();
        let mut conf = stretched_pair();
        let sim = SimulationState::new(0.002, 10);
        let mut alg = LincsAlgorithm::new(NtcMode::AllBonds, 4, 4, true, false);
        assert!(alg.apply(&topo, &mut conf, &sim).is_err());
    }
}
