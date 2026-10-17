//! Put charge groups back into the periodic box, as gromosXX's `Lattice_Shift_Tracker` does.
//!
//! gromosXX runs this once per step, right after the centre-of-mass removal and before the force
//! field (`create_md_sequence.cc`: `Lattice_Shift_Tracker` for every non-vacuum boundary), so a
//! molecule that drifts out of the box is translated back as a whole — the charge group is the
//! unit, which is what keeps the charge-group cutoff and the written configuration meaningful.
//!
//! The shift is a whole box vector, so no interaction changes: every distance goes through the
//! minimum-image convention either way. What it does change is the coordinates that get written,
//! and which image of a molecule the charge-group cutoff sees.
//!
//! Source: md++/src/math/periodicity.cc `put_chargegroups_into_box`,
//!         md++/src/algorithm/integration/lattice_shift.cc

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::math::Periodicity;
use newen_core::topology::Topology;

use crate::constraints::box_periodicity;

/// Translates each charge group so that its reference point lies inside the box.
///
/// gromosXX uses the **centre of geometry** for solute charge groups and the **first atom** for
/// solvent ones; both are reproduced here.
#[derive(Debug, Clone, Default)]
pub struct LatticeShift;

impl LatticeShift {
    pub fn new() -> Self {
        Self
    }
}

impl Algorithm for LatticeShift {
    /// gromosXX applies the same shift before the first force evaluation (the tracker is part of
    /// the sequence from step 0), so the initial configuration is wrapped too.
    fn init(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        self.apply(topo, conf, sim)
    }

    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        let pbc = box_periodicity(conf);
        if matches!(pbc, Periodicity::Vacuum(_)) {
            return Ok(());
        }
        let n_solute_cg = topo.num_solute_chargegroups;
        for (k, cg) in topo.chargegroups.iter().enumerate() {
            if cg.atoms.is_empty() {
                continue;
            }
            // solute: centre of geometry; solvent: the first atom (gromosXX's two loops)
            let reference = if k < n_solute_cg {
                cg.center_of_geometry(&conf.current().pos)
            } else {
                conf.current().pos[cg.atoms[0]]
            };
            let shift = pbc.put_into_box(reference) - reference;
            if shift == newen_core::math::Vec3::ZERO {
                continue;
            }
            for &a in &cg.atoms {
                conf.current_mut().pos[a] += shift;
            }
        }
        Ok(())
    }

    fn name(&self) -> &str {
        "Lattice_Shift_Tracker"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::configuration::{Box as SimBox, BoxType};
    use newen_core::math::{Mat3, Vec3};
    use newen_core::topology::ChargeGroup;

    fn cubic(l: f64) -> SimBox {
        let vectors = Mat3::from_diagonal(Vec3::new(l, l, l));
        SimBox {
            box_type: BoxType::Rectangular,
            vectors,
            inv_vectors: vectors.inverse(),
        }
    }

    /// A two-atom solute charge group whose centre has drifted past the box edge moves back as a
    /// whole: the atoms keep their relative geometry, and the shift is exactly one box vector.
    #[test]
    fn a_charge_group_is_translated_as_a_whole() {
        let mut topo = Topology::new();
        topo.chargegroups = vec![ChargeGroup { atoms: vec![0, 1] }];
        topo.num_solute_chargegroups = 1;
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().box_config = cubic(3.0);
        conf.current_mut().pos[0] = Vec3::new(3.4, 0.5, 0.5);
        conf.current_mut().pos[1] = Vec3::new(3.6, 0.5, 0.5);
        let sim = SimulationState::new(0.002, 1);
        LatticeShift::new().apply(&topo, &mut conf, &sim).unwrap();
        let (a, b) = (conf.current().pos[0], conf.current().pos[1]);
        assert!(
            (b - a - Vec3::new(0.2, 0.0, 0.0)).length() < 1e-12,
            "geometry kept"
        );
        assert!(a.x < 3.0 && a.x >= 0.0, "moved into the box: {a:?}");
    }

    /// A solvent charge group is placed by its *first atom*, not its centre — gromosXX's rule, and
    /// it matters for a water straddling the boundary.
    #[test]
    fn solvent_is_placed_by_its_first_atom() {
        let mut topo = Topology::new();
        topo.chargegroups = vec![ChargeGroup {
            atoms: vec![0, 1, 2],
        }];
        topo.num_solute_chargegroups = 0; // this one is solvent
        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().box_config = cubic(3.0);
        conf.current_mut().pos[0] = Vec3::new(2.95, 0.5, 0.5); // inside
        conf.current_mut().pos[1] = Vec3::new(3.05, 0.5, 0.5); // outside
        conf.current_mut().pos[2] = Vec3::new(3.02, 0.6, 0.5); // outside
        let before = conf.current().pos.clone();
        let sim = SimulationState::new(0.002, 1);
        LatticeShift::new().apply(&topo, &mut conf, &sim).unwrap();
        // the first atom is already inside, so nothing moves
        for (a, b) in conf.current().pos.iter().zip(&before) {
            assert!((*a - *b).length() < 1e-12);
        }
    }

    /// Vacuum has no box: the algorithm is a no-op.
    #[test]
    fn vacuum_is_untouched() {
        let mut topo = Topology::new();
        topo.chargegroups = vec![ChargeGroup { atoms: vec![0] }];
        topo.num_solute_chargegroups = 1;
        let mut conf = Configuration::new(1, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(99.0, 0.0, 0.0);
        let sim = SimulationState::new(0.002, 1);
        LatticeShift::new().apply(&topo, &mut conf, &sim).unwrap();
        assert_eq!(conf.current().pos[0], Vec3::new(99.0, 0.0, 0.0));
    }
}
