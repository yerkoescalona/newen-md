//! Roto-translational constraints (`ROTTRANS`).
//!
//! Port of gromosXX `algorithm::Rottrans_Constraints`
//! (`md++/src/algorithm/constraints/rottrans.cc`). Removes the six rigid-body degrees of
//! freedom of atoms `1..=RTCLAST` — three translations and three rotations — by *constraining*
//! them rather than by subtracting the centre-of-mass motion every so often
//! ([`RemoveCOMMotion`](super::remove_com_motion::RemoveCOMMotion), `COMTRANSROT`). The
//! constraint keeps the group's total momentum and its angular momentum about the reference
//! orientation at zero for every step, which is what a run with position restraints switched
//! off needs if the solute is not to drift and tumble.
//!
//! Runs immediately after SHAKE, as in `create_constraints.cc:294-305`.
//!
//! **Reference state.** `init` centres the system on the constrained group's centre of mass,
//! stores the resulting positions of `1..=RTCLAST` as the reference orientation, and inverts
//! the two `theta` matrices (mass and inertia). gromosXX can instead read that state from the
//! configuration's `ROTTRANSREFPOS` block (`NTIRTC = 0`); this port always resets it, which is
//! `NTIRTC = 1`.

use crate::constraints::box_periodicity;
use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::Configuration;
use newen_core::gather::put_chargegroups_into_box;
use newen_core::math::{Mat3, Vec3};
use newen_core::topology::Topology;

/// Constrains the rigid-body motion of the first `last` atoms.
#[derive(Debug, Clone)]
pub struct RottransConstraints {
    /// RTCLAST — the number of leading atoms that are constrained.
    last: usize,
    /// The reference orientation: the constrained atoms' positions after centring, at `init`.
    reference: Vec<Vec3>,
    /// `theta_inv_trans` — diagonal, `1 / total mass`.
    inv_mass_total: f64,
    /// `theta_inv_rot` — the inverse of the reference inertia tensor.
    inv_theta_rot: Mat3,
}

impl RottransConstraints {
    /// `last` is RTCLAST, the 1-based index of the last constrained atom.
    pub fn new(last: usize) -> Self {
        Self {
            last,
            reference: Vec::new(),
            inv_mass_total: 0.0,
            inv_theta_rot: Mat3::ZERO,
        }
    }

    pub fn last(&self) -> usize {
        self.last
    }
}

impl Algorithm for RottransConstraints {
    /// gromosXX `_init`: centre on the group's centre of mass, keep the centred positions as the
    /// reference orientation, and invert the mass and inertia matrices.
    fn init(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        let last = self.last.min(topo.mass.len());
        if last == 0 {
            return Err("ROTTRANS: RTCLAST is 0, so nothing is constrained".to_string());
        }

        let total_mass: f64 = topo.mass[..last].iter().sum();
        if total_mass <= 0.0 {
            return Err(format!(
                "ROTTRANS: atoms 1..{last} have no mass, so they cannot be constrained"
            ));
        }
        // gromosXX `_init` (`rottrans.cc:201-232`): the whole system, both states, is shifted by
        // the centre of mass of the chain-gathered group and put back into the box, so the
        // reference is centred (its `ROTTRANSREFPOS` and trajectory, `aladip_multibath_rottrans`).
        // Uncentred, the rotational correction λ_rot × r moves the group's centre of mass.
        let pbc = box_periodicity(conf);
        let (cgs, n_solute_cgs) = (&topo.chargegroups, topo.num_solute_chargegroups);
        put_chargegroups_into_box(&mut conf.current_mut().pos, cgs, n_solute_cgs, &pbc);
        let com = {
            let pos = &conf.current().pos;
            let mut prev = pos[0];
            let mut sum = Vec3::ZERO;
            for (i, &x) in pos[..last].iter().enumerate() {
                let p = pbc.nearest_image(x, prev);
                sum += (p + prev) * topo.mass[i];
                prev += p;
            }
            sum / total_mass
        };
        for state in [true, false] {
            let pos = if state {
                &mut conf.current_mut().pos
            } else {
                &mut conf.old_mut().pos
            };
            for x in pos.iter_mut() {
                *x -= com;
            }
            put_chargegroups_into_box(pos, cgs, n_solute_cgs, &pbc);
        }

        self.reference = conf.current().pos[..last].to_vec();
        self.inv_mass_total = 1.0 / total_mass;

        // theta[3..6] of `rottrans.cc`: the inertia tensor of the reference orientation.
        let mut theta = Mat3::ZERO;
        for (i, r) in self.reference.iter().enumerate() {
            let m = topo.mass[i];
            theta.x_axis.x += m * (r.z * r.z + r.y * r.y);
            theta.x_axis.y += m * (-r.y * r.x);
            theta.x_axis.z += m * (-r.x * r.z);
            theta.y_axis.x += m * (-r.x * r.y);
            theta.y_axis.y += m * (r.x * r.x + r.z * r.z);
            theta.y_axis.z += m * (-r.y * r.z);
            theta.z_axis.x += m * (-r.x * r.z);
            theta.z_axis.y += m * (-r.y * r.z);
            theta.z_axis.z += m * (r.x * r.x + r.y * r.y);
        }

        // gromosXX inverts it by hand from the three rows (`d` is the determinant); a linear
        // molecule or a single atom makes `d` zero and the rotational constraint meaningless.
        let (t0, t1, t2) = (theta.x_axis, theta.y_axis, theta.z_axis);
        let d = t0.cross(t1).dot(t2);
        if d.abs() < 1e-12 {
            return Err(format!(
                "ROTTRANS: atoms 1..{last} have a singular inertia tensor (collinear or a \
                 single atom); there is no rotation to constrain"
            ));
        }
        self.inv_theta_rot = Mat3::from_cols(t1.cross(t2) / d, t0.cross(t2) / -d, t0.cross(t1) / d);
        Ok(())
    }

    /// gromosXX `_apply`: the constraint multipliers from this step's displacement, then the
    /// position correction and a velocity recomputed from it for *every* atom.
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        let last = self.last.min(topo.mass.len());
        let mut c_trans = Vec3::ZERO;
        let mut c_rot = Vec3::ZERO;
        for i in 0..last {
            let diff = conf.current().pos[i] - conf.old().pos[i];
            let m = topo.mass[i];
            let r = self.reference[i];
            c_trans += diff * m;
            c_rot += r.cross(diff) * m;
        }

        let lambda_trans = -c_trans * self.inv_mass_total;
        let lambda_rot = -(self.inv_theta_rot * c_rot);

        for i in 0..last {
            let r = self.reference[i];
            conf.current_mut().pos[i] += lambda_trans + lambda_rot.cross(r);
        }

        let inv_dt = 1.0 / sim.dt;
        for i in 0..conf.current().pos.len() {
            conf.current_mut().vel[i] = (conf.current().pos[i] - conf.old().pos[i]) * inv_dt;
        }
        Ok(())
    }

    fn name(&self) -> &str {
        "Rottrans_Constraints"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::configuration::Box as SimBox;

    /// Four atoms in a square, given a rigid-body kick: the constraint must take all of it back.
    fn square() -> (Topology, Configuration) {
        let mut topo = Topology::new();
        for _ in 0..4 {
            topo.mass.push(12.0);
            topo.inverse_mass.push(1.0 / 12.0);
        }
        let mut conf = Configuration::new(4, 1, 1);
        let square = [
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(-1.0, 1.0, 0.0),
            Vec3::new(-1.0, -1.0, 0.2),
            Vec3::new(1.0, -1.0, -0.2),
        ];
        conf.current_mut().pos = square.to_vec();
        conf.current_mut().box_config = SimBox::rectangular(10.0, 10.0, 10.0);
        conf.copy_current_to_old();
        (topo, conf)
    }

    fn momenta(topo: &Topology, conf: &Configuration, reference: &[Vec3], dt: f64) -> (Vec3, Vec3) {
        let mut p = Vec3::ZERO;
        let mut l = Vec3::ZERO;
        for i in 0..reference.len() {
            let v = (conf.current().pos[i] - conf.old().pos[i]) / dt;
            p += v * topo.mass[i];
            l += reference[i].cross(v) * topo.mass[i];
        }
        (p, l)
    }

    #[test]
    fn a_rigid_body_kick_is_taken_back() {
        let (topo, mut conf) = square();
        let sim = SimulationState::new(0.002, 1);
        let mut rtc = RottransConstraints::new(4);
        rtc.init(&topo, &mut conf, &sim).unwrap();

        // translate and rotate the whole group between old and current
        let shift = Vec3::new(0.03, -0.02, 0.01);
        let omega = Vec3::new(0.0, 0.0, 0.05);
        for i in 0..4 {
            let r = conf.current().pos[i];
            conf.current_mut().pos[i] = r + shift + omega.cross(r);
        }
        let reference = rtc.reference.clone();
        let (p0, l0) = momenta(&topo, &conf, &reference, sim.dt);
        assert!(p0.length() > 1.0 && l0.length() > 1.0, "kick was too small");

        rtc.apply(&topo, &mut conf, &sim).unwrap();

        let (p, l) = momenta(&topo, &conf, &reference, sim.dt);
        assert!(p.length() < 1e-9, "linear momentum left: {p:?}");
        assert!(l.length() < 1e-9, "angular momentum left: {l:?}");
    }

    #[test]
    fn internal_motion_survives() {
        let (topo, mut conf) = square();
        let sim = SimulationState::new(0.002, 1);
        let mut rtc = RottransConstraints::new(4);
        rtc.init(&topo, &mut conf, &sim).unwrap();

        // a breathing mode: no net momentum, no net angular momentum
        for i in 0..4 {
            let r = conf.current().pos[i];
            conf.current_mut().pos[i] = r * 1.01;
        }
        let before = conf.current().pos.clone();
        rtc.apply(&topo, &mut conf, &sim).unwrap();
        for i in 0..4 {
            assert!(
                (conf.current().pos[i] - before[i]).length() < 1e-12,
                "atom {i} was moved: {:?} -> {:?}",
                before[i],
                conf.current().pos[i]
            );
        }
    }

    #[test]
    fn a_single_atom_has_no_rotation_to_constrain() {
        let (topo, mut conf) = square();
        let sim = SimulationState::new(0.002, 1);
        let err = RottransConstraints::new(1)
            .init(&topo, &mut conf, &sim)
            .unwrap_err();
        assert!(err.contains("singular inertia tensor"), "{err}");
    }
}
