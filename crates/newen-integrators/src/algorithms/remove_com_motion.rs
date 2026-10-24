//! Centre of mass motion removal algorithm.
//!
//! Equivalent to GROMOS `algorithm::Remove_COM_Motion`
//! (md++/src/algorithm/constraints/remove_com_motion.cc).
//! Removes translational and/or rotational centre-of-mass motion from velocities.
//!
//! At step 0 (initial): controlled by NTICOM from the INITIALISE block
//! (`start.remove_com_translation`/`remove_com_rotation`, `in_parameter.cc:1242-1258`).
//! At step > 0 (periodic): controlled by the *sign* of NSCM from the COMTRANSROT
//! block (`in_parameter.cc:1305-1319`): positive removes translation only every
//! `NSCM` steps, negative removes translation *and* rotation every `|NSCM|` steps,
//! zero turns periodic removal off.
//!
//! **PBC position convention.** GROMOS calls `gather_chargegroups` during
//! configuration init (`periodicity.cc:175`), which uses `put_into_box` to wrap
//! each chargegroup's reference atom to the `[-L/2, L/2]` range (nearest image
//! relative to the origin). This leaves `conf.current().pos` in the `[-L/2, L/2]`
//! convention rather than the `[0, L]` range written by the coordinate file. The
//! angular momentum and inertia tensor must use the same convention; otherwise the
//! computed COM and angular momentum diverge from GROMOS. We replicate this by
//! applying minimum-image wrapping (GROMOS `put_into_box`, i.e. nearest image
//! relative to origin) to each position *locally* inside `remove_com_rotation`.
//! The stored `conf.current().pos` is NOT modified.
//!
//! **PBC periodic-rotation suppression.** GROMOS also disables *periodic* COM
//! rotation removal when the boundary condition is not vacuum
//! (`configuration.cc:555-560`, `param.centreofmass.remove_rot = false`). The
//! *initial* removal at step 0 (controlled by NTICOM) is unaffected. We mirror
//! this: when the box is periodic, no rotation is removed after step 0, whatever
//! the sign of NSCM. Which step removes what is `newen_logic::schedule::com::com_removal`.

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::{BoxType, Configuration};
use newen_core::math::{Mat3, Vec3};
use newen_core::topology::Topology;
use newen_logic::schedule::com::com_removal;

/// GROMOS `math::epsilon` (gmath.h:689) — used for the inertia-tensor
/// determinant singularity check (signed comparison, not `abs()`).
const EPSILON: f64 = 0.000000000001;

/// Removes centre-of-mass translational and/or rotational motion from velocities.
///
/// GROMOS convention: only modifies `conf.current().vel`.
#[derive(Debug, Clone)]
pub struct RemoveCOMMotion {
    /// NTICOM from INITIALISE block: 0=off, 1=remove translation, 2=remove translation+rotation
    pub nticom: i32,
    /// NSCM from the COMTRANSROT block, with its sign (`in_parameter.cc:1305-1319`)
    nscm: i32,
}

impl RemoveCOMMotion {
    /// `nscm` carries its sign from the COMTRANSROT block (`in_parameter.cc:1305-1319`):
    /// `> 0` removes translation every `nscm` steps, `< 0` removes translation and
    /// rotation every `|nscm|` steps, `0` disables periodic removal.
    pub fn new(nticom: i32, nscm: i32) -> Self {
        Self { nticom, nscm }
    }
}

/// Calculate (and optionally remove) translational COM motion.
///
/// Mirrors GROMOS `Remove_COM_Motion::remove_com_translation`
/// (`remove_com_motion.cc:90-119`): `v_com = Σ(m_i·v_i) / Σ(m_i)`,
/// `E_kin_trans = 0.5 * M * |v_com|²`, then `v_i -= v_com` if `remove`.
fn remove_com_translation(topo: &Topology, conf: &mut Configuration, remove: bool) -> f64 {
    let n_atoms = topo.num_atoms();
    let mut com_v = Vec3::ZERO;
    let mut com_mass = 0.0;

    for i in 0..n_atoms {
        com_mass += topo.mass[i];
        com_v += topo.mass[i] * conf.current().vel[i];
    }
    com_v /= com_mass;

    let ekin_trans = 0.5 * com_mass * com_v.length_squared();

    if remove {
        for i in 0..n_atoms {
            conf.current_mut().vel[i] -= com_v;
        }
    }

    ekin_trans
}

/// Apply minimum-image wrapping to a position, equivalent to GROMOS
/// `Periodicity<rectangular>::put_into_box` (`periodicity.cc:38-42`):
/// `nearest_image(v, origin=(0,0,0), v)` which maps each component to
/// `(-L/2, L/2]` via `v - round(v/L)*L`.
#[inline]
fn put_into_box_rect(pos: Vec3, box_dims: Vec3) -> Vec3 {
    Vec3::new(
        pos.x - (pos.x / box_dims.x).round() * box_dims.x,
        pos.y - (pos.y / box_dims.y).round() * box_dims.y,
        pos.z - (pos.z / box_dims.z).round() * box_dims.z,
    )
}

/// Calculate (and optionally remove) rotational COM motion.
///
/// Mirrors GROMOS `Remove_COM_Motion::remove_com_rotation`
/// (`remove_com_motion.cc:121-222`): builds the COM position/velocity at the
/// velocity time point (`r - 0.5*dt*v`), the angular momentum `L = Σ m (r×v)`
/// and inertia tensor `I`, then inverts `I` via the explicit cofactor-expansion
/// formula (note: GROMOS checks `denom < math::epsilon`, a *signed* comparison,
/// not `denom.abs() < epsilon` — this is intentional and must be preserved
/// bit-for-bit), and computes `ω = I⁻¹L / det(I)`.
///
/// **Position convention.** GROMOS wraps positions to `[-L/2, L/2]` via
/// `gather_chargegroups` during init (equivalent to `put_into_box` for each atom).
/// We replicate this locally: each `pos[i]` is wrapped before use so that the
/// COM and angular momentum match GROMOS bit-for-bit.
fn remove_com_rotation(topo: &Topology, conf: &mut Configuration, dt: f64, remove: bool) -> f64 {
    let n_atoms = topo.num_atoms();

    // Determine box dimensions for wrapping (vacuum = no wrapping).
    let box_dims = match conf.current().box_config.box_type {
        BoxType::Rectangular => Some(conf.current().box_config.dimensions()),
        _ => None,
    };

    // Helper: wrap a position to [-L/2, L/2] if periodic (GROMOS put_into_box).
    let wrap = |pos: Vec3| -> Vec3 {
        if let Some(l) = box_dims {
            put_into_box_rect(pos, l)
        } else {
            pos
        }
    };

    let mut com_v = Vec3::ZERO;
    let mut com_r = Vec3::ZERO;
    let mut com_mass = 0.0;

    for i in 0..n_atoms {
        let m = topo.mass[i];
        let r = wrap(conf.current().pos[i]);
        com_mass += m;
        com_v += m * conf.current().vel[i];
        com_r += m * r - 0.5 * m * conf.current().vel[i] * dt;
    }
    com_v /= com_mass;
    com_r /= com_mass;

    let mut com_l = Vec3::ZERO;
    let mut com_i = Mat3::ZERO;

    for i in 0..n_atoms {
        let m = topo.mass[i];
        let r = wrap(conf.current().pos[i]) - 0.5 * dt * conf.current().vel[i] - com_r;

        com_l += m * r.cross(conf.current().vel[i] - com_v);

        // Inertia tensor (remove_com_motion.cc:185-193); `Mat3` is column-major
        // (glam `DMat3`), so `x_axis.x` = row 0 col 0, `x_axis.y` = row 1 col 0, etc.
        com_i.x_axis.x += m * (r.y * r.y + r.z * r.z); // I(0,0)
        com_i.y_axis.y += m * (r.x * r.x + r.z * r.z); // I(1,1)
        com_i.z_axis.z += m * (r.x * r.x + r.y * r.y); // I(2,2)
        com_i.x_axis.y -= m * r.x * r.y; // I(1,0) = I(0,1)
        com_i.y_axis.x -= m * r.x * r.y; // I(0,1)
        com_i.x_axis.z -= m * r.x * r.z; // I(2,0) = I(0,2)
        com_i.z_axis.x -= m * r.x * r.z; // I(0,2)
        com_i.y_axis.z -= m * r.y * r.z; // I(2,1) = I(1,2)
        com_i.z_axis.y -= m * r.y * r.z; // I(1,2)
    }

    // Explicit cofactor-expansion inversion (remove_com_motion.cc:197-216).
    // Variable names: i00 = I(0,0), i01 = I(0,1) = I(1,0), etc.
    let i00 = com_i.x_axis.x; // I(0,0)
    let i01 = com_i.x_axis.y; // I(1,0) = I(0,1) (symmetric)
    let i02 = com_i.x_axis.z; // I(2,0) = I(0,2) (symmetric)
    let i11 = com_i.y_axis.y; // I(1,1)
    let i12 = com_i.y_axis.z; // I(2,1) = I(1,2) (symmetric)
    let i20 = com_i.z_axis.x; // I(0,2) == i02 (symmetric matrix)
    let i22 = com_i.z_axis.z; // I(2,2)

    let denom = -i20 * i20 * i11 + 2.0 * i01 * i02 * i12 - i00 * i12 * i12 - i01 * i01 * i22
        + i00 * i11 * i22;

    let ii00 = -i12 * i12 + i11 * i22;
    let ii01 = i02 * i12 - i01 * i22;
    let ii02 = -i02 * i11 + i01 * i12;
    let ii11 = -i02 * i02 + i00 * i22;
    let ii12 = i01 * i02 - i00 * i12;
    let ii22 = -i01 * i01 + i00 * i11;

    // GROMOS compares the signed determinant against `math::epsilon`
    // (remove_com_motion.cc:212), not its absolute value.
    let com_o = if denom < EPSILON {
        Vec3::ZERO
    } else {
        Vec3::new(
            ii00 * com_l.x + ii01 * com_l.y + ii02 * com_l.z,
            ii01 * com_l.x + ii11 * com_l.y + ii12 * com_l.z,
            ii02 * com_l.x + ii12 * com_l.y + ii22 * com_l.z,
        ) / denom
    };

    let ekin_rot = 0.5 * com_o.dot(com_l);
    log::debug!(
        "    com_r={:.12e} {:.12e} {:.12e}  com_l={:.12e} {:.12e} {:.12e}  denom={:.6e}  com_o={:.12e} {:.12e} {:.12e}",
        com_r.x, com_r.y, com_r.z,
        com_l.x, com_l.y, com_l.z,
        denom,
        com_o.x, com_o.y, com_o.z,
    );

    if remove {
        for i in 0..n_atoms {
            let r = wrap(conf.current().pos[i]) - 0.5 * dt * conf.current().vel[i] - com_r;
            conf.current_mut().vel[i] -= com_o.cross(r);
        }

        // sanity check: residual angular momentum about COM should be ~zero
        if log::log_enabled!(log::Level::Debug) {
            let mut residual_l = Vec3::ZERO;
            let mut residual_v = Vec3::ZERO;
            for i in 0..n_atoms {
                let m = topo.mass[i];
                let r = wrap(conf.current().pos[i]) - 0.5 * dt * conf.current().vel[i] - com_r;
                residual_l += m * r.cross(conf.current().vel[i]);
                residual_v += m * conf.current().vel[i];
            }
            log::debug!(
                "    residual_L={:?} residual_Mv={:?}",
                residual_l,
                residual_v
            );
        }
    }

    ekin_rot
}

impl Algorithm for RemoveCOMMotion {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        sim: &SimulationState,
    ) -> Result<(), String> {
        // `Remove_COM_Motion::apply` (remove_com_motion.cc:240-264): NTICOM at step 0, NSCM
        // after it, no periodic rotation removal in PBC (configuration.cc:555-560)
        let periodic = !matches!(
            conf.current().box_config.box_type,
            newen_core::configuration::BoxType::Vacuum
        );
        let (remove_trans, remove_rot) = com_removal(
            sim.step,
            i64::from(self.nticom),
            i64::from(self.nscm),
            periodic,
        );

        if !remove_trans && !remove_rot {
            return Ok(());
        }

        if remove_trans {
            let ekin_trans = remove_com_translation(topo, conf, true);
            log::debug!("  COM translation removed: E_kin_trans={:.6e}", ekin_trans);
        }
        if remove_rot {
            let ekin_rot = remove_com_rotation(topo, conf, sim.dt, true);
            log::debug!("  COM rotation removed: E_kin_rot={:.6e}", ekin_rot);
        }

        Ok(())
    }

    fn name(&self) -> &str {
        "Remove_COM_Motion"
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

    fn two_atoms() -> (Topology, Configuration) {
        let topo = bare(&[1.0, 3.0]);
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().vel[0] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().vel[1] = Vec3::new(1.0, 0.0, 0.0);
        (topo, conf)
    }

    fn momentum(topo: &Topology, conf: &Configuration) -> Vec3 {
        topo.mass[0] * conf.current().vel[0] + topo.mass[1] * conf.current().vel[1]
    }

    #[test]
    fn nticom_1_removes_the_initial_translation() {
        let (topo, mut conf) = two_atoms();
        let sim = SimulationState::new(0.002, 10); // step 0
        RemoveCOMMotion::new(1, 0)
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        let p = momentum(&topo, &conf);
        assert!(p.length() < 1e-14, "COM momentum left: {p:?}");
    }

    #[test]
    fn nticom_0_leaves_the_initial_velocities() {
        let (topo, mut conf) = two_atoms();
        let sim = SimulationState::new(0.002, 10);
        RemoveCOMMotion::new(0, 0)
            .apply(&topo, &mut conf, &sim)
            .unwrap();
        assert_eq!(conf.current().vel[0], Vec3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn nscm_sets_the_period() {
        let (topo, mut conf) = two_atoms();
        let mut sim = SimulationState::new(0.002, 10);
        let mut alg = RemoveCOMMotion::new(0, 5);
        sim.step = 3;
        alg.apply(&topo, &mut conf, &sim).unwrap();
        assert_eq!(
            conf.current().vel[0],
            Vec3::new(1.0, 0.0, 0.0),
            "no removal off-period"
        );
        sim.step = 5;
        alg.apply(&topo, &mut conf, &sim).unwrap();
        assert!(
            momentum(&topo, &conf).length() < 1e-14,
            "removal expected at step 5"
        );
    }
}
