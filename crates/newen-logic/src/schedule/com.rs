//! Centre-of-mass motion removal: on which steps the translation and the rotation of the whole
//! system are removed from the velocities.
//!
//! Step 0 follows INITIALISE's NTICOM (≥ 1 translation, ≥ 2 also rotation); every later step
//! follows COMTRANSROT's NSCM: every `|NSCM|` steps, translation when `NSCM ≠ 0`, rotation too
//! when `NSCM < 0`, but never rotation in a periodic box (gromosXX `remove_com_motion.cc:240-262`,
//! `in_parameter.cc:1306-1317`, `configuration.cc:555-560`).
//!
//! Theory: `spec/NewenTheory/Dynamics/Schedule.lean` (section "Centre-of-mass motion removal");
//! model and bridge: `spec/NewenSpec/Schedule.lean`.

/// `(translation, rotation)` at `step` for NTICOM `nticom`, NSCM `nscm` and `periodic`, exactly
/// the rule of `RemoveCOMMotion`.
///
/// The Lean refinement assumes `nscm > i64::MIN` (callers pass `i32` values).
///
/// spec: schedule.com.com_removal
#[allow(clippy::manual_is_multiple_of)]
pub fn com_removal(step: usize, nticom: i64, nscm: i64, periodic: bool) -> (bool, bool) {
    if step == 0 {
        return (nticom >= 1, nticom >= 2);
    }
    let mut every = nscm;
    if nscm < 0 {
        every = -nscm;
    }
    if every == 0 {
        return (false, false);
    }
    if step % (every as usize) != 0 {
        return (false, false);
    }
    (true, nscm < 0 && !periodic)
}

#[cfg(test)]
mod tests {
    use super::com_removal;

    #[test]
    fn step_zero_follows_nticom() {
        assert_eq!(com_removal(0, 0, -5, false), (false, false));
        assert_eq!(com_removal(0, 1, 0, false), (true, false));
        assert_eq!(com_removal(0, 2, 0, true), (true, true));
    }

    #[test]
    fn later_steps_follow_nscm() {
        assert_eq!(com_removal(10, 0, 10, false), (true, false));
        assert_eq!(com_removal(10, 0, -10, false), (true, true));
        assert_eq!(com_removal(10, 0, -10, true), (true, false));
        assert_eq!(com_removal(5, 0, -10, false), (false, false));
        assert_eq!(com_removal(5, 2, 0, false), (false, false));
    }
}
