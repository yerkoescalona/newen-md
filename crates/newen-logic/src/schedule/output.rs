//! Output steps: which steps write a trajectory frame, how many frames a run writes, and what a
//! converged minimisation writes at the step it stops on.
//!
//! Steps are 0-based: step `s` is the `s`-th step of the run, and the frame written for it
//! describes the configuration at its start (gromosXX `sim.steps()`). WRITETRAJ's `NTWX`, `NTWF`,
//! `NTWE`, `NTWG` are the `every` arguments; 0 turns that output off.
//!
//! Theory: `spec/NewenTheory/Dynamics/Schedule.lean` (section "Output steps"); model and bridge:
//! `spec/NewenSpec/Schedule.lean`.

/// A frame is due at `step` when `every > 0` and `step` is a multiple of `every` (0-based steps;
/// gromosXX `steps % n == 0`).
///
/// spec: schedule.output.due
// `%` rather than `is_multiple_of`: plain arithmetic keeps the Lean extraction self-explanatory.
#[allow(clippy::manual_is_multiple_of)]
pub fn due(step: usize, every: usize) -> bool {
    if every == 0 {
        return false;
    }
    step % every == 0
}

/// How many of the steps `0..n_steps` are due: 0 when `every == 0`, else ⌈n_steps / every⌉.
///
/// spec: schedule.output.frames_in
#[allow(clippy::manual_is_multiple_of)]
pub fn frames_in(n_steps: usize, every: usize) -> usize {
    if every == 0 {
        return 0;
    }
    let whole = n_steps / every;
    if n_steps % every == 0 {
        whole
    } else {
        whole + 1
    }
}

/// The frames a converged minimisation writes at the step it stops on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StopFrames {
    /// the positions of the stop step (NTWX)
    pub positions: bool,
    /// the forces of the step before (NTWF)
    pub forces: bool,
    /// the energies of the step before (NTWE)
    pub energies: bool,
    /// the free-energy frame of the step before (NTWG)
    pub free_energy: bool,
}

/// Whether a minimisation stops at `step` (`step > nmin && converged`, `converged` decided by the
/// caller from |ΔE| < DELE) and, if so, which frames it writes: the positions of `step` when due,
/// and the forces, energies and free-energy frame of `step − 1` when due (gromosXX writes the
/// top-of-loop positions of `step` and its final frames of `step − 1`).
///
/// There are no velocity frames: gromosXX refuses NTWV ≠ 0 in a minimisation, and so does
/// newen's plan validation.
///
/// spec: schedule.output.minimisation_stop
pub fn minimisation_stop(
    step: usize,
    nmin: usize,
    converged: bool,
    ntwx: usize,
    ntwf: usize,
    ntwe: usize,
    ntwg: usize,
) -> Option<StopFrames> {
    if !converged || step <= nmin {
        return None;
    }
    // step > nmin >= 0, so the step before exists
    let last = step - 1;
    Some(StopFrames {
        positions: due(step, ntwx),
        forces: due(last, ntwf),
        energies: due(last, ntwe),
        free_energy: due(last, ntwg),
    })
}

/// The step the final configuration carries: the stop step of a converged minimisation, else
/// `n_steps` (gromosXX `sim.steps()` after the loop).
///
/// spec: schedule.output.final_step
pub fn final_step(n_steps: usize, stopped_at: Option<usize>) -> usize {
    match stopped_at {
        Some(step) => step,
        None => n_steps,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_turns_output_off() {
        assert!(!due(0, 0));
        assert!(!due(5, 0));
        assert_eq!(frames_in(10, 0), 0);
    }

    #[test]
    fn frames_on_multiples() {
        assert!(due(0, 5) && due(10, 5) && !due(6, 5));
        assert_eq!(frames_in(10, 5), 2);
        assert_eq!(frames_in(11, 5), 3);
        assert_eq!(frames_in(0, 5), 0);
        assert_eq!(frames_in(usize::MAX, 2), usize::MAX / 2 + 1);
    }

    #[test]
    fn a_minimisation_stops_after_nmin() {
        assert_eq!(minimisation_stop(1, 1, true, 1, 1, 1, 1), None);
        assert_eq!(minimisation_stop(6, 1, false, 1, 1, 1, 1), None);
        // aladip_vacuum_em_ntw5: NTWX = NTWF = NTWE = 5, stop at step 6
        assert_eq!(
            minimisation_stop(6, 1, true, 5, 5, 5, 0),
            Some(StopFrames {
                positions: false,
                forces: true,
                energies: true,
                free_energy: false
            })
        );
    }

    #[test]
    fn the_final_step() {
        assert_eq!(final_step(500, Some(6)), 6);
        assert_eq!(final_step(500, None), 500);
    }
}
