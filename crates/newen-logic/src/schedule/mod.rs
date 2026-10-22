//! Step schedules: which MD steps do what.
//!
//! Theory: `spec/NewenTheory/Dynamics/Schedule.lean`; models and bridges:
//! `spec/NewenSpec/Schedule.lean`.

pub mod com;
pub mod output;

/// Whether some free-energy output step (every `ntwg` steps) follows a pairlist step (every
/// `nsnb` steps) that was not itself an output step.
///
/// With a twin-range pairlist the long-range interactions are evaluated only on pairlist steps
/// and reused until the next one. Under PRECALCLAM gromosXX adds the long-range λ_k energies to
/// the written rows only when that pairlist step was also an output step, so such an output step
/// writes λ_k rows without their long-range part. That happens on some output step exactly when
/// `nsnb > 1`, `ntwg > 0`, and neither of `nsnb` and `ntwg` divides the other. `ntwg == 0` means
/// no output; `nsnb <= 1` means every step is a pairlist step, so nothing is reused.
///
/// spec: schedule.output_misses_longrange
// `%` rather than `is_multiple_of`: plain arithmetic keeps the Lean extraction self-explanatory.
#[allow(clippy::manual_is_multiple_of)]
pub fn output_misses_longrange(nsnb: usize, ntwg: usize) -> bool {
    if nsnb <= 1 || ntwg == 0 {
        return false;
    }
    ntwg % nsnb != 0 && nsnb % ntwg != 0
}

#[cfg(test)]
mod tests {
    use super::output_misses_longrange;

    #[test]
    fn neither_divides_the_other() {
        assert!(output_misses_longrange(5, 2));
        assert!(output_misses_longrange(4, 6));
    }

    #[test]
    fn one_divides_the_other() {
        assert!(!output_misses_longrange(5, 10));
        assert!(!output_misses_longrange(5, 1));
        assert!(!output_misses_longrange(5, 5));
    }

    #[test]
    fn every_step_a_pairlist_step_or_no_output() {
        assert!(!output_misses_longrange(1, 2));
        assert!(!output_misses_longrange(0, 2));
        assert!(!output_misses_longrange(5, 0));
    }
}
