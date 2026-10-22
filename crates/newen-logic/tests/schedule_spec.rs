//! The step schedules against their meaning and their Lean models.
//!
//! `output_misses_longrange`, on every `1 <= nsnb, ntwg <= GRID`:
//! - `misses_by_simulation` below walks the MD steps one by one, the way a run does, without
//!   using the closed form;
//! - `spec/vectors/schedule.json` holds the Lean model `missesModel` on the same grid.
//!
//! `due`, `frames_in`, `minimisation_stop`, `final_step` and `com_removal`:
//! - brute force: multiples walked by addition, frames counted step by step, gromosXX's
//!   top-of-loop writer simulated against an end-of-step writer that adds the stop frames, and
//!   `RemoveCOMMotion`'s flag-based rule as it was written before `com_removal`;
//! - `spec/vectors/schedule2.json` holds the Lean models on the grids below.
//!
//! The vector files are written by `lake exe spec-vectors` and checked in, so these tests need
//! no Lean.

use newen_logic::schedule::com::com_removal;
use newen_logic::schedule::output::{due, final_step, frames_in, minimisation_stop, StopFrames};
use newen_logic::schedule::output_misses_longrange;

/// Grid bound, equal to `grid` in `spec/SpecVectors.lean`.
const GRID: usize = 200;

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Walks steps `0..lcm(nsnb, ntwg)` (the schedule repeats after that): remember the last pairlist
/// step and whether it was an output step; report an output step that reuses one that was not.
fn misses_by_simulation(nsnb: usize, ntwg: usize) -> bool {
    let period = nsnb / gcd(nsnb, ntwg) * ntwg;
    let mut pairlist_was_output = true;
    for step in 0..period {
        let is_output = step % ntwg == 0;
        if step % nsnb == 0 {
            pairlist_was_output = is_output;
        }
        if is_output && !pairlist_was_output {
            return true;
        }
    }
    false
}

#[test]
fn rust_equals_simulation_on_grid() {
    let mut mismatches = Vec::new();
    for nsnb in 1..=GRID {
        for ntwg in 1..=GRID {
            if output_misses_longrange(nsnb, ntwg) != misses_by_simulation(nsnb, ntwg) {
                mismatches.push((nsnb, ntwg));
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "(nsnb, ntwg) mismatches: {mismatches:?}"
    );
}

/// The `0`/`1` rows of the vector file, `rows[nsnb - 1][ntwg - 1]`.
fn lean_rows() -> Vec<Vec<bool>> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/schedule.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    text.lines()
        .map(|l| l.trim().trim_end_matches(',').trim_matches('"'))
        .filter(|l| !l.is_empty() && l.bytes().all(|b| b == b'0' || b == b'1'))
        .map(|l| l.bytes().map(|b| b == b'1').collect())
        .collect()
}

#[test]
fn rust_equals_lean_model_on_grid() {
    let rows = lean_rows();
    assert_eq!(rows.len(), GRID, "vector file has {} rows", rows.len());
    let mut mismatches = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(
            row.len(),
            GRID,
            "row nsnb = {} has {} entries",
            i + 1,
            row.len()
        );
        for (j, &lean) in row.iter().enumerate() {
            if output_misses_longrange(i + 1, j + 1) != lean {
                mismatches.push((i + 1, j + 1));
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "(nsnb, ntwg) mismatches: {mismatches:?}"
    );
}

// Grid bounds, equal to `outputEvery`, `outputSteps`, `stopSteps`, `stopNmin`, `stopEvery`,
// `finalSteps`, `comSteps`, `comNscm` in `spec/SpecVectors.lean`.
const OUTPUT_EVERY: usize = 64;
const OUTPUT_STEPS: usize = 256;
const STOP_STEPS: usize = 16;
const STOP_NMIN: usize = 4;
const STOP_EVERY: usize = 4;
const FINAL_STEPS: usize = 20;
const COM_STEPS: usize = 60;
const COM_NSCM: i64 = 6;

/// `step` is `k · every` for some `k`, found by adding `every` up (never for `every == 0`).
fn due_by_walking(step: usize, every: usize) -> bool {
    if every == 0 {
        return false;
    }
    let mut multiple = 0;
    while multiple < step {
        multiple += every;
    }
    multiple == step
}

#[test]
fn due_and_frames_in_equal_walking_the_steps() {
    for every in 0..=OUTPUT_EVERY {
        let mut frames = 0;
        for step in 0..=OUTPUT_STEPS {
            assert_eq!(frames_in(step, every), frames, "n {step} every {every}");
            let d = due(step, every);
            assert_eq!(d, due_by_walking(step, every), "step {step} every {every}");
            if d {
                frames += 1;
            }
        }
    }
}

/// The frames of each output of a minimisation stopping at `k`: positions, forces, energies,
/// free energy.
type Frames = [Vec<usize>; 4];

/// gromosXX: at the top of each step `s = 0..=k` the positions of `s` and, for `s >= 1`, the
/// forces, energies and free energy of `s − 1`; then the final write of `k − 1`'s.
fn top_of_loop(k: usize, every: [usize; 4]) -> Frames {
    let mut f: Frames = Default::default();
    for s in 0..=k {
        if due_by_walking(s, every[0]) {
            f[0].push(s);
        }
        for o in 1..4 {
            if s >= 1 && due_by_walking(s - 1, every[o]) {
                f[o].push(s - 1);
            }
        }
    }
    for o in 1..4 {
        if due_by_walking(k - 1, every[o]) {
            f[o].push(k - 1);
        }
    }
    f
}

/// newen's loop: at the end of each completed step `s < k` all four outputs of `s`, then the
/// frames `minimisation_stop` names.
fn end_of_step(k: usize, every: [usize; 4], stop: StopFrames) -> Frames {
    let mut f: Frames = Default::default();
    for s in 0..k {
        for o in 0..4 {
            if due_by_walking(s, every[o]) {
                f[o].push(s);
            }
        }
    }
    if stop.positions {
        f[0].push(k);
    }
    for (o, on) in [(1, stop.forces), (2, stop.energies), (3, stop.free_energy)] {
        if on {
            f[o].push(k - 1);
        }
    }
    f
}

/// The four write frequencies of grid column `code`: NTWX, NTWF, NTWE, NTWG.
fn stop_every(code: usize) -> [usize; 4] {
    [3, 2, 1, 0].map(|k| code / STOP_EVERY.pow(k) % STOP_EVERY)
}

#[test]
fn the_minimisation_stop_writes_gromosxx_frames() {
    let mut stops = 0;
    for step in 0..STOP_STEPS {
        for nmin in 0..STOP_NMIN {
            for converged in [false, true] {
                for code in 0..STOP_EVERY.pow(4) {
                    let [x, f, e, g] = stop_every(code);
                    let got = minimisation_stop(step, nmin, converged, x, f, e, g);
                    assert_eq!(got.is_some(), converged && step > nmin);
                    if let Some(stop) = got {
                        assert_eq!(
                            end_of_step(step, [x, f, e, g], stop),
                            top_of_loop(step, [x, f, e, g]),
                            "step {step} nmin {nmin} every {:?}",
                            [x, f, e, g]
                        );
                        stops += 1;
                    }
                }
            }
        }
    }
    assert!(stops > 1000);
}

#[test]
fn final_step_is_the_stop_or_the_run_length() {
    for n in 0..=FINAL_STEPS {
        assert_eq!(final_step(n, None), n);
        for k in 0..=FINAL_STEPS {
            assert_eq!(final_step(n, Some(k)), k);
        }
    }
}

/// gromosXX's centre-of-mass removal, case by case, as data: `(step, NTICOM, NSCM, periodic,
/// (translation, rotation))`. Step 0 reads INITIALISE NTICOM only (`in_parameter.cc:1243-1258`:
/// 0 nothing, 1 translation, 2 both; `remove_com_motion.cc:240-243`); later steps remove when
/// `step % |NSCM| == 0` (`remove_com_motion.cc:245-248`), translation for any NSCM ≠ 0 and
/// rotation for NSCM < 0 (`in_parameter.cc:1308-1319`), never rotation under PBC
/// (`configuration.cc:555-560`, the periodic flag only, `remove_com_motion.cc:258-261`).
/// NTICOM outside 0..=2 is refused when the run is built.
const GROMOSXX_COM_CASES: &[(usize, i64, i64, bool, (bool, bool))] = &[
    // step 0: NTICOM, whatever NSCM and the boundary
    (0, 0, 0, false, (false, false)),
    (0, 0, -4, false, (false, false)),
    (0, 0, 4, true, (false, false)),
    (0, 1, 0, false, (true, false)),
    (0, 1, -4, true, (true, false)),
    (0, 2, 0, false, (true, true)),
    (0, 2, 4, false, (true, true)),
    (0, 2, -4, true, (true, true)),
    // NSCM = 0: never after step 0
    (1, 2, 0, false, (false, false)),
    (12, 1, 0, false, (false, false)),
    // NSCM > 0: translation every NSCM steps
    (4, 0, 4, false, (true, false)),
    (8, 2, 4, true, (true, false)),
    (6, 0, 4, false, (false, false)),
    (1, 0, 1, false, (true, false)),
    (7, 0, 1, true, (true, false)),
    // NSCM < 0: translation and rotation every |NSCM| steps, rotation not under PBC
    (3, 0, -3, false, (true, true)),
    (9, 1, -3, false, (true, true)),
    (3, 0, -3, true, (true, false)),
    (4, 0, -3, false, (false, false)),
    (5, 2, -3, true, (false, false)),
    (1, 0, -1, false, (true, true)),
    (60, 0, -6, true, (true, false)),
    (59, 0, -6, false, (false, false)),
];

#[test]
fn com_removal_equals_the_gromosxx_cases() {
    for &(step, nticom, nscm, periodic, expected) in GROMOSXX_COM_CASES {
        assert_eq!(
            com_removal(step, nticom, nscm, periodic),
            expected,
            "step {step} nticom {nticom} nscm {nscm} periodic {periodic}"
        );
    }
    assert_eq!(com_removal(7, 0, i64::MAX, false), (false, false));
    assert_eq!(com_removal(7, 0, -7, false), (true, true));
}

/// The rows of section `name` in `schedule2.json`, one per line: quoted strings or `[..]` lists.
fn schedule2_section(name: &str) -> Vec<String> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/schedule2.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    let start = format!("\"{name}\": [");
    text.lines()
        .skip_while(|l| l.trim() != start)
        .skip(1)
        .take_while(|l| l.trim() != "]" && l.trim() != "],")
        .map(|l| l.trim().trim_end_matches(',').trim_matches('"').to_string())
        .collect()
}

fn ints(row: &str) -> Vec<usize> {
    row.trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .map(|v| v.trim().parse().expect("integer in schedule2.json"))
        .collect()
}

#[test]
fn rust_equals_lean_models_on_schedule2() {
    let due_rows = schedule2_section("due");
    assert_eq!(due_rows.len(), OUTPUT_EVERY + 1);
    for (every, row) in due_rows.iter().enumerate() {
        assert_eq!(row.len(), OUTPUT_STEPS + 1);
        for (step, c) in row.bytes().enumerate() {
            assert_eq!(due(step, every), c == b'1', "due({step}, {every})");
        }
    }
    let frames_rows = schedule2_section("frames_in");
    assert_eq!(frames_rows.len(), OUTPUT_EVERY + 1);
    for (every, row) in frames_rows.iter().enumerate() {
        let lean = ints(row);
        assert_eq!(lean.len(), OUTPUT_STEPS + 1);
        for (n, &f) in lean.iter().enumerate() {
            assert_eq!(frames_in(n, every), f, "frames_in({n}, {every})");
        }
    }
    let stop_rows = schedule2_section("minimisation_stop");
    assert_eq!(stop_rows.len(), STOP_STEPS * STOP_NMIN * 2);
    for (i, row) in stop_rows.iter().enumerate() {
        let (step, nmin, converged) = (i / (2 * STOP_NMIN), i / 2 % STOP_NMIN, i % 2 == 1);
        assert_eq!(row.len(), STOP_EVERY.pow(4));
        for (code, c) in row.chars().enumerate() {
            let [x, f, e, g] = stop_every(code);
            let got = minimisation_stop(step, nmin, converged, x, f, e, g);
            let want = c.to_digit(16).map(|d| StopFrames {
                positions: d & 1 != 0,
                forces: d & 2 != 0,
                energies: d & 4 != 0,
                free_energy: d & 8 != 0,
            });
            assert_eq!(got, want, "row {i} column {code}");
        }
    }
    let final_rows = schedule2_section("final_step");
    assert_eq!(final_rows.len(), FINAL_STEPS + 1);
    for (n, row) in final_rows.iter().enumerate() {
        let lean = ints(row);
        assert_eq!(final_step(n, None), lean[0]);
        for k in 0..=FINAL_STEPS {
            assert_eq!(final_step(n, Some(k)), lean[k + 1]);
        }
    }
    let com_rows = schedule2_section("com_removal");
    let n_nscm = (2 * COM_NSCM + 1) as usize;
    assert_eq!(com_rows.len(), 5 * n_nscm * 2);
    for (i, row) in com_rows.iter().enumerate() {
        let nticom = (i / (2 * n_nscm)) as i64 - 1;
        let nscm = (i / 2 % n_nscm) as i64 - COM_NSCM;
        let periodic = i % 2 == 1;
        assert_eq!(row.len(), COM_STEPS + 1);
        for (step, c) in row.chars().enumerate() {
            let d = c.to_digit(16).expect("digit");
            assert_eq!(
                com_removal(step, nticom, nscm, periodic),
                (d & 1 != 0, d & 2 != 0),
                "row {i} step {step}"
            );
        }
    }
}
