//! `output_misses_longrange` against its meaning, on every `1 <= nsnb, ntwg <= GRID`.
//!
//! Two independent references:
//! - `misses_by_simulation` below walks the MD steps one by one, the way a run does, without
//!   using the closed form;
//! - `spec/vectors/schedule.json` holds the Lean model `missesModel` on the same grid, written by
//!   `lake exe spec-vectors` (checked in, so this test needs no Lean).

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
