//! `total_dof` against its meaning and its Lean model.
//!
//! Two independent references:
//! - `dof_by_molecules` below builds the system the way `newen-run` sees it (one solute, `k`
//!   copies of a solvent molecule) and counts it molecule by molecule, atom by atom and
//!   constraint by constraint, without the closed form;
//! - `spec/vectors/dof.json` holds the Lean model `totalDofModel` on a grid of its four inputs,
//!   written by `lake exe spec-vectors` (checked in, so this test needs no Lean).

use newen_logic::thermo::dof::total_dof;

/// Count bound of the vector grid, equal to `dofGrid` in `spec/SpecVectors.lean`.
const GRID: u64 = 10;
/// The `ndfmin` values of the vector grid, equal to `ndfmins` in `spec/SpecVectors.lean`.
const NDFMINS: std::ops::RangeInclusive<i64> = -3..=6;

/// A molecule as the theory sees it: atoms and the constraints imposed on them.
struct Molecule {
    atoms: u64,
    constraints: u64,
}

/// Every atom brings 3 coordinates, every constraint takes one away, and `ndfmin` is removed
/// once from the whole system.
fn dof_by_molecules(molecules: &[Molecule], ndfmin: i64) -> i64 {
    let mut dof: i64 = 0;
    for m in molecules {
        for _ in 0..m.atoms {
            dof += 3;
        }
        for _ in 0..m.constraints {
            dof -= 1;
        }
    }
    dof - ndfmin
}

/// Solute plus `k` solvent molecules, through the arguments `newen-run::total_dof` passes:
/// all atoms, `k × (constraints per solvent molecule)`, and the solute constraint pairs.
#[test]
fn rust_equals_molecule_count_on_grid() {
    let mut mismatches = Vec::new();
    for solute_atoms in 0..=8u64 {
        for solute_constraints in 0..=8u64 {
            for k in 0..=6u64 {
                for solvent_atoms in 1..=4u64 {
                    for solvent_constraints in 0..=4u64 {
                        for ndfmin in [0, 3, 6, 7] {
                            let mut system = vec![Molecule {
                                atoms: solute_atoms,
                                constraints: solute_constraints,
                            }];
                            for _ in 0..k {
                                system.push(Molecule {
                                    atoms: solvent_atoms,
                                    constraints: solvent_constraints,
                                });
                            }
                            let expected = dof_by_molecules(&system, ndfmin);
                            let got = total_dof(
                                solute_atoms + k * solvent_atoms,
                                k * solvent_constraints,
                                solute_constraints,
                                ndfmin,
                            );
                            if got != expected {
                                mismatches.push((
                                    solute_atoms,
                                    solute_constraints,
                                    k,
                                    solvent_atoms,
                                    solvent_constraints,
                                    ndfmin,
                                    got,
                                    expected,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(mismatches.is_empty(), "mismatches: {mismatches:?}");
}

/// Large counts, up to the bound of the Lean refinement, against 128-bit arithmetic.
#[test]
fn rust_equals_wide_arithmetic_on_large_counts() {
    let bound = i64::MAX as i128;
    for &n in &[0u64, 1, 1_000_000, 1 << 40, (i64::MAX / 3) as u64] {
        for &s in &[0u64, 7, 1 << 20, n] {
            for &c in &[0u64, 3, n / 2] {
                for &f in &[-6i64, 0, 3, 6] {
                    let wide = 3 * n as i128 + s as i128 + c as i128 + (f as i128).abs();
                    if wide > bound {
                        continue;
                    }
                    let expected = 3 * n as i128 - s as i128 - c as i128 - f as i128;
                    assert_eq!(total_dof(n, s, c, f) as i128, expected, "{n} {s} {c} {f}");
                }
            }
        }
    }
}

/// The rows of the vector file, `rows[(n_atoms * 11 + solvent) * 11 + solute][ndfmin + 3]`.
fn lean_rows() -> Vec<Vec<i64>> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../spec/vectors/dof.json");
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    text.lines()
        .map(|l| l.trim().trim_end_matches(','))
        .filter(|l| l.starts_with('[') && l.ends_with(']'))
        .map(|l| {
            l[1..l.len() - 1]
                .split(',')
                .map(|v| v.trim().parse().expect("integer in dof.json"))
                .collect()
        })
        .collect()
}

#[test]
fn rust_equals_lean_model_on_grid() {
    let rows = lean_rows();
    let side = GRID + 1;
    assert_eq!(
        rows.len() as u64,
        side * side * side,
        "vector file has {} rows",
        rows.len()
    );
    let mut mismatches = Vec::new();
    for n in 0..side {
        for s in 0..side {
            for c in 0..side {
                let row = &rows[((n * side + s) * side + c) as usize];
                assert_eq!(row.len(), NDFMINS.count(), "row ({n}, {s}, {c})");
                for (&lean, f) in row.iter().zip(NDFMINS) {
                    if total_dof(n, s, c, f) != lean {
                        mismatches.push((n, s, c, f));
                    }
                }
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "(n_atoms, solvent, solute, ndfmin) mismatches: {mismatches:?}"
    );
}
