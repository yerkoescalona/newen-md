//! Pair kinds and group lookups against their meaning and their Lean model.
//!
//! Three independent references:
//! - brute force below: plain membership (`slice::contains`) in priority order for the pair
//!   kind, the least index (`Iterator::position`) for the group, adjacent windows for the NRE
//!   rule, and every cell of the matrix for `energy_cell`;
//! - `spec/vectors/topology.json` holds the Lean models on enumerated small lists, written by
//!   `lake exe spec-vectors` (checked in, so this test needs no Lean);
//! - hand cases from the gromosXX references: the butane topology's exclusions, 1-4 pairs and
//!   charge groups, and three FORCE `NRE` lines.

use newen_logic::topology::groups::{check_last_atoms, energy_cell, group_of, LastAtomsCheck};
use newen_logic::topology::pairs::{contains_sorted, pair_kind, PairKind};

/// Every ascending subset of `0..n`, by bitmask.
fn subsets(n: usize) -> Vec<Vec<usize>> {
    (0..1usize << n)
        .map(|m| (0..n).filter(|b| m >> b & 1 == 1).collect())
        .collect()
}

/// Every list of at most `k` entries from `0..vals`, any order, repeats allowed.
fn lists(vals: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = vec![vec![]];
    let mut last = vec![vec![]];
    for _ in 0..k {
        last = last
            .iter()
            .flat_map(|t: &Vec<usize>| {
                (0..vals).map(move |v| {
                    let mut t = t.clone();
                    t.push(v);
                    t
                })
            })
            .collect();
        out.extend(last.iter().cloned());
    }
    out
}

fn kind_by_membership(e: &[usize], f: &[usize], a: &[usize], j: usize) -> PairKind {
    if e.contains(&j) {
        PairKind::Excluded
    } else if f.contains(&j) {
        PairKind::OneFour
    } else if a.contains(&j) {
        PairKind::AtomPair
    } else {
        PairKind::Normal
    }
}

fn group_by_position(last_atoms: &[usize], atom: usize) -> Option<usize> {
    last_atoms.iter().position(|&l| atom <= l)
}

fn check_by_windows(last_atoms: &[usize], n_atoms: usize) -> LastAtomsCheck {
    if last_atoms.is_empty() {
        return LastAtomsCheck::Empty;
    }
    if let Some(w) = last_atoms.windows(2).position(|w| w[0] >= w[1]) {
        return LastAtomsCheck::NotIncreasing(w + 1);
    }
    if n_atoms == 0 || *last_atoms.last().unwrap() != n_atoms - 1 {
        return LastAtomsCheck::LastNotFinal;
    }
    LastAtomsCheck::Ok
}

#[test]
fn contains_sorted_is_membership_on_every_ascending_list() {
    for list in subsets(10) {
        for j in 0..=11 {
            assert_eq!(contains_sorted(&list, j), list.contains(&j), "{list:?} {j}");
        }
    }
}

#[test]
fn pair_kind_is_membership_by_priority() {
    let subs = subsets(5);
    let mut seen = [0usize; 4];
    for e in &subs {
        for f in &subs {
            for a in &subs {
                for j in 0..=6 {
                    let k = pair_kind(e, f, a, j);
                    assert_eq!(k, kind_by_membership(e, f, a, j), "{e:?} {f:?} {a:?} {j}");
                    seen[k as usize] += 1;
                }
            }
        }
    }
    assert!(seen.iter().all(|&c| c > 0), "{seen:?}");
}

#[test]
fn group_of_is_the_least_index_and_the_nre_rule_partitions() {
    let mut ok = 0usize;
    for last_atoms in lists(7, 4) {
        for n_atoms in 0..=8 {
            let check = check_last_atoms(&last_atoms, n_atoms);
            assert_eq!(
                check,
                check_by_windows(&last_atoms, n_atoms),
                "{last_atoms:?} {n_atoms}"
            );
            for atom in 0..=8 {
                assert_eq!(
                    group_of(&last_atoms, atom),
                    group_by_position(&last_atoms, atom)
                );
            }
            if check == LastAtomsCheck::Ok {
                ok += 1;
                // every atom in exactly one non-empty group, in order; nothing past the atoms
                let groups: Vec<usize> = (0..n_atoms)
                    .map(|a| group_of(&last_atoms, a).expect("an atom without a group"))
                    .collect();
                assert_eq!(groups[0], 0);
                assert!(groups.windows(2).all(|w| w[1] == w[0] || w[1] == w[0] + 1));
                assert_eq!(*groups.last().unwrap(), last_atoms.len() - 1);
                assert_eq!(group_of(&last_atoms, n_atoms), None);
            }
        }
    }
    assert!(ok > 90, "only {ok} partitions");
}

#[test]
fn energy_cells_are_row_major_and_distinct() {
    for n in 0..=6usize {
        let mut cells = vec![false; n * n];
        for i in 0..=7 {
            for j in 0..=7 {
                let cell = energy_cell(i, j, n);
                if i < n && j < n {
                    let c = cell.expect("a cell in range");
                    assert_eq!(c, i * n + j);
                    assert!(!std::mem::replace(&mut cells[c], true), "cell {c} twice");
                } else {
                    assert_eq!(cell, None, "({i}, {j}) of {n}");
                }
            }
        }
        assert!(cells.iter().all(|&c| c));
    }
}

fn lean_rows() -> Vec<Vec<i64>> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../spec/vectors/topology.json"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    text.lines()
        .map(|l| l.trim().trim_end_matches(','))
        .filter(|l| l.starts_with('[') && l.ends_with(']'))
        .map(|l| {
            l[1..l.len() - 1]
                .split(',')
                .map(|v| v.trim().parse().expect("integer in topology.json"))
                .collect()
        })
        .collect()
}

/// Reads a length-prefixed list from the front of `it`.
fn take_list(it: &mut std::slice::Iter<'_, i64>) -> Vec<usize> {
    let n = *it.next().expect("row too short") as usize;
    (0..n)
        .map(|_| *it.next().expect("row too short") as usize)
        .collect()
}

fn kind_code(k: PairKind) -> i64 {
    match k {
        PairKind::Normal => 0,
        PairKind::Excluded => 1,
        PairKind::OneFour => 2,
        PairKind::AtomPair => 3,
    }
}

fn check_code(c: LastAtomsCheck) -> [i64; 2] {
    match c {
        LastAtomsCheck::Ok => [0, 0],
        LastAtomsCheck::Empty => [1, 0],
        LastAtomsCheck::NotIncreasing(g) => [2, g as i64],
        LastAtomsCheck::LastNotFinal => [3, 0],
    }
}

fn opt_code(o: Option<usize>) -> i64 {
    o.map_or(-1, |g| g as i64)
}

#[test]
fn rust_equals_lean_model_on_vectors() {
    let rows = lean_rows();
    let mut per_tag = [0usize; 4];
    let mut not_increasing = 0usize;
    for row in &rows {
        let mut it = row[1..].iter();
        match row[0] {
            0 => {
                let list = take_list(&mut it);
                let got: Vec<i64> = (0..it.len())
                    .map(|j| contains_sorted(&list, j) as i64)
                    .collect();
                assert_eq!(got, it.as_slice(), "row {row:?}");
            },
            1 => {
                let e = take_list(&mut it);
                let f = take_list(&mut it);
                let a = take_list(&mut it);
                let got: Vec<i64> = (0..it.len())
                    .map(|j| kind_code(pair_kind(&e, &f, &a, j)))
                    .collect();
                assert_eq!(got, it.as_slice(), "row {row:?}");
            },
            2 => {
                let n_atoms = *it.next().unwrap() as usize;
                let last_atoms = take_list(&mut it);
                let check = check_last_atoms(&last_atoms, n_atoms);
                if matches!(check, LastAtomsCheck::NotIncreasing(_)) {
                    not_increasing += 1;
                }
                let mut got = check_code(check).to_vec();
                got.extend((0..=n_atoms).map(|a| opt_code(group_of(&last_atoms, a))));
                assert_eq!(got, it.as_slice(), "row {row:?}");
            },
            3 => {
                let v: Vec<usize> = it.by_ref().take(3).map(|&x| x as usize).collect();
                let want = *it.next().unwrap();
                assert_eq!(opt_code(energy_cell(v[0], v[1], v[2])), want, "row {row:?}");
                assert!(it.next().is_none(), "row {row:?} too long");
            },
            t => panic!("unknown tag {t}"),
        }
        per_tag[row[0] as usize] += 1;
    }
    assert!(
        per_tag.iter().all(|&c| c >= 100),
        "rows per tag {per_tag:?}"
    );
    assert!(not_increasing > 0, "no NotIncreasing case in the vectors");
}

/// `butane_vacuum.topo`, 0-based and symmetric as newen stores it: atom 1 excludes 2 and 3 and
/// has the 1-4 partner 4; atom 2 excludes 3 and 4; atom 3 excludes 4.
#[test]
fn butane_pairs() {
    let exclusions: [&[usize]; 4] = [&[1, 2], &[0, 2, 3], &[0, 1, 3], &[1, 2]];
    let one_four: [&[usize]; 4] = [&[3], &[], &[], &[0]];
    for i in 0..4 {
        for j in 0..4 {
            if i == j {
                continue;
            }
            let k = pair_kind(exclusions[i], one_four[i], &[], j);
            assert_eq!(k, pair_kind(exclusions[j], one_four[j], &[], i));
            let want = if (i, j) == (0, 3) || (i, j) == (3, 0) {
                PairKind::OneFour
            } else {
                PairKind::Excluded
            };
            assert_eq!(k, want, "({i}, {j})");
        }
    }
    // a partner in another molecule (atom 4 of a second butane) is normal
    assert_eq!(
        pair_kind(exclusions[0], one_four[0], &[], 4),
        PairKind::Normal
    );
}

/// Butane's charge groups, CGC `1 0 0 1`: atom 1 alone, atoms 2-4 together.
#[test]
fn butane_charge_groups() {
    let last_atoms = [0, 3];
    assert_eq!(check_last_atoms(&last_atoms, 4), LastAtomsCheck::Ok);
    let groups: Vec<Option<usize>> = (0..5).map(|a| group_of(&last_atoms, a)).collect();
    assert_eq!(groups, [Some(0), Some(1), Some(1), Some(1), None]);
}

/// FORCE `NRE` of three gromosXX references (1-based in the file).
#[test]
fn force_nre_lines() {
    // pair_lj: `2 1 2`, two atoms
    assert_eq!(check_last_atoms(&[0, 1], 2), LastAtomsCheck::Ok);
    assert_eq!(group_of(&[0, 1], 1), Some(1));
    // aladip_harmonic_covalent_groups: `3 4 12 72`, 12 solute atoms and 20 waters
    let nre = [3, 11, 71];
    assert_eq!(check_last_atoms(&nre, 72), LastAtomsCheck::Ok);
    let eg: Vec<usize> = (0..72).map(|a| group_of(&nre, a).unwrap()).collect();
    assert_eq!((eg[3], eg[4], eg[11], eg[12], eg[71]), (0, 1, 1, 2, 2));
    assert_eq!(energy_cell(eg[3], eg[71], 3), Some(2));
    assert_eq!(energy_cell(eg[71], eg[3], 3), Some(6));
    // aladip_vacuum_fep_groups: `2 6 12`, 12 atoms
    assert_eq!(check_last_atoms(&[5, 11], 12), LastAtomsCheck::Ok);
    assert_eq!(
        (group_of(&[5, 11], 5), group_of(&[5, 11], 6)),
        (Some(0), Some(1))
    );
    // the same line against too few or too many atoms
    assert_eq!(check_last_atoms(&[5, 11], 11), LastAtomsCheck::LastNotFinal);
    assert_eq!(check_last_atoms(&[5, 11], 13), LastAtomsCheck::LastNotFinal);
}
