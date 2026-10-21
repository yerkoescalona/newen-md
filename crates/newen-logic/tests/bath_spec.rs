//! The per-bath degrees of freedom against their meaning and their Lean model.
//!
//! Three independent references:
//! - `by_atoms` below counts atom by atom without the range formula: it scans the ranges for each
//!   atom's baths, gives 6 halves to the COM bath of a group's first atom and to the IR bath of
//!   every other atom, takes 1 per constraint end and 6 from each bath of atom 0 for
//!   roto-translational constraints; the spread is checked by 128-bit cross-multiplication;
//! - `spec/vectors/bath.json` holds the Lean models on a family of small systems, written by
//!   `lake exe spec-vectors` (checked in, so this test needs no Lean);
//! - the two DOFSET splits of the gromosXX reference `aladip_vacuum_multibath_ndfmin3`.

use newen_logic::thermo::bath::{
    bath_half_dofs, baths_of, ndfmin_spread, Constraint, DofRange, HalfDof, Spread,
};

/// The baths of `atom`: those of the first range reaching it.
fn bath_of_atom(ranges: &[DofRange], atom: usize) -> Option<(usize, usize)> {
    for r in ranges {
        if r.last_atom >= atom {
            return Some((r.com_bath, r.ir_bath));
        }
    }
    None
}

/// Atom-by-atom half degrees of freedom of every bath.
fn by_atoms(
    n_atoms: usize,
    n_baths: usize,
    ranges: &[DofRange],
    group_firsts: &[usize],
    constraints: &[Constraint],
    rottrans: bool,
) -> Vec<HalfDof> {
    let mut halves = vec![HalfDof { com: 0, ir: 0 }; n_baths];
    for atom in 0..n_atoms {
        let (com, ir) = bath_of_atom(ranges, atom).expect("ranges cover the atoms");
        if group_firsts.contains(&atom) {
            halves[com].com += 6;
        } else {
            halves[ir].ir += 6;
        }
    }
    for c in constraints {
        for atom in [c.i, c.j] {
            if let Some((_, ir)) = bath_of_atom(ranges, atom) {
                halves[ir].ir -= 1;
            }
        }
    }
    if rottrans && n_atoms > 0 {
        let (com, ir) = bath_of_atom(ranges, 0).expect("atom 0 has baths");
        halves[com].com -= 6;
        halves[ir].ir -= 6;
    }
    halves
}

/// The spread means: bath b keeps `(h_b / 2) · (1 − 2·ndfmin / H)`, all baths together
/// `H/2 − ndfmin`; with `H ≤ 0` each keeps `h_b / 2`.
fn check_spread(h: &[i64], ndfmin: i64, s: &Spread) -> Result<(), String> {
    let total: i128 = h.iter().map(|&x| x as i128).sum();
    let den = s.denominator as i128;
    if s.numerators.len() != h.len() || den <= 0 {
        return Err(format!("shape {s:?}"));
    }
    if total <= 0 {
        let unchanged = s.numerators == h && s.denominator == 2;
        return if unchanged {
            Ok(())
        } else {
            Err(format!("H = {total}: {s:?}"))
        };
    }
    let kept = total - 2 * ndfmin as i128;
    for (&hb, &num) in h.iter().zip(&s.numerators) {
        // num / den == hb · kept / (2 H)
        if num as i128 * 2 * total != hb as i128 * kept * den {
            return Err(format!("bath {hb}: {s:?}"));
        }
    }
    let sum: i128 = s.numerators.iter().map(|&x| x as i128).sum();
    // Σ num / den == H/2 − ndfmin
    if 2 * sum != (total - 2 * ndfmin as i128) * den {
        return Err(format!("sum: {s:?}"));
    }
    Ok(())
}

/// Ascending subsets of `xs` with at most `k` elements.
fn subsets<T: Clone>(xs: &[T], k: usize) -> Vec<Vec<T>> {
    let mut out = vec![vec![]];
    if k == 0 || xs.is_empty() {
        return out;
    }
    for rest in subsets(&xs[1..], k - 1) {
        let mut with = vec![xs[0].clone()];
        with.extend(rest);
        out.push(with);
    }
    out.extend(subsets(&xs[1..], k).into_iter().skip(1));
    out
}

/// Every bath assignment of `n_ranges` ranges over `n_baths` baths, at most about `cap` of them
/// (evenly strided).
fn assignments(n_ranges: usize, n_baths: usize, cap: usize) -> Vec<Vec<(usize, usize)>> {
    let total = (n_baths * n_baths).pow(n_ranges as u32);
    let stride = total.div_ceil(cap).max(1);
    (0..total)
        .step_by(stride)
        .map(|mut code| {
            (0..n_ranges)
                .map(|_| {
                    let pair = (code % n_baths, (code / n_baths) % n_baths);
                    code /= n_baths * n_baths;
                    pair
                })
                .collect()
        })
        .collect()
}

#[test]
fn rust_equals_atom_count_on_grid() {
    let mut checked = 0usize;
    let mut half_integer = 0usize;
    let mut nonpositive = 0usize;
    for n in 1..=6usize {
        let inner: Vec<usize> = (0..n - 1).collect();
        let pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
            .collect();
        let constraint_sets = subsets(&pairs, 2);
        let stride = constraint_sets.len().div_ceil(25).max(1);
        let group_tails: Vec<usize> = (1..n).collect();
        for split in subsets(&inner, 2) {
            let mut lasts = split.clone();
            lasts.push(n - 1);
            for n_baths in 1..=3 {
                for assignment in assignments(lasts.len(), n_baths, 81) {
                    let ranges: Vec<DofRange> = lasts
                        .iter()
                        .zip(&assignment)
                        .map(|(&last_atom, &(com_bath, ir_bath))| DofRange {
                            last_atom,
                            com_bath,
                            ir_bath,
                        })
                        .collect();
                    for tail in subsets(&group_tails, 2) {
                        let mut groups = vec![0];
                        groups.extend(tail);
                        for cs in constraint_sets.iter().step_by(stride) {
                            let constraints: Vec<Constraint> =
                                cs.iter().map(|&(i, j)| Constraint { i, j }).collect();
                            for rottrans in [false, true] {
                                let got = bath_half_dofs(
                                    n_baths,
                                    &ranges,
                                    &groups,
                                    &constraints,
                                    rottrans,
                                );
                                let want =
                                    by_atoms(n, n_baths, &ranges, &groups, &constraints, rottrans);
                                assert_eq!(
                                    got, want,
                                    "n {n} ranges {ranges:?} groups {groups:?} \
                                     constraints {cs:?} rottrans {rottrans}"
                                );
                                let h: Vec<i64> = got.iter().map(|b| b.com + b.ir).collect();
                                if h.iter().any(|x| x % 2 != 0) {
                                    half_integer += 1;
                                }
                                if h.iter().sum::<i64>() <= 0 {
                                    nonpositive += 1;
                                }
                                for ndfmin in [0, 3, 6] {
                                    let s = ndfmin_spread(&h, ndfmin);
                                    if let Err(e) = check_spread(&h, ndfmin, &s) {
                                        panic!("{h:?} ndfmin {ndfmin}: {e}");
                                    }
                                }
                                checked += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(checked > 100_000, "only {checked} systems");
    assert!(half_integer > 0 && nonpositive > 0);
}

/// A vector row, split into its fields by the layout in `bath.json`.
struct Row {
    n_baths: usize,
    rottrans: bool,
    ndfmin: i64,
    ranges: Vec<DofRange>,
    groups: Vec<usize>,
    constraints: Vec<Constraint>,
    halves: Vec<HalfDof>,
    spread: Spread,
}

fn parse_row(v: &[i64]) -> Row {
    let mut it = v.iter().copied();
    let mut next = || it.next().expect("row too short");
    let _n_atoms = next();
    let n_baths = next() as usize;
    let rottrans = next() == 1;
    let ndfmin = next();
    let n_ranges = next() as usize;
    let ranges = (0..n_ranges)
        .map(|_| DofRange {
            last_atom: next() as usize,
            com_bath: next() as usize,
            ir_bath: next() as usize,
        })
        .collect();
    let n_groups = next() as usize;
    let groups = (0..n_groups).map(|_| next() as usize).collect();
    let n_constraints = next() as usize;
    let constraints = (0..n_constraints)
        .map(|_| Constraint {
            i: next() as usize,
            j: next() as usize,
        })
        .collect();
    let halves = (0..n_baths)
        .map(|_| HalfDof {
            com: next(),
            ir: next(),
        })
        .collect();
    let numerators = (0..n_baths).map(|_| next()).collect();
    let denominator = next();
    assert!(it.next().is_none(), "row too long");
    Row {
        n_baths,
        rottrans,
        ndfmin,
        ranges,
        groups,
        constraints,
        halves,
        spread: Spread {
            numerators,
            denominator,
        },
    }
}

fn lean_rows() -> Vec<Vec<i64>> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../spec/vectors/bath.json");
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("{path}: {e} (regenerate with `cd spec && lake exe spec-vectors`)")
    });
    text.lines()
        .map(|l| l.trim().trim_end_matches(','))
        .filter(|l| l.starts_with('[') && l.ends_with(']'))
        .map(|l| {
            l[1..l.len() - 1]
                .split(',')
                .map(|v| v.trim().parse().expect("integer in bath.json"))
                .collect()
        })
        .collect()
}

#[test]
fn rust_equals_lean_model_on_vectors() {
    let rows = lean_rows();
    assert!(rows.len() >= 1000, "vector file has {} rows", rows.len());
    let mut half_integer = 0usize;
    let mut nonpositive = 0usize;
    for v in &rows {
        let row = parse_row(v);
        let got = bath_half_dofs(
            row.n_baths,
            &row.ranges,
            &row.groups,
            &row.constraints,
            row.rottrans,
        );
        assert_eq!(got, row.halves, "row {v:?}");
        let h: Vec<i64> = got.iter().map(|b| b.com + b.ir).collect();
        assert_eq!(ndfmin_spread(&h, row.ndfmin), row.spread, "row {v:?}");
        if h.iter().any(|x| x % 2 != 0) {
            half_integer += 1;
        }
        if h.iter().sum::<i64>() <= 0 {
            nonpositive += 1;
        }
    }
    assert!(half_integer > 0, "no half-integer bath in the vectors");
    assert!(nonpositive > 0, "no H <= 0 case in the vectors");
}

/// The 11 bond constraints of alanine dipeptide (`aladip_tg2.topo`, BONDH and BOND, 0-based).
const ALADIP_BONDS: [(usize, usize); 11] = [
    (3, 4),
    (9, 10),
    (0, 1),
    (1, 2),
    (1, 3),
    (3, 5),
    (5, 6),
    (5, 7),
    (7, 8),
    (7, 9),
    (9, 11),
];

/// 12 atoms, NTC 3, two baths split after atom `split` (1-based last atom of bath 1), one
/// temperature group per bath, NDFMIN 3: the halves and the kept dof as hundredths.
fn aladip(split: usize) -> (Vec<i64>, Vec<i64>) {
    let ranges = [
        DofRange {
            last_atom: split - 1,
            com_bath: 0,
            ir_bath: 0,
        },
        DofRange {
            last_atom: 11,
            com_bath: 1,
            ir_bath: 1,
        },
    ];
    let constraints: Vec<Constraint> = ALADIP_BONDS
        .iter()
        .map(|&(i, j)| Constraint { i, j })
        .collect();
    let halves = bath_half_dofs(2, &ranges, &[0, split], &constraints, false);
    let h: Vec<i64> = halves.iter().map(|b| b.com + b.ir).collect();
    let s = ndfmin_spread(&h, 3);
    let hundredths = s
        .numerators
        .iter()
        .map(|&num| {
            assert_eq!(num * 100 % s.denominator, 0);
            num * 100 / s.denominator
        })
        .collect();
    (h, hundredths)
}

#[test]
fn gromosxx_aladip_split_6_7() {
    assert_eq!(aladip(6), (vec![24, 26], vec![1056, 1144]));
}

#[test]
fn gromosxx_aladip_split_7_8() {
    assert_eq!(aladip(7), (vec![29, 21], vec![1276, 924]));
}

#[test]
fn nothing_to_spread_over_keeps_the_halves() {
    for h in [vec![0, 0], vec![-6, 0, 4], vec![-1]] {
        assert_eq!(
            ndfmin_spread(&h, 3),
            Spread {
                numerators: h.clone(),
                denominator: 2
            }
        );
    }
}

#[test]
fn an_atom_beyond_the_last_range_has_no_bath() {
    let ranges = [DofRange {
        last_atom: 4,
        com_bath: 0,
        ir_bath: 1,
    }];
    assert_eq!(baths_of(&ranges, 4), Some((0, 1)));
    assert_eq!(baths_of(&ranges, 5), None);
}
