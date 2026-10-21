//! Degrees of freedom per temperature bath, and which bath an atom belongs to.
//!
//! A `MULTIBATH` block couples every atom to two baths: the centre-of-mass translation of its
//! temperature group to a COM bath, everything else (internal motion and rotation) to an IR bath.
//! The counts are kept in half degrees of freedom, so that a constraint between two baths (½ to
//! each) stays an integer, and the NDFMIN spread is an exact rational.
//!
//! Theory: `spec/NewenTheory/Dynamics/DegreesOfFreedom.lean` (section "Temperature baths"); model,
//! bridge and Aeneas refinement: `spec/NewenSpec/Thermo/Bath.lean`, `spec/NewenSpec/Refinement.lean`.
//!
//! Not modelled: virtual atoms (gromosXX `multibath.cc:95-104`), position, angle and dihedral
//! constraints (`topology.cc:1062-1070, 1230-1325`); newen refuses NTPOR 3 and has no angle or
//! dihedral constraints.

/// One DOFSET line of a `MULTIBATH` block (gromosXX `bath_index_struct`): the atoms after the
/// previous line's `last_atom` up to this one couple their centre-of-mass motion to `com_bath`
/// and their internal and rotational motion to `ir_bath`. Lines are in file order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct DofRange {
    /// last atom of the range (0-based, inclusive)
    pub last_atom: usize,
    /// bath index that takes the centre-of-mass degrees of freedom of the range
    pub com_bath: usize,
    /// bath index that takes the internal/rotational degrees of freedom
    pub ir_bath: usize,
}

/// One distance constraint between two atoms (global 0-based indices).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Constraint {
    pub i: usize,
    pub j: usize,
}

/// Half degrees of freedom of one bath before NDFMIN: `com + ir` is twice its dof.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HalfDof {
    /// twice the centre-of-mass translational degrees of freedom
    pub com: i64,
    /// twice the internal and rotational degrees of freedom
    pub ir: i64,
}

/// NDFMIN spread as exact rationals: bath `b` keeps `numerators[b] / denominator` degrees of
/// freedom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Spread {
    pub numerators: Vec<i64>,
    pub denominator: i64,
}

/// The `(com, ir)` baths of `atom`: those of the first range whose `last_atom >= atom`; `None`
/// when no range reaches it (gromosXX asserts there; plan validation guarantees coverage).
///
/// spec: thermo.bath.baths_of
pub fn baths_of(ranges: &[DofRange], atom: usize) -> Option<(usize, usize)> {
    let mut r = 0;
    while r < ranges.len() {
        if atom <= ranges[r].last_atom {
            return Some((ranges[r].com_bath, ranges[r].ir_bath));
        }
        r += 1;
    }
    None
}

/// How many of `group_firsts` (the first atom of every temperature group) lie in `first..=last`.
/// 0 when `first > last`.
///
/// spec: thermo.bath.groups_starting_in
pub fn groups_starting_in(group_firsts: &[usize], first: usize, last: usize) -> u64 {
    let mut count: u64 = 0;
    let mut g = 0;
    while g < group_firsts.len() {
        if first <= group_firsts[g] && group_firsts[g] <= last {
            count += 1;
        }
        g += 1;
    }
    count
}

/// Half degrees of freedom of bath `bath`.
///
/// Every atom has 3 degrees of freedom, 6 halves. A temperature group's centre-of-mass
/// translation (6 halves) goes to the COM bath of its range, the rest of the range's atoms to
/// its IR bath. A distance constraint removes one degree of freedom, half from the IR bath of
/// each end. Roto-translational constraints remove 3 from each bath of atom 0.
///
/// With `start_r` = previous range's `last_atom + 1` (0 for the first range), `n_r` = atoms in
/// range r (`last_atom − start_r + 1`, 0 when the range is empty or backwards) and
/// `g_r = groups_starting_in(group_firsts, start_r, last_atom)`:
///
/// ```text
/// com = Σ_{r : com_bath == bath} 6·g_r            − 6·[rottrans ∧ baths_of(0).com == bath]
/// ir  = Σ_{r : ir_bath == bath} (6·n_r − 6·g_r)   − Σ_c ([baths_of(c.i).ir == bath] + [baths_of(c.j).ir == bath])
///                                                 − 6·[rottrans ∧ baths_of(0).ir == bath]
/// ```
///
/// A constraint end whose atom no range reaches contributes nothing. The Lean refinement holds
/// while the counts fit in `i64` and no `last_atom` is `usize::MAX` (`BathHalfDofInRange`).
///
/// spec: thermo.bath.bath_half_dof
pub fn bath_half_dof(
    bath: usize,
    ranges: &[DofRange],
    group_firsts: &[usize],
    constraints: &[Constraint],
    rottrans: bool,
) -> HalfDof {
    let mut com: i64 = 0;
    let mut ir: i64 = 0;
    let mut start: usize = 0;
    let mut r = 0;
    while r < ranges.len() {
        let last = ranges[r].last_atom;
        let mut atoms: i64 = 0;
        if start <= last {
            atoms = (last - start + 1) as i64;
        }
        let groups = groups_starting_in(group_firsts, start, last) as i64;
        if ranges[r].com_bath == bath {
            com += 6 * groups;
        }
        if ranges[r].ir_bath == bath {
            ir += 6 * atoms - 6 * groups;
        }
        start = last + 1;
        r += 1;
    }
    let mut ends: i64 = 0;
    let mut c = 0;
    while c < constraints.len() {
        ends += end_in_bath(ranges, constraints[c].i, bath);
        ends += end_in_bath(ranges, constraints[c].j, bath);
        c += 1;
    }
    ir -= ends;
    if rottrans {
        if let Some((com_0, ir_0)) = baths_of(ranges, 0) {
            if com_0 == bath {
                com -= 6;
            }
            if ir_0 == bath {
                ir -= 6;
            }
        }
    }
    HalfDof { com, ir }
}

/// 1 if `atom` is reached by a range whose IR bath is `bath`, else 0: one constraint end.
fn end_in_bath(ranges: &[DofRange], atom: usize, bath: usize) -> i64 {
    if let Some((_, ir)) = baths_of(ranges, atom) {
        if ir == bath {
            return 1;
        }
    }
    0
}

/// [`bath_half_dof`] for every bath `0..n_baths`, in order.
///
/// spec: thermo.bath.bath_half_dofs
pub fn bath_half_dofs(
    n_baths: usize,
    ranges: &[DofRange],
    group_firsts: &[usize],
    constraints: &[Constraint],
    rottrans: bool,
) -> Vec<HalfDof> {
    let mut halves = Vec::new();
    let mut b = 0;
    while b < n_baths {
        halves.push(bath_half_dof(
            b,
            ranges,
            group_firsts,
            constraints,
            rottrans,
        ));
        b += 1;
    }
    halves
}

/// NDFMIN removed from every bath in proportion to its degrees of freedom.
///
/// With `H = Σ half_dofs`: if `H > 0`, `numerators[b] = half_dofs[b] · (H − 2·ndfmin)` and
/// `denominator = 2·H`, so bath b keeps `(h_b/2)·(1 − 2·ndfmin/H)` and all baths together lose
/// exactly `ndfmin`. Otherwise there is nothing to spread over: `numerators = half_dofs`,
/// `denominator = 2`. gromosXX sums H in an `int`, truncating half-integer baths (D40); here it
/// is exact. The Lean refinement holds while `(Σ|h_b| + 2·|ndfmin| + 2)² ≤ i64::MAX`.
///
/// spec: thermo.bath.ndfmin_spread
pub fn ndfmin_spread(half_dofs: &[i64], ndfmin: i64) -> Spread {
    let mut total: i64 = 0;
    let mut b = 0;
    while b < half_dofs.len() {
        total += half_dofs[b];
        b += 1;
    }
    let mut factor: i64 = 1;
    let mut denominator: i64 = 2;
    if total > 0 {
        factor = total - 2 * ndfmin;
        denominator = 2 * total;
    }
    let mut numerators = Vec::new();
    let mut k = 0;
    while k < half_dofs.len() {
        numerators.push(half_dofs[k] * factor);
        k += 1;
    }
    Spread {
        numerators,
        denominator,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn range(last_atom: usize, com_bath: usize, ir_bath: usize) -> DofRange {
        DofRange {
            last_atom,
            com_bath,
            ir_bath,
        }
    }

    #[test]
    fn an_atom_belongs_to_the_first_range_that_reaches_it() {
        let ranges = [range(1, 0, 1), range(4, 2, 2)];
        assert_eq!(baths_of(&ranges, 0), Some((0, 1)));
        assert_eq!(baths_of(&ranges, 2), Some((2, 2)));
        assert_eq!(baths_of(&ranges, 5), None);
        assert_eq!(baths_of(&[], 0), None);
    }

    #[test]
    fn groups_are_counted_inside_the_interval() {
        assert_eq!(groups_starting_in(&[0, 2, 5], 1, 5), 2);
        assert_eq!(groups_starting_in(&[0, 2, 5], 3, 2), 0);
    }

    #[test]
    fn one_bath_has_six_halves_per_atom_less_constraints() {
        // three atoms, one group, one constraint: 3·3 − 1 = 8 dof = 16 halves, 6 of them COM
        let h = bath_half_dof(
            0,
            &[range(2, 0, 0)],
            &[0],
            &[Constraint { i: 0, j: 1 }],
            false,
        );
        assert_eq!(h, HalfDof { com: 6, ir: 10 });
        let r = bath_half_dof(0, &[range(2, 0, 0)], &[0], &[], true);
        assert_eq!(r, HalfDof { com: 0, ir: 6 });
    }

    #[test]
    fn a_constraint_across_two_baths_gives_half_to_each() {
        let ranges = [range(0, 0, 0), range(1, 1, 1)];
        let halves = bath_half_dofs(2, &ranges, &[0, 1], &[Constraint { i: 0, j: 1 }], false);
        assert_eq!(
            halves,
            vec![HalfDof { com: 6, ir: -1 }, HalfDof { com: 6, ir: -1 }]
        );
    }

    #[test]
    fn ndfmin_is_spread_in_proportion() {
        // halves 24 and 26 (12 and 13 dof), NDFMIN 3: 12 − 3·12/25 and 13 − 3·13/25
        let s = ndfmin_spread(&[24, 26], 3);
        assert_eq!(s.denominator, 100);
        assert_eq!(s.numerators, vec![24 * 44, 26 * 44]);
    }

    #[test]
    fn nothing_to_spread_over() {
        let s = ndfmin_spread(&[0, -2], 3);
        assert_eq!(
            s,
            Spread {
                numerators: vec![0, -2],
                denominator: 2
            }
        );
    }
}
