//! What kind of nonbonded pair two atoms form.
//!
//! Atoms one or two bonds apart (1-2, 1-3) are excluded: no Lennard-Jones or Coulomb term from
//! the pairlist, only the reaction-field correction of excluded pairs. Atoms three bonds apart
//! (1-4) interact with the special 1-4 parameters, also off the pairlist. A `PERTATOMPAIR` or
//! EDS pair belongs to the perturbed atom's own loop and stays off the pairlist too. Every other
//! pair is a normal pair.
//!
//! Each atom `i` carries three partner lists, ascending and symmetric (`j` is in a list of `i`
//! exactly when `i` is in the same list of `j`), so one atom's lists decide the pair.
//!
//! Theory: `spec/NewenTheory/Topology/Pairs.lean` (`pairKind`, by membership); model, bridge and
//! Aeneas refinement: `spec/NewenSpec/Topology/Pairs.lean`, `spec/NewenSpec/Refinement.lean`.

/// The kind of the pair `(i, j)`, in priority order: an excluded pair is never also booked as a
/// 1-4 or an owned pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairKind {
    /// on the pairlist: full Lennard-Jones and Coulomb (or reaction-field) interaction
    Normal,
    /// 1-2 or 1-3 neighbour (gromosXX `exclusion`)
    Excluded,
    /// 1-4 neighbour (gromosXX `one_four_pair`), computed with the 1-4 parameters
    OneFour,
    /// `PERTATOMPAIR` or EDS pair, computed by the perturbed atom's own loop
    AtomPair,
}

/// Whether `j` is in the ascending list `list`. An empty list, or one whose last (largest) entry
/// is below `j`, answers at once (gromosXX `Exclusions::is_excluded`, `exclusions.h:79-90`):
/// exclusions are intramolecular, so a partner in another molecule rarely needs the scan.
///
/// spec: topology.pairs.contains_sorted
#[inline]
pub fn contains_sorted(list: &[usize], j: usize) -> bool {
    let n = list.len();
    if n == 0 {
        return false;
    }
    if list[n - 1] < j {
        return false;
    }
    let mut k = 0;
    while k < n {
        if list[k] == j {
            return true;
        }
        k += 1;
    }
    false
}

/// The kind of the pair `(i, j)` from atom `i`'s three ascending partner lists: `Excluded` if
/// `j` is in `exclusions_i`, else `OneFour` if in `one_four_i`, else `AtomPair` if in
/// `atom_pair_i`, else `Normal`.
///
/// spec: topology.pairs.pair_kind
#[inline]
pub fn pair_kind(
    exclusions_i: &[usize],
    one_four_i: &[usize],
    atom_pair_i: &[usize],
    j: usize,
) -> PairKind {
    if contains_sorted(exclusions_i, j) {
        return PairKind::Excluded;
    }
    if contains_sorted(one_four_i, j) {
        return PairKind::OneFour;
    }
    if contains_sorted(atom_pair_i, j) {
        return PairKind::AtomPair;
    }
    PairKind::Normal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn membership_with_the_early_out() {
        assert!(!contains_sorted(&[], 0));
        assert!(contains_sorted(&[1, 2], 2));
        assert!(!contains_sorted(&[1, 2], 3));
        assert!(!contains_sorted(&[1, 3], 2));
        assert!(!contains_sorted(&[1, 3], 0));
    }

    #[test]
    fn exclusion_wins_over_one_four_and_atom_pair() {
        assert_eq!(pair_kind(&[2], &[2], &[2], 2), PairKind::Excluded);
        assert_eq!(pair_kind(&[1], &[2], &[2], 2), PairKind::OneFour);
        assert_eq!(pair_kind(&[1], &[3], &[2], 2), PairKind::AtomPair);
        assert_eq!(pair_kind(&[1], &[3], &[4], 2), PairKind::Normal);
    }
}
