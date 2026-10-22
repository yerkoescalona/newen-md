//! Which group an atom belongs to, from a list of last atoms.
//!
//! GROMOS describes a partition of the atoms into consecutive groups by the last atom of each
//! group (FORCE `NRE` for energy groups, `MULTIBATH` DOFSET for temperature baths): group `g`
//! holds the atoms after group `g − 1`'s last atom up to its own. Energies are booked per pair of
//! energy groups, in an `n × n` matrix.
//!
//! Theory: `spec/NewenTheory/Topology/Groups.lean` (`groupOf`, `LastAtomsOk`); model, bridge and
//! Aeneas refinement: `spec/NewenSpec/Topology/Groups.lean`, `spec/NewenSpec/Refinement.lean`.

/// The group of `atom`: the first group whose last atom (0-based) reaches it; `None` past the
/// last group (gromosXX `topo.atom_energy_group(i)`, `topology.h:738-743`, filled by
/// `in_topology.cc:449-458`).
///
/// spec: topology.groups.group_of
pub fn group_of(last_atoms: &[usize], atom: usize) -> Option<usize> {
    let mut g = 0;
    while g < last_atoms.len() {
        if atom <= last_atoms[g] {
            return Some(g);
        }
        g += 1;
    }
    None
}

/// Outcome of [`check_last_atoms`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LastAtomsCheck {
    /// the groups partition the atoms
    Ok,
    /// no group at all
    Empty,
    /// entry `g` (0-based) does not exceed entry `g − 1`, so group `g` would be empty
    NotIncreasing(usize),
    /// the last group does not end at the last atom
    LastNotFinal,
}

/// The rule that makes a list of 0-based last atoms partition `n_atoms` atoms into non-empty
/// consecutive groups: at least one group, strictly increasing, the last entry `n_atoms − 1`
/// (so `n_atoms ≥ 1`). Checked in that order.
///
/// spec: topology.groups.check_last_atoms
pub fn check_last_atoms(last_atoms: &[usize], n_atoms: usize) -> LastAtomsCheck {
    let n = last_atoms.len();
    if n == 0 {
        return LastAtomsCheck::Empty;
    }
    // the first entry that does not exceed its predecessor, or `n`
    let mut g = 1;
    while g < n && last_atoms[g - 1] < last_atoms[g] {
        g += 1;
    }
    if g < n {
        return LastAtomsCheck::NotIncreasing(g);
    }
    if n_atoms == 0 {
        return LastAtomsCheck::LastNotFinal;
    }
    if last_atoms[n - 1] != n_atoms - 1 {
        return LastAtomsCheck::LastNotFinal;
    }
    LastAtomsCheck::Ok
}

/// Row-major cell `eg_i · n_groups + eg_j` of the `n_groups × n_groups` energy matrix, `None`
/// when a group index is out of range.
///
/// spec: topology.groups.energy_cell
#[inline]
pub fn energy_cell(eg_i: usize, eg_j: usize, n_groups: usize) -> Option<usize> {
    if eg_i < n_groups && eg_j < n_groups {
        return Some(eg_i * n_groups + eg_j);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_atom_belongs_to_the_first_group_that_reaches_it() {
        assert_eq!(group_of(&[1, 4], 0), Some(0));
        assert_eq!(group_of(&[1, 4], 2), Some(1));
        assert_eq!(group_of(&[1, 4], 5), None);
        assert_eq!(group_of(&[], 0), None);
    }

    #[test]
    fn the_nre_rule() {
        assert_eq!(check_last_atoms(&[1, 4], 5), LastAtomsCheck::Ok);
        assert_eq!(check_last_atoms(&[], 5), LastAtomsCheck::Empty);
        assert_eq!(
            check_last_atoms(&[1, 1, 4], 5),
            LastAtomsCheck::NotIncreasing(1)
        );
        assert_eq!(check_last_atoms(&[1, 3], 5), LastAtomsCheck::LastNotFinal);
        assert_eq!(check_last_atoms(&[0], 0), LastAtomsCheck::LastNotFinal);
    }

    #[test]
    fn cells_are_row_major() {
        assert_eq!(energy_cell(1, 2, 3), Some(5));
        assert_eq!(energy_cell(3, 0, 3), None);
        assert_eq!(energy_cell(0, 3, 3), None);
    }
}
