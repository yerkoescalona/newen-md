//! The atom → temperature-bath mapping of a `MULTIBATH` block (gromosXX `simulation/multibath.h`).
//!
//! The kinetic energy per bath, its centre-of-mass / internal split, the kinetic dH/dλ, every
//! thermostat's scaling and the degrees of freedom all read this one mapping.

/// One DOFSET line of a `MULTIBATH` block (gromosXX `bath_index_struct`): the atoms after the
/// previous line's `last_atom` up to this one couple their centre-of-mass motion to `com_bath`
/// and their internal and rotational motion to `ir_bath`. Lines are in file order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct BathRange {
    /// last atom of the range (0-based, inclusive)
    pub last_atom: usize,
    /// bath index that takes the centre-of-mass degrees of freedom of the range
    pub com_bath: usize,
    /// bath index that takes the internal/rotational degrees of freedom
    pub ir_bath: usize,
}

/// gromosXX `Multibath::in_bath` (`multibath.h:159-171`): the `(com, ir)` baths of `atom`, those
/// of the first range whose `last_atom` reaches it. Bath 0 when no range does (gromosXX asserts;
/// plan validation guarantees the last range ends at the last atom).
pub fn in_bath(ranges: &[BathRange], atom: usize) -> (usize, usize) {
    ranges
        .iter()
        .find(|r| r.last_atom >= atom)
        .map_or((0, 0), |r| (r.com_bath, r.ir_bath))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_atom_belongs_to_the_first_range_that_reaches_it() {
        // DOFSET `12 2 2  72 1 1`: the solute couples to bath 2, the solvent to bath 1
        let ranges = [
            BathRange {
                last_atom: 11,
                com_bath: 1,
                ir_bath: 1,
            },
            BathRange {
                last_atom: 71,
                com_bath: 0,
                ir_bath: 0,
            },
        ];
        assert_eq!(in_bath(&ranges, 0), (1, 1));
        assert_eq!(in_bath(&ranges, 11), (1, 1));
        assert_eq!(in_bath(&ranges, 12), (0, 0));
        assert_eq!(in_bath(&ranges, 71), (0, 0));
    }
}
