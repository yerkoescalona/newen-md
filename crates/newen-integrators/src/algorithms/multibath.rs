//! The atom → temperature-bath mapping of a `MULTIBATH` block (gromosXX `simulation/multibath.h`).
//!
//! The kinetic energy per bath, its centre-of-mass / internal split, the kinetic dH/dλ, every
//! thermostat's scaling and the degrees of freedom all read this one mapping. It is defined, with
//! its Lean proof, in `newen_logic::thermo::bath`.

pub use newen_logic::thermo::bath::DofRange as BathRange;

/// gromosXX `Multibath::in_bath` (`multibath.h:159-171`): the `(com, ir)` baths of `atom`, those
/// of the first range whose `last_atom` reaches it. Bath 0 when no range does (gromosXX asserts;
/// plan validation guarantees the last range ends at the last atom).
pub fn in_bath(ranges: &[BathRange], atom: usize) -> (usize, usize) {
    newen_logic::thermo::bath::baths_of(ranges, atom).unwrap_or((0, 0))
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
