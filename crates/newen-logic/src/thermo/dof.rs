//! Kinetic degrees of freedom of the whole system.
//!
//! Theory: `spec/NewenTheory/Dynamics/DegreesOfFreedom.lean`; model, bridge and Aeneas
//! refinement: `spec/NewenSpec/Thermo/Dof.lean`.

/// Total kinetic degrees of freedom: `3 · n_atoms − solvent_constraints − solute_constraints −
/// ndfmin`.
///
/// `n_atoms` point particles in three dimensions have `3 · n_atoms` coordinates; every
/// independent holonomic constraint the run imposes removes one, and `ndfmin` (NDFMIN, BOUNDCOND
/// block) removes overall motion, e.g. 3 for the centre-of-mass translation. This is the
/// `N_df` of equipartition, `⟨K⟩ = ½ N_df k_B T`, so it sets the instantaneous temperature and
/// what the thermostat couples to.
///
/// The arithmetic is signed: more constraints than coordinates is not a physical system, and it
/// gives a negative count instead of a wrapped one. The Lean proof covers every input with
/// `3 · n_atoms + solvent_constraints + solute_constraints + |ndfmin| ≤ i64::MAX`. Beyond that an
/// operation can overflow (a panic with overflow checks on) and a count `≥ 2^63` would wrap in
/// its `as i64`; counts taken from Rust collections reach neither.
///
/// spec: thermo.dof.total_dof
pub fn total_dof(
    n_atoms: u64,
    solvent_constraints: u64,
    solute_constraints: u64,
    ndfmin: i64,
) -> i64 {
    let coordinates = 3 * (n_atoms as i64);
    coordinates - (solvent_constraints as i64) - (solute_constraints as i64) - ndfmin
}

#[cfg(test)]
mod tests {
    use super::total_dof;

    #[test]
    fn unconstrained_atoms() {
        assert_eq!(total_dof(5, 0, 0, 0), 15);
    }

    #[test]
    fn constraints_and_ndfmin_are_subtracted() {
        // two bonded solute atoms and one rigid water: 15 − 3 − 1 − 3
        assert_eq!(total_dof(5, 3, 1, 3), 8);
    }

    #[test]
    fn more_constraints_than_coordinates_is_negative() {
        assert_eq!(total_dof(1, 2, 2, 0), -1);
        assert_eq!(total_dof(0, 0, 0, 3), -3);
    }
}
