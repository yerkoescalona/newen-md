//! Kinetic degrees of freedom — one formula for the thermostat and the reported temperature.
//!
//! gromosXX (`Multibath::calculate_degrees_of_freedom`): `3N − N_constraints − NDFMIN`, where
//! `N_constraints` counts the solvent constraints of every solvent molecule (when the solvent
//! is constrained) and the solute distance constraints (when the solute is constrained).
//!
//! Before this crate, three copies of this computation existed and none counted the solute
//! constraints (`md.rs` carried a `TODO`); they also subtracted the solvent constraints whenever
//! *any* constraint algorithm was on, even one acting on the solute only. 3.9 A11.

use newen_core::Topology;
use newen_integrators::constraints::{NtcMode, ShakeBuffers};

use newen_integrators::algorithms::{in_bath, BathRange};

use crate::ConstraintSelection;

/// The degrees of freedom of one temperature bath, as gromosXX's `bath_struct` keeps them.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BathDof {
    /// All of the bath's kinetic degrees of freedom (`dof`): what the thermostat couples to.
    pub dof: f64,
    /// Centre-of-mass translation (`com_dof`).
    pub com: f64,
    /// Internal and rotational (`ir_dof`).
    pub ir: f64,
}

/// Degrees of freedom per temperature bath, as gromosXX's
/// `Multibath::calculate_degrees_of_freedom` counts them:
///
/// * every DOFSET range contributes `3 × (number of temperature groups starting in it)` to the
///   COM bath and the remaining `3 × (atoms in the range)` minus that to the IR bath;
/// * constraints are subtracted from the bath that owns their atoms;
/// * roto-translational constraints (`rottrans`) take 3 from the COM and 3 from the IR degrees
///   of freedom of atom 1's baths (`topology.cc:1370-1402`);
/// * NDFMIN is spread over the baths in proportion to their degrees of freedom, from `dof`
///   only (`multibath.cc:139-149`).
///
/// With one bath over the whole system `dof` is exactly [`total_dof`] less [`rottrans_dof`].
pub fn bath_dofs(
    topo: &Topology,
    sel: &ConstraintSelection,
    ntc: NtcMode,
    ndfmin: i32,
    n_baths: usize,
    ranges: &[BathRange],
    rottrans: bool,
) -> Vec<BathDof> {
    let n_atoms = topo.num_atoms();
    let mut com = vec![0.0; n_baths];
    let mut ir = vec![0.0; n_baths];
    let mut last: i64 = -1;
    for r in ranges {
        let group_starts = topo
            .temperature_groups
            .iter()
            .filter_map(|g| g.first().copied())
            .filter(|&first| first as i64 > last && first <= r.last_atom)
            .count();
        let atoms_in_range = r.last_atom as i64 - last;
        com[r.com_bath] += 3.0 * group_starts as f64;
        ir[r.ir_bath] += 3.0 * atoms_in_range as f64 - 3.0 * group_starts as f64;
        last = r.last_atom as i64;
    }
    // constraints go to the IR bath of their atoms, half to each end
    // (gromosXX Topology::calculate_constraint_dof, topology.cc:1076-1110)
    let ir_bath = |atom: usize| in_bath(ranges, atom).1;
    if sel.solute_constrained() {
        for (i, j) in ShakeBuffers::new(topo, ntc, false).solute_constraint_pairs() {
            ir[ir_bath(i)] -= 0.5;
            ir[ir_bath(j)] -= 0.5;
        }
    }
    if sel.solvent_constrained() {
        let per_molecule = topo.solvent_constraint_template.len();
        let n_solute = topo.num_solute_atoms();
        let atoms_per_molecule = topo.solvent_atom_template.len().max(1);
        for m in 0..topo.num_solvent_molecules() {
            let first = n_solute + m * atoms_per_molecule;
            if first < n_atoms {
                ir[ir_bath(first)] -= per_molecule as f64;
            }
        }
    }
    if rottrans && topo.num_solute_atoms() > 0 {
        let (com_bath, ir_bath) = in_bath(ranges, 0);
        com[com_bath] -= 3.0;
        ir[ir_bath] -= 3.0;
    }
    let mut baths: Vec<BathDof> = com
        .iter()
        .zip(&ir)
        .map(|(&com, &ir)| BathDof {
            dof: com + ir,
            com,
            ir,
        })
        .collect();
    // NDFMIN in proportion to each bath's degrees of freedom (gromosXX's "Martin Stroet edit")
    let total: f64 = baths.iter().map(|b| b.dof).sum();
    if ndfmin != 0 && total > 0.0 {
        for b in baths.iter_mut() {
            b.dof -= ndfmin as f64 * (b.dof / total);
        }
    }
    baths
}

/// [`bath_dofs`]' total degrees of freedom per bath.
pub fn bath_dof(
    topo: &Topology,
    sel: &ConstraintSelection,
    ntc: NtcMode,
    ndfmin: i32,
    n_baths: usize,
    ranges: &[BathRange],
    rottrans: bool,
) -> Vec<f64> {
    bath_dofs(topo, sel, ntc, ndfmin, n_baths, ranges, rottrans)
        .iter()
        .map(|b| b.dof)
        .collect()
}

/// The degrees of freedom roto-translational constraints remove from a single bath: 3 COM and
/// 3 internal-rotational when there is a solute (`topology.cc:1370-1402`).
pub fn rottrans_dof(topo: &Topology, rottrans: bool) -> f64 {
    if rottrans && topo.num_solute_atoms() > 0 {
        6.0
    } else {
        0.0
    }
}

/// Total kinetic degrees of freedom of the system.
pub fn total_dof(topo: &Topology, sel: &ConstraintSelection, ntc: NtcMode, ndfmin: i32) -> f64 {
    let n_atoms = topo.num_atoms();
    let solvent_constraints = if sel.solvent_constrained() {
        topo.num_solvent_molecules() * topo.solvent_constraint_template.len()
    } else {
        0
    };
    let solute_constraints = if sel.solute_constrained() {
        ShakeBuffers::new(topo, ntc, false)
            .solute_constraint_pairs()
            .count()
    } else {
        0
    };
    (3 * n_atoms - solvent_constraints - solute_constraints) as f64 - ndfmin as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::topology::{
        Bond, BondParameters, MolTypeAtom, SolventAtomTemplate, SolventConstraintTemplate,
    };

    /// Two bonded solute atoms plus one rigid SPC water: 3·5 = 15 degrees of freedom.
    fn dipeptide_and_a_water() -> Topology {
        let mut topo = Topology::new();
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0; 2];
        topo.moltypes[0].atoms = vec![
            MolTypeAtom {
                name: "C1".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
            MolTypeAtom {
                name: "C2".into(),
                residue_nr: 0,
                residue_name: String::new(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            },
        ];
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 0.0,
            r0: 0.15,
        });
        topo.solvent_atom_template = vec![
            SolventAtomTemplate {
                iac: 0,
                name: "OW".into(),
                mass: 15.9994,
                charge: -0.82,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW1".into(),
                mass: 1.008,
                charge: 0.41,
            },
            SolventAtomTemplate {
                iac: 1,
                name: "HW2".into(),
                mass: 1.008,
                charge: 0.41,
            },
        ];
        topo.solvent_constraint_template = vec![
            SolventConstraintTemplate {
                i: 0,
                j: 1,
                length: 0.1,
            },
            SolventConstraintTemplate {
                i: 0,
                j: 2,
                length: 0.1,
            },
            SolventConstraintTemplate {
                i: 1,
                j: 2,
                length: 0.1633,
            },
        ];
        topo.solvate(1);
        topo
    }

    #[test]
    fn everything_constrained() {
        let topo = dipeptide_and_a_water();
        assert_eq!(topo.num_atoms(), 5);
        let sel = ConstraintSelection::from_codes(3, 1, 1, true);
        // 15 − 3 (water) − 1 (solute bond) − 3 (NDFMIN)
        assert_eq!(total_dof(&topo, &sel, NtcMode::AllBonds, 3), 8.0);
    }

    #[test]
    fn solute_constraints_are_not_subtracted_when_only_the_solvent_is_constrained() {
        let topo = dipeptide_and_a_water();
        let sel = ConstraintSelection::from_codes(1, 1, 1, true);
        assert_eq!(total_dof(&topo, &sel, NtcMode::SolventOnly, 0), 12.0);
    }

    #[test]
    fn nothing_constrained() {
        let topo = dipeptide_and_a_water();
        let sel = ConstraintSelection::from_codes(1, 1, 0, true);
        assert_eq!(total_dof(&topo, &sel, NtcMode::SolventOnly, 0), 15.0);
    }
}
