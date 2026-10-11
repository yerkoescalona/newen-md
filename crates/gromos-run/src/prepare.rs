//! Everything that happens *before* the algorithm sequence — lifted from `md.rs` so the
//! Python binding does it the same way (it did not: NSM came from the parameter file instead
//! of the coordinate file, validation was skipped, the topology's PHYSICALCONSTANTS were
//! ignored — 3.8 §2).
//!
//! GROMOS order (`in_topology` → `in_configuration` → `check_configuration`), kept:
//! perturbation topology → truncated-octahedron transform → NSM from the coordinates and
//! `solvate` → topology validation → atom-count check → coordinate validation →
//! configuration (initial velocities per `NTIVEL`) → configuration validation.

use std::path::PathBuf;

use gromos_core::{
    configuration::{Box as SimBox, Configuration},
    gather::gather_chargegroups,
    math::{
        truncoct_triclinic, truncoct_triclinic_box, Mat3, Periodicity, Rectangular, Triclinic, Vec3,
    },
    random::generate_velocities,
    units::PhysicalConstants,
    validation::{
        validate_configuration, validate_coordinates, validate_topology, ValidationReport,
    },
    Topology,
};
use gromos_io::{coordinate::CoordinateData, imd::ImdParameters, ptp::read_pttopo};

use crate::{RunError, RunInputs};

/// Positions, velocities and box as read from a coordinate file (or built in memory).
#[derive(Debug, Clone)]
pub struct Coordinates {
    pub positions: Vec<Vec3>,
    /// May be empty or of the wrong length: then velocities start at zero (or are generated
    /// when `NTIVEL=1`), exactly as the `md` binary does.
    pub velocities: Vec<Vec3>,
    /// Zero for a vacuum system.
    pub box_dims: Vec3,
}

impl From<CoordinateData> for Coordinates {
    fn from(c: CoordinateData) -> Self {
        Self {
            positions: c.positions,
            velocities: c.velocities,
            box_dims: c.box_dims,
        }
    }
}

/// Something the caller may want to report; never an error.
pub enum PrepareNote {
    /// A perturbation topology was given but `NTG=0`: ignored, as gromosXX does.
    PerturbationIgnored { path: PathBuf },
    /// A perturbation topology was merged (`NTG != 0`).
    PerturbationLoaded {
        atoms: usize,
        bonds: usize,
        angles: usize,
        impropers: usize,
        dihedrals: usize,
        soft_bonds: usize,
        soft_angles: usize,
        soft_impropers: usize,
    },
    /// The number of solvent molecules in the coordinate file differs from `NSM`; the
    /// coordinate file wins.
    NsmAdjusted { imd: usize, coordinates: usize },
    /// The last energy group (`FORCE` NRE) ends before the last atom: gromosXX appends one more
    /// group ending at the last atom, with a warning ("added an additional energy group").
    EnergyGroupAdded {
        /// last atom of the appended group, 1-based
        last_atom: usize,
    },
    /// Non-fatal validation findings (errors or warnings) for `stage`.
    Validation {
        stage: &'static str,
        report: ValidationReport,
    },
}

/// The system, ready for [`crate::build_sequence_from_imd`].
pub struct Prepared {
    pub topology: Topology,
    pub configuration: Configuration,
    /// From the topology's PHYSICALCONSTANTS block (gromosXX honours it; so does `Forcefield`).
    pub physical_constants: PhysicalConstants,
    /// Box edge lengths as read (before any truncated-octahedron transform); zero for vacuum.
    pub box_dims: Vec3,
    /// The lower-triangular triclinic box when `NTB=-1`, else `None`.
    pub truncoct_box: Option<Mat3>,
    pub notes: Vec<PrepareNote>,
}

/// Prepare topology and configuration for a run.
///
/// `topology` is a built topology, solvated or not (the Python object path may already have
/// called `solvate`; the binary never has). `inputs.pttopo` is merged only when `NTG != 0`.
pub fn prepare_system(
    imd: &ImdParameters,
    mut topology: Topology,
    physical_constants: PhysicalConstants,
    coords: Coordinates,
    inputs: &RunInputs,
) -> Result<Prepared, RunError> {
    let mut notes = Vec::new();

    // 1. Perturbation topology (FEP), only when NTG != 0 — gromosXX ignores it otherwise.
    if let Some(path) = &inputs.pttopo {
        if imd.ntg != 0 {
            let pt = read_pttopo(path).map_err(|e| RunError::Io {
                what: format!("perturbation topology '{}'", path.display()),
                source: e,
            })?;
            topology.apply_perturbation(pt.perturbed_solute, pt.is_perturbed);
            let p = &topology.perturbed_solute;
            notes.push(PrepareNote::PerturbationLoaded {
                atoms: p.atoms.len(),
                bonds: p.bonds.len(),
                angles: p.angles.len(),
                impropers: p.improper_dihedrals.len(),
                dihedrals: p.proper_dihedrals.len(),
                soft_bonds: p.soft_bonds.len(),
                soft_angles: p.soft_angles.len(),
                soft_impropers: p.soft_impropers.len(),
            });
        } else {
            notes.push(PrepareNote::PerturbationIgnored { path: path.clone() });
        }
    }
    if imd.ntg == 0 && !topology.perturbed_solute.atoms.is_empty() {
        // Only reachable from the object path (`Topology.apply_perturbation` in Python): the
        // regular solute already lost its perturbed terms, so running with NTG=0 would not be
        // what the binary does with the same files (it ignores the .ptp entirely).
        return Err(RunError::Inconsistent(
            "the topology carries a perturbation but NTG=0; set NTG != 0 or use an \
             unperturbed topology"
                .to_string(),
        ));
    }

    // gromosXX `Topology::update_for_lambda`: a perturbed atom's mass is λ-mixed,
    // m(λ) = (1 − λ)·m_A + λ·m_B (the mass λ has no NLAM exponent). It changes the kinetic
    // energy and the integration, so it is part of the prepared topology — whichever way the
    // perturbation arrived (the `.ptp` file above, or `Topology.apply_perturbation` from Python).
    if imd.ntg != 0 && !topology.perturbed_solute.atoms.is_empty() {
        let lambda = imd.rlam;
        let mixed: Vec<(usize, f64)> = topology
            .perturbed_solute
            .atoms
            .iter()
            .map(|a| (a.seq, (1.0 - lambda) * a.a_mass + lambda * a.b_mass))
            .collect();
        for (seq, mass) in mixed {
            if seq < topology.mass.len() && mass > 0.0 {
                topology.mass[seq] = mass;
                topology.inverse_mass[seq] = 1.0 / mass;
            }
        }
    }
    // 2. NTB=-1: truncated octahedron. GROMOS converts the legacy "cube edge length L" BOX
    //    block into the lower-triangular triclinic box vectors and rotates positions and
    //    velocities into that frame on read (io/configuration/in_configuration.cc).
    let Coordinates {
        mut positions,
        mut velocities,
        box_dims,
    } = coords;
    let truncoct_box = if imd.ntb == -1 {
        let cubic = Mat3::from_diagonal(box_dims);
        let triclinic_box = truncoct_triclinic_box(cubic, true);
        truncoct_triclinic(&mut positions, true);
        truncoct_triclinic(&mut velocities, true);
        Some(triclinic_box)
    } else {
        None
    };

    // 3. NSM from the coordinate file (the binary's rule; the parameter file's NSM is only a
    //    hint), then solvate — unless the caller already did.
    let unsolvated = topology.num_atoms() == topology.num_solute_atoms();
    if unsolvated {
        let actual_nsm = if imd.nsm > 0 && !topology.solvent_atom_template.is_empty() {
            let atoms_per_solvent = topology.solvent_atom_template.len();
            let n_solute = topology.num_solute_atoms();
            let remaining = positions.len().saturating_sub(n_solute);
            if remaining % atoms_per_solvent != 0 {
                return Err(RunError::SolventCount {
                    coordinates: positions.len(),
                    solute: n_solute,
                    atoms_per_solvent,
                });
            }
            let from_coords = remaining / atoms_per_solvent;
            if from_coords != imd.nsm {
                notes.push(PrepareNote::NsmAdjusted {
                    imd: imd.nsm,
                    coordinates: from_coords,
                });
            }
            from_coords
        } else {
            imd.nsm
        };
        topology.solvate(actual_nsm);
    }

    // 3a. gromosXX in_topology.cc:278-290 (with solute molecules only): the last temperature
    //     group of the solute ends with the last solute atom.
    check_solute_temperature_groups(&topology, imd.npm)?;

    // 3b. Energy groups over every atom, solvent included (gromosXX assigns them after
    //     `solvate`, in_topology.cc:433-459).
    if let Some(note) = assign_energy_groups(&mut topology, &imd.nre)? {
        notes.push(note);
    }

    // 3c. Hand the bonds that a constraint algorithm owns to the topology, so the bonded force
    //     loop leaves them alone — gromosXX does this when it creates the constraints
    //     (`create_constraints.cc`): a constrained bond has no potential energy and no force.
    {
        use gromos_integrators::constraints::ShakeBuffers;
        let sel = crate::ConstraintSelection::from_imd(
            imd,
            topology.num_atoms() > topology.num_solute_atoms(),
        );
        topology.constrained_bonds.clear();
        if sel.solute_constrained() {
            let ntc = crate::ConstraintSelection::ntc_mode_of(imd.ntc);
            {
                for (i, j, _) in &ShakeBuffers::new(&topology, ntc, false).solute_constraints {
                    topology.constrained_bonds.insert((*i, *j));
                    topology.constrained_bonds.insert((*j, *i));
                }
            }
        }
    }

    // 3d. Make every charge group whole around its first atom, as gromosXX does on reading a
    //     configuration. A `.cnf` written with the molecules put back into the box has charge
    //     groups split across the boundary; the pairlist works on their centres of geometry.
    if imd.ntb != 0 && box_dims.x > 0.0 {
        let pbc = match truncoct_box {
            Some(b) => Periodicity::Triclinic(Triclinic::new(b)),
            None => Periodicity::Rectangular(Rectangular::new(box_dims)),
        };
        gather_chargegroups(&mut positions, &topology.chargegroups, &pbc);
    }

    // 4. Topology validation (fatal stops; errors and warnings are reported).
    validate(&mut notes, "topology", validate_topology(&topology))?;

    // 5. check_configuration: atom counts must match.
    let n_atoms = topology.num_atoms();
    if positions.len() != n_atoms {
        return Err(RunError::AtomCountMismatch {
            topology: n_atoms,
            coordinates: positions.len(),
        });
    }

    // 6. Coordinate validation.
    validate(
        &mut notes,
        "coordinates",
        validate_coordinates(&positions, Some(box_dims)),
    )?;

    // 7. Configuration (double-buffered state). Velocities: generated (NTIVEL=1), read, or zero.
    let n_baths = imd.num_temp_baths.max(1);
    let mut conf = Configuration::new(n_atoms, n_baths, topology.energy_groups.len().max(1));
    conf.current_mut().pos = positions;
    conf.current_mut().vel = if imd.ntivel == 1 {
        generate_velocities(imd.tempi, imd.ig as u32, &topology.mass)
    } else if velocities.len() == n_atoms {
        velocities
    } else {
        vec![Vec3::ZERO; n_atoms]
    };
    conf.current_mut().box_config = match truncoct_box {
        Some(triclinic_box) => SimBox::truncated_octahedral(triclinic_box),
        None => SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z),
    };
    conf.copy_current_to_old();

    // 8. Configuration validation (topology + coordinates together).
    validate(
        &mut notes,
        "configuration",
        validate_configuration(&topology, &conf),
    )?;

    Ok(Prepared {
        topology,
        configuration: conf,
        physical_constants,
        box_dims,
        truncoct_box,
        notes,
    })
}

/// Partition the atoms into energy groups from the `FORCE` block's NRE (1-based last atom of
/// each group), as gromosXX's `In_Topology` does (in_topology.cc:433-459): none → one group
/// over everything; a last NRE beyond the last atom is an error; one before the last atom gets
/// a trailing group appended (and a note). NRE must be positive and strictly increasing
/// (in_parameter.cc:978-993).
fn assign_energy_groups(
    topology: &mut Topology,
    nre: &[usize],
) -> Result<Option<PrepareNote>, RunError> {
    let n_atoms = topology.num_atoms();
    let mut last_atoms = Vec::with_capacity(nre.len() + 1);
    let mut previous = 0usize;
    for &e in nre {
        if e == 0 || e <= previous {
            return Err(RunError::Recipe(format!(
                "FORCE block: energy group last atoms {nre:?} must be positive and strictly                  increasing"
            )));
        }
        previous = e;
        last_atoms.push(e - 1);
    }
    let mut note = None;
    match last_atoms.last().copied() {
        None => last_atoms.push(n_atoms.saturating_sub(1)),
        Some(last) if last + 1 > n_atoms => {
            return Err(RunError::Recipe(format!(
                "FORCE block: the last energy group ends at atom {} but the system has {n_atoms}                  atoms",
                last + 1
            )));
        },
        Some(last) if last + 1 < n_atoms => {
            last_atoms.push(n_atoms - 1);
            note = Some(PrepareNote::EnergyGroupAdded { last_atom: n_atoms });
        },
        Some(_) => {},
    }
    topology.energy_groups.clear();
    topology.atom_to_energy_group.clear();
    topology.atom_to_energy_group.reserve(n_atoms);
    let mut atom = 0usize;
    for (g, &last) in last_atoms.iter().enumerate() {
        let mut members = Vec::with_capacity((last + 1).saturating_sub(atom));
        while atom <= last && atom < n_atoms {
            members.push(atom);
            topology.atom_to_energy_group.push(g);
            atom += 1;
        }
        topology.energy_groups.push(members);
    }
    Ok(note)
}

/// gromosXX `In_Topology::read` (in_topology.cc:278-290): a TEMPERATUREGROUPS block whose last
/// group stops before the last solute atom is an error for MD (gromos++ reads such a topology).
fn check_solute_temperature_groups(topo: &Topology, npm: usize) -> Result<(), RunError> {
    let n_solute = topo.num_solute_atoms();
    if npm == 0 || n_solute == 0 {
        return Ok(());
    }
    let end = topo
        .temperature_groups
        .iter()
        .filter_map(|g| g.last().copied())
        .filter(|&last| last < n_solute)
        .max();
    match end {
        Some(last) if last + 1 != n_solute => Err(RunError::Inconsistent(format!(
            "Error in TEMPERATUREGROUPS block: last temperature group has to end with last \
             solute atom (it ends at atom {}, the solute at atom {n_solute})",
            last + 1
        ))),
        _ => Ok(()),
    }
}

fn validate(
    notes: &mut Vec<PrepareNote>,
    stage: &'static str,
    report: ValidationReport,
) -> Result<(), RunError> {
    if report.has_fatal() {
        return Err(RunError::Validation { stage, report });
    }
    if report.has_errors() || !report.warnings.is_empty() {
        notes.push(PrepareNote::Validation { stage, report });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_io::coordinate::read_coordinates;
    use gromos_io::topology::{build_topology, read_topology_file};
    use std::path::PathBuf;

    fn shared() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../gromos-md/tests/gromosXX_references/shared")
    }

    fn aladip(imd: &ImdParameters, pttopo: bool) -> Prepared {
        let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
        let constants = parsed.physical_constants;
        let topo = build_topology(parsed);
        let coords: Coordinates = read_coordinates(shared().join("aladip.conf"))
            .unwrap()
            .into();
        let inputs = RunInputs {
            pttopo: pttopo.then(|| shared().join("aladip.pttopo")),
            ..Default::default()
        };
        prepare_system(imd, topo, constants, coords, &inputs).unwrap()
    }

    #[test]
    fn energy_groups_cover_every_atom_as_gromosxx_assigns_them() {
        // NRE stops at the last solute atom: the solvent gets a group of its own
        let imd = ImdParameters {
            nsm: 20,
            negr: 1,
            nre: vec![12],
            ..Default::default()
        };
        let p = aladip(&imd, false);
        let t = &p.topology;
        assert_eq!(t.energy_groups.len(), 2);
        assert_eq!(t.atom_to_energy_group.len(), 72);
        assert_eq!((t.energy_group(11), t.energy_group(12)), (0, 1));
        assert!(p
            .notes
            .iter()
            .any(|n| matches!(n, PrepareNote::EnergyGroupAdded { last_atom: 72 })));
        assert_eq!(p.configuration.current().energies.num_energy_groups(), 2);
        assert_eq!(
            p.configuration
                .current()
                .perturbed_energy_derivatives
                .lj_energy
                .len(),
            2
        );

        // no NRE: one group over everything
        let imd = ImdParameters {
            nsm: 20,
            ..Default::default()
        };
        let p = aladip(&imd, false);
        assert_eq!(p.topology.energy_groups, vec![(0..72).collect::<Vec<_>>()]);
    }

    #[test]
    fn energy_groups_beyond_the_last_atom_are_an_error() {
        let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
        let constants = parsed.physical_constants;
        let topo = build_topology(parsed);
        let coords: Coordinates = read_coordinates(shared().join("aladip.conf"))
            .unwrap()
            .into();
        let imd = ImdParameters {
            nsm: 20,
            nre: vec![6, 73],
            ..Default::default()
        };
        let err = prepare_system(&imd, topo, constants, coords, &RunInputs::default());
        assert!(matches!(err, Err(RunError::Recipe(m)) if m.contains("ends at atom 73")));
    }

    #[test]
    fn solvates_to_the_configuration_size() {
        let imd = ImdParameters {
            nsm: 20,
            ..Default::default()
        };
        let p = aladip(&imd, false);
        assert_eq!(p.topology.num_atoms(), 72);
        assert_eq!(p.topology.num_solute_atoms(), 12);
        assert_eq!(p.configuration.current().pos.len(), 72);
    }

    #[test]
    fn perturbation_mixes_the_masses_only_when_ntg_is_set() {
        let mut imd = ImdParameters {
            nsm: 20,
            ..Default::default()
        };
        let plain = aladip(&imd, false);
        let ignored = aladip(&imd, true); // NTG = 0: gromosXX ignores the perturbation topology
        assert_eq!(ignored.topology.mass, plain.topology.mass);
        assert!(ignored.topology.perturbed_solute.is_empty());
        imd.ntg = 1;
        imd.rlam = 0.5;
        let mixed = aladip(&imd, true);
        assert!(!mixed.topology.perturbed_solute.is_empty());
        assert_ne!(
            mixed.topology.mass, plain.topology.mass,
            "λ-mixed masses expected"
        );
    }
}
