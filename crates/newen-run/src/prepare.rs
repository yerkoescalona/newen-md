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

use newen_core::{
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
    IndividualLambdas, NetCharge, PrecalcLambdas, Topology,
};
use newen_io::{coordinate::CoordinateData, imd::ImdParameters, ptp::read_pttopo};

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
    /// λ of the configuration's `PERTDATA` block, if it has one; a run with NRDGL=1 starts
    /// from it.
    pub pertdata: Option<f64>,
    /// The data lines of the configuration's `AEDSSEARCH` block, if it has one; an A-EDS search
    /// with NTIAEDSS = 0 continues from it.
    pub aedssearch: Option<Vec<String>>,
    /// The data lines of the configuration's `NHCVARIABLES` block, if it has one; a
    /// Nosé–Hoover run with NTINHT = 0 starts from it.
    pub nhcvariables: Option<Vec<String>>,
}

impl From<CoordinateData> for Coordinates {
    fn from(c: CoordinateData) -> Self {
        Self {
            positions: c.positions,
            velocities: c.velocities,
            box_dims: c.box_dims,
            pertdata: c.pertdata,
            aedssearch: c.aedssearch,
            nhcvariables: c.nhcvariables,
        }
    }
}

/// Something the caller may want to report; never an error.
pub enum PrepareNote {
    /// A perturbation topology was given but `NTG=0`: ignored, as gromosXX does.
    PerturbationIgnored { path: PathBuf },
    /// The EDS end states of a perturbation topology's `MPERTATOM` block were merged.
    EdsLoaded {
        /// number of EDS atoms
        atoms: usize,
        /// number of end states
        states: usize,
    },
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
    /// Net charge of the solute in the two states of the perturbation (`NTG != 0`).
    PerturbationNetCharge { charge: NetCharge },
    /// The perturbation changes the solute's net charge under a reaction field (NLRELE −1, 1
    /// or 4): the ΔG of a charge change carries a finite-size artefact and needs a correction
    /// (e.g. `dGslv_pbsolv`). A warning, not a refusal; gromosXX says nothing.
    NetChargeChanged {
        /// ΔQ = Q_B − Q_A (e)
        delta: f64,
        nlrele: i32,
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
    /// A perturbed atom pair whose energy groups have different individual λ for Lennard-Jones
    /// and reaction field: gromosXX runs it but reports a critical message ("Perturbed atom
    /// pairs with different individual lambdas for lj and crf not implemented",
    /// `perturbed_nonbonded_pair.cc:309-312`) — its force uses the LJ λ for both.
    AtomPairLambdasDiffer {
        /// the pair's atoms, 1-based
        atoms: (usize, usize),
    },
    /// Perturbed SHAKE together with PRECALCLAM: the constraints add to the `AB_bond`
    /// derivatives at every λ_k but not to the energies, so BAR or reweighting in extended TI
    /// is inappropriate (gromosXX warns, `in_perturbation.cc:327-334`).
    PerturbedShakePrecalc {
        /// number of perturbed distance constraints
        constraints: usize,
    },
    /// PRECALCLAM with a twin-range pairlist whose update (NSNB) and free-energy output (NTWG)
    /// steps interleave: the λ_k rows written after a pairlist step that was not an output
    /// step miss the long-range terms (gromosXX behaviour, reproduced; porting.md), so BAR or
    /// MBAR from them is biased. A warning, not a refusal; gromosXX says nothing
    /// ([`precalclam_longrange_gaps`], [`precalclam_longrange_warning`]).
    PrecalclamLongRangeGaps { nsnb: usize, ntwg: usize },
    /// The configuration has an `AEDSSEARCH` block but the A-EDS search starts afresh
    /// (NTIAEDSS = 1): ignored, as gromosXX does with a warning.
    AedsSearchIgnored,
    /// The configuration has an `NHCVARIABLES` block but the Nosé–Hoover variables start at
    /// zero (NTINHT = 1): ignored, as gromosXX does with a warning ("Nose-Hoover-Chains
    /// variables found but not read.", `in_configuration.cc:1166-1169`).
    NoseHooverVariablesIgnored,
    /// A Nosé–Hoover run asks for its variables from the configuration (NTINHT = 0) but the
    /// configuration has no `NHCVARIABLES` block: they start at zero. A warning, not a refusal
    /// (gromosXX refuses it for a chain and reads nothing for plain Nosé–Hoover); callers report
    /// it under the intended `nh_restart` only ([`nose_hoover_missing_warning`]).
    NoseHooverVariablesMissing,
    /// The run starts at a λ outside [0, 1] (a `PERTDATA` λ; RLAM is checked to [0, 1]): every
    /// perturbed term and mass is extrapolated. A warning, not a refusal; gromosXX says nothing
    /// (porting.md, slow growth). Slow growth that leaves [0, 1] during the run is reported by
    /// [`crate::advance_lambda`].
    LambdaOutOfRange { lambda: f64 },
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

    // 1. Perturbation topology: its PERT* blocks when NTG != 0 (FEP), its MPERTATOM block with
    //    EDS or AEDS on (gromosXX ignores the file with neither, and requires it for EDS).
    let aeds = imd.aeds.as_ref().filter(|a| a.enabled);
    let eds = imd
        .eds
        .as_ref()
        .filter(|e| e.enabled)
        .map(|e| e.num_states)
        .or(aeds.map(|a| a.num_states));
    if eds.is_some() && inputs.pttopo.is_none() && topology.eds_solute.is_empty() {
        return Err(RunError::MissingInput {
            flag: "pttopo",
            reason: "EDS = 1 or AEDS = 1 (the end states come from its MPERTATOM block)"
                .to_string(),
        });
    }
    if let Some(path) = &inputs.pttopo {
        if imd.ntg != 0 || eds.is_some() {
            let pt = read_pttopo(path).map_err(|e| RunError::Io {
                what: format!("perturbation topology '{}'", path.display()),
                source: e,
            })?;
            check_eds_perturbation(eds, &pt, imd.ntg != 0)?;
            if imd.ntg != 0 {
                topology
                    .apply_perturbation(pt.perturbed_solute.clone(), pt.is_perturbed.clone())
                    .map_err(RunError::Inconsistent)?;
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
            }
            // (the object path may have merged it already: `Topology.apply_perturbation`)
            let merged = !topology.eds_solute.is_empty();
            if let Some(m) = pt.eds.as_ref().filter(|_| eds.is_some() && !merged) {
                topology
                    .apply_eds_perturbation(m.num_states, m.state_names.clone(), m.atoms.clone())
                    .map_err(RunError::Inconsistent)?;
                notes.push(PrepareNote::EdsLoaded {
                    atoms: m.atoms.len(),
                    states: m.num_states,
                });
            }
        } else {
            notes.push(PrepareNote::PerturbationIgnored { path: path.clone() });
        }
    }
    match (eds, topology.eds_solute.is_empty()) {
        (Some(n), false) if n != topology.eds_solute.num_states => {
            return Err(RunError::Inconsistent(format!(
                "Number of perturbed states given in perturbation topology ({}) and input file \
                 ({}) do not match",
                topology.eds_solute.num_states, n
            )))
        },
        (None, false) => {
            return Err(RunError::Inconsistent(
                "the topology carries EDS end states (MPERTATOM) but the input has no EDS \
                 block with EDS = 1"
                    .to_string(),
            ))
        },
        _ => {},
    }
    check_ti_eds(&topology)?;
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
    if imd.ntg != 0 && !topology.perturbed_solute.atoms.is_empty() {
        let charge = topology.perturbation_net_charge();
        notes.push(PrepareNote::PerturbationNetCharge { charge });
        if charge.delta() != 0.0 && matches!(imd.nlrele, -1 | 1 | 4) {
            notes.push(PrepareNote::NetChargeChanged {
                delta: charge.delta(),
                nlrele: imd.nlrele,
            });
        }
    }

    let Coordinates {
        mut positions,
        mut velocities,
        box_dims,
        pertdata,
        aedssearch,
        nhcvariables,
    } = coords;
    // gromosXX `read_nose_hoover_chains` (`in_configuration.cc:1148-1178`): with NTINHT = 0 the
    // variables come from the configuration's NHCVARIABLES block. Kept whole here; the
    // thermostat restores them per `nh_restart` (`restored_zeta`).
    let nose_hoover_run = imd.temp_bath.first().is_some_and(|b| b.algorithm >= 1);
    let nose_hoover = match nhcvariables.as_ref().filter(|_| nose_hoover_run) {
        Some(lines) if imd.ntinht == 0 => Some(
            newen_io::coordinate::parse_nhcvariables(lines).map_err(|e| {
                RunError::Inconsistent(format!("configuration NHCVARIABLES block: {e}"))
            })?,
        ),
        Some(_) => {
            notes.push(PrepareNote::NoseHooverVariablesIgnored);
            None
        },
        None => {
            if nose_hoover_run && imd.ntinht == 0 {
                notes.push(PrepareNote::NoseHooverVariablesMissing);
            }
            None
        },
    };
    // gromosXX `read_aedssearch` (`in_configuration.cc:1399-1455`): a search form (2-6) with
    // NTIAEDSS = 0 continues the configuration's AEDSSEARCH block (it overrides the input's
    // EMAX, EMIN and EIR); with NTIAEDSS = 1 the block is ignored with a warning.
    let aeds_search = match aeds.filter(|a| a.form != 1) {
        Some(a) if !a.ntiaedss => {
            let lines = aedssearch.as_ref().ok_or_else(|| {
                RunError::Inconsistent(
                    "AEDS NTIAEDSS = 0: no AEDSSEARCH block in configuration".to_string(),
                )
            })?;
            Some(
                newen_io::coordinate::parse_aedssearch(lines, a.num_states).map_err(|e| {
                    RunError::Inconsistent(format!("configuration AEDSSEARCH block: {e}"))
                })?,
            )
        },
        Some(_) => {
            if aedssearch.is_some() {
                notes.push(PrepareNote::AedsSearchIgnored);
            }
            None
        },
        None => None,
    };
    // gromosXX `in_topology.cc:367-370`: the run's λ and NLAM live on the topology. With NTG=1
    // and NRDGL=1, λ comes from the configuration's PERTDATA block instead, which must be
    // there (`in_configuration.cc:876-897`, `_read_pertdata` :2351-2371).
    topology.lambda = if imd.ntg != 0 && imd.nrdgl == 1 {
        pertdata.ok_or_else(|| {
            RunError::Inconsistent(
                "PERTURBATION NRDGL = 1 reads λ from the configuration, but it has no \
                 (readable) PERTDATA block"
                    .to_string(),
            )
        })?
    } else {
        imd.rlam
    };
    topology.lambda_exp = imd.nlam;
    if imd.ntg != 0 && !crate::build::lambda_in_range(topology.lambda) {
        notes.push(PrepareNote::LambdaOutOfRange {
            lambda: topology.lambda,
        });
    }
    // 2. NTB=-1: truncated octahedron. GROMOS converts the legacy "cube edge length L" BOX
    //    block into the lower-triangular triclinic box vectors and rotates positions and
    //    velocities into that frame on read (io/configuration/in_configuration.cc).
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

    // 3b'. The individual λ over those groups (LAMBDAS, identity where no row applies), then
    //      gromosXX `Topology::init` → `update_for_lambda`: evaluate them at RLAM and mix the
    //      perturbed masses with the mass λ of each atom's group. The masses change the kinetic
    //      energy and the integration, so they are part of the prepared topology — whichever
    //      way the perturbation arrived (`.ptp` above, or `Topology.apply_perturbation`).
    topology.individual_lambdas = individual_lambdas(imd, topology.energy_groups.len())?;
    topology.update_for_lambda();
    {
        use newen_core::lambdas::InteractionLambda::{Crf, Lj};
        let t = &topology.individual_lambdas;
        let group = |i: usize| topology.atom_to_energy_group.get(i).copied().unwrap_or(0);
        for ap in &topology.perturbed_solute.atom_pairs {
            let (n1, n2) = (group(ap.i), group(ap.j));
            if t.lambda(Lj, n1, n2) != t.lambda(Crf, n1, n2) {
                notes.push(PrepareNote::AtomPairLambdasDiffer {
                    atoms: (ap.i + 1, ap.j + 1),
                });
            }
        }
    }

    // 3c. Hand the bonds that a constraint algorithm owns to the topology, so the bonded force
    //     loop leaves them alone — gromosXX does this when it creates the constraints
    //     (`create_constraints.cc`): a constrained bond has no potential energy and no force.
    {
        use newen_integrators::constraints::ShakeBuffers;
        let sel = crate::ConstraintSelection::from_imd(
            imd,
            topology.num_atoms() > topology.num_solute_atoms(),
        );
        topology.constrained_bonds.clear();
        hand_perturbed_bonds_to_constraints(&mut topology, imd, &sel)?;
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

    // 3c'. PRECALCLAM (extended TI) on a perturbed run: the λ grid and the NTWG steps it is
    //      evaluated on travel with the topology, as λ and NLAM do.
    topology.precalc_lambdas = imd
        .precalclam
        .filter(|b| imd.ntg != 0 && b.nr_lambdas > 0)
        .map(|b| PrecalcLambdas {
            nr_lambdas: b.nr_lambdas,
            min_lambda: b.min_lambda,
            max_lambda: b.max_lambda,
            every: imd.ntwg,
        });
    let n_constraints = topology.perturbed_solute.distance_constraints.len();
    if topology.precalc_lambdas.is_some() && n_constraints > 0 {
        notes.push(PrepareNote::PerturbedShakePrecalc {
            constraints: n_constraints,
        });
    }
    if let Some((nsnb, ntwg)) = precalclam_longrange_gaps(imd) {
        notes.push(PrepareNote::PrecalclamLongRangeGaps { nsnb, ntwg });
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
    let n_groups = topology.energy_groups.len().max(1);
    let mut conf = Configuration::new(n_atoms, n_baths, n_groups);
    if let Some(pc) = &topology.precalc_lambdas {
        // gromosXX `Configuration::init` sizes both records for NRLAM (`configuration.cc:371-379`);
        // `copy_current_to_old` below sizes the old state alike
        let state = conf.current_mut();
        state.energies.precalc.resize(pc.nr_lambdas, n_groups);
        state
            .perturbed_energy_derivatives
            .precalc
            .resize(pc.nr_lambdas, n_groups);
    }
    if !topology.eds_solute.is_empty() {
        // gromosXX `Configuration::init` (`configuration.cc:403-419`): the end states in both
        // records and the end-state forces and virials
        let n = topology.eds_solute.num_states;
        let state = conf.current_mut();
        state.energies.resize_eds(n);
        state.perturbed_energy_derivatives.resize_eds(n);
        conf.eds.resize(n, n_atoms);
    }
    conf.aeds = aeds_search;
    conf.nose_hoover = nose_hoover;
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

/// Whether the `MPERTATOM` block of a perturbation topology can be merged (gromosXX
/// `in_perturbation.cc`): required with EDS on, which also forbids the PERT* blocks of a run without perturbation
/// ("block X not supported for eds!", `:1805-1820`) and needs NPTB = NUMSTATES (`:1650`);
/// refused without EDS (gromosXX reads it anyway and finds NPTB ≠ 0 states). With
/// perturbation the PERT* blocks run next to it (TI-EDS; [`check_ti_eds`]).
fn check_eds_perturbation(
    eds_states: Option<usize>,
    pt: &newen_io::ptp::PerturbedTopology,
    perturbed: bool,
) -> Result<(), RunError> {
    let err = |m: String| Err(RunError::Inconsistent(m));
    let Some(n) = eds_states else {
        if pt.eds.is_some() {
            return err(
                "the perturbation topology has an MPERTATOM block (EDS end states) but the \
                 input has no EDS block with EDS = 1"
                    .to_string(),
            );
        }
        return Ok(());
    };
    let Some(m) = &pt.eds else {
        return err(
            "EDS = 1 but the perturbation topology has no MPERTATOM block (the end states)"
                .to_string(),
        );
    };
    if m.num_states != n {
        return err(format!(
            "Number of perturbed states given in perturbation topology ({}) and input file ({}) \
             do not match",
            m.num_states, n
        ));
    }
    if !perturbed {
        if let Some(b) = pt.blocks.iter().find(|b| b.as_str() != "MPERTATOM") {
            return err(format!(
                "block {b} not supported for eds! (the perturbation topology of an EDS run \
                 without PERTURBATION NTG = 1 holds MPERTATOM only)"
            ));
        }
    }
    Ok(())
}

/// TI-EDS (EDS with perturbed atoms, gromosXX 869b8a24): what gromosXX accepts but does not
/// compute, refused on every path to a topology with EDS end states and perturbed terms. An
/// atom in both `PERTATOMPARAM` and `MPERTATOM` falls through gromosXX's pair classification
/// (code 6, no case: its pairs are dropped silently, `perturbed_nonbonded_innerloop.cc:54-65`);
/// a `PERTATOMPAIR` pair has no EDS branch (`perturbed_nonbonded_pair.cc` uses the topology
/// type and charge in every state).
fn check_ti_eds(topology: &Topology) -> Result<(), RunError> {
    let eds = &topology.eds_solute;
    if eds.is_empty() {
        return Ok(());
    }
    let eds_atom = |i: usize| eds.atoms.iter().any(|a| a.seq == i);
    let pert = &topology.perturbed_solute;
    if let Some(a) = pert.atoms.iter().find(|a| eds_atom(a.seq)) {
        return Err(RunError::Inconsistent(format!(
            "atom {} is in both PERTATOMPARAM and MPERTATOM: an atom is either perturbed or \
             has EDS end states (gromosXX drops its pairs silently)",
            a.seq + 1
        )));
    }
    if let Some(p) = pert
        .atom_pairs
        .iter()
        .find(|p| eds_atom(p.i) || eds_atom(p.j))
    {
        return Err(RunError::Inconsistent(format!(
            "PERTATOMPAIR {}-{} involves an EDS atom: not supported (gromosXX computes it with \
             the topology type and charge in every end state)",
            p.i + 1,
            p.j + 1
        )));
    }
    Ok(())
}

/// Partition the atoms into energy groups from the `FORCE` block's NRE (1-based last atom of
/// each group), as gromosXX's `In_Topology` does (in_topology.cc:433-459): none → one group
/// over everything; a last NRE beyond the last atom is an error; one before the last atom gets
/// a trailing group appended (and a note). NRE must be positive and strictly increasing
/// (in_parameter.cc:978-993).
/// The individual-λ table of a run over `n_groups` energy groups: the identity, with the rows
/// of a `LAMBDAS` block that applies (NTIL on, NTG ≠ 0) written in file order (a later row for
/// the same kind and pair overwrites, `in_parameter.cc:4459-4468`).
fn individual_lambdas(imd: &ImdParameters, n_groups: usize) -> Result<IndividualLambdas, RunError> {
    let mut t = IndividualLambdas::identity(n_groups);
    let Some(block) = imd.lambdas.as_ref().filter(|l| l.applies(imd.ntg)) else {
        return Ok(t);
    };
    for r in &block.rows {
        if r.g1 == 0 || r.g2 < r.g1 || r.g2 > imd.nre.len().max(1) {
            return Err(RunError::Recipe(format!(
                "LAMBDAS row {} {} {}: the groups must satisfy 1 <= NILG1 <= NILG2 <= NEGR ({})",
                r.kind.gromos_name(),
                r.g1,
                r.g2,
                imd.nre.len().max(1)
            )));
        }
        t.set_row(r.kind, r.g1 - 1, r.g2 - 1, r.coefficients);
    }
    Ok(t)
}

/// gromosXX `in_perturbation.cc:183-300`: with NTC = 3 every perturbed bond becomes a perturbed
/// distance constraint (perturbed SHAKE holds it at its λ-mixed length); with NTC = 2 perturbed
/// bonds stay flexible, and one that perturbs a BONDH bond (which NTC = 2 constrains) is an
/// error there ("Perturbation of a non-existing bond"). LINCS has no perturbed variant: a perturbed bond it
/// would own is refused instead of being left unconstrained and without a force, which is what
/// gromosXX does (it only warns, `create_constraints.cc:139-143`).
fn hand_perturbed_bonds_to_constraints(
    topology: &mut Topology,
    imd: &ImdParameters,
    sel: &crate::ConstraintSelection,
) -> Result<(), RunError> {
    let p = &topology.perturbed_solute;
    if p.bonds.is_empty() || !sel.solute_constrained() {
        return Ok(());
    }
    match imd.ntc {
        3 if sel.solute_lincs => Err(RunError::Recipe(format!(
            "NTC = 3 with LINCS for the solute and {} perturbed bonds: LINCS has no perturbed \
             variant (no λ-dependent lengths, no dH/dλ); use SHAKE",
            p.bonds.len()
        ))),
        3 => {
            let bonds = std::mem::take(&mut topology.perturbed_solute.bonds);
            topology.perturbed_solute.distance_constraints = bonds;
            Ok(())
        },
        2 => {
            // the bonds `ShakeBuffers` constrains for NTC = 2 (NtcMode::HydrogenBonds)
            if let Some(b) = p.bonds.iter().find(|b| b.hydrogen) {
                return Err(RunError::Recipe(format!(
                    "NTC = 2 constrains the bond {}-{} to a hydrogen, which the perturbation \
                     topology perturbs: gromosXX has no perturbed constraint for NTC = 2 \
                     (\"Perturbation of a non-existing bond\")",
                    b.i + 1,
                    b.j + 1
                )));
            }
            Ok(())
        },
        _ => Ok(()),
    }
}

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

/// `(NSNB, NTWG)` when a perturbed run with PRECALCLAM writes λ_k rows without the long-range
/// terms on some output steps: a twin-range pairlist (RCUTP < RCUTL, NSNB > 1) keeps the
/// long-range pairs from its last update, and their λ_k terms are evaluated only when that
/// update is also an NTWG step (`perturbed_nonbonded_set.cc:87-96`). That happens on some
/// output step exactly when neither of NSNB and NTWG divides the other.
pub fn precalclam_longrange_gaps(imd: &ImdParameters) -> Option<(usize, usize)> {
    let precalc = imd.ntg != 0 && imd.precalclam.is_some_and(|b| b.nr_lambdas > 0);
    let twin_range = imd.rcutp < imd.rcutl - 1e-10;
    let (nsnb, ntwg) = (imd.nsnb, imd.ntwg);
    (precalc && twin_range && newen_logic::schedule::output_misses_longrange(nsnb, ntwg))
        .then_some((nsnb, ntwg))
}

/// The warning for [`PrepareNote::NoseHooverVariablesMissing`] (`md` and Python print the
/// same), or `None` when `recipe` does not ask for it: a thermostat other than Nosé–Hoover, or
/// gromosXX's restart (`nh_restart = "gromosxx"`), whose plain Nosé–Hoover never reads them.
pub fn nose_hoover_missing_warning(recipe: &crate::RunRecipe) -> Option<String> {
    use crate::recipe::{NhRestart, ThermostatAlgorithm};
    let t = recipe.ensemble.thermostat.as_ref()?;
    (t.algorithm != ThermostatAlgorithm::Berendsen && t.nh_restart == NhRestart::Intended).then(
        || {
            "INITIALISE NTINHT = 0 asks for the Nose-Hoover thermostat variables from the \
             configuration, but it has no NHCVARIABLES block: they start at zero"
                .to_string()
        },
    )
}

/// The warning for [`PrepareNote::PrecalclamLongRangeGaps`] (`md` and Python print the same).
pub fn precalclam_longrange_warning(nsnb: usize, ntwg: usize) -> String {
    format!(
        "PRECALCLAM with a twin-range pairlist updated every {nsnb} steps (NSNB) and \
         free-energy output every {ntwg} steps (NTWG): as in gromosXX, the lambda_k energies \
         and derivatives written on an output step that follows a pairlist step which was not \
         an output step miss the long-range contribution (the dH/dlambda totals and the forces \
         keep it), so BAR and MBAR from these rows are biased; set NTWG to a multiple of NSNB"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_io::coordinate::read_coordinates;
    use newen_io::topology::{build_topology, read_topology_file};
    use std::path::PathBuf;

    fn shared() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../newen-md/tests/gromosXX_references/shared")
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

    /// NTINHT = 0 without an NHCVARIABLES block is noted (variables start at zero), reported
    /// under the intended restart only; NTINHT = 1, a block, or no Nosé–Hoover say nothing.
    #[test]
    fn missing_nose_hoover_variables_are_noted() {
        let refs = shared().join("..");
        let imd = newen_io::imd::read_imd_file(
            refs.join("water_216_nvt_nosehoover/water_216_nvt_nosehoover.in"),
        )
        .unwrap();
        let run = |imd: &ImdParameters, conf: &str| {
            let parsed = read_topology_file(refs.join("water_216_box/water_216_box.topo")).unwrap();
            let constants = parsed.physical_constants;
            let coords: Coordinates = read_coordinates(refs.join(conf)).unwrap().into();
            let p = prepare_system(
                imd,
                build_topology(parsed),
                constants,
                coords,
                &RunInputs::default(),
            )
            .unwrap();
            p.notes
                .iter()
                .any(|n| matches!(n, PrepareNote::NoseHooverVariablesMissing))
        };
        let fresh = "water_216_box/water_216_box.conf";
        assert!(run(&imd, fresh));
        let mut reset = imd.clone();
        reset.ntinht = 1;
        assert!(!run(&reset, fresh));
        assert!(!run(
            &imd,
            "water_216_nvt_nosehoover_cont_fixed/water_216_nvt_nosehoover_cont_fixed.conf"
        ));
        let mut weak = imd.clone();
        weak.temp_bath[0].algorithm = 0;
        assert!(!run(&weak, fresh));

        let (mut recipe, _) = crate::RunRecipe::from_imd(&imd).unwrap();
        assert!(nose_hoover_missing_warning(&recipe)
            .unwrap()
            .contains("no NHCVARIABLES block: they start at zero"));
        recipe
            .set("ensemble.thermostat.nh_restart", "gromosxx")
            .unwrap();
        assert!(nose_hoover_missing_warning(&recipe).is_none());
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

    #[test]
    fn a_pertdata_lambda_outside_zero_one_is_a_warning() {
        let imd = |nrdgl: i32| ImdParameters {
            nsm: 20,
            ntg: 1,
            nrdgl,
            rlam: 0.5,
            ..Default::default()
        };
        let note = |imd: &ImdParameters, pertdata: Option<f64>| {
            let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
            let constants = parsed.physical_constants;
            let mut coords: Coordinates = read_coordinates(shared().join("aladip.conf"))
                .unwrap()
                .into();
            coords.pertdata = pertdata;
            let inputs = RunInputs {
                pttopo: Some(shared().join("aladip.pttopo")),
                ..Default::default()
            };
            let p = prepare_system(imd, build_topology(parsed), constants, coords, &inputs);
            p.unwrap().notes.iter().find_map(|n| match n {
                PrepareNote::LambdaOutOfRange { lambda } => Some(*lambda),
                _ => None,
            })
        };
        // aladip_vacuum_fep_slowgrowth_cont's continuation would start at its last λ, 1.125
        assert_eq!(note(&imd(1), Some(1.125)), Some(1.125));
        assert_eq!(note(&imd(1), Some(-0.5)), Some(-0.5));
        assert_eq!(note(&imd(1), Some(1.0 + 1e-12)), None);
        assert_eq!(note(&imd(1), Some(0.525)), None);
        assert_eq!(
            note(&imd(0), Some(1.125)),
            None,
            "NRDGL = 0: RLAM, not PERTDATA"
        );
    }

    #[test]
    fn a_net_charge_change_under_reaction_field_is_a_warning() {
        // CB of the neutral dipeptide: +1 in state A, 0 in state B
        let ptp = |qa: f64, qb: f64| {
            format!(
                "TITLE\ncharge\nEND\nPERTATOMPARAM\n1\n    1 1 CB 14 15.0350 {qa:6.3} 14 15.0350 \
                 {qb:6.3} 0.0 0.0\nEND\n"
            )
        };
        let run = |text: String, nlrele: i32, tag: &str| {
            let path = tmp_ptp(tag, &text);
            let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
            let constants = parsed.physical_constants;
            let coords: Coordinates = read_coordinates(shared().join("aladip.conf"))
                .unwrap()
                .into();
            let imd = ImdParameters {
                nsm: 20,
                ntg: 1,
                nlrele,
                ..Default::default()
            };
            let inputs = RunInputs {
                pttopo: Some(path.clone()),
                ..Default::default()
            };
            let p = prepare_system(&imd, build_topology(parsed), constants, coords, &inputs);
            std::fs::remove_file(path).ok();
            p.unwrap().notes
        };
        let charge = |notes: &[PrepareNote]| {
            notes.iter().find_map(|n| match n {
                PrepareNote::PerturbationNetCharge { charge } => Some(*charge),
                _ => None,
            })
        };
        let changed = |notes: &[PrepareNote]| {
            notes.iter().find_map(|n| match n {
                PrepareNote::NetChargeChanged { delta, nlrele } => Some((*delta, *nlrele)),
                _ => None,
            })
        };

        let notes = run(ptp(1.0, 0.0), 1, "plus1");
        let q = charge(&notes).expect("net charge reported");
        assert!((q.a - 1.0).abs() < 1e-9 && q.b.abs() < 1e-9, "{q:?}");
        let (delta, nlrele) = changed(&notes).expect("warned");
        assert!((delta + 1.0).abs() < 1e-9 && nlrele == 1);
        // the other reaction-field modes warn too; no electrostatics or a lattice sum does not
        assert!(changed(&run(ptp(1.0, 0.0), -1, "plus1_lserf")).is_some());
        assert!(changed(&run(ptp(1.0, 0.0), 0, "plus1_none")).is_none());
        assert!(changed(&run(ptp(1.0, 0.0), 2, "plus1_ewald")).is_none());
        // equal end states: ΔQ = 0, no warning
        let notes = run(ptp(0.0, 0.0), 1, "neutral");
        assert_eq!(charge(&notes).unwrap().delta(), 0.0);
        assert!(changed(&notes).is_none());
    }

    #[test]
    fn precalclam_with_interleaved_pairlist_and_output_steps_is_a_warning() {
        let gaps = |imd: &ImdParameters| {
            let parsed = read_topology_file(shared().join("ch4_spc.top")).unwrap();
            let constants = parsed.physical_constants;
            let coords: Coordinates = read_coordinates(shared().join("ch4_spc.cnf"))
                .unwrap()
                .into();
            let inputs = RunInputs {
                pttopo: Some(shared().join("ch4_spc_dummy.ptp")),
                ..Default::default()
            };
            let p = prepare_system(imd, build_topology(parsed), constants, coords, &inputs);
            p.unwrap().notes.iter().find_map(|n| match n {
                PrepareNote::PrecalclamLongRangeGaps { nsnb, ntwg } => Some((*nsnb, *ntwg)),
                _ => None,
            })
        };
        let twinrange = newen_io::imd::read_imd_file(
            shared().join("../ch4_water_fep_precalc_twinrange/ch4_water_fep_precalc_twinrange.in"),
        )
        .unwrap();
        assert_eq!((twinrange.nsnb, twinrange.ntwg), (5, 2));
        assert_eq!(gaps(&twinrange), Some((5, 2)));
        let with = |nsnb: usize, ntwg: usize| ImdParameters {
            nsnb,
            ntwg,
            ..twinrange.clone()
        };
        assert_eq!(gaps(&with(5, 10)), None, "NTWG a multiple of NSNB");
        assert_eq!(
            gaps(&with(5, 1)),
            None,
            "every pairlist step is an output step"
        );
        assert_eq!(
            gaps(&with(1, 2)),
            None,
            "NSNB 1: the pairlist is rebuilt every step"
        );
        assert_eq!(gaps(&with(4, 6)), Some((4, 6)));
        let single_range = ImdParameters {
            rcutl: twinrange.rcutp,
            ..twinrange.clone()
        };
        assert_eq!(gaps(&single_range), None);
        let no_precalc = ImdParameters {
            precalclam: None,
            ..twinrange.clone()
        };
        assert_eq!(gaps(&no_precalc), None);
        let text = precalclam_longrange_warning(5, 2);
        assert!(text.contains("every 5 steps (NSNB)") && text.contains("every 2 steps (NTWG)"));
    }

    fn try_aladip(imd: &ImdParameters) -> Result<Prepared, RunError> {
        let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
        let constants = parsed.physical_constants;
        let topo = build_topology(parsed);
        let coords: Coordinates = read_coordinates(shared().join("aladip.conf"))
            .unwrap()
            .into();
        let inputs = RunInputs {
            pttopo: Some(shared().join("aladip.pttopo")),
            ..Default::default()
        };
        prepare_system(imd, topo, constants, coords, &inputs)
    }

    fn perturbed(ntc: i32, ntcp: i32) -> ImdParameters {
        ImdParameters {
            nsm: 20,
            ntg: 1,
            rlam: 0.125,
            ntc,
            ntcp,
            ..Default::default()
        }
    }

    #[test]
    fn ntc_3_hands_the_perturbed_bonds_to_perturbed_shake() {
        let p = aladip(&perturbed(3, 1), true);
        let ps = &p.topology.perturbed_solute;
        assert!(ps.bonds.is_empty(), "no flexible perturbed bond left");
        assert_eq!(ps.distance_constraints.len(), 2);
        // NTC = 1: they stay flexible
        let p = aladip(&perturbed(1, 1), true);
        assert_eq!(p.topology.perturbed_solute.bonds.len(), 2);
        assert!(p.topology.perturbed_solute.distance_constraints.is_empty());
    }

    #[test]
    fn lincs_cannot_own_a_perturbed_bond() {
        let err = try_aladip(&perturbed(3, 2)).err().expect("refused");
        assert!(err.to_string().contains("LINCS has no perturbed"), "{err}");
    }

    #[test]
    fn ntc_2_keeps_heavy_perturbed_bonds_flexible() {
        // aladip.pttopo perturbs two bonds between heavy atoms: gromosXX keeps them as perturbed
        // bonds under NTC = 2 (with a warning), and so does prepare_system
        let p = aladip(&perturbed(2, 1), true);
        assert_eq!(p.topology.perturbed_solute.bonds.len(), 2);
        assert!(p.topology.perturbed_solute.distance_constraints.is_empty());
    }

    #[test]
    fn ntc_2_constrains_the_bondh_bonds_whatever_the_mixed_masses() {
        // aladip.pttopo takes H11 from 1.008 to 2.008: at RLAM = 1 both atoms of BONDH 10-11 are
        // heavier than 2, and gromosXX still constrains it (BONDH membership, not masses)
        let mut imd = perturbed(2, 1);
        imd.rlam = 1.0;
        let p = aladip(&imd, true);
        assert!(p.topology.mass[10] > 2.0 && p.topology.mass[9] > 2.0);
        let bondh: Vec<(usize, usize)> = p.topology.moltypes[0]
            .bonds
            .iter()
            .filter(|b| b.hydrogen)
            .map(|b| (b.i, b.j))
            .collect();
        assert_eq!(bondh, vec![(3, 4), (9, 10)]);
        for (i, j) in bondh {
            assert!(p.topology.constrained_bonds.contains(&(i, j)));
        }
        assert_eq!(
            p.topology.constrained_bonds.len(),
            4,
            "the two BONDH bonds, both ways"
        );
    }

    #[test]
    fn ntc_2_refuses_a_perturbed_bondh_bond() {
        // aladip_hblock.pttopo perturbs BONDH 4-5: "Perturbation of a non-existing bond"
        let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
        let constants = parsed.physical_constants;
        let coords: Coordinates = read_coordinates(shared().join("aladip.conf"))
            .unwrap()
            .into();
        let inputs = RunInputs {
            pttopo: Some(shared().join("aladip_hblock.pttopo")),
            ..Default::default()
        };
        let err = prepare_system(
            &perturbed(2, 1),
            build_topology(parsed),
            constants,
            coords,
            &inputs,
        )
        .err()
        .expect("refused");
        assert!(err.to_string().contains("bond 4-5"), "{err}");
    }

    /// An EDS run of the vacuum aladip with `pttopo` (in `shared/`, or a file written here).
    fn eds_prepare(
        ntg: i32,
        n_states: usize,
        pttopo: Option<PathBuf>,
    ) -> Result<Prepared, RunError> {
        let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
        let constants = parsed.physical_constants;
        let coords: Coordinates =
            read_coordinates(shared().join("../aladip_vacuum/aladip_vacuum.conf"))
                .unwrap()
                .into();
        let imd = ImdParameters {
            ntg,
            eds: Some(newen_io::imd::EdsBlock {
                enabled: true,
                form: 1,
                num_states: n_states,
                s: vec![0.5],
                pairs: Vec::new(),
                eir: vec![0.0; n_states],
            }),
            ..Default::default()
        };
        let inputs = RunInputs {
            pttopo,
            ..Default::default()
        };
        prepare_system(&imd, build_topology(parsed), constants, coords, &inputs)
    }

    fn tmp_ptp(tag: &str, text: &str) -> PathBuf {
        let p =
            std::env::temp_dir().join(format!("newen_prepare_{tag}_{}.ptp", std::process::id()));
        std::fs::write(&p, text).unwrap();
        p
    }

    #[test]
    fn eds_merges_mpertatom_and_sizes_the_end_states() {
        let p = eds_prepare(0, 2, Some(shared().join("aladip_eds2.mpt"))).unwrap();
        let eds = &p.topology.eds_solute;
        assert_eq!(eds.num_states, 2);
        assert_eq!(
            eds.atoms.iter().map(|a| a.seq).collect::<Vec<_>>(),
            [0, 3, 10]
        );
        // CB (1) owns its excluded partners 2, 3 and N (4), an EDS atom listed after it
        assert!(eds.atoms[0].exclusions.contains(&3));
        assert!(!eds.atoms[1].exclusions.contains(&0));
        assert!(p.topology.exclusions[0].is_empty() && p.topology.is_excluded_or_14(0, 3));
        let c = &p.configuration;
        assert_eq!(c.current().energies.eds_vi.len(), 2);
        assert_eq!(c.old().perturbed_energy_derivatives.eds_vi.len(), 2);
        assert_eq!(c.eds.num_states(), 2);
        assert!(p.notes.iter().any(|n| matches!(
            n,
            PrepareNote::EdsLoaded {
                atoms: 3,
                states: 2
            }
        )));
    }

    #[test]
    fn eds_refuses_inconsistent_perturbation_topologies() {
        let err = |r: Result<Prepared, RunError>| r.err().expect("refused").to_string();
        // no @pttopo
        let e = err(eds_prepare(0, 2, None));
        assert!(e.contains("pttopo"), "{e}");
        // NPTB ≠ NUMSTATES
        let e = err(eds_prepare(0, 3, Some(shared().join("aladip_eds2.mpt"))));
        assert!(e.contains("do not match"), "{e}");
        // no MPERTATOM
        let e = err(eds_prepare(1, 2, Some(shared().join("aladip.pttopo"))));
        assert!(e.contains("no MPERTATOM"), "{e}");
        // PERT* blocks without perturbation
        let mpert = std::fs::read_to_string(shared().join("aladip_eds2.mpt")).unwrap();
        let bonded = "PERTBONDSTRETCH\n1\n2 4 26 20\nEND\n";
        let path = tmp_ptp("pertblock", &format!("{mpert}{bonded}"));
        let e = err(eds_prepare(0, 2, Some(path.clone())));
        assert!(e.contains("PERTBONDSTRETCH not supported for eds"), "{e}");
        // …which a perturbed run reads
        assert!(eds_prepare(1, 2, Some(path.clone())).is_ok());
        std::fs::remove_file(path).ok();
        // perturbed atoms next to EDS atoms (TI-EDS) run…
        let atom = "PERTATOMPARAM\n1\n7 1 C 12 12.011 0.38 12 12.011 0.38 0.5 0.5\nEND\n";
        let path = tmp_ptp("pertatom", &format!("{mpert}{atom}"));
        assert!(eds_prepare(1, 2, Some(path.clone())).is_ok());
        std::fs::remove_file(path).ok();
        // …but not an atom in both PERTATOMPARAM and MPERTATOM (N, atom 4)…
        let atom = "PERTATOMPARAM\n1\n4 1 N 5 14.0067 -0.28 5 14.0067 0.1 0.5 0.5\nEND\n";
        let path = tmp_ptp("pertatom_eds", &format!("{mpert}{atom}"));
        let e = err(eds_prepare(1, 2, Some(path.clone())));
        assert!(e.contains("both PERTATOMPARAM and MPERTATOM"), "{e}");
        std::fs::remove_file(path).ok();
        // …nor a perturbed atom pair of an EDS atom
        let pair = "PERTATOMPAIR\n1\n4 10 1 0\nEND\n";
        let path = tmp_ptp("pertpair_eds", &format!("{mpert}{pair}"));
        let e = err(eds_prepare(1, 2, Some(path.clone())));
        assert!(e.contains("PERTATOMPAIR 4-10 involves an EDS atom"), "{e}");
        std::fs::remove_file(path).ok();
        // MPERTATOM without EDS
        let mut imd = perturbed(1, 1);
        imd.nsm = 0;
        let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
        let constants = parsed.physical_constants;
        let coords: Coordinates =
            read_coordinates(shared().join("../aladip_vacuum/aladip_vacuum.conf"))
                .unwrap()
                .into();
        let inputs = RunInputs {
            pttopo: Some(shared().join("aladip_eds2.mpt")),
            ..Default::default()
        };
        let e = err(prepare_system(
            &imd,
            build_topology(parsed),
            constants,
            coords,
            &inputs,
        ));
        assert!(e.contains("MPERTATOM") && e.contains("no EDS block"), "{e}");
    }

    #[test]
    fn nrdgl_1_takes_lambda_from_pertdata() {
        let run = |nrdgl: i32, pertdata: Option<f64>| {
            let parsed = read_topology_file(shared().join("aladip.topo")).unwrap();
            let constants = parsed.physical_constants;
            let mut coords: Coordinates = read_coordinates(shared().join("aladip.conf"))
                .unwrap()
                .into();
            coords.pertdata = pertdata;
            let mut imd = perturbed(1, 1);
            imd.nrdgl = nrdgl;
            let inputs = RunInputs {
                pttopo: Some(shared().join("aladip.pttopo")),
                ..Default::default()
            };
            prepare_system(&imd, build_topology(parsed), constants, coords, &inputs)
        };
        assert_eq!(run(1, Some(0.5)).unwrap().topology.lambda, 0.5);
        assert_eq!(run(0, Some(0.5)).unwrap().topology.lambda, 0.125);
        let err = run(1, None).err().expect("PERTDATA required");
        assert!(err.to_string().contains("PERTDATA"), "{err}");
    }

    #[test]
    fn lambdas_mass_row_mixes_with_its_groups_lambda() {
        use newen_io::imd::{LambdasBlock, LambdasRow};
        let mut imd = perturbed(1, 1);
        imd.negr = 2;
        imd.nre = vec![6, 72];
        imd.lambdas = Some(LambdasBlock {
            enabled: true,
            rows: vec![LambdasRow {
                kind: newen_core::lambdas::InteractionLambda::Mass,
                g1: 1,
                g2: 1,
                coefficients: [0.0, 1.0, -1.0, 1.0, 0.0],
            }],
        });
        let p = aladip(&imd, true);
        let t = &p.topology;
        let lam: f64 = 0.125;
        let (b, c, d) = (1.0, -1.0, 1.0);
        let l = b * lam * lam * lam + c * lam * lam + d * lam;
        for a in &t.perturbed_solute.atoms {
            let l_atom = if t.energy_group(a.seq) == 0 { l } else { lam };
            assert_eq!(t.mass[a.seq], (1.0 - l_atom) * a.a_mass + l_atom * a.b_mass);
        }
        // off: the identity, every mass at RLAM
        imd.lambdas.as_mut().unwrap().enabled = false;
        let p = aladip(&imd, true);
        for a in &p.topology.perturbed_solute.atoms {
            assert_eq!(
                p.topology.mass[a.seq],
                (1.0 - lam) * a.a_mass + lam * a.b_mass
            );
        }
    }
}
