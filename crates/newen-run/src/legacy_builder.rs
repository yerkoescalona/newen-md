//! The step-1 builder, frozen as a **test oracle**: `build_sequence_from_imd` must produce a
//! sequence whose energies, forces and positions are bit-identical to this one on every
//! reference system (`build::tests::plan_matches_legacy_builder`). Compiled only under
//! `cfg(test)`; deleted in 3.9 step 4 once the recipe path has carried the suites for a
//! release.
#![allow(dead_code, unused_imports)]

use newen_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::{BoxType, Configuration},
    math::{Periodicity, Rectangular, Triclinic, Vacuum, Vec3},
    pairlist::{PairlistAlgorithm, PairlistContainer},
    Topology,
};
use newen_forces::{
    nonbonded::CRFParameters,
    restraints::{
        DistanceRestraint, DistanceRestraints, PerturbedDistanceRestraint,
        PerturbedDistanceRestraints, PositionRestraint, PositionRestraints,
    },
};
use newen_integrators::{
    algorithms::{
        BerendsenBarostat, BerendsenBarostatParams, BerendsenThermostat, EnergyCalculation,
        Forcefield, LatticeShift, LeapFrogPosition, LeapFrogVelocity, LincsAlgorithm,
        NoseHooverThermostat, PressureCalculation, RemoveCOMMotion, SettleAlgorithm,
        ShakeAlgorithm, SteepestDescentAlgorithm, TemperatureCalculation, VirialType,
    },
    constraints::{NtcMode, ShakeParameters},
};
use newen_io::{
    distanceres::read_distanceres,
    imd::ImdParameters,
    posres::{build_posres_entries, read_posresspec, read_refpos},
};

use crate::{total_dof, ConstraintSelection, Prepared, RunError, RunInputs, RunOptions};

/// What the thermostat was built as — for the binary's report and the Python getters.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermostatSummary {
    /// `"Berendsen"`, `"Nose-Hoover"` or `"Nose-Hoover-Chain(n)"`.
    pub label: String,
    pub temperature: f64,
    pub tau: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BarostatSummary {
    pub pressure0: f64,
    pub tau: f64,
    pub compressibility: f64,
}

/// Facts about the built sequence that the caller may report. Never consulted by the
/// sequence itself.
#[derive(Debug, Clone)]
pub struct BuildSummary {
    pub is_minimization: bool,
    pub constraints: ConstraintSelection,
    /// SHAKE parameters when SHAKE is in the sequence.
    pub shake: Option<ShakeParameters>,
    /// LINCS expansion orders (solute, solvent).
    pub lincs_orders: (usize, usize),
    pub thermostat: Option<ThermostatSummary>,
    pub barostat: Option<BarostatSummary>,
    /// Number of position-restrained atoms.
    pub position_restraints: usize,
    /// (unperturbed, perturbed) distance restraints.
    pub distance_restraints: (usize, usize),
    /// Pairs in the initial pairlist.
    pub initial_pairs: usize,
    pub parallel_nonbonded: bool,
    /// Kinetic degrees of freedom (see [`crate::total_dof`]); what the thermostat couples to
    /// and what the temperature must be reported with.
    pub total_dof: f64,
}

pub struct Built {
    pub sequence: AlgorithmSequence,
    pub summary: BuildSummary,
}

struct Restraints {
    position: Option<PositionRestraints>,
    distance: Option<(DistanceRestraints, PerturbedDistanceRestraints)>,
}

/// Load position/distance restraints exactly as the binary's `@posresspec`/`@refpos`/
/// `@distrest` handling does (same file formats, same NTPOR/NTPORB/NTDIR dispatch).
fn load_restraints(
    imd: &ImdParameters,
    inputs: &RunInputs,
    positions: &[Vec3],
) -> Result<Restraints, RunError> {
    let position = if imd.ntpor > 0 {
        let por_file = inputs
            .posresspec
            .as_ref()
            .ok_or_else(|| RunError::MissingInput {
                flag: "posresspec",
                reason: format!("NTPOR={}", imd.ntpor),
            })?;
        let restrained_atoms = read_posresspec(por_file).map_err(|e| RunError::Io {
            what: "position restraint spec".into(),
            source: e,
        })?;
        // NTPORB=1: reference positions from @refpos; NTPORB=0: the startup configuration.
        let ref_positions = if imd.ntporb >= 1 {
            let rpr_file = inputs
                .refpos
                .as_ref()
                .ok_or_else(|| RunError::MissingInput {
                    flag: "refpos",
                    reason: format!("NTPORB={}", imd.ntporb),
                })?;
            read_refpos(rpr_file).map_err(|e| RunError::Io {
                what: "reference positions".into(),
                source: e,
            })?
        } else {
            positions.to_vec()
        };
        let entries = build_posres_entries(&restrained_atoms, &ref_positions);
        let mut pr = PositionRestraints::new();
        for entry in &entries {
            pr.add_restraint(PositionRestraint::new(
                entry.atom,
                entry.reference_pos,
                imd.cpor,
            ));
        }
        Some(pr)
    } else {
        None
    };

    let distance = if imd.ntdir != 0 {
        let dr_file = inputs
            .distrest
            .as_ref()
            .ok_or_else(|| RunError::MissingInput {
                flag: "distrest",
                reason: format!("NTDIR={}", imd.ntdir),
            })?;
        let (unperturbed, perturbed) = read_distanceres(dr_file).map_err(|e| RunError::Io {
            what: "distance restraint file".into(),
            source: e,
        })?;
        let mode = imd.ntdir.abs();
        let k = imd.cdir;
        let r_linear = imd.dir0;

        let mut dr = DistanceRestraints::new();
        for spec in &unperturbed {
            let mut r = DistanceRestraint::new(
                spec.atom1, spec.atom2, spec.r0, spec.w0, spec.rah, k, r_linear, mode,
            );
            if imd.ntdir < 0 {
                r = r.with_time_averaging(imd.taudir, imd.dt);
            }
            dr.add(r);
        }
        let mut pdr = PerturbedDistanceRestraints::new();
        for spec in perturbed.iter().filter(|_| imd.ntg != 0) {
            let mut r = PerturbedDistanceRestraint::new(
                spec.atom1, spec.atom2, spec.n, spec.m, spec.a_r0, spec.b_r0, spec.a_w0, spec.b_w0,
                spec.rah, k, r_linear, mode,
            );
            if imd.ntdir < 0 {
                r = r.with_time_averaging(imd.taudir, imd.dt);
            }
            pdr.add(r);
        }
        Some((dr, pdr))
    } else {
        None
    };

    Ok(Restraints { position, distance })
}

fn shake_parameters(imd: &ImdParameters, sel: &ConstraintSelection) -> ShakeParameters {
    let ntc = if sel.solute_shake {
        ConstraintSelection::ntc_mode(imd)
    } else {
        NtcMode::SolventOnly
    };
    ShakeParameters {
        tolerance: imd.shake_tol,
        max_iterations: 1000,
        ntc,
    }
}

/// SHAKE → SETTLE → LINCS, in the binary's order.
fn push_constraint_algorithms(
    seq: &mut AlgorithmSequence,
    imd: &ImdParameters,
    sel: &ConstraintSelection,
) -> Option<ShakeParameters> {
    let mut shake = None;
    if sel.shake_enabled() {
        let params = shake_parameters(imd, sel);
        let mut alg = ShakeAlgorithm::new(params.clone());
        alg.include_solvent = sel.solvent_shake;
        // NTISHK=1: shake initial positions; NTISHK=2: positions + velocities.
        if imd.ntishk >= 1 {
            alg.shake_initial_positions = true;
        }
        if imd.ntishk >= 2 {
            alg.shake_initial_velocities = true;
        }
        seq.push(Box::new(alg));
        shake = Some(params);
    }
    if sel.settle_enabled {
        seq.push(Box::new(SettleAlgorithm::new()));
    }
    if sel.lincs_enabled {
        seq.push(Box::new(LincsAlgorithm::new(
            ConstraintSelection::ntc_mode(imd),
            imd.lincs_order_solute,
            imd.lincs_order_solvent,
            sel.solute_lincs,
            sel.solvent_lincs,
        )));
    }
    shake
}

/// The periodicity a prepared system runs under: triclinic for a truncated-octahedron box,
/// vacuum for an all-zero box, rectangular otherwise. Exposed for callers that attach extra
/// algorithms needing the same boundary conditions as `Forcefield` (e.g. an ML term).
fn periodicity_of(prepared: &Prepared) -> Periodicity {
    let b = prepared.box_dims;
    if let Some(triclinic_box) = prepared.truncoct_box {
        Periodicity::Triclinic(Triclinic::new(triclinic_box))
    } else if b.x == 0.0 && b.y == 0.0 && b.z == 0.0 {
        Periodicity::Vacuum(Vacuum)
    } else {
        Periodicity::Rectangular(Rectangular::new(b))
    }
}

fn virial_type(imd: &ImdParameters) -> VirialType {
    match imd
        .pressure_parameters
        .as_ref()
        .map(|p| p.virial)
        .unwrap_or(0)
    {
        2 => VirialType::Molecular,
        1 => VirialType::Atomic,
        _ => VirialType::None,
    }
}

/// Build the algorithm sequence for `imd` on a prepared system.
pub(crate) fn legacy_build_sequence(
    imd: &ImdParameters,
    prepared: &Prepared,
    inputs: &RunInputs,
    options: &RunOptions,
) -> Result<Built, RunError> {
    let topo = &prepared.topology;
    let conf = &prepared.configuration;
    let n_atoms = topo.num_atoms();

    // Nonbonded interactions: CRF interior dielectric is always 1 in GROMOS.
    let crf_params = CRFParameters::new(imd.rcutl, 1.0, imd.epsrf, imd.appak);
    let lj_params = Forcefield::convert_lj_parameters(topo);

    // Pairlist. PAIRLIST block SIZE; 0.0 means "auto" and is resolved by the cell-list
    // algorithm to half the short cutoff, as gromosXX does.
    let mut pairlist = PairlistContainer::new(imd.rcutp, imd.rcutl, 0.0);
    pairlist.update_frequency = imd.nsnb;
    pairlist.grid_size = imd.size;

    let periodicity = periodicity_of(prepared);
    let box_type = match &periodicity {
        Periodicity::Rectangular(_) => BoxType::Rectangular,
        Periodicity::Triclinic(_) => BoxType::Triclinic,
        Periodicity::Vacuum(_) => BoxType::Vacuum,
    };
    let pairlist_algorithm = PairlistAlgorithm::from_imd(
        imd.algorithm,
        n_atoms,
        box_type,
        !topo.chargegroups.is_empty(),
        imd.type_,
    );
    pairlist_algorithm.update(topo, conf, &mut pairlist, &periodicity);
    let initial_pairs = pairlist.total_pairs();

    let restraints = load_restraints(imd, inputs, &conf.current().pos)?;
    let position_restraints = restraints
        .position
        .as_ref()
        .map(|p| p.restraints.len())
        .unwrap_or(0);
    let distance_restraints = restraints
        .distance
        .as_ref()
        .map(|(d, p)| (d.restraints.len(), p.restraints.len()))
        .unwrap_or((0, 0));

    let is_minimization = imd.ntem > 0;
    let has_solvent = n_atoms > topo.num_solute_atoms();
    let constraints = ConstraintSelection::from_imd(imd, has_solvent);
    let ntc_mode = ConstraintSelection::ntc_mode(imd);
    let dof = total_dof(topo, &constraints, ntc_mode, imd.ndfmin);
    let parallel_nonbonded = options.parallel.resolve(n_atoms);

    // Forcefield (bonded + nonbonded), shared by both branches.
    let mut forcefield = Forcefield::new(
        lj_params,
        crf_params,
        periodicity,
        pairlist,
        pairlist_algorithm,
    );
    // Constants from the topology's PHYSICALCONSTANTS block (mirrors gromosXX).
    forcefield.four_pi_eps_i = prepared.physical_constants.four_pi_eps_i;
    forcefield.ntf_bond = imd.ntf[0] != 0;
    forcefield.ntf_angle = imd.ntf[1] != 0;
    forcefield.ntf_improper = imd.ntf[2] != 0;
    forcefield.ntf_dihedral = imd.ntf[3] != 0;
    forcefield.parallel_nonbonded = parallel_nonbonded;
    if !topo.solvent_atom_template.is_empty() {
        forcefield.atoms_per_solvent = topo.solvent_atom_template.len();
    }
    // The binary only configures the virial for dynamics with pressure coupling.
    if !is_minimization && imd.couple_pressure {
        forcefield.virial_type = virial_type(imd);
    }
    if let Some(pr) = restraints.position {
        forcefield.position_restraints = pr;
    }
    if let Some((dr, pdr)) = restraints.distance {
        forcefield.distance_restraints = dr;
        forcefield.perturbed_distance_restraints = pdr;
    }
    if imd.ntg != 0 {
        forcefield.global_alphlj = imd.alphlj;
        forcefield.global_alphc = imd.alphc;
    }

    let mut seq = AlgorithmSequence::new();
    let shake;
    let mut thermostat = None;
    let mut barostat = None;

    if is_minimization {
        // Forcefield → SteepestDescent → constraints → EnergyCalculation. GROMOS applies
        // constraints during minimisation too; no TemperatureCalculation (E_total = E_pot).
        seq.push(Box::new(forcefield));
        seq.push(Box::new(
            SteepestDescentAlgorithm::new()
                .with_tolerance(imd.dele)
                .with_step_sizes(imd.dx0, imd.dxm)
                .with_min_steps(imd.nmin)
                .with_force_limit(imd.flim),
        ));
        shake = push_constraint_algorithms(&mut seq, imd, &constraints);
        seq.push(Box::new(EnergyCalculation::new()));
    } else {
        // 1. COM motion removal (GROMOS: first in sequence, before the forcefield).
        if imd.nticom >= 1 || imd.nscm != 0 {
            seq.push(Box::new(RemoveCOMMotion::new(imd.nticom, imd.nscm)));
        }
        // 1b. Lattice-shift tracker: charge groups back into the box before the forces
        //     (gromosXX `create_md_sequence.cc`, every non-vacuum boundary).
        if imd.ntb != 0 {
            seq.push(Box::new(LatticeShift::new()));
        }
        // 2. Forcefield.
        seq.push(Box::new(forcefield));
        // 3. Leap-frog velocity step.
        seq.push(Box::new(LeapFrogVelocity::new()));
        // 3b. Thermostat between the velocity and position updates (GROMOS convention).
        //     MULTIBATH algorithm: 0 → Berendsen, 1 → Nosé-Hoover, N ≥ 2 → NH chain of length N.
        let bath = imd.temp_bath.first();
        let temperature = bath.and_then(|b| b.temp0.first().copied()).unwrap_or(300.0);
        let tau = bath.and_then(|b| b.tau.first().copied()).unwrap_or(-1.0);
        if tau > 0.0 {
            let algorithm = bath.map(|b| b.algorithm).unwrap_or(0);
            let nhc_chain = bath.map(|b| b.nhc_chain).unwrap_or(1);
            let label = match algorithm {
                0 => {
                    seq.push(Box::new(BerendsenThermostat::new_single_bath(
                        temperature,
                        tau,
                        dof,
                        n_atoms,
                    )));
                    "Berendsen".to_string()
                },
                1 => {
                    seq.push(Box::new(NoseHooverThermostat::new_single_bath(
                        temperature,
                        tau,
                        dof,
                        n_atoms,
                    )));
                    "Nose-Hoover".to_string()
                },
                n => {
                    let chain = (n as usize).max(nhc_chain).max(2);
                    seq.push(Box::new(NoseHooverThermostat::new_chain_bath(
                        temperature,
                        tau,
                        dof,
                        n_atoms,
                        chain,
                    )));
                    format!("Nose-Hoover-Chain({n})")
                },
            };
            thermostat = Some(ThermostatSummary {
                label,
                temperature,
                tau,
            });
        }
        // 4. Leap-frog position step.
        seq.push(Box::new(LeapFrogPosition::new()));
        // 5. Constraints.
        shake = push_constraint_algorithms(&mut seq, imd, &constraints);
        // 6. Temperature / kinetic energy.
        seq.push(Box::new(TemperatureCalculation::new()));
        // 7. Pressure calculation and barostat (PRESSURESCALE).
        if imd.couple_pressure {
            seq.push(Box::new(PressureCalculation::new(virial_type(imd))));
            let pp = imd.pressure_parameters.as_ref();
            let params = BerendsenBarostatParams {
                pressure0: pp.map(|p| p.pressure0[0][0]).unwrap_or(1.0),
                compressibility: pp.map(|p| p.compressibility[0][0]).unwrap_or(4.575e-4),
                tau: pp.map(|p| p.tau_p).unwrap_or(0.5),
            };
            barostat = Some(BarostatSummary {
                pressure0: params.pressure0,
                tau: params.tau,
                compressibility: params.compressibility,
            });
            seq.push(Box::new(BerendsenBarostat::new(params)));
        }
        // 8. Energy finalisation.
        seq.push(Box::new(EnergyCalculation::new()));
    }

    Ok(Built {
        sequence: seq,
        summary: BuildSummary {
            is_minimization,
            constraints,
            shake,
            lincs_orders: (imd.lincs_order_solute, imd.lincs_order_solvent),
            thermostat,
            barostat,
            position_restraints,
            distance_restraints,
            initial_pairs,
            parallel_nonbonded,
            total_dof: dof,
        },
    })
}
