//! Stages 2 and 3 of the builder: `instantiate` turns a validated plan into runnable
//! algorithms, `start` primes the run (3.9).
//!
//! ```text
//! RunRecipe ── build_plan ──▶ Vec<AlgorithmSpec> ── validate_plan ── instantiate ──▶ AlgorithmSequence ── start
//! ```
//!
//! `instantiate` reads **only** the plan (plus the prepared topology/configuration and the
//! boundary conditions): every physical value it needs was resolved by `build_plan`, so an
//! edited plan is exactly what runs (drift guard G8). The recipe-aware wrappers below add the
//! reporting summary the binary prints.

use std::path::Path;

use gromos_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::{BoxType, Configuration},
    math::{Periodicity, Rectangular, Triclinic, Vacuum},
    pairlist::{PairlistAlgorithm, PairlistContainer},
    Topology,
};
use gromos_forces::{
    nonbonded::CRFParameters,
    restraints::{
        DistanceRestraint, DistanceRestraints, PerturbedDistanceRestraint,
        PerturbedDistanceRestraints, PositionRestraint, PositionRestraints,
    },
};
use gromos_integrators::{
    algorithms::{
        BerendsenBarostat, BerendsenBarostatParams, BerendsenThermostat, BerendsenThermostatParams,
        EnergyCalculation, Forcefield, LatticeShift, LeapFrogPosition, LeapFrogVelocity,
        LincsAlgorithm, NoseHooverThermostat, PressureCalculation, RemoveCOMMotion,
        RottransConstraints, SettleAlgorithm, ShakeAlgorithm, SteepestDescentAlgorithm,
        TemperatureCalculation, VirialType,
    },
    constraints::ShakeParameters,
};
use gromos_io::{
    distanceres::read_distanceres,
    imd::ImdParameters,
    posres::{build_posres_entries, read_posresspec, read_refpos},
};

use crate::plan::{
    build_plan, validate_plan, AlgorithmSpec, DistanceRestraintsPlan, ForcefieldPlan,
    PositionRestraintsPlan,
};
use crate::recipe::{
    ConstraintAlgorithm, Diagnostics, RunRecipe, TermSpec, ThermostatAlgorithm, VirialKind,
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

/// The result of `instantiate`: the sequence plus the counts the caller may report.
pub struct Instantiated {
    pub sequence: AlgorithmSequence,
    pub initial_pairs: usize,
    pub position_restraints: usize,
    pub distance_restraints: (usize, usize),
}

/// A built run: the sequence, the data it came from, and the report.
pub struct Built {
    pub sequence: AlgorithmSequence,
    pub recipe: RunRecipe,
    pub plan: Vec<AlgorithmSpec>,
    pub diagnostics: Diagnostics,
    pub summary: BuildSummary,
}

/// The periodicity a prepared system runs under: triclinic for a truncated-octahedron box,
/// vacuum for an all-zero box, rectangular otherwise.
pub fn periodicity_of(prepared: &Prepared) -> Periodicity {
    let b = prepared.box_dims;
    if let Some(triclinic_box) = prepared.truncoct_box {
        Periodicity::Triclinic(Triclinic::new(triclinic_box))
    } else if b.x == 0.0 && b.y == 0.0 && b.z == 0.0 {
        Periodicity::Vacuum(Vacuum)
    } else {
        Periodicity::Rectangular(Rectangular::new(b))
    }
}

fn virial_type(kind: VirialKind) -> VirialType {
    match kind {
        VirialKind::Molecular => VirialType::Molecular,
        VirialKind::Atomic => VirialType::Atomic,
        VirialKind::None => VirialType::None,
    }
}

fn io_err(what: &str, path: &Path) -> impl FnOnce(gromos_io::IoError) -> RunError {
    let what = format!("{what} '{}'", path.display());
    move |source| RunError::Io { what, source }
}

fn load_position_restraints(
    plan: &PositionRestraintsPlan,
    startup_positions: &[gromos_core::math::Vec3],
) -> Result<PositionRestraints, RunError> {
    let restrained_atoms =
        read_posresspec(&plan.spec).map_err(io_err("position restraint spec", &plan.spec))?;
    // NTPORB=1: reference positions from a file; NTPORB=0: the start-up configuration.
    let ref_positions = match &plan.reference {
        Some(path) => read_refpos(path).map_err(io_err("reference positions", path))?,
        None => startup_positions.to_vec(),
    };
    let entries = build_posres_entries(&restrained_atoms, &ref_positions);
    let mut pr = PositionRestraints::new();
    for entry in &entries {
        pr.add_restraint(PositionRestraint::new(
            entry.atom,
            entry.reference_pos,
            plan.k,
        ));
    }
    Ok(pr)
}

fn load_distance_restraints(
    plan: &DistanceRestraintsPlan,
) -> Result<(DistanceRestraints, PerturbedDistanceRestraints), RunError> {
    let (unperturbed, perturbed) =
        read_distanceres(&plan.file).map_err(io_err("distance restraint file", &plan.file))?;
    let mut dr = DistanceRestraints::new();
    for spec in &unperturbed {
        let mut r = DistanceRestraint::new(
            spec.atom1,
            spec.atom2,
            spec.r0,
            spec.w0,
            spec.rah,
            plan.k,
            plan.r_linear,
            plan.mode,
        );
        if let Some(tau) = plan.time_averaging {
            r = r.with_time_averaging(tau, plan.dt);
        }
        dr.add(r);
    }
    // gromosXX creates the perturbed restraints only for a perturbed run
    // (`create_special.cc:125-133`; `in_distanceres.cc:520-524` ignores PERTDISRESSPEC otherwise)
    let mut pdr = PerturbedDistanceRestraints::new();
    for spec in perturbed.iter().filter(|_| plan.perturbed) {
        let mut r = PerturbedDistanceRestraint::new(
            spec.atom1,
            spec.atom2,
            spec.n,
            spec.m,
            spec.a_r0,
            spec.b_r0,
            spec.a_w0,
            spec.b_w0,
            spec.rah,
            plan.k,
            plan.r_linear,
            plan.mode,
        );
        if let Some(tau) = plan.time_averaging {
            r = r.with_time_averaging(tau, plan.dt);
        }
        pdr.add(r);
    }
    Ok((dr, pdr))
}

struct ForcefieldBuilt {
    forcefield: Forcefield,
    initial_pairs: usize,
    position_restraints: usize,
    distance_restraints: (usize, usize),
}

fn instantiate_forcefield(
    ff: &ForcefieldPlan,
    topo: &Topology,
    conf: &Configuration,
    periodicity: &Periodicity,
) -> Result<ForcefieldBuilt, RunError> {
    let n_atoms = topo.num_atoms();
    // CRF interior dielectric is always 1 in GROMOS; the constants come from RCRF, not RCUTL
    // (gromosXX `nonbonded_term.cc:51-63` reads `nonbonded.rf_cutoff`).
    let crf_params = CRFParameters::new(ff.rf_cutoff, 1.0, ff.epsilon_rf, ff.kappa);
    let lj_params = Forcefield::convert_lj_parameters(topo);

    let mut pairlist = PairlistContainer::new(ff.cutoff_short, ff.cutoff_long, 0.0);
    pairlist.update_frequency = ff.pairlist_every;
    pairlist.grid_size = ff.grid_size;
    let box_type = match periodicity {
        Periodicity::Rectangular(_) => BoxType::Rectangular,
        Periodicity::Triclinic(_) => BoxType::Triclinic,
        Periodicity::Vacuum(_) => BoxType::Vacuum,
    };
    let pairlist_algorithm = PairlistAlgorithm::from_imd(
        ff.pairlist_algorithm,
        n_atoms,
        box_type,
        !topo.chargegroups.is_empty(),
        ff.pairlist_type,
    );
    pairlist_algorithm.update(topo, conf, &mut pairlist, periodicity);
    let initial_pairs = pairlist.total_pairs();

    let mut forcefield = Forcefield::new(
        lj_params,
        crf_params,
        periodicity.clone(),
        pairlist,
        pairlist_algorithm,
    );
    forcefield.four_pi_eps_i = ff.four_pi_eps_i;
    forcefield.rf_excluded = ff.rf_excluded;
    forcefield.covalent_form = gromos_forces::bonded::CovalentForm {
        quartic_bonds: !ff.covalent_form.0,
        cosine_harmonic_angles: !ff.covalent_form.1,
        arbitrary_phase_shifts: !ff.covalent_form.2,
    };
    forcefield.ntf_bond = ff.bonds;
    forcefield.ntf_angle = ff.angles;
    forcefield.ntf_improper = ff.impropers;
    forcefield.ntf_dihedral = ff.dihedrals;
    forcefield.parallel_nonbonded = ff.parallel;
    forcefield.atoms_per_solvent = ff.atoms_per_solvent;
    forcefield.virial_type = virial_type(ff.virial);

    let mut position_restraints = 0;
    if let Some(plan) = &ff.position_restraints {
        let pr = load_position_restraints(plan, &conf.current().pos)?;
        position_restraints = pr.restraints.len();
        forcefield.position_restraints = pr;
    }
    let mut distance_restraints = (0, 0);
    if let Some(plan) = &ff.distance_restraints {
        let (dr, pdr) = load_distance_restraints(plan)?;
        distance_restraints = (dr.restraints.len(), pdr.restraints.len());
        forcefield.distance_restraints = dr;
        forcefield.perturbed_distance_restraints = pdr;
    }
    if let Some(p) = &ff.perturbation {
        // λ, NLAM and the individual λ come with the topology (`prepare_system`)
        forcefield.global_alphlj = p.alpha_lj;
        forcefield.global_alphc = p.alpha_c;
    }

    Ok(ForcefieldBuilt {
        forcefield,
        initial_pairs,
        position_restraints,
        distance_restraints,
    })
}

/// The `orchestrator` entry: one provider per term, each over its own region. Adding a term
/// kind = a `TermSpec` variant (recipe.rs), its registry lines (plan.rs) and one arm here
/// (3.9 step 5).
fn instantiate_orchestrator(
    terms: &[TermSpec],
    topo: &Topology,
    periodicity: &Periodicity,
) -> Result<Box<dyn gromos_core::algorithm::Algorithm>, RunError> {
    use gromos_core::selection::AtomSelection;
    use gromos_forces::orchestrator::ProviderOrchestrator;
    use gromos_forces::orchestrator_algorithm::ProviderOrchestratorAlgorithm;
    use gromos_forces::zones::{Zone, ZonePartition};

    if terms.is_empty() {
        return Err(RunError::InvalidPlan(
            "orchestrator without terms".to_string(),
        ));
    }
    let n_atoms = topo.num_atoms();
    let region_of = |name: &str, region: &str, buffer: Option<&str>| {
        let partition = ZonePartition::from_selections(topo, region, buffer)
            .map_err(|e| RunError::Recipe(format!("{name} region: {e}")))?;
        AtomSelection::from_indices(partition.atoms_in(Zone::Inner), n_atoms)
            .map_err(|e| RunError::Recipe(format!("{name} region: {e}")))
    };

    // Energies are reported per term under the registry name; a repeated kind is indexed
    // (`xtb`, `xtb:1`) so two terms of one kind stay distinguishable (G10).
    let label_of = |index: usize| {
        let name = terms[index].name();
        if terms.iter().filter(|t| t.name() == name).count() > 1 {
            format!("{name}:{index}")
        } else {
            name.to_string()
        }
    };

    let mut orchestrator = ProviderOrchestrator::new();
    for (index, term) in terms.iter().enumerate() {
        match term {
            #[cfg(feature = "ml")]
            TermSpec::Schnet {
                model,
                cutoff,
                elements,
                region,
                buffer,
                ..
            } => {
                use gromos_forces::nonbonded::schnet::SchNetInteraction;
                let selection = region_of("schnet", region, buffer.as_deref())?;
                let provider = SchNetInteraction::load(model, *cutoff, elements.clone())
                    .map_err(|e| RunError::Recipe(format!("schnet model '{model}': {e}")))?;
                orchestrator.register_labelled(label_of(index), selection, Box::new(provider));
            },
            #[cfg(not(feature = "ml"))]
            TermSpec::Schnet { .. } => {
                return Err(RunError::MissingFeature {
                    term: term.name().to_string(),
                    feature: "ml",
                })
            },
            TermSpec::Xtb {
                region,
                elements,
                gfn,
                charge,
                multiplicity,
                work_dir,
                timeout_s,
                ..
            } => {
                use gromos_forces::nonbonded::xtb::XtbInteraction;
                let selection = region_of("xtb", region, None)?;
                if let Some(&beyond) = selection.indices().iter().find(|&&i| i >= elements.len()) {
                    return Err(RunError::Recipe(format!(
                        "xtb term: `elements` has {} entries but the region includes atom {} \
                         (one atomic number per global atom index)",
                        elements.len(),
                        beyond + 1
                    )));
                }
                let dir = match work_dir {
                    Some(dir) => std::path::PathBuf::from(dir),
                    None => std::env::temp_dir().join(format!("gromos-rs-xtb-term-{index}")),
                };
                let provider =
                    XtbInteraction::new(dir, *gfn, *charge, *multiplicity, elements.clone())
                        .map_err(|e| RunError::Recipe(format!("xtb term: {e}")))?
                        .with_timeout(std::time::Duration::from_secs(*timeout_s));
                orchestrator.register_labelled(label_of(index), selection, Box::new(provider));
            },
        }
    }
    Ok(Box::new(ProviderOrchestratorAlgorithm::new(
        orchestrator,
        periodicity.clone(),
    )))
}

/// Stage 3: construct the algorithms of a **validated** plan.
///
/// Reads only the plan, the prepared topology/configuration and the boundary conditions —
/// no recipe, no parameter file (drift guard G8).
pub fn instantiate(
    plan: &[AlgorithmSpec],
    topo: &Topology,
    conf: &Configuration,
    periodicity: &Periodicity,
) -> Result<Instantiated, RunError> {
    let n_atoms = topo.num_atoms();
    let mut seq = AlgorithmSequence::new();
    let mut initial_pairs = 0;
    let mut position_restraints = 0;
    let mut distance_restraints = (0, 0);

    for spec in plan {
        match spec {
            AlgorithmSpec::RemoveCom { initial, every } => {
                seq.push(Box::new(RemoveCOMMotion::new(*initial, *every)));
            },
            AlgorithmSpec::Forcefield(ff) => {
                let built = instantiate_forcefield(ff, topo, conf, periodicity)?;
                initial_pairs = built.initial_pairs;
                position_restraints = built.position_restraints;
                distance_restraints = built.distance_restraints;
                seq.push(Box::new(built.forcefield));
            },
            AlgorithmSpec::Orchestrator { terms } => {
                seq.push(instantiate_orchestrator(terms, topo, periodicity)?);
            },
            AlgorithmSpec::LatticeShift => seq.push(Box::new(LatticeShift::new())),
            AlgorithmSpec::LeapFrogVelocity => seq.push(Box::new(LeapFrogVelocity::new())),
            AlgorithmSpec::Thermostat {
                algorithm,
                baths,
                ranges,
            } => {
                let first = baths
                    .first()
                    .ok_or_else(|| RunError::Recipe("thermostat without a bath".to_string()))?;
                if ranges.last().map(|r| r.last_atom + 1) != Some(n_atoms)
                    || ranges
                        .iter()
                        .any(|r| r.com_bath.max(r.ir_bath) >= baths.len())
                {
                    return Err(RunError::Recipe(format!(
                        "thermostat ranges {ranges:?} do not map the {n_atoms} atoms onto {} \
                         bath(s)",
                        baths.len()
                    )));
                }
                match algorithm {
                    ThermostatAlgorithm::Berendsen => {
                        let params = baths
                            .iter()
                            .map(|b| BerendsenThermostatParams {
                                temperature: b.temperature,
                                tau: b.tau,
                                dof: b.dof,
                            })
                            .collect();
                        seq.push(Box::new(BerendsenThermostat::new(params, ranges.clone())))
                    },
                    ThermostatAlgorithm::NoseHoover => {
                        let mut alg = NoseHooverThermostat::new_single_bath(
                            first.temperature,
                            first.tau,
                            first.dof,
                            n_atoms,
                        );
                        alg.ranges.clone_from(ranges);
                        seq.push(Box::new(alg))
                    },
                    ThermostatAlgorithm::NoseHooverChain(n) => {
                        let mut alg = NoseHooverThermostat::new_chain_bath(
                            first.temperature,
                            first.tau,
                            first.dof,
                            n_atoms,
                            (*n).max(2),
                        );
                        alg.ranges.clone_from(ranges);
                        seq.push(Box::new(alg))
                    },
                }
            },
            AlgorithmSpec::LeapFrogPosition => seq.push(Box::new(LeapFrogPosition::new())),
            AlgorithmSpec::Shake {
                tolerance,
                max_iterations,
                solute,
                include_solvent,
                initial_positions,
                initial_velocities,
            } => {
                let mut alg = ShakeAlgorithm::new(ShakeParameters {
                    tolerance: *tolerance,
                    max_iterations: *max_iterations,
                    ntc: ConstraintSelection::ntc_mode_of(solute.ntc()),
                });
                alg.include_solvent = *include_solvent;
                alg.shake_initial_positions = *initial_positions;
                alg.shake_initial_velocities = *initial_velocities;
                seq.push(Box::new(alg));
            },
            AlgorithmSpec::Settle => seq.push(Box::new(SettleAlgorithm::new())),
            AlgorithmSpec::Lincs {
                solute,
                order_solute,
                order_solvent,
                constrain_solute,
                constrain_solvent,
            } => seq.push(Box::new(LincsAlgorithm::new(
                ConstraintSelection::ntc_mode_of(solute.ntc()),
                *order_solute,
                *order_solvent,
                *constrain_solute,
                *constrain_solvent,
            ))),
            AlgorithmSpec::Rottrans { last } => seq.push(Box::new(RottransConstraints::new(*last))),
            AlgorithmSpec::SteepestDescent {
                tolerance,
                step0,
                step_max,
                min_steps,
                force_limit,
            } => seq.push(Box::new(
                SteepestDescentAlgorithm::new()
                    .with_tolerance(*tolerance)
                    .with_step_sizes(*step0, *step_max)
                    .with_min_steps(*min_steps)
                    .with_force_limit(*force_limit),
            )),
            AlgorithmSpec::TemperatureCalculation { n_baths, ranges } => seq.push(Box::new(
                TemperatureCalculation::with_baths(*n_baths, ranges),
            )),
            AlgorithmSpec::PressureCalculation { virial } => {
                seq.push(Box::new(PressureCalculation::new(virial_type(*virial))))
            },
            AlgorithmSpec::Barostat {
                pressure0,
                compressibility,
                tau,
            } => seq.push(Box::new(BerendsenBarostat::new(BerendsenBarostatParams {
                pressure0: *pressure0,
                compressibility: *compressibility,
                tau: *tau,
            }))),
            AlgorithmSpec::EnergyCalculation => seq.push(Box::new(EnergyCalculation::new())),
        }
    }

    Ok(Instantiated {
        sequence: seq,
        initial_pairs,
        position_restraints,
        distance_restraints,
    })
}

fn summarize(
    recipe: &RunRecipe,
    plan: &[AlgorithmSpec],
    topo: &Topology,
    inst: &Instantiated,
) -> BuildSummary {
    let has_solvent = topo.num_atoms() > topo.num_solute_atoms();
    let c = &recipe.constraints;
    let sel = ConstraintSelection::from_codes(
        c.solute.ntc(),
        ConstraintAlgorithm::code(c.solute_algorithm),
        ConstraintAlgorithm::code(c.solvent_algorithm),
        has_solvent,
    );
    let mut summary = BuildSummary {
        is_minimization: recipe.is_minimization(),
        constraints: sel,
        shake: None,
        lincs_orders: (c.lincs_order_solute, c.lincs_order_solvent),
        thermostat: None,
        barostat: None,
        position_restraints: inst.position_restraints,
        distance_restraints: inst.distance_restraints,
        initial_pairs: inst.initial_pairs,
        parallel_nonbonded: false,
        total_dof: total_dof(
            topo,
            &sel,
            ConstraintSelection::ntc_mode_of(c.solute.ntc()),
            recipe.boundary.ndfmin,
        ),
    };
    for spec in plan {
        match spec {
            AlgorithmSpec::Forcefield(ff) => summary.parallel_nonbonded = ff.parallel,
            AlgorithmSpec::Shake {
                tolerance,
                max_iterations,
                solute,
                ..
            } => {
                summary.shake = Some(ShakeParameters {
                    tolerance: *tolerance,
                    max_iterations: *max_iterations,
                    ntc: ConstraintSelection::ntc_mode_of(solute.ntc()),
                })
            },
            AlgorithmSpec::Thermostat {
                algorithm, baths, ..
            } => {
                let label = match algorithm {
                    ThermostatAlgorithm::Berendsen => "Berendsen".to_string(),
                    ThermostatAlgorithm::NoseHoover => "Nose-Hoover".to_string(),
                    ThermostatAlgorithm::NoseHooverChain(n) => format!("Nose-Hoover-Chain({n})"),
                };
                // the summary names the first bath; a multi-bath run says how many there are
                summary.thermostat = Some(ThermostatSummary {
                    label: if baths.len() > 1 {
                        format!("{label} ({} baths)", baths.len())
                    } else {
                        label
                    },
                    temperature: baths.first().map(|b| b.temperature).unwrap_or(0.0),
                    tau: baths.first().map(|b| b.tau).unwrap_or(0.0),
                })
            },
            AlgorithmSpec::Barostat {
                pressure0,
                compressibility,
                tau,
            } => {
                summary.barostat = Some(BarostatSummary {
                    pressure0: *pressure0,
                    tau: *tau,
                    compressibility: *compressibility,
                })
            },
            _ => {},
        }
    }
    summary
}

/// Build a run from a recipe: plan → validate → instantiate, plus the report.
pub fn build_sequence_from_recipe(
    recipe: &RunRecipe,
    prepared: &Prepared,
) -> Result<Built, RunError> {
    let topo = &prepared.topology;
    let plan = build_plan(recipe, topo, prepared.physical_constants.four_pi_eps_i)?;
    validate_plan(&plan)?;
    let periodicity = periodicity_of(prepared);
    let inst = instantiate(&plan, topo, &prepared.configuration, &periodicity)?;
    let summary = summarize(recipe, &plan, topo, &inst);
    Ok(Built {
        sequence: inst.sequence,
        recipe: recipe.clone(),
        plan,
        diagnostics: Diagnostics::default(),
        summary,
    })
}

/// Build a run from a recipe and an already-built (possibly edited) plan: validate →
/// instantiate, plus the report. `Simulation(system, recipe, plan=...)` goes through here.
pub fn build_sequence_from_plan(
    recipe: &RunRecipe,
    plan: Vec<AlgorithmSpec>,
    prepared: &Prepared,
) -> Result<Built, RunError> {
    // the recipe's refusals hold whichever way the plan was made (`build_plan` runs them too)
    crate::plan::check_perturbation(recipe)?;
    crate::plan::check_plan_against_recipe(recipe, &plan)?;
    validate_plan(&plan)?;
    let periodicity = periodicity_of(prepared);
    let inst = instantiate(
        &plan,
        &prepared.topology,
        &prepared.configuration,
        &periodicity,
    )?;
    let summary = summarize(recipe, &plan, &prepared.topology, &inst);
    Ok(Built {
        sequence: inst.sequence,
        recipe: recipe.clone(),
        plan,
        diagnostics: Diagnostics::default(),
        summary,
    })
}

/// Build a run from GROMOS parameters: the `.imd` front-end of [`build_sequence_from_recipe`].
pub fn build_sequence_from_imd(
    imd: &ImdParameters,
    prepared: &Prepared,
    inputs: &RunInputs,
    options: &RunOptions,
) -> Result<Built, RunError> {
    let (mut recipe, diagnostics) = RunRecipe::from_imd_with(imd, &options.passthrough)?;
    recipe.execution.parallel = options.parallel;
    recipe.inputs = inputs.clone();
    let mut built = build_sequence_from_recipe(&recipe, prepared)?;
    built.diagnostics = diagnostics;
    if let Some(partition) = options.pair_partition {
        if let Some(ff) = built
            .sequence
            .find_mut::<gromos_integrators::algorithms::Forcefield>()
        {
            ff.pair_partition = Some(partition);
        }
    }
    Ok(built)
}

/// Initialise the sequence and run step 0 (the initial force evaluation — GROMOS convention).
///
/// Shared so that the binary and the binding start a run identically: `init` is physics
/// (NTISHK shakes the initial positions/velocities), and step 0 primes the double-buffered
/// state. The caller advances `state` afterwards.
pub fn start(
    sequence: &mut AlgorithmSequence,
    topo: &Topology,
    conf: &mut Configuration,
    state: &SimulationState,
) -> Result<(), RunError> {
    sequence
        .init(topo, conf, state)
        .map_err(|e| RunError::Init(e.to_string()))?;
    sequence
        .run_step(topo, conf, state)
        .map_err(|e| RunError::Step {
            step: state.step,
            message: e.to_string(),
        })
}

#[cfg(test)]
mod tests {
    //! The recipe path must be **bit-identical** to the step-1 builder it replaced, on every
    //! reference system (the Rust-side A/B 3.9 step 2 promised).

    use std::path::{Path, PathBuf};

    use gromos_core::algorithm::SimulationState;
    use gromos_io::{
        coordinate::read_coordinates,
        imd::read_imd_file,
        topology::{build_topology, read_topology_file},
    };

    use super::*;
    use crate::{legacy_builder::legacy_build_sequence, prepare_system, Coordinates};
    use gromos_integrators::algorithms::BathRange;

    fn refs() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
    }

    fn toml_str(text: &str, key: &str) -> Option<String> {
        text.lines().find_map(|l| {
            let l = l.trim();
            l.strip_prefix(&format!("{key} = "))
                .map(|v| v.trim().trim_matches('"').to_string())
        })
    }

    struct Loaded {
        imd: ImdParameters,
        inputs: RunInputs,
        prepare: Box<dyn Fn() -> Prepared>,
    }

    fn load(dir: &Path) -> Option<Loaded> {
        let toml = std::fs::read_to_string(dir.join("input.toml")).ok()?;
        let path = |key: &str| toml_str(&toml, key).map(|v| dir.join(v));
        let topo_path = path("topology")?;
        let conf_path = path("configuration")?;
        let imd = read_imd_file(path("parameters")?).ok()?;
        let inputs = RunInputs {
            pttopo: path("pttopo"),
            posresspec: path("posresspec"),
            refpos: path("refpos"),
            distrest: path("distrest"),
        };
        let imd2 = imd.clone();
        let inputs2 = inputs.clone();
        let prepare = Box::new(move || {
            let parsed = read_topology_file(&topo_path).expect("topology");
            let pc = parsed.physical_constants;
            let topo = build_topology(parsed);
            let coords: Coordinates = read_coordinates(&conf_path).expect("conf").into();
            prepare_system(&imd2, topo, pc, coords, &inputs2).expect("prepare_system")
        });
        Some(Loaded {
            imd,
            inputs,
            prepare,
        })
    }

    /// Run `steps` steps and return every number a parity comparison cares about.
    fn run(mut prepared: Prepared, mut seq: AlgorithmSequence, dt: f64, steps: usize) -> Vec<f64> {
        let topo = prepared.topology;
        let conf = &mut prepared.configuration;
        let mut state = SimulationState::new(dt, steps);
        start(&mut seq, &topo, conf, &state).expect("start");
        state.advance();
        for _ in 0..steps {
            seq.run_step(&topo, conf, &state).expect("run_step");
            state.advance();
        }
        let old = conf.old();
        let mut out = vec![
            old.energies.total(),
            old.energies.kinetic_total,
            old.energies.potential_total,
        ];
        out.extend(old.pos.iter().flat_map(|v| [v.x, v.y, v.z]));
        out.extend(old.force.iter().flat_map(|v| [v.x, v.y, v.z]));
        out
    }

    /// `.imd` blocks the test-only legacy builder never modelled.
    const LEGACY_UNMODELLED_BLOCKS: [&str; 2] = ["LAMBDAS", "PRECALCLAM"];

    #[test]
    fn plan_matches_legacy_builder_bit_for_bit() {
        let mut dirs: Vec<PathBuf> = std::fs::read_dir(refs())
            .expect("reference dir")
            .flatten()
            .map(|e| e.path())
            .filter(|p| p.is_dir() && p.join("input.toml").exists())
            .collect();
        dirs.sort();
        assert!(dirs.len() >= 40);
        // The same passthrough policy the binary uses: blocks that carry no physics for the
        // recipe (GAMD and EDS are applied by the binary itself, INNERLOOP selects a kernel).
        let options = RunOptions {
            passthrough: crate::recipe::PassthroughPolicy::allow(["GAMD", "EDS", "INNERLOOP"]),
            ..RunOptions::default()
        };
        let mut compared = 0;
        for dir in dirs {
            let name = dir.file_name().unwrap().to_string_lossy().to_string();
            let Some(loaded) = load(&dir) else { continue };
            // The legacy builder is the pre-3.9 path this test uses as an oracle, and it only ever
            // built a single temperature bath (`imd.temp_bath.first()`). For an input with several
            // baths there is nothing to compare against — the recipe path is checked directly
            // against gromosXX by `test_gromosXX_references.rs` instead.
            if loaded
                .imd
                .temp_bath
                .first()
                .is_some_and(|b| b.temp0.len() > 1)
            {
                continue;
            }
            // Likewise COVALENTFORM: the legacy builder always used the quartic/cosine-harmonic
            // forms, so an input asking for the harmonic ones has no counterpart to compare with.
            if loaded.imd.ntbbh != 0 || loaded.imd.ntbah != 0 || loaded.imd.ntbdn != 0 {
                continue;
            }
            // And it built the reaction field from RCUTL, so an RCRF that differs from it has no
            // legacy counterpart (item 1.8).
            if loaded.imd.rcrf != loaded.imd.rcutl {
                continue;
            }
            // Nor did it model individual λ (LAMBDAS) or extended TI (PRECALCLAM): an input
            // carrying either has no legacy counterpart either.
            let unmodelled: Vec<&str> = LEGACY_UNMODELLED_BLOCKS
                .iter()
                .copied()
                .filter(|b| loaded.imd.raw_blocks.contains_key(*b))
                .collect();
            if !unmodelled.is_empty() {
                println!(
                    "{name}: skipped, the legacy builder does not model {}",
                    unmodelled.join(", ")
                );
                continue;
            }
            let steps = loaded.imd.nstlim.min(10);
            let dt = loaded.imd.dt;

            let p_legacy = (loaded.prepare)();
            let legacy = legacy_build_sequence(&loaded.imd, &p_legacy, &loaded.inputs, &options)
                .unwrap_or_else(|e| panic!("{name}: legacy builder: {e}"));
            let legacy_names: Vec<String> = legacy
                .sequence
                .algorithm_names()
                .iter()
                .map(|s| s.to_string())
                .collect();
            let a = run(p_legacy, legacy.sequence, dt, steps);

            let p_new = (loaded.prepare)();
            let built = build_sequence_from_imd(&loaded.imd, &p_new, &loaded.inputs, &options)
                .unwrap_or_else(|e| panic!("{name}: recipe builder: {e}"));
            let new_names: Vec<String> = built
                .sequence
                .algorithm_names()
                .iter()
                .map(|s| s.to_string())
                .collect();
            assert_eq!(new_names, legacy_names, "{name}: sequences differ");
            let b = run(p_new, built.sequence, dt, steps);

            assert_eq!(a.len(), b.len(), "{name}");
            for (i, (x, y)) in a.iter().zip(&b).enumerate() {
                assert!(
                    x.to_bits() == y.to_bits(),
                    "{name}: value {i} differs after {steps} steps: {x:e} vs {y:e}"
                );
            }
            compared += 1;
        }
        assert!(compared >= 40, "only {compared} systems compared");
    }

    /// The DOFSET ranges are the one atom → bath mapping (gromosXX `Multibath::in_bath`), also
    /// when the bath numbers do not ascend with the atoms.
    #[test]
    fn dofset_ranges_map_atoms_to_baths_in_file_order() {
        let loaded = load(&refs().join("aladip_multibath_swapped")).expect("reference");
        let built = build_sequence_from_imd(
            &loaded.imd,
            &(loaded.prepare)(),
            &loaded.inputs,
            &RunOptions::default(),
        )
        .unwrap();
        let swapped = vec![
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
        let thermostat = built.plan.iter().find_map(|a| match a {
            AlgorithmSpec::Thermostat { ranges, .. } => Some(ranges.clone()),
            _ => None,
        });
        assert_eq!(thermostat, Some(swapped.clone()));
        let kinetic = built.plan.iter().find_map(|a| match a {
            AlgorithmSpec::TemperatureCalculation { n_baths, ranges } => {
                Some((*n_baths, ranges.clone()))
            },
            _ => None,
        });
        assert_eq!(kinetic, Some((2, swapped)));
    }

    #[test]
    fn a_temperature_group_split_over_two_baths_is_refused() {
        // gromosXX multibath.cc:161-181: the solute (atoms 1-12) is one temperature group
        let loaded = load(&refs().join("aladip_multibath_swapped")).expect("reference");
        let mut imd = loaded.imd.clone();
        imd.temp_bath[0].dof_sets = vec![[6, 1, 1], [72, 2, 2]];
        let err = build_sequence_from_imd(
            &imd,
            &(loaded.prepare)(),
            &loaded.inputs,
            &RunOptions::default(),
        )
        .err()
        .expect("refused");
        assert!(
            err.to_string()
                .contains("Temperature group distributed over multiple baths"),
            "{err}"
        );

        imd.temp_bath[0].dof_sets = vec![[42, 1, 1], [12, 2, 2], [72, 1, 1]];
        let err = build_sequence_from_imd(
            &imd,
            &(loaded.prepare)(),
            &loaded.inputs,
            &RunOptions::default(),
        )
        .err()
        .expect("refused");
        assert!(err.to_string().contains("Multibath not sorted"), "{err}");
    }

    #[test]
    fn perturbed_distance_restraints_only_in_a_perturbed_run() {
        // gromosXX ignores PERTDISRESSPEC with NTG = 0 (in_distanceres.cc:520-524)
        let dir = refs().join("aladip_vacuum_fep_distres");
        let mut loaded = load(&dir).expect("reference");
        let built = build_sequence_from_imd(
            &loaded.imd,
            &(loaded.prepare)(),
            &loaded.inputs,
            &RunOptions::default(),
        )
        .unwrap();
        assert_eq!(built.summary.distance_restraints, (6, 6));
        loaded.imd.ntg = 0;
        let imd = loaded.imd.clone();
        let inputs = loaded.inputs.clone();
        let parsed = read_topology_file(dir.join("../shared/aladip.topo")).unwrap();
        let pc = parsed.physical_constants;
        let coords: Coordinates = read_coordinates(dir.join("../aladip_vacuum/aladip_vacuum.conf"))
            .unwrap()
            .into();
        let prepared = prepare_system(&imd, build_topology(parsed), pc, coords, &inputs).unwrap();
        let built =
            build_sequence_from_imd(&imd, &prepared, &inputs, &RunOptions::default()).unwrap();
        assert_eq!(built.summary.distance_restraints, (6, 0));
    }

    #[test]
    fn a_plan_from_another_recipe_still_meets_the_recipes_refusals() {
        // `Simulation(system, recipe, plan=…)`: the plan was made from a recipe without DLAMT
        let dir = refs().join("aladip_vacuum_fep");
        let loaded = load(&dir).expect("reference");
        let prepared = (loaded.prepare)();
        let (mut recipe, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        recipe.inputs = loaded.inputs.clone();
        let plan = build_plan(
            &recipe,
            &prepared.topology,
            prepared.physical_constants.four_pi_eps_i,
        )
        .unwrap();
        assert!(build_sequence_from_plan(&recipe, plan.clone(), &prepared).is_ok());
        let mut slow_growth = recipe.clone();
        slow_growth.perturbation.dlamt = 0.01;
        let err = build_sequence_from_plan(&slow_growth, plan, &prepared)
            .err()
            .expect("refused");
        assert!(err.to_string().contains("slow growth"), "{err}");
    }

    #[test]
    fn solute_temperature_groups_must_reach_the_last_solute_atom() {
        // gromosXX in_topology.cc:278-290
        let dir = refs().join("aladip_multibath_swapped");
        let loaded = load(&dir).expect("reference");
        let parsed = read_topology_file(dir.join("../shared/aladip.topo")).unwrap();
        let pc = parsed.physical_constants;
        let mut topo = build_topology(parsed);
        topo.temperature_groups = vec![(0..6).collect()];
        let coords: Coordinates = read_coordinates(dir.join("../shared/aladip.conf"))
            .unwrap()
            .into();
        let err = prepare_system(&loaded.imd, topo, pc, coords, &loaded.inputs)
            .err()
            .expect("refused");
        assert!(err.to_string().contains("TEMPERATUREGROUPS"), "{err}");
    }
}
