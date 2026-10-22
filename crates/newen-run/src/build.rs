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

use newen_core::{
    algorithm::{AlgorithmSequence, SimulationState},
    configuration::{BoxType, Configuration},
    math::{Periodicity, Rectangular, Triclinic, Vacuum},
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
        BathRange, BerendsenBarostat, BerendsenBarostatParams, BerendsenThermostat,
        BerendsenThermostatParams, EnergyCalculation, Forcefield, LatticeShift, LeapFrogPosition,
        LeapFrogVelocity, LincsAlgorithm, NoseHooverThermostat, NoseHooverThermostatParams,
        PressureCalculation, RemoveCOMMotion, RottransConstraints, SettleAlgorithm, ShakeAlgorithm,
        SteepestDescentAlgorithm, TemperatureCalculation, VirialType,
    },
    constraints::ShakeParameters,
};
use newen_io::{
    distanceres::read_distanceres,
    imd::ImdParameters,
    posres::{build_posres_entries, read_posresspec, read_refpos},
};

use crate::dof::{bath_dofs, rottrans_dof, BathDof};
use crate::plan::{
    build_plan, validate_plan, AlgorithmSpec, DistanceRestraintsPlan, ForcefieldPlan,
    PositionRestraintsPlan,
};
use crate::recipe::{
    ConstraintAlgorithm, Diagnostics, RunRecipe, TermSpec, ThermostatAlgorithm, VirialKind,
};
use crate::{total_dof, ConstraintSelection, Prepared, RunError, RunInputs, RunOptions};
use newen_forces::nonbonded::ewald::{EwaldSum, LatticeSum, LatticeSumParams};
use newen_forces::nonbonded::lattice_sum::ChargeShape;
use newen_forces::nonbonded::p3m::{P3m, P3mParams};
use newen_forces::nonbonded::spme::{Spme, SpmeParams};

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
    /// The temperature baths the run reports (`VOLUMEPRESSURE03`): their total, COM and
    /// internal-rotational degrees of freedom. One bath over the whole system without a
    /// multi-bath thermostat, as gromosXX's default `Multibath`.
    pub baths: Vec<BathDof>,
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

/// The periodicity a prepared system runs under: triclinic for a truncated-octahedron box or a
/// GENBOX box with an angle other than 90°,
/// vacuum for an all-zero box, rectangular otherwise.
pub fn periodicity_of(prepared: &Prepared) -> Periodicity {
    let b = prepared.box_dims;
    if let Some(triclinic_box) = prepared.truncoct_box.or(prepared.triclinic_box) {
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

fn io_err(what: &str, path: &Path) -> impl FnOnce(newen_io::IoError) -> RunError {
    let what = format!("{what} '{}'", path.display());
    move |source| RunError::Io { what, source }
}

fn load_position_restraints(
    plan: &PositionRestraintsPlan,
    startup_positions: &[newen_core::math::Vec3],
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
    forcefield.covalent_form = newen_forces::bonded::CovalentForm {
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
    if let Some(ls) = &ff.lattice_sum {
        let shape = ChargeShape::from_nshape(ls.nshape).ok_or_else(|| {
            RunError::Recipe(format!(
                "NONBONDED block: NSHAPE parameter is {}, should be: >=-1 && <=10",
                ls.nshape
            ))
        })?;
        let nk = ls.nk.map(|n| n as i64);
        let params = LatticeSumParams {
            shape,
            width: ls.ashape,
            nk,
            kcut: ls.kcut,
            na2clc: ls.na2clc,
            tola2: ls.tola2,
            eps_ls: ls.epsls,
            // gromosXX sets `eps0_i = four_pi_eps0_i * 4.0 * Pi` from the topology's
            // PHYSICALCONSTANTS (`in_topology.cc:1405-1409`)
            eps0_inv: ff.four_pi_eps_i * 4.0 * std::f64::consts::PI,
            four_pi_eps_i: ff.four_pi_eps_i,
        };
        forcefield.lattice_sum = Some(match (&ls.p3m, &ls.spme) {
            (None, Some(m)) => {
                let order = m.nspord;
                if !(order % 2 == 0 && (4..=12).contains(&order))
                    || !m.grid.iter().all(|n| n.is_power_of_two())
                {
                    return Err(RunError::Recipe(format!(
                        "NONBONDED SPME settings out of range (grid {:?} powers of two, NSPORD \
                         even in 4..=12): {m:?}",
                        m.grid
                    )));
                }
                LatticeSum::Spme(Box::new(Spme::new(
                    params,
                    SpmeParams {
                        grid: m.grid,
                        order: order as usize,
                    },
                )))
            },
            (None, None) => LatticeSum::Ewald(Box::new(EwaldSum::new(params))),
            (Some(m), _) => {
                let order = |name: &str, v: i32, range: std::ops::RangeInclusive<i32>| {
                    if range.contains(&v) {
                        Ok(v as u8)
                    } else {
                        Err(RunError::Recipe(format!(
                            "NONBONDED block: {name} parameter is {v}, should be in {range:?}"
                        )))
                    }
                };
                let positive = m.nalias > 0 && m.nqeval >= 0 && m.faccur > 0.0;
                if !positive || !m.grid.iter().all(|n| n.is_power_of_two()) {
                    return Err(RunError::Recipe(format!(
                        "NONBONDED P3M settings out of range (grid {:?} powers of two, NALIAS > \
                         0, NQEVAL >= 0, FACCUR > 0): {m:?}",
                        m.grid
                    )));
                }
                LatticeSum::P3m(Box::new(P3m::new(
                    params,
                    P3mParams {
                        grid: m.grid,
                        order: order("NASORD", m.nasord, 1..=5)?,
                        nfdord: order("NFDORD", m.nfdord, 0..=5)?,
                        nalias: i64::from(m.nalias),
                        nqeval: m.nqeval as usize,
                        faccur: m.faccur,
                    },
                )))
            },
        });
    }

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
) -> Result<Box<dyn newen_core::algorithm::Algorithm>, RunError> {
    use newen_core::selection::AtomSelection;
    use newen_forces::orchestrator::ProviderOrchestrator;
    use newen_forces::orchestrator_algorithm::ProviderOrchestratorAlgorithm;
    use newen_forces::zones::{Zone, ZonePartition};

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
                use newen_forces::nonbonded::schnet::SchNetInteraction;
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
                use newen_forces::nonbonded::xtb::XtbInteraction;
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
                    None => std::env::temp_dir().join(format!("newen-xtb-term-{index}")),
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
    // the constraint algorithms add their virial only for the atomic virial
    let mut atomic_virial = false;
    // the force field computes end states exactly when the topology has EDS atoms; only the
    // EDS algorithm uses them
    let eds_states = plan.iter().find_map(|a| match a {
        AlgorithmSpec::Eds(e) => Some(e.num_states),
        _ => None,
    });
    match (eds_states, topo.eds_solute.is_empty()) {
        (Some(n), false) if n != topo.eds_solute.num_states => {
            return Err(RunError::InvalidPlan(format!(
                "the EDS step mixes {n} states, the topology's MPERTATOM block has {}",
                topo.eds_solute.num_states
            )))
        },
        (Some(_), true) => {
            return Err(RunError::InvalidPlan(
                "an EDS step without EDS atoms (no MPERTATOM in the topology)".to_string(),
            ))
        },
        (None, false) => {
            return Err(RunError::InvalidPlan(
                "the topology has EDS atoms (MPERTATOM) but the plan has no EDS step".to_string(),
            ))
        },
        _ => {},
    }

    for spec in plan {
        match spec {
            AlgorithmSpec::RemoveCom { initial, every } => {
                seq.push(Box::new(RemoveCOMMotion::new(*initial, *every)));
            },
            AlgorithmSpec::Forcefield(ff) => {
                atomic_virial = ff.virial == VirialKind::Atomic;
                let built = instantiate_forcefield(ff, topo, conf, periodicity)?;
                initial_pairs = built.initial_pairs;
                position_restraints = built.position_restraints;
                distance_restraints = built.distance_restraints;
                seq.push(Box::new(built.forcefield));
            },
            AlgorithmSpec::Orchestrator { terms } => {
                seq.push(instantiate_orchestrator(terms, topo, periodicity)?);
            },
            AlgorithmSpec::Eds(e) => {
                use newen_integrators::eds::{AedsForm, AedsParams, Eds, EdsForm};
                if let Some(a) = &e.aeds {
                    let form = AedsForm::from_gromos(a.form).ok_or_else(|| {
                        RunError::InvalidPlan(format!("AEDS FORM {}: must be 1-6", a.form))
                    })?;
                    let params = AedsParams {
                        form,
                        bmax_type: a.bmax_type,
                        bmax: a.bmax,
                        asteps: a.asteps,
                        bsteps: a.bsteps,
                    };
                    seq.push(Box::new(Eds::accelerated(
                        params,
                        a.emax,
                        a.emin,
                        e.eir.clone(),
                        a.fullemin,
                        a.restart,
                        e.beta,
                    )));
                    continue;
                }
                let form = match e.form {
                    1 => EdsForm::SingleS {
                        s: e.s.first().copied().unwrap_or(1.0),
                    },
                    2 => EdsForm::MultiS { s: e.s.clone() },
                    _ => EdsForm::PairS {
                        pairs: e.pairs.clone(),
                    },
                };
                seq.push(Box::new(Eds::new(form, e.eir.clone(), e.beta)));
            },
            AlgorithmSpec::LatticeShift => seq.push(Box::new(LatticeShift::new())),
            AlgorithmSpec::LeapFrogVelocity => seq.push(Box::new(LeapFrogVelocity::new())),
            AlgorithmSpec::Thermostat {
                algorithm,
                baths,
                ranges,
                read_variables,
            } => {
                if baths.is_empty() {
                    return Err(RunError::Recipe("thermostat without a bath".to_string()));
                }
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
                    ThermostatAlgorithm::NoseHoover | ThermostatAlgorithm::NoseHooverChain(_) => {
                        // one chain length for every bath, as gromosXX's MULTIBATH ALGORITHM
                        let nhc = match algorithm {
                            ThermostatAlgorithm::NoseHooverChain(n) => (*n).max(2),
                            _ => 1,
                        };
                        let params = baths
                            .iter()
                            .map(|b| NoseHooverThermostatParams {
                                temperature: b.temperature,
                                tau: b.tau,
                                dof: b.dof,
                                nhc,
                            })
                            .collect();
                        let mut alg = NoseHooverThermostat::new(params, ranges.clone());
                        alg.read_variables = *read_variables;
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
                alg.atomic_virial = atomic_virial;
                seq.push(Box::new(alg));
            },
            AlgorithmSpec::Settle => {
                let mut alg = SettleAlgorithm::new();
                alg.atomic_virial = atomic_virial;
                seq.push(Box::new(alg))
            },
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
    let ntc_mode = ConstraintSelection::ntc_mode_of(c.solute.ntc());
    let rottrans = c.rottrans_last.is_some();
    let (n_baths, ranges) = plan
        .iter()
        .find_map(|spec| match spec {
            AlgorithmSpec::TemperatureCalculation { n_baths, ranges } if !ranges.is_empty() => {
                Some((*n_baths, ranges.clone()))
            },
            _ => None,
        })
        .unwrap_or_else(|| {
            let whole = BathRange {
                last_atom: topo.num_atoms().saturating_sub(1),
                com_bath: 0,
                ir_bath: 0,
            };
            (1, vec![whole])
        });
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
        total_dof: total_dof(topo, &sel, ntc_mode, recipe.boundary.ndfmin)
            - rottrans_dof(topo, rottrans),
        baths: bath_dofs(
            topo,
            &sel,
            ntc_mode,
            recipe.boundary.ndfmin,
            n_baths,
            &ranges,
            rottrans,
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
    let plan = build_plan(recipe, topo, &prepared.physical_constants)?;
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
    crate::plan::check_start_time(recipe)?;
    crate::plan::check_com_motion(recipe)?;
    crate::plan::check_electrostatics(recipe)?;
    crate::plan::check_pairlist(recipe)?;
    crate::plan::check_outputs(recipe)?;
    crate::plan::check_dihedral_phases(recipe, &prepared.topology)?;
    crate::plan::check_perturbation(recipe)?;
    crate::plan::check_eds(recipe)?;
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
    for (key, value) in &options.recipe_overrides {
        recipe.set(key, value)?;
    }
    let mut built = build_sequence_from_recipe(&recipe, prepared)?;
    built.diagnostics = diagnostics;
    if let Some(partition) = options.pair_partition {
        if let Some(ff) = built
            .sequence
            .find_mut::<newen_integrators::algorithms::Forcefield>()
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

/// Slow growth: after every step ([`start`]'s included), λ += DLAMT·Δt and the individual λ
/// and perturbed masses follow ([`Topology::update_for_lambda`]). A no-op unless the run is
/// perturbed (NTG = 1) with DLAMT ≠ 0.
///
/// gromosXX `Slow_Growth::apply` (`slow_growth.cc:43-53`), pushed after the barostat and before
/// the energy totals (`create_md_sequence.cc:279-286`): the step's forces, integration,
/// constraints and kinetic energy used the old λ, the next step uses the new one. λ is not
/// bounded: it passes 1 when RLAM + DLAMT·NSTLIM·Δt > 1, as in gromosXX. The step's outputs
/// carry the old λ (`.trg`, gromosXX `topo.old_lambda()`), except the `.tre` mass, summed over
/// the masses already at the new λ (`out_configuration.cc:1786`); the final configuration's
/// `PERTDATA` holds the new λ.
///
/// Returns the new λ on the step that takes it out of [0, 1] (from inside), for the caller to
/// warn ([`lambda_warning`]); gromosXX neither clamps nor warns.
pub fn advance_lambda(recipe: &RunRecipe, topo: &mut Topology) -> Option<f64> {
    let p = &recipe.perturbation;
    if p.enabled == 0 || p.dlamt == 0.0 {
        return None;
    }
    let before = topo.lambda;
    topo.lambda += p.dlamt * recipe.control.dt;
    topo.update_for_lambda();
    (lambda_in_range(before) && !lambda_in_range(topo.lambda)).then_some(topo.lambda)
}

/// λ in [0, 1], up to 10⁻⁹: a sum of DLAMT·Δt steps meant to end on 1 may overshoot it by
/// rounding.
pub fn lambda_in_range(lambda: f64) -> bool {
    (-1e-9..=1.0 + 1e-9).contains(&lambda)
}

/// The warning for a λ outside [0, 1]: reached by slow growth after step `after_step`
/// ([`advance_lambda`]), or, with `None`, the λ a run starts at
/// ([`crate::PrepareNote::LambdaOutOfRange`]).
pub fn lambda_warning(after_step: Option<usize>, lambda: f64) -> String {
    let what = match after_step {
        Some(step) => format!("slow growth takes lambda out of [0, 1] after step {step}"),
        None => "the run starts outside [0, 1]".to_string(),
    };
    format!(
        "{what}: lambda = {lambda:.6}; every perturbed term and mass is extrapolated past the \
         end state"
    )
}

#[cfg(test)]
mod tests {
    //! The recipe path must be **bit-identical** to the step-1 builder it replaced, on every
    //! reference system (the Rust-side A/B 3.9 step 2 promised).

    use std::path::{Path, PathBuf};

    use newen_core::algorithm::SimulationState;
    use newen_io::{
        coordinate::read_coordinates,
        imd::read_imd_file,
        topology::{build_topology, read_topology_file},
    };

    use super::*;
    use crate::{legacy_builder::legacy_build_sequence, prepare_system, Coordinates};
    use newen_integrators::algorithms::BathRange;

    fn refs() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../newen-md/tests/gromosXX_references")
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
    const LEGACY_UNMODELLED_BLOCKS: [&str; 4] = ["LAMBDAS", "PRECALCLAM", "EDS", "AEDS"];

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
        // recipe (GAMD is applied by the binary itself, INNERLOOP selects a kernel).
        let options = RunOptions {
            passthrough: crate::recipe::PassthroughPolicy::allow(["GAMD", "INNERLOOP"]),
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
            // Nor did it model individual λ (LAMBDAS), extended TI (PRECALCLAM) or (A)EDS: an input
            // carrying any of them has no legacy counterpart either.
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
            if matches!(loaded.imd.nlrele, 2 | 3 | 5) {
                println!("{name}: skipped, the legacy builder does not model lattice sums");
                continue;
            }
            let steps = loaded.imd.nstlim.min(10);
            let dt = loaded.imd.dt;

            let p_legacy = (loaded.prepare)();
            // Nor did it restore Nosé–Hoover variables (a continuation's NHCVARIABLES, D41).
            if p_legacy.configuration.nose_hoover.is_some() {
                continue;
            }
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
    fn slow_growth_moves_lambda_and_the_masses_without_a_bound() {
        // gromosXX slow_growth.cc:49-52; atom 1 (CB) is 15.035 in A, 12.035 in B
        let dir = refs().join("aladip_vacuum_fep");
        let loaded = load(&dir).expect("reference");
        let mut topo = (loaded.prepare)().topology;
        let (mut recipe, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        assert_eq!(
            (recipe.perturbation.lambda, recipe.control.dt),
            (0.125, 0.002)
        );
        advance_lambda(&recipe, &mut topo);
        assert_eq!(topo.lambda, 0.125, "DLAMT = 0: nothing moves");
        recipe.perturbation.dlamt = 2.0;
        let mut expected = 0.125;
        for _ in 0..3 {
            advance_lambda(&recipe, &mut topo);
            expected += 2.0 * 0.002;
        }
        assert_eq!(topo.lambda, expected);
        assert_eq!(topo.mass[0], (1.0 - expected) * 15.035 + expected * 12.035);
        topo.lambda = 0.995;
        assert_eq!(
            advance_lambda(&recipe, &mut topo),
            None,
            "λ = 0.999: inside"
        );
        assert_eq!(
            advance_lambda(&recipe, &mut topo),
            Some(0.995 + 2.0 * 0.002 + 2.0 * 0.002),
            "the step that leaves [0, 1] is reported"
        );
        assert_eq!(
            topo.lambda,
            0.995 + 2.0 * 0.002 + 2.0 * 0.002,
            "no clamp at 1"
        );
        assert_eq!(advance_lambda(&recipe, &mut topo), None, "reported once");
        assert!(lambda_warning(Some(2), 1.007).contains("after step 2: lambda = 1.007000"));
        recipe.perturbation.enabled = 0;
        let before = topo.lambda;
        advance_lambda(&recipe, &mut topo);
        assert_eq!(topo.lambda, before, "NTG = 0: inert");
    }

    #[test]
    fn a_plan_from_another_recipe_still_meets_the_recipes_refusals() {
        // `Simulation(system, recipe, plan=…)`: the plan was made from a recipe without NSCALE
        let dir = refs().join("aladip_vacuum_fep");
        let loaded = load(&dir).expect("reference");
        let prepared = (loaded.prepare)();
        let (mut recipe, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        recipe.inputs = loaded.inputs.clone();
        let plan = build_plan(&recipe, &prepared.topology, &prepared.physical_constants).unwrap();
        assert!(build_sequence_from_plan(&recipe, plan.clone(), &prepared).is_ok());
        let mut scaled = recipe.clone();
        scaled.perturbation.nscale = 1;
        let err = build_sequence_from_plan(&scaled, plan, &prepared)
            .err()
            .expect("refused");
        assert!(err.to_string().contains("NSCALE"), "{err}");
    }

    #[test]
    fn lattice_sums_run_with_an_atomic_pairlist() {
        let dir = refs().join("nacl_water_box");
        let loaded = load(&dir).expect("reference");
        let prepared = (loaded.prepare)();
        let (mut recipe, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        recipe.inputs = loaded.inputs.clone();
        assert_eq!(recipe.perturbation.enabled, 0);
        let rf_plan =
            build_plan(&recipe, &prepared.topology, &prepared.physical_constants).unwrap();
        let refused = |r: &RunRecipe, message: &str| {
            let err = build_plan(r, &prepared.topology, &prepared.physical_constants)
                .expect_err("refused");
            assert!(err.to_string().contains(message), "{err}");
            let err = build_sequence_from_plan(r, rf_plan.clone(), &prepared)
                .err()
                .expect("refused");
            assert!(err.to_string().contains(message), "{err}");
        };
        // P3M needs the atomic cutoff too
        let mut p3m = recipe.clone();
        p3m.forcefield.electrostatics.electrostatics = 3;
        refused(&p3m, "chargegroup based cutoff does not work with P3M");
        // Ewald needs the atomic cutoff (gromosXX FeatureChecker)
        let mut ewald = recipe.clone();
        ewald.forcefield.electrostatics.electrostatics = 2;
        refused(&ewald, "chargegroup based cutoff does not work with Ewald");
        ewald.forcefield.pairlist.kind = crate::recipe::PairlistKind::Atomic;
        refused(&ewald, "multiple energy groups does not work with Ewald");

        // the Ewald reference builds; a plan made from its reaction-field variant does not run
        // the reaction field in its place
        let loaded = load(&refs().join("nacl_pair_box_ewald")).expect("reference");
        let prepared = (loaded.prepare)();
        let (mut ewald, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        ewald.inputs = loaded.inputs.clone();
        let plan = build_plan(&ewald, &prepared.topology, &prepared.physical_constants).unwrap();
        let mut rf = ewald.clone();
        rf.forcefield.electrostatics.electrostatics = 1;
        let rf_plan = build_plan(&rf, &prepared.topology, &prepared.physical_constants).unwrap();
        let err = build_sequence_from_plan(&ewald, rf_plan, &prepared)
            .err()
            .expect("refused");
        assert!(err.to_string().contains("lattice-sum"), "{err}");
        let err = build_sequence_from_plan(&rf, plan.clone(), &prepared)
            .err()
            .expect("refused");
        assert!(err.to_string().contains("lattice-sum"), "{err}");
        assert!(build_sequence_from_plan(&ewald, plan, &prepared).is_ok());
        // nor one lattice-sum method in place of the other
        let mut p3m = ewald.clone();
        p3m.forcefield.electrostatics.electrostatics = 3;
        let p3m_plan = build_plan(&p3m, &prepared.topology, &prepared.physical_constants).unwrap();
        let err = build_sequence_from_plan(&ewald, p3m_plan.clone(), &prepared)
            .err()
            .expect("refused");
        assert!(err.to_string().contains("has P3M lattice-sum"), "{err}");
        assert!(build_sequence_from_plan(&p3m, p3m_plan, &prepared).is_ok());
        for nlrele in [-1, 0, 1] {
            let mut rf = recipe.clone();
            rf.forcefield.electrostatics.electrostatics = nlrele;
            assert!(
                crate::plan::check_electrostatics(&rf).is_ok(),
                "NLRELE = {nlrele}"
            );
        }
    }

    /// Both build paths refuse what `refuse` makes of `recipe`, with `message`.
    fn assert_refused_both_ways(
        recipe: &RunRecipe,
        prepare: &dyn Fn() -> Prepared,
        refuse: &dyn Fn(&mut RunRecipe, &mut Prepared),
        message: &str,
    ) {
        let prepared = prepare();
        let plan = build_plan(recipe, &prepared.topology, &prepared.physical_constants).unwrap();
        let (mut r, mut p) = (recipe.clone(), prepare());
        refuse(&mut r, &mut p);
        let err = build_plan(&r, &p.topology, &p.physical_constants).expect_err("refused");
        assert!(err.to_string().contains(message), "{err}");
        let err = build_sequence_from_plan(&r, plan, &p)
            .err()
            .expect("refused");
        assert!(err.to_string().contains(message), "{err}");
    }

    #[test]
    fn grid_pairlist_is_refused_in_a_truncated_octahedron() {
        // gromosXX FeatureChecker: `pbc_t` is never unlocked against `pairlist_grid`
        let dir = refs().join("aladip_trunc_oct");
        let loaded = load(&dir).expect("reference");
        let (mut recipe, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        recipe.inputs = loaded.inputs.clone();
        assert_eq!(
            recipe.boundary.kind,
            crate::recipe::BoundaryKind::TruncatedOctahedron
        );
        assert_refused_both_ways(
            &recipe,
            &*loaded.prepare,
            &|r, _| r.forcefield.pairlist.algorithm = crate::recipe::PairlistAlgorithm::Grid,
            "truncated octahedral PBC does not work with grid-based pairlist algorithm",
        );
        let mut cell = recipe.clone();
        cell.forcefield.pairlist.algorithm = crate::recipe::PairlistAlgorithm::GridCell;
        assert!(
            crate::plan::check_pairlist(&cell).is_ok(),
            "grid_cell is unlocked"
        );
        let mut rect = recipe.clone();
        rect.forcefield.pairlist.algorithm = crate::recipe::PairlistAlgorithm::Grid;
        rect.boundary.kind = crate::recipe::BoundaryKind::Rectangular;
        assert!(crate::plan::check_pairlist(&rect).is_ok());
    }

    #[test]
    fn ntbdn_1_refuses_a_dihedral_phase_other_than_0_or_180() {
        // gromosXX in_topology.cc:960-962
        let dir = refs().join("butane_vacuum");
        let loaded = load(&dir).expect("reference");
        let (mut recipe, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        recipe.inputs = loaded.inputs.clone();
        recipe.forcefield.covalent_form.limited_phase_shifts = true;
        let phase = |p: &mut Prepared, deg: f64| {
            let d = &mut p.topology.dihedral_parameters[0];
            d.pd = deg * std::f64::consts::PI / 180.0;
            d.cospd = d.pd.cos();
        };
        assert_refused_both_ways(
            &recipe,
            &*loaded.prepare,
            &|_, p| phase(p, 90.0),
            "dihedral neither 0 nor 180, this is not allowed if NTBDN=1",
        );
        let mut p = (loaded.prepare)();
        for deg in [0.0, 180.0] {
            phase(&mut p, deg);
            assert!(
                crate::plan::check_dihedral_phases(&recipe, &p.topology).is_ok(),
                "{deg}°"
            );
        }
        phase(&mut p, 90.0);
        let mut arbitrary = recipe.clone();
        arbitrary.forcefield.covalent_form.limited_phase_shifts = false;
        assert!(crate::plan::check_dihedral_phases(&arbitrary, &p.topology).is_ok());
        let mut off = recipe.clone();
        off.forcefield.dihedrals = false;
        assert!(
            crate::plan::check_dihedral_phases(&off, &p.topology).is_ok(),
            "gromosXX checks only with the dihedral term on"
        );
    }

    #[test]
    fn step_t_outside_its_range_is_refused() {
        // gromosXX in_parameter.cc:298: T >= 0, or -1 (time from the TIMESTEP block)
        let dir = refs().join("aladip_vacuum");
        let loaded = load(&dir).expect("reference");
        let (mut recipe, _) = RunRecipe::from_imd(&loaded.imd).unwrap();
        recipe.inputs = loaded.inputs.clone();
        assert_refused_both_ways(
            &recipe,
            &*loaded.prepare,
            &|r, _| r.control.t0 = -1.0,
            "STEP T = -1 (the start time from the configuration's TIMESTEP block) is not \
             implemented",
        );
        assert_refused_both_ways(
            &recipe,
            &*loaded.prepare,
            &|r, _| r.control.t0 = -0.5,
            "STEP T = -0.5: must be >= 0",
        );
        let mut later = recipe.clone();
        later.control.t0 = 10.0;
        assert!(crate::plan::check_start_time(&later).is_ok());
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
