//! `AlgorithmSpec` — the MD step as data (stage 1 of the builder, 3.9).
//!
//! `build_plan` turns a [`RunRecipe`] into the ordered list of algorithms the GROMOS step runs,
//! with **every value resolved** (degrees of freedom, virial kind, `four_pi_eps_i`, restraint
//! file paths, the parallel decision, …). `instantiate` (in `build.rs`) reads *only* this list —
//! never the recipe — so an edited plan is exactly what runs (drift guard G8).
//!
//! `validate_plan` enforces the GROMOS ordering invariants, declared per kind in
//! [`AlgorithmSpec::rules`] (drift guard G9); every plan goes through it, edited or not.
//!
//! Registries: [`AlgorithmSpec::name`] and [`TermSpec`]'s are exhaustive `match`es, so a new
//! variant without a name is a compile error; [`AlgorithmSpec::examples`] must list one value
//! per variant (`registry_is_complete` test), which is what the Python side introspects.

use std::path::PathBuf;

use newen_core::{units::PhysicalConstants, Topology};
use serde::{Deserialize, Serialize};

use crate::dof::{bath_dof, rottrans_dof};
use crate::recipe::{
    BoundaryKind, ConstraintAlgorithm, Coupling, SoluteConstraints, ThermostatAlgorithm,
};
use crate::recipe::{RunRecipe, TermSpec, Thermostat, VirialKind};
use crate::{total_dof, ConstraintSelection, RunError};
use newen_integrators::algorithms::{in_bath, BathRange};
use newen_integrators::constraints::NtcMode;

/// Position restraints as the plan carries them: the input files, the force constant, and
/// whether the reference positions are the start-up configuration (`reference: None`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PositionRestraintsPlan {
    pub spec: PathBuf,
    /// `None` = the start-up positions (NTPORB=0), resolved at instantiation from the
    /// configuration, so the plan stays the same run after run.
    pub reference: Option<PathBuf>,
    pub k: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DistanceRestraintsPlan {
    pub file: PathBuf,
    /// |NTDIR|
    pub mode: i32,
    pub k: f64,
    pub r_linear: f64,
    /// `Some(tau)` for time-averaged restraints (NTDIR < 0); needs `dt`.
    pub time_averaging: Option<f64>,
    pub dt: f64,
    /// A perturbed run (NTG = 1): only then are the file's PERTDISRESSPEC rows restraints.
    #[serde(default)]
    pub perturbed: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PerturbationPlan {
    pub alpha_lj: f64,
    pub alpha_c: f64,
    /// A `LAMBDAS` block sets individual λ (NTIL on). λ itself, NLAM and the individual-λ table
    /// travel with the prepared topology (`Topology::lambda`, `lambda_exp`,
    /// `individual_lambdas`), as in gromosXX.
    #[serde(default)]
    pub individual_lambdas: bool,
    /// `Perturbation::atom_pair_kernel`.
    #[serde(default, skip_serializing_if = "is_intended")]
    pub atom_pair_kernel: crate::recipe::AtomPairKernel,
}

fn is_intended(k: &crate::recipe::AtomPairKernel) -> bool {
    *k == crate::recipe::AtomPairKernel::Intended
}

/// The EDS reference Hamiltonian, resolved (gromosXX `algorithm::EDS`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdsPlan {
    /// FORM: 1 single s, 2 one s per state pair, 3 one s per listed pair.
    pub form: i32,
    pub num_states: usize,
    /// FORM 1: one value; FORM 2: N(N−1)/2 in pair order (0,1), (0,2) … (1,2) …
    pub s: Vec<f64>,
    /// FORM 3: the listed pairs, states 0-based.
    pub pairs: Vec<(usize, usize, f64)>,
    /// E_i^R per state (kJ/mol).
    pub eir: Vec<f64>,
    /// 1/(k_B T): T of the first MULTIBATH bath, k_B from the topology's PHYSICALCONSTANTS
    /// (gromosXX `eds.cc:385-391`, `in_topology.cc:1407`).
    pub beta: f64,
    /// Accelerated EDS (an `AEDS` block): FORM 1, s = 1, and these settings; `eir` holds the
    /// starting offsets.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aeds: Option<AedsPlan>,
}

/// The settings of accelerated EDS, resolved (gromosXX `read_AEDS`, `eds.cc:57-378`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AedsPlan {
    /// FORM 1-6.
    pub form: i32,
    /// starting E_max and E_min (kJ/mol)
    pub emax: f64,
    pub emin: f64,
    /// NTIAEDSS = 0 in FORMS 2-6: continue the configuration's AEDSSEARCH state.
    pub restart: bool,
    /// RESTREMIN = 0: E_min unrestricted from the start (gromosXX `fullemin`).
    pub fullemin: bool,
    pub bmax_type: u32,
    pub bmax: f64,
    pub asteps: u32,
    pub bsteps: u32,
    /// `Aeds::energy_columns`.
    #[serde(default, skip_serializing_if = "is_default")]
    pub energy_columns: crate::recipe::AedsColumns,
    /// `Aeds::prevalence`.
    #[serde(default, skip_serializing_if = "is_default")]
    pub prevalence: crate::recipe::AedsPrevalence,
    /// `Aeds::target_state`.
    #[serde(default, skip_serializing_if = "is_default")]
    pub target_state: crate::recipe::AedsTarget,
}

fn is_default<T: Default + PartialEq>(v: &T) -> bool {
    *v == T::default()
}

/// Everything `Forcefield` needs, resolved.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ForcefieldPlan {
    pub cutoff_short: f64,
    pub cutoff_long: f64,
    pub pairlist_every: usize,
    pub grid_size: f64,
    /// PAIRLIST ALGORITHM code (0 standard, 1 grid, 2 grid_cell) — `PairlistAlgorithm::from_imd`
    pub pairlist_algorithm: i32,
    /// PAIRLIST TYPE code (0 chargegroup, 1 atomic)
    pub pairlist_type: i32,
    /// RCRF: the reaction-field radius the RF constants are built from (gromosXX
    /// `nonbonded.rf_cutoff`), which need not be the long-range cutoff.
    pub rf_cutoff: f64,
    pub epsilon_rf: f64,
    pub kappa: f64,
    /// NSLFEXCL: reaction-field contribution of the excluded pairs and the self term
    /// (gromosXX `nonbonded.rf_excluded`; 0 is the GROMOS96 behaviour, 1 the current one).
    pub rf_excluded: bool,
    pub four_pi_eps_i: f64,
    pub bonds: bool,
    /// COVALENTFORM: NTBBH / NTBAH / NTBDN, as `(harmonic_bonds, harmonic_angles,
    /// limited_phase_shifts)` — the GROMOS defaults are all false.
    #[serde(default)]
    pub covalent_form: (bool, bool, bool),
    pub angles: bool,
    pub impropers: bool,
    pub dihedrals: bool,
    pub parallel: bool,
    pub atoms_per_solvent: usize,
    pub virial: VirialKind,
    pub position_restraints: Option<PositionRestraintsPlan>,
    pub distance_restraints: Option<DistanceRestraintsPlan>,
    pub perturbation: Option<PerturbationPlan>,
}

/// One algorithm of the MD step, fully resolved. One spec = one `Algorithm`.
///
/// `Forcefield` is much larger than the unit variants; that is the nature of the data (it is
/// what the force field needs) and a plan has one of them, so it is not boxed.
/// One temperature bath of a `MULTIBATH` block, resolved onto the topology.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThermostatBath {
    pub temperature: f64,
    /// TAU (ps); < 0 = this bath is not coupled
    pub tau: f64,
    /// Kinetic degrees of freedom of this bath (resolved from DOFSET and the constraints).
    pub dof: f64,
}

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum AlgorithmSpec {
    RemoveCom {
        initial: i32,
        every: i32,
    },
    Forcefield(ForcefieldPlan),
    Orchestrator {
        terms: Vec<TermSpec>,
    },
    /// Mix the EDS end states into the forces, right after the force field.
    Eds(EdsPlan),
    /// Put every charge group back into the box (gromosXX `Lattice_Shift_Tracker`), once per step
    /// before the force field. Only for periodic boundaries.
    LatticeShift,
    LeapFrogVelocity,
    Thermostat {
        algorithm: ThermostatAlgorithm,
        /// One entry per MULTIBATH bath, in file order.
        baths: Vec<ThermostatBath>,
        /// The DOFSET ranges that assign the atoms to `baths` (gromosXX `bath_index`), in file
        /// order; the last one ends at the last atom.
        ranges: Vec<BathRange>,
    },
    LeapFrogPosition,
    Shake {
        tolerance: f64,
        max_iterations: usize,
        solute: SoluteConstraints,
        include_solvent: bool,
        initial_positions: bool,
        initial_velocities: bool,
    },
    Settle,
    Lincs {
        solute: SoluteConstraints,
        order_solute: usize,
        order_solvent: usize,
        constrain_solute: bool,
        constrain_solvent: bool,
    },
    /// ROTTRANS — constrain the six rigid-body degrees of freedom of the first `last` atoms.
    /// Runs immediately after the distance constraints, as gromosXX creates it
    /// (`create_constraints.cc:294`).
    Rottrans {
        last: usize,
    },
    SteepestDescent {
        tolerance: f64,
        step0: f64,
        step_max: f64,
        min_steps: usize,
        force_limit: f64,
    },
    TemperatureCalculation {
        /// Number of temperature baths the kinetic energy is booked to.
        #[serde(default)]
        n_baths: usize,
        /// The DOFSET ranges that assign the atoms to the baths, in file order; empty = one bath
        /// over the whole system. The kinetic energy is accumulated per bath, as gromosXX does.
        #[serde(default)]
        ranges: Vec<BathRange>,
    },
    PressureCalculation {
        virial: VirialKind,
    },
    Barostat {
        pressure0: f64,
        compressibility: f64,
        tau: f64,
    },
    EnergyCalculation,
}

/// Ordering rules of one kind, consulted by `validate_plan`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KindRules {
    /// At most one of this kind.
    pub unique: bool,
    /// Must appear exactly once.
    pub required: bool,
    /// Must be the last algorithm.
    pub last: bool,
    /// Must be the first algorithm, if present.
    pub first: bool,
    /// Must come immediately after this kind.
    pub must_follow: Option<&'static str>,
    /// Must appear somewhere after every one of these kinds.
    pub after: &'static [&'static str],
    /// Must appear somewhere before every one of these kinds.
    pub before: &'static [&'static str],
    /// Cannot coexist with these kinds.
    pub excludes: &'static [&'static str],
    /// Requires these kinds to be present.
    pub requires: &'static [&'static str],
}

const NONE: KindRules = KindRules {
    unique: false,
    required: false,
    last: false,
    first: false,
    must_follow: None,
    after: &[],
    before: &[],
    excludes: &[],
    requires: &[],
};

const DYNAMICS_ONLY: &[&str] = &[
    "leap_frog_velocity",
    "leap_frog_position",
    "thermostat",
    "temperature_calculation",
    "pressure_calculation",
    "barostat",
    "remove_com",
];

impl AlgorithmSpec {
    /// Every kind, in canonical order. `examples()` and the tests keep this list honest.
    pub const KINDS: &'static [&'static str] = &[
        "remove_com",
        "forcefield",
        "orchestrator",
        "eds",
        "leap_frog_velocity",
        "thermostat",
        "leap_frog_position",
        "shake",
        "settle",
        "lincs",
        "rottrans",
        "steepest_descent",
        "temperature_calculation",
        "pressure_calculation",
        "barostat",
        "energy_calculation",
    ];

    /// The serde tag of this spec — exhaustive, so a new variant cannot be nameless.
    pub fn name(&self) -> &'static str {
        match self {
            AlgorithmSpec::RemoveCom { .. } => "remove_com",
            AlgorithmSpec::Forcefield(_) => "forcefield",
            AlgorithmSpec::Orchestrator { .. } => "orchestrator",
            AlgorithmSpec::Eds(_) => "eds",
            AlgorithmSpec::LatticeShift => "lattice_shift",
            AlgorithmSpec::LeapFrogVelocity => "leap_frog_velocity",
            AlgorithmSpec::Thermostat { .. } => "thermostat",
            AlgorithmSpec::LeapFrogPosition => "leap_frog_position",
            AlgorithmSpec::Shake { .. } => "shake",
            AlgorithmSpec::Settle => "settle",
            AlgorithmSpec::Lincs { .. } => "lincs",
            AlgorithmSpec::Rottrans { .. } => "rottrans",
            AlgorithmSpec::SteepestDescent { .. } => "steepest_descent",
            AlgorithmSpec::TemperatureCalculation { .. } => "temperature_calculation",
            AlgorithmSpec::PressureCalculation { .. } => "pressure_calculation",
            AlgorithmSpec::Barostat { .. } => "barostat",
            AlgorithmSpec::EnergyCalculation => "energy_calculation",
        }
    }

    /// One example value per variant — the registry Python introspects (`newen.algorithms()`)
    /// and the test that proves `KINDS` is complete.
    pub fn examples() -> Vec<AlgorithmSpec> {
        vec![
            AlgorithmSpec::RemoveCom {
                initial: 1,
                every: 10,
            },
            AlgorithmSpec::Forcefield(ForcefieldPlan {
                cutoff_short: 0.8,
                cutoff_long: 1.4,
                pairlist_every: 5,
                grid_size: 0.0,
                pairlist_algorithm: 0,
                pairlist_type: 0,
                rf_cutoff: 1.4,
                epsilon_rf: 61.0,
                kappa: 0.0,
                rf_excluded: true,
                four_pi_eps_i: newen_core::units::four_pi_eps_i,
                bonds: true,
                covalent_form: (false, false, false),
                angles: true,
                impropers: true,
                dihedrals: true,
                parallel: false,
                atoms_per_solvent: 3,
                virial: VirialKind::None,
                position_restraints: None,
                distance_restraints: None,
                perturbation: None,
            }),
            AlgorithmSpec::Orchestrator { terms: Vec::new() },
            AlgorithmSpec::Eds(EdsPlan {
                form: 1,
                num_states: 2,
                s: vec![0.5],
                pairs: Vec::new(),
                eir: vec![0.0, 10.0],
                beta: 1.0 / (newen_core::units::kB * 300.0),
                aeds: None,
            }),
            AlgorithmSpec::LeapFrogVelocity,
            AlgorithmSpec::Thermostat {
                algorithm: ThermostatAlgorithm::Berendsen,
                baths: vec![ThermostatBath {
                    temperature: 300.0,
                    tau: 0.1,
                    dof: 1944.0,
                }],
                ranges: vec![BathRange {
                    last_atom: 647,
                    com_bath: 0,
                    ir_bath: 0,
                }],
            },
            AlgorithmSpec::LeapFrogPosition,
            AlgorithmSpec::Shake {
                tolerance: 1e-4,
                max_iterations: 1000,
                solute: SoluteConstraints::AllBonds,
                include_solvent: true,
                initial_positions: false,
                initial_velocities: false,
            },
            AlgorithmSpec::Settle,
            AlgorithmSpec::Lincs {
                solute: SoluteConstraints::AllBonds,
                order_solute: 4,
                order_solvent: 4,
                constrain_solute: true,
                constrain_solvent: false,
            },
            AlgorithmSpec::Rottrans { last: 1 },
            AlgorithmSpec::SteepestDescent {
                tolerance: 0.1,
                step0: 0.01,
                step_max: 0.05,
                min_steps: 1,
                force_limit: 0.0,
            },
            AlgorithmSpec::TemperatureCalculation {
                n_baths: 1,
                ranges: Vec::new(),
            },
            AlgorithmSpec::PressureCalculation {
                virial: VirialKind::Molecular,
            },
            AlgorithmSpec::Barostat {
                pressure0: 1.0,
                compressibility: 4.575e-4,
                tau: 0.5,
            },
            AlgorithmSpec::EnergyCalculation,
        ]
    }

    /// Ordering invariants per kind (GROMOS `md.cc` step order).
    pub fn rules(kind: &str) -> KindRules {
        match kind {
            "remove_com" => KindRules {
                unique: true,
                first: true,
                excludes: &["steepest_descent"],
                ..NONE
            },
            "forcefield" => KindRules {
                unique: true,
                required: true,
                ..NONE
            },
            "orchestrator" => KindRules {
                unique: true,
                must_follow: Some("forcefield"),
                ..NONE
            },
            // gromosXX `create_md_sequence.cc:134-141`: after the force field, before anything
            // that uses the forces
            "eds" => KindRules {
                unique: true,
                after: &["forcefield", "orchestrator"],
                before: &[
                    "leap_frog_velocity",
                    "steepest_descent",
                    "shake",
                    "settle",
                    "lincs",
                    "rottrans",
                ],
                requires: &["forcefield"],
                ..NONE
            },
            "leap_frog_velocity" => KindRules {
                unique: true,
                after: &["forcefield"],
                before: &["leap_frog_position"],
                requires: &["leap_frog_position", "temperature_calculation"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "thermostat" => KindRules {
                unique: true,
                after: &["leap_frog_velocity"],
                before: &["leap_frog_position"],
                requires: &["leap_frog_velocity"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "leap_frog_position" => KindRules {
                unique: true,
                after: &["leap_frog_velocity"],
                requires: &["leap_frog_velocity"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "shake" | "settle" | "lincs" => KindRules {
                unique: true,
                after: &["forcefield"],
                before: &["temperature_calculation", "energy_calculation"],
                ..NONE
            },
            "rottrans" => KindRules {
                unique: true,
                after: &["forcefield"],
                before: &["temperature_calculation", "energy_calculation"],
                ..NONE
            },
            "steepest_descent" => KindRules {
                unique: true,
                after: &["forcefield"],
                excludes: DYNAMICS_ONLY,
                ..NONE
            },
            "temperature_calculation" => KindRules {
                unique: true,
                after: &["leap_frog_position"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "pressure_calculation" => KindRules {
                unique: true,
                after: &["temperature_calculation"],
                before: &["barostat"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "barostat" => KindRules {
                unique: true,
                after: &["pressure_calculation"],
                requires: &["pressure_calculation"],
                excludes: &["steepest_descent"],
                ..NONE
            },
            "energy_calculation" => KindRules {
                unique: true,
                required: true,
                last: true,
                ..NONE
            },
            _ => NONE,
        }
    }
}

impl TermSpec {
    pub const KINDS: &'static [&'static str] = &["schnet", "xtb"];

    pub fn name(&self) -> &'static str {
        match self {
            TermSpec::Schnet { .. } => "schnet",
            TermSpec::Xtb { .. } => "xtb",
        }
    }

    pub fn examples() -> Vec<TermSpec> {
        vec![
            TermSpec::Schnet {
                model: "model.pt".into(),
                cutoff: 0.5,
                elements: vec![8, 1, 1],
                region: "1:res(SOL:a)".into(),
                buffer: None,
                coupling: Coupling::Delta,
            },
            TermSpec::Xtb {
                region: "1:a".into(),
                elements: vec![8, 1, 1],
                gfn: 2,
                charge: 0,
                multiplicity: 1,
                work_dir: None,
                timeout_s: 600,
                coupling: Coupling::Delta,
            },
        ]
    }

    /// The cargo feature this term needs, if any.
    pub fn feature(&self) -> Option<&'static str> {
        match self {
            TermSpec::Schnet { .. } => Some("ml"),
            TermSpec::Xtb { .. } => None,
        }
    }

    /// Whether this term's provider reports a virial (a barostat needs one).
    pub fn provides_virial(&self) -> bool {
        match self {
            TermSpec::Schnet { .. } => false,
            TermSpec::Xtb { .. } => false,
        }
    }

    pub fn coupling(&self) -> Coupling {
        match self {
            TermSpec::Schnet { coupling, .. } | TermSpec::Xtb { coupling, .. } => *coupling,
        }
    }
}

/// Stage 1: the recipe as an ordered, fully resolved list of algorithms.
///
/// Needs the topology (constraint counts → degrees of freedom, solvent shape); the restraint
/// file paths come from `recipe.inputs` and become plan values.
/// Resolve a recipe's `MULTIBATH` onto the topology: one [`ThermostatBath`] per bath, with the
/// DOFSET ranges turned into per-bath degrees of freedom, and the ranges themselves — the one
/// atom → bath mapping (gromosXX `Multibath::in_bath`) the thermostat and the kinetic energy use.
///
/// gromosXX allows a DOFSET line to send a range's centre-of-mass and internal degrees of freedom
/// to *different* baths; the velocity scaling of such a split range is not implemented here, so a
/// COM-bath ≠ IR-bath line is refused by name rather than silently coupled to one of them. Every
/// GROMOS input we have seen (including all of the LiveCoMS tutorials) uses COMBATH = IRBATH.
fn thermostat_baths(
    t: &Thermostat,
    topo: &Topology,
    sel: &ConstraintSelection,
    ntc_mode: NtcMode,
    ndfmin: i32,
    total: f64,
    rottrans: bool,
) -> Result<(Vec<ThermostatBath>, Vec<BathRange>), RunError> {
    let n_baths = t.baths.len();
    let last_atom = topo.num_atoms().saturating_sub(1);
    if n_baths <= 1 {
        let baths = t
            .baths
            .iter()
            .map(|b| ThermostatBath {
                temperature: b.temperature,
                tau: b.tau,
                dof: total,
            })
            .collect();
        let whole = BathRange {
            last_atom,
            com_bath: 0,
            ir_bath: 0,
        };
        return Ok((baths, vec![whole]));
    }
    if t.algorithm != ThermostatAlgorithm::Berendsen {
        return Err(RunError::Recipe(format!(
            "{n_baths} temperature baths with {:?} coupling: only weak coupling (Berendsen) is \
             implemented for more than one bath",
            t.algorithm
        )));
    }
    if t.dof_sets.is_empty() {
        return Err(RunError::Recipe(format!(
            "{n_baths} temperature baths but no DOFSET lines: gromosXX needs one range per bath"
        )));
    }
    let mut ranges = Vec::with_capacity(t.dof_sets.len());
    for set in &t.dof_sets {
        let [last, com, ir] = *set;
        if com != ir {
            return Err(RunError::Recipe(format!(
                "DOFSET line with COMBATH {com} ≠ IRBATH {ir}: splitting a range's centre-of-mass \
                 and internal degrees of freedom over two baths is not implemented"
            )));
        }
        if com == 0 || com > n_baths {
            return Err(RunError::Recipe(format!(
                "DOFSET line refers to bath {com}, but MULTIBATH declares {n_baths}"
            )));
        }
        if last == 0 || last > topo.num_atoms() {
            return Err(RunError::Recipe(format!(
                "DOFSET line's LAST atom {last} is outside the topology ({} atoms)",
                topo.num_atoms()
            )));
        }
        // gromosXX `Multibath::check_state` (multibath.cc:189-196) means to refuse this, but its
        // running maximum is never updated; the DOF count and the thermostat's ranges assume it
        if ranges
            .last()
            .is_some_and(|r: &BathRange| last - 1 < r.last_atom)
        {
            return Err(RunError::Recipe(format!(
                "Multibath not sorted: DOFSET line with LAST {last} follows LAST {}",
                ranges.last().map_or(0, |r| r.last_atom + 1)
            )));
        }
        ranges.push(BathRange {
            last_atom: last - 1,
            com_bath: com - 1,
            ir_bath: ir - 1,
        });
    }
    if ranges.last().map(|r| r.last_atom) != Some(last_atom) {
        return Err(RunError::Recipe(format!(
            "the last DOFSET range ends at atom {} but the topology has {} atoms: gromosXX \
             requires the last bath to reach the last atom",
            ranges.last().map(|r| r.last_atom + 1).unwrap_or(0),
            topo.num_atoms()
        )));
    }
    // gromosXX `Multibath::calculate_degrees_of_freedom` (multibath.cc:161-181)
    for g in &topo.temperature_groups {
        let Some(&first) = g.first() else { continue };
        let baths = in_bath(&ranges, first);
        if let Some(&atom) = g.iter().find(|&&a| in_bath(&ranges, a) != baths) {
            return Err(RunError::Recipe(format!(
                "Multibath: Temperature group distributed over multiple baths (the group of \
                 atoms {}-{} is split at atom {} by the DOFSET ranges)",
                first + 1,
                g.last().map_or(first, |&l| l) + 1,
                atom + 1
            )));
        }
    }
    let dof = bath_dof(topo, sel, ntc_mode, ndfmin, n_baths, &ranges, rottrans);
    let baths = t
        .baths
        .iter()
        .enumerate()
        .map(|(i, b)| ThermostatBath {
            temperature: b.temperature,
            tau: b.tau,
            dof: dof[i],
        })
        .collect();
    Ok((baths, ranges))
}

/// Whether a `LAMBDAS` block sets individual λ: NTIL on in a perturbed run (gromosXX ignores
/// the block otherwise).
fn individual_lambdas_on(recipe: &RunRecipe) -> bool {
    recipe.perturbation.enabled != 0 && ntil_on(recipe)
}

/// NTIL on, perturbed run or not: gromosXX raises its `individual_lambdas` feature flag before
/// it looks at NTG (`in_parameter.cc:4367-4379`), so the feature locks apply either way.
fn ntil_on(recipe: &RunRecipe) -> bool {
    recipe
        .perturbation
        .lambdas
        .as_ref()
        .is_some_and(|l| l.enabled)
}

/// Whether the recipe runs EDS (an `EDS` block with EDS = 1).
pub(crate) fn eds_on(recipe: &RunRecipe) -> bool {
    recipe.eds.as_ref().is_some_and(|e| e.enabled)
}

/// The EDS settings gromosXX refuses (`in_parameter.cc:3843-3961`, `check_parameter.cc:146-175`
/// and the feature locks it leaves in place: unconstrained solvent, flexible SHAKE, triclinic
/// and truncated-octahedral boxes, lattice sums, QM/MM) and the ones newen does not port
/// (shifted reaction field, NLRELE = 4; the coarse-grained and virtual-atom variants need
/// blocks the recipe does not model, refused as such). Checked again on a recipe built in
/// memory, so the Python front-end meets the same rules as an `.imd` file.
pub(crate) fn check_eds(recipe: &RunRecipe) -> Result<(), RunError> {
    check_aeds(recipe)?;
    let Some(e) = recipe.eds.as_ref().filter(|e| e.enabled) else {
        return Ok(());
    };
    let err = |m: String| Err(RunError::Recipe(format!("EDS: {m}")));
    let n = e.num_states;
    if n < 2 {
        return err(format!("NUMSTATES = {n}: must be >= 2"));
    }
    match e.form {
        1 if e.s.len() != 1 => return err(format!("FORM 1 needs one S, has {}", e.s.len())),
        2 if e.s.len() != n * (n - 1) / 2 => {
            return err(format!(
                "FORM 2 with {n} states needs {} S values, has {}",
                n * (n - 1) / 2,
                e.s.len()
            ))
        },
        3 if e.pairs.len() != n - 1 => {
            return err(format!(
                "FORM 3 with {n} states needs {} pairs, has {}",
                n - 1,
                e.pairs.len()
            ))
        },
        1..=3 => {},
        f => return err(format!("FORM = {f}: must be 1, 2 or 3")),
    }
    if let Some(s) =
        e.s.iter()
            .chain(e.pairs.iter().map(|p| &p.2))
            .find(|&&s| s.is_nan() || s <= 0.0)
    {
        return err(format!("smoothness parameter S = {s}: must be > 0"));
    }
    // gromosXX does not check the pair states (it indexes with them, `eds.cc:635, 660`)
    for &(i, j, _) in &e.pairs {
        let ok = |k: i32| k >= 1 && k as usize <= n;
        if !ok(i) || !ok(j) || i == j {
            return err(format!(
                "FORM 3 pair {i} {j}: the states must be two different numbers in 1..={n}"
            ));
        }
    }
    if e.eir.len() != n {
        return err(format!("{n} states but {} EIR values", e.eir.len()));
    }
    let Some(t) = recipe
        .ensemble
        .thermostat
        .as_ref()
        .filter(|t| t.baths.iter().any(|b| b.tau != -1.0))
    else {
        return err(
            "EDS requires temperature coupling (a MULTIBATH bath with TAU != -1; gromosXX also \
             accepts stochastic dynamics, which newen does not implement)"
                .to_string(),
        );
    };
    if t.baths
        .iter()
        .any(|b| b.temperature != t.baths[0].temperature)
    {
        return err("all baths must have the same temperature".to_string());
    }
    if t.baths[0].temperature == 0.0 {
        return err("baths must have non zero temperature with EDS".to_string());
    }
    if recipe.minimisation.method >= 2 {
        return err(format!(
            "NTEM = {}: EDS cannot run with conjugate gradients",
            recipe.minimisation.method
        ));
    }
    if recipe.constraints.solvent_algorithm.is_none() {
        return err(
            "the solvent must be constrained (CONSTRAINT NTCS = SHAKE, LINCS or SETTLE): \
             gromosXX keeps unconstrained solvent locked with EDS, in vacuum too"
                .to_string(),
        );
    }
    if recipe.constraints.solute != SoluteConstraints::None
        && recipe.constraints.solute_algorithm == Some(ConstraintAlgorithm::Settle)
    {
        return err("flexible SHAKE (NTCP = 3) cannot be combined with EDS".to_string());
    }
    if matches!(
        recipe.boundary.kind,
        BoundaryKind::Triclinic | BoundaryKind::TruncatedOctahedron
    ) {
        return err(
            "only vacuum and rectangular boxes (NTB 0, 1): gromosXX keeps the triclinic and \
             truncated-octahedral boxes locked with EDS"
                .to_string(),
        );
    }
    match recipe.forcefield.electrostatics.electrostatics {
        2 | 3 => {
            return err("lattice-sum electrostatics (NLRELE 2/3) cannot be combined".to_string())
        },
        4 => {
            return err(
                "the shifted reaction field (NLRELE = 4) has EDS kernels in gromosXX that \
                 newen does not port"
                    .to_string(),
            )
        },
        _ => {},
    }
    if !recipe.forcefield.terms.is_empty() {
        return err("QM/ML terms cannot be combined with EDS".to_string());
    }
    if recipe.passthrough.contains_key("GAMD") {
        return err("Gaussian accelerated MD (GAMD) cannot be combined with EDS".to_string());
    }
    Ok(())
}

/// The AEDS settings gromosXX refuses (`in_parameter.cc:3968-4140`, `check_parameter.cc:146-175`,
/// the `aeds` feature locks, `check_parameter.cc:2571-2657`), the ones on which it crashes or
/// computes nonsense (refused here instead: EMIN > EMAX, item 1.8 D14; BSTEPS < ASTEPS, D13), and
/// the combinations newen has no reference for (triclinic and truncated-octahedral boxes).
pub(crate) fn check_aeds(recipe: &RunRecipe) -> Result<(), RunError> {
    let Some(a) = recipe.aeds.as_ref() else {
        return Ok(());
    };
    let err = |m: String| Err(RunError::Recipe(format!("AEDS: {m}")));
    if eds_on(recipe) {
        // gromosXX's error fires whatever the AEDS block says (`in_parameter.cc:4029-4045`)
        return err("AEDS cannot be used in combination with EDS (EDS block with EDS = 1)".into());
    }
    if !a.enabled {
        return Ok(());
    }
    let n = a.num_states;
    if n < 2 {
        return err(format!("NUMSTATES = {n}: must be >= 2"));
    }
    let Some(form) = newen_integrators::eds::AedsForm::from_gromos(a.form) else {
        return err(format!("FORM = {}: must be 1, 2, 3, 4, 5 or 6", a.form));
    };
    if a.eir.len() != n {
        return err(format!("{n} states but {} EIR values", a.eir.len()));
    }
    if a.emin > a.emax {
        return err(format!(
            "EMIN = {} is larger than EMAX = {}: gromosXX only warns, then leaves EIR unread \
             and crashes (item 1.8 D14)",
            a.emin, a.emax
        ));
    }
    if !matches!(a.bmaxtype, 1 | 2) {
        return err(format!("BMAXTYPE = {}: must be 1 or 2", a.bmaxtype));
    }
    if a.bmax.is_nan() || a.bmax <= 0.0 {
        return err(format!("BMAX = {}: must be > 0", a.bmax));
    }
    if a.asteps == 0 || a.bsteps == 0 {
        return err(format!(
            "ASTEPS = {}, BSTEPS = {}: both must be > 0",
            a.asteps, a.bsteps
        ));
    }
    if form.searches_offsets() && a.bsteps < a.asteps {
        return err(format!(
            "BSTEPS = {} < ASTEPS = {}: gromosXX computes BSTEPS − ASTEPS unsigned, so the \
             memory time wraps to about 4·10⁹ steps after the first step and the offsets freeze \
             (item 1.8 D13); use BSTEPS >= ASTEPS",
            a.bsteps, a.asteps
        ));
    }
    let Some(t) = recipe
        .ensemble
        .thermostat
        .as_ref()
        .filter(|t| t.baths.iter().any(|b| b.tau != -1.0))
    else {
        return err(
            "AEDS requires temperature coupling (a MULTIBATH bath with TAU != -1; gromosXX also \
             accepts stochastic dynamics, which newen does not implement)"
                .to_string(),
        );
    };
    if t.baths
        .iter()
        .any(|b| b.temperature != t.baths[0].temperature)
    {
        return err("all baths must have the same temperature".to_string());
    }
    if t.baths[0].temperature == 0.0 {
        return err("baths must have non zero temperature with AEDS".to_string());
    }
    if recipe.minimisation.method != 0 {
        return err(format!(
            "NTEM = {}: gromosXX keeps energy minimisation locked with AEDS",
            recipe.minimisation.method
        ));
    }
    if matches!(
        recipe.boundary.kind,
        BoundaryKind::Triclinic | BoundaryKind::TruncatedOctahedron
    ) {
        return err(
            "only vacuum and rectangular boxes (NTB 0, 1): gromosXX allows the triclinic and \
             truncated-octahedral boxes with AEDS, but newen has no EDS reference in them"
                .to_string(),
        );
    }
    match recipe.forcefield.electrostatics.electrostatics {
        2 | 3 => {
            return err("lattice-sum electrostatics (NLRELE 2/3) cannot be combined".to_string())
        },
        4 => {
            return err(
                "the shifted reaction field (NLRELE = 4) has EDS kernels in gromosXX that \
                 newen does not port"
                    .to_string(),
            )
        },
        _ => {},
    }
    if !recipe.forcefield.terms.is_empty() {
        return err("QM/ML terms cannot be combined with AEDS".to_string());
    }
    if recipe.passthrough.contains_key("GAMD") {
        return err("Gaussian accelerated MD (GAMD) cannot be combined with AEDS".to_string());
    }
    if ntil_on(recipe) {
        return err(
            "LAMBDAS on (individual λ): gromosXX keeps individual lambdas locked with AEDS"
                .to_string(),
        );
    }
    if recipe
        .perturbation
        .precalclam
        .is_some_and(|pc| pc.nr_lambdas > 0)
    {
        return err("PRECALCLAM: gromosXX keeps precalculated lambdas locked with AEDS".into());
    }
    Ok(())
}

/// The plan of an EDS or AEDS run (`None` without either); `k_boltzmann` from the topology.
fn eds_plan(recipe: &RunRecipe, k_boltzmann: f64) -> Option<EdsPlan> {
    let temperature = recipe
        .ensemble
        .thermostat
        .as_ref()?
        .baths
        .first()?
        .temperature;
    if let Some(a) = recipe.aeds.as_ref().filter(|a| a.enabled) {
        return Some(EdsPlan {
            form: 1,
            num_states: a.num_states,
            s: vec![1.0],
            pairs: Vec::new(),
            eir: a.eir.clone(),
            beta: 1.0 / (k_boltzmann * temperature),
            aeds: Some(AedsPlan {
                form: a.form,
                emax: a.emax,
                emin: a.emin,
                restart: !a.ntiaedss && a.form != 1,
                fullemin: !a.restremin,
                bmax_type: a.bmaxtype as u32,
                bmax: a.bmax,
                asteps: a.asteps,
                bsteps: a.bsteps,
                energy_columns: a.energy_columns,
                prevalence: a.prevalence,
                target_state: a.target_state,
            }),
        });
    }
    let e = recipe.eds.as_ref().filter(|e| e.enabled)?;
    Some(EdsPlan {
        form: e.form,
        num_states: e.num_states,
        s: e.s.clone(),
        pairs: e
            .pairs
            .iter()
            .map(|&(i, j, s)| (i as usize - 1, j as usize - 1, s))
            .collect(),
        eir: e.eir.clone(),
        beta: 1.0 / (k_boltzmann * temperature),
        aeds: None,
    })
}

/// Lattice-sum electrostatics (NLRELE 2 Ewald, 3 P3M) for every run, perturbed or not: nothing
/// in the build path calls newen-forces' `electrostatics.rs`/`pme.rs`, and the force field
/// would run the reaction field in their place. Every build path runs it.
pub(crate) fn check_electrostatics(recipe: &RunRecipe) -> Result<(), RunError> {
    let nlrele = recipe.forcefield.electrostatics.electrostatics;
    if matches!(nlrele, 2 | 3) {
        return Err(RunError::Recipe(format!(
            "NONBONDED NLRELE = {nlrele}: lattice-sum electrostatics (NLRELE 2/3) are not \
             implemented (the reaction field would run in their place); use NLRELE 1"
        )));
    }
    Ok(())
}

/// The PERTURBATION settings newen does not implement, refused instead of ignored, and the
/// features gromosXX refuses together with individual λ (`check_parameter.cc:1581-1644`: slow
/// growth, lattice sum, the conjugate-gradient minimiser, QM/MM, GaMD and AEDS stay locked
/// against `individual_lambdas`; plain EDS is unlocked). Every build path runs it.
pub(crate) fn check_perturbation(recipe: &RunRecipe) -> Result<(), RunError> {
    let p = &recipe.perturbation;
    if ntil_on(recipe) {
        let refuse = |what: &str| {
            Err(RunError::Recipe(format!(
                "LAMBDAS on (individual λ) together with {what}: gromosXX refuses this \
                 combination"
            )))
        };
        if matches!(recipe.forcefield.electrostatics.electrostatics, 2 | 3) {
            return refuse("lattice-sum electrostatics (NLRELE 2/3)");
        }
        if recipe.minimisation.method >= 2 {
            return refuse("the conjugate-gradient minimiser (NTEM >= 2)");
        }
        if !recipe.forcefield.terms.is_empty() {
            return refuse("QM/ML terms");
        }
        if recipe.passthrough.contains_key("GAMD") {
            return refuse("Gaussian accelerated MD (GAMD)");
        }
    }
    if let Some(pc) = p.precalclam.filter(|pc| pc.nr_lambdas > 0) {
        check_precalclam(recipe, &pc)?;
    }
    // gromosXX reads DLAMT with the range ">=0" whatever NTG says (`in_parameter.cc:879`)
    if p.dlamt.is_nan() || p.dlamt < 0.0 {
        return Err(RunError::Recipe(format!(
            "PERTURBATION DLAMT = {}: must be >= 0",
            p.dlamt
        )));
    }
    if p.enabled == 0 {
        return Ok(());
    }
    if p.dlamt != 0.0 {
        check_slow_growth(recipe)?;
    }
    if p.nscale != 0 {
        return Err(RunError::Recipe(format!(
            "PERTURBATION NSCALE = {}: scaled interactions (SCALEDINTERACTIONS) are not \
             implemented; set NSCALE = 0",
            p.nscale
        )));
    }
    if individual_lambdas_on(recipe) {
        let negr = recipe.forcefield.energy_groups.count.max(1);
        for r in p.lambdas.iter().flat_map(|l| &l.rows) {
            if r.g1 == 0 || r.g2 < r.g1 || r.g2 > negr {
                return Err(RunError::Recipe(format!(
                    "LAMBDAS row {} {} {}: the groups must satisfy 1 <= NILG1 <= NILG2 <= NEGR \
                     ({negr})",
                    r.kind.0.gromos_name(),
                    r.g1,
                    r.g2,
                )));
            }
        }
    }
    Ok(())
}

/// Slow growth (NTG = 1, DLAMT > 0): the features gromosXX never unlocks against its
/// `slow_growth` flag (`check_parameter.cc:387-388`; the unlocks at `:1603`, `:2530`, `:2594`,
/// `:2748`, `:2819`, `:4217` and `:4489` are commented out; GAMD, the shifted reaction field and
/// position constraints have none) that a recipe can express.
fn check_slow_growth(recipe: &RunRecipe) -> Result<(), RunError> {
    let refuse = |what: &str| {
        Err(RunError::Recipe(format!(
            "slow growth (PERTURBATION DLAMT = {}) together with {what}: gromosXX refuses this \
             combination",
            recipe.perturbation.dlamt
        )))
    };
    if ntil_on(recipe) {
        return refuse("individual λ (LAMBDAS with NTIL on)");
    }
    if recipe
        .perturbation
        .precalclam
        .is_some_and(|pc| pc.nr_lambdas > 0)
    {
        return refuse("PRECALCLAM (NRLAM > 0)");
    }
    if eds_on(recipe) {
        return refuse("EDS");
    }
    if recipe.aeds.as_ref().is_some_and(|a| a.enabled) {
        return refuse("A-EDS");
    }
    match recipe.forcefield.electrostatics.electrostatics {
        2 | 3 => return refuse("lattice-sum electrostatics (NLRELE 2/3)"),
        4 => return refuse("the shifted reaction field (NLRELE 4)"),
        _ => {},
    }
    if recipe.minimisation.method >= 2 {
        return refuse("the conjugate-gradient minimiser (NTEM >= 2)");
    }
    if recipe.passthrough.contains_key("GAMD") {
        return refuse("Gaussian accelerated MD (GAMD)");
    }
    if recipe.forcefield.restraints.position.mode == 3 {
        return refuse("position constraints (NTPOR = 3)");
    }
    Ok(())
}

/// PRECALCLAM with NRLAM > 0: gromosXX's requirements (`check_parameter.cc:261-266`) and the
/// features it locks together with `precalculate_lambdas` (`:4194-4280`) that a recipe can
/// express. NRLAM = 1 is refused: gromosXX divides by NRLAM − 1 and computes NaN.
fn check_precalclam(recipe: &RunRecipe, pc: &crate::recipe::PrecalcLam) -> Result<(), RunError> {
    let err = |m: String| Err(RunError::Recipe(m));
    // The IMD parser's range checks (`parse_precalclam`, gromosXX `in_parameter.cc:4478-4513`),
    // repeated for recipes not built from an IMD file; NRLAM >= 0 is the type's.
    for (name, v) in [("MINLAM", pc.min_lambda), ("MAXLAM", pc.max_lambda)] {
        if !(0.0..=1.0).contains(&v) {
            return err(format!("PRECALCLAM {name} = {v}: must be >= 0 and <= 1"));
        }
    }
    if pc.min_lambda >= pc.max_lambda {
        return err("PRECALCLAM: MINLAM should be smaller than MAXLAM".to_string());
    }
    if recipe.perturbation.enabled == 0 {
        return err("PRECALCLAM cannot be on without perturbation (NTG = 0)".to_string());
    }
    let (ntwe, ntwg) = (recipe.outputs.energy_every, recipe.outputs.ntwg);
    if ntwe == 0 || ntwg == 0 || ntwe != ntwg {
        return err(format!(
            "PRECALCLAM requires NTWE = NTWG > 0 (WRITETRAJ has NTWE = {ntwe}, NTWG = {ntwg})"
        ));
    }
    if pc.nr_lambdas == 1 {
        return err(
            "PRECALCLAM NRLAM = 1: the λ step is (MAXLAM − MINLAM)/(NRLAM − 1), a division by \
             zero (gromosXX computes NaN); use NRLAM = 0 (off) or NRLAM >= 2"
                .to_string(),
        );
    }
    let refuse = |what: &str| {
        err(format!(
            "PRECALCLAM together with {what}: gromosXX refuses this combination"
        ))
    };
    if matches!(recipe.forcefield.electrostatics.electrostatics, 2 | 3) {
        return refuse("lattice-sum electrostatics (NLRELE 2/3)");
    }
    if recipe.minimisation.method >= 2 {
        return refuse("the conjugate-gradient minimiser (NTEM >= 2)");
    }
    if !recipe.forcefield.terms.is_empty() {
        return refuse("QM/ML terms");
    }
    if recipe.passthrough.contains_key("GAMD") {
        return refuse("Gaussian accelerated MD (GAMD)");
    }
    if eds_on(recipe) {
        return refuse("EDS");
    }
    Ok(())
}

/// The individual-λ lock of a plan built elsewhere (`build_sequence_from_plan`): whether
/// LAMBDAS applies is the recipe's, not the plan's editable flag.
pub(crate) fn check_plan_against_recipe(
    recipe: &RunRecipe,
    plan: &[AlgorithmSpec],
) -> Result<(), RunError> {
    let has_terms = plan
        .iter()
        .any(|a| matches!(a, AlgorithmSpec::Orchestrator { terms } if !terms.is_empty()));
    if ntil_on(recipe) && has_terms {
        return Err(RunError::Recipe(
            "LAMBDAS on (individual λ) together with QM/ML terms: gromosXX refuses this \
             combination"
                .to_string(),
        ));
    }
    Ok(())
}

pub fn build_plan(
    recipe: &RunRecipe,
    topo: &Topology,
    constants: &PhysicalConstants,
) -> Result<Vec<AlgorithmSpec>, RunError> {
    let four_pi_eps_i = constants.four_pi_eps_i;
    let inputs = &recipe.inputs;
    let n_atoms = topo.num_atoms();
    let has_solvent = n_atoms > topo.num_solute_atoms();
    let c = &recipe.constraints;
    let sel = ConstraintSelection::from_codes(
        c.solute.ntc(),
        ConstraintAlgorithm::code(c.solute_algorithm),
        ConstraintAlgorithm::code(c.solvent_algorithm),
        has_solvent,
    );
    let ntc_mode = ConstraintSelection::ntc_mode_of(c.solute.ntc());
    let rottrans = recipe.constraints.rottrans_last.is_some();
    let dof =
        total_dof(topo, &sel, ntc_mode, recipe.boundary.ndfmin) - rottrans_dof(topo, rottrans);

    let ff = &recipe.forcefield;
    let position_restraints = if ff.restraints.position.mode > 0 {
        let spec = inputs
            .posresspec
            .clone()
            .ok_or_else(|| RunError::MissingInput {
                flag: "posresspec",
                reason: format!("NTPOR={}", ff.restraints.position.mode),
            })?;
        let reference = if ff.restraints.position.reference >= 1 {
            Some(
                inputs
                    .refpos
                    .clone()
                    .ok_or_else(|| RunError::MissingInput {
                        flag: "refpos",
                        reason: format!("NTPORB={}", ff.restraints.position.reference),
                    })?,
            )
        } else {
            None
        };
        Some(PositionRestraintsPlan {
            spec,
            reference,
            k: ff.restraints.position.k,
        })
    } else {
        None
    };
    let distance_restraints = if ff.restraints.distance.mode != 0 {
        let d = &ff.restraints.distance;
        if d.mode < 0 && d.forcescale != 0 {
            return Err(RunError::Recipe(format!(
                "DISTANCERES FORCESCALE = {}: only FORCESCALE = 0 (d<r>/dr = 1) is implemented \
                 for time-averaged distance restraints",
                d.forcescale
            )));
        }
        if d.mode < 0 && d.ntdira != 0 {
            return Err(RunError::Recipe(
                "DISTANCERES NTDIRA = 1: reading the initial averages (DISRESEXPAVE) from the \
                 configuration is not implemented; use NTDIRA = 0"
                    .to_string(),
            ));
        }
        Some(DistanceRestraintsPlan {
            file: inputs
                .distrest
                .clone()
                .ok_or_else(|| RunError::MissingInput {
                    flag: "distrest",
                    reason: format!("NTDIR={}", d.mode),
                })?,
            mode: d.mode.abs(),
            k: d.k,
            r_linear: d.r_linear,
            time_averaging: (d.mode < 0).then_some(d.tau),
            dt: recipe.control.dt,
            perturbed: recipe.perturbation.enabled != 0,
        })
    } else {
        None
    };
    check_electrostatics(recipe)?;
    check_perturbation(recipe)?;
    check_eds(recipe)?;
    let eds = eds_plan(recipe, constants.kB);
    let perturbation = (recipe.perturbation.enabled != 0).then_some(PerturbationPlan {
        alpha_lj: recipe.perturbation.alpha_lj,
        alpha_c: recipe.perturbation.alpha_c,
        individual_lambdas: individual_lambdas_on(recipe),
        atom_pair_kernel: recipe.perturbation.atom_pair_kernel,
    });

    let is_minimization = recipe.is_minimization();
    let barostat = recipe.ensemble.barostat.as_ref();
    let forcefield = AlgorithmSpec::Forcefield(ForcefieldPlan {
        cutoff_short: ff.pairlist.cutoff_short,
        cutoff_long: ff.pairlist.cutoff_long,
        pairlist_every: ff.pairlist.every,
        grid_size: ff.pairlist.grid_size,
        pairlist_algorithm: recipe.to_imd().algorithm,
        pairlist_type: recipe.to_imd().type_,
        rf_cutoff: ff.electrostatics.rf_cutoff,
        epsilon_rf: ff.electrostatics.epsilon_rf,
        kappa: ff.electrostatics.kappa,
        rf_excluded: ff.electrostatics.self_exclusion != 0,
        four_pi_eps_i,
        bonds: ff.bonds,
        covalent_form: (
            ff.covalent_form.harmonic_bonds,
            ff.covalent_form.harmonic_angles,
            ff.covalent_form.limited_phase_shifts,
        ),
        angles: ff.angles,
        impropers: ff.impropers,
        dihedrals: ff.dihedrals,
        parallel: recipe.execution.parallel.resolve(n_atoms),
        atoms_per_solvent: if topo.solvent_atom_template.is_empty() {
            3
        } else {
            topo.solvent_atom_template.len()
        },
        // The binary only configures the virial for dynamics with pressure coupling.
        virial: match (is_minimization, barostat) {
            (false, Some(b)) => b.virial,
            _ => VirialKind::None,
        },
        position_restraints,
        distance_restraints,
        perturbation,
    });

    let mut plan = Vec::new();
    let constraints = |plan: &mut Vec<AlgorithmSpec>| {
        if sel.shake_enabled() {
            plan.push(AlgorithmSpec::Shake {
                tolerance: c.shake_tolerance,
                max_iterations: 1000,
                solute: if sel.solute_shake {
                    c.solute
                } else {
                    SoluteConstraints::None
                },
                include_solvent: sel.solvent_shake,
                initial_positions: recipe.control.initial_shake >= 1,
                initial_velocities: recipe.control.initial_shake >= 2,
            });
        }
        if sel.settle_enabled {
            plan.push(AlgorithmSpec::Settle);
        }
        if sel.lincs_enabled {
            plan.push(AlgorithmSpec::Lincs {
                solute: c.solute,
                order_solute: c.lincs_order_solute,
                order_solvent: c.lincs_order_solvent,
                constrain_solute: sel.solute_lincs,
                constrain_solvent: sel.solvent_lincs,
            });
        }
        if let Some(last) = recipe.constraints.rottrans_last {
            plan.push(AlgorithmSpec::Rottrans { last });
        }
    };

    if is_minimization {
        plan.push(forcefield);
        if !ff.terms.is_empty() {
            plan.push(AlgorithmSpec::Orchestrator {
                terms: ff.terms.clone(),
            });
        }
        if let Some(e) = eds.clone() {
            plan.push(AlgorithmSpec::Eds(e));
        }
        let m = &recipe.minimisation;
        plan.push(AlgorithmSpec::SteepestDescent {
            tolerance: m.tolerance,
            step0: m.step0,
            step_max: m.step_max,
            min_steps: m.min_steps,
            force_limit: m.force_limit,
        });
        constraints(&mut plan);
        plan.push(AlgorithmSpec::EnergyCalculation);
    } else {
        let com = &recipe.control.com_motion;
        if com.initial >= 1 || com.every != 0 {
            plan.push(AlgorithmSpec::RemoveCom {
                initial: com.initial,
                every: com.every,
            });
        }
        // gromosXX inserts the lattice-shift tracker here for every non-vacuum boundary
        // (`create_md_sequence.cc`), so molecules are wrapped before the forces are evaluated.
        if recipe.boundary.kind != BoundaryKind::Vacuum {
            plan.push(AlgorithmSpec::LatticeShift);
        }
        plan.push(forcefield);
        if !ff.terms.is_empty() {
            plan.push(AlgorithmSpec::Orchestrator {
                terms: ff.terms.clone(),
            });
        }
        if let Some(e) = eds {
            plan.push(AlgorithmSpec::Eds(e));
        }
        plan.push(AlgorithmSpec::LeapFrogVelocity);
        // with one bath the kinetic energy needs no mapping: only its total feeds the thermostat
        let mut kinetic = AlgorithmSpec::TemperatureCalculation {
            n_baths: 1,
            ranges: Vec::new(),
        };
        if let Some(t) = &recipe.ensemble.thermostat {
            let (baths, ranges) = thermostat_baths(
                t,
                topo,
                &sel,
                ntc_mode,
                recipe.boundary.ndfmin,
                dof,
                rottrans,
            )?;
            if baths.len() > 1 {
                kinetic = AlgorithmSpec::TemperatureCalculation {
                    n_baths: baths.len(),
                    ranges: ranges.clone(),
                };
            }
            if baths.iter().any(|b| b.tau > 0.0) {
                plan.push(AlgorithmSpec::Thermostat {
                    algorithm: t.algorithm,
                    baths,
                    ranges,
                });
            }
        }
        plan.push(AlgorithmSpec::LeapFrogPosition);
        constraints(&mut plan);
        plan.push(kinetic);
        if let Some(b) = barostat {
            plan.push(AlgorithmSpec::PressureCalculation { virial: b.virial });
            plan.push(AlgorithmSpec::Barostat {
                pressure0: b.pressure0[0][0],
                compressibility: b.compressibility[0][0],
                tau: b.tau,
            });
        }
        plan.push(AlgorithmSpec::EnergyCalculation);
    }
    Ok(plan)
}

/// Check the GROMOS ordering invariants of a plan (edited or not).
pub fn validate_plan(plan: &[AlgorithmSpec]) -> Result<(), RunError> {
    let names: Vec<&str> = plan.iter().map(|a| a.name()).collect();
    let err = |msg: String| Err(RunError::InvalidPlan(msg));
    let position = |kind: &str| names.iter().position(|n| *n == kind);
    let count = |kind: &str| names.iter().filter(|n| **n == kind).count();

    if plan.is_empty() {
        return err("the plan is empty".into());
    }
    for kind in AlgorithmSpec::KINDS {
        let rules = AlgorithmSpec::rules(kind);
        let n = count(kind);
        if rules.required && n != 1 {
            return err(format!("{kind} must appear exactly once (found {n})"));
        }
        if rules.unique && n > 1 {
            return err(format!("{kind} may appear only once (found {n})"));
        }
        let Some(i) = position(kind) else { continue };
        if rules.first && i != 0 {
            return err(format!("{kind} must be the first algorithm"));
        }
        if rules.last && i != names.len() - 1 {
            return err(format!("{kind} must be the last algorithm"));
        }
        if let Some(prev) = rules.must_follow {
            if i == 0 || names[i - 1] != prev {
                return err(format!("{kind} must come immediately after {prev}"));
            }
        }
        for other in rules.after {
            if let Some(j) = position(other) {
                if j > i {
                    return err(format!("{kind} must come after {other}"));
                }
            }
        }
        for other in rules.before {
            if let Some(j) = position(other) {
                if j < i {
                    return err(format!("{kind} must come before {other}"));
                }
            }
        }
        for other in rules.excludes {
            if position(other).is_some() {
                return err(format!("{kind} cannot be combined with {other}"));
            }
        }
        for other in rules.requires {
            if position(other).is_none() {
                return err(format!("{kind} requires {other}"));
            }
        }
    }
    // gromosXX refuses QM/MM together with individual λ (`check_parameter.cc`)
    let individual_lambdas = plan.iter().any(|a| {
        matches!(a, AlgorithmSpec::Forcefield(ff)
            if ff.perturbation.as_ref().is_some_and(|p| p.individual_lambdas))
    });
    if individual_lambdas && position("orchestrator").is_some() {
        return err("individual λ (LAMBDAS on) cannot be combined with QM/ML terms".into());
    }
    // Terms: a barostat needs every provider to report a virial; `Replace` coupling needs a
    // zone-aware Forcefield (2.8) — declared, never silently additive.
    for spec in plan {
        if let AlgorithmSpec::Orchestrator { terms } = spec {
            for term in terms {
                if term.coupling() == Coupling::Replace {
                    return err(format!(
                        "term {} with coupling=replace requires the zone-aware Forcefield \
                         (2.8); only coupling=delta is available",
                        term.name()
                    ));
                }
                if position("barostat").is_some() && !term.provides_virial() {
                    return err(format!(
                        "term {} reports no virial and cannot be combined with a barostat",
                        term.name()
                    ));
                }
            }
        }
    }
    Ok(())
}

/// A plan as JSON (the machine form; `md @dump`).
pub fn plan_to_json(plan: &[AlgorithmSpec]) -> Result<String, RunError> {
    serde_json::to_string_pretty(plan).map_err(|e| RunError::Serde(e.to_string()))
}

pub fn plan_from_json(text: &str) -> Result<Vec<AlgorithmSpec>, RunError> {
    serde_json::from_str(text).map_err(|e| RunError::Serde(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_is_complete() {
        let names: Vec<&str> = AlgorithmSpec::examples().iter().map(|a| a.name()).collect();
        assert_eq!(names, AlgorithmSpec::KINDS, "examples() and KINDS disagree");
        let term_names: Vec<&str> = TermSpec::examples().iter().map(|t| t.name()).collect();
        assert_eq!(term_names, TermSpec::KINDS);
    }

    #[test]
    fn examples_round_trip_through_json() {
        for spec in AlgorithmSpec::examples() {
            let text = serde_json::to_string(&spec).unwrap();
            let back: AlgorithmSpec = serde_json::from_str(&text).unwrap();
            assert_eq!(back, spec);
        }
        for term in TermSpec::examples() {
            let text = serde_json::to_string(&term).unwrap();
            let back: TermSpec = serde_json::from_str(&text).unwrap();
            assert_eq!(back, term);
        }
    }

    #[test]
    fn unknown_term_parameter_is_an_error() {
        let text = r#"{"kind":"schnet","model":"m.pt","cutof":0.5,"elements":[8],"region":"a"}"#;
        assert!(serde_json::from_str::<TermSpec>(text).is_err());
        let text = r#"{"kind":"xtb","region":"1:a","elements":[8,1,1],"gfn2":true}"#;
        assert!(serde_json::from_str::<TermSpec>(text).is_err());
    }

    #[test]
    fn xtb_term_defaults_and_kinds() {
        let term: TermSpec =
            serde_json::from_str(r#"{"kind":"xtb","region":"1:a","elements":[8,1,1]}"#).unwrap();
        assert_eq!(
            term,
            TermSpec::Xtb {
                region: "1:a".into(),
                elements: vec![8, 1, 1],
                gfn: 2,
                charge: 0,
                multiplicity: 1,
                work_dir: None,
                timeout_s: 600,
                coupling: Coupling::Delta,
            }
        );
        assert_eq!(term.feature(), None, "xtb needs no cargo feature");
        assert!(!term.provides_virial());
    }

    fn xtb_term(coupling: Coupling) -> TermSpec {
        TermSpec::Xtb {
            region: "1:a".into(),
            elements: vec![8, 1, 1],
            gfn: 2,
            charge: 0,
            multiplicity: 1,
            work_dir: None,
            timeout_s: 600,
            coupling,
        }
    }

    fn nve_plan_with(term: TermSpec) -> Vec<AlgorithmSpec> {
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        vec![
            pick("forcefield"),
            AlgorithmSpec::Orchestrator { terms: vec![term] },
            pick("leap_frog_velocity"),
            pick("leap_frog_position"),
            pick("temperature_calculation"),
            pick("energy_calculation"),
        ]
    }

    #[test]
    fn xtb_delta_coupling_is_a_valid_plan() {
        validate_plan(&nve_plan_with(xtb_term(Coupling::Delta))).unwrap();
    }

    #[test]
    fn xtb_replace_coupling_is_rejected_with_the_2_8_message() {
        let err = validate_plan(&nve_plan_with(xtb_term(Coupling::Replace))).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("coupling=replace") && msg.contains("2.8"),
            "{msg}"
        );
    }

    #[test]
    fn xtb_with_a_barostat_is_rejected_no_virial() {
        let mut plan = dynamics_plan();
        let after_ff = plan.iter().position(|a| a.name() == "forcefield").unwrap() + 1;
        plan.insert(
            after_ff,
            AlgorithmSpec::Orchestrator {
                terms: vec![xtb_term(Coupling::Delta)],
            },
        );
        let msg = validate_plan(&plan).unwrap_err().to_string();
        assert!(msg.contains("xtb") && msg.contains("virial"), "{msg}");
    }

    fn dynamics_plan() -> Vec<AlgorithmSpec> {
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        vec![
            pick("remove_com"),
            pick("forcefield"),
            pick("leap_frog_velocity"),
            pick("thermostat"),
            pick("leap_frog_position"),
            pick("shake"),
            pick("temperature_calculation"),
            pick("pressure_calculation"),
            pick("barostat"),
            pick("energy_calculation"),
        ]
    }

    #[test]
    fn canonical_plans_validate() {
        validate_plan(&dynamics_plan()).unwrap();
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        validate_plan(&[
            pick("forcefield"),
            pick("steepest_descent"),
            pick("shake"),
            pick("energy_calculation"),
        ])
        .unwrap();
    }

    #[test]
    fn violations_are_rejected_with_a_named_reason() {
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        let mut p = dynamics_plan();
        p.remove(9); // no energy_calculation
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("energy_calculation"));

        let mut p = dynamics_plan();
        p.swap(7, 8); // barostat before pressure_calculation
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("pressure_calculation"));

        let mut p = dynamics_plan();
        p.insert(3, pick("orchestrator")); // orchestrator not right after forcefield
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("immediately after"));

        let mut p = dynamics_plan();
        p.insert(2, pick("steepest_descent")); // SD in a dynamics plan
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("cannot be combined"));

        let mut p = dynamics_plan();
        p.push(pick("forcefield")); // second forcefield, and energy_calculation no longer last
        assert!(validate_plan(&p).is_err());
    }

    // ── PERTURBATION / LAMBDAS refusals (1.8) ───────────────────────────────────────────

    fn perturbed_recipe(lambdas_on: bool) -> RunRecipe {
        use crate::recipe::{LambdaKind, LambdaRow, Lambdas};
        let mut r = RunRecipe::default();
        r.perturbation.enabled = 1;
        r.perturbation.lambda = 0.25;
        r.forcefield.energy_groups.count = 2;
        r.perturbation.lambdas = Some(Lambdas {
            enabled: lambdas_on,
            rows: vec![LambdaRow {
                kind: LambdaKind(newen_core::lambdas::InteractionLambda::Lj),
                g1: 1,
                g2: 2,
                a: 0.0,
                b: 0.0,
                c: 0.1,
                d: 0.9,
                e: 0.0,
            }],
        });
        r
    }

    fn refusal(r: &RunRecipe) -> String {
        check_perturbation(r).unwrap_err().to_string()
    }

    #[test]
    fn precalclam_needs_perturbation_equal_output_steps_and_two_points() {
        use crate::recipe::PrecalcLam;
        let mut r = perturbed_recipe(false);
        r.outputs.energy_every = 2;
        r.outputs.ntwg = 2;
        r.perturbation.precalclam = Some(PrecalcLam {
            nr_lambdas: 11,
            min_lambda: 0.0,
            max_lambda: 1.0,
        });
        assert!(check_perturbation(&r).is_ok());
        r.outputs.ntwg = 1;
        assert!(refusal(&r).contains("NTWE = NTWG > 0"), "{}", refusal(&r));
        r.outputs.ntwg = 2;
        r.perturbation.precalclam.as_mut().unwrap().nr_lambdas = 1;
        assert!(refusal(&r).contains("NRLAM = 1"), "{}", refusal(&r));
        r.perturbation.precalclam.as_mut().unwrap().nr_lambdas = 0;
        r.outputs.ntwg = 0;
        assert!(check_perturbation(&r).is_ok(), "NRLAM = 0 is off");
        r.perturbation.precalclam.as_mut().unwrap().nr_lambdas = 3;
        r.outputs.ntwg = 2;
        r.eds = Some(eds_recipe(1, 2));
        assert!(refusal(&r).contains("EDS"), "{}", refusal(&r));
        r.eds.as_mut().unwrap().enabled = false;
        assert!(
            check_perturbation(&r).is_ok(),
            "an EDS block with EDS = 0 locks nothing"
        );
        r.eds = None;
        r.perturbation.enabled = 0;
        assert!(
            refusal(&r).contains("without perturbation"),
            "{}",
            refusal(&r)
        );
    }

    #[test]
    fn precalclam_range_is_checked_on_a_recipe() {
        use crate::recipe::PrecalcLam;
        let mut r = perturbed_recipe(false);
        r.outputs.energy_every = 2;
        r.outputs.ntwg = 2;
        r.perturbation.precalclam = Some(PrecalcLam {
            nr_lambdas: 5,
            min_lambda: 0.0,
            max_lambda: 1.5,
        });
        assert!(refusal(&r).contains("MAXLAM = 1.5"), "{}", refusal(&r));
        r.perturbation.precalclam = Some(PrecalcLam {
            nr_lambdas: 5,
            min_lambda: 0.6,
            max_lambda: 0.4,
        });
        assert!(
            refusal(&r).contains("smaller than MAXLAM"),
            "{}",
            refusal(&r)
        );
    }

    #[test]
    fn slow_growth_refuses_what_gromosxx_refuses() {
        use crate::recipe::PrecalcLam;
        let slow = |r: &RunRecipe, what: &str| {
            let m = refusal(r);
            m.contains("slow growth") && m.contains(what)
        };
        let mut r = perturbed_recipe(false);
        r.perturbation.dlamt = 0.001;
        assert!(check_perturbation(&r).is_ok());
        r.minimisation.method = 1;
        assert!(
            check_perturbation(&r).is_ok(),
            "steepest descent is unlocked"
        );
        r.minimisation.method = 2;
        assert!(slow(&r, "conjugate-gradient"), "{}", refusal(&r));
        r.minimisation.method = 0;
        for (nlrele, what) in [
            (2, "lattice-sum"),
            (3, "lattice-sum"),
            (4, "shifted reaction"),
        ] {
            let mut e = r.clone();
            e.forcefield.electrostatics.electrostatics = nlrele;
            assert!(slow(&e, what), "{}", refusal(&e));
        }
        let mut e = r.clone();
        e.forcefield.restraints.position.mode = 3;
        assert!(slow(&e, "NTPOR = 3"), "{}", refusal(&e));
        let mut e = r.clone();
        e.passthrough.insert("GAMD".into(), vec![]);
        assert!(slow(&e, "GAMD"), "{}", refusal(&e));
        let mut e = r.clone();
        e.eds = Some(eds_recipe(1, 2));
        assert!(slow(&e, "EDS"), "{}", refusal(&e));
        let mut e = r.clone();
        e.aeds = aeds_run(1).aeds;
        assert!(slow(&e, "A-EDS"), "{}", refusal(&e));
        let mut e = r.clone();
        e.outputs.energy_every = 1;
        e.outputs.ntwg = 1;
        e.perturbation.precalclam = Some(PrecalcLam {
            nr_lambdas: 3,
            min_lambda: 0.0,
            max_lambda: 1.0,
        });
        assert!(slow(&e, "PRECALCLAM"), "{}", refusal(&e));
        let mut e = perturbed_recipe(true);
        e.perturbation.dlamt = 0.001;
        assert!(slow(&e, "individual λ"), "{}", refusal(&e));
        // gromosXX's range check on DLAMT holds whatever NTG says
        r.perturbation.dlamt = -0.001;
        assert!(refusal(&r).contains(">= 0"), "{}", refusal(&r));
        r.perturbation.enabled = 0;
        assert!(refusal(&r).contains(">= 0"), "{}", refusal(&r));
        // and without perturbation slow growth is inert
        r.perturbation.dlamt = 0.001;
        r.forcefield.electrostatics.electrostatics = 2;
        assert!(check_perturbation(&r).is_ok());
    }

    #[test]
    fn scaled_interactions_are_refused() {
        let mut r = perturbed_recipe(false);
        assert!(check_perturbation(&r).is_ok());
        r.perturbation.nscale = 1;
        assert!(
            refusal(&r).contains("SCALEDINTERACTIONS"),
            "{}",
            refusal(&r)
        );
        // without perturbation gromosXX never reads them
        r.perturbation.enabled = 0;
        assert!(check_perturbation(&r).is_ok());
    }

    #[test]
    fn individual_lambdas_refuse_what_gromosxx_refuses() {
        let mut r = perturbed_recipe(true);
        assert!(check_perturbation(&r).is_ok());
        r.forcefield.electrostatics.electrostatics = 2;
        assert!(refusal(&r).contains("lattice-sum"), "{}", refusal(&r));
        r.forcefield.electrostatics.electrostatics = 1;
        r.minimisation.method = 2;
        assert!(
            refusal(&r).contains("conjugate-gradient"),
            "{}",
            refusal(&r)
        );
        r.minimisation.method = 1;
        assert!(
            check_perturbation(&r).is_ok(),
            "steepest descent is allowed"
        );
        r.minimisation.method = 0;
        r.forcefield.terms = vec![xtb_term(Coupling::Delta)];
        assert!(refusal(&r).contains("QM/ML"), "{}", refusal(&r));
        // the same block switched off refuses nothing
        let mut off = perturbed_recipe(false);
        off.forcefield.terms = vec![xtb_term(Coupling::Delta)];
        off.forcefield.electrostatics.electrostatics = 2;
        assert!(check_perturbation(&off).is_ok());
    }

    #[test]
    fn ntil_on_locks_features_without_perturbation_too() {
        // gromosXX sets `individual_lambdas` before it checks NTG (in_parameter.cc:4367-4379)
        let mut r = perturbed_recipe(true);
        r.perturbation.enabled = 0;
        assert!(check_perturbation(&r).is_ok());
        r.forcefield.electrostatics.electrostatics = 2;
        assert!(refusal(&r).contains("lattice-sum"), "{}", refusal(&r));
    }

    #[test]
    fn individual_lambdas_refuse_gamd_but_not_eds() {
        let mut r = perturbed_recipe(true);
        r.passthrough
            .insert("GAMD".into(), vec!["1 0 0".into(), "0 0 0".into()]);
        assert!(refusal(&r).contains("GAMD"), "{}", refusal(&r));
        r.passthrough.remove("GAMD");
        r.eds = Some(eds_recipe(1, 2));
        assert!(check_perturbation(&r).is_ok(), "plain EDS is unlocked");
        // LAMBDAS off: nothing to lock
        let mut off = perturbed_recipe(false);
        off.passthrough
            .insert("GAMD".into(), vec!["1 0 0".into(), "0 0 0".into()]);
        assert!(check_perturbation(&off).is_ok());
    }

    // ── EDS refusals (gromosXX read_EDS, check_parameter.cc:146-175, feature locks) ─────

    fn eds_recipe(form: i32, n: usize) -> crate::recipe::Eds {
        crate::recipe::Eds {
            enabled: true,
            form,
            num_states: n,
            s: match form {
                1 => vec![0.5],
                2 => vec![0.3; n * (n - 1) / 2],
                _ => Vec::new(),
            },
            pairs: if form == 3 {
                (1..n as i32).map(|i| (i, i + 1, 0.2)).collect()
            } else {
                Vec::new()
            },
            eir: vec![0.0; n],
        }
    }

    /// An EDS run gromosXX accepts: one coupled bath at 300 K, SHAKE on the solvent.
    fn eds_run(form: i32, n: usize) -> RunRecipe {
        use crate::recipe::{Bath, Thermostat};
        let mut r = RunRecipe {
            eds: Some(eds_recipe(form, n)),
            ..Default::default()
        };
        r.ensemble.thermostat = Some(Thermostat {
            algorithm: ThermostatAlgorithm::Berendsen,
            baths: vec![Bath {
                temperature: 300.0,
                tau: 0.1,
            }],
            dof_sets: Vec::new(),
        });
        r.constraints.solvent_algorithm = Some(ConstraintAlgorithm::Shake);
        r
    }

    fn eds_refusal(r: &RunRecipe) -> String {
        check_eds(r).unwrap_err().to_string()
    }

    #[test]
    fn eds_accepts_the_three_forms() {
        for (form, n) in [(1, 2), (1, 4), (2, 3), (3, 3)] {
            check_eds(&eds_run(form, n)).unwrap();
        }
        let mut off = RunRecipe {
            eds: Some(eds_recipe(1, 2)),
            ..Default::default()
        };
        off.eds.as_mut().unwrap().enabled = false;
        check_eds(&off).unwrap();
    }

    #[test]
    fn eds_refuses_what_gromosxx_refuses() {
        use crate::recipe::Bath;
        let case = |f: &dyn Fn(&mut RunRecipe), what: &str| {
            let mut r = eds_run(1, 2);
            f(&mut r);
            let e = eds_refusal(&r);
            assert!(e.contains(what), "{what}: {e}");
        };
        case(&|r| r.ensemble.thermostat = None, "temperature coupling");
        case(
            &|r| {
                for b in &mut r.ensemble.thermostat.as_mut().unwrap().baths {
                    b.tau = -1.0;
                }
            },
            "temperature coupling",
        );
        case(
            &|r| {
                r.ensemble.thermostat.as_mut().unwrap().baths.push(Bath {
                    temperature: 310.0,
                    tau: 0.1,
                })
            },
            "same temperature",
        );
        case(
            &|r| r.ensemble.thermostat.as_mut().unwrap().baths[0].temperature = 0.0,
            "non zero temperature",
        );
        case(&|r| r.minimisation.method = 2, "conjugate gradients");
        case(&|r| r.constraints.solvent_algorithm = None, "NTCS");
        case(
            &|r| {
                r.constraints.solute = SoluteConstraints::AllBonds;
                r.constraints.solute_algorithm = Some(ConstraintAlgorithm::Settle);
            },
            "flexible SHAKE",
        );
        case(
            &|r| r.boundary.kind = BoundaryKind::TruncatedOctahedron,
            "truncated-octahedral",
        );
        case(&|r| r.boundary.kind = BoundaryKind::Triclinic, "triclinic");
        case(
            &|r| r.forcefield.electrostatics.electrostatics = 2,
            "lattice-sum",
        );
        case(
            &|r| r.forcefield.electrostatics.electrostatics = 4,
            "shifted reaction field",
        );
        case(
            &|r| r.forcefield.terms = vec![xtb_term(Coupling::Delta)],
            "QM/ML",
        );
        case(
            &|r| {
                r.passthrough.insert("GAMD".into(), vec![]);
            },
            "GAMD",
        );
        // the recipe's own consistency (a Python recipe never went through the .imd parser)
        case(&|r| r.eds.as_mut().unwrap().s.clear(), "FORM 1 needs one S");
        case(&|r| r.eds.as_mut().unwrap().s = vec![0.0], "must be > 0");
        case(
            &|r| {
                r.eds.as_mut().unwrap().eir.pop();
            },
            "EIR",
        );
        case(&|r| r.eds.as_mut().unwrap().num_states = 1, "NUMSTATES = 1");
        case(&|r| r.eds.as_mut().unwrap().form = 4, "FORM = 4");
        let mut r = eds_run(2, 3);
        r.eds.as_mut().unwrap().s.pop();
        assert!(eds_refusal(&r).contains("needs 3 S"), "{}", eds_refusal(&r));
        let mut r = eds_run(3, 3);
        r.eds.as_mut().unwrap().pairs[1] = (2, 4, 0.1);
        assert!(eds_refusal(&r).contains("pair 2 4"), "{}", eds_refusal(&r));
        r.eds.as_mut().unwrap().pairs[1] = (2, 2, 0.1);
        assert!(eds_refusal(&r).contains("pair 2 2"), "{}", eds_refusal(&r));
        r.eds.as_mut().unwrap().pairs.pop();
        assert!(
            eds_refusal(&r).contains("needs 2 pairs"),
            "{}",
            eds_refusal(&r)
        );
    }

    // ── AEDS refusals (gromosXX read_AEDS, check_parameter.cc:146-175, the aeds locks) ───────
    /// An AEDS run gromosXX accepts (`aladip_vacuum_aeds_f4`'s settings): no solvent
    /// constraints needed, unlike EDS.
    fn aeds_run(form: i32) -> RunRecipe {
        let mut r = eds_run(1, 3);
        r.eds = None;
        r.constraints.solvent_algorithm = None;
        r.aeds = Some(crate::recipe::Aeds {
            enabled: true,
            form,
            num_states: 3,
            emax: 0.0,
            emin: 0.0,
            eir: vec![0.0, 0.0, -3.0],
            ntiaedss: true,
            restremin: true,
            bmaxtype: 1,
            bmax: 0.1,
            asteps: 20,
            bsteps: 20,
            energy_columns: Default::default(),
            prevalence: Default::default(),
            target_state: Default::default(),
        });
        r
    }

    #[test]
    fn aeds_accepts_the_six_forms_and_plans_single_s() {
        for form in 1..=6 {
            check_eds(&aeds_run(form)).unwrap();
        }
        let plan = eds_plan(&aeds_run(4), 0.00831441).unwrap();
        assert_eq!((plan.form, plan.s.clone()), (1, vec![1.0]));
        let a = plan.aeds.unwrap();
        assert_eq!((a.form, a.restart, a.fullemin), (4, false, false));
        // NTIAEDSS = 0 continues the configuration's state, except in FORM 1
        let mut r = aeds_run(2);
        r.aeds.as_mut().unwrap().ntiaedss = false;
        assert!(eds_plan(&r, 0.00831441).unwrap().aeds.unwrap().restart);
        let mut r = aeds_run(1);
        r.aeds.as_mut().unwrap().ntiaedss = false;
        assert!(!eds_plan(&r, 0.00831441).unwrap().aeds.unwrap().restart);
        // BSTEPS < ASTEPS is harmless where no memory time is used (FORMS 1 and 3)
        for form in [1, 3] {
            let mut r = aeds_run(form);
            r.aeds.as_mut().unwrap().bsteps = 10;
            check_eds(&r).unwrap();
        }
        // AEDS = 0: nothing to check, no plan
        let mut off = aeds_run(4);
        off.aeds.as_mut().unwrap().enabled = false;
        check_eds(&off).unwrap();
        assert!(eds_plan(&off, 0.00831441).is_none());
    }

    #[test]
    fn aeds_refuses_what_gromosxx_refuses_or_crashes_on() {
        let case = |f: &dyn Fn(&mut RunRecipe), what: &str| {
            let mut r = aeds_run(4);
            f(&mut r);
            let e = check_eds(&r).unwrap_err().to_string();
            assert!(e.contains(what) && e.contains("AEDS"), "{what}: {e}");
        };
        fn a(r: &mut RunRecipe) -> &mut crate::recipe::Aeds {
            r.aeds.as_mut().unwrap()
        }
        case(&|r| r.eds = Some(eds_recipe(1, 3)), "combination with EDS");
        case(
            &|r| {
                r.eds = Some(eds_recipe(1, 3));
                a(r).enabled = false;
            },
            "combination with EDS",
        );
        case(&|r| a(r).emin = 1.0, "EMIN = 1 is larger than EMAX = 0");
        case(&|r| a(r).bsteps = 10, "BSTEPS = 10 < ASTEPS = 20");
        case(&|r| a(r).form = 7, "FORM = 7");
        case(&|r| a(r).num_states = 1, "NUMSTATES = 1");
        case(
            &|r| a(r).eir.pop().map(|_| ()).unwrap(),
            "3 states but 2 EIR",
        );
        case(&|r| a(r).bmaxtype = 3, "BMAXTYPE = 3");
        case(&|r| a(r).bmax = 0.0, "BMAX = 0");
        case(&|r| a(r).asteps = 0, "both must be > 0");
        case(&|r| r.ensemble.thermostat = None, "temperature coupling");
        case(
            &|r| r.ensemble.thermostat.as_mut().unwrap().baths[0].temperature = 0.0,
            "non zero temperature",
        );
        case(&|r| r.minimisation.method = 1, "NTEM = 1");
        case(&|r| r.boundary.kind = BoundaryKind::Triclinic, "triclinic");
        case(
            &|r| r.forcefield.electrostatics.electrostatics = 2,
            "lattice-sum",
        );
        case(
            &|r| {
                r.passthrough.insert("GAMD".into(), Vec::new());
            },
            "GAMD",
        );
        case(
            &|r| {
                r.perturbation.precalclam = Some(crate::recipe::PrecalcLam {
                    nr_lambdas: 3,
                    min_lambda: 0.0,
                    max_lambda: 1.0,
                })
            },
            "PRECALCLAM",
        );
    }

    #[test]
    fn eds_plan_resolves_beta_and_zero_based_pairs_and_sits_after_the_force_field() {
        let mut r = eds_run(3, 3);
        r.ensemble.thermostat.as_mut().unwrap().baths[0].temperature = 250.0;
        let kb = 0.0083;
        let e = eds_plan(&r, kb).unwrap();
        assert_eq!(e.beta, 1.0 / (kb * 250.0));
        assert_eq!(e.pairs, [(0, 1, 0.2), (1, 2, 0.2)]);
        let ex = AlgorithmSpec::examples();
        let pick = |kind: &str| ex.iter().find(|a| a.name() == kind).unwrap().clone();
        let mut p = dynamics_plan();
        p.insert(2, pick("eds"));
        validate_plan(&p).unwrap();
        let mut p = dynamics_plan();
        p.insert(1, pick("eds"));
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("eds must come after"));
        let mut p = dynamics_plan();
        p.insert(3, pick("eds"));
        assert!(validate_plan(&p)
            .unwrap_err()
            .to_string()
            .contains("before leap_frog_velocity"));
    }

    #[test]
    fn a_plan_with_terms_is_refused_by_the_recipe_not_the_plan_flag() {
        // a plan built with LAMBDAS off (flag false) and edited to carry a term
        let plan = nve_plan_with(xtb_term(Coupling::Delta));
        validate_plan(&plan).unwrap();
        let err = check_plan_against_recipe(&perturbed_recipe(true), &plan)
            .unwrap_err()
            .to_string();
        assert!(err.contains("QM/ML"), "{err}");
        assert!(check_plan_against_recipe(&perturbed_recipe(false), &plan).is_ok());
    }

    #[test]
    fn lambdas_rows_must_name_existing_groups() {
        let mut r = perturbed_recipe(true);
        r.forcefield.energy_groups.count = 1;
        assert!(
            refusal(&r).contains("NILG1 <= NILG2 <= NEGR"),
            "{}",
            refusal(&r)
        );
    }

    #[test]
    fn an_edited_plan_cannot_add_terms_to_individual_lambdas() {
        let mut plan = nve_plan_with(xtb_term(Coupling::Delta));
        let AlgorithmSpec::Forcefield(ff) = &mut plan[0] else {
            panic!("forcefield first")
        };
        ff.perturbation = Some(PerturbationPlan {
            alpha_lj: 0.5,
            alpha_c: 0.5,
            individual_lambdas: true,
            atom_pair_kernel: Default::default(),
        });
        let err = validate_plan(&plan).unwrap_err().to_string();
        assert!(err.contains("individual λ"), "{err}");
    }

    #[test]
    fn lambdas_kind_is_a_gromos_name_or_number() {
        use crate::recipe::LambdaKind;
        let k: LambdaKind = serde_json::from_str(r#""vdw_soft""#).unwrap();
        assert_eq!(k.0, newen_core::lambdas::InteractionLambda::LjSoftness);
        let k: LambdaKind = serde_json::from_str(r#""12""#).unwrap();
        assert_eq!(serde_json::to_string(&k).unwrap(), r#""mass""#);
        assert!(serde_json::from_str::<LambdaKind>(r#""angleres""#).is_err());
        // a bare number, as in a .imd, from JSON (Python dicts) and TOML
        let k: LambdaKind = serde_json::from_str("5").unwrap();
        assert_eq!(k.0, newen_core::lambdas::InteractionLambda::Lj);
        assert!(serde_json::from_str::<LambdaKind>("13").is_err());
        #[derive(serde::Deserialize)]
        struct Row {
            kind: LambdaKind,
        }
        let row: Row = toml::from_str("kind = 12").unwrap();
        assert_eq!(row.kind.0, newen_core::lambdas::InteractionLambda::Mass);
    }
}
