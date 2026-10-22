//! `RunRecipe` — the one description of "what to run", as plain data.
//!
//! The GROMOS `.imd` file and the Python objects are both *front-ends* to this struct
//! (3.9): `RunRecipe::from_imd` parses one, `to_imd`/`to_imd_string` write one back,
//! serde makes it a dict/JSON/TOML. Groups are orthogonal concerns; field names stay close to
//! GROMOS where the semantics are GROMOS's (this is recipe version 1 — renames are a version
//! bump, not a silent change).
//!
//! Two rules keep it honest:
//! - **Lossless.** Every field of every modelled `.imd` block is here, so
//!   `from_imd(to_imd(from_imd(x))) == from_imd(x)` and gromosXX accepts what `to_imd` writes.
//! - **Absent means absent.** A missing optional block is what gromosXX does without it (no
//!   bath, no COM removal, no pressure coupling …) — never a "default" value smuggled in; the
//! parameter defaults in `newen-io` were aligned to that (3.9 A18).
//!
//! Blocks the recipe does not model (GAMD, QMMM, …) are refused unless the caller lists
//! them in a [`PassthroughPolicy`]; allowed ones ride along verbatim in `passthrough`.
//! REPLICA_EDS is refused whatever the policy: it acts together with the modelled EDS block.

use std::collections::{BTreeMap, BTreeSet};

use newen_core::lambdas::InteractionLambda;
use newen_io::imd::{
    AedsBlock, EdsBlock, EdsStatePair, ImdParameters, LambdasBlock, LambdasRow, NonbondedExtra,
    PrecalcLamBlock, PressureParameters, TempBathParameters,
};
use newen_io::imd_write::{write_imd, MODELLED_BLOCKS};
use serde::{Deserialize, Serialize};

use crate::{inputs::ParallelPolicy, RunError, RunInputs};

/// Schema version written into every recipe (`to_dict`/TOML/JSON); bump on any renaming.
pub const RECIPE_VERSION: u32 = 1;

/// Notes from `from_imd`: what was absent and defaulted, what passed through. Never an error.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Diagnostics {
    pub notes: Vec<String>,
}

/// Which unmodelled blocks a caller accepts as inert passthrough.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct PassthroughPolicy {
    pub allowed: BTreeSet<String>,
}

impl PassthroughPolicy {
    pub fn allow<I: IntoIterator<Item = S>, S: Into<String>>(names: I) -> Self {
        Self {
            allowed: names.into_iter().map(Into::into).collect(),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Groups
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct SystemSpec {
    /// NPM — number of solute molecules (GROMOS keeps it at 1).
    pub npm: usize,
    /// NSM — solvent molecules. A *hint*: the coordinate file decides (`prepare_system`).
    pub nsm: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct ComMotion {
    /// NTICOM: 0 = none, 1 = translation, 2 = translation + rotation at start-up.
    pub initial: i32,
    /// NSCM: > 0 remove translation every NSCM steps, < 0 translation + rotation every |NSCM|,
    /// 0 = never (gromosXX's meaning of an absent COMTRANSROT block).
    pub every: i32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Control {
    /// NSTLIM
    pub steps: usize,
    /// T — initial time (ps)
    pub t0: f64,
    /// DT (ps)
    pub dt: f64,
    /// IG — random seed
    pub seed: i64,
    /// NTIVEL = 1: generate Maxwell–Boltzmann velocities at `initial_temperature`.
    pub generate_velocities: bool,
    /// TEMPI (K) — used when `generate_velocities`; carried otherwise.
    pub initial_temperature: f64,
    /// NTISHK: 0 none, 1 shake initial positions, 2 positions + velocities.
    pub initial_shake: i32,
    pub ntinht: i32,
    pub ntinhb: i32,
    pub ntishi: i32,
    pub ntirtc: i32,
    pub ntisti: i32,
    pub com_motion: ComMotion,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundaryKind {
    Vacuum,
    Rectangular,
    Triclinic,
    TruncatedOctahedron,
}

impl BoundaryKind {
    pub fn from_ntb(ntb: i32) -> Result<Self, RunError> {
        Ok(match ntb {
            0 => Self::Vacuum,
            1 => Self::Rectangular,
            2 => Self::Triclinic,
            -1 => Self::TruncatedOctahedron,
            other => return Err(RunError::Recipe(format!("unknown NTB={other}"))),
        })
    }
    pub fn ntb(self) -> i32 {
        match self {
            Self::Vacuum => 0,
            Self::Rectangular => 1,
            Self::Triclinic => 2,
            Self::TruncatedOctahedron => -1,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Boundary {
    pub kind: BoundaryKind,
    /// NDFMIN — degrees of freedom subtracted from the first bath.
    pub ndfmin: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PairlistAlgorithm {
    Standard,
    Grid,
    GridCell,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PairlistKind {
    ChargeGroup,
    Atomic,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct PairlistSpec {
    pub algorithm: PairlistAlgorithm,
    /// NSNB — rebuild every N steps
    pub every: usize,
    /// RCUTP (nm)
    pub cutoff_short: f64,
    /// RCUTL (nm)
    pub cutoff_long: f64,
    /// SIZE (nm); 0 = auto
    pub grid_size: f64,
    pub kind: PairlistKind,
}

/// NONBONDED lines 3–7 — lattice-sum settings, inert for reaction-field runs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct LatticeSum {
    pub nshape: i32,
    pub ashape: f64,
    pub na2clc: i32,
    pub tola2: f64,
    pub epsls: f64,
    pub nk: [usize; 3],
    pub nk2: f64,
    pub ng: [usize; 3],
    pub nasord: i32,
    pub nfdord: i32,
    pub nalias: i32,
    pub nspord: i32,
    pub nqeval: i32,
    pub faccur: f64,
    pub nrdgrd: i32,
    pub nwrgdr: i32,
    pub nlrlj: i32,
    pub slvdns: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct NonbondedSpec {
    /// NLRELE: 0 cutoff, 1 reaction field, 2 Ewald, 3 P3M
    pub electrostatics: i32,
    /// APPAK — reaction-field inverse Debye length (nm⁻¹)
    pub kappa: f64,
    /// RCRF (nm)
    pub rf_cutoff: f64,
    /// EPSRF
    pub epsilon_rf: f64,
    /// NSLFEXCL — 0: no reaction-field term for the excluded pairs and no self term (GROMOS96);
    /// 1: include them (gromosXX default). Consumed by the force field (`ForcefieldPlan::rf_excluded`).
    pub self_exclusion: i32,
    pub lattice_sum: LatticeSum,
    /// Unused PME leftovers carried for the round trip.
    pub pme_order: usize,
    pub pme_alpha: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct EnergyGroups {
    /// NEGR
    pub count: usize,
    /// NRE(1..NEGR) — last atom of each group (1-indexed); empty = one group over everything.
    pub last_atoms: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct PositionRestraintsSpec {
    /// NTPOR: 0 off, 1 restrain (CPOR), 2 B-factor weighted, 3 constrain
    pub mode: i32,
    /// NTPORB: 0 reference = start-up positions, 1 read reference positions
    pub reference: i32,
    /// NTPORS
    pub scaling: i32,
    /// CPOR (kJ mol⁻¹ nm⁻²)
    pub k: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct DistanceRestraintsSpec {
    /// NTDIR: 0 off, ±1 instantaneous, ±2 with W0; negative = time-averaged
    pub mode: i32,
    pub ntdira: i32,
    /// CDIR
    pub k: f64,
    /// DIR0 — linear region threshold (nm)
    pub r_linear: f64,
    /// TAUDIR (ps)
    pub tau: f64,
    pub forcescale: i32,
    pub vdir: i32,
    /// NTWDIR
    pub write_every: i32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct RestraintsSpec {
    pub position: PositionRestraintsSpec,
    pub distance: DistanceRestraintsSpec,
}

/// How an engine term couples to the classical force field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Coupling {
    /// Added on top of the classical forces (what exists today: an honest correction term).
    #[default]
    Delta,
    /// Replaces the classical treatment of the term's own pairs — needs a zone-aware
    /// `Forcefield` (2.8) and is rejected at instantiation until then.
    Replace,
}

/// An additive force provider registered through the `ProviderOrchestrator`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum TermSpec {
    /// A TorchScript SchNetPack model over `region` (an atom-selection string).
    Schnet {
        model: String,
        cutoff: f64,
        elements: Vec<i64>,
        region: String,
        #[serde(default)]
        buffer: Option<String>,
        #[serde(default)]
        coupling: Coupling,
    },
    /// A real `xtb` (GFN-xTB) subprocess over `region`, additive on top of the classical terms
    /// (`coupling: Delta`; no cargo feature). `elements`: atomic number per global atom index
    /// (the `XtbInteraction` convention — it must cover every atom of `region`). `work_dir`
    /// defaults to a per-term directory under the system temp dir, so two xtb terms never
    /// collide; `timeout_s` bounds every xtb call.
    Xtb {
        region: String,
        elements: Vec<i64>,
        #[serde(default = "default_gfn")]
        gfn: u8,
        #[serde(default)]
        charge: i32,
        #[serde(default = "default_multiplicity")]
        multiplicity: i32,
        #[serde(default)]
        work_dir: Option<String>,
        #[serde(default = "default_xtb_timeout_s")]
        timeout_s: u64,
        #[serde(default)]
        coupling: Coupling,
    },
}

fn default_gfn() -> u8 {
    2
}

fn default_multiplicity() -> i32 {
    1
}

fn default_xtb_timeout_s() -> u64 {
    600
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct ForcefieldSpec {
    /// NTF flags
    pub bonds: bool,
    pub angles: bool,
    pub impropers: bool,
    pub dihedrals: bool,
    pub charges: bool,
    pub nonbonded: bool,
    /// COVALENTFORM: which functional form the covalent terms take.
    #[serde(default)]
    pub covalent_form: CovalentForm,
    pub energy_groups: EnergyGroups,
    pub pairlist: PairlistSpec,
    pub electrostatics: NonbondedSpec,
    pub restraints: RestraintsSpec,
    /// Additive providers (none in an `.imd`; Python adds them).
    pub terms: Vec<TermSpec>,
}

/// gromosXX `COVALENTFORM` (NTBBH, NTBAH, NTBDN): which functional form each covalent term takes.
/// The defaults are the GROMOS ones, which is what an input without the block means.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct CovalentForm {
    /// NTBBH = 1: harmonic bonds instead of the quartic form
    pub harmonic_bonds: bool,
    /// NTBAH = 1: harmonic angles instead of the cosine-harmonic form
    pub harmonic_angles: bool,
    /// NTBDN = 1: dihedral phase shifts limited to 0° and 180° instead of arbitrary ones
    pub limited_phase_shifts: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SoluteConstraints {
    None,
    HydrogenBonds,
    AllBonds,
}

impl SoluteConstraints {
    pub fn from_ntc(ntc: i32) -> Result<Self, RunError> {
        Ok(match ntc {
            1 => Self::None,
            2 => Self::HydrogenBonds,
            3 | 4 => Self::AllBonds,
            other => return Err(RunError::Recipe(format!("unknown NTC={other}"))),
        })
    }
    pub fn ntc(self) -> i32 {
        match self {
            Self::None => 1,
            Self::HydrogenBonds => 2,
            Self::AllBonds => 3,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConstraintAlgorithm {
    Shake,
    Lincs,
    Settle,
}

impl ConstraintAlgorithm {
    pub fn from_code(code: i32, what: &str) -> Result<Option<Self>, RunError> {
        Ok(match code {
            0 => None,
            1 => Some(Self::Shake),
            2 => Some(Self::Lincs),
            3 => Some(Self::Settle),
            other => return Err(RunError::Recipe(format!("unknown {what}={other}"))),
        })
    }
    pub fn code(algorithm: Option<Self>) -> i32 {
        match algorithm {
            None => 0,
            Some(Self::Shake) => 1,
            Some(Self::Lincs) => 2,
            Some(Self::Settle) => 3,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct ConstraintsSpec {
    /// NTC
    pub solute: SoluteConstraints,
    /// NTCP — algorithm for the solute constraints (when `solute != None`)
    pub solute_algorithm: Option<ConstraintAlgorithm>,
    /// NTCP0 for SHAKE
    pub shake_tolerance: f64,
    /// NTCP0 for LINCS
    pub lincs_order_solute: usize,
    /// NTCS — solvent algorithm; `None` = solvent unconstrained
    pub solvent_algorithm: Option<ConstraintAlgorithm>,
    /// NTCS0 for SHAKE
    pub solvent_shake_tolerance: f64,
    /// NTCS0 for LINCS
    pub lincs_order_solvent: usize,
    /// ROTTRANS RTCLAST when RTC = 1: the six rigid-body degrees of freedom of atoms
    /// `1..=rottrans_last` are constrained. `None` = no ROTTRANS block, or RTC = 0.
    #[serde(default)]
    pub rottrans_last: Option<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ThermostatAlgorithm {
    Berendsen,
    NoseHoover,
    /// Nosé–Hoover chain of the given length (≥ 2)
    NoseHooverChain(usize),
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bath {
    pub temperature: f64,
    /// TAU (ps); ≤ 0 = this bath is not coupled
    pub tau: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Thermostat {
    pub algorithm: ThermostatAlgorithm,
    pub baths: Vec<Bath>,
    /// DOFSET lines: [LAST, COM-BATH, IR-BATH]; empty = one set over all atoms.
    #[serde(default)]
    pub dof_sets: Vec<[usize; 3]>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VirialKind {
    None,
    Atomic,
    Molecular,
}

impl VirialKind {
    fn from_code(code: i32) -> Result<Self, RunError> {
        Ok(match code {
            0 => Self::None,
            1 => Self::Atomic,
            2 => Self::Molecular,
            other => return Err(RunError::Recipe(format!("unknown VIRIAL={other}"))),
        })
    }
    pub fn code(self) -> i32 {
        match self {
            Self::None => 0,
            Self::Atomic => 1,
            Self::Molecular => 2,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Barostat {
    /// COUPLE: 1 calc (virial only), 2 scale
    pub couple: i32,
    /// SCALE: 1 iso, 2 aniso, 3 full, 4 semianiso
    pub scale: i32,
    pub semi: [i32; 3],
    pub compressibility: [[f64; 3]; 3],
    /// TAUP (ps)
    pub tau: f64,
    pub virial: VirialKind,
    pub pressure0: [[f64; 3]; 3],
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields, default)]
pub struct Ensemble {
    /// `None` = no MULTIBATH block = no temperature coupling.
    pub thermostat: Option<Thermostat>,
    /// `None` = no PRESSURESCALE block (or COUPLE off).
    pub barostat: Option<Barostat>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Minimisation {
    /// NTEM: 0 = dynamics (leap-frog), 1 = steepest descent, 2 = conjugate gradient
    pub method: i32,
    /// NCYC
    pub cycles: usize,
    /// DELE (kJ/mol)
    pub tolerance: f64,
    /// DX0 (nm)
    pub step0: f64,
    /// DXM (nm)
    pub step_max: f64,
    /// NMIN
    pub min_steps: usize,
    /// FLIM
    pub force_limit: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Perturbation {
    /// NTG: 0 off, 1 on
    pub enabled: i32,
    pub nrdgl: i32,
    /// RLAM
    pub lambda: f64,
    /// DLAMT, ps⁻¹: slow growth, λ += DLAMT·Δt after every step (`advance_lambda`)
    pub dlamt: f64,
    /// ALPHLJ
    pub alpha_lj: f64,
    /// ALPHC
    pub alpha_c: f64,
    /// NLAM
    pub nlam: i32,
    pub nscale: i32,
    /// LAMBDAS block: individual λ per interaction kind and energy-group pair; `None` = absent.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lambdas: Option<Lambdas>,
    /// PRECALCLAM block: energies at `nr_lambdas` extra λ points (extended TI); `None` = absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub precalclam: Option<PrecalcLam>,
}

/// The `PRECALCLAM` block: the energies of every perturbed interaction (and their λ
/// derivatives) at `nr_lambdas` points from `min_lambda` to `max_lambda`, written to the energy
/// and free-energy trajectories on every NTWG step (gromosXX extended TI).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrecalcLam {
    /// NRLAM: 0 = off
    pub nr_lambdas: usize,
    /// MINLAM
    pub min_lambda: f64,
    /// MAXLAM
    pub max_lambda: f64,
}

/// The `LAMBDAS` block. Its rows act only when `enabled` and the run is perturbed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Lambdas {
    /// NTIL
    pub enabled: bool,
    #[serde(default)]
    pub rows: Vec<LambdaRow>,
}

/// One `LAMBDAS` row: `L = a·λ⁴ + b·λ³ + c·λ² + d·λ + e` for `kind` between energy groups
/// `g1 ≤ g2` (1-based, as in the file).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LambdaRow {
    pub kind: LambdaKind,
    pub g1: usize,
    pub g2: usize,
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub e: f64,
}

/// A `LAMBDAS` NTLI: one of gromosXX's twelve names (`bond` … `mass`, `vdw` for Lennard-Jones)
/// or its number 1-12 (as a number or a string); serialised as the name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(into = "String")]
pub struct LambdaKind(pub InteractionLambda);

impl<'de> Deserialize<'de> for LambdaKind {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct Visitor;
        impl serde::de::Visitor<'_> for Visitor {
            type Value = LambdaKind;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("a LAMBDAS kind: a gromosXX name or its number 1-12")
            }
            fn visit_str<E: serde::de::Error>(self, s: &str) -> Result<LambdaKind, E> {
                LambdaKind::try_from(s.to_string()).map_err(E::custom)
            }
            fn visit_i64<E: serde::de::Error>(self, n: i64) -> Result<LambdaKind, E> {
                self.visit_str(&n.to_string())
            }
            fn visit_u64<E: serde::de::Error>(self, n: u64) -> Result<LambdaKind, E> {
                self.visit_str(&n.to_string())
            }
        }
        d.deserialize_any(Visitor)
    }
}

impl TryFrom<String> for LambdaKind {
    type Error = String;
    fn try_from(s: String) -> Result<Self, String> {
        InteractionLambda::from_ntli(&s)
            .map(LambdaKind)
            .ok_or_else(|| {
                format!(
                    "unknown LAMBDAS kind {s:?}: expected bond, angle, dihedral, improper, vdw, \
                 vdw_soft, crf, crf_soft, distanceres, distancefield, dihedralres or mass"
                )
            })
    }
}

impl From<LambdaKind> for String {
    fn from(k: LambdaKind) -> String {
        k.0.gromos_name().to_string()
    }
}

impl Lambdas {
    fn from_block(b: &LambdasBlock) -> Self {
        Self {
            enabled: b.enabled,
            rows: b
                .rows
                .iter()
                .map(|r| {
                    let [a, b, c, d, e] = r.coefficients;
                    LambdaRow {
                        kind: LambdaKind(r.kind),
                        g1: r.g1,
                        g2: r.g2,
                        a,
                        b,
                        c,
                        d,
                        e,
                    }
                })
                .collect(),
        }
    }

    fn to_block(&self) -> LambdasBlock {
        LambdasBlock {
            enabled: self.enabled,
            rows: self
                .rows
                .iter()
                .map(|r| LambdasRow {
                    kind: r.kind.0,
                    g1: r.g1,
                    g2: r.g2,
                    coefficients: [r.a, r.b, r.c, r.d, r.e],
                })
                .collect(),
        }
    }
}

/// The `EDS` block: enveloping distribution sampling over `num_states` end states, whose types
/// and charges come from the perturbation topology's `MPERTATOM` block (gromosXX `read_EDS`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Eds {
    /// EDS = 1 (`false`: the block is present and off).
    pub enabled: bool,
    /// FORM: 1 one smoothness parameter, 2 one per state pair, 3 one per listed pair.
    pub form: i32,
    /// NUMSTATES
    pub num_states: usize,
    /// S: FORM 1 one value, FORM 2 N(N−1)/2 in pair order (1,2), (1,3) … (2,3) …
    #[serde(default)]
    pub s: Vec<f64>,
    /// FORM 3: N−1 rows `(i, j, s)`, states 1-based as in the file.
    #[serde(default)]
    pub pairs: Vec<(i32, i32, f64)>,
    /// EIR: the energy offset of every state (kJ/mol).
    #[serde(default)]
    pub eir: Vec<f64>,
}

impl Eds {
    fn from_block(b: &EdsBlock) -> Self {
        Self {
            enabled: b.enabled,
            form: b.form,
            num_states: b.num_states,
            s: b.s.clone(),
            pairs: b.pairs.iter().map(|p| (p.i, p.j, p.s)).collect(),
            eir: b.eir.clone(),
        }
    }

    fn to_block(&self) -> EdsBlock {
        EdsBlock {
            enabled: self.enabled,
            form: self.form,
            num_states: self.num_states,
            s: self.s.clone(),
            pairs: self
                .pairs
                .iter()
                .map(|&(i, j, s)| EdsStatePair { i, j, s })
                .collect(),
            eir: self.eir.clone(),
        }
    }
}

/// The `AEDS` block: accelerated EDS over `num_states` end states (single s = 1, the boost
/// between E_min and E_max, and in FORMS 2-6 the parameter search; gromosXX `read_AEDS`). The
/// end states come from the perturbation topology's `MPERTATOM` block, as for EDS.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Aeds {
    /// AEDS = 1 (`false`: the block is present and off).
    pub enabled: bool,
    /// FORM: 1 fixed; 2 search E^R; 3 search E_max/E_min; 4 search all; 5 E_max/E_min search
    /// with the "advanced" offset update; 6 the advanced offset update alone.
    pub form: i32,
    /// NUMSTATES
    pub num_states: usize,
    /// EMAX (kJ/mol); FORMS 3-5 start the search at V_mix instead.
    #[serde(default)]
    pub emax: f64,
    /// EMIN (kJ/mol)
    #[serde(default)]
    pub emin: f64,
    /// EIR: the energy offset of every state (kJ/mol).
    #[serde(default)]
    pub eir: Vec<f64>,
    /// NTIAEDSS: `true` start a fresh search; `false` continue the configuration's
    /// `AEDSSEARCH` block, which then overrides EMAX, EMIN and EIR.
    #[serde(default)]
    pub ntiaedss: bool,
    /// RESTREMIN: hold E_min at or above the lowest average end-state energy until every state
    /// has been visited.
    #[serde(default)]
    pub restremin: bool,
    /// BMAXTYPE: 1 BMAX in kJ/mol, 2 in standard deviations of the target state's energy.
    #[serde(default)]
    pub bmaxtype: i32,
    /// BMAX: the maximum energy barrier parameter.
    #[serde(default)]
    pub bmax: f64,
    /// ASTEPS: memory time (steps) at the start of each run; FORMS 5/6 offset step ASTEPS·Δt.
    #[serde(default)]
    pub asteps: u32,
    /// BSTEPS: memory time (steps) at the end of each run.
    #[serde(default)]
    pub bsteps: u32,
}

impl Aeds {
    fn from_block(b: &AedsBlock) -> Self {
        Self {
            enabled: b.enabled,
            form: b.form,
            num_states: b.num_states,
            emax: b.emax,
            emin: b.emin,
            eir: b.eir.clone(),
            ntiaedss: b.ntiaedss,
            restremin: b.restremin,
            bmaxtype: b.bmaxtype,
            bmax: b.bmax,
            asteps: b.asteps,
            bsteps: b.bsteps,
        }
    }

    fn to_block(&self) -> AedsBlock {
        AedsBlock {
            enabled: self.enabled,
            form: self.form,
            num_states: self.num_states,
            emax: self.emax,
            emin: self.emin,
            eir: self.eir.clone(),
            ntiaedss: self.ntiaedss,
            restremin: self.restremin,
            bmaxtype: self.bmaxtype,
            bmax: self.bmax,
            asteps: self.asteps,
            bsteps: self.bsteps,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Outputs {
    /// NTWX
    pub trajectory_every: usize,
    pub ntwse: usize,
    /// NTWV — velocity-trajectory write frequency (0 = off)
    pub velocities_every: usize,
    /// NTWF — force-trajectory write frequency (0 = off)
    pub forces_every: usize,
    /// NTWE
    pub energy_every: usize,
    pub ntwg: usize,
    pub ntwb: usize,
    /// NTPR
    pub print_every: usize,
    pub ntpp: usize,
    pub energy_special: bool,
}

/// Not physics, not in the `.imd` (3.9 A15): how the run executes.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields, default)]
pub struct Execution {
    pub parallel: ParallelPolicy,
}

// ---------------------------------------------------------------------------------------------
// The recipe
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct RunRecipe {
    pub version: u32,
    pub title: String,
    pub system: SystemSpec,
    pub control: Control,
    pub boundary: Boundary,
    pub forcefield: ForcefieldSpec,
    pub constraints: ConstraintsSpec,
    pub ensemble: Ensemble,
    pub minimisation: Minimisation,
    pub perturbation: Perturbation,
    /// EDS block; `None` = absent.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub eds: Option<Eds>,
    /// AEDS block; `None` = absent.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aeds: Option<Aeds>,
    pub outputs: Outputs,
    pub execution: Execution,
    /// Auxiliary input files (perturbation topology, restraint specifications) — part of the run
    /// description, so a recipe built in memory is self-describing (3.9 A5).
    pub inputs: RunInputs,
    /// Unmodelled blocks the caller allowed through, verbatim.
    pub passthrough: BTreeMap<String, Vec<String>>,
}

// Group defaults are *derived* from `ImdParameters::default` (one table, 3.9 G7);
// the `Default` impls below exist for serde's `default` and go through that single source.
macro_rules! default_from_imd {
    ($($ty:ident => $field:ident),* $(,)?) => {
        $(impl Default for $ty {
            fn default() -> Self {
                RunRecipe::default().$field
            }
        })*
    };
}
default_from_imd!(
    SystemSpec => system,
    Control => control,
    Boundary => boundary,
    ForcefieldSpec => forcefield,
    ConstraintsSpec => constraints,
    Minimisation => minimisation,
    Perturbation => perturbation,
    Outputs => outputs,
);
impl Default for ComMotion {
    fn default() -> Self {
        RunRecipe::default().control.com_motion
    }
}
impl Default for PairlistSpec {
    fn default() -> Self {
        RunRecipe::default().forcefield.pairlist
    }
}
impl Default for LatticeSum {
    fn default() -> Self {
        RunRecipe::default().forcefield.electrostatics.lattice_sum
    }
}
impl Default for NonbondedSpec {
    fn default() -> Self {
        RunRecipe::default().forcefield.electrostatics
    }
}
impl Default for EnergyGroups {
    fn default() -> Self {
        RunRecipe::default().forcefield.energy_groups
    }
}
impl Default for PositionRestraintsSpec {
    fn default() -> Self {
        RunRecipe::default().forcefield.restraints.position
    }
}
impl Default for DistanceRestraintsSpec {
    fn default() -> Self {
        RunRecipe::default().forcefield.restraints.distance
    }
}
impl Default for RestraintsSpec {
    fn default() -> Self {
        RunRecipe::default().forcefield.restraints
    }
}

impl Default for RunRecipe {
    /// What gromosXX does with every optional block absent — derived from
    /// `ImdParameters::default()`, the single defaults table.
    fn default() -> Self {
        Self::from_imd(&ImdParameters::default())
            .expect("ImdParameters::default() must always map to a recipe")
            .0
    }
}

impl RunRecipe {
    /// Parse a recipe out of GROMOS parameters. Unmodelled blocks are errors (see
    /// [`from_imd_with`](Self::from_imd_with)).
    pub fn from_imd(imd: &ImdParameters) -> Result<(Self, Diagnostics), RunError> {
        Self::from_imd_with(imd, &PassthroughPolicy::default())
    }

    /// Parse a recipe out of GROMOS parameters, allowing the listed unmodelled blocks to pass
    /// through verbatim.
    pub fn from_imd_with(
        imd: &ImdParameters,
        policy: &PassthroughPolicy,
    ) -> Result<(Self, Diagnostics), RunError> {
        let mut diag = Diagnostics::default();
        let present = |block: &str| imd.raw_blocks.contains_key(block);
        let parsed_from_file = !imd.raw_blocks.is_empty();

        if parsed_from_file {
            for required in [
                "SYSTEM",
                "STEP",
                "BOUNDCOND",
                "FORCE",
                "INITIALISE",
                "PAIRLIST",
                "NONBONDED",
            ] {
                if !present(required) {
                    diag.notes.push(format!(
                        "{required} block absent: gromosXX would refuse this input; defaults used"
                    ));
                }
            }
            for (block, meaning) in [
                ("MULTIBATH", "no temperature coupling"),
                ("COMTRANSROT", "no centre-of-mass motion removal"),
                ("PRESSURESCALE", "no pressure coupling"),
                (
                    "CONSTRAINT",
                    "no constraints (NTC=1, NTCS=1 with SHAKE tolerance 1e-4)",
                ),
                ("PERTURBATION", "no free-energy perturbation"),
                ("POSITIONRES", "no position restraints"),
                ("DISTANCERES", "no distance restraints"),
                ("ENERGYMIN", "dynamics, not minimisation"),
                ("WRITETRAJ", "no trajectory/energy output"),
                ("PRINTOUT", "no step printout"),
            ] {
                if !present(block) {
                    diag.notes.push(format!("{block} block absent: {meaning}"));
                }
            }
        }

        if let Some(l) = &imd.lambdas {
            for w in l.warnings(imd.ntg) {
                diag.notes
                    .push(format!("LAMBDAS: {w} (gromosXX warns and runs)"));
            }
        }

        let mut passthrough = BTreeMap::new();
        for (name, lines) in &imd.raw_blocks {
            if MODELLED_BLOCKS.contains(&name.as_str()) {
                continue;
            }
            // refused whatever the policy: it drives the EDS machinery
            if name == "REPLICA_EDS" {
                return Err(RunError::Recipe(
                    "REPLICA_EDS block: replica-exchange EDS is not implemented".to_string(),
                ));
            }
            if !policy.allowed.contains(name) {
                return Err(RunError::UnknownBlock {
                    block: name.clone(),
                });
            }
            diag.notes.push(format!(
                "{name} block passed through verbatim (not modelled)"
            ));
            passthrough.insert(name.clone(), lines.clone());
        }

        let thermostat = match imd.temp_bath.first() {
            None => None,
            Some(b) => {
                let algorithm = match b.algorithm {
                    0 => ThermostatAlgorithm::Berendsen,
                    1 => ThermostatAlgorithm::NoseHoover,
                    n if n >= 2 => {
                        ThermostatAlgorithm::NoseHooverChain((n as usize).max(b.nhc_chain))
                    },
                    other => {
                        return Err(RunError::Recipe(format!(
                            "unknown MULTIBATH algorithm {other}"
                        )))
                    },
                };
                let baths = b
                    .temp0
                    .iter()
                    .zip(&b.tau)
                    .map(|(&temperature, &tau)| Bath { temperature, tau })
                    .collect();
                Some(Thermostat {
                    algorithm,
                    baths,
                    dof_sets: b.dof_sets.clone(),
                })
            },
        };

        let barostat = match (&imd.pressure_parameters, imd.couple_pressure) {
            (Some(pp), true) => Some(Barostat {
                couple: pp.couple.max(1),
                scale: pp.algorithm,
                semi: pp.semi,
                compressibility: pp.compressibility,
                tau: pp.tau_p,
                virial: VirialKind::from_code(pp.virial)?,
                pressure0: pp.pressure0,
            }),
            _ => None,
        };

        let x = &imd.nonbonded_extra;
        let recipe = RunRecipe {
            version: RECIPE_VERSION,
            title: imd.title.clone(),
            system: SystemSpec {
                npm: imd.npm,
                nsm: imd.nsm,
            },
            control: Control {
                steps: imd.nstlim,
                t0: imd.t0,
                dt: imd.dt,
                seed: imd.ig,
                generate_velocities: imd.ntivel == 1,
                initial_temperature: imd.tempi,
                initial_shake: imd.ntishk,
                ntinht: imd.ntinht,
                ntinhb: imd.ntinhb,
                ntishi: imd.ntishi,
                ntirtc: imd.ntirtc,
                ntisti: imd.ntisti,
                com_motion: ComMotion {
                    initial: imd.nticom,
                    every: imd.nscm,
                },
            },
            boundary: Boundary {
                kind: BoundaryKind::from_ntb(imd.ntb)?,
                ndfmin: imd.ndfmin,
            },
            forcefield: ForcefieldSpec {
                bonds: imd.ntf[0] != 0,
                angles: imd.ntf[1] != 0,
                impropers: imd.ntf[2] != 0,
                dihedrals: imd.ntf[3] != 0,
                charges: imd.ntf[4] != 0,
                nonbonded: imd.ntf[5] != 0,
                covalent_form: CovalentForm {
                    harmonic_bonds: imd.ntbbh == 1,
                    harmonic_angles: imd.ntbah == 1,
                    limited_phase_shifts: imd.ntbdn == 1,
                },
                energy_groups: EnergyGroups {
                    count: imd.negr,
                    last_atoms: imd.nre.clone(),
                },
                pairlist: PairlistSpec {
                    algorithm: match imd.algorithm {
                        0 => PairlistAlgorithm::Standard,
                        1 => PairlistAlgorithm::Grid,
                        2 => PairlistAlgorithm::GridCell,
                        other => {
                            return Err(RunError::Recipe(format!(
                                "unknown PAIRLIST algorithm {other}"
                            )))
                        },
                    },
                    every: imd.nsnb,
                    cutoff_short: imd.rcutp,
                    cutoff_long: imd.rcutl,
                    grid_size: imd.size,
                    kind: match imd.type_ {
                        0 => PairlistKind::ChargeGroup,
                        1 => PairlistKind::Atomic,
                        other => {
                            return Err(RunError::Recipe(format!("unknown PAIRLIST type {other}")))
                        },
                    },
                },
                electrostatics: NonbondedSpec {
                    electrostatics: imd.nlrele,
                    kappa: imd.appak,
                    rf_cutoff: imd.rcrf,
                    epsilon_rf: imd.epsrf,
                    self_exclusion: imd.nslfexcl,
                    lattice_sum: LatticeSum {
                        nshape: x.nshape,
                        ashape: x.ashape,
                        na2clc: x.na2clc,
                        tola2: x.tola2,
                        epsls: x.epsls,
                        nk: [imd.grid_x, imd.grid_y, imd.grid_z],
                        nk2: x.nk2,
                        ng: x.ng,
                        nasord: x.nasord,
                        nfdord: x.nfdord,
                        nalias: x.nalias,
                        nspord: x.nspord,
                        nqeval: x.nqeval,
                        faccur: x.faccur,
                        nrdgrd: x.nrdgrd,
                        nwrgdr: x.nwrgdr,
                        nlrlj: x.nlrlj,
                        slvdns: x.slvdns,
                    },
                    pme_order: imd.pme_order,
                    pme_alpha: imd.pme_alpha,
                },
                restraints: RestraintsSpec {
                    position: PositionRestraintsSpec {
                        mode: imd.ntpor,
                        reference: imd.ntporb,
                        scaling: imd.ntpors,
                        k: imd.cpor,
                    },
                    distance: DistanceRestraintsSpec {
                        mode: imd.ntdir,
                        ntdira: imd.ntdira,
                        k: imd.cdir,
                        r_linear: imd.dir0,
                        tau: imd.taudir,
                        forcescale: imd.forcescale,
                        vdir: imd.vdir,
                        write_every: imd.ntwdir,
                    },
                },
                terms: Vec::new(),
            },
            constraints: ConstraintsSpec {
                solute: SoluteConstraints::from_ntc(imd.ntc)?,
                solute_algorithm: ConstraintAlgorithm::from_code(imd.ntcp, "NTCP")?,
                shake_tolerance: imd.shake_tol,
                lincs_order_solute: imd.lincs_order_solute,
                solvent_algorithm: ConstraintAlgorithm::from_code(imd.ntcs, "NTCS")?,
                solvent_shake_tolerance: imd.ntcs0,
                lincs_order_solvent: imd.lincs_order_solvent,
                rottrans_last: (imd.rtc != 0 && imd.rtclast > 0).then_some(imd.rtclast),
            },
            ensemble: Ensemble {
                thermostat,
                barostat,
            },
            minimisation: Minimisation {
                method: imd.ntem,
                cycles: imd.ncyc,
                tolerance: imd.dele,
                step0: imd.dx0,
                step_max: imd.dxm,
                min_steps: imd.nmin,
                force_limit: imd.flim,
            },
            perturbation: Perturbation {
                enabled: imd.ntg,
                nrdgl: imd.nrdgl,
                lambda: imd.rlam,
                dlamt: imd.dlamt,
                alpha_lj: imd.alphlj,
                alpha_c: imd.alphc,
                nlam: imd.nlam,
                nscale: imd.nscale,
                lambdas: imd.lambdas.as_ref().map(Lambdas::from_block),
                precalclam: imd.precalclam.map(|b| PrecalcLam {
                    nr_lambdas: b.nr_lambdas,
                    min_lambda: b.min_lambda,
                    max_lambda: b.max_lambda,
                }),
            },
            eds: imd.eds.as_ref().map(Eds::from_block),
            aeds: imd.aeds.as_ref().map(Aeds::from_block),
            outputs: Outputs {
                trajectory_every: imd.ntwx,
                ntwse: imd.ntwse,
                velocities_every: imd.ntwv,
                forces_every: imd.ntwf,
                energy_every: imd.ntwe,
                ntwg: imd.ntwg,
                ntwb: imd.ntwb,
                print_every: imd.ntpr,
                ntpp: imd.ntpp,
                energy_special: imd.ntwe_special,
            },
            execution: Execution::default(),
            inputs: RunInputs::default(),
            passthrough,
        };
        Ok((recipe, diag))
    }

    /// The GROMOS parameters this recipe describes (the inverse of `from_imd`; `passthrough`
    /// blocks become `raw_blocks` so the writer re-emits them).
    pub fn to_imd(&self) -> ImdParameters {
        let r = self;
        let ls = &r.forcefield.electrostatics.lattice_sum;
        let temp_bath = match &r.ensemble.thermostat {
            None => Vec::new(),
            Some(t) => {
                let (algorithm, nhc_chain) = match t.algorithm {
                    ThermostatAlgorithm::Berendsen => (0, 1),
                    ThermostatAlgorithm::NoseHoover => (1, 1),
                    ThermostatAlgorithm::NoseHooverChain(n) => (n as i32, n),
                };
                vec![TempBathParameters {
                    algorithm,
                    nhc_chain,
                    num_bath_groups: t.baths.len().max(1),
                    temp0: t.baths.iter().map(|b| b.temperature).collect(),
                    tau: t.baths.iter().map(|b| b.tau).collect(),
                    dof: vec![0; t.baths.len().max(1)],
                    dof_sets: t.dof_sets.clone(),
                }]
            },
        };
        let (couple_pressure, pressure_parameters) = match &r.ensemble.barostat {
            None => (false, None),
            Some(b) => (
                true,
                Some(PressureParameters {
                    couple: b.couple,
                    algorithm: b.scale,
                    semi: b.semi,
                    pressure0: b.pressure0,
                    compressibility: b.compressibility,
                    tau_p: b.tau,
                    virial: b.virial.code(),
                }),
            ),
        };
        ImdParameters {
            title: r.title.clone(),
            npm: r.system.npm,
            nsm: r.system.nsm,
            nstlim: r.control.steps,
            t0: r.control.t0,
            dt: r.control.dt,
            ntb: r.boundary.kind.ntb(),
            rtc: i32::from(r.constraints.rottrans_last.is_some()),
            rtclast: r.constraints.rottrans_last.unwrap_or(0),
            ndfmin: r.boundary.ndfmin,
            num_temp_baths: temp_bath.len(),
            temp_bath,
            couple_pressure,
            pressure_parameters,
            force_groups: Vec::new(),
            ntc: r.constraints.solute.ntc(),
            ntcp: ConstraintAlgorithm::code(r.constraints.solute_algorithm),
            ntcs: ConstraintAlgorithm::code(r.constraints.solvent_algorithm),
            shake_tol: r.constraints.shake_tolerance,
            lincs_order_solute: r.constraints.lincs_order_solute,
            lincs_order_solvent: r.constraints.lincs_order_solvent,
            algorithm: match r.forcefield.pairlist.algorithm {
                PairlistAlgorithm::Standard => 0,
                PairlistAlgorithm::Grid => 1,
                PairlistAlgorithm::GridCell => 2,
            },
            nsnb: r.forcefield.pairlist.every,
            rcutp: r.forcefield.pairlist.cutoff_short,
            rcutl: r.forcefield.pairlist.cutoff_long,
            size: r.forcefield.pairlist.grid_size,
            type_: match r.forcefield.pairlist.kind {
                PairlistKind::ChargeGroup => 0,
                PairlistKind::Atomic => 1,
            },
            nlrele: r.forcefield.electrostatics.electrostatics,
            appak: r.forcefield.electrostatics.kappa,
            rcrf: r.forcefield.electrostatics.rf_cutoff,
            epsrf: r.forcefield.electrostatics.epsilon_rf,
            nslfexcl: r.forcefield.electrostatics.self_exclusion,
            ntbbh: i32::from(r.forcefield.covalent_form.harmonic_bonds),
            ntbah: i32::from(r.forcefield.covalent_form.harmonic_angles),
            ntbdn: i32::from(r.forcefield.covalent_form.limited_phase_shifts),
            ntf: [
                r.forcefield.bonds as i32,
                r.forcefield.angles as i32,
                r.forcefield.impropers as i32,
                r.forcefield.dihedrals as i32,
                r.forcefield.charges as i32,
                r.forcefield.nonbonded as i32,
            ],
            negr: r.forcefield.energy_groups.count,
            nre: r.forcefield.energy_groups.last_atoms.clone(),
            nscm: r.control.com_motion.every,
            grid_x: ls.nk[0],
            grid_y: ls.nk[1],
            grid_z: ls.nk[2],
            pme_order: r.forcefield.electrostatics.pme_order,
            pme_alpha: r.forcefield.electrostatics.pme_alpha,
            ntivel: r.control.generate_velocities as i32,
            ntishk: r.control.initial_shake,
            ntinht: r.control.ntinht,
            ntinhb: r.control.ntinhb,
            ntishi: r.control.ntishi,
            nticom: r.control.com_motion.initial,
            ig: r.control.seed,
            tempi: r.control.initial_temperature,
            ntwx: r.outputs.trajectory_every,
            ntwe: r.outputs.energy_every,
            ntwv: r.outputs.velocities_every,
            ntwf: r.outputs.forces_every,
            ntwe_special: r.outputs.energy_special,
            ntpr: r.outputs.print_every,
            ntpor: r.forcefield.restraints.position.mode,
            ntporb: r.forcefield.restraints.position.reference,
            ntpors: r.forcefield.restraints.position.scaling,
            cpor: r.forcefield.restraints.position.k,
            ntdir: r.forcefield.restraints.distance.mode,
            ntdira: r.forcefield.restraints.distance.ntdira,
            cdir: r.forcefield.restraints.distance.k,
            dir0: r.forcefield.restraints.distance.r_linear,
            taudir: r.forcefield.restraints.distance.tau,
            forcescale: r.forcefield.restraints.distance.forcescale,
            vdir: r.forcefield.restraints.distance.vdir,
            ntwdir: r.forcefield.restraints.distance.write_every,
            ntg: r.perturbation.enabled,
            nrdgl: r.perturbation.nrdgl,
            rlam: r.perturbation.lambda,
            dlamt: r.perturbation.dlamt,
            alphlj: r.perturbation.alpha_lj,
            alphc: r.perturbation.alpha_c,
            nlam: r.perturbation.nlam,
            nscale: r.perturbation.nscale,
            lambdas: r.perturbation.lambdas.as_ref().map(Lambdas::to_block),
            precalclam: r.perturbation.precalclam.map(|p| PrecalcLamBlock {
                nr_lambdas: p.nr_lambdas,
                min_lambda: p.min_lambda,
                max_lambda: p.max_lambda,
            }),
            eds: r.eds.as_ref().map(Eds::to_block),
            aeds: r.aeds.as_ref().map(Aeds::to_block),
            ntem: r.minimisation.method,
            nmin: r.minimisation.min_steps,
            dele: r.minimisation.tolerance,
            dx0: r.minimisation.step0,
            dxm: r.minimisation.step_max,
            flim: r.minimisation.force_limit,
            ntirtc: r.control.ntirtc,
            ntisti: r.control.ntisti,
            ntwse: r.outputs.ntwse,
            ntwg: r.outputs.ntwg,
            ntwb: r.outputs.ntwb,
            ntpp: r.outputs.ntpp,
            ncyc: r.minimisation.cycles,
            ntcs0: r.constraints.solvent_shake_tolerance,
            nonbonded_extra: NonbondedExtra {
                nshape: ls.nshape,
                ashape: ls.ashape,
                na2clc: ls.na2clc,
                tola2: ls.tola2,
                epsls: ls.epsls,
                nk2: ls.nk2,
                ng: ls.ng,
                nasord: ls.nasord,
                nfdord: ls.nfdord,
                nalias: ls.nalias,
                nspord: ls.nspord,
                nqeval: ls.nqeval,
                faccur: ls.faccur,
                nrdgrd: ls.nrdgrd,
                nwrgdr: ls.nwrgdr,
                nlrlj: ls.nlrlj,
                slvdns: ls.slvdns,
            },
            raw_blocks: r
                .passthrough
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect(),
        }
    }

    /// The `.imd` text for this recipe (gromosXX-readable). `n_atoms` fills the atom counts a
    /// recipe built in memory does not know (MULTIBATH DOFSET, FORCE NRE).
    pub fn to_imd_string(&self, n_atoms: Option<usize>) -> String {
        write_imd(&self.to_imd(), n_atoms)
    }

    pub fn is_minimization(&self) -> bool {
        self.minimisation.method > 0
    }

    /// Recipe as TOML (the human-readable echo written next to a run's outputs).
    pub fn to_toml(&self) -> Result<String, RunError> {
        toml::to_string_pretty(self).map_err(|e| RunError::Serde(e.to_string()))
    }

    pub fn from_toml(text: &str) -> Result<Self, RunError> {
        toml::from_str(text).map_err(|e| RunError::Serde(e.to_string()))
    }

    /// Recipe as JSON (the machine form; `md @dump`).
    pub fn to_json(&self) -> Result<String, RunError> {
        serde_json::to_string_pretty(self).map_err(|e| RunError::Serde(e.to_string()))
    }

    pub fn from_json(text: &str) -> Result<Self, RunError> {
        serde_json::from_str(text).map_err(|e| RunError::Serde(e.to_string()))
    }

    /// Set one field by its dotted path (`execution.parallel`), for the fields no
    /// `.imd` block carries (`md @recipe_set`). `value` is read as JSON, or as a string when it
    /// is not JSON; a path or value the recipe does not accept is an error.
    pub fn set(&mut self, key: &str, value: &str) -> Result<(), RunError> {
        let err = |m: String| RunError::Recipe(format!("recipe field '{key}': {m}"));
        let mut root = serde_json::to_value(&*self).map_err(|e| err(e.to_string()))?;
        let mut slot = &mut root;
        let mut parts = key.split('.').peekable();
        while let Some(part) = parts.next() {
            let serde_json::Value::Object(table) = slot else {
                return Err(err(format!("'{part}' is not inside a group")));
            };
            slot = if parts.peek().is_some() {
                table
                    .get_mut(part)
                    .ok_or_else(|| err(format!("no group '{part}'")))?
            } else {
                table.entry(part).or_insert(serde_json::Value::Null)
            };
        }
        *slot = serde_json::from_str(value)
            .unwrap_or_else(|_| serde_json::Value::String(value.to_string()));
        *self = serde_json::from_value(root).map_err(|e| err(e.to_string()))?;
        Ok(())
    }

    // --- Factories: the same runs `ImdParameters::{nve,nvt,npt,steepest_descent}` build ---

    pub fn nve(dt: f64, steps: usize, constraints: &str) -> Result<Self, RunError> {
        let imd = ImdParameters::nve(dt, steps, constraints).map_err(RunError::Recipe)?;
        Ok(Self::from_imd(&imd)?.0)
    }

    pub fn nvt(
        dt: f64,
        steps: usize,
        temperature: f64,
        constraints: &str,
    ) -> Result<Self, RunError> {
        let imd =
            ImdParameters::nvt(dt, steps, temperature, constraints).map_err(RunError::Recipe)?;
        Ok(Self::from_imd(&imd)?.0)
    }

    pub fn npt(
        dt: f64,
        steps: usize,
        temperature: f64,
        pressure: f64,
        constraints: &str,
    ) -> Result<Self, RunError> {
        let imd = ImdParameters::npt(dt, steps, temperature, pressure, constraints)
            .map_err(RunError::Recipe)?;
        Ok(Self::from_imd(&imd)?.0)
    }

    pub fn minimize(steps: usize) -> Self {
        Self::from_imd(&ImdParameters::steepest_descent(steps))
            .expect("factory parameters always map")
            .0
    }
}
