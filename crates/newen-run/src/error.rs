//! The one error type of the run assembly.
//!
//! Every failure the `md` binary used to report with `eprintln!` + `process::exit(1)` is a
//! variant here, so the binary and the Python binding render the same condition the same way
//! (3.9 C10). Python maps variants onto builtin exception types in `pyo3-newen`.

use std::fmt;

use newen_core::validation::ValidationReport;
use newen_io::IoError;

pub enum RunError {
    /// Reading an input file failed.
    Io { what: String, source: IoError },
    /// The parameters ask for an input file that was not given
    /// (e.g. `NTPOR=1` without a position-restraint specification).
    MissingInput { flag: &'static str, reason: String },
    /// Topology and coordinate file disagree on the number of atoms.
    AtomCountMismatch { topology: usize, coordinates: usize },
    /// The coordinate file's non-solute atoms are not a whole number of solvent molecules.
    SolventCount {
        coordinates: usize,
        solute: usize,
        atoms_per_solvent: usize,
    },
    /// A validation stage found fatal problems; the report carries them.
    Validation {
        stage: &'static str,
        report: ValidationReport,
    },
    /// Inputs contradict each other (e.g. a perturbed topology with `NTG=0`).
    Inconsistent(String),
    /// The parameter file carries a block the recipe does not model and the caller did not
    /// allow to pass through (3.9 A17: a physics-bearing block must never be ignored).
    UnknownBlock { block: String },
    /// A value the recipe cannot represent (e.g. an unknown NTB code, several temperature baths).
    Recipe(String),
    /// A plan violates the GROMOS step-order invariants (`validate_plan`).
    InvalidPlan(String),
    /// A term/algorithm needs a feature this build lacks (e.g. `schnet` without `--features ml`).
    MissingFeature { term: String, feature: &'static str },
    /// Recipe or plan text could not be (de)serialised.
    Serde(String),
    /// `AlgorithmSequence::init` failed.
    Init(String),
    /// `AlgorithmSequence::run_step` failed at `step`.
    Step { step: usize, message: String },
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RunError::Io { what, source } => write!(f, "error reading {what}: {source}"),
            RunError::MissingInput { flag, reason } => {
                write!(f, "{reason} (no @{flag} file specified)")
            },
            RunError::AtomCountMismatch {
                topology,
                coordinates,
            } => write!(
                f,
                "number of atoms in topology ({topology}) != coordinates ({coordinates})"
            ),
            RunError::SolventCount {
                coordinates,
                solute,
                atoms_per_solvent,
            } => write!(
                f,
                "({coordinates} coords - {solute} solute) is not divisible by \
                 {atoms_per_solvent} atoms/solvent"
            ),
            RunError::Validation { stage, report } => write!(
                f,
                "fatal errors in {stage} validation ({} finding(s)) - cannot continue",
                report.errors.len()
            ),
            RunError::Inconsistent(msg) => write!(f, "inconsistent inputs: {msg}"),
            RunError::UnknownBlock { block } => write!(
                f,
                "input block {block} is not modelled by newen and was not allowed to pass \
                 through — it would be silently ignored (allow it explicitly if it is inert)"
            ),
            RunError::Recipe(msg) => write!(f, "cannot represent the run as a recipe: {msg}"),
            RunError::InvalidPlan(msg) => write!(f, "invalid algorithm plan: {msg}"),
            RunError::MissingFeature { term, feature } => write!(
                f,
                "{term} requires the `{feature}` feature, which this build does not include"
            ),
            RunError::Serde(msg) => write!(f, "recipe/plan (de)serialisation failed: {msg}"),
            RunError::Init(msg) => write!(f, "error initializing algorithm sequence: {msg}"),
            RunError::Step { step, message } => write!(f, "error at step {step}: {message}"),
        }
    }
}

impl fmt::Debug for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RunError({self})")
    }
}

impl std::error::Error for RunError {}
