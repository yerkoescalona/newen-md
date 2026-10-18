//! ext_ti_ana — extrapolate ⟨∂H/∂λ⟩ from one simulation to other λ values and other
//! λ-dependences (gromos++ `ext_ti_ana`, `programs/ext_ti_ana.cc`).
//!
//! gromosXX's `PRECALCLAM` writes, per frame, the energies (`.tre`, `PRECALCLAM`) and their
//! derivatives (`.trg`, `FREEPRECALCLAM`) of every perturbed term at `nrlambdas` grid values of
//! λ, plus the A and B dihedral energies (`ABDIH`, `FREEABDIH`). From these the energy and
//! ∂H/∂λ at any λ (and for any polynomial λ-dependence `@lam<X>`) is an interpolation, and
//! the ⟨∂H/∂λ⟩ a simulation at `@slam` predicts for it a reweighted average
//! `Σ dH·exp(−ΔH/kT) / Σ exp(−ΔH/kT)` (gromos++ `Stat::lnXexpave`). One `predicted_l<slam>.dat`
//! per softness combination; `@bar_data`/`@dHdl_data` also write the unweighted energies and
//! derivatives per frame for `bar` and `mbar`.
//!
//! Deliberate differences from gromos++ (item 1.8):
//! - The dihedral term of the predicted ∂H/∂λ (D5) carries the factor dΛ_dih/dλ like every
//!   other term. The oracle is gromos++ 9cd9f6b3 patched with `D5_ext_ti_ana_dlam.patch`
//!   (`crates/newen-md/tests/gromosXX_references/patches/`).
//! - `@imd` goes through newen's IMD reader, which knows gromosXX's `LAMBDAS` rows by name
//!   as well as by number and maps NTLI 10/11 as gromosXX does (distancefield/dihedralres;
//!   gromos++'s own type table swaps them, `ext_ti_ana.cc:404` used at `:1371-1385`; its
//!   mk_script reader ignores named rows and uses the rows even when NTIL is off). mk_script's
//!   parser diagnostics on stdout are not reproduced.
//! - The `PRECALCLAM` columns follow the trajectory's layout: with 12 columns (md++ 2023-04-15,
//!   `out_configuration.cc:3230-3260`) disres/dihres/disfld are columns 8/10/11 (gromos++ reads
//!   9 and 10, `ext_ti_ana.cc:781-782`, which are angres and dihres there); with 11 they are 8/9/10.
//! - With `@pmin` > 0, `@bar_data`/`@dHdl_data` values are stored at the λ they name
//!   (gromos++ keys by absolute index but looks up `p − pmin`, `ext_ti_ana.cc:513, 865`).
//! - `@dHdl_data` files carry one row per frame without `@bar_data` too (gromos++ loops over the
//!   `@bar_data` energies and writes none, `ext_ti_ana.cc:1025`).
//! - Values gromos++ leaves undefined (a parameter neither in `@imd` nor given, `@pmax` past the
//!   grid, a grid of one λ, `@bar_data` outside `@pmin`/`@pmax`) are errors; the bootstrap seed
//!   can be given (library), gromos++ takes `time(NULL)`; an `@slj_scrf_file` whose last `END`
//!   has no newline is read (gromos++ drops that block).

use super::bar::BarFile;
use super::common::{
    cpow, draw_index, get_value, get_values, population_variance, wall_clock_seed,
};
use crate::args::Arguments;
use crate::energy_sequence::EnergySequence;
use crate::lnexp::{ln_x_exp_ave, lnexp_ave};
use newen_core::lambdas::InteractionLambda;
use newen_io::args::{Program, Status};
use newen_io::energy_traj::{EnergyTraj, OFFICIAL_LIBRARY};
use newen_io::table::{cpp_fixed, cpp_g};
use std::io::Write;

pub const USAGE: &str = "# ext_ti_ana
\t@en_files       <energy files>
\t@fr_files       <free energy files>
\t@library        <library for block information>
\t@temp           <simulation temperature>
\t@nrlambdas      <no. of precalculated lambdas>
\t@minlam         <minimum precalculated lambda>
\t@maxlam         <maximum precalculated lambda>
\t@slam           <lambda value of simulation>
\t@NLAMs          <lambda exponent of simulation>
\t[@lam<X>_sim    <a b c d e; coefficients for individual lambda dependence, default: 0 0 0 1 0>]
\t                (where <X> is one of: slj, scrf, lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)
\t[@imd           <gromos input parameter file>]
\t                (flags temp, nrlambdas, minlam, maxlam, slam, NLAMs and lamX_sim can be
\t                skipped if they can be read from the provided imd file, from PRECALCLAM,
\t                PERTURBATION and MULTIBATH blocks.)
\t[@NLAMp         <lambda exponent value to predict for, default: @NLAMs >]
\t[@lam<X>          <a b c d e, coefficients for individual lambda dependence, default: @lam<X>_sim>]
\t                (where <X> is one of: slj, scrf, lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)
\t[@slj_scrf_file < file with sets of slj and scrf lambda coefficients>]
\t[@no_<X>        <exclude <X> free energy derivative contribution>]
\t                (where <X> is one of: lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)
\t[@countframes   <count nr of contributing frames for each plam>]
\t[@pmin          <min index of prediction>]
\t[@pmax          <max index of prediction>]
\t[@bootstrap     <no. of bootstraps>]
\t[@outdir        <directory to write output to>]
\t[@cpus           <number of threads> Default: 1]
\t[@lam_precision <lambda value precision in outfiles, default: 2>]
\t[@verbose       <print used parameters to file header>]
\t[@bar_data      <print energies to be used for BAR (not reweighted)>]
\t[@dHdl_data     <print dHdl values to be used for MBAR (not reweighted)>]";

pub const KNOWN: &[&str] = &[
    "en_files",
    "fr_files",
    "library",
    "nrlambdas",
    "minlam",
    "maxlam",
    "slam",
    "NLAMs",
    "NLAMp",
    "temp",
    "slj_scrf_file",
    "outdir",
    "imd",
    "lamslj",
    "lamlj",
    "lamscrf",
    "lamcrf",
    "lamkin",
    "lambond",
    "lamang",
    "lamimpr",
    "lamdih",
    "lamdisres",
    "lamdihres",
    "lamdisfld",
    "lamslj_sim",
    "lamlj_sim",
    "lamscrf_sim",
    "lamcrf_sim",
    "lamkin_sim",
    "lambond_sim",
    "lamang_sim",
    "lamimpr_sim",
    "lamdih_sim",
    "lamdisres_sim",
    "lamdihres_sim",
    "lamdisfld_sim",
    "no_lj",
    "no_crf",
    "no_kin",
    "no_bond",
    "no_ang",
    "no_dih",
    "no_disres",
    "no_dihres",
    "no_disfld",
    "no_impr",
    "pmin",
    "pmax",
    "bootstrap",
    "countframes",
    "lam_precision",
    "verbose",
    "cpus",
    "bar_data",
    "dHdl_data",
];

pub const PROGRAM: Program = Program {
    name: "ext_ti_ana",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// The λ-dependent terms, in gromos++'s order (`ext_ti_ana.cc:404`): the index of a type here
/// is the index into [`Options::coefficients`].
pub const TYPES: [&str; 12] = [
    "kin", "bond", "ang", "dih", "impr", "lj", "slj", "crf", "scrf", "disres", "dihres", "disfld",
];
const KIN: usize = 0;
const BOND: usize = 1;
const ANG: usize = 2;
const DIH: usize = 3;
const IMPR: usize = 4;
const LJ: usize = 5;
const SLJ: usize = 6;
const CRF: usize = 7;
const SCRF: usize = 8;
const DISRES: usize = 9;
const DIHRES: usize = 10;
const DISFLD: usize = 11;

/// The identity polynomial, gromos++'s default for every `@lam<X>_sim`.
pub const IDENTITY: [f64; 5] = [0.0, 0.0, 0.0, 1.0, 0.0];

/// The `@no_<X>` switches: each drops that term's contribution to the predicted ∂H/∂λ (only
/// the derivative; the energies used for reweighting keep it, as in gromos++).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Exclude {
    pub lj: bool,
    pub crf: bool,
    pub kin: bool,
    pub bond: bool,
    pub ang: bool,
    pub impr: bool,
    pub dih: bool,
    pub disres: bool,
    pub dihres: bool,
    pub disfld: bool,
}

/// Everything the program's arguments say, resolved (the `@imd` values already merged in).
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    pub en_files: Vec<String>,
    pub fr_files: Vec<String>,
    /// The `ene_ana` library file; `None` is newen-io's built-in one (library callers only:
    /// the program requires `@library`, as gromos++ does).
    pub library: Option<String>,
    pub temp: f64,
    pub nrlambdas: usize,
    pub minlam: f64,
    pub maxlam: f64,
    pub slam: f64,
    pub nlams: f64,
    pub nlamp: f64,
    /// `@lam<X>_sim`, per [`TYPES`] entry: the λ-dependence the simulation used.
    pub coefficients_sim: [[f64; 5]; 12],
    /// `@lam<X>`: the λ-dependence to predict for.
    pub coefficients: [[f64; 5]; 12],
    pub slj_scrf_file: Option<String>,
    pub exclude: Exclude,
    pub countframes: bool,
    pub pmin: usize,
    pub pmax: usize,
    pub bootstrap: usize,
    /// The bootstrap seed base (the stream of SLJ set `i` starts at `seed + i`); `None` is the
    /// wall clock, as gromos++ does.
    pub seed: Option<u32>,
    pub outdir: String,
    pub lam_precision: usize,
    pub verbose: bool,
    /// `@bar_data`: `None` off, `Some(empty)` every predicted λ, else those λ (sorted).
    pub bar_data: Option<Vec<f64>>,
    /// `@dHdl_data`, the same way.
    pub dhdl_data: Option<Vec<f64>>,
    pub cpus: u32,
}

/// One interpolated λ: its value and derivative from the polynomial, and the two grid points
/// and weights that interpolate a precalculated table at it (`LambdaStruct`).
#[derive(Debug, Clone, Copy)]
struct Lam {
    lam: f64,
    dlam: f64,
    lo: usize,
    hi: usize,
    wlow: f64,
    whigh: f64,
}

impl Lam {
    #[inline]
    fn at(&self, table: &[Vec<f64>], col: usize) -> f64 {
        self.wlow * table[self.lo][col] + self.whigh * table[self.hi][col]
    }
}

/// `[ a b c d e ]` in the C++ default float format (`ext_ti_ana.cc:221-231`).
fn coeff_vec(c: &[f64]) -> String {
    let mut s = String::from("[");
    for &v in c {
        s.push(' ');
        s.push_str(&g6(v));
    }
    s.push_str(" ]");
    s
}

fn g6(x: f64) -> String {
    if x == 0.0 && x.is_sign_negative() {
        return "-0".into();
    }
    cpp_g(x, 6)
}

/// The λ grid of the precalculated values: `p·step + minlam` (`ext_ti_ana.cc:487-494`).
struct Grid {
    minlam: f64,
    step: f64,
    p_lam: Vec<f64>,
}

impl Grid {
    fn new(n: usize, minlam: f64, maxlam: f64) -> Grid {
        let step = (maxlam - minlam) / (n - 1) as f64;
        Grid {
            minlam,
            step,
            p_lam: (0..n).map(|p| p as f64 * step + minlam).collect(),
        }
    }

    /// gromos++ `getSingleLambda` (`ext_ti_ana.cc:1123-1166`). `checked` is false for lj, crf
    /// and dih, whose λ only enters as a factor and may leave the grid.
    fn single(&self, c: &[f64; 5], slam: f64, checked: bool) -> Result<Lam, String> {
        let lam = c[0] * cpow(slam, 4.0)
            + c[1] * cpow(slam, 3.0)
            + c[2] * cpow(slam, 2.0)
            + c[3] * slam
            + c[4];
        let dlam =
            4.0 * c[0] * cpow(slam, 3.0) + 3.0 * c[1] * cpow(slam, 2.0) + 2.0 * c[2] * slam + c[3];
        let n = self.p_lam.len();
        let maxlam = self.minlam + (n - 1) as f64 * self.step;
        if checked && (lam > maxlam || lam < self.minlam) {
            return Err(format!(
                "at lambda={}: with lambda coefficients {}\n            we get final lambda \
                 values outside the precalculated\n            range [{}-{}], interpolation not \
                 possible.",
                g6(slam),
                coeff_vec(c),
                g6(self.minlam),
                g6(maxlam)
            ));
        }
        let q = (lam - self.minlam) / self.step;
        let (lo, hi) = (q.floor(), q.ceil());
        let in_grid = lo >= 0.0 && hi < n as f64;
        if !in_grid {
            if checked {
                return Err(format!(
                    "at lambda={}: lambda {} rounds outside the precalculated grid",
                    g6(slam),
                    g6(lam)
                ));
            }
            // lj, crf and dih never read the table through these
            return Ok(Lam {
                lam,
                dlam,
                lo: 0,
                hi: 0,
                wlow: 0.0,
                whigh: 1.0,
            });
        }
        let (lo, hi) = (lo as usize, hi as usize);
        let (wlow, whigh) = if lo == hi {
            (0.0, 1.0)
        } else {
            let p = &self.p_lam;
            (
                (lam - p[hi]) / (p[lo] - p[hi]),
                (lam - p[lo]) / (p[hi] - p[lo]),
            )
        };
        Ok(Lam {
            lam,
            dlam,
            lo,
            hi,
            wlow,
            whigh,
        })
    }

    /// gromos++ `getLambdas`: [`Grid::single`] at every grid λ.
    fn all(&self, c: &[f64; 5], checked: bool) -> Result<Vec<Lam>, String> {
        self.p_lam
            .iter()
            .map(|&l| self.single(c, l, checked))
            .collect()
    }
}

/// The unweighted per-frame energies `@bar_data` writes: `e_s` at the simulated λ and
/// `e_p[k]` at `plam[k]` (kJ/mol).
#[derive(Debug, Clone, PartialEq)]
pub struct BarData {
    pub plam: Vec<f64>,
    pub e_s: Vec<f64>,
    pub e_p: Vec<Vec<f64>>,
}

/// The unweighted per-frame ∂H/∂λ `@dHdl_data` writes, `de[k]` at `plam[k]` (kJ/mol).
#[derive(Debug, Clone, PartialEq)]
pub struct DhdlData {
    pub plam: Vec<f64>,
    pub de: Vec<Vec<f64>>,
}

/// The predictions for one combination of softness polynomials (the single `@lamslj`/`@lamscrf`
/// one, or one SLJ × SCRF pair of `@slj_scrf_file`).
#[derive(Debug, Clone, PartialEq)]
pub struct Prediction {
    /// The `@slj_scrf_file` labels, `None` without the file.
    pub labels: Option<(i64, i64)>,
    pub slj: [f64; 5],
    pub scrf: [f64; 5],
    /// The predicted λ (`p_LAM[pmin..pmax]`).
    pub plam: Vec<f64>,
    /// ⟨∂H/∂λ⟩ predicted at each `plam` (kJ/mol).
    pub dhdl: Vec<f64>,
    /// The bootstrap standard deviation (`@bootstrap`).
    pub error: Option<Vec<f64>>,
    /// Frames with ΔH ≤ ΔG + 0.001 kJ/mol (`@countframes`).
    pub contributing_frames: Option<Vec<usize>>,
    pub bar: Option<BarData>,
    pub dhdl_data: Option<DhdlData>,
}

/// Everything one run computes.
#[derive(Debug, Clone, PartialEq)]
pub struct Analysis {
    pub n_frames: usize,
    pub predictions: Vec<Prediction>,
}

/// The program: parse, compute, write the files. Nothing goes to `out`.
pub fn run(argv: &[String], _out: &mut dyn Write) -> Result<(), String> {
    let opts = Options::from_argv(argv)?;
    let analysis = analyse(&opts)?;
    write_files(&opts, &analysis)?;
    Ok(())
}

/// [`run`] that also returns what it computed and the files it wrote.
pub fn evaluate(argv: &[String]) -> Result<(Options, Analysis, Vec<String>), String> {
    let opts = Options::from_argv(argv)?;
    let analysis = analyse(&opts)?;
    let files = write_files(&opts, &analysis)?;
    Ok((opts, analysis, files))
}

impl Options {
    /// gromos++'s argument handling (`ext_ti_ana.cc:313-611`), `@imd` included.
    pub fn from_argv(argv: &[String]) -> Result<Options, String> {
        let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
        Self::from_args(&args, true)
    }

    /// The options of a call on frames held in memory ([`analyse_frames`], [`bar_file`]): a
    /// simulation at `slam` with `NLAM = nlam` and the identity λ-dependence of every term,
    /// predicting for the same; `@bar_data` at `bar_data`; every other value gromos++'s
    /// default. No files are named.
    pub fn for_frames(
        temp: f64,
        nrlambdas: usize,
        minlam: f64,
        maxlam: f64,
        slam: f64,
        nlam: f64,
        mut bar_data: Vec<f64>,
    ) -> Result<Options, String> {
        if nrlambdas < 2 {
            return Err(format!(
                "nrlambdas: {nrlambdas} precalculated lambdas leave no grid to interpolate on"
            ));
        }
        bar_data.sort_by(|a, b| a.total_cmp(b));
        Ok(Options {
            en_files: Vec::new(),
            fr_files: Vec::new(),
            library: None,
            temp,
            nrlambdas,
            minlam,
            maxlam,
            slam,
            nlams: nlam,
            nlamp: nlam,
            coefficients_sim: [IDENTITY; 12],
            coefficients: [IDENTITY; 12],
            slj_scrf_file: None,
            exclude: Exclude::default(),
            countframes: false,
            pmin: 0,
            pmax: nrlambdas,
            bootstrap: 0,
            seed: None,
            outdir: ".".into(),
            lam_precision: 2,
            verbose: false,
            bar_data: Some(bar_data),
            dhdl_data: None,
            cpus: 1,
        })
    }

    /// `library_required = false` lets a library caller leave `@library` out for newen-io's
    /// built-in one ([`Options::library`] `None`).
    pub fn from_args(args: &Arguments, library_required: bool) -> Result<Options, String> {
        if args.count("en_files") <= 0 || args.count("fr_files") <= 0 {
            return Err(format!(
                "en_files and fr_files need to be specified:\n{USAGE}"
            ));
        }
        if library_required && args.count("library") <= 0 {
            return Err(format!("no library file specified:\n{USAGE}"));
        }
        let library = args.values("library").first().cloned();
        let bootstrap: i64 = get_value(args, "bootstrap", Some(0))?;
        let lam_precision: i64 = get_value(args, "lam_precision", Some(2))?;
        let cpus: i64 = get_value(args, "cpus", Some(1))?;
        let present = |k: &str| args.has(k);
        let exclude = Exclude {
            lj: present("no_lj"),
            crf: present("no_crf"),
            kin: present("no_kin"),
            bond: present("no_bond"),
            ang: present("no_ang"),
            impr: present("no_impr"),
            dih: present("no_dih"),
            disres: present("no_disres"),
            dihres: present("no_dihres"),
            disfld: present("no_disfld"),
        };
        let countframes = args.count("countframes") == 0;
        let outdir: String = get_value(args, "outdir", Some(".".to_string()))?;

        let mut coefficients_sim = [IDENTITY; 12];
        for (t, name) in TYPES.iter().enumerate() {
            coefficients_sim[t] = to5(get_values(
                args,
                &format!("lam{name}_sim"),
                5,
                Some(IDENTITY.to_vec()),
            )?);
        }

        // values the imd file supplies; a flag given as well overrides it
        let (mut temp, mut nrlambdas, mut minlam, mut maxlam, mut slam, mut nlams) =
            (None, None, None, None, None, None);
        if args.has("imd") {
            let path: String = get_value(args, "imd", None)?;
            if !std::path::Path::new(&path).is_file() {
                return Err(format!("@imd: file does not exist: {path}\n"));
            }
            let imd = newen_run::read_imd(&path).map_err(|e| format!("@imd: {e}"))?;
            if imd.raw_blocks.contains_key("PERTURBATION") {
                slam = Some(imd.rlam);
                nlams = Some(imd.nlam as f64);
            }
            if imd.raw_blocks.contains_key("MULTIBATH") {
                temp = imd.temp_bath.first().and_then(|b| b.temp0.first().copied());
            }
            if let Some(pre) = imd.precalclam().map_err(|e| format!("@imd: {e}"))? {
                minlam = Some(pre.minlam);
                maxlam = Some(pre.maxlam);
                nrlambdas = Some(pre.nrlam as i64);
            }
            if let Some(block) = &imd.lambdas {
                eprintln!(
                    "# WARNING: LAMBDAS block: if there is more than one set of coefficients\n\
                     #      for the same property (different lambda dependence for different\n\
                     #      energy groups), ext_ti_ana will use only the first set!"
                );
                if block.applies(imd.ntg) {
                    // the first row of a type wins (gromos++ `read_lambdas`)
                    let mut first: [Option<[f64; 5]>; 12] = [None; 12];
                    let mut warn = false;
                    for row in &block.rows {
                        let Some(t) = lambdas_type(row.kind) else {
                            continue;
                        };
                        match first[t] {
                            Some(c) => warn |= c != row.coefficients,
                            None => {
                                first[t] = Some(row.coefficients);
                                if !args.has(&format!("lam{}_sim", TYPES[t])) {
                                    coefficients_sim[t] = row.coefficients;
                                }
                            },
                        }
                    }
                    if warn {
                        eprintln!(
                            "# WARNING: LAMBDAS block: different coefficients for different \
                             energy groups"
                        );
                    }
                }
            }
            eprintln!("# Read parameters from: {path}");
            let flags = ["slam", "NLAMs", "temp", "nrlambdas", "minlam", "maxlam"];
            let given: Vec<String> = flags
                .iter()
                .filter(|f| args.has(f))
                .map(|f| format!(" @{f}"))
                .collect();
            if !given.is_empty() {
                eprintln!(
                    "# WARNING: {} flag(s) overwrite values from the imd file",
                    given.concat()
                );
            }
        }
        let need = |key: &str, from_imd: bool| -> String {
            if from_imd {
                format!("@{key}: not given and not in the @imd file")
            } else {
                format!("Argument '{key}' is required.")
            }
        };
        let imd = args.has("imd");
        let temp: f64 = get_value(args, "temp", temp).map_err(|e| or_need(e, need("temp", imd)))?;
        let nrlambdas: i64 = get_value(args, "nrlambdas", nrlambdas)
            .map_err(|e| or_need(e, need("nrlambdas", imd)))?;
        let minlam: f64 =
            get_value(args, "minlam", minlam).map_err(|e| or_need(e, need("minlam", imd)))?;
        let maxlam: f64 =
            get_value(args, "maxlam", maxlam).map_err(|e| or_need(e, need("maxlam", imd)))?;
        let slam: f64 = get_value(args, "slam", slam).map_err(|e| or_need(e, need("slam", imd)))?;
        let nlams: f64 =
            get_value(args, "NLAMs", nlams).map_err(|e| or_need(e, need("NLAMs", imd)))?;
        let nlamp: f64 = get_value(args, "NLAMp", Some(nlams))?;
        if nrlambdas < 2 {
            return Err(format!(
                "@nrlambdas: {nrlambdas} precalculated lambdas leave no grid to interpolate on"
            ));
        }
        let nrlambdas = nrlambdas as usize;

        let list = |key: &str| -> Result<Option<Vec<f64>>, String> {
            if !args.has(key) {
                return Ok(None);
            }
            let mut v: Vec<f64> = get_values(args, key, args.values(key).len(), None)?;
            v.sort_by(|a, b| a.total_cmp(b));
            Ok(Some(v))
        };
        let bar_data = list("bar_data")?;
        let dhdl_data = list("dHdl_data")?;

        let mut coefficients = coefficients_sim;
        for (t, name) in TYPES.iter().enumerate() {
            coefficients[t] = to5(get_values(
                args,
                &format!("lam{name}"),
                5,
                Some(coefficients_sim[t].to_vec()),
            )?);
        }
        let slj_scrf_file = if args.count("slj_scrf_file") > 0 {
            if args.has("lamslj") || args.has("lamscrf") {
                return Err("specify slj_scrf_file OR lamslj/lamscrf".into());
            }
            Some(args.values("slj_scrf_file")[0].clone())
        } else {
            None
        };
        let pmin: i64 = get_value(args, "pmin", Some(0))?;
        let pmax: i64 = get_value(args, "pmax", Some(nrlambdas as i64))?;
        if pmin < 0 || pmax < pmin || pmax > nrlambdas as i64 {
            return Err(format!(
                "@pmin {pmin} / @pmax {pmax}: need 0 <= pmin <= pmax <= nrlambdas ({nrlambdas})"
            ));
        }
        if bootstrap < 0 || lam_precision < 0 || cpus < 0 {
            return Err("@bootstrap, @lam_precision and @cpus cannot be negative".into());
        }
        Ok(Options {
            en_files: args.values("en_files").to_vec(),
            fr_files: args.values("fr_files").to_vec(),
            library,
            temp,
            nrlambdas,
            minlam,
            maxlam,
            slam,
            nlams,
            nlamp,
            coefficients_sim,
            coefficients,
            slj_scrf_file,
            exclude,
            countframes,
            pmin: pmin as usize,
            pmax: pmax as usize,
            bootstrap: bootstrap as usize,
            seed: None,
            outdir,
            lam_precision: lam_precision as usize,
            verbose: args.has("verbose"),
            bar_data,
            dhdl_data,
            cpus: cpus as u32,
        })
    }
}

/// A "required" error from `get_value` reworded when the value could have come from `@imd`.
fn or_need(e: String, need: String) -> String {
    if e.contains("is required") {
        need
    } else {
        e
    }
}

/// The `@bar_data`/`@dHdl_data` column grid point `p` feeds: every predicted λ in order when
/// the flag has no values, else the listed λ that matched `p`.
#[inline]
fn column(
    list: &Option<Vec<f64>>,
    index: &[Option<usize>],
    p: usize,
    pmin: usize,
) -> Option<usize> {
    match list {
        None => None,
        Some(v) if v.is_empty() => Some(p - pmin),
        Some(_) => index[p],
    }
}

/// The ext_ti_ana type a gromosXX `LAMBDAS` kind sets (gromosXX's own numbering: NTLI 10 is
/// distancefield, 11 dihedralres, `in_parameter.cc:4339-4341`). Mass and the restraint kinds
/// ext_ti_ana does not model set nothing.
fn lambdas_type(kind: InteractionLambda) -> Option<usize> {
    Some(match kind {
        InteractionLambda::Bond => BOND,
        InteractionLambda::Angle => ANG,
        InteractionLambda::Dihedral => DIH,
        InteractionLambda::Improper => IMPR,
        InteractionLambda::Lj => LJ,
        InteractionLambda::LjSoftness => SLJ,
        InteractionLambda::Crf => CRF,
        InteractionLambda::CrfSoftness => SCRF,
        InteractionLambda::Disres => DISRES,
        InteractionLambda::Disfield => DISFLD,
        InteractionLambda::Dihres => DIHRES,
        InteractionLambda::Mass | InteractionLambda::Angres | InteractionLambda::Colvarres => {
            return None
        },
    })
}

fn to5(v: Vec<f64>) -> [f64; 5] {
    [v[0], v[1], v[2], v[3], v[4]]
}

/// The SLJ and SCRF sets of an `@slj_scrf_file` (`read_slj_scrf_file`, `ext_ti_ana.cc:1185-1245`):
/// a TITLE, then `SLJ`/`SCRF` blocks of `label a b c d e` lines.
/// Labelled softness polynomials: `(label, [a, b, c, d, e])`.
type SoftnessSets = Vec<(i64, [f64; 5])>;

fn read_slj_scrf(path: &str) -> Result<(SoftnessSets, SoftnessSets), String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("failed to open slj_scrf file {path}: {e}"))?;
    let mut blocks: Vec<Vec<String>> = Vec::new();
    let mut cur: Option<Vec<String>> = None;
    for raw in text.lines() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        match cur.as_mut() {
            None => cur = Some(vec![line.to_string()]),
            Some(b) => {
                b.push(line.to_string());
                if line.starts_with("END") {
                    blocks.push(cur.take().unwrap());
                }
            },
        }
    }
    if let Some(b) = cur {
        return Err(format!(
            "file {path} is corrupted. No END in {} block.",
            b[0]
        ));
    }
    if blocks.first().map(|b| b[0].as_str()) != Some("TITLE") {
        return Err(format!("{path}: TITLE block expected"));
    }
    let (mut slj, mut scrf) = (Vec::new(), Vec::new());
    for b in &blocks[1..] {
        let target = match b[0].as_str() {
            "SLJ" => &mut slj,
            "SCRF" => &mut scrf,
            other => {
                return Err(format!(
                    "read unknown block:{other}\n known are: SLJ, SCRF\n"
                ))
            },
        };
        for line in &b[1..b.len() - 1] {
            let tok: Vec<&str> = line.split_whitespace().collect();
            let label = tok.first().and_then(|t| t.parse::<i64>().ok());
            let c: Option<Vec<f64>> = tok.get(1..6).map(|t| {
                t.iter()
                    .filter_map(|v| super::common::parse_finite(v))
                    .collect()
            });
            match (label, c) {
                (Some(l), Some(c)) if c.len() == 5 => target.push((l, to5(c))),
                _ => {
                    return Err(format!(
                        "{path}: block {}: expected 'label a b c d e', got '{line}'",
                        b[0]
                    ))
                },
            }
        }
    }
    Ok((slj, scrf))
}

/// Read the trajectories and compute every prediction (`ext_ti_ana.cc:540-1092`), without
/// writing anything.
pub fn analyse(o: &Options) -> Result<Analysis, String> {
    analyse_source(o, None)
}

/// [`analyse`] over frames held in memory instead of `@en_files`/`@fr_files` (`o.en_files`,
/// `o.fr_files` and `o.library` are not read): what a run that kept its PRECALCLAM energies,
/// e.g. a Python `Simulation`, gives the program.
pub fn analyse_frames(o: &Options, frames: &[PrecalcFrame]) -> Result<Analysis, String> {
    analyse_source(o, Some(frames))
}

/// One frame's extended-TI data: the `PRECALCLAM` rows of the energy trajectory and the
/// `FREEPRECALCLAM` rows of the free-energy trajectory (one per grid λ, 11 or 12 columns), and
/// their `ABDIH`/`FREEABDIH` rows (A, B).
#[derive(Debug, Clone, PartialEq)]
pub struct PrecalcFrame {
    pub energies: Vec<Vec<f64>>,
    pub derivatives: Vec<Vec<f64>>,
    pub abdih: [f64; 2],
    pub free_abdih: [f64; 2],
}

/// The tables of one frame, borrowed from where they are kept.
struct FrameTables<'a> {
    se: &'a [Vec<f64>],
    sfe: &'a [Vec<f64>],
    se_dih: [f64; 2],
    sfe_dih: [f64; 2],
}

/// The `@en_files`/`@fr_files` frame stream.
struct FileFrames {
    en: EnergySequence,
    fr: EnergySequence,
    etrj: EnergyTraj,
}

impl FileFrames {
    fn next(&mut self) -> Result<Option<FrameTables<'_>>, String> {
        if !self.en.read_frame(&mut self.etrj)? || !self.fr.read_frame(&mut self.etrj)? {
            return Ok(None);
        }
        let se = self
            .etrj
            .table("PRECALCLAM")
            .ok_or("the library declares no PRECALCLAM subblock")?;
        let sfe = self
            .etrj
            .table("FREEPRECALCLAM")
            .ok_or("the library declares no FREEPRECALCLAM subblock")?;
        Ok(Some(FrameTables {
            se,
            sfe,
            se_dih: abdih(&self.etrj, "ABDIH")?,
            sfe_dih: abdih(&self.etrj, "FREEABDIH")?,
        }))
    }
}

enum Frames<'a> {
    Files(Box<FileFrames>),
    Memory(std::slice::Iter<'a, PrecalcFrame>),
}

impl Frames<'_> {
    fn next(&mut self) -> Result<Option<FrameTables<'_>>, String> {
        match self {
            Frames::Files(f) => f.next(),
            Frames::Memory(it) => Ok(it.next().map(|f| FrameTables {
                se: &f.energies,
                sfe: &f.derivatives,
                se_dih: f.abdih,
                sfe_dih: f.free_abdih,
            })),
        }
    }
}

/// `@bar_data` of one simulation from frames held in memory: `o.bar_data`'s λ, the energies
/// [`analyse_frames`] gives them, as the [`BarFile`] `bar` would read from
/// `bar_data_l<slam>.dat` (at full precision, the simulated λ `o.slam` unrounded). Every
/// `o.bar_data` λ must be a grid point.
pub fn bar_file(o: &Options, frames: &[PrecalcFrame]) -> Result<BarFile, String> {
    let wanted = o.bar_data.as_deref().unwrap_or_default();
    let grid = Grid::new(o.nrlambdas, o.minlam, o.maxlam);
    let eps = lam_epsilon(o.lam_precision);
    if let Some(&l) = wanted
        .iter()
        .find(|&&l| !grid.p_lam.iter().any(|&p| (l - p).abs() < eps))
    {
        return Err(format!(
            "lambda {} is not a point of the PRECALCLAM grid ({} points on [{}, {}]): BAR and \
             MBAR need the energies at every simulated lambda",
            g6(l),
            o.nrlambdas,
            g6(o.minlam),
            g6(o.maxlam)
        ));
    }
    let a = analyse_frames(o, frames)?;
    let b = a
        .predictions
        .into_iter()
        .next()
        .and_then(|p| p.bar)
        .ok_or("bar_file: Options::bar_data is not set")?;
    Ok(BarFile {
        slam: o.slam,
        plam: b.plam,
        e_s: b.e_s,
        e_p: b.e_p,
    })
}

/// `10^-(lam_precision + 2)`: how close a `@bar_data`/`@dHdl_data` λ must be to its grid point.
fn lam_epsilon(lam_precision: usize) -> f64 {
    cpow(10.0, -((lam_precision + 2) as f64))
}

fn analyse_source(o: &Options, frames: Option<&[PrecalcFrame]>) -> Result<Analysis, String> {
    let n = o.nrlambdas;
    let grid = Grid::new(n, o.minlam, o.maxlam);
    let lam_epsilon = lam_epsilon(o.lam_precision);
    // @bar_data / @dHdl_data λ → grid index; every value must match one grid point
    let index_of = |vals: &Option<Vec<f64>>, what: &str| -> Result<Vec<Option<usize>>, String> {
        let mut idx = vec![None; n];
        if let Some(v) = vals {
            if v.is_empty() {
                return Ok(idx);
            }
            for (i, &l) in v.iter().enumerate() {
                for (p, &pl) in grid.p_lam.iter().enumerate() {
                    if (l - pl).abs() < lam_epsilon {
                        idx[p] = Some(i);
                    }
                }
            }
            if idx.iter().flatten().count() != v.len() {
                return Err(format!(
                    "not all {what} values are available - check precision?"
                ));
            }
            for i in 0..v.len() {
                let p = idx.iter().position(|x| *x == Some(i));
                if !p.is_some_and(|p| (o.pmin..o.pmax).contains(&p)) {
                    return Err(format!(
                        "{what} value {} is not among the predicted lambdas (@pmin/@pmax)",
                        g6(v[i])
                    ));
                }
            }
        }
        Ok(idx)
    };
    let bar_index = index_of(&o.bar_data, "bar_data")?;
    let dhdl_index = index_of(&o.dhdl_data, "dHdl_data")?;

    // the reference λ at slam and the predicted ones on the grid, type by type
    let mut reference = Vec::with_capacity(12);
    let mut predicted = Vec::with_capacity(12);
    for t in 0..12 {
        let checked = !matches!(t, LJ | CRF | DIH);
        reference.push(grid.single(&o.coefficients_sim[t], o.slam, checked)?);
        predicted.push(grid.all(&o.coefficients[t], checked)?);
    }
    let (slj_sets, scrf_sets, labelled) = match &o.slj_scrf_file {
        Some(path) => {
            eprintln!("# slj_scrf_file specified");
            let (a, b) = read_slj_scrf(path)?;
            (a, b, true)
        },
        None => {
            eprintln!("# single slj and scrf value");
            (
                vec![(0, o.coefficients[SLJ])],
                vec![(0, o.coefficients[SCRF])],
                false,
            )
        },
    };
    let slj: Vec<Vec<Lam>> = slj_sets
        .iter()
        .map(|(_, c)| grid.all(c, true))
        .collect::<Result<_, _>>()?;
    let scrf: Vec<Vec<Lam>> = scrf_sets
        .iter()
        .map(|(_, c)| grid.all(c, true))
        .collect::<Result<_, _>>()?;

    let kbt = crate::gromospp_kb_warn() * o.temp;
    let mut source = match frames {
        Some(f) => Frames::Memory(f.iter()),
        None => {
            let lib_text = match &o.library {
                Some(path) => std::io::read_to_string(
                    newen_io::gz::open_text(path)
                        .map_err(|_| format!("failed to open library file {path}"))?,
                )
                .map_err(|e| format!("{path}: {e}"))?,
                None => OFFICIAL_LIBRARY.to_string(),
            };
            let lib_name = o.library.as_deref().unwrap_or("<built-in library>");
            let mut etrj =
                EnergyTraj::from_library(&lib_text).map_err(|e| format!("{lib_name}: {e}"))?;
            etrj.add_constant("MASS", 0.0);
            Frames::Files(Box::new(FileFrames {
                en: EnergySequence::new(o.en_files.clone(), "ENERTRJ", lib_name),
                fr: EnergySequence::new(o.fr_files.clone(), "FRENERTRJ", lib_name),
                etrj,
            }))
        },
    };

    let nr_plam = o.pmax - o.pmin;
    let (ni, nj) = (slj.len(), scrf.len());
    let combos = ni * nj;
    let nr_bar = match &o.bar_data {
        None => 0,
        Some(v) if v.is_empty() => nr_plam,
        Some(v) => v.len(),
    };
    let nr_dhdl = match &o.dhdl_data {
        None => 0,
        Some(v) if v.is_empty() => nr_plam,
        Some(v) => v.len(),
    };
    let mut x: Vec<Vec<f64>> = vec![Vec::new(); combos * nr_plam];
    let mut vyvr: Vec<Vec<f64>> = vec![Vec::new(); combos * nr_plam];
    let mut e_sim: Vec<f64> = Vec::new();
    let mut e_bar: Vec<Vec<f64>> = vec![Vec::new(); combos * nr_bar];
    let mut e_dhdl: Vec<Vec<f64>> = vec![Vec::new(); combos * nr_dhdl];

    let ex = o.exclude;
    let r = &reference;
    let mut n_frames = 0usize;
    let mut columns: Option<(usize, usize, usize)> = None;
    while let Some(FrameTables {
        se,
        sfe,
        se_dih,
        sfe_dih,
    }) = source.next()?
    {
        n_frames += 1;
        if se.len() != n {
            return Err(format!(
                "Expected {n} precalculated values but found {} in the energy trajectory\n",
                se.len()
            ));
        }
        if sfe.len() != n {
            return Err(format!(
                "Expected {n} precalculated values but found {} in the free energy trajectory\n",
                sfe.len()
            ));
        }
        let (c_disres, c_dihres, c_disfld) = match columns {
            Some(c) => c,
            None => {
                let width = se[0].len().min(sfe[0].len());
                let c = match width {
                    w if w >= 12 => (8, 10, 11),
                    11 => (8, 9, 10),
                    w => {
                        return Err(format!(
                            "PRECALCLAM rows have {w} columns; ext_ti_ana reads 11 (md++ before \
                             2023) or 12"
                        ))
                    },
                };
                columns = Some(c);
                c
            },
        };
        // the reference energy, in gromos++'s term order (ext_ti_ana.cc:762-784)
        let pow1_rlj = cpow(1.0 - r[LJ].lam, o.nlams);
        let pow_rlj = cpow(r[LJ].lam, o.nlams);
        let pow1_rcrf = cpow(1.0 - r[CRF].lam, o.nlams);
        let pow_rcrf = cpow(r[CRF].lam, o.nlams);
        let mut e_s = 0.0;
        e_s += pow1_rlj * r[SLJ].at(se, 0) + pow_rlj * r[SLJ].at(se, 1);
        e_s += pow1_rcrf * r[SCRF].at(se, 2) + pow_rcrf * r[SCRF].at(se, 3);
        e_s += r[KIN].at(se, 4);
        e_s += r[BOND].at(se, 5);
        e_s += r[ANG].at(se, 6);
        e_s += r[IMPR].at(se, 7);
        e_s += r[DISRES].at(se, c_disres);
        e_s += r[DIHRES].at(se, c_dihres);
        e_s += r[DISFLD].at(se, c_disfld);
        e_s += (1.0 - r[DIH].lam) * se_dih[0] + r[DIH].lam * se_dih[1];
        if o.bar_data.is_some() {
            e_sim.push(e_s);
        }

        for p in o.pmin..o.pmax {
            let l = |t: usize| &predicted[t][p];
            let pow1_lj = cpow(1.0 - l(LJ).lam, o.nlamp);
            let pow_lj = cpow(l(LJ).lam, o.nlamp);
            let pow1_crf = cpow(1.0 - l(CRF).lam, o.nlamp);
            let pow_crf = cpow(l(CRF).lam, o.nlamp);

            let mut e = 0.0;
            e += l(KIN).at(se, 4);
            e += l(BOND).at(se, 5);
            e += l(ANG).at(se, 6);
            e += l(IMPR).at(se, 7);
            e += l(DISRES).at(se, c_disres);
            e += l(DIHRES).at(se, c_dihres);
            e += l(DISFLD).at(se, c_disfld);
            e += (1.0 - l(DIH).lam) * se_dih[0] + l(DIH).lam * se_dih[1];

            let mut de = 0.0;
            let dpow1_lj = o.nlamp * cpow(1.0 - l(LJ).lam, o.nlamp - 1.0);
            let dpow_lj = o.nlamp * cpow(l(LJ).lam, o.nlamp - 1.0);
            let dpow1_crf = o.nlamp * cpow(1.0 - l(CRF).lam, o.nlamp - 1.0);
            let dpow_crf = o.nlamp * cpow(l(CRF).lam, o.nlamp - 1.0);
            if !ex.kin {
                de += l(KIN).at(sfe, 4) * l(KIN).dlam;
            }
            if !ex.bond {
                de += l(BOND).at(sfe, 5) * l(BOND).dlam;
            }
            if !ex.ang {
                de += l(ANG).at(sfe, 6) * l(ANG).dlam;
            }
            if !ex.impr {
                de += l(IMPR).at(sfe, 7) * l(IMPR).dlam;
            }
            if !ex.disres {
                de += l(DISRES).at(sfe, c_disres) * l(DISRES).dlam;
            }
            if !ex.dihres {
                de += l(DIHRES).at(sfe, c_dihres) * l(DIHRES).dlam;
            }
            if !ex.disfld {
                de += l(DISFLD).at(sfe, c_disfld) * l(DISFLD).dlam;
            }
            if !ex.dih {
                de += sfe_dih[0] * l(DIH).dlam;
            }

            for (i, slj_i) in slj.iter().enumerate() {
                let s = &slj_i[p];
                for (j, scrf_j) in scrf.iter().enumerate() {
                    let c = &scrf_j[p];
                    let mut e_mult = e;
                    e_mult += pow1_lj * s.at(se, 0) + pow_lj * s.at(se, 1);
                    e_mult += pow1_crf * c.at(se, 2) + pow_crf * c.at(se, 3);
                    let mut de_mult = de;
                    if !ex.lj {
                        de_mult += (-dpow1_lj * s.at(se, 0) + dpow_lj * s.at(se, 1)) * l(LJ).dlam
                            + (pow1_lj * s.at(sfe, 0) - pow_lj * s.at(sfe, 1)) * s.dlam;
                    }
                    if !ex.crf {
                        de_mult += (-dpow1_crf * c.at(se, 2) + dpow_crf * c.at(se, 3))
                            * l(CRF).dlam
                            + (pow1_crf * c.at(sfe, 2) - pow_crf * c.at(sfe, 3)) * c.dlam;
                    }
                    let diff = -(e_mult - e_s) / kbt;
                    let k = (i * nj + j) * nr_plam + (p - o.pmin);
                    x[k].push(de_mult);
                    vyvr[k].push(diff);
                    let combo = i * nj + j;
                    if let Some(col) = column(&o.bar_data, &bar_index, p, o.pmin) {
                        e_bar[combo * nr_bar + col].push(e_mult);
                    }
                    if let Some(col) = column(&o.dhdl_data, &dhdl_index, p, o.pmin) {
                        e_dhdl[combo * nr_dhdl + col].push(de_mult);
                    }
                }
            }
        }
    }
    if n_frames == 0 {
        return Err(match frames {
            None => "no frames read from @en_files/@fr_files".into(),
            Some(_) => "no frames".into(),
        });
    }
    if o.cpus != 1 {
        return Err(format!(
            "@cpus: newen ext_ti_ana runs on one thread (gromos++ without OpenMP refuses \
             @cpus != 1 too).\n\n{USAGE}"
        ));
    }

    let seed0 = o.seed.unwrap_or_else(wall_clock_seed);
    let mut predictions = Vec::with_capacity(combos);
    for (i, (slj_label, slj_c)) in slj_sets.iter().enumerate() {
        let mut seed = seed0.wrapping_add(i as u32);
        for (j, (scrf_label, scrf_c)) in scrf_sets.iter().enumerate() {
            let combo = i * nj + j;
            let mut plam = Vec::with_capacity(nr_plam);
            let mut dhdl = Vec::with_capacity(nr_plam);
            let mut error = Vec::new();
            let mut counts = Vec::new();
            for p in 0..nr_plam {
                let (xs, ys) = (&x[combo * nr_plam + p], &vyvr[combo * nr_plam + p]);
                plam.push(grid.p_lam[p + o.pmin]);
                let (lnx, sign) = ln_x_exp_ave(xs, ys)?;
                let lnexpave = lnexp_ave(ys);
                let final_v = (lnx - lnexpave).exp() * sign as f64;
                dhdl.push(final_v);
                if o.bootstrap > 0 {
                    let length = xs.len();
                    let mut boot = Vec::with_capacity(o.bootstrap);
                    for _ in 0..o.bootstrap {
                        let mut xn = Vec::with_capacity(length);
                        let mut yn = Vec::with_capacity(length);
                        for _ in 0..length {
                            let k = draw_index(&mut seed, length);
                            xn.push(xs[k]);
                            yn.push(ys[k]);
                        }
                        let (lnx, sign) = ln_x_exp_ave(&xn, &yn)?;
                        boot.push((lnx - lnexp_ave(&yn)).exp() * sign as f64);
                    }
                    error.push(population_variance(&boot).sqrt());
                }
                if o.countframes {
                    let deltag = -kbt * lnexpave;
                    counts.push(ys.iter().filter(|&&v| -kbt * v <= deltag + 0.001).count());
                }
            }
            let bar = o.bar_data.as_ref().map(|v| BarData {
                plam: if v.is_empty() {
                    grid.p_lam[o.pmin..o.pmax].to_vec()
                } else {
                    v.clone()
                },
                e_s: e_sim.clone(),
                e_p: e_bar[combo * nr_bar..(combo + 1) * nr_bar].to_vec(),
            });
            let dhdl_data = o.dhdl_data.as_ref().map(|v| DhdlData {
                plam: if v.is_empty() {
                    grid.p_lam[o.pmin..o.pmax].to_vec()
                } else {
                    v.clone()
                },
                de: e_dhdl[combo * nr_dhdl..(combo + 1) * nr_dhdl].to_vec(),
            });
            predictions.push(Prediction {
                labels: labelled.then_some((*slj_label, *scrf_label)),
                slj: *slj_c,
                scrf: *scrf_c,
                plam,
                dhdl,
                error: (o.bootstrap > 0).then_some(error),
                contributing_frames: o.countframes.then_some(counts),
                bar,
                dhdl_data,
            });
        }
    }
    Ok(Analysis {
        n_frames,
        predictions,
    })
}

/// Row 0 of `ABDIH`/`FREEABDIH`: the A and B dihedral energies (derivatives).
fn abdih(etrj: &EnergyTraj, name: &str) -> Result<[f64; 2], String> {
    let row = etrj
        .table(name)
        .and_then(|t| t.first())
        .ok_or_else(|| format!("the trajectory has no {name} subblock"))?;
    if row.len() < 2 {
        return Err(format!("{name}: expected 2 values, found {}", row.len()));
    }
    Ok([row[0], row[1]])
}

/// `<outdir>/<stem>_l<slam>[_slj<a>_scrf<b>].dat` (`ext_ti_ana.cc:965-967, 1094-1099`).
fn file_name(o: &Options, stem: &str, labels: Option<(i64, i64)>) -> String {
    let lam = cpp_fixed(o.slam, o.lam_precision);
    match labels {
        None => format!("{}/{stem}_l{lam}.dat", o.outdir),
        Some((a, b)) => format!("{}/{stem}_l{lam}_slj{a}_scrf{b}.dat", o.outdir),
    }
}

/// The text of `predicted_l<slam>.dat` (`ext_ti_ana.cc:929-1092`).
pub fn predicted_text(o: &Options, pr: &Prediction) -> String {
    let mut s = format!("#SLAM {}\n", g6(o.slam));
    let mut coeffs = o.coefficients;
    coeffs[SLJ] = pr.slj;
    coeffs[SCRF] = pr.scrf;
    if o.verbose {
        s.push_str(&format!(
            "# NLAMs={}, NLAMp={}, T={}\n",
            g6(o.nlams),
            g6(o.nlamp),
            g6(o.temp)
        ));
        s.push_str("#----- lambda coefficients: ------------\n");
        s.push_str("#        coeffs_sim    coeffs_prediction\n");
        for (t, name) in TYPES.iter().enumerate() {
            s.push_str(&format!(
                "# {name:>4}: {} | {}\n",
                coeff_vec(&o.coefficients_sim[t]),
                coeff_vec(&coeffs[t])
            ));
        }
        s.push_str("#---------------------------------------\n");
    } else {
        for t in [SLJ, SCRF] {
            s.push_str(&format!(
                "# {:>4}: {} | {}\n",
                TYPES[t],
                coeff_vec(&o.coefficients_sim[t]),
                coeff_vec(&coeffs[t])
            ));
        }
    }
    s.push_str(&format!("#plam {:>15}", "pred_dHdl"));
    if pr.error.is_some() {
        s.push_str(&format!("{:>15}", "err"));
    }
    if pr.contributing_frames.is_some() {
        s.push_str(&format!(" {:>15}", "contrib.frames"));
    }
    s.push('\n');
    for (k, (&l, &v)) in pr.plam.iter().zip(&pr.dhdl).enumerate() {
        s.push_str(&cpp_fixed(l, o.lam_precision));
        s.push_str(&format!(" {:>15}", cpp_fixed(v, 6)));
        if let Some(e) = &pr.error {
            s.push_str(&format!(" {:>15}", cpp_fixed(e[k], 6)));
        }
        if let Some(c) = &pr.contributing_frames {
            s.push_str(&format!(" {:>15}", c[k]));
        }
        s.push('\n');
    }
    s
}

/// The text of `bar_data_l<slam>.dat` (`ext_ti_ana.cc:958-994`), the input of `bar`.
pub fn bar_text(o: &Options, b: &BarData) -> String {
    let prec = o.lam_precision;
    let mut s = format!("#      nr lam_p\n{:>15}\n", b.plam.len());
    s.push_str(&format!("#         lam_s{:>19}\n", "lam_p ..."));
    s.push_str(&format!("{:>15}", cpp_fixed(o.slam, prec)));
    for &l in &b.plam {
        s.push_str(&format!("{:>15}", cpp_fixed(l, prec)));
    }
    s.push('\n');
    s.push_str(&format!("#           E_s{:>19}\n", "E_p(lam_p) ..."));
    for (n, &e) in b.e_s.iter().enumerate() {
        s.push_str(&format!("{:>15} ", cpp_fixed(e, 8)));
        for col in &b.e_p {
            s.push_str(&format!("{:>15} ", cpp_fixed(col[n], 8)));
        }
        s.push('\n');
    }
    s
}

/// The text of `dHdl_data_l<slam>.dat` (`ext_ti_ana.cc:997-1032`).
pub fn dhdl_text(o: &Options, d: &DhdlData, n_frames: usize) -> String {
    let prec = o.lam_precision;
    let mut s = format!("# lam_s {:>15}\n", cpp_fixed(o.slam, prec));
    s.push_str(&format!("#      nr lam_p\n{:>15}\n", d.plam.len()));
    s.push_str("#         lam_p\n");
    for &l in &d.plam {
        s.push_str(&format!("{:>15}", cpp_fixed(l, prec)));
    }
    s.push('\n');
    s.push_str("#          dE_p ...\n");
    for n in 0..n_frames {
        for col in &d.de {
            s.push_str(&format!("{:>15} ", cpp_fixed(col[n], 8)));
        }
        s.push('\n');
    }
    s
}

/// Write every file of `a` into `o.outdir` (created when missing, as gromos++ does); returns
/// their paths.
pub fn write_files(o: &Options, a: &Analysis) -> Result<Vec<String>, String> {
    let dir = std::path::Path::new(&o.outdir);
    if !dir.exists() {
        std::fs::create_dir(dir).map_err(|_| "Error creating directory!\n".to_string())?;
    } else if !dir.is_dir() {
        return Err("Error creating directory, file exists!\n".into());
    }
    let mut written = Vec::new();
    let mut put = |path: String, text: String| -> Result<(), String> {
        std::fs::write(&path, text).map_err(|e| format!("{path}: {e}"))?;
        written.push(path);
        Ok(())
    };
    for pr in &a.predictions {
        if let Some(b) = &pr.bar {
            put(file_name(o, "bar_data", pr.labels), bar_text(o, b))?;
        }
        if let Some(d) = &pr.dhdl_data {
            put(
                file_name(o, "dHdl_data", pr.labels),
                dhdl_text(o, d, a.n_frames),
            )?;
        }
        put(file_name(o, "predicted", pr.labels), predicted_text(o, pr))?;
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn data(name: &str) -> String {
        format!(
            "{}/tests/data/ext_ti_ana/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    /// The five-window CH4 chain's λ = 0.25 window, with `@bar_data` at the five window λ.
    fn window_options() -> Options {
        let argv: Vec<String> = [
            "ext_ti_ana",
            "@en_files",
            &data("windows/ch4_l0.25.tre"),
            "@fr_files",
            &data("windows/ch4_l0.25.trg"),
            "@library",
            &data("ene_ana.md++.lib"),
            "@temp",
            "300",
            "@nrlambdas",
            "21",
            "@minlam",
            "0",
            "@maxlam",
            "1",
            "@slam",
            "0.25",
            "@NLAMs",
            "1",
            "@bar_data",
            "0",
            "0.25",
            "0.5",
            "0.75",
            "1",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let argv = newen_io::normalize(argv).unwrap();
        Options::from_argv(&argv).unwrap()
    }

    fn frames_of(o: &Options) -> Vec<PrecalcFrame> {
        let lib = std::fs::read_to_string(o.library.as_ref().unwrap()).unwrap();
        let mut f = FileFrames {
            en: EnergySequence::new(o.en_files.clone(), "ENERTRJ", "lib"),
            fr: EnergySequence::new(o.fr_files.clone(), "FRENERTRJ", "lib"),
            etrj: EnergyTraj::from_library(&lib).unwrap(),
        };
        f.etrj.add_constant("MASS", 0.0);
        let mut out = Vec::new();
        while let Some(t) = f.next().unwrap() {
            out.push(PrecalcFrame {
                energies: t.se.to_vec(),
                derivatives: t.sfe.to_vec(),
                abdih: t.se_dih,
                free_abdih: t.sfe_dih,
            });
        }
        out
    }

    #[test]
    fn frames_in_memory_give_what_the_files_give() {
        let o = window_options();
        let frames = frames_of(&o);
        assert_eq!(frames.len(), 20);
        let from_files = analyse(&o).unwrap();
        assert_eq!(analyse_frames(&o, &frames).unwrap(), from_files);

        // the same window through `for_frames`: gromos++'s defaults, no files
        let m = Options::for_frames(
            300.0,
            21,
            0.0,
            1.0,
            0.25,
            1.0,
            vec![1.0, 0.0, 0.5, 0.25, 0.75],
        )
        .unwrap();
        let b = bar_file(&m, &frames).unwrap();
        let expect = from_files.predictions[0].bar.clone().unwrap();
        assert_eq!(b.slam, 0.25);
        assert_eq!(b.plam, expect.plam);
        assert_eq!(b.e_s, expect.e_s);
        assert_eq!(b.e_p, expect.e_p);
    }

    #[test]
    fn bar_file_needs_every_lambda_on_the_grid() {
        let frames = frames_of(&window_options());
        let m = Options::for_frames(300.0, 21, 0.0, 1.0, 0.25, 1.0, vec![0.25, 0.33]).unwrap();
        let e = bar_file(&m, &frames).unwrap_err();
        assert!(
            e.contains("0.33 is not a point of the PRECALCLAM grid"),
            "{e}"
        );
        let ok = Options::for_frames(300.0, 21, 0.0, 1.0, 0.25, 1.0, vec![0.25, 0.3]).unwrap();
        assert!(analyse_frames(&ok, &[]).unwrap_err().contains("no frames"));
        assert!(Options::for_frames(300.0, 1, 0.0, 1.0, 0.25, 1.0, vec![]).is_err());
    }
}
