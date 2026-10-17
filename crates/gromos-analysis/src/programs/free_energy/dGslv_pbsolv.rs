//! dGslv_pbsolv — the Poisson–Boltzmann correction to the electrostatic solvation free energy
//! (gromos++ `dGslv_pbsolv`, `programs/dGslv_pbsolv.cc`).
//!
//! Four finite-difference solves on the solute: non-periodic in vacuum and in `@epsNPBC`, periodic
//! in vacuum and in `@epsSOLV`. Their differences are the CB/NPBC and LS/PBC solvation free
//! energies, written to the `@verbose` log (`# DGRESULT NPBC|PBC`). `@schemeELEC RF` adds two
//! periodic FFT solves of the reaction potential, lattice-sum and reaction-field, whose free
//! energies end their log sections (`# Convergence after N steps with dG`). stdout carries the
//! per-atom potentials of every solve. gromos++ documents the correction (LS: CB/NPBC[FD] −
//! LS/PBC[FD]; RF: that + LS/PBC[FFT] − RF/PBC[FFT]) without computing it;
//! [`DgslvPbsolv::documented_correction`] does.
//!
//! Deliberate departures from gromos++: `@epsRF 1` is refused (gromos++ solves it without keeping
//! the potentials, then aborts on the empty column, `dGslv_pbsolv.cc:299`), an FFT read-out atom
//! whose interpolation corners leave the grid is an error (gromos++ reads outside it,
//! `FFTPoissonIterator.cc:728-737`), FFT grids need more than 100 points (gromos++ logs
//! element 100), `@gridspacing` must be positive in coordinate mode and
//! `@nogridpoints` in PQR mode (gromos++ divides by zero), `@increasegrid` takes exactly three
//! values (gromos++ reads the neighbouring arguments when given fewer), numbers must parse
//! completely (gromos++ uses `atof`/`atoi`), solvent atoms are refused in `@atoms` and
//! `@atomsTOcharge`, a grid index outside the grid is an error (undefined behaviour in gromos++),
//! and PQR mode ends with success (gromos++ calls `exit(1)`, `dGslv_pbsolv.cc:1060`).

use crate::args::Arguments;
use crate::pb::{
    atomic_radii_vdw, dgslv_fd, dgslv_fft, g, pb_radius, DgslvFd, DgslvFft, FdSetup, FftSetup,
    PbParameters, PbSystem, RadiusDefinition, SolveOrder, GROMOSPP_FOUR_PI_EPS_I,
};
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::ordered_atoms;
use gromos_core::topology::MolTypeAtom;
use gromos_core::{Topology, Vec3};
use gromos_io::args::{Program, Status};
use gromos_io::coordinate::read_coordinates;
use gromos_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use std::io::Write;

pub const USAGE: &str = "# dGslv_pbsolv

# USAGE

# -----------------------------------------------------------------------------------------
# if you have a gromos coordinate file:
\t@topo            <molecular topology file>
\t@pbc             <boundary type> [<gather method>]
\t@coord           <gromos96 coordinates>
\t@probeIAC        <integer atom code for radius calculation
\t                  (for water oxygen it would be 4 or 5, depending on the ff)>
\t@atoms           <atoms to include for the pb calculations
\t                  expected in gromos format>
\t@atomsTOcharge   <atoms to charge; expected in gromos format>
\t@rminORsigma     <how to calculate the radii - rmin (0) or sigma (1); default: 0>

# -----------------------------------------------------------------------------------------
# if you have a pqr file:
# short description - the last two elements in each ATOM or HETATOM line
# are the atom charge and the atom radii (in Angstrom);
# hydrogen atoms can have zero radius (they get 0.05 nm)
\t@pqr             <pqr file>
\t@coordinates     <box coordinates in X Y Z direction (in nm) that were used
\t                  in the simulation>
\t@atoms           <atoms to include for the pb calculations>
\t                  typically, all atoms should be included
\t                  (simply type 'a' to include all atoms)
\t                  if not, atoms are indexed from 1 and individual atoms
\t                  can be picked by
\t                  using a comma delimiter; ranges can be specified using a hyphen
\t                  e.g. '1,2,3,4,5,8,9' is the same as '1-5,8,9'>
\t@atomsTOcharge   <atoms to charge; e.g. the perturbed atoms in the
\t                  free energy calculations;
\t                  same rules as above to pick individual atoms>

# -----------------------------------------------------------------------------------------
# general input:
\t@schemeELEC      <electrostatics scheme: LS or RF>
\t if RF: @rcut    <cutoff distance in nm>
\t if RF: @epsRF   <dielectric permittivity of the reaction field (0: conducting; 1 is not supported)>
\t@epsSOLV         <relative dielectric permittivity of the employed solvent model>
\t@gridspacing     <grid spacing in nm; should not be much higher than 0.02 but
\t                  lower numbers are computationally much more expensive in
\t                  terms of time and memory usage>
\t[@nogridpoints    <number of grid points per axis; PQR input only, where it sets the
\t                  grid spacing to the X box edge / this number>]
\t[@NPBCsize        <optional; if argument is given change the NPBC grid size, based
\t                  on the box dimension given in nm; default the npbc box is
\t                  extended by 4 nm>]
\t[@epsNPBC        <relative dielectrict permittivity for calculation of macroscopic,
\t                  non-periodic boundary conditions; default 78.4 (for water)>]
\t[@maxiter        <maximum number of iteration steps; default 600>]
\t[@cubesFFT       <number of cubes in the fast Fourier transformation for
\t                  boundary smoothing; default 4>]
\t[@probeRAD       <probe radius in nm; default 0.14 (for water)>]
\t[@HRAD           <accepted and ignored, as in gromos++>]
\t[@radscal        <scale non-H radii with this factor (use this only in case you
\t                  want to play with radii); default 1.0>]
\t[@increasegrid   <takes three integer values for X Y Z; grid for PBC calculations gets increased by the number of given gridpoints;
\t                  may be usefull if atoms close to the border of the box extend the grid!>]
\t[@verbose        <path to log file to document status and errors>]";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "atoms",
    "atomsTOcharge",
    "coord",
    "pqr",
    "schemeELEC",
    "epsSOLV",
    "epsRF",
    "rcut",
    "gridspacing",
    "coordinates",
    "maxiter",
    "nogridpoints",
    "NPBCsize",
    "cubesFFT",
    "probeIAC",
    "probeRAD",
    "HRAD",
    "epsNPBC",
    "radscal",
    "rminORsigma",
    "increasegrid",
    "verbose",
];

pub const PROGRAM: Program = Program {
    name: "dGslv_pbsolv",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// gromos++ `throw gromos::Exception("dGslv_pbsolv", msg)`.
fn fail(msg: &str) -> String {
    format!("dGslv_pbsolv: {msg}")
}

/// The first value of `key` as a number, `None` when the key has no value (gromos++
/// `args.count(key) > 0`).
fn number(args: &Arguments, key: &str) -> Result<Option<f64>, String> {
    match args.values(key).first() {
        None => Ok(None),
        Some(v) => v
            .parse::<f64>()
            .map(Some)
            .map_err(|_| fail(&format!("@{key}: cannot parse '{v}' as a number"))),
    }
}

fn integer(args: &Arguments, key: &str) -> Result<Option<i64>, String> {
    match args.values(key).first() {
        None => Ok(None),
        Some(v) => v
            .parse::<i64>()
            .map(Some)
            .map_err(|_| fail(&format!("@{key}: cannot parse '{v}' as an integer"))),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Scheme {
    Ls,
    Rf,
}

/// The arguments gromos++ reads before it looks at the system (`dGslv_pbsolv.cc:420-528`).
struct Settings {
    scheme: Scheme,
    rcut: f64,
    eps_rf: f64,
    fftcub: usize,
    probe_rad: f64,
    hydrogen_rad: f64,
    radscal: f64,
    increasegrid: [i64; 3],
    epssolvent: f64,
    eps_npbc: f64,
    maxiter: usize,
    definition: RadiusDefinition,
}

fn read_settings(args: &Arguments, log: &mut String) -> Result<Settings, String> {
    let scheme_text = args
        .values("schemeELEC")
        .first()
        .map(|s| s.to_uppercase())
        .unwrap_or_default();
    let scheme = match scheme_text.as_str() {
        "LS" => Scheme::Ls,
        "RF" => Scheme::Rf,
        _ => {
            return Err(fail(&format!(
                "schemeELEC format {scheme_text} unknown. Exiting ..."
            )))
        },
    };
    log.push_str(&format!("# READ: schemeELEC {scheme_text}\n"));

    let mut rcut = 0.0;
    if scheme == Scheme::Rf {
        rcut = number(args, "rcut")?.unwrap_or(-1.0);
        if rcut < 0.0 {
            return Err(fail("No or negative RF cutoff given. Exiting ..."));
        }
        log.push_str(&format!("# READ: rcut {}\n", g(rcut)));
    }

    let probe_rad = number(args, "probeRAD")?.unwrap_or(0.14);
    if probe_rad < 0.0 {
        return Err(fail(
            "The probe radius (probeRAD) must not be negative. Exiting ...",
        ));
    }
    log.push_str(&format!("# READ: probe_rad {}\n", g(probe_rad)));

    // @radH is read by gromos++ but is not one of its known arguments, so it is always 0.05.
    let hydrogen_rad = 0.05;
    log.push_str(&format!("# READ: hydrogen_rad {}\n", g(hydrogen_rad)));

    let radscal = number(args, "radscal")?.unwrap_or(1.0);
    log.push_str(&format!("# READ: radscal {}\n", g(radscal)));

    let mut increasegrid = [0i64; 3];
    if args.count("increasegrid") > 0 {
        let v = args.values("increasegrid");
        if v.len() != 3 {
            return Err(fail(&format!(
                "@increasegrid takes three integers (X Y Z), got {}",
                v.len()
            )));
        }
        for (d, s) in v.iter().enumerate() {
            increasegrid[d] = s
                .parse::<i64>()
                .map_err(|_| fail(&format!("@increasegrid: cannot parse '{s}' as an integer")))?;
        }
        log.push_str(&format!(
            "# READ: increasegrid {} {} {}\n",
            increasegrid[0], increasegrid[1], increasegrid[2]
        ));
    }

    let epssolvent = number(args, "epsSOLV")?.unwrap_or(0.0);
    if epssolvent < 1.0 {
        return Err(fail("The solvent permittivity (epsSOLV) not given or given value is smaller than 1. Exiting ..."));
    }
    log.push_str(&format!("# READ: epssolvent {}\n", g(epssolvent)));

    let mut eps_rf = -1.0;
    if scheme == Scheme::Rf {
        eps_rf = number(args, "epsRF")?.unwrap_or(-1.0);
        if eps_rf < 1.0 && eps_rf != 0.0 {
            return Err(fail("Reaction field permittivity (epsRF) not given or given value not allowed. Exiting ..."));
        }
        log.push_str(&format!("# READ: epsRF {}\n", g(eps_rf)));
    }

    let eps_npbc = number(args, "epsNPBC")?.unwrap_or(78.4);
    if eps_npbc < 1.0 {
        return Err(fail("The permittivity for the non-periodic boundary conditions (epsNPBC) must not be smaller than 1. Exiting ..."));
    }
    log.push_str(&format!("# READ: epsNPBC {}\n", g(eps_npbc)));

    let maxiter = integer(args, "maxiter")?.unwrap_or(600);
    if maxiter <= 0 {
        return Err(fail(
            "The maximum number of iterations (maxiter) must be positive. Exiting ...",
        ));
    }
    log.push_str(&format!("# READ: maxiter {maxiter}\n"));

    let mut definition = RadiusDefinition::Rmin;
    if args.count("coord") > 0 {
        let r = integer(args, "rminORsigma")?.unwrap_or(0);
        definition = match r {
            0 => RadiusDefinition::Rmin,
            1 => RadiusDefinition::Sigma,
            _ => {
                return Err(fail(
                    "The rminORsigma flag should be 0 (rmin) or 1 (sigma). Exiting ...",
                ))
            },
        };
        log.push_str(&format!("# READ: rminORsigma {r}\n"));
    }

    let fftcub = integer(args, "cubesFFT")?.unwrap_or(4);
    if fftcub <= 0 {
        return Err(fail("The FFT cubelet number must be positive. Exiting ..."));
    }
    log.push_str(&format!("# READ: fftcub {fftcub}\n"));

    // gromos++ takes epsRF within tiny_real of 1 as "no reaction-field correction" (SC), solves it
    // without keeping the potentials and aborts on the empty table column (dGslv_pbsolv.cc:280-299).
    if scheme == Scheme::Rf && (eps_rf - 1.0).abs() <= 0.0000000001 {
        return Err(fail(
            "@epsRF 1 (no reaction-field correction) is not supported: gromos++ does not keep the \
             potentials of that solve and aborts. Use @epsRF 0 (conducting) or a value above 1",
        ));
    }

    Ok(Settings {
        scheme,
        rcut,
        eps_rf,
        fftcub: fftcub as usize,
        probe_rad,
        hydrogen_rad,
        radscal,
        increasegrid,
        epssolvent,
        eps_npbc,
        maxiter: maxiter as usize,
        definition,
    })
}

/// The atoms of every value of `key`, in the order written, each once (gromos++ `AtomSpecifier`
/// keeps insertion order and skips repeats). `prefix` is prepended to values without a molecule.
fn select(
    args: &Arguments,
    key: &str,
    topo: &Topology,
    prefix: Option<&str>,
) -> Result<Vec<usize>, String> {
    let mut out: Vec<usize> = Vec::new();
    for v in args.values(key) {
        let spec = match prefix {
            Some(p) if !v.contains(':') => format!("{p}{v}"),
            _ => v.clone(),
        };
        for a in ordered_atoms(&spec, topo).map_err(|e| fail(&format!("@{key}: {e}")))? {
            if !out.contains(&a) {
                out.push(a);
            }
        }
    }
    Ok(out)
}

/// `ceil(len / h)` as gromos++ turns it into an int.
fn grid_points(len: f64, h: f64) -> i64 {
    (len / h).ceil() as i64
}

/// Grid sizes after the program's own checks; gromos++ overflows `int` on absurd grids.
fn grid(n: [i64; 3], what: &str) -> Result<[usize; 3], String> {
    let total = n.iter().try_fold(1i64, |acc, &x| acc.checked_mul(x));
    match total {
        Some(t) if n.iter().all(|&x| x > 0) && t <= i32::MAX as i64 => {
            Ok([n[0] as usize, n[1] as usize, n[2] as usize])
        },
        _ => Err(fail(&format!(
            "the {what} grid {} x {} x {} is empty or too large",
            n[0], n[1], n[2]
        ))),
    }
}

struct Grids {
    pbc: [usize; 3],
    npbc: [usize; 3],
    h: f64,
}

/// `@nogridpoints`, `@gridspacing`, `@NPBCsize` and `@increasegrid` (`dGslv_pbsolv.cc:594-660`,
/// PQR `:896-955`). `box_edges` is the X, Y, Z edge; `pqr` selects the PQR rule where the
/// spacing is `X edge / @nogridpoints`.
fn read_grids(
    args: &Arguments,
    settings: &Settings,
    box_edges: [f64; 3],
    pqr: bool,
    log: &mut String,
) -> Result<Grids, String> {
    let [a, b, c] = box_edges;
    let mut ngrid = [0i64; 3];
    let mut nogridpoints = 0i64;
    if let Some(n) = number(args, "nogridpoints")? {
        // atof into an int
        nogridpoints = n as i64;
        if nogridpoints < 0 {
            return Err(fail(
                "The number of grids must not be 0 or negative. Exiting ...",
            ));
        }
        log.push_str(&format!("# READ: nogridpoints {nogridpoints}\n"));
        ngrid = [nogridpoints; 3];
    }

    let mut h = 0.0;
    if let Some(s) = number(args, "gridspacing")? {
        h = s;
        if h < 0.0 {
            return Err(fail("The grid spacing must not be negative. Exiting ..."));
        }
        log.push_str(&format!("# READ: gridspacing {}\n", g(h)));
        if h > 0.0 {
            ngrid = [grid_points(a, h), grid_points(b, h), grid_points(c, h)];
        }
    }
    if pqr {
        if nogridpoints <= 0 {
            return Err(fail(
                "@nogridpoints must be given and positive with @pqr: the grid spacing is the X box \
                 edge divided by it",
            ));
        }
        if h == 0.0 && args.count("gridspacing") > 0 {
            return Err(fail("@gridspacing must be positive"));
        }
        h = a / nogridpoints as f64;
    } else if h <= 0.0 {
        return Err(fail(
            "@gridspacing must be given and positive: the grids are built from it",
        ));
    }

    let npbc = if let Some(size) = number(args, "NPBCsize")? {
        if size <= 0.0 {
            return Err(fail(
                "The box size for NPBC must not be 0 or negative. Exiting ...",
            ));
        }
        log.push_str(&format!("# READ: NPBCsize {}\n", g(size)));
        [grid_points(size, h); 3]
    } else {
        [
            grid_points(a + 4.0, h),
            grid_points(b + 4.0, h),
            grid_points(c + 4.0, h),
        ]
    };

    for (d, axis) in ["x", "y", "z"].iter().enumerate() {
        let inc = settings.increasegrid[d];
        if inc > 0 {
            ngrid[d] += inc;
            log.push_str(&format!(
                "Increased the grid on the {axis} axis by {inc} according to @increasegrid flag!\n"
            ));
        }
    }

    if !pqr {
        log.push_str(&format!("# Boxsize: {}\n", g(a)));
        log.push_str(&format!("# Adapted gridspacing: {}\n", g(h)));
    }
    log.push_str(&format!(
        "# Calculated number of gridpoints PBC (x,y,z): {} {} {}\n",
        ngrid[0], ngrid[1], ngrid[2]
    ));
    log.push_str(&format!(
        "# Calculated number of gridpoints NPBC (x,y,z): {} {} {}\n",
        npbc[0], npbc[1], npbc[2]
    ));
    Ok(Grids {
        pbc: grid(ngrid, "PBC")?,
        npbc: grid(npbc, "NPBC")?,
        h,
    })
}

/// One read-out atom (`@atomsTOcharge`), a row of the stdout table.
#[derive(Debug, Clone)]
pub struct ReadoutAtom {
    /// Index among the solute atoms (0-based).
    pub atom: usize,
    /// Residue number within its molecule, from 0 (printed + 1).
    pub resnum: usize,
    pub name: String,
    pub charge: f64,
}

/// Everything one `dGslv_pbsolv` run computes.
#[derive(Debug, Clone)]
pub struct DgslvPbsolv {
    /// The read-out atoms, in the order of the potentials.
    pub atoms: Vec<ReadoutAtom>,
    /// The four finite-difference solves and their two free energies.
    pub fd: DgslvFd,
    /// The LS and RF FFT solves (`@schemeELEC RF` only).
    pub fft: Option<DgslvFft>,
}

impl DgslvPbsolv {
    /// The correction gromos++ documents but does not compute (`dGslv_pbsolv.cc:54-60`), kJ/mol:
    /// CB/NPBC[FD] − LS/PBC[FD], plus LS/PBC[FFT] − RF/PBC[FFT] with `@schemeELEC RF`.
    pub fn documented_correction(&self) -> f64 {
        let fd = self.fd.dg_npbc - self.fd.dg_pbc;
        match &self.fft {
            None => fd,
            Some(fft) => fd + fft.ls.dg - fft.rf.dg,
        }
    }
}

/// `writeout` (`dGslv_pbsolv.cc:306-331`).
fn writeout(out: &mut dyn Write, r: &DgslvPbsolv) -> Result<(), String> {
    let fd = &r.fd;
    let mut text = format!(
        "{:>6}{:>7}{:>7}{:>10}{:>16}{:>16}{:>16}{:>16}",
        "# atom", "nres", "name", "charge", "NPBC_SLV", "NPBC_VAC", "PBC_SLV", "PBC_VAC"
    );
    if r.fft.is_some() {
        text.push_str(&format!("{:>16}{:>16}", "FFT_LS_PBC", "FFT_RF_PBC"));
    }
    text.push('\n');
    for (i, a) in r.atoms.iter().enumerate() {
        text.push_str(&format!(
            "{:>6}{:>7}{:>7}{:>10.4}{:>16.7}{:>16.7}{:>16.7}{:>16.7}",
            i + 1,
            a.resnum + 1,
            a.name,
            a.charge,
            fd.npbc_slv.potentials[i],
            fd.npbc_vac.potentials[i],
            fd.pbc_slv.potentials[i],
            fd.pbc_vac.potentials[i],
        ));
        if let Some(fft) = &r.fft {
            text.push_str(&format!(
                "{:>16.7}{:>16.7}",
                fft.ls.potentials[i], fft.rf.potentials[i]
            ));
        }
        text.push('\n');
    }
    text.push('\n');
    out.write_all(text.as_bytes()).map_err(io_err)
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    evaluate(argv, out).map(|_| ())
}

/// The program as a library function: runs `argv` (program name first, `@key` form), writes the
/// table (and any gather warning) to `out` and the `@verbose` log file, and returns the results
/// with the log text.
pub fn evaluate(argv: &[String], out: &mut dyn Write) -> Result<(DgslvPbsolv, String), String> {
    if argv.len() < 2 {
        return Err(USAGE.to_string());
    }
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let log_path = if args.has("verbose") {
        match args.values("verbose").first() {
            Some(p) if !p.is_empty() => Some(p.clone()),
            _ => return Err(fail("verbose - no file name given")),
        }
    } else {
        None
    };
    let mut log = String::new();
    let result = run_logged(&args, out, &mut log);
    if let Some(path) = log_path {
        std::fs::write(&path, &log).map_err(|e| fail(&format!("@verbose {path}: {e}")))?;
    }
    result.map(|r| (r, log))
}

fn run_logged(
    args: &Arguments,
    out: &mut dyn Write,
    log: &mut String,
) -> Result<DgslvPbsolv, String> {
    let settings = read_settings(args, log)?;
    if args.count("coord") > 0 {
        coordinate_mode(args, &settings, out, log)
    } else if args.count("pqr") > 0 {
        pqr_mode(args, &settings, out, log)
    } else {
        Err(fail("No input file. Use @coord or @pqr"))
    }
}

fn log_readout(log: &mut String, readout: &[usize], topo: &Topology) {
    log.push_str("# READ: atomsTOcharge \n");
    for &a in readout {
        log.push_str(&format!(
            "# READ: mol {} {}\n",
            topo.molecule_nr(a).unwrap_or(0),
            topo.atom_name(a).unwrap_or("")
        ));
    }
}

fn coordinate_mode(
    args: &Arguments,
    settings: &Settings,
    out: &mut dyn Write,
    log: &mut String,
) -> Result<DgslvPbsolv, String> {
    let topo_path = args
        .value("topo")
        .map_err(|_| fail("@topo is required with @coord"))?;
    let parsed = read_topology_file(topo_path).map_err(|e| fail(&format!("@topo: {e}")))?;
    let fpepsi = parsed.physical_constants.four_pi_eps_i;
    let mut topo = build_topology(parsed);
    let n_solute = topo.num_solute_atoms();

    let atoms = select(args, "atoms", &topo, None)?;
    if atoms.is_empty() {
        return Err(fail("No atoms specified. Exiting ..."));
    }
    let readout = select(args, "atomsTOcharge", &topo, None)?;
    if readout.is_empty() {
        return Err(fail("No atoms to charge specified. Exiting ..."));
    }
    for (key, set) in [("atoms", &atoms), ("atomsTOcharge", &readout)] {
        if let Some(a) = set.iter().find(|&&a| a >= n_solute) {
            return Err(fail(&format!(
                "@{key}: atom {} is a solvent atom; only solute atoms are supported",
                a + 1
            )));
        }
    }
    log_readout(log, &readout, &topo);

    let pbc = Pbc::from_args(args)?;
    // GatherParser.cc:79-85: a boundary without a gather method announces the COG default.
    let pbc_values = args.values("pbc");
    if pbc_values.len() == 1 && pbc_values[0] != "v" {
        out.write_all(
            b"###### GATHER WARNING ######\n\
              # NO Gathering method specified ! \n\
              # Thus : if the system is requested to be gathered, \n\
              # the gathering will be done using the COG gather method\n",
        )
        .map_err(io_err)?;
    }

    let coord_path = args.value("coord")?;
    let conf = read_coordinates(coord_path).map_err(|e| fail(&format!("@coord: {e}")))?;
    let mut positions = conf.positions;
    if positions.is_empty() {
        return Err(fail("No position block in coordinate file. Exiting ..."));
    }
    solvate_to_atoms(&mut topo, positions.len()).map_err(|e| fail(&format!("@coord: {e}")))?;
    pbc.gather(&topo, &mut positions, &pbc.periodicity(conf.box_dims));

    let edges = [conf.box_dims.x, conf.box_dims.y, conf.box_dims.z];
    if edges == [0.0; 3] {
        return Err(fail("No box block in coordinate file. Exiting ..."));
    }
    if edges.iter().any(|&e| e <= 0.0) {
        return Err(fail("At least one coordinate is zero. Exiting ..."));
    }

    let grids = read_grids(args, settings, edges, false, log)?;
    let ppp = PbParameters::new_logged(settings.epssolvent, fpepsi, log);

    let probe_iac = integer(args, "probeIAC")?.unwrap_or(0);
    if probe_iac <= 0 {
        return Err(fail(
            "The probe integer atom code (probeIAC) must not be negative or 0. Exiting ...",
        ));
    }
    log.push_str(&format!("# READ: probe_iac {probe_iac}\n"));
    log.push_str(&format!("# RADII: based on probe_iac {probe_iac}\n"));
    log.push_str(match settings.definition {
        RadiusDefinition::Sigma => "# RADII: based on sigma\n",
        RadiusDefinition::Rmin => "# RADII: based on rmin\n",
    });

    let mut radii = atomic_radii_vdw(&topo, (probe_iac - 1) as usize, settings.probe_rad)
        .map_err(|e| fail(&e))?;
    for &a in &atoms {
        radii[a] = pb_radius(
            radii[a],
            settings.definition,
            settings.probe_rad,
            settings.hydrogen_rad,
            settings.radscal,
            ppp.tiny_real,
        );
    }

    positions.truncate(n_solute);
    let sys = PbSystem {
        positions,
        charges: topo.charge[..n_solute].to_vec(),
        radii,
        atoms,
        readout,
        fpepsi,
    };
    let atoms: Vec<ReadoutAtom> = sys
        .readout
        .iter()
        .map(|&a| ReadoutAtom {
            atom: a,
            resnum: residue_in_molecule(&topo, a),
            name: topo.atom_name(a).unwrap_or("").to_string(),
            charge: sys.charges[a],
        })
        .collect();
    solve_and_write(
        &sys,
        settings,
        &grids,
        SolveOrder::VacuumFirst,
        atoms,
        out,
        log,
    )
}

/// The FD solves, the FFT pair with `@schemeELEC RF`, and the table.
fn solve_and_write(
    sys: &PbSystem,
    settings: &Settings,
    grids: &Grids,
    order: SolveOrder,
    atoms: Vec<ReadoutAtom>,
    out: &mut dyn Write,
    log: &mut String,
) -> Result<DgslvPbsolv, String> {
    let fd = dgslv_fd(sys, &fd_setup(settings, grids), order, log)?;
    let fft = match settings.scheme {
        Scheme::Ls => None,
        Scheme::Rf => Some(dgslv_fft(
            sys,
            &FftSetup {
                grid: grids.pbc,
                h: grids.h,
                eps_solv: settings.epssolvent,
                eps_rf: settings.eps_rf,
                rcut: settings.rcut,
                ncubes: settings.fftcub,
                maxiter: settings.maxiter,
            },
            log,
        )?),
    };
    let result = DgslvPbsolv { atoms, fd, fft };
    writeout(out, &result)?;
    Ok(result)
}

fn fd_setup(settings: &Settings, grids: &Grids) -> FdSetup {
    FdSetup {
        npbc_grid: grids.npbc,
        pbc_grid: grids.pbc,
        h: grids.h,
        eps_solv: settings.epssolvent,
        eps_npbc: settings.eps_npbc,
        maxiter: settings.maxiter,
    }
}

/// gromos++ numbers residues from 0 within each molecule (`LinearTopology.cc:254-259`).
fn residue_in_molecule(topo: &Topology, a: usize) -> usize {
    let first = topo
        .molecule_nr(a)
        .and_then(|m| topo.molecules.get(m))
        .map_or(a, |r| r.start);
    let r = topo.residue_nr(a).unwrap_or(1);
    let r0 = topo.residue_nr(first).unwrap_or(1);
    r.saturating_sub(r0)
}

/// One ATOM/HETATM record of a PQR file.
struct PqrAtom {
    name: String,
    pos: Vec3,
    charge: f64,
    radius: f64,
}

/// gromos++ `InPDB::readPQR` (`InPDB.cc:506-627`): ATOM and HETATM lines of 10 or 11
/// whitespace-separated fields (the 11th is a chain identifier after the residue name), charge
/// and radius last; positions and radius in Å.
fn read_pqr(path: &str) -> Result<Vec<PqrAtom>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| fail(&format!("Could not open the PQR file {path}: {e}")))?;
    let mut atoms = Vec::new();
    for line in text.lines() {
        let record: String = line.chars().take(6).collect();
        let record = record.split_whitespace().next().unwrap_or("");
        if record != "ATOM" && record != "HETATM" {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 10 {
            return Err(fail("At least one ATOM or HETATOM line has less than 10 elements. Don't know what to do with it!"));
        }
        if f.len() > 11 {
            return Err(fail("At least one ATOM or HETATOM line has more than 11 elements. Don't know what to do with it!"));
        }
        let o = f.len() - 10; // the chain identifier shifts the numbers by one
        let bad = || fail(&format!("bad line in PDB file:\n{line}"));
        let num = |s: &str| s.parse::<f64>().map_err(|_| bad());
        f[1].parse::<i64>().map_err(|_| bad())?;
        f[4 + o].parse::<i64>().map_err(|_| bad())?;
        atoms.push(PqrAtom {
            name: f[2].to_string(),
            pos: Vec3::new(num(f[5 + o])?, num(f[6 + o])?, num(f[7 + o])?),
            charge: num(f[8 + o])?,
            radius: num(f[9 + o])?,
        });
    }
    Ok(atoms)
}

fn pqr_mode(
    args: &Arguments,
    settings: &Settings,
    out: &mut dyn Write,
    log: &mut String,
) -> Result<DgslvPbsolv, String> {
    let pqr = read_pqr(args.value("pqr")?)?;

    // One molecule of the PQR atoms (gromos++ makes them argon, IAC 31; the type is never read).
    let mut topo = Topology::new();
    for a in &pqr {
        topo.moltypes[0].atoms.push(MolTypeAtom {
            name: a.name.clone(),
            residue_nr: 1,
            residue_name: String::new(),
            iac: 30,
            mass: 0.0,
            charge: a.charge,
            ..Default::default()
        });
        topo.iac.push(30);
        topo.mass.push(0.0);
        topo.charge.push(a.charge);
    }
    topo.molecules.push(0..pqr.len());

    // A specifier without a molecule refers to the one molecule (dGslv_pbsolv.cc:837-861).
    let atoms = select(args, "atoms", &topo, Some("1:"))?;
    if atoms.is_empty() {
        return Err(fail("No atoms specified. Exiting ..."));
    }
    let readout = select(args, "atomsTOcharge", &topo, Some("1:"))?;
    if readout.is_empty() {
        return Err(fail("No atoms to charge specified. Exiting ..."));
    }
    log_readout(log, &readout, &topo);

    let mut edges = [0.0; 3];
    for (d, v) in args.values("coordinates").iter().take(3).enumerate() {
        edges[d] = v
            .parse::<f64>()
            .map_err(|_| fail(&format!("@coordinates: cannot parse '{v}' as a number")))?;
    }
    if edges.iter().any(|&e| e <= 0.0) {
        return Err(fail("Coordinates must be positive and given in X Y Z (space-separated, no comma inbetween). Exiting ..."));
    }

    let grids = read_grids(args, settings, edges, true, log)?;

    // No topology: gromos++ warns once, when PB_Parameters first asks for 1/(4πε₀).
    eprintln!(
        "WARNING: No topology given or no value for four_pi_eps_i in topology file.\n         \
         Setting four_pi_eps_i to hardcoded value from gmath/Physics.cc:\n         \
         four_pi_eps_i = {}",
        g(GROMOSPP_FOUR_PI_EPS_I)
    );
    let ppp = PbParameters::new_logged(settings.epssolvent, GROMOSPP_FOUR_PI_EPS_I, log);

    let mut radii: Vec<f64> = pqr.iter().map(|a| a.radius / 10.0).collect();
    for (i, &a) in atoms.iter().enumerate() {
        // PQR input never reads @rminORsigma (dGslv_pbsolv.cc:511).
        let r = radii[a];
        radii[a] = pb_radius(
            r,
            RadiusDefinition::Rmin,
            settings.probe_rad,
            settings.hydrogen_rad,
            settings.radscal,
            ppp.tiny_real,
        );
        if r.abs() < ppp.tiny_real {
            log.push_str(&format!(
                "# !!!! H atom number (0 radius) {i} : assign rad(i) = {}\n",
                g(radii[a])
            ));
        }
        log.push_str(&format!(
            "# radius of atom {i} : rad(i) = {}\n",
            g(radii[a])
        ));
    }
    for (i, &a) in readout.iter().enumerate() {
        log.push_str(&format!(
            "# radius of atom_to_charge {i} : rad(i) = {}\n",
            g(radii[a])
        ));
    }

    let sys = PbSystem {
        positions: pqr.iter().map(|a| a.pos / 10.0).collect(),
        charges: pqr.iter().map(|a| a.charge).collect(),
        radii,
        atoms,
        readout,
        fpepsi: GROMOSPP_FOUR_PI_EPS_I,
    };
    // gromos++ prints residue 1 for every PQR atom (its residue map is empty).
    let atoms: Vec<ReadoutAtom> = sys
        .readout
        .iter()
        .map(|&a| ReadoutAtom {
            atom: a,
            resnum: 0,
            name: pqr[a].name.clone(),
            charge: pqr[a].charge,
        })
        .collect();
    solve_and_write(
        &sys,
        settings,
        &grids,
        SolveOrder::SolventFirst,
        atoms,
        out,
        log,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(a: &[&str]) -> Vec<String> {
        let mut v = vec!["dGslv_pbsolv".to_string()];
        v.extend(a.iter().map(|s| s.to_string()));
        gromos_io::normalize(v).unwrap()
    }

    fn err(a: &[&str]) -> String {
        let mut out = Vec::new();
        run(&argv(a), &mut out).expect_err("should be refused")
    }

    #[test]
    fn refuses_unknown_and_malformed_arguments() {
        // gromos++ reads @radH but does not list it among its known arguments.
        assert!(err(&["@radH", "0.1"]).contains("unknown argument @radH"));
        assert!(err(&["@gridpointsXYZ", "10", "10", "10"]).contains("unknown argument"));
        assert!(err(&["@schemeELEC", "XX"]).contains("schemeELEC format XX unknown"));
        assert!(err(&["@schemeELEC", "ls"]).contains("solvent permittivity"));
        assert!(err(&[
            "@schemeELEC",
            "LS",
            "@epsSOLV",
            "66.6",
            "@increasegrid",
            "1",
            "2"
        ])
        .contains("three integers"));
        assert!(
            err(&["@schemeELEC", "LS", "@epsSOLV", "66.6", "@maxiter", "0"]).contains("maxiter")
        );
        assert!(err(&["@schemeELEC", "LS", "@epsSOLV", "abc"]).contains("cannot parse"));
        assert!(err(&["@schemeELEC", "LS", "@epsSOLV", "66.6"]).contains("No input file"));
        assert!(err(&[
            "@schemeELEC",
            "RF",
            "@rcut",
            "1.4",
            "@epsSOLV",
            "66.6",
            "@epsRF",
            "1"
        ])
        .contains("@epsRF 1"));
        assert!(err(&["@verbose"]).contains("no file name"));
        // @HRAD is accepted and ignored.
        assert!(
            err(&["@HRAD", "0.1", "@schemeELEC", "LS", "@epsSOLV", "66.6"])
                .contains("No input file")
        );
    }

    #[test]
    fn grid_sizes_follow_gromospp() {
        let settings = Settings {
            scheme: Scheme::Ls,
            rcut: 0.0,
            eps_rf: -1.0,
            fftcub: 4,
            probe_rad: 0.14,
            hydrogen_rad: 0.05,
            radscal: 1.0,
            increasegrid: [0; 3],
            epssolvent: 66.6,
            eps_npbc: 78.4,
            maxiter: 600,
            definition: RadiusDefinition::Rmin,
        };
        let edge = 3.767055681;
        let args = Arguments::from_vec(argv(&["@gridspacing", "0.1"]), KNOWN, USAGE).unwrap();
        let mut log = String::new();
        let g = read_grids(&args, &settings, [edge; 3], false, &mut log).unwrap();
        assert_eq!(g.pbc, [38; 3]);
        assert_eq!(g.npbc, [78; 3]);
        assert!(log.contains("# Calculated number of gridpoints PBC (x,y,z): 38 38 38\n"));

        // @increasegrid grows the periodic grid only; @NPBCsize sets the other.
        let args = Arguments::from_vec(
            argv(&[
                "@gridspacing",
                "0.1",
                "@NPBCsize",
                "5.0",
                "@increasegrid",
                "2",
                "2",
                "2",
            ]),
            KNOWN,
            USAGE,
        )
        .unwrap();
        let s = Settings {
            increasegrid: [2; 3],
            ..settings
        };
        let g = read_grids(&args, &s, [edge; 3], false, &mut log).unwrap();
        assert_eq!(g.pbc, [40; 3]);
        assert_eq!(g.npbc, [50; 3]);

        // Without @gridspacing gromos++ would divide by zero.
        let args = Arguments::from_vec(argv(&["@nogridpoints", "40"]), KNOWN, USAGE).unwrap();
        assert!(read_grids(&args, &settings, [edge; 3], false, &mut log).is_err());
    }
}
