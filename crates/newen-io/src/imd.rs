//! Parser for GROMOS input parameter files (.imd)
//!
//! The IMD format uses block-based keywords to specify simulation parameters.
//!
//! # Format Example
//! ```text
//! TITLE
//!   My simulation
//! END
//! SYSTEM
//!   NPM     1
//!   NSM     0
//! END
//! STEP
//!   NSTLIM  1000
//!   T       0.0
//!   DT      0.002
//! END
//! ```
//!
//! # Major Blocks
//! - TITLE: Simulation description
//! - SYSTEM: System composition (NPM, NSM)
//! - STEP: Integration parameters (NSTLIM, DT, T)
//! - BOUNDCOND: Boundary conditions (NTB, NDFMIN)
//! - MULTIBATH: Temperature coupling groups
//! - PRESSURESCALE: Pressure coupling
//! - FORCE: Force field terms to compute
//! - CONSTRAINT: Constraint algorithm settings
//! - PAIRLIST: Neighbor list parameters
//! - NONBONDED: Cutoffs and long-range electrostatics
//! - INITIALISE: Initial velocities
//! - WRITETRAJ: Trajectory output settings
//! - PRINTOUT: Energy output settings

use crate::IoError;
use newen_core::lambdas::InteractionLambda;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// GROMOS simulation parameters from .imd file
#[derive(Debug, Clone, PartialEq)]
pub struct ImdParameters {
    /// Simulation title
    pub title: String,

    // SYSTEM block
    pub npm: usize, // Number of (identical) protein molecules
    pub nsm: usize, // Number of (identical) solvent molecules

    // STEP block
    pub nstlim: usize, // Number of MD steps
    pub t0: f64,       // Initial time (ps)
    pub dt: f64,       // Time step (ps)

    // BOUNDCOND block
    pub ntb: i32, // Boundary type (0=vacuum, 1=rectangular, 2=triclinic, -1=truncated octahedron)
    pub ndfmin: i32, // Minimum number of degrees of freedom

    // MULTIBATH block
    pub num_temp_baths: usize,
    pub temp_bath: Vec<TempBathParameters>,

    // PRESSURESCALE block
    pub couple_pressure: bool,
    pub pressure_parameters: Option<PressureParameters>,

    // FORCE block
    pub force_groups: Vec<Vec<(usize, usize)>>, // Energy group pairs

    // CONSTRAINT block
    pub ntc: i32,       // SHAKE constraints (1=none, 2=H-bonds, 3=all bonds, 4=all)
    pub ntcp: i32,      // P-SHAKE (pressure-SHAKE)
    pub ntcs: i32,      // Solvent SHAKE/SETTLE
    pub shake_tol: f64, // SHAKE tolerance
    pub lincs_order_solute: usize, // LINCS expansion order for solute (NTCP0 when NTCP=lincs)
    pub lincs_order_solvent: usize, // LINCS expansion order for solvent (NTCS0 when NTCS=lincs)

    // PAIRLIST block
    pub algorithm: i32, // Pairlist algorithm
    pub nsnb: usize,    // Update frequency
    pub rcutp: f64,     // Short-range cutoff (nm)
    pub rcutl: f64,     // Long-range cutoff (nm)
    pub size: f64,      // Grid cell size (nm)
    pub type_: i32,     // Pairlist type

    // COVALENTFORM block: which functional form each covalent term takes.
    /// NTBBH: 0 = quartic bonds (default), 1 = harmonic
    pub ntbbh: i32,
    /// NTBAH: 0 = cosine-harmonic angles (default), 1 = harmonic
    pub ntbah: i32,
    /// NTBDN: 0 = dihedrals with arbitrary phase shifts (default), 1 = phase shifts 0° or 180°
    pub ntbdn: i32,

    // NONBONDED block
    pub nlrele: i32,   // Long-range electrostatics (0=cutoff, 1=RF, 2=Ewald, 3=P3M)
    pub appak: f64,    // Reaction field κ (nm⁻¹)
    pub rcrf: f64,     // Reaction field cutoff (nm)
    pub epsrf: f64,    // Reaction field permittivity
    pub nslfexcl: i32, // Exclusions

    // FORCE block
    /// NTF force flags: [bonds, angles, improper, dihedral, charge, nonbonded]
    /// 0=off, 1=on for each term
    pub ntf: [i32; 6],
    /// Number of energy groups
    pub negr: usize,
    /// Last atom of each energy group (1-indexed)
    pub nre: Vec<usize>,

    // COMTRANSROT
    pub nscm: i32, // COM motion removal: >0 translation every NSCM steps, <0 translation+rotation every |NSCM| steps, 0=off

    // PME-specific parameters
    pub grid_x: usize,    // PME grid size X
    pub grid_y: usize,    // PME grid size Y
    pub grid_z: usize,    // PME grid size Z
    pub pme_order: usize, // PME spline order (typically 4)
    pub pme_alpha: f64,   // PME Ewald parameter

    // INITIALISE block
    pub ntivel: i32, // Initial velocities (0=read, 1=generate)
    pub ntishk: i32, // SHAKE initial configuration
    pub ntinht: i32, // Initial temperature assignment
    pub ntinhb: i32, // Initial bond constraints
    pub ntishi: i32, // Initial SHAKE iterations
    pub nticom: i32, // Initial COM motion removal (0=off, 1=translation, 2=trans+rot)
    pub ig: i64,     // Random seed
    pub tempi: f64,  // Initial temperature (K)

    // WRITETRAJ block
    pub ntwx: usize,        // Trajectory write frequency
    pub ntwe: usize,        // Energy write frequency
    pub ntwv: usize,        // Velocity-trajectory write frequency (0 = off)
    pub ntwf: usize,        // Force-trajectory write frequency (0 = off)
    pub ntwe_special: bool, // Special energy format

    // PRINTOUT block
    pub ntpr: usize, // Print frequency

    // POSITIONRES block
    /// Position restraints (0=off, 1=restrain with CPOR, 2=use B-factors, 3=constrain)
    pub ntpor: i32,
    /// Position restraint read B-factors (0=no, 1=yes)
    pub ntporb: i32,
    /// Position restraint scaling (0=off, 1=on)
    pub ntpors: i32,
    /// Position restraint force constant (kJ mol⁻¹ nm⁻²)
    pub cpor: f64,
    /// ROTTRANS: RTC — roto-translational constraints on (1) or off (0)
    pub rtc: i32,
    /// ROTTRANS: RTCLAST — last atom (1-based) of the constrained group
    pub rtclast: usize,

    // DISTANCERES block
    /// NTDIR: -2=time-avg+force, -1=time-avg, 0=off, 1=instantaneous, 2=instantaneous+w0
    pub ntdir: i32,
    /// NTDIRA: time-averaging flag
    pub ntdira: i32,
    /// CDIR: force constant for distance restraints (kJ mol⁻¹ nm⁻²)
    pub cdir: f64,
    /// DIR0: linear region threshold (nm)
    pub dir0: f64,
    /// TAUDIR: time constant for time averaging (ps)
    pub taudir: f64,
    /// FORCESCALE: force scaling mode for time-averaged restraints (0,1,2)
    pub forcescale: i32,
    /// VDIR: averaging exponent (unused in plain instantaneous mode)
    pub vdir: i32,
    /// NTWDIR: trajectory write frequency for restraints (0=off)
    pub ntwdir: i32,

    // PERTURBATION block (for FEP/lambda calculations)
    /// NTG: perturbation flag (0=off, 1=on)
    pub ntg: i32,
    /// NRDGL: read lambda from startup (0=no)
    pub nrdgl: i32,
    /// RLAM: starting lambda value
    pub rlam: f64,
    /// DLAMT: rate of λ change in ps⁻¹ (slow growth adds DLAMT·Δt per step)
    pub dlamt: f64,
    /// ALPHLJ: LJ soft-core alpha
    pub alphlj: f64,
    /// ALPHC: Coulomb soft-core alpha
    pub alphc: f64,
    /// NLAM: lambda exponent
    pub nlam: i32,
    /// NSCALE: energy group scaling (0=off)
    pub nscale: i32,
    /// LAMBDAS block (individual λ per interaction kind and energy-group pair); `None` = absent.
    pub lambdas: Option<LambdasBlock>,
    /// PRECALCLAM block (extended TI: energies at extra λ points); `None` = absent.
    pub precalclam: Option<PrecalcLamBlock>,
    /// EDS block (enveloping distribution sampling); `None` = absent.
    pub eds: Option<EdsBlock>,
    /// AEDS block (accelerated EDS); `None` = absent.
    pub aeds: Option<AedsBlock>,

    // ENERGYMIN block
    /// Energy minimization method (0=off, 1=steepest descent, 2=Fletcher-Reeves CG)
    pub ntem: i32,
    /// Minimum steps before convergence check
    pub nmin: usize,
    /// Energy convergence tolerance (kJ/mol)
    pub dele: f64,
    /// Initial step size (nm)
    pub dx0: f64,
    /// Maximum step size (nm)
    pub dxm: f64,
    /// Maximum force magnitude limit (0=unlimited)
    pub flim: f64,

    // Fields earlier parser versions dropped — kept so a read → write round trip is lossless
    // (the engine ignores most of them; gromosXX does not, and the writer must reproduce them).
    /// INITIALISE line 2: NTIRTC (read COM rotation/translation from the configuration)
    pub ntirtc: i32,
    /// INITIALISE line 3: NTISTI (stochastic-dynamics integrals)
    pub ntisti: i32,
    /// WRITETRAJ: NTWSE (selection trajectory), NTWG (free-energy), NTWB (block averages)
    pub ntwse: usize,
    pub ntwg: usize,
    pub ntwb: usize,
    /// PRINTOUT: NTPP (dihedral-angle transition printing)
    pub ntpp: usize,
    /// ENERGYMIN: NCYC (conjugate-gradient reset cycles)
    pub ncyc: usize,
    /// CONSTRAINT: NTCS0 when the solvent algorithm is SHAKE (its tolerance)
    pub ntcs0: f64,
    /// NONBONDED lines 3–7: lattice-sum settings (reaction-field runs ignore them)
    pub nonbonded_extra: NonbondedExtra,

    /// Raw blocks for custom parsing
    pub raw_blocks: HashMap<String, Vec<String>>,
}

/// NONBONDED block, lines 3–7 (NSHAPE … SLVDNS). Only relevant for lattice-sum
/// electrostatics; carried so the block round-trips.
#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedExtra {
    pub nshape: i32,
    pub ashape: f64,
    pub na2clc: i32,
    pub tola2: f64,
    pub epsls: f64,
    /// NK2 (the NKX/NKY/NKZ grid lives in `grid_x`/`grid_y`/`grid_z`)
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

impl Default for NonbondedExtra {
    fn default() -> Self {
        Self {
            nshape: -1,
            ashape: 1.4,
            na2clc: 2,
            tola2: 1e-10,
            epsls: 0.0,
            nk2: 100.0,
            ng: [32, 32, 32],
            nasord: 3,
            nfdord: 2,
            nalias: 3,
            nspord: 4,
            nqeval: 100000,
            faccur: 1.6,
            nrdgrd: 0,
            nwrgdr: 0,
            nlrlj: 0,
            slvdns: 33.3,
        }
    }
}

/// Temperature bath parameters (MULTIBATH block)
#[derive(Debug, Clone, PartialEq)]
pub struct TempBathParameters {
    pub algorithm: i32, // Coupling algorithm: 0=Berendsen, 1=NHC single, N>=2=NHC chain length N
    pub nhc_chain: usize, // NHC chain length (only used when algorithm >= 2)
    pub num_bath_groups: usize,
    pub temp0: Vec<f64>, // Reference temperatures (K)
    pub tau: Vec<f64>,   // Coupling times (ps)
    pub dof: Vec<usize>, // Degrees of freedom per group
    /// DOFSET lines: [LAST, COM-BATH, IR-BATH] per set (empty = synthesise one covering all atoms)
    pub dof_sets: Vec<[usize; 3]>,
}

/// Pressure bath parameters (PRESSURESCALE block)
#[derive(Debug, Clone, PartialEq)]
pub struct PressureParameters {
    /// COUPLE: 0=off, 1=calc (virial only), 2=scale (couple)
    pub couple: i32,
    pub algorithm: i32, // SCALE mode (0=off, 1=iso, 2=aniso, 3=full, 4=semianiso)
    /// SEMI: semi-anisotropic coupling flags per axis
    pub semi: [i32; 3],
    pub pressure0: [[f64; 3]; 3], // Reference pressure tensor (GPa/bar)
    pub compressibility: [[f64; 3]; 3], // Isothermal compressibility
    pub tau_p: f64,               // Coupling time (ps)
    pub virial: i32,              // Virial calculation method
}

impl Default for ImdParameters {
    fn default() -> Self {
        Self {
            title: String::from("GROMOS simulation"),
            npm: 1,
            nsm: 0,
            nstlim: 1000,
            t0: 0.0,
            dt: 0.002,
            ntb: 1,
            ndfmin: 0,
            // No MULTIBATH block = no temperature coupling (gromosXX adds an uncoupled bath).
            num_temp_baths: 0,
            temp_bath: Vec::new(),
            couple_pressure: false,
            pressure_parameters: None,
            force_groups: Vec::new(),
            ntc: 1,
            ntcp: 0,
            ntcs: 1,
            shake_tol: 1e-4,
            lincs_order_solute: 4,
            lincs_order_solvent: 4,
            algorithm: 0,
            nsnb: 5,
            rcutp: 0.8,
            rcutl: 1.4,
            size: 0.4,
            type_: 0,
            ntbbh: 0,
            ntbah: 0,
            ntbdn: 0,
            nlrele: 1,
            appak: 0.0,
            rcrf: 1.4,
            epsrf: 0.0,
            nslfexcl: 1,             // gromosXX parameter.h: `rf_excluded(true)`
            ntf: [1, 1, 1, 1, 1, 1], // All force terms on by default
            negr: 1,
            nre: Vec::new(),
            // No COMTRANSROT block = no COM-motion removal (gromosXX parameter.h: skip_step 0).
            nscm: 0,
            grid_x: 64,
            grid_y: 64,
            grid_z: 64,
            pme_order: 4,
            pme_alpha: 0.0,
            ntivel: 0,
            ntishk: 0,
            ntinht: 0,
            ntinhb: 0,
            ntishi: 1000,
            nticom: 0,
            ig: 12345,
            tempi: 300.0,
            ntwx: 100,
            ntwe: 100,
            ntwv: 0,
            ntwf: 0,
            ntwe_special: false,
            ntpr: 100,
            ntpor: 0,
            ntporb: 0,
            ntpors: 0,
            cpor: 0.0,
            ntdir: 0,
            ntdira: 0,
            cdir: 0.0,
            dir0: 0.0,
            taudir: 0.0,
            forcescale: 0,
            vdir: 0,
            ntwdir: 0,
            ntg: 0,
            nrdgl: 0,
            rlam: 0.0,
            dlamt: 0.0,
            alphlj: 0.0,
            alphc: 0.0,
            nlam: 1,
            nscale: 0,
            ntem: 0,
            nmin: 1,
            dele: 0.1,
            dx0: 0.01,
            dxm: 0.05,
            flim: 0.0,
            ntirtc: 0,
            rtc: 0,
            rtclast: 0,
            ntisti: 0,
            ntwse: 0,
            ntwg: 0,
            ntwb: 0,
            ntpp: 0,
            ncyc: 0,
            ntcs0: 1e-4,
            nonbonded_extra: NonbondedExtra::default(),
            lambdas: None,
            precalclam: None,
            eds: None,
            aeds: None,
            raw_blocks: HashMap::new(),
        }
    }
}

impl Default for TempBathParameters {
    fn default() -> Self {
        Self {
            algorithm: 0,
            nhc_chain: 1,
            num_bath_groups: 1,
            temp0: vec![300.0],
            tau: vec![0.1],
            dof: vec![0],
            dof_sets: Vec::new(),
        }
    }
}

/// Parse GROMOS .imd/.in parameter file.
///
/// Handles both key-value format and GROMOS positional format.
/// In positional format, comment lines starting with `#` describe fields,
/// and the following data line contains values in that order.
pub fn read_imd_file<P: AsRef<Path>>(path: P) -> Result<ImdParameters, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    parse_imd(reader)
}

/// Parse `.imd` text held in memory (same grammar as [`read_imd_file`]).
pub fn parse_imd_str(text: &str) -> Result<ImdParameters, IoError> {
    parse_imd(text.as_bytes())
}

/// Parse an `.imd` stream. A block whose numbers do not parse is an error naming the block
/// and the offending token — never a silent zero.
pub fn parse_imd<R: BufRead>(reader: R) -> Result<ImdParameters, IoError> {
    let mut params = ImdParameters::default();
    let mut current_block = String::new();
    let mut block_lines: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim().to_string();

        if trimmed.is_empty() {
            continue;
        }

        // Check for block start/end
        if trimmed == "END" {
            if !current_block.is_empty() {
                // Filter out comment lines for positional parsing
                let data_lines: Vec<String> = block_lines
                    .iter()
                    .filter(|l| !l.starts_with('#'))
                    .cloned()
                    .collect();
                parse_block(&mut params, &current_block, &data_lines).map_err(|e| match e {
                    IoError::ParseError(m) => {
                        IoError::ParseError(format!("block {current_block}: {m}"))
                    },
                    other => other,
                })?;
                params
                    .raw_blocks
                    .insert(current_block.clone(), block_lines.clone());
                block_lines.clear();
            }
            current_block.clear();
        } else if !current_block.is_empty() {
            // Inside a block - collect ALL lines (including comments for raw_blocks)
            block_lines.push(trimmed);
        } else if !trimmed.starts_with('#') {
            // New block starting (not a comment)
            current_block = trimmed;
        }
    }

    // LAMBDAS depends on PERTURBATION (NTG) and FORCE (NEGR), which may come later in the file;
    // gromosXX reads it after both (`in_parameter.cc:121`).
    if let Some(lines) = params.raw_blocks.get("LAMBDAS") {
        let data: Vec<String> = lines
            .iter()
            .filter(|l| !l.starts_with('#'))
            .cloned()
            .collect();
        params.lambdas = Some(
            parse_lambdas(&data, params.ntg, params.nre.len()).map_err(|e| match e {
                IoError::ParseError(m) => IoError::ParseError(format!("block LAMBDAS: {m}")),
                other => other,
            })?,
        );
    }

    Ok(params)
}

/// The `PRECALCLAM` block (gromosXX `read_PRECALCLAM`, `in_parameter.cc:4478-4513`): the
/// energies of the perturbed interactions at NRLAM extra λ points from MINLAM to MAXLAM.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrecalcLamBlock {
    /// NRLAM: 0 = off.
    pub nr_lambdas: usize,
    /// MINLAM
    pub min_lambda: f64,
    /// MAXLAM
    pub max_lambda: f64,
}

/// Parse the data lines of a `PRECALCLAM` block, with gromosXX's range checks.
pub fn parse_precalclam(lines: &[String]) -> Result<PrecalcLamBlock, IoError> {
    let mut t = Tokens::new("PRECALCLAM", lines);
    // gromosXX defaults (`parameter.h:2978`): 0, 0.0, 1.0
    let mut b = PrecalcLamBlock {
        nr_lambdas: 0,
        min_lambda: 0.0,
        max_lambda: 1.0,
    };
    if let Some(v) = t.next_i32("NRLAM")? {
        if v < 0 {
            return Err(IoError::ParseError(format!("NRLAM = {v}: must be >= 0")));
        }
        b.nr_lambdas = v as usize;
    }
    for (name, slot) in [("MINLAM", &mut b.min_lambda), ("MAXLAM", &mut b.max_lambda)] {
        if let Some(v) = t.next_f64(name)? {
            if !(0.0..=1.0).contains(&v) {
                return Err(IoError::ParseError(format!(
                    "{name} = {v}: must be >= 0 and <= 1"
                )));
            }
            *slot = v;
        }
    }
    if b.min_lambda >= b.max_lambda {
        return Err(IoError::ParseError(
            "MINLAM should be smaller than MAXLAM".to_string(),
        ));
    }
    Ok(b)
}

/// One FORM 3 pair of an `EDS` block: end states `i` and `j` (1-based, as in the file) and
/// their smoothness parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdsStatePair {
    pub i: i32,
    pub j: i32,
    pub s: f64,
}

/// The `EDS` block (gromosXX `read_EDS`, `in_parameter.cc:3843-3961`).
#[derive(Debug, Clone, PartialEq)]
pub struct EdsBlock {
    /// EDS: 0 off, 1 on.
    pub enabled: bool,
    /// FORM: 1 one smoothness parameter, 2 one per state pair, 3 one per listed pair.
    pub form: i32,
    /// NUMSTATES (as read; gromosXX sets it to 0 when EDS is off).
    pub num_states: usize,
    /// S: one value (FORM 1) or N(N−1)/2 in upper-triangle row-major pair order (FORM 2).
    pub s: Vec<f64>,
    /// FORM 3: the N−1 rows `i j S`.
    pub pairs: Vec<EdsStatePair>,
    /// EIR: the energy offset of every state.
    pub eir: Vec<f64>,
}

/// Parse the data lines of an `EDS` block as gromosXX does: a stream of tokens EDS, FORM,
/// NUMSTATES (FORM in 1..3 and NUMSTATES ≥ 2 even with EDS = 0), then for EDS = 1 the S values
/// of the form and NUMSTATES offsets EIR. Left-over tokens are an error (`get_final_messages`,
/// which with EDS = 0 runs right after NUMSTATES).
pub fn parse_eds(lines: &[String]) -> Result<EdsBlock, IoError> {
    // gromosXX's `io::getline` drops everything after a `#` (blockinput.cc:37-62)
    let lines: Vec<String> = lines
        .iter()
        .map(|l| l.split('#').next().unwrap_or_default().trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    fn int(t: &mut Tokens<'_>, name: &str) -> Result<i32, IoError> {
        t.next_i32(name)?
            .ok_or_else(|| IoError::ParseError(format!("missing {name}")))
    }
    fn real(t: &mut Tokens<'_>, name: &str) -> Result<f64, IoError> {
        t.next_f64(name)?
            .ok_or_else(|| IoError::ParseError(format!("missing {name}")))
    }
    let mut t = Tokens::new("EDS", &lines);
    let err = IoError::ParseError;
    let enabled = match int(&mut t, "EDS")? {
        0 => false,
        1 => true,
        v => return Err(err(format!("EDS = {v}: must be 0 or 1"))),
    };
    let form = int(&mut t, "FORM")?;
    if !(1..=3).contains(&form) {
        return Err(err(format!("FORM = {form}: must be 1, 2 or 3")));
    }
    let n = int(&mut t, "NUMSTATES")?;
    if n < 2 {
        return Err(err(format!("NUMSTATES = {n}: must be >= 2")));
    }
    let n = n as usize;
    let leftover = |t: &mut Tokens<'_>| match t.next_raw() {
        Some(extra) => Err(IoError::ParseError(format!(
            "Left-over parameters in EDS block (from {extra:?})"
        ))),
        None => Ok(()),
    };
    let mut b = EdsBlock {
        enabled,
        form,
        num_states: n,
        s: Vec::new(),
        pairs: Vec::new(),
        eir: Vec::new(),
    };
    if !enabled {
        leftover(&mut t)?;
        return Ok(b);
    }
    let positive = |v: f64, name: &str| {
        if v > 0.0 {
            Ok(v)
        } else {
            Err(IoError::ParseError(format!("{name} = {v}: must be > 0")))
        }
    };
    match form {
        1 => b.s.push(positive(real(&mut t, "S[0]")?, "S[0]")?),
        2 => {
            for k in 0..n * (n - 1) / 2 {
                let name = format!("S[{k}]");
                b.s.push(positive(real(&mut t, &name)?, &name)?);
            }
        },
        _ => {
            for k in 0..n - 1 {
                let i = int(&mut t, &format!("i[{k}]"))?;
                let j = int(&mut t, &format!("j[{k}]"))?;
                let name = format!("S[{k}]");
                let s = positive(real(&mut t, &name)?, &name)?;
                b.pairs.push(EdsStatePair { i, j, s });
            }
        },
    }
    for k in 0..n {
        b.eir.push(real(&mut t, &format!("EIR[{k}]"))?);
    }
    leftover(&mut t)?;
    Ok(b)
}

/// The `AEDS` block (gromosXX `read_AEDS`, `in_parameter.cc:3968-4140`). With AEDS = 0 only
/// FORM and NUMSTATES are read; the other fields are then zero.
#[derive(Debug, Clone, PartialEq)]
pub struct AedsBlock {
    /// AEDS: 0 off, 1 on.
    pub enabled: bool,
    /// FORM 1-6: 1 fixed parameters; 2 search E^R; 3 search E_max/E_min; 4 search all; 5 and 6
    /// the "advanced" offset search with (5) or without (6) the E_max/E_min search.
    pub form: i32,
    /// NUMSTATES
    pub num_states: usize,
    /// EMAX (kJ/mol)
    pub emax: f64,
    /// EMIN (kJ/mol)
    pub emin: f64,
    /// EIR: the energy offset of every state.
    pub eir: Vec<f64>,
    /// NTIAEDSS: 1 start a fresh search; 0 continue the configuration's AEDSSEARCH block.
    pub ntiaedss: bool,
    /// RESTREMIN: 1 keep E_min at or above the lowest average end-state energy until every
    /// state has been visited.
    pub restremin: bool,
    /// BMAXTYPE: 1 BMAX in kJ/mol, 2 BMAX in standard deviations.
    pub bmaxtype: i32,
    /// BMAX: the maximum energy barrier parameter.
    pub bmax: f64,
    /// ASTEPS: memory time (steps) at the start of the run.
    pub asteps: u32,
    /// BSTEPS: memory time (steps) at the end of the run.
    pub bsteps: u32,
}

/// Parse the data lines of an `AEDS` block as gromosXX does: the token stream AEDS, FORM,
/// NUMSTATES (FORM in 1..6 and NUMSTATES ≥ 2 even with AEDS = 0, which stops there), then EMAX,
/// EMIN, NUMSTATES offsets EIR, NTIAEDSS, RESTREMIN, BMAXTYPE, BMAX, ASTEPS, BSTEPS. Left-over
/// tokens are an error. EMIN > EMAX parses: gromosXX only warns and then reads nothing more
/// (`:4084-4088`); the recipe refuses it.
pub fn parse_aeds(lines: &[String]) -> Result<AedsBlock, IoError> {
    let lines: Vec<String> = lines
        .iter()
        .map(|l| l.split('#').next().unwrap_or_default().trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    fn int(t: &mut Tokens<'_>, name: &str) -> Result<i32, IoError> {
        t.next_i32(name)?
            .ok_or_else(|| IoError::ParseError(format!("missing {name}")))
    }
    fn real(t: &mut Tokens<'_>, name: &str) -> Result<f64, IoError> {
        t.next_f64(name)?
            .ok_or_else(|| IoError::ParseError(format!("missing {name}")))
    }
    fn one_of(v: i32, name: &str, allowed: &[i32]) -> Result<i32, IoError> {
        if allowed.contains(&v) {
            Ok(v)
        } else {
            Err(IoError::ParseError(format!(
                "AEDS {name} = {v}: must be one of {allowed:?}"
            )))
        }
    }
    let mut t = Tokens::new("AEDS", &lines);
    let enabled = one_of(int(&mut t, "AEDS")?, "AEDS", &[0, 1])? == 1;
    let form = one_of(int(&mut t, "FORM")?, "FORM", &[1, 2, 3, 4, 5, 6])?;
    let n = int(&mut t, "NUMSTATES")?;
    if n < 2 {
        return Err(IoError::ParseError(format!(
            "AEDS NUMSTATES = {n}: must be >= 2"
        )));
    }
    let leftover = |t: &mut Tokens<'_>| match t.next_raw() {
        Some(extra) => Err(IoError::ParseError(format!(
            "Left-over parameters in AEDS block (from {extra:?})"
        ))),
        None => Ok(()),
    };
    let mut b = AedsBlock {
        enabled,
        form,
        num_states: n as usize,
        emax: 0.0,
        emin: 0.0,
        eir: Vec::new(),
        ntiaedss: false,
        restremin: false,
        bmaxtype: 0,
        bmax: 0.0,
        asteps: 0,
        bsteps: 0,
    };
    if !enabled {
        leftover(&mut t)?;
        return Ok(b);
    }
    b.emax = real(&mut t, "EMAX")?;
    b.emin = real(&mut t, "EMIN")?;
    for k in 0..b.num_states {
        b.eir.push(real(&mut t, &format!("EIR[{k}]"))?);
    }
    b.ntiaedss = one_of(int(&mut t, "NTIAEDSS")?, "NTIAEDSS", &[0, 1])? == 1;
    b.restremin = one_of(int(&mut t, "RESTREMIN")?, "RESTREMIN", &[0, 1])? == 1;
    b.bmaxtype = one_of(int(&mut t, "BMAXTYPE")?, "BMAXTYPE", &[1, 2])?;
    b.bmax = real(&mut t, "BMAX")?;
    if b.bmax.is_nan() || b.bmax <= 0.0 {
        return Err(IoError::ParseError(format!(
            "AEDS BMAX = {}: must be > 0",
            b.bmax
        )));
    }
    for (name, slot) in [("ASTEPS", &mut b.asteps), ("BSTEPS", &mut b.bsteps)] {
        let v = int(&mut t, name)?;
        if v <= 0 {
            return Err(IoError::ParseError(format!(
                "AEDS {name} = {v}: must be > 0"
            )));
        }
        *slot = v as u32;
    }
    leftover(&mut t)?;
    Ok(b)
}

/// One `LAMBDAS` row: `L = a·λ⁴ + b·λ³ + c·λ² + d·λ + e` for `kind` between energy groups
/// `g1` and `g2` (1-based, as in the file, `g1 ≤ g2`).
#[derive(Debug, Clone, PartialEq)]
pub struct LambdasRow {
    pub kind: InteractionLambda,
    pub g1: usize,
    pub g2: usize,
    /// `[ALI, BLI, CLI, DLI, ELI]`
    pub coefficients: [f64; 5],
}

/// The `LAMBDAS` block (gromosXX `read_LAMBDAS`, `in_parameter.cc:4326-4472`).
///
/// The rows act only when the block is on *and* the run is perturbed (NTG ≠ 0); otherwise
/// gromosXX never reads them and every individual λ stays the identity. They are kept anyway
/// (when they parse) so the block round-trips.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LambdasBlock {
    /// NTIL: `on`/`1` or `off`/`0`.
    pub enabled: bool,
    pub rows: Vec<LambdasRow>,
}

impl LambdasBlock {
    /// Whether the rows set the individual λ of a run with this NTG.
    pub fn applies(&self, ntg: i32) -> bool {
        self.enabled && ntg != 0
    }

    /// The warnings gromosXX prints for this block (it runs anyway).
    pub fn warnings(&self, ntg: i32) -> Vec<String> {
        let mut out = Vec::new();
        if self.enabled && ntg == 0 {
            out.push("LAMBDAS block without perturbation is ignored".to_string());
            return out;
        }
        if !self.enabled {
            return out;
        }
        for r in &self.rows {
            if !r.kind.is_nonbonded() && r.g1 != r.g2 {
                out.push(format!(
                    "NILG1 != NILG2 in LAMBDAS block only allowed for nonbonded interactions \
                     ({} {} {}: never read, the term uses its own group's diagonal)",
                    r.kind.gromos_name(),
                    r.g1,
                    r.g2
                ));
            }
        }
        out
    }
}

/// The `PRECALCLAM` block (gromosXX `read_PRECALCLAM`, `in_parameter.cc:4478-4515`): energies
/// and dH/dλ precalculated at `nrlam` λ values evenly spaced over `[minlam, maxlam]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Precalclam {
    pub nrlam: usize,
    pub minlam: f64,
    pub maxlam: f64,
}

impl ImdParameters {
    /// The `PRECALCLAM` block, `None` when the file has none. Read from the stored block text
    /// with gromosXX's checks (NRLAM ≥ 0, MINLAM and MAXLAM in [0, 1], MINLAM < MAXLAM).
    pub fn precalclam(&self) -> Result<Option<Precalclam>, IoError> {
        let Some(lines) = self.raw_blocks.get("PRECALCLAM") else {
            return Ok(None);
        };
        let data: Vec<String> = lines
            .iter()
            .map(|l| l.split('#').next().unwrap_or_default().trim().to_string())
            .filter(|l| !l.is_empty())
            .collect();
        let mut t = Tokens::new("PRECALCLAM", &data);
        let missing = |name: &str| IoError::ParseError(format!("block PRECALCLAM: missing {name}"));
        let nrlam = t.next_i32("NRLAM")?.ok_or_else(|| missing("NRLAM"))?;
        let minlam = t.next_f64("MINLAM")?.ok_or_else(|| missing("MINLAM"))?;
        let maxlam = t.next_f64("MAXLAM")?.ok_or_else(|| missing("MAXLAM"))?;
        let bad = |m: String| Err(IoError::ParseError(format!("block PRECALCLAM: {m}")));
        if nrlam < 0 {
            return bad(format!("NRLAM must be >= 0, found {nrlam}"));
        }
        for (name, v) in [("MINLAM", minlam), ("MAXLAM", maxlam)] {
            if !(0.0..=1.0).contains(&v) {
                return bad(format!("{name} must be in [0, 1], found {v}"));
            }
        }
        if minlam >= maxlam {
            return bad("MINLAM should be smaller than MAXLAM".to_string());
        }
        Ok(Some(Precalclam {
            nrlam: nrlam as usize,
            minlam,
            maxlam,
        }))
    }
}

/// Parse the data lines of a `LAMBDAS` block. `ntg` is PERTURBATION's NTG and `n_groups` the
/// number of energy groups FORCE declares (NEGR, before gromosXX's implicit trailing group).
pub fn parse_lambdas(lines: &[String], ntg: i32, n_groups: usize) -> Result<LambdasBlock, IoError> {
    // gromosXX's `io::getline` drops everything after a `#` (blockinput.cc:37-62)
    let lines: Vec<String> = lines
        .iter()
        .map(|l| l.split('#').next().unwrap_or_default().trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    let mut t = Tokens::new("LAMBDAS", &lines);
    let enabled = match t.next_raw().as_deref() {
        Some("on") | Some("1") => true,
        Some("off") | Some("0") => false,
        Some(other) => {
            return Err(IoError::ParseError(format!(
                "NTIL must be on, off, 1 or 0, found {other:?}"
            )))
        },
        None => return Err(IoError::ParseError("missing NTIL".to_string())),
    };
    // gromosXX: `num = block.numlines() - 3`, one row per line after NTIL's
    let num = lines.len().saturating_sub(1);
    let rows = read_lambdas_rows(&mut t, num, n_groups);
    if enabled && ntg != 0 {
        let rows = rows?;
        // `block.get_final_messages()` (in_parameter.cc:4470, blockinput.cc:281-286)
        if let Some(extra) = t.next_raw() {
            return Err(IoError::ParseError(format!(
                "Left-over parameters in LAMBDAS block (from {extra:?}): one row per line"
            )));
        }
        return Ok(LambdasBlock { enabled, rows });
    }
    // off, or on without perturbation: gromosXX returns before reading a row
    Ok(LambdasBlock {
        enabled,
        rows: rows.unwrap_or_default(),
    })
}

fn read_lambdas_rows(
    t: &mut Tokens<'_>,
    num: usize,
    n_groups: usize,
) -> Result<Vec<LambdasRow>, IoError> {
    let mut rows = Vec::with_capacity(num);
    for i in 0..num {
        let missing = |name: &str| IoError::ParseError(format!("missing {name}[{i}]"));
        let ntli = t.next_raw().ok_or_else(|| missing("NTLI"))?;
        let kind = InteractionLambda::from_ntli(&ntli).ok_or_else(|| {
            IoError::ParseError(format!(
                "NTLI[{i}] = {ntli:?} is not one of 1..12, bond, angle, dihedral, improper, \
                 vdw, vdw_soft, crf, crf_soft, distanceres, distancefield, dihedralres, mass"
            ))
        })?;
        let mut group = |name: &str| -> Result<usize, IoError> {
            let v = t
                .next_i32(&format!("{name}[{i}]"))?
                .ok_or_else(|| missing(name))?;
            if v <= 0 {
                return Err(IoError::ParseError(format!(
                    "{name}[{i}] must be > 0, found {v}"
                )));
            }
            Ok(v as usize)
        };
        let g1 = group("NILG1")?;
        let g2 = group("NILG2")?;
        let mut coefficients = [0.0; 5];
        for (c, name) in coefficients
            .iter_mut()
            .zip(["ALI", "BLI", "CLI", "DLI", "ELI"])
        {
            *c = t
                .next_f64(&format!("{name}[{i}]"))?
                .ok_or_else(|| missing(name))?;
        }
        if g1 > n_groups || g2 > n_groups {
            return Err(IoError::ParseError(format!(
                "NILG1 and NILG2 need to be smaller than the number of energy groups \
                 (row {kind} {g1} {g2}, NEGR = {n_groups})",
                kind = kind.gromos_name()
            )));
        }
        if g2 < g1 {
            return Err(IoError::ParseError(format!(
                "only give NILG2 >= NILG1 in LAMBDAS BLOCK (row {} {g1} {g2})",
                kind.gromos_name()
            )));
        }
        rows.push(LambdasRow {
            kind,
            g1,
            g2,
            coefficients,
        });
    }
    Ok(rows)
}

/// Parse a specific IMD block using GROMOS positional format.
///
/// Data lines are the non-comment lines within the block, in order.
/// The format follows GROMOS conventions where each data line
/// corresponds to a group of parameters described by preceding comments.
fn parse_block(
    params: &mut ImdParameters,
    block_name: &str,
    data_lines: &[String],
) -> Result<(), IoError> {
    match block_name {
        "TITLE" => {
            // Keep the lines: the title is echoed verbatim into every output file.
            params.title = data_lines.join("\n");
        },
        "SYSTEM" => {
            let mut t = Tokens::new("SYSTEM", data_lines);
            if let Some(v) = t.next_usize("NPM")? {
                params.npm = v;
            }
            if let Some(v) = t.next_usize("NSM")? {
                params.nsm = v;
            }
        },
        "STEP" => {
            let mut t = Tokens::new("STEP", data_lines);
            if let Some(v) = t.next_usize("NSTLIM")? {
                params.nstlim = v;
            }
            if let Some(v) = t.next_f64("T")? {
                params.t0 = v;
            }
            if let Some(v) = t.next_f64("DT")? {
                params.dt = v;
            }
        },
        "BOUNDCOND" => {
            let mut t = Tokens::new("BOUNDCOND", data_lines);
            if let Some(v) = t.next_i32("NTB")? {
                params.ntb = v;
            }
            if let Some(v) = t.next_i32("NDFMIN")? {
                params.ndfmin = v;
            }
        },
        "MULTIBATH" => {
            // gromosXX reads this block as a value stream (`read_MULTIBATH`): the algorithm name
            // (with the chain length right after it for nose-hoover-chains), then NBATHS, then
            // TEMP0/TAU per bath, then DOFSET and LAST/COMBATH/IRBATH per set. Line breaks carry
            // no meaning, so a file may put the whole block on one line.
            let mut bath = TempBathParameters::default();
            let mut t = Tokens::new("MULTIBATH", data_lines);
            let Some(algorithm) = t.next_raw() else {
                return Ok(());
            };
            match algorithm.as_str() {
                "weak-coupling" => bath.algorithm = 0,
                "nose-hoover" => bath.algorithm = 1,
                "nose-hoover-chains" => {
                    // the chain length follows the name
                    let n = t.next_usize("NUM")?.unwrap_or(3).max(2);
                    bath.algorithm = n as i32;
                    bath.nhc_chain = n;
                },
                other => bath.algorithm = parse_i32(other)?,
            }
            let Some(nbaths) = t.next_usize("NBATHS")? else {
                return Err(IoError::ParseError(
                    "block MULTIBATH: missing value for NBATHS".to_string(),
                ));
            };
            bath.num_bath_groups = nbaths;
            bath.temp0.clear();
            bath.tau.clear();
            for _ in 0..nbaths {
                match (t.next_f64("TEMP0")?, t.next_f64("TAU")?) {
                    (Some(temp0), Some(tau)) => {
                        bath.temp0.push(temp0);
                        bath.tau.push(tau);
                    },
                    _ => {
                        return Err(IoError::ParseError(format!(
                            "block MULTIBATH: NBATHS = {nbaths} but fewer TEMP0/TAU pairs follow"
                        )))
                    },
                }
            }
            let n_sets = t.next_usize("DOFSET")?.unwrap_or(0);
            for _ in 0..n_sets {
                match (
                    t.next_usize("LAST")?,
                    t.next_usize("COMBATH")?,
                    t.next_usize("IRBATH")?,
                ) {
                    (Some(last), Some(com), Some(ir)) => bath.dof_sets.push([last, com, ir]),
                    _ => {
                        return Err(IoError::ParseError(format!(
                            "block MULTIBATH: DOFSET = {n_sets} but fewer LAST/COMBATH/IRBATH \
                             triplets follow"
                        )))
                    },
                }
            }
            bath.dof = vec![0; bath.num_bath_groups.max(1)];
            params.num_temp_baths = bath.num_bath_groups;
            params.temp_bath = vec![bath];
        },
        "PRESSURESCALE" => {
            // GROMOS PRESSURESCALE block format:
            //   Line 0: COUPLE SCALE COMP TAUP VIRIAL
            //     COUPLE: off/calc/scale (or 0/1/2)
            //     SCALE: off/iso/aniso/full/semianiso (or 0/1/2/3/4)
            //     COMP: compressibility value
            //     TAUP: pressure coupling time
            //     VIRIAL: none/atomic/molecular (or 0/1/2)
            //   Line 1: SEMI (3 values, for semianisotropic)
            //   Lines 2-4: reference pressure (3x3)
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if !v.is_empty() {
                    // Parse COUPLE keyword/number
                    let couple = match v[0].to_lowercase().as_str() {
                        "off" => 0,
                        "calc" => 1,
                        "scale" => 2,
                        other => parse_i32(other)?,
                    };
                    params.couple_pressure = couple > 0;
                    if couple > 0 {
                        // Parse SCALE keyword
                        let scale = if v.len() >= 2 {
                            match v[1].to_lowercase().as_str() {
                                "off" => 0,
                                "iso" | "isotropic" => 1,
                                "aniso" | "anisotropic" => 2,
                                "full" | "full_anisotropic" => 3,
                                "semianiso" | "semi_anisotropic" => 4,
                                other => parse_i32(other)?,
                            }
                        } else {
                            0
                        };
                        // Parse COMP (compressibility)
                        let comp = if v.len() >= 3 {
                            parse_f64(&v[2])?
                        } else {
                            4.575e-4
                        };
                        // Parse TAUP
                        let tau_p = if v.len() >= 4 { parse_f64(&v[3])? } else { 0.5 };
                        // Parse VIRIAL keyword
                        let virial = if v.len() >= 5 {
                            match v[4].to_lowercase().as_str() {
                                "none" => 0,
                                "atomic" => 1,
                                "molecular" => 2,
                                other => parse_i32(other)?,
                            }
                        } else {
                            0
                        };

                        let mut pp = PressureParameters {
                            couple,
                            algorithm: scale, // store SCALE mode (iso/aniso/etc.)
                            semi: [1, 1, 1],
                            pressure0: [[0.0; 3]; 3],
                            compressibility: [[comp; 3]; 3],
                            tau_p,
                            virial,
                        };

                        // Line 1: SEMI (semianisotropic coupling flags)
                        let mut dl = 1;
                        if dl < data_lines.len() {
                            let sv = parse_values(&data_lines[dl]);
                            if sv.len() >= 3 {
                                pp.semi =
                                    [parse_i32(&sv[0])?, parse_i32(&sv[1])?, parse_i32(&sv[2])?];
                            }
                            dl += 1;
                        }

                        // Lines 2-4: reference pressure (3x3)
                        for row in 0..3 {
                            if dl < data_lines.len() {
                                let pv = parse_values(&data_lines[dl]);
                                for (col, tok) in pv.iter().take(3).enumerate() {
                                    pp.pressure0[row][col] = parse_f64(tok)?;
                                }
                                dl += 1;
                            }
                        }
                        params.pressure_parameters = Some(pp);
                    }
                }
            }
        },
        "CONSTRAINT" => {
            // gromosXX `read_CONSTRAINT`: NTC, then the solute algorithm and *its* options, then
            // the solvent algorithm and its options — a value stream, so the file may put all of
            // them on one line (the GROMOS tutorials) or one per line (our reference inputs).
            // Which options follow depends on the algorithm: SETTLE reads none, LINCS reads an
            // expansion order, SHAKE reads a tolerance.
            let mut t = Tokens::new("CONSTRAINT", data_lines);
            if let Some(ntc) = t.next_raw() {
                params.ntc = match ntc.as_str() {
                    "off" => 0,
                    "solvent" => 1,
                    "hydrogen" => 2,
                    "all" => 3,
                    "specified" => 4,
                    _ => parse_i32(&ntc)?,
                };
            }
            if let Some(alg) = t.next_raw() {
                params.ntcp = constraint_algorithm_code("NTCP", &alg)?;
                match params.ntcp {
                    0 | 3 => {},
                    2 => {
                        if let Some(v) = t.next_usize("NTCP[0]")? {
                            params.lincs_order_solute = v;
                        }
                    },
                    _ => {
                        if let Some(v) = t.next_f64("NTCP[0]")? {
                            params.shake_tol = v;
                        }
                    },
                }
            }
            if let Some(alg) = t.next_raw() {
                params.ntcs = constraint_algorithm_code("NTCS", &alg)?;
                match params.ntcs {
                    0 | 3 => {},
                    2 => {
                        if let Some(v) = t.next_usize("NTCS[0]")? {
                            params.lincs_order_solvent = v;
                        }
                    },
                    _ => {
                        if let Some(v) = t.next_f64("NTCS[0]")? {
                            params.ntcs0 = v;
                        }
                    },
                }
            }
        },
        "FORCE" => {
            // NTF(1..6), then NEGR and NEGR energy-group boundaries — a value stream, so the
            // six flags and the group list may be wrapped anywhere.
            let mut t = Tokens::new("FORCE", data_lines);
            for i in 0..6 {
                if let Some(v) = t.next_i32("NTF")? {
                    params.ntf[i] = v;
                }
            }
            if let Some(negr) = t.next_usize("NEGR")? {
                params.negr = negr;
                let mut nre = Vec::with_capacity(negr);
                for _ in 0..negr {
                    match t.next_usize("NRE")? {
                        // gromosXX in_parameter.cc:978-993: each NRE > 0, strictly increasing
                        Some(v) if v == 0 || nre.last().is_some_and(|&p| v <= p) => {
                            return Err(IoError::ParseError(format!(
                                "block FORCE: NRE values must be positive and strictly \
                                 increasing, found {v} after {:?}",
                                nre
                            )));
                        },
                        Some(v) => nre.push(v),
                        None => break,
                    }
                }
                params.nre = nre;
            }
        },
        "PAIRLIST" => {
            let mut t = Tokens::new("PAIRLIST", data_lines);
            if let Some(raw) = t.next_raw() {
                // ALGORITHM: faithful to gromosXX in_parameter.cc:1419-1422
                //   "standard" / 0 → Standard_Pairlist_Algorithm
                //   "grid"     / 1 → Extended_Grid_Pairlist_Algorithm (production default)
                //   "grid_cell"/ 2 → Grid_Cell_Pairlist (Heinz & Hünenberger 2004)
                params.algorithm = match raw.as_str() {
                    "standard" => 0,
                    "grid" => 1,
                    "grid_cell" => 2,
                    _ => parse_i32(&raw)?,
                };
            }
            if let Some(v) = t.next_usize("NSNB")? {
                params.nsnb = v;
            }
            if let Some(v) = t.next_f64("RCUTP")? {
                params.rcutp = v;
            }
            if let Some(v) = t.next_f64("RCUTL")? {
                params.rcutl = v;
            }
            if let Some(raw) = t.next_raw() {
                // SIZE can be "auto" or a number
                params.size = match raw.as_str() {
                    "auto" => 0.0,
                    _ => parse_f64(&raw)?,
                };
            }
            if let Some(raw) = t.next_raw() {
                params.type_ = match raw.as_str() {
                    "chargegroup" => 0,
                    "atomic" => 1,
                    _ => parse_i32(&raw)?,
                };
            }
        },
        "NONBONDED" => {
            // gromosXX reads this block as a stream of values (`In_Parameter::read_NONBONDED`,
            // one `get_next_parameter` per name), so the line wrapping of the file is irrelevant:
            // real inputs put NLRELE on the same line as APPAK/RCRF/EPSRF/NSLFEXCL, ours put it
            // on its own. Values missing at the end of the block keep their defaults.
            let mut t = Tokens::new("NONBONDED", data_lines);
            if let Some(v) = t.next_i32("NLRELE")? {
                params.nlrele = v;
            }
            if let Some(v) = t.next_f64("APPAK")? {
                params.appak = v;
            }
            if let Some(v) = t.next_f64("RCRF")? {
                params.rcrf = v;
            }
            if let Some(v) = t.next_f64("EPSRF")? {
                params.epsrf = v;
            }
            if let Some(v) = t.next_i32("NSLFEXCL")? {
                params.nslfexcl = v;
            }
            // Lattice-sum settings: unused by a reaction-field run, carried for the round trip.
            let x = &mut params.nonbonded_extra;
            if let Some(v) = t.next_i32("NSHAPE")? {
                x.nshape = v;
            }
            if let Some(v) = t.next_f64("ASHAPE")? {
                x.ashape = v;
            }
            if let Some(v) = t.next_i32("NA2CLC")? {
                x.na2clc = v;
            }
            if let Some(v) = t.next_f64("TOLA2")? {
                x.tola2 = v;
            }
            if let Some(v) = t.next_f64("EPSLS")? {
                x.epsls = v;
            }
            for (name, slot) in [
                ("NKX", &mut params.grid_x),
                ("NKY", &mut params.grid_y),
                ("NKZ", &mut params.grid_z),
            ] {
                if let Some(v) = t.next_usize(name)? {
                    *slot = v;
                }
            }
            let x = &mut params.nonbonded_extra;
            if let Some(v) = t.next_f64("KCUT")? {
                x.nk2 = v;
            }
            for (name, k) in [("NGX", 0), ("NGY", 1), ("NGZ", 2)] {
                if let Some(v) = t.next_usize(name)? {
                    x.ng[k] = v;
                }
            }
            for (name, slot) in [
                ("NASORD", &mut x.nasord),
                ("NFDORD", &mut x.nfdord),
                ("NALIAS", &mut x.nalias),
                ("NSPORD", &mut x.nspord),
                ("NQEVAL", &mut x.nqeval),
            ] {
                if let Some(v) = t.next_i32(name)? {
                    *slot = v;
                }
            }
            if let Some(v) = t.next_f64("FACCUR")? {
                x.faccur = v;
            }
            for (name, slot) in [
                ("NRDGRD", &mut x.nrdgrd),
                ("NWRGRD", &mut x.nwrgdr),
                ("NLRLJ", &mut x.nlrlj),
            ] {
                if let Some(v) = t.next_i32(name)? {
                    *slot = v;
                }
            }
            if let Some(v) = t.next_f64("SLVDNS")? {
                x.slvdns = v;
            }
        },
        "INITIALISE" => {
            let mut t = Tokens::new("INITIALISE", data_lines);
            for (name, slot) in [
                ("NTIVEL", &mut params.ntivel),
                ("NTISHK", &mut params.ntishk),
                ("NTINHT", &mut params.ntinht),
                ("NTINHB", &mut params.ntinhb),
                ("NTISHI", &mut params.ntishi),
                ("NTIRTC", &mut params.ntirtc),
                ("NTICOM", &mut params.nticom),
                ("NTISTI", &mut params.ntisti),
            ] {
                if let Some(v) = t.next_i32(name)? {
                    *slot = v;
                }
            }
            if let Some(raw) = t.next_raw() {
                params.ig = raw.parse::<i64>().map_err(|_| {
                    IoError::ParseError(format!(
                        "block INITIALISE: IG expects an integer seed, found {raw:?}"
                    ))
                })?;
            }
            if let Some(v) = t.next_f64("TEMPI")? {
                params.tempi = v;
            }
        },
        "WRITETRAJ" => {
            let mut t = Tokens::new("WRITETRAJ", data_lines);
            for (name, slot) in [
                ("NTWX", &mut params.ntwx),
                ("NTWSE", &mut params.ntwse),
                ("NTWV", &mut params.ntwv),
                ("NTWF", &mut params.ntwf),
                ("NTWE", &mut params.ntwe),
                ("NTWG", &mut params.ntwg),
                ("NTWB", &mut params.ntwb),
            ] {
                if let Some(v) = t.next_usize(name)? {
                    *slot = v;
                }
            }
        },
        "PRINTOUT" => {
            let mut t = Tokens::new("PRINTOUT", data_lines);
            if let Some(v) = t.next_usize("NTPR")? {
                params.ntpr = v;
            }
            if let Some(v) = t.next_usize("NTPP")? {
                params.ntpp = v;
            }
        },
        "COMTRANSROT" => {
            let mut t = Tokens::new("COMTRANSROT", data_lines);
            if let Some(v) = t.next_i32("NSCM")? {
                params.nscm = v;
            }
        },
        "ROTTRANS" => {
            let mut t = Tokens::new("ROTTRANS", data_lines);
            if let Some(v) = t.next_i32("RTC")? {
                params.rtc = v;
            }
            if let Some(v) = t.next_usize("RTCLAST")? {
                params.rtclast = v;
            }
        },
        "POSITIONRES" => {
            let mut t = Tokens::new("POSITIONRES", data_lines);
            for (name, slot) in [
                ("NTPOR", &mut params.ntpor),
                ("NTPORB", &mut params.ntporb),
                ("NTPORS", &mut params.ntpors),
            ] {
                if let Some(v) = t.next_i32(name)? {
                    *slot = v;
                }
            }
            if let Some(v) = t.next_f64("CPOR")? {
                params.cpor = v;
            }
        },
        "DISTANCERES" => {
            let mut t = Tokens::new("DISTANCERES", data_lines);
            if let Some(v) = t.next_i32("NTDIR")? {
                params.ntdir = v;
            }
            if let Some(v) = t.next_i32("NTDIRA")? {
                params.ntdira = v;
            }
            for (name, slot) in [
                ("CDIR", &mut params.cdir),
                ("DIR0", &mut params.dir0),
                ("TAUDIR", &mut params.taudir),
            ] {
                if let Some(v) = t.next_f64(name)? {
                    *slot = v;
                }
            }
            for (name, slot) in [
                ("FORCESCALE", &mut params.forcescale),
                ("VDIR", &mut params.vdir),
                ("NTWDIR", &mut params.ntwdir),
            ] {
                if let Some(v) = t.next_i32(name)? {
                    *slot = v;
                }
            }
        },
        "PRECALCLAM" => {
            params.precalclam = Some(parse_precalclam(data_lines)?);
        },
        "EDS" => {
            params.eds = Some(parse_eds(data_lines)?);
        },
        "AEDS" => {
            params.aeds = Some(parse_aeds(data_lines)?);
        },
        "PERTURBATION" => {
            let mut t = Tokens::new("PERTURBATION", data_lines);
            if let Some(v) = t.next_i32("NTG")? {
                params.ntg = v;
            }
            if let Some(v) = t.next_i32("NRDGL")? {
                params.nrdgl = v;
            }
            for (name, slot) in [
                ("RLAM", &mut params.rlam),
                ("DLAMT", &mut params.dlamt),
                ("ALPHLJ", &mut params.alphlj),
                ("ALPHC", &mut params.alphc),
            ] {
                if let Some(v) = t.next_f64(name)? {
                    *slot = v;
                }
            }
            if let Some(v) = t.next_i32("NLAM")? {
                params.nlam = v;
            }
            if let Some(v) = t.next_i32("NSCALE")? {
                params.nscale = v;
            }
        },
        "COVALENTFORM" => {
            // gromosXX `read_COVALENTFORM`: NTBBH, NTBAH, NTBDN — which functional form the bond,
            // angle and dihedral terms take. The defaults (0 0 0) are the GROMOS ones.
            let mut t = Tokens::new("COVALENTFORM", data_lines);
            for (name, slot) in [
                ("NTBBH", &mut params.ntbbh),
                ("NTBAH", &mut params.ntbah),
                ("NTBDN", &mut params.ntbdn),
            ] {
                if let Some(v) = t.next_i32(name)? {
                    if v != 0 && v != 1 {
                        return Err(IoError::ParseError(format!(
                            "block COVALENTFORM: {name} must be 0 or 1, found {v}"
                        )));
                    }
                    *slot = v;
                }
            }
        },
        "ENERGYMIN" => {
            let mut t = Tokens::new("ENERGYMIN", data_lines);
            if let Some(v) = t.next_i32("NTEM")? {
                params.ntem = v;
            }
            if let Some(v) = t.next_usize("NCYC")? {
                params.ncyc = v;
            }
            for (name, slot) in [
                ("DELE", &mut params.dele),
                ("DX0", &mut params.dx0),
                ("DXM", &mut params.dxm),
            ] {
                if let Some(v) = t.next_f64(name)? {
                    *slot = v;
                }
            }
            if let Some(v) = t.next_usize("NMIN")? {
                params.nmin = v;
            }
            if let Some(v) = t.next_f64("FLIM")? {
                params.flim = v;
            }
        },
        _ => {
            // Unknown block - stored in raw_blocks by caller
        },
    }

    Ok(())
}

/// Maps the `constraints` convenience string used by `ImdParameters::nve/nvt/npt`
/// to the GROMOS NTC constraint mode.
///
/// - `"none"` → NTC=1 (no solute constraints — matches the pre-M1 default; a
///   factory-built run of a constrained system, e.g. one with solute H-bonds,
///   will diverge unless this is set to `"hbonds"`/`"allbonds"`)
/// - `"hbonds"` → NTC=2 (SHAKE-constrain solute bonds involving hydrogen)
/// - `"allbonds"` → NTC=3 (SHAKE-constrain all solute bonds)
///
/// This only sets NTC (solute constraints). Solvent rigidity (SETTLE, `NTCS`)
/// is orthogonal and unaffected — see P3.5.
pub fn ntc_from_constraints(constraints: &str) -> Result<i32, String> {
    match constraints {
        "none" => Ok(1),
        "hbonds" => Ok(2),
        "allbonds" => Ok(3),
        other => Err(format!(
            "unknown constraints mode {other:?}: expected \"none\", \"hbonds\", or \"allbonds\""
        )),
    }
}

impl ImdParameters {
    /// NVE ensemble: microcanonical, no thermostat or barostat.
    ///
    /// `constraints`: `"none"` | `"hbonds"` | `"allbonds"` — see [`ntc_from_constraints`].
    pub fn nve(dt: f64, steps: usize, constraints: &str) -> Result<Self, String> {
        let ntc = ntc_from_constraints(constraints)?;
        Ok(Self {
            title: String::from("NVE simulation"),
            nstlim: steps,
            dt,
            ntc,
            num_temp_baths: 1,
            temp_bath: vec![TempBathParameters {
                tau: vec![-1.0], // tau=-1 → no coupling (NVE)
                ..TempBathParameters::default()
            }],
            ..Self::default()
        })
    }

    /// NVT ensemble: canonical, Berendsen thermostat.
    ///
    /// `constraints`: `"none"` | `"hbonds"` | `"allbonds"` — see [`ntc_from_constraints`].
    pub fn nvt(dt: f64, steps: usize, temperature: f64, constraints: &str) -> Result<Self, String> {
        let ntc = ntc_from_constraints(constraints)?;
        Ok(Self {
            title: String::from("NVT simulation"),
            nstlim: steps,
            dt,
            ntc,
            num_temp_baths: 1,
            temp_bath: vec![TempBathParameters {
                algorithm: 0,
                temp0: vec![temperature],
                tau: vec![0.1],
                ..TempBathParameters::default()
            }],
            tempi: temperature,
            ..Self::default()
        })
    }

    /// NPT ensemble: isothermal-isobaric, Berendsen thermostat + barostat.
    ///
    /// `constraints`: `"none"` | `"hbonds"` | `"allbonds"` — see [`ntc_from_constraints`].
    pub fn npt(
        dt: f64,
        steps: usize,
        temperature: f64,
        pressure: f64,
        constraints: &str,
    ) -> Result<Self, String> {
        let ntc = ntc_from_constraints(constraints)?;
        let p_tensor = [[pressure; 3]; 3];
        let comp_tensor = [[4.575e-4; 3]; 3]; // water compressibility (nm² kJ⁻¹ mol)
        Ok(Self {
            title: String::from("NPT simulation"),
            nstlim: steps,
            dt,
            ntc,
            num_temp_baths: 1,
            temp_bath: vec![TempBathParameters {
                algorithm: 0,
                temp0: vec![temperature],
                tau: vec![0.1],
                ..TempBathParameters::default()
            }],
            tempi: temperature,
            couple_pressure: true,
            pressure_parameters: Some(PressureParameters {
                couple: 2,
                algorithm: 1,
                semi: [1, 1, 1],
                pressure0: p_tensor,
                compressibility: comp_tensor,
                tau_p: 0.5,
                virial: 1,
            }),
            ..Self::default()
        })
    }

    /// Steepest-descent energy minimization.
    pub fn steepest_descent(steps: usize) -> Self {
        Self {
            title: String::from("Energy minimization"),
            nstlim: steps,
            dt: 0.0,
            ntem: 1,
            nmin: 1,
            dele: 0.1,
            dx0: 0.01,
            dxm: 0.05,
            ntb: 1,
            ntc: 1,
            ..Self::default()
        }
    }
}

/// Split a line into whitespace-separated tokens
/// The constraint algorithm codes of a `CONSTRAINT` block, by gromosXX's names and numbers
/// (`in_parameter.cc`: shake 1, lincs 2, flexshake 3, settle 4, m_shake 5, gpu_shake 6), mapped
/// onto the codes this crate uses (1 shake, 2 lincs, 3 settle). The algorithms newen does not
/// implement are refused by name instead of being mapped onto a different one.
fn constraint_algorithm_code(what: &str, token: &str) -> Result<i32, IoError> {
    Ok(match token {
        "off" | "0" => 0,
        "shake" | "1" => 1,
        "lincs" | "2" => 2,
        "settle" | "4" => 3,
        "flexshake" | "3" => {
            return Err(IoError::ParseError(format!(
                "{what} = flexshake is not implemented by newen"
            )))
        },
        "m_shake" | "5" | "gpu_shake" | "6" => {
            return Err(IoError::ParseError(format!(
                "{what} = {token} is not implemented by newen"
            )))
        },
        other => {
            return Err(IoError::ParseError(format!(
                "{what} expects shake/lincs/settle (or 1/2/4), found {other:?}"
            )))
        },
    })
}

/// A block read as a stream of whitespace-separated values, the way gromosXX reads every input
/// block (`io::Block::get_next_parameter`): the parameters come in a fixed order and the file may
/// wrap them across lines however it likes. `next_*` returns `None` once the block runs out, so a
/// short block leaves the remaining parameters at their defaults.
struct Tokens<'a> {
    block: &'a str,
    values: Vec<String>,
    pos: usize,
}

impl<'a> Tokens<'a> {
    fn new(block: &'a str, lines: &[String]) -> Self {
        Tokens {
            block,
            values: lines.iter().flat_map(|l| parse_values(l)).collect(),
            pos: 0,
        }
    }

    fn next_raw(&mut self) -> Option<String> {
        let v = self.values.get(self.pos).cloned();
        if v.is_some() {
            self.pos += 1;
        }
        v
    }

    fn next_i32(&mut self, name: &str) -> Result<Option<i32>, IoError> {
        match self.next_raw() {
            None => Ok(None),
            Some(v) => v.parse::<i32>().map(Some).map_err(|_| {
                IoError::ParseError(format!(
                    "block {}: {name} expects an integer, found {v:?}",
                    self.block
                ))
            }),
        }
    }

    fn next_usize(&mut self, name: &str) -> Result<Option<usize>, IoError> {
        match self.next_raw() {
            None => Ok(None),
            Some(v) => v.parse::<usize>().map(Some).map_err(|_| {
                IoError::ParseError(format!(
                    "block {}: {name} expects a non-negative integer, found {v:?}",
                    self.block
                ))
            }),
        }
    }

    fn next_f64(&mut self, name: &str) -> Result<Option<f64>, IoError> {
        match self.next_raw() {
            None => Ok(None),
            Some(v) => v.parse::<f64>().map(Some).map_err(|_| {
                IoError::ParseError(format!(
                    "block {}: {name} expects a number, found {v:?}",
                    self.block
                ))
            }),
        }
    }
}

fn parse_values(line: &str) -> Vec<String> {
    line.split_whitespace().map(|s| s.to_string()).collect()
}

fn parse_f64(s: &str) -> Result<f64, IoError> {
    s.parse::<f64>()
        .map_err(|_| IoError::ParseError(format!("expected a number, found {s:?}")))
}

fn parse_i32(s: &str) -> Result<i32, IoError> {
    s.parse::<i32>()
        .map_err(|_| IoError::ParseError(format!("expected an integer, found {s:?}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(content: &str, suffix: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("newen_imd_test_{suffix}.tmp"));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_default_parameters() {
        let params = ImdParameters::default();
        assert_eq!(params.nstlim, 1000);
        assert_eq!(params.dt, 0.002);
        assert_eq!(params.ntc, 1); // No constraints
    }

    #[test]
    fn test_parse_gromosxx_format() {
        // Actual GROMOS .in file format (positional)
        let content = "\
TITLE
Two Argon atoms - LJ pair reference
END
SYSTEM
#      NPM      NSM
         1        0
END
STEP
#   NSTLIM         T        DT
        10       0.0     0.002
END
BOUNDCOND
#      NTB    NDFMIN
         0         0
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
     0        0         0        0     0  1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     2     1      2
END
PAIRLIST
#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE
        standard        1       0.8     1.4     auto    chargegroup
END
NONBONDED
#   NLRELE
         0
#    APPAK      RCRF     EPSRF  NSLFEXCL
       0.0       1.4       1.0         1
END
CONSTRAINT
#       NTC
        1
#       NTCP
        shake
#       NTCP0(1)
        0.0001
#       NTCS
        shake
#       NTCS0(1)
        0.0001
END
WRITETRAJ
#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
         1         0         0         1         1         0         0
END
PRINTOUT
#     NTPR      NTPP
         1         0
END
MULTIBATH
    weak-coupling
#   NBATHS
    1
#   TEMP0  TAU
    300.0  -1.0
#   DOFSET
    1
#   LAST   COM-BATH  IR-BATH
    2      1         1
END
INITIALISE
#  NTIVEL  NTISHK  NTINHT  NTINHB
        0       0       0       0
#  NTISHI  NTIRTC  NTICOM
        1       0       0
#  NTISTI
        0
#      IG   TEMPI
   210185     0.0
END
";
        let path = write_tmp(content, "gxx_format");
        let params = read_imd_file(&path).expect("Failed to parse GROMOS format");

        assert_eq!(params.npm, 1);
        assert_eq!(params.nsm, 0);
        assert_eq!(params.nstlim, 10);
        assert_eq!(params.t0, 0.0);
        assert_eq!(params.dt, 0.002);
        assert_eq!(params.ntb, 0);
        assert_eq!(params.ndfmin, 0);
        assert_eq!(params.ntf, [0, 0, 0, 0, 0, 1]);
        assert_eq!(params.negr, 2);
        assert_eq!(params.nre, vec![1, 2]);
        assert_eq!(params.nlrele, 0);
        assert_eq!(params.appak, 0.0);
        assert_eq!(params.rcrf, 1.4);
        assert_eq!(params.epsrf, 1.0);
        assert_eq!(params.nslfexcl, 1);
        assert_eq!(params.nsnb, 1);
        assert_eq!(params.rcutp, 0.8);
        assert_eq!(params.rcutl, 1.4);
        assert_eq!(params.ntc, 1);
        assert_eq!(params.shake_tol, 0.0001);
        assert_eq!(params.ntwx, 1);
        assert_eq!(params.ntwe, 1);
        assert_ne!(params.ntwf, 0);
        assert_eq!(params.ntpr, 1);
        assert_eq!(params.temp_bath[0].temp0, vec![300.0]);
        assert_eq!(params.temp_bath[0].tau, vec![-1.0]);
        assert_eq!(params.ntivel, 0);
        assert_eq!(params.ig, 210185);
        assert_eq!(params.tempi, 0.0);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_parse_water_box_format() {
        let content = "\
TITLE
216 SPC waters in rectangular box
END
SYSTEM
#      NPM      NSM
         1        0
END
STEP
#   NSTLIM         T        DT
       100       0.0     0.002
END
BOUNDCOND
#      NTB    NDFMIN
         1         0
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
     1        1         1        1     1  1
# NEGR    NRE(1)
     1     648
END
NONBONDED
#   NLRELE
         1
#    APPAK      RCRF     EPSRF  NSLFEXCL
       0.0       0.9      62.0         1
END
PAIRLIST
#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE
        standard        5       0.8     0.9     auto    chargegroup
END
";
        let path = write_tmp(content, "water_box");
        let params = read_imd_file(&path).expect("Failed to parse water box");

        assert_eq!(params.nstlim, 100);
        assert_eq!(params.ntb, 1);
        assert_eq!(params.ntf, [1, 1, 1, 1, 1, 1]);
        assert_eq!(params.negr, 1);
        assert_eq!(params.nre, vec![648]);
        assert_eq!(params.nlrele, 1);
        assert_eq!(params.rcrf, 0.9);
        assert_eq!(params.epsrf, 62.0);
        assert_eq!(params.nsnb, 5);
        assert_eq!(params.rcutl, 0.9);

        std::fs::remove_file(path).ok();
    }

    /// PERTURBATION block: multi-line format (GROMOS splits 8 values across 2 non-comment
    /// lines). This was a real bug: `data_lines.first()` silently dropped ALPHLJ/ALPHC/NLAM/NSCALE.
    #[test]
    fn test_perturbation_block_multiline() {
        let content = "\
TITLE
CH4 FEP test
END
SYSTEM
#      NPM      NSM
         1      999
END
PERTURBATION
#     NTG   NRDGL    RLAM   DLAMT
        1       0   0.500     0.0
#  ALPHLJ   ALPHC    NLAM  NSCALE
      0.5     0.5       1       0
END
";
        let path = write_tmp(content, "perturbation_multiline");
        let p = read_imd_file(&path).expect("should parse multi-line PERTURBATION");

        assert_eq!(p.ntg, 1, "NTG");
        assert_eq!(p.nrdgl, 0, "NRDGL");
        assert!((p.rlam - 0.5).abs() < 1e-10, "RLAM");
        assert!((p.dlamt - 0.0).abs() < 1e-10, "DLAMT");
        assert!(
            (p.alphlj - 0.5).abs() < 1e-10,
            "ALPHLJ — was silently dropped before fix"
        );
        assert!(
            (p.alphc - 0.5).abs() < 1e-10,
            "ALPHC  — was silently dropped before fix"
        );
        assert_eq!(p.nlam, 1, "NLAM");
        assert_eq!(p.nscale, 0, "NSCALE");

        std::fs::remove_file(path).ok();
    }

    /// PERTURBATION block: single-line format (all 8 values on one line — also valid).
    #[test]
    fn test_perturbation_block_singleline() {
        let content = "\
TITLE
test
END
SYSTEM
         1        0
END
PERTURBATION
        1    0    0.125    0.0    1.5    1.5    2    0
END
";
        let path = write_tmp(content, "perturbation_singleline");
        let p = read_imd_file(&path).expect("should parse single-line PERTURBATION");
        assert_eq!(p.ntg, 1);
        assert!((p.rlam - 0.125).abs() < 1e-10, "RLAM");
        assert!((p.alphlj - 1.5).abs() < 1e-10, "ALPHLJ");
        assert!((p.alphc - 1.5).abs() < 1e-10, "ALPHC");
        assert_eq!(p.nlam, 2, "NLAM");
        std::fs::remove_file(path).ok();
    }

    /// Verify MULTIBATH algorithm mapping matches GROMOS exactly:
    ///   "weak-coupling" → 0, "nose-hoover" → 1, "nose-hoover-chains N" → N (N ≥ 2)
    #[test]
    fn test_multibath_algorithm_mapping() {
        fn make_imd(algo_line: &str) -> String {
            format!("\
TITLE\ntest\nEND\nSYSTEM\n#      NPM      NSM\n         1        0\nEND\n\
MULTIBATH\n    {algo_line}\n#   NBATHS\n    1\n#   TEMP0  TAU\n    300.0  0.1\n#   DOFSET\n    1\n#   LAST   COM-BATH  IR-BATH\n    648    1         1\nEND\n")
        }

        // 0 = Berendsen / weak-coupling
        let p = read_imd_file(write_tmp(&make_imd("weak-coupling"), "nhc_wc")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 0, "weak-coupling should map to 0");

        // 1 = NHC single
        let p = read_imd_file(write_tmp(&make_imd("nose-hoover"), "nhc_nh")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 1, "nose-hoover should map to 1");

        // N = NHC chain of length N (N >= 2)
        let p = read_imd_file(write_tmp(&make_imd("nose-hoover-chains  3"), "nhc_chain3")).unwrap();
        assert_eq!(
            p.temp_bath[0].algorithm, 3,
            "nose-hoover-chains 3 should map to 3"
        );
        assert_eq!(p.temp_bath[0].nhc_chain, 3, "nhc_chain should be 3");

        // Numeric fallback: "0" → 0, "1" → 1
        let p = read_imd_file(write_tmp(&make_imd("0"), "nhc_num0")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 0, "numeric 0 should map to 0");

        let p = read_imd_file(write_tmp(&make_imd("1"), "nhc_num1")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 1, "numeric 1 should map to 1");
    }

    #[test]
    fn test_factory_nve() {
        let p = ImdParameters::nve(0.001, 500, "none").unwrap();
        assert_eq!(p.nstlim, 500);
        assert!((p.dt - 0.001).abs() < 1e-15);
        assert!(
            (p.temp_bath[0].tau[0] - (-1.0)).abs() < 1e-15,
            "NVE needs tau=-1"
        );
        assert_eq!(p.ntc, 1, "\"none\" -> NTC=1");
    }

    #[test]
    fn test_factory_nvt() {
        let p = ImdParameters::nvt(0.002, 1000, 298.0, "none").unwrap();
        assert_eq!(p.nstlim, 1000);
        assert!((p.dt - 0.002).abs() < 1e-15);
        assert!((p.temp_bath[0].temp0[0] - 298.0).abs() < 1e-10);
        assert_eq!(p.temp_bath[0].algorithm, 0, "NVT uses Berendsen");
        assert!(!p.couple_pressure);
        assert_eq!(p.ntc, 1);
    }

    #[test]
    fn test_factory_npt() {
        let p = ImdParameters::npt(0.002, 1000, 300.0, 1.0, "none").unwrap();
        assert_eq!(p.nstlim, 1000);
        assert!((p.temp_bath[0].temp0[0] - 300.0).abs() < 1e-10);
        assert!(p.couple_pressure, "NPT must couple pressure");
        let pp = p.pressure_parameters.unwrap();
        assert!((pp.pressure0[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_factory_steepest_descent() {
        let p = ImdParameters::steepest_descent(200);
        assert_eq!(p.nstlim, 200);
        assert_eq!(p.ntem, 1);
    }

    #[test]
    fn test_factory_constraints_knob() {
        assert_eq!(ImdParameters::nve(0.001, 10, "hbonds").unwrap().ntc, 2);
        assert_eq!(ImdParameters::nve(0.001, 10, "allbonds").unwrap().ntc, 3);
        assert_eq!(
            ImdParameters::nvt(0.001, 10, 300.0, "hbonds").unwrap().ntc,
            2
        );
        assert_eq!(
            ImdParameters::npt(0.001, 10, 300.0, 1.0, "allbonds")
                .unwrap()
                .ntc,
            3
        );
        assert!(ImdParameters::nve(0.001, 10, "bogus").is_err());
    }

    #[test]
    fn energy_groups_out_of_order_are_an_error() {
        // gromosXX in_parameter.cc:986-993: NRE strictly increasing
        let content = "TITLE\nx\nEND\nFORCE\n 1 1 1 1 1 1\n 2 5 5\nEND\n";
        let path = write_tmp(content, "nre_order");
        let err = read_imd_file(&path).unwrap_err().to_string();
        assert!(err.contains("strictly increasing"), "{err}");
    }

    // ── EDS (gromosXX read_EDS, in_parameter.cc:3843-3961) ───────────────────────────────────

    fn eds_block(body: &str) -> Result<Option<EdsBlock>, IoError> {
        parse_imd_str(&format!("EDS\n{body}\nEND\n")).map(|p| p.eds)
    }

    #[test]
    fn eds_parses_the_three_forms_and_round_trips() {
        let one = eds_block("# EDS\n 1\n# FORM NUMSTATES\n 1 2\n# S\n 0.5\n# EIR\n 0.0 10.0")
            .unwrap()
            .unwrap();
        assert!(one.enabled);
        assert_eq!((one.form, one.num_states), (1, 2));
        assert_eq!(
            (one.s.clone(), one.eir.clone()),
            (vec![0.5], vec![0.0, 10.0])
        );
        let multi = eds_block("1\n2 3\n0.3 0.2 0.1\n0 10 -5").unwrap().unwrap();
        assert_eq!(multi.s, [0.3, 0.2, 0.1]);
        // FORM 3: n−1 rows `i j S`, a token stream (one row may span lines)
        let pairs = eds_block("1\n3 3\n1 2 0.3\n2\n3 0.1\n0 10 -5")
            .unwrap()
            .unwrap();
        assert_eq!(
            pairs.pairs,
            [
                EdsStatePair { i: 1, j: 2, s: 0.3 },
                EdsStatePair { i: 2, j: 3, s: 0.1 }
            ]
        );
        assert!(pairs.s.is_empty());
        let off = eds_block("0\n1 2").unwrap().unwrap();
        assert!(!off.enabled && off.eir.is_empty());
        for b in [one, multi, pairs, off] {
            let p = ImdParameters {
                eds: Some(b.clone()),
                ..Default::default()
            };
            let again = parse_imd_str(&crate::imd_write::write_imd(&p, None)).unwrap();
            assert_eq!(again.eds, Some(b));
        }
        assert_eq!(parse_imd_str("TITLE\nx\nEND\n").unwrap().eds, None);
    }

    #[test]
    fn eds_refuses_what_gromosxx_refuses() {
        for (body, what) in [
            // EDS = 0 stops after NUMSTATES: anything more is left over
            ("0\n1 2\n0.5\n0 0", "Left-over"),
            ("1\n1 2\n0.5\n0 0 7", "Left-over"),
            ("2\n1 2\n0.5\n0 0", "EDS = 2"),
            ("1\n4 2\n0.5\n0 0", "FORM = 4"),
            ("0\n4 2", "FORM = 4"),
            ("0\n1 1", "NUMSTATES = 1"),
            ("1\n1 2\n0.0\n0 0", "S[0] = 0"),
            // FORM 2 with three states needs three S: the stream takes the first EIR as S[2]
            ("1\n2 3\n0.3 0.2\n0 10 -5", "S[2] = 0"),
            ("1\n2 3\n0.3 0.2 0.1\n0 10", "missing EIR[2]"),
            ("1\n3 3\n1 2 0.3\n0 10 -5", "S[1]"),
            ("1\n1 2\nx\n0 0", "S[0] expects a number"),
        ] {
            let e = eds_block(body).unwrap_err().to_string();
            assert!(e.contains(what) && e.contains("EDS"), "{body:?}: {e}");
        }
    }

    // ── AEDS (gromosXX read_AEDS, in_parameter.cc:3968-4140) ─────────────────────────────────
    fn aeds_block(body: &str) -> Result<Option<AedsBlock>, IoError> {
        parse_imd_str(&format!("AEDS\n{body}\nEND\n")).map(|p| p.aeds)
    }

    #[test]
    fn aeds_parses_the_token_stream() {
        // gromosXX's own example block
        let b = aeds_block(
            "# AEDS\n  1\n# FORM  NUMSTATES\n  4  5\n# EMAX  EMIN\n  10  -50\n# EIR\n  0 -5 -140 \
             -560 -74\n# NTIAEDSS  RESTREMIN  BMAXTYPE  BMAX  ASTEPS  BSTEPS\n  1 1 2 3 500 50000",
        )
        .unwrap()
        .unwrap();
        assert!(b.enabled && b.ntiaedss && b.restremin);
        assert_eq!((b.form, b.num_states, b.bmaxtype), (4, 5, 2));
        assert_eq!((b.emax, b.emin, b.bmax), (10.0, -50.0, 3.0));
        assert_eq!(b.eir, vec![0.0, -5.0, -140.0, -560.0, -74.0]);
        assert_eq!((b.asteps, b.bsteps), (500, 50000));
        // EMIN > EMAX parses (gromosXX warns); the recipe refuses it
        assert_eq!(
            aeds_block("1\n1 2\n0 10\n0 0\n1 1 1 2 10 10")
                .unwrap()
                .unwrap()
                .emin,
            10.0
        );
        let off = aeds_block("0\n6 2").unwrap().unwrap();
        assert!(!off.enabled && off.eir.is_empty());
    }

    #[test]
    fn aeds_refuses_what_gromosxx_refuses() {
        for (body, what) in [
            ("0\n1 2\n0 0", "Left-over"),
            ("1\n1 2\n0 0\n0 0\n1 1 1 2 10 10 3", "Left-over"),
            ("2\n1 2", "AEDS = 2"),
            ("1\n7 2", "FORM = 7"),
            ("0\n1 1", "NUMSTATES = 1"),
            ("1\n1 2\n0 0\n0 0\n2 1 1 2 10 10", "NTIAEDSS = 2"),
            ("1\n1 2\n0 0\n0 0\n1 2 1 2 10 10", "RESTREMIN = 2"),
            ("1\n1 2\n0 0\n0 0\n1 1 3 2 10 10", "BMAXTYPE = 3"),
            ("1\n1 2\n0 0\n0 0\n1 1 1 0 10 10", "BMAX = 0"),
            ("1\n1 2\n0 0\n0 0\n1 1 1 2 0 10", "ASTEPS = 0"),
            ("1\n1 2\n0 0\n0 0\n1 1 1 2 10 -1", "BSTEPS = -1"),
            ("1\n1 2\n0 0\n0\n1 1 1 2 10", "missing ASTEPS"),
        ] {
            let e = aeds_block(body).unwrap_err().to_string();
            assert!(e.contains(what) && e.contains("AEDS"), "{body:?}: {e}");
        }
    }

    // ── PRECALCLAM (gromosXX read_PRECALCLAM, in_parameter.cc:4478-4513) ───────────────────

    #[test]
    fn precalclam_parses_checks_its_ranges_and_round_trips() {
        let p = parse_imd_str("PRECALCLAM\n# NRLAM MINLAM MAXLAM\n 11 0.2 0.7\nEND\n").unwrap();
        let b = PrecalcLamBlock {
            nr_lambdas: 11,
            min_lambda: 0.2,
            max_lambda: 0.7,
        };
        assert_eq!(p.precalclam, Some(b));
        let again = parse_imd_str(&crate::imd_write::write_imd(&p, None)).unwrap();
        assert_eq!(again.precalclam, Some(b));
        assert_eq!(parse_imd_str("TITLE\nx\nEND\n").unwrap().precalclam, None);
        for (body, what) in [
            (" -1 0.0 1.0", "NRLAM"),
            (" 11 0.0 1.5", "MAXLAM"),
            (" 11 -0.1 1.0", "MINLAM"),
            (" 11 0.5 0.5", "MINLAM should be smaller than MAXLAM"),
        ] {
            let err = parse_imd_str(&format!("PRECALCLAM\n{body}\nEND\n"))
                .unwrap_err()
                .to_string();
            assert!(
                err.contains(what) && err.contains("PRECALCLAM"),
                "{body}: {err}"
            );
        }
    }

    // ── LAMBDAS (gromosXX read_LAMBDAS, in_parameter.cc:4326-4472) ──────────────────────────

    /// A perturbed input with two energy groups and the given LAMBDAS body; `ntg` sets NTG.
    fn lambdas_input(ntg: i32, body: &str) -> String {
        format!(
            "FORCE\n 1 1 1 1 1 1\n 2 6 12\nEND\nPERTURBATION\n {ntg} 0 0.125 0.0\n 0.5 0.5 1 0\nEND\n\
             LAMBDAS\n{body}END\n"
        )
    }

    fn parse_lambdas_text(ntg: i32, body: &str) -> Result<ImdParameters, IoError> {
        parse_imd_str(&lambdas_input(ntg, body))
    }

    fn lambdas_err(ntg: i32, body: &str) -> String {
        parse_lambdas_text(ntg, body).unwrap_err().to_string()
    }

    #[test]
    fn precalclam_reads_and_checks_its_three_values() {
        let p = parse_imd_str("PRECALCLAM\n# NRLAM MINLAM MAXLAM\n  21 0.0 1.0\nEND\n").unwrap();
        let pre = p.precalclam().unwrap().unwrap();
        assert_eq!((pre.nrlam, pre.minlam, pre.maxlam), (21, 0.0, 1.0));
        assert!(parse_imd_str("TITLE\nx\nEND\n")
            .unwrap()
            .precalclam()
            .unwrap()
            .is_none());
        // the range checks already fire in parse_imd_str; a missing MAXLAM only in precalclam()
        for bad in ["11 0.7 0.2", "11 0 1.5", "-1 0 1", "11 0"] {
            let r =
                parse_imd_str(&format!("PRECALCLAM\n{bad}\nEND\n")).and_then(|p| p.precalclam());
            assert!(r.is_err(), "{bad}");
        }
    }

    #[test]
    fn lambdas_ntil_takes_on_off_one_zero() {
        for (ntil, on) in [("on", true), ("1", true), ("off", false), ("0", false)] {
            let p = parse_lambdas_text(1, &format!(" {ntil}\n")).unwrap();
            assert_eq!(p.lambdas.as_ref().unwrap().enabled, on, "{ntil}");
        }
        let err = lambdas_err(1, " yes\n");
        assert!(err.contains("NTIL"), "{err}");
    }

    #[test]
    fn lambdas_off_ignores_the_rows() {
        // gromosXX returns right after NTIL: a row it would refuse is not even read
        let p = parse_lambdas_text(1, " off\n angleres 1 1 0 0 1 0 0\n").unwrap();
        let l = p.lambdas.unwrap();
        assert!(!l.enabled && !l.applies(1));
        assert!(l.warnings(1).is_empty());
        // rows that parse are kept for the round trip, inert
        let p = parse_lambdas_text(1, " off\n vdw 1 2 0 0 1 0 0\n").unwrap();
        let l = p.lambdas.unwrap();
        assert_eq!(l.rows.len(), 1);
        assert!(!l.applies(1));
    }

    #[test]
    fn lambdas_on_without_perturbation_is_ignored_with_a_warning() {
        let p = parse_lambdas_text(0, " on\n angleres 1 1 0 0 1 0 0\n").unwrap();
        let l = p.lambdas.unwrap();
        assert!(l.enabled && !l.applies(0));
        assert_eq!(
            l.warnings(0),
            vec!["LAMBDAS block without perturbation is ignored".to_string()]
        );
    }

    #[test]
    fn lambdas_ntli_is_a_name_or_its_number() {
        let p = parse_lambdas_text(
            1,
            " on\n vdw 1 2 0 0 0.1 0.9 0\n 5 1 1 0 1 0 0 0\n 12 2 2 0.2 0 0 0.8 0\n",
        )
        .unwrap();
        let rows = p.lambdas.unwrap().rows;
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].kind, InteractionLambda::Lj);
        assert_eq!((rows[0].g1, rows[0].g2), (1, 2));
        assert_eq!(rows[0].coefficients, [0.0, 0.0, 0.1, 0.9, 0.0]);
        assert_eq!(rows[1].kind, InteractionLambda::Lj);
        assert_eq!(rows[2].kind, InteractionLambda::Mass);
    }

    #[test]
    fn lambdas_refuses_angleres_and_unknown_kinds() {
        for ntli in ["angleres", "13", "0", "colvarres", "lj"] {
            let err = lambdas_err(1, &format!(" on\n {ntli} 1 1 0 0 1 0 0\n"));
            assert!(err.contains("NTLI"), "{ntli}: {err}");
        }
    }

    #[test]
    fn lambdas_groups_must_be_positive() {
        let err = lambdas_err(1, " on\n bond 0 1 0 0 1 0 0\n");
        assert!(err.contains("NILG1[0] must be > 0"), "{err}");
    }

    #[test]
    fn lambdas_groups_must_exist_before_the_implicit_one() {
        // NEGR = 2: group 3 would be gromosXX's implicit trailing group, which LAMBDAS cannot name
        let err = lambdas_err(1, " on\n vdw 1 3 0 0 1 0 0\n");
        assert!(
            err.contains("smaller than the number of energy groups"),
            "{err}"
        );
    }

    #[test]
    fn lambdas_second_group_must_not_precede_the_first() {
        let err = lambdas_err(1, " on\n vdw 2 1 0 0 1 0 0\n");
        assert!(err.contains("NILG2 >= NILG1"), "{err}");
    }

    #[test]
    fn lambdas_off_diagonal_bonded_row_is_kept_with_a_warning() {
        let p = parse_lambdas_text(1, " on\n bond 1 2 0 0 1 0 0\n crf 1 2 0 0 1 0 0\n").unwrap();
        let l = p.lambdas.unwrap();
        assert_eq!(l.rows.len(), 2);
        let w = l.warnings(1);
        assert_eq!(w.len(), 1, "{w:?}");
        assert!(w[0].contains("only allowed for nonbonded"), "{w:?}");
    }

    #[test]
    fn lambdas_rows_are_checked_against_blocks_that_come_later() {
        // gromosXX reads LAMBDAS after FORCE and PERTURBATION whatever the file order
        let text = "LAMBDAS\n on\n vdw 1 3 0 0 1 0 0\nEND\nFORCE\n 1 1 1 1 1 1\n 2 6 12\nEND\n\
                    PERTURBATION\n 1 0 0.125 0.0\n 0.5 0.5 1 0\nEND\n";
        let err = parse_imd_str(text).unwrap_err().to_string();
        assert!(
            err.contains("smaller than the number of energy groups"),
            "{err}"
        );
    }

    #[test]
    fn lambdas_missing_value_is_an_error() {
        let err = lambdas_err(1, " on\n vdw 1 2 0 0 1 0\n");
        assert!(err.contains("missing ELI[0]"), "{err}");
    }

    #[test]
    fn lambdas_left_over_values_are_an_error() {
        // two rows on one line: gromosXX reads one row per line and reports the rest
        let err = lambdas_err(1, " on\n vdw 1 2 0 0 1 0 0  crf 1 2 0 0 1 0 0\n");
        assert!(err.contains("Left-over parameters"), "{err}");
        let err = lambdas_err(1, " on\n vdw 1 2 0 0 1 0 0 7\n");
        assert!(err.contains("Left-over parameters"), "{err}");
        // a trailing comment is not a value
        let p = parse_lambdas_text(1, " on # NTIL\n vdw 1 2 0 0 1 0 0 # the LJ row\n").unwrap();
        assert_eq!(p.lambdas.unwrap().rows.len(), 1);
        // a block that does not apply is never read to the end
        assert!(parse_lambdas_text(0, " on\n vdw 1 2 0 0 1 0 0 7\n").is_ok());
    }

    #[test]
    fn lambdas_round_trip_through_the_writer() {
        let p = parse_lambdas_text(
            1,
            " on\n bond 1 1 0 0 1 0 0\n vdw 1 2 0 0 0.1 0.9 0\n crf_soft 2 2 1 1 -1 0 0\n",
        )
        .unwrap();
        let text = crate::imd_write::write_imd(&p, None);
        let back = parse_imd_str(&text).unwrap();
        assert_eq!(back.lambdas, p.lambdas, "{text}");
    }
}
