//! Poisson–Boltzmann electrostatics: the gromos++ `pb` library (`src/pb/`), ported for
//! `dGslv_pbsolv`.
//!
//! A solve sees a [`PbSystem`]: the solute positions and charges, the PB radii, the atoms that are
//! charged on the grid and carve out the low-dielectric cavity (`@atoms`) and the atoms whose
//! potentials and energy are read out (`@atomsTOcharge`). [`solve_fd`] is one finite-difference
//! solve (gromos++ `FDPoissonBoltzmann` with the ICCG solvers); [`dgslv_fd`] is the sequence of
//! four solves `dGslv_pbsolv` runs for every scheme. [`solve_fft`] is one periodic FFT solve of
//! the reaction potential (gromos++ `FFTPoisson`), [`dgslv_fft`] the LS and RF pair
//! `@schemeELEC RF` adds.
//!
//! Every routine reproduces gromos++'s floating-point expression order: the iteration count of the
//! L1-converged ICCG solver, and with it the printed potentials, depend on it (item 1.8). The FFT
//! solves differ from gromos++'s FFTW in the last bits.
//!
//! The log (`@verbose`) is gromos++'s, line for line, in C++ stream format.

pub mod fd;
pub mod fft;
pub mod fft_fields;
pub mod fft_poisson;
pub mod iccg;
pub mod params;
pub mod radii;

pub use fd::{solve_fd, FdBoundary, FdSolution};
pub use fft::FftGrid;
pub use fft_poisson::{
    dgslv_fft, solve_fft, DgslvFft, FftScheme, FftSetup, FftSolution, LogFormat, StreamLog,
};
pub use params::{PbParameters, GROMOSPP_FOUR_PI_EPS_I};
pub use radii::{atomic_radii_vdw, pb_radius, RadiusDefinition};

use newen_core::Vec3;

/// The atoms a PB solve sees, all in nm and e.
///
/// Every solute atom is listed: the grid is centred on the bounding box of all of them
/// (gromos++ `PositionUtils::getmin/maxcoordinates` over the solute molecules), whichever subset is
/// charged.
#[derive(Debug, Clone)]
pub struct PbSystem {
    /// Position of every solute atom.
    pub positions: Vec<Vec3>,
    /// Charge of every solute atom.
    pub charges: Vec<f64>,
    /// PB radius of every solute atom; only those of [`PbSystem::atoms`] are read.
    pub radii: Vec<f64>,
    /// The atoms charged on the grid and treated as solute dielectric (`@atoms`), in order.
    pub atoms: Vec<usize>,
    /// The atoms whose potential and energy are read out (`@atomsTOcharge`), in order.
    pub readout: Vec<usize>,
    /// 1/(4πε₀) in kJ mol⁻¹ nm e⁻²: the topology's PHYSICALCONSTANTS value, or
    /// [`GROMOSPP_FOUR_PI_EPS_I`] without a topology.
    pub fpepsi: f64,
}

impl PbSystem {
    /// The lowest and highest coordinate over all solute atoms, per axis, starting from ±1e8 as
    /// gromos++ does (`PositionUtils.cc:252-290`).
    pub fn extent(&self) -> (Vec3, Vec3) {
        extent(&self.positions)
    }
}

/// The lowest and highest coordinate per axis, starting from ±1e8 (`PositionUtils.cc:252-290`).
pub fn extent(positions: &[Vec3]) -> (Vec3, Vec3) {
    let mut lo = [100000000.0f64; 3];
    let mut hi = [-100000000.0f64; 3];
    for p in positions {
        for d in 0..3 {
            if p[d] < lo[d] {
                lo[d] = p[d];
            }
            if p[d] > hi[d] {
                hi[d] = p[d];
            }
        }
    }
    (Vec3::from_array(lo), Vec3::from_array(hi))
}

/// The four finite-difference solves of `dGslv_pbsolv` and the two free energies formed from them
/// (`dGslv_pbsolv.cc:740-752`).
#[derive(Debug, Clone)]
pub struct DgslvFd {
    pub npbc_vac: FdSolution,
    pub npbc_slv: FdSolution,
    pub pbc_vac: FdSolution,
    pub pbc_slv: FdSolution,
    /// `npbc_slv.dg − npbc_vac.dg`: the electrostatic solvation free energy with macroscopic
    /// (Coulombic, non-periodic) boundary conditions, kJ/mol.
    pub dg_npbc: f64,
    /// `pbc_slv.dg − pbc_vac.dg`: the same under periodic (lattice-sum) boundary conditions.
    pub dg_pbc: f64,
}

/// Grids and permittivities of the finite-difference solves.
#[derive(Debug, Clone, Copy)]
pub struct FdSetup {
    /// Grid points of the non-periodic grid.
    pub npbc_grid: [usize; 3],
    /// Grid points of the periodic grid.
    pub pbc_grid: [usize; 3],
    /// Grid spacing, nm.
    pub h: f64,
    /// Relative permittivity of the solvent in the periodic solve (`@epsSOLV`).
    pub eps_solv: f64,
    /// Relative permittivity of the solvent in the non-periodic solve (`@epsNPBC`).
    pub eps_npbc: f64,
    /// Iteration limit of the ICCG solver.
    pub maxiter: usize,
}

/// Which order the program runs the two solves of a pair in: the coordinate mode solves the
/// vacuum first, the PQR mode the solvent first (`dGslv_pbsolv.cc:740-758` vs `:1030-1047`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveOrder {
    VacuumFirst,
    SolventFirst,
}

/// The NPBC pair, then the PBC pair, each with its banner and `# DGRESULT` line in the log.
pub fn dgslv_fd(
    sys: &PbSystem,
    setup: &FdSetup,
    order: SolveOrder,
    log: &mut String,
) -> Result<DgslvFd, String> {
    let npbc = |eps: f64, banner: &str, log: &mut String| {
        banner_lines(log, banner);
        solve_fd(
            sys,
            setup.npbc_grid,
            setup.h,
            FdBoundary::Npbc,
            eps,
            setup.maxiter,
            log,
        )
    };
    let pbc = |eps: f64, banner: &str, log: &mut String| {
        banner_lines(log, banner);
        solve_fd(
            sys,
            setup.pbc_grid,
            setup.h,
            FdBoundary::Pbc,
            eps,
            setup.maxiter,
            log,
        )
    };
    // fd_ls_pbc_slv closes its section with one more star line (dGslv_pbsolv.cc:210).
    let pbc_slv = |log: &mut String| {
        let s = pbc(setup.eps_solv, "FD LS PBC EPSSOLV", log)?;
        log.push_str(STARS);
        Ok::<_, String>(s)
    };

    let (npbc_vac, npbc_slv) = match order {
        SolveOrder::VacuumFirst => {
            let v = npbc(1.0, "FD LS NPBC VAC", log)?;
            let s = npbc(setup.eps_npbc, "FD LS NPBC EPSSOLV", log)?;
            (v, s)
        },
        SolveOrder::SolventFirst => {
            let s = npbc(setup.eps_npbc, "FD LS NPBC EPSSOLV", log)?;
            let v = npbc(1.0, "FD LS NPBC VAC", log)?;
            (v, s)
        },
    };
    let dg_npbc = npbc_slv.dg - npbc_vac.dg;
    log.push_str(&format!("# DGRESULT NPBC {}\n", g(dg_npbc)));

    let (pbc_vac, pbc_slv) = match order {
        SolveOrder::VacuumFirst => {
            let v = pbc(1.0, "FD LS PBC VAC", log)?;
            let s = pbc_slv(log)?;
            (v, s)
        },
        SolveOrder::SolventFirst => {
            let s = pbc_slv(log)?;
            let v = pbc(1.0, "FD LS PBC VAC", log)?;
            (v, s)
        },
    };
    let dg_pbc = pbc_slv.dg - pbc_vac.dg;
    log.push_str(&format!("# DGRESULT PBC {}\n", g(dg_pbc)));

    Ok(DgslvFd {
        npbc_vac,
        npbc_slv,
        pbc_vac,
        pbc_slv,
        dg_npbc,
        dg_pbc,
    })
}

const STARS: &str = "# ************************************************** \n";

fn banner_lines(log: &mut String, banner: &str) {
    log.push_str(STARS);
    log.push_str(&format!("# *** {banner} *** \n"));
    log.push_str(STARS);
}

/// A double the way a C++ stream prints it by default (`%g`, 6 significant digits, `-0` kept).
pub(crate) fn g(x: f64) -> String {
    if x == 0.0 && x.is_sign_negative() {
        return "-0".to_string();
    }
    newen_io::table::cpp_g(x, 6)
}
