//! gromos++'s periodic Poisson solver, used by `dGslv_pbsolv` for the lattice-sum (LS) and
//! reaction-field (RF) reaction potentials: `pb::FFTPoisson` and `pb::FFTPoissonIterator` with
//! their interaction tensors (`FFTDipoleDipole_{LS,RF}`, `FFTChargeDipole_{LS,RF}`).
//!
//! The solute is a cavity of permittivity 1 (the occupancy of [`inside_sol`]) in a continuum of
//! permittivity εs. Starting from the vacuum field of the charges (zero inside the cavity), the
//! iteration applies the k-space dipole–dipole tensor of the boundary condition, returns to real
//! space, and corrects the field inside the cavity with an adaptive over-relaxation until the RMS
//! field inside drops below `convergence_fft`. The reaction potential `φ̂ = −i·P(k)(εs−1)/k²·k·Ê`
//! is then interpolated at the read-out atoms; `dG = ½ Σ qφ`. The potentials are reaction-field
//! potentials only, in kJ mol⁻¹ e⁻¹.
//!
//! The log is gromos++'s, line for line, including the eighteen `PB_Parameters` blocks its member
//! objects print and the float format `FFTPoisson::atomshift` leaves on the stream (`fixed`,
//! three decimals) for the rest of the log ([`LogFormat`]).

use super::fft::{slot, Fft3d, FftGrid};
use super::fft_fields::{
    inside_sol, integr_inside, ls_vacuum_field, rf_vacuum_field, IonLists, TINY,
};
use super::params::PbParameters;
use super::{extent, g, PbSystem};
use newen_core::Vec3;
use rustfft::num_complex::Complex64;

/// The float format of gromos++'s log stream: C++ default (`%g`), or `fixed` with a precision.
/// `FFTPoisson::atomshift` sets `fixed`/3 and nothing resets it (`FFTPoisson.cc:422-423`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LogFormat {
    fixed: Option<usize>,
}

impl LogFormat {
    /// A double as the stream prints it.
    pub fn f(self, x: f64) -> String {
        match self.fixed {
            None => g(x),
            Some(p) => {
                if x.is_nan() {
                    if x.is_sign_negative() { "-nan" } else { "nan" }.to_string()
                } else if x.is_infinite() {
                    if x > 0.0 { "inf" } else { "-inf" }.to_string()
                } else {
                    format!("{x:.p$}")
                }
            },
        }
    }
}

/// The `@verbose` text together with its stream's float format.
pub struct StreamLog<'a> {
    pub text: &'a mut String,
    pub fmt: LogFormat,
}

impl StreamLog<'_> {
    pub fn line(&mut self, s: impl AsRef<str>) {
        self.text.push_str(s.as_ref());
        self.text.push('\n');
    }

    /// `PB_Parameters(os)` (`None`, the default ε of 78.4) or `PB_Parameters(eps, os)`: every
    /// gromos++ PB object owns one and logs its constants on construction.
    fn pb_parameters(&mut self, eps: Option<f64>) {
        let f = self.fmt;
        match eps {
            None => self.line(
                "# initialising PB_Parameters using default epsilon of 78.4 for adjusted boundary conditions",
            ),
            Some(e) => self.line(format!(
                "# initialising PB_Parameters using epsilon of {} for adjusted boundary conditions",
                f.f(e)
            )),
        }
        let p = PbParameters::new(eps.unwrap_or(78.4), 0.0);
        let floats = |name: &str, v: f64| format!("# PB_PARAMS: {name} {}", f.f(v));
        let lines = [
            floats("epssolute", p.epssolute),
            floats("tiny_real", p.tiny_real),
            floats("tinfoil_boundary", p.tinfoil_boundary),
            floats("vacuum_boundary", p.vacuum_boundary),
            floats("adjusted_boundary", p.adjusted_boundary),
            floats("ddg_conv", p.ddg_conv),
            floats("hradius", p.hradius),
            floats("alpha1", p.alpha1),
            floats("alpha2", p.alpha2),
            format!("# PB_PARAMS: nalias1 {}", p.nalias1),
            format!("# PB_PARAMS: nalias2 {}", p.nalias2),
            floats("FFTlambda", p.fft_lambda),
            floats("kappa", p.kappa),
            floats("convergence_fd", p.convergence_fd),
            floats("convergence_fft", p.convergence_fft),
            floats("quadr", p.quadr),
            floats("tole", p.tole),
            floats("rcut", p.rcut),
            format!("# PB_PARAMS: kmax {}", p.kmax),
            floats("xiew", p.xiew),
        ];
        for l in lines {
            self.line(l);
        }
    }
}

/// The boundary condition of one FFT solve (gromos++ `FFTInteractionTypeCodes`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FftScheme {
    /// Type 0: lattice sum, tinfoil boundary, Ewald-type vacuum field with hat-shaped charges.
    Ls,
    /// Type 1: reaction field beyond `rcut` with permittivity `eps_rf` in the vacuum field.
    Rf,
}

/// Grid and permittivities of the FFT solves.
#[derive(Debug, Clone, Copy)]
pub struct FftSetup {
    /// Grid points per axis (the periodic FD grid, `@increasegrid` included).
    pub grid: [usize; 3],
    /// Grid spacing the program was given; the grid spans `grid·h`.
    pub h: f64,
    /// Relative permittivity of the solvent (`@epsSOLV`), the continuum of both solves and the
    /// boundary permittivity of the RF tensors.
    pub eps_solv: f64,
    /// Reaction-field permittivity (`@epsRF`, 0 = conducting) of the RF vacuum field.
    pub eps_rf: f64,
    /// Reaction-field cutoff (`@rcut`), nm.
    pub rcut: f64,
    /// Cubelets per half axis for the occupancy and the cutoff sphere (`@cubesFFT`).
    pub ncubes: usize,
    /// Iteration limit (`@maxiter`).
    pub maxiter: usize,
}

/// The result of one FFT solve.
#[derive(Debug, Clone)]
pub struct FftSolution {
    /// Reaction-field potential at every read-out atom, kJ mol⁻¹ e⁻¹.
    pub potentials: Vec<f64>,
    /// ½ Σ qᵢ φᵢ over the read-out atoms, kJ/mol.
    pub dg: f64,
    /// Iterations taken.
    pub steps: usize,
    /// Whether the RMS field inside the solute fell to `convergence_fft`.
    pub converged: bool,
}

/// The two FFT solves of `dGslv_pbsolv @schemeELEC RF` (`dGslv_pbsolv.cc:754-757`).
#[derive(Debug, Clone)]
pub struct DgslvFft {
    pub ls: FftSolution,
    pub rf: FftSolution,
    /// Solute positions after the LS solve centred `@atoms` on the grid and both wrapped them
    /// into it; the RF solve ran on these.
    pub positions: Vec<Vec3>,
}

/// `fft_ls_pbc` then `fft_rf_pbc` (`dGslv_pbsolv.cc:235-304`): the LS solve centres `@atoms` on the
/// grid, and the RF solve runs on the shifted positions. `epsRF = 1` (gromos++'s "SC" branch,
/// which never returns its potentials and aborts) is refused.
pub fn dgslv_fft(sys: &PbSystem, setup: &FftSetup, log: &mut String) -> Result<DgslvFft, String> {
    check_eps_rf(setup.eps_rf)?;
    let mut positions = sys.positions.clone();
    let mut stream = StreamLog {
        text: log,
        fmt: LogFormat::default(),
    };
    let ls = solve_fft(sys, &mut positions, setup, FftScheme::Ls, true, &mut stream)?;
    let rf = solve_fft(
        sys,
        &mut positions,
        setup,
        FftScheme::Rf,
        false,
        &mut stream,
    )?;
    Ok(DgslvFft { ls, rf, positions })
}

/// gromos++ treats `|epsRF − 1| ≤ tiny_real` as "no reaction-field correction" (SC, type 2) and
/// solves it without keeping the potentials (`dGslv_pbsolv.cc:280-299`); not supported.
fn check_eps_rf(eps_rf: f64) -> Result<(), String> {
    if (eps_rf - 1.0).abs() <= TINY {
        return Err(
            "@epsRF 1 (no reaction-field correction, gromos++'s SC branch) is not supported: \
             gromos++ solves it without returning the potentials and aborts"
                .into(),
        );
    }
    Ok(())
}

/// The interaction tensors of one boundary condition.
#[derive(Debug, Clone, Copy)]
enum Tensors {
    /// `FFTDipoleDipole_LS(0, εs)` and `FFTChargeDipole_LS`.
    Ls { es_1_es: f64, k0: f64 },
    /// `FFTDipoleDipole_RF(εB, rc, εs)` and `FFTChargeDipole_RF(εB, rc, εs)`.
    Rf {
        eps_b: f64,
        cut: f64,
        two_eps_prime_plus_one: f64,
        eps_s_minus_one: f64,
        upper: f64,
        two_es1_eprime: f64,
        es_two_eprime_plus_one: f64,
        k0: f64,
        first: f64,
        second: f64,
    },
}

/// The k = 0 factor of `FFTDipoleDipole_LS(εB, εs)` (`FFTDipoleDipole_LS.cc:55-61`).
fn k0_factor(eps_b: f64, eps_s: f64) -> f64 {
    if eps_b.abs() < TINY {
        1.0
    } else {
        (2.0 * eps_b + 1.0) / (2.0 * eps_b + eps_s)
    }
}

/// `3/x²·(−cos x + sin x/x)` as `FFTDipoleDipole_RF::fkr` writes it.
#[inline]
pub(crate) fn fkr(kr: f64) -> f64 {
    let kri = 1.0 / kr;
    let kr2i = kri * kri;
    3.0 * kr2i * (-kr.cos() + kr.sin() * kri)
}

impl Tensors {
    fn ls(eps_s: f64) -> Self {
        Tensors::Ls {
            es_1_es: (eps_s - 1.0) / eps_s,
            k0: k0_factor(0.0, eps_s),
        }
    }

    fn rf(eps_b: f64, cut: f64, eps_s: f64) -> Self {
        let two_eps_prime_plus_one = 2.0 * eps_b + 1.0;
        let eps_s_minus_one = eps_s - 1.0;
        let (first, second) = if eps_b.abs() < TINY {
            (1.0, 0.0)
        } else {
            let tepo = 2.0 * eps_b + 1.0;
            let tepm = 2.0 * (eps_b - 1.0);
            (tepm / tepo, 3.0 / tepo)
        };
        Tensors::Rf {
            eps_b,
            cut,
            two_eps_prime_plus_one,
            eps_s_minus_one,
            upper: eps_s_minus_one * two_eps_prime_plus_one * two_eps_prime_plus_one,
            two_es1_eprime: 2.0 * (eps_s - 1.0) * eps_b,
            es_two_eprime_plus_one: eps_s * two_eps_prime_plus_one,
            k0: k0_factor(eps_b, eps_s),
            first,
            second,
        }
    }

    /// `FFTDipoleDipole::updateTensor` (`FFTDipoleDipole.cc:67-106`) on `tt = k kᵀ`.
    fn update(&self, k2: f64, tt: &mut [[f64; 3]; 3]) {
        let k0 = match *self {
            Tensors::Ls { k0, .. } | Tensors::Rf { k0, .. } => k0,
        };
        if k2 < TINY {
            *tt = [[0.0; 3]; 3];
            for (ii, row) in tt.iter_mut().enumerate() {
                row[ii] = k0;
            }
            return;
        }
        let (a, b) = match *self {
            Tensors::Ls { es_1_es, .. } => (1.0, -es_1_es),
            Tensors::Rf {
                eps_b,
                cut,
                two_eps_prime_plus_one: tepo,
                eps_s_minus_one: esm1,
                upper,
                two_es1_eprime,
                es_two_eprime_plus_one,
                ..
            } => {
                // computeAFactor caches f(kr) for computeBFactor.
                let fckr = fkr(k2.sqrt() * cut);
                let a = if 0.0 == eps_b {
                    1.0
                } else {
                    tepo / (tepo + esm1 * fckr)
                };
                let b = if 0.0 == eps_b {
                    -esm1 * (1.0 - fckr) / (1.0 + esm1 * (1.0 - fckr))
                } else {
                    -upper * (fckr - 1.0)
                        / ((tepo + esm1 * fckr) * (two_es1_eprime * fckr - es_two_eprime_plus_one))
                };
                (a, b)
            },
        };
        let b = b / k2;
        for row in tt.iter_mut() {
            for v in row.iter_mut() {
                *v *= b;
            }
        }
        for (ii, row) in tt.iter_mut().enumerate() {
            row[ii] += a;
        }
    }

    /// `FFTChargeDipole::polarization` (`FFTChargeDipole_RF.cc:98-102`; 1 for LS).
    fn polarization(&self, k2: f64) -> f64 {
        match *self {
            Tensors::Ls { .. } => 1.0,
            Tensors::Rf {
                cut, first, second, ..
            } => {
                let kr = k2.sqrt() * cut;
                1.0 - first * fkr(kr) - second * kr.sin() / kr
            },
        }
    }
}

const STARS: &str = "# ************************************************** ";

/// The `FFTBoundaryCondition::dumpparameters` lines (`FFTBoundaryCondition.cc:51-61`).
fn dump_boundary(log: &mut StreamLog, ty: i32, name: &str, ppp: &PbParameters) {
    let f = log.fmt;
    log.line("# BOUNDARY PARAMETERS");
    log.line("# -------------------");
    log.line(format!("# TYPE {ty}"));
    log.line(format!("# STYPE {name}"));
    log.line("# HAT CHARGE SHAPING FUNCTION");
    log.line(format!("# ALPHA1: {}", f.f(ppp.alpha1)));
    log.line(format!("# ALPHA2: {}", f.f(ppp.alpha2)));
    log.line(format!("# NALIAS1: {}", ppp.nalias1));
    log.line(format!("# NALIAS2: {}", ppp.nalias2));
}

/// One FFT solve as `dGslv_pbsolv`'s `fft_ls_pbc` / `fft_rf_pbc` run it: the grid and boundary
/// objects, `FFTPoisson` (centring `@atoms` on the grid when `shift_atoms`) and `solve_poisson`.
/// `positions` holds every solute atom and is changed the way gromos++ changes the shared
/// `System`: centring and wrapping move the `@atoms` entries.
pub fn solve_fft(
    sys: &PbSystem,
    positions: &mut [Vec3],
    setup: &FftSetup,
    scheme: FftScheme,
    shift_atoms: bool,
    log: &mut StreamLog,
) -> Result<FftSolution, String> {
    if setup.grid.contains(&0) || !(setup.h > 0.0 && setup.h.is_finite()) {
        return Err(format!(
            "FFT grid {:?} with spacing {}: every axis needs points and a positive spacing",
            setup.grid, setup.h
        ));
    }
    if setup.grid.iter().product::<usize>() <= 100 {
        // gromos++ logs element 100 of the grids (FFTPoisson.cc:341-346).
        return Err(format!(
            "FFT grid {:?}: gromos++ needs more than 100 grid points",
            setup.grid
        ));
    }
    if setup.eps_solv < 1.0 {
        return Err("the solvent permittivity of an FFT solve must be at least 1".into());
    }
    if setup.ncubes == 0 || setup.maxiter == 0 {
        return Err("@cubesFFT and @maxiter must be positive".into());
    }
    if scheme == FftScheme::Rf {
        check_eps_rf(setup.eps_rf)?;
    }
    let ppp = PbParameters::new(setup.eps_solv, sys.fpepsi);
    let eps = setup.eps_solv;
    let len = setup.grid.map(|n| n as f64 * setup.h);
    let gt = FftGrid::new(setup.grid, len, setup.ncubes, ppp.pi);

    // fft_ls_pbc / fft_rf_pbc (dGslv_pbsolv.cc:235-304)
    log.pb_parameters(None); // FFTGridType
    log.line(STARS);
    let (ty, name) = match scheme {
        FftScheme::Ls => {
            log.line("# *** FFT LS (PBC) *** ");
            (0, "LS")
        },
        FftScheme::Rf => {
            log.line("# *** FFT RF (PBC) *** ");
            (1, "RF")
        },
    };
    log.line(STARS);
    if scheme == FftScheme::Rf {
        log.line("# WITH REACTION FIELD CORRECTION : RF ... ");
    }
    log.pb_parameters(None); // FFTBoundaryCondition
    if scheme == FftScheme::Ls {
        dump_boundary(log, ty, name, &ppp);
    }

    // FFTPoisson: members ppp, pbiterator, bc, gt in declaration order (FFTPoisson.h).
    log.pb_parameters(Some(eps));
    let tensors = iterator_setup(log, scheme, eps, setup);
    log.pb_parameters(None);
    log.pb_parameters(None);
    let f = log.fmt;
    log.line(format!("# epssolvent {}", f.f(eps)));
    for (i, &a) in sys.readout.iter().enumerate() {
        log.line(format!("# atom {i} charge {}", f.f(sys.charges[a])));
        log.line(format!("# atom {i} radius {}", f.f(sys.radii[a])));
    }
    if shift_atoms {
        log.line("# FFTPoisson setup ... center atoms ");
        center_atoms_on_grid(sys, positions, gt.center, log);
    }
    if scheme == FftScheme::Rf {
        dump_boundary(log, ty, name, &ppp);
    }
    log.line("# call solve_poisson ...");

    // solve_poisson (FFTPoisson.cc:164-251)
    let mut radii: Vec<f64> = sys.atoms.iter().map(|&a| sys.radii[a]).collect();
    atomshift(sys, positions, &mut radii, &gt, setup.h, log);
    gridcheck(sys, positions, &radii, &gt, log)?;

    log.pb_parameters(None); // FFTInsideOutside
    log.line("# FFTPoisson setup ...  insideoutside");
    // inside_sol reads the topology radii, not the ones atomshift shrank.
    let mut inside = inside_sol(&gt, gt.ncubes, positions, &sys.atoms, &sys.radii);
    log.line("# FFTPoisson setup ... done insideoutside");
    log.line("# FFTPoisson setup ... grid");
    let mut fft = Fft3d::new(gt);
    let zero = vec![Complex64::default(); gt.ngr3];
    let mut v = [zero.clone(), zero.clone(), zero.clone()];
    log.line("# FFTPoisson setup ... done gridresizing");
    log.line("# FFTPoisson setup ... setup vacfield ");
    setup_vacuum_field(
        sys, positions, &gt, &ppp, scheme, setup, &inside, &mut v, &mut fft, log,
    )?;
    integr_inside(&mut inside, log);

    iterate(
        sys, positions, &gt, &ppp, &tensors, setup, &inside, &mut v, &mut fft, log,
    )
}

/// The `FFTPoissonIterator` constructor's log and its tensors (`FFTPoissonIterator.cc:56-116`).
fn iterator_setup(log: &mut StreamLog, scheme: FftScheme, eps: f64, setup: &FftSetup) -> Tensors {
    // members ppp(epssolvent), bc, gt
    log.pb_parameters(Some(eps));
    log.pb_parameters(None);
    log.pb_parameters(None);
    let ppp = PbParameters::new(eps, 0.0);
    let f = log.fmt;
    log.line(format!(
        "# FFTPoissonIterator: lambda {}",
        f.f(ppp.fft_lambda)
    ));
    log.line(format!("# FFTPoissonIterator: maxsteps {}", setup.maxiter));
    log.line(format!(
        "# FFTPoissonIterator: ls_boundary_eps {}",
        f.f(ppp.tinfoil_boundary)
    ));
    log.line(format!(
        "# FFTPoissonIterator: rf_boundary_eps {}",
        f.f(ppp.adjusted_boundary)
    ));
    log.line(format!(
        "# FFTPoissonIterator: sc_boundary_eps {}",
        f.f(ppp.vacuum_boundary)
    ));
    log.line(format!("# FFTPoissonIterator: epssolvent {}", f.f(eps)));
    let dipole_dipole = |log: &mut StreamLog, eps_b: f64| {
        log.pb_parameters(Some(eps));
        log.line(format!("# FFTDipoleDipole: epsS {}", f.f(eps)));
        log.line(format!("# FFTDipoleDipole: epsB {}", f.f(eps_b)));
    };
    let tensors = match scheme {
        FftScheme::Ls => {
            log.line("# FFTPoissonIterator: get InterxTensor for LS ");
            dipole_dipole(log, ppp.tinfoil_boundary);
            log.pb_parameters(Some(eps)); // FFTChargeDipole_LS
            Tensors::ls(eps)
        },
        FftScheme::Rf => {
            // The RF tensors take the solvent permittivity as boundary (rf_boundary_eps), not
            // @epsRF, which only enters the vacuum field (FFTPoissonIterator.cc:91-94).
            let eps_b = ppp.adjusted_boundary;
            log.line("# FFTPoissonIterator: get InterxTensor for RF ");
            dipole_dipole(log, eps_b);
            dipole_dipole(log, eps_b); // the embedded FFTDipoleDipole_LS for k = 0
            log.line(format!("# FFTDipoleDipole_RF : epsB {}", f.f(eps_b)));
            log.line(format!("# FFTDipoleDipole_RF : epsS {}", f.f(eps)));
            log.pb_parameters(Some(eps));
            log.line(format!("# FFTChargeDipole_RF: epssolvent {}", f.f(eps)));
            log.line(format!("# FFTChargeDipole_RF: epsB {}", f.f(eps_b)));
            Tensors::rf(eps_b, setup.rcut, eps)
        },
    };
    log.line("# FFTPoissonIterator setup ... make fftplans ");
    tensors
}

/// `FFTPoisson::center_atoms_on_grid` (`FFTPoisson.cc:635-683`): moves `@atoms` by the grid centre
/// minus the centre of the bounding box of all solute atoms.
fn center_atoms_on_grid(
    sys: &PbSystem,
    positions: &mut [Vec3],
    center: [f64; 3],
    log: &mut StreamLog,
) {
    let mmcenter = |positions: &[Vec3]| {
        let (lo, hi) = extent(positions);
        [
            (lo.x + hi.x) * 0.5,
            (lo.y + hi.y) * 0.5,
            (lo.z + hi.z) * 0.5,
        ]
    };
    let f = log.fmt;
    let m = mmcenter(positions);
    log.line(format!(
        "# initial min max center x,y,z = {} {} {}",
        f.f(m[0]),
        f.f(m[1]),
        f.f(m[2])
    ));
    let tvec = [center[0] - m[0], center[1] - m[1], center[2] - m[2]];
    log.line(format!(
        "translation vector x,y,z = {} {} {}",
        f.f(tvec[0]),
        f.f(tvec[1]),
        f.f(tvec[2])
    ));
    for &a in &sys.atoms {
        positions[a].x += tvec[0];
        positions[a].y += tvec[1];
        positions[a].z += tvec[2];
    }
    log.line(format!(
        "# gridcenter x,y,z = {} {} {}",
        f.f(center[0]),
        f.f(center[1]),
        f.f(center[2])
    ));
    let m = mmcenter(positions);
    log.line(format!(
        "# final mmcenterAtoms (after centering on grid) x,y,z = {} {} {}",
        f.f(m[0]),
        f.f(m[1]),
        f.f(m[2])
    ));
}

/// `FFTPoisson::atomshift` (`FFTPoisson.cc:418-590`): wraps `@atoms` positions lying outside the
/// grid back into it and shrinks radii protruding from it. The shrunk radii only reach
/// `gridcheck`; the occupancy uses the topology radii (gromos++ does the same).
fn atomshift(
    sys: &PbSystem,
    positions: &mut [Vec3],
    radii: &mut [f64],
    gt: &FftGrid,
    h: f64,
    log: &mut StreamLog,
) {
    for (i, &a) in sys.atoms.iter().enumerate() {
        let p = positions[a];
        log.line(format!(
            "# GRIDCHECK for atom {} radius {:.7} charge {:>6.3} pos {:.3} {:.3} {:.3}",
            i + 1,
            radii[i],
            sys.charges[a],
            p.x,
            p.y,
            p.z
        ));
        log.fmt = LogFormat { fixed: Some(3) };
        let f = log.fmt;
        for (d, axis) in ['X', 'Y', 'Z'].into_iter().enumerate() {
            let len = gt.len[d];
            let shifted = |log: &mut StreamLog, p: Vec3| {
                log.line(format!(
                    "# Needed to shift atom {} on the {axis} axis to the other side of the box ",
                    i + 1
                ));
                log.line(format!(
                    "# New coordinates atom {} {} {} {}",
                    i + 1,
                    f.f(p.x),
                    f.f(p.y),
                    f.f(p.z)
                ));
            };
            let shrunk = |log: &mut StreamLog, r: f64| {
                log.line(format!(
                    "# Needed to shrink the radius of atom {} on the {axis} axis since it exceeded the grid ",
                    i + 1
                ));
                log.line(format!("# New radius atom {} {}", i + 1, f.f(r)));
            };
            let upper_radius = |log: &mut StreamLog, pos: f64, r: &mut f64| {
                if pos + *r > len {
                    let protruding = pos + *r - len + h / 100.0;
                    *r -= protruding;
                    shrunk(log, *r);
                }
            };
            if positions[a][d] > len {
                let protruding = positions[a][d] - len;
                positions[a][d] = 0.0 + protruding;
                shifted(log, positions[a]);
            }
            upper_radius(log, positions[a][d], &mut radii[i]);
            if positions[a][d] < 0.0 {
                let protruding = positions[a][d] - 0.0;
                positions[a][d] = len + protruding;
                shifted(log, positions[a]);
            }
            if positions[a][d] - radii[i] < 0.0 {
                let protruding = positions[a][d] - radii[i] - 0.0 - h / 100.0;
                radii[i] += protruding;
                shrunk(log, radii[i]);
            }
            upper_radius(log, positions[a][d], &mut radii[i]);
        }
    }
}

/// `FFTPoisson::gridcheck` (`FFTPoisson.cc:592-630`).
fn gridcheck(
    sys: &PbSystem,
    positions: &[Vec3],
    radii: &[f64],
    gt: &FftGrid,
    log: &mut StreamLog,
) -> Result<(), String> {
    for (i, &a) in sys.atoms.iter().enumerate() {
        let (p, rad) = (positions[a], radii[i]);
        if p.x + rad > gt.len[0]
            || p.x - rad < 0.0
            || p.y + rad > gt.len[1]
            || p.y - rad < 0.0
            || p.z + rad > gt.len[2]
            || p.z - rad < 0.0
        {
            let f = log.fmt;
            log.line("# ATOM EXTENDING GRID!");
            log.line(format!("# ATOM ID:     {}", i + 1));
            log.line(format!(
                "# ATOM POS:    {} {} {}",
                f.f(p.x),
                f.f(p.y),
                f.f(p.z)
            ));
            log.line(format!("# ATOM RADIUS: {}", f.f(rad)));
            log.line(format!(
                "# GRID START   {} {} {}",
                f.f(0.0),
                f.f(0.0),
                f.f(0.0)
            ));
            log.line(format!(
                "# GRID BOX     {} {} {}",
                f.f(gt.len[0]),
                f.f(gt.len[1]),
                f.f(gt.len[2])
            ));
            return Err("FDPoissonBolzmann: Atom extending grid. Exiting ...".into());
        }
    }
    Ok(())
}

/// `FFTPoisson::setupVacuumField` (`FFTPoisson.cc:263-349`): the vacuum field in real space, set to
/// zero inside the solute by the fractional occupancy.
#[allow(clippy::too_many_arguments)]
fn setup_vacuum_field(
    sys: &PbSystem,
    positions: &[Vec3],
    gt: &FftGrid,
    ppp: &PbParameters,
    scheme: FftScheme,
    setup: &FftSetup,
    inside: &[f64],
    v: &mut [Vec<Complex64>; 3],
    fft: &mut Fft3d,
    log: &mut StreamLog,
) -> Result<(), String> {
    log.line("# setting up vacuum field ... ");
    let ty = match scheme {
        FftScheme::Ls => 0,
        FftScheme::Rf => 1,
    };
    log.line(format!("# bc.type is: {ty}"));
    match scheme {
        FftScheme::Ls => {
            log.line("# need LS field ");
            // FFTVacuumField members gt, bc, ppp; FFTVacuumField_LS's own ppp; the generic
            // FFTChargeShapingFunction; then the hat function.
            for _ in 0..5 {
                log.pb_parameters(None);
            }
            log.line(format!("# FFTVacuumField : cstype {}", ppp.default_chshape));
            log.line("# FFTVacuumField : csfunc->HatType 0");
            log.line("# FFTVacuumField : csfunc->ParabolaType 1");
            log.pb_parameters(None);
            let lists = IonLists::new(sys, ppp.alpha1, ppp.alpha2, log)?;
            let eps0 = 1.0 / (sys.fpepsi * 4.0 * ppp.pi);
            ls_vacuum_field(
                gt,
                sys,
                positions,
                &lists,
                [ppp.alpha1, ppp.alpha2],
                [ppp.nalias1 as i64, ppp.nalias2 as i64],
                eps0,
                ppp.pi,
                v,
                log,
            );
            for c in v.iter_mut() {
                fft.backward(c);
            }
        },
        FftScheme::Rf => {
            log.line("# need RF field ");
            for _ in 0..3 {
                log.pb_parameters(None);
            }
            rf_vacuum_field(gt, sys, positions, setup.eps_rf, setup.rcut, v, log)?;
        },
    }
    for c in v.iter_mut() {
        for (x, &ins) in c.iter_mut().zip(inside) {
            let heaviside = 1.0 - ins;
            *x = Complex64::new(x.re * heaviside, x.im * heaviside);
        }
    }
    let f = log.fmt;
    for m in [0, 100] {
        for (c, name) in ["vx", "vy", "vz"].iter().enumerate() {
            log.line(format!(
                "# after setupfield: {name}[{m}] {}",
                f.f(slot(&v[c], m))
            ));
        }
    }
    Ok(())
}

/// `FFTPoissonIterator::iterate_poisson` and `postIteration` (`FFTPoissonIterator.cc:120-290`).
#[allow(clippy::too_many_arguments)]
fn iterate(
    sys: &PbSystem,
    positions: &[Vec3],
    gt: &FftGrid,
    ppp: &PbParameters,
    tensors: &Tensors,
    setup: &FftSetup,
    inside: &[f64],
    v: &mut [Vec<Complex64>; 3],
    fft: &mut Fft3d,
    log: &mut StreamLog,
) -> Result<FftSolution, String> {
    let [nx, ny, nz] = gt.n;
    let (kx, ky, kz) = (gt.k_vectors(0), gt.k_vectors(1), gt.k_vectors(2));

    // The tensors depend on k only: built once instead of once per iteration.
    let mut tt = Vec::with_capacity(gt.ngr3);
    for &kxv in &kx {
        for &kyv in &ky {
            for &kzv in &kz {
                let k = [kxv, kyv, kzv];
                let k2 = kxv * kxv + kyv * kyv + kzv * kzv;
                let mut t = [[0.0; 3]; 3];
                for r in 0..3 {
                    for c in 0..3 {
                        t[r][c] = k[r] * k[c];
                    }
                }
                tensors.update(k2, &mut t);
                tt.push(t);
            }
        }
    }

    let zero = vec![Complex64::default(); gt.ngr3];
    let mut e = [zero.clone(), zero.clone(), zero.clone()];
    let mut lambda = ppp.fft_lambda;
    let convergence = ppp.convergence_fft;
    let mut steps = 1usize;
    let mut converged = false;
    let mut e_inside = 0.0;
    let trio = |log: &mut StreamLog, tag: &str, g: &[Vec<Complex64>; 3]| {
        let f = log.fmt;
        for (c, (axis, n)) in [('X', nx), ('Y', ny), ('Z', nz)].into_iter().enumerate() {
            log.line(format!(
                "# {tag}{axis} {} {} {}",
                f.f(slot(&g[c], 0)),
                f.f(slot(&g[c], n / 2)),
                f.f(slot(&g[c], n - 1))
            ));
        }
    };

    log.line("# START ITERATION ... ");
    loop {
        let f = log.fmt;
        log.line(format!("# ITER {steps}"));
        log.line(format!("# lambda {}", f.f(lambda)));
        for c in 0..3 {
            fft.forward(&v[c], &mut e[c]);
        }
        log.line("# done FT vacuum field ");
        trio(log, "FT_VF", v);
        for (l, t) in tt.iter().enumerate() {
            let er = [e[0][l].re, e[1][l].re, e[2][l].re];
            let ei = [e[0][l].im, e[1][l].im, e[2][l].im];
            for r in 0..3 {
                let mut re = 0.0;
                let mut im = 0.0;
                for u in 0..3 {
                    re += t[r][u] * er[u];
                    im += t[r][u] * ei[u];
                }
                e[r][l] = Complex64::new(re, im);
            }
        }
        log.line("# done compute Efield from vacuum field ");
        trio(log, "FT_EF", &e);
        for c in e.iter_mut() {
            fft.backward(c);
        }
        trio(log, "RS_EF", &e);

        let prev = e_inside;
        e_inside = residual(&e, inside);
        let delta = if steps == 1 { 0.0 } else { prev - e_inside };
        if e_inside <= convergence {
            converged = true;
            break;
        }
        if steps == setup.maxiter {
            break;
        }
        // updateVacuumField (FFTPoissonIterator.cc:295-330)
        let f = log.fmt;
        log.line(format!("# deltaSigma = {}", f.f(delta)));
        if delta >= 0.0 {
            lambda *= 1.2;
        } else {
            lambda *= 0.5;
        }
        log.line(format!("# Continuing with lambda: {}", f.f(lambda)));
        for (l, &ins) in inside.iter().enumerate() {
            if 0.0 != ins {
                let lam_fact = -lambda * ins;
                for c in 0..3 {
                    let (vv, ee) = (v[c][l], e[c][l]);
                    v[c][l] = Complex64::new(vv.re + lam_fact * ee.re, vv.im + lam_fact * ee.im);
                }
            }
        }
        trio(log, "UP_VF", v);
        steps += 1;
    }

    // postIteration
    for c in e.iter_mut() {
        fft.forward_in_place(c);
    }
    trio(log, "POST_FT_EF", &e);
    let es_1 = setup.eps_solv - 1.0;
    let mut rpot = vec![Complex64::default(); gt.ngr3];
    let mut l = 0;
    for &kxv in &kx {
        let kx2 = kxv * kxv;
        for &kyv in &ky {
            let ky2 = kyv * kyv;
            for &kzv in &kz {
                let k2 = kx2 + ky2 + kzv * kzv;
                if k2 > TINY {
                    let pp = tensors.polarization(k2) * es_1 / k2;
                    rpot[l] = Complex64::new(
                        pp * (kxv * e[0][l].im + kyv * e[1][l].im + kzv * e[2][l].im),
                        -pp * (kxv * e[0][l].re + kyv * e[1][l].re + kzv * e[2][l].re),
                    );
                }
                l += 1;
            }
        }
    }
    for c in e.iter_mut() {
        fft.backward(c);
    }
    trio(log, "POST_RS_EF", &e);
    fft.backward(&mut rpot);

    let head = if converged {
        "# Convergence after"
    } else {
        "# Terminated after"
    };
    // gromos++ writes the head, then free_energy_restricted logs the atoms, then the value.
    log.text
        .push_str(&format!("{head} {steps} steps with dG [kJ/mol]: "));
    let (potentials, dg) = free_energy_restricted(sys, positions, gt, &rpot, log)?;
    let f = log.fmt;
    log.line(f.f(dg));
    Ok(FftSolution {
        potentials,
        dg,
        steps,
        converged,
    })
}

/// `computeResidualField` (`FFTPoissonIterator.cc:559-618`): the RMS field over the solute points,
/// `sqrt(Σ inside·|E|² / count)` with gromos++'s `int` count.
fn residual(e: &[Vec<Complex64>; 3], inside: &[f64]) -> f64 {
    let mut e_ave = 0.0;
    let mut count: i32 = 0;
    for (l, &ins) in inside.iter().enumerate() {
        if ins > TINY {
            let (x, y, z) = (e[0][l], e[1][l], e[2][l]);
            let e_abs =
                x.re * x.re + x.im * x.im + y.re * y.re + y.im * y.im + z.re * z.re + z.im * z.im;
            e_ave += e_abs * ins;
            count = (count as f64 + ins) as i32;
        }
    }
    (e_ave / count as f64).sqrt()
}

/// `free_energy_restricted` (`FFTPoissonIterator.cc:695-769`): the real part of `pot` trilinearly
/// interpolated at each read-out atom (z-fastest corners, the `+1` corners not wrapped), and
/// `½ Σ q φ`.
fn free_energy_restricted(
    sys: &PbSystem,
    positions: &[Vec3],
    gt: &FftGrid,
    pot: &[Complex64],
    log: &mut StreamLog,
) -> Result<(Vec<f64>, f64), String> {
    let [_, ny, nz] = gt.n;
    let mut dg = 0.0;
    let mut potentials = Vec::with_capacity(sys.readout.len());
    for (ion, &a) in sys.readout.iter().enumerate() {
        let p = positions[a];
        let cell = |x: f64, len: f64, d: f64| {
            let c = (x - (x / len).floor() * len) / d;
            let i = c as i64;
            (i, c - i as f64)
        };
        let (i, fx) = cell(p.x, gt.len[0], gt.dr[0]);
        let (j, fy) = cell(p.y, gt.len[1], gt.dr[1]);
        let (k, fz) = cell(p.z, gt.len[2], gt.dr[2]);
        let (ny, nz) = (ny as i64, nz as i64);
        let s = |di: i64, dj: i64, dk: i64| k + dk + nz * (j + dj + ny * (i + di));
        let last = s(1, 1, 1);
        if s(0, 0, 0) < 0 || last >= gt.ngr3 as i64 {
            return Err(format!(
                "FFT read-out: atom {} at ({}, {}, {}) is too close to the grid edge to \
                 interpolate (gromos++ reads outside its grid)",
                a + 1,
                p.x,
                p.y,
                p.z
            ));
        }
        let at = |di, dj, dk| pot[s(di, dj, dk) as usize].re;
        let pv = (1.0 - fx) * (1.0 - fy) * (1.0 - fz) * at(0, 0, 0)
            + (1.0 - fx) * (1.0 - fy) * fz * at(0, 0, 1)
            + (1.0 - fx) * fy * (1.0 - fz) * at(0, 1, 0)
            + (1.0 - fx) * fy * fz * at(0, 1, 1)
            + fx * (1.0 - fy) * (1.0 - fz) * at(1, 0, 0)
            + fx * (1.0 - fy) * fz * at(1, 0, 1)
            + fx * fy * (1.0 - fz) * at(1, 1, 0)
            + fx * fy * fz * at(1, 1, 1);
        potentials.push(pv);
        let f = log.fmt;
        log.line(format!(
            "# atom  {ion} charge {} ele. pot. {}",
            f.f(sys.charges[a]),
            f.f(pv)
        ));
        dg += sys.charges[a] * pv;
    }
    Ok((potentials, 0.5 * dg))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rf_fkr_tends_to_one_at_small_kr() {
        // f(x) = 3(sin x − x cos x)/x³ = 1 − x²/10 + O(x⁴).
        for &x in &[1e-2f64, 3e-2, 1e-1] {
            let expected = 1.0 - x * x / 10.0 + x.powi(4) / 280.0;
            assert!((fkr(x) - expected).abs() < 1e-6, "{x}: {}", fkr(x));
        }
        // At large x it decays like −3cos(x)/x².
        let x = 200.0f64;
        assert!((fkr(x) + 3.0 * x.cos() / (x * x)).abs() < 1e-6);
    }

    #[test]
    fn tensors_reduce_to_their_limits() {
        // LS: A = 1, B = −(εs−1)/εs; the tensor on k̂ is 1/εs.
        let t = Tensors::ls(66.6);
        let k = [0.3, -0.4, 1.2];
        let k2 = k.iter().map(|x| x * x).sum::<f64>();
        let mut tt = [[0.0; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                tt[r][c] = k[r] * k[c];
            }
        }
        t.update(k2, &mut tt);
        let khat: Vec<f64> = k.iter().map(|x| x / k2.sqrt()).collect();
        let along: f64 = (0..3)
            .map(|r| (0..3).map(|c| khat[r] * tt[r][c] * khat[c]).sum::<f64>())
            .sum();
        assert!((along - 1.0 / 66.6).abs() < 1e-14);
        // k = 0: the conducting K0 factor.
        let mut tt = [[5.0; 3]; 3];
        t.update(0.0, &mut tt);
        assert_eq!(tt, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        // RF with εB = εs: the K0 factor (2εB+1)/(2εB+εs).
        let t = Tensors::rf(66.6, 1.4, 66.6);
        let mut tt = [[0.0; 3]; 3];
        t.update(0.0, &mut tt);
        assert_eq!(tt[0][0], (2.0 * 66.6 + 1.0) / (2.0 * 66.6 + 66.6));
        // As kr → 0 the RF polarization tends to 1 − 2(εB−1)/(2εB+1) − 3/(2εB+1) = 0.
        assert!(t.polarization(1e-8).abs() < 1e-6);
        assert_eq!(Tensors::ls(66.6).polarization(4.0), 1.0);
    }

    #[test]
    fn fixed_format_follows_the_cpp_stream() {
        let f = LogFormat { fixed: Some(3) };
        assert_eq!(f.f(-0.0), "-0.000");
        assert_eq!(f.f(-1e-5), "-0.000");
        assert_eq!(f.f(66.6), "66.600");
        assert_eq!(LogFormat::default().f(66.6), "66.6");
        assert_eq!(f.f(f64::NAN), "nan");
    }
}
