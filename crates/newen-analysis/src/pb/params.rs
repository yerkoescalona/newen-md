//! gromos++ `pb::PB_Parameters` (`src/pb/PB_Parameters.cc`): the fixed constants of the PB
//! solvers. None of them can be set from the command line.

use super::g;

/// gromos++'s 1/(4πε₀) when no topology supplies one (`gmath/Physics.cc:53`), kJ mol⁻¹ nm e⁻².
pub const GROMOSPP_FOUR_PI_EPS_I: f64 = 138.935484611;

/// The constants of one gromos++ `PB_Parameters` object (`PB_Parameters.cc:52-84`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PbParameters {
    pub epssolute: f64,
    pub tiny_real: f64,
    pub fpepsi: f64,
    pub pi: f64,
    pub tinfoil_boundary: f64,
    pub vacuum_boundary: f64,
    /// The solvent permittivity the object was built with.
    pub adjusted_boundary: f64,
    pub ddg_conv: f64,
    pub hradius: f64,
    pub nmax_surfpoints: i32,
    pub default_chshape: i32,
    pub alpha1: f64,
    pub alpha2: f64,
    pub nalias1: i32,
    pub nalias2: i32,
    pub threadnum: i32,
    pub fft_lambda: f64,
    /// Inverse Debye length of the NPBC boundary potential; 0, so the boundary is Coulombic.
    pub kappa: f64,
    pub debugvar: i32,
    /// Relative L1 residual at which the ICCG solver stops.
    pub convergence_fd: f64,
    pub convergence_fft: f64,
    pub tole: f64,
    pub rcut: f64,
    pub kmax: i32,
    pub quadr: f64,
    pub xiew: f64,
}

impl PbParameters {
    /// `PB_Parameters(epssolvent, os)` without the log lines; `fpepsi` is the value gromos++'s
    /// `physConst` holds at that moment (the topology's, or [`GROMOSPP_FOUR_PI_EPS_I`]).
    pub fn new(epssolvent: f64, fpepsi: f64) -> Self {
        PbParameters {
            epssolute: 1.0,
            tiny_real: 0.0000000001,
            fpepsi,
            // gmath::physConst.get_pi()
            pi: 4.0 * 1.0f64.atan(),
            tinfoil_boundary: 0.0,
            vacuum_boundary: 1.0,
            adjusted_boundary: epssolvent,
            ddg_conv: -0.001,
            hradius: 0.05,
            nmax_surfpoints: 15000,
            default_chshape: 0,
            alpha1: 0.15,
            alpha2: 0.05,
            nalias1: 2,
            nalias2: 8,
            threadnum: 1,
            fft_lambda: 1.5,
            kappa: 0.0,
            debugvar: 0,
            convergence_fd: 0.000001,
            convergence_fft: 0.3,
            tole: 0.000001,
            rcut: 0.9,
            kmax: 64,
            quadr: 0.0082,
            xiew: -2.837297,
        }
    }

    /// `PB_Parameters(epssolvent, os)`: the object and the 21 lines its constructor logs.
    pub fn new_logged(epssolvent: f64, fpepsi: f64, log: &mut String) -> Self {
        log.push_str(&format!(
            "# initialising PB_Parameters using epsilon of {} for adjusted boundary conditions\n",
            g(epssolvent)
        ));
        let p = Self::new(epssolvent, fpepsi);
        let lines: [(&str, f64); 20] = [
            ("epssolute", p.epssolute),
            ("tiny_real", p.tiny_real),
            ("tinfoil_boundary", p.tinfoil_boundary),
            ("vacuum_boundary", p.vacuum_boundary),
            ("adjusted_boundary", p.adjusted_boundary),
            ("ddg_conv", p.ddg_conv),
            ("hradius", p.hradius),
            ("alpha1", p.alpha1),
            ("alpha2", p.alpha2),
            ("nalias1", p.nalias1 as f64),
            ("nalias2", p.nalias2 as f64),
            ("FFTlambda", p.fft_lambda),
            ("kappa", p.kappa),
            ("convergence_fd", p.convergence_fd),
            ("convergence_fft", p.convergence_fft),
            ("quadr", p.quadr),
            ("tole", p.tole),
            ("rcut", p.rcut),
            ("kmax", p.kmax as f64),
            ("xiew", p.xiew),
        ];
        for (name, v) in lines {
            log.push_str(&format!("# PB_PARAMS: {name} {}\n", g(v)));
        }
        p
    }
}
