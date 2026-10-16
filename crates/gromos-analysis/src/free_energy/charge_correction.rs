//! The analytic correction terms for the charging free energy of an ion simulated under periodic
//! boundary conditions, after Kastenholz and Hünenberger[^kh2006]: the companions of
//! `dGslv_pbsolv`.
//!
//! The raw charging free energy ΔG_chg of an explicit-solvent simulation (TI, BAR or MBAR over
//! the charge) depends on the box and on the electrostatic scheme. The paper corrects it to the
//! value for a macroscopic, non-periodic, Coulombic system (Eqs. 9, 16, 30, 48, 56):
//!
//! ΔG_chg^cor = ΔG_chg + ΔG_A1 + ΔG_A2 + ΔG_B + ΔG_C1 + ΔG_C2
//!
//! for a spherical ion of charge q_I and radius R_I in a cubic box of edge L filled with a
//! solvent of permittivity ε_S:
//!
//! | term | what it corrects | LS | CM, BM, BA |
//! |------|------------------|----|------------|
//! | A1 | no ion–solvent interaction beyond the cut-off | 0 (Eq. 24) | Born, Eq. 25 |
//! | A2 | wrong solvent polarisation inside the cut-off | 0 (Eq. 24) | integral equation, Eqs. 26–29 |
//! | B  | finite size and periodicity | Eq. 32 | empirical fit, Eq. 33 |
//! | C1 | potential summed per atom, not per molecule (the solvent's quadrupole-moment trace) | Eq. 37 | Eqs. 39, 41, 42 |
//! | C2 | interface potential (finite-system effect) | SPC fit, Eq. 47 | 0 (Eq. 46) |
//!
//! A1, A2 and B are quadratic in the charge, C1 and C2 linear. Every term is **added** to the raw
//! ΔG_chg. For a change of the charge from q_A to q_B each term is term(q_B) − term(q_A), so the
//! quadratic ones scale with q_B² − q_A², not with (q_B − q_A)². The corrected value is the
//! intrinsic one: it excludes the surface term ΔG_srf = q_I Δχ° (Eqs. 11, 49), the cavity term
//! and the standard-state term (Eq. 15).
//!
//! **What `dGslv_pbsolv` already contains.** Its correction, CB/NPBC − LS/PBC (+ LS/PBC[FFT] −
//! RF/PBC[FFT] under a reaction field), is the continuum-electrostatics difference between the
//! macroscopic Coulombic system and the simulated one for the real solute shape: the polarisation
//! terms, B for a lattice sum and A1 + A2 + B for a reaction field. The analytic A1, A2 and B here
//! are the same quantities for a sphere and must not be added to it. What is left on top of
//! `dGslv_pbsolv` is C1 + C2 ([`ChargeCorrection::beyond_pb`]), which no continuum model sees:
//! it comes from the discrete solvent molecules.
//!
//! **The reaction-field self term comes off first.** GROMOS books a reaction-field self energy
//! E_self = −½ q_i² (4πε₀)⁻¹ (1 − C_rf/2)/R_C for every charged atom (gromosXX
//! `nonbonded_term.cc:79-99`, gromos-forces `rf_excluded_interactions` and
//! `perturbed_self_energy_correction`), and the perturbed loop books its λ-derivative, so a raw
//! GROMOS TI/BAR/MBAR ΔG of a charge change contains ΔG_self = −½ (4πε₀)⁻¹ (1 − C_rf/2)/R_C
//! Σ_i (q_B,i² − q_A,i²) exactly ([`rf_self_energy`]). The scheme here is built for a raw ΔG
//! without it (Reif, Hünenberger and Oostenbrink, J. Chem. Theory Comput. 8, 3705 (2012); Reif
//! and Oostenbrink, J. Comput. Chem. 35, 227 (2014)), so the corrected value is
//!
//! ΔG_chg^cor = [ΔG_raw − ΔG_self] + A1 + A2 + B + C1 + C2 ([`ChargeCorrection::corrected`])
//!
//! or, with `dGslv_pbsolv`, [ΔG_raw − ΔG_self] + PB(B) − PB(A) + C1 + C2
//! ([`ChargeCorrection::corrected_with_pb`]). ΔG_self is not part of A1: for q 0 → 1 at
//! R_C = 1.4 nm and EPSRF 61 it is −73.8 kJ/mol next to A1 = −48.8, and leaving it in shifts
//! the corrected value by that much.
//!
//! Units are gromos-core's: nm, e, kJ/mol, with 1/(4πε₀) = `four_pi_eps_i` (kJ mol⁻¹ nm e⁻²),
//! so ε₀⁻¹ = 4π·`four_pi_eps_i`.
//!
//! Not covered: the generalised reaction field (GM, GA; its A2 needs a three-dimensional solve,
//! Sec. II D), lattice sums with vacuum or adjusted boundaries (LV, LA) and finite droplets (FS).
//!
//! [^kh2006]: M. A. Kastenholz and P. H. Hünenberger, "Computation of methodology-independent
//!     ionic solvation free energies from molecular simulations. II. The hydration free energy
//!     of the sodium cation", J. Chem. Phys. 124, 224501 (2006), doi:10.1063/1.2201698. Equation
//!     and table numbers in this module are that paper's.

use gromos_forces::CRFParameters;
use nalgebra::{DMatrix, DVector};
use std::f64::consts::PI;
use std::str::FromStr;

/// ξ_LS, the Wigner self-term of a cubic lattice (Eq. 32).
pub const XI_LS: f64 = -2.837297;

/// Number density of water at 1 g cm⁻³ (nm⁻³), the paper's ρ (Sec. III).
pub const WATER_NUMBER_DENSITY: f64 = 33.427;

/// Solvent permittivity the paper uses for A1 and B: water at 298 K (Sec. III).
pub const WATER_PERMITTIVITY: f64 = 78.0;

/// Ionic radius the paper uses for every term (nm, Sec. III).
pub const PAPER_ION_RADIUS: f64 = 0.2;

/// The electrostatic scheme of the simulation (Sec. II C).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scheme {
    /// `LS`: lattice sum with tinfoil boundary conditions.
    LatticeSum,
    /// `CM`: Coulomb interactions truncated at a molecular cut-off.
    Cutoff,
    /// `BM`: Barker–Watts reaction field with a molecular cut-off; GROMOS's reaction field with
    /// charge groups (κ = 0).
    ReactionField,
    /// `BA`: Barker–Watts reaction field with an atomic cut-off.
    ReactionFieldAtomic,
}

impl FromStr for Scheme {
    type Err = String;

    /// `LS`, `CM`, `BM` (or `RF`) and `BA`, in any case.
    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_uppercase().as_str() {
            "LS" => Ok(Scheme::LatticeSum),
            "CM" => Ok(Scheme::Cutoff),
            "BM" | "RF" => Ok(Scheme::ReactionField),
            "BA" => Ok(Scheme::ReactionFieldAtomic),
            _ => Err(format!(
                "unknown electrostatics scheme '{s}' (LS, CM, BM or RF, BA)"
            )),
        }
    }
}

/// A rigid solvent model with one van der Waals site M, as the C terms need it (Sec. II F).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolventModel {
    pub name: &'static str,
    /// (charge in e, distance from M in nm) of every charged site.
    pub sites: &'static [(f64, f64)],
    /// Relative permittivity of the model, where a value is quoted with it.
    pub permittivity: Option<f64>,
}

impl SolventModel {
    /// The quadrupole-moment trace relative to M, γ = Σ_i q_i r_i² (e nm², Eq. 35).
    pub fn quadrupole_trace(&self) -> f64 {
        self.sites.iter().map(|&(q, r)| q * r * r).sum()
    }

    /// The model called `name` (case-insensitive; `SPC/E` also as `SPCE`).
    pub fn find(name: &str) -> Option<&'static SolventModel> {
        let key = name.to_ascii_uppercase().replace('/', "");
        SOLVENTS
            .iter()
            .copied()
            .find(|m| m.name.to_ascii_uppercase().replace('/', "") == key)
    }
}

/// SPC (Berendsen, Postma, van Gunsteren and Hermans, in *Intermolecular Forces*, Reidel 1981,
/// pp. 331–342): O −0.82 e, H +0.41 e at 0.1 nm; γ = 0.0082 e nm². The permittivity 66.6 is the
/// one the paper uses for SPC (its ref. 250: Glättli, Daura and van Gunsteren, J. Chem. Phys.
/// 116, 9811 (2002)). With ρ = 33.427 nm⁻³ it gives the paper's ΔΦ_ODL = 79.76 kJ mol⁻¹ e⁻¹
/// (C1 of Table V).
pub const SPC: SolventModel = SolventModel {
    name: "SPC",
    sites: &[(-0.82, 0.0), (0.41, 0.1), (0.41, 0.1)],
    permittivity: Some(66.6),
};

/// SPC/E (Berendsen, Grigera and Straatsma, J. Phys. Chem. 91, 6269 (1987)): H +0.4238 e at
/// 0.1 nm; γ = 0.008476 e nm².
pub const SPC_E: SolventModel = SolventModel {
    name: "SPC/E",
    sites: &[(-0.8476, 0.0), (0.4238, 0.1), (0.4238, 0.1)],
    permittivity: None,
};

/// TIP3P (Jorgensen, Chandrasekhar, Madura, Impey and Klein, J. Chem. Phys. 79, 926 (1983)):
/// H +0.417 e at 0.09572 nm; γ = 0.0076414 e nm².
pub const TIP3P: SolventModel = SolventModel {
    name: "TIP3P",
    sites: &[(-0.834, 0.0), (0.417, 0.09572), (0.417, 0.09572)],
    permittivity: None,
};

/// TIP4P (Jorgensen et al. 1983, as [`TIP3P`]): H +0.52 e at 0.09572 nm, the charge −1.04 e on
/// the virtual site 0.015 nm from O; γ = 0.0092948 e nm².
pub const TIP4P: SolventModel = SolventModel {
    name: "TIP4P",
    sites: &[(-1.04, 0.015), (0.52, 0.09572), (0.52, 0.09572)],
    permittivity: None,
};

/// The tabulated solvent models.
pub const SOLVENTS: [&SolventModel; 4] = [&SPC, &SPC_E, &TIP3P, &TIP4P];

/// What the correction is computed for.
#[derive(Debug, Clone, PartialEq)]
pub struct ChargeCorrectionInput {
    pub scheme: Scheme,
    /// Charge of the ion before and after the change (e).
    pub charge_a: f64,
    pub charge_b: f64,
    /// Edge L of the cubic box (nm).
    pub box_edge: f64,
    /// Cut-off R_C (nm); not read for [`Scheme::LatticeSum`].
    pub cutoff: f64,
    /// Reaction-field permittivity ε_BW (GROMOS EPSRF; 0 is a conducting boundary, as in
    /// GROMOS); read for [`Scheme::ReactionField`] and [`Scheme::ReactionFieldAtomic`] only.
    pub eps_rf: f64,
    /// Solvent permittivity ε_S of A1 and B (the paper: [`WATER_PERMITTIVITY`]).
    pub eps_solvent: f64,
    /// Permittivity of the solvent model, ε_S of A2 (the paper: 66.6 for SPC); `None` takes the
    /// model's.
    pub eps_model: Option<f64>,
    /// Ionic radius R_I (nm).
    pub ion_radius: f64,
    /// Solvent number density ρ (nm⁻³).
    pub density: f64,
    pub solvent: SolventModel,
    /// γ in place of the model's [`SolventModel::quadrupole_trace`] (e nm²).
    pub quadrupole_trace: Option<f64>,
    /// 1/(4πε₀) in kJ mol⁻¹ nm e⁻².
    pub four_pi_eps_i: f64,
}

impl ChargeCorrectionInput {
    /// Charging from 0 to `charge` with the paper's choices: R_I = 0.2 nm, ε_S = 78 for A1 and
    /// B, the model's permittivity for A2, ρ = 33.427 nm⁻³, gromos-core's 1/(4πε₀).
    pub fn new(
        scheme: Scheme,
        charge: f64,
        box_edge: f64,
        cutoff: f64,
        eps_rf: f64,
        solvent: SolventModel,
    ) -> Self {
        Self {
            scheme,
            charge_a: 0.0,
            charge_b: charge,
            box_edge,
            cutoff,
            eps_rf,
            eps_solvent: WATER_PERMITTIVITY,
            eps_model: None,
            ion_radius: PAPER_ION_RADIUS,
            density: WATER_NUMBER_DENSITY,
            solvent,
            quadrupole_trace: None,
            four_pi_eps_i: gromos_core::units::four_pi_eps_i,
        }
    }
}

/// The correction terms, kJ/mol, each to be added to the raw ΔG of the charge change.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChargeCorrection {
    /// ΔG_A1: the Born energy of the solvent beyond the cut-off (Eq. 25); 0 for LS.
    pub a1: f64,
    /// ΔG_A2: the polarisation error inside the cut-off (Eqs. 26–29); 0 for LS.
    pub a2: f64,
    /// ΔG_B: finite size and periodicity (Eq. 32 for LS, Eq. 33 otherwise).
    pub b: f64,
    /// ΔG_C1: the summation scheme, from the solvent's quadrupole-moment trace (Eqs. 37–42).
    pub c1: f64,
    /// ΔG_C2: the interface potential; 0 for cut-off schemes (Eq. 46), the SPC fit for LS
    /// (Eq. 47).
    pub c2: f64,
    /// ΔΦ_ODL = (6ε₀)⁻¹ ρ γ (kJ mol⁻¹ e⁻¹, Eq. 34).
    pub delta_phi_odl: f64,
    /// The quadrupole-moment trace γ used (e nm²).
    pub quadrupole_trace: f64,
}

impl ChargeCorrection {
    /// A1 + A2 + B + C1 + C2: the whole correction for a spherical ion, without `dGslv_pbsolv`.
    pub fn total(&self) -> f64 {
        self.a1 + self.a2 + self.b + self.c1 + self.c2
    }

    /// C1 + C2: what to add on top of `dGslv_pbsolv`'s correction, which already holds the
    /// polarisation terms A1 + A2 + B.
    pub fn beyond_pb(&self) -> f64 {
        self.c1 + self.c2
    }

    /// The corrected ΔG of the charge change for a spherical ion,
    /// [raw − ΔG_self] + A1 + A2 + B + C1 + C2, with `rf_self` from [`rf_self_energy`] (0 for a
    /// lattice sum).
    pub fn corrected(&self, raw: f64, rf_self: f64) -> f64 {
        raw - rf_self + self.total()
    }

    /// The corrected ΔG with `dGslv_pbsolv`: [raw − ΔG_self] + PB(B) − PB(A) + C1 + C2, where
    /// `pb_a` and `pb_b` are the corrections `dGslv_pbsolv` reports for the two end states.
    pub fn corrected_with_pb(&self, raw: f64, rf_self: f64, pb_a: f64, pb_b: f64) -> f64 {
        raw - rf_self + (pb_b - pb_a) + self.beyond_pb()
    }
}

/// The correction terms of [`ChargeCorrectionInput`] (module documentation).
pub fn charge_correction(input: &ChargeCorrectionInput) -> Result<ChargeCorrection, String> {
    let ChargeCorrectionInput {
        scheme,
        charge_a: qa,
        charge_b: qb,
        box_edge,
        cutoff,
        eps_rf,
        eps_solvent,
        ion_radius,
        density,
        four_pi_eps_i: fpepsi,
        ..
    } = *input;
    positive("box edge", box_edge)?;
    positive("ion radius", ion_radius)?;
    positive("solvent density", density)?;
    positive("1/(4 pi eps0)", fpepsi)?;
    at_least_one("solvent permittivity", eps_solvent)?;
    let gamma = input
        .quadrupole_trace
        .unwrap_or_else(|| input.solvent.quadrupole_trace());
    let phi = delta_phi_odl(density, gamma, fpepsi);
    let square = qb * qb - qa * qa;
    let linear = qb - qa;

    let (a1, a2, b, c1, c2) = if scheme == Scheme::LatticeSum {
        let f = liquid_fraction(box_edge, ion_radius)?;
        let c2 = if input.solvent.name == SPC.name && input.quadrupole_trace.is_none() {
            c2_lattice_sum_spc(1.0, f, ion_radius)
        } else {
            return Err(format!(
                "C2 for a lattice sum is the SPC fit of Eq. 47; there is none for {}",
                if input.quadrupole_trace.is_some() {
                    "a solvent given by its quadrupole-moment trace"
                } else {
                    input.solvent.name
                }
            ));
        };
        (
            0.0,
            0.0,
            square * b_lattice_sum(1.0, box_edge, ion_radius, eps_solvent, fpepsi),
            linear * c1(scheme, 1.0, phi, f, ion_radius, 0.0, 1.0),
            linear * c2,
        )
    } else {
        positive("cut-off", cutoff)?;
        if cutoff <= ion_radius {
            return Err(format!(
                "cut-off {cutoff} nm must exceed the ion radius {ion_radius} nm"
            ));
        }
        let eps_rf = match scheme {
            Scheme::Cutoff => 1.0,
            _ => {
                if eps_rf != 0.0 {
                    at_least_one("reaction-field permittivity", eps_rf)?;
                }
                eps_rf
            },
        };
        let eps_model = input
            .eps_model
            .or(input.solvent.permittivity)
            .ok_or_else(|| {
                format!(
                    "A2 needs the permittivity of {}: give eps_model",
                    input.solvent.name
                )
            })?;
        at_least_one("model permittivity", eps_model)?;
        let a1 = a1_cutoff(1.0, cutoff, eps_solvent, fpepsi);
        let a2 = a2_cutoff(1.0, cutoff, ion_radius, eps_model, eps_rf, fpepsi);
        let b = b_cutoff(
            scheme,
            1.0,
            box_edge,
            cutoff,
            ion_radius,
            eps_solvent,
            a1 + a2,
            fpepsi,
        );
        (
            square * a1,
            square * a2,
            square * b,
            linear * c1(scheme, 1.0, phi, 1.0, ion_radius, cutoff, eps_rf),
            0.0,
        )
    };
    Ok(ChargeCorrection {
        a1,
        a2,
        b,
        c1,
        c2,
        delta_phi_odl: phi,
        quadrupole_trace: gamma,
    })
}

fn positive(what: &str, x: f64) -> Result<(), String> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(format!("{what} must be positive, got {x}"))
    }
}

fn at_least_one(what: &str, x: f64) -> Result<(), String> {
    if x.is_finite() && x >= 1.0 {
        Ok(())
    } else {
        Err(format!("{what} must be at least 1, got {x}"))
    }
}

/// (8πε₀)⁻¹ q² (1 − ε⁻¹) / R: the Born free energy of a charge `q` in a cavity of radius R,
/// with the sign removed.
fn born(q: f64, radius: f64, eps: f64, fpepsi: f64) -> f64 {
    0.5 * fpepsi * q * q * (1.0 - 1.0 / eps) / radius
}

/// The Barker–Watts coefficient a₂ = (ε_BW − 1)/(2ε_BW + 1) of Eq. A4; ε_BW = 0 is conducting
/// (GROMOS's EPSRF convention), a₂ = 1/2.
fn bw_a2(eps_rf: f64) -> f64 {
    if eps_rf == 0.0 {
        0.5
    } else {
        (eps_rf - 1.0) / (2.0 * eps_rf + 1.0)
    }
}

/// ΔG_A1 = −(8πε₀)⁻¹ q² (1 − ε_S⁻¹) R_C⁻¹ (Eq. 25): the ion–solvent interactions beyond the
/// cut-off, as a Born continuum. Negative.
pub fn a1_cutoff(q: f64, cutoff: f64, eps_solvent: f64, fpepsi: f64) -> f64 {
    -born(q, cutoff, eps_solvent, fpepsi)
}

/// ΔG_A2 (Eq. 29): the error in the solvent polarisation inside the cut-off sphere for truncated
/// Coulomb (`eps_rf` = 1) or Barker–Watts reaction-field interactions, against the Born
/// polarisation.
///
/// The reduced radial polarisation p(r) = 4π q⁻¹ (1 − ε_S⁻¹)⁻¹ P(r) (Born: r⁻²) solves the
/// Fredholm equation of the second kind (Eq. 26)
///
/// p(r) = g(r) + ∫ K(r, r′) p(r′) dr′ over max(R_I, |r − R_C|) ≤ r′ ≤ r + R_C, r > R_I,
///
/// with (Eqs. A2, A3, A5, A6) g(r) = H(R_C − r)[r⁻² − 2a₂ R_C⁻³ r] and
/// K(r, r′) = (1 + a₂)(1 − ε_S⁻¹)/(8R_C³r²) [2R_C²(r² + r′²) − (r² − r′²)² − R_C⁴],
/// a₂ = (ε_BW − 1)/(2ε_BW + 1) (0 for truncated Coulomb). Then
///
/// ΔG_A2 = −(8πε₀)⁻¹ q² (1 − ε_S⁻¹) [R_I⁻¹ − R_C⁻¹ − ∫ p dr + 2a₂R_C⁻³ ∫ r³ p dr],
///
/// both integrals over R_I…R_C (Eq. 29 with a₋₁ = 1, a₁ = 0). Solved by product integration:
/// p piecewise linear on a grid with a node on each side of R_C (g, and so p, jumps there),
/// uniform up to 3R_C and growing by 5 % per element to 13R_C, beyond which p is taken as 0; the
/// kernel × hat-function integrals by three-point Gauss–Legendre on each element's part inside
/// the window (exact: the integrand is a polynomial of degree 5); the bracket at the spacings
/// R_C/80 and R_C/160, Richardson-extrapolated (the error is second order). The extrapolated
/// value moves by less than 0.003 kJ/mol when both spacings are halved, and by less than
/// 10⁻⁴ kJ/mol when the outer edge goes from 13R_C to 25R_C.
pub fn a2_cutoff(
    q: f64,
    cutoff: f64,
    ion_radius: f64,
    eps_model: f64,
    eps_rf: f64,
    fpepsi: f64,
) -> f64 {
    let a2 = bw_a2(eps_rf);
    let coarse = polarisation_bracket(cutoff, ion_radius, eps_model, a2, cutoff / 80.0);
    let fine = polarisation_bracket(cutoff, ion_radius, eps_model, a2, cutoff / 160.0);
    let bracket = (4.0 * fine - coarse) / 3.0;
    -0.5 * fpepsi * q * q * (1.0 - 1.0 / eps_model) * bracket
}

/// Three-point Gauss–Legendre nodes and weights on [−1, 1].
const GAUSS3: [(f64, f64); 3] = [
    (-0.774_596_669_241_483_4, 5.0 / 9.0),
    (0.0, 8.0 / 9.0),
    (0.774_596_669_241_483_4, 5.0 / 9.0),
];

/// The bracket of Eq. 29 on a grid of spacing `h` ([`a2_cutoff`]).
fn polarisation_bracket(rc: f64, ri: f64, eps: f64, a2: f64, h: f64) -> f64 {
    // Nodes: [R_I, R_C] (index 0..=n1, the last at R_C⁻), then [R_C⁺, ...].
    let n1 = ((rc - ri) / h).ceil() as usize;
    let mut nodes: Vec<f64> = (0..=n1)
        .map(|i| ri + (rc - ri) * i as f64 / n1 as f64)
        .collect();
    let n2 = (2.0 * rc / h).ceil() as usize;
    nodes.extend((0..=n2).map(|i| rc + 2.0 * rc * i as f64 / n2 as f64));
    let mut step = h;
    while *nodes.last().unwrap() < 13.0 * rc {
        step *= 1.05;
        let next = nodes.last().unwrap() + step;
        nodes.push(next);
    }
    let n = nodes.len();
    let elements: Vec<(usize, usize)> = (0..n1)
        .map(|e| (e, e + 1))
        .chain((n1 + 1..n - 1).map(|e| (e, e + 1)))
        .collect();

    let b = 2.0 * a2 / rc.powi(3);
    let g = DVector::from_fn(n, |i, _| {
        if i <= n1 {
            let r = nodes[i];
            1.0 / (r * r) - b * r
        } else {
            0.0
        }
    });
    let prefactor = (1.0 + a2) * (1.0 - 1.0 / eps) / (8.0 * rc.powi(3));
    let mut a = DMatrix::<f64>::identity(n, n);
    for i in 0..n {
        let r = nodes[i];
        let (lo, hi) = (ri.max((r - rc).abs()), r + rc);
        for &(ea, eb) in &elements {
            let (xa, xb) = (nodes[ea], nodes[eb]);
            let (l, u) = (lo.max(xa), hi.min(xb));
            if u <= l {
                continue;
            }
            let (mut wa, mut wb) = (0.0, 0.0);
            for (x, w) in GAUSS3 {
                let t = 0.5 * (u - l) * x + 0.5 * (u + l);
                let wt = 0.5 * (u - l) * w;
                let k = prefactor / (r * r)
                    * (2.0 * rc * rc * (r * r + t * t) - (r * r - t * t).powi(2) - rc.powi(4));
                wa += wt * k * (xb - t) / (xb - xa);
                wb += wt * k * (t - xa) / (xb - xa);
            }
            a[(i, ea)] -= wa;
            a[(i, eb)] -= wb;
        }
    }
    let p = a
        .lu()
        .solve(&g)
        .expect("the polarisation equation (I − K) p = g is regular");

    let (mut i0, mut i3) = (0.0, 0.0);
    for e in 0..n1 {
        let (xa, xb) = (nodes[e], nodes[e + 1]);
        for (x, w) in GAUSS3 {
            let t = 0.5 * (xb - xa) * x + 0.5 * (xb + xa);
            let wt = 0.5 * (xb - xa) * w;
            let pt = (p[e] * (xb - t) + p[e + 1] * (t - xa)) / (xb - xa);
            i0 += wt * pt;
            i3 += wt * t.powi(3) * pt;
        }
    }
    1.0 / ri - 1.0 / rc - i0 + b * i3
}

/// ΔG_B for a lattice sum with tinfoil boundaries (Eq. 32):
/// (8πε₀)⁻¹ q² (1 − ε_S⁻¹) L⁻¹ [ξ_LS + (4π/3)(R_I/L)² − (16π²/45)(R_I/L)⁵]. Negative: the
/// periodic copies of the ion under-polarise the solvent.
pub fn b_lattice_sum(q: f64, box_edge: f64, ion_radius: f64, eps_solvent: f64, fpepsi: f64) -> f64 {
    let x = ion_radius / box_edge;
    born(q, box_edge, eps_solvent, fpepsi)
        * (XI_LS + 4.0 * PI / 3.0 * x * x - 16.0 * PI * PI / 45.0 * x.powi(5))
}

/// ΔG_B for a cut-off scheme (Eq. 33), the paper's fit to continuum electrostatics for
/// ε_S = 78 over L/(2R_C) = 0…1.7 (its ref. 41):
/// −(R_I/R_C) 10^{μ L/(2R_C) + ν} [(8πε₀)⁻¹ q² (1 − ε_S⁻¹) R_I⁻¹ + ΔG_A^CT], with
/// ΔG_A^CT = ΔG_A1 + ΔG_A2 (`a_ct`, for the same charge) and μ = −2.20, ν = 1.29 (CM) or
/// μ = −1.47, ν = 0.72 (BM, BA). Not for LS.
#[allow(clippy::too_many_arguments)]
pub fn b_cutoff(
    scheme: Scheme,
    q: f64,
    box_edge: f64,
    cutoff: f64,
    ion_radius: f64,
    eps_solvent: f64,
    a_ct: f64,
    fpepsi: f64,
) -> f64 {
    let (mu, nu) = match scheme {
        Scheme::Cutoff => (-2.20, 1.29),
        Scheme::ReactionField | Scheme::ReactionFieldAtomic => (-1.47, 0.72),
        Scheme::LatticeSum => return 0.0,
    };
    -(ion_radius / cutoff)
        * 10f64.powf(mu * box_edge / (2.0 * cutoff) + nu)
        * (born(q, ion_radius, eps_solvent, fpepsi) + a_ct)
}

/// ΔΦ_ODL = (6ε₀)⁻¹ ρ γ (kJ mol⁻¹ e⁻¹, Eq. 34), ε₀⁻¹ = 4π·`fpepsi`.
pub fn delta_phi_odl(density: f64, quadrupole_trace: f64, fpepsi: f64) -> f64 {
    4.0 * PI * fpepsi / 6.0 * density * quadrupole_trace
}

/// f = 1 − 4πR_I³/(3L³), the fraction of the box the solvent fills (Eq. 38).
pub fn liquid_fraction(box_edge: f64, ion_radius: f64) -> Result<f64, String> {
    let f = 1.0 - 4.0 * PI * ion_radius.powi(3) / (3.0 * box_edge.powi(3));
    if f > 0.0 {
        Ok(f)
    } else {
        Err(format!(
            "an ion of radius {ion_radius} nm does not fit a box of edge {box_edge} nm"
        ))
    }
}

/// ΔG_C1 for a charge `q`, from ΔΦ_ODL = `phi`:
/// LS −q f ΔΦ_ODL (Eq. 37, `f` from [`liquid_fraction`]); CM 0 (Eq. 39);
/// BM −q [2(ε_BW − 1)/(2ε_BW + 1)] (1 − R_I³/R_C³) ΔΦ_ODL (Eq. 41);
/// BA that − q [3/(2ε_BW + 1)] ΔΦ_ODL (Eq. 42). ε_BW = 0 is conducting. Negative for a
/// cation in water (γ > 0).
pub fn c1(
    scheme: Scheme,
    q: f64,
    phi: f64,
    f: f64,
    ion_radius: f64,
    cutoff: f64,
    eps_rf: f64,
) -> f64 {
    match scheme {
        Scheme::LatticeSum => -q * f * phi,
        Scheme::Cutoff => 0.0,
        Scheme::ReactionField | Scheme::ReactionFieldAtomic => {
            let a2 = bw_a2(eps_rf);
            let bm = -q * 2.0 * a2 * (1.0 - (ion_radius / cutoff).powi(3)) * phi;
            if scheme == Scheme::ReactionField {
                bm
            } else {
                // 3/(2ε + 1) = 1 − 2a₂, which is 0 for a conducting boundary.
                bm - q * (1.0 - 2.0 * a2) * phi
            }
        },
    }
}

/// ΔG_C2 for a lattice sum, q (f − 1)(70.3 − 10.9 R_I⁻¹) kJ/mol with q in e and R_I in nm
/// (Eq. 47): the paper's fit for SPC at 300 K.
pub fn c2_lattice_sum_spc(q: f64, f: f64, ion_radius: f64) -> f64 {
    q * (f - 1.0) * (70.3 - 10.9 / ion_radius)
}

/// ΔG_self, the change of the GROMOS reaction-field self energy when the charges go from
/// `charges_a` to `charges_b` (e): −½ (4πε₀)⁻¹ (1 − C_rf/2)/R_C Σ_i (q_B,i² − q_A,i²) kJ/mol, with
/// C_rf of gromos-forces' [`CRFParameters`] for the reaction-field radius `cutoff` (RCRF, nm),
/// `eps_rf` (EPSRF; 0 is conducting), `kappa` (RCRF's κ, nm⁻¹) and `eps_cs` (EPSCS). A raw
/// GROMOS ΔG holds it exactly; subtract it before the correction terms (module documentation).
/// `cutoff` 0 is no reaction field (gromos-forces' convention) and gives 0.
pub fn rf_self_energy(
    charges_a: &[f64],
    charges_b: &[f64],
    cutoff: f64,
    eps_rf: f64,
    kappa: f64,
    eps_cs: f64,
    fpepsi: f64,
) -> Result<f64, String> {
    if charges_a.len() != charges_b.len() {
        return Err(format!(
            "{} charges in state A, {} in state B",
            charges_a.len(),
            charges_b.len()
        ));
    }
    if cutoff < 0.0 || eps_cs <= 0.0 || eps_rf < 0.0 || kappa < 0.0 {
        return Err(format!(
            "reaction field needs cutoff ≥ 0, eps_cs > 0, eps_rf ≥ 0 and kappa ≥ 0 \
             (got {cutoff}, {eps_cs}, {eps_rf}, {kappa})"
        ));
    }
    let crf_cut = CRFParameters::new(cutoff, eps_cs, eps_rf, kappa).crf_cut;
    let dq2: f64 = charges_a
        .iter()
        .zip(charges_b)
        .map(|(qa, qb)| qb * qb - qa * qa)
        .sum();
    Ok(-0.5 * fpepsi * crf_cut * dq2)
}

#[cfg(test)]
mod tests {
    //! Against Tables III and IV of the paper: Na⁺ (R_I = 0.2 nm) in SPC, ε_S = 78.0 for A1
    //! and B, 66.6 for A2, ε_BW = 66.6, ρ = 33.427 nm⁻³, L = [(N_w + 1)/ρ]^{1/3} (Sec. III; the
    //! tables print L rounded, and the LS row for N_w = 2048 prints 3.13 where the box is
    //! 3.94 nm). The tables give two decimals.

    use super::*;

    const F: f64 = gromos_core::units::four_pi_eps_i;
    const EPS_S: f64 = 78.0;
    const EPS_SPC: f64 = 66.6;
    const RI: f64 = 0.2;

    fn edge(nw: u32) -> f64 {
        ((nw as f64 + 1.0) / WATER_NUMBER_DENSITY).cbrt()
    }

    fn close(what: &str, ours: f64, paper: f64, tol: f64) {
        assert!(
            (ours - paper).abs() <= tol,
            "{what}: ours {ours:.4}, paper {paper:.2} (tolerance {tol})"
        );
    }

    /// Table III, the LS rows: N_w, ΔG_B, ΔG_C1, ΔG_C2.
    const LS: [(u32, f64, f64, f64); 10] = [
        (4, -293.17, -61.89, -3.54),
        (8, -259.71, -69.83, -1.97),
        (16, -221.47, -74.50, -1.04),
        (32, -183.85, -77.05, -0.54),
        (64, -150.00, -78.39, -0.27),
        (128, -121.08, -79.07, -0.14),
        (256, -97.09, -79.41, -0.07),
        (512, -77.55, -79.59, -0.03),
        (1024, -61.79, -79.67, -0.02),
        (2048, -49.16, -79.72, -0.01),
    ];

    #[test]
    fn lattice_sum_terms_match_table_iii() {
        let phi = delta_phi_odl(WATER_NUMBER_DENSITY, SPC.quadrupole_trace(), F);
        for (nw, b, c1_, c2) in LS {
            let l = edge(nw);
            let f = liquid_fraction(l, RI).unwrap();
            close(
                &format!("B LS N_w={nw}"),
                b_lattice_sum(1.0, l, RI, EPS_S, F),
                b,
                0.01,
            );
            close(
                &format!("C1 LS N_w={nw}"),
                c1(Scheme::LatticeSum, 1.0, phi, f, RI, 0.0, 0.0),
                c1_,
                0.01,
            );
            close(
                &format!("C2 LS N_w={nw}"),
                c2_lattice_sum_spc(1.0, f, RI),
                c2,
                0.01,
            );
        }
    }

    #[test]
    fn spc_quadrupole_trace_gives_the_papers_potential_jump() {
        // Table V lists C1 = −q ΔΦ_ODL = −79.76 kJ/mol for every GM/GA row.
        close("gamma SPC", SPC.quadrupole_trace(), 0.0082, 1e-12);
        close(
            "dPhi_ODL",
            delta_phi_odl(WATER_NUMBER_DENSITY, SPC.quadrupole_trace(), F),
            79.76,
            0.005,
        );
        close("gamma SPC/E", SPC_E.quadrupole_trace(), 0.008476, 1e-12);
        close("gamma TIP3P", TIP3P.quadrupole_trace(), 0.0076414, 1e-7);
        close("gamma TIP4P", TIP4P.quadrupole_trace(), 0.0092948, 1e-7);
        assert_eq!(SolventModel::find("spc/e"), Some(&SPC_E));
        assert_eq!(SolventModel::find("SPCE"), Some(&SPC_E));
        assert_eq!(SolventModel::find("tip4p"), Some(&TIP4P));
        assert!(SolventModel::find("TIP5P").is_none());
    }

    const CUTOFFS: [f64; 5] = [0.8, 1.0, 1.2, 1.4, 1.6];

    /// A1 depends on R_C alone; every CM, BM and BA row of Tables III and IV carries one of
    /// these five (Table IV's BA row N_w = 1024, R_C = 1.2 prints −557.15 for −57.15).
    #[test]
    fn a1_matches_tables_iii_and_iv() {
        for (rc, paper) in CUTOFFS.iter().zip([-85.72, -68.58, -57.15, -48.98, -42.86]) {
            close(
                &format!("A1 R_C={rc}"),
                a1_cutoff(1.0, *rc, EPS_S, F),
                paper,
                0.01,
            );
        }
    }

    /// Our converged A2 for CM and BM at the five cut-offs (Richardson-extrapolated solve).
    fn a2_table(eps_rf: f64) -> Vec<f64> {
        CUTOFFS
            .iter()
            .map(|&rc| a2_cutoff(1.0, rc, RI, EPS_SPC, eps_rf, F))
            .collect()
    }

    #[test]
    fn a2_reaction_field_matches_table_iv() {
        let paper = [-42.36, -34.14, -28.55, -24.54, -21.52];
        for ((rc, ours), paper) in CUTOFFS.iter().zip(a2_table(EPS_SPC)).zip(paper) {
            close(&format!("A2 BM R_C={rc}"), ours, paper, 0.1);
        }
        // A conducting boundary (GROMOS EPSRF 0) is the limit of a large ε_BW.
        close(
            "A2 BM conducting",
            a2_cutoff(1.0, 1.4, RI, EPS_SPC, 0.0, F),
            a2_cutoff(1.0, 1.4, RI, EPS_SPC, 1e9, F),
            1e-6,
        );
    }

    #[test]
    fn a2_cutoff_coulomb_matches_table_iii_to_a_quarter_kj() {
        // Not reproduced to 0.1 kJ/mol: our solve of Eqs. 26-29 (with the kernel of Eq. A3,
        // which is Eq. 28 at N_a = 0) converges 0.07-0.25 kJ/mol below Table III, the gap
        // shrinking as R_C grows, while the same solver gives Table IV's BM column within
        // 0.06. An unconverged solve at spacing R_C/80 alone lands within 0.12 of the table,
        // so the table's CM column is most likely its own discretisation error. Achieved
        // tolerance: 0.26 kJ/mol.
        let paper = [24.48, 20.58, 17.61, 15.33, 13.54];
        for ((rc, ours), paper) in CUTOFFS.iter().zip(a2_table(1.0)).zip(paper) {
            close(&format!("A2 CM R_C={rc}"), ours, paper, 0.26);
            assert!(ours < paper, "A2 CM R_C={rc} lies below the table");
        }
    }

    /// Tables III and IV, ΔG_B of the CM, BM and BA rows: (scheme, N_w, R_C, ΔG_B).
    const B_CUTOFF: [(Scheme, u32, f64, f64); 34] = [
        (Scheme::Cutoff, 512, 0.8, -0.53),
        (Scheme::Cutoff, 512, 1.0, -2.12),
        (Scheme::Cutoff, 512, 1.2, -5.20),
        (Scheme::Cutoff, 1024, 0.8, -0.07),
        (Scheme::Cutoff, 1024, 1.0, -0.41),
        (Scheme::Cutoff, 1024, 1.2, -1.33),
        (Scheme::Cutoff, 1024, 1.4, -2.99),
        (Scheme::Cutoff, 2048, 0.8, -0.01),
        (Scheme::Cutoff, 2048, 1.0, -0.05),
        (Scheme::Cutoff, 2048, 1.2, -0.24),
        (Scheme::Cutoff, 2048, 1.4, -0.69),
        (Scheme::Cutoff, 2048, 1.6, -1.49),
        (Scheme::ReactionField, 512, 0.8, -1.47),
        (Scheme::ReactionField, 512, 1.0, -3.76),
        (Scheme::ReactionField, 512, 1.2, -6.76),
        (Scheme::ReactionField, 1024, 0.8, -0.38),
        (Scheme::ReactionField, 1024, 1.0, -1.26),
        (Scheme::ReactionField, 1024, 1.2, -2.72),
        (Scheme::ReactionField, 1024, 1.4, -4.59),
        (Scheme::ReactionField, 2048, 0.8, -0.07),
        (Scheme::ReactionField, 2048, 1.0, -0.32),
        (Scheme::ReactionField, 2048, 1.2, -0.87),
        (Scheme::ReactionField, 2048, 1.4, -1.72),
        (Scheme::ReactionField, 2048, 1.6, -2.82),
        (Scheme::ReactionFieldAtomic, 512, 0.8, -1.47),
        (Scheme::ReactionFieldAtomic, 512, 1.0, -3.76),
        (Scheme::ReactionFieldAtomic, 1024, 0.8, -0.38),
        (Scheme::ReactionFieldAtomic, 1024, 1.0, -1.26),
        (Scheme::ReactionFieldAtomic, 1024, 1.2, -2.72),
        (Scheme::ReactionFieldAtomic, 2048, 0.8, -0.07),
        (Scheme::ReactionFieldAtomic, 2048, 1.0, -0.32),
        (Scheme::ReactionFieldAtomic, 2048, 1.2, -0.87),
        (Scheme::ReactionFieldAtomic, 2048, 1.4, -1.72),
        (Scheme::ReactionFieldAtomic, 2048, 1.6, -2.82),
    ];

    #[test]
    fn b_cutoff_matches_tables_iii_and_iv() {
        let (cm, bm) = (a2_table(1.0), a2_table(EPS_SPC));
        for (scheme, nw, rc, paper) in B_CUTOFF {
            let k = CUTOFFS.iter().position(|&c| c == rc).unwrap();
            let a2 = if scheme == Scheme::Cutoff {
                cm[k]
            } else {
                bm[k]
            };
            let a_ct = a1_cutoff(1.0, rc, EPS_S, F) + a2;
            let ours = b_cutoff(scheme, 1.0, edge(nw), rc, RI, EPS_S, a_ct, F);
            close(
                &format!("B {scheme:?} N_w={nw} R_C={rc}"),
                ours,
                paper,
                0.02,
            );
        }
    }

    #[test]
    fn c1_cutoff_schemes_match_tables_iii_and_iv() {
        let phi = delta_phi_odl(WATER_NUMBER_DENSITY, SPC.quadrupole_trace(), F);
        let bm = [-76.76, -77.35, -77.62, -77.75, -77.82];
        let ba = [-78.54, -79.13, -79.40, -79.53, -79.60];
        for (k, &rc) in CUTOFFS.iter().enumerate() {
            let c = |s| c1(s, 1.0, phi, 1.0, RI, rc, EPS_SPC);
            assert_eq!(c(Scheme::Cutoff), 0.0);
            close(
                &format!("C1 BM R_C={rc}"),
                c(Scheme::ReactionField),
                bm[k],
                0.01,
            );
            close(
                &format!("C1 BA R_C={rc}"),
                c(Scheme::ReactionFieldAtomic),
                ba[k],
                0.01,
            );
        }
    }

    /// The corrected hydration free energy of the tables, ΔG°_slv = ΔG_chg + the five terms +
    /// (ΔG_cav + ΔG_cmp = 13.62), from [`charge_correction`] on the paper's inputs: one row of
    /// each scheme. (scheme, N_w, R_C, ΔG_chg, ΔG°_slv, tolerance)
    #[test]
    fn corrected_free_energies_match_the_tables() {
        let rows = [
            // Table III; the CM row carries A2's 0.25 kJ/mol (see above).
            (Scheme::LatticeSum, 1024, 0.0, -340.49, -468.35, 0.02),
            (Scheme::LatticeSum, 32, 0.0, -209.85, -457.67, 0.02),
            (Scheme::Cutoff, 1024, 0.8, -424.47, -472.16, 0.27),
            (Scheme::Cutoff, 2048, 1.6, -447.19, -464.37, 0.1),
            // Table IV.
            (Scheme::ReactionField, 2048, 1.4, -329.10, -468.47, 0.1),
            (Scheme::ReactionFieldAtomic, 512, 0.8, -271.22, -465.68, 0.1),
        ];
        for (scheme, nw, rc, chg, slv, tol) in rows {
            let mut input = ChargeCorrectionInput::new(scheme, 1.0, edge(nw), rc, EPS_SPC, SPC);
            input.four_pi_eps_i = F;
            let c = charge_correction(&input).unwrap();
            close(
                &format!("dG_slv {scheme:?} N_w={nw} R_C={rc}"),
                chg + c.total() + 13.62,
                slv,
                tol,
            );
            close("beyond_pb", c.beyond_pb(), c.c1 + c.c2, 0.0);
        }
    }

    #[test]
    fn a_charge_change_scales_the_quadratic_and_linear_terms_apart() {
        let unit = ChargeCorrectionInput::new(Scheme::ReactionField, 1.0, 3.0, 1.4, 0.0, SPC);
        let half = ChargeCorrectionInput {
            charge_a: 0.5,
            ..unit.clone()
        };
        let (u, h) = (
            charge_correction(&unit).unwrap(),
            charge_correction(&half).unwrap(),
        );
        for (what, uu, hh) in [("a1", u.a1, h.a1), ("a2", u.a2, h.a2), ("b", u.b, h.b)] {
            close(what, hh, 0.75 * uu, 1e-9);
        }
        close("c1", h.c1, 0.5 * u.c1, 1e-12);
        assert_eq!(h.c2, 0.0);
        // A conducting reaction field: BM's C1 factor is 1, BA adds nothing.
        close(
            "c1 conducting",
            u.c1,
            -(1.0 - (0.2f64 / 1.4).powi(3)) * u.delta_phi_odl,
            1e-12,
        );
    }

    #[test]
    fn inputs_without_a_formula_are_refused() {
        let ls =
            |solvent| ChargeCorrectionInput::new(Scheme::LatticeSum, 1.0, 3.0, 0.0, 0.0, solvent);
        assert!(charge_correction(&ls(SPC)).is_ok());
        assert!(charge_correction(&ls(TIP3P))
            .unwrap_err()
            .contains("Eq. 47"));
        let rf = ChargeCorrectionInput::new(Scheme::ReactionField, 1.0, 3.0, 1.4, 61.0, TIP3P);
        assert!(charge_correction(&rf).unwrap_err().contains("eps_model"));
        assert!(charge_correction(&ChargeCorrectionInput {
            eps_model: Some(94.0),
            ..rf.clone()
        })
        .is_ok());
        for bad in [
            ChargeCorrectionInput {
                cutoff: 0.1,
                ..rf.clone()
            },
            ChargeCorrectionInput {
                eps_rf: 0.5,
                eps_model: Some(94.0),
                ..rf.clone()
            },
            ChargeCorrectionInput {
                box_edge: 0.2,
                ..ls(SPC)
            },
        ] {
            assert!(charge_correction(&bad).is_err(), "{bad:?}");
        }
        assert_eq!("rf".parse::<Scheme>(), Ok(Scheme::ReactionField));
        assert!("GM".parse::<Scheme>().is_err());
    }

    /// ΔG_self for q 0 → 1, R_C = 1.4 nm, κ = 0, ε_CS = 1 at three EPSRF (GROMOS 61, the
    /// paper's 66.6 for SPC, 78.5).
    #[test]
    fn rf_self_energy_of_a_unit_charge() {
        for (eps_rf, want) in [(61.0, -73.825), (66.6, -73.875), (78.5, -73.959)] {
            let got = rf_self_energy(&[0.0], &[1.0], 1.4, eps_rf, 0.0, 1.0, F).unwrap();
            close(&format!("ΔG_self at EPSRF {eps_rf}"), got, want, 0.01);
        }
        // q_B² − q_A², summed over atoms; conducting boundary C_rf = −1
        let two = rf_self_energy(&[0.5, -1.0], &[1.0, 0.0], 1.4, 61.0, 0.0, 1.0, F).unwrap();
        let one = rf_self_energy(&[0.0], &[1.0], 1.4, 61.0, 0.0, 1.0, F).unwrap();
        assert!((two - one * (1.0 + 0.0 - 0.25 - 1.0)).abs() < 1e-12);
        let conducting = rf_self_energy(&[0.0], &[1.0], 1.4, 0.0, 0.0, 1.0, F).unwrap();
        assert!((conducting + 0.5 * F * 1.5 / 1.4).abs() < 1e-12);
        assert_eq!(
            rf_self_energy(&[1.0], &[0.0], 0.0, 61.0, 0.0, 1.0, F),
            Ok(0.0)
        );
        assert!(rf_self_energy(&[1.0], &[0.0, 0.0], 1.4, 61.0, 0.0, 1.0, F).is_err());
        assert!(rf_self_energy(&[1.0], &[0.0], 1.4, 61.0, 0.0, 0.0, F).is_err());
    }

    /// ΔG_self is what TI over the engine's own perturbed self term gives: the integral of
    /// gromos-forces' `perturbed_self_energy_correction` dH/dλ on a one-atom charge change,
    /// trapezoid over λ (exact for n = 1, a fine grid for n = 2), and the energy difference of
    /// its end states.
    #[test]
    fn rf_self_energy_is_the_ti_integral_of_the_perturbed_self_term() {
        use gromos_forces::nonbonded::{
            perturbed_self_energy_correction, EnergyGroups, GroupEnergies, PertAtomInfo,
            PerturbedLambdaParams,
        };
        let (qa, qb, rc, eps_rf, kappa) = (1.0, -0.3, 1.4, 61.0, 0.7);
        let crf = CRFParameters::new(rc, 1.0, eps_rf, kappa);
        let pert = [Some(PertAtomInfo {
            a_iac: 0,
            b_iac: 0,
            a_charge: qa,
            b_charge: qb,
            alpha_lj: 0.5,
            alpha_crf: 0.5,
        })];
        let at = |lam: f64, n: i32| {
            let lp = PerturbedLambdaParams::from_lambda(lam, lam, lam, lam, 1.0, 1.0, 1.0, 1.0, n);
            let mut sink = GroupEnergies::new(1);
            perturbed_self_energy_correction(
                &pert,
                &crf,
                &[lp],
                F,
                EnergyGroups::single(),
                &mut sink,
                None,
            )
        };
        let want = rf_self_energy(&[qa], &[qb], rc, eps_rf, kappa, 1.0, F).unwrap();
        for (n, points, tol) in [(1, 11, 1e-12), (2, 1001, 1e-5)] {
            let h = 1.0 / (points - 1) as f64;
            let dhdl: Vec<f64> = (0..points).map(|i| at(i as f64 * h, n).1).collect();
            let ti = h * (dhdl.iter().sum::<f64>() - 0.5 * (dhdl[0] + dhdl[points - 1]));
            assert!((ti - want).abs() < tol, "n = {n}: TI {ti}, ΔG_self {want}");
            let de = at(1.0, n).0 - at(0.0, n).0;
            assert!(
                (de - want).abs() < 1e-12,
                "n = {n}: ΔE {de}, ΔG_self {want}"
            );
        }
    }

    #[test]
    fn corrected_dg_takes_the_self_term_off_the_raw_value() {
        let c = ChargeCorrection {
            a1: -48.0,
            a2: 10.0,
            b: -5.0,
            c1: -77.0,
            c2: 0.0,
            delta_phi_odl: 79.76,
            quadrupole_trace: 0.0082,
        };
        assert!((c.corrected(-400.0, -73.8) - (-326.2 - 120.0)).abs() < 1e-12);
        assert!(
            (c.corrected_with_pb(-400.0, -73.8, 0.0, -77.7) - (-326.2 - 77.7 - 77.0)).abs() < 1e-12
        );
    }
}
