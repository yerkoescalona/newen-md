//! Individual λ per interaction kind and energy-group pair (gromosXX `LAMBDAS`, item 1.8).
//!
//! gromosXX links every perturbed interaction to the global λ through a polynomial
//! `L = a·λ⁴ + b·λ³ + c·λ² + d·λ + e` chosen per interaction kind and energy-group pair
//! (`simulation/parameter.h:740-800`, `topology/topology.cc:1419-1478`). The defaults
//! `a = b = c = e = 0`, `d = 1` give `L = λ` and `L′ = 1` exactly. One table, one evaluator:
//! the engine and the analysis programs read λ from [`IndividualLambdas`].

/// The interaction kinds that carry their own λ, in gromosXX's `interaction_lambda_enum` order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InteractionLambda {
    /// Bond stretching (and perturbed SHAKE constraint lengths).
    Bond,
    /// Bond angles.
    Angle,
    /// Proper dihedrals.
    Dihedral,
    /// Improper dihedrals.
    Improper,
    /// Lennard-Jones (`vdw`).
    Lj,
    /// Lennard-Jones soft-core (`vdw_soft`).
    LjSoftness,
    /// Coulomb reaction field (`crf`).
    Crf,
    /// Coulomb reaction field soft-core (`crf_soft`).
    CrfSoftness,
    /// Distance restraints (`distanceres`).
    Disres,
    /// Distance-field restraints (`distancefield`).
    Disfield,
    /// Dihedral restraints (`dihedralres`).
    Dihres,
    /// Perturbed masses.
    Mass,
    /// Angle restraints (no `LAMBDAS` row can set it: gromosXX rejects `angleres`).
    Angres,
    /// Collective-variable restraints (no `LAMBDAS` row can set it either).
    Colvarres,
}

impl InteractionLambda {
    /// Number of kinds (gromosXX `last_interaction_lambda`).
    pub const COUNT: usize = 14;

    /// Every kind, in table order.
    pub const ALL: [InteractionLambda; Self::COUNT] = [
        Self::Bond,
        Self::Angle,
        Self::Dihedral,
        Self::Improper,
        Self::Lj,
        Self::LjSoftness,
        Self::Crf,
        Self::CrfSoftness,
        Self::Disres,
        Self::Disfield,
        Self::Dihres,
        Self::Mass,
        Self::Angres,
        Self::Colvarres,
    ];

    /// Position in the table (the gromosXX enum value).
    #[inline(always)]
    pub fn index(self) -> usize {
        self as usize
    }

    /// The kind a numeric `LAMBDAS` NTLI selects: `k` in 1..=12 is kind `k − 1`
    /// (`in_parameter.cc:4419-4442`); anything else, 13 (`angleres`) included, is `None`.
    pub fn from_index(k: usize) -> Option<Self> {
        (1..=12).contains(&k).then(|| Self::ALL[k - 1])
    }

    /// The kind a named `LAMBDAS` NTLI selects (`in_parameter.cc:4390`); `None` for any other
    /// name, `angleres` included.
    pub fn from_gromos_name(name: &str) -> Option<Self> {
        Some(match name {
            "bond" => Self::Bond,
            "angle" => Self::Angle,
            "dihedral" => Self::Dihedral,
            "improper" => Self::Improper,
            "vdw" => Self::Lj,
            "vdw_soft" => Self::LjSoftness,
            "crf" => Self::Crf,
            "crf_soft" => Self::CrfSoftness,
            "distanceres" => Self::Disres,
            "distancefield" => Self::Disfield,
            "dihedralres" => Self::Dihres,
            "mass" => Self::Mass,
            _ => return None,
        })
    }

    /// An NTLI token as gromosXX accepts it: one of the twelve names or the numbers 1-12.
    pub fn from_ntli(token: &str) -> Option<Self> {
        Self::from_gromos_name(token).or_else(|| token.parse().ok().and_then(Self::from_index))
    }

    /// The `LAMBDAS` name of this kind (`angleres`/`colvarres` for the two kinds no row sets).
    pub fn gromos_name(self) -> &'static str {
        match self {
            Self::Bond => "bond",
            Self::Angle => "angle",
            Self::Dihedral => "dihedral",
            Self::Improper => "improper",
            Self::Lj => "vdw",
            Self::LjSoftness => "vdw_soft",
            Self::Crf => "crf",
            Self::CrfSoftness => "crf_soft",
            Self::Disres => "distanceres",
            Self::Disfield => "distancefield",
            Self::Dihres => "dihedralres",
            Self::Mass => "mass",
            Self::Angres => "angleres",
            Self::Colvarres => "colvarres",
        }
    }

    /// The four nonbonded kinds, the only ones gromosXX reads off the diagonal of the table.
    pub fn is_nonbonded(self) -> bool {
        matches!(
            self,
            Self::Lj | Self::LjSoftness | Self::Crf | Self::CrfSoftness
        )
    }
}

/// The identity polynomial: `L = λ`, `L′ = 1` (`in_parameter.cc:1013-1022`).
pub const IDENTITY_COEFFICIENTS: [f64; 5] = [0.0, 0.0, 0.0, 1.0, 0.0];

/// The individual-λ table: coefficients `[a, b, c, d, e]` and the evaluated `L`, `L′` for every
/// interaction kind and energy-group pair, flat `[kind][g1][g2]`.
#[derive(Debug, Clone, PartialEq)]
pub struct IndividualLambdas {
    n_groups: usize,
    coefficients: Vec<[f64; 5]>,
    lambda: Vec<f64>,
    derivative: Vec<f64>,
}

impl Default for IndividualLambdas {
    fn default() -> Self {
        Self::identity(1)
    }
}

impl IndividualLambdas {
    /// Every kind and pair at the identity polynomial, evaluated at λ = 0.
    pub fn identity(n_groups: usize) -> Self {
        let n_groups = n_groups.max(1);
        let cells = InteractionLambda::COUNT * n_groups * n_groups;
        let mut t = Self {
            n_groups,
            coefficients: vec![IDENTITY_COEFFICIENTS; cells],
            lambda: vec![0.0; cells],
            derivative: vec![0.0; cells],
        };
        t.update(0.0);
        t
    }

    /// Number of energy groups the table covers.
    pub fn n_groups(&self) -> usize {
        self.n_groups
    }

    #[inline(always)]
    fn cell(&self, kind: InteractionLambda, g1: usize, g2: usize) -> usize {
        (kind.index() * self.n_groups + g1) * self.n_groups + g2
    }

    /// Set the polynomial of `kind` for the pair `(g1, g2)` and `(g2, g1)` (0-based groups);
    /// a later call for the same pair overwrites. Call [`Self::update`] afterwards.
    pub fn set_row(&mut self, kind: InteractionLambda, g1: usize, g2: usize, coeff: [f64; 5]) {
        let (a, b) = (self.cell(kind, g1, g2), self.cell(kind, g2, g1));
        self.coefficients[a] = coeff;
        self.coefficients[b] = coeff;
    }

    /// The coefficients `[a, b, c, d, e]` of `kind` for `(g1, g2)`.
    pub fn coefficients(&self, kind: InteractionLambda, g1: usize, g2: usize) -> [f64; 5] {
        self.coefficients[self.cell(kind, g1, g2)]
    }

    /// Grow (or shrink) to `n_groups`, keeping the rows of the groups both sizes share and
    /// giving the new ones the identity (gromosXX pads the implicit trailing energy group so,
    /// `in_topology.cc:461-490`). Call [`Self::update`] afterwards.
    pub fn resize(&mut self, n_groups: usize) {
        let n_groups = n_groups.max(1);
        if n_groups == self.n_groups {
            return;
        }
        let mut t = Self::identity(n_groups);
        let shared = n_groups.min(self.n_groups);
        for kind in InteractionLambda::ALL {
            for g1 in 0..shared {
                for g2 in 0..shared {
                    let c = t.cell(kind, g1, g2);
                    t.coefficients[c] = self.coefficients(kind, g1, g2);
                }
            }
        }
        *self = t;
    }

    /// Evaluate every polynomial and its derivative at the global `lambda`, in gromosXX's
    /// expression and order (`Topology::update_for_lambda`, `topology.cc:1442-1455`): no Horner
    /// form, no clamping, the full square visited with both `[n1][n2]` and `[n2][n1]` written
    /// from the coefficients of `(n1, n2)`.
    pub fn update(&mut self, lambda: f64) {
        let lam = lambda;
        let n = self.n_groups;
        for kind in InteractionLambda::ALL {
            for n1 in 0..n {
                for n2 in 0..n {
                    let [a, b, c, d, e] = self.coefficients(kind, n1, n2);
                    let l = a * lam * lam * lam * lam
                        + b * lam * lam * lam
                        + c * lam * lam
                        + d * lam
                        + e;
                    let dl = 4.0 * a * lam * lam * lam + 3.0 * b * lam * lam + 2.0 * c * lam + d;
                    let (p, q) = (self.cell(kind, n1, n2), self.cell(kind, n2, n1));
                    self.lambda[p] = l;
                    self.derivative[p] = dl;
                    self.lambda[q] = l;
                    self.derivative[q] = dl;
                }
            }
        }
    }

    /// `L` of `kind` for the pair `(g1, g2)`.
    #[inline(always)]
    pub fn lambda(&self, kind: InteractionLambda, g1: usize, g2: usize) -> f64 {
        self.lambda[self.cell(kind, g1, g2)]
    }

    /// `L′ = dL/dλ` of `kind` for the pair `(g1, g2)`.
    #[inline(always)]
    pub fn derivative(&self, kind: InteractionLambda, g1: usize, g2: usize) -> f64 {
        self.derivative[self.cell(kind, g1, g2)]
    }

    /// `(L, L′)` of `kind` for the pair `(g1, g2)`.
    #[inline(always)]
    pub fn get(&self, kind: InteractionLambda, g1: usize, g2: usize) -> (f64, f64) {
        let c = self.cell(kind, g1, g2);
        (self.lambda[c], self.derivative[c])
    }
}

/// gromosXX `PRECALCLAM` (extended TI, item 1.8): the energies of every perturbed interaction
/// and their λ derivatives on the grid `λ_k = k·(max − min)/(n − 1) + min`, `k = 0 … n − 1`,
/// evaluated on the steps that write the free-energy trajectory (`step % every == 0`, NTWG).
/// The grid is the raw λ_k: no NLAM exponent, no individual-λ polynomial.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrecalcLambdas {
    /// NRLAM (≥ 2).
    pub nr_lambdas: usize,
    /// MINLAM
    pub min_lambda: f64,
    /// MAXLAM
    pub max_lambda: f64,
    /// NTWG: the energies are computed on the steps this divides.
    pub every: usize,
}

impl PrecalcLambdas {
    /// λ_k, as gromosXX computes it at every site (`(lam_index * lambda_step) + min_lam`).
    #[inline(always)]
    pub fn lambda(&self, k: usize) -> f64 {
        let step = (self.max_lambda - self.min_lambda) / (self.nr_lambdas - 1) as f64;
        (k as f64 * step) + self.min_lambda
    }

    /// Whether `step` is one the energies are computed on
    /// (`sim.steps() % sim.param().write.free_energy == 0`).
    #[inline(always)]
    pub fn active(&self, step: usize) -> bool {
        self.every > 0 && step.is_multiple_of(self.every)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_gives_lambda_and_a_unit_derivative_exactly() {
        let mut t = IndividualLambdas::identity(3);
        for lam in [0.0, 0.125, 0.3, 1.0 / 3.0, 0.9, 1.0] {
            t.update(lam);
            for kind in InteractionLambda::ALL {
                for g1 in 0..3 {
                    for g2 in 0..3 {
                        assert_eq!(t.lambda(kind, g1, g2).to_bits(), lam.to_bits());
                        assert_eq!(t.derivative(kind, g1, g2), 1.0);
                    }
                }
            }
        }
    }

    /// The individual λ of `X/check/lambdas.t.cc` at RLAM = 0.125 (group 1), as gromosXX's
    /// expression evaluates them.
    #[test]
    fn polynomials_at_one_eighth_match_the_check_suite() {
        let mut t = IndividualLambdas::identity(2);
        t.set_row(InteractionLambda::Bond, 0, 0, [0.0, 0.0, 1.0, 0.0, 0.0]);
        t.set_row(InteractionLambda::Angle, 0, 0, [1.0, 0.0, 0.0, 0.0, 0.0]);
        t.set_row(InteractionLambda::Dihedral, 0, 0, [0.5, 0.0, 0.5, 0.0, 0.0]);
        t.set_row(
            InteractionLambda::Improper,
            0,
            0,
            [0.33, 0.33, 0.34, 0.0, 0.0],
        );
        t.set_row(InteractionLambda::Mass, 0, 0, [0.0, 1.0, -1.0, 1.0, 0.0]);
        t.update(0.125);
        let at = |k| t.get(k, 0, 0);
        assert_eq!(at(InteractionLambda::Bond), (0.015625, 0.25));
        assert_eq!(at(InteractionLambda::Angle), (0.000244140625, 0.0078125));
        assert_eq!(
            at(InteractionLambda::Dihedral),
            (0.0079345703125, 0.12890625)
        );
        assert_eq!(
            at(InteractionLambda::Improper),
            (0.00603759765625, 0.10304687500000001)
        );
        assert_eq!(at(InteractionLambda::Mass), (0.111328125, 0.796875));
        // group 2 untouched
        assert_eq!(t.get(InteractionLambda::Bond, 1, 1), (0.125, 1.0));
    }

    #[test]
    fn rows_are_symmetric_and_the_last_one_wins() {
        let mut t = IndividualLambdas::identity(2);
        t.set_row(InteractionLambda::Lj, 0, 1, [0.0, 0.0, 0.1, 0.9, 0.0]);
        t.set_row(InteractionLambda::Lj, 0, 1, [0.0, 0.0, 0.0, 0.5, 0.0]);
        t.update(0.5);
        assert_eq!(t.get(InteractionLambda::Lj, 0, 1), (0.25, 0.5));
        assert_eq!(t.get(InteractionLambda::Lj, 1, 0), (0.25, 0.5));
        assert_eq!(t.get(InteractionLambda::Lj, 0, 0), (0.5, 1.0));
    }

    #[test]
    fn resize_pads_new_groups_with_the_identity() {
        let mut t = IndividualLambdas::identity(1);
        t.set_row(InteractionLambda::Crf, 0, 0, [0.0, 0.0, 1.0, 0.0, 0.0]);
        t.resize(2);
        t.update(0.5);
        assert_eq!(t.get(InteractionLambda::Crf, 0, 0), (0.25, 1.0));
        assert_eq!(t.get(InteractionLambda::Crf, 0, 1), (0.5, 1.0));
        assert_eq!(t.get(InteractionLambda::Crf, 1, 1), (0.5, 1.0));
    }

    #[test]
    fn ntli_accepts_the_twelve_names_and_numbers_only() {
        assert_eq!(
            InteractionLambda::from_ntli("vdw"),
            Some(InteractionLambda::Lj)
        );
        assert_eq!(
            InteractionLambda::from_ntli("5"),
            Some(InteractionLambda::Lj)
        );
        assert_eq!(
            InteractionLambda::from_ntli("12"),
            Some(InteractionLambda::Mass)
        );
        assert_eq!(InteractionLambda::from_ntli("13"), None);
        assert_eq!(InteractionLambda::from_ntli("0"), None);
        assert_eq!(InteractionLambda::from_ntli("angleres"), None);
        for k in 1..=12 {
            let kind = InteractionLambda::from_index(k).unwrap();
            assert_eq!(kind.index(), k - 1);
            assert_eq!(
                InteractionLambda::from_gromos_name(kind.gromos_name()),
                Some(kind)
            );
        }
    }
}
