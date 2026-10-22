//! LJ/CRF parameter types and force storage

use newen_core::math::Vec3;
use newen_logic::topology::groups::energy_cell;

/// Lennard-Jones parameters for atom pair
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct LJParameters {
    pub c6: f64,   // LJ C6 coefficient (attraction)
    pub c12: f64,  // LJ C12 coefficient (repulsion)
    pub cs6: f64,  // 1-4 C6 coefficient
    pub cs12: f64, // 1-4 C12 coefficient
}

impl From<&newen_core::topology::LJParameters> for LJParameters {
    fn from(p: &newen_core::topology::LJParameters) -> Self {
        Self {
            c6: p.c6,
            c12: p.c12,
            cs6: p.cs6,
            cs12: p.cs12,
        }
    }
}

/// Flat LJ parameter matrix for cache-friendly access.
///
/// Stores parameters in a contiguous array indexed as `[type_i * n_types + type_j]`.
/// Eliminates the double indirection of `Vec<Vec<LJParameters>>`.
#[derive(Debug, Clone)]
pub struct LJParamMatrix {
    data: Vec<LJParameters>,
    n_types: usize,
}

impl LJParamMatrix {
    /// Create from a nested Vec (converts from legacy format).
    pub fn from_nested(nested: &[Vec<LJParameters>]) -> Self {
        let n_types = nested.len();
        let mut data = Vec::with_capacity(n_types * n_types);
        for row in nested {
            assert_eq!(row.len(), n_types, "LJ parameter matrix must be square");
            data.extend_from_slice(row);
        }
        Self { data, n_types }
    }

    /// Get parameters for atom type pair (i, j).
    #[inline(always)]
    pub fn get(&self, type_i: usize, type_j: usize) -> LJParameters {
        unsafe { *self.data.get_unchecked(type_i * self.n_types + type_j) }
    }

    /// Number of atom types.
    pub fn n_types(&self) -> usize {
        self.n_types
    }
}

/// Coulomb Reaction Field parameters
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CRFParameters {
    pub crf_cut: f64,    // (1 - crf/2) / cutoff  (energy shift constant)
    pub crf_2cut3i: f64, // crf / (2 * cutoff^3)  (energy quadratic term)
    pub crf_cut3i: f64,  // crf / cutoff^3         (force correction term)
    pub cutoff_sq: f64,  // cutoff^2 for Heaviside truncation (GROMOS XXHEAVISIDE)
}

impl CRFParameters {
    /// Calculate CRF parameters from physical constants, as gromosXX's
    /// `Nonbonded_Term::init` and `Eds_Nonbonded_Term::init` do (`nonbonded_term.cc:51-121`,
    /// `perturbed_nonbonded_term.cc:51-111`).
    ///
    /// # Arguments
    /// * `cutoff` - Reaction-field radius RCRF in nm (gromosXX `nonbonded.rf_cutoff`; not RCUTL)
    /// * `epsilon` - Dielectric permittivity inside the cutoff
    /// * `rf_epsilon` - Dielectric permittivity outside the cutoff (reaction field); 0 means a
    ///   conducting boundary, `C_RF = −1`
    /// * `rf_kappa` - Debye-Hückel screening parameter in nm⁻¹
    ///
    /// RCRF = 0 switches the reaction field off: every constant is zero (plain Coulomb).
    ///
    /// spec: nonbonded.crf_parameters
    pub fn new(cutoff: f64, epsilon: f64, rf_epsilon: f64, rf_kappa: f64) -> Self {
        if cutoff <= 0.0 {
            return Self {
                crf_cut: 0.0,
                crf_2cut3i: 0.0,
                crf_cut3i: 0.0,
                cutoff_sq: 0.0,
            };
        }
        let kappa_cut = rf_kappa * cutoff;
        let kappa_cut2 = kappa_cut * kappa_cut;

        let crf = if rf_epsilon == 0.0 {
            -1.0
        } else {
            (2.0 * (epsilon - rf_epsilon) * (1.0 + kappa_cut) - rf_epsilon * kappa_cut2)
                / ((epsilon + 2.0 * rf_epsilon) * (1.0 + kappa_cut) + rf_epsilon * kappa_cut2)
        };

        let cut3i = 1.0 / (cutoff * cutoff * cutoff);
        let crf_cut3i = crf * cut3i;
        let crf_2cut3i = crf_cut3i / 2.0;
        let crf_cut = (1.0 - crf / 2.0) / cutoff;

        Self {
            crf_cut,
            crf_2cut3i,
            crf_cut3i,
            cutoff_sq: cutoff * cutoff,
        }
    }
}

/// Which energy group each atom belongs to (gromosXX `topo.atom_energy_group()`) and how many
/// groups there are. With one group nothing is booked per group, and the kernels run their
/// ungrouped instantiation.
#[derive(Debug, Clone, Copy)]
pub struct EnergyGroups<'a> {
    atom_group: &'a [usize],
    n_groups: usize,
}

impl<'a> EnergyGroups<'a> {
    /// One energy group: no per-group booking.
    pub const fn single() -> EnergyGroups<'static> {
        EnergyGroups {
            atom_group: &[],
            n_groups: 1,
        }
    }

    /// `atom_group[i]` is the energy group of atom `i`; it must cover every atom a kernel sees
    /// when `n_groups > 1`.
    pub fn new(atom_group: &'a [usize], n_groups: usize) -> Self {
        Self {
            atom_group,
            n_groups,
        }
    }

    /// More than one group: energies are booked per group pair.
    #[inline(always)]
    pub fn is_grouped(&self) -> bool {
        self.n_groups > 1
    }

    /// Number of energy groups.
    pub fn n_groups(&self) -> usize {
        self.n_groups
    }

    /// Row-major cell `[eg(i)][eg(j)]` of an `n_groups × n_groups` matrix.
    #[inline(always)]
    pub fn cell(&self, i: usize, j: usize) -> usize {
        energy_cell(self.atom_group[i], self.atom_group[j], self.n_groups)
            .expect("energy groups cover every atom")
    }
}

/// Pair energies and their λ derivatives per energy-group pair, flat row-major
/// `n_groups × n_groups` (gromosXX `lj_energy[eg(i)][eg(j)]`, `crf_energy`, and the same two of
/// `perturbed_energy_derivatives`). Both triangles are booked; the writer folds them.
/// Empty with one energy group.
#[derive(Debug, Clone, Default)]
pub struct GroupEnergies {
    pub n_groups: usize,
    pub lj: Vec<f64>,
    pub crf: Vec<f64>,
    pub dhdl_lj: Vec<f64>,
    pub dhdl_crf: Vec<f64>,
}

impl GroupEnergies {
    /// Zeroed matrices for `n_groups` groups; empty for one group.
    pub fn new(n_groups: usize) -> Self {
        let cells = if n_groups > 1 { n_groups * n_groups } else { 0 };
        Self {
            n_groups: if n_groups > 1 { n_groups } else { 1 },
            lj: vec![0.0; cells],
            crf: vec![0.0; cells],
            dhdl_lj: vec![0.0; cells],
            dhdl_crf: vec![0.0; cells],
        }
    }

    /// Whether the matrices are sized for more than one group.
    #[inline(always)]
    pub fn is_grouped(&self) -> bool {
        !self.lj.is_empty()
    }

    pub fn clear(&mut self) {
        self.lj.fill(0.0);
        self.crf.fill(0.0);
        self.dhdl_lj.fill(0.0);
        self.dhdl_crf.fill(0.0);
    }

    /// Element-wise sum of `other` into `self` (same group count, or `other` ungrouped).
    pub fn merge(&mut self, other: &GroupEnergies) {
        if !other.is_grouped() {
            return;
        }
        assert_eq!(self.lj.len(), other.lj.len(), "energy-group count mismatch");
        for (dst, src) in [
            (&mut self.lj, &other.lj),
            (&mut self.crf, &other.crf),
            (&mut self.dhdl_lj, &other.dhdl_lj),
            (&mut self.dhdl_crf, &other.dhdl_crf),
        ] {
            for (d, s) in dst.iter_mut().zip(src) {
                *d += *s;
            }
        }
    }

    /// Add a pair's energies and λ derivatives to flat cell `k`.
    #[inline(always)]
    pub fn book(&mut self, k: usize, e_lj: f64, e_crf: f64, de_lj: f64, de_crf: f64) {
        self.lj[k] += e_lj;
        self.crf[k] += e_crf;
        self.dhdl_lj[k] += de_lj;
        self.dhdl_crf[k] += de_crf;
    }
}

/// Extended-TI (PRECALCLAM, item 1.8) energies of the perturbed nonbonded terms at every λ_k
/// of the grid: A- and B-state Lennard-Jones and reaction field, energies and derivatives, flat
/// `[k][cell]` with `cell` the row-major energy-group pair (gromosXX `A_lj_energy[k][g1][g2]` and
/// its three siblings, on `energies` and `perturbed_energy_derivatives`). Empty — and never
/// touched by a kernel — when PRECALCLAM is off.
#[derive(Debug, Clone, Default)]
pub struct PrecalcNonbonded {
    /// λ_k of the grid.
    pub lambdas: Vec<f64>,
    n_cells: usize,
    /// `A_lj`, `B_lj`, `A_crf`, `B_crf` energies.
    pub energy: [Vec<f64>; 4],
    /// The same four on the derivative record (`A_de`, `B_de`).
    pub deriv: [Vec<f64>; 4],
}

impl PrecalcNonbonded {
    /// Zeroed arrays for the grid `grid` and `n_groups` energy groups.
    pub fn new(grid: &newen_core::lambdas::PrecalcLambdas, n_groups: usize) -> Self {
        let n_cells = n_groups.max(1) * n_groups.max(1);
        let len = grid.nr_lambdas * n_cells;
        Self {
            lambdas: (0..grid.nr_lambdas).map(|k| grid.lambda(k)).collect(),
            n_cells,
            energy: std::array::from_fn(|_| vec![0.0; len]),
            deriv: std::array::from_fn(|_| vec![0.0; len]),
        }
    }

    /// Whether the arrays are sized (PRECALCLAM on).
    #[inline(always)]
    pub fn is_active(&self) -> bool {
        !self.lambdas.is_empty()
    }

    /// Number of energy-group pairs per λ_k.
    pub fn n_cells(&self) -> usize {
        self.n_cells
    }

    pub fn clear(&mut self) {
        for v in self.energy.iter_mut().chain(self.deriv.iter_mut()) {
            v.fill(0.0);
        }
    }

    /// Element-wise sum of `other` into `self`; a no-op when `other` is empty.
    pub fn merge(&mut self, other: &PrecalcNonbonded) {
        if !other.is_active() {
            return;
        }
        assert_eq!(
            self.energy[0].len(),
            other.energy[0].len(),
            "PRECALCLAM grid mismatch"
        );
        for (dst, src) in self
            .energy
            .iter_mut()
            .chain(self.deriv.iter_mut())
            .zip(other.energy.iter().chain(other.deriv.iter()))
        {
            for (d, s) in dst.iter_mut().zip(src) {
                *d += *s;
            }
        }
    }

    /// Add one λ_k's `[A_lj, B_lj, A_crf, B_crf]` energies and derivatives to `cell`.
    #[inline(always)]
    pub fn book(&mut self, k: usize, cell: usize, e: [f64; 4], de: [f64; 4]) {
        let idx = k * self.n_cells + cell;
        for c in 0..4 {
            self.energy[c][idx] += e[c];
            self.deriv[c][idx] += de[c];
        }
    }
}

/// Storage for forces and energies
#[repr(C)]
#[derive(Debug, Clone)]
pub struct ForceStorage {
    pub forces: Vec<Vec3>,
    pub e_lj: f64,
    pub e_crf: f64,
    pub virial: [[f64; 3]; 3],
    /// dH/dλ of the Lennard-Jones part of these pairs (perturbed runs; zero otherwise)
    pub dhdl_lj: f64,
    /// dH/dλ of the reaction-field part of these pairs
    pub dhdl_crf: f64,
    /// The same energies and derivatives per energy-group pair (empty with one group); the
    /// scalars above stay the totals.
    pub groups: GroupEnergies,
    /// PRECALCLAM energies of these pairs (empty when off).
    pub precalc: PrecalcNonbonded,
}

impl ForceStorage {
    pub fn new(n_atoms: usize) -> Self {
        Self::with_groups(n_atoms, 1)
    }

    /// Storage that also books energies per energy-group pair when `n_groups > 1`.
    pub fn with_groups(n_atoms: usize, n_groups: usize) -> Self {
        Self {
            forces: vec![Vec3::ZERO; n_atoms],
            e_lj: 0.0,
            e_crf: 0.0,
            virial: [[0.0; 3]; 3],
            dhdl_lj: 0.0,
            dhdl_crf: 0.0,
            groups: GroupEnergies::new(n_groups),
            precalc: PrecalcNonbonded::default(),
        }
    }

    pub fn clear(&mut self) {
        self.forces.fill(Vec3::ZERO);
        self.e_lj = 0.0;
        self.e_crf = 0.0;
        self.virial = [[0.0; 3]; 3];
        self.dhdl_lj = 0.0;
        self.dhdl_crf = 0.0;
        self.groups.clear();
        self.precalc.clear();
    }

    /// Accumulate another ForceStorage into this one (thread-local reduction).
    pub fn merge(&mut self, other: &ForceStorage) {
        for (dst, src) in self.forces.iter_mut().zip(other.forces.iter()) {
            *dst += *src;
        }
        self.e_lj += other.e_lj;
        self.e_crf += other.e_crf;
        for a in 0..3 {
            for b in 0..3 {
                self.virial[a][b] += other.virial[a][b];
            }
        }
        self.dhdl_lj += other.dhdl_lj;
        self.dhdl_crf += other.dhdl_crf;
        self.groups.merge(&other.groups);
        self.precalc.merge(&other.precalc);
    }

    /// Panic unless this storage books the group pairs `groups` describes (checked once per
    /// kernel call, outside the pair loop).
    #[inline]
    pub(crate) fn check_groups(&self, groups: &EnergyGroups<'_>) {
        if groups.is_grouped() {
            assert_eq!(
                self.groups.lj.len(),
                groups.n_groups() * groups.n_groups(),
                "ForceStorage is not sized for {} energy groups",
                groups.n_groups()
            );
        }
    }
}

/// A group of atom pairs sharing the same PBC shift (same charge-group pair).
///
/// The `ref_atom_i` and `ref_atom_j` are used to compute nearest_image once;
/// all atom pairs in the group then use the resulting shift offset.
#[derive(Debug, Clone, Copy)]
pub struct CGPairGroup {
    /// Reference atom from CG i (typically first atom in the CG)
    pub ref_atom_i: u32,
    /// Reference atom from CG j (typically first atom in the CG)
    pub ref_atom_j: u32,
    /// Start index into the flat pair array
    pub start: u32,
    /// End index (exclusive) into the flat pair array
    pub end: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_reaction_field_of_the_interior_dielectric_is_a_plain_cutoff() {
        // ε_RF = ε → C_RF = 0: no quadratic term, no force correction, shift 1/R_c
        let crf = CRFParameters::new(1.4, 1.0, 1.0, 0.0);
        assert_eq!(crf.crf_2cut3i, 0.0);
        assert_eq!(crf.crf_cut3i, 0.0);
        assert!((crf.crf_cut - 1.0 / 1.4).abs() < 1e-15);
        assert!((crf.cutoff_sq - 1.96).abs() < 1e-12);
    }

    #[test]
    fn a_conducting_boundary_has_crf_minus_one() {
        // gromosXX nonbonded_term.cc:58-71: EPSRF = 0 → crf = −1, not the ε formula's 2
        let crf = CRFParameters::new(1.2, 1.0, 0.0, 0.0);
        let cut3i = 1.0 / (1.2 * 1.2 * 1.2);
        assert_eq!(crf.crf_cut3i, -cut3i);
        assert_eq!(crf.crf_2cut3i, -cut3i / 2.0);
        assert_eq!(crf.crf_cut, 1.5 / 1.2);
        assert_eq!(crf.cutoff_sq, 1.2 * 1.2);
    }

    #[test]
    fn spc_reaction_field_gives_the_gromos_constants() {
        // ε_RF = 61: C_RF = 2(1−61)/(1+122) = −120/123
        let crf = CRFParameters::new(1.4, 1.0, 61.0, 0.0);
        let c_rf = 2.0 * (1.0 - 61.0) / (1.0 + 2.0 * 61.0);
        assert!((crf.crf_cut3i - c_rf / 1.4_f64.powi(3)).abs() < 1e-15);
        assert!((crf.crf_2cut3i - 0.5 * crf.crf_cut3i).abs() < 1e-15);
        assert!((crf.crf_cut - (1.0 - c_rf / 2.0) / 1.4).abs() < 1e-15);
    }

    #[test]
    fn lj_matrix_is_flat_and_symmetric() {
        let mut a = LJParameters {
            c6: 0.0,
            c12: 0.0,
            cs6: 0.0,
            cs12: 0.0,
        };
        a.c6 = 1.0;
        let mut b = LJParameters {
            c6: 0.0,
            c12: 0.0,
            cs6: 0.0,
            cs12: 0.0,
        };
        b.c6 = 2.0;
        let mut ab = LJParameters {
            c6: 0.0,
            c12: 0.0,
            cs6: 0.0,
            cs12: 0.0,
        };
        ab.c6 = 3.0;
        let m = LJParamMatrix::from_nested(&[vec![a, ab], vec![ab, b]]);
        assert_eq!(m.n_types(), 2);
        assert_eq!(m.get(0, 0).c6, 1.0);
        assert_eq!(m.get(1, 1).c6, 2.0);
        assert_eq!(m.get(0, 1).c6, m.get(1, 0).c6);
    }
}
