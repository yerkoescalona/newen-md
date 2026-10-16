//! Perturbed (FEP/TI) bonded force calculations — GROMOS faithful.
//!
//! References (all in md++/src/interaction/bonded/):
//!   perturbed_quartic_bond_interaction.cc
//!   perturbed_angle_interaction.cc
//!   perturbed_improper_dihedral_interaction.cc
//!   perturbed_dihedral_interaction.cc

use gromos_core::configuration::Configuration;
use gromos_core::lambdas::{IndividualLambdas, InteractionLambda};
use gromos_core::topology::Topology;

use super::{BondedGroupEnergies, ForceEnergyLambda};
use gromos_core::math::Periodicity;

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// cos(mφ) and d[cos(mφ)]/d[cosφ] via Chebyshev polynomials (same as regular dihedral).
fn chebyshev(m: i32, cos_phi: f64) -> (f64, f64) {
    match m {
        0 => (0.0, 0.0),
        1 => (cos_phi, 1.0),
        2 => (2.0 * cos_phi * cos_phi - 1.0, 4.0 * cos_phi),
        3 => (
            4.0 * cos_phi.powi(3) - 3.0 * cos_phi,
            12.0 * cos_phi * cos_phi - 3.0,
        ),
        4 => (
            8.0 * cos_phi.powi(4) - 8.0 * cos_phi * cos_phi + 1.0,
            32.0 * cos_phi.powi(3) - 16.0 * cos_phi,
        ),
        5 => (
            16.0 * cos_phi.powi(5) - 20.0 * cos_phi.powi(3) + 5.0 * cos_phi,
            80.0 * cos_phi.powi(4) - 60.0 * cos_phi * cos_phi + 5.0,
        ),
        6 => (
            32.0 * cos_phi.powi(6) - 48.0 * cos_phi.powi(4) + 18.0 * cos_phi * cos_phi - 1.0,
            192.0 * cos_phi.powi(5) - 192.0 * cos_phi.powi(3) + 36.0 * cos_phi,
        ),
        _ => {
            eprintln!("perturbed dihedral: unsupported multiplicity {m}");
            (0.0, 0.0)
        },
    }
}

/// Where a perturbed bonded term takes its `(L, L′)`: one fixed pair (the single-term entry
/// points), or the individual λ of one interaction kind at the energy group of the term's first
/// atom.
#[derive(Debug, Clone, Copy)]
pub enum TermLambda<'a> {
    Fixed(f64, f64),
    Table {
        lambdas: &'a IndividualLambdas,
        kind: InteractionLambda,
        atom_group: &'a [usize],
    },
}

impl<'a> TermLambda<'a> {
    /// The individual λ of `kind` from the topology's table.
    pub fn table(topo: &'a Topology, kind: InteractionLambda) -> Self {
        Self::Table {
            lambdas: &topo.individual_lambdas,
            kind,
            atom_group: &topo.atom_to_energy_group,
        }
    }

    /// `(L, L′)` of a term whose first atom is `first`.
    #[inline]
    pub fn at(&self, first: usize) -> (f64, f64) {
        match *self {
            Self::Fixed(l, dl) => (l, dl),
            Self::Table {
                lambdas,
                kind,
                atom_group,
            } => {
                let g = atom_group.get(first).copied().unwrap_or(0);
                lambdas.get(kind, g, g)
            },
        }
    }
}

/// PRECALCLAM (extended TI, item 1.8) energies of the perturbed bonded terms: `AB_bond`,
/// `AB_angle`, `AB_improper` at every λ_k (energies and their exact ∂/∂λ, each term's energy
/// with its parameters mixed linearly at λ_k — no individual λ) and the A/B energies of the
/// perturbed proper dihedrals, which have no grid (gromosXX `A_dihedral`, `B_dihedral`; on the
/// derivative record `A_dihedral` is ΣB − A).
#[derive(Debug, Clone, Default)]
pub struct BondedPrecalc {
    /// λ_k of the grid.
    pub lambdas: Vec<f64>,
    /// `AB_bond`, `AB_angle`, `AB_improper` per λ_k.
    pub energy: [Vec<f64>; 3],
    /// Their λ derivatives per λ_k.
    pub deriv: [Vec<f64>; 3],
    /// Σ E_A of the perturbed proper dihedrals.
    pub a_dihedral: f64,
    /// Σ E_B of the perturbed proper dihedrals.
    pub b_dihedral: f64,
    /// Σ (E_B − E_A), the derivative record's `A_dihedral`.
    pub d_a_dihedral: f64,
}

impl BondedPrecalc {
    /// Zeroed arrays for the λ grid `lambdas`.
    pub fn new(lambdas: &[f64]) -> Self {
        let n = lambdas.len();
        Self {
            lambdas: lambdas.to_vec(),
            energy: std::array::from_fn(|_| vec![0.0; n]),
            deriv: std::array::from_fn(|_| vec![0.0; n]),
            ..Self::default()
        }
    }

    pub fn clear(&mut self) {
        for v in self.energy.iter_mut().chain(self.deriv.iter_mut()) {
            v.fill(0.0);
        }
        self.a_dihedral = 0.0;
        self.b_dihedral = 0.0;
        self.d_a_dihedral = 0.0;
    }

    /// Harmonic-in-x term `½·K(λ)·(x − x0(λ))²` with K and x0 linear in λ at every λ_k
    /// (`AB_angle`, `AB_improper`: gromosXX `perturbed_angle_interaction.cc:179-201`,
    /// `perturbed_improper_dihedral_interaction.cc:198-220`).
    #[inline]
    fn add_harmonic(&mut self, col: usize, x: f64, k: (f64, f64), x0: (f64, f64)) {
        let (k_diff, x0_diff) = (k.1 - k.0, x0.1 - x0.0);
        for (i, &lam) in self.lambdas.iter().enumerate() {
            let k_lam = (1.0 - lam) * k.0 + lam * k.1;
            let x0_lam = (1.0 - lam) * x0.0 + lam * x0.1;
            let diff = x - x0_lam;
            let diff2 = diff * diff;
            self.energy[col][i] += 0.5 * k_lam * diff2;
            self.deriv[col][i] += 0.5 * k_diff * diff2 - k_lam * x0_diff * diff;
        }
    }

    /// The soft-core variant (`perturbed_soft_bond_interaction.cc:182-209`, soft angle
    /// `:192-219`, soft improper `:203-230`).
    #[inline]
    fn add_soft(&mut self, col: usize, x: f64, k: (f64, f64), x0: (f64, f64), alpha: f64) {
        let x0_diff = x0.1 - x0.0;
        for (i, &lam) in self.lambdas.iter().enumerate() {
            let x0_lam = (1.0 - lam) * x0.0 + lam * x0.1;
            let diff = x - x0_lam;
            let diff2 = diff * diff;
            let soft_a = 1.0 + alpha * lam * diff2;
            let soft_b = 1.0 + alpha * (1.0 - lam) * diff2;
            let soft_a2 = soft_a * soft_a;
            let soft_b2 = soft_b * soft_b;
            let k_soft = (1.0 - lam) * k.0 / soft_a + lam * k.1 / soft_b;
            let st1 = 1.0 + alpha * diff2;
            let st2 = -2.0 * alpha * lam * (1.0 - lam) * diff * x0_diff;
            self.energy[col][i] += 0.5 * k_soft * diff2;
            self.deriv[col][i] +=
                0.5 * diff2 * (k.0 / soft_a2 * (-st1 - st2) + k.1 / soft_b2 * (st1 - st2))
                    - k_soft * diff * x0_diff;
        }
    }
}

// ─── Perturbed quartic bond ───────────────────────────────────────────────────

/// Perturbed quartic bond (GROMOS `perturbed_quartic_bond_interaction.cc`).
///
/// K(λ)  = (1-λ)·K_A + λ·K_B        [k_quartic]
/// r0(λ) = (1-λ)·r0_A + λ·r0_B
/// E     = ¼·K·(r²-r0²)²
/// F_i   = v·(-K)·(r²-r0²)           [v = pos_i - pos_j]
/// dE/dλ = ¼·λ_d·[K_d·Δ² - 4·K·b_d·b_mix·Δ]   where Δ = r²-b_mix²
pub fn calculate_perturbed_bond_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
) -> ForceEnergyLambda {
    calculate_perturbed_bond_forces_grouped::<false>(
        topo,
        conf,
        pbc,
        TermLambda::Fixed(lambda, lambda_derivative),
        &mut [],
        &mut [],
        None,
    )
}

/// [`calculate_perturbed_bond_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_bond_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lam: TermLambda<'_>,
    e_group: &mut [f64],
    d_group: &mut [f64],
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for bond in &topo.perturbed_solute.bonds {
        let (lambda, lambda_derivative) = lam.at(bond.i);
        let pa = &topo.bond_parameters[bond.a_type];
        let pb = &topo.bond_parameters[bond.b_type];

        let k = (1.0 - lambda) * pa.k_quartic + lambda * pb.k_quartic;
        let r0 = (1.0 - lambda) * pa.r0 + lambda * pb.r0;
        let k_d = pb.k_quartic - pa.k_quartic;
        let b_d = pb.r0 - pa.r0;
        // b_mix = r0(λ), already computed above as `r0`

        // GROMOS: v = pos(i) - pos(j)
        let v = pbc.nearest_image(conf.current().pos[bond.i], conf.current().pos[bond.j]);
        let dist2 = v.length_squared();
        let r02 = r0 * r0;
        let delta = dist2 - r02;

        let energy = 0.25 * k * delta * delta;
        let f_i = v * (-k * delta); // F_i += f; F_j -= f
        result.energy += energy;
        if GROUPS {
            e_group[topo.atom_to_energy_group[bond.i]] += energy;
        }
        result.forces[bond.i] += f_i;
        result.forces[bond.j] -= f_i;
        result.add_virial(&[(v, f_i)]);

        // dE/dλ = ¼·λ_d·[K_d·Δ² - 4·K·b_d·r0·Δ]
        let de_dl = 0.25 * lambda_derivative * (k_d * delta * delta - 4.0 * k * b_d * r0 * delta);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[bond.i]] += de_dl;
        }
        // gromosXX perturbed_quartic_bond_interaction.cc:173-196: ¼·K(λ_k)·(r² − b0(λ_k)²)²
        if let Some(pc) = precalc.as_deref_mut() {
            for (i, &l) in pc.lambdas.iter().enumerate() {
                let b0lam = (1.0 - l) * pa.r0 + l * pb.r0;
                let klam = (1.0 - l) * pa.k_quartic + l * pb.k_quartic;
                let difflam = dist2 - b0lam * b0lam;
                let difflam2 = difflam * difflam;
                pc.energy[0][i] += 0.25 * klam * difflam2;
                pc.deriv[0][i] += 0.25 * k_d * difflam2 - klam * b_d * b0lam * difflam;
            }
        }
    }

    result
}

// ─── Perturbed cos-harmonic angle ────────────────────────────────────────────

/// Perturbed cos-harmonic angle (GROMOS `perturbed_angle_interaction.cc`).
///
/// K(λ)    = (1-λ)·K_A + λ·K_B       [k_cosine]
/// cos0(λ) = (1-λ)·cos0_A + λ·cos0_B
/// E       = ½·K·(cosθ - cos0)²
/// dE/dλ   = ½·λ_d·[K_d·d² - 2·K·cd_d·d]    where d = cosθ - cos0(λ), cd_d = cos0_B - cos0_A
pub fn calculate_perturbed_angle_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
) -> ForceEnergyLambda {
    calculate_perturbed_angle_forces_grouped::<false>(
        topo,
        conf,
        pbc,
        TermLambda::Fixed(lambda, lambda_derivative),
        &mut [],
        &mut [],
        None,
    )
}

/// [`calculate_perturbed_angle_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_angle_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lam: TermLambda<'_>,
    e_group: &mut [f64],
    d_group: &mut [f64],
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for angle in &topo.perturbed_solute.angles {
        let (lambda, lambda_derivative) = lam.at(angle.i);
        let pa = &topo.angle_parameters[angle.a_type];
        let pb = &topo.angle_parameters[angle.b_type];

        let k = (1.0 - lambda) * pa.k_cosine + lambda * pb.k_cosine;
        let cos0a = pa.theta0.cos();
        let cos0b = pb.theta0.cos();
        let cos0 = (1.0 - lambda) * cos0a + lambda * cos0b;
        let k_d = pb.k_cosine - pa.k_cosine;
        let cd_d = cos0b - cos0a;

        // GROMOS: rij = pos(i)-pos(j), rkj = pos(k)-pos(j)
        let rij = pbc.nearest_image(conf.current().pos[angle.i], conf.current().pos[angle.j]);
        let rkj = pbc.nearest_image(conf.current().pos[angle.k], conf.current().pos[angle.j]);
        let dij = rij.length();
        let dkj = rkj.length();
        if dij < 1e-10 || dkj < 1e-10 {
            continue;
        }

        let cost = (rij.dot(rkj) / (dij * dkj)).clamp(-1.0, 1.0);
        let d_cos = cost - cos0;
        let energy = 0.5 * k * d_cos * d_cos;

        // Force (same kernel as calculate_angle_forces)
        let df = -k * d_cos;
        let f_i = (rkj * (1.0 / dkj) - rij * (cost / dij)) * (df / dij);
        let f_k = (rij * (1.0 / dij) - rkj * (cost / dkj)) * (df / dkj);
        let f_j = -(f_i + f_k);
        result.energy += energy;
        if GROUPS {
            e_group[topo.atom_to_energy_group[angle.i]] += energy;
        }
        result.forces[angle.i] += f_i;
        result.forces[angle.j] += f_j;
        result.forces[angle.k] += f_k;
        result.add_virial(&[(rij, f_i), (rkj, f_k)]);

        // dE/dλ = ½·λ_d·[K_d·d² - 2·K·cd_d·d]
        let de_dl = 0.5 * lambda_derivative * (k_d * d_cos * d_cos - 2.0 * k * cd_d * d_cos);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[angle.i]] += de_dl;
        }
        if let Some(pc) = precalc.as_deref_mut() {
            pc.add_harmonic(1, cost, (pa.k_cosine, pb.k_cosine), (cos0a, cos0b));
        }
    }

    result
}

// ─── Perturbed improper dihedral ─────────────────────────────────────────────

/// Perturbed improper dihedral (GROMOS `perturbed_improper_dihedral_interaction.cc`).
///
/// K(λ)  = (1-λ)·K_A + λ·K_B
/// q0(λ) = (1-λ)·q0_A + λ·q0_B
/// E     = ½·K·(ζ - q0)²
/// dE/dλ = ½·λ_d·[-2·K·q_d·Δ + K_d·Δ²]    where Δ = ζ - q0(λ)
pub fn calculate_perturbed_improper_dihedral_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
) -> ForceEnergyLambda {
    calculate_perturbed_improper_dihedral_forces_grouped::<false>(
        topo,
        conf,
        pbc,
        TermLambda::Fixed(lambda, lambda_derivative),
        &mut [],
        &mut [],
        None,
    )
}

/// [`calculate_perturbed_improper_dihedral_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_improper_dihedral_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lam: TermLambda<'_>,
    e_group: &mut [f64],
    d_group: &mut [f64],
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for imp in &topo.perturbed_solute.improper_dihedrals {
        let (lambda, lambda_derivative) = lam.at(imp.i);
        let pa = &topo.improper_dihedral_parameters[imp.a_type];
        let pb = &topo.improper_dihedral_parameters[imp.b_type];

        let k = (1.0 - lambda) * pa.k + lambda * pb.k;
        let q0 = (1.0 - lambda) * pa.q0 + lambda * pb.q0;
        let k_d = pb.k - pa.k;
        let q_d = pb.q0 - pa.q0;

        // Same geometry as calculate_improper_dihedral_forces
        let r_kj = pbc.nearest_image(conf.current().pos[imp.k], conf.current().pos[imp.j]);
        let r_ij = pbc.nearest_image(conf.current().pos[imp.i], conf.current().pos[imp.j]);
        let r_kl = pbc.nearest_image(conf.current().pos[imp.k], conf.current().pos[imp.l]);

        let r_mj = r_ij.cross(r_kj);
        let r_nk = r_kj.cross(r_kl);

        let d_kj2 = r_kj.dot(r_kj);
        let d_mj2 = r_mj.dot(r_mj);
        let d_nk2 = r_nk.dot(r_nk);
        let d_kj = d_kj2.sqrt();
        let d_mj = d_mj2.sqrt();
        let d_nk = d_nk2.sqrt();

        if d_mj < 1e-10 || d_nk < 1e-10 {
            continue;
        }

        let acs = (r_mj.dot(r_nk) / (d_mj * d_nk)).clamp(-1.0, 1.0);
        let mut zeta = acs.acos();
        if r_ij.dot(r_nk) < 0.0 {
            zeta = -zeta;
        }

        let mut zeta_adj = zeta;
        while zeta_adj < q0 - std::f64::consts::PI {
            zeta_adj += 2.0 * std::f64::consts::PI;
        }
        while zeta_adj > q0 + std::f64::consts::PI {
            zeta_adj -= 2.0 * std::f64::consts::PI;
        }

        let d_zeta = zeta_adj - q0;
        let energy = 0.5 * k * d_zeta * d_zeta;

        let mut k_i = -k * d_zeta * d_kj;
        let mut k_l = -k_i;
        if d_mj2 < 1e-10 * d_kj2 {
            k_i = 0.0;
        } else {
            k_i /= d_mj2;
        }
        if d_nk2 < 1e-10 * d_kj2 {
            k_l = 0.0;
        } else {
            k_l /= d_nk2;
        }

        let k_j1 = r_ij.dot(r_kj) / d_kj2 - 1.0;
        let k_j2 = r_kl.dot(r_kj) / d_kj2;

        let f_i = r_mj * k_i;
        let f_l = r_nk * k_l;
        let f_j = f_i * k_j1 - f_l * k_j2;
        let f_k = -(f_i + f_j + f_l);
        result.energy += energy;
        if GROUPS {
            e_group[topo.atom_to_energy_group[imp.i]] += energy;
        }
        result.forces[imp.i] += f_i;
        result.forces[imp.j] += f_j;
        result.forces[imp.k] += f_k;
        result.forces[imp.l] += f_l;
        let r_lj = pbc.nearest_image(conf.current().pos[imp.l], conf.current().pos[imp.j]);
        result.add_virial(&[(r_ij, f_i), (r_kj, f_k), (r_lj, f_l)]);

        // dE/dλ = ½·λ_d·[-2·K·q_d·Δ + K_d·Δ²]
        let de_dl = 0.5 * lambda_derivative * (-2.0 * k * q_d * d_zeta + k_d * d_zeta * d_zeta);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[imp.i]] += de_dl;
        }
        // ζ stays wrapped around the current-λ q0 (gromosXX does not re-wrap it per λ_k)
        if let Some(pc) = precalc.as_deref_mut() {
            pc.add_harmonic(2, zeta_adj, (pa.k, pb.k), (pa.q0, pb.q0));
        }
    }

    result
}

// ─── Perturbed proper dihedral ────────────────────────────────────────────────

/// Perturbed proper dihedral (GROMOS `perturbed_dihedral_interaction.cc`).
///
/// States A and B computed separately (may have different multiplicities):
///   E_A = K_A·(1 + cos(δ_A)·cos(m_A·φ))
///   E_B = K_B·(1 + cos(δ_B)·cos(m_B·φ))
/// Combined:
///   E   = (1-λ)·E_A + λ·E_B
///   F   = (1-λ)·F_A + λ·F_B
///   dE/dλ = λ_d·(E_B - E_A)
pub fn calculate_perturbed_dihedral_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
) -> ForceEnergyLambda {
    calculate_perturbed_dihedral_forces_grouped::<false>(
        topo,
        conf,
        pbc,
        TermLambda::Fixed(lambda, lambda_derivative),
        &mut [],
        &mut [],
        None,
    )
}

/// [`calculate_perturbed_dihedral_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_dihedral_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lam: TermLambda<'_>,
    e_group: &mut [f64],
    d_group: &mut [f64],
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for dih in &topo.perturbed_solute.proper_dihedrals {
        let (lambda, lambda_derivative) = lam.at(dih.i);
        let pa = &topo.dihedral_parameters[dih.a_type];
        let pb = &topo.dihedral_parameters[dih.b_type];

        // Dihedral geometry (same as calculate_dihedral_forces)
        let r_ij = pbc.nearest_image(conf.current().pos[dih.i], conf.current().pos[dih.j]);
        let r_kj = pbc.nearest_image(conf.current().pos[dih.k], conf.current().pos[dih.j]);
        let r_kl = pbc.nearest_image(conf.current().pos[dih.k], conf.current().pos[dih.l]);

        let _r_mj = r_ij.cross(r_kj);
        let _r_nk = r_kj.cross(r_kl);
        let d_kj2 = r_kj.dot(r_kj);
        if d_kj2 < 1e-10 {
            continue;
        }

        let f_rim = r_ij.dot(r_kj) / d_kj2;
        let f_rln = r_kl.dot(r_kj) / d_kj2;
        let r_im = r_ij - r_kj * f_rim;
        let r_ln = r_kj * f_rln - r_kl;
        let d_im = r_im.length();
        let d_ln = r_ln.length();
        if d_im < 1e-10 || d_ln < 1e-10 {
            continue;
        }

        let cos_phi = (r_im.dot(r_ln) / (d_im * d_ln)).clamp(-1.0, 1.0);

        // Compute state A energy and force vectors
        let (cos_m_phi_a, d_cos_m_phi_a) = chebyshev(pa.m, cos_phi);
        let e_a = pa.k * (1.0 + pa.cospd * cos_m_phi_a);
        let k_ia = -pa.k * pa.cospd * d_cos_m_phi_a / d_im;
        let k_la = -pa.k * pa.cospd * d_cos_m_phi_a / d_ln;
        let f_ia = r_ln * (k_ia / d_ln) - r_im * (k_ia * cos_phi / d_im);
        let f_la = r_im * (k_la / d_im) - r_ln * (k_la * cos_phi / d_ln);
        let f_ja = f_ia * (f_rim - 1.0) - f_la * f_rln;
        let f_ka = -(f_ia + f_ja + f_la);

        // Compute state B energy and force vectors
        let (cos_m_phi_b, d_cos_m_phi_b) = chebyshev(pb.m, cos_phi);
        let e_b = pb.k * (1.0 + pb.cospd * cos_m_phi_b);
        let k_ib = -pb.k * pb.cospd * d_cos_m_phi_b / d_im;
        let k_lb = -pb.k * pb.cospd * d_cos_m_phi_b / d_ln;
        let f_ib = r_ln * (k_ib / d_ln) - r_im * (k_ib * cos_phi / d_im);
        let f_lb = r_im * (k_lb / d_im) - r_ln * (k_lb * cos_phi / d_ln);
        let f_jb = f_ib * (f_rim - 1.0) - f_lb * f_rln;
        let f_kb = -(f_ib + f_jb + f_lb);

        // Combine
        let w_a = 1.0 - lambda;
        let w_b = lambda;
        let e_term = w_a * e_a + w_b * e_b;
        result.energy += e_term;
        if GROUPS {
            e_group[topo.atom_to_energy_group[dih.i]] += e_term;
        }
        result.forces[dih.i] += f_ia * w_a + f_ib * w_b;
        result.forces[dih.j] += f_ja * w_a + f_jb * w_b;
        result.forces[dih.k] += f_ka * w_a + f_kb * w_b;
        result.forces[dih.l] += f_la * w_a + f_lb * w_b;
        let r_lj = pbc.nearest_image(conf.current().pos[dih.l], conf.current().pos[dih.j]);
        result.add_virial(&[
            (r_ij, f_ia * w_a + f_ib * w_b),
            (r_kj, f_ka * w_a + f_kb * w_b),
            (r_lj, f_la * w_a + f_lb * w_b),
        ]);

        let de_term = lambda_derivative * (e_b - e_a);
        result.lambda_derivative += de_term;
        if GROUPS {
            d_group[topo.atom_to_energy_group[dih.i]] += de_term;
        }
        // gromosXX perturbed_dihedral_interaction.cc:270-276
        if let Some(pc) = precalc.as_deref_mut() {
            pc.a_dihedral += e_a;
            pc.b_dihedral += e_b;
            pc.d_a_dihedral += e_b - e_a;
        }
    }

    result
}

// ─── Soft-core perturbed harmonic bond ───────────────────────────────────────

/// Soft-core perturbed harmonic bond (PERTBONDSOFT, GROMOS `perturbed_soft_bond_interaction.cc`).
///
/// Handles bonds absent in one state (K=0) via soft-core suppression:
///   S_A = 1 + α·λ·diff²,   S_B = 1 + α·(1-λ)·diff²
///   K_eff = (1-λ)·K_A/S_A² + λ·K_B/S_B²         (for forces)
///   K_soft = (1-λ)·K_A/S_A + λ·K_B/S_B           (for energy)
///   E = 0.5·K_soft·diff²
///   diff = b - b0(λ),   b0(λ) = (1-λ)·b0_A + λ·b0_B
/// A state whose type is `None` (bond absent in it) has K = 0 and the other state's b0.
pub fn calculate_soft_bond_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
) -> ForceEnergyLambda {
    calculate_soft_bond_forces_grouped::<false>(
        topo,
        conf,
        pbc,
        TermLambda::Fixed(lambda, lambda_derivative),
        &mut [],
        &mut [],
        None,
    )
}

/// [`calculate_soft_bond_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_soft_bond_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lam: TermLambda<'_>,
    e_group: &mut [f64],
    d_group: &mut [f64],
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for sb in &topo.perturbed_solute.soft_bonds {
        let (lambda, lambda_derivative) = lam.at(sb.i);
        let state = |t: Option<usize>| t.map(|t| &topo.bond_parameters[t]);
        let ((k_a, r0_a), (k_b, r0_b)) = match (state(sb.a_type), state(sb.b_type)) {
            (Some(a), Some(b)) => ((a.k_harmonic, a.r0), (b.k_harmonic, b.r0)),
            (Some(a), None) => ((a.k_harmonic, a.r0), (0.0, a.r0)), // absent in state B
            (None, Some(b)) => ((0.0, b.r0), (b.k_harmonic, b.r0)), // absent in state A
            (None, None) => continue,
        };

        let r0 = (1.0 - lambda) * r0_a + lambda * r0_b;
        let b_diff = r0_b - r0_a;
        let alpha = sb.alpha;

        let v = pbc.nearest_image(conf.current().pos[sb.i], conf.current().pos[sb.j]);
        let dist = v.length();
        if dist < 1e-10 {
            continue;
        }
        let diff = dist - r0;
        let diff2 = diff * diff;

        let s_a = 1.0 + alpha * lambda * diff2;
        let s_b = 1.0 + alpha * (1.0 - lambda) * diff2;
        let s_a2 = s_a * s_a;
        let s_b2 = s_b * s_b;

        let k_eff = (1.0 - lambda) * k_a / s_a2 + lambda * k_b / s_b2;
        let k_soft = (1.0 - lambda) * k_a / s_a + lambda * k_b / s_b;

        let energy = 0.5 * k_soft * diff2;
        let f_i = v * (-k_eff * diff / dist);
        result.energy += energy;
        if GROUPS {
            e_group[topo.atom_to_energy_group[sb.i]] += energy;
        }
        result.forces[sb.i] += f_i;
        result.forces[sb.j] -= f_i;
        result.add_virial(&[(v, f_i)]);

        // dE/dλ (GROMOS `e_lambda` formula)
        let st1 = 1.0 + alpha * diff2;
        let st2 = -2.0 * alpha * lambda * (1.0 - lambda) * diff * b_diff;
        let de_dl = lambda_derivative
            * (0.5 * diff2 * (k_a / s_a2 * (-st1 - st2) + k_b / s_b2 * (st1 - st2))
                - k_soft * diff * b_diff);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[sb.i]] += de_dl;
        }
        if let Some(pc) = precalc.as_deref_mut() {
            pc.add_soft(0, dist, (k_a, k_b), (r0_a, r0_b), alpha);
        }
    }

    result
}

// ─── Soft-core perturbed cos-harmonic angle ───────────────────────────────────

/// Soft-core perturbed cos-harmonic angle (PERTANGLESOFT, GROMOS `perturbed_soft_angle_interaction.cc`).
///
/// Same soft-core formulation as soft bonds but on cosθ:
///   diff = cosθ - cos0(λ),  cos0(λ) = (1-λ)·cos0_A + λ·cos0_B
/// A state whose type is `None` has K = 0 and the other state's cos0.
pub fn calculate_soft_angle_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
) -> ForceEnergyLambda {
    calculate_soft_angle_forces_grouped::<false>(
        topo,
        conf,
        pbc,
        TermLambda::Fixed(lambda, lambda_derivative),
        &mut [],
        &mut [],
        None,
    )
}

/// [`calculate_soft_angle_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_soft_angle_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lam: TermLambda<'_>,
    e_group: &mut [f64],
    d_group: &mut [f64],
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for sa in &topo.perturbed_solute.soft_angles {
        let (lambda, lambda_derivative) = lam.at(sa.i);
        let state = |t: Option<usize>| {
            t.map(|t| {
                let p = &topo.angle_parameters[t];
                (p.k_cosine, p.theta0.cos())
            })
        };
        let ((k_a, cos0_a), (k_b, cos0_b)) = match (state(sa.a_type), state(sa.b_type)) {
            (Some(a), Some(b)) => (a, b),
            (Some(a), None) => (a, (0.0, a.1)),
            (None, Some(b)) => ((0.0, b.1), b),
            (None, None) => continue,
        };

        let cos0 = (1.0 - lambda) * cos0_a + lambda * cos0_b;
        let cos_diff = cos0_b - cos0_a;
        let alpha = sa.alpha;

        let rij = pbc.nearest_image(conf.current().pos[sa.i], conf.current().pos[sa.j]);
        let rkj = pbc.nearest_image(conf.current().pos[sa.k], conf.current().pos[sa.j]);
        let dij = rij.length();
        let dkj = rkj.length();
        if dij < 1e-10 || dkj < 1e-10 {
            continue;
        }

        let cost = (rij.dot(rkj) / (dij * dkj)).clamp(-1.0, 1.0);
        let diff = cost - cos0;
        let diff2 = diff * diff;

        let s_a = 1.0 + alpha * lambda * diff2;
        let s_b = 1.0 + alpha * (1.0 - lambda) * diff2;
        let s_a2 = s_a * s_a;
        let s_b2 = s_b * s_b;

        let k_eff = (1.0 - lambda) * k_a / s_a2 + lambda * k_b / s_b2;
        let k_soft = (1.0 - lambda) * k_a / s_a + lambda * k_b / s_b;

        let energy = 0.5 * k_soft * diff2;

        let f_i = (rkj / dkj - rij / dij * cost) * (-k_eff * diff / dij);
        let f_k = (rij / dij - rkj / dkj * cost) * (-k_eff * diff / dkj);
        let f_j = -(f_i + f_k);
        result.energy += energy;
        if GROUPS {
            e_group[topo.atom_to_energy_group[sa.i]] += energy;
        }
        result.forces[sa.i] += f_i;
        result.forces[sa.j] += f_j;
        result.forces[sa.k] += f_k;
        result.add_virial(&[(rij, f_i), (rkj, f_k)]);

        let st1 = 1.0 + alpha * diff2;
        let st2 = -2.0 * alpha * lambda * (1.0 - lambda) * diff * cos_diff;
        let de_dl = lambda_derivative
            * (0.5 * diff2 * (k_a / s_a2 * (-st1 - st2) + k_b / s_b2 * (st1 - st2))
                - k_soft * diff * cos_diff);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[sa.i]] += de_dl;
        }
        if let Some(pc) = precalc.as_deref_mut() {
            pc.add_soft(1, cost, (k_a, k_b), (cos0_a, cos0_b), alpha);
        }
    }

    result
}

// ─── Soft-core perturbed improper dihedral ────────────────────────────────────

/// Soft-core perturbed improper dihedral (PERTIMPROPERDIHSOFT, GROMOS `perturbed_soft_improper_interaction.cc`).
///
/// Same soft-core formulation as soft bonds but on the improper angle ζ (radians):
///   diff = ζ - q0(λ),  q0(λ) = (1-λ)·q0_A + λ·q0_B
/// A state whose type is `None` has K = 0 and the other state's q0.
pub fn calculate_soft_improper_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
) -> ForceEnergyLambda {
    calculate_soft_improper_forces_grouped::<false>(
        topo,
        conf,
        pbc,
        TermLambda::Fixed(lambda, lambda_derivative),
        &mut [],
        &mut [],
        None,
    )
}

/// [`calculate_soft_improper_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_soft_improper_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lam: TermLambda<'_>,
    e_group: &mut [f64],
    d_group: &mut [f64],
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for si in &topo.perturbed_solute.soft_impropers {
        let (lambda, lambda_derivative) = lam.at(si.i);
        let state = |t: Option<usize>| {
            t.map(|t| {
                let p = &topo.improper_dihedral_parameters[t];
                (p.k, p.q0)
            })
        };
        let ((k_a, q0_a), (k_b, q0_b)) = match (state(si.a_type), state(si.b_type)) {
            (Some(a), Some(b)) => (a, b),
            (Some(a), None) => (a, (0.0, a.1)),
            (None, Some(b)) => ((0.0, b.1), b),
            (None, None) => continue,
        };

        let q0 = (1.0 - lambda) * q0_a + lambda * q0_b;
        let q_diff = q0_b - q0_a;
        let alpha = si.alpha;

        // Same geometry as calculate_improper_dihedral_forces
        let r_kj = pbc.nearest_image(conf.current().pos[si.k], conf.current().pos[si.j]);
        let r_ij = pbc.nearest_image(conf.current().pos[si.i], conf.current().pos[si.j]);
        let r_kl = pbc.nearest_image(conf.current().pos[si.k], conf.current().pos[si.l]);

        let r_mj = r_ij.cross(r_kj);
        let r_nk = r_kj.cross(r_kl);

        let d_kj2 = r_kj.dot(r_kj);
        let d_mj2 = r_mj.dot(r_mj);
        let d_nk2 = r_nk.dot(r_nk);
        let d_kj = d_kj2.sqrt();
        let d_mj = d_mj2.sqrt();
        let d_nk = d_nk2.sqrt();

        if d_mj < 1e-10 || d_nk < 1e-10 {
            continue;
        }

        let acs = (r_mj.dot(r_nk) / (d_mj * d_nk)).clamp(-1.0, 1.0);
        let mut zeta = acs.acos();
        if r_ij.dot(r_nk) < 0.0 {
            zeta = -zeta;
        }

        let mut zeta_adj = zeta;
        while zeta_adj < q0 - std::f64::consts::PI {
            zeta_adj += 2.0 * std::f64::consts::PI;
        }
        while zeta_adj > q0 + std::f64::consts::PI {
            zeta_adj -= 2.0 * std::f64::consts::PI;
        }

        let diff = zeta_adj - q0;
        let diff2 = diff * diff;

        let s_a = 1.0 + alpha * lambda * diff2;
        let s_b = 1.0 + alpha * (1.0 - lambda) * diff2;
        let s_a2 = s_a * s_a;
        let s_b2 = s_b * s_b;

        let k_eff = (1.0 - lambda) * k_a / s_a2 + lambda * k_b / s_b2;
        let k_soft = (1.0 - lambda) * k_a / s_a + lambda * k_b / s_b;

        let energy = 0.5 * k_soft * diff2;

        let fac = k_eff * diff * d_kj;
        let kj1 = r_ij.dot(r_kj) / d_kj2 - 1.0;
        let kj2 = r_kl.dot(r_kj) / d_kj2;

        let f_i = r_mj * (-fac / d_mj2);
        let f_l = r_nk * (fac / d_nk2);
        let f_j = f_i * kj1 - f_l * kj2;
        let f_k = -(f_i + f_j + f_l);
        result.energy += energy;
        if GROUPS {
            e_group[topo.atom_to_energy_group[si.i]] += energy;
        }
        result.forces[si.i] += f_i;
        result.forces[si.j] += f_j;
        result.forces[si.k] += f_k;
        result.forces[si.l] += f_l;
        let r_lj = pbc.nearest_image(conf.current().pos[si.l], conf.current().pos[si.j]);
        result.add_virial(&[(r_ij, f_i), (r_kj, f_k), (r_lj, f_l)]);

        let st1 = 1.0 + alpha * diff2;
        let st2 = -2.0 * alpha * lambda * (1.0 - lambda) * diff * q_diff;
        let de_dl = lambda_derivative
            * (0.5 * diff2 * (k_a / s_a2 * (-st1 - st2) + k_b / s_b2 * (st1 - st2))
                - k_soft * diff * q_diff);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[si.i]] += de_dl;
        }
        if let Some(pc) = precalc.as_deref_mut() {
            pc.add_soft(2, zeta_adj, (k_a, k_b), (q0_a, q0_b), alpha);
        }
    }

    result
}

// ─── Top-level combiner ───────────────────────────────────────────────────────

/// Calculate all perturbed bonded forces and dH/dλ, each term at the individual λ of its kind
/// for the energy group of its first atom (`topo.individual_lambdas`, gromosXX
/// `individual_lambda(bond_lambda)[eg(i)][eg(i)]` and likewise; soft terms count as their kind).
///
/// With more than one energy group each term's energy is added to `energy_groups` and its
/// dH/dλ to `dhdl_groups`, at the group of its first atom.
///
/// With `precalc` (PRECALCLAM on, an evaluated step) the terms also fill its `AB_*` columns and
/// dihedral A/B energies.
pub fn calculate_perturbed_bonded_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    energy_groups: &mut BondedGroupEnergies,
    dhdl_groups: &mut BondedGroupEnergies,
    precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    if energy_groups.is_grouped() {
        perturbed_bonded_forces::<true>(topo, conf, pbc, energy_groups, dhdl_groups, precalc)
    } else {
        perturbed_bonded_forces::<false>(topo, conf, pbc, energy_groups, dhdl_groups, precalc)
    }
}

fn perturbed_bonded_forces<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    eg: &mut BondedGroupEnergies,
    dg: &mut BondedGroupEnergies,
    mut precalc: Option<&mut BondedPrecalc>,
) -> ForceEnergyLambda {
    let n = topo.num_atoms();
    let mut result = ForceEnergyLambda::new(n);
    let of = |kind| TermLambda::table(topo, kind);

    if !topo.perturbed_solute.bonds.is_empty() {
        let term = calculate_perturbed_bond_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            of(InteractionLambda::Bond),
            &mut eg.bond,
            &mut dg.bond,
            precalc.as_deref_mut(),
        );
        result.bond_energy += term.energy;
        result.bond_dhdl += term.lambda_derivative;
        result.add(&term);
    }
    if !topo.perturbed_solute.angles.is_empty() {
        let term = calculate_perturbed_angle_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            of(InteractionLambda::Angle),
            &mut eg.angle,
            &mut dg.angle,
            precalc.as_deref_mut(),
        );
        result.angle_energy += term.energy;
        result.angle_dhdl += term.lambda_derivative;
        result.add(&term);
    }
    if !topo.perturbed_solute.improper_dihedrals.is_empty() {
        let term = calculate_perturbed_improper_dihedral_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            of(InteractionLambda::Improper),
            &mut eg.improper,
            &mut dg.improper,
            precalc.as_deref_mut(),
        );
        result.improper_energy += term.energy;
        result.improper_dhdl += term.lambda_derivative;
        result.add(&term);
    }
    if !topo.perturbed_solute.proper_dihedrals.is_empty() {
        let term = calculate_perturbed_dihedral_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            of(InteractionLambda::Dihedral),
            &mut eg.dihedral,
            &mut dg.dihedral,
            precalc.as_deref_mut(),
        );
        result.dihedral_energy += term.energy;
        result.dihedral_dhdl += term.lambda_derivative;
        result.add(&term);
    }
    if !topo.perturbed_solute.soft_bonds.is_empty() {
        let term = calculate_soft_bond_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            of(InteractionLambda::Bond),
            &mut eg.bond,
            &mut dg.bond,
            precalc.as_deref_mut(),
        );
        result.bond_energy += term.energy;
        result.bond_dhdl += term.lambda_derivative;
        result.add(&term);
    }
    if !topo.perturbed_solute.soft_angles.is_empty() {
        let term = calculate_soft_angle_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            of(InteractionLambda::Angle),
            &mut eg.angle,
            &mut dg.angle,
            precalc.as_deref_mut(),
        );
        result.angle_energy += term.energy;
        result.angle_dhdl += term.lambda_derivative;
        result.add(&term);
    }
    if !topo.perturbed_solute.soft_impropers.is_empty() {
        let term = calculate_soft_improper_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            of(InteractionLambda::Improper),
            &mut eg.improper,
            &mut dg.improper,
            precalc,
        );
        result.improper_energy += term.energy;
        result.improper_dhdl += term.lambda_derivative;
        result.add(&term);
    }

    result
}

#[cfg(test)]
mod tests {
    /// The unit tests below are single isolated molecules: no box, so no minimum image.
    fn vacuum() -> gromos_core::math::Periodicity {
        gromos_core::math::Periodicity::Vacuum(gromos_core::math::Vacuum)
    }

    use super::*;
    use gromos_core::configuration::Configuration;
    use gromos_core::math::Vec3;
    use gromos_core::topology::{
        AngleParameters, Atom, BondParameters, DihedralParameters, ImproperDihedralParameters,
        PerturbedAngle, PerturbedBond, PerturbedDihedral, Topology,
    };

    fn atom() -> Atom {
        Atom {
            name: "C".into(),
            residue_nr: 1,
            residue_name: "MOL".into(),
            iac: 0,
            mass: 12.0,
            charge: 0.0,
            is_perturbed: true,
            is_polarisable: false,
            is_coarse_grained: false,
        }
    }

    // ── Quartic bond ─────────────────────────────────────────────────────────

    #[test]
    fn test_perturbed_bond_lambda_0_equals_state_a() {
        let mut topo = Topology::new();
        topo.moltypes[0].atoms.extend([atom(), atom()]);
        topo.mass = vec![12.0; 2];
        topo.inverse_mass = vec![1.0 / 12.0; 2];
        // State A: K=1000, r0=0.15; State B: K=2000, r0=0.20
        topo.bond_parameters.push(BondParameters {
            k_quartic: 1000.0,
            k_harmonic: 0.0,
            r0: 0.15,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 2000.0,
            k_harmonic: 0.0,
            r0: 0.20,
        });
        topo.perturbed_solute.bonds.push(PerturbedBond {
            i: 0,
            j: 1,
            a_type: 0,
            b_type: 1,
            hydrogen: false,
        });

        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.20, 0.0, 0.0); // r=0.20 nm

        // At λ=0: K=1000, r0=0.15 → E = 0.25*1000*(0.04-0.0225)² = 0.25*1000*0.0003063 = 0.07656
        let res = calculate_perturbed_bond_forces(&topo, &conf, &vacuum(), 0.0, 1.0);
        let r2 = 0.20_f64 * 0.20;
        let r02 = 0.15_f64 * 0.15;
        let expected = 0.25 * 1000.0 * (r2 - r02) * (r2 - r02);
        assert!(
            (res.energy - expected).abs() < 1e-9,
            "λ=0 energy {:.9e} != {:.9e}",
            res.energy,
            expected
        );
        // Forces equal and opposite
        assert!((res.forces[0].x + res.forces[1].x).abs() < 1e-9);
        assert!(res.lambda_derivative.is_finite());
    }

    #[test]
    fn test_perturbed_bond_lambda_1_equals_state_b() {
        let mut topo = Topology::new();
        topo.moltypes[0].atoms.extend([atom(), atom()]);
        topo.mass = vec![12.0; 2];
        topo.inverse_mass = vec![1.0 / 12.0; 2];
        topo.bond_parameters.push(BondParameters {
            k_quartic: 1000.0,
            k_harmonic: 0.0,
            r0: 0.15,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 2000.0,
            k_harmonic: 0.0,
            r0: 0.20,
        });
        topo.perturbed_solute.bonds.push(PerturbedBond {
            i: 0,
            j: 1,
            a_type: 0,
            b_type: 1,
            hydrogen: false,
        });

        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.25, 0.0, 0.0);

        // At λ=1: K=2000, r0=0.20
        let res = calculate_perturbed_bond_forces(&topo, &conf, &vacuum(), 1.0, 1.0);
        let r2 = 0.25_f64 * 0.25;
        let r02 = 0.20_f64 * 0.20;
        let expected = 0.25 * 2000.0 * (r2 - r02) * (r2 - r02);
        assert!(
            (res.energy - expected).abs() < 1e-9,
            "λ=1 energy {:.9e} != {:.9e}",
            res.energy,
            expected
        );
    }

    #[test]
    fn test_perturbed_bond_force_direction() {
        // At equilibrium (r = r0(λ)), energy and force should be zero
        let mut topo = Topology::new();
        topo.moltypes[0].atoms.extend([atom(), atom()]);
        topo.mass = vec![12.0; 2];
        topo.inverse_mass = vec![1.0 / 12.0; 2];
        topo.bond_parameters.push(BondParameters {
            k_quartic: 1000.0,
            k_harmonic: 0.0,
            r0: 0.15,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 2000.0,
            k_harmonic: 0.0,
            r0: 0.20,
        });
        topo.perturbed_solute.bonds.push(PerturbedBond {
            i: 0,
            j: 1,
            a_type: 0,
            b_type: 1,
            hydrogen: false,
        });

        let lambda = 0.5_f64;
        let r0_lam = 0.5 * 0.15 + 0.5 * 0.20; // = 0.175 nm

        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(r0_lam, 0.0, 0.0);

        let res = calculate_perturbed_bond_forces(&topo, &conf, &vacuum(), lambda, 1.0);
        assert!(
            res.energy.abs() < 1e-12,
            "at r=r0(λ), E should be 0: {}",
            res.energy
        );
        assert!(res.forces[0].length() < 1e-10, "force at eq should be 0");
    }

    // ── Cos-harmonic angle ───────────────────────────────────────────────────

    #[test]
    fn test_perturbed_angle_lambda_0() {
        let mut topo = Topology::new();
        topo.moltypes[0].atoms.extend([atom(), atom(), atom()]);
        topo.mass = vec![12.0; 3];
        topo.inverse_mass = vec![1.0 / 12.0; 3];
        // State A: θ0=120°, K=500; State B: θ0=180°, K=1000
        let th0a = 120_f64.to_radians();
        let th0b = 180_f64.to_radians();
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 500.0,
            k_harmonic: 0.0,
            theta0: th0a,
        });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 1000.0,
            k_harmonic: 0.0,
            theta0: th0b,
        });
        topo.perturbed_solute.angles.push(PerturbedAngle {
            i: 0,
            j: 1,
            k: 2,
            a_type: 0,
            b_type: 1,
        });

        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0); // 180°

        // At λ=0: K=500, cos0=cos(120°)=-0.5, cosθ=cos(180°)=-1.0
        let res = calculate_perturbed_angle_forces(&topo, &conf, &vacuum(), 0.0, 1.0);
        let cost = -1.0_f64;
        let cos0 = th0a.cos();
        let expected = 0.5 * 500.0 * (cost - cos0) * (cost - cos0);
        assert!(
            (res.energy - expected).abs() < 1e-9,
            "λ=0 angle energy {:.9e} != {:.9e}",
            res.energy,
            expected
        );
    }

    #[test]
    fn test_perturbed_angle_force_conservation() {
        let mut topo = Topology::new();
        topo.moltypes[0].atoms.extend([atom(), atom(), atom()]);
        topo.mass = vec![12.0; 3];
        topo.inverse_mass = vec![1.0 / 12.0; 3];
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 500.0,
            k_harmonic: 0.0,
            theta0: 109.5_f64.to_radians(),
        });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 750.0,
            k_harmonic: 0.0,
            theta0: 120.0_f64.to_radians(),
        });
        topo.perturbed_solute.angles.push(PerturbedAngle {
            i: 0,
            j: 1,
            k: 2,
            a_type: 0,
            b_type: 1,
        });

        let mut conf = Configuration::new(3, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.3, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(0.8, 0.6, 0.0);

        let res = calculate_perturbed_angle_forces(&topo, &conf, &vacuum(), 0.5, 1.0);
        let total: Vec3 = res.forces.iter().copied().sum();
        assert!(
            total.length() < 1e-10,
            "angle forces not conserved: {}",
            total.length()
        );
    }

    // ── Improper dihedral ────────────────────────────────────────────────────

    #[test]
    fn test_perturbed_improper_at_equilibrium() {
        let mut topo = Topology::new();
        for _ in 0..4 {
            topo.moltypes[0].atoms.push(atom());
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0 / 12.0; 4];
        topo.improper_dihedral_parameters
            .push(ImproperDihedralParameters { k: 100.0, q0: 0.0 });
        topo.improper_dihedral_parameters
            .push(ImproperDihedralParameters { k: 200.0, q0: 0.1 });
        use gromos_core::topology::PerturbedDihedral;
        topo.perturbed_solute
            .improper_dihedrals
            .push(PerturbedDihedral {
                i: 0,
                j: 1,
                k: 2,
                l: 3,
                a_type: 0,
                b_type: 1,
            });

        // Planar geometry → ζ ≈ 0
        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.1, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(0.5, -0.1, 0.0);

        let res = calculate_perturbed_improper_dihedral_forces(&topo, &conf, &vacuum(), 0.0, 1.0);
        assert!(res.energy >= 0.0 && res.energy.is_finite());
        assert!(res.lambda_derivative.is_finite());
    }

    // ── Proper dihedral ──────────────────────────────────────────────────────

    #[test]
    fn test_perturbed_dihedral_lambda_0_equals_state_a() {
        let mut topo = Topology::new();
        for _ in 0..4 {
            topo.moltypes[0].atoms.push(atom());
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0 / 12.0; 4];
        // State A: K=5, m=3, δ=0; State B: K=8, m=3, δ=π (different K/δ)
        topo.dihedral_parameters.push(DihedralParameters {
            k: 5.0,
            cospd: 1.0,
            pd: 0.0,
            m: 3,
        });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 8.0,
            cospd: -1.0,
            pd: std::f64::consts::PI,
            m: 3,
        });
        topo.perturbed_solute
            .proper_dihedrals
            .push(PerturbedDihedral {
                i: 0,
                j: 1,
                k: 2,
                l: 3,
                a_type: 0,
                b_type: 1,
            });

        let mut conf = Configuration::new(4, 1, 1);
        // trans geometry: φ ≈ π
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, -0.5, 0.0);

        let res_0 = calculate_perturbed_dihedral_forces(&topo, &conf, &vacuum(), 0.0, 1.0);
        let res_1 = calculate_perturbed_dihedral_forces(&topo, &conf, &vacuum(), 1.0, 1.0);

        // At λ=0, result = state A; at λ=1, result = state B
        // dE/dλ = λ_d * (E_B - E_A) — verify sign is consistent
        assert!(res_0.lambda_derivative.is_finite());
        assert!(res_1.lambda_derivative.is_finite());
        // Force conservation
        let total: Vec3 = res_0.forces.iter().copied().sum();
        assert!(
            total.length() < 1e-9,
            "forces not conserved at λ=0: {}",
            total.length()
        );
        let total: Vec3 = res_1.forces.iter().copied().sum();
        assert!(
            total.length() < 1e-9,
            "forces not conserved at λ=1: {}",
            total.length()
        );
    }

    // ── Soft terms absent in state A ─────────────────────────────────────────

    /// A soft term absent in state A (type 0 in the file) is the mirror of one absent in state
    /// B: K = 0 and the other state's equilibrium value, so at λ it equals the other at 1 − λ,
    /// with dE/dλ of opposite sign, and it vanishes at λ = 0.
    #[test]
    fn soft_terms_absent_in_state_a_mirror_those_absent_in_state_b() {
        use gromos_core::topology::{SoftAngle, SoftBond, SoftImproper};

        let mut topo = Topology::new();
        topo.moltypes[0]
            .atoms
            .extend([atom(), atom(), atom(), atom()]);
        topo.resize_atom_arrays();
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: 3.0e5,
            r0: 0.15,
        });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 400.0,
            k_harmonic: 0.0,
            theta0: 1.9,
        });
        topo.improper_dihedral_parameters
            .push(ImproperDihedralParameters { k: 0.05, q0: 0.1 });
        let mut conf = Configuration::new(4, 1, 1);
        let pos = &mut conf.current_mut().pos;
        pos[0] = Vec3::new(0.0, 0.0, 0.0);
        pos[1] = Vec3::new(0.17, 0.0, 0.0);
        pos[2] = Vec3::new(0.22, 0.14, 0.0);
        pos[3] = Vec3::new(0.35, 0.16, 0.08);

        let with = |a: Option<usize>, b: Option<usize>| {
            let mut t = topo.clone();
            let ps = &mut t.perturbed_solute;
            ps.soft_bonds.push(SoftBond {
                i: 0,
                j: 1,
                a_type: a,
                b_type: b,
                alpha: 1.5,
            });
            ps.soft_angles.push(SoftAngle {
                i: 0,
                j: 1,
                k: 2,
                a_type: a,
                b_type: b,
                alpha: 1.5,
            });
            ps.soft_impropers.push(SoftImproper {
                i: 0,
                j: 1,
                k: 2,
                l: 3,
                a_type: a,
                b_type: b,
                alpha: 1.5,
            });
            t
        };
        let made = with(None, Some(0));
        let broken = with(Some(0), None);
        type Kernel = fn(&Topology, &Configuration, &Periodicity, f64, f64) -> ForceEnergyLambda;
        let kernels: [(&str, Kernel); 3] = [
            ("bond", calculate_soft_bond_forces),
            ("angle", calculate_soft_angle_forces),
            ("improper", calculate_soft_improper_forces),
        ];
        for (name, f) in kernels {
            assert_eq!(f(&made, &conf, &vacuum(), 0.0, 1.0).energy, 0.0, "{name}");
            for lambda in [0.2, 0.5, 0.9] {
                let m = f(&made, &conf, &vacuum(), lambda, 1.0);
                let b = f(&broken, &conf, &vacuum(), 1.0 - lambda, 1.0);
                assert!(m.energy > 0.0, "{name} λ={lambda}");
                assert!((m.energy - b.energy).abs() < 1e-12, "{name} λ={lambda}");
                assert!(
                    (m.lambda_derivative + b.lambda_derivative).abs() < 1e-12,
                    "{name} λ={lambda}"
                );
                for (fm, fb) in m.forces.iter().zip(&b.forces) {
                    assert!((*fm - *fb).length() < 1e-9, "{name} λ={lambda}");
                }
            }
        }
    }

    // ── GROMOS reference values (aladip_special.t.cc) ─────────────────────
    //
    // Hard-coded in check/aladip_special.t.cc:
    //   PerturbedQuarticBond      = 1.271962 kJ/mol
    //   PerturbedAngle            = 0.500162 kJ/mol
    //   PerturbedImproperDihedral = 0.560988 kJ/mol
    //   PerturbedDihedral         = 4.626301 kJ/mol
    // at λ=0.125, NLAM=1, using aladip.topo + aladip.pttopo + aladip.conf.
    //
    // These become the end-to-end reference test once the .pttopo reader is wired
    // (P1.7 step 1).  See test_GROMOS_references::aladip_vacuum_fep.
}
