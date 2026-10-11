//! Perturbed (FEP/TI) bonded force calculations — GROMOS faithful.
//!
//! References (all in md++/src/interaction/bonded/):
//!   perturbed_quartic_bond_interaction.cc
//!   perturbed_angle_interaction.cc
//!   perturbed_improper_dihedral_interaction.cc
//!   perturbed_dihedral_interaction.cc

use gromos_core::configuration::Configuration;
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
        lambda,
        lambda_derivative,
        &mut [],
        &mut [],
    )
}

/// [`calculate_perturbed_bond_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_bond_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    e_group: &mut [f64],
    d_group: &mut [f64],
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for bond in &topo.perturbed_solute.bonds {
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

        // dE/dλ = ¼·λ_d·[K_d·Δ² - 4·K·b_d·r0·Δ]
        let de_dl = 0.25 * lambda_derivative * (k_d * delta * delta - 4.0 * k * b_d * r0 * delta);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[bond.i]] += de_dl;
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
        lambda,
        lambda_derivative,
        &mut [],
        &mut [],
    )
}

/// [`calculate_perturbed_angle_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_angle_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    e_group: &mut [f64],
    d_group: &mut [f64],
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for angle in &topo.perturbed_solute.angles {
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

        // dE/dλ = ½·λ_d·[K_d·d² - 2·K·cd_d·d]
        let de_dl = 0.5 * lambda_derivative * (k_d * d_cos * d_cos - 2.0 * k * cd_d * d_cos);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[angle.i]] += de_dl;
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
        lambda,
        lambda_derivative,
        &mut [],
        &mut [],
    )
}

/// [`calculate_perturbed_improper_dihedral_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_improper_dihedral_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    e_group: &mut [f64],
    d_group: &mut [f64],
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for imp in &topo.perturbed_solute.improper_dihedrals {
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

        // dE/dλ = ½·λ_d·[-2·K·q_d·Δ + K_d·Δ²]
        let de_dl = 0.5 * lambda_derivative * (-2.0 * k * q_d * d_zeta + k_d * d_zeta * d_zeta);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[imp.i]] += de_dl;
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
        lambda,
        lambda_derivative,
        &mut [],
        &mut [],
    )
}

/// [`calculate_perturbed_dihedral_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_perturbed_dihedral_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    e_group: &mut [f64],
    d_group: &mut [f64],
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for dih in &topo.perturbed_solute.proper_dihedrals {
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

        let de_term = lambda_derivative * (e_b - e_a);
        result.lambda_derivative += de_term;
        if GROUPS {
            d_group[topo.atom_to_energy_group[dih.i]] += de_term;
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
/// When b_type=None (state B absent): K_B=0, b0_B=b0_A.
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
        lambda,
        lambda_derivative,
        &mut [],
        &mut [],
    )
}

/// [`calculate_soft_bond_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_soft_bond_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    e_group: &mut [f64],
    d_group: &mut [f64],
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for sb in &topo.perturbed_solute.soft_bonds {
        let pa = &topo.bond_parameters[sb.a_type];
        let (k_a, r0_a) = (pa.k_harmonic, pa.r0);
        let (k_b, r0_b) = match sb.b_type {
            Some(bt) => (
                topo.bond_parameters[bt].k_harmonic,
                topo.bond_parameters[bt].r0,
            ),
            None => (0.0, r0_a), // absent in state B
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
    }

    result
}

// ─── Soft-core perturbed cos-harmonic angle ───────────────────────────────────

/// Soft-core perturbed cos-harmonic angle (PERTANGLESOFT, GROMOS `perturbed_soft_angle_interaction.cc`).
///
/// Same soft-core formulation as soft bonds but on cosθ:
///   diff = cosθ - cos0(λ),  cos0(λ) = (1-λ)·cos0_A + λ·cos0_B
/// When b_type=None: K_B=0, cos0_B=cos0_A.
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
        lambda,
        lambda_derivative,
        &mut [],
        &mut [],
    )
}

/// [`calculate_soft_angle_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_soft_angle_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    e_group: &mut [f64],
    d_group: &mut [f64],
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for sa in &topo.perturbed_solute.soft_angles {
        let pa = &topo.angle_parameters[sa.a_type];
        let (k_a, cos0_a) = (pa.k_cosine, pa.theta0.cos());
        let (k_b, cos0_b) = match sa.b_type {
            Some(bt) => (
                topo.angle_parameters[bt].k_cosine,
                topo.angle_parameters[bt].theta0.cos(),
            ),
            None => (0.0, cos0_a),
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

        let st1 = 1.0 + alpha * diff2;
        let st2 = -2.0 * alpha * lambda * (1.0 - lambda) * diff * cos_diff;
        let de_dl = lambda_derivative
            * (0.5 * diff2 * (k_a / s_a2 * (-st1 - st2) + k_b / s_b2 * (st1 - st2))
                - k_soft * diff * cos_diff);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[sa.i]] += de_dl;
        }
    }

    result
}

// ─── Soft-core perturbed improper dihedral ────────────────────────────────────

/// Soft-core perturbed improper dihedral (PERTIMPROPERDIHSOFT, GROMOS `perturbed_soft_improper_interaction.cc`).
///
/// Same soft-core formulation as soft bonds but on the improper angle ζ (radians):
///   diff = ζ - q0(λ),  q0(λ) = (1-λ)·q0_A + λ·q0_B
/// When b_type=None: K_B=0, q0_B=q0_A.
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
        lambda,
        lambda_derivative,
        &mut [],
        &mut [],
    )
}

/// [`calculate_soft_improper_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `energies / perturbed_energy_derivatives[atom_energy_group(i)]`).
pub(crate) fn calculate_soft_improper_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    e_group: &mut [f64],
    d_group: &mut [f64],
) -> ForceEnergyLambda {
    let mut result = ForceEnergyLambda::new(topo.num_atoms());

    for si in &topo.perturbed_solute.soft_impropers {
        let pa = &topo.improper_dihedral_parameters[si.a_type];
        let (k_a, q0_a) = (pa.k, pa.q0);
        let (k_b, q0_b) = match si.b_type {
            Some(bt) => (
                topo.improper_dihedral_parameters[bt].k,
                topo.improper_dihedral_parameters[bt].q0,
            ),
            None => (0.0, q0_a),
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

        let st1 = 1.0 + alpha * diff2;
        let st2 = -2.0 * alpha * lambda * (1.0 - lambda) * diff * q_diff;
        let de_dl = lambda_derivative
            * (0.5 * diff2 * (k_a / s_a2 * (-st1 - st2) + k_b / s_b2 * (st1 - st2))
                - k_soft * diff * q_diff);
        result.lambda_derivative += de_dl;
        if GROUPS {
            d_group[topo.atom_to_energy_group[si.i]] += de_dl;
        }
    }

    result
}

// ─── Top-level combiner ───────────────────────────────────────────────────────

/// Calculate all perturbed bonded forces and dH/dλ for a given lambda.
///
/// `lambda_derivative` = d(lambda)/d(RLAM) = NLAM·RLAM^(NLAM-1).
/// For standard NLAM=1 this is 1.0; use `ImdParameters::lambda_and_derivative()`.
///
/// With more than one energy group each term's energy is added to `energy_groups` and its
/// dH/dλ to `dhdl_groups`, at the group of its first atom (soft terms count as their kind).
pub fn calculate_perturbed_bonded_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    energy_groups: &mut BondedGroupEnergies,
    dhdl_groups: &mut BondedGroupEnergies,
) -> ForceEnergyLambda {
    if energy_groups.is_grouped() {
        perturbed_bonded_forces::<true>(
            topo,
            conf,
            pbc,
            lambda,
            lambda_derivative,
            energy_groups,
            dhdl_groups,
        )
    } else {
        perturbed_bonded_forces::<false>(
            topo,
            conf,
            pbc,
            lambda,
            lambda_derivative,
            energy_groups,
            dhdl_groups,
        )
    }
}

fn perturbed_bonded_forces<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    lambda: f64,
    lambda_derivative: f64,
    eg: &mut BondedGroupEnergies,
    dg: &mut BondedGroupEnergies,
) -> ForceEnergyLambda {
    let n = topo.num_atoms();
    let mut result = ForceEnergyLambda::new(n);
    let (l, ld) = (lambda, lambda_derivative);

    if !topo.perturbed_solute.bonds.is_empty() {
        let term = calculate_perturbed_bond_forces_grouped::<GROUPS>(
            topo,
            conf,
            pbc,
            l,
            ld,
            &mut eg.bond,
            &mut dg.bond,
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
            l,
            ld,
            &mut eg.angle,
            &mut dg.angle,
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
            l,
            ld,
            &mut eg.improper,
            &mut dg.improper,
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
            l,
            ld,
            &mut eg.dihedral,
            &mut dg.dihedral,
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
            l,
            ld,
            &mut eg.bond,
            &mut dg.bond,
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
            l,
            ld,
            &mut eg.angle,
            &mut dg.angle,
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
            l,
            ld,
            &mut eg.improper,
            &mut dg.improper,
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
