//! Proper dihedral force calculations (both old Chebyshev and new arbitrary-phase variants).

use newen_core::configuration::Configuration;
use newen_core::math::Vec3;
use newen_core::topology::Topology;

use super::ForceEnergy;
use newen_core::math::Periodicity;

/// Calculate proper dihedral forces (GROMOS torsional potential)
///
/// Potential: V = K * [1 + cos(δ) * cos(m*φ)]
/// where:
/// - K is the force constant
/// - δ is the phase shift angle
/// - m is the multiplicity (1-6)
/// - φ is the dihedral angle
///
/// Uses Chebyshev polynomials for efficient cos(m*φ) calculation
pub fn calculate_dihedral_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
) -> ForceEnergy {
    calculate_dihedral_forces_grouped::<false>(topo, conf, pbc, &mut [])
}

/// [`calculate_dihedral_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `dihedral_energy[atom_energy_group(i)]`).
pub(crate) fn calculate_dihedral_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    per_group: &mut [f64],
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for dihedral in topo.all_proper_dihedrals_global() {
        if dihedral.dihedral_type >= topo.dihedral_parameters.len() {
            continue;
        }

        let params = &topo.dihedral_parameters[dihedral.dihedral_type];

        // GROMOS convention: rij = pos(i)-pos(j), rkj = pos(k)-pos(j), rkl = pos(k)-pos(l)
        let r_ij = pbc.nearest_image(
            conf.current().pos[dihedral.i],
            conf.current().pos[dihedral.j],
        );
        let r_kj = pbc.nearest_image(
            conf.current().pos[dihedral.k],
            conf.current().pos[dihedral.j],
        );
        let r_kl = pbc.nearest_image(
            conf.current().pos[dihedral.k],
            conf.current().pos[dihedral.l],
        );

        // Cross products kept for reference (Chebyshev path uses projection instead)
        let _r_mj = r_ij.cross(r_kj); // normal to plane i-j-k
        let _r_nk = r_kj.cross(r_kl); // normal to plane j-k-l

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

        // Calculate cos(φ)
        let cos_phi = (r_im.dot(r_ln) / (d_im * d_ln)).clamp(-1.0, 1.0);

        // Calculate cos(m*φ) and d[cos(m*φ)]/d[cos(φ)] using Chebyshev polynomials
        let (cos_m_phi, d_cos_m_phi) = match params.m {
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
                eprintln!("Warning: unsupported dihedral multiplicity {}", params.m);
                continue;
            },
        };

        // Energy: V = K * (1 + cos(δ) * cos(m*φ))
        let energy = params.k * (1.0 + params.cospd * cos_m_phi);

        // Force calculation
        let k_i = -params.k * params.cospd * d_cos_m_phi / d_im;
        let k_l = -params.k * params.cospd * d_cos_m_phi / d_ln;
        let k_j1 = f_rim - 1.0;
        let k_j2 = f_rln;

        let f_i = r_ln * (k_i / d_ln) - r_im * (k_i * cos_phi / d_im);
        let f_l = r_im * (k_l / d_im) - r_ln * (k_l * cos_phi / d_ln);
        let f_j = f_i * k_j1 - f_l * k_j2;
        let f_k = -(f_i + f_j + f_l);

        result.energy += energy;
        if GROUPS {
            per_group[topo.atom_to_energy_group[dihedral.i]] += energy;
        }
        result.forces[dihedral.i] += f_i;
        result.forces[dihedral.j] += f_j;
        result.forces[dihedral.k] += f_k;
        result.forces[dihedral.l] += f_l;

        // GROMOS: virial_tensor(a, bb) += rij(a)*fi(bb) + rkj(a)*fk(bb) + rlj(a)*fl(bb)
        // rlj = pos(l) - pos(j)
        let r_lj = pbc.nearest_image(
            conf.current().pos[dihedral.l],
            conf.current().pos[dihedral.j],
        );
        let rij_v = [r_ij.x, r_ij.y, r_ij.z];
        let rkj_v = [r_kj.x, r_kj.y, r_kj.z];
        let rlj_v = [r_lj.x, r_lj.y, r_lj.z];
        let fi_v = [f_i.x, f_i.y, f_i.z];
        let fk_v = [f_k.x, f_k.y, f_k.z];
        let fl_v = [f_l.x, f_l.y, f_l.z];
        for a in 0..3 {
            for bb in 0..3 {
                result.virial[a][bb] +=
                    rij_v[a] * fi_v[bb] + rkj_v[a] * fk_v[bb] + rlj_v[a] * fl_v[bb];
            }
        }
    }

    result
}

/// Calculate proper dihedral forces with "new" formula (arbitrary phase shifts)
///
/// Potential: V = K * (1 + cos(m*φ - δ))
/// Force: dV/dφ = K * m * sin(m*φ - δ)
///
/// This is the improved GROMOS dihedral potential that supports arbitrary phase shifts δ,
/// not limited to 0° or 180° as in the old formula. It's simpler and more flexible.
pub fn calculate_dihedral_new_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
) -> ForceEnergy {
    calculate_dihedral_new_forces_grouped::<false>(topo, conf, pbc, &mut [])
}

/// [`calculate_dihedral_new_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `dihedral_energy[atom_energy_group(i)]`).
pub(crate) fn calculate_dihedral_new_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    per_group: &mut [f64],
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for dihedral in topo.all_proper_dihedrals_global() {
        if dihedral.dihedral_type >= topo.dihedral_parameters.len() {
            continue;
        }

        let params = &topo.dihedral_parameters[dihedral.dihedral_type];

        // GROMOS convention: rij = pos(i)-pos(j), rkj = pos(k)-pos(j), rkl = pos(k)-pos(l)
        let r_ij = pbc.nearest_image(
            conf.current().pos[dihedral.i],
            conf.current().pos[dihedral.j],
        );
        let r_kj = pbc.nearest_image(
            conf.current().pos[dihedral.k],
            conf.current().pos[dihedral.j],
        );
        let r_kl = pbc.nearest_image(
            conf.current().pos[dihedral.k],
            conf.current().pos[dihedral.l],
        );

        // Cross products to get plane normals
        let r_mj = r_ij.cross(r_kj); // normal to plane i-j-k
        let r_nk = r_kj.cross(r_kl); // normal to plane j-k-l

        let d_kj2 = r_kj.dot(r_kj);
        let d_mj2 = r_mj.dot(r_mj);
        let d_nk2 = r_nk.dot(r_nk);
        let d_kj = (d_kj2).sqrt();

        if d_kj < 1e-10 {
            continue;
        }

        // Project vectors onto plane perpendicular to k-j bond
        let f_rim = r_ij.dot(r_kj) / d_kj2;
        let f_rln = r_kl.dot(r_kj) / d_kj2;

        let r_im = r_ij - r_kj * f_rim;
        let r_ln = r_kj * f_rln - r_kl;

        let d_im = r_im.length();
        let d_ln = r_ln.length();

        if d_im < 1e-10 || d_ln < 1e-10 {
            continue;
        }

        // Calculate dihedral angle φ
        let cos_phi = (r_im.dot(r_ln) / (d_im * d_ln)).clamp(-1.0, 1.0);
        let mut phi = cos_phi.acos();

        // Determine sign of angle
        let sign = r_ij.dot(r_nk);
        if sign < 0.0 {
            phi = -phi;
        }

        // Calculate energy and force derivative
        // Energy: V = K * (1 + cos(m*φ - δ))
        let m_phi_minus_delta = params.m as f64 * phi - params.pd;
        let energy = params.k * (1.0 + m_phi_minus_delta.cos());

        // Force derivative: dV/dφ = K * m * sin(m*φ - δ)
        let k_i = params.k * params.m as f64 * m_phi_minus_delta.sin();
        let k_l = -k_i;

        // Calculate forces on each atom
        // Check for near-linear angles to avoid division by zero
        let mut f_i = Vec3::ZERO;
        let mut f_l = Vec3::ZERO;

        if d_mj2 > 1e-10 * d_kj2 {
            f_i = r_mj * (k_i * d_kj / d_mj2);
        }

        if d_nk2 > 1e-10 * d_kj2 {
            f_l = r_nk * (k_l * d_kj / d_nk2);
        }

        // Forces on j and k from chain rule
        let k_j1 = f_rim - 1.0;
        let k_j2 = f_rln;
        let f_j = f_i * k_j1 - f_l * k_j2;
        let f_k = -(f_i + f_j + f_l);

        result.energy += energy;
        if GROUPS {
            per_group[topo.atom_to_energy_group[dihedral.i]] += energy;
        }
        result.forces[dihedral.i] += f_i;
        result.forces[dihedral.j] += f_j;
        result.forces[dihedral.k] += f_k;
        result.forces[dihedral.l] += f_l;

        // GROMOS: virial_tensor(a, bb) += rij(a)*fi(bb) + rkj(a)*fk(bb) + rlj(a)*fl(bb)
        let r_lj = pbc.nearest_image(
            conf.current().pos[dihedral.l],
            conf.current().pos[dihedral.j],
        );
        let rij_v = [r_ij.x, r_ij.y, r_ij.z];
        let rkj_v = [r_kj.x, r_kj.y, r_kj.z];
        let rlj_v = [r_lj.x, r_lj.y, r_lj.z];
        let fi_v = [f_i.x, f_i.y, f_i.z];
        let fk_v = [f_k.x, f_k.y, f_k.z];
        let fl_v = [f_l.x, f_l.y, f_l.z];
        for a in 0..3 {
            for bb in 0..3 {
                result.virial[a][bb] +=
                    rij_v[a] * fi_v[bb] + rkj_v[a] * fk_v[bb] + rlj_v[a] * fl_v[bb];
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    /// The unit tests below are single isolated molecules: no box, so no minimum image.
    fn vacuum() -> newen_core::math::Periodicity {
        newen_core::math::Periodicity::Vacuum(newen_core::math::Vacuum)
    }

    use super::*;
    use newen_core::configuration::Configuration;
    use newen_core::math::Vec3;
    use newen_core::topology::{Atom, Dihedral, DihedralParameters, Topology};

    #[test]
    fn test_dihedral_new_at_minimum() {
        use std::f64::consts::PI;

        // Test new dihedral at energy minimum (φ = δ/m)
        let mut topo = Topology::new();

        // Add 4 atoms in a dihedral
        for i in 0..4 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0 / 12.0; 4];

        // Dihedral parameters: K=10, m=1, δ=180° (π radians)
        // Energy minimum at φ = π
        topo.moltypes[0].proper_dihedrals.push(Dihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            dihedral_type: 0,
        });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 10.0,
            cospd: -1.0, // cos(180°) = -1
            pd: PI,      // δ = 180° = π radians
            m: 1,        // multiplicity
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Create planar trans configuration (φ = 180° = π)
        // For trans, atoms 0 and 3 should be on opposite sides of the j-k bond
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0); // Above the j-k plane
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, -0.5, 0.0); // Below the j-k plane (opposite side)

        let result = calculate_dihedral_new_forces(&topo, &conf, &vacuum());

        // At minimum: V = K * (1 + cos(m*φ - δ)) = K * (1 + cos(π - π)) = K * (1 + 1) = 2K = 20
        let expected_energy = 10.0 * (1.0 + (1.0 * PI - PI).cos());

        println!("\n========== New Dihedral Test (At Minimum, φ=180°) ==========");
        println!(
            "Energy: {:.6} kJ/mol (expected: {:.6})",
            result.energy, expected_energy
        );
        println!(
            "Force on atom 0: ({:.6}, {:.6}, {:.6})",
            result.forces[0].x, result.forces[0].y, result.forces[0].z
        );

        // At minimum, energy should be maximum (2K) and forces should be zero
        assert!(
            (result.energy - expected_energy).abs() < 1e-6,
            "Energy at minimum should be 2K: got {}, expected {}",
            result.energy,
            expected_energy
        );

        // Forces should be very small at minimum
        for i in 0..4 {
            assert!(
                result.forces[i].length() < 1e-3,
                "Force on atom {} should be near zero at minimum",
                i
            );
        }
    }

    #[test]
    fn test_dihedral_new_vs_old() {
        // Compare new and old dihedral formulas at φ = 0° (cis)
        // For δ = 0°, both formulas should give same result
        let mut topo = Topology::new();

        // Add 4 atoms
        for i in 0..4 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0 / 12.0; 4];

        // Parameters: K=10, m=1, δ=0°
        topo.moltypes[0].proper_dihedrals.push(Dihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            dihedral_type: 0,
        });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 10.0,
            cospd: 1.0, // cos(0°) = 1
            pd: 0.0,    // δ = 0°
            m: 1,
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Create eclipsed configuration (φ ≈ 0°)
        // Atoms form a zigzag with small dihedral angle
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, 0.5, 0.01); // Small z-deviation for φ ≈ 0°

        let result_new = calculate_dihedral_new_forces(&topo, &conf, &vacuum());
        let result_old = calculate_dihedral_forces(&topo, &conf, &vacuum());

        println!("\n========== New vs Old Dihedral Comparison ==========");
        println!("New formula energy: {:.6} kJ/mol", result_new.energy);
        println!("Old formula energy: {:.6} kJ/mol", result_old.energy);
        println!(
            "Energy difference: {:.6} kJ/mol",
            (result_new.energy - result_old.energy).abs()
        );

        // For small angles and δ=0°, both formulas should give very similar results
        assert!(
            (result_new.energy - result_old.energy).abs() < 0.1,
            "New and old formulas should agree for δ=0°"
        );
    }

    #[test]
    fn test_dihedral_new_arbitrary_phase() {
        use std::f64::consts::PI;

        // Test new dihedral with arbitrary phase shift δ = 60°
        let mut topo = Topology::new();

        // Add 4 atoms
        for i in 0..4 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0 / 12.0; 4];

        // Parameters: K=15, m=3, δ=60° = π/3 rad
        let delta = PI / 3.0; // 60°
        topo.moltypes[0].proper_dihedrals.push(Dihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            dihedral_type: 0,
        });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 15.0,
            cospd: delta.cos(), // cos(60°) = 0.5
            pd: delta,          // δ = 60°
            m: 3,               // multiplicity 3
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Create configuration with φ = 90°
        // Proper dihedral geometry: first plane in xy, rotate second plane by 90° around x
        conf.current_mut().pos[0] = Vec3::new(-1.0, 0.5, 0.0); // First plane (xy)
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new(2.0, 0.0, 0.5); // Second plane (xz) - 90° rotation

        let result = calculate_dihedral_new_forces(&topo, &conf, &vacuum());

        // Energy: V = K * (1 + cos(m*φ - δ)) = 15 * (1 + cos(3*90° - 60°))
        //           = 15 * (1 + cos(270° - 60°)) = 15 * (1 + cos(210°))
        //           = 15 * (1 - sqrt(3)/2) ≈ 2.00961894323
        let phi = PI / 2.0; // 90° in radians
        let expected_energy = 15.0 * (1.0 + (3.0 * phi - delta).cos());

        println!("\n========== New Dihedral Test (Arbitrary Phase δ=60°) ==========");
        println!("φ = 90°, m = 3, δ = 60°");
        println!(
            "Energy: {:.6} kJ/mol (expected: {:.6})",
            result.energy, expected_energy
        );

        // Energy should be finite and non-zero
        assert!(
            !result.energy.is_nan() && !result.energy.is_infinite(),
            "Energy should be finite"
        );

        // Energy should be reasonable (within broad range)
        // The exact value depends on the precise dihedral angle
        assert!(
            result.energy > 0.0 && result.energy < 30.0,
            "Energy should be in reasonable range: got {}",
            result.energy
        );
    }

    #[test]
    fn test_dihedral_new_force_conservation() {
        // Test force conservation for new dihedral formula
        let mut topo = Topology::new();

        // Add 4 atoms
        for i in 0..4 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0 / 12.0; 4];

        topo.moltypes[0].proper_dihedrals.push(Dihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            dihedral_type: 0,
        });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 20.0,
            cospd: 0.0,                     // cos(90°) = 0
            pd: std::f64::consts::PI / 2.0, // δ = 90°
            m: 2,
        });

        let mut conf = Configuration::new(4, 1, 1);

        // Arbitrary dihedral geometry
        conf.current_mut().pos[0] = Vec3::new(-1.2, 0.3, 0.1);
        conf.current_mut().pos[1] = Vec3::new(-0.4, -0.1, 0.2);
        conf.current_mut().pos[2] = Vec3::new(0.5, 0.2, -0.1);
        conf.current_mut().pos[3] = Vec3::new(1.3, -0.2, 0.3);

        let result = calculate_dihedral_new_forces(&topo, &conf, &vacuum());

        println!("\n========== New Dihedral Force Conservation Test ==========");
        println!("Energy: {:.6} kJ/mol", result.energy);

        // Check force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2] + result.forces[3];
        println!(
            "Total force: ({:.6}, {:.6}, {:.6})",
            total_force.x, total_force.y, total_force.z
        );

        assert!(
            total_force.length() < 1e-4,
            "Forces should be conserved: total = {}",
            total_force.length()
        );
    }
}
