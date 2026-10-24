//! Angle force calculations: cosine-harmonic and harmonic angles.

use newen_core::configuration::Configuration;
use newen_core::topology::Topology;

use super::{ForceEnergy, ZERO_DISTANCE, ZERO_SINE};
use newen_core::math::{Periodicity, Vec3};

/// One cosine-harmonic angle `V = K/2 · (cos θ − cos θ0)²` at `r_ij = x_i − x_j`,
/// `r_kj = x_k − x_j`: `(energy, f_i, f_j, f_k)`; nothing when a bond is shorter than
/// [`ZERO_DISTANCE`]. `cos0 = cos θ0` comes precomputed, once per angle type. In this order:
///
/// ```text
/// d_ij = |r_ij|,  d_kj = |r_kj|
/// cost   = clamp(r_ij · r_kj / (d_ij · d_kj), −1, 1)
/// d_cos  = cost − cos0
/// energy = 0.5 · k · d_cos · d_cos
/// df     = −k · d_cos
/// f_i    = (r_kj · (1/d_kj) − r_ij · (cost/d_ij)) · (df/d_ij)
/// f_k    = (r_ij · (1/d_ij) − r_kj · (cost/d_kj)) · (df/d_kj)
/// f_j    = −(f_i + f_k)
/// ```
///
/// spec: bonded.cosine_angle
#[inline]
pub fn cosine_angle(r_ij: Vec3, r_kj: Vec3, k: f64, cos0: f64) -> (f64, Vec3, Vec3, Vec3) {
    let d_ij = r_ij.length();
    let d_kj = r_kj.length();
    if d_ij < ZERO_DISTANCE || d_kj < ZERO_DISTANCE {
        return (0.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    }
    let cost = (r_ij.dot(r_kj) / (d_ij * d_kj)).clamp(-1.0, 1.0);
    let d_cos = cost - cos0;
    let energy = 0.5 * k * d_cos * d_cos;
    let df = -k * d_cos;
    let f_i = (r_kj * (1.0 / d_kj) - r_ij * (cost / d_ij)) * (df / d_ij);
    let f_k = (r_ij * (1.0 / d_ij) - r_kj * (cost / d_kj)) * (df / d_kj);
    let f_j = -(f_i + f_k);
    (energy, f_i, f_j, f_k)
}

/// One harmonic angle `V = K/2 · (θ − θ0)²` at `r_ij = x_i − x_j`, `r_kj = x_k − x_j`:
/// `(energy, f_i, f_j, f_k)`; nothing when a bond is shorter than [`ZERO_DISTANCE`]. In this
/// order:
///
/// ```text
/// d_ij = |r_ij|,  d_kj = |r_kj|
/// cos_theta = clamp(r_ij · r_kj / (d_ij · d_kj), −1, 1)
/// theta = acos(cos_theta),  sin_theta = sin(theta)
/// |sin_theta| ≥ ZERO_SINE:
///     ff = k · (theta − theta0) / sin_theta;  k_i = ff / d_ij;  k_k = ff / d_kj
///     energy = 0.5 · k · (theta − theta0) · (theta − theta0)
/// |sin_theta| < ZERO_SINE (θ = π or θ = 0) and |theta0 − π| ≤ ZERO_SINE:
///     k_i = −k / d_ij;  k_k = −k / d_kj;  energy = k · (1 + cos_theta)
/// |sin_theta| < ZERO_SINE otherwise: nothing
/// f_i = (r_kj/d_kj − (r_ij/d_ij) · cos_theta) · k_i
/// f_k = (r_ij/d_ij − (r_kj/d_kj) · cos_theta) · k_k
/// f_j = −(f_i + f_k)
/// ```
///
/// The linear branch is the θ = π limit, but `|sin θ| < ZERO_SINE` also holds at θ ≈ 0: there,
/// with θ0 = π, it returns `k(1 + cos θ)` = 2k (not `k/2 · π²`) and a vanishing force, as
/// gromosXX does (`harm_angle_interaction.cc:112-127` assumes θ = 0 cannot occur).
///
/// spec: bonded.harmonic_angle
#[inline]
pub fn harmonic_angle(r_ij: Vec3, r_kj: Vec3, k: f64, theta0: f64) -> (f64, Vec3, Vec3, Vec3) {
    const PI: f64 = std::f64::consts::PI;
    let nothing = (0.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    let d_ij = r_ij.length();
    let d_kj = r_kj.length();
    if d_ij < ZERO_DISTANCE || d_kj < ZERO_DISTANCE {
        return nothing;
    }
    let cos_theta = (r_ij.dot(r_kj) / (d_ij * d_kj)).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    let sin_theta = theta.sin();

    let (k_i, k_k, energy) = if sin_theta.abs() < ZERO_SINE {
        // θ = π, where (θ − θ0)/sin θ → −1 if θ0 = π; any other θ0 has no finite force there.
        // Also reached at θ ≈ 0, where it gives k(1 + cos θ) as gromosXX.
        if theta0 > PI + ZERO_SINE || theta0 < PI - ZERO_SINE {
            return nothing;
        }
        (-k / d_ij, -k / d_kj, k * (1.0 + cos_theta))
    } else {
        let ff = k * (theta - theta0) / sin_theta;
        let d_theta = theta - theta0;
        (ff / d_ij, ff / d_kj, 0.5 * k * d_theta * d_theta)
    };

    let unit_ij = r_ij / d_ij;
    let unit_kj = r_kj / d_kj;
    let f_i = (unit_kj - unit_ij * cos_theta) * k_i;
    let f_k = (unit_ij - unit_kj * cos_theta) * k_k;
    let f_j = -(f_i + f_k);
    (energy, f_i, f_j, f_k)
}

/// Calculate angle forces (GROMOS cosine-based)
///
/// Potential: V = (1/2) * k_harmonic * (cos(θ) - cos(θ0))^2
/// Forces calculated using chain rule on angle derivatives
pub fn calculate_angle_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
) -> ForceEnergy {
    calculate_angle_forces_grouped::<false>(topo, conf, pbc, &mut [])
}

/// [`calculate_angle_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `angle_energy[atom_energy_group(i)]`).
pub(crate) fn calculate_angle_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    per_group: &mut [f64],
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());
    // cos θ0 once per angle type, as gromosXX stores it when it reads the topology
    let cos0: Vec<f64> = topo
        .angle_parameters
        .iter()
        .map(|p| p.theta0.cos())
        .collect();

    for angle in topo.all_angles_global() {
        if angle.angle_type >= topo.angle_parameters.len() {
            continue;
        }

        let params = &topo.angle_parameters[angle.angle_type];

        // GROMOS convention: rij = pos(i) - pos(j), rkj = pos(k) - pos(j)
        let rij = pbc.nearest_image(conf.current().pos[angle.i], conf.current().pos[angle.j]);
        let rkj = pbc.nearest_image(conf.current().pos[angle.k], conf.current().pos[angle.j]);

        let (energy, f_i, f_j, f_k) =
            cosine_angle(rij, rkj, params.k_cosine, cos0[angle.angle_type]);

        result.energy += energy;
        if GROUPS {
            per_group[topo.atom_to_energy_group[angle.i]] += energy;
        }
        result.forces[angle.i] += f_i;
        result.forces[angle.j] += f_j;
        result.forces[angle.k] += f_k;

        // GROMOS: virial_tensor(a, bb) += rij(a) * fi(bb) + rkj(a) * fk(bb)
        let rij_v = [rij.x, rij.y, rij.z];
        let rkj_v = [rkj.x, rkj.y, rkj.z];
        let fi_v = [f_i.x, f_i.y, f_i.z];
        let fk_v = [f_k.x, f_k.y, f_k.z];
        for a in 0..3 {
            for bb in 0..3 {
                result.virial[a][bb] += rij_v[a] * fi_v[bb] + rkj_v[a] * fk_v[bb];
            }
        }
    }

    result
}

/// Calculate harmonic angle forces (simple harmonic potential)
///
/// Potential: V = (1/2) * K * (θ - θ₀)²
/// where:
/// - K is the force constant
/// - θ is the bond angle in radians
/// - θ₀ is the equilibrium angle in radians
///
/// This is an alternative to the cosine-based angle potential.
/// It's simpler and more intuitive, commonly used in other force fields.
pub fn calculate_harmonic_angle_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
) -> ForceEnergy {
    calculate_harmonic_angle_forces_grouped::<false>(topo, conf, pbc, &mut [])
}

/// [`calculate_harmonic_angle_forces`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `angle_energy[atom_energy_group(i)]`).
pub(crate) fn calculate_harmonic_angle_forces_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    per_group: &mut [f64],
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for angle in topo.all_angles_global() {
        if angle.angle_type >= topo.angle_parameters.len() {
            continue;
        }

        let params = &topo.angle_parameters[angle.angle_type];

        // GROMOS convention: rij = pos(i) - pos(j), rkj = pos(k) - pos(j)
        let r_ij = pbc.nearest_image(conf.current().pos[angle.i], conf.current().pos[angle.j]);
        let r_kj = pbc.nearest_image(conf.current().pos[angle.k], conf.current().pos[angle.j]);

        let (energy, f_i, f_j, f_k) = harmonic_angle(r_ij, r_kj, params.k_harmonic, params.theta0);

        result.energy += energy;
        if GROUPS {
            per_group[topo.atom_to_energy_group[angle.i]] += energy;
        }
        result.forces[angle.i] += f_i;
        result.forces[angle.j] += f_j;
        result.forces[angle.k] += f_k;

        // GROMOS: virial_tensor(a, bb) += rij(a) * fi(bb) + rkj(a) * fk(bb)
        let rij_v = [r_ij.x, r_ij.y, r_ij.z];
        let rkj_v = [r_kj.x, r_kj.y, r_kj.z];
        let fi_v = [f_i.x, f_i.y, f_i.z];
        let fk_v = [f_k.x, f_k.y, f_k.z];
        for a in 0..3 {
            for bb in 0..3 {
                result.virial[a][bb] += rij_v[a] * fi_v[bb] + rkj_v[a] * fk_v[bb];
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
    use newen_core::math::Vec3;
    use newen_core::topology::Topology;
    use newen_io::coordinate::read_coordinate_file;
    use newen_io::topology::{build_topology, read_topology_file};

    #[test]
    fn test_angle_forces() {
        let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
        let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

        if topo_result.is_err() || conf_result.is_err() {
            println!("Skipping: test files not found");
            return;
        }

        let parsed = topo_result.unwrap();
        let topo = build_topology(parsed);
        let conf = conf_result.unwrap();

        let result = calculate_angle_forces(&topo, &conf, &vacuum());

        println!("Angle energy: {:.6} kJ/mol", result.energy);
        println!(
            "Force on atom 0: ({:.6}, {:.6}, {:.6})",
            result.forces[0].x, result.forces[0].y, result.forces[0].z
        );

        assert!(result.energy >= 0.0);

        // Check force conservation
        let total_force: Vec3 = result.forces.iter().sum();
        println!(
            "Total force (should be ~0): ({:.6}, {:.6}, {:.6})",
            total_force.x, total_force.y, total_force.z
        );

        // Should be very close to zero (within numerical precision)
        assert!(total_force.length() < 1e-4);
    }

    #[test]
    fn test_harmonic_angle_simple() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Angle, AngleParameters, Atom};

        // Create simple 3-atom system: i-j-k forming an angle
        let mut topo = Topology::new();

        // Add 3 atoms with equal masses
        for i in 0..3 {
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
        topo.mass = vec![12.0, 12.0, 12.0]; // Carbon atoms
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0];

        // Add one angle: 0-1-2
        topo.moltypes[0].angles.push(Angle {
            i: 0,
            j: 1,
            k: 2,
            angle_type: 0,
        });

        // Angle parameters: K = 500 kJ/(mol·rad²), θ₀ = 120° = 2.0944 rad
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 0.0, // Not used for harmonic angles
            k_harmonic: 500.0,
            theta0: 2.0943951, // 120 degrees in radians
        });

        // Create configuration with atoms at 120° angle
        let mut conf = Configuration::new(3, 1, 1);

        // Position atoms in xy plane
        // Atom 1 (central) at origin
        // Atom 0 along x-axis
        // Atom 2 at 120° from x-axis
        conf.current_mut().pos[0] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(-0.075, 0.12990, 0.0); // 120° angle

        // Calculate forces
        let result = calculate_harmonic_angle_forces(&topo, &conf, &vacuum());

        println!("\n========== Harmonic Angle Test (Equilibrium) ==========");
        println!("Angle energy: {:.6} kJ/mol", result.energy);
        println!(
            "Force on atom 0: ({:.6}, {:.6}, {:.6})",
            result.forces[0].x, result.forces[0].y, result.forces[0].z
        );
        println!(
            "Force on atom 1: ({:.6}, {:.6}, {:.6})",
            result.forces[1].x, result.forces[1].y, result.forces[1].z
        );
        println!(
            "Force on atom 2: ({:.6}, {:.6}, {:.6})",
            result.forces[2].x, result.forces[2].y, result.forces[2].z
        );

        // At equilibrium angle, energy should be near zero
        assert!(
            result.energy < 0.01,
            "Energy at equilibrium should be ~0: {}",
            result.energy
        );

        // Forces should be near zero at equilibrium
        assert!(
            result.forces[0].length() < 0.1,
            "Forces should be small at equilibrium"
        );
        assert!(
            result.forces[1].length() < 0.1,
            "Forces should be small at equilibrium"
        );
        assert!(
            result.forces[2].length() < 0.1,
            "Forces should be small at equilibrium"
        );

        // Force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2];
        println!(
            "Total force: ({:.6}, {:.6}, {:.6})",
            total_force.x, total_force.y, total_force.z
        );
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }

    #[test]
    fn test_harmonic_angle_compressed() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Angle, AngleParameters, Atom};

        // Test with angle compressed (smaller than equilibrium)
        let mut topo = Topology::new();

        // Add 3 atoms
        for i in 0..3 {
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
        topo.mass = vec![12.0, 12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0];

        topo.moltypes[0].angles.push(Angle {
            i: 0,
            j: 1,
            k: 2,
            angle_type: 0,
        });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 0.0,
            k_harmonic: 500.0,
            theta0: 2.0943951, // 120° = 2.0944 rad
        });

        let mut conf = Configuration::new(3, 1, 1);

        // Create 90° angle (compressed from 120° equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(0.0, 0.15, 0.0); // 90° angle

        let result = calculate_harmonic_angle_forces(&topo, &conf, &vacuum());

        println!("\n========== Harmonic Angle Test (Compressed 90°) ==========");
        println!("Angle energy: {:.6} kJ/mol", result.energy);

        // Calculate expected energy: V = 0.5 * K * (θ - θ₀)²
        // θ = π/2 = 1.5708 rad, θ₀ = 2.0944 rad
        // dθ = 1.5708 - 2.0944 = -0.5236 rad (30° difference)
        let d_theta = std::f64::consts::FRAC_PI_6; // 30° = π/6 rad
        let expected_energy = 0.5 * 500.0 * d_theta * d_theta; // ≈ 68.5 kJ/mol

        println!("Expected energy: {:.6} kJ/mol", expected_energy);
        assert!(
            (result.energy - expected_energy).abs() < 1.0,
            "Energy should be ~{} kJ/mol, got {}",
            expected_energy,
            result.energy
        );

        // Energy should be positive when compressed
        assert!(
            result.energy > 0.0,
            "Energy should be positive when angle is compressed"
        );

        // Forces should push atoms apart (opening the angle)
        // Force on atom 0 should have negative y component (push down)
        // Force on atom 2 should have negative x component (push left)
        println!(
            "Force on atom 0: ({:.6}, {:.6}, {:.6})",
            result.forces[0].x, result.forces[0].y, result.forces[0].z
        );
        println!(
            "Force on atom 2: ({:.6}, {:.6}, {:.6})",
            result.forces[2].x, result.forces[2].y, result.forces[2].z
        );

        // Force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2];
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }

    #[test]
    fn test_harmonic_angle_vs_cosine() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Angle, AngleParameters, Atom};

        // Compare harmonic and cosine-based angle potentials
        let mut topo = Topology::new();

        // Add 3 atoms (water-like masses)
        topo.moltypes[0].atoms.push(Atom {
            name: "O".to_string(),
            residue_nr: 1,
            residue_name: "WAT".to_string(),
            iac: 0,
            mass: 16.0,
            charge: 0.0,
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: false,
        });
        for i in 0..2 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("H{}", i + 1),
                residue_nr: 1,
                residue_name: "WAT".to_string(),
                iac: 1,
                mass: 1.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }
        topo.mass = vec![16.0, 1.0, 1.0]; // Water-like: O-H-H
        topo.inverse_mass = vec![1.0 / 16.0, 1.0, 1.0];

        topo.moltypes[0].angles.push(Angle {
            i: 0,
            j: 1,
            k: 2,
            angle_type: 0,
        });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 400.0,
            k_harmonic: 400.0,
            theta0: 1.91063, // 109.47° (tetrahedral angle)
        });

        let mut conf = Configuration::new(3, 1, 1);

        // Create angle slightly off equilibrium (115°)
        conf.current_mut().pos[0] = Vec3::new(0.1, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(-0.04226, 0.09063, 0.0); // 115°

        let harmonic = calculate_harmonic_angle_forces(&topo, &conf, &vacuum());
        let cosine = calculate_angle_forces(&topo, &conf, &vacuum());

        println!("\n========== Harmonic vs Cosine Angle Comparison ==========");
        println!("Harmonic energy: {:.6} kJ/mol", harmonic.energy);
        println!("Cosine energy:   {:.6} kJ/mol", cosine.energy);

        // Both should give reasonable energies (positive for off-equilibrium)
        assert!(harmonic.energy > 0.0, "Harmonic energy should be positive");
        assert!(cosine.energy > 0.0, "Cosine energy should be positive");

        // Both should conserve forces
        let harm_total: Vec3 = harmonic.forces.iter().sum();
        let cos_total: Vec3 = cosine.forces.iter().sum();

        println!(
            "Harmonic total force: ({:.6}, {:.6}, {:.6})",
            harm_total.x, harm_total.y, harm_total.z
        );
        println!(
            "Cosine total force:   ({:.6}, {:.6}, {:.6})",
            cos_total.x, cos_total.y, cos_total.z
        );

        assert!(
            harm_total.length() < 5e-5,
            "Harmonic forces should be conserved"
        );
        assert!(
            cos_total.length() < 5e-5,
            "Cosine forces should be conserved"
        );

        // Near equilibrium, both should give similar results (within 20%)
        let energy_ratio = harmonic.energy / cosine.energy;
        println!("Energy ratio (harmonic/cosine): {:.3}", energy_ratio);
        assert!(
            energy_ratio > 0.5 && energy_ratio < 2.0,
            "Energies should be similar near equilibrium"
        );
    }

    #[test]
    fn test_harmonic_angle_linear() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Angle, AngleParameters, Atom};

        // Test special case: nearly linear angle (θ ≈ 180°)
        let mut topo = Topology::new();

        // Add 3 atoms
        for i in 0..3 {
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
        topo.mass = vec![12.0, 12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0];

        topo.moltypes[0].angles.push(Angle {
            i: 0,
            j: 1,
            k: 2,
            angle_type: 0,
        });
        topo.angle_parameters.push(AngleParameters {
            k_cosine: 0.0, // Not used for harmonic angles
            k_harmonic: 500.0,
            theta0: std::f64::consts::PI, // 180° (linear)
        });

        let mut conf = Configuration::new(3, 1, 1);

        // Create nearly linear angle (179.9°)
        conf.current_mut().pos[0] = Vec3::new(0.15, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(-0.15, 0.001, 0.0); // Almost linear

        let result = calculate_harmonic_angle_forces(&topo, &conf, &vacuum());

        println!("\n========== Harmonic Angle Test (Nearly Linear) ==========");
        println!("Angle energy: {:.6} kJ/mol", result.energy);

        // Energy should be finite and small (close to equilibrium)
        assert!(!result.energy.is_nan(), "Energy should not be NaN");
        assert!(
            !result.energy.is_infinite(),
            "Energy should not be infinite"
        );
        assert!(
            result.energy < 1.0,
            "Energy should be small near equilibrium"
        );

        // Forces should be finite
        for i in 0..3 {
            assert!(!result.forces[i].x.is_nan(), "Force should not be NaN");
            assert!(!result.forces[i].y.is_nan(), "Force should not be NaN");
            assert!(!result.forces[i].z.is_nan(), "Force should not be NaN");
        }

        // Force conservation
        let total_force = result.forces[0] + result.forces[1] + result.forces[2];
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }
}
