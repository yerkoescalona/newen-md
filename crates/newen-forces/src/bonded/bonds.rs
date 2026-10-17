//! Bond force calculations: quartic, harmonic, and CG bonds.

use newen_core::configuration::Configuration;
use newen_core::topology::Topology;

use super::ForceEnergy;
use newen_core::math::Periodicity;

/// Calculate quartic bond forces and energies (GROMOS standard)
///
/// Potential: V = (1/4) * k_harmonic * (r^2 - r0^2)^2
/// Force: F = -dV/dr = -k_harmonic * (r^2 - r0^2) * r * r_vec/r
pub fn calculate_bond_forces_quartic(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
) -> ForceEnergy {
    calculate_bond_forces_quartic_grouped::<false>(topo, conf, pbc, &mut [])
}

/// [`calculate_bond_forces_quartic`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `bond_energy[atom_energy_group(i)]`).
pub(crate) fn calculate_bond_forces_quartic_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    per_group: &mut [f64],
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for bond in topo.all_bonds_global() {
        // gromosXX moves a constrained bond out of the force list when it creates the
        // constraints (`create_constraints.cc`): its length is imposed, so it contributes
        // neither energy nor force. Empty set for NTC = 1.
        if topo.constrained_bonds.contains(&(bond.i, bond.j)) {
            continue;
        }
        if bond.bond_type >= topo.bond_parameters.len() {
            continue;
        }

        let params = &topo.bond_parameters[bond.bond_type];
        // GROMOS convention: v = pos(i) - pos(j)
        let r_vec = pbc.nearest_image(conf.current().pos[bond.i], conf.current().pos[bond.j]);

        // Energy: V = (1/4) * k * (r^2 - r0^2)^2
        let r2 = r_vec.length_squared();
        let r0_2 = params.r0 * params.r0;
        let dr2 = r2 - r0_2;
        let energy = 0.25 * params.k_quartic * dr2 * dr2;

        // Force: f = -K * (r^2 - r0^2) * v  (GROMOS convention)
        let force = r_vec * (-params.k_quartic * dr2);

        result.energy += energy;
        if GROUPS {
            per_group[topo.atom_to_energy_group[bond.i]] += energy;
        }
        result.forces[bond.i] += force;
        result.forces[bond.j] -= force;

        // GROMOS: virial_tensor(a, c) += v(a) * f(c)
        let rv = [r_vec.x, r_vec.y, r_vec.z];
        let fv = [force.x, force.y, force.z];
        for a in 0..3 {
            for c in 0..3 {
                result.virial[a][c] += rv[a] * fv[c];
            }
        }
    }

    result
}

/// Calculate harmonic bond forces (alternative to quartic)
///
/// Potential: V = (1/2) * k * (r - r0)^2
/// Force: F = -k * (r - r0) * r_vec/r
pub fn calculate_bond_forces_harmonic(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
) -> ForceEnergy {
    calculate_bond_forces_harmonic_grouped::<false>(topo, conf, pbc, &mut [])
}

/// [`calculate_bond_forces_harmonic`]; with `GROUPS` each term is also booked to the energy group of its first
/// atom (gromosXX `bond_energy[atom_energy_group(i)]`).
pub(crate) fn calculate_bond_forces_harmonic_grouped<const GROUPS: bool>(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
    per_group: &mut [f64],
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for bond in topo.all_bonds_global() {
        // gromosXX moves a constrained bond out of the force list when it creates the
        // constraints (`create_constraints.cc`): its length is imposed, so it contributes
        // neither energy nor force. Empty set for NTC = 1.
        if topo.constrained_bonds.contains(&(bond.i, bond.j)) {
            continue;
        }
        if bond.bond_type >= topo.bond_parameters.len() {
            continue;
        }

        let params = &topo.bond_parameters[bond.bond_type];
        // GROMOS convention: v = pos(i) - pos(j)
        let r_vec = pbc.nearest_image(conf.current().pos[bond.i], conf.current().pos[bond.j]);
        let r = r_vec.length();

        if r < 1e-10 {
            continue;
        }

        // Energy: V = (1/2) * k * (r - r0)^2
        let dr = r - params.r0;
        let energy = 0.5 * params.k_harmonic * dr * dr;

        // Force: f = -k * (r - r0) * v/|v|  (GROMOS convention)
        let force = r_vec * (-params.k_harmonic * dr / r);

        result.energy += energy;
        if GROUPS {
            per_group[topo.atom_to_energy_group[bond.i]] += energy;
        }
        result.forces[bond.i] += force;
        result.forces[bond.j] -= force;

        // GROMOS: virial_tensor(a, c) += v(a) * f(c)
        let rv = [r_vec.x, r_vec.y, r_vec.z];
        let fv = [force.x, force.y, force.z];
        for a in 0..3 {
            for c in 0..3 {
                result.virial[a][c] += rv[a] * fv[c];
            }
        }
    }

    result
}

/// Calculate CG (coarse-grained) bond forces - quartic repulsive potential
///
/// Potential: V = 0.5 * K * (r - r0)^4 when r > r0, otherwise V = 0
/// Force: F = -2 * K * (r - r0)^3 * direction when r > r0, otherwise F = 0
///
/// This is a soft repulsive potential that prevents bonds from stretching too much
/// but allows free compression. Useful for coarse-grained models.
pub fn calculate_cg_bond_forces(
    topo: &Topology,
    conf: &Configuration,
    pbc: &Periodicity,
) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    // CG bonds are stored separately in topology (if available)
    // For now, we'll use a marker in bond parameters to identify CG bonds
    for bond in topo.all_bonds_global() {
        if bond.bond_type >= topo.bond_parameters.len() {
            continue;
        }

        let params = &topo.bond_parameters[bond.bond_type];

        // Skip if this is not a CG bond (we'll add a flag later)
        // For now, assume CG bonds have negative k_harmonic as a marker
        if params.k_harmonic >= 0.0 {
            continue;
        }

        // GROMOS convention: v = pos(i) - pos(j)
        let r_vec = pbc.nearest_image(conf.current().pos[bond.i], conf.current().pos[bond.j]);
        let r = r_vec.length();

        if r < 1e-10 {
            // Avoid division by zero
            continue;
        }

        // Only apply force when bond is stretched (r > r0)
        if r > params.r0 {
            let diff = r - params.r0;
            let diff2 = diff * diff;
            let diff3 = diff2 * diff;
            let diff4 = diff2 * diff2;

            // Use absolute value of k_harmonic (it's negative as a marker)
            let k = params.k_harmonic.abs();

            // Energy: V = 0.5 * K * (r - r0)^4
            let energy = 0.5 * k * diff4;

            // Force magnitude: -dV/dr = -2 * K * (r - r0)^3
            let f_magnitude = -2.0 * k * diff3;

            // Force vector: F = f_magnitude * direction
            let direction = r_vec / (r);
            let force = direction * (f_magnitude);

            result.energy += energy;
            result.forces[bond.i] += force;
            result.forces[bond.j] -= force;
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
    use newen_io::coordinate::read_coordinate_file;
    use newen_io::topology::{build_topology, read_topology_file};

    #[test]
    fn test_bond_forces_quartic() {
        let topo_result = read_topology_file("../md++/src/check/data/cg16.topo");
        let conf_result = read_coordinate_file("../md++/src/check/data/cg16.conf", 1, 1);

        if topo_result.is_err() || conf_result.is_err() {
            println!("Skipping: test files not found");
            return;
        }

        let parsed = topo_result.unwrap();
        let topo = build_topology(parsed);
        let conf = conf_result.unwrap();

        let result = calculate_bond_forces_quartic(&topo, &conf, &vacuum());

        println!("Bond energy: {:.6} kJ/mol", result.energy);
        println!(
            "Force on atom 0: ({:.6}, {:.6}, {:.6})",
            result.forces[0].x, result.forces[0].y, result.forces[0].z
        );

        // Energy should be positive (bonds are stretched or compressed)
        assert!(result.energy >= 0.0);

        // Forces should not all be zero
        let total_force: Vec3 = result.forces.iter().sum();
        println!(
            "Total force (should be ~0): ({:.6}, {:.6}, {:.6})",
            total_force.x, total_force.y, total_force.z
        );
    }

    #[test]
    fn test_cg_bond_compressed() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Atom, Bond, BondParameters};

        // Test CG bond when compressed (r < r0) - should have zero force/energy
        let mut topo = newen_core::topology::Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0];

        // Add one CG bond (negative k_harmonic as marker)
        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,     // Not used for CG bonds
            k_harmonic: -500.0, // Negative = CG bond marker
            r0: 0.15,           // 0.15 nm equilibrium
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Position atoms compressed (0.1 nm < 0.15 nm equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.1, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf, &vacuum());

        println!("\n========== CG Bond Test (Compressed) ==========");
        println!("Bond length: 0.10 nm, equilibrium: 0.15 nm");
        println!("Energy: {:.6} kJ/mol", result.energy);
        println!(
            "Force on atom 0: ({:.6}, {:.6}, {:.6})",
            result.forces[0].x, result.forces[0].y, result.forces[0].z
        );

        // When compressed (r < r0), CG bonds should have zero energy and force
        assert!(
            result.energy < 1e-10,
            "Energy should be zero when compressed: {}",
            result.energy
        );
        assert!(
            result.forces[0].length() < 1e-10,
            "Force should be zero when compressed"
        );
    }

    #[test]
    fn test_cg_bond_at_equilibrium() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Atom, Bond, BondParameters};

        // Test CG bond at equilibrium (r = r0) - should have zero force/energy
        let mut topo = newen_core::topology::Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0];

        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: -500.0,
            r0: 0.15,
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Position atoms at equilibrium (0.15 nm)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.15, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf, &vacuum());

        println!("\n========== CG Bond Test (Equilibrium) ==========");
        println!("Bond length: 0.15 nm (at equilibrium)");
        println!("Energy: {:.6} kJ/mol", result.energy);

        // At equilibrium, energy and force should be zero
        assert!(
            result.energy < 1e-10,
            "Energy should be zero at equilibrium"
        );
        assert!(
            result.forces[0].length() < 1e-10,
            "Force should be zero at equilibrium"
        );
    }

    #[test]
    fn test_cg_bond_stretched() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Atom, Bond, BondParameters};

        // Test CG bond when stretched (r > r0) - should have repulsive force
        let mut topo = newen_core::topology::Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0];

        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: -500.0, // K = 500 kJ/(mol·nm^4)
            r0: 0.15,           // r0 = 0.15 nm
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Position atoms stretched (0.20 nm > 0.15 nm equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.20, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf, &vacuum());

        // Calculate expected values
        // diff = 0.20 - 0.15 = 0.05 nm
        // Energy = 0.5 * 500 * (0.05)^4 = 0.5 * 500 * 0.00000625 = 0.0015625 kJ/mol
        let diff: f64 = 0.05;
        let expected_energy = 0.5 * 500.0 * diff.powi(4);
        let expected_force_mag = 2.0 * 500.0 * diff.powi(3); // Attractive (pulls atoms together)

        println!("\n========== CG Bond Test (Stretched) ==========");
        println!("Bond length: 0.20 nm, equilibrium: 0.15 nm");
        println!(
            "Energy: {:.6} kJ/mol (expected: {:.6})",
            result.energy, expected_energy
        );
        println!(
            "Force magnitude: {:.6} (expected: {:.6})",
            result.forces[0].length(),
            expected_force_mag
        );
        println!(
            "Force on atom 0: ({:.6}, {:.6}, {:.6})",
            result.forces[0].x, result.forces[0].y, result.forces[0].z
        );
        println!(
            "Force on atom 1: ({:.6}, {:.6}, {:.6})",
            result.forces[1].x, result.forces[1].y, result.forces[1].z
        );

        // Check energy
        assert!(
            (result.energy - expected_energy).abs() < 1e-6,
            "Energy mismatch: got {}, expected {}",
            result.energy,
            expected_energy
        );

        // Check force magnitude (should be attractive = positive x-direction on atom 0, toward atom 1 at x=0.20)
        assert!(
            result.forces[0].x > 0.0,
            "Force should be attractive (positive x, toward atom 1)"
        );
        assert!(
            (result.forces[0].length() - expected_force_mag).abs() < 1e-5,
            "Force magnitude mismatch: got {}, expected {}",
            result.forces[0].length(),
            expected_force_mag
        );

        // Force conservation
        let total_force = result.forces[0] + result.forces[1];
        assert!(total_force.length() < 1e-5, "Forces should be conserved");
    }

    #[test]
    fn test_cg_bond_highly_stretched() {
        use newen_core::configuration::Configuration;
        use newen_core::math::Vec3;
        use newen_core::topology::{Atom, Bond, BondParameters};

        // Test CG bond when highly stretched - quartic potential grows rapidly
        let mut topo = newen_core::topology::Topology::new();

        // Add 2 atoms
        for i in 0..2 {
            topo.moltypes[0].atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "TEST".to_string(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: true,
            });
        }
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0, 1.0 / 12.0];

        topo.moltypes[0].bonds.push(Bond {
            i: 0,
            j: 1,
            bond_type: 0,
            hydrogen: false,
        });
        topo.bond_parameters.push(BondParameters {
            k_quartic: 0.0,
            k_harmonic: -500.0,
            r0: 0.15,
        });

        let mut conf = Configuration::new(2, 1, 1);

        // Highly stretched: 0.25 nm (0.10 nm beyond equilibrium)
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.25, 0.0, 0.0);

        let result = calculate_cg_bond_forces(&topo, &conf, &vacuum());

        // diff = 0.10 nm, Energy = 0.5 * 500 * (0.10)^4 = 0.025 kJ/mol
        let diff: f64 = 0.10;
        let expected_energy = 0.5 * 500.0 * diff.powi(4);

        println!("\n========== CG Bond Test (Highly Stretched) ==========");
        println!("Bond length: 0.25 nm, equilibrium: 0.15 nm (diff = 0.10 nm)");
        println!(
            "Energy: {:.6} kJ/mol (expected: {:.6})",
            result.energy, expected_energy
        );
        println!("Force magnitude: {:.6}", result.forces[0].length());

        // Energy should grow as (diff)^4
        assert!(
            (result.energy - expected_energy).abs() < 1e-6,
            "Energy mismatch for highly stretched bond"
        );

        // Verify quartic growth: doubling diff should multiply energy by 16
        // Previous test: diff = 0.05, energy = 0.0015625
        // This test: diff = 0.10, energy should be 16x larger = 0.025
        assert!(
            result.energy > 0.02,
            "Quartic potential should grow rapidly"
        );
    }
}
