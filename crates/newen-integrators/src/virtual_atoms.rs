//! Virtual Atoms for Coarse-Grained and United-Atom Models
//!
//! Direct translation of md++/src/algorithm/virtualatoms/
//!
//! **Purpose**: Handle virtual (dummy) atoms whose positions are constructed
//! geometrically from real atoms
//!
//! Virtual atoms are used for:
//! - **Coarse-grained models** (Martini, SIRAH): Virtual sites for interaction centers
//! - **United-atom models**: Explicit hydrogens on heavy atoms
//! - **Rigid groups**: TIP4P/TIP5P water models with virtual sites
//! - **Extended force fields**: CHARMM Drude polarizable models
//!
//! **Algorithm**:
//! 1. Before force calculation: Construct virtual atom positions from parent atoms
//! 2. After force calculation: Redistribute forces from virtual atoms to parents
//!
//! **Virtual atom types** (GROMOS):
//! - Type 1: COM (center of mass) of N atoms
//! - Type 2: On line between two atoms
//! - Type 3: In plane of three atoms
//! - Type 4: Out of plane (tetrahedral, e.g., aliphatic H)

use newen_core::configuration::Configuration;
use newen_core::math::Vec3;

/// Virtual atom definition
///
/// Defines how a virtual atom's position is constructed from parent atoms
#[derive(Debug, Clone)]
pub struct VirtualAtom {
    /// Index of the virtual atom
    pub atom_index: usize,
    /// Type of virtual atom (1=COM, 2=linear, 3=planar, 4=out-of-plane)
    pub virt_type: i32,
    /// Parent atom indices (1-4 atoms depending on type)
    pub parent_atoms: Vec<usize>,
    /// Geometric parameters (distances, angles)
    pub parameters: Vec<f64>,
    /// Masses for COM calculation (type 1 only)
    pub masses: Vec<f64>,
}

impl VirtualAtom {
    /// Create COM virtual atom (type 1)
    ///
    /// Position: r_virt = Σ(mᵢ * rᵢ) / Σmᵢ
    pub fn com(atom_index: usize, parents: Vec<usize>, masses: Vec<f64>) -> Self {
        assert_eq!(
            parents.len(),
            masses.len(),
            "Number of parents must match masses"
        );
        Self {
            atom_index,
            virt_type: 1,
            parent_atoms: parents,
            parameters: Vec::new(),
            masses,
        }
    }

    /// Create linear virtual atom (type 2)
    ///
    /// Position: r_virt = r1 + a * (r2 - r1)
    /// Parameter a: distance fraction along bond
    pub fn linear(atom_index: usize, parent1: usize, parent2: usize, fraction: f64) -> Self {
        Self {
            atom_index,
            virt_type: 2,
            parent_atoms: vec![parent1, parent2],
            parameters: vec![fraction],
            masses: Vec::new(),
        }
    }

    /// Create planar virtual atom (type 3)
    ///
    /// Position in plane of three atoms
    /// r_virt = r1 + a*(r2-r1) + b*(r3-r1)
    pub fn planar(
        atom_index: usize,
        parent1: usize,
        parent2: usize,
        parent3: usize,
        a: f64,
        b: f64,
    ) -> Self {
        Self {
            atom_index,
            virt_type: 3,
            parent_atoms: vec![parent1, parent2, parent3],
            parameters: vec![a, b],
            masses: Vec::new(),
        }
    }

    /// Create out-of-plane virtual atom (type 4)
    ///
    /// Position: tetrahedral geometry (e.g., aliphatic H)
    /// Uses four parent atoms to define 3D position
    pub fn out_of_plane(
        atom_index: usize,
        parent1: usize,
        parent2: usize,
        parent3: usize,
        parent4: usize,
        params: Vec<f64>,
    ) -> Self {
        Self {
            atom_index,
            virt_type: 4,
            parent_atoms: vec![parent1, parent2, parent3, parent4],
            parameters: params,
            masses: Vec::new(),
        }
    }

    /// Construct virtual atom position from parent atoms
    ///
    /// # Algorithm
    /// Depends on virtual atom type:
    /// - Type 1 (COM): Weighted average of parent positions
    /// - Type 2 (Linear): Linear interpolation between two atoms
    /// - Type 3 (Planar): Linear combination in plane
    /// - Type 4 (Out-of-plane): 3D geometric construction
    pub fn construct_position(&self, conf: &Configuration) -> Vec3 {
        match self.virt_type {
            1 => self.construct_com(conf),
            2 => self.construct_linear(conf),
            3 => self.construct_planar(conf),
            4 => self.construct_out_of_plane(conf),
            _ => panic!("Unknown virtual atom type: {}", self.virt_type),
        }
    }

    /// Type 1: Center of mass
    fn construct_com(&self, conf: &Configuration) -> Vec3 {
        let mut total_mass = 0.0;
        let mut com = Vec3::ZERO;

        for (i, &parent_idx) in self.parent_atoms.iter().enumerate() {
            let mass = self.masses[i];
            total_mass += mass;
            com += conf.current().pos[parent_idx] * (mass);
        }

        com / (total_mass)
    }

    /// Type 2: Linear interpolation
    fn construct_linear(&self, conf: &Configuration) -> Vec3 {
        let r1 = conf.current().pos[self.parent_atoms[0]];
        let r2 = conf.current().pos[self.parent_atoms[1]];
        let a = self.parameters[0];

        r1 + (r2 - r1) * a
    }

    /// Type 3: Planar
    fn construct_planar(&self, conf: &Configuration) -> Vec3 {
        let r1 = conf.current().pos[self.parent_atoms[0]];
        let r2 = conf.current().pos[self.parent_atoms[1]];
        let r3 = conf.current().pos[self.parent_atoms[2]];
        let a = self.parameters[0];
        let b = self.parameters[1];

        r1 + (r2 - r1) * a + (r3 - r1) * b
    }

    /// Type 4: Out of plane (tetrahedral)
    fn construct_out_of_plane(&self, conf: &Configuration) -> Vec3 {
        // TODO: Implement full tetrahedral geometry
        // For now, return parent1 position as placeholder
        conf.current().pos[self.parent_atoms[0]]
    }

    /// Redistribute force from virtual atom to parent atoms
    ///
    /// # Algorithm
    /// Force redistribution follows chain rule:
    /// F_parent_i = -∂V/∂r_parent_i = -∂V/∂r_virt * ∂r_virt/∂r_parent_i
    ///
    /// The derivative ∂r_virt/∂r_parent depends on virtual atom type
    pub fn redistribute_force(&self, force_on_virt: Vec3, conf: &mut Configuration) {
        match self.virt_type {
            1 => self.redistribute_com(force_on_virt, conf),
            2 => self.redistribute_linear(force_on_virt, conf),
            3 => self.redistribute_planar(force_on_virt, conf),
            4 => self.redistribute_out_of_plane(force_on_virt, conf),
            _ => panic!("Unknown virtual atom type: {}", self.virt_type),
        }
    }

    /// Type 1: Redistribute COM force
    fn redistribute_com(&self, force_on_virt: Vec3, conf: &mut Configuration) {
        let total_mass: f64 = self.masses.iter().sum();

        for (i, &parent_idx) in self.parent_atoms.iter().enumerate() {
            let weight = self.masses[i] / total_mass;
            conf.current_mut().force[parent_idx] += force_on_virt * weight;
        }
    }

    /// Type 2: Redistribute linear force
    fn redistribute_linear(&self, force_on_virt: Vec3, conf: &mut Configuration) {
        let a = self.parameters[0];
        conf.current_mut().force[self.parent_atoms[0]] += force_on_virt * (1.0 - a);
        conf.current_mut().force[self.parent_atoms[1]] += force_on_virt * a;
    }

    /// Type 3: Redistribute planar force
    fn redistribute_planar(&self, force_on_virt: Vec3, conf: &mut Configuration) {
        let a = self.parameters[0];
        let b = self.parameters[1];

        conf.current_mut().force[self.parent_atoms[0]] += force_on_virt * (1.0 - a - b);
        conf.current_mut().force[self.parent_atoms[1]] += force_on_virt * a;
        conf.current_mut().force[self.parent_atoms[2]] += force_on_virt * b;
    }

    /// Type 4: Redistribute out-of-plane force
    fn redistribute_out_of_plane(&self, force_on_virt: Vec3, conf: &mut Configuration) {
        // TODO: Implement full tetrahedral force redistribution
        // For now, apply all force to first parent
        conf.current_mut().force[self.parent_atoms[0]] += force_on_virt;
    }
}

/// Virtual atoms manager
///
/// Manages all virtual atoms in the system
#[derive(Debug, Clone)]
pub struct VirtualAtomsManager {
    /// List of all virtual atoms
    pub virtual_atoms: Vec<VirtualAtom>,
    /// Enable/disable virtual atoms
    pub enabled: bool,
}

impl VirtualAtomsManager {
    pub fn new() -> Self {
        Self {
            virtual_atoms: Vec::new(),
            enabled: true,
        }
    }

    pub fn add_virtual_atom(&mut self, virt_atom: VirtualAtom) {
        self.virtual_atoms.push(virt_atom);
    }

    /// Construct all virtual atom positions before force calculation
    ///
    /// Call this BEFORE calculating forces
    pub fn construct_positions(&self, conf: &mut Configuration) {
        if !self.enabled {
            return;
        }

        for virt in &self.virtual_atoms {
            let pos = virt.construct_position(conf);
            conf.current_mut().pos[virt.atom_index] = pos;
        }
    }

    /// Redistribute forces from virtual atoms to parent atoms
    ///
    /// Call this AFTER calculating forces
    pub fn redistribute_forces(&self, conf: &mut Configuration) {
        if !self.enabled {
            return;
        }

        for virt in &self.virtual_atoms {
            let force = conf.current().force[virt.atom_index];
            // Zero the force on virtual atom (it's been redistributed)
            conf.current_mut().force[virt.atom_index] = Vec3::ZERO;
            // Redistribute to parents
            virt.redistribute_force(force, conf);
        }
    }

    /// Construct positions and velocities for virtual atoms
    ///
    /// Velocities are derived from parent atom velocities
    pub fn construct_velocities(&self, conf: &mut Configuration) {
        if !self.enabled {
            return;
        }

        for virt in &self.virtual_atoms {
            // Velocity construction follows same rules as position
            // but applied to velocity field instead
            let vel = match virt.virt_type {
                1 => self.construct_com_velocity(virt, conf),
                2 => self.construct_linear_velocity(virt, conf),
                3 => self.construct_planar_velocity(virt, conf),
                4 => self.construct_out_of_plane_velocity(virt, conf),
                _ => Vec3::ZERO,
            };
            conf.current_mut().vel[virt.atom_index] = vel;
        }
    }

    fn construct_com_velocity(&self, virt: &VirtualAtom, conf: &Configuration) -> Vec3 {
        let mut total_mass = 0.0;
        let mut com_vel = Vec3::ZERO;

        for (i, &parent_idx) in virt.parent_atoms.iter().enumerate() {
            let mass = virt.masses[i];
            total_mass += mass;
            com_vel += conf.current().vel[parent_idx] * (mass);
        }

        com_vel / (total_mass)
    }

    fn construct_linear_velocity(&self, virt: &VirtualAtom, conf: &Configuration) -> Vec3 {
        let v1 = conf.current().vel[virt.parent_atoms[0]];
        let v2 = conf.current().vel[virt.parent_atoms[1]];
        let a = virt.parameters[0];

        v1 * (1.0 - a) + v2 * a
    }

    fn construct_planar_velocity(&self, virt: &VirtualAtom, conf: &Configuration) -> Vec3 {
        let v1 = conf.current().vel[virt.parent_atoms[0]];
        let v2 = conf.current().vel[virt.parent_atoms[1]];
        let v3 = conf.current().vel[virt.parent_atoms[2]];
        let a = virt.parameters[0];
        let b = virt.parameters[1];

        v1 * (1.0 - a - b) + v2 * a + v3 * b
    }

    fn construct_out_of_plane_velocity(&self, virt: &VirtualAtom, conf: &Configuration) -> Vec3 {
        // TODO: Implement full tetrahedral velocity
        conf.current().vel[virt.parent_atoms[0]]
    }
}

impl Default for VirtualAtomsManager {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_com_virtual_atom() {
        let mut conf = Configuration::new(4, 1, 1);

        // Set up three real atoms forming a triangle
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new(0.0, 1.0, 0.0);

        // Virtual atom at COM (equal masses)
        let virt = VirtualAtom::com(3, vec![0, 1, 2], vec![1.0, 1.0, 1.0]);
        let pos = virt.construct_position(&conf);

        // COM should be at (1/3, 1/3, 0)
        assert!((pos.x - 1.0 / 3.0).abs() < 1e-6);
        assert!((pos.y - 1.0 / 3.0).abs() < 1e-6);
        assert!(pos.z.abs() < 1e-6);
    }

    #[test]
    fn test_linear_virtual_atom() {
        let mut conf = Configuration::new(3, 1, 1);

        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);

        // Virtual atom at midpoint (fraction = 0.5)
        let virt = VirtualAtom::linear(2, 0, 1, 0.5);
        let pos = virt.construct_position(&conf);

        // Should be at (0.5, 0, 0)
        assert!((pos.x - 0.5).abs() < 1e-6);
        assert!(pos.y.abs() < 1e-6);
        assert!(pos.z.abs() < 1e-6);
    }

    #[test]
    fn test_force_redistribution() {
        let mut conf = Configuration::new(3, 1, 1);

        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(1.0, 0.0, 0.0);

        // Linear virtual atom
        let virt = VirtualAtom::linear(2, 0, 1, 0.5);

        // Apply force of 10 N in x-direction to virtual atom
        let force = Vec3::new(10.0, 0.0, 0.0);
        virt.redistribute_force(force, &mut conf);

        // Each parent should get 5 N (50% each)
        assert!((conf.current().force[0].x - 5.0).abs() < 1e-6);
        assert!((conf.current().force[1].x - 5.0).abs() < 1e-6);
    }
}
