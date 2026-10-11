//! Configuration module - system state and energy management
//!
//! Direct Rust translation of GROMOS configuration structures from:
//! - md++/src/configuration/configuration.h
//! - md++/src/configuration/energy.h

use crate::math::{Mat3, Vec3};

/// Simulation box representation
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoxType {
    /// No periodic boundary (gas-phase / vacuum simulations).
    Vacuum,
    /// Orthorhombic (rectangular) periodic box.
    Rectangular,
    /// General triclinic periodic box.
    Triclinic,
    /// Truncated-octahedral box stored as a triclinic cell.
    TruncatedOctahedral,
}

/// Simulation box with periodic boundary conditions
#[derive(Debug, Clone)]
pub struct Box {
    /// Periodicity type (vacuum, rectangular, triclinic, truncated-octahedral).
    pub box_type: BoxType,
    /// Box vectors as columns of a 3×3 matrix (nm).
    pub vectors: Mat3,
    /// Inverse of `vectors`; used for wrapping positions into the primary cell.
    pub inv_vectors: Mat3,
}

impl Box {
    /// Create vacuum box (no periodicity)
    pub fn vacuum() -> Self {
        Self {
            box_type: BoxType::Vacuum,
            vectors: Mat3::IDENTITY,
            inv_vectors: Mat3::IDENTITY,
        }
    }

    /// Create rectangular box with dimensions (lx, ly, lz)
    pub fn rectangular(lx: f64, ly: f64, lz: f64) -> Self {
        let vectors = Mat3::from_cols(
            Vec3::new(lx, 0.0, 0.0),
            Vec3::new(0.0, ly, 0.0),
            Vec3::new(0.0, 0.0, lz),
        );

        let inv_vectors = Mat3::from_cols(
            Vec3::new(1.0 / lx, 0.0, 0.0),
            Vec3::new(0.0, 1.0 / ly, 0.0),
            Vec3::new(0.0, 0.0, 1.0 / lz),
        );

        Self {
            box_type: BoxType::Rectangular,
            vectors,
            inv_vectors,
        }
    }

    /// Create triclinic box from box vectors
    pub fn triclinic(vectors: Mat3) -> Self {
        Self {
            box_type: BoxType::Triclinic,
            vectors,
            inv_vectors: vectors.inverse(),
        }
    }

    /// Create a truncated-octahedron box from its (lower-triangular)
    /// triclinic box vectors, as produced by
    /// `math::truncoct_triclinic_box(cubic, true)` from the legacy NTB=-1
    /// "cube edge length" BOX block (P1.4).
    pub fn truncated_octahedral(vectors: Mat3) -> Self {
        Self {
            box_type: BoxType::TruncatedOctahedral,
            vectors,
            inv_vectors: vectors.inverse(),
        }
    }

    /// Get box volume
    pub fn volume(&self) -> f64 {
        // Volume = determinant of box matrix
        self.vectors.determinant()
    }

    /// Get box dimensions (for rectangular box)
    pub fn dimensions(&self) -> Vec3 {
        Vec3::new(
            self.vectors.x_axis.x,
            self.vectors.y_axis.y,
            self.vectors.z_axis.z,
        )
    }
}

/// Energy storage for the system, laid out like gromosXX's `configuration::Energy`
/// (`configuration/energy.h`): totals, per-energy-group arrays, energy-group-pair matrices and
/// per-bath kinetic energies. A `State` carries two of them: the energies and their λ
/// derivatives (`State::perturbed_energy_derivatives`).
///
/// The `*_total` scalars are the engine's totals; the per-group arrays and the matrices are what
/// the `# bonded`, `# nonbonded` and `# special` sections of `ENERGY03` / `FREEENERDERIVS03`
/// print. Neither object is cleared each step: whoever fills a field assigns it every step.
#[derive(Debug, Clone)]
pub struct Energy {
    /// Total kinetic energy (kJ mol⁻¹).
    pub kinetic_total: f64,
    /// Total potential energy (kJ mol⁻¹).
    pub potential_total: f64,

    /// Sum of all covalent bond stretch energies.
    pub bond_total: f64,
    /// Sum of all bond-angle bending energies.
    pub angle_total: f64,
    /// Sum of all proper dihedral torsion energies.
    pub dihedral_total: f64,
    /// Sum of all improper dihedral (planarity / chirality) energies.
    pub improper_total: f64,
    /// Sum of all cross-dihedral coupling energies (GROMOS CMAP-equivalent).
    pub cross_dihedral_total: f64,

    /// Total Lennard-Jones (van der Waals) nonbonded energy.
    pub lj_total: f64,
    /// Total Coulomb reaction-field electrostatic energy.
    pub crf_total: f64,
    /// Total lattice-sum (Ewald / P3M / PME) electrostatic energy.
    pub ls_total: f64,
    /// Real-space contribution to the lattice-sum energy.
    pub ls_realspace_total: f64,
    /// k-space contribution to the lattice-sum energy.
    pub ls_kspace_total: f64,

    /// Sum of special interaction energies (position restraints, etc.).
    pub special_total: f64,
    /// Energy of each additive term of the run (`Term("xtb", …)`, `Term("schnet", …)`), keyed by
    /// the term's registry name (`xtb`, or `xtb:1` when a kind occurs more than once), in plan
    /// order. Their sum is booked in `special_total`; this keeps them visible one by one
    /// (3.9 G10) so a compensating error between two terms cannot hide.
    pub term_energies: Vec<(String, f64)>,
    /// Sum of position-restraint energies (also part of `special_total`).
    pub posrest_total: f64,
    /// Sum of distance-restraint energies.
    pub distanceres_total: f64,
    /// Solvent-accessible surface-area (SASA) solvation energy.
    pub sasa_total: f64,
    /// Energy contribution from holonomic constraints (SHAKE/LINCS).
    pub constraint_total: f64,
    /// gromosXX `entropy_term`: on the derivative object, dH/dλ · H of the potential
    /// (`energy_calculation.cc:87-89`); zero on the energies.
    pub entropy_term: f64,

    /// Kinetic energy from the current (new) velocities, used for thermostat scaling.
    pub kinetic_energy_new: f64,

    /// Kinetic energy per temperature bath: `com_kinetic_energy + ir_kinetic_energy`.
    pub kinetic_energy: Vec<f64>,
    /// Centre-of-mass translational kinetic energy of the temperature groups, per bath.
    pub com_kinetic_energy: Vec<f64>,
    /// Internal and rotational kinetic energy of the temperature groups, per bath.
    pub ir_kinetic_energy: Vec<f64>,

    /// Per-bath kinetic energy from the current (new) velocities — what a thermostat scales on.
    /// Empty for a single-bath run, where `kinetic_energy_new` carries the same number.
    pub kinetic_energy_new_bath: Vec<f64>,

    /// Bond energy per energy group (booked to the group of the bond's first atom).
    pub bond_energy: Vec<f64>,
    /// Angle energy per energy group.
    pub angle_energy: Vec<f64>,
    /// Improper-dihedral energy per energy group.
    pub improper_energy: Vec<f64>,
    /// Proper-dihedral energy per energy group.
    pub dihedral_energy: Vec<f64>,
    /// Cross-dihedral energy per energy group.
    pub crossdihedral_energy: Vec<f64>,
    /// Position-restraint energy per energy group.
    pub posrest_energy: Vec<f64>,
    /// Distance-restraint energy per energy group.
    pub distanceres_energy: Vec<f64>,
    /// Distance-field restraint energy per energy group.
    pub disfieldres_energy: Vec<f64>,
    /// Angle-restraint energy per energy group.
    pub angrest_energy: Vec<f64>,
    /// Dihedral-restraint energy per energy group.
    pub dihrest_energy: Vec<f64>,
    /// Constraint energy per energy group.
    pub constraints_energy: Vec<f64>,
    /// SASA energy per energy group.
    pub sasa_energy: Vec<f64>,
    /// SASA volume energy per energy group.
    pub sasa_volume_energy: Vec<f64>,
    /// RDC restraint energy per energy group.
    pub rdc_energy: Vec<f64>,

    /// Per-energy-group-pair LJ energies (energy_groups × energy_groups matrix).
    pub lj_energy: Vec<Vec<f64>>,
    /// Per-energy-group-pair CRF energies (energy_groups × energy_groups matrix).
    pub crf_energy: Vec<Vec<f64>>,
    /// Per-energy-group-pair lattice-sum k-space energies.
    pub ls_k_energy: Vec<Vec<f64>>,
    /// Per-energy-group-pair shifted-RF extra energy, original form.
    pub shift_extra_orig: Vec<Vec<f64>>,
    /// Per-energy-group-pair shifted-RF extra energy, physical form.
    pub shift_extra_phys: Vec<Vec<f64>>,

    /// Scalar virial (trace of the virial tensor), used for pressure coupling.
    pub virial_total: f64,
}

impl Energy {
    /// Create a zeroed energy accumulator for the given number of temperature and energy groups.
    pub fn new(num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        let mut e = Self {
            kinetic_total: 0.0,
            potential_total: 0.0,
            bond_total: 0.0,
            angle_total: 0.0,
            dihedral_total: 0.0,
            improper_total: 0.0,
            cross_dihedral_total: 0.0,
            lj_total: 0.0,
            crf_total: 0.0,
            ls_total: 0.0,
            ls_realspace_total: 0.0,
            ls_kspace_total: 0.0,
            special_total: 0.0,
            term_energies: Vec::new(),
            posrest_total: 0.0,
            distanceres_total: 0.0,
            sasa_total: 0.0,
            constraint_total: 0.0,
            entropy_term: 0.0,
            kinetic_energy_new: 0.0,
            kinetic_energy: Vec::new(),
            com_kinetic_energy: Vec::new(),
            ir_kinetic_energy: Vec::new(),
            kinetic_energy_new_bath: Vec::new(),
            bond_energy: Vec::new(),
            angle_energy: Vec::new(),
            improper_energy: Vec::new(),
            dihedral_energy: Vec::new(),
            crossdihedral_energy: Vec::new(),
            posrest_energy: Vec::new(),
            distanceres_energy: Vec::new(),
            disfieldres_energy: Vec::new(),
            angrest_energy: Vec::new(),
            dihrest_energy: Vec::new(),
            constraints_energy: Vec::new(),
            sasa_energy: Vec::new(),
            sasa_volume_energy: Vec::new(),
            rdc_energy: Vec::new(),
            lj_energy: Vec::new(),
            crf_energy: Vec::new(),
            ls_k_energy: Vec::new(),
            shift_extra_orig: Vec::new(),
            shift_extra_phys: Vec::new(),
            virial_total: 0.0,
        };
        e.resize(num_energy_groups, num_temperature_groups);
        e
    }

    /// Size the per-group arrays to `n_groups`, the group-pair matrices to `n_groups × n_groups`
    /// and the per-bath arrays to `n_baths` (gromosXX `Energy::resize(num, numb)`), zeroing them.
    pub fn resize(&mut self, n_groups: usize, n_baths: usize) {
        for v in self.group_arrays_mut() {
            v.clear();
            v.resize(n_groups, 0.0);
        }
        for m in self.group_matrices_mut() {
            m.clear();
            m.resize(n_groups, vec![0.0; n_groups]);
        }
        for v in [
            &mut self.kinetic_energy,
            &mut self.com_kinetic_energy,
            &mut self.ir_kinetic_energy,
        ] {
            v.clear();
            v.resize(n_baths, 0.0);
        }
    }

    /// Number of energy groups the per-group arrays are sized for.
    pub fn num_energy_groups(&self) -> usize {
        self.bond_energy.len()
    }

    fn group_arrays_mut(&mut self) -> [&mut Vec<f64>; 14] {
        [
            &mut self.bond_energy,
            &mut self.angle_energy,
            &mut self.improper_energy,
            &mut self.dihedral_energy,
            &mut self.crossdihedral_energy,
            &mut self.posrest_energy,
            &mut self.distanceres_energy,
            &mut self.disfieldres_energy,
            &mut self.angrest_energy,
            &mut self.dihrest_energy,
            &mut self.constraints_energy,
            &mut self.sasa_energy,
            &mut self.sasa_volume_energy,
            &mut self.rdc_energy,
        ]
    }

    fn group_matrices_mut(&mut self) -> [&mut Vec<Vec<f64>>; 5] {
        [
            &mut self.lj_energy,
            &mut self.crf_energy,
            &mut self.ls_k_energy,
            &mut self.shift_extra_orig,
            &mut self.shift_extra_phys,
        ]
    }

    /// With a single energy group every group cell is the total of its term, as gromosXX's
    /// group sums are: book the totals into group 0 (and pair (0, 0)). No-op for any other
    /// group count — those cells are booked by the terms themselves.
    pub fn book_totals_into_single_group(&mut self) {
        if self.bond_energy.len() != 1 || self.lj_energy.len() != 1 {
            return;
        }
        self.bond_energy[0] = self.bond_total;
        self.angle_energy[0] = self.angle_total;
        self.improper_energy[0] = self.improper_total;
        self.dihedral_energy[0] = self.dihedral_total;
        self.crossdihedral_energy[0] = self.cross_dihedral_total;
        self.posrest_energy[0] = self.posrest_total;
        self.distanceres_energy[0] = self.distanceres_total;
        self.constraints_energy[0] = self.constraint_total;
        self.sasa_energy[0] = self.sasa_total;
        self.lj_energy[0][0] = self.lj_total;
        self.crf_energy[0][0] = self.crf_total;
    }

    /// Total energy (kinetic + potential + special), mirroring GROMOS:
    /// total = potential_total + kinetic_total + special_total
    /// where special_total includes distanceres_total, posrest, angrest, etc.
    pub fn total(&self) -> f64 {
        self.kinetic_total + self.potential_total + self.special_total + self.distanceres_total
    }

    /// Update potential total from components (bonded + nonbonded, without special)
    pub fn update_potential_total(&mut self) {
        self.potential_total = self.bond_total
            + self.angle_total
            + self.dihedral_total
            + self.improper_total
            + self.cross_dihedral_total
            + self.lj_total
            + self.crf_total
            + self.ls_total
            + self.sasa_total;
    }

    /// Clear all energies (sizes are kept).
    pub fn clear(&mut self) {
        self.kinetic_total = 0.0;
        self.potential_total = 0.0;
        self.bond_total = 0.0;
        self.angle_total = 0.0;
        self.dihedral_total = 0.0;
        self.improper_total = 0.0;
        self.cross_dihedral_total = 0.0;
        self.lj_total = 0.0;
        self.crf_total = 0.0;
        self.ls_total = 0.0;
        self.ls_realspace_total = 0.0;
        self.ls_kspace_total = 0.0;
        self.special_total = 0.0;
        self.term_energies.clear();
        self.posrest_total = 0.0;
        self.distanceres_total = 0.0;
        self.sasa_total = 0.0;
        self.constraint_total = 0.0;
        self.entropy_term = 0.0;
        self.kinetic_energy_new = 0.0;
        self.kinetic_energy_new_bath.fill(0.0);
        self.virial_total = 0.0;

        self.kinetic_energy.fill(0.0);
        self.com_kinetic_energy.fill(0.0);
        self.ir_kinetic_energy.fill(0.0);
        for v in self.group_arrays_mut() {
            v.fill(0.0);
        }
        for m in self.group_matrices_mut() {
            for row in m.iter_mut() {
                row.fill(0.0);
            }
        }
    }
}

/// State of the system at a single time point
///
/// Contains positions, velocities, forces, and all associated properties
#[derive(Debug, Clone)]
pub struct State {
    /// Atomic positions (nm).
    pub pos: Vec<Vec3>,
    /// Atomic velocities (nm ps⁻¹).
    pub vel: Vec<Vec3>,
    /// Total forces on each atom (kJ mol⁻¹ nm⁻¹).
    pub force: Vec<Vec3>,
    /// Holonomic constraint forces (kJ mol⁻¹ nm⁻¹).
    pub constraint_force: Vec<Vec3>,
    /// Simulation box geometry.
    pub box_config: Box,
    /// Virial tensor (kJ mol⁻¹); used for pressure calculation.
    pub virial_tensor: Mat3,
    /// Kinetic-energy tensor (kJ mol⁻¹); used for pressure calculation.
    pub kinetic_energy_tensor: Mat3,
    /// Pressure tensor (bar); computed by `PressureCalculation`.
    pub pressure_tensor: Mat3,
    /// All energy components for this state.
    pub energies: Energy,
    /// dH/dλ of every energy component, laid out like `energies` (gromosXX
    /// `perturbed_energy_derivatives`); zero unless the run is perturbed.
    pub perturbed_energy_derivatives: Energy,
}

impl State {
    /// Create a zeroed state for `n_atoms` atoms with the given group counts.
    pub fn new(n_atoms: usize, num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            pos: vec![Vec3::ZERO; n_atoms],
            vel: vec![Vec3::ZERO; n_atoms],
            force: vec![Vec3::ZERO; n_atoms],
            constraint_force: vec![Vec3::ZERO; n_atoms],
            box_config: Box::vacuum(),
            virial_tensor: Mat3::ZERO,
            kinetic_energy_tensor: Mat3::ZERO,
            pressure_tensor: Mat3::ZERO,
            energies: Energy::new(num_temperature_groups, num_energy_groups),
            perturbed_energy_derivatives: Energy::new(num_temperature_groups, num_energy_groups),
        }
    }

    /// Clear forces and virial
    pub fn clear_forces(&mut self) {
        self.force.fill(Vec3::ZERO);
        self.virial_tensor = Mat3::ZERO;
        self.kinetic_energy_tensor = Mat3::ZERO;
    }

    /// Clear constraint forces
    pub fn clear_constraint_forces(&mut self) {
        self.constraint_force.fill(Vec3::ZERO);
    }

    /// Calculate kinetic energy
    ///
    /// KE = 0.5 * sum_i(m_i * v_i^2)
    pub fn calculate_kinetic_energy(&mut self, masses: &[f64]) {
        self.energies.kinetic_total = self
            .vel
            .iter()
            .zip(masses.iter())
            .map(|(v, &m)| 0.5 * m * v.length_squared())
            .sum();
    }

    /// Calculate kinetic energy per group
    pub fn calculate_kinetic_energy_groups(
        &mut self,
        masses: &[f64],
        temperature_groups: &[Vec<usize>],
    ) {
        for (group_idx, group) in temperature_groups.iter().enumerate() {
            self.energies.kinetic_energy[group_idx] = group
                .iter()
                .map(|&i| 0.5 * masses[i] * self.vel[i].length_squared())
                .sum();
        }

        // Total kinetic energy is sum of all groups
        self.energies.kinetic_total = self.energies.kinetic_energy.iter().sum();
    }

    /// Calculate temperature from kinetic energy
    ///
    /// T = 2 * KE / (k_B * N_dof)
    /// where N_dof = 3N - N_constraints
    pub fn temperature(&self, num_degrees_of_freedom: usize) -> f64 {
        2.0 * self.energies.kinetic_total / (crate::units::kB * num_degrees_of_freedom as f64)
    }

    /// Calculate pressure from virial
    ///
    /// P = (2*KE - Virial) / (3*V)
    pub fn pressure(&self) -> f64 {
        let volume = self.box_config.volume();
        let virial_trace =
            self.virial_tensor.x_axis.x + self.virial_tensor.y_axis.y + self.virial_tensor.z_axis.z;

        (2.0 * self.energies.kinetic_total - virial_trace) / (3.0 * volume)
    }
}

/// Configuration - manages current and previous states with double buffering
///
/// GROMOS uses double buffering (pointer swapping) for efficient state management.
/// This is translated to Rust using two State structs and an index.
#[derive(Debug, Clone)]
pub struct Configuration {
    state1: State,
    state2: State,
    current_idx: usize, // 0 or 1

    /// Lattice shifts for periodic boundary tracking (FEP)
    /// Stores how many times each atom crossed periodic boundaries
    /// Dimension: [n_atoms][3] for x, y, z shifts
    pub lattice_shifts: Vec<[i32; 3]>,
}

impl Configuration {
    /// Create a new double-buffered configuration for `n_atoms` atoms.
    pub fn new(n_atoms: usize, num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            state1: State::new(n_atoms, num_temperature_groups, num_energy_groups),
            state2: State::new(n_atoms, num_temperature_groups, num_energy_groups),
            current_idx: 0,
            lattice_shifts: vec![[0, 0, 0]; n_atoms],
        }
    }

    /// Get reference to current state
    #[inline]
    pub fn current(&self) -> &State {
        if self.current_idx == 0 {
            &self.state1
        } else {
            &self.state2
        }
    }

    /// Get mutable reference to current state
    #[inline]
    pub fn current_mut(&mut self) -> &mut State {
        if self.current_idx == 0 {
            &mut self.state1
        } else {
            &mut self.state2
        }
    }

    /// Get reference to old (previous) state
    #[inline]
    pub fn old(&self) -> &State {
        if self.current_idx == 0 {
            &self.state2
        } else {
            &self.state1
        }
    }

    /// Get mutable reference to old state
    #[inline]
    pub fn old_mut(&mut self) -> &mut State {
        if self.current_idx == 0 {
            &mut self.state2
        } else {
            &mut self.state1
        }
    }

    /// Exchange current and old states (zero-cost pointer swap)
    ///
    /// This is the Rust equivalent of GROMOS C++ pointer swapping
    #[inline]
    pub fn exchange_state(&mut self) {
        self.current_idx = 1 - self.current_idx;
    }

    /// Copy current state to old state
    pub fn copy_current_to_old(&mut self) {
        if self.current_idx == 0 {
            self.state2 = self.state1.clone();
        } else {
            self.state1 = self.state2.clone();
        }
    }
}

/// Stochastic dynamics additional variables
#[derive(Debug, Clone)]
pub struct StochasticVariables {
    /// Per-atom friction coefficients γ (ps⁻¹) for Langevin dynamics.
    pub gamma: Vec<f64>,
    /// Per-atom random forces (kJ mol⁻¹ nm⁻¹) drawn each step.
    pub random_force: Vec<Vec3>,
}

impl StochasticVariables {
    /// Create zeroed stochastic variables for `n_atoms` atoms.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            gamma: vec![0.0; n_atoms],
            random_force: vec![Vec3::ZERO; n_atoms],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_box_rectangular() {
        let box_config = Box::rectangular(3.0, 4.0, 5.0);

        assert_eq!(box_config.box_type, BoxType::Rectangular);
        assert_eq!(box_config.dimensions(), Vec3::new(3.0, 4.0, 5.0));
        assert!((box_config.volume() - 60.0).abs() < 1e-6);
    }

    #[test]
    fn test_energy_total() {
        let mut energy = Energy::new(1, 1);

        energy.kinetic_total = 100.0;
        energy.bond_total = 50.0;
        energy.lj_total = 30.0;
        energy.crf_total = 20.0;

        energy.update_potential_total();

        assert_eq!(energy.potential_total, 100.0);
        assert_eq!(energy.total(), 200.0);
    }

    #[test]
    fn energy_is_laid_out_per_group_pair_and_bath() {
        let mut e = Energy::new(2, 3);
        assert_eq!(e.num_energy_groups(), 3);
        assert_eq!((e.lj_energy.len(), e.ls_k_energy[2].len()), (3, 3));
        assert_eq!((e.kinetic_energy.len(), e.ir_kinetic_energy.len()), (2, 2));
        e.resize(1, 4);
        assert_eq!((e.rdc_energy.len(), e.shift_extra_phys.len()), (1, 1));
        assert_eq!(e.com_kinetic_energy.len(), 4);
        // a single group carries the totals
        e.bond_total = 1.5;
        e.lj_total = -2.0;
        e.posrest_total = 0.25;
        e.book_totals_into_single_group();
        assert_eq!(
            (e.bond_energy[0], e.lj_energy[0][0], e.posrest_energy[0]),
            (1.5, -2.0, 0.25)
        );
        e.clear();
        assert_eq!((e.bond_energy[0], e.lj_energy[0][0]), (0.0, 0.0));
        // more than one group: the cells are the terms' to book
        let mut two = Energy::new(1, 2);
        two.bond_total = 1.0;
        two.book_totals_into_single_group();
        assert_eq!(two.bond_energy, vec![0.0, 0.0]);
    }

    #[test]
    fn test_state_kinetic_energy() {
        let mut state = State::new(2, 1, 1);
        let masses = vec![1.0, 2.0];

        // Set velocities
        state.vel[0] = Vec3::new(1.0, 0.0, 0.0);
        state.vel[1] = Vec3::new(0.0, 2.0, 0.0);

        state.calculate_kinetic_energy(&masses);

        // KE = 0.5 * (1.0 * 1.0^2 + 2.0 * 2.0^2) = 0.5 * (1 + 8) = 4.5
        assert!((state.energies.kinetic_total - 4.5).abs() < 1e-6);
    }

    #[test]
    fn test_temperature_calculation() {
        let mut state = State::new(100, 1, 1);
        let masses = vec![1.0; 100];

        // Set random velocities for testing
        for vel in &mut state.vel {
            *vel = Vec3::new(0.1, 0.1, 0.1);
        }

        state.calculate_kinetic_energy(&masses);

        // 3N degrees of freedom (no constraints)
        let temp = state.temperature(300);

        assert!(temp > 0.0);
        assert!(temp < 100.0); // Sanity check
    }

    #[test]
    fn test_state_exchange() {
        let mut conf = Configuration::new(2, 1, 1);

        // Set some values in current state
        conf.current_mut().pos[0] = Vec3::new(1.0, 2.0, 3.0);

        // Exchange
        conf.exchange_state();

        // Old state should now have the values
        assert_eq!(conf.old().pos[0], Vec3::new(1.0, 2.0, 3.0));

        // Current state should be fresh
        assert_eq!(conf.current().pos[0], Vec3::ZERO);
    }

    #[test]
    fn test_clear_forces() {
        let mut state = State::new(10, 1, 1);

        // Set some forces
        state.force[0] = Vec3::new(1.0, 2.0, 3.0);
        state.virial_tensor = Mat3::from_cols(
            Vec3::new(10.0, 0.0, 0.0),
            Vec3::new(0.0, 10.0, 0.0),
            Vec3::new(0.0, 0.0, 10.0),
        );

        state.clear_forces();

        assert_eq!(state.force[0], Vec3::ZERO);
        assert_eq!(state.virial_tensor, Mat3::ZERO);
    }
}
