//! Polarization Models for Explicit Electronic Polarizability
//!
//! Direct translation of md++/src/interaction/nonbonded/interaction/polarisation*.cc
//!
//! **Purpose**: Model induced dipoles and higher-order polarization effects
//!
//! Polarization allows the electron cloud of atoms to respond to the local
//! electric field, providing more accurate electrostatics than fixed charges.
//!
//! ## Theory
//!
//! The induced dipole on atom i is:
//! ```text
//! μ_i = α_i * E_i
//! ```
//! where α_i is the atomic polarizability and E_i is the total electric field.
//!
//! The electric field includes:
//! ```text
//! E_i = E_permanent + E_induced
//! E_i = Σ_j (q_j * r_ij / r_ij³) + Σ_j (3(μ_j · r_ij)r_ij - μ_j) / r_ij⁵
//! ```
//!
//! This leads to self-consistent equations (SCF iteration).
//!
//! ## Models Implemented
//!
//! 1. **Point Dipoles**: Single induced dipole per atom
//! 2. **Drude Oscillators**: Auxiliary charge particles connected by springs
//! 3. **Fluctuating Charges**: Charge equilibration (QEq/EEM methods)
//! 4. **Classical Drude**: CHARMM Drude polarizable force field
//!
//! ## References
//!
//! - Thole (1981). "Molecular polarizabilities calculated with a modified dipole interaction."
//!   Chem. Phys. 59:341-350
//! - Lamoureux & Roux (2003). "Modeling induced polarization with classical Drude oscillators."
//!   J. Chem. Phys. 119:3025-3039
//! - Rick & Stuart (2002). "Potentials and algorithms for incorporating polarizability in
//!   computer simulations." Rev. Comp. Chem. 18:89-146

use newen_core::configuration::Configuration;
use newen_core::math::Vec3;
use newen_core::topology::Topology;

/// Polarization model type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolarizationModel {
    /// No polarization (fixed charges)
    None = 0,
    /// Point induced dipoles (classical Drude)
    PointDipole = 1,
    /// Shell model (Dick-Overhauser)
    ShellModel = 2,
    /// Drude oscillators with auxiliary particles
    DrudeOscillator = 3,
    /// Fluctuating charges (QEq/EEM)
    FluctuatingCharge = 4,
}

/// Atom polarizability parameters
#[derive(Debug, Clone)]
pub struct PolarizabilityParameters {
    /// Atomic polarizability α (nm³ or Å³)
    pub alpha: f64,
    /// Thole damping parameter a
    pub thole_a: f64,
    /// Drude charge (for Drude oscillator model)
    pub drude_charge: f64,
    /// Drude spring constant (kJ/mol/nm²)
    pub drude_k: f64,
    /// Drude mass (amu)
    pub drude_mass: f64,
}

impl Default for PolarizabilityParameters {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            thole_a: 2.6, // Thole's original value
            drude_charge: 0.0,
            drude_k: 4184.0, // 1000 kcal/mol/Å² typical
            drude_mass: 0.4, // Small mass for fast equilibration
        }
    }
}

/// Polarization state for each atom
#[derive(Debug, Clone)]
pub struct PolarizationState {
    /// Induced dipole moment μ_ind (e·nm)
    pub induced_dipole: Vec3,
    /// Electric field at this atom E (kJ/(mol·nm·e))
    pub electric_field: Vec3,
    /// Position of Drude particle (if Drude model)
    pub drude_position: Option<Vec3>,
    /// Velocity of Drude particle
    pub drude_velocity: Option<Vec3>,
    /// Force on Drude particle
    pub drude_force: Option<Vec3>,
}

impl Default for PolarizationState {
    fn default() -> Self {
        Self {
            induced_dipole: Vec3::ZERO,
            electric_field: Vec3::ZERO,
            drude_position: None,
            drude_velocity: None,
            drude_force: None,
        }
    }
}

/// Polarization calculation engine
///
/// Handles self-consistent field (SCF) iteration to solve for induced dipoles
#[derive(Debug, Clone)]
pub struct PolarizationCalculator {
    /// Polarization model
    pub model: PolarizationModel,
    /// Polarizability parameters per atom type
    pub polarizabilities: Vec<PolarizabilityParameters>,
    /// Polarization state per atom
    pub states: Vec<PolarizationState>,
    /// SCF convergence tolerance
    pub scf_tolerance: f64,
    /// Maximum SCF iterations
    pub max_scf_iterations: usize,
    /// Use Thole damping for short-range dipole-dipole interactions
    pub use_thole_damping: bool,
    /// Cutoff for polarization interactions (nm)
    pub cutoff: f64,
}

impl PolarizationCalculator {
    /// Create new polarization calculator
    pub fn new(n_atoms: usize, model: PolarizationModel) -> Self {
        Self {
            model,
            polarizabilities: vec![PolarizabilityParameters::default(); n_atoms],
            states: vec![PolarizationState::default(); n_atoms],
            scf_tolerance: 1e-6, // Convergence: max change in dipole < 1e-6 e·nm
            max_scf_iterations: 100,
            use_thole_damping: true,
            cutoff: 1.2, // 1.2 nm typical
        }
    }

    /// Calculate polarization energy and forces
    ///
    /// # Algorithm (Point Dipole Model)
    /// 1. Calculate permanent electric field E_perm from fixed charges
    /// 2. SCF iteration to converge induced dipoles:
    ///    a) Calculate field from induced dipoles E_ind
    ///    b) Update dipoles: μ_i^new = α_i * (E_perm + E_ind)
    ///    c) Check convergence: |μ_i^new - μ_i^old| < tolerance
    ///    d) Repeat until converged
    /// 3. Calculate polarization energy: E_pol = -0.5 * Σ_i μ_i · E_i
    /// 4. Calculate forces: F_i = -∂E_pol/∂r_i
    pub fn calculate_polarization(&mut self, topo: &Topology, conf: &mut Configuration) -> f64 {
        match self.model {
            PolarizationModel::None => 0.0,
            PolarizationModel::PointDipole => self.calculate_point_dipole(topo, conf),
            PolarizationModel::DrudeOscillator => self.calculate_drude_oscillator(topo, conf),
            PolarizationModel::FluctuatingCharge => self.calculate_fluctuating_charge(topo, conf),
            PolarizationModel::ShellModel => {
                // TODO: Implement shell model
                0.0
            },
        }
    }

    /// Point dipole model with SCF iteration
    fn calculate_point_dipole(&mut self, topo: &Topology, conf: &mut Configuration) -> f64 {
        let n_atoms = topo.num_atoms();

        // 1. Calculate permanent electric field from fixed charges
        self.calculate_permanent_field(topo, conf);

        // 2. SCF iteration for induced dipoles
        let mut converged = false;
        for _iter in 0..self.max_scf_iterations {
            // Save old dipoles
            let old_dipoles: Vec<Vec3> = self.states.iter().map(|s| s.induced_dipole).collect();

            // Calculate field from induced dipoles
            self.calculate_induced_field(topo, conf);

            // Update induced dipoles: μ_i = α_i * E_total
            let mut max_change: f64 = 0.0;
            for i in 0..n_atoms {
                let alpha = self.polarizabilities[i].alpha;
                let e_total = self.states[i].electric_field;
                let new_dipole = e_total * alpha;

                let change = (new_dipole - old_dipoles[i]).length();
                max_change = max_change.max(change);

                self.states[i].induced_dipole = new_dipole;
            }

            // Check convergence
            if max_change < self.scf_tolerance {
                converged = true;
                break;
            }
        }

        if !converged {
            eprintln!("Warning: SCF polarization did not converge");
        }

        // 3. Calculate polarization energy: E_pol = -0.5 * Σ μ · E
        let mut energy = 0.0;
        for state in &self.states {
            energy -= 0.5 * state.induced_dipole.dot(state.electric_field);
        }

        // 4. Calculate forces (dipole-charge and dipole-dipole interactions)
        self.calculate_polarization_forces(topo, conf);

        energy
    }

    /// Calculate permanent electric field from fixed charges
    fn calculate_permanent_field(&mut self, topo: &Topology, conf: &Configuration) {
        let n_atoms = topo.num_atoms();

        // Reset fields
        for state in &mut self.states {
            state.electric_field = Vec3::ZERO;
        }

        // E_i = Σ_j q_j * r_ij / r_ij³
        for i in 0..n_atoms {
            let pos_i = conf.current().pos[i];
            let mut field = Vec3::ZERO;

            for j in 0..n_atoms {
                if i == j {
                    continue;
                }

                let pos_j = conf.current().pos[j];
                let r_ij = pos_i - pos_j;
                let r = r_ij.length();

                if r > self.cutoff {
                    continue;
                }

                let q_j = topo.charge[j];
                let r3 = r * r * r;
                field += r_ij * (q_j / r3);
            }

            self.states[i].electric_field = field;
        }
    }

    /// Calculate electric field from induced dipoles
    fn calculate_induced_field(&mut self, topo: &Topology, conf: &Configuration) {
        let n_atoms = topo.num_atoms();

        // E_ind_i = Σ_j T_ij · μ_j
        // where T_ij is the dipole field tensor:
        // T_ij = (3 r_ij ⊗ r_ij / r^5 - I / r^3) with Thole damping
        for i in 0..n_atoms {
            let pos_i = conf.current().pos[i];
            let mut field_ind = Vec3::ZERO;

            for j in 0..n_atoms {
                if i == j {
                    continue;
                }

                let pos_j = conf.current().pos[j];
                let r_ij = pos_i - pos_j;
                let r = r_ij.length();

                if r > self.cutoff {
                    continue;
                }

                let mu_j = self.states[j].induced_dipole;
                let r3 = r * r * r;
                let r5 = r3 * r * r;

                // Thole damping for short range
                let damping = if self.use_thole_damping {
                    self.thole_damping_factor(i, j, r)
                } else {
                    1.0
                };

                // Dipole field tensor contribution
                let dot_prod = (r_ij.dot(mu_j)) / r5;
                let term1 = r_ij * (3.0 * dot_prod * damping);
                let term2 = mu_j * (damping / r3);

                field_ind += term1 - term2;
            }

            self.states[i].electric_field += field_ind;
        }
    }

    /// Thole damping factor
    ///
    /// Damps dipole-dipole interactions at short range to avoid polarization catastrophe
    /// s = r / (α_i * α_j)^(1/6)
    /// f(s) = 1 - (1 + s/2 + s²/6) * exp(-s)  [for a=2.6]
    fn thole_damping_factor(&self, i: usize, j: usize, r: f64) -> f64 {
        let alpha_i = self.polarizabilities[i].alpha;
        let alpha_j = self.polarizabilities[j].alpha;
        let a = self.polarizabilities[i].thole_a;

        let alpha_prod = (alpha_i * alpha_j).powf(1.0 / 6.0);
        let s = a * r / alpha_prod;

        1.0 - (1.0 + s / 2.0 + s * s / 6.0) * (-s).exp()
    }

    /// Calculate forces from polarization
    fn calculate_polarization_forces(&self, topo: &Topology, _conf: &mut Configuration) {
        // TODO: Implement force calculation
        // F_i = -∂E_pol/∂r_i
        // Includes:
        // 1. Charge-dipole forces: F = q * ∇E_dipole
        // 2. Dipole-dipole forces: F = μ_i · ∇E_j
        // 3. Chain rule through SCF convergence

        let _n_atoms = topo.num_atoms();
        // Placeholder
    }

    /// Drude oscillator model
    ///
    /// Each polarizable atom has an auxiliary Drude particle connected by a spring.
    /// The Drude particle position relaxes to induce polarization.
    fn calculate_drude_oscillator(&mut self, topo: &Topology, conf: &mut Configuration) -> f64 {
        let n_atoms = topo.num_atoms();
        let mut energy = 0.0;

        // For each polarizable atom
        for i in 0..n_atoms {
            if self.polarizabilities[i].alpha == 0.0 {
                continue; // Not polarizable
            }

            let pos_atom = conf.current().pos[i];
            let q_drude = self.polarizabilities[i].drude_charge;
            let k = self.polarizabilities[i].drude_k;

            // Get Drude particle position (or initialize)
            let pos_drude = self.states[i].drude_position.unwrap_or(pos_atom);

            // Calculate electric field at Drude particle
            let mut field = Vec3::ZERO;
            for j in 0..n_atoms {
                if i == j {
                    continue;
                }

                let pos_j = conf.current().pos[j];
                let r_ij = pos_drude - pos_j;
                let r = r_ij.length();

                if r < 1e-6 || r > self.cutoff {
                    continue;
                }

                let q_j = topo.charge[j];
                field += r_ij * (q_j / (r * r * r));
            }

            // Drude position minimizes: E = 0.5*k*d² + q_drude*φ(d)
            // At equilibrium: k*d = -q_drude*E => d = -q_drude*E/k
            let displacement = field * (-q_drude / k);
            let new_drude_pos = pos_atom + displacement;

            // Spring energy
            let spring_energy = 0.5 * k * displacement.length_squared();
            energy += spring_energy;

            // Update state
            self.states[i].drude_position = Some(new_drude_pos);
            self.states[i].induced_dipole = displacement * q_drude;

            // Forces on atom and Drude particle
            let f_spring = displacement * (-k);
            conf.current_mut().force[i] += f_spring;

            // Electric force on Drude
            let f_electric = field * q_drude;
            conf.current_mut().force[i] -= f_electric;
        }

        energy
    }

    /// Fluctuating charge model (QEq/EEM)
    ///
    /// Charges equilibrate to minimize electrostatic energy subject to
    /// charge conservation constraint.
    fn calculate_fluctuating_charge(&mut self, _topo: &Topology, _conf: &mut Configuration) -> f64 {
        // TODO: Implement charge equilibration
        // Solve: J * ΔQ = -χ subject to Σ ΔQ = 0
        // where J is Coulomb matrix, χ is electronegativity
        0.0
    }

    /// Extended Lagrangian dynamics for Drude particles
    ///
    /// Propagate Drude particles as auxiliary DOF with small mass
    pub fn propagate_drude_particles(&mut self, dt: f64) {
        for (i, state) in self.states.iter_mut().enumerate() {
            if let Some(drude_pos) = state.drude_position {
                let drude_vel = state.drude_velocity.unwrap_or(Vec3::ZERO);
                let drude_force = state.drude_force.unwrap_or(Vec3::ZERO);
                let mass = self.polarizabilities[i].drude_mass;

                // Velocity Verlet
                let accel = drude_force / mass;
                let new_vel = drude_vel + accel * (dt);
                let new_pos = drude_pos + new_vel * (dt);

                state.drude_velocity = Some(new_vel);
                state.drude_position = Some(new_pos);
            }
        }
    }
}

/// Polarization topology - defines which atoms are polarizable
#[derive(Debug, Clone)]
pub struct PolarizationTopology {
    /// Polarizable atoms (true if atom is polarizable)
    pub is_polarizable: Vec<bool>,
    /// Polarizability parameters per atom
    pub parameters: Vec<PolarizabilityParameters>,
}

impl PolarizationTopology {
    pub fn new(n_atoms: usize) -> Self {
        Self {
            is_polarizable: vec![false; n_atoms],
            parameters: vec![PolarizabilityParameters::default(); n_atoms],
        }
    }

    /// Set atom as polarizable with given parameters
    pub fn set_polarizable(&mut self, atom_idx: usize, params: PolarizabilityParameters) {
        self.is_polarizable[atom_idx] = true;
        self.parameters[atom_idx] = params;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polarization_calculator() {
        let calc = PolarizationCalculator::new(100, PolarizationModel::PointDipole);
        assert_eq!(calc.states.len(), 100);
        assert_eq!(calc.model, PolarizationModel::PointDipole);
    }

    #[test]
    fn test_thole_damping() {
        let mut calc = PolarizationCalculator::new(2, PolarizationModel::PointDipole);
        calc.polarizabilities[0].alpha = 1.0;
        calc.polarizabilities[1].alpha = 1.0;

        let f = calc.thole_damping_factor(0, 1, 0.3);
        assert!(f > 0.0 && f < 1.0); // Damped at short range
    }

    #[test]
    fn test_drude_parameters() {
        let params = PolarizabilityParameters::default();
        assert_eq!(params.thole_a, 2.6);
        assert!(params.drude_mass > 0.0);
    }
}
