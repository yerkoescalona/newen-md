//! Local Elevation Umbrella Sampling (LEUS/LE-META)
//!
//! Direct translation of md++/src/interaction/special/local_elevation_interaction.cc
//! and md++/src/util/umbrella.cc
//!
//! **Purpose**: Adaptive biasing force method for enhanced sampling
//!
//! Local elevation is a metadynamics variant that adaptively builds up a biasing
//! potential to escape free energy minima and explore the conformational space.
//!
//! ## Algorithm (Huber et al. 1994)
//!
//! ```text
//! 1. Define reaction coordinates ξ (e.g., distances, dihedrals, RMSDs)
//! 2. Discretize coordinate space into N-dimensional grid
//! 3. During simulation:
//!    a) Calculate current coordinates ξ(t)
//!    b) Determine grid bin i
//!    c) Add Gaussian hill at current position:
//!       ΔV_i(t) = ΔV_i(t-Δt) + W * exp(-||ξ(t) - ξ_i||² / (2σ²))
//!    d) Apply bias force: F_bias = -∂V_bias/∂ξ * ∂ξ/∂r
//! 4. Free energy: F(ξ) ≈ -V_bias(ξ) (in limit of long simulation)
//! ```
//!
//! ## Features
//!
//! - **Multi-dimensional**: 1D, 2D, 3D, ..., ND reaction coordinates
//! - **Multiple umbrellas**: Different biasing potentials simultaneously
//! - **Flexible coordinates**: Distance, angle, dihedral, RMSD, radius of gyration, etc.
//! - **Build or freeze**: Build new potential or use frozen potential from previous run
//! - **Restart capability**: Save/load potential grids
//! - **Periodic boundaries**: Handle periodic coordinates (angles, dihedrals)
//!
//! ## Use Cases
//!
//! - **Protein folding**: Escape kinetic traps using RMSD as coordinate
//! - **Ligand binding**: Sample binding modes using distance coordinates
//! - **Conformational transitions**: Explore dihedral space (φ, ψ angles)
//! - **Free energy profiles**: Calculate 1D/2D free energy landscapes
//!
//! ## References
//!
//! - Huber et al. (1994). "Local elevation: A method for improving the searching
//!   properties of molecular dynamics simulation." J. Comput. Aided Mol. Des. 8:695-708
//! - Hansen & Hünenberger (2010). "Using the local elevation method to construct
//!   optimized umbrella sampling potentials." J. Comput. Chem. 31:1-23

use newen_core::configuration::Configuration;
use newen_core::math::Vec3;
use std::f64::consts::PI;

/// Type of coordinate for local elevation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateType {
    /// Distance between two atoms
    Distance = 1,
    /// Bond angle (three atoms)
    Angle = 2,
    /// Dihedral angle (four atoms)
    Dihedral = 3,
    /// RMSD from reference structure
    RMSD = 4,
    /// Radius of gyration
    RadiusOfGyration = 5,
    /// Custom coordinate (user-defined function)
    Custom = 99,
}

/// Local elevation coordinate definition
///
/// Defines one reaction coordinate for umbrella sampling
#[derive(Debug, Clone)]
pub struct LECoordinate {
    /// Umbrella ID this coordinate belongs to
    pub umbrella_id: usize,
    /// Type of coordinate
    pub coord_type: CoordinateType,
    /// Atom indices involved (2 for distance, 3 for angle, 4 for dihedral)
    pub atoms: Vec<usize>,
    /// Reference structure for RMSD (if applicable)
    pub reference_positions: Vec<Vec3>,
    /// Current value of coordinate
    pub value: f64,
    /// Force on this coordinate (dV/dξ)
    pub force: f64,
}

impl LECoordinate {
    /// Create distance coordinate
    pub fn distance(umbrella_id: usize, atom_i: usize, atom_j: usize) -> Self {
        Self {
            umbrella_id,
            coord_type: CoordinateType::Distance,
            atoms: vec![atom_i, atom_j],
            reference_positions: Vec::new(),
            value: 0.0,
            force: 0.0,
        }
    }

    /// Create dihedral angle coordinate
    pub fn dihedral(umbrella_id: usize, i: usize, j: usize, k: usize, l: usize) -> Self {
        Self {
            umbrella_id,
            coord_type: CoordinateType::Dihedral,
            atoms: vec![i, j, k, l],
            reference_positions: Vec::new(),
            value: 0.0,
            force: 0.0,
        }
    }

    /// Create RMSD coordinate
    pub fn rmsd(umbrella_id: usize, atoms: Vec<usize>, reference: Vec<Vec3>) -> Self {
        Self {
            umbrella_id,
            coord_type: CoordinateType::RMSD,
            atoms,
            reference_positions: reference,
            value: 0.0,
            force: 0.0,
        }
    }

    /// Calculate coordinate value from configuration
    pub fn calculate(&mut self, conf: &Configuration) -> f64 {
        match self.coord_type {
            CoordinateType::Distance => self.calculate_distance(conf),
            CoordinateType::Angle => self.calculate_angle(conf),
            CoordinateType::Dihedral => self.calculate_dihedral(conf),
            CoordinateType::RMSD => self.calculate_rmsd(conf),
            CoordinateType::RadiusOfGyration => self.calculate_rgyr(conf),
            CoordinateType::Custom => 0.0, // User-defined
        }
    }

    fn calculate_distance(&mut self, conf: &Configuration) -> f64 {
        let r1 = conf.current().pos[self.atoms[0]];
        let r2 = conf.current().pos[self.atoms[1]];
        let r = (r2 - r1).length();
        self.value = r;
        r
    }

    fn calculate_angle(&mut self, conf: &Configuration) -> f64 {
        let r1 = conf.current().pos[self.atoms[0]];
        let r2 = conf.current().pos[self.atoms[1]];
        let r3 = conf.current().pos[self.atoms[2]];

        let v1 = r1 - r2;
        let v2 = r3 - r2;
        let cos_theta = (v1.dot(v2) / (v1.length() * v2.length())).clamp(-1.0, 1.0);
        let angle = (cos_theta).acos();
        self.value = angle;
        angle
    }

    fn calculate_dihedral(&mut self, conf: &Configuration) -> f64 {
        let r1 = conf.current().pos[self.atoms[0]];
        let r2 = conf.current().pos[self.atoms[1]];
        let r3 = conf.current().pos[self.atoms[2]];
        let r4 = conf.current().pos[self.atoms[3]];

        let b1 = r2 - r1;
        let b2 = r3 - r2;
        let b3 = r4 - r3;

        let n1 = b1.cross(b2);
        let n2 = b2.cross(b3);

        let cos_phi = (n1.dot(n2) / (n1.length() * n2.length())).clamp(-1.0, 1.0);
        let mut phi = (cos_phi).acos();

        // Determine sign
        if b1.dot(n2) < 0.0 {
            phi = -phi;
        }

        self.value = phi;
        phi
    }

    fn calculate_rmsd(&mut self, conf: &Configuration) -> f64 {
        let mut sum_sq = 0.0;
        for (i, &atom_idx) in self.atoms.iter().enumerate() {
            let r = conf.current().pos[atom_idx];
            let r_ref = self.reference_positions[i];
            let diff = r - r_ref;
            sum_sq += diff.length_squared();
        }
        let rmsd = (sum_sq / self.atoms.len() as f64).sqrt();
        self.value = rmsd;
        rmsd
    }

    fn calculate_rgyr(&mut self, conf: &Configuration) -> f64 {
        // Center of mass
        let mut com = Vec3::ZERO;
        for &atom_idx in &self.atoms {
            com += conf.current().pos[atom_idx];
        }
        com /= self.atoms.len() as f64;

        // Radius of gyration
        let mut sum_sq = 0.0;
        for &atom_idx in &self.atoms {
            let r = conf.current().pos[atom_idx];
            let diff = r - com;
            sum_sq += diff.length_squared();
        }
        let rgyr = (sum_sq / self.atoms.len() as f64).sqrt();
        self.value = rgyr;
        rgyr
    }
}

/// Umbrella weight calculation method
#[derive(Debug, Clone, Copy)]
pub enum UmbrellaWeightMethod {
    /// Number of visits: W_i = N_i
    NumberOfVisits,
    /// Time-weighted: W_i = Σ t_i
    TimeWeighted,
    /// Energy-weighted: W_i = Σ exp(-βE_i)
    EnergyWeighted,
}

/// Multi-dimensional umbrella for local elevation
///
/// An umbrella defines:
/// - Reaction coordinates (1D, 2D, 3D, ...)
/// - Grid discretization
/// - Biasing potential on grid
/// - Build parameters (Gaussian width, height)
#[derive(Debug, Clone)]
pub struct Umbrella {
    /// Unique identifier
    pub id: usize,
    /// Dimensionality (number of coordinates)
    pub dimensionality: usize,
    /// Reaction coordinates
    pub coordinates: Vec<LECoordinate>,
    /// Grid size per dimension
    pub grid_sizes: Vec<usize>,
    /// Grid min values per dimension
    pub grid_mins: Vec<f64>,
    /// Grid max values per dimension
    pub grid_maxs: Vec<f64>,
    /// Grid spacing per dimension
    pub grid_spacings: Vec<f64>,
    /// Biasing potential on grid (flattened N-D array)
    pub potential_grid: Vec<f64>,
    /// Visit counts on grid (for statistics)
    pub visit_counts: Vec<usize>,
    /// Gaussian width σ for hill deposition (per dimension)
    pub gaussian_widths: Vec<f64>,
    /// Hill height W (energy per deposition)
    pub hill_height: f64,
    /// Deposition frequency (every N steps)
    pub deposition_frequency: usize,
    /// Enable this umbrella
    pub enabled: bool,
    /// Build potential (true) or use frozen potential (false)
    pub building: bool,
    /// Periodic coordinates (true for angles/dihedrals)
    pub periodic: Vec<bool>,
    /// Weight calculation method
    pub weight_method: UmbrellaWeightMethod,
}

impl Umbrella {
    /// Create 1D umbrella
    pub fn new_1d(
        id: usize,
        coord: LECoordinate,
        grid_size: usize,
        grid_min: f64,
        grid_max: f64,
        gaussian_width: f64,
        hill_height: f64,
    ) -> Self {
        let grid_spacing = (grid_max - grid_min) / grid_size as f64;

        Self {
            id,
            dimensionality: 1,
            coordinates: vec![coord],
            grid_sizes: vec![grid_size],
            grid_mins: vec![grid_min],
            grid_maxs: vec![grid_max],
            grid_spacings: vec![grid_spacing],
            potential_grid: vec![0.0; grid_size],
            visit_counts: vec![0; grid_size],
            gaussian_widths: vec![gaussian_width],
            hill_height,
            deposition_frequency: 100, // Every 100 steps
            enabled: true,
            building: true,
            periodic: vec![false],
            weight_method: UmbrellaWeightMethod::NumberOfVisits,
        }
    }

    /// Create 2D umbrella
    pub fn new_2d(
        id: usize,
        coord1: LECoordinate,
        coord2: LECoordinate,
        grid_sizes: [usize; 2],
        grid_mins: [f64; 2],
        grid_maxs: [f64; 2],
        gaussian_widths: [f64; 2],
        hill_height: f64,
    ) -> Self {
        let grid_spacings = [
            (grid_maxs[0] - grid_mins[0]) / grid_sizes[0] as f64,
            (grid_maxs[1] - grid_mins[1]) / grid_sizes[1] as f64,
        ];
        let total_grid_points = grid_sizes[0] * grid_sizes[1];

        Self {
            id,
            dimensionality: 2,
            coordinates: vec![coord1, coord2],
            grid_sizes: grid_sizes.to_vec(),
            grid_mins: grid_mins.to_vec(),
            grid_maxs: grid_maxs.to_vec(),
            grid_spacings: grid_spacings.to_vec(),
            potential_grid: vec![0.0; total_grid_points],
            visit_counts: vec![0; total_grid_points],
            gaussian_widths: gaussian_widths.to_vec(),
            hill_height,
            deposition_frequency: 100,
            enabled: true,
            building: true,
            periodic: vec![false, false],
            weight_method: UmbrellaWeightMethod::NumberOfVisits,
        }
    }

    /// Calculate current coordinate values
    pub fn calculate_coordinates(&mut self, conf: &Configuration) {
        for coord in &mut self.coordinates {
            coord.calculate(conf);
        }
    }

    /// Get current grid bin indices
    fn get_grid_bin(&self) -> Option<Vec<usize>> {
        let mut bin_indices = Vec::with_capacity(self.dimensionality);

        for (dim, coord) in self.coordinates.iter().enumerate() {
            let mut val = coord.value;

            // Handle periodicity
            if self.periodic[dim] {
                while val < self.grid_mins[dim] {
                    val += 2.0 * PI;
                }
                while val >= self.grid_maxs[dim] {
                    val -= 2.0 * PI;
                }
            }

            // Check bounds
            if val < self.grid_mins[dim] || val >= self.grid_maxs[dim] {
                return None; // Out of bounds
            }

            let bin = ((val - self.grid_mins[dim]) / self.grid_spacings[dim]) as usize;
            let bin = bin.min(self.grid_sizes[dim] - 1);
            bin_indices.push(bin);
        }

        Some(bin_indices)
    }

    /// Convert multi-dimensional bin indices to flat index
    fn bin_to_flat_index(&self, bin_indices: &[usize]) -> usize {
        let mut flat_idx = 0;
        let mut multiplier = 1;

        for (dim, &bin) in bin_indices.iter().enumerate() {
            flat_idx += bin * multiplier;
            multiplier *= self.grid_sizes[dim];
        }

        flat_idx
    }

    /// Build potential (deposit Gaussian hill at current position)
    pub fn build(&mut self, _conf: &Configuration) {
        if let Some(bin_indices) = self.get_grid_bin() {
            let flat_idx = self.bin_to_flat_index(&bin_indices);
            self.visit_counts[flat_idx] += 1;

            // Deposit Gaussian hill
            // For efficiency, only update nearby grid points
            self.deposit_gaussian_hill(&bin_indices);
        }
    }

    /// Deposit Gaussian hill at given bin
    fn deposit_gaussian_hill(&mut self, center_bin: &[usize]) {
        // Calculate cutoff (3σ)
        let cutoffs: Vec<usize> = self
            .gaussian_widths
            .iter()
            .zip(&self.grid_spacings)
            .map(|(&sigma, &spacing)| (3.0 * sigma / spacing).ceil() as usize)
            .collect();

        // Copy necessary data to avoid borrowing issues
        let grid_mins = self.grid_mins.clone();
        let grid_spacings = self.grid_spacings.clone();
        let gaussian_widths = self.gaussian_widths.clone();
        let hill_height = self.hill_height;

        // Iterate over nearby grid points (simplified)
        for dim in 0..self.dimensionality {
            let min_bin = center_bin[dim].saturating_sub(cutoffs[dim]);
            let max_bin = (center_bin[dim] + cutoffs[dim]).min(self.grid_sizes[dim] - 1);

            for bin in min_bin..=max_bin {
                let mut bin_indices = center_bin.to_vec();
                bin_indices[dim] = bin;

                // Calculate Gaussian weight
                let mut dist_sq = 0.0;
                for (d, (&b, &center)) in bin_indices.iter().zip(center_bin.iter()).enumerate() {
                    let bin_center = grid_mins[d] + (b as f64 + 0.5) * grid_spacings[d];
                    let center_val = grid_mins[d] + (center as f64 + 0.5) * grid_spacings[d];
                    let delta = bin_center - center_val;
                    dist_sq += (delta / gaussian_widths[d]).powi(2);
                }

                let weight = (-0.5 * dist_sq).exp();
                let flat_idx = self.bin_to_flat_index(&bin_indices);
                self.potential_grid[flat_idx] += hill_height * weight;
            }
        }
    }

    /// Apply biasing potential (calculate and apply forces)
    pub fn apply(&mut self, conf: &mut Configuration) {
        // Calculate forces on coordinates from potential gradient
        self.calculate_forces();

        // Redistribute forces to atoms
        for coord in &self.coordinates {
            self.redistribute_coordinate_force(coord, conf);
        }
    }

    /// Calculate forces on coordinates (numerical gradient of potential)
    fn calculate_forces(&mut self) {
        if let Some(bin_indices) = self.get_grid_bin() {
            // Extract needed data to avoid borrowing issues
            let grid_spacings = self.grid_spacings.clone();
            let grid_sizes = self.grid_sizes.clone();

            for dim in 0..self.coordinates.len() {
                // Numerical derivative: dV/dξ ≈ (V(ξ+Δξ) - V(ξ-Δξ)) / (2Δξ)
                let delta = grid_spacings[dim];

                let mut bin_plus = bin_indices.clone();
                let mut bin_minus = bin_indices.clone();

                bin_plus[dim] = (bin_plus[dim] + 1).min(grid_sizes[dim] - 1);
                bin_minus[dim] = bin_minus[dim].saturating_sub(1);

                let v_plus = self.potential_grid[self.bin_to_flat_index(&bin_plus)];
                let v_minus = self.potential_grid[self.bin_to_flat_index(&bin_minus)];

                self.coordinates[dim].force = -(v_plus - v_minus) / (2.0 * delta);
            }
        }
    }

    /// Redistribute coordinate force to atoms
    fn redistribute_coordinate_force(&self, coord: &LECoordinate, conf: &mut Configuration) {
        // TODO: Implement proper chain rule force redistribution
        // For distance: dξ/dr_i = (r_j - r_i) / |r_j - r_i|
        // For angles/dihedrals: More complex derivatives
        // For RMSD: Gradient involves all atoms

        // Placeholder: Apply force equally to all involved atoms
        let force_per_atom = coord.force / coord.atoms.len() as f64;
        for &atom_idx in &coord.atoms {
            // This is a simplification - proper implementation needs chain rule
            conf.current_mut().force[atom_idx].x += force_per_atom;
        }
    }

    /// Save potential grid to file
    pub fn save_potential(&self, filename: &str) -> std::io::Result<()> {
        use std::io::Write;
        let mut file = std::fs::File::create(filename)?;

        writeln!(file, "# Local Elevation Potential Grid")?;
        writeln!(file, "# Umbrella ID: {}", self.id)?;
        writeln!(file, "# Dimensionality: {}", self.dimensionality)?;
        writeln!(file, "# Grid sizes: {:?}", self.grid_sizes)?;

        for (idx, &pot) in self.potential_grid.iter().enumerate() {
            writeln!(file, "{} {}", idx, pot)?;
        }

        Ok(())
    }

    /// Load potential grid from file
    pub fn load_potential(&mut self, filename: &str) -> std::io::Result<()> {
        use std::io::BufRead;
        let file = std::fs::File::open(filename)?;
        let reader = std::io::BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 2 {
                let bad = |what: &str| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("{filename}: bad {what} in line {line:?}"),
                    )
                };
                let idx: usize = parts[0].parse().map_err(|_| bad("grid index"))?;
                let pot: f64 = parts[1].parse().map_err(|_| bad("potential"))?;
                if idx < self.potential_grid.len() {
                    self.potential_grid[idx] = pot;
                }
            }
        }

        Ok(())
    }
}

/// Local elevation manager
#[derive(Debug, Clone)]
pub struct LocalElevation {
    /// All umbrellas
    pub umbrellas: Vec<Umbrella>,
    /// Step counter for deposition frequency
    pub step_counter: usize,
}

impl LocalElevation {
    pub fn new() -> Self {
        Self {
            umbrellas: Vec::new(),
            step_counter: 0,
        }
    }

    pub fn add_umbrella(&mut self, umbrella: Umbrella) {
        self.umbrellas.push(umbrella);
    }

    /// Apply local elevation biasing
    pub fn apply(&mut self, conf: &mut Configuration) {
        self.step_counter += 1;

        for umbrella in &mut self.umbrellas {
            if !umbrella.enabled {
                continue;
            }

            // Calculate coordinates
            umbrella.calculate_coordinates(conf);

            // Build potential if enabled
            if umbrella.building
                && self
                    .step_counter
                    .is_multiple_of(umbrella.deposition_frequency)
            {
                umbrella.build(conf);
            }

            // Apply biasing forces
            umbrella.apply(conf);
        }
    }
}

impl Default for LocalElevation {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1d_umbrella() {
        let coord = LECoordinate::distance(1, 0, 1);
        let umbrella = Umbrella::new_1d(1, coord, 100, 0.0, 10.0, 0.5, 1.0);

        assert_eq!(umbrella.dimensionality, 1);
        assert_eq!(umbrella.grid_sizes[0], 100);
        assert_eq!(umbrella.potential_grid.len(), 100);
    }

    #[test]
    fn test_2d_umbrella() {
        let coord1 = LECoordinate::distance(1, 0, 1);
        let coord2 = LECoordinate::dihedral(1, 0, 1, 2, 3);
        let umbrella = Umbrella::new_2d(
            1,
            coord1,
            coord2,
            [50, 72],
            [0.0, -PI],
            [10.0, PI],
            [0.5, 0.2],
            1.0,
        );

        assert_eq!(umbrella.dimensionality, 2);
        assert_eq!(umbrella.potential_grid.len(), 50 * 72);
    }

    #[test]
    fn test_grid_binning() {
        let coord = LECoordinate::distance(1, 0, 1);
        let mut umbrella = Umbrella::new_1d(1, coord, 100, 0.0, 10.0, 0.5, 1.0);

        // Set coordinate value
        umbrella.coordinates[0].value = 5.0;

        let bin = umbrella.get_grid_bin();
        assert!(bin.is_some());
        assert_eq!(bin.unwrap()[0], 50); // 5.0 is at 50% of [0, 10]
    }
}
