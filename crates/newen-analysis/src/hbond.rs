//! Hydrogen bond analysis

/// Hydrogen bond criteria
pub struct HBondCriteria {
    /// Maximum donor-acceptor distance (nm)
    pub max_distance: f64,
    /// Maximum donor-H-acceptor angle (degrees)
    pub max_angle: f64,
}

impl Default for HBondCriteria {
    fn default() -> Self {
        Self {
            max_distance: 0.35,
            max_angle: 30.0,
        }
    }
}

/// Detected hydrogen bond
pub struct HBond {
    pub donor_index: usize,
    pub hydrogen_index: usize,
    pub acceptor_index: usize,
    pub distance: f64,
    pub angle: f64,
}
