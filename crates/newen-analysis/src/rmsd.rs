//! Root Mean Square Deviation (RMSD) calculation

use newen_core::Vec3;

/// Calculate RMSD between two structures
pub fn calculate_rmsd(positions: &[Vec3], reference: &[Vec3]) -> f64 {
    if positions.len() != reference.len() {
        return f64::NAN;
    }

    let sum_sq: f64 = positions
        .iter()
        .zip(reference.iter())
        .map(|(p, r)| {
            let dx = p.x - r.x;
            let dy = p.y - r.y;
            let dz = p.z - r.z;
            dx * dx + dy * dy + dz * dz
        })
        .sum();

    (sum_sq / positions.len() as f64).sqrt()
}
