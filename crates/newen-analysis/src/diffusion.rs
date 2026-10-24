//! Diffusion coefficient calculation (MSD)

use newen_core::Vec3;

/// Calculate mean square displacement
pub fn calculate_msd(trajectory: &[Vec<Vec3>], tau: usize) -> f64 {
    if trajectory.len() <= tau {
        return 0.0;
    }

    let n_frames = trajectory.len() - tau;
    let n_atoms = trajectory[0].len();
    let mut msd = 0.0;

    for t in 0..n_frames {
        for i in 0..n_atoms {
            let r0 = &trajectory[t][i];
            let rt = &trajectory[t + tau][i];
            let dx = rt.x - r0.x;
            let dy = rt.y - r0.y;
            let dz = rt.z - r0.z;
            msd += dx * dx + dy * dy + dz * dz;
        }
    }

    (msd) / (n_frames as f64 * n_atoms as f64)
}

/// Calculate diffusion coefficient from MSD
pub fn calculate_diffusion_coefficient(msd: f64, tau: f64, dimensions: usize) -> f64 {
    msd / (2.0 * dimensions as f64 * tau)
}
