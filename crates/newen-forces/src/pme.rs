//! Particle Mesh Ewald implementation
//!
//! Long-range electrostatics using PME/SPME methods

use newen_core::Vec3;
use num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;

/// PME calculator with FFT
pub struct Pme {
    /// Grid dimensions
    grid_size: [usize; 3],
    /// Ewald splitting parameter
    alpha: f64,
    /// B-spline order
    // TODO: use in B-spline charge spreading once PME innerloop is implemented
    #[allow(dead_code)]
    order: usize,
    /// Charge grid
    charge_grid: Vec<Complex64>,
    /// FFT planners
    fft_forward: Arc<dyn Fft<f64>>,
    // TODO: use in back-transform once PME force calculation is implemented
    #[allow(dead_code)]
    fft_inverse: Arc<dyn Fft<f64>>,
}

impl Pme {
    pub fn new(grid_size: [usize; 3], alpha: f64, order: usize) -> Self {
        let total_size = grid_size[0] * grid_size[1] * grid_size[2];

        let mut planner = FftPlanner::new();
        let fft_forward = planner.plan_fft_forward(total_size);
        let fft_inverse = planner.plan_fft_inverse(total_size);

        Self {
            grid_size,
            alpha,
            order,
            charge_grid: vec![Complex64::new(0.0, 0.0); total_size],
            fft_forward,
            fft_inverse,
        }
    }

    /// Spread charges to grid using B-spline interpolation
    pub fn spread_charges(&mut self, positions: &[Vec3], charges: &[f64], box_size: &[f64; 3]) {
        // Clear grid
        for c in &mut self.charge_grid {
            *c = Complex64::new(0.0, 0.0);
        }

        for (pos, &charge) in positions.iter().zip(charges.iter()) {
            // Convert position to grid coordinates
            let gx = (pos.x / box_size[0]) * self.grid_size[0] as f64;
            let gy = (pos.y / box_size[1]) * self.grid_size[1] as f64;
            let gz = (pos.z / box_size[2]) * self.grid_size[2] as f64;

            // B-spline spreading (simplified - 1st order)
            let ix = gx as usize % self.grid_size[0];
            let iy = gy as usize % self.grid_size[1];
            let iz = gz as usize % self.grid_size[2];

            let idx = ix + iy * self.grid_size[0] + iz * self.grid_size[0] * self.grid_size[1];
            self.charge_grid[idx] += Complex64::new(charge, 0.0);
        }
    }

    /// Calculate reciprocal space energy
    pub fn reciprocal_energy(&mut self, box_size: &[f64; 3]) -> f64 {
        // Forward FFT
        self.fft_forward.process(&mut self.charge_grid);

        let volume = box_size[0] * box_size[1] * box_size[2];
        let prefactor = 1.0 / (2.0 * std::f64::consts::PI * volume);

        let mut energy = 0.0;

        for kz in 0..self.grid_size[2] {
            for ky in 0..self.grid_size[1] {
                for kx in 0..self.grid_size[0] {
                    if kx == 0 && ky == 0 && kz == 0 {
                        continue; // Skip k=0
                    }

                    // Wave vectors
                    let mx = if kx <= self.grid_size[0] / 2 {
                        kx
                    } else {
                        kx - self.grid_size[0]
                    };
                    let my = if ky <= self.grid_size[1] / 2 {
                        ky
                    } else {
                        ky - self.grid_size[1]
                    };
                    let mz = if kz <= self.grid_size[2] / 2 {
                        kz
                    } else {
                        kz - self.grid_size[2]
                    };

                    let kx_real = 2.0 * std::f64::consts::PI * mx as f64 / box_size[0];
                    let ky_real = 2.0 * std::f64::consts::PI * my as f64 / box_size[1];
                    let kz_real = 2.0 * std::f64::consts::PI * mz as f64 / box_size[2];

                    let k2 = kx_real * kx_real + ky_real * ky_real + kz_real * kz_real;

                    // Gaussian screening
                    let exp_factor = (-k2 / (4.0 * self.alpha * self.alpha)).exp();

                    let idx =
                        kx + ky * self.grid_size[0] + kz * self.grid_size[0] * self.grid_size[1];
                    let rho_k = &self.charge_grid[idx];
                    let rho_k_sq = rho_k.norm_sqr();

                    energy += prefactor * exp_factor / k2 * rho_k_sq;
                }
            }
        }

        energy
    }
}

/// B-spline basis function
pub fn bspline(order: usize, x: f64) -> f64 {
    match order {
        1 => {
            if (0.0..1.0).contains(&x) {
                1.0
            } else {
                0.0
            }
        },
        2 => {
            if (0.0..1.0).contains(&x) {
                x
            } else if (1.0..2.0).contains(&x) {
                2.0 - x
            } else {
                0.0
            }
        },
        n => {
            let n_f = n as f64;
            (x / (n_f - 1.0)) * bspline(n - 1, x)
                + ((n_f - x) / (n_f - 1.0)) * bspline(n - 1, x - 1.0)
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bspline_normalization() {
        // B-splines should integrate to 1
        let order = 4;
        let mut sum = 0.0;
        let dx = 0.01;
        let mut x = 0.0;
        while x < order as f64 {
            sum += bspline(order, x) * dx;
            x += dx;
        }
        assert!((sum - 1.0).abs() < 0.05);
    }

    #[test]
    fn test_pme_creation() {
        let pme = Pme::new([32, 32, 32], 3.5, 4);
        assert_eq!(pme.grid_size, [32, 32, 32]);
    }
}
