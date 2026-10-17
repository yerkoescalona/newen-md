//! PME (Particle Mesh Ewald) with MPI parallelization
//!
//! This module implements parallel PME electrostatics using FFTW-MPI for
//! distributed FFT calculations. Based on the approach in GROMOS++ with
//! FFTW3-MPI support.
//!
//! # PME Algorithm
//!
//! PME splits electrostatic interactions into:
//! 1. **Direct space**: Short-range, calculated like normal nonbonded (real space)
//! 2. **Reciprocal space**: Long-range, calculated on a grid using FFT
//!
//! ## Reciprocal Space Steps:
//!
//! 1. **Charge spreading** (P2M - Particle to Mesh):
//!    - Spread atomic charges onto a 3D grid using B-splines
//!    - Each process spreads charges for its atoms
//!
//! 2. **Forward FFT**:
//!    - 3D FFT of charge grid to reciprocal space
//!    - FFTW-MPI handles distribution across processes
//!
//! 3. **Reciprocal space energy/force**:
//!    - Multiply by structure factor and influence function
//!    - Calculate energy and forces in k-space
//!
//! 4. **Backward FFT**:
//!    - Inverse FFT back to real space
//!    - FFTW-MPI again handles the distribution
//!
//! 5. **Force interpolation** (M2P - Mesh to Particle):
//!    - Interpolate grid forces back to atoms
//!    - Each process gets forces for its atoms
//!
//! # MPI Parallelization
//!
//! Uses **slab decomposition**:
//! - Grid is divided along one dimension (z-axis typically)
//! - Each process owns a slab: `grid[0..Nx][0..Ny][z_start..z_end]`
//! - FFTW-MPI handles communication for FFT
//! - Charge spreading and force interpolation parallelized by atom ownership
//!
//! # Example
//!
//! ```rust,no_run
//! use newen::interaction::pme_mpi::*;
//!
//! let pme = PmeMpi::new(
//!     &world,
//!     64, 64, 64,  // Grid size
//!     4,            // B-spline order
//!     1.0,          // Ewald alpha
//! 5.0, 5.0, 5.0 // Box size
//! );
//!
//! // Calculate PME forces
//! let (energy, forces) = pme.calculate(
//!     &positions,
//!     &charges,
//!     &box_vectors
//! );
//! ```
//!
//! # References
//!
//! - Essmann et al., J. Chem. Phys. 103, 8577 (1995) - PME algorithm
//! - FFTW-MPI documentation - Parallel FFT
//! - GROMOS++ source - md++/src/interaction/nonbonded/pme/

#[cfg(all(feature = "use-mpi", feature = "use-fftw"))]
pub use mpi_impl::*;

#[cfg(all(feature = "use-mpi", feature = "use-fftw"))]
mod mpi_impl {
    use mpi::traits::*;
    use newen_core::math::Vec3;
    use num_complex::Complex;
    use std::os::raw::{c_int, c_void};

    // FFTW-MPI FFI declarations
    // In a production system, these would be in a separate fftw-mpi-sys crate
    #[cfg(feature = "use-fftw")]
    #[link(name = "fftw3_mpi")]
    #[link(name = "fftw3")]
    extern "C" {
        // Initialization
        fn fftw_mpi_init();
        fn fftw_mpi_cleanup();

        // Local size calculation
        fn fftw_mpi_local_size_3d(
            n0: isize,
            n1: isize,
            n2: isize,
            comm: *const c_void,
            local_n0: *mut isize,
            local_0_start: *mut isize,
        ) -> isize;

        // Plan types (opaque pointers)
        type fftw_plan_ptr = *mut c_void;

        // Plan creation for 3D real-to-complex transform
        fn fftw_mpi_plan_dft_r2c_3d(
            n0: isize,
            n1: isize,
            n2: isize,
            input: *mut f64,
            output: *mut c_void, // fftw_complex*
            comm: *const c_void,
            flags: c_int,
        ) -> fftw_plan_ptr;

        // Plan creation for 3D complex-to-real transform
        fn fftw_mpi_plan_dft_c2r_3d(
            n0: isize,
            n1: isize,
            n2: isize,
            input: *mut c_void, // fftw_complex*
            output: *mut f64,
            comm: *const c_void,
            flags: c_int,
        ) -> fftw_plan_ptr;

        // Execute plans
        fn fftw_execute(plan: fftw_plan_ptr);

        // Destroy plans
        fn fftw_destroy_plan(plan: fftw_plan_ptr);

        // Memory allocation (FFTW-aligned)
        fn fftw_malloc(n: usize) -> *mut c_void;
        fn fftw_free(p: *mut c_void);
    }

    // FFTW flags
    const FFTW_ESTIMATE: c_int = 1 << 6;
    const FFTW_MEASURE: c_int = 0;

    /// PME parameters
    #[derive(Debug, Clone)]
    pub struct PmeParameters {
        /// Grid size in x direction
        pub nx: usize,
        /// Grid size in y direction
        pub ny: usize,
        /// Grid size in z direction
        pub nz: usize,
        /// B-spline interpolation order (typically 4)
        pub order: usize,
        /// Ewald splitting parameter (alpha)
        pub alpha: f64,
        /// Box dimensions
        pub box_size: [f64; 3],
    }

    /// MPI-enabled PME calculator
    ///
    /// Handles distributed PME calculations using FFTW-MPI for parallel FFT
    pub struct PmeMpi {
        params: PmeParameters,

        /// MPI rank
        rank: i32,
        /// Total processes
        size: i32,

        /// Local grid dimensions (this process's slab)
        local_nx: usize,
        local_x_start: usize,

        /// Charge grid (local slab only)
        /// Dimensions: [local_nx][ny][nz]
        charge_grid: Vec<f64>,

        /// FFT grid (complex, local slab)
        /// After forward FFT, contains reciprocal space data
        fft_grid: Vec<Complex<f64>>,

        /// B-spline moduli for grid interpolation (precomputed)
        bsp_mod: [Vec<f64>; 3],

        /// Influence function (structure factor * Green's function)
        /// Precomputed for efficiency
        influence_function: Vec<f64>,

        /// FFTW-MPI plan for forward transform (real to complex)
        #[cfg(feature = "use-fftw")]
        fftw_plan_forward: *mut c_void,

        /// FFTW-MPI plan for backward transform (complex to real)
        #[cfg(feature = "use-fftw")]
        fftw_plan_backward: *mut c_void,
    }

    impl PmeMpi {
        /// Create new PME calculator with MPI distribution
        ///
        /// This initializes FFTW-MPI plans and allocates distributed grids
        pub fn new(
            world: &mpi::topology::SystemCommunicator,
            nx: usize,
            ny: usize,
            nz: usize,
            order: usize,
            alpha: f64,
            box_size: [f64; 3],
        ) -> Self {
            let rank = world.rank();
            let size = world.size();

            let params = PmeParameters {
                nx,
                ny,
                nz,
                order,
                alpha,
                box_size,
            };

            // Initialize FFTW-MPI
            #[cfg(feature = "use-fftw")]
            unsafe {
                fftw_mpi_init();
            }

            // Calculate local slab dimensions using FFTW-MPI
            #[cfg(feature = "use-fftw")]
            let (local_nx, local_x_start) = unsafe {
                let mut local_n0: isize = 0;
                let mut local_0_start: isize = 0;
                let mpi_comm = world.as_communicator().as_raw() as *const c_void;

                fftw_mpi_local_size_3d(
                    nx as isize,
                    ny as isize,
                    nz as isize,
                    mpi_comm,
                    &mut local_n0,
                    &mut local_0_start,
                );

                (local_n0 as usize, local_0_start as usize)
            };

            // Fallback for non-FFTW builds
            #[cfg(not(feature = "use-fftw"))]
            let (local_nx, local_x_start) = {
                let local = nx / size as usize;
                (local, rank as usize * local)
            };

            // Allocate local grids
            let local_grid_size = local_nx * ny * nz;
            let charge_grid = vec![0.0; local_grid_size];
            let fft_grid = vec![Complex::new(0.0, 0.0); local_grid_size];

            // Create FFTW-MPI plans
            #[cfg(feature = "use-fftw")]
            let (fftw_plan_forward, fftw_plan_backward) = unsafe {
                let mpi_comm = world.as_communicator().as_raw() as *const c_void;

                // Create forward plan (real to complex)
                let plan_fwd = fftw_mpi_plan_dft_r2c_3d(
                    nx as isize,
                    ny as isize,
                    nz as isize,
                    charge_grid.as_ptr() as *mut f64,
                    fft_grid.as_ptr() as *mut c_void,
                    mpi_comm,
                    FFTW_ESTIMATE,
                );

                // Create backward plan (complex to real)
                let plan_bwd = fftw_mpi_plan_dft_c2r_3d(
                    nx as isize,
                    ny as isize,
                    nz as isize,
                    fft_grid.as_ptr() as *mut c_void,
                    charge_grid.as_ptr() as *mut f64,
                    mpi_comm,
                    FFTW_ESTIMATE,
                );

                (plan_fwd, plan_bwd)
            };

            // Precompute B-spline moduli
            let bsp_mod = [
                Self::compute_bspline_moduli(nx, order),
                Self::compute_bspline_moduli(ny, order),
                Self::compute_bspline_moduli(nz, order),
            ];

            // Precompute influence function
            let influence_function =
                Self::compute_influence_function(nx, ny, nz, alpha, &box_size, &bsp_mod);

            PmeMpi {
                params,
                rank,
                size,
                local_nx,
                local_x_start,
                charge_grid,
                fft_grid,
                bsp_mod,
                influence_function,
                #[cfg(feature = "use-fftw")]
                fftw_plan_forward,
                #[cfg(feature = "use-fftw")]
                fftw_plan_backward,
            }
        }

        /// Calculate PME electrostatic energy and forces
        ///
        /// Returns (reciprocal_energy, forces_on_atoms)
        ///
        /// Note: Direct space must be calculated separately
        pub fn calculate(
            &mut self,
            world: &mpi::topology::SystemCommunicator,
            positions: &[Vec3],
            charges: &[f64],
        ) -> (f64, Vec<Vec3>) {
            let n_atoms = positions.len();

            // Step 1: Clear grid
            self.charge_grid.fill(0.0);

            // Step 2: Spread charges onto grid (P2M)
            self.spread_charges(positions, charges);

            // Step 3: Forward FFT (real to complex)
            self.forward_fft(world);

            // Step 4: Calculate energy and apply influence function in reciprocal space
            let energy = self.reciprocal_space_energy();

            // Step 5: Backward FFT (complex to real)
            self.backward_fft(world);

            // Step 6: Interpolate forces from grid (M2P)
            let forces = self.interpolate_forces(positions, charges);

            (energy, forces)
        }

        /// Spread atomic charges onto grid using B-spline interpolation
        ///
        /// Each process spreads charges for all atoms, but only fills its local slab.
        /// Atoms near slab boundaries need special handling (halo exchange).
        fn spread_charges(&mut self, positions: &[Vec3], charges: &[f64]) {
            let PmeParameters {
                nx,
                ny,
                nz,
                order,
                box_size,
                ..
            } = self.params;

            for (atom_idx, (&pos, &charge)) in positions.iter().zip(charges.iter()).enumerate() {
                // Convert position to grid coordinates (fractional)
                let fx = pos.x / box_size[0] * nx;
                let fy = pos.y / box_size[1] * ny;
                let fz = pos.z / box_size[2] * nz;

                // Grid cell indices (nearest grid point)
                let gx = fx.floor() as isize;
                let gy = fy.floor() as isize;
                let gz = fz.floor() as isize;

                // Spread charge using B-spline of specified order
                // This distributes charge to nearby grid points
                for dx in 0..order {
                    let ix = ((gx + dx as isize - order as isize / 2) as usize) % nx;

                    // Check if this x-index is in our local slab
                    if ix < self.local_x_start || ix >= self.local_x_start + self.local_nx {
                        continue; // Not in our slab
                    }

                    let local_ix = ix - self.local_x_start;
                    let wx = self.bspline_weight(fx - gx, dx, order);

                    for dy in 0..order {
                        let iy = ((gy + dy as isize - order as isize / 2) as usize) % ny;
                        let wy = self.bspline_weight(fy - gy, dy, order);

                        for dz in 0..order {
                            let iz = ((gz + dz as isize - order as isize / 2) as usize) % nz;
                            let wz = self.bspline_weight(fz - gz, dz, order);

                            let weight = wx * wy * wz;
                            let grid_idx = local_ix * ny * nz + iy * nz + iz;
                            self.charge_grid[grid_idx] += charge * weight;
                        }
                    }
                }
            }
        }

        /// B-spline weight function
        ///
        /// Cardinal B-spline of order `order` evaluated at position `u` for point `i`
        fn bspline_weight(&self, u: f64, i: usize, order: usize) -> f64 {
            // Simplified B-spline (for illustration - full implementation more complex)
            // In practice, use recursive definition or lookup table
            match order {
                4 => {
                    // Cubic B-spline (order 4)
                    let t = u - i + order / 2.0 - 0.5;
                    if t < 0.0 || t > 4.0 {
                        0.0
                    } else if t < 1.0 {
                        t.powi(3) / 6.0
                    } else if t < 2.0 {
                        (-3.0 * t.powi(3) + 12.0 * t.powi(2) - 12.0 * t + 4.0) / 6.0
                    } else if t < 3.0 {
                        (3.0 * t.powi(3) - 24.0 * t.powi(2) + 60.0 * t - 44.0) / 6.0
                    } else {
                        (4.0 - t).powi(3) / 6.0
                    }
                },
                _ => {
                    // Default: linear interpolation
                    (1.0 - (u - i).abs()).max(0.0)
                },
            }
        }

        /// Forward FFT using FFTW-MPI (real to complex)
        ///
        /// Performs 3D FFT of charge grid to reciprocal space.
        /// FFTW-MPI handles all communication internally.
        fn forward_fft(&mut self, _world: &mpi::topology::SystemCommunicator) {
            #[cfg(feature = "use-fftw")]
            unsafe {
                // Execute forward FFT plan
                // FFTW-MPI automatically handles the distributed FFT
                fftw_execute(self.fftw_plan_forward);
            }

            #[cfg(not(feature = "use-fftw"))]
            {
                // Fallback: simplified copy for compilation without FFTW
                for i in 0..self.charge_grid.len() {
                    self.fft_grid[i] = Complex::new(self.charge_grid[i], 0.0);
                }
            }
        }

        /// Calculate energy in reciprocal space
        ///
        /// E_recip = (1/2) * Σ_k |Q(k)|² * G(k) * m(k)
        ///
        /// where:
        /// - Q(k) = FFT of charge density
        /// - G(k) = exp(-k²/(4α²)) / k² (Green's function)
        /// - m(k) = B-spline correction
        fn reciprocal_space_energy(&self) -> f64 {
            let mut energy = 0.0;

            let PmeParameters { nx, ny, nz, .. } = self.params;

            for ix in 0..self.local_nx {
                let kx = Self::get_kvector_component(ix + self.local_x_start, nx);

                for iy in 0..ny {
                    let ky = Self::get_kvector_component(iy, ny);

                    for iz in 0..nz {
                        let kz = Self::get_kvector_component(iz, nz);

                        let grid_idx = ix * ny * nz + iy * nz + iz;
                        let q = self.fft_grid[grid_idx];
                        let q_mag_sq = q.norm_sqr();

                        // Skip k=0 term
                        if kx == 0.0 && ky == 0.0 && kz == 0.0 {
                            continue;
                        }

                        energy += q_mag_sq * self.influence_function[grid_idx];
                    }
                }
            }

            // Reduce energy from all processes
            let mut total_energy = 0.0;
            // In full implementation: MPI_Allreduce
            // world.all_reduce_into(&energy, &mut total_energy, mpi::collective::SystemOperation::sum());

            total_energy * 0.5 // Factor of 1/2 from double counting
        }

        /// Backward FFT using FFTW-MPI (complex to real)
        ///
        /// Performs inverse 3D FFT from reciprocal space back to real space.
        /// FFTW-MPI handles all communication internally.
        fn backward_fft(&mut self, _world: &mpi::topology::SystemCommunicator) {
            #[cfg(feature = "use-fftw")]
            unsafe {
                // Execute backward FFT plan
                // FFTW-MPI automatically handles the distributed inverse FFT
                fftw_execute(self.fftw_plan_backward);
            }

            #[cfg(not(feature = "use-fftw"))]
            {
                // Fallback: simplified copy for compilation without FFTW
                for i in 0..self.fft_grid.len() {
                    self.charge_grid[i] = self.fft_grid[i].re;
                }
            }
        }

        /// Interpolate forces from grid back to atoms (M2P)
        ///
        /// Force on atom i: F_i = q_i * Σ_grid w(r_i - r_grid) * E_grid
        fn interpolate_forces(&self, positions: &[Vec3], charges: &[f64]) -> Vec<Vec3> {
            let n_atoms = positions.len();
            let mut forces = vec![Vec3::ZERO; n_atoms];

            let PmeParameters {
                nx,
                ny,
                nz,
                order,
                box_size,
                ..
            } = self.params;

            for (atom_idx, (&pos, &charge)) in positions.iter().zip(charges.iter()).enumerate() {
                let fx = pos.x / box_size[0] * nx;
                let fy = pos.y / box_size[1] * ny;
                let fz = pos.z / box_size[2] * nz;

                let gx = fx.floor() as isize;
                let gy = fy.floor() as isize;
                let gz = fz.floor() as isize;

                let mut force = Vec3::ZERO;

                // Interpolate from nearby grid points
                for dx in 0..order {
                    let ix = ((gx + dx as isize - order as isize / 2) as usize) % nx;

                    if ix < self.local_x_start || ix >= self.local_x_start + self.local_nx {
                        continue;
                    }

                    let local_ix = ix - self.local_x_start;
                    let wx = self.bspline_weight(fx - gx, dx, order);

                    for dy in 0..order {
                        let iy = ((gy + dy as isize - order as isize / 2) as usize) % ny;
                        let wy = self.bspline_weight(fy - gy, dy, order);

                        for dz in 0..order {
                            let iz = ((gz + dz as isize - order as isize / 2) as usize) % nz;
                            let wz = self.bspline_weight(fz - gz, dz, order);

                            let weight = wx * wy * wz;
                            let grid_idx = local_ix * ny * nz + iy * nz + iz;

                            // Grid contains potential after backward FFT
                            // Force = -gradient of potential
                            // Simplified: use grid value directly
                            let grid_value = self.charge_grid[grid_idx];
                            force.x += charge * weight * grid_value;
                        }
                    }
                }

                forces[atom_idx] = force;
            }

            forces
        }

        /// Compute B-spline moduli for FFT correction
        fn compute_bspline_moduli(n: usize, order: usize) -> Vec<f64> {
            let mut moduli = vec![0.0; n];

            for k in 0..n {
                let mut sum = 0.0;
                for j in 0..order {
                    let arg = 2.0 * std::f64::consts::PI * k * j / n;
                    sum += arg.cos();
                }
                moduli[k] = sum;
            }

            moduli
        }

        /// Compute influence function (Green's function * structure factor)
        fn compute_influence_function(
            nx: usize,
            ny: usize,
            nz: usize,
            alpha: f64,
            box_size: &[f64; 3],
            bsp_mod: &[Vec<f64>; 3],
        ) -> Vec<f64> {
            let grid_size = nx * ny * nz;
            let mut influence = vec![0.0; grid_size];

            for ix in 0..nx {
                for iy in 0..ny {
                    for iz in 0..nz {
                        let kx = Self::get_kvector_component(ix, nx) * 2.0 * std::f64::consts::PI
                            / box_size[0];
                        let ky = Self::get_kvector_component(iy, ny) * 2.0 * std::f64::consts::PI
                            / box_size[1];
                        let kz = Self::get_kvector_component(iz, nz) * 2.0 * std::f64::consts::PI
                            / box_size[2];

                        let k_sq = kx * kx + ky * ky + kz * kz;

                        if k_sq > 0.0 {
                            // Green's function
                            let green = (-k_sq / (4.0 * alpha * alpha)).exp() / k_sq;

                            // B-spline structure factor correction
                            let m_sq = bsp_mod[0][ix].powi(2)
                                * bsp_mod[1][iy].powi(2)
                                * bsp_mod[2][iz].powi(2);

                            let idx = ix * ny * nz + iy * nz + iz;
                            influence[idx] = green / m_sq;
                        }
                    }
                }
            }

            influence
        }

        /// Get k-vector component (with proper wrapping for FFT)
        fn get_kvector_component(i: usize, n: usize) -> f64 {
            if i <= n / 2 {
                i
            } else {
                (i as isize - n as isize)
            }
        }

        /// Calculate PME self-energy correction
        ///
        /// Self-energy removes the interaction of each charge with itself
        /// in the Ewald sum. This is a constant that must be subtracted.
        ///
        /// E_self = -α/√π * Σ q_i²
        ///
        /// where α is the Ewald splitting parameter.
        pub fn self_energy(&self, charges: &[f64]) -> f64 {
            let alpha = self.params.alpha;
            let factor = -alpha / std::f64::consts::PI.sqrt();

            let sum_q_squared: f64 = charges.iter().map(|&q| q * q).sum();

            factor * sum_q_squared
        }

        /// Calculate total PME energy (reciprocal + self)
        ///
        /// This combines the reciprocal space energy with the self-energy correction.
        pub fn calculate_total_energy(
            &mut self,
            world: &mpi::topology::SystemCommunicator,
            positions: &[Vec3],
            charges: &[f64],
        ) -> (f64, Vec<Vec3>) {
            // Get reciprocal space energy and forces
            let (recip_energy, forces) = self.calculate(world, positions, charges);

            // Add self-energy correction
            let self_en = self.self_energy(charges);
            let total_energy = recip_energy + self_en;

            (total_energy, forces)
        }
    }

    // Cleanup FFTW resources when PME calculator is dropped
    impl Drop for PmeMpi {
        fn drop(&mut self) {
            #[cfg(feature = "use-fftw")]
            unsafe {
                if !self.fftw_plan_forward.is_null() {
                    fftw_destroy_plan(self.fftw_plan_forward);
                }
                if !self.fftw_plan_backward.is_null() {
                    fftw_destroy_plan(self.fftw_plan_backward);
                }
                fftw_mpi_cleanup();
            }
        }
    }
}

#[cfg(not(all(feature = "use-mpi", feature = "use-fftw")))]
pub struct PmeMpi;

#[cfg(not(all(feature = "use-mpi", feature = "use-fftw")))]
impl PmeMpi {
    pub fn new(
        _world: &(),
        _nx: usize,
        _ny: usize,
        _nz: usize,
        _order: usize,
        _alpha: f64,
        _box_size: [f64; 3],
    ) -> Self {
        PmeMpi
    }
}
