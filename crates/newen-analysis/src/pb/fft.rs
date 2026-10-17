//! The FFT grid of gromos++'s periodic Poisson solver: `pb::FFTGridType` (`src/pb/FFTGridType.cc`)
//! and the two transforms of `FFTPoissonIterator` (`fft_grid`, `bft_grid`), with rustfft in place
//! of FFTW.
//!
//! Grid conventions, all from gromos++: point `(i, j, k)` sits at `(i·drx, j·dry, k·drz)`, the flat
//! index is `k + nz·(j + ny·i)` (z fastest), and the k-vector component of index `i` is
//! `((i + n/2) % n − n/2)·dk` with integer `n/2`. A complex grid is one `Complex64` per point (the
//! C++ interleaves real and imaginary parts in a `vector<double>`; [`slot`] reads that layout).
//! The forward transform has FFTW's sign (−i) and is scaled by `vol/ngr3`, the backward one (+i)
//! by `1/vol`, so a round trip is the identity.

use rustfft::num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;

/// gromos++ `FFTGridType`: `n` points per axis spanning `len` (the driver passes `len = n·h`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FftGrid {
    pub n: [usize; 3],
    pub ngr3: usize,
    pub len: [f64; 3],
    /// `0.5·len`.
    pub center: [f64; 3],
    pub vol: f64,
    /// Real-space spacing `len/n`.
    pub dr: [f64; 3],
    /// k-space spacing `2π/len`.
    pub dk: [f64; 3],
    /// Cubelets per half axis of a grid cell (`@cubesFFT`): a cell is split into `(2·ncubes+1)³`.
    pub ncubes: usize,
}

impl FftGrid {
    /// `FFTGridType(ngrdx, …, xlen, …, ncubes)` (`FFTGridType.cc:34-67`); `pi` is gromos++'s
    /// `4·atan(1)`.
    pub fn new(n: [usize; 3], len: [f64; 3], ncubes: usize, pi: f64) -> Self {
        let pi2 = 2.0 * pi;
        FftGrid {
            n,
            ngr3: n[0] * n[1] * n[2],
            len,
            center: [0.5 * len[0], 0.5 * len[1], 0.5 * len[2]],
            vol: len[0] * len[1] * len[2],
            dr: [
                len[0] / n[0] as f64,
                len[1] / n[1] as f64,
                len[2] / n[2] as f64,
            ],
            dk: [pi2 / len[0], pi2 / len[1], pi2 / len[2]],
            ncubes,
        }
    }

    /// Flat index of point `(i, j, k)`.
    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k + self.n[2] * (j + self.n[1] * i)
    }

    /// The k-vector component of every index along `axis`.
    pub fn k_vectors(&self, axis: usize) -> Vec<f64> {
        (0..self.n[axis])
            .map(|i| k_component(i, self.n[axis], self.dk[axis]))
            .collect()
    }
}

/// `((i + n/2) % n − n/2)·dk` with integer `n/2`: for even `n`, index `n/2` maps to `−n/2`.
#[inline]
pub fn k_component(i: usize, n: usize, dk: f64) -> f64 {
    let (i, n) = (i as i64, n as i64);
    (((i + n / 2) % n) - n / 2) as f64 * dk
}

/// Element `m` of a complex grid seen as gromos++'s interleaved `vector<double>`: the real part of
/// point `m/2` for even `m`, its imaginary part for odd `m`.
#[inline]
pub fn slot(v: &[Complex64], m: usize) -> f64 {
    let c = v[m / 2];
    if m.is_multiple_of(2) {
        c.re
    } else {
        c.im
    }
}

/// The 3D complex transforms of one grid, as sequential 1D passes over z, y and x.
pub struct Fft3d {
    grid: FftGrid,
    forward: [Arc<dyn Fft<f64>>; 3],
    inverse: [Arc<dyn Fft<f64>>; 3],
    scratch: Vec<Complex64>,
    lines: Vec<Complex64>,
}

impl Fft3d {
    pub fn new(grid: FftGrid) -> Self {
        let mut planner = FftPlanner::<f64>::new();
        let forward = grid.n.map(|n| planner.plan_fft_forward(n));
        let inverse = grid.n.map(|n| planner.plan_fft_inverse(n));
        let scratch_len = forward
            .iter()
            .chain(inverse.iter())
            .map(|p| p.get_inplace_scratch_len())
            .max()
            .unwrap_or(0);
        Fft3d {
            grid,
            forward,
            inverse,
            scratch: vec![Complex64::default(); scratch_len],
            lines: vec![Complex64::default(); grid.ngr3],
        }
    }

    /// The unscaled 3D transform of `data` in place.
    fn transform(&mut self, data: &mut [Complex64], forward: bool) {
        let [nx, ny, nz] = self.grid.n;
        let plans = if forward {
            &self.forward
        } else {
            &self.inverse
        };
        // z: contiguous lines.
        plans[2].process_with_scratch(data, &mut self.scratch);
        // y: per x plane, gather the lines k-major, transform, scatter back.
        let plane = ny * nz;
        for i in 0..nx {
            let p = &mut data[i * plane..(i + 1) * plane];
            for j in 0..ny {
                for k in 0..nz {
                    self.lines[k * ny + j] = p[j * nz + k];
                }
            }
            plans[1].process_with_scratch(&mut self.lines[..plane], &mut self.scratch);
            for j in 0..ny {
                for k in 0..nz {
                    p[j * nz + k] = self.lines[k * ny + j];
                }
            }
        }
        // x: gather every (j, k) line.
        for i in 0..nx {
            for jk in 0..plane {
                self.lines[jk * nx + i] = data[i * plane + jk];
            }
        }
        plans[0].process_with_scratch(&mut self.lines, &mut self.scratch);
        for i in 0..nx {
            for jk in 0..plane {
                data[i * plane + jk] = self.lines[jk * nx + i];
            }
        }
    }

    /// `fft_grid(r, k)`: `k = FFT₋(r)·vol/ngr3`, `r` unchanged (`FFTPoissonIterator.cc:954-1007`).
    pub fn forward(&mut self, r: &[Complex64], k: &mut [Complex64]) {
        k.copy_from_slice(r);
        self.forward_in_place(k);
    }

    /// `fft_grid(x, x)`: the k-space result overwrites the input.
    pub fn forward_in_place(&mut self, x: &mut [Complex64]) {
        self.transform(x, true);
        let fact = self.grid.vol / self.grid.ngr3 as f64;
        for c in x.iter_mut() {
            *c = Complex64::new(c.re * fact, c.im * fact);
        }
    }

    /// `bft_grid(x)`: `x = FFT₊(x)/vol` (`FFTPoissonIterator.cc:1010-1053`).
    pub fn backward(&mut self, x: &mut [Complex64]) {
        self.transform(x, false);
        let fact = 1.0 / self.grid.vol;
        for c in x.iter_mut() {
            *c = Complex64::new(c.re * fact, c.im * fact);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn k_index_wraps_at_half_the_grid() {
        // n = 6: 0 1 2 −3 −2 −1; n = 5: 0 1 2 −2 −1.
        let k: Vec<f64> = (0..6).map(|i| k_component(i, 6, 1.0)).collect();
        assert_eq!(k, [0.0, 1.0, 2.0, -3.0, -2.0, -1.0]);
        let k: Vec<f64> = (0..5).map(|i| k_component(i, 5, 0.5)).collect();
        assert_eq!(k, [0.0, 0.5, 1.0, -1.0, -0.5]);
    }

    #[test]
    fn grid_follows_fftgridtype() {
        let g = FftGrid::new([40, 40, 40], [40.0 * 0.1; 3], 4, std::f64::consts::PI);
        assert_eq!(g.ngr3, 64000);
        assert_eq!(g.dr[0], 4.0 / 40.0);
        assert_eq!(g.dk[1], 2.0 * std::f64::consts::PI / 4.0);
        assert_eq!(g.center, [2.0; 3]);
        assert_eq!(g.index(1, 2, 3), 3 + 40 * (2 + 40));
        assert_eq!(slot(&[Complex64::new(1.0, 2.0)], 1), 2.0);
    }

    #[test]
    fn transforms_match_a_direct_dft_and_round_trip() {
        let g = FftGrid::new([3, 4, 5], [0.9, 1.2, 1.5], 1, std::f64::consts::PI);
        let mut fft = Fft3d::new(g);
        let r: Vec<Complex64> = (0..g.ngr3)
            .map(|m| Complex64::new((m as f64 * 0.37).sin(), (m as f64 * 0.11).cos()))
            .collect();
        let mut k = vec![Complex64::default(); g.ngr3];
        fft.forward(&r, &mut k);
        let fact = g.vol / g.ngr3 as f64;
        for (a, b, c) in [(0, 0, 0), (1, 2, 3), (2, 3, 4)] {
            let mut sum = Complex64::default();
            for i in 0..3 {
                for j in 0..4 {
                    for l in 0..5 {
                        let ph = -2.0
                            * std::f64::consts::PI
                            * ((a * i) as f64 / 3.0 + (b * j) as f64 / 4.0 + (c * l) as f64 / 5.0);
                        sum += r[g.index(i, j, l)] * Complex64::new(ph.cos(), ph.sin());
                    }
                }
            }
            let got = k[g.index(a, b, c)];
            assert!((got - sum * fact).norm() < 1e-12, "{got} vs {}", sum * fact);
        }
        fft.backward(&mut k);
        for (x, y) in k.iter().zip(&r) {
            assert!((x - y).norm() < 1e-13);
        }
    }
}
