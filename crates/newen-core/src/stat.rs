//! Running statistics with block-averaging error estimate.
//!
//! Faithful port of newen `gmath/Stat` (Oostenbrink et al.).
//! Implements the Allen–Tildesley block-averaging `ee()` described in
//! "Computer Simulation of Liquids", Oxford (1987), appendix D.
//!
//! # Usage
//! ```
//! use newen_core::stat::Stat;
//!
//! let mut s = Stat::new();
//! for x in [1.0_f64, 2.0, 3.0, 4.0, 5.0] { s.add(x); }
//! assert!((s.ave() - 3.0).abs() < 1e-12);
//! ```

/// Running statistics accumulator with block-averaging error estimate.
///
/// Stores all values in memory (same as newen Stat<T>).
/// Call `add()` to accumulate data, then query `ave()`, `rmsd()`, `ee()`.
#[derive(Debug, Clone, Default)]
pub struct Stat {
    vals: Vec<f64>,
}

impl Stat {
    /// Create an empty accumulator.
    pub fn new() -> Self {
        Self { vals: Vec::new() }
    }

    /// Append one observation.
    /// The values in insertion order.
    pub fn values(&self) -> &[f64] {
        &self.vals
    }

    /// gromos++ `Stat::ee` exactly: block averaging from block size 50 up; with fewer than 200
    /// values there is no block and gromos++ prints `nan`. [`Self::ee`] falls back to the naive
    /// standard error instead; programs that reproduce gromos++ output use this one.
    pub fn ee_strict(&self) -> f64 {
        let n = self.vals.len();
        if 4.0 * 50.0 >= n as f64 {
            return f64::NAN;
        }
        self.ee()
    }

    /// Record one value.
    pub fn add(&mut self, v: f64) {
        self.vals.push(v);
    }

    /// Number of observations.
    pub fn n(&self) -> usize {
        self.vals.len()
    }

    /// Arithmetic mean of all observations. Panics if empty.
    pub fn ave(&self) -> f64 {
        assert!(!self.vals.is_empty(), "Stat::ave() called on empty series");
        self.vals.iter().sum::<f64>() / self.vals.len() as f64
    }

    /// Mean of slice `[b, e)`.
    fn subave(&self, b: usize, e: usize) -> f64 {
        self.vals[b..e].iter().sum::<f64>() / (e - b) as f64
    }

    /// Mean-square deviation (biased variance): ⟨x²⟩ − ⟨x⟩².
    pub fn msd(&self) -> f64 {
        let n = self.vals.len() as f64;
        let mean = self.ave();
        let sq_mean: f64 = self.vals.iter().map(|&x| x * x).sum::<f64>() / n;
        (sq_mean - mean * mean).max(0.0)
    }

    /// Root-mean-square deviation: √(msd()).
    pub fn rmsd(&self) -> f64 {
        self.msd().sqrt()
    }

    /// gromos++ `Stat::msd` exactly: `<x²> − <x>²` with no clamp, so a constant series comes
    /// out very slightly negative and [`Self::rmsd_strict`] is then `NaN` — which is what
    /// gromos++ prints. [`Self::msd`] clamps at zero instead; programs that reproduce
    /// gromos++ output use these.
    pub fn msd_strict(&self) -> f64 {
        let n = self.vals.len() as f64;
        let mean = self.ave();
        let sq_mean: f64 = self.vals.iter().map(|&x| x * x).sum::<f64>() / n;
        sq_mean - mean * mean
    }

    /// √(`msd_strict`) — `NaN` for a constant series, as gromos++ prints it.
    pub fn rmsd_strict(&self) -> f64 {
        self.msd_strict().sqrt()
    }

    /// Minimum value.
    pub fn min(&self) -> f64 {
        self.vals.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Maximum value.
    pub fn max(&self) -> f64 {
        self.vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Block-averaging error estimate (Allen & Tildesley / newen algorithm).
    ///
    /// Partitions the series into blocks of increasing size, fits the RMSD of
    /// block averages linearly vs 1/block_size, and extrapolates to infinite
    /// block size.  Returns `rmsd / sqrt(N)` when the series is too short to
    /// form even one valid block (N < 200).
    pub fn ee(&self) -> f64 {
        let n = self.vals.len();
        if n < 2 {
            return 0.0;
        }

        let run_ave = self.ave();
        let run_rmsd = self.rmsd();
        if run_rmsd == 0.0 {
            return 0.0;
        }

        // Build block sizes: start at 50, multiply by 1.07177 until 4*blksz >= N
        let mut block_sizes: Vec<usize> = Vec::new();
        let mut blksz: f64 = 50.0;
        while 4.0 * blksz < n as f64 {
            let bs = blksz as usize;
            block_sizes.push(bs);
            let old = bs;
            while old == blksz as usize {
                blksz *= 1.07177;
            }
        }

        if block_sizes.is_empty() {
            // Too few data to do block averaging — fall back to naive SEM
            return run_rmsd / (n as f64).sqrt();
        }

        let nb = block_sizes.len();
        let mut fit = vec![0.0f64; nb];
        let mut x = vec![0.0f64; nb];

        for (j, &bs) in block_sizes.iter().enumerate() {
            let n_blocks = n / bs;
            let mut rmsd2 = 0.0f64;
            for i in 0..n_blocks {
                let bave = self.subave(i * bs, (i + 1) * bs);
                let d = bave - run_ave;
                rmsd2 += d * d;
            }
            rmsd2 /= n_blocks as f64;
            fit[j] = (bs as f64 * rmsd2) / (run_rmsd * run_rmsd);
            x[j] = 1.0 / bs as f64;
        }

        // Linear regression: fit = a*(1/bs) + b  →  extrapolate to 1/bs=0 → b
        let nb_f = nb as f64;
        let sx: f64 = x.iter().sum();
        let sf: f64 = fit.iter().sum();
        let sfx: f64 = x.iter().zip(fit.iter()).map(|(xi, fi)| xi * fi).sum();
        let sxx: f64 = x.iter().map(|xi| xi * xi).sum();

        let denom = sx * sx / nb_f - sxx;
        if denom.abs() < 1e-30 {
            return run_rmsd / (n as f64).sqrt();
        }
        let a = (sf * sx / nb_f - sfx) / denom;
        let b = (sf - a * sx) / nb_f;

        if b <= 0.0 {
            return run_rmsd / (n as f64).sqrt();
        }
        (b / n as f64).sqrt() * run_rmsd
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ave_rmsd() {
        let mut s = Stat::new();
        for x in [2.0_f64, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0] {
            s.add(x);
        }
        assert!((s.ave() - 5.0).abs() < 1e-12, "ave={}", s.ave());
        // msd = E[x²] - E[x]² = (4+16+16+16+25+25+49+81)/8 - 25 = 29 - 25 = 4
        assert!((s.msd() - 4.0).abs() < 1e-12, "msd={}", s.msd());
        assert!((s.rmsd() - 2.0).abs() < 1e-12, "rmsd={}", s.rmsd());
    }

    #[test]
    fn test_ee_constant_series() {
        let mut s = Stat::new();
        for _ in 0..1000 {
            s.add(3.0); // arbitrary constant — just needs to be non-zero and repeated
        }
        // All values identical → rmsd=0 → ee=0
        assert_eq!(s.ee(), 0.0);
    }

    #[test]
    fn test_ee_iid_normal_approx() {
        // i.i.d. standard normal: ee ≈ 1/sqrt(N)
        // Use a deterministic series with mean=0, std≈1 for reproducibility
        let mut s = Stat::new();
        let n = 500usize;
        for i in 0..n {
            // Simple deterministic pseudo-random-ish values with mean ~0
            let v = ((i as f64 * 1.6180339887) % 1.0) * 2.0 - 1.0;
            s.add(v);
        }
        let ee = s.ee();
        // Should be O(1/sqrt(N)) ≈ 0.045; just check it's finite and positive
        assert!(ee > 0.0, "ee={ee}");
        assert!(ee < 1.0, "ee={ee} unexpectedly large");
    }

    #[test]
    fn test_min_max() {
        let mut s = Stat::new();
        for x in [3.0_f64, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0] {
            s.add(x);
        }
        assert_eq!(s.min(), 1.0);
        assert_eq!(s.max(), 9.0);
    }
}
