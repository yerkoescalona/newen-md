//! GSL-compatible random number generation.
//!
//! GROMOS generates initial velocities (NTIVEL=1) using GSL's default
//! generator, `mt19937`, combined with `gsl_ran_gaussian` (the Marsaglia
//! polar / Box-Muller method). To reproduce GROMOS trajectories bit-for-bit
//! we re-implement those two algorithms exactly rather than relying on a
//! generic RNG crate.
//!
//! Refs: GROMOS `math/random.{h,cc}` (`RandomGeneratorGSL`),
//! `util/generate_velocities.cc`.

use crate::math::Vec3;

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

/// Mersenne Twister (MT19937), seeded and stepped exactly like GSL's
/// `gsl_rng_mt19937` (`rng/mt.c`).
#[derive(Debug, Clone)]
pub struct GslMt19937 {
    mt: [u32; N],
    mti: usize,
}

impl GslMt19937 {
    /// Seed the generator. GSL maps a seed of 0 to its default seed (4357),
    /// matching `mt19937_set`.
    pub fn new(seed: u32) -> Self {
        let s = if seed == 0 { 4357 } else { seed };
        let mut mt = [0u32; N];
        mt[0] = s;
        for i in 1..N {
            mt[i] = (1_812_433_253u32.wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30)))
                .wrapping_add(i as u32);
        }
        GslMt19937 { mt, mti: N }
    }

    fn next_u32(&mut self) -> u32 {
        if self.mti >= N {
            for kk in 0..N - M {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] = self.mt[kk + M] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            for kk in N - M..N - 1 {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] =
                    self.mt[kk - (N - M)] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            let y = (self.mt[N - 1] & UPPER_MASK) | (self.mt[0] & LOWER_MASK);
            self.mt[N - 1] = self.mt[M - 1] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };

            self.mti = 0;
        }

        let mut k = self.mt[self.mti];
        k ^= k >> 11;
        k ^= (k << 7) & 0x9d2c_5680;
        k ^= (k << 15) & 0xefc6_0000;
        k ^= k >> 18;
        self.mti += 1;
        k
    }

    /// Uniform deviate in [0, 1), matching `gsl_rng_uniform` for `mt19937`
    /// (`mt19937_get_double`: 32-bit draw scaled by 1/2^32).
    pub fn uniform(&mut self) -> f64 {
        self.next_u32() as f64 * (1.0 / 4_294_967_296.0)
    }

    /// Uniform deviate in (0, 1), matching `gsl_rng_uniform_pos` (re-draws on 0).
    pub fn uniform_pos(&mut self) -> f64 {
        loop {
            let x = self.uniform();
            if x != 0.0 {
                return x;
            }
        }
    }
}

/// Gaussian deviate with mean 0 and standard deviation `sigma`, matching
/// `gsl_ran_gaussian` (Marsaglia polar method): draw `(x, y)` uniformly in the
/// unit disk and return `sigma * y * sqrt(-2 ln(r^2) / r^2)`, discarding `x`.
pub fn gsl_ran_gaussian(rng: &mut GslMt19937, sigma: f64) -> f64 {
    loop {
        let x = -1.0 + 2.0 * rng.uniform_pos();
        let y = -1.0 + 2.0 * rng.uniform_pos();
        let r2 = x * x + y * y;
        if r2 <= 1.0 && r2 != 0.0 {
            return sigma * y * (-2.0 * r2.ln() / r2).sqrt();
        }
    }
}

/// Generate initial velocities from the Maxwell-Boltzmann distribution at
/// temperature `temperature`, one Gaussian draw per Cartesian component with
/// `sigma = sqrt(k_B * T / m)`.
///
/// Mirrors GROMOS `util::generate_velocities`: both the "current" and "old"
/// velocity arrays are filled with the same freshly generated values (GROMOS
/// stores the same draw into `old_vel` and `vel`).
///
/// Ref: `util/generate_velocities.cc`. Uses GROMOS's `math::k_Boltzmann`
/// (`math/math.cc`), not the CODATA value used elsewhere in this codebase, so
/// that the generated trajectory matches the reference bit-for-bit.
pub fn generate_velocities(temperature: f64, seed: u32, masses: &[f64]) -> Vec<Vec3> {
    let mut rng = GslMt19937::new(seed);
    masses
        .iter()
        .map(|&mass| {
            let sd = (crate::units::kB * temperature / mass).sqrt();
            Vec3::new(
                gsl_ran_gaussian(&mut rng, sd),
                gsl_ran_gaussian(&mut rng, sd),
                gsl_ran_gaussian(&mut rng, sd),
            )
        })
        .collect()
}
