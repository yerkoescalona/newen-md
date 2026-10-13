//! What `ext_ti_ana`, `bar`, `ext_ti_merge` and `mbar` share: gromos++'s argument getters, its
//! random numbers, its `std::max`/`std::min`, and the [`Parity`] switch.

use gromos_io::args::Arguments;
use std::str::FromStr;

/// What `bar`, `ext_ti_merge` and `ext_ti_ana` print where gromos++ 9cd9f6b3 has a defect
/// (item 1.8): the corrected value (the default), or the released gromos++ value
/// (`@parity gromospp`), kept for the byte-for-byte parity tests against that release.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Parity {
    /// The defect corrected, as the patched gromos++ prints it.
    #[default]
    Corrected,
    /// The defect reproduced, as gromos++ 9cd9f6b3 prints it.
    GromosppReleased,
}

impl Parity {
    /// The program argument (`@parity gromospp`).
    pub const KEY: &'static str = "parity";

    /// `gromospp` or `corrected`.
    pub fn parse(value: &str) -> Result<Parity, String> {
        match value {
            "gromospp" => Ok(Parity::GromosppReleased),
            "corrected" => Ok(Parity::Corrected),
            v => Err(format!(
                "@parity: '{v}' is neither 'gromospp' (reproduce gromos++ 9cd9f6b3) nor \
                 'corrected' (the default)"
            )),
        }
    }

    /// `@parity` of a program's arguments: absent is [`Parity::Corrected`]; exactly one value.
    pub(crate) fn from_args(args: &Arguments) -> Result<Parity, String> {
        if !args.has(Self::KEY) {
            return Ok(Parity::Corrected);
        }
        match args.values(Self::KEY) {
            [v] => Parity::parse(v),
            _ => Err("@parity takes one value: gromospp or corrected".into()),
        }
    }
}

/// gromos++ `Arguments::getValue<T>(key, required, default)` (`Arguments.h:196-221`): absent →
/// the default (an error when there is none); given without a value, or with one that does not
/// parse → "invalid value"; further values are ignored.
pub(crate) fn get_value<T: FromStr>(
    args: &Arguments,
    key: &str,
    default: Option<T>,
) -> Result<T, String> {
    if !args.has(key) {
        return default.ok_or_else(|| format!("Argument '{key}' is required."));
    }
    let v = args.values(key).first().map_or("", String::as_str);
    v.trim()
        .parse::<T>()
        .map_err(|_| format!("Argument '{key}' has an invalid value: cannot interpret '{v}'."))
}

/// gromos++ `Arguments::getValues<T>(key, num, required, default)` (`Arguments.h:224-258`):
/// absent → the default; given → exactly `num` values.
pub(crate) fn get_values<T: FromStr>(
    args: &Arguments,
    key: &str,
    num: usize,
    default: Option<Vec<T>>,
) -> Result<Vec<T>, String> {
    if !args.has(key) {
        return default.ok_or_else(|| format!("Argument '{key}' is required ({num} values)."));
    }
    let vals = args.values(key);
    if vals.len() < num {
        return Err(format!("Not enough values for argument '{key}' given."));
    }
    if vals.len() > num {
        return Err(format!("Too many arguments for argument '{key}'."));
    }
    vals.iter()
        .enumerate()
        .map(|(i, v)| {
            v.parse::<T>().map_err(|_| {
                format!(
                    "Argument '{key}' has an invalid value (value number {}): cannot interpret '{v}'.",
                    i + 1
                )
            })
        })
        .collect()
}

/// `std::max(a, b)`: `a < b ? b : a` (differs from `f64::max` only on NaN).
#[inline]
pub(crate) fn cmax(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

/// `std::min(a, b)`: `b < a ? b : a`.
#[inline]
pub(crate) fn cmin(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

/// C `pow(x, y)` through libm, never rewritten by the optimiser into multiplications (the
/// gromos++ oracle build calls `pow` for `pow(slam, 4)` too).
#[inline]
pub(crate) fn cpow(x: f64, y: f64) -> f64 {
    x.powf(std::hint::black_box(y))
}

/// glibc `RAND_MAX`.
pub(crate) const RAND_MAX: f64 = 2147483647.0;

/// glibc `rand_r` (`stdlib/rand_r.c`): three steps of the ANSI LCG combined into 31 bits, so a
/// seed gives gromos++'s bootstrap draws exactly.
pub(crate) fn rand_r(seed: &mut u32) -> u32 {
    let mut next = *seed;
    next = next.wrapping_mul(1103515245).wrapping_add(12345);
    let mut result = (next / 65536) % 2048;
    next = next.wrapping_mul(1103515245).wrapping_add(12345);
    result <<= 10;
    result ^= (next / 65536) % 1024;
    next = next.wrapping_mul(1103515245).wrapping_add(12345);
    result <<= 10;
    result ^= (next / 65536) % 1024;
    *seed = next;
    result
}

/// gromos++'s draw `int((r * length) / RAND_MAX)`; `r == RAND_MAX` (index `length`, out of
/// bounds there) is kept in range.
pub(crate) fn draw_index(seed: &mut u32, length: usize) -> usize {
    let r = rand_r(seed) as f64;
    (((r * length as f64) / RAND_MAX) as usize).min(length.saturating_sub(1))
}

/// gromos++'s bootstrap seed base, `time(NULL)`.
pub(crate) fn wall_clock_seed() -> u32 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs() as u32)
}

/// gromos++ `Stat::covariance(X, X)` (`Stat.cc:181-196`): the population variance, in its
/// operation order.
pub(crate) fn population_variance(x: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mut sum, mut ssum) = (0.0, 0.0);
    for &v in x {
        sum += v;
        ssum += v * v;
    }
    let ave = sum / n;
    ssum /= n;
    ssum - ave * (sum / n)
}

/// One `double` of a C++ `istringstream >> double`: libstdc++ refuses `nan` and `inf`, which
/// Rust's parser accepts.
pub(crate) fn parse_finite(token: &str) -> Option<f64> {
    token.parse::<f64>().ok().filter(|v| v.is_finite())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rand_r_is_glibcs() {
        // glibc: unsigned s = 1; rand_r(&s) three times
        let mut s = 1u32;
        let got: Vec<u32> = (0..3).map(|_| rand_r(&mut s)).collect();
        assert_eq!(got, vec![476707713, 1186278907, 505671508]);
    }
}
