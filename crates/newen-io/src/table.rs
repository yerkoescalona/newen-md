//! Numeric column files — the time series `ene_ana`, `tser` and friends write and the free-energy
//! programs read back: whitespace-separated numbers, `#` starts a comment, blank lines skipped.
//! gromos++ reads these with `istringstream` line by line; a line that does not parse is an error.

use std::path::Path;

use crate::IoError;

/// Every data row of `path` as numbers. Rows may differ in length; the caller checks what it needs.
pub fn read_columns<P: AsRef<Path>>(path: P) -> Result<Vec<Vec<f64>>, IoError> {
    let text = std::fs::read_to_string(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    parse_columns(&text, &path.as_ref().display().to_string())
}

/// [`read_columns`] on text already in memory.
pub fn parse_columns(text: &str, what: &str) -> Result<Vec<Vec<f64>>, IoError> {
    let mut rows = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let data = line.split('#').next().unwrap_or("").trim();
        if data.is_empty() {
            continue;
        }
        let row: Result<Vec<f64>, _> = data.split_whitespace().map(|t| t.parse::<f64>()).collect();
        match row {
            Ok(r) => rows.push(r),
            Err(_) => {
                return Err(IoError::ParseError(format!(
                    "{what}, line {}: not a numeric row: {data:?}",
                    n + 1
                )))
            },
        }
    }
    Ok(rows)
}

/// C++ default stream formatting of a double with `precision` significant digits (`%g`):
/// fixed when the exponent is in [−5, precision), scientific otherwise, trailing zeros dropped.
pub fn cpp_g(x: f64, precision: usize) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    if !x.is_finite() {
        return if x.is_nan() {
            "nan".into()
        } else if x > 0.0 {
            "inf".into()
        } else {
            "-inf".into()
        };
    }
    let p = precision.max(1);
    // exponent after rounding to p significant digits
    let sci = format!("{:.*e}", p - 1, x);
    let (mant, ex) = sci.split_once('e').unwrap();
    let exp: i32 = ex.parse().unwrap();
    if exp < -4 || exp >= p as i32 {
        let m = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        return format!("{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs());
    }
    let decimals = (p as i32 - 1 - exp).max(0) as usize;
    let t = format!("{:.*}", decimals, x);
    if t.contains('.') {
        t.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        t
    }
}

/// C++ `std::fixed` with `precision` decimals (`%.*f`): ties round half to even, as glibc does
/// on the exact binary value (`0.125` → `0.12`); NaN prints `nan`/`-nan` by its sign bit.
pub fn cpp_fixed(x: f64, precision: usize) -> String {
    if x.is_nan() {
        return if x.is_sign_negative() { "-nan" } else { "nan" }.to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    format!("{x:.precision$}")
}

/// C++ `std::scientific` with `precision` digits: mantissa, `e`, signed two-digit exponent
/// (`1.57000e+07`).
pub fn cpp_e(x: f64, precision: usize) -> String {
    let s = format!("{:.*e}", precision, x);
    let (mantissa, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    format!(
        "{mantissa}e{}{:02}",
        if exp < 0 { '-' } else { '+' },
        exp.abs()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comments_and_blank_lines_are_skipped() {
        let rows = parse_columns("# t  x\n\n 0.0 1.5 # trailing\n0.002 2.5\n", "test").unwrap();
        assert_eq!(rows, vec![vec![0.0, 1.5], vec![0.002, 2.5]]);
    }

    #[test]
    fn fixed_rounds_ties_to_even_like_glibc() {
        assert_eq!(cpp_fixed(0.125, 2), "0.12");
        assert_eq!(cpp_fixed(0.375, 2), "0.38");
        assert_eq!(cpp_fixed(2.5, 0), "2");
        assert_eq!(cpp_fixed(-0.125, 2), "-0.12");
        assert_eq!(cpp_fixed(0.125, 3), "0.125");
        assert_eq!(cpp_fixed(-f64::NAN, 6), "-nan");
        assert_eq!(cpp_fixed(f64::INFINITY, 6), "inf");
    }

    #[test]
    fn a_word_is_an_error() {
        assert!(parse_columns("0.0 abc\n", "test").is_err());
    }
}
