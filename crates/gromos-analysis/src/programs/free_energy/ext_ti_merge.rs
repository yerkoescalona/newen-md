//! ext_ti_merge — combine the `ext_ti_ana` predictions of several simulations into one
//! ⟨∂H/∂λ⟩ curve and integrate it (gromos++ `ext_ti_merge`, `programs/ext_ti_merge.cc`).
//!
//! Between two neighbouring simulated λ, the prediction of each at the same λ_P is weighted
//! linearly by its distance from λ_P; the last simulated λ is appended as is. The curve is
//! integrated with the trapezoid rule, the error by gromos++'s propagation.
//!
//! Reproduced from gromos++ as it is: one `no_errors` flag for all files (the first two-column
//! line switches errors off everywhere), and a third column read as the error even when
//! `ext_ti_ana @countframes` wrote a frame count there.
//!
//! Deliberate differences from gromos++ (item 1.8):
//! - The error of the integral. gromos++ 9cd9f6b3's sum skips the point N−2
//!   (`ext_ti_merge.cc:282`, `i < N-2`). Corrected by default since 0.0.57; `@parity gromospp`
//!   ([`Parity::GromosppReleased`]) reproduces the defect. The oracle of the corrected output is
//!   gromos++ 9cd9f6b3 patched with `D8_ext_ti_merge_loop.patch`
//!   (`crates/gromos-md/tests/gromosXX_references/patches/`), that of the reproduced one the
//!   released build. [`Merged::integral_error_corrected`] includes every point in both modes.
//! - With fewer than two points the error is 0 (gromos++ reads past the end); no points at all
//!   is an error (gromos++ crashes).

use super::common::{get_value, parse_finite};
use super::Parity;
use crate::args::Arguments;
use crate::programs::io_err;
use gromos_io::args::{Program, Status};
use gromos_io::table::{cpp_fixed, cpp_g};
use std::collections::BTreeMap;
use std::io::Write;

pub const USAGE: &str = "# ext_ti_merge
\t@files            <data files (ext_ti_ana output)>
\t@noerrors         <do not read and use the error column which might be in the files>
\t[@slam            <lambdas of the simulation, optional if found in the header of the data files after '#SLAM '>]
\t[@lam_precision   <lambda value precision in outfiles, default: 2>]
\t[@parity          <gromospp: the integral's error as gromos++ 9cd9f6b3 prints it (default: corrected)>]";

pub const KNOWN: &[&str] = &[
    "files",
    "slam",
    "noerrors",
    "weightfunction",
    "lam_precision",
    "parity",
];

pub const PROGRAM: Program = Program {
    name: "ext_ti_merge",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// One `predicted_l<slam>.dat`: λ_P, ⟨∂H/∂λ⟩ and (third column) its error.
#[derive(Debug, Clone, PartialEq)]
pub struct Predicted {
    pub slam: f64,
    pub plam: Vec<f64>,
    pub dhdl: Vec<f64>,
    pub error: Vec<f64>,
}

/// The merged curve and its integral (kJ/mol).
#[derive(Debug, Clone, PartialEq)]
pub struct Merged {
    pub lambda: Vec<f64>,
    pub dhdl: Vec<f64>,
    pub error: Option<Vec<f64>>,
    pub integral: f64,
    /// The printed error of the integral: [`Merged::integral_error_corrected`], or under
    /// [`Parity::GromosppReleased`] gromos++'s (the point N−2 left out).
    pub integral_error: Option<f64>,
    /// The error of the integral with every point, whatever the parity.
    pub integral_error_corrected: Option<f64>,
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    evaluate(argv, out).map(|_| ())
}

/// [`run`], returning the curve it printed.
pub fn evaluate(argv: &[String], out: &mut dyn Write) -> Result<Merged, String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let mut no_errors = args.has("noerrors");
    let parity = Parity::from_args(&args)?;
    if args.count("files") <= 1 {
        return Err(format!("At least two data files required!\n{USAGE}"));
    }
    let weightfunction: String = get_value(&args, "weightfunction", Some("default".into()))?;
    let given_slam: Option<Vec<f64>> = if args.count("slam") > 0 {
        if args.count("slam") != args.count("files") {
            return Err("@slam needs to provide one lambda value for each given file".into());
        }
        Some(
            args.values("slam")
                .iter()
                .map(|v| parse_finite(v).ok_or_else(|| format!("@slam: cannot parse '{v}'")))
                .collect::<Result<_, _>>()?,
        )
    } else {
        None
    };
    let lam_precision: i64 = get_value(&args, "lam_precision", Some(2))?;
    if lam_precision < 0 {
        return Err("@lam_precision cannot be negative".into());
    }
    let mut files = Vec::new();
    for (k, path) in args.values("files").iter().enumerate() {
        let text =
            std::fs::read_to_string(path).map_err(|_| format!("Could not open file '{path}'"))?;
        let slam = given_slam.as_ref().map(|s| s[k]);
        files.push(
            parse_predicted(&text, slam, &mut no_errors).map_err(|e| format!("{path}: {e}"))?,
        );
    }
    let slam: Vec<f64> = files.iter().map(|f| f.slam).collect();
    let mut head = String::from("# Interpolating from data simulated at lambdas: ");
    for &s in &slam {
        head.push_str(&cpp_g(s, 6));
        head.push(' ');
    }
    head.push('\n');
    out.write_all(head.as_bytes()).map_err(io_err)?;
    if weightfunction != "default" {
        return Err(format!("Unknown weightfunction {weightfunction}"));
    }
    let merged = merge(&files, no_errors, parity)?;
    let prec = lam_precision as usize;
    let mut text = String::new();
    for (k, (&l, &v)) in merged.lambda.iter().zip(&merged.dhdl).enumerate() {
        text.push_str(&format!("{:>4} {}", cpp_fixed(l, prec), cpp_fixed(v, 6)));
        if let Some(e) = &merged.error {
            text.push_str(&format!(" {}", cpp_fixed(e[k], 6)));
        }
        text.push('\n');
    }
    text.push_str(&format!("# integral: {}", cpp_fixed(merged.integral, 6)));
    if let Some(e) = merged.integral_error {
        text.push_str(&format!(" +/- {}", cpp_fixed(e, 6)));
    }
    text.push('\n');
    out.write_all(text.as_bytes()).map_err(io_err)?;
    Ok(merged)
}

/// Read one `predicted_l<slam>.dat` (`ext_ti_merge.cc:180-224`). `slam` from `@slam`, else the
/// file's first `#SLAM` line. `no_errors` is shared by all files: a line with only two columns
/// sets it. Every line but a last one without a newline is read; a blank line is an error.
pub fn parse_predicted(
    text: &str,
    slam: Option<f64>,
    no_errors: &mut bool,
) -> Result<Predicted, String> {
    let read_slam = slam.is_some();
    let mut file_slam = slam;
    let mut f = Predicted {
        slam: 0.0,
        plam: Vec::new(),
        dhdl: Vec::new(),
        error: Vec::new(),
    };
    let parts: Vec<&str> = text.split('\n').collect();
    let mut read_slam_file = false;
    for raw in &parts[..parts.len() - 1] {
        let line = raw.trim();
        if !read_slam && !read_slam_file && line.contains("#SLAM") {
            let v = line
                .split_whitespace()
                .nth(1)
                .and_then(parse_finite)
                .ok_or_else(|| format!("cannot read the lambda after '#SLAM' in '{line}'"))?;
            file_slam = Some(v);
            read_slam_file = true;
        } else if !line.starts_with('#') {
            let tok: Vec<&str> = line.split_whitespace().collect();
            let (Some(p), Some(d)) = (
                tok.first().and_then(|t| parse_finite(t)),
                tok.get(1).and_then(|t| parse_finite(t)),
            ) else {
                return Err(
                    "can't read this file format, all lines not starting with a '#' \
                            should contain 2 or (with errors) 3 columns of floats"
                        .into(),
                );
            };
            f.plam.push(p);
            f.dhdl.push(d);
            if tok.len() > 2 && !*no_errors {
                match parse_finite(tok[2]) {
                    Some(e) => f.error.push(e),
                    None => *no_errors = true,
                }
            } else {
                *no_errors = true;
            }
        }
    }
    f.slam = file_slam.ok_or(
        "You have to provide the simulated lambda values either with @slam or in the file's \
         header after '#SLAM '",
    )?;
    Ok(f)
}

/// gromos++ `weight_function` and the integrals (`ext_ti_merge.cc:266-344`); `parity` picks
/// [`Merged::integral_error`].
pub fn merge(files: &[Predicted], no_errors: bool, parity: Parity) -> Result<Merged, String> {
    // one entry per simulated λ (a later file with the same λ replaces an earlier one)
    let mut by_slam: BTreeMap<u64, &Predicted> = BTreeMap::new();
    let key = |x: f64| x.to_bits();
    for f in files {
        by_slam.insert(key(f.slam), f);
    }
    let mut slam: Vec<f64> = files.iter().map(|f| f.slam).collect();
    slam.sort_by(|a, b| a.total_cmp(b));
    let (mut lam, mut val, mut err) = (Vec::new(), Vec::new(), Vec::new());
    for w in slam.windows(2) {
        let (sl1, sl2) = (w[0], w[1]);
        if sl1 == sl2 {
            return Err("Two slam values are the same!\n".into());
        }
        let pick = |f: &Predicted| -> Vec<(f64, f64, f64)> {
            (0..f.plam.len())
                .filter(|&j| f.plam[j] >= sl1 && f.plam[j] < sl2)
                .map(|j| (f.plam[j], f.dhdl[j], f.error.get(j).copied().unwrap_or(0.0)))
                .collect()
        };
        let a = pick(by_slam[&key(sl1)]);
        let b = pick(by_slam[&key(sl2)]);
        if a.iter().map(|t| t.0).ne(b.iter().map(|t| t.0)) {
            eprintln!(
                "ext_ti_merge: predicted lambdas between {} and {} are not the same, skipping.",
                cpp_g(sl1, 6),
                cpp_g(sl2, 6)
            );
            continue;
        }
        if a.is_empty() {
            eprintln!(
                "ext_ti_merge: no predicted lambdas between {} and {}.",
                cpp_g(sl1, 6),
                cpp_g(sl2, 6)
            );
            continue;
        }
        for (&(p, d0, e0), &(_, d1, e1)) in a.iter().zip(&b) {
            let ws1 = (p - sl2) / (sl1 - sl2);
            let ws2 = (p - sl1) / (sl2 - sl1);
            lam.push(p);
            val.push(d0 * ws1 + d1 * ws2);
            err.push((e0 * e0 * ws1 * ws1 + e1 * e1 * ws2 * ws2).sqrt());
        }
    }
    let last = *slam.last().unwrap();
    let lf = by_slam[&key(last)];
    if let Some(j) = lf.plam.iter().position(|&p| p == last) {
        lam.push(lf.plam[j]);
        val.push(lf.dhdl[j]);
        err.push(lf.error.get(j).copied().unwrap_or(0.0));
    }
    if lam.is_empty() {
        return Err("no predicted lambdas to integrate".into());
    }
    let integral = integrate(&lam, &val);
    let (error, integral_error, integral_error_corrected) = if no_errors {
        (None, None, None)
    } else {
        let c = integrate_error(&lam, &err, true);
        let e = match parity {
            Parity::Corrected => c,
            Parity::GromosppReleased => integrate_error(&lam, &err, false),
        };
        (Some(err), Some(e), Some(c))
    };
    Ok(Merged {
        lambda: lam,
        dhdl: val,
        error,
        integral,
        integral_error,
        integral_error_corrected,
    })
}

/// The trapezoid rule, in gromos++'s operation order.
pub fn integrate(x: &[f64], y: &[f64]) -> f64 {
    let mut integral = 0.0;
    for i in 0..x.len().saturating_sub(1) {
        integral += 0.5 * (y[i + 1] + y[i]) * (x[i + 1] - x[i]);
    }
    integral
}

/// The error of [`integrate`] from independent point errors `y`:
/// √(Σ w_i² y_i²) with the trapezoid weights. `all_points = false` is gromos++ 9cd9f6b3's loop
/// (`i < N-2`, so the point N−2 is missing); fewer than two points give 0.
pub fn integrate_error(x: &[f64], y: &[f64], all_points: bool) -> f64 {
    let n = x.len();
    if n < 2 {
        return 0.0;
    }
    let mut integral = 0.25 * (x[1] - x[0]) * (x[1] - x[0]) * y[0] * y[0];
    let end = if all_points {
        n - 1
    } else {
        n.saturating_sub(2)
    };
    for i in 1..end.max(1) {
        integral += y[i] * y[i] * (x[i + 1] - x[i - 1]) * (x[i + 1] - x[i - 1]) / 4.0;
    }
    integral += 0.25 * (x[n - 1] - x[n - 2]) * (x[n - 1] - x[n - 2]) * y[n - 1] * y[n - 1];
    integral.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_integral_skips_the_second_to_last_point_like_gromospp() {
        let x = [0.0, 0.5, 1.0];
        let y = [1.0, 2.0, 1.0];
        // ends only: 0.25·0.25·1 twice
        assert!((integrate_error(&x, &y, false) - (2.0f64 * 0.0625).sqrt()).abs() < 1e-15);
        // with the middle point: + 4·1/4
        assert!((integrate_error(&x, &y, true) - (0.125f64 + 1.0).sqrt()).abs() < 1e-15);
        assert_eq!(integrate_error(&[0.0], &[1.0], false), 0.0);
        assert!((integrate(&x, &[0.0, 1.0, 2.0]) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn two_columns_switch_errors_off_for_every_file() {
        let mut no_errors = false;
        let a = parse_predicted("#SLAM 0\n0.0 1.0 0.1\n0.5 2.0 0.2\n", None, &mut no_errors);
        assert_eq!(a.unwrap().error, vec![0.1, 0.2]);
        assert!(!no_errors);
        let b = parse_predicted("#SLAM 1\n0.5 2.0\n1.0 3.0 0.3\n", None, &mut no_errors).unwrap();
        assert!(no_errors);
        assert!(b.error.is_empty());
        // an unterminated last line is not read
        let c = parse_predicted("#SLAM 1\n0.5 2.0\n1.0 3.0", None, &mut no_errors).unwrap();
        assert_eq!(c.plam, vec![0.5]);
    }
}
