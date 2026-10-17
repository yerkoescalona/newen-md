//! bar — free-energy differences between neighbouring λ by Bennett's acceptance ratio
//! (gromos++ `bar`, `programs/bar.cc`).
//!
//! One file per simulated λ, as `ext_ti_ana @bar_data` writes them: the frame energies at the
//! simulated λ and at other λ. The files are sorted by their λ and each neighbouring pair gives
//! ΔG by the self-consistent BAR iteration in kT units (Shirts et al., PRL 91 (2003) 140601),
//! with gromos++'s analytical error, an optional bootstrap error and the overlap integral of the
//! two ΔE distributions (`@printdist` writes them to `dist_<k>_<k+1>.dat`).
//!
//! Deliberate differences from gromos++ (item 1.8):
//! - The analytical error. gromos++ 9cd9f6b3 divides each side's relative variance by an extra
//!   (ln⟨f⟩)² (`bar.cc:639`; dfmult has no such factor, `dfmult.cc:212-217`). Corrected by
//!   default since 0.0.57;`@parity gromospp` ([`Parity::GromosppReleased`]) reproduces the
//!   defect. The oracle of the corrected output is gromos++ 9cd9f6b3 patched with
//!   `D4_bar_error.patch` (`crates/newen-md/tests/gromosXX_references/patches/`), that of the
//!   reproduced one the released build. [`Interval::error_corrected`] is the corrected error in
//!   both modes.
//! - The bootstrap seed can be given (library); gromos++ takes `time(NULL)`.

use super::common::{
    cmax, cmin, draw_index, get_value, parse_finite, population_variance, wall_clock_seed,
};
use super::Parity;
use crate::args::Arguments;
use crate::distribution::Distribution;
use crate::lnexp::{lnexp_ave, lnexp_covariance, lnexp_stat_ineff};
use crate::programs::io_err;
use newen_io::args::{Program, Status};
use newen_io::table::{cpp_fixed, cpp_g};
use std::io::Write;
use std::path::{Path, PathBuf};

pub const USAGE: &str = "# bar
\t@files            <data files (ext_ti_ana bar_data output)>
\t@temp             <temperature>
\t[@maxiterations   <maximum number of iterations (default: 500)>]
\t[@convergence     <relative change in free energy for convergence (default: 1E-5)>]
\t[@printdist       <write out distributions>]
\t[@bootstrap       <number of bootstrap estimates for error estimates (default: 0)>]
\t[@parity          <gromospp: the error as gromos++ 9cd9f6b3 prints it (default: corrected)>]";

pub const KNOWN: &[&str] = &[
    "files",
    "maxiterations",
    "convergence",
    "bootstrap",
    "temp",
    "printdist",
    "parity",
];

pub const PROGRAM: Program = Program {
    name: "bar",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// One `ext_ti_ana @bar_data` file: the simulated λ, the other λ, and per frame the energy at
/// the simulated λ (`e_s`) and at each other λ (`e_p[k]`), kJ/mol.
#[derive(Debug, Clone, PartialEq)]
pub struct BarFile {
    pub slam: f64,
    pub plam: Vec<f64>,
    pub e_s: Vec<f64>,
    pub e_p: Vec<Vec<f64>>,
}

/// How to run the estimator.
#[derive(Debug, Clone, PartialEq)]
pub struct BarOptions {
    pub temp: f64,
    pub max_iterations: usize,
    /// Relative change of ΔG between iterations that ends them.
    pub convergence: f64,
    pub bootstrap: usize,
    /// Bootstrap seed base (interval `k` starts at `seed + k`); `None` is the wall clock.
    pub seed: Option<u32>,
    /// Where `dist_<k>_<k+1>.dat` go (`@printdist`: the working directory); `None` writes none.
    pub printdist: Option<PathBuf>,
    /// `@parity`: which analytical error [`Interval::error`] is.
    pub parity: Parity,
}

/// One neighbouring pair, kJ/mol.
#[derive(Debug, Clone, PartialEq)]
pub struct Interval {
    pub lam_i: f64,
    pub lam_j: f64,
    /// G(λ_j) − G(λ_i).
    pub dg: f64,
    /// The analytical error printed: [`Interval::error_corrected`], or under
    /// [`Parity::GromosppReleased`] gromos++'s (with the extra 1/(ln⟨f⟩)², `bar.cc:639`).
    pub error: f64,
    /// The analytical error √(σ²⟨f_i⟩/⟨f_i⟩² + σ²⟨f_j⟩/⟨f_j⟩²)·kT, whatever the parity.
    pub error_corrected: f64,
    /// The bootstrap standard deviation (`@bootstrap`).
    pub bootstrap: Option<f64>,
    /// Overlap integral of the two ΔE distributions (100 bins).
    pub overlap: f64,
    pub iterations: usize,
}

/// All intervals and their sums (errors added in quadrature).
#[derive(Debug, Clone, PartialEq)]
pub struct BarResult {
    pub intervals: Vec<Interval>,
    pub dg_total: f64,
    /// The printed total error (from [`Interval::error`]).
    pub error_total: f64,
    /// The total from [`Interval::error_corrected`].
    pub error_total_corrected: f64,
    pub bootstrap_total: Option<f64>,
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    evaluate(argv, out).map(|_| ())
}

/// [`run`], returning the numbers it printed.
pub fn evaluate(argv: &[String], out: &mut dyn Write) -> Result<BarResult, String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let max_iterations: i64 = get_value(&args, "maxiterations", Some(500))?;
    let convergence: f64 = get_value(&args, "convergence", Some(1e-5))?;
    let temp: f64 = get_value(&args, "temp", None)?;
    let printdist = args.has("printdist");
    let bootstrap: i64 = get_value(&args, "bootstrap", Some(0))?;
    let parity = Parity::from_args(&args)?;
    if args.count("files") <= 1 {
        return Err(format!("At least two data files required!\n{USAGE}"));
    }
    if max_iterations < 0 || bootstrap < 0 {
        return Err("@maxiterations and @bootstrap cannot be negative".into());
    }
    let files = args
        .values("files")
        .iter()
        .map(|p| read_bar_file(p))
        .collect::<Result<Vec<_>, _>>()?;
    let opts = BarOptions {
        temp,
        max_iterations: max_iterations as usize,
        convergence,
        bootstrap: bootstrap as usize,
        seed: None,
        printdist: printdist.then(|| PathBuf::from(".")),
        parity,
    };
    bar(&files, &opts, out)
}

/// Read a `bar_data` file the way gromos++ does (`bar.cc:241-288`): lines with `#` in column 0
/// are skipped; the first data line is the number of other λ, the second the simulated λ and
/// those λ, every following line a frame. A last line without a newline is not read; a short,
/// blank or non-numeric line (`nan`/`inf` included) is an error.
pub fn read_bar_file(path: &str) -> Result<BarFile, String> {
    let text =
        std::fs::read_to_string(path).map_err(|_| format!("Could not open file '{path}'"))?;
    let f = parse_bar_file(&text).map_err(|e| format!("{path}: {e}"))?;
    eprintln!(
        "read file {path} with {} plam values and {} data_points",
        f.plam.len(),
        f.e_s.len()
    );
    Ok(f)
}

/// C++ `std::getline` over a whole text: each line and whether the stream hit EOF reading it.
struct Lines<'a> {
    parts: Vec<&'a str>,
    pos: usize,
}

impl<'a> Lines<'a> {
    fn new(text: &'a str) -> Self {
        Lines {
            parts: text.split('\n').collect(),
            pos: 0,
        }
    }

    fn getline(&mut self) -> (&'a str, bool) {
        if self.pos >= self.parts.len() {
            return ("", true);
        }
        let line = self.parts[self.pos];
        self.pos += 1;
        (line, self.pos == self.parts.len())
    }

    /// `do getline while (line[0] == '#')`; a comment as the unterminated last line ends the
    /// text (gromos++ spins on it forever).
    fn data(&mut self) -> (&'a str, bool) {
        loop {
            let (line, eof) = self.getline();
            if !line.starts_with('#') {
                return (line, eof);
            }
            if eof {
                return ("", true);
            }
        }
    }
}

/// [`read_bar_file`] on text in memory.
pub fn parse_bar_file(text: &str) -> Result<BarFile, String> {
    let mut lines = Lines::new(text);
    let (line, _) = lines.data();
    let nr: usize = line
        .split_whitespace()
        .next()
        .and_then(|t| t.parse::<i64>().ok())
        .filter(|&n| n >= 0)
        .ok_or("failed to read number of predicted lambda values")? as usize;
    let (line, _) = lines.data();
    let lam: Vec<f64> = line
        .split_whitespace()
        .take(nr + 1)
        .map(parse_finite)
        .collect::<Option<_>>()
        .filter(|v: &Vec<f64>| v.len() == nr + 1)
        .ok_or("failed to read slam and or plam values")?;
    let mut f = BarFile {
        slam: lam[0],
        plam: lam[1..].to_vec(),
        e_s: Vec::new(),
        e_p: vec![Vec::new(); nr],
    };
    let (mut line, mut eof) = lines.data();
    while !eof {
        let row: Vec<f64> = line
            .split_whitespace()
            .take(nr + 1)
            .map(parse_finite)
            .collect::<Option<_>>()
            .filter(|v: &Vec<f64>| v.len() == nr + 1)
            .ok_or("failed to read E_s and or E_p values")?;
        f.e_s.push(row[0]);
        for (col, v) in f.e_p.iter_mut().zip(&row[1..]) {
            col.push(*v);
        }
        (line, eof) = lines.data();
    }
    Ok(f)
}

/// The free energies of every neighbouring pair of `files` (sorted by their simulated λ),
/// printing gromos++'s table to `out` as it goes (`bar.cc:295-442`).
pub fn bar(files: &[BarFile], o: &BarOptions, out: &mut dyn Write) -> Result<BarResult, String> {
    if files.len() < 2 {
        return Err(format!("At least two data files required!\n{USAGE}"));
    }
    let kbt = crate::gromospp_kb_warn() * o.temp;
    // gromos++ divides every energy by kT as it reads it
    let reduced: Vec<(Vec<f64>, Vec<Vec<f64>>)> = files
        .iter()
        .map(|f| {
            (
                f.e_s.iter().map(|e| e / kbt).collect(),
                f.e_p
                    .iter()
                    .map(|c| c.iter().map(|e| e / kbt).collect())
                    .collect(),
            )
        })
        .collect();
    let mut lam_s: Vec<f64> = files.iter().map(|f| f.slam).collect();
    lam_s.sort_by(|a, b| a.total_cmp(b));

    let mut head = format!("#{:>7}{:>8}{:>15}{:>15}", "lam_i", "lam_j", "DG", "ee");
    if o.bootstrap > 0 {
        head.push_str(&format!("{:>15}", "bootstrap ee"));
    }
    head.push_str(&format!("{:>15}\n", "overlap int."));
    out.write_all(head.as_bytes()).map_err(io_err)?;

    let seed0 = o.seed.unwrap_or_else(wall_clock_seed);
    let (mut dg_tot, mut err_tot, mut err_tot_c, mut boot_tot) = (0.0, 0.0, 0.0, 0.0);
    let mut intervals = Vec::with_capacity(files.len() - 1);
    for k in 0..files.len() - 1 {
        // literal gromos++ matching: the last file and the last column that fit win
        let (mut i, mut j, mut i_j, mut j_i) = (None, None, None, None);
        for (l, f) in files.iter().enumerate() {
            if f.slam == lam_s[k] {
                i = Some(l);
                for (m, &p) in f.plam.iter().enumerate() {
                    if p == lam_s[k + 1] {
                        i_j = Some(m);
                    }
                }
            }
            if f.slam == lam_s[k + 1] {
                j = Some(l);
                for (m, &p) in f.plam.iter().enumerate() {
                    if p == lam_s[k] {
                        j_i = Some(m);
                    }
                }
            }
        }
        let (Some(i), Some(j)) = (i, j) else {
            return Err("could not find lambda_s value in data".into());
        };
        let Some(i_j) = i_j else {
            return Err(format!(
                "Could not find data for plam = {} in data file for slam = {}\n",
                cpp_g(lam_s[k + 1], 6),
                cpp_g(lam_s[k], 6)
            ));
        };
        let Some(j_i) = j_i else {
            return Err(format!(
                "Could not find data for plam = {} in data file for slam = {}\n",
                cpp_g(lam_s[k], 6),
                cpp_g(lam_s[k + 1], 6)
            ));
        };
        let (e_ii, e_ij) = (&reduced[i].0, &reduced[i].1[i_j]);
        let (e_jj, e_ji) = (&reduced[j].0, &reduced[j].1[j_i]);
        let r = do_bar(
            e_ii,
            e_ij,
            e_jj,
            e_ji,
            o.max_iterations,
            o.convergence,
            true,
        )?;
        let (error_gpp, error_c) = r.error.expect("asked for the error");
        let error = match o.parity {
            Parity::Corrected => error_c,
            Parity::GromosppReleased => error_gpp,
        };
        let (dg, error, error_c) = (r.dg * kbt, error * kbt, error_c * kbt);
        dg_tot += dg;
        err_tot += error * error;
        err_tot_c += error_c * error_c;

        let boot = if o.bootstrap > 0 {
            let mut seed = seed0.wrapping_add(k as u32);
            let (ni, nj) = (e_ii.len(), e_jj.len());
            let mut fb = Vec::with_capacity(o.bootstrap);
            for _ in 0..o.bootstrap {
                let (mut s_i, mut p_i, mut s_j, mut p_j) = (
                    Vec::with_capacity(ni),
                    Vec::with_capacity(ni),
                    Vec::with_capacity(nj),
                    Vec::with_capacity(nj),
                );
                for _ in 0..ni {
                    let x = draw_index(&mut seed, ni);
                    s_i.push(e_ii[x]);
                    p_i.push(e_ij[x]);
                }
                for _ in 0..nj {
                    let x = draw_index(&mut seed, nj);
                    s_j.push(e_jj[x]);
                    p_j.push(e_ji[x]);
                }
                let b = do_bar(
                    &s_i,
                    &p_i,
                    &s_j,
                    &p_j,
                    o.max_iterations,
                    o.convergence,
                    false,
                )?;
                fb.push(b.dg * kbt);
            }
            let var = population_variance(&fb);
            boot_tot += var;
            Some(var)
        } else {
            None
        };

        let print = o.printdist.as_deref().map(|d| (d, k));
        let overlap = compute_overlap(e_ii, e_ij, e_jj, e_ji, print)?;

        let mut row = format!(
            "{:>8}{:>8}{:>15} +/- {:>10}",
            cpp_fixed(lam_s[k], 3),
            cpp_fixed(lam_s[k + 1], 3),
            cpp_fixed(dg, 5),
            cpp_fixed(error, 5)
        );
        if let Some(var) = boot {
            row.push_str(&format!("{:>15}", cpp_fixed(var.sqrt(), 5)));
        }
        row.push_str(&format!("{:>15}\n", cpp_fixed(overlap, 5)));
        out.write_all(row.as_bytes()).map_err(io_err)?;
        intervals.push(Interval {
            lam_i: lam_s[k],
            lam_j: lam_s[k + 1],
            dg,
            error,
            error_corrected: error_c,
            bootstrap: boot.map(f64::sqrt),
            overlap,
            iterations: r.iterations,
        });
    }
    let mut total = format!(
        "\n# total{:>24} +/- {:>10}",
        cpp_fixed(dg_tot, 5),
        cpp_fixed(err_tot.sqrt(), 5)
    );
    if o.bootstrap > 0 {
        total.push_str(&format!("{:>15}", cpp_fixed(boot_tot.sqrt(), 5)));
    }
    total.push('\n');
    out.write_all(total.as_bytes()).map_err(io_err)?;
    Ok(BarResult {
        intervals,
        dg_total: dg_tot,
        error_total: err_tot.sqrt(),
        error_total_corrected: err_tot_c.sqrt(),
        bootstrap_total: (o.bootstrap > 0).then(|| boot_tot.sqrt()),
    })
}

/// What [`do_bar`] returns, in kT.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BarStep {
    /// G(j) − G(i) / kT.
    pub dg: f64,
    /// gromos++ 9cd9f6b3's analytical error and the corrected one, when asked for.
    pub error: Option<(f64, f64)>,
    pub iterations: usize,
}

/// gromos++ `do_bar` (`bar.cc:453-644`): reduced energies (÷kT) of the frames of state i at i
/// and j and of state j at j and i. ΔG starts at 0 and is iterated until its relative change is
/// at most `convergence` (always at least two iterations) or `max_iterations` is reached.
pub fn do_bar(
    e_ii: &[f64],
    e_ij: &[f64],
    e_jj: &[f64],
    e_ji: &[f64],
    max_iterations: usize,
    convergence: f64,
    with_error: bool,
) -> Result<BarStep, String> {
    let (n_i, n_j) = (e_ii.len(), e_jj.len());
    if n_i == 0 || n_j == 0 {
        return Err("bar: a state without frames".into());
    }
    let m = (n_i as f64 / n_j as f64).ln();
    let mut reldiff = 1e15;
    let mut dg = 0.0;
    let mut max_ij = -1e15;
    let mut max_ji = -1e15;
    for n in 0..n_i {
        if e_ij[n] - e_ii[n] > max_ij {
            max_ij = e_ij[n] - e_ii[n];
        }
    }
    for n in 0..n_j {
        if e_jj[n] - e_ji[n] > max_ji {
            max_ji = e_jj[n] - e_ji[n];
        }
    }
    let exp_ij: Vec<f64> = (0..n_i)
        .map(|n| (e_ij[n] - e_ii[n] - max_ij).exp())
        .collect();
    let exp_ji: Vec<f64> = (0..n_j)
        .map(|n| (e_jj[n] - e_ji[n] - max_ji).exp())
        .collect();
    let berg = |s: f64, v: f64| cmax(s, v) + (1.0 + (cmin(s, v) - cmax(s, v)).exp()).ln();

    let mut k = 0usize;
    while k < max_iterations && reldiff > convergence {
        let dg_prev = dg;
        let mx_ij = m + max_ij - dg;
        let exp_mx_ij = (-mx_ij).exp();
        let mx_ji = m + max_ji - dg;
        let exp_mx_ji = (-mx_ji).exp();
        let mut lnsum_ij = -mx_ij - (exp_mx_ij + exp_ij[0]).ln();
        for n in 1..n_i {
            let val = -mx_ij - (exp_mx_ij + exp_ij[n]).ln();
            lnsum_ij = berg(lnsum_ij, val);
        }
        let mut lnsum_ji = -max_ji - (exp_mx_ji + exp_ji[0]).ln() + e_jj[0] - e_ji[0];
        for n in 1..n_j {
            let val = -max_ji - (exp_mx_ji + exp_ji[n]).ln() + e_jj[n] - e_ji[n];
            lnsum_ji = berg(lnsum_ji, val);
        }
        dg += lnsum_ji - lnsum_ij;
        if k > 0 {
            reldiff = ((dg - dg_prev) / dg_prev).abs();
        }
        k += 1;
    }
    if k == max_iterations {
        eprintln!("WARNING: Iteration finished because maximum number of iterations reached");
    }
    if !with_error {
        return Ok(BarStep {
            dg,
            error: None,
            iterations: k,
        });
    }
    if n_i < 2 || n_j < 2 {
        return Err("bar: the error estimate needs at least two frames per state".into());
    }
    let mx_ij = m + max_ij - dg;
    let exp_mx_ij = (-mx_ij).exp();
    let mx_ji = m + max_ji - dg;
    let exp_mx_ji = (-mx_ji).exp();
    let delta_ij: Vec<f64> = (0..n_i)
        .map(|n| -mx_ij - (exp_mx_ij + exp_ij[n]).ln())
        .collect();
    let delta_ji: Vec<f64> = (0..n_j)
        .map(|n| -max_ji - (exp_mx_ji + exp_ji[n]).ln() + e_jj[n] - e_ji[n])
        .collect();
    let side = |d: &[f64], n: usize| -> Result<(f64, f64), String> {
        let ave = lnexp_ave(d);
        let (var, sign) = lnexp_covariance(d, d)?;
        if sign < 0 {
            return Err("Got negative variance!".into());
        }
        let si = lnexp_stat_ineff(d, d)?;
        let d2 = var - (n as f64).ln() + si;
        Ok((ave, (d2 - 2.0 * ave).exp()))
    };
    let (d_ij, error_i2) = side(&delta_ij, n_i)?;
    let (d_ji, error_j2) = side(&delta_ji, n_j)?;
    // gromos++ 9cd9f6b3 divides by d² a second time (bar.cc:639)
    let error = (error_i2 / (d_ij * d_ij) + error_j2 / (d_ji * d_ji)).sqrt();
    let corrected = (error_i2 + error_j2).sqrt();
    Ok(BarStep {
        dg,
        error: Some((error, corrected)),
        iterations: k,
    })
}

/// gromos++ `compute_overlap` (`bar.cc:646-696`): the ΔE distributions of both states on 100
/// bins over their joint [min, max) and `2 Σ P_i P_j / (P_i + P_j)`; with `print = (dir, k)` the
/// two normalised distributions go to `dir/dist_<k>_<k+1>.dat`.
pub fn compute_overlap(
    e_ii: &[f64],
    e_ij: &[f64],
    e_jj: &[f64],
    e_ji: &[f64],
    print: Option<(&Path, usize)>,
) -> Result<f64, String> {
    let d_ij: Vec<f64> = e_ij.iter().zip(e_ii).map(|(a, b)| a - b).collect();
    let d_ji: Vec<f64> = e_jj.iter().zip(e_ji).map(|(a, b)| a - b).collect();
    let smin = |v: &[f64]| v[1..].iter().fold(v[0], |m, &x| if x < m { x } else { m });
    let smax = |v: &[f64]| v[1..].iter().fold(v[0], |m, &x| if x > m { x } else { m });
    let (min_i, min_j) = (smin(&d_ij), smin(&d_ji));
    let (max_i, max_j) = (smax(&d_ij), smax(&d_ji));
    let min = if min_i < min_j { min_i } else { min_j };
    let max = if max_i > max_j { max_i } else { max_j };
    let nbins = 100;
    let di = Distribution::over_values_in(&d_ij, min, max, nbins, false)?;
    let dj = Distribution::over_values_in(&d_ji, min, max, nbins, false)?;
    let mut oi = 0.0;
    for n in 0..nbins {
        let pi = di.count(n) as f64 / di.n_values() as f64;
        let pj = dj.count(n) as f64 / dj.n_values() as f64;
        if pi + pj != 0.0 {
            oi += (pi * pj) / (pi + pj);
        }
    }
    if let Some((dir, k)) = print {
        let g = |x: f64| cpp_g(x, 6);
        let mut s = format!(
            "# Distribution from state {k} to state {}\n#\n# DeltaE    P(DeltaE)\n",
            k + 1
        );
        s.push_str(&di.write_normalized(&g));
        s.push_str(&format!(
            "# Distribution from state {} to state {k}\n#\n# DeltaE    P(-DeltaE)\n",
            k + 1
        ));
        s.push_str(&dj.write_normalized(&g));
        let path = dir.join(format!("dist_{k}_{}.dat", k + 1));
        std::fs::write(&path, s).map_err(|e| format!("{}: {e}", path.display()))?;
    }
    Ok(2.0 * oi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reader_skips_comments_and_drops_an_unterminated_last_line() {
        let text = "#      nr lam_p\n              1\n# x\n           0.00           0.50\n\
                    #           E_s\n    1.0 2.0 \n    3.0 4.0 \n    5.0 6.0";
        let f = parse_bar_file(text).unwrap();
        assert_eq!(f.slam, 0.0);
        assert_eq!(f.plam, vec![0.5]);
        assert_eq!(f.e_s, vec![1.0, 3.0]);
        assert_eq!(f.e_p, vec![vec![2.0, 4.0]]);
        // blank and short lines, and non-finite numbers, are errors
        for bad in [
            "1\n0 0.5\n1.0 2.0\n\n",
            "1\n0 0.5\n1.0\n",
            "1\n0 0.5\n1.0 nan\n",
        ] {
            assert!(parse_bar_file(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn bar_of_identical_states_is_zero_and_symmetric() {
        let e: Vec<f64> = (0..50).map(|i| (i as f64 * 0.37).sin()).collect();
        let r = do_bar(&e, &e, &e, &e, 500, 1e-10, false).unwrap();
        assert!(r.dg.abs() < 1e-12, "{}", r.dg);
        assert!(r.iterations >= 2);
    }
}
