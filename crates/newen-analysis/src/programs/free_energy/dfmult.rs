//! dfmult — free-energy differences from a reference state to several end states, and between
//! the end states, with the log-exponential estimators of gromos++ `dfmult`.
//!
//! Usage: dfmult @temp <T> @stateR <file> @endstates <files…>
//!
//! Every file is a time series `time  H`; ΔF_iR = −kT ln <exp(−β(H_i − H_R))>_R with the error
//! from the variance and the statistical inefficiency (gromos++ `Stat::lnexp*`).

use crate::args::Arguments;
use crate::lnexp::{lnexp_ave, lnexp_covariance, lnexp_stat_ineff};
use crate::programs::io_err;
use newen_io::args::{Program, Status};
use newen_io::table::read_columns;
use std::io::Write;

pub const USAGE: &str = "# dfmult
\t@temp      <temperature for perturbation>
\t@stateR    <energy file for state R>
\t@endstates <energy files of endstates>";

pub const KNOWN: &[&str] = &["temp", "stateR", "endstates"];

pub const PROGRAM: Program = Program {
    name: "dfmult",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// C++ `scientific` with 7 decimals in an 18-wide field.
fn sci(x: f64) -> String {
    let s = format!("{:.7e}", x);
    let (m, e) = s.split_once('e').unwrap();
    let exp: i32 = e.parse().unwrap();
    format!(
        "{:>18}",
        format!("{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    )
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let temp: f64 = args.require("temp")?;
    let kt = crate::GROMOSPP_BOLTZMANN * temp;
    let r = read_columns(
        args.value("stateR")
            .map_err(|_| "@stateR is required".to_string())?,
    )
    .map_err(|e| format!("state R: {e}"))?;
    if args.count("endstates") < 1 {
        return Err("@endstates is required".into());
    }
    let mut all: Vec<Vec<f64>> = Vec::new();
    for f in args.values("endstates") {
        let x = read_columns(f).map_err(|e| format!("end state {f}: {e}"))?;
        if x.len() != r.len() {
            return Err("Could not read from energy file. Check lengths!".into());
        }
        let mut series = Vec::with_capacity(x.len());
        for (rr, rx) in r.iter().zip(&x) {
            if rr.len() < 2 || rx.len() < 2 {
                return Err("Could not read from energy file. Check lengths!".into());
            }
            series.push(-(rx[1] - rr[1]) / kt);
        }
        all.push(series);
    }
    let mut text = String::new();
    text.push_str(&format!("{:>36}{:>18}\n", "#DF (kJ/mol)", "err"));
    for i in 0..all.len() {
        let n = all[i].len() as f64;
        let df_ir = lnexp_ave(&all[i]);
        let (var_ii, sign) = lnexp_covariance(&all[i], &all[i])?;
        if sign < 0 {
            return Err("Got negative variance!".into());
        }
        let si_ii = lnexp_stat_ineff(&all[i], &all[i])?;
        let d2i = var_ii - n.ln() + si_ii;
        text.push_str(&format!(
            "{:>18}{}{}\n",
            format!("DF_{}_R", i + 1),
            sci(-kt * df_ir),
            sci(kt * (d2i - 2.0 * df_ir).exp().sqrt())
        ));
        for j in i + 1..all.len() {
            let df_jr = lnexp_ave(&all[j]);
            let (var_jj, sign_jj) = lnexp_covariance(&all[j], &all[j])?;
            if sign_jj < 0 {
                return Err("Got negative variance!".into());
            }
            let si_jj = lnexp_stat_ineff(&all[j], &all[j])?;
            let d2j = var_jj - n.ln() + si_jj;
            let df_ji = df_jr - df_ir;
            let (var_ji, sign_ji) = lnexp_covariance(&all[j], &all[i])?;
            let mut si_ji = lnexp_stat_ineff(&all[j], &all[i])?;
            let max_si_ji = ((si_ii + si_jj).exp().sqrt()
                * (var_ii + var_jj - 2.0 * var_ji).exp().sqrt().abs())
            .ln();
            if si_ji > max_si_ji && max_si_ji > 0.0 {
                eprintln!("setting si_ji to max_si_ji");
                si_ji = max_si_ji;
            }
            let djdi = var_ji - n.ln() + si_ji;
            let err = kt
                * ((d2i - 2.0 * df_ir).exp() + (d2j - 2.0 * df_jr).exp()
                    - 2.0 * (djdi - (df_jr + df_ir)).exp() * sign_ji as f64)
                    .sqrt();
            text.push_str(&format!(
                "{:>18}{}{}\n",
                format!("DF_{}_{}", j + 1, i + 1),
                sci(-kt * df_ji),
                sci(err)
            ));
        }
    }
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
