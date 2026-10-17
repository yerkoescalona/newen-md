//! EDS parameter updates shared by `eds_update_1` and `eds_update_2` (gromos++ programs of the
//! same names): the energy-offset update from reference-state and end-state energy time series,
//! and the per-pair statistics the s-updates are built from. Log-sum-exp is written exactly as
//! gromos++ writes it so that the iterated updates agree digit for digit.

use crate::lnexp::lnexp_ave;

/// One pair of EDS states and the parameters gromos++ keeps per pair.
#[derive(Debug, Clone, Default)]
pub struct Pair {
    pub i: usize,
    pub j: usize,
    pub s_old: f64,
    pub s_new: f64,
    pub eir: f64,
    pub ejr: f64,
    pub term: f64,
    pub ln_i: f64,
    pub ln_j: f64,
    /// V_j − V_i time series
    pub x: Vec<f64>,
}

/// ln(eᵃ + eᵇ) as gromos++ writes it: max + ln(1 + e^(min − max)).
pub fn ln_sum_exp(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a < b { (b, a) } else { (a, b) };
    hi + (1.0 + (lo - hi).exp()).ln()
}

/// gromos++ `read_data`: the second column of a `time value` file; `#` starts a comment and blank
/// lines are skipped. (gromos++ also drops the last line when the file has no final newline; we
/// read it.)
pub fn read_series(name: &str, path: &str) -> Result<Vec<f64>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|_| format!("Could not open time series file for {name}.\n"))?;
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("");
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        let t = it.next().and_then(|s| s.parse::<f64>().ok());
        let q = it.next().and_then(|s| s.parse::<f64>().ok());
        match (t, q) {
            (Some(_), Some(q)) => out.push(q),
            _ => {
                return Err(format!(
                    "Error when reading from {name} time series file.\n"
                ))
            },
        }
    }
    Ok(out)
}

/// gromos++ `calculate_EiR`: the new energy offset of state `p_i` from the reference-state
/// energies `vr` and the end-state energies `vy`. `pair_s` selects the s of a pair that enters the
/// reference-state energy (`s_old` in eds_update_1, `s_new` in eds_update_2); `factor` is the
/// functional-form prefactor.
#[allow(clippy::too_many_arguments)]
pub fn calculate_eir(
    pairs: &[Pair],
    pair_s: &dyn Fn(&Pair) -> f64,
    vr: &[f64],
    vy: &[Vec<f64>],
    eir_old: &[f64],
    beta: f64,
    p_i: usize,
    factor: f64,
) -> f64 {
    let numstat = vy.len();
    let mut statesum = Vec::with_capacity(vr.len());
    let mut weight = Vec::with_capacity(vr.len());
    for k in 0..vr.len() {
        let mut sum_ln = 1.0;
        for j in 0..numstat {
            if j != p_i {
                let deir = eir_old[j] - eir_old[p_i];
                let diff_ln = -beta * ((vy[j][k] - vy[p_i][k]) - deir);
                sum_ln = ln_sum_exp(sum_ln, diff_ln);
            }
        }
        let s0 = pair_s(&pairs[0]);
        let beta_s = -beta * s0;
        let part_a = beta_s * (vy[pairs[0].i][k] - pairs[0].eir);
        let part_b = beta_s * (vy[pairs[0].j][k] - pairs[0].ejr);
        let mut vrnew = ln_sum_exp(part_a, part_b) / s0;
        for p in &pairs[1..] {
            let s = pair_s(p);
            let beta_s = -beta * s;
            let part_a = beta_s * (vy[p.i][k] - p.eir);
            let part_b = beta_s * (vy[p.j][k] - p.ejr);
            let elem = ln_sum_exp(part_a, part_b) / s;
            vrnew = ln_sum_exp(vrnew, elem);
        }
        vrnew = -(1.0 / beta) * (vrnew + factor.ln());
        let diff = -beta * (vrnew - vr[k]);
        statesum.push(diff - sum_ln);
        weight.push(diff);
    }
    let ln_x_y = lnexp_ave(&statesum) - lnexp_ave(&weight);
    -(ln_x_y / beta) + eir_old[p_i]
}

/// gromos++ eds_update_1 pair statistics: `x` = V_j − V_i and the log-average terms `ln_i`,
/// `ln_j` of the pair.
pub fn pair_statistics(pair: &mut Pair, vr: &[f64], vy: &[Vec<f64>], beta: f64) {
    pair.x = (0..vy[pair.j].len())
        .map(|k| vy[pair.j][k] - vy[pair.i][k])
        .collect();
    let (mut weight_a, mut weight_b, mut sum_a, mut sum_b) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for k in 0..vr.len() {
        let tempref_a = -beta * (vy[pair.i][k] - vr[k]);
        let tempref_b = -beta * (vy[pair.j][k] - vr[k]);
        sum_a.push(-beta * pair.x[k].abs() + tempref_a);
        sum_b.push(-beta * pair.x[k].abs() + tempref_b);
        weight_a.push(tempref_a);
        weight_b.push(tempref_b);
    }
    pair.ln_i = lnexp_ave(&sum_a) - lnexp_ave(&weight_a);
    pair.ln_j = lnexp_ave(&sum_b) - lnexp_ave(&weight_b);
}

/// gromos++ eds_update_1 `calculate_s`: s = min(1/|−ln_i − term|, 1/|−ln_j + term|), capped at 1.
pub fn update_s(pair: &mut Pair) {
    let s_a = 1.0 / (-pair.ln_i - pair.term).abs();
    let s_b = 1.0 / (-pair.ln_j + pair.term).abs();
    pair.s_new = s_a.min(s_b).min(1.0);
}
