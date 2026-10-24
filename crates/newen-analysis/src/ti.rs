//! Plain multi-window thermodynamic integration: ⟨∂H/∂λ⟩ per window from the total dH/dλ of
//! each `.trg` frame, integrated over λ with the trapezoid rule.
//!
//! This is what newen's `ext_ti_ana` did before it became the gromos++ program of that name
//! (which reweights one simulation onto other λ values instead, item 1.8). It has no gromos++
//! counterpart; `scripts/ti_ch4.py --estimator ti` and `newen.analysis.ti_trapezoid` use it.

use newen_core::stat::Stat;
use newen_io::read_free_energy_trajectory;

/// One λ window: its ⟨∂H/∂λ⟩ (kJ/mol), the error estimate and the number of frames used.
#[derive(Debug, Clone, PartialEq)]
pub struct TiWindow {
    pub lambda: f64,
    pub mean: f64,
    /// `Stat::ee` (block averaging) with 4 or more frames, the rmsd with fewer.
    pub ee: f64,
    pub n: usize,
}

impl TiWindow {
    /// ⟨∂H/∂λ⟩ of the dH/dλ samples of one window; `None` without samples.
    pub fn from_samples(lambda: f64, samples: impl IntoIterator<Item = f64>) -> Option<Self> {
        let mut stat = Stat::new();
        for x in samples {
            stat.add(x);
        }
        let n = stat.n();
        if n == 0 {
            return None;
        }
        let ee = if n >= 4 { stat.ee() } else { stat.rmsd() };
        Some(TiWindow {
            lambda,
            mean: stat.ave(),
            ee,
            n,
        })
    }
}

/// The windows in ascending λ and, with two or more, ΔG = ∫⟨∂H/∂λ⟩dλ by the trapezoid rule and
/// its error from the per-window errors with the trapezoid weights (windows independent).
#[derive(Debug, Clone, PartialEq)]
pub struct TiResult {
    pub windows: Vec<TiWindow>,
    pub dg: Option<f64>,
    pub dg_error: Option<f64>,
}

/// ⟨∂H/∂λ⟩ of each `.trg` file (`FREEENER[1]`, the total dH/dλ of a frame) after dropping its
/// first `skip` frames. λ is `lambdas[i]` when given, else `RLAM` of the file's first frame. A
/// file without free-energy frames is skipped with a warning on stderr.
pub fn trapezoid_from_trg(
    files: &[String],
    lambdas: Option<&[f64]>,
    skip: usize,
) -> Result<TiResult, String> {
    if let Some(l) = lambdas {
        if l.len() != files.len() {
            return Err(format!(
                "{} lambda values for {} files: give one per file",
                l.len(),
                files.len()
            ));
        }
    }
    let mut windows = Vec::new();
    for (i, path) in files.iter().enumerate() {
        let traj = read_free_energy_trajectory(path).map_err(|e| format!("{path}: {e}"))?;
        for w in traj.warnings {
            eprintln!("{w}");
        }
        if traj.frames.is_empty() {
            eprintln!("Warning: {path} contains no free-energy frames — skipping");
            continue;
        }
        let lambda = lambdas.map_or(traj.frames[0].lambda, |l| l[i]);
        match TiWindow::from_samples(lambda, traj.frames.iter().skip(skip).map(|f| f.dhdl_total)) {
            Some(w) => windows.push(w),
            None => eprintln!("Warning: window λ={lambda} has 0 samples after skipping {skip}"),
        }
    }
    Ok(trapezoid(windows))
}

/// Sort the windows by λ and integrate: ΔG = Σ ½(y_k + y_{k+1})(λ_{k+1} − λ_k),
/// σ²(ΔG) = Σ w_k² σ_k² with w_k = ½(λ_{k+1} − λ_{k−1}) (one-sided at the ends).
pub fn trapezoid(mut windows: Vec<TiWindow>) -> TiResult {
    windows.sort_by(|a, b| a.lambda.total_cmp(&b.lambda));
    let k = windows.len();
    if k < 2 {
        return TiResult {
            windows,
            dg: None,
            dg_error: None,
        };
    }
    let mut dg = 0.0;
    for p in windows.windows(2) {
        dg += 0.5 * (p[0].mean + p[1].mean) * (p[1].lambda - p[0].lambda);
    }
    let mut var = 0.0;
    for i in 0..k {
        let lo = windows[i.saturating_sub(1)].lambda;
        let hi = windows[(i + 1).min(k - 1)].lambda;
        let w = 0.5 * (hi - lo);
        var += w * w * windows[i].ee * windows[i].ee;
    }
    TiResult {
        windows,
        dg: Some(dg),
        dg_error: Some(var.sqrt()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn w(lambda: f64, mean: f64, ee: f64) -> TiWindow {
        TiWindow {
            lambda,
            mean,
            ee,
            n: 10,
        }
    }

    #[test]
    fn trapezoid_integrates_a_line_exactly_and_sorts() {
        // y = 2λ + 1 on [0, 1]: ∫ = 2
        let r = trapezoid(vec![w(1.0, 3.0, 0.1), w(0.0, 1.0, 0.1), w(0.5, 2.0, 0.2)]);
        assert_eq!(r.windows[0].lambda, 0.0);
        assert!((r.dg.unwrap() - 2.0).abs() < 1e-12);
        // weights 0.25, 0.5, 0.25
        let e = (0.25f64 * 0.25 * 0.01 + 0.25 * 0.04 + 0.25 * 0.25 * 0.01).sqrt();
        assert!((r.dg_error.unwrap() - e).abs() < 1e-12);
        assert_eq!(trapezoid(vec![w(0.0, 1.0, 0.1)]).dg, None);
    }
}
