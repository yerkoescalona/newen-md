//! mbar — free energies of all simulated λ at once by the multistate Bennett acceptance ratio
//! (Shirts & Chodera, J. Chem. Phys. 129 (2008) 124105). No gromos++ counterpart.
//!
//! Input: the `ext_ti_ana @bar_data` files of K simulations, each listing its frames' energies at
//! the other simulated λ (`@bar_data` with every simulated λ). The reduced potentials
//! u_k(x_n) = E_k(x_n)/kT of every frame at every state give the free energies f_k from the
//! self-consistent equations f_i = −ln Σ_n exp(−u_i(x_n)) / Σ_k N_k exp(f_k − u_k(x_n)),
//! iterated, then finished by Newton–Raphson on the convex MBAR objective. The uncertainty is
//! the asymptotic covariance Θ = V Σ (I − Σ Vᵀ N V Σ)⁺ Σ Vᵀ of the weight matrix W = U Σ Vᵀ
//! (eq. D8 of Shirts & Chodera; pymbar's `svd-ew` form, pseudo-inverse cut at 10⁻¹⁰ of the
//! largest singular value), and δ(f_j − f_i)² = Θ_ii + Θ_jj − 2Θ_ij. A frame's energy at its own
//! λ is `E_s`, as `bar` uses it, so two states give `bar`'s ΔG.

use super::bar::{read_bar_file, BarFile};
use super::common::get_value;
use crate::args::Arguments;
use crate::programs::io_err;
use gromos_io::args::{Program, Status};
use gromos_io::table::{cpp_fixed, cpp_g};
use nalgebra::{DMatrix, DVector};
use std::io::Write;

pub const USAGE: &str = "# mbar
\t@files            <data files (ext_ti_ana bar_data output, with every simulated lambda)>
\t@temp             <temperature>
\t[@maxiterations   <maximum number of self-consistent iterations (default: 10000)>]
\t[@convergence     <largest change of a reduced free energy that ends them (default: 1E-10)>]";

pub const KNOWN: &[&str] = &["files", "temp", "maxiterations", "convergence"];

pub const PROGRAM: Program = Program {
    name: "mbar",
    known: KNOWN,
    usage: USAGE,
    status: Status::NoCounterpart,
    run,
};

/// The MBAR solution in kT units.
#[derive(Debug, Clone, PartialEq)]
pub struct MbarResult {
    /// f_k − f_0.
    pub f: Vec<f64>,
    /// `delta_f[i][j]` = f_j − f_i.
    pub delta_f: Vec<Vec<f64>>,
    /// Its asymptotic standard error.
    pub d_delta_f: Vec<Vec<f64>>,
    pub self_consistent_iterations: usize,
    pub newton_iterations: usize,
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    evaluate(argv, out).map(|_| ())
}

/// [`run`], returning the simulated λ and the solution (kT units; the printed numbers are
/// these times kT).
pub fn evaluate(argv: &[String], out: &mut dyn Write) -> Result<(Vec<f64>, MbarResult), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let temp: f64 = get_value(&args, "temp", None)?;
    let max_iterations: i64 = get_value(&args, "maxiterations", Some(10000))?;
    let convergence: f64 = get_value(&args, "convergence", Some(1e-10))?;
    if args.count("files") <= 1 {
        return Err(format!("At least two data files required!\n{USAGE}"));
    }
    if max_iterations < 0 {
        return Err("@maxiterations cannot be negative".into());
    }
    let files = args
        .values("files")
        .iter()
        .map(|p| read_bar_file(p))
        .collect::<Result<Vec<_>, _>>()?;
    let kbt = crate::gromospp_kb_warn() * temp;
    let (lam, u_kn, n_k) = reduced_potentials(&files, kbt)?;
    let r = mbar(&u_kn, &n_k, max_iterations as usize, convergence)?;
    eprintln!(
        "# MBAR: {} self-consistent and {} Newton-Raphson iterations",
        r.self_consistent_iterations, r.newton_iterations
    );
    let mut s = format!("#{:>7}{:>8}{:>15}{:>15}\n", "lam_i", "lam_j", "DG", "ee");
    for k in 0..lam.len() - 1 {
        s.push_str(&format!(
            "{:>8}{:>8}{:>15} +/- {:>10}\n",
            cpp_fixed(lam[k], 3),
            cpp_fixed(lam[k + 1], 3),
            cpp_fixed(r.delta_f[k][k + 1] * kbt, 5),
            cpp_fixed(r.d_delta_f[k][k + 1] * kbt, 5)
        ));
    }
    let last = lam.len() - 1;
    s.push_str(&format!(
        "\n# total{:>24} +/- {:>10}\n",
        cpp_fixed(r.delta_f[0][last] * kbt, 5),
        cpp_fixed(r.d_delta_f[0][last] * kbt, 5)
    ));
    out.write_all(s.as_bytes()).map_err(io_err)?;
    Ok((lam, r))
}

/// The states (the files' simulated λ, ascending) and the reduced potentials `u_kn[k][n]` of
/// every frame `n` (the frames of the lowest λ first) at every state `k`, with `N_k`. A file
/// must list every other simulated λ among its `plam` (exact match, as `bar` pairs them).
pub fn reduced_potentials(
    files: &[BarFile],
    kbt: f64,
) -> Result<(Vec<f64>, Vec<Vec<f64>>, Vec<usize>), String> {
    let mut order: Vec<usize> = (0..files.len()).collect();
    order.sort_by(|&a, &b| files[a].slam.total_cmp(&files[b].slam));
    let lam: Vec<f64> = order.iter().map(|&l| files[l].slam).collect();
    if lam.windows(2).any(|w| w[0] == w[1]) {
        return Err("two files with the same simulated lambda".into());
    }
    let k_states = lam.len();
    let mut u_kn = vec![Vec::new(); k_states];
    let mut n_k = Vec::with_capacity(k_states);
    for &l in &order {
        let f = &files[l];
        n_k.push(f.e_s.len());
        for (k, &lk) in lam.iter().enumerate() {
            if lk == f.slam {
                u_kn[k].extend(f.e_s.iter().map(|e| e / kbt));
                continue;
            }
            let col = f.plam.iter().rposition(|&p| p == lk).ok_or_else(|| {
                format!(
                    "the data file for slam = {} has no energies at plam = {}: \
                     give ext_ti_ana @bar_data every simulated lambda",
                    cpp_g(f.slam, 6),
                    cpp_g(lk, 6)
                )
            })?;
            u_kn[k].extend(f.e_p[col].iter().map(|e| e / kbt));
        }
    }
    Ok((lam, u_kn, n_k))
}

/// ln Σ exp(v), stable.
fn logsumexp(v: impl Iterator<Item = f64> + Clone) -> f64 {
    let m = v.clone().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return m;
    }
    m + v.map(|x| (x - m).exp()).sum::<f64>().ln()
}

/// MBAR on reduced potentials `u_kn` (K rows of N = Σ N_k values; the samples of state 0 first,
/// then state 1, …) with `n_k` samples per state. `tolerance` bounds the largest change of an
/// f_k in the self-consistent iteration before the Newton–Raphson finish.
pub fn mbar(
    u_kn: &[Vec<f64>],
    n_k: &[usize],
    max_iterations: usize,
    tolerance: f64,
) -> Result<MbarResult, String> {
    let k_states = n_k.len();
    let n: usize = n_k.iter().sum();
    if k_states < 2 || u_kn.len() != k_states {
        return Err("mbar: need at least two states and one row of u_kn per state".into());
    }
    if u_kn.iter().any(|r| r.len() != n) {
        return Err(format!(
            "mbar: every row of u_kn needs N = sum(N_k) = {n} values"
        ));
    }
    if n_k.contains(&0) {
        return Err("mbar: a state without samples".into());
    }
    if u_kn.iter().flatten().any(|v| !v.is_finite()) {
        return Err("mbar: u_kn has a non-finite value".into());
    }
    let ln_n: Vec<f64> = n_k.iter().map(|&c| (c as f64).ln()).collect();
    let log_denominators = |f: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|s| logsumexp((0..k_states).map(|k| ln_n[k] + f[k] - u_kn[k][s])))
            .collect()
    };
    let mut f = vec![0.0; k_states];
    let mut sc = 0;
    while sc < max_iterations {
        let ld = log_denominators(&f);
        let mut f_new: Vec<f64> = (0..k_states)
            .map(|i| -logsumexp((0..n).map(|s| -u_kn[i][s] - ld[s])))
            .collect();
        let f0 = f_new[0];
        f_new.iter_mut().for_each(|x| *x -= f0);
        let change = f_new
            .iter()
            .zip(&f)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        f = f_new;
        sc += 1;
        if change < tolerance {
            break;
        }
    }

    // W[n][k] = exp(f_k − u_kn − ln D_n): Σ_n W_nk = 1 at the solution
    let weights = |f: &[f64]| -> DMatrix<f64> {
        let ld = log_denominators(f);
        DMatrix::from_fn(n, k_states, |s, k| (f[k] - u_kn[k][s] - ld[s]).exp())
    };
    let gradient = |w: &DMatrix<f64>| -> Vec<f64> {
        (0..k_states)
            .map(|k| n_k[k] as f64 * (w.column(k).sum() - 1.0))
            .collect()
    };
    let norm = |g: &[f64]| g.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mut w = weights(&f);
    let mut g = gradient(&w);
    let mut newton = 0;
    for _ in 0..50 {
        if norm(&g) < 1e-12 * n as f64 {
            break;
        }
        // Hessian of the MBAR objective, f_0 held at 0
        let m = k_states - 1;
        let h = DMatrix::from_fn(m, m, |a, b| {
            let (i, j) = (a + 1, b + 1);
            let (ni, nj) = (n_k[i] as f64, n_k[j] as f64);
            let cross: f64 = w.column(i).dot(&w.column(j));
            let diag = if i == j { ni * w.column(i).sum() } else { 0.0 };
            diag - ni * nj * cross
        });
        let rhs = DVector::from_fn(m, |a, _| -g[a + 1]);
        let Some(step) = h.lu().solve(&rhs) else {
            break;
        };
        let trial: Vec<f64> = (0..k_states)
            .map(|k| if k == 0 { f[0] } else { f[k] + step[k - 1] })
            .collect();
        let w_trial = weights(&trial);
        let g_trial = gradient(&w_trial);
        if norm(&g_trial) >= norm(&g) {
            break;
        }
        f = trial;
        w = w_trial;
        g = g_trial;
        newton += 1;
    }

    let theta = asymptotic_covariance(&w, n_k)?;
    let mut delta_f = vec![vec![0.0; k_states]; k_states];
    let mut d_delta_f = vec![vec![0.0; k_states]; k_states];
    for i in 0..k_states {
        for j in 0..k_states {
            delta_f[i][j] = f[j] - f[i];
            let mut d2 = theta[(i, i)] + theta[(j, j)] - 2.0 * theta[(i, j)];
            if d2 < 0.0 && d2 > -1e-10 {
                d2 = 0.0;
            }
            d_delta_f[i][j] = d2.sqrt();
        }
    }
    Ok(MbarResult {
        f,
        delta_f,
        d_delta_f,
        self_consistent_iterations: sc,
        newton_iterations: newton,
    })
}

/// Θ = V Σ (I − Σ Vᵀ N V Σ)⁺ Σ Vᵀ for W = U Σ Vᵀ (Shirts & Chodera 2008, appendix D).
fn asymptotic_covariance(w: &DMatrix<f64>, n_k: &[usize]) -> Result<DMatrix<f64>, String> {
    let k_states = n_k.len();
    let svd = w.clone().svd(false, true);
    let v_t = svd.v_t.ok_or("mbar: SVD of the weight matrix failed")?;
    let v = v_t.transpose();
    let sigma = DMatrix::from_diagonal(&svd.singular_values);
    let nd = DMatrix::from_fn(
        k_states,
        k_states,
        |i, j| {
            if i == j {
                n_k[i] as f64
            } else {
                0.0
            }
        },
    );
    let inner =
        DMatrix::<f64>::identity(k_states, k_states) - &sigma * v.transpose() * &nd * &v * &sigma;
    let pinv = pseudo_inverse(&inner, 1e-10)?;
    Ok(&v * &sigma * pinv * &sigma * v.transpose())
}

/// Moore–Penrose pseudo-inverse, singular values below `rcond` × the largest dropped
/// (`numpy.linalg.pinv`).
fn pseudo_inverse(a: &DMatrix<f64>, rcond: f64) -> Result<DMatrix<f64>, String> {
    let svd = a.clone().svd(true, true);
    let (u, v_t) = match (svd.u, svd.v_t) {
        (Some(u), Some(v_t)) => (u, v_t),
        _ => return Err("mbar: SVD for the pseudo-inverse failed".into()),
    };
    let s = &svd.singular_values;
    let cut = rcond * s.iter().fold(0.0f64, |m, &x| m.max(x));
    let inv = DMatrix::from_fn(s.len(), s.len(), |i, j| {
        if i == j && s[i] > cut {
            1.0 / s[i]
        } else {
            0.0
        }
    });
    Ok(v_t.transpose() * inv * u.transpose())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two harmonic states u_k(x) = (x − μ_k)²/2 with exact samples of the first one's
    /// distribution quantiles: the free energies are equal (same width), so Δf → 0, and
    /// shifting a state's potential by c shifts f by c.
    #[test]
    fn shift_of_a_potential_is_its_free_energy() {
        let xs: Vec<f64> = (0..200).map(|i| -2.0 + 4.0 * i as f64 / 199.0).collect();
        let u0: Vec<f64> = xs.iter().map(|x| 0.5 * x * x).collect();
        let u1: Vec<f64> = u0.iter().map(|u| u + 1.25).collect();
        let n_k = [100, 100];
        let r = mbar(&[u0, u1], &n_k, 10000, 1e-12).unwrap();
        assert!((r.delta_f[0][1] - 1.25).abs() < 1e-10, "{:?}", r.delta_f);
        assert!(r.d_delta_f[0][1].abs() < 1e-6);
    }
}
