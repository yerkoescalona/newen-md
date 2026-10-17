//! eds_update_1 — update the EDS parameters s and EᵢR of N states from the energy time series of
//! a reference-state simulation (gromos++ `eds_update_1`).
//!
//! Usage: eds_update_1 @temp <T> @numstat <N> @form <1|2|3> @vr <file> @vy <N files>
//!        @s <s values> @EiR <N values> [@update_tree <0|1>] [@tree <file>]
//!
//! Functional forms: 1 = a single s, 2 = one s per pair (N(N−1)/2), 3 = one s per edge of a
//! maximum spanning tree (N−1 edges read from `@tree`; rebuilt during the iteration with
//! `@update_tree 1`). `tree_new.dat` is written in the working directory as gromos++ does (empty
//! unless form 3). The time series are `time energy` columns; `#` starts a comment.

use crate::args::Arguments;
use crate::distribution::cpp_g;
use crate::eds::{calculate_eir, pair_statistics, read_series, update_s, Pair};
use crate::programs::io_err;
use crate::GROMOSPP_BOLTZMANN;
use newen_io::args::{Program, Status};
use std::io::Write;

pub const USAGE: &str = "# eds_update_1
\t@temp          <temperature for perturbation>
\t@numstat        <number of EDS states (= N)>
\t@form           <functional form (i.e. 1 s, N(N-1)/2 s or (N-1) s)>
\t@vr             <energy time series of state R>
\t@vy             <energy time series of states Y (N files)>
\t@s              <old s parameters (number according to functional form)>
\t@EiR            <old energy offset parameters (N values)>
\t[@update_tree   <updating of maximum spanning tree>]
\t[@tree          <maximum spanning tree (required for @form = 3)>]";

pub const KNOWN: &[&str] = &[
    "temp",
    "numstat",
    "form",
    "vr",
    "vy",
    "s",
    "EiR",
    "tree",
    "update_tree",
];

pub const PROGRAM: Program = Program {
    name: "eds_update_1",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

const MAX_ITER: usize = 300;
const CRIT: f64 = 0.0000001;

fn g(x: f64) -> String {
    cpp_g(x, 6)
}

fn parse_list(values: &[String], what: &str) -> Result<Vec<f64>, String> {
    values
        .iter()
        .map(|v| {
            v.parse::<f64>()
                .map_err(|_| format!("{what} not numeric or too many values"))
        })
        .collect()
}

/// The energy-offset half of one iteration (gromos++'s loop over `calculate_EiR` followed by the
/// shift that keeps E₀R = 0); returns Σ|ΔEᵢR|.
fn offsets_step(
    pairs: &[Pair],
    vr: &[f64],
    vy: &[Vec<f64>],
    eir_old: &mut [f64],
    eir_new: &mut [f64],
    beta: f64,
    factor: f64,
) -> f64 {
    let numstat = eir_old.len();
    for j in 0..numstat {
        eir_new[j] = calculate_eir(pairs, &|p| p.s_old, vr, vy, eir_old, beta, j, factor);
    }
    let e0 = eir_new[0];
    let mut diff = 0.0;
    for j in 0..numstat {
        eir_new[j] -= e0;
        diff += (eir_new[j] - eir_old[j]).abs();
        eir_old[j] = eir_new[j];
    }
    diff
}

/// gromos++ `calculate_tree`: grow a maximum spanning tree over the pair s values, one edge per
/// step, always the largest s that connects a tree state with a state outside it.
fn maximum_spanning_tree(
    pairs: &[Pair],
    treepairs: &mut [usize],
    numstat: usize,
    out: &mut String,
) {
    let mut indicator = vec![1; numstat];
    indicator[0] = 0;
    let max_iter = 500;
    let mut iter = 0;
    let mut p = 0;
    loop {
        let (mut i_max, mut j_max, mut mem_max) = (-5i64, -5i64, 0usize);
        let mut s_max = -100.0;
        for i in 0..numstat {
            for j in 0..numstat {
                if j != i {
                    let mem = if j < i {
                        j * (2 * numstat - j - 1) / 2 + (i - j - 1)
                    } else {
                        i * (2 * numstat - i - 1) / 2 + (j - i - 1)
                    };
                    if s_max < pairs[mem].s_old && indicator[i] != indicator[j] {
                        s_max = pairs[mem].s_old;
                        i_max = i as i64;
                        j_max = j as i64;
                        mem_max = mem;
                    }
                }
            }
        }
        if j_max > 0 && p < treepairs.len() {
            indicator[j_max as usize] = 0;
            indicator[i_max as usize] = 0;
            out.push_str(&format!(
                "i_max = {i_max} , j_max = {j_max}, s_max = {}\n",
                g(s_max)
            ));
            treepairs[p] = mem_max;
            p += 1;
            iter += 1;
        } else {
            break;
        }
        if iter >= max_iter {
            break;
        }
    }
    if iter == max_iter {
        out.push_str("Warning: tree is not converged!\n");
    }
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let temp: f64 = args.require("temp")?;
    let beta = 1.0 / (GROMOSPP_BOLTZMANN * temp);
    let numstat: usize = args.require("numstat")?;
    if numstat < 2 {
        return Err("numstat must be at least 2".into());
    }
    let numpairs1 = numstat * (numstat - 1) / 2;
    let numpairs3 = numstat - 1;
    let form: u32 = args.require("form")?;
    let npairs = match form {
        1 | 2 => numpairs1,
        3 => numpairs3,
        _ => return Err("functional form must be 1,2 or 3".into()),
    };
    let mut pairs = vec![Pair::default(); npairs];
    let mut eir_old = parse_list(args.values("EiR"), "EiR")?;
    if eir_old.len() < numstat {
        return Err("not enough EiR values (must be N)".into());
    }
    if eir_old.len() > numstat {
        return Err("EiR not numeric or too many values".into());
    }
    let mut smallest_s_old = 1.0;
    if form > 1 {
        let sv = parse_list(args.values("s"), "s")?;
        if sv.len() > npairs {
            return Err("s not numeric or too many values".into());
        }
        if form == 2 && sv.len() < numpairs1 {
            return Err("not enough s values (must be N(N-1)/2)".into());
        }
        if form == 3 && sv.len() < numpairs3 {
            return Err("not enough s values (must be (N-1))".into());
        }
        for (p, s) in pairs.iter_mut().zip(&sv) {
            p.s_old = *s;
        }
    } else {
        smallest_s_old = args
            .require("s")
            .map_err(|_| "functional form not numeric".to_string())?;
    }
    // gromos++ opens tree_new.dat before anything else, so it exists (empty) for forms 1 and 2
    std::fs::write("tree_new.dat", "").map_err(|e| format!("tree_new.dat: {e}"))?;
    let mut tree_out = String::new();
    let mut text = String::new();
    let mut tree = 0u32;
    if form == 3 {
        tree = args
            .require("update_tree")
            .map_err(|_| "updating tree switch not numeric".to_string())?;
        if tree != 0 && tree != 1 {
            return Err("update_tree must be 0 or 1".into());
        }
        let path = args.value("tree")?;
        let tree_text =
            std::fs::read_to_string(path).map_err(|_| format!("Could not open file '{path}'"))?;
        let mut i = 0;
        for line in tree_text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let v: Vec<usize> = line
                .split_whitespace()
                .take(2)
                .map(|t| {
                    t.parse::<usize>()
                        .map_err(|_| format!("failed to read values from line\n{line}\ngot\n"))
                })
                .collect::<Result<_, _>>()?;
            if v.len() < 2 {
                return Err(format!("failed to read values from line\n{line}\ngot\n"));
            }
            if i >= numpairs3 {
                return Err("too many pairs in the tree file (must be (N-1))".into());
            }
            pairs[i].i = v[0] - 1;
            pairs[i].j = v[1] - 1;
            text.push_str(&format!("tree: pair = {} and {}\n", v[0], v[1]));
            i += 1;
        }
        if i < numpairs3 {
            return Err("not enough pairs (must be (N-1))".into());
        }
    } else {
        let mut p = 0;
        for j in 0..numstat - 1 {
            for k in j + 1..numstat {
                pairs[p].i = j;
                pairs[p].j = k;
                p += 1;
            }
        }
    }
    let vr = read_series("vr", args.value("vr")?)?;
    let vy: Vec<Vec<f64>> = args
        .values("vy")
        .iter()
        .map(|p| read_series("vy", p))
        .collect::<Result<_, _>>()?;
    if vy.len() < numstat {
        return Err("not enough energy time series vy".into());
    }
    if vy.len() > numstat {
        return Err("too many energy time series vy (must be N)".into());
    }
    if vy[0].len() != vr.len() {
        return Err("Time series files differ in length!\n".into());
    }
    text.push_str("# ############### New EDS Parameters: #################\n");
    let mut eir_new = vec![0.0; numstat];
    match (form, tree) {
        (1, _) => {
            let factor = 1.0 / (numstat - 1) as f64;
            for p in pairs.iter_mut() {
                p.eir = eir_old[p.i];
                p.ejr = eir_old[p.j];
                pair_statistics(p, &vr, &vy, beta);
                p.s_old = smallest_s_old;
            }
            let mut iterations = 0usize;
            let smallest_s_new = loop {
                let diff_eir =
                    offsets_step(&pairs, &vr, &vy, &mut eir_old, &mut eir_new, beta, factor);
                let mut s_min = 1.0;
                for p in pairs.iter_mut() {
                    p.eir = eir_old[p.i];
                    p.ejr = eir_old[p.j];
                    p.term = beta * (p.ejr - p.eir);
                    update_s(p);
                    if p.s_new < s_min {
                        s_min = p.s_new;
                    }
                }
                for p in pairs.iter_mut() {
                    p.s_old = s_min;
                }
                let diff_s = (s_min - smallest_s_old).abs();
                smallest_s_old = s_min;
                if diff_eir < CRIT && diff_s < CRIT {
                    break s_min;
                }
                iterations += 1;
                text.push_str(&format!("iteration: {iterations}\n"));
                if iterations >= MAX_ITER {
                    break s_min;
                }
            };
            text.push_str(&format!(
                "# iterations: {iterations} of {MAX_ITER}\n\n# new s parameter\n{}\n\n",
                g(smallest_s_new)
            ));
        },
        (2, _) | (3, 0) => {
            let factor = if form == 2 {
                1.0 / (numstat - 1) as f64
            } else {
                numstat as f64 / (2 * (numstat - 1)) as f64
            };
            for p in pairs.iter_mut() {
                p.eir = eir_old[p.i];
                p.ejr = eir_old[p.j];
                pair_statistics(p, &vr, &vy, beta);
            }
            let mut iterations = 0usize;
            loop {
                let diff_eir =
                    offsets_step(&pairs, &vr, &vy, &mut eir_old, &mut eir_new, beta, factor);
                let mut diff_s = 0.0;
                for p in pairs.iter_mut() {
                    p.eir = eir_old[p.i];
                    p.ejr = eir_old[p.j];
                    p.term = beta * (p.ejr - p.eir);
                    update_s(p);
                    diff_s += (p.s_new - p.s_old).abs();
                    p.s_old = p.s_new;
                }
                if diff_eir < CRIT && diff_s < CRIT {
                    break;
                }
                iterations += 1;
                text.push_str(&format!("iteration: {iterations}\n"));
                if iterations >= MAX_ITER {
                    break;
                }
            }
            text.push_str(&format!(
                "# iterations: {iterations} of {MAX_ITER}\n\n# new s parameters\n"
            ));
            for p in &pairs {
                text.push_str(&format!("{} ", g(p.s_new)));
                if form == 3 {
                    tree_out.push_str(&format!("{} {} {}\n", p.i + 1, p.j + 1, g(p.s_new)));
                }
            }
            text.push_str(if form == 2 { "\n\n" } else { "\n" });
        },
        _ => {
            // form 3 with the tree rebuilt during the iteration
            let factor = numstat as f64 / (2 * (numstat - 1)) as f64;
            let mut allpairs = vec![Pair::default(); numpairs1];
            let mut treepairs = vec![0usize; numpairs3];
            let (mut p, mut t) = (0, 0);
            for j in 0..numstat - 1 {
                for k in j + 1..numstat {
                    allpairs[p].i = j;
                    allpairs[p].j = k;
                    allpairs[p].s_old = 1.0;
                    for l in 0..numpairs3 {
                        if j == pairs[l].i && k == pairs[l].j && t < numpairs3 {
                            treepairs[t] = p;
                            allpairs[p].s_old = pairs[l].s_old;
                            t += 1;
                        }
                    }
                    p += 1;
                }
            }
            for ap in allpairs.iter_mut() {
                ap.eir = eir_old[ap.i];
                ap.ejr = eir_old[ap.j];
                pair_statistics(ap, &vr, &vy, beta);
            }
            let mut iterations = 0usize;
            loop {
                if iterations == 1 || (iterations.is_multiple_of(MAX_ITER / 10) && iterations != 0)
                {
                    text.push_str(&format!("update tree at iteration: {iterations}\n"));
                    maximum_spanning_tree(&allpairs, &mut treepairs, numstat, &mut text);
                    for j in 0..numpairs3 {
                        pairs[j] = allpairs[treepairs[j]].clone();
                    }
                }
                let diff_eir =
                    offsets_step(&pairs, &vr, &vy, &mut eir_old, &mut eir_new, beta, factor);
                let mut diff_s = 0.0;
                for ap in allpairs.iter_mut() {
                    ap.eir = eir_old[ap.i];
                    ap.ejr = eir_old[ap.j];
                    ap.term = beta * (ap.ejr - ap.eir);
                    update_s(ap);
                    diff_s += (ap.s_new - ap.s_old).abs();
                    ap.s_old = ap.s_new;
                }
                for j in 0..numpairs3 {
                    let ap = &allpairs[treepairs[j]];
                    pairs[j].s_old = ap.s_old;
                    pairs[j].eir = ap.eir;
                    pairs[j].ejr = ap.ejr;
                }
                if diff_eir < CRIT && diff_s < CRIT {
                    break;
                }
                iterations += 1;
                text.push_str(&format!("iteration: {iterations}\n"));
                if iterations >= MAX_ITER {
                    break;
                }
            }
            text.push_str(&format!(
                "# iterations: {iterations} of {MAX_ITER}\n\n# new tree pairs\n"
            ));
            for p in &pairs {
                text.push_str(&format!("{}  {}\n", p.i + 1, p.j + 1));
                tree_out.push_str(&format!("{}  {} {}\n", p.i + 1, p.j + 1, g(p.s_old)));
            }
            text.push_str("# new s parameters\n");
            for p in &pairs {
                text.push_str(&format!("{} ", g(p.s_old)));
            }
            text.push('\n');
        },
    }
    text.push_str("# new energy offsets: \n");
    for e in &eir_old {
        text.push_str(&format!("{} ", g(*e)));
    }
    text.push_str("\n\n");
    std::fs::write("tree_new.dat", tree_out).map_err(|e| format!("tree_new.dat: {e}"))?;
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
