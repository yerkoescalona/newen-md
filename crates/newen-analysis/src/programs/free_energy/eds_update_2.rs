//! eds_update_2 — update the EDS parameters s and EᵢR of two states from the energy time series
//! of a reference-state simulation (gromos++ `eds_update_2`).
//!
//! Usage: eds_update_2 @temp <T> @vr <file> @vy <2 files> @s <current s> @s_old <old s>
//!        @EiR <2 values> @update <1|2> @eunder <E> [@etrans <E>] [@scale <f>]
//!
//! Update 1 scales s by the occupancy ratio of the two states (or up when both average energies
//! lie below `@eunder`); update 2 by the populations of the two states and the transition region
//! `@eunder ± @etrans` of V_B − V_A. The time series are `time energy` columns; `#` starts a
//! comment.

use crate::args::Arguments;
use crate::distribution::cpp_g;
use crate::eds::{calculate_eir, read_series, Pair};
use crate::programs::io_err;
use crate::GROMOSPP_BOLTZMANN;
use newen_io::args::{Program, Status};
use std::io::Write;

pub const USAGE: &str = "# eds_update_2
\t@temp          <temperature for perturbation>
\t@vr             <energy time series of state R>
\t@vy             <energy time series of states Y (2 files)>
\t@s              <current s parameter>
\t@s_old          <old s parameter>
\t@EiR            <energy offset parameters (2 values)>
\t@update         <which update scheme should be used>
\t@eunder         <energy threshold if update=1; separation energy if update=2>
\t@etrans         <ignored if update=1; size of transition region if update=2>
\t[@scale         <scaling factor to modify default factors>]";

pub const KNOWN: &[&str] = &[
    "temp", "vr", "vy", "s", "s_old", "EiR", "update", "eunder", "etrans", "scale",
];

pub const PROGRAM: Program = Program {
    name: "eds_update_2",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

fn g(x: f64) -> String {
    cpp_g(x, 6)
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let temp: f64 = args.require("temp")?;
    let beta = 1.0 / (GROMOSPP_BOLTZMANN * temp);
    let scale: f64 = args.get("scale", 1.0)?;
    let update: i32 = args.require("update")?;
    if update != 1 && update != 2 {
        return Err("update must be 1 or 2".into());
    }
    let eunder: f64 = args.require("eunder")?;
    let etrans: f64 = if update == 2 {
        args.require("etrans")?
    } else {
        0.0
    };
    let eir_vals = args.values("EiR");
    if eir_vals.len() < 2 {
        return Err("not enough EiR values (must be 2)".into());
    }
    if eir_vals.len() > 2 {
        return Err("EiR not numeric or too many values".into());
    }
    let mut eir_old: Vec<f64> = eir_vals
        .iter()
        .map(|v| {
            v.parse::<f64>()
                .map_err(|_| "EiR not numeric or too many values".to_string())
        })
        .collect::<Result<_, _>>()?;
    let s_cur: f64 = args
        .require("s")
        .map_err(|_| "current s parameter not numeric".to_string())?;
    let s_old: f64 = args
        .require("s_old")
        .map_err(|_| "old s parameter not numeric".to_string())?;
    let vr = read_series("vr", args.value("vr")?)?;
    let vy: Vec<Vec<f64>> = args
        .values("vy")
        .iter()
        .map(|p| read_series("vy", p))
        .collect::<Result<_, _>>()?;
    if vy.len() < 2 {
        return Err("not enough energy time series vy".into());
    }
    if vy[0].len() != vr.len() {
        return Err("Time series files differ in length!\n".into());
    }
    let n = vr.len();
    let mut text = String::from("# ############### New EDS Parameters: #################\n");
    // gromos++ calculate_s
    let s_new = if update == 1 {
        let mut occup = [0usize; 2];
        for k in 0..n {
            let state = if vy[1][k] < vy[0][k] { 1 } else { 0 };
            occup[state] += 1;
        }
        text.push_str(&format!(
            "# state0 = {} # state1 = {} \n",
            occup[0], occup[1]
        ));
        let downscale = 1.0 - ((scale * 1.1 - 1.0) / 2.0);
        let ave_a = vy[0].iter().fold(0.0, |a, v| a + v) / n as f64;
        let ave_b = vy[1].iter().fold(0.0, |a, v| a + v) / n as f64;
        text.push_str(&format!(
            "# s_old = {}, current s = {}\n",
            g(s_old),
            g(s_cur)
        ));
        if ave_a < eunder && ave_b < eunder {
            s_cur * scale * 1.1
        } else {
            // the state with the higher offset is the numerator of the occupancy ratio
            let (a, b) = if eir_old[0] > eir_old[1] {
                (occup[0], occup[1])
            } else {
                (occup[1], occup[0])
            };
            let ratio = if b > 0 { a as f64 / b as f64 } else { a as f64 };
            if ratio * 1.1 < 1.0 && ratio * 0.9 < 1.0 && a != 0 {
                s_cur * scale * 1.1
            } else if (ratio * 1.1 > 1.0 && ratio * 0.9 > 1.0) || a == 0 {
                if 2.0 * s_cur == s_old {
                    s_cur / 2.0
                } else {
                    s_cur * downscale
                }
            } else {
                s_cur
            }
        }
    } else {
        let (mut state_a, mut state_b, mut state_t, mut numtrans) =
            (0usize, 0usize, 0usize, 0usize);
        let (mut visit, mut visitold) = (0, 0);
        for i in 0..n {
            let d = vy[1][i] - vy[0][i];
            if d > eunder + etrans {
                state_a += 1;
                visit = 0;
            } else if d < eunder - etrans {
                state_b += 1;
                visit = 1;
            } else {
                state_t += 1;
            }
            if i > 0 && visitold != visit {
                numtrans += 1;
            }
            visitold = visit;
        }
        let mintrans = 10;
        let s = if state_t > state_a && state_t > state_b {
            s_cur * 1.1 * scale
        } else if (state_t >= state_a || state_t >= state_b) && numtrans > mintrans {
            s_cur * 1.05 * scale
        } else if (state_t < state_a || state_t < state_b) && numtrans < mintrans {
            s_cur * 0.95 * scale
        } else if state_t as f64 * 3.0 > state_a as f64
            && state_t as f64 * 3.0 > state_b as f64
            && numtrans > mintrans
        {
            s_cur * 1.02 * scale
        } else {
            s_cur
        };
        text.push_str(&format!(
            "# state A = {state_a}, state B = {state_b}, state T = {state_t}\n# number of A-B and B-A transitions = {numtrans}\n# s_old = {}, current s = {}\n",
            g(s_old),
            g(s_cur)
        ));
        s
    };
    let pair = Pair {
        i: 0,
        j: 1,
        s_old,
        s_new,
        eir: eir_old[0],
        ejr: eir_old[1],
        ..Default::default()
    };
    let pairs = [pair];
    let mut eir_new = [0.0; 2];
    for j in 0..2 {
        eir_new[j] = calculate_eir(&pairs, &|p| p.s_new, &vr, &vy, &eir_old, beta, j, 1.0);
    }
    let e0 = eir_new[0];
    for e in eir_new.iter_mut() {
        *e -= e0;
    }
    eir_old.copy_from_slice(&eir_new);
    text.push_str(&format!("# new s parameter\n{}\n\n", g(s_new)));
    text.push_str(&format!(
        "# new energy offsets: \n{} {} \n\n",
        g(eir_old[0]),
        g(eir_old[1])
    ));
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
