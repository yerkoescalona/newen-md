//! tser — time series of properties (gromos++ `tser`).
//!
//! Usage: tser @topo <top> @pbc <r|v> [cog|nog] [@time <t0 dt>] @prop <properties> @traj <trc…>
//!             [@nots] [@dist <steps [min max] [periodic]>] [@norm] [@solv] [@skip <n>] [@stride <n>]
//!
//! Properties are gromos++ property specifiers (`crate::property`): `d%1:1,2`,
//! `a%1:1,2,3`, `t%1:1,2,3,4`, `tp%…`, `o%atom(…)%cart(…)`, `op%…`, `pr%…`, `pa%…`.

use crate::args::Arguments;
use crate::distribution::{trim_float, Distribution};
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::PropertyContainer;
use crate::time::Time;
use newen_io::args::{Program, Status};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# tser
\t@topo      <molecular topology file>
\t@pbc       <boundary type> [<gathermethod>]
\t@time      <time and dt>
\t@prop      <properties>
\t[@nots     (do not write time series)]
\t[@dist     <steps [min max] [periodic]>]
\t[@norm     (normalise distribution)]
\t[@solv     (read in solvent)]
\t@traj      <trajectory files>
\t[@skip     <skip n first frames>]
\t[@stride   <take every n-th frame>]";

pub const KNOWN: &[&str] = &[
    "topo", "pbc", "time", "prop", "traj", "skip", "stride", "nots", "dist", "norm", "solv",
];

pub const PROGRAM: Program = Program {
    name: "tser",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let mut time = Time::from_args(&args)?;
    let do_tser = args.count("nots") < 0;
    let normalize = args.count("norm") >= 0;
    let dist = args.values("dist");
    let do_dist = !dist.is_empty();
    let mut dist_steps = 0usize;
    let mut bounds: Option<(f64, f64)> = None;
    let mut periodic = false;
    if do_dist {
        dist_steps = dist[0]
            .parse::<usize>()
            .map_err(|_| format!("distribution: wrong number of steps specified {}", dist[0]))?;
        if dist_steps == 0 {
            return Err(format!(
                "distribution: wrong number of steps specified {}",
                dist[0]
            ));
        }
        if dist.len() >= 3 {
            let lo: f64 = dist[1]
                .parse()
                .map_err(|_| "distribution: bad minimum".to_string())?;
            let hi: f64 = dist[2]
                .parse()
                .map_err(|_| "distribution: bad maximum".to_string())?;
            if hi < lo {
                return Err("distribution: maximum value smaller than minimum value".into());
            }
            bounds = Some((lo, hi));
        }
        if let Some(kw) = dist.get(3) {
            if kw == "periodic" {
                periodic = true;
            } else {
                return Err(format!("distribution: keyword {kw} not known"));
            }
        }
    }
    let mut topo =
        build_topology(read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?);
    // gromos++ reads every atom of the frames (`select("ALL")`): solvate to the first frame
    if let Some(first) = args.values("traj").first() {
        if let Ok(Some(frame)) = TrajectoryReader::new(first)
            .map_err(|e| format!("@traj {first}: {e}"))?
            .read_frame()
        {
            solvate_to_atoms(&mut topo, frame.positions.len()).map_err(|e| e.to_string())?;
        }
    }
    let pbc = Pbc::from_args(&args)?;
    let mut props = PropertyContainer::parse(args.values("prop"), &topo)?;
    let skip: usize = args.get("skip", 0)?;
    let stride: usize = args.get("stride", 1)?.max(1);
    let mut text = String::new();
    if do_tser {
        text.push_str(&format!("#\n#{:>9}\t\t{}\n", "time", props.title(&topo)));
    }
    let mut frame_no = 0usize;
    for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let use_frame = frame_no >= skip && (frame_no - skip).is_multiple_of(stride);
            frame_no += 1;
            if !use_frame {
                continue;
            }
            let t = time.next(frame.time);
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            let values = props.calc(&pos, &periodicity);
            if do_tser {
                text.push_str(&Time::fmt(t));
                text.push_str("\t\t");
                for v in &values {
                    text.push_str(&trim_float(*v));
                    text.push_str("\t\t");
                }
                text.push('\n');
            }
        }
    }
    if do_tser {
        text.push_str("# Averages over run: (<average> <rmsd> <error estimate>)\n");
        for p in &props.props {
            text.push_str(&format!("# {}\t{}\n", p.title(&topo), p.average()));
        }
    }
    if do_dist {
        for p in &props.props {
            let values = p.stat.values();
            let d = match bounds {
                Some((lo, hi)) => {
                    Distribution::over_values_in(values, lo, hi, dist_steps, periodic)?
                },
                None => Distribution::over_values(values, dist_steps)?,
            };
            text.push_str("\n#\n");
            text.push_str(&format!("# Distribution of      {}\n", p.title(&topo)));
            text.push_str(&format!("# values:              {}\n", p.stat.n()));
            text.push_str(&format!("# average value:       {:.6}\n", p.stat.ave()));
            text.push_str(&format!("# rmsd:                {:.6}\n", p.stat.rmsd()));
            let ee = p.stat.ee_strict();
            text.push_str(&format!(
                "# error estimate:      {}\n",
                if ee.is_nan() {
                    "-nan".to_string()
                } else {
                    format!("{:.6}", ee)
                }
            ));
            text.push_str(&format!("# maximum value at:    {:.6}\n", d.max_val_at()));
            let fixed6 = |x: f64| format!("{:.6}", x);
            text.push_str(&if normalize {
                d.write_normalized(&fixed6)
            } else {
                d.write(&fixed6)
            });
        }
    }
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
