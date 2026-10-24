//! bilayer_dist — distribution of selected atoms along z relative to the centre of mass of a
//! set of atoms, e.g. head groups relative to a bilayer's centre (gromos++ `bilayer_dist`).
//!
//! Usage: bilayer_dist @topo <top> @pbc <r|v> [cog|nog] [@time <t0 dt>] @atoms <spec>
//!        @selection <spec> [@translate] [@grid <n=100>] [@mult <f=1>] [@density] @traj <trc…>

use crate::args::Arguments;
use crate::distribution::{trim_float, Distribution};
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::ordered_atoms;
use crate::time::Time;
use newen_core::Vec3;
use newen_io::args::{Program, Status};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# bilayer_dist
\t@topo           <molecular topology file>
\t@pbc            <boundary type> [<gathermethod>]
\t[@time           <time and dt>]
\t@atoms          <AtomSpecifier>
\t@selection      <AtomSpecifier>
\t[@translate      (only in case com is not in between layers)]
\t[@grid          <integer; default: 100>]
\t@traj           <trajectory files>
\t[@mult          <double; default: 1>]
\t[@density       (in this case the density is calculated)]";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "time",
    "atoms",
    "selection",
    "translate",
    "grid",
    "traj",
    "mult",
    "density",
];

pub const PROGRAM: Program = Program {
    name: "bilayer_dist",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
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
    let mut time = Time::from_args(&args)?;
    let atoms = ordered_atoms(&args.values("atoms").join(";"), &topo)?;
    if atoms.is_empty() {
        return Err("@atoms: no atoms selected".into());
    }
    let grid: usize = args.get("grid", 100)?;
    let selection = ordered_atoms(&args.values("selection").join(";"), &topo)?;
    if selection.is_empty() {
        return Err("@selection: argument is empty!".into());
    }
    let translate = args.count("translate") >= 0;
    let density = args.count("density") >= 0;
    let mult: f64 = args.get("mult", 1.0)?;
    let mut zcoord = Vec::new();
    let (mut min_d, mut max_d) = (0.0, 0.0);
    let mut frames = 0usize;
    let mut area = 0.0;
    writeln!(out, "# DISTRIBUTION: ").map_err(io_err)?;
    for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            time.next(frame.time);
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            let mut cm = Vec3::ZERO;
            let mut mass = 0.0;
            for &a in &atoms {
                cm += pos[a] * topo.mass[a];
                mass += topo.mass[a];
            }
            cm /= mass;
            if translate {
                cm.z -= frame.box_dims.z / 2.0;
            }
            min_d = -frame.box_dims.z / 2.0;
            max_d = frame.box_dims.z / 2.0;
            area = frame.box_dims.x * frame.box_dims.y;
            for &s in &selection {
                // gromos++: pos − nearestImage(pos, cm) = the minimum-image vector pos − cm
                zcoord.push(periodicity.nearest_image(pos[s], cm).z);
            }
            frames += 1;
        }
    }
    let mut dist = Distribution::new(min_d, max_d, grid)?;
    for &z in &zcoord {
        dist.add(z);
    }
    let binwidth = (max_d - min_d) / grid as f64;
    if !density {
        for k in 0..grid {
            writeln!(
                out,
                "{}\t{}",
                trim_float(dist.centre(k)),
                trim_float(dist.count(k) as f64 * mult / (binwidth * dist.n_values() as f64))
            )
            .map_err(io_err)?;
        }
    } else {
        for k in 0..grid {
            writeln!(
                out,
                "{}{:>15}",
                trim_float(dist.centre(k)),
                trim_float(dist.count(k) as f64 * mult / (area * binwidth * frames as f64))
            )
            .map_err(io_err)?;
        }
    }
    Ok(())
}
