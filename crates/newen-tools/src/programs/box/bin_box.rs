//! bin_box - Bin coordinates into grid cells
//!
//! Usage: bin_box @topo <file> @traj <file> @bins <nx ny nz>
//!
//! Divides simulation box into bins and counts atoms per bin

use crate::programs::io_err;
use newen::io::topology::{build_topology, read_topology_file};
use newen::io::trajectory::TrajectoryReader;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "bin_box - Count atoms per grid cell of the box\n\nUsage: bin_box @topo <file> @traj <file> @bins <nx> <ny> <nz>\n\nArguments:\n  @topo  Molecular topology file\n  @traj  Trajectory file\n  @bins  Number of bins per dimension (three values)\n\nNote: gromos++'s bin_box builds a binary mixture; this one bins a trajectory.";

pub const KNOWN: &[&str] = &["topo", "traj", "bins"];

pub const PROGRAM: Program = Program {
    name: "bin_box",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let traj_file = args.value("traj")?;

    let bins = args.values("bins");
    if bins.len() != 3 {
        return Err(format!("@bins takes three values, got {}", bins.len()));
    }
    let parse = |v: &String| -> Result<usize, String> {
        match v.parse::<usize>() {
            Ok(0) | Err(_) => Err(format!("@bins: '{v}' is not a positive integer")),
            Ok(n) => Ok(n),
        }
    };
    let (nx, ny, nz) = (parse(&bins[0])?, parse(&bins[1])?, parse(&bins[2])?);

    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;
    let _topo = build_topology(topo_data);
    let mut traj =
        TrajectoryReader::new(traj_file).map_err(|e| format!("@traj '{traj_file}': {e}"))?;

    eprintln!("# Binning into {} x {} x {} grid", nx, ny, nz);

    writeln!(out, "# Time (ps)    Bin_X  Bin_Y  Bin_Z  Count").map_err(io_err)?;

    if let Ok(Some(frame)) = traj.read_frame() {
        let box_dims = frame.box_dims;
        let dx = box_dims.x / nx as f64;
        let dy = box_dims.y / ny as f64;
        let dz = box_dims.z / nz as f64;

        // Count atoms in each bin
        let mut bins = vec![vec![vec![0; nz]; ny]; nx];

        for pos in &frame.positions {
            let ix = ((pos.x / dx).floor() as usize).min(nx - 1);
            let iy = ((pos.y / dy).floor() as usize).min(ny - 1);
            let iz = ((pos.z / dz).floor() as usize).min(nz - 1);
            bins[ix][iy][iz] += 1;
        }

        // Print non-empty bins
        #[allow(clippy::needless_range_loop)] // one index over several arrays
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    if bins[ix][iy][iz] > 0 {
                        writeln!(
                            out,
                            "{:12.4} {:6} {:6} {:6} {:6}",
                            frame.time, ix, iy, iz, bins[ix][iy][iz]
                        )
                        .map_err(io_err)?;
                    }
                }
            }
        }
    }
    Ok(())
}
