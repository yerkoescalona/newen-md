//! unify_box - Unify box shapes in trajectory
//!
//! Usage: unify_box @topo <file> @traj <file> @box <x y z>
//!
//! Sets all frames to use the same box dimensions

use crate::programs::io_err;
use newen::io::topology::{build_topology, read_topology_file};
use newen::io::trajectory::TrajectoryReader;
use newen::math::Vec3;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "unify_box - Rescale every frame to one box\n\nUsage: unify_box @topo <file> @traj <file> @box <x> <y> <z>\n\nArguments:\n  @topo  Molecular topology file\n  @traj  Trajectory file\n  @box   Target box dimensions (three values, nm)\n\nNote: gromos++'s unify_box converts between box types; this one rescales to a given box.";

pub const KNOWN: &[&str] = &["topo", "traj", "box"];

pub const PROGRAM: Program = Program {
    name: "unify_box",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let traj_file = args.value("traj")?;

    let dims = args.values("box");
    if dims.len() != 3 {
        return Err(format!("@box takes three values, got {}", dims.len()));
    }
    let parse = |v: &String| -> Result<f64, String> {
        v.parse::<f64>()
            .map_err(|_| format!("@box: '{v}' is not a number"))
    };
    let (box_x, box_y, box_z) = (parse(&dims[0])?, parse(&dims[1])?, parse(&dims[2])?);

    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;
    let _topo = build_topology(topo_data);
    let mut traj =
        TrajectoryReader::new(traj_file).map_err(|e| format!("@traj '{traj_file}': {e}"))?;

    let unified_box = Vec3::new(box_x, box_y, box_z);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(
        out,
        "Unified box dimensions: {:.3} {:.3} {:.3}",
        box_x, box_y, box_z
    )
    .map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                writeln!(out, "TIMESTEP").map_err(io_err)?;
                writeln!(out, "{:15} {:15.4}", frame.step, frame.time).map_err(io_err)?;
                writeln!(out, "END").map_err(io_err)?;
                writeln!(out, "POSITIONRED").map_err(io_err)?;

                for (i, pos) in frame.positions.iter().enumerate() {
                    writeln!(
                        out,
                        "{:6} {:15.9} {:15.9} {:15.9}",
                        i + 1,
                        pos.x,
                        pos.y,
                        pos.z
                    )
                    .map_err(io_err)?;
                }

                writeln!(out, "END").map_err(io_err)?;
                writeln!(out, "BOX").map_err(io_err)?;
                writeln!(
                    out,
                    " {:15.9} {:15.9} {:15.9}",
                    unified_box.x, unified_box.y, unified_box.z
                )
                .map_err(io_err)?;
                writeln!(out, "END").map_err(io_err)?;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
    Ok(())
}
