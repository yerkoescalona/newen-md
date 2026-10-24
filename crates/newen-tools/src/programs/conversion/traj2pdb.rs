//! traj2pdb - Convert trajectory to PDB format
//!
//! Usage: traj2pdb @traj <file> @frame <n>
//!
//! Converts trajectory frame to PDB format

use crate::programs::io_err;
use newen::io::trajectory::TrajectoryReader;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "traj2pdb - Write a trajectory frame as PDB\n\nUsage: traj2pdb @traj <file> [@frame <n>]\n\nArguments:\n  @traj   Trajectory file\n  @frame  Frame to write (1-based, default 1)";

pub const KNOWN: &[&str] = &["traj", "frame"];

pub const PROGRAM: Program = Program {
    name: "traj2pdb",
    known: KNOWN,
    usage: USAGE,
    status: Status::NoCounterpart,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let traj_file = args.value("traj")?;
    let target_frame: usize = args.get("frame", 1)?;

    let mut traj =
        TrajectoryReader::new(traj_file).map_err(|e| format!("@traj '{traj_file}': {e}"))?;

    eprintln!("# Converting trajectory to PDB");
    eprintln!("# Target frame: {}", target_frame);

    let mut frame_idx = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                frame_idx += 1;

                if frame_idx == target_frame {
                    writeln!(
                        out,
                        "REMARK   Frame {} at time {:.2} ps",
                        frame_idx, frame.time
                    )
                    .map_err(io_err)?;
                    writeln!(
                        out,
                        "CRYST1{:9.3}{:9.3}{:9.3}  90.00  90.00  90.00 P 1           1",
                        frame.box_dims.x * 10.0, // Convert nm to Å
                        frame.box_dims.y * 10.0,
                        frame.box_dims.z * 10.0
                    )
                    .map_err(io_err)?;

                    for (idx, pos) in frame.positions.iter().enumerate() {
                        writeln!(out, "ATOM  {:5}  CA  ALA A{:4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C",
                                 idx + 1,
                                 (idx / 100) + 1,
                                 pos.x * 10.0,  // Convert nm to Å
                                 pos.y * 10.0,
                                 pos.z * 10.0).map_err(io_err)?;
                    }

                    writeln!(out, "END").map_err(io_err)?;
                    break;
                }
            },
            Ok(None) => {
                if frame_idx < target_frame {
                    eprintln!(
                        "Warning: Only {} frames available, requested frame {}",
                        frame_idx, target_frame
                    );
                }
                break;
            },
            Err(_) => break,
        }
    }

    eprintln!("# Conversion complete");
    Ok(())
}
