//! check_box - Check simulation box dimensions
//!
//! Usage: check_box @traj <file>
//!
//! Validates box dimensions throughout trajectory

use crate::programs::io_err;
use newen::io::trajectory::TrajectoryReader;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "check_box - Box statistics over a trajectory\n\nUsage: check_box @traj <file>\n\nArguments:\n  @traj  Trajectory file\n\nNote: gromos++'s check_box reports close atom pairs; this one reports box dimensions and volume.";

pub const KNOWN: &[&str] = &["traj"];

pub const PROGRAM: Program = Program {
    name: "check_box",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let traj_file = args.value("traj")?;

    let mut traj =
        TrajectoryReader::new(traj_file).map_err(|e| format!("@traj '{traj_file}': {e}"))?;

    writeln!(
        out,
        "# Frame    Time (ps)    Box_X (nm)    Box_Y (nm)    Box_Z (nm)    Volume (nm³)"
    )
    .map_err(io_err)?;

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let box_dims = frame.box_dims;
                let volume = box_dims.x * box_dims.y * box_dims.z;

                writeln!(
                    out,
                    "{:8} {:12.4} {:12.6} {:12.6} {:12.6} {:12.4}",
                    frame_idx + 1,
                    frame.time,
                    box_dims.x,
                    box_dims.y,
                    box_dims.z,
                    volume
                )
                .map_err(io_err)?;
                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Checked {} frames", frame_idx);
    Ok(())
}
