//! gch - Generate constraint history
//!
//! Usage: gch @traj <file> @constraints <file>
//!
//! Generates history of constraint violations

use crate::programs::io_err;
use newen::io::trajectory::TrajectoryReader;
use newen_io::args::{Arguments, Program, Status};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

pub const USAGE: &str = "gch - Constraint-distance history over a trajectory\n\nUsage: gch @traj <file> [@constraints <file>]\n\nArguments:\n  @traj         Trajectory file\n  @constraints  Constraint list: one \"i j distance\" per line (1-based atoms)\n\nNote: gromos++'s gch regenerates hydrogen coordinates; this one checks constraint distances.";

pub const KNOWN: &[&str] = &["traj", "constraints"];

pub const PROGRAM: Program = Program {
    name: "gch",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let traj_file = args.value("traj")?;
    let const_file = args.has("constraints").then(|| args.value("constraints"));
    let const_file = match const_file {
        Some(r) => Some(r?.to_string()),
        None => None,
    };

    let mut constraints: Vec<(usize, usize, f64)> = Vec::new();

    if let Some(path) = const_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                let parts: Vec<&str> = line_str.split_whitespace().collect();
                if parts.len() >= 3 {
                    if let (Ok(i), Ok(j), Ok(d)) = (
                        parts[0].parse::<usize>(),
                        parts[1].parse::<usize>(),
                        parts[2].parse::<f64>(),
                    ) {
                        constraints.push((i - 1, j - 1, d));
                    }
                }
            }
        }
    }

    if constraints.is_empty() {
        eprintln!("# Warning: No constraints loaded, using default (atoms 1-2, 0.1 nm)");
        constraints.push((0, 1, 0.1));
    }

    let mut traj =
        TrajectoryReader::new(traj_file).map_err(|e| format!("@traj '{traj_file}': {e}"))?;

    eprintln!("# Constraint history");
    eprintln!("# Constraints: {}", constraints.len());

    writeln!(
        out,
        "# Frame    Time (ps)    Constraint    Distance (nm)    Violation (nm)"
    )
    .map_err(io_err)?;

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                for (cidx, &(i, j, target)) in constraints.iter().enumerate() {
                    if i < frame.positions.len() && j < frame.positions.len() {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                        let viol = (dist - target).abs();

                        writeln!(
                            out,
                            "{:8} {:12.4} {:12} {:16.6} {:16.6}",
                            frame_idx + 1,
                            frame.time,
                            cidx + 1,
                            dist,
                            viol
                        )
                        .map_err(io_err)?;
                    }
                }
                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Processed {} frames", frame_idx);
    Ok(())
}
