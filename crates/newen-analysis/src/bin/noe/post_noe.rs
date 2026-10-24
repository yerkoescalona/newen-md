//! post_noe - Post-process NOE data
//!
//! Usage: post_noe @noe <file> @traj <file>
//!
//! Post-processes NOE restraint violations

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: post_noe @noe <file> @traj <file>");
        process::exit(1);
    }

    let mut noe_file = None;
    let mut traj_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@noe" => {
                i += 1;
                noe_file = Some(args[i].clone());
            },
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    // Read NOE restraints
    let mut noe_pairs: Vec<(usize, usize, f64)> = Vec::new();

    if let Some(path) = noe_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                let parts: Vec<&str> = line_str.split_whitespace().collect();
                if parts.len() >= 3 {
                    if let (Ok(i), Ok(j), Ok(target)) = (
                        parts[0].parse::<usize>(),
                        parts[1].parse::<usize>(),
                        parts[2].parse::<f64>(),
                    ) {
                        noe_pairs.push((i - 1, j - 1, target));
                    }
                }
            }
        }
    }

    if noe_pairs.is_empty() {
        eprintln!("Warning: No NOE restraints loaded");
        process::exit(1);
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# NOE post-processing");
    eprintln!("# Restraints: {}", noe_pairs.len());

    println!("# Frame    Time (ps)    Violations    Max_Viol (nm)    Avg_Viol (nm)");

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut violations = 0;
                let mut max_viol = 0.0f64;
                let mut total_viol = 0.0f64;

                for &(i, j, target) in &noe_pairs {
                    if i < frame.positions.len() && j < frame.positions.len() {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                        let viol = (dist - target).abs();
                        if viol > 0.05 {
                            violations += 1;
                            max_viol = max_viol.max(viol);
                            total_viol += viol;
                        }
                    }
                }

                let avg_viol = if violations > 0 {
                    total_viol / violations as f64
                } else {
                    0.0
                };

                println!(
                    "{:8} {:12.4} {:12} {:15.6} {:15.6}",
                    frame_idx + 1,
                    frame.time,
                    violations,
                    max_viol,
                    avg_viol
                );

                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Analyzed {} frames", frame_idx);
}
