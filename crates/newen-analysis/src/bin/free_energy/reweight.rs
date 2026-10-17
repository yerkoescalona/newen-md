//! reweight - Reweight trajectory frames
//!
//! Usage: reweight @traj <file> @weights <file>
//!
//! Applies statistical reweighting to trajectory frames

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: reweight @traj <file> @weights <file>");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut weights_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@weights" => {
                i += 1;
                weights_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    // Read weights
    let file = File::open(weights_file.unwrap()).unwrap();
    let reader = BufReader::new(file);
    let weights: Vec<f64> = reader
        .lines()
        .map_while(Result::ok)
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .filter_map(|line| line.trim().parse().ok())
        .collect();

    eprintln!("# Loaded {} weights", weights.len());

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    println!("# Frame    Time (ps)    Weight");

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let weight = if frame_idx < weights.len() {
                    weights[frame_idx]
                } else {
                    1.0
                };

                println!("{:8} {:12.4} {:12.6}", frame_idx + 1, frame.time, weight);
                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Processed {} frames", frame_idx);
}
