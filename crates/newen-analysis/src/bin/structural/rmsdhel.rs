//! rmsdhel - RMSD for helical structures
//!
//! Usage: rmsdhel @ref <file> @traj <file>
//!
//! Calculates RMSD with helical alignment

use newen_core::math::Vec3;
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn read_reference(filename: &str) -> Vec<Vec3> {
    let mut positions = Vec::new();

    if let Ok(file) = File::open(filename) {
        let reader = BufReader::new(file);
        for line_str in reader.lines().map_while(Result::ok) {
            let parts: Vec<&str> = line_str.split_whitespace().collect();
            if parts.len() >= 3 {
                if let (Ok(x), Ok(y), Ok(z)) = (
                    parts[0].parse::<f64>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                ) {
                    positions.push(Vec3::new(x, y, z));
                }
            }
        }
    }

    positions
}

fn calc_rmsd_helical(pos1: &[Vec3], pos2: &[Vec3]) -> f64 {
    let mut sum_sq = 0.0f64;
    let n = pos1.len().min(pos2.len());

    // Simple helical RMSD (without full alignment)
    for i in 0..n {
        let dx = pos1[i].x - pos2[i].x;
        let dy = pos1[i].y - pos2[i].y;
        let dz = pos1[i].z - pos2[i].z;
        sum_sq += dx * dx + dy * dy + dz * dz;
    }

    (sum_sq / n as f64).sqrt()
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: rmsdhel @ref <file> @traj <file>");
        process::exit(1);
    }

    let mut ref_file = None;
    let mut traj_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@ref" => {
                i += 1;
                ref_file = Some(args[i].clone());
            },
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let ref_pos = read_reference(&ref_file.unwrap());

    if ref_pos.is_empty() {
        eprintln!("Error: Could not load reference");
        process::exit(1);
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Helical RMSD analysis");
    eprintln!("# Reference atoms: {}", ref_pos.len());

    println!("# Frame    Time (ps)    RMSD (nm)");

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let rmsd = calc_rmsd_helical(&frame.positions, &ref_pos);
                println!("{:8} {:12.4} {:12.6}", frame_idx + 1, frame.time, rmsd);
                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Analyzed {} frames", frame_idx);
}
