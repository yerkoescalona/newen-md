//! epath - Path ensemble analysis
//!
//! Usage: epath @traj <file> @ref1 <file> @ref2 <file>
//!
//! Analyzes transition paths between two reference structures

use newen_core::math::Vec3;
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn read_reference_structure(filename: &str) -> Vec<Vec3> {
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

fn calc_rmsd(pos1: &[Vec3], pos2: &[Vec3]) -> f64 {
    let mut sum_sq = 0.0f64;
    let n = pos1.len().min(pos2.len());

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

    if args.len() < 7 {
        eprintln!("Usage: epath @traj <file> @ref1 <file> @ref2 <file>");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut ref1_file = None;
    let mut ref2_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@ref1" => {
                i += 1;
                ref1_file = Some(args[i].clone());
            },
            "@ref2" => {
                i += 1;
                ref2_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let ref1_pos = read_reference_structure(&ref1_file.unwrap());
    let ref2_pos = read_reference_structure(&ref2_file.unwrap());

    if ref1_pos.is_empty() || ref2_pos.is_empty() {
        eprintln!("Error: Could not load reference structures");
        process::exit(1);
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Path ensemble analysis");
    eprintln!("# Reference 1: {} atoms", ref1_pos.len());
    eprintln!("# Reference 2: {} atoms", ref2_pos.len());

    println!("# Frame    Time (ps)    RMSD_ref1 (nm)    RMSD_ref2 (nm)    Progress");

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let rmsd1 = calc_rmsd(&frame.positions, &ref1_pos);
                let rmsd2 = calc_rmsd(&frame.positions, &ref2_pos);

                // Progress coordinate: 0 at ref1, 1 at ref2
                let progress = rmsd1 / (rmsd1 + rmsd2);

                println!(
                    "{:8} {:12.4} {:17.6} {:17.6} {:12.6}",
                    frame_idx + 1,
                    frame.time,
                    rmsd1,
                    rmsd2,
                    progress
                );

                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Analyzed {} frames", frame_idx);
}
