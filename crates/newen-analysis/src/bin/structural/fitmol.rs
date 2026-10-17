//! fitmol - Fit molecule to reference structure
//!
//! Usage: fitmol @ref <file> @traj <file>
//!
//! Performs least-squares fitting to reference

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

fn calc_centroid(positions: &[Vec3]) -> Vec3 {
    let mut sum = Vec3::ZERO;
    for pos in positions {
        sum += *pos;
    }
    sum / (positions.len() as f64)
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: fitmol @ref <file> @traj <file>");
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

    let ref_centroid = calc_centroid(&ref_pos);

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Molecular fitting");
    eprintln!("# Reference atoms: {}", ref_pos.len());

    println!("TITLE");
    println!("  Fitted structures");
    println!("END");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let n = frame.positions.len().min(ref_pos.len());
                let frame_centroid = calc_centroid(&frame.positions[0..n]);

                println!("POSITION");
                println!(
                    "  {:15.9} {:15.9} {:15.9}",
                    frame.box_dims.x, frame.box_dims.y, frame.box_dims.z
                );

                // Translate to match centroids
                for i in 0..n {
                    let shifted = frame.positions[i] - frame_centroid + ref_centroid;
                    println!("{:15.9} {:15.9} {:15.9}", shifted.x, shifted.y, shifted.z);
                }

                println!("END");
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
