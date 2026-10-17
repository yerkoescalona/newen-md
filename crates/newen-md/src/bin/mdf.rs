//! mdf - Mean force analysis
//!
//! Usage: mdf @traj <file> @coord <dimension>
//!
//! Calculates mean force along a reaction coordinate

use newen::io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: mdf @traj <file> [@coord <x|y|z>] [@bins <n>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut coord = "x";
    let mut n_bins = 50;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@coord" => {
                i += 1;
                coord = &args[i];
            },
            "@bins" => {
                i += 1;
                n_bins = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Mean force analysis");
    eprintln!("# Coordinate: {}", coord);
    eprintln!("# Bins: {}", n_bins);

    let mut histogram = vec![0usize; n_bins];
    let mut min_coord = f64::MAX;
    let mut max_coord = f64::MIN;

    // First pass: determine range
    let mut all_coords = Vec::new();
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                if let Some(first_pos) = frame.positions.first() {
                    let coord_val = match coord {
                        "y" => first_pos.y,
                        "z" => first_pos.z,
                        _ => first_pos.x,
                    };
                    all_coords.push(coord_val);
                    min_coord = min_coord.min(coord_val);
                    max_coord = max_coord.max(coord_val);
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    let range = max_coord - min_coord;
    let bin_width = range / n_bins as f64;

    // Build histogram
    for &coord_val in &all_coords {
        let bin = ((coord_val - min_coord) / bin_width) as usize;
        if bin < n_bins {
            histogram[bin] += 1;
        }
    }

    println!("# Coordinate    Count    PMF (kJ/mol)");

    let kb_t = 2.5; // kJ/mol at 300K
    let max_count = *histogram.iter().max().unwrap_or(&1);

    for (idx, &count) in histogram.iter().enumerate() {
        let coord_val = min_coord + (idx as f64 + 0.5) * bin_width;
        let pmf = if count > 0 {
            -kb_t * (count as f64 / max_count as f64).ln()
        } else {
            0.0
        };

        println!("{:12.6} {:8} {:12.4}", coord_val, count, pmf);
    }
}
