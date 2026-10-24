//! propertyaver - Average molecular properties
//!
//! Usage: propertyaver @traj <file> @property <name>
//!
//! Averages molecular properties over trajectory

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: propertyaver @traj <file> [@property <volume|box_x|box_y|box_z>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut property = "volume";

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@property" => {
                i += 1;
                property = &args[i];
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Property averaging");
    eprintln!("# Property: {}", property);

    let mut values = Vec::new();

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let value = match property {
                    "box_x" => frame.box_dims.x,
                    "box_y" => frame.box_dims.y,
                    "box_z" => frame.box_dims.z,
                    "natoms" => frame.positions.len() as f64,
                    _ => frame.box_dims.x * frame.box_dims.y * frame.box_dims.z,
                };

                values.push(value);
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    if values.is_empty() {
        eprintln!("Error: No data");
        process::exit(1);
    }

    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
    let std_dev = variance.sqrt();
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    println!("# Property: {}", property);
    println!("# Frames: {}", values.len());
    println!("# Mean: {:.6}", mean);
    println!("# Std Dev: {:.6}", std_dev);
    println!("# Min: {:.6}", min);
    println!("# Max: {:.6}", max);
}
