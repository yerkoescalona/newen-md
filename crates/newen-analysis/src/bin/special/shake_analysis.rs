//! shake_analysis - Analyze SHAKE constraint violations
//!
//! Usage: shake_analysis @topo <file> @traj <file>
//!
//! Analyzes bond constraint violations

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: shake_analysis @topo <file> @traj <file> [@tol <tolerance>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut tolerance = 0.001f64; // nm

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@tol" => {
                i += 1;
                tolerance = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# SHAKE constraint analysis");
    eprintln!("# Tolerance: {} nm", tolerance);
    eprintln!("# Bonds: {}", topo.moltypes[0].bonds.len());

    println!("# Time (ps)    Violations    Max_Error (nm)    Avg_Error (nm)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut violations = 0;
                let mut max_error = 0.0f64;
                let mut total_error = 0.0f64;
                let mut bond_count = 0;

                for bond in &topo.moltypes[0].bonds {
                    let i = bond.i;
                    let j = bond.j;
                    if i < frame.positions.len() && j < frame.positions.len() {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                        // Assume target bond length is 0.1 nm (C-C bond)
                        let target = 0.1f64;
                        let error = (dist - target).abs();

                        if error > tolerance {
                            violations += 1;
                        }

                        max_error = max_error.max(error);
                        total_error += error;
                        bond_count += 1;
                    }
                }

                let avg_error = if bond_count > 0 {
                    total_error / bond_count as f64
                } else {
                    0.0
                };

                println!(
                    "{:12.4} {:12} {:17.6} {:17.6}",
                    frame.time, violations, max_error, avg_error
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
