//! sasa_hasel - SASA using Hasel algorithm
//!
//! Usage: sasa_hasel @topo <file> @traj <file>
//!
//! Calculates solvent accessible surface area using Hasel method

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: sasa_hasel @topo <file> @traj <file> [@probe <r>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut probe_radius = 0.14f64; // 1.4 Å

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
            "@probe" => {
                i += 1;
                probe_radius = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let num_atoms = topo.mass.len();

    eprintln!("# SASA (Hasel method)");
    eprintln!("# Probe radius: {} nm", probe_radius);
    eprintln!("# Atoms: {}", num_atoms);

    println!("# Time (ps)    SASA (nm²)");

    // Default atomic radii (Hasel)
    let default_radius = 0.15f64;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut total_sasa = 0.0f64;

                // Hasel approximation: pi * sum(p_i * r_i²)
                // where p_i is atomic accessibility parameter
                for i in 0..num_atoms {
                    let r_i = default_radius + probe_radius;

                    // Calculate accessibility (simplified)
                    let mut neighbors = 0;
                    for j in 0..num_atoms {
                        if i != j {
                            let dx = frame.positions[i].x - frame.positions[j].x;
                            let dy = frame.positions[i].y - frame.positions[j].y;
                            let dz = frame.positions[i].z - frame.positions[j].z;
                            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                            if dist < 0.5 {
                                // Nearby atom
                                neighbors += 1;
                            }
                        }
                    }

                    // Accessibility parameter (decreases with neighbors)
                    let p_i = if neighbors > 0 {
                        1.0 / (1.0 + 0.1 * neighbors as f64)
                    } else {
                        1.0
                    };

                    total_sasa += std::f64::consts::PI * p_i * r_i * r_i;
                }

                println!("{:12.4} {:12.4}", frame.time, total_sasa);
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
