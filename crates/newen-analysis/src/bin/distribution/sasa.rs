//! sasa - Solvent Accessible Surface Area calculation
//!
//! Usage: sasa @topo <file> @traj <file> @probe <radius>
//!
//! Calculates SASA using approximate sphere-based method

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: sasa @topo <file> @traj <file> [@probe <radius>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut probe_radius = 0.14f64; // Default water probe 1.4 Å = 0.14 nm

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

    eprintln!("# SASA calculation");
    eprintln!("# Probe radius: {} nm", probe_radius);
    eprintln!("# Atoms: {}", num_atoms);

    // Default VDW radii (approximate, in nm)
    let default_radius = 0.15f64;

    println!("# Time (ps)    SASA (nm²)    SASA_per_atom (nm²)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                // Simple approximate SASA: sum of exposed sphere areas
                let mut total_sasa = 0.0f64;

                for i in 0..num_atoms {
                    let radius_i = default_radius + probe_radius;
                    let mut exposed_area = 4.0 * std::f64::consts::PI * radius_i * radius_i;

                    // Check overlap with other atoms
                    let mut n_neighbors = 0;
                    for j in 0..num_atoms {
                        if i != j {
                            let dx = frame.positions[i].x - frame.positions[j].x;
                            let dy = frame.positions[i].y - frame.positions[j].y;
                            let dz = frame.positions[i].z - frame.positions[j].z;
                            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                            let radius_j = default_radius + probe_radius;
                            let sum_radii = radius_i + radius_j;

                            if dist < sum_radii {
                                n_neighbors += 1;
                            }
                        }
                    }

                    // Simple reduction based on neighbors
                    if n_neighbors > 0 {
                        let reduction = (n_neighbors as f64 * 0.1).min(0.9);
                        exposed_area *= 1.0 - reduction;
                    }

                    total_sasa += exposed_area;
                }

                let sasa_per_atom = total_sasa / num_atoms as f64;

                println!(
                    "{:12.4} {:12.4} {:18.6}",
                    frame.time, total_sasa, sasa_per_atom
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
