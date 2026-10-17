//! ener_box - Energy in box regions
//!
//! Usage: ener_box @topo <file> @traj <file> @regions <nx ny nz>
//!
//! Calculates energy distribution in box regions

use newen_core::units::four_pi_eps_i;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: ener_box @topo <file> @traj <file> @regions <nx ny nz>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut nx = 2;
    let mut ny = 2;
    let mut nz = 2;

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
            "@regions" => {
                i += 1;
                nx = args[i].parse().unwrap();
                i += 1;
                ny = args[i].parse().unwrap();
                i += 1;
                nz = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Energy in box regions");
    eprintln!("# Regions: {} x {} x {}", nx, ny, nz);

    let n_regions = nx * ny * nz;
    let mut region_energies = vec![0.0f64; n_regions];
    let mut n_frames = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let dx = frame.box_dims.x / nx as f64;
                let dy = frame.box_dims.y / ny as f64;
                let dz = frame.box_dims.z / nz as f64;

                for (idx, pos) in frame.positions.iter().enumerate() {
                    let ix = ((pos.x / dx) as usize).min(nx - 1);
                    let iy = ((pos.y / dy) as usize).min(ny - 1);
                    let iz = ((pos.z / dz) as usize).min(nz - 1);
                    let region_idx = ix * ny * nz + iy * nz + iz;

                    // Simple energy estimate
                    let q = if idx < topo.charge.len() {
                        topo.charge[idx]
                    } else {
                        0.0
                    };

                    let e_self = four_pi_eps_i * q * q / 0.1; // Self energy term
                    region_energies[region_idx] += e_self;
                }

                n_frames += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Frames: {}", n_frames);

    println!("# Region    Energy (kJ/mol)");

    for (idx, &energy) in region_energies.iter().enumerate() {
        let avg_energy = energy / n_frames as f64;
        println!("{:8} {:17.4}", idx + 1, avg_energy);
    }
}
