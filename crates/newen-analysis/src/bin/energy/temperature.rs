//! temperature - Calculate temperature from velocities
//!
//! Usage: temperature @topo <file> @traj <file> @dof <n>
//!
//! Calculates instantaneous temperature from velocity trajectory

use newen_core::units::kB;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: temperature @topo <file> @traj <file> @dof <n>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut dof = 0;

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
            "@dof" => {
                i += 1;
                dof = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    println!("# Time (ps)    Temperature (K)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut ke = 0.0f64;

                // Calculate kinetic energy from velocities
                if let Some(ref velocities) = frame.velocities {
                    for (i, vel) in velocities.iter().enumerate() {
                        if i < topo.mass.len() {
                            let v_sq = vel.x * vel.x + vel.y * vel.y + vel.z * vel.z;
                            ke += 0.5 * topo.mass[i] * v_sq;
                        }
                    }
                } else {
                    // Use positions as velocities (for position trajectories)
                    for (i, pos) in frame.positions.iter().enumerate() {
                        if i < topo.mass.len() {
                            let v_sq = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
                            ke += 0.5 * topo.mass[i] * v_sq;
                        }
                    }
                }

                let temp = 2.0 * ke / (kB * dof as f64);
                println!("{:12.4} {:12.4}", frame.time, temp);
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
