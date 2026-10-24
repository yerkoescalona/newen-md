//! follow - Follow molecules across periodic boundaries
//!
//! Usage: follow @topo <file> @traj <file>
//!
//! Tracks molecules continuously across PBC jumps

use newen_core::math::Vec3;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: follow @topo <file> @traj <file>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;

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
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let natoms = topo.mass.len();
    let mut prev_positions = vec![Vec3::ZERO; natoms];
    let mut shifts = vec![Vec3::ZERO; natoms];
    let mut first_frame = true;

    println!("TITLE");
    println!("Followed trajectory");
    println!("END");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let box_dims = frame.box_dims;

                if !first_frame {
                    // Detect and correct PBC jumps
                    for i in 0..natoms.min(frame.positions.len()) {
                        let pos = frame.positions[i];
                        let delta = pos - prev_positions[i];

                        // Check if jump occurred
                        if delta.x.abs() > box_dims.x * 0.5 {
                            shifts[i].x -= delta.x.signum() * box_dims.x;
                        }
                        if delta.y.abs() > box_dims.y * 0.5 {
                            shifts[i].y -= delta.y.signum() * box_dims.y;
                        }
                        if delta.z.abs() > box_dims.z * 0.5 {
                            shifts[i].z -= delta.z.signum() * box_dims.z;
                        }

                        prev_positions[i] = pos;
                    }
                }

                println!("TIMESTEP");
                println!("{:15} {:15.4}", frame.step, frame.time);
                println!("END");
                println!("POSITIONRED");

                for (i, pos) in frame.positions.iter().enumerate() {
                    if i < natoms {
                        let corrected = *pos + shifts[i];
                        println!(
                            "{:6} {:15.9} {:15.9} {:15.9}",
                            i + 1,
                            corrected.x,
                            corrected.y,
                            corrected.z
                        );
                        if first_frame {
                            prev_positions[i] = *pos;
                        }
                    }
                }

                println!("END");
                println!("BOX");
                println!(
                    " {:15.9} {:15.9} {:15.9}",
                    box_dims.x, box_dims.y, box_dims.z
                );
                println!("END");

                first_frame = false;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
