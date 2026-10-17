//! cog - Calculate center of geometry
//!
//! Usage: cog @topo <file> @traj <file> [@mass]
//!
//! Calculates center of geometry (or mass) for each frame

use newen_core::math::Vec3;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: cog @topo <file> @traj <file> [@mass]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut use_mass = false;

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
            "@mass" | "@com" => {
                use_mass = true;
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    println!("# Time (ps)    COG/COM_x (nm)  COG/COM_y (nm)  COG/COM_z (nm)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let center = if use_mass {
                    // Center of mass
                    let mut com = Vec3::ZERO;
                    let mut total_mass = 0.0f64;

                    for (i, pos) in frame.positions.iter().enumerate() {
                        if i < topo.mass.len() {
                            let mass = topo.mass[i];
                            com += *pos * mass;
                            total_mass += topo.mass[i];
                        }
                    }

                    if total_mass > 0.0 {
                        com /= total_mass;
                    }
                    com
                } else {
                    // Center of geometry
                    let mut cog = Vec3::ZERO;
                    for pos in &frame.positions {
                        cog += *pos;
                    }
                    cog /= frame.positions.len() as f64;
                    cog
                };

                println!(
                    "{:12.4} {:15.9} {:15.9} {:15.9}",
                    frame.time, center.x, center.y, center.z
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
