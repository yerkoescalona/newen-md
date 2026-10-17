//! gathtraj - Gather trajectory with PBC
//!
//! Usage: gathtraj @topo <file> @traj <file> @pbc <type>
//!
//! Applies periodic boundary conditions to trajectory

use newen_core::math::Vec3;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn apply_pbc(mut pos: Vec3, box_dims: Vec3) -> Vec3 {
    if box_dims.x > 0.0 {
        pos.x -= (pos.x / box_dims.x).floor() * box_dims.x;
    }
    if box_dims.y > 0.0 {
        pos.y -= (pos.y / box_dims.y).floor() * box_dims.y;
    }
    if box_dims.z > 0.0 {
        pos.z -= (pos.z / box_dims.z).floor() * box_dims.z;
    }
    pos
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: gathtraj @topo <file> @traj <file> @pbc <type>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut pbc_type = String::new();

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
            "@pbc" => {
                i += 1;
                pbc_type = args[i].clone();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    println!("TITLE");
    println!("Gathered trajectory with PBC");
    println!("END");

    loop {
        match traj.read_frame() {
            Ok(Some(mut frame)) => {
                if pbc_type != "v" {
                    let box_dims = frame.box_dims;
                    for pos in &mut frame.positions {
                        *pos = apply_pbc(*pos, box_dims);
                    }
                }

                println!("TIMESTEP");
                println!("{:15} {:15.4}", frame.step, frame.time);
                println!("END");
                println!("POSITIONRED");

                for (i, pos) in frame.positions.iter().enumerate() {
                    println!("{:6} {:15.9} {:15.9} {:15.9}", i + 1, pos.x, pos.y, pos.z);
                }

                println!("END");
                println!("BOX");
                println!(
                    " {:15.9} {:15.9} {:15.9}",
                    frame.box_dims.x, frame.box_dims.y, frame.box_dims.z
                );
                println!("END");
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
