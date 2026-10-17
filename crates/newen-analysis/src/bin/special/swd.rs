//! swd - Sliding window distance analysis
//!
//! Usage: swd @topo <file> @traj <file> @atoms1 <i> @atoms2 <j>
//!
//! Calculates distance between two atoms over time

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 9 {
        eprintln!("Usage: swd @topo <file> @traj <file> @atoms1 <i> @atoms2 <j>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut atom1 = 0;
    let mut atom2 = 0;

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
            "@atoms1" => {
                i += 1;
                atom1 = args[i].parse::<usize>().unwrap() - 1;
            },
            "@atoms2" => {
                i += 1;
                atom2 = args[i].parse::<usize>().unwrap() - 1;
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    println!("# Time (ps)    Distance (nm)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                if atom1 < frame.positions.len() && atom2 < frame.positions.len() {
                    let dx = frame.positions[atom1].x - frame.positions[atom2].x;
                    let dy = frame.positions[atom1].y - frame.positions[atom2].y;
                    let dz = frame.positions[atom1].z - frame.positions[atom2].z;
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                    println!("{:12.4} {:12.6}", frame.time, dist);
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
