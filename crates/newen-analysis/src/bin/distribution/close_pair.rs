//! close_pair - Find close atom pairs
//!
//! Usage: close_pair @topo <file> @traj <file> @cutoff <d>
//!
//! Reports atom pairs within cutoff distance

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: close_pair @topo <file> @traj <file> @cutoff <dist>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut cutoff = 0.5;

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
            "@cutoff" => {
                i += 1;
                cutoff = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let cutoff_sq = cutoff * cutoff;

    println!("# Time (ps)    Atom_i  Atom_j  Distance (nm)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                for i in 0..frame.positions.len() {
                    for j in (i + 1)..frame.positions.len() {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist_sq = dx * dx + dy * dy + dz * dz;

                        if dist_sq <= cutoff_sq {
                            let dist = dist_sq.sqrt();
                            println!("{:12.4} {:8} {:8} {:12.6}", frame.time, i + 1, j + 1, dist);
                        }
                    }
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
