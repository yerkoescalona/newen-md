//! pairlist - Generate pairlist for nonbonded interactions
//!
//! Usage: pairlist @topo <file> @traj <file> @cutoff <d>
//!
//! Generates list of atom pairs within cutoff

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: pairlist @topo <file> @traj <file> @cutoff <dist>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut cutoff = 1.4;

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

    eprintln!("# Generating pairlist with cutoff {} nm", cutoff);

    if let Ok(Some(frame)) = traj.read_frame() {
        let mut pair_count = 0;

        println!("# Frame at time {} ps", frame.time);

        for i in 0..frame.positions.len() {
            for j in (i + 1)..frame.positions.len() {
                let dx = frame.positions[i].x - frame.positions[j].x;
                let dy = frame.positions[i].y - frame.positions[j].y;
                let dz = frame.positions[i].z - frame.positions[j].z;
                let dist_sq = dx * dx + dy * dy + dz * dz;

                if dist_sq <= cutoff_sq {
                    println!("{:6} {:6}", i + 1, j + 1);
                    pair_count += 1;
                }
            }
        }

        eprintln!("  {} pairs within cutoff", pair_count);
    }
}
