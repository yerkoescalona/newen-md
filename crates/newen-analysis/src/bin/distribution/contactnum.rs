//! contactnum - Calculate contact numbers
//!
//! Usage: contactnum @topo <file> @traj <file> @cutoff <distance>
//!
//! Counts number of contacts within cutoff distance

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: contactnum @topo <file> @traj <file> @cutoff <distance>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut cutoff = 0.6f64; // Default 6 Å = 0.6 nm

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
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let cutoff_sq = cutoff * cutoff;
    let num_atoms = topo.mass.len();

    eprintln!("# Contact number analysis");
    eprintln!("# Cutoff: {} nm", cutoff);
    eprintln!("# Atoms: {}", num_atoms);

    println!("# Time (ps)    Avg_Contacts    Total_Contacts");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut total_contacts = 0;

                for i in 0..num_atoms {
                    for j in (i + 1)..num_atoms {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist_sq = dx * dx + dy * dy + dz * dz;

                        if dist_sq <= cutoff_sq {
                            total_contacts += 1;
                        }
                    }
                }

                let avg_contacts = (2.0 * total_contacts as f64) / (num_atoms as f64);

                println!(
                    "{:12.4} {:15.4} {:18}",
                    frame.time, avg_contacts, total_contacts
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
