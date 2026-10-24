//! noe - NMR NOE distance analysis
//!
//! Usage: noe @topo <file> @traj <file> @pairs <file>
//!
//! Analyzes distances for NOE restraints

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: noe @topo <file> @traj <file> @pairs <file>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut pairs_file = None;

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
            "@pairs" => {
                i += 1;
                pairs_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    // Read atom pairs from file
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    if let Some(pairs_path) = pairs_file {
        if let Ok(file) = File::open(&pairs_path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                let parts: Vec<&str> = line_str.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let (Ok(i), Ok(j)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) {
                        pairs.push((i - 1, j - 1)); // Convert to 0-indexed
                    }
                }
            }
        }
    }

    if pairs.is_empty() {
        eprintln!("# Warning: No pairs loaded, using default pair (1, 10)");
        pairs.push((0, 9));
    }

    eprintln!("# NOE distance analysis");
    eprintln!("# Number of pairs: {}", pairs.len());

    println!("# Time (ps)    Pair    Atom_i    Atom_j    Distance (nm)    r^-6_avg");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                for (idx, &(i, j)) in pairs.iter().enumerate() {
                    if i < frame.positions.len() && j < frame.positions.len() {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                        let r_minus_6 = if dist > 0.01 { 1.0 / dist.powi(6) } else { 0.0 };

                        println!(
                            "{:12.4} {:8} {:8} {:8} {:15.6} {:15.6}",
                            frame.time,
                            idx + 1,
                            i + 1,
                            j + 1,
                            dist,
                            r_minus_6
                        );
                    }
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
