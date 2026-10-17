//! rdf_matrix - Radial distribution function matrix
//!
//! Usage: rdf_matrix @topo <file> @traj <file>
//!
//! Calculates RDF matrix between all atom types

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: rdf_matrix @topo <file> @traj <file> [@bins <n>] [@rmax <r>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut n_bins = 100;
    let mut r_max = 1.5f64; // nm

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
            "@bins" => {
                i += 1;
                n_bins = args[i].parse().unwrap();
            },
            "@rmax" => {
                i += 1;
                r_max = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let n_atoms = topo.mass.len();

    eprintln!("# RDF matrix calculation");
    eprintln!("# Atoms: {}", n_atoms);
    eprintln!("# Bins: {}", n_bins);
    eprintln!("# r_max: {} nm", r_max);

    let bin_width = r_max / n_bins as f64;
    let mut histograms = vec![vec![0usize; n_bins]; n_atoms];
    let mut n_frames = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                for i in 0..n_atoms.min(frame.positions.len()) {
                    for j in (i + 1)..n_atoms.min(frame.positions.len()) {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                        if dist < r_max {
                            let bin = (dist / bin_width) as usize;
                            if bin < n_bins {
                                histograms[i][bin] += 1;
                            }
                        }
                    }
                }

                n_frames += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Frames: {}", n_frames);

    println!("# r (nm)    RDF (averaged over all pairs)");

    for bin in 0..n_bins {
        let r = (bin as f64 + 0.5) * bin_width;
        let total_count: usize = histograms.iter().map(|h| h[bin]).sum();
        let avg_count = total_count as f64 / n_atoms as f64;

        println!("{:8.4} {:12.6}", r, avg_count);
    }
}
