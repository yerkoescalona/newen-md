//! predict_noe - Predict NOE intensities from structure
//!
//! Usage: predict_noe @topo <file> @traj <file>
//!
//! Predicts NOE cross-peaks from trajectory

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: predict_noe @topo <file> @traj <file> [@cutoff <r>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut cutoff = 0.5f64; // nm

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

    eprintln!("# NOE prediction");
    eprintln!("# Distance cutoff: {} nm", cutoff);

    println!("# Atom_i    Atom_j    <r^-6> (nm^-6)    <r> (nm)    Predicted_NOE");

    let cutoff_sq = cutoff * cutoff;

    // Collect distances
    let mut pair_distances: std::collections::HashMap<(usize, usize), Vec<f64>> =
        std::collections::HashMap::new();

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let n_atoms = frame.positions.len();

                for i in 0..n_atoms {
                    for j in (i + 1)..n_atoms.min(i + 20) {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist_sq = dx * dx + dy * dy + dz * dz;

                        if dist_sq <= cutoff_sq {
                            let dist = dist_sq.sqrt();
                            pair_distances.entry((i, j)).or_default().push(dist);
                        }
                    }
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    // Calculate NOE predictions
    for ((i, j), distances) in &pair_distances {
        if !distances.is_empty() {
            let avg_r: f64 = distances.iter().sum::<f64>() / distances.len() as f64;
            let avg_r_minus_6: f64 =
                distances.iter().map(|&r| 1.0 / (r).powi(6)).sum::<f64>() / distances.len() as f64;

            let noe_intensity = avg_r_minus_6.powf(-1.0 / 6.0);

            println!(
                "{:8} {:8} {:17.6e} {:11.6} {:16.6}",
                i + 1,
                j + 1,
                avg_r_minus_6,
                avg_r,
                noe_intensity
            );
        }
    }

    eprintln!("# Predicted {} NOE cross-peaks", pair_distances.len());
}
