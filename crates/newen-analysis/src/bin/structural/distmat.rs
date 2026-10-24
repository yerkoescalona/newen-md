//! distmat - Distance matrix calculation
//!
//! Usage: distmat @traj <file> @atoms <start> <end>
//!
//! Calculates pairwise distance matrix

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: distmat @traj <file> [@atoms <start> <end>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut atom_start = 0;
    let mut atom_end = 10;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@atoms" => {
                i += 1;
                atom_start = args[i].parse::<usize>().unwrap() - 1;
                i += 1;
                atom_end = args[i].parse::<usize>().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Distance matrix");
    eprintln!("# Atoms: {} to {}", atom_start + 1, atom_end);

    let n_atoms = atom_end - atom_start;
    let mut dist_sums = vec![vec![0.0f64; n_atoms]; n_atoms];
    let mut n_frames = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                for i in 0..n_atoms {
                    let idx_i = atom_start + i;
                    if idx_i >= frame.positions.len() {
                        break;
                    }

                    for j in 0..n_atoms {
                        let idx_j = atom_start + j;
                        if idx_j >= frame.positions.len() {
                            break;
                        }

                        let dx = frame.positions[idx_i].x - frame.positions[idx_j].x;
                        let dy = frame.positions[idx_i].y - frame.positions[idx_j].y;
                        let dz = frame.positions[idx_i].z - frame.positions[idx_j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                        dist_sums[i][j] += dist;
                    }
                }

                n_frames += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Frames: {}", n_frames);

    println!("# Average distance matrix (nm)");
    println!("# Row: atom_i, Column: atom_j");

    for i in 0..n_atoms {
        for j in 0..n_atoms {
            let avg_dist = dist_sums[i][j] / n_frames as f64;
            print!("{:8.4} ", avg_dist);
        }
        println!();
    }
}
