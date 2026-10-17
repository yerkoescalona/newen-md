//! rmsdmat - Calculate RMSD matrix between frames
//!
//! Usage: rmsdmat @topo <file> @traj <file>
//!
//! Computes pairwise RMSD between all trajectory frames

use newen_core::math::Vec3;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn calc_rmsd(pos1: &[Vec3], pos2: &[Vec3]) -> f64 {
    let mut sum_sq = 0.0f64;
    let n = pos1.len().min(pos2.len());

    for i in 0..n {
        let dx = pos1[i].x - pos2[i].x;
        let dy = pos1[i].y - pos2[i].y;
        let dz = pos1[i].z - pos2[i].z;
        sum_sq += dx * dx + dy * dy + dz * dz;
    }

    (sum_sq / n as f64).sqrt()
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: rmsdmat @topo <file> @traj <file>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;

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
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    // Read all frames
    let mut frames = Vec::new();
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                frames.push(frame.positions.clone());
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    let n_frames = frames.len();
    eprintln!("# Loaded {} frames", n_frames);

    println!("# RMSD Matrix ({} x {})", n_frames, n_frames);
    println!("# Row: Frame_i, Columns: RMSD to Frame_j");

    for i in 0..n_frames {
        print!("{:6}", i + 1);
        for j in 0..n_frames {
            let rmsd = calc_rmsd(&frames[i], &frames[j]);
            print!(" {:8.4}", rmsd);
        }
        println!();
    }
}
