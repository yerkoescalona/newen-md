//! cluster - Trajectory clustering based on RMSD
//!
//! Usage: cluster @topo <file> @traj <file> @cutoff <rmsd>
//!
//! Performs RMSD-based clustering of trajectory frames

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

    if args.len() < 7 {
        eprintln!("Usage: cluster @topo <file> @traj <file> @cutoff <rmsd>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut cutoff = 0.1f64; // Default 0.1 nm

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

    eprintln!("# Clustering with RMSD cutoff: {} nm", cutoff);

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

    // Simple single-linkage clustering
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut assigned = vec![false; n_frames];

    for i in 0..n_frames {
        if !assigned[i] {
            let mut cluster = vec![i];
            assigned[i] = true;

            for j in (i + 1)..n_frames {
                if !assigned[j] {
                    let rmsd = calc_rmsd(&frames[i], &frames[j]);
                    if rmsd <= cutoff {
                        cluster.push(j);
                        assigned[j] = true;
                    }
                }
            }

            clusters.push(cluster);
        }
    }

    eprintln!("# Found {} clusters", clusters.len());

    println!("# Cluster    Size    Representative    Members");
    for (idx, cluster) in clusters.iter().enumerate() {
        print!("{:8} {:8} {:16}", idx + 1, cluster.len(), cluster[0] + 1);
        for &member in cluster.iter().take(10) {
            print!(" {}", member + 1);
        }
        if cluster.len() > 10 {
            print!(" ...");
        }
        println!();
    }
}
