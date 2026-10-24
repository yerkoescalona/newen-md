//! postcluster - Post-process clustering results
//!
//! Usage: postcluster @clusters <file> @traj <file>
//!
//! Analyzes and post-processes clustering output

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: postcluster @clusters <file> @traj <file>");
        process::exit(1);
    }

    let mut cluster_file = None;
    let mut traj_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@clusters" => {
                i += 1;
                cluster_file = Some(args[i].clone());
            },
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    // Read cluster assignments
    let mut assignments = Vec::new();
    if let Some(path) = cluster_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                if let Ok(cluster_id) = line_str.trim().parse::<usize>() {
                    assignments.push(cluster_id);
                }
            }
        }
    }

    if assignments.is_empty() {
        eprintln!("Warning: No cluster assignments loaded");
        process::exit(1);
    }

    let n_clusters = *assignments.iter().max().unwrap_or(&0) + 1;

    eprintln!("# Post-cluster analysis");
    eprintln!("# Frames: {}", assignments.len());
    eprintln!("# Clusters: {}", n_clusters);

    // Count cluster populations
    let mut populations = vec![0usize; n_clusters];
    for &cluster_id in &assignments {
        if cluster_id < n_clusters {
            populations[cluster_id] += 1;
        }
    }

    println!("# Cluster    Population    Percentage");
    for (cluster_id, &pop) in populations.iter().enumerate() {
        let percentage = 100.0 * pop as f64 / assignments.len() as f64;
        println!("{:8} {:13} {:12.2}", cluster_id + 1, pop, percentage);
    }

    // Analyze trajectory if provided
    if let Some(traj_path) = traj_file {
        if let Ok(mut traj) = TrajectoryReader::new(&traj_path) {
            eprintln!("\n# Analyzing cluster properties from trajectory");

            let mut cluster_times: Vec<Vec<f64>> = vec![Vec::new(); n_clusters];
            let mut frame_idx = 0;

            loop {
                match traj.read_frame() {
                    Ok(Some(frame)) => {
                        if frame_idx < assignments.len() {
                            let cluster_id = assignments[frame_idx];
                            if cluster_id < n_clusters {
                                cluster_times[cluster_id].push(frame.time);
                            }
                        }
                        frame_idx += 1;
                    },
                    Ok(None) => break,
                    Err(_) => break,
                }
            }

            println!("\n# Cluster    First_Time    Last_Time     Lifetime");
            for (cluster_id, times) in cluster_times.iter().enumerate() {
                if !times.is_empty() {
                    let first = times.first().unwrap();
                    let last = times.last().unwrap();
                    let lifetime = last - first;
                    println!(
                        "{:8} {:13.2} {:13.2} {:13.2}",
                        cluster_id + 1,
                        first,
                        last,
                        lifetime
                    );
                }
            }
        }
    }
}
