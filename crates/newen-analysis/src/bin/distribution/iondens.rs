//! iondens - Ion density distribution
//!
//! Usage: iondens @topo <file> @traj <file> @ion <atom> @bins <n>
//!
//! Calculates radial ion density distribution

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: iondens @topo <file> @traj <file> @ion <atom> [@bins <n>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut ion_idx = 0;
    let mut n_bins = 50;

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
            "@ion" => {
                i += 1;
                ion_idx = args[i].parse::<usize>().unwrap() - 1;
            },
            "@bins" => {
                i += 1;
                n_bins = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let num_atoms = topo.mass.len();
    let max_dist = 2.0f64; // nm
    let bin_width = max_dist / n_bins as f64;

    let mut histogram = vec![0usize; n_bins];
    let mut n_frames = 0;

    eprintln!("# Ion density distribution");
    eprintln!("# Ion atom: {}", ion_idx + 1);
    eprintln!("# Bins: {}", n_bins);
    eprintln!("# Max distance: {} nm", max_dist);

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let ion_pos = frame.positions[ion_idx];

                for j in 0..num_atoms {
                    if j != ion_idx {
                        let dx = frame.positions[j].x - ion_pos.x;
                        let dy = frame.positions[j].y - ion_pos.y;
                        let dz = frame.positions[j].z - ion_pos.z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                        let bin = (dist / bin_width) as usize;
                        if bin < n_bins {
                            histogram[bin] += 1;
                        }
                    }
                }

                n_frames += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Analyzed {} frames", n_frames);

    println!("# Distance (nm)    Density");

    for (idx, &count) in histogram.iter().enumerate() {
        let r = (idx as f64 + 0.5) * bin_width;
        let volume = 4.0 / 3.0
            * std::f64::consts::PI
            * ((r + bin_width / 2.0).powi(3) - (r - bin_width / 2.0).powi(3));
        let density = if volume > 0.0 && n_frames > 0 {
            count as f64 / (volume * n_frames as f64)
        } else {
            0.0
        };

        println!("{:14.6} {:12.6}", r, density);
    }
}
