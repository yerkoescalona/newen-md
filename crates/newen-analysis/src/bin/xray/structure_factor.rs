//! structure_factor - Calculate structure factors
//!
//! Usage: structure_factor @topo <file> @traj <file>
//!
//! Calculates structure factors for X-ray scattering

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: structure_factor @topo <file> @traj <file> [@qmax <q>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut q_max = 2.0f64; // nm^-1

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
            "@qmax" => {
                i += 1;
                q_max = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Structure factor calculation");
    eprintln!("# q_max: {} nm^-1", q_max);

    let n_q = 100;
    let dq = q_max / n_q as f64;

    let mut s_q_sum = vec![0.0f64; n_q];
    let mut n_frames = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let n_atoms = frame.positions.len();

                for iq in 0..n_q {
                    let q = (iq as f64 + 0.5) * dq;
                    let mut sum_real = 0.0f64;
                    let mut sum_imag = 0.0f64;

                    for i in 0..n_atoms {
                        // Simplified: q along x-direction
                        let phase = q * frame.positions[i].x;
                        sum_real += phase.cos();
                        sum_imag += phase.sin();
                    }

                    let s_q = (sum_real * sum_real + sum_imag * sum_imag) / (n_atoms as f64);
                    s_q_sum[iq] += s_q;
                }

                n_frames += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Analyzed {} frames", n_frames);

    println!("# q (nm^-1)    S(q)");

    for iq in 0..n_q {
        let q = (iq as f64 + 0.5) * dq;
        let s_q_avg = s_q_sum[iq] / n_frames as f64;
        println!("{:12.6} {:12.6}", q, s_q_avg);
    }
}
