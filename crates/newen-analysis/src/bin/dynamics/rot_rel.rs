//! rot_rel - Rotational relaxation time calculation
//!
//! Usage: rot_rel @traj <file> @vector <i j>
//!
//! Calculates rotational correlation time

use newen_core::math::Vec3;
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn normalize(v: Vec3) -> Vec3 {
    let len = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
    if len > 0.0 {
        Vec3::new(v.x / len, v.y / len, v.z / len)
    } else {
        v
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: rot_rel @traj <file> @vector <atom_i> <atom_j>");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut atom_i = 0;
    let mut atom_j = 1;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@vector" => {
                i += 1;
                atom_i = args[i].parse::<usize>().unwrap() - 1;
                i += 1;
                atom_j = args[i].parse::<usize>().unwrap() - 1;
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Rotational relaxation analysis");
    eprintln!("# Vector: atom {} to atom {}", atom_i + 1, atom_j + 1);

    // Collect all vectors
    let mut vectors = Vec::new();

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                if atom_i < frame.positions.len() && atom_j < frame.positions.len() {
                    let v = Vec3::new(
                        frame.positions[atom_j].x - frame.positions[atom_i].x,
                        frame.positions[atom_j].y - frame.positions[atom_i].y,
                        frame.positions[atom_j].z - frame.positions[atom_i].z,
                    );
                    vectors.push(normalize(v));
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    let n_frames = vectors.len();
    eprintln!("# Frames: {}", n_frames);

    if n_frames < 2 {
        eprintln!("Error: Need at least 2 frames");
        process::exit(1);
    }

    println!("# Lag (frames)    C(t)    P2(t)");

    // Calculate correlation functions
    let max_lag = n_frames / 4;

    for lag in 0..max_lag {
        let mut c1 = 0.0f64;
        let mut c2 = 0.0f64;
        let mut count = 0;

        for t in 0..(n_frames - lag) {
            let v0 = vectors[t];
            let vt = vectors[t + lag];

            // C1(t) = <v(0)·v(t)>
            let dot = v0.x * vt.x + v0.y * vt.y + v0.z * vt.z;
            c1 += dot;

            // P2(t) = <(3*cos²θ - 1)/2>
            let p2 = 0.5 * (3.0 * dot * dot - 1.0);
            c2 += p2;

            count += 1;
        }

        if count > 0 {
            c1 /= count as f64;
            c2 /= count as f64;
        }

        println!("{:14} {:8.6} {:8.6}", lag, c1, c2);
    }
}
