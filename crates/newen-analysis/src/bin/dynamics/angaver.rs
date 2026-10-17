//! angaver - Average angles from trajectory
//!
//! Usage: angaver @topo <file> @traj <file>
//!
//! Calculates average bond angles

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn calc_angle(pos1: (f64, f64, f64), pos2: (f64, f64, f64), pos3: (f64, f64, f64)) -> f64 {
    let v1x = pos1.0 - pos2.0;
    let v1y = pos1.1 - pos2.1;
    let v1z = pos1.2 - pos2.2;

    let v2x = pos3.0 - pos2.0;
    let v2y = pos3.1 - pos2.1;
    let v2z = pos3.2 - pos2.2;

    let dot = v1x * v2x + v1y * v2y + v1z * v2z;
    let len1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt();
    let len2 = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt();

    if len1 > 0.0 && len2 > 0.0 {
        let cos_theta = (dot / (len1 * len2)).clamp(-1.0, 1.0);
        cos_theta.acos().to_degrees()
    } else {
        0.0
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: angaver @topo <file> @traj <file>");
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
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let n_angles = topo.moltypes[0].angles.len();

    eprintln!("# Bond angle averaging");
    eprintln!("# Angles: {}", n_angles);

    let mut angle_sums = vec![0.0f64; n_angles];
    let mut n_frames = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                for (idx, angle) in topo.moltypes[0].angles.iter().enumerate() {
                    if angle.i < frame.positions.len()
                        && angle.j < frame.positions.len()
                        && angle.k < frame.positions.len()
                    {
                        let pos1 = (
                            frame.positions[angle.i].x,
                            frame.positions[angle.i].y,
                            frame.positions[angle.i].z,
                        );
                        let pos2 = (
                            frame.positions[angle.j].x,
                            frame.positions[angle.j].y,
                            frame.positions[angle.j].z,
                        );
                        let pos3 = (
                            frame.positions[angle.k].x,
                            frame.positions[angle.k].y,
                            frame.positions[angle.k].z,
                        );

                        let angle_val = calc_angle(pos1, pos2, pos3);
                        angle_sums[idx] += angle_val;
                    }
                }

                n_frames += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    if n_frames == 0 {
        eprintln!("Error: No frames read");
        process::exit(1);
    }

    eprintln!("# Frames: {}", n_frames);

    println!("# Angle    Atom_i    Atom_j    Atom_k    Avg_Angle (deg)");

    for (idx, angle) in topo.moltypes[0].angles.iter().enumerate() {
        let avg = angle_sums[idx] / n_frames as f64;
        println!(
            "{:8} {:8} {:8} {:8} {:17.4}",
            idx + 1,
            angle.i + 1,
            angle.j + 1,
            angle.k + 1,
            avg
        );
    }
}
