//! ditrans - Distance in transition analysis
//!
//! Usage: ditrans @topo <file> @traj <file> @ref <i> @target <j> @cutoff <d>
//!
//! Tracks when atoms transition across distance threshold

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 11 {
        eprintln!("Usage: ditrans @topo <file> @traj <file> @ref <i> @target <j> @cutoff <d>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut ref_atom = 0;
    let mut target_atom = 0;
    let mut cutoff = 0.5;

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
            "@ref" => {
                i += 1;
                ref_atom = args[i].parse::<usize>().unwrap() - 1;
            },
            "@target" => {
                i += 1;
                target_atom = args[i].parse::<usize>().unwrap() - 1;
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

    let mut was_inside = false;
    let cutoff_sq = cutoff * cutoff;

    println!("# Time (ps)    State    Distance (nm)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                if ref_atom < frame.positions.len() && target_atom < frame.positions.len() {
                    let dx = frame.positions[ref_atom].x - frame.positions[target_atom].x;
                    let dy = frame.positions[ref_atom].y - frame.positions[target_atom].y;
                    let dz = frame.positions[ref_atom].z - frame.positions[target_atom].z;
                    let dist_sq = dx * dx + dy * dy + dz * dz;
                    let dist = dist_sq.sqrt();

                    let is_inside = dist_sq <= cutoff_sq;
                    let state = if is_inside { "IN" } else { "OUT" };

                    if is_inside != was_inside {
                        println!(
                            "{:12.4} {:>6} {:12.6}  # TRANSITION",
                            frame.time, state, dist
                        );
                        was_inside = is_inside;
                    } else {
                        println!("{:12.4} {:>6} {:12.6}", frame.time, state, dist);
                    }
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
