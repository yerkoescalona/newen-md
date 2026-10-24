//! eps_field - Electric field analysis
//!
//! Usage: eps_field @topo <file> @traj <file> @atom <idx>
//!
//! Calculates electric field at specified atom

use newen_core::units::four_pi_eps_i;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: eps_field @topo <file> @traj <file> @atom <idx>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut target_atom = 0;

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
            "@atom" => {
                i += 1;
                target_atom = args[i].parse::<usize>().unwrap() - 1;
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    let num_atoms = topo.mass.len();

    if target_atom >= num_atoms {
        eprintln!(
            "Error: Atom {} out of range (max: {})",
            target_atom + 1,
            num_atoms
        );
        process::exit(1);
    }

    eprintln!("# Electric field analysis");
    eprintln!("# Target atom: {}", target_atom + 1);

    println!("# Time (ps)    E_x (V/nm)    E_y (V/nm)    E_z (V/nm)    |E| (V/nm)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let target_pos = frame.positions[target_atom];
                let mut ex = 0.0f64;
                let mut ey = 0.0f64;
                let mut ez = 0.0f64;

                for j in 0..num_atoms {
                    if j != target_atom {
                        let dx = target_pos.x - frame.positions[j].x;
                        let dy = target_pos.y - frame.positions[j].y;
                        let dz = target_pos.z - frame.positions[j].z;
                        let r_sq = dx * dx + dy * dy + dz * dz;
                        let r = r_sq.sqrt();

                        if r > 0.01 {
                            let q_j = topo.charge[j];
                            let factor = four_pi_eps_i * q_j / (r * r * r);

                            ex += factor * dx;
                            ey += factor * dy;
                            ez += factor * dz;
                        }
                    }
                }

                let e_mag = (ex * ex + ey * ey + ez * ez).sqrt();

                println!(
                    "{:12.4} {:13.6} {:13.6} {:13.6} {:13.6}",
                    frame.time, ex, ey, ez, e_mag
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
