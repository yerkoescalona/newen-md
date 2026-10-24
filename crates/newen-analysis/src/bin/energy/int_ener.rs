//! int_ener - Interaction energies between groups
//!
//! Usage: int_ener @topo <file> @traj <file> @group1 <start> <end> @group2 <start> <end>
//!
//! Calculates interaction energies between two atom groups

use newen_core::units::four_pi_eps_i;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 11 {
        eprintln!(
            "Usage: int_ener @topo <file> @traj <file> @group1 <start> <end> @group2 <start> <end>"
        );
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut g1_start = 0;
    let mut g1_end = 10;
    let mut g2_start = 10;
    let mut g2_end = 20;

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
            "@group1" => {
                i += 1;
                g1_start = args[i].parse::<usize>().unwrap() - 1;
                i += 1;
                g1_end = args[i].parse::<usize>().unwrap();
            },
            "@group2" => {
                i += 1;
                g2_start = args[i].parse::<usize>().unwrap() - 1;
                i += 1;
                g2_end = args[i].parse::<usize>().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Interaction energy analysis");
    eprintln!("# Group 1: atoms {} to {}", g1_start + 1, g1_end);
    eprintln!("# Group 2: atoms {} to {}", g2_start + 1, g2_end);

    println!("# Time (ps)    E_elec (kJ/mol)    E_vdw (kJ/mol)    E_tot (kJ/mol)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut e_elec = 0.0f64;
                let mut e_vdw = 0.0f64;

                for i in g1_start..g1_end.min(frame.positions.len()) {
                    for j in g2_start..g2_end.min(frame.positions.len()) {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                        if dist > 0.01 {
                            // Electrostatic
                            let qi = topo.charge[i];
                            let qj = topo.charge[j];
                            e_elec += four_pi_eps_i * qi * qj / dist;

                            // VDW (simple LJ)
                            let sigma = 0.3f64;
                            let epsilon = 1.0f64;
                            let sr = sigma / dist;
                            let sr6 = sr.powi(6);
                            let sr12 = sr6 * sr6;
                            e_vdw += 4.0 * epsilon * (sr12 - sr6);
                        }
                    }
                }

                let e_tot = e_elec + e_vdw;

                println!(
                    "{:12.4} {:18.4} {:17.4} {:17.4}",
                    frame.time, e_elec, e_vdw, e_tot
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
