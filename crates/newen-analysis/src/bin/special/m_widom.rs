//! m_widom - Widom insertion for chemical potential
//!
//! Usage: m_widom @topo <file> @traj <file> @insertions <n>
//!
//! Calculates excess chemical potential using Widom insertion

use newen_core::units::{four_pi_eps_i, kB};
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn calc_insertion_energy(
    insert_pos: (f64, f64, f64),
    positions: &[(f64, f64, f64)],
    charges: &[f64],
) -> f64 {
    let mut energy = 0.0f64;
    let insert_charge = 1.0f64; // Test particle charge

    for (i, &(x, y, z)) in positions.iter().enumerate() {
        let dx = insert_pos.0 - x;
        let dy = insert_pos.1 - y;
        let dz = insert_pos.2 - z;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist > 0.01 {
            let q = charges.get(i).copied().unwrap_or(0.0);
            energy += four_pi_eps_i * insert_charge * q / dist;

            // Simple LJ
            let sigma = 0.3f64;
            let epsilon = 1.0f64;
            let sr = sigma / dist;
            let sr6 = sr.powi(6);
            let sr12 = sr6 * sr6;
            energy += 4.0 * epsilon * (sr12 - sr6);
        }
    }

    energy
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: m_widom @topo <file> @traj <file> [@insertions <n>] [@temp <T>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut n_insertions = 100;
    let mut temperature = 300.0;

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
            "@insertions" => {
                i += 1;
                n_insertions = args[i].parse().unwrap();
            },
            "@temp" => {
                i += 1;
                temperature = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Widom insertion analysis");
    eprintln!("# Insertions per frame: {}", n_insertions);
    eprintln!("# Temperature: {} K", temperature);

    let beta = 1.0 / (kB * temperature);
    let mut total_exp_u = 0.0f64;
    let mut n_frames = 0;

    println!("# Frame    <exp(-βU)>    μ_ex (kJ/mol)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let positions: Vec<(f64, f64, f64)> =
                    frame.positions.iter().map(|p| (p.x, p.y, p.z)).collect();

                let mut frame_exp_u = 0.0f64;

                for _ in 0..n_insertions {
                    // Random insertion point
                    let x = (n_frames as f64 * 0.123).fract() * frame.box_dims.x;
                    let y = (n_frames as f64 * 0.456).fract() * frame.box_dims.y;
                    let z = (n_frames as f64 * 0.789).fract() * frame.box_dims.z;

                    let u = calc_insertion_energy((x, y, z), &positions, &topo.charge);
                    frame_exp_u += (-beta * u).exp();
                }

                frame_exp_u /= n_insertions as f64;
                total_exp_u += frame_exp_u;
                n_frames += 1;

                let mu_ex = -(1.0 / beta) * frame_exp_u.ln();
                println!("{:8} {:14.6} {:14.4}", n_frames, frame_exp_u, mu_ex);
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    if n_frames > 0 {
        let avg_exp_u = total_exp_u / n_frames as f64;
        let mu_ex = -(1.0 / beta) * avg_exp_u.ln();
        eprintln!("\n# Average excess chemical potential: {:.4} kJ/mol", mu_ex);
    }
}
