//! solute_entropy - Calculate solute entropy
//!
//! Usage: solute_entropy @traj <file> @temp <T>
//!
//! Estimates configurational entropy of solute

use newen_core::units::kB;
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: solute_entropy @traj <file> [@temp <T>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut temperature = 300.0;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@temp" => {
                i += 1;
                temperature = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Solute entropy calculation");
    eprintln!("# Temperature: {} K", temperature);

    // Collect positions for quasi-harmonic analysis
    let mut all_positions = Vec::new();
    let mut n_atoms = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                n_atoms = frame.positions.len();
                all_positions.push(frame.positions.clone());
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    let n_frames = all_positions.len();
    eprintln!("# Frames: {}", n_frames);
    eprintln!("# Atoms: {}", n_atoms);

    if n_frames < 2 {
        eprintln!("Error: Need at least 2 frames");
        process::exit(1);
    }

    // Calculate positional fluctuations (simplified)
    let mut fluctuations = vec![0.0f64; n_atoms];

    for atom_idx in 0..n_atoms {
        let mut mean_x = 0.0f64;
        let mut mean_y = 0.0f64;
        let mut mean_z = 0.0f64;

        for frame_pos in &all_positions {
            mean_x += frame_pos[atom_idx].x;
            mean_y += frame_pos[atom_idx].y;
            mean_z += frame_pos[atom_idx].z;
        }

        mean_x /= n_frames as f64;
        mean_y /= n_frames as f64;
        mean_z /= n_frames as f64;

        let mut variance = 0.0f64;
        for frame_pos in &all_positions {
            let dx = frame_pos[atom_idx].x - mean_x;
            let dy = frame_pos[atom_idx].y - mean_y;
            let dz = frame_pos[atom_idx].z - mean_z;
            variance += dx * dx + dy * dy + dz * dz;
        }

        fluctuations[atom_idx] = variance / n_frames as f64;
    }

    // Estimate entropy from fluctuations (quasi-harmonic approximation)
    let mut total_entropy = 0.0f64;

    for &fluct in &fluctuations {
        if fluct > 1e-10 {
            // S = k_B * ln(fluct) + constant
            let s_atom = kB * (fluct * 2.0 * std::f64::consts::PI * std::f64::consts::E).ln();
            total_entropy += s_atom;
        }
    }

    println!("# Configurational Entropy Estimation");
    println!("# Temperature: {} K", temperature);
    println!("# Atoms: {}", n_atoms);
    println!("# Total entropy: {:.4} kJ/(mol·K)", total_entropy);
    println!(
        "# Entropy per atom: {:.4} kJ/(mol·K)",
        total_entropy / n_atoms as f64
    );
}
