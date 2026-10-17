//! cry - Crystallography analysis
//!
//! Usage: cry @traj <file> @cell <a b c alpha beta gamma>
//!
//! Analyzes crystal structures from trajectories

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: cry @traj <file> [@cell <a> <b> <c> <alpha> <beta> <gamma>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut cell_a = 1.0f64;
    let mut cell_b = 1.0f64;
    let mut cell_c = 1.0f64;
    let mut alpha = 90.0f64;
    let mut beta = 90.0f64;
    let mut gamma = 90.0f64;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@cell" => {
                i += 1;
                cell_a = args[i].parse().unwrap();
                i += 1;
                cell_b = args[i].parse().unwrap();
                i += 1;
                cell_c = args[i].parse().unwrap();
                i += 1;
                alpha = args[i].parse().unwrap();
                i += 1;
                beta = args[i].parse().unwrap();
                i += 1;
                gamma = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Crystallography analysis");
    eprintln!("# Unit cell: a={} b={} c={} nm", cell_a, cell_b, cell_c);
    eprintln!("# Angles: α={} β={} γ={} degrees", alpha, beta, gamma);

    // Calculate cell volume
    let alpha_rad = alpha.to_radians();
    let beta_rad = beta.to_radians();
    let gamma_rad = gamma.to_radians();

    let volume = cell_a
        * cell_b
        * cell_c
        * (1.0 - alpha_rad.cos().powi(2) - beta_rad.cos().powi(2) - gamma_rad.cos().powi(2)
            + 2.0 * alpha_rad.cos() * beta_rad.cos() * gamma_rad.cos())
        .sqrt();

    println!("# Unit Cell Parameters");
    println!("# a = {:.6} nm", cell_a);
    println!("# b = {:.6} nm", cell_b);
    println!("# c = {:.6} nm", cell_c);
    println!("# α = {:.2}°", alpha);
    println!("# β = {:.2}°", beta);
    println!("# γ = {:.2}°", gamma);
    println!("# Volume = {:.6} nm³", volume);

    println!("\n# Frame    Atoms    Density (g/cm³)");

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let n_atoms = frame.positions.len();
                let box_volume = frame.box_dims.x * frame.box_dims.y * frame.box_dims.z;

                // Estimate density (assuming average atomic mass ~12 g/mol)
                let avg_mass = 12.0; // g/mol
                let na = 6.022e23; // Avogadro's number
                let total_mass_g = (n_atoms as f64 * avg_mass) / na;
                let volume_cm3 = box_volume * 1e-21; // nm³ to cm³
                let density = total_mass_g / volume_cm3;

                println!("{:8} {:8} {:16.4}", frame_idx + 1, n_atoms, density);

                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
