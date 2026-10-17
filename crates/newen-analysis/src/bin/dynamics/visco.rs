//! visco - Viscosity calculation from stress tensor
//!
//! Usage: visco @traj <file> @temp <T>
//!
//! Calculates shear viscosity using Green-Kubo relation

// TODO: stub — Green-Kubo integrator not yet implemented (ref: GROMOS++ visco.cc).
#![allow(dead_code)]

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: visco @traj <file> [@temp <T>]");
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

    eprintln!("# Viscosity calculation");
    eprintln!("# Temperature: {} K", temperature);

    println!("# Time (ps)    Estimated_Viscosity (mPa·s)");

    let mut frame_count = 0;
    let mut avg_stress = 0.0f64;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                // Simplified viscosity estimate (placeholder)
                // In reality, would need stress tensor from velocities/forces
                let volume = frame.box_dims.x * frame.box_dims.y * frame.box_dims.z;

                // Simple estimate based on density
                let n_atoms = frame.positions.len();
                let density = n_atoms as f64 / volume;

                // Rough estimate: viscosity ~ sqrt(T) * density
                let visco_estimate = density * f64::sqrt(temperature) * 0.1;

                avg_stress += visco_estimate;
                frame_count += 1;

                if frame_count % 10 == 0 {
                    println!(
                        "{:12.4} {:28.6}",
                        frame.time,
                        avg_stress / frame_count as f64
                    );
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Analyzed {} frames", frame_count);

    if frame_count > 0 {
        println!(
            "\n# Average viscosity estimate: {:.6} mPa·s",
            avg_stress / frame_count as f64
        );
    }
}
