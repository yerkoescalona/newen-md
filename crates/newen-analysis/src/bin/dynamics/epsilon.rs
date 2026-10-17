// TODO: stub — CONV_FACTOR/EPSILON_0 used by dielectric formula not yet implemented (ref: GROMOS++ epsilon.cc).
#![allow(dead_code)]
//! epsilon - Dielectric constant calculation
//!
//! Usage: epsilon @topo <file> @traj <file> @temp <T>
//!
//! Calculates static dielectric constant from dipole fluctuations

// TODO: stub — constants defined but the Kirkwood-Fröhlich formula that uses
// them is not yet implemented. Remove this once the dielectric calculation body
// is wired up (ref: GROMOS++ epsilon.cc, ~line 180).

use newen_core::math::Vec3;
use newen_core::units::kB;
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;
const EPSILON_0: f64 = 8.854187817e-12; // F/m (in SI units)
const CONV_FACTOR: f64 = 1389.35459; // Conversion factor for GROMOS units

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: epsilon @topo <file> @traj <file> [@temp <T>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
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

    eprintln!("# Dielectric constant calculation");
    eprintln!("# Temperature: {} K", temperature);

    let mut dipole_moments = Vec::new();

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut dipole = Vec3::ZERO;

                for (idx, pos) in frame.positions.iter().enumerate() {
                    if idx < topo.charge.len() {
                        let charge = topo.charge[idx];
                        dipole += *pos * charge;
                    }
                }

                let dipole_magnitude =
                    (dipole.x * dipole.x + dipole.y * dipole.y + dipole.z * dipole.z).sqrt();
                dipole_moments.push(dipole_magnitude);
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    if dipole_moments.is_empty() {
        eprintln!("Error: No frames read");
        process::exit(1);
    }

    let n_frames = dipole_moments.len();
    let avg_dipole: f64 = dipole_moments.iter().sum::<f64>() / n_frames as f64;
    let avg_dipole_sq: f64 = dipole_moments.iter().map(|&d| d * d).sum::<f64>() / n_frames as f64;

    let variance = avg_dipole_sq - avg_dipole * avg_dipole;

    eprintln!("# Analyzed {} frames", n_frames);
    eprintln!("# Average dipole: {:.4} e·nm", avg_dipole);

    println!("# Temperature (K)    <M> (e·nm)    <M²> (e²·nm²)    Variance    Epsilon");
    println!(
        "{:16.2} {:14.6} {:17.6} {:12.6} {:12.4}",
        temperature,
        avg_dipole,
        avg_dipole_sq,
        variance,
        1.0 + variance / (3.0 * kB * temperature)
    );
}
