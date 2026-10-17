//! tstrip - Strip solvent from trajectory
//!
//! Usage: tstrip @topo <file> @traj <file> [@natoms <n>] [@stride <n>]
//!
//! Removes solvent molecules, keeping only first N atoms (solute)

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn print_usage() {
    eprintln!("tstrip - Strip solvent from trajectory");
    eprintln!();
    eprintln!("Usage: tstrip @topo <file> @traj <file> [@natoms <n>] [@stride <n>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo     Molecular topology file");
    eprintln!("  @traj     Input trajectory file");
    eprintln!("  @natoms   Number of solute atoms to keep (default: from topology)");
    eprintln!("  @stride   Write every nth frame (default: 1)");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  tstrip @topo system.top @traj output.trc");
    eprintln!("  tstrip @topo system.top @traj output.trc @natoms 1000 @stride 10");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        print_usage();
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut natoms = None;
    let mut stride = 1;

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
            "@natoms" => {
                i += 1;
                natoms = Some(args[i].parse().unwrap());
            },
            "@stride" => {
                i += 1;
                stride = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    // Read topology
    let topo_data = read_topology_file(topo_file.unwrap()).unwrap_or_else(|e| {
        eprintln!("Error reading topology: {}", e);
        process::exit(1);
    });

    let topo = build_topology(topo_data);
    let num_solute = natoms.unwrap_or(topo.num_solute_atoms());

    // Open trajectory
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap_or_else(|e| {
        eprintln!("Error opening trajectory: {}", e);
        process::exit(1);
    });

    eprintln!("Stripping to {} solute atoms", num_solute);
    eprintln!("Stride: {}", stride);

    println!("TITLE");
    println!("Stripped trajectory");
    println!("END");

    let mut frame_count = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                if frame_count % stride == 0 {
                    println!("TIMESTEP");
                    println!("{:15} {:15.4}", frame.step, frame.time);
                    println!("END");
                    println!("POSITIONRED");

                    for (i, pos) in frame.positions.iter().take(num_solute).enumerate() {
                        println!("{:6} {:15.9} {:15.9} {:15.9}", i + 1, pos.x, pos.y, pos.z);
                    }

                    println!("END");
                    println!("BOX");
                    println!(
                        " {:15.9} {:15.9} {:15.9}",
                        frame.box_dims.x, frame.box_dims.y, frame.box_dims.z
                    );
                    println!("END");
                }
                frame_count += 1;
            },
            Ok(None) => break,
            Err(e) => {
                eprintln!("Error reading frame: {}", e);
                process::exit(1);
            },
        }
    }

    eprintln!("Processed {} frames", frame_count);
}
