//! xray_map - Generate X-ray electron density maps
//!
//! Usage: xray_map @topo <file> @traj <file>
//!
//! Generates electron density maps from trajectory

use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: xray_map @topo <file> @traj <file> [@grid <n>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut grid_size = 30;

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
            "@grid" => {
                i += 1;
                grid_size = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# X-ray electron density map");
    eprintln!("# Grid size: {}³", grid_size);

    let mut density_map = vec![vec![vec![0.0f64; grid_size]; grid_size]; grid_size];
    let mut n_frames = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let grid_spacing_x = frame.box_dims.x / grid_size as f64;
                let grid_spacing_y = frame.box_dims.y / grid_size as f64;
                let grid_spacing_z = frame.box_dims.z / grid_size as f64;

                for pos in &frame.positions {
                    let ix = ((pos.x / grid_spacing_x) as usize).min(grid_size - 1);
                    let iy = ((pos.y / grid_spacing_y) as usize).min(grid_size - 1);
                    let iz = ((pos.z / grid_spacing_z) as usize).min(grid_size - 1);

                    density_map[ix][iy][iz] += 1.0;
                }

                n_frames += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Frames: {}", n_frames);

    // Normalize and output
    println!("# X-ray density map");
    println!("# Grid: {} x {} x {}", grid_size, grid_size, grid_size);
    println!("# Frames: {}", n_frames);

    for ix in 0..grid_size {
        for iy in 0..grid_size {
            for iz in 0..grid_size {
                let density = density_map[ix][iy][iz] / n_frames as f64;
                if density > 0.1 {
                    println!("{:4} {:4} {:4} {:10.4}", ix, iy, iz, density);
                }
            }
        }
    }
}
