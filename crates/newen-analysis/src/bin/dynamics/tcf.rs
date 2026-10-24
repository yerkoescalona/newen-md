//! tcf - Time correlation function analysis
//!
//! Usage: tcf @traj <file> @property <name>
//!
//! Calculates time correlation functions

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: tcf @traj <file> [@property <velocity|position>] [@maxlag <n>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut property = "velocity";
    let mut max_lag = 100;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@property" => {
                i += 1;
                property = &args[i];
            },
            "@maxlag" => {
                i += 1;
                max_lag = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!("# Time correlation function");
    eprintln!("# Property: {}", property);
    eprintln!("# Max lag: {}", max_lag);

    // Collect all frame data
    let mut frame_data = Vec::new();
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                // For simplicity, use center of mass position
                let mut com_x = 0.0f64;
                let mut com_y = 0.0f64;
                let mut com_z = 0.0f64;
                let n = frame.positions.len() as f64;

                for pos in &frame.positions {
                    com_x += pos.x;
                    com_y += pos.y;
                    com_z += pos.z;
                }

                frame_data.push((com_x / n, com_y / n, com_z / n));
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    let n_frames = frame_data.len();
    eprintln!("# Loaded {} frames", n_frames);

    if n_frames < 2 {
        eprintln!("Error: Need at least 2 frames");
        process::exit(1);
    }

    println!("# Lag (frames)    Correlation");

    // Calculate autocorrelation
    for lag in 0..max_lag.min(n_frames / 2) {
        let mut correlation = 0.0f64;
        let mut count = 0usize;

        for i in 0..(n_frames - lag) {
            let (x0, y0, z0) = frame_data[i];
            let (x1, y1, z1) = frame_data[i + lag];

            correlation += x0 * x1 + y0 * y1 + z0 * z1;
            count += 1;
        }

        let avg_corr = if count > 0 {
            correlation / count as f64
        } else {
            0.0
        };

        println!("{:14} {:16.8}", lag, avg_corr);
    }
}
