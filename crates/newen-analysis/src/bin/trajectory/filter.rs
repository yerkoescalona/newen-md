//! filter - Filter trajectory frames based on criteria
//!
//! Usage: filter @traj <file> @every <n> @start <frame> @end <frame>
//!
//! Selects frames from trajectory based on filters

use newen_io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: filter @traj <file> [@every <n>] [@start <frame>] [@end <frame>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut every = 1;
    let mut start_frame = 0;
    let mut end_frame = usize::MAX;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@every" => {
                i += 1;
                every = args[i].parse().unwrap();
            },
            "@start" => {
                i += 1;
                start_frame = args[i].parse::<usize>().unwrap() - 1; // Convert to 0-indexed
            },
            "@end" => {
                i += 1;
                end_frame = args[i].parse::<usize>().unwrap() - 1; // Convert to 0-indexed
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(traj_file.unwrap()).unwrap();

    eprintln!(
        "# Filtering frames: start={}, end={}, every={}",
        start_frame + 1,
        if end_frame == usize::MAX {
            "end".to_string()
        } else {
            (end_frame + 1).to_string()
        },
        every
    );

    let mut frame_idx = 0;
    let mut output_count = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                if frame_idx >= start_frame
                    && frame_idx <= end_frame
                    && (frame_idx - start_frame) % every == 0
                {
                    println!("TIMESTEP");
                    println!("  {}", frame_idx + 1);
                    println!("  {}", frame.time);
                    println!("POSITION");
                    println!(
                        "  {:15.9} {:15.9} {:15.9}",
                        frame.box_dims.x, frame.box_dims.y, frame.box_dims.z
                    );

                    for pos in &frame.positions {
                        println!("{:15.9} {:15.9} {:15.9}", pos.x, pos.y, pos.z);
                    }
                    println!("END");
                    output_count += 1;
                }

                frame_idx += 1;
                if frame_idx > end_frame {
                    break;
                }
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!(
        "# Processed {} frames, output {} frames",
        frame_idx, output_count
    );
}
