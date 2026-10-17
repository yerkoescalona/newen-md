//! xrayts - X-ray time series analysis
//!
//! Usage: xrayts @obs <file> @calc <file>
//!
//! Analyzes time series of X-ray observables

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: xrayts @obs <file> [@calc <file>]");
        process::exit(1);
    }

    let mut obs_file = None;
    let mut calc_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@obs" => {
                i += 1;
                obs_file = Some(args[i].clone());
            },
            "@calc" => {
                i += 1;
                calc_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let mut obs_data = Vec::new();
    if let Some(path) = obs_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                if let Ok(val) = line_str.trim().parse::<f64>() {
                    obs_data.push(val);
                }
            }
        }
    }

    let mut calc_data = Vec::new();
    if let Some(path) = calc_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                if let Ok(val) = line_str.trim().parse::<f64>() {
                    calc_data.push(val);
                }
            }
        }
    }

    if obs_data.is_empty() {
        eprintln!("Error: No observed data");
        process::exit(1);
    }

    eprintln!("# X-ray time series analysis");
    eprintln!("# Observed points: {}", obs_data.len());
    eprintln!("# Calculated points: {}", calc_data.len());

    println!("# Point    Observed    Calculated    Deviation");

    let n = obs_data.len().min(if calc_data.is_empty() {
        obs_data.len()
    } else {
        calc_data.len()
    });

    let mut sum_dev_sq = 0.0f64;

    for i in 0..n {
        let obs = obs_data[i];
        let calc = if i < calc_data.len() {
            calc_data[i]
        } else {
            0.0
        };
        let dev = obs - calc;

        sum_dev_sq += dev * dev;

        println!("{:8} {:12.6} {:12.6} {:12.6}", i + 1, obs, calc, dev);
    }

    let rmsd = (sum_dev_sq / n as f64).sqrt();
    eprintln!("\n# RMSD: {:.6}", rmsd);
}
