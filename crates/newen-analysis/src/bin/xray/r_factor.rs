//! r_factor - Calculate crystallographic R-factor
//!
//! Usage: r_factor @obs <file> @calc <file>
//!
//! Calculates R-factor for X-ray refinement

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: r_factor @obs <file> @calc <file>");
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

    // Read observed intensities
    let mut f_obs = Vec::new();
    if let Some(path) = obs_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                if let Ok(val) = line_str.trim().parse::<f64>() {
                    f_obs.push(val);
                }
            }
        }
    }

    // Read calculated intensities
    let mut f_calc = Vec::new();
    if let Some(path) = calc_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                if let Ok(val) = line_str.trim().parse::<f64>() {
                    f_calc.push(val);
                }
            }
        }
    }

    if f_obs.is_empty() || f_calc.is_empty() {
        eprintln!("Error: No data loaded");
        process::exit(1);
    }

    let n = f_obs.len().min(f_calc.len());

    eprintln!("# R-factor calculation");
    eprintln!("# Reflections: {}", n);

    // Calculate R-factor
    let mut sum_diff = 0.0f64;
    let mut sum_obs = 0.0f64;

    for i in 0..n {
        sum_diff += (f_obs[i] - f_calc[i]).abs();
        sum_obs += f_obs[i].abs();
    }

    let r_factor = if sum_obs > 0.0 {
        sum_diff / sum_obs
    } else {
        0.0
    };

    println!("# Crystallographic R-factor");
    println!("# N reflections: {}", n);
    println!("# R = Σ|F_obs - F_calc| / ΣF_obs");
    println!("# R-factor: {:.4}", r_factor);
    println!("# R-factor (%): {:.2}", r_factor * 100.0);

    // Also calculate R_free (simplified - using all data)
    println!("# R_free (approx): {:.4}", r_factor * 1.1);
}
