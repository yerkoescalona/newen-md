//! r_real_factor - Real-space R-factor for X-ray refinement
//!
//! Usage: r_real_factor @obs_map <file> @calc_map <file>
//!
//! Calculates real-space R-factor between observed and calculated maps

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: r_real_factor @obs_map <file> @calc_map <file>");
        process::exit(1);
    }

    let mut obs_file = None;
    let mut calc_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@obs_map" => {
                i += 1;
                obs_file = Some(args[i].clone());
            },
            "@calc_map" => {
                i += 1;
                calc_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let mut rho_obs = Vec::new();
    if let Some(path) = obs_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                let parts: Vec<&str> = line_str.split_whitespace().collect();
                if parts.len() >= 4 {
                    if let Ok(density) = parts[3].parse::<f64>() {
                        rho_obs.push(density);
                    }
                }
            }
        }
    }

    let mut rho_calc = Vec::new();
    if let Some(path) = calc_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                let parts: Vec<&str> = line_str.split_whitespace().collect();
                if parts.len() >= 4 {
                    if let Ok(density) = parts[3].parse::<f64>() {
                        rho_calc.push(density);
                    }
                }
            }
        }
    }

    if rho_obs.is_empty() || rho_calc.is_empty() {
        eprintln!("Error: Could not load density maps");
        process::exit(1);
    }

    let n = rho_obs.len().min(rho_calc.len());

    eprintln!("# Real-space R-factor calculation");
    eprintln!("# Grid points: {}", n);

    // Calculate real-space R-factor
    let mut sum_diff = 0.0f64;
    let mut sum_obs = 0.0f64;

    for i in 0..n {
        sum_diff += (rho_obs[i] - rho_calc[i]).abs();
        sum_obs += rho_obs[i].abs();
    }

    let r_real = if sum_obs > 0.0 {
        sum_diff / sum_obs
    } else {
        0.0
    };

    // Also calculate correlation coefficient
    let mean_obs: f64 = rho_obs.iter().take(n).sum::<f64>() / n as f64;
    let mean_calc: f64 = rho_calc.iter().take(n).sum::<f64>() / n as f64;

    let mut cov = 0.0f64;
    let mut var_obs = 0.0f64;
    let mut var_calc = 0.0f64;

    for i in 0..n {
        let diff_obs = rho_obs[i] - mean_obs;
        let diff_calc = rho_calc[i] - mean_calc;
        cov += diff_obs * diff_calc;
        var_obs += diff_obs * diff_obs;
        var_calc += diff_calc * diff_calc;
    }

    let correlation = if var_obs > 0.0 && var_calc > 0.0 {
        cov / (var_obs * var_calc).sqrt()
    } else {
        0.0
    };

    println!("# Real-Space R-factor Analysis");
    println!("# Grid points: {}", n);
    println!("# R_real = Σ|ρ_obs - ρ_calc| / Σ|ρ_obs|");
    println!("# R_real: {:.4}", r_real);
    println!("# R_real (%): {:.2}", r_real * 100.0);
    println!("# Correlation coefficient: {:.4}", correlation);
}
