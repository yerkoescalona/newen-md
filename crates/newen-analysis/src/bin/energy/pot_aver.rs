//! pot_aver - Average potential energy surfaces
//!
//! Usage: pot_aver @ene <file>
//!
//! Averages potential energy components

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: pot_aver @ene <file>");
        process::exit(1);
    }

    let mut ene_file = None;

    let mut i = 1;
    while i < args.len() {
        if args[i].as_str() == "@ene" {
            i += 1;
            ene_file = Some(args[i].clone());
        }
        i += 1;
    }

    let mut energies = Vec::new();

    if let Some(path) = ene_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                if let Ok(val) = line_str.trim().parse::<f64>() {
                    energies.push(val);
                }
            }
        }
    }

    if energies.is_empty() {
        eprintln!("Error: No energy data loaded");
        process::exit(1);
    }

    eprintln!("# Potential energy averaging");
    eprintln!("# Data points: {}", energies.len());

    let mean = energies.iter().sum::<f64>() / energies.len() as f64;
    let variance =
        energies.iter().map(|&e| (e - mean).powi(2)).sum::<f64>() / energies.len() as f64;
    let std_dev = variance.sqrt();
    let min = energies.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    println!("# Potential Energy Statistics");
    println!("# Points: {}", energies.len());
    println!("# Mean: {:.6} kJ/mol", mean);
    println!("# Std Dev: {:.6} kJ/mol", std_dev);
    println!("# Min: {:.6} kJ/mol", min);
    println!("# Max: {:.6} kJ/mol", max);
    println!("# Range: {:.6} kJ/mol", max - min);
}
