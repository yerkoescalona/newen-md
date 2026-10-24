//! prep_noe - Prepare NOE restraint files
//!
//! Usage: prep_noe @pairs <file> @force <k>
//!
//! Prepares NOE restraint topology

use crate::programs::io_err;
use newen_io::args::{Arguments, Program, Status};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

pub const USAGE: &str = "prep_noe - Write NOE distance restraints\n\nUsage: prep_noe @pairs <file> [@dist <nm>] [@force <kJ/mol/nm^2>]\n\nArguments:\n  @pairs  Atom pairs, one \"i j\" per line\n  @dist   Default restraint distance (nm)\n  @force  Force constant\n\nNote: gromos++'s prep_noe reads an NOE specification file; this one takes plain pairs.";

pub const KNOWN: &[&str] = &["pairs", "dist", "force"];

pub const PROGRAM: Program = Program {
    name: "prep_noe",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let pairs_file = args.value("pairs")?.to_string();
    let force_const: f64 = args.get("force", 1000.0)?; // kJ/(mol.nm^2)
    let target_dist: f64 = args.get("dist", 0.3)?; // nm

    let mut pairs = Vec::new();
    {
        let path = &pairs_file;
        if let Ok(file) = File::open(path) {
            let reader = BufReader::new(file);
            for line_str in reader.lines().map_while(Result::ok) {
                let parts: Vec<&str> = line_str.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let (Ok(i), Ok(j)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) {
                        pairs.push((i, j));
                    }
                }
            }
        }
    }

    eprintln!("# Preparing NOE restraints");
    eprintln!("# Pairs: {}", pairs.len());
    eprintln!("# Force constant: {} kJ/(mol·nm²)", force_const);
    eprintln!("# Target distance: {} nm", target_dist);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  NOE distance restraints").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "DISTANCERES").map_err(io_err)?;
    writeln!(out, "  # Number of restraints").map_err(io_err)?;
    writeln!(out, "    {}", pairs.len()).map_err(io_err)?;

    for (idx, &(i, j)) in pairs.iter().enumerate() {
        writeln!(out, "  # Restraint {}: atoms {} - {}", idx + 1, i, j).map_err(io_err)?;
        writeln!(
            out,
            "    {:6} {:6} {:12.4} {:12.4}",
            i, j, target_dist, force_const
        )
        .map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;

    eprintln!("# NOE preparation complete");
    Ok(())
}
