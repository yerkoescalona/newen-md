//! prep_posres - Generate position restraint files (.por and .rpr)
//!
//! Usage: prep_posres @topo <topology> @pos <coordinates.cnf>
//!                    [@por <output.por>] [@rpr <output.rpr>]
//!
//! Creates two files from a solvated coordinate file:
//!   .por (POSRESSPEC) — lists solute atoms to be restrained
//!   .rpr (REFPOSITION) — reference positions for all atoms
//!
//! The number of solute atoms is determined from the topology.

use crate::programs::parse_clap;
use clap::Parser;
use newen_io::args::{Program, Status};
use newen_io::topology::read_topology_file;
use newen_io::{read_g96_labeled, write_por, write_rpr};
use std::io::Write;

pub const USAGE: &str = "prep_posres - Write position-restraint and reference-position files\n\nUsage: prep_posres @topo <file> @pos <file> [@por <file>] [@rpr <file>]\n\nArguments:\n  @topo  Molecular topology file\n  @pos   Coordinate file\n  @por   Position-restraint specification to write (default <pos>.por)\n  @rpr   Reference positions to write (default <pos>.rpr)\n\nBoth output files are written to disk; stdout carries nothing.";

pub const KNOWN: &[&str] = &["topo", "pos", "por", "rpr"];

pub const PROGRAM: Program = Program {
    name: "prep_posres",
    known: KNOWN,
    usage: USAGE,
    status: Status::NoCounterpart,
    run,
};

#[derive(Parser)]
#[command(
    name = "prep_posres",
    version,
    about = "Generate position restraint files (.por/.rpr)"
)]
struct Args {
    /// Molecular topology file
    #[arg(long)]
    topo: String,

    /// Input coordinate file (solvated system)
    #[arg(long)]
    pos: String,

    /// Output .por file (POSRESSPEC, solute atoms only)
    #[arg(long)]
    por: Option<String>,

    /// Output .rpr file (REFPOSITION, all atoms)
    #[arg(long)]
    rpr: Option<String>,
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let Some(args) = parse_clap::<Args>(argv, out)? else {
        return Ok(());
    };

    let topo =
        read_topology_file(&args.topo).map_err(|e| format!("Error reading topology: {e}"))?;
    let n_solute = topo.n_atoms;

    let coord =
        read_g96_labeled(&args.pos).map_err(|e| format!("Error reading coordinates: {e}"))?;

    eprintln!("prep_posres:");
    eprintln!("  Topology:     {}", args.topo);
    eprintln!("  Coordinates:  {}", args.pos);
    eprintln!("  Total atoms:  {}", coord.atoms.len());
    eprintln!("  Solute atoms: {}", n_solute);

    if n_solute > coord.atoms.len() {
        return Err(format!(
            "Error: topology has {n_solute} solute atoms but coordinate file has only {} atoms",
            coord.atoms.len()
        ));
    }

    // Derive default names from input
    let base = args.pos.trim_end_matches(".cnf").trim_end_matches(".g96");

    // Write .por
    let por_path = args.por.unwrap_or_else(|| format!("{}.por", base));
    write_por(&por_path, &coord.atoms, n_solute).map_err(|e| format!("Error writing .por: {e}"))?;
    eprintln!("  Wrote .por:   {} ({} solute atoms)", por_path, n_solute);

    // Write .rpr
    let rpr_path = args.rpr.unwrap_or_else(|| format!("{}.rpr", base));
    write_rpr(&rpr_path, &coord.atoms).map_err(|e| format!("Error writing .rpr: {e}"))?;
    eprintln!("  Wrote .rpr:   {} ({} atoms)", rpr_path, coord.atoms.len());
    Ok(())
}
