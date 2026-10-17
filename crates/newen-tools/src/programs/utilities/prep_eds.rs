//! prep_eds - Prepare EDS simulation files
//!
//! Usage: prep_eds @topo <file> @states <n>
//!
//! Prepares topology and parameters for EDS simulations

use crate::programs::io_err;
use newen::io::topology::{build_topology, read_topology_file};
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "prep_eds - Report an EDS setup\n\nUsage: prep_eds @topo <file> @states <n> [@emin <value>]\n\nArguments:\n  @topo    Molecular topology file\n  @states  Number of EDS states\n  @emin    Energy offset\n\nNote: gromos++'s prep_eds writes the EDS perturbation topology; this one reports.";

pub const KNOWN: &[&str] = &["topo", "states", "emin"];

pub const PROGRAM: Program = Program {
    name: "prep_eds",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let n_states: usize = args.get("states", 2)?;
    let emin: f64 = args.get("emin", -100.0)?;

    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;
    let _topo = build_topology(topo_data);

    eprintln!("# Preparing EDS simulation");
    eprintln!("# Number of states: {}", n_states);
    eprintln!("# E_min: {} kJ/mol", emin);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  EDS simulation parameters").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "EDS").map_err(io_err)?;
    writeln!(out, "  # Number of states").map_err(io_err)?;
    writeln!(out, "    {}", n_states).map_err(io_err)?;
    writeln!(out, "  # E_min (kJ/mol)").map_err(io_err)?;
    writeln!(out, "    {:.4}", emin).map_err(io_err)?;
    writeln!(out, "  # State energies (initial)").map_err(io_err)?;

    for i in 0..n_states {
        writeln!(out, "    # State {}", i + 1).map_err(io_err)?;
        writeln!(out, "      0.0").map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;

    eprintln!("# EDS preparation complete");
    Ok(())
}
