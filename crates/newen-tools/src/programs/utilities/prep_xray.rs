//! prep_xray - Prepare X-ray refinement files
//!
//! Usage: prep_xray @topo <file> @resolution <d>
//!
//! Prepares X-ray restraint topology

use crate::programs::io_err;
use newen::io::topology::{build_topology, read_topology_file};
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "prep_xray - Report an X-ray restraint setup\n\nUsage: prep_xray @topo <file> [@resolution <nm>] [@weight <w>]\n\nArguments:\n  @topo        Molecular topology file\n  @resolution  Resolution (nm)\n  @weight      Restraint weight\n\nNote: gromos++'s prep_xray needs clipper; this one reports a skeleton.";

pub const KNOWN: &[&str] = &["topo", "resolution", "weight"];

pub const PROGRAM: Program = Program {
    name: "prep_xray",
    known: KNOWN,
    usage: USAGE,
    status: Status::NoCounterpart,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let resolution: f64 = args.get("resolution", 2.0)?;
    let weight: f64 = args.get("weight", 1.0)?;

    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;
    let topo = build_topology(topo_data);

    let n_atoms = topo.mass.len();

    eprintln!("# Preparing X-ray refinement");
    eprintln!("# Atoms: {}", n_atoms);
    eprintln!("# Resolution: {} Å", resolution);
    eprintln!("# Weight: {}", weight);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  X-ray refinement parameters").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "XRAYRES").map_err(io_err)?;
    writeln!(out, "  # Resolution (Å)").map_err(io_err)?;
    writeln!(out, "    {:.4}", resolution).map_err(io_err)?;
    writeln!(out, "  # Weight").map_err(io_err)?;
    writeln!(out, "    {:.4}", weight).map_err(io_err)?;
    writeln!(out, "  # Number of atoms").map_err(io_err)?;
    writeln!(out, "    {}", n_atoms).map_err(io_err)?;
    writeln!(out, "  # B-factors (initial)").map_err(io_err)?;

    for i in 0..n_atoms {
        writeln!(out, "    {:6} {:8.2}", i + 1, 20.0).map_err(io_err)?; // Default B-factor
    }

    writeln!(out, "END").map_err(io_err)?;

    eprintln!("# X-ray preparation complete");
    Ok(())
}
