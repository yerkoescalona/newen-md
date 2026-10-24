//! top2build - Convert topology to build file
//!
//! Usage: top2build @topo <file>
//!
//! Converts GROMOS topology to building block format

use crate::programs::io_err;
use newen::io::topology::read_topology_file;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "top2build - Write a topology's solute as an MTB building block\n\nUsage: top2build @topo <file>\n\nArguments:\n  @topo  Molecular topology file";

pub const KNOWN: &[&str] = &["topo"];

pub const PROGRAM: Program = Program {
    name: "top2build",
    known: KNOWN,
    usage: USAGE,
    status: Status::NoCounterpart,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;

    eprintln!("# Converting topology to build file");
    eprintln!("# Atoms: {}", topo_data.masses.len());

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  Building blocks from topology").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "BUILDBLOCK").map_err(io_err)?;
    writeln!(out, "  # Number of blocks: 1").map_err(io_err)?;
    writeln!(out, "  # Block 1").map_err(io_err)?;
    writeln!(out, "  MOL").map_err(io_err)?;
    writeln!(out, "  # Number of atoms").map_err(io_err)?;
    writeln!(out, "    {}", topo_data.masses.len()).map_err(io_err)?;

    for (i, mass) in topo_data.masses.iter().enumerate() {
        let charge = topo_data.charges.get(i).copied().unwrap_or(0.0);
        let iac = topo_data.iac.get(i).copied().unwrap_or(1);

        writeln!(out, "  # Atom {}", i + 1).map_err(io_err)?;
        writeln!(
            out,
            "    {:6} {:10.4} {:12.6} {:4}",
            i + 1,
            mass,
            charge,
            iac
        )
        .map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;

    eprintln!("# Conversion complete");
    Ok(())
}
