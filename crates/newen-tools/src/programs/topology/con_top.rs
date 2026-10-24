//! con_top - Concatenate topology files
//!
//! Usage: con_top @topo <file1> @topo <file2> ...
//!
//! Combines multiple topology files into one

use crate::programs::io_err;
use newen::io::topology::{build_topology, read_topology_file};
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "con_top - Concatenate topologies\n\nUsage: con_top @topo <file> <file> [...]\n\nArguments:\n  @topo  Molecular topology files, in order\n\nNote: gromos++'s con_top converts force-field parameters; this one concatenates.";

pub const KNOWN: &[&str] = &["topo"];

pub const PROGRAM: Program = Program {
    name: "con_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_files = args.values("topo").to_vec();

    if topo_files.len() < 2 {
        return Err("Error: Need at least 2 topology files".to_string());
    }

    let mut all_topos = Vec::new();

    for topo_file in &topo_files {
        let topo_data =
            read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;
        let topo = build_topology(topo_data);
        all_topos.push(topo);
        eprintln!(
            "Read: {} ({} atoms)",
            topo_file,
            all_topos.last().unwrap().mass.len()
        );
    }

    let total_atoms: usize = all_topos.iter().map(|t| t.mass.len()).sum();
    let total_bonds: usize = all_topos.iter().map(|t| t.moltypes[0].bonds.len()).sum();

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "Combined topology from {} files", topo_files.len()).map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out).map_err(io_err)?;
    writeln!(out, "# Total atoms: {}", total_atoms).map_err(io_err)?;
    writeln!(out, "# Total bonds: {}", total_bonds).map_err(io_err)?;
    writeln!(out).map_err(io_err)?;

    let mut atom_offset = 0;
    for (idx, topo) in all_topos.iter().enumerate() {
        writeln!(out, "# From: {}", topo_files[idx]).map_err(io_err)?;
        for (i, &mass) in topo.mass.iter().enumerate() {
            writeln!(
                out,
                "# Atom {:6}: Mass {:8.4}, Charge {:10.6}",
                atom_offset + i + 1,
                mass,
                topo.charge[i]
            )
            .map_err(io_err)?;
        }
        atom_offset += topo.mass.len();
    }

    eprintln!();
    eprintln!(
        "Combined {} atoms from {} topologies",
        total_atoms,
        topo_files.len()
    );
    Ok(())
}
