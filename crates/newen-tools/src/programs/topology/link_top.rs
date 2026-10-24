//! link_top - Link multiple topology files
//!
//! Usage: link_top @topo1 <file> @topo2 <file>
//!
//! Combines topologies into a single file

use crate::programs::io_err;
use newen::io::topology::read_topology_file;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "link_top - Report how topologies would link\n\nUsage: link_top @topo1 <file> @topo2 <file> [@topo3 <file>]\n\nArguments:\n  @topo1  First molecular topology file\n  @topo2  Second molecular topology file\n  @topo3  Third molecular topology file\n\nNote: gromos++'s link_top applies link building blocks; this one reports.";

pub const KNOWN: &[&str] = &["topo1", "topo2", "topo3", "topo"];

pub const PROGRAM: Program = Program {
    name: "link_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    // gromos++ numbers the keys; @topo takes them all at once
    let mut topo_files: Vec<String> = Vec::new();
    for key in ["topo", "topo1", "topo2", "topo3"] {
        topo_files.extend(args.values(key).iter().cloned());
    }

    if topo_files.is_empty() {
        return Err("Error: No topology files specified".to_string());
    }

    eprintln!("# Linking {} topology files", topo_files.len());

    let mut all_masses = Vec::new();
    let mut all_charges = Vec::new();
    let mut all_iacs = Vec::new();

    for (idx, topo_file) in topo_files.iter().enumerate() {
        eprintln!("# Reading topology {}: {}", idx + 1, topo_file);

        match read_topology_file(topo_file) {
            Ok(topo_data) => {
                all_masses.extend(topo_data.masses);
                all_charges.extend(topo_data.charges);
                all_iacs.extend(topo_data.iac);
            },
            Err(e) => {
                eprintln!("Warning: Could not read {}: {:?}", topo_file, e);
            },
        }
    }

    let total_atoms = all_masses.len();
    eprintln!("# Total atoms in linked topology: {}", total_atoms);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  Linked topology").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "TOPVERSION").map_err(io_err)?;
    writeln!(out, "    2.0").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "ATOMTYPENAME").map_err(io_err)?;
    writeln!(out, "    C").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "RESNAME").map_err(io_err)?;
    writeln!(out, "    MOL").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "SOLUTEATOM").map_err(io_err)?;
    writeln!(out, "  # NRP: number of solute atoms").map_err(io_err)?;
    writeln!(out, "    {}", total_atoms).map_err(io_err)?;

    for i in 0..total_atoms {
        writeln!(out, "  # {:6} RES  ATOM  IAC   MASS    CHARGE", i + 1).map_err(io_err)?;
        writeln!(
            out,
            "    {:6}   1     1  {:3}  {:8.4}  {:10.6}",
            i + 1,
            all_iacs[i],
            all_masses[i],
            all_charges[i]
        )
        .map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "SOLUTEMOLECULES").map_err(io_err)?;
    writeln!(out, "  # NSPM: number of separate molecules").map_err(io_err)?;
    writeln!(out, "    1").map_err(io_err)?;
    writeln!(out, "  # NSP[1...NSPM]: atom sequence number of last atom").map_err(io_err)?;
    writeln!(out, "    {}", total_atoms).map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    Ok(())
}
