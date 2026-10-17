//! addvirt_top - Add virtual atoms to topology
//!
//! Usage: addvirt_top @topo <file> @nvirt <n>
//!
//! Adds virtual (massless) atoms to topology

use crate::programs::io_err;
use newen::io::topology::read_topology_file;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "addvirt_top - Report virtual atoms that would be added\n\nUsage: addvirt_top @topo <file> @nvirt <n>\n\nArguments:\n  @topo   Molecular topology file\n  @nvirt  Number of virtual atoms\n\nNote: gromos++'s addvirt_top writes the modified topology; this one reports.";

pub const KNOWN: &[&str] = &["topo", "nvirt"];

pub const PROGRAM: Program = Program {
    name: "addvirt_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let n_virtual: usize = args.get("nvirt", 10)?;

    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;

    let original_atoms = topo_data.masses.len();
    let total_atoms = original_atoms + n_virtual;

    eprintln!("# Adding {} virtual atoms to topology", n_virtual);
    eprintln!("# Original atoms: {}", original_atoms);
    eprintln!("# Total atoms: {}", total_atoms);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  Topology with virtual atoms").map_err(io_err)?;
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

    // Print original atoms
    for i in 0..original_atoms {
        writeln!(out, "  # {:6} RES  ATOM  IAC   MASS    CHARGE", i + 1).map_err(io_err)?;
        writeln!(
            out,
            "    {:6}   1     1  {:3}  {:8.4}  {:10.6}",
            i + 1,
            topo_data.iac[i],
            topo_data.masses[i],
            topo_data.charges[i]
        )
        .map_err(io_err)?;
    }

    // Print virtual atoms (mass = 0)
    for i in 0..n_virtual {
        let atom_idx = original_atoms + i + 1;
        writeln!(
            out,
            "  # {:6} RES  ATOM  IAC   MASS    CHARGE (virtual)",
            atom_idx
        )
        .map_err(io_err)?;
        writeln!(
            out,
            "    {:6}   1     1    1    0.0000      0.000000",
            atom_idx
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
