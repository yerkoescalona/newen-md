//! pert_top - Create perturbation topology
//!
//! Usage: pert_top @topo1 <file> @topo2 <file>
//!
//! Creates perturbation topology for free energy calculations

use crate::programs::io_err;
use newen::io::topology::read_topology_file;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "pert_top - Compare two topologies\n\nUsage: pert_top @topo1 <file> @topo2 <file>\n\nArguments:\n  @topo1  State A topology\n  @topo2  State B topology\n\nNote: gromos++'s pert_top modifies selected atoms; this one reports the differences.";

pub const KNOWN: &[&str] = &["topo1", "topo2"];

pub const PROGRAM: Program = Program {
    name: "pert_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo1_file = args.value("topo1")?;
    let topo2_file = args.value("topo2")?;

    let topo1_data =
        read_topology_file(topo1_file).map_err(|e| format!("@topo1 '{topo1_file}': {e}"))?;
    let topo2_data =
        read_topology_file(topo2_file).map_err(|e| format!("@topo2 '{topo2_file}': {e}"))?;

    let n1 = topo1_data.masses.len();
    let n2 = topo2_data.masses.len();
    let n_atoms = n1.max(n2);

    eprintln!("# Creating perturbation topology");
    eprintln!("# State A atoms: {}", n1);
    eprintln!("# State B atoms: {}", n2);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  Perturbation topology").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "PERTATOM").map_err(io_err)?;
    writeln!(out, "  # Number of perturbed atoms").map_err(io_err)?;
    writeln!(out, "    {}", n_atoms).map_err(io_err)?;

    for i in 0..n_atoms {
        let mass_a = if i < n1 { topo1_data.masses[i] } else { 0.0 };
        let charge_a = if i < n1 { topo1_data.charges[i] } else { 0.0 };
        let iac_a = if i < n1 { topo1_data.iac[i] } else { 0 };

        let mass_b = if i < n2 { topo2_data.masses[i] } else { 0.0 };
        let charge_b = if i < n2 { topo2_data.charges[i] } else { 0.0 };
        let iac_b = if i < n2 { topo2_data.iac[i] } else { 0 };

        writeln!(out, "  # Atom {}", i + 1).map_err(io_err)?;
        writeln!(
            out,
            "    # State A: mass={:.4} charge={:.6} IAC={}",
            mass_a, charge_a, iac_a
        )
        .map_err(io_err)?;
        writeln!(
            out,
            "    # State B: mass={:.4} charge={:.6} IAC={}",
            mass_b, charge_b, iac_b
        )
        .map_err(io_err)?;
        writeln!(
            out,
            "    {:6} {:10.4} {:12.6} {:4} {:10.4} {:12.6} {:4}",
            i + 1,
            mass_a,
            charge_a,
            iac_a,
            mass_b,
            charge_b,
            iac_b
        )
        .map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;

    eprintln!("# Perturbation topology created");
    Ok(())
}
