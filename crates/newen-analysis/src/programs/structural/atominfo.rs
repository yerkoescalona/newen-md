//! atominfo — Report information about selected atoms (newen compatible).
//!
//! Usage:
//!   atominfo @topo <topology> @atomspec <spec> [@sort]
//!
//! Reads a GROMOS topology, resolves the AtomSpecifier, and prints the TITLE +
//! ATOMS block in newen format for every selected atom.
//!
//! Output columns (ATOMS block):
//!   mol:atom  gromos_nr  residue_nr  residue_name  atom_name  iac(1-based)  charge

use crate::programs::io_err;
use newen_core::selection::AtomSelection;
use newen_io::args::{Arguments, Program, Status};
use newen_io::topology::{build_topology, read_topology_file};
use std::io::Write;

pub const USAGE: &str = "# atominfo
\t@topo          <topology>
\t@atomspec      <atom specifier>
\t[@sort         <sort the selected atoms by global index (default: yes)>]";

pub const KNOWN: &[&str] = &["topo", "atomspec", "sort"];

pub const PROGRAM: Program = Program {
    name: "atominfo",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    if argv.len() < 2 {
        return Err(USAGE.to_string());
    }
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_path = args.value("topo")?;
    let spec = args.value("atomspec")?; // @sort is accepted and ignored: the selection is sorted

    let topo_data =
        read_topology_file(topo_path).map_err(|e| format!("Error reading topology: {e}"))?;
    let topo = build_topology(topo_data);

    let selection = AtomSelection::from_string(spec, &topo)
        .map_err(|e| format!("Error in atom specifier '{spec}': {e}"))?;

    // Print TITLE block
    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "\tatominfo").map_err(io_err)?;
    writeln!(out, "\t{spec}").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;

    // Print ATOMS block
    writeln!(out, "ATOMS").map_err(io_err)?;
    for &global_idx in selection.indices() {
        let mol_nr = topo.molecule_nr(global_idx).map(|m| m + 1).unwrap_or(1); // 1-based
        let local_nr = if topo.molecules.is_empty() {
            global_idx + 1
        } else {
            let start = topo.molecules.get(mol_nr - 1).map(|r| r.start).unwrap_or(0);
            global_idx - start + 1
        };
        let gromos_nr = global_idx + 1;
        let res_nr = topo.residue_nr(global_idx).unwrap_or(0);
        let res_name = topo.residue_name(global_idx).unwrap_or("???");
        let atom_name = topo.atom_name(global_idx).unwrap_or("???");
        let iac = topo.iac.get(global_idx).map(|&v| v + 1).unwrap_or(0); // 1-based
        let charge = topo.charge.get(global_idx).copied().unwrap_or(0.0);

        writeln!(
            out,
            "{:>11}  {:>9}  {:>9}  {:>10}  {:>9}  {:>9}  {:>9}",
            format!("{mol_nr}:{local_nr}"),
            gromos_nr,
            res_nr,
            res_name,
            atom_name,
            iac,
            charge,
        )
        .map_err(io_err)?;
    }
    writeln!(out, "END").map_err(io_err)?;
    Ok(())
}
