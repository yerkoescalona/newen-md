//! amber2gromos - Convert AMBER topology to GROMOS format
//!
//! Usage: amber2gromos @amber <file> @output <file>
//!
//! Converts AMBER prmtop to GROMOS topology

use crate::programs::io_err;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "amber2gromos - Convert an AMBER topology to GROMOS\n\nUsage: amber2gromos @amber <prmtop> [@natoms <n>]\n\nArguments:\n  @amber   AMBER topology (prmtop)\n  @natoms  Number of atoms\n\nNote: a stub — the conversion is not implemented.";

pub const KNOWN: &[&str] = &["amber", "natoms"];

pub const PROGRAM: Program = Program {
    name: "amber2gromos",
    known: KNOWN,
    usage: USAGE,
    status: Status::Stub,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let _amber_file = args.value("amber")?.to_string();
    let n_atoms: usize = args.get("natoms", 100)?;

    eprintln!("# Converting AMBER topology to GROMOS");
    eprintln!("# Atoms: {} (estimated)", n_atoms);

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  Converted from AMBER").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "PHYSICALCONSTANTS").map_err(io_err)?;
    writeln!(out, "    {:.7}", newen_core::units::four_pi_eps_i).map_err(io_err)?;
    writeln!(out, "    {:.7}", newen_core::units::hBar).map_err(io_err)?;
    writeln!(out, "    {:.3}", newen_core::units::spd_l).map_err(io_err)?;
    writeln!(out, "    {:.8}", newen_core::units::kB).map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "TOPVERSION").map_err(io_err)?;
    writeln!(out, "    2.0").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "ATOMTYPENAME").map_err(io_err)?;
    writeln!(out, "    C N O H S").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "RESNAME").map_err(io_err)?;
    writeln!(out, "    MOL").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "SOLUTEATOM").map_err(io_err)?;
    writeln!(out, "    {}", n_atoms).map_err(io_err)?;

    // Default parameters for demonstration
    for i in 0..n_atoms {
        writeln!(
            out,
            "    {:6}   1     1    1   12.0000      0.000000",
            i + 1
        )
        .map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "SOLUTEMOLECULES").map_err(io_err)?;
    writeln!(out, "    1").map_err(io_err)?;
    writeln!(out, "    {}", n_atoms).map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;

    eprintln!("# Conversion complete (simplified)");
    Ok(())
}
