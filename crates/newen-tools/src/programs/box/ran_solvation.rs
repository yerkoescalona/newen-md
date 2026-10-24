//! ran_solvation - Random solvation of solute
//!
//! Usage: ran_solvation @solute <file> @solvent <file> @box <x y z> @nsolv <n>
//!
//! Randomly places solvent molecules around solute

use crate::programs::io_err;
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "ran_solvation - Random solvent placement around a solute\n\nUsage: ran_solvation @solute <file> @solvent <file> @nsolv <n> @box <x> <y> <z> [@seed <n>]\n\nArguments:\n  @solute   Solute coordinate file\n  @solvent  Solvent coordinate file\n  @nsolv    Number of solvent molecules\n  @box      Box dimensions (three values, nm)\n  @seed     Random seed (default 12345)\n\nNote: this is not gromos++'s ran_solvation, which takes a different argument model.";

pub const KNOWN: &[&str] = &["solute", "solvent", "nsolv", "box", "seed"];

pub const PROGRAM: Program = Program {
    name: "ran_solvation",
    known: KNOWN,
    usage: USAGE,
    status: Status::Stub,
    run,
};

fn lcg_random(state: &mut u64) -> f64 {
    const A: u64 = 1103515245;
    const C: u64 = 12345;
    const M: u64 = 2u64.pow(31);

    *state = (A.wrapping_mul(*state).wrapping_add(C)) % M;
    *state as f64 / M as f64
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let _solute_file = args.value("solute")?.to_string();
    let _solvent_file = args.value("solvent")?.to_string();
    let n_solvent: usize = args.get("nsolv", 1000)?;
    let seed: u64 = args.get("seed", 12345)?;

    let dims = args.values("box");
    let (box_x, box_y, box_z) = match dims.len() {
        0 => (5.0, 5.0, 5.0),
        3 => {
            let p = |s: &String| -> Result<f64, String> {
                s.parse::<f64>()
                    .map_err(|_| format!("@box: '{s}' is not a number"))
            };
            (p(&dims[0])?, p(&dims[1])?, p(&dims[2])?)
        },
        n => return Err(format!("@box takes three values, got {n}")),
    };

    eprintln!("# Random solvation");
    eprintln!("# Box: {} x {} x {} nm", box_x, box_y, box_z);
    eprintln!("# Solvent molecules: {}", n_solvent);

    let mut rng = seed;

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  Random solvation").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "POSITION").map_err(io_err)?;
    writeln!(out, "  {:15.9} {:15.9} {:15.9}", box_x, box_y, box_z).map_err(io_err)?;

    // Simplified: place 3-atom water molecules
    for _ in 0..n_solvent {
        let cx = lcg_random(&mut rng) * box_x;
        let cy = lcg_random(&mut rng) * box_y;
        let cz = lcg_random(&mut rng) * box_z;

        // Oxygen
        writeln!(out, "{:15.9} {:15.9} {:15.9}", cx, cy, cz).map_err(io_err)?;

        // H1
        writeln!(out, "{:15.9} {:15.9} {:15.9}", cx + 0.01, cy, cz).map_err(io_err)?;

        // H2
        writeln!(out, "{:15.9} {:15.9} {:15.9}", cx, cy + 0.01, cz).map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;

    eprintln!("# Solvation complete: {} atoms", n_solvent * 3);
    Ok(())
}
