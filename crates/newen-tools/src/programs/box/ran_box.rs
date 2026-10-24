//! ran_box - Generate random box configuration
//!
//! Usage: ran_box @topo <file> @box <x y z> @seed <n>
//!
//! Creates random positions within box dimensions

use crate::programs::io_err;
use newen::io::topology::{build_topology, read_topology_file};
use newen_io::args::{Arguments, Program, Status};
use std::io::Write;

pub const USAGE: &str = "ran_box - Random positions in a box\n\nUsage: ran_box @topo <file> @box <x> <y> <z> [@seed <n>]\n\nArguments:\n  @topo  Molecular topology file (for the atom count)\n  @box   Box dimensions (three values, nm)\n  @seed  Random seed (default 12345)\n\nNote: this is not gromos++'s ran_box, which places whole molecules from a building block.";

pub const KNOWN: &[&str] = &["topo", "box", "seed"];

pub const PROGRAM: Program = Program {
    name: "ran_box",
    known: KNOWN,
    usage: USAGE,
    status: Status::Stub,
    run,
};

fn lcg_random(state: &mut u64) -> f64 {
    // Linear congruential generator
    const A: u64 = 1103515245;
    const C: u64 = 12345;
    const M: u64 = 2u64.pow(31);

    *state = (A.wrapping_mul(*state).wrapping_add(C)) % M;
    *state as f64 / M as f64
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let seed: u64 = args.get("seed", 12345)?;

    let dims = args.values("box");
    let (box_x, box_y, box_z) = match dims.len() {
        0 => (3.0, 3.0, 3.0),
        3 => {
            let p = |s: &String| -> Result<f64, String> {
                s.parse::<f64>()
                    .map_err(|_| format!("@box: '{s}' is not a number"))
            };
            (p(&dims[0])?, p(&dims[1])?, p(&dims[2])?)
        },
        n => return Err(format!("@box takes three values, got {n}")),
    };

    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("@topo '{topo_file}': {e}"))?;
    let topo = build_topology(topo_data);
    let num_atoms = topo.mass.len();

    eprintln!("# Generating random positions for {} atoms", num_atoms);
    eprintln!("# Box: {} x {} x {} nm", box_x, box_y, box_z);
    eprintln!("# Seed: {}", seed);

    let mut rng_state = seed;

    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "  Random box configuration").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "POSITION").map_err(io_err)?;
    writeln!(out, "  {:15.9} {:15.9} {:15.9}", box_x, box_y, box_z).map_err(io_err)?;

    for _ in 0..num_atoms {
        let x = lcg_random(&mut rng_state) * box_x;
        let y = lcg_random(&mut rng_state) * box_y;
        let z = lcg_random(&mut rng_state) * box_z;
        writeln!(out, "{:15.9} {:15.9} {:15.9}", x, y, z).map_err(io_err)?;
    }

    writeln!(out, "END").map_err(io_err)?;
    Ok(())
}
