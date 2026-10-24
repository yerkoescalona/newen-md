//! inbox — put every atom into the periodic box (gromos++ `inbox`).
//!
//! Usage: inbox @topo <file> @pos <file> @pbc <r|v> [@shift <x> <y> <z>]
//!    or: inbox @topo <file> @traj <file> @pbc <r|v>
//!
//! gromos++: each atom is replaced by its periodic image nearest to the box centre
//! (K/2, L/2, M/2), after an optional shift; vacuum leaves the coordinates alone.

use crate::programs::io_err;
use newen::math::Vec3;
use newen_io::args::{Arguments, Program, Status};
use newen_io::read_g96_labeled;
use std::io::Write;

pub const USAGE: &str = "inbox - Put every atom into the central box\n\nUsage: inbox @topo <file> @pos <file> [@pbc <r|v>] [@shift <x> <y> <z>]\n\nArguments:\n  @topo   Molecular topology file (read, but every atom is treated alike)\n  @pos    Input coordinate file (also accepted as @traj)\n  @pbc    Periodic boundary conditions (r = rectangular, v = vacuum)\n  @shift  Shift applied before boxing (three values, nm)";

pub const KNOWN: &[&str] = &["topo", "pos", "traj", "pbc", "shift"];

pub const PROGRAM: Program = Program {
    name: "inbox",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// gromos++ rectangular `nearestImage(r, centre)`: the image of `r` nearest to `centre`.
fn nearest_image_to(r: Vec3, centre: Vec3, b: Vec3) -> Vec3 {
    let d = r - centre;
    Vec3::new(
        r.x - b.x * (d.x / b.x).round(),
        r.y - b.y * (d.y / b.y).round(),
        r.z - b.z * (d.z / b.z).round(),
    )
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    // @topo is accepted and unused: every atom is treated alike
    let pos_file = if args.has("pos") {
        args.value("pos")?
    } else {
        args.value("traj").map_err(|_| USAGE.to_string())?
    };
    let pbc = match args.has("pbc") {
        false => 'r',
        true => args
            .value("pbc")?
            .chars()
            .next()
            .ok_or_else(|| USAGE.to_string())?,
    };
    let shift = if args.has("shift") {
        let v = args.values("shift");
        if v.len() != 3 {
            return Err(format!("@shift takes three values, got {}", v.len()));
        }
        let p = |s: &String| -> Result<f64, String> {
            s.parse::<f64>()
                .map_err(|_| format!("@shift: '{s}' is not a number"))
        };
        Vec3::new(p(&v[0])?, p(&v[1])?, p(&v[2])?)
    } else {
        Vec3::ZERO
    };

    let data =
        read_g96_labeled(pos_file).map_err(|e| format!("Error reading '{pos_file}': {e}"))?;
    let b = data.box_dims;
    if pbc != 'v' && b == Vec3::ZERO {
        return Err(format!("Error: no box in '{pos_file}'"));
    }
    let centre = b * 0.5;
    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(out, "inbox: {}", pos_file).map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "POSITION").map_err(io_err)?;
    for (k, a) in data.atoms.iter().enumerate() {
        let p = a.pos + shift;
        let p = if pbc == 'v' {
            p
        } else {
            nearest_image_to(p, centre, b)
        };
        writeln!(
            out,
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}",
            a.res_num,
            a.res_name,
            a.atom_name,
            k + 1,
            p.x,
            p.y,
            p.z
        )
        .map_err(io_err)?;
    }
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "GENBOX").map_err(io_err)?;
    writeln!(out, "{:>8}", if pbc == 'v' { 0 } else { 1 }).map_err(io_err)?;
    writeln!(out, "{:15.9}{:15.9}{:15.9}", b.x, b.y, b.z).map_err(io_err)?;
    writeln!(out, "{:15.9}{:15.9}{:15.9}", 90.0, 90.0, 90.0).map_err(io_err)?;
    writeln!(out, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0).map_err(io_err)?;
    writeln!(out, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0).map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    Ok(())
}
