//! copy_box — duplicate a periodic box along its axes (gromos++ `copy_box`).
//!
//! Usage: copy_box @topo <file> @pos <file> @dir <x[:N],y[:N],z[:N]>
//!    or: copy_box @topo <file> @traj <file> @replicate <nx> <ny> <nz>
//!
//! As gromos++ does: the images are built direction by direction (each existing image is
//! shifted by 1..N box lengths), the solute molecules of every image come first, then the
//! solvent of every image, and the box grows by N lengths along each duplicated axis.

use crate::programs::io_err;
use newen::math::Vec3;
use newen_io::args::{Arguments, Program, Status};
use newen_io::topology::read_topology_file;
use newen_io::{read_g96_labeled, G96Atom};
use std::fmt::Write as FmtWrite;
use std::io::Write;

pub const USAGE: &str = "copy_box - Replicate a box along the axes\n\nUsage: copy_box @topo <file> @pos <file> @dir <x:N,y:M,...>\n   or: copy_box @topo <file> @pos <file> @replicate <nx> <ny> <nz>\n\nArguments:\n  @topo       Molecular topology file\n  @pos        Input coordinate file (also accepted as @traj)\n  @dir        Directions and multiplicities, e.g. \"x:1,y:2\"\n  @replicate  Copies per dimension (three values; 1 means no copy)";

pub const KNOWN: &[&str] = &["topo", "pos", "traj", "dir", "replicate"];

pub const PROGRAM: Program = Program {
    name: "copy_box",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    // gromos++ accepts the coordinates under either name
    let pos_file = if args.has("pos") {
        args.value("pos")?
    } else {
        args.value("traj").map_err(|_| USAGE.to_string())?
    };

    let mut dirs: Vec<(char, usize)> = Vec::new();
    for item in args.values("dir").iter().flat_map(|v| v.split(',')) {
        let (axis, n) = item.split_once(':').unwrap_or((item, "1"));
        let axis = axis
            .trim()
            .chars()
            .next()
            .ok_or_else(|| USAGE.to_string())?;
        if !"xyz".contains(axis) {
            return Err(
                "Error: only x, y and z directions are supported (rectangular boxes)".to_string(),
            );
        }
        let n: usize = n.trim().parse().map_err(|_| USAGE.to_string())?;
        if n == 0 {
            return Err("Error: multiplicity in @dir must be >= 1".to_string());
        }
        dirs.push((axis, n));
    }
    if args.has("replicate") {
        let v = args.values("replicate");
        if v.len() != 3 {
            return Err(format!("@replicate takes three values, got {}", v.len()));
        }
        for (k, axis) in ['x', 'y', 'z'].into_iter().enumerate() {
            let n: usize = v[k].parse().map_err(|_| USAGE.to_string())?;
            if n > 1 {
                dirs.push((axis, n - 1));
            }
        }
    }

    let topo = read_topology_file(topo_file)
        .map_err(|e| format!("Error reading topology '{topo_file}': {e}"))?;
    let data =
        read_g96_labeled(pos_file).map_err(|e| format!("Error reading '{pos_file}': {e}"))?;
    if data.box_dims == Vec3::ZERO {
        return Err(format!("Error: '{pos_file}' has no box"));
    }
    let n_solute = topo.n_atoms;
    if data.atoms.len() < n_solute {
        return Err(format!(
            "Error: {} atoms in the coordinates, {n_solute} in the topology",
            data.atoms.len()
        ));
    }
    let per_solvent = topo.solvent_atoms.len().max(1);
    let base = data.box_dims;
    let axis_shift = |axis: char| match axis {
        'x' => Vec3::new(base.x, 0.0, 0.0),
        'y' => Vec3::new(0.0, base.y, 0.0),
        _ => Vec3::new(0.0, 0.0, base.z),
    };
    // images, as gromos++ builds them: for each direction, every existing image shifted 1..N times
    let mut images: Vec<Vec3> = vec![Vec3::ZERO];
    for &(axis, n) in &dirs {
        let mut new_images = Vec::new();
        for k in 1..=n {
            let shift = axis_shift(axis) * k as f64;
            for &im in &images {
                new_images.push(im + shift);
            }
        }
        images.extend(new_images);
    }
    let mut new_box = base;
    for &(axis, n) in &dirs {
        match axis {
            'x' => new_box.x += base.x * n as f64,
            'y' => new_box.y += base.y * n as f64,
            _ => new_box.z += base.z * n as f64,
        }
    }
    let n_res = topo.residue_names.len();
    let (solute, solvent) = data.atoms.split_at(n_solute);
    let n_solvent_mols = solvent.len() / per_solvent;
    eprintln!(
        "# {} images of {} solute atoms and {} solvent molecules; box ({:.6}, {:.6}, {:.6})",
        images.len(),
        n_solute,
        n_solvent_mols,
        new_box.x,
        new_box.y,
        new_box.z
    );
    writeln!(out, "TITLE").map_err(io_err)?;
    writeln!(
        out,
        "Copy_box: {} duplicated with dir={} ({} copies)",
        pos_file,
        dirs.iter()
            .map(|(a, n)| format!("{a}:{n}"))
            .collect::<Vec<_>>()
            .join(","),
        images.len()
    )
    .map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "POSITION").map_err(io_err)?;
    let mut positions = String::new();
    let mut line = |res: usize, res_name: &str, atom: &G96Atom, serial: usize, p: Vec3| {
        let _ = writeln!(
            positions,
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}",
            res, res_name, atom.atom_name, serial, p.x, p.y, p.z
        );
    };
    let mut serial = 0;
    for (k, &shift) in images.iter().enumerate() {
        for a in solute {
            serial += 1;
            line(a.res_num + k * n_res, &a.res_name, a, serial, a.pos + shift);
        }
    }
    let mut mol = 0;
    for &shift in &images {
        for (j, a) in solvent.iter().enumerate() {
            if j % per_solvent == 0 {
                mol += 1;
            }
            serial += 1;
            line(mol, "SOLV", a, serial, a.pos + shift);
        }
    }
    write!(out, "{positions}").map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    writeln!(out, "GENBOX").map_err(io_err)?;
    writeln!(out, "{:>8}", 1).map_err(io_err)?;
    writeln!(
        out,
        "{:15.9}{:15.9}{:15.9}",
        new_box.x, new_box.y, new_box.z
    )
    .map_err(io_err)?;
    writeln!(out, "{:15.9}{:15.9}{:15.9}", 90.0, 90.0, 90.0).map_err(io_err)?;
    writeln!(out, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0).map_err(io_err)?;
    writeln!(out, "{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0).map_err(io_err)?;
    writeln!(out, "END").map_err(io_err)?;
    Ok(())
}
