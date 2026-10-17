//! explode — put solute molecules on a cubic grid at a given spacing (gromos++ `explode`).
//!
//! Usage: explode @topo <top> @pos <cnf> @nsm <n molecules> @dist <spacing nm>
//!
//! The grid has ceil(n^(1/3)) points per side; molecule i goes to grid point i with its centre
//! of mass at the point's centre; the box is the grid's size. Solvent is dropped.

use crate::programs::io_err;
use newen_core::Vec3;
use newen_io::args::{Arguments, Program, Status};
use newen_io::coordinate::format_g96;
use newen_io::read_g96_labeled;
use newen_io::topology::read_topology_file;
use std::io::Write;

pub const USAGE: &str = "# explode
\t@topo <molecular topology file>
\t@pos  <coordinates for the molecules>
\t@nsm  <number of solute molecules
\t@dist <distance to put between molecules>";

pub const KNOWN: &[&str] = &["topo", "pos", "nsm", "dist"];

pub const PROGRAM: Program = Program {
    name: "explode",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let nsm: usize = args.require("nsm")?;
    let dist: f64 = args.require("dist")?;
    let nsm3 = (nsm as f64).powf(1.0 / 3.0).ceil() as usize;
    let box_len = nsm3 as f64 * dist;
    let half = Vec3::splat(dist / 2.0);
    let topo = read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?;
    let pos_file = args.value("pos")?;
    let data = read_g96_labeled(pos_file).map_err(|e| format!("@pos: {e}"))?;
    let mols: Vec<std::ops::Range<usize>> = {
        let mut v = Vec::new();
        let mut start = 0;
        for &end in &topo.solute_molecules {
            v.push(start..end);
            start = end;
        }
        if v.is_empty() {
            v.push(0..topo.n_atoms);
        }
        v
    };
    if mols.len() < nsm {
        return Err(format!(
            "@nsm {nsm}: the topology has {} solute molecules",
            mols.len()
        ));
    }
    if data.atoms.len() < topo.n_atoms {
        return Err(format!(
            "{pos_file}: {} atoms, topology has {}",
            data.atoms.len(),
            topo.n_atoms
        ));
    }
    let mut placed = Vec::new();
    let mut count = 0usize;
    for i in 0..nsm3 {
        for j in 0..nsm3 {
            for k in 0..nsm3 {
                if count >= nsm {
                    count += 1;
                    continue;
                }
                let shift = Vec3::new(i as f64 * dist, j as f64 * dist, k as f64 * dist);
                let m = &mols[count];
                let mut com = Vec3::ZERO;
                let mut mass = 0.0;
                for a in m.clone() {
                    com += data.atoms[a].pos * topo.masses[a];
                    mass += topo.masses[a];
                }
                com /= mass;
                let rc = com - half;
                for a in m.clone() {
                    let mut atom = data.atoms[a].clone();
                    atom.pos = atom.pos - rc + shift;
                    placed.push(atom);
                }
                count += 1;
            }
        }
    }
    if count > nsm {
        eprintln!("NOTE: {} empty positions on the grid ({} grid points) since only {} molecules positioned", count - nsm, nsm3 * nsm3 * nsm3, nsm);
    }
    let title = format!("Explode : {nsm} molecules put at intermolecular distance {dist} nm\nTaken from: {pos_file}");
    write!(
        out,
        "{}",
        format_g96(&title, &placed, Some(Vec3::splat(box_len)))
    )
    .map_err(io_err)?;
    Ok(())
}
