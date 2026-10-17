//! bilayer_oparam — order parameters of lipid tails from the local molecular frame of every
//! selected atom and its two bonded neighbours (gromos++ `bilayer_oparam`).
//!
//! Usage: bilayer_oparam @topo <top> @pbc <r|v> [cog|nog] @atoms <spec> @refvec <x y z> @traj <trc…>
//!
//! `@atoms` selects the atoms of the first lipid; the same atom numbers are taken in every
//! molecule up to the last selected one. First and last atoms of a molecule have no two
//! neighbours and are rejected, as in gromos++.

use crate::args::Arguments;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::ordered_atoms;
use newen_core::Vec3;
use newen_io::args::{Program, Status};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# bilayer_oparam
\t@topo           <molecular topology file>
\t@pbc            <boundary type> [<gathermethod>]
\t[@time           <time and dt>]
\t@atoms          <atomspecifier>
\t@refvec        <x y z>]
\t@traj           <trajectory files>";

pub const KNOWN: &[&str] = &["topo", "pbc", "atoms", "refvec", "time", "traj"];

pub const PROGRAM: Program = Program {
    name: "bilayer_oparam",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// C++ prints a NaN as `nan` or `-nan` by its sign bit; the fixed/precision format otherwise.
fn fixed_or_nan(x: f64, decimals: usize) -> String {
    if x.is_nan() {
        if x.is_sign_negative() {
            "-nan".to_string()
        } else {
            "nan".to_string()
        }
    } else {
        format!("{:.*}", decimals, x)
    }
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let rv = args.values("refvec");
    let mut refvec = Vec3::ZERO;
    for (k, v) in rv.iter().take(3).enumerate() {
        let x: f64 = v.parse().map_err(|_| format!("@refvec: bad number {v}"))?;
        match k {
            0 => refvec.x = x,
            1 => refvec.y = x,
            _ => refvec.z = x,
        }
    }
    let refvec = refvec.normalize();
    let mut topo =
        build_topology(read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?);
    // gromos++ reads every atom of the frames (`select("ALL")`): solvate to the first frame
    if let Some(first) = args.values("traj").first() {
        if let Ok(Some(frame)) = TrajectoryReader::new(first)
            .map_err(|e| format!("@traj {first}: {e}"))?
            .read_frame()
        {
            solvate_to_atoms(&mut topo, frame.positions.len()).map_err(|e| e.to_string())?;
        }
    }
    let bilayer_atoms = ordered_atoms(&args.values("atoms").join(";"), &topo)?;
    if bilayer_atoms.is_empty() {
        return Err("@atoms: no atoms selected".into());
    }
    let mols = &topo.molecules;
    let mol_of = |a: usize| mols.iter().position(|r| r.contains(&a)).unwrap_or(0);
    let moln = mol_of(*bilayer_atoms.last().unwrap()) + 1;
    let num_atperlip = bilayer_atoms.len() / moln;
    let last_mol = &mols[moln - 1];
    let mut atoms = Vec::new();
    for &a in bilayer_atoms.iter().take(num_atperlip) {
        let local = a - mols[mol_of(a)].start;
        if local == 0 {
            return Err("cannot calculate the oparam for the first atom!".into());
        }
        if local >= last_mol.len() - 1 {
            return Err("cannot calculate the oparam for the last atom!".into());
        }
        atoms.push(local);
    }
    if atoms.is_empty() {
        return Err("at least one atom needs to be defined!".into());
    }
    // gromos++ Neighbours: the bonded partners of an atom in the order of the BOND blocks
    let first_mol = mols[mol_of(bilayer_atoms[0])].clone();
    let mut bonds: Vec<(usize, usize)> = topo
        .all_bonds_global()
        .filter(|b| first_mol.contains(&b.i))
        .map(|b| (b.i - first_mol.start, b.j - first_mol.start))
        .collect();
    bonds.sort(); // gromos++ keeps bonds in a std::set ordered by (i, j)
    let neighbours = |a: usize| -> Vec<usize> {
        let mut n = Vec::new();
        for &(i, j) in &bonds {
            if i == a {
                n.push(j);
            }
            if j == a {
                n.push(i);
            }
        }
        n
    };
    let mut at: Vec<usize> = Vec::new();
    for &a in &atoms {
        let n = neighbours(a);
        let first = *n
            .first()
            .ok_or_else(|| format!("atom {} has no bonded neighbour", a + 1))?;
        at.push(first);
        at.push(a);
        at.push(*n.last().unwrap());
    }
    let pbc = Pbc::from_args(&args)?;
    let n_sel = at.len() / 3;
    let mut avcos2 = vec![Vec3::ZERO; n_sel];
    let mut num_frames = 0usize;
    for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            num_frames += 1;
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            for m in 0..moln {
                let off = mols[m].start;
                for j in 0..n_sel {
                    let (p0, p1, p2) = (
                        pos[off + at[3 * j]],
                        pos[off + at[3 * j + 1]],
                        pos[off + at[3 * j + 2]],
                    );
                    let z = (p2 - p0).normalize();
                    let y = p1 - p0;
                    let y = (y - (p1 - p0).dot(z) * z).normalize();
                    let x = z.cross(y).normalize();
                    let c = Vec3::new(refvec.dot(x), refvec.dot(y), refvec.dot(z));
                    avcos2[j] += c * c;
                }
            }
        }
    }
    for v in avcos2.iter_mut() {
        *v /= (num_frames * moln) as f64;
    }
    let s: Vec<Vec3> = avcos2
        .iter()
        .map(|c| {
            Vec3::new(
                (3.0 * c.x - 1.0) / 2.0,
                (3.0 * c.y - 1.0) / 2.0,
                (3.0 * c.z - 1.0) / 2.0,
            )
        })
        .collect();
    writeln!(
        out,
        "{:>4}{:>8}{:>8}{:>8}{:>8}{:>8}",
        "Atom", "SX", "SY", "SZ", "SCHOP", "|SCD|"
    )
    .map_err(io_err)?;
    writeln!(out).map_err(io_err)?;
    for j in 0..n_sel {
        let c = avcos2[j];
        let scdop = -(c.x + (c.y - 1.0) / 2.0);
        let scdeb = (2.0 * s[j].x / 3.0) + s[j].y / 3.0;
        writeln!(
            out,
            "{:>4}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            at[3 * j + 1] + 1,
            fixed_or_nan(s[j].x, 4),
            fixed_or_nan(s[j].y, 4),
            fixed_or_nan(s[j].z, 4),
            fixed_or_nan(scdop, 4),
            fixed_or_nan(s[j].y.abs(), 4),
            fixed_or_nan(scdeb, 4)
        )
        .map_err(io_err)?;
    }
    Ok(())
}
