//! rmsf — atom-positional root-mean-square fluctuation about the average structure
//! (gromos++ `rmsf`, `programs/rmsf.cc`).
//!
//! Every frame is superimposed on the reference using `@atomsfit` (defaulting to `@atomsrmsf`);
//! ⟨r⟩ and ⟨r²⟩ are then accumulated over `@atomsrmsf` and the fluctuation is
//! `sqrt(⟨r²⟩ − ⟨r⟩²)`.

use crate::args::Arguments;
use crate::fit::superimpose;
use crate::pbc::Pbc;
use crate::programs::io_err;
use newen_core::{math::Vec3, selection::AtomSelection, Topology};
use newen_io::args::{Program, Status};
use newen_io::coordinate::read_coordinates;
use newen_io::pdb::PDBAtom;
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::fmt::Write as _;
use std::io::Write;

pub const USAGE: &str = "# rmsf
\t@topo        <molecular topology file>
\t@pbc         <boundary type> [<gathermethod>]
\t@atomsrmsf   <atoms to consider for rmsf>
\t[@atomsfit   <atoms to consider for fit>]
\t[@ref        <reference coordinates (if absent, the first frame of @traj is reference)>]
\t[@outpdb     <write average structure in pdb format with rmsf values in the b-factor column>]
\t@traj        <trajectory files>";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "atomsrmsf",
    "atomsfit",
    "ref",
    "outpdb",
    "traj",
];

pub const PROGRAM: Program = Program {
    name: "rmsf",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// The atoms named by every value of `key`, unioned (gromos++ adds one specifier per value).
fn select(args: &Arguments, key: &str, topo: &Topology) -> Result<Vec<usize>, String> {
    let spec = args.values(key).join(";");
    AtomSelection::from_string(&spec, topo)
        .map(|s| s.indices().to_vec())
        .map_err(|e| format!("@{key} '{spec}': {e}"))
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;

    let mut topo =
        build_topology(read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?);
    let trajectories = args.values("traj");
    if trajectories.is_empty() {
        return Err("@traj: no trajectory files given".into());
    }
    if let Some(frame) = TrajectoryReader::new(&trajectories[0])
        .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
        .read_frame()
        .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
    {
        solvate_to_atoms(&mut topo, frame.positions.len()).map_err(|e| e.to_string())?;
    }

    let pbc = Pbc::from_args(&args)?;

    let (mut reference, ref_box) = if args.has("ref") {
        let path = args.value("ref")?;
        let conf = read_coordinates(path).map_err(|e| format!("@ref {path}: {e}"))?;
        (conf.positions, conf.box_dims)
    } else {
        let frame = TrajectoryReader::new(&trajectories[0])
            .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
            .read_frame()
            .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
            .ok_or_else(|| format!("@traj {}: no frames", trajectories[0]))?;
        (frame.positions, frame.box_dims)
    };
    let ref_periodicity = pbc.periodicity(ref_box);
    pbc.gather(&topo, &mut reference, &ref_periodicity);

    let rmsf_atoms = select(&args, "atomsrmsf", &topo)?;
    if rmsf_atoms.is_empty() {
        return Err("No atoms specified for RMSF calculation".into());
    }
    let fit_atoms = if args.has("atomsfit") {
        select(&args, "atomsfit", &topo)?
    } else {
        writeln!(out, "# @atomsrmsf atoms are taken for fit.").map_err(io_err)?;
        rmsf_atoms.clone()
    };

    let mut sum_pos = vec![Vec3::ZERO; rmsf_atoms.len()];
    let mut sum_sq = vec![0.0f64; rmsf_atoms.len()];
    let mut n_frames = 0usize;

    for traj in trajectories {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            if !fit_atoms.is_empty() {
                superimpose(&mut pos, &reference, &fit_atoms, None);
            }
            for (k, &a) in rmsf_atoms.iter().enumerate() {
                sum_pos[k] += pos[a];
                sum_sq[k] += pos[a].length_squared();
            }
            n_frames += 1;
        }
    }
    if n_frames == 0 {
        return Err("@traj: no frames".into());
    }

    let n = n_frames as f64;
    let mut text = String::from("#\n#  at          rmsf name\n");
    let mut average = Vec::with_capacity(rmsf_atoms.len());
    let mut fluctuation = Vec::with_capacity(rmsf_atoms.len());
    for (k, &a) in rmsf_atoms.iter().enumerate() {
        let mean = sum_pos[k] / n;
        let f = (sum_sq[k] / n - mean.length_squared()).max(0.0).sqrt();
        let _ = writeln!(
            text,
            "{:5}{:14.8} {:>4}",
            k + 1,
            f,
            topo.atom_name(a).unwrap_or("?")
        );
        average.push(mean);
        fluctuation.push(f);
    }
    write!(out, "{text}").map_err(io_err)?;

    if let Some(path) = args
        .values("outpdb")
        .first()
        .cloned()
        .or_else(|| args.has("outpdb").then(|| "aver.pdb".to_string()))
    {
        let mut pdb = String::from("REMARK Average structure\n");
        for (k, &a) in rmsf_atoms.iter().enumerate() {
            let atom = PDBAtom {
                serial: k + 1,
                name: topo.atom_name(a).unwrap_or("X").to_string(),
                resname: topo.residue_name(a).unwrap_or("UNK").to_string(),
                chain: " ".to_string(),
                resnum: topo.residue_nr(a).unwrap_or(1) as i32,
                coord: average[k],
                occupancy: 1.0,
                bfactor: fluctuation[k],
                element: String::new(),
            };
            pdb.push_str(&atom.to_pdb_line(false));
            pdb.push('\n');
        }
        pdb.push_str("END\n");
        std::fs::write(&path, pdb).map_err(|e| format!("@outpdb {path}: {e}"))?;
    }
    Ok(())
}
