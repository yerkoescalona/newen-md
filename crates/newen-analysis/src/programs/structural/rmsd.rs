//! rmsd — atom-positional root-mean-square deviation from a reference structure
//! (gromos++ `rmsd`, `programs/rmsd.cc`).
//!
//! The fit and the RMSD are taken over two independent atom sets: the frame is superimposed on
//! the reference using `@atomsfit` (defaulting to `@atomsrmsd`), and the deviation is then
//! measured over `@atomsrmsd`. Both sets are shifted by the *fit* set's centre of geometry
//! before the rotation, as gromos++ does (`rmsd.cc:337`).

use crate::args::Arguments;
use crate::fit::rmsd_after_fit;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::PropertyContainer;
use crate::time::Time;
use newen_core::{selection::AtomSelection, Topology};
use newen_io::args::{Program, Status};
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# rmsd
\t@topo       <molecular topology file>
\t@pbc        <boundary type> [<gathermethod>]
\t@time       <time and dt>
\t@atomsrmsd  <atoms to consider for rmsd>
\t[@atomsfit  <atoms to consider for fit>]
\t[@prop      <properties>]
\t[@ref       <reference coordinates (if absent, the first frame of @traj is reference)>]
\t[@printatoms  <print list of selected atoms>]
\t@traj       <trajectory files>";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "time",
    "atomsrmsd",
    "atomsfit",
    "prop",
    "ref",
    "traj",
    "printatoms",
    "reftopo",
    "refpbc",
];

pub const PROGRAM: Program = Program {
    name: "rmsd",
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
    for unported in ["reftopo", "refpbc"] {
        if args.has(unported) {
            return Err(format!(
                "@{unported} is not ported — the reference must use @topo"
            ));
        }
    }
    if !args.has("atomsrmsd") && !args.has("prop") {
        return Err("No rmsd atoms or property specified!".into());
    }

    let mut time = Time::from_args(&args)?;
    let mut topo =
        build_topology(read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?);

    let trajectories = args.values("traj");
    if trajectories.is_empty() {
        return Err("@traj: no trajectory files given".into());
    }

    // gromos++ reads every atom of every frame (`select("ALL")`).
    if let Some(frame) = TrajectoryReader::new(&trajectories[0])
        .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
        .read_frame()
        .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
    {
        solvate_to_atoms(&mut topo, frame.positions.len()).map_err(|e| e.to_string())?;
    }

    let pbc = Pbc::from_args(&args)?;

    // Reference: @ref, else the first frame of the first trajectory.
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

    let rmsd_atoms = if args.has("atomsrmsd") {
        select(&args, "atomsrmsd", &topo)?
    } else {
        Vec::new()
    };
    let fit_atoms = if args.has("atomsfit") {
        select(&args, "atomsfit", &topo)?
    } else {
        if args.has("atomsrmsd") {
            writeln!(out, "# @atomsrmsd atoms are taken for fit.").map_err(io_err)?;
        }
        rmsd_atoms.clone()
    };
    if args.has("printatoms") {
        for (i, a) in fit_atoms.iter().enumerate() {
            let r = rmsd_atoms.get(i).copied();
            eprintln!(
                "# fit {:>6} {:<5} | rmsd {}",
                a + 1,
                topo.atom_name(*a).unwrap_or("?"),
                r.map_or("-".to_string(), |r| format!(
                    "{:>6} {}",
                    r + 1,
                    topo.atom_name(r).unwrap_or("?")
                ))
            );
        }
    }

    let mut props = if args.has("prop") {
        Some((
            PropertyContainer::parse(args.values("prop"), &topo)?,
            PropertyContainer::parse(args.values("prop"), &topo)?,
        ))
    } else {
        None
    };
    let ref_prop = props
        .as_mut()
        .map(|(r, _)| r.calc(&reference, &ref_periodicity));

    let mut text = String::new();
    for traj in trajectories {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let t = time.next(frame.time);
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);

            let r_prop = match (&mut props, &ref_prop) {
                (Some((_, sys)), Some(refv)) => {
                    let v = sys.calc(&pos, &periodicity);
                    let n = refv.len().max(1) as f64;
                    let sum: f64 = refv.iter().zip(&v).map(|(a, b)| (a - b) * (a - b)).sum();
                    Some((sum / n).sqrt())
                },
                _ => None,
            };
            // rmsd_after_fit superimposes on `fit_atoms` and measures over `rmsd_atoms`,
            // which is gromos++'s FastRotationalFit::fit + ::rmsd pair.
            let r = (!rmsd_atoms.is_empty())
                .then(|| rmsd_after_fit(&mut pos, &reference, &fit_atoms, &rmsd_atoms, None));

            text.push_str(&Time::fmt(t));
            match (r, r_prop) {
                (Some(r), Some(rp)) => text.push_str(&format!("{r:15.9}{rp:15.9}")),
                (Some(r), None) => text.push_str(&format!("{r:10.5}")),
                (None, Some(rp)) => text.push_str(&format!("{rp:15.9}")),
                (None, None) => {},
            }
            text.push('\n');
        }
    }
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
