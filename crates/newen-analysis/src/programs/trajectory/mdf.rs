//! mdf — minimum distance from each centre atom to a set of atoms, per frame (gromos++ `mdf`).
//!
//! Usage: mdf @topo <top> @pbc <r|v> [cog|nog] [@time <t0 dt>] @centre <atoms> @with <atoms> @traj <trc…>
//!
//! One file `MIN_<mol>:<atom>.dat` per centre atom, as gromos++ writes them:
//! `time  min distance  # closest atom`.

use crate::args::Arguments;
use crate::pbc::Pbc;

use crate::property::ordered_atoms;
use crate::time::Time;
use newen_io::args::{Program, Status};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# mdf
\t@topo   <topology>
\t@pbc    <boundary type> [<gather method>]
\t@time   <time and dt>
\t@centre <atoms to take as centre>
\t@with   <atoms to calculate the distance for>
\t@traj   <trajectory files>";

pub const KNOWN: &[&str] = &["topo", "pbc", "time", "centre", "with", "nsm", "traj"];

pub const PROGRAM: Program = Program {
    name: "mdf",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// gromos++ `AtomSpecifier::toString(i)`: `mol:atom`, 1-based; solvent as `s:atom`.
fn atom_label(topo: &newen_core::Topology, atom: usize) -> String {
    let n_solute = topo.num_solute_atoms();
    if atom >= n_solute {
        return format!("s:{}", atom - n_solute + 1);
    }
    for (m, r) in topo.molecules.iter().enumerate() {
        if r.contains(&atom) {
            return format!("{}:{}", m + 1, atom - r.start + 1);
        }
    }
    format!("1:{}", atom + 1)
}

pub fn run(argv: &[String], _out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
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
    let mut time = Time::from_args(&args)?;
    let centre = ordered_atoms(&args.values("centre").join(";"), &topo)?;
    let with = ordered_atoms(&args.values("with").join(";"), &topo)?;
    if centre.is_empty() || with.is_empty() {
        return Err("@centre and @with must each select atoms".into());
    }
    let pbc = Pbc::from_args(&args)?;
    let mut outputs: Vec<String> = vec![String::new(); centre.len()];
    for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let t = time.next(frame.time);
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            for (ci, &c) in centre.iter().enumerate() {
                let mut min2 = f64::MAX;
                let mut minat = 0;
                for (wi, &w) in with.iter().enumerate() {
                    if w == c {
                        continue;
                    }
                    let d = periodicity.nearest_image(pos[c], pos[w]);
                    if d.length_squared() < min2 {
                        min2 = d.length_squared();
                        minat = wi;
                    }
                }
                // the stream keeps Time's fixed nine decimals for the distance
                outputs[ci].push_str(&format!(
                    "{}\t{:.9}\t# {}\n",
                    Time::fmt(t),
                    min2.sqrt(),
                    atom_label(&topo, with[minat])
                ));
            }
        }
    }
    for (ci, &c) in centre.iter().enumerate() {
        let name = format!("MIN_{}.dat", atom_label(&topo, c));
        std::fs::write(&name, &outputs[ci]).map_err(|e| format!("{name}: {e}"))?;
    }
    Ok(())
}
