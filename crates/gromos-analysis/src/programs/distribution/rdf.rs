//! rdf — radial distribution function of `@with` atoms around `@centre` atoms over a trajectory
//! (gromos++ `rdf`, `programs/rdf.cc`). The physics is [`crate::rdf::Rdf`].
//!
//! Usage: rdf @topo <top> @pbc <r|v> [cog|nog] @centre <spec…> @with <spec…> @cut <r>
//!        @grid <n> [@nointra] @traj <trc…>

use crate::args::Arguments;
use crate::distribution::cpp_g;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::rdf::Rdf;
use gromos_core::{selection::AtomSelection, Topology};
use gromos_io::args::{Program, Status};
use gromos_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use gromos_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# rdf
\t@topo   <molecular topology file>
\t@pbc        <boundary type> [<gather method>]
\t@centre <atoms to take as centre>
\t@with   <atoms to calculate distances for>
\t@cut    <maximum distance>
\t@grid   <number of points>
\t[@nointra   <skip intramolecular atoms>]
\t[@doDCF     <compute dipole-dipole correlation> <norm>]
\t@traj   <trajectory files>";

pub const KNOWN: &[&str] = &[
    "topo", "centre", "with", "cut", "grid", "nointra", "traj", "pbc", "doDCF",
];

pub const PROGRAM: Program = Program {
    name: "rdf",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// The atoms named by every value of `key`, unioned (gromos++ adds one specifier per value).
fn select(args: &Arguments, key: &str, topo: &Topology) -> Result<Vec<usize>, String> {
    let spec = args.values(key).join(";");
    let atoms = AtomSelection::from_string(&spec, topo)
        .map(|s| s.indices().to_vec())
        .map_err(|e| format!("@{key} '{spec}': {e}"))?;
    if atoms.is_empty() {
        return Err(format!("No atoms specified for {key} atoms!"));
    }
    Ok(atoms)
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    if args.has("doDCF") {
        return Err("@doDCF is not supported".into());
    }
    let mut topo =
        build_topology(read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?);
    let trajectories = args.values("traj");
    if trajectories.is_empty() {
        return Err("@traj: no trajectory files given".into());
    }
    // gromos++ reads every atom of the frames (`select("ALL")`): solvate to the first frame
    if let Some(frame) = TrajectoryReader::new(&trajectories[0])
        .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
        .read_frame()
        .map_err(|e| format!("@traj {}: {e}", trajectories[0]))?
    {
        solvate_to_atoms(&mut topo, frame.positions.len()).map_err(|e| e.to_string())?;
    }
    let pbc = Pbc::from_args(&args)?;
    let centres = select(&args, "centre", &topo)?;
    let withs = select(&args, "with", &topo)?;
    let cut: f64 = args.require("cut")?;
    let grid: usize = args.require("grid")?;
    let nointra = args.has("nointra");

    // gromos++ molecule numbers: solute molecules by index, all solvent as one molecule (−1)
    let n_solute = topo.num_solute_atoms();
    let mut molecule_of = vec![-1i64; topo.num_atoms()];
    for (m, range) in topo.molecules.iter().enumerate() {
        if range.end <= n_solute {
            molecule_of[range.clone()].fill(m as i64);
        }
    }

    let mut rdf = Rdf::new(cut, grid)?;
    for traj in trajectories {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            if frame.positions.len() != molecule_of.len() {
                return Err(format!(
                    "{traj}: a frame has {} atoms, the topology {}",
                    frame.positions.len(),
                    molecule_of.len()
                ));
            }
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            let volume = frame.box_dims.x * frame.box_dims.y * frame.box_dims.z;
            rdf.add_frame(
                &pos,
                &centres,
                &withs,
                &molecule_of,
                &periodicity,
                volume,
                nointra,
            );
        }
    }
    for (k, g) in rdf.g().into_iter().enumerate() {
        writeln!(out, "{:>15}{:>15}", cpp_g(rdf.r(k), 9), cpp_g(g, 9)).map_err(io_err)?;
    }
    Ok(())
}
