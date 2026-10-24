//! ener — Recalculate interaction energies from a trajectory.
//!
//! Usage:
//!   ener @topo <top> @traj <trc> [@cutoff <nm>] [@epsrf <eps>] [@pbc]

use crate::programs::io_err;
use newen_core::{
    gather::gather_molecules,
    math::{Periodicity, Rectangular, Vec3},
};
use newen_forces::energy::{single_point_energy, EnergyParams};
use newen_io::args::{Arguments, Program, Status};
use newen_io::topology::{build_topology, read_topology_file};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# ener
\t@topo          <topology>
\t@traj          <trajectory>
\t[@cutoff       <nonbonded cutoff in nm (default: 1.4)>]
\t[@epsrf        <reaction-field dielectric (default: 78.5)>]
\t[@pbc          <gather the molecules before computing the energy>]";

pub const KNOWN: &[&str] = &["topo", "traj", "cutoff", "epsrf", "pbc"];

pub const PROGRAM: Program = Program {
    name: "ener",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    if argv.len() < 2 {
        return Err(USAGE.to_string());
    }
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let traj_file = args.value("traj")?;
    let do_pbc = args.has("pbc");
    let defaults = EnergyParams::default();
    let mut params = EnergyParams {
        cutoff: args.get("cutoff", defaults.cutoff)?,
        rf_epsilon: args.get("epsrf", defaults.rf_epsilon)?,
        ..defaults
    };

    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("Error reading topology: {e}"))?;
    let topo = build_topology(topo_data);
    params.atoms_per_solvent = topo.solvent_atom_template.len().max(3);

    let mut reader =
        TrajectoryReader::new(traj_file).map_err(|e| format!("Error opening trajectory: {e}"))?;

    writeln!(out, "# ener: topo={topo_file} traj={traj_file}").map_err(io_err)?;
    writeln!(
        out,
        "# cutoff={:.3} nm  epsrf={:.1}",
        params.cutoff, params.rf_epsilon
    )
    .map_err(io_err)?;
    writeln!(
        out,
        "# {:>12}  {:>14}  {:>14}  {:>14}  {:>14}",
        "time(ps)", "E_bond(kJ/mol)", "E_lj(kJ/mol)", "E_crf(kJ/mol)", "E_pot(kJ/mol)"
    )
    .map_err(io_err)?;

    while let Ok(Some(frame)) = reader.read_frame() {
        let mut pos = frame.positions.clone();

        if do_pbc && frame.box_dims.x > 0.0 {
            let periodicity = Periodicity::Rectangular(Rectangular::new(frame.box_dims));
            gather_molecules(&mut pos, &topo.molecules, &periodicity, Vec3::ZERO);
        }

        let e = single_point_energy(&topo, &pos, frame.box_dims, &params);

        writeln!(
            out,
            "{:14.6}  {:14.6}  {:14.6}  {:14.6}  {:14.6}",
            frame.time, e.bond, e.lj, e.crf, e.potential
        )
        .map_err(io_err)?;
    }
    Ok(())
}
