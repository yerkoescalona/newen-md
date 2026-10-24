//! jepot — the J-value local-elevation restraint potential (gromos++ `jepot`).
//!
//! Usage: jepot @jval <JVALRESSPEC file> @K <force constant> @ngrid <bins> [@angles ALL|CURR]
//!        [@fin <configuration with a JVALUERESEPS block>] [@time <t0 dt>]
//!        [@timespec ALL|EVERY|SPEC] [@timepts <n…>] [@topo <top>] [@pbc <r>]
//!        [@postraj <trc…>] [@restraj <trs…>]
//!
//! V(φ) = K Σ_bins w ε_bin exp(−Δφ²/2δ²), δ the bin width, with the elevation weights ε from a
//! restraint trajectory (`@restraj`, every selected frame) or from a final configuration
//! (`@fin`). `ALL` scans φ over 0…360° in 1° steps; `CURR` evaluates V at the dihedral of the
//! corresponding position-trajectory frame (a single restraint that is a dihedral of the
//! topology).

use crate::args::Arguments;
use crate::distribution::cpp_g;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::Property;
use crate::time::{Time, TimeSpec};
use newen_io::args::{Program, Status};
use newen_io::jvalue::{read_jvalrespec, read_jvalue_eps_block, read_restraint_trajectory};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::io::Write;

pub const USAGE: &str = "# jepot
\t@jval        <jvalue restraint specifications>
\t@K           <force constant>
\t@ngrid       <number of grid points>
\t[@angles     <angle values over which to compute the LE potential: ALL (default) or CURR>]
\t[@fin         <file containing final coordinates (if not time-series)>]
\t[@time       <time dt (optional and only if time-series)>]
\t[@timespec   <timepoints at which to compute the LE potential: ALL (default), EVERY or SPEC (if time-series)>]
\t[@timepts    <timepoints at which to compute the LE potential (if time-series and timespec EVERY or SPEC)>]
\t[@topo       <molecular topology file (if CURR)>]
\t[@pbc        <boundary type (if CURR)>]
\t[@postraj    <position trajectory files (if CURR)>]
\t[@restraj    <restraint trajectory files (if time-series)>]";

pub const KNOWN: &[&str] = &[
    "jval", "K", "ngrid", "angles", "fin", "time", "timespec", "timepts", "topo", "pbc", "postraj",
    "restraj",
];

pub const PROGRAM: Program = Program {
    name: "jepot",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// the local-elevation potential at φ (rad)
fn le_potential(phi: f64, eps: &[f64], k: f64, weight: f64, ngrid: usize) -> f64 {
    let binwidth = 2.0 * PI / ngrid as f64;
    let mut v = 0.0;
    for (bin, e) in eps.iter().enumerate().take(ngrid) {
        let phi0 = (bin as f64 + 0.5) * binwidth;
        let mut phi_bin = phi;
        while phi_bin < phi0 - PI {
            phi_bin += 2.0 * PI;
        }
        while phi_bin > phi0 + PI {
            phi_bin -= 2.0 * PI;
        }
        let delta_phi = phi_bin - phi0;
        v += k * weight * e * (-delta_phi * delta_phi / (2.0 * binwidth * binwidth)).exp();
    }
    v
}

/// gromos++'s scan `for (phi = 0; phi < 2π; phi += π/180)`; `fixed` once a time has been printed
/// (the stream then stays in fixed notation with nine decimals)
fn scan(eps: &[f64], k: f64, weight: f64, ngrid: usize, fixed: bool) -> String {
    let mut out = String::new();
    let mut phi = 0.0;
    while phi < 2.0 * PI {
        let v = le_potential(phi, eps, k, weight, ngrid);
        let deg = 180.0 * phi / PI;
        if fixed {
            out.push_str(&format!("{deg:>18.9}{v:>18.9}\n"));
        } else {
            out.push_str(&format!("{:>18}{:>18}\n", cpp_g(deg, 6), cpp_g(v, 6)));
        }
        phi += PI / 180.0;
    }
    out
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let angles = if args.count("angles") > 0 {
        args.value("angles")?.to_string()
    } else {
        "ALL".to_string()
    };
    if angles != "ALL" && angles != "CURR" {
        return Err(format!("angle specification {angles} unknown.\n"));
    }
    let kps = read_jvalrespec(args.value("jval")?).map_err(|e| e.to_string())?;
    if angles == "CURR" && kps.len() > 1 {
        return Err(
            "To compute current-angle jepot, jval specification file should only contain one restraint\n".into(),
        );
    }
    let k: f64 = args.require("K")?;
    let ngrid: usize = args.require("ngrid")?;
    if ngrid == 0 {
        return Err("@ngrid must be positive".into());
    }
    let mut time = Time::from_args(&args)?;
    let mut ts_j = TimeSpec::from_args(&args)?;
    let je_ts = args.count("restraj") > 0;
    if !je_ts {
        if args.count("fin") < 0 {
            return Err("either specify @restraj or @fin:".into());
        } else if angles == "CURR" {
            return Err("CURR does not work with @fin: use @postraj and @restraj instead and use @timespec and @timepts to select the last frame".into());
        }
    }
    let mut text = String::new();
    let mut dihedral_angles: Vec<f64> = Vec::new();
    if angles == "CURR" {
        if !(args.count("topo") > 0 && args.count("pbc") > 0 && args.count("postraj") > 0) {
            return Err("To compute current-angle jepot, you must specify a topology file and the periodic boundary conditions\n".into());
        }
        let parsed = read_topology_file(args.value("topo")?).map_err(|e| e.to_string())?;
        let mut topo = build_topology(parsed);
        if let Some(first) = args.values("postraj").first() {
            if let Ok(Some(frame)) = TrajectoryReader::new(first)
                .map_err(|e| format!("@postraj {first}: {e}"))?
                .read_frame()
            {
                solvate_to_atoms(&mut topo, frame.positions.len()).map_err(|e| e.to_string())?;
            }
        }
        let pbc = Pbc::from_args(&args)?;
        let spec = format!(
            "t%1:{}",
            kps.iter()
                .map(|kp| kp
                    .atoms
                    .iter()
                    .map(|a| (a + 1).to_string())
                    .collect::<Vec<_>>()
                    .join(","))
                .collect::<Vec<_>>()
                .join("")
        );
        let mut prop = Property::parse(&spec, &topo)?;
        if prop.type_name() != "Torsion" {
            return Err("Only dihedral (torsion) properties allowed".into());
        }
        // gromos++ getTopologyType: the dihedral has to exist in the topology
        let dihedrals: BTreeSet<(usize, usize, usize, usize)> = topo
            .all_proper_dihedrals_global()
            .map(|d| (d.i, d.j, d.k, d.l))
            .chain(
                topo.all_improper_dihedrals_global()
                    .map(|d| (d.i, d.j, d.k, d.l)),
            )
            .collect();
        let a = &prop.atoms;
        if !(dihedrals.contains(&(a[0], a[1], a[2], a[3]))
            || dihedrals.contains(&(a[3], a[2], a[1], a[0])))
        {
            return Err(format!("Property {} not found", prop.title(&topo)));
        }
        let mut ts_a = TimeSpec::from_args(&args)?;
        'post: for traj in args.values("postraj") {
            let mut reader =
                TrajectoryReader::new(traj).map_err(|e| format!("@postraj {traj}: {e}"))?;
            while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
                let periodicity = pbc.periodicity(frame.box_dims);
                let mut pos = frame.positions;
                pbc.gather(&topo, &mut pos, &periodicity);
                if ts_a.select() {
                    dihedral_angles.push(prop.calc(&pos, &periodicity));
                }
                if ts_a.done {
                    break 'post;
                }
            }
        }
        text.push_str(&format!("{:>15}{:>18}{:>18}\n", "#time", "phi", "V_J_LE"));
    }
    if je_ts {
        'res: for path in args.values("restraj") {
            let frames = read_restraint_trajectory(path).map_err(|e| e.to_string())?;
            for frame in frames {
                let t = time.next(frame.time);
                if ts_j.select() {
                    for e in &frame.jvalue_eps {
                        for kp in &kps {
                            let atoms = [
                                kp.atoms[0] + 1,
                                kp.atoms[1] + 1,
                                kp.atoms[2] + 1,
                                kp.atoms[3] + 1,
                            ];
                            if e.atoms != atoms {
                                continue;
                            }
                            if e.epsilon.len() < ngrid {
                                return Err(format!(
                                    "{path}: a JVALUERESEPS entry has {} values, @ngrid is {ngrid}",
                                    e.epsilon.len()
                                ));
                            }
                            if angles == "CURR" {
                                let deg = dihedral_angles
                                    .get(ts_j.written - 1)
                                    .ok_or("@postraj has fewer selected frames than @restraj")?;
                                let phi = (deg / 180.0) * PI;
                                let v = le_potential(phi, &e.epsilon, k, kp.weight, ngrid);
                                text.push_str(&format!(
                                    "{:>15.9}{:>18.9}{:>18.9}\n",
                                    t,
                                    180.0 * phi / PI,
                                    v
                                ));
                            } else {
                                text.push_str(&format!(
                                    "#time:\t{}\t\tatoms:\t{}\t{}\t{}\t{}\n\n",
                                    Time::fmt(t),
                                    e.atoms[0],
                                    e.atoms[1],
                                    e.atoms[2],
                                    e.atoms[3]
                                ));
                                text.push_str(&scan(&e.epsilon, k, kp.weight, ngrid, true));
                                text.push_str("\n\n");
                            }
                        }
                    }
                }
                if ts_j.done {
                    break 'res;
                }
            }
        }
    } else {
        let path = args.value("fin")?;
        let eps = read_jvalue_eps_block(path).map_err(|e| e.to_string())?;
        if eps.len() != kps.len() {
            return Err(format!(
                "J-value specification file has different number of restraints to final jepot file {path}"
            ));
        }
        for (kp, e) in kps.iter().zip(&eps) {
            if e.len() < ngrid {
                return Err(format!(
                    "{path}: a JVALUERESEPS line has {} values, @ngrid is {ngrid}",
                    e.len()
                ));
            }
            text.push_str(&format!(
                "\n\n# atoms: {} {} {} {} \n",
                kp.atoms[0] + 1,
                kp.atoms[1] + 1,
                kp.atoms[2] + 1,
                kp.atoms[3] + 1
            ));
            text.push_str(&scan(e, k, kp.weight, ngrid, false));
        }
    }
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
