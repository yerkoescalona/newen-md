//! jval — ³J-coupling constants from dihedral angles through the Karplus relation
//! (gromos++ `jval`).
//!
//! Usage: jval @topo <top> @pbc <r|v> [cog|nog] @jval <JVALRESSPEC file> [@timeseries] [@rmsd]
//!        [@time <t0 dt>] [@timespec ALL|EVERY|SPEC] [@timepts <n…>] [@avg] @traj <trc…>
//!
//! The specification lines are `i j k l fdum J0 delta A B C` (global atom numbers, 1-based);
//! J = A cos²(φ+δ) + B cos(φ+δ) + C.

use crate::args::Arguments;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::Property;
use crate::time::{Time, TimeSpec};
use newen_core::stat::Stat;
use newen_io::args::{Program, Status};
use newen_io::jvalue::read_jvalrespec;
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# jval
\t@topo         <molecular topology file>
\t@pbc          <boundary type [<gather method>]>
\t@jval         <j-value specification file>
\t[@timeseries  <write time-series of j-values>
\t[@rmsd        <write the rmsd over all j-values as a time-series>
\t[@time        <time and dt> (optional and only if time-series)]
\t[@timespec    <timepoints at which to compute the j-values: ALL (default), EVERY or SPEC (if time-series)]
\t[@timepts     <timepoints at which to compute the j-values (if time-series and timespec EVERY or SPEC)]
\t@traj         <position trajectory file(s)>";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "jval",
    "timeseries",
    "rmsd",
    "time",
    "timespec",
    "timepts",
    "avg",
    "traj",
];

pub const PROGRAM: Program = Program {
    name: "jval",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
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
    let pbc = Pbc::from_args(&args)?;
    let kps = read_jvalrespec(args.value("jval")?).map_err(|e| e.to_string())?;
    let mol_of = |a: usize| {
        topo.molecules
            .iter()
            .position(|r| r.contains(&a))
            .unwrap_or(0)
    };
    let mut props: Vec<Property> = Vec::new();
    for kp in &kps {
        let m = mol_of(kp.atoms[0]);
        let off = topo.molecules.get(m).map(|r| r.start).unwrap_or(0);
        let spec = format!(
            "t%{}:{},{},{},{}",
            m + 1,
            kp.atoms[0] - off + 1,
            kp.atoms[1] - off + 1,
            kp.atoms[2] - off + 1,
            kp.atoms[3] - off + 1
        );
        props.push(Property::parse(&spec, &topo)?);
    }
    let mut time = Time::from_args(&args)?;
    let time_series_given = args.count("time") >= 0;
    let mut timespec = TimeSpec::from_args(&args)?;
    let jv_ts = args.count("timeseries") >= 0;
    let rmsd_ts = args.count("rmsd") >= 0;
    if rmsd_ts && jv_ts {
        return Err("cannot write a time-series of the j-values and the rmsd: use either @rmsd or @timeseries".into());
    }
    let avg_ts = args.count("avg") >= 0;
    let mut text = String::new();
    if rmsd_ts {
        text.push_str(&format!("{:>8}{:>15}\n", "#time", "rmsd"));
    }
    let mut stat_j: Vec<Stat> = (0..kps.len()).map(|_| Stat::new()).collect();
    let mut stat_phi: Vec<Stat> = (0..kps.len()).map(|_| Stat::new()).collect();
    'traj: for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let t = time.next(frame.time);
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            let compute = timespec.select();
            if compute {
                let values: Vec<f64> = props
                    .iter_mut()
                    .map(|p| p.calc(&pos, &periodicity))
                    .collect();
                if time_series_given && !avg_ts {
                    text.push_str(&format!("#\n#\n# TIME\t{}\n", Time::fmt(t)));
                    text.push_str(&format!("#{:>4}{:>12}\n", "num", "j"));
                }
                if jv_ts && !rmsd_ts {
                    text.push_str(&format!("#\n#\n# TIME\t{}\n", Time::fmt(t)));
                    text.push_str(&format!("#{:>4}{:>16}\n", "num", "j-value"));
                }
                let mut rmsd = 0.0;
                for (i, kp) in kps.iter().enumerate() {
                    let cosphi = ((values[i] + kp.delta).to_radians()).cos();
                    let j = kp.a * cosphi * cosphi + kp.b * cosphi + kp.c;
                    stat_j[i].add(j);
                    stat_phi[i].add(values[i]);
                    if jv_ts && !rmsd_ts {
                        // the stream still carries Time's fixed nine decimals
                        text.push_str(&format!("{:>5}{:>16.9}\n", i + 1, j));
                    }
                    rmsd += (j - kp.j0) * (j - kp.j0);
                }
                rmsd /= kps.len() as f64;
                if rmsd_ts {
                    text.push_str(&format!("{:>8}{:>15.9}\n", Time::fmt(t), rmsd.sqrt()));
                }
            }
            if timespec.done {
                break 'traj;
            }
        }
    }
    text.push_str("#\n#\n# JVALUE AVERAGES\n#\n# (print with a2ps --landscape --columns=1 --font-size=9)\n#\n");
    text.push_str(&format!(
        "#{:>4}{:>4}{:>9}{:>15}{:>21}{:>11}{:>5}{:>5}{:>9}{:>7}{:>10}{:>9}{:>9}{:>8}{:>10}\n",
        "num",
        "mol",
        "residue",
        "atom names",
        "atom numbers",
        "A",
        "B",
        "C",
        "delta",
        "J0",
        "phi ave",
        "phi rmsd",
        "J ave",
        "J rmsd",
        "|Jave-J0|"
    ));
    let (mut sum, mut abssum, mut ssum) = (0.0, 0.0, 0.0);
    for (i, kp) in kps.iter().enumerate() {
        let m = mol_of(kp.atoms[0]);
        let off = topo.molecules.get(m).map(|r| r.start).unwrap_or(0);
        let j_atom = kp.atoms[1];
        let mut line = format!(
            "{:>5}{:>4}{:>4}{:>5}{:>5}{:>5}{:>5}{:>5}{:>7}{:>5}{:>5}{:>5}",
            i + 1,
            m + 1,
            topo.residue_nr(j_atom).unwrap_or(0),
            topo.residue_name(j_atom).unwrap_or(""),
            topo.atom_name(kp.atoms[0]).unwrap_or(""),
            topo.atom_name(kp.atoms[1]).unwrap_or(""),
            topo.atom_name(kp.atoms[2]).unwrap_or(""),
            topo.atom_name(kp.atoms[3]).unwrap_or(""),
            kp.atoms[0] - off + 1,
            kp.atoms[1] - off + 1,
            kp.atoms[2] - off + 1,
            kp.atoms[3] - off + 1
        );
        line.push_str(&format!("{:>7.1}{:>5.1}{:>5.1}", kp.a, kp.b, kp.c));
        line.push_str(&format!("{:>7.0}", kp.delta));
        line.push_str(&format!("{:>7.2}", kp.j0));
        line.push_str(&format!("{:>10.1}", stat_phi[i].ave()));
        line.push_str(&format!(
            "{:>9.2}{:>9.2}",
            stat_phi[i].rmsd(),
            stat_j[i].ave()
        ));
        line.push_str(&format!("{:>8.3}", stat_j[i].rmsd()));
        line.push_str(&format!("{:>10.1}\n", (stat_j[i].ave() - kp.j0).abs()));
        text.push_str(&line);
        sum += stat_j[i].ave() - kp.j0;
        abssum += (stat_j[i].ave() - kp.j0).abs();
        ssum += (kp.j0 - stat_j[i].ave()).powi(2);
    }
    let n = kps.len() as f64;
    // the stream is left at fixed precision 1 by the table
    text.push_str(&format!("\n#{:>30}{:.1}\n", "average deviation ", sum / n));
    text.push_str(&format!(
        "#{:>30}{:.1}\n",
        "average absolute deviation ",
        abssum / n
    ));
    text.push_str(&format!(
        "#{:>30}{:.1}\n",
        "root-mean-square deviation ",
        (ssum / n).sqrt()
    ));
    text.push_str(&format!(
        "#{:>30}{:.1}\n",
        "rmsd over deviations ",
        ((ssum - sum * sum / n) / n).sqrt()
    ));
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
