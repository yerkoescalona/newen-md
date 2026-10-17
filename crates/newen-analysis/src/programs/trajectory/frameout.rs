//! frameout — Write individual snapshots from a trajectory.
//!
//! Full port of GROMOS frameout. All I/O goes through newen-io.
//!
//! Usage:
//!   frameout @topo <top> @traj <trc>
//!            [@pbc <r|t|v> [gather]]
//!            [@spec <ALL|EVERY <n>|SPEC>] [@frames <list>]
//!            [@include <SOLUTE|SOLVENT|ALL>]
//!            [@ref <cnf>] [@atomsfit <spec>]
//!            [@out <cnf|pdb|trc>] [@name <prefix>] [@single]
//!            [@notimeblock] [@time <t_start> [t_end] [dt]]

use crate::fit::superimpose;
use newen_core::{
    gather::gather_molecules,
    math::{Mat3, Periodicity, Rectangular, Triclinic, Vacuum, Vec3},
    selection::AtomSelection,
};
use newen_io::args::{Arguments, Program, Status};
use newen_io::{
    coordinate::read_coordinates,
    g96::write_g96,
    pdb::write_pdb_positions,
    topology::{build_topology, read_topology_file},
    trajectory::{TrajectoryReader, TrajectoryWriter},
};
use std::io::Write;

pub const USAGE: &str = "# frameout
\t@topo          <topology>
\t@traj          <trajectory>
\t[@pbc          <r|t|v> [cog|nog]]
\t[@spec         <ALL|EVERY <n>|SPEC>]
\t[@frames       <frame list, e.g. 1,4,7-9>]
\t[@include      <SOLUTE|SOLVENT|ALL>]
\t[@ref          <reference coordinates for the fit>]
\t[@atomsfit     <atoms used for the rotational fit>]
\t[@out          <cnf|pdb|trc>]
\t[@name         <output file prefix>]
\t[@single       <one file instead of one per frame>]
\t[@notimeblock  <omit the TIMESTEP block>]
\t[@time         <t_start> [t_end] [dt]]";

pub const KNOWN: &[&str] = &[
    "topo",
    "traj",
    "ref",
    "atomsfit",
    "out",
    "name",
    "include",
    "single",
    "notimeblock",
    "pbc",
    "spec",
    "frames",
    "time",
];

pub const PROGRAM: Program = Program {
    name: "frameout",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

enum FrameSpec {
    All,
    Every(usize),
    Spec(Vec<usize>),
}

impl FrameSpec {
    fn matches(&self, idx: usize) -> bool {
        match self {
            FrameSpec::All => true,
            FrameSpec::Every(n) => idx.is_multiple_of(*n),
            FrameSpec::Spec(list) => list.contains(&idx),
        }
    }
}

fn parse_index_list(s: &str) -> Vec<usize> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let p = part.trim();
        if let Some((a, b)) = p.split_once('-') {
            if let (Ok(a), Ok(b)) = (a.parse::<usize>(), b.parse::<usize>()) {
                out.extend(a..=b);
                continue;
            }
        }
        if let Ok(n) = p.parse::<usize>() {
            out.push(n);
        }
    }
    out
}

pub fn run(argv: &[String], _out: &mut dyn Write) -> Result<(), String> {
    if argv.len() < 2 {
        return Err(USAGE.to_string());
    }
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.values("topo").first().cloned();
    let traj_file = args.value("traj")?;
    let ref_file = args.values("ref").first().cloned();
    let fit_spec = args.values("atomsfit").first().cloned();
    let outfmt = args.values("out").first().map_or("cnf", |s| s.as_str());
    let name_pfx = args.values("name").first().map_or("frame", |s| s.as_str());
    let include_str = args
        .values("include")
        .first()
        .map_or("ALL".to_string(), |s| s.to_uppercase());
    let single = args.has("single");
    let notimeblock = args.has("notimeblock");

    // @pbc <r|t|v> [cog|nog]: with no second value, everything but vacuum gathers
    let pbc = args.values("pbc");
    let pbc_type = pbc.first().map_or("v".to_string(), |s| s.to_lowercase());
    let do_gather = match pbc.get(1) {
        Some(g) => !g.eq_ignore_ascii_case("nog"),
        None => pbc_type != "v",
    };

    // @frames names the frames @spec SPEC (or a plain @spec ALL) writes
    let extra_frames = args
        .values("frames")
        .first()
        .map_or_else(Vec::new, |s| parse_index_list(s));
    let spec = args.values("spec");
    let frame_spec = match spec.first().map(|s| s.to_uppercase()).as_deref() {
        Some("EVERY") => FrameSpec::Every(spec.get(1).and_then(|n| n.parse().ok()).unwrap_or(1)),
        _ if !extra_frames.is_empty() => FrameSpec::Spec(extra_frames),
        _ => FrameSpec::All,
    };

    let tv: Vec<f64> = args
        .values("time")
        .iter()
        .map_while(|v| v.parse().ok())
        .collect();
    let time_range = (
        tv.first().copied().unwrap_or(f64::NEG_INFINITY),
        tv.get(1).copied().unwrap_or(f64::INFINITY),
        tv.get(2).copied().unwrap_or(0.0),
    );

    // Load topology
    let topo = topo_file
        .as_ref()
        .map(|p| {
            read_topology_file(p)
                .map(build_topology)
                .map_err(|e| format!("Error reading topology: {e}"))
        })
        .transpose()?;
    let topo_ref = topo.as_ref();

    // Reference positions for fit
    let reference = ref_file
        .as_ref()
        .map(|p| {
            read_coordinates(p)
                .map(|c| c.positions)
                .map_err(|e| format!("Error reading @ref: {e}"))
        })
        .transpose()?;

    // Fit atom indices
    let fit_indices: Option<Vec<usize>> = match (&topo_ref, &fit_spec, &reference) {
        (Some(t), Some(spec), _) => Some(
            AtomSelection::from_string(spec, t)
                .map_err(|e| format!("Error in @atomsfit '{spec}': {e}"))?
                .indices()
                .to_vec(),
        ),
        (Some(t), None, Some(_)) => Some((0..t.num_solute_atoms()).collect()),
        _ => None,
    };

    // Atom index filter for @include
    let write_indices: Option<Vec<usize>> = topo_ref.map(|t| match include_str.as_str() {
        "SOLUTE" => (0..t.num_solute_atoms()).collect(),
        "SOLVENT" => (t.num_solute_atoms()..t.num_atoms()).collect(),
        _ => (0..t.num_atoms()).collect(),
    });

    // Periodicity builder
    let make_periodicity = |b: Vec3| -> Periodicity {
        match pbc_type.as_str() {
            "r" if b.x > 0.0 => Periodicity::Rectangular(Rectangular::new(b)),
            "t" if b.x > 0.0 => {
                let m = Mat3::from_cols(
                    Vec3::new(b.x, 0., 0.),
                    Vec3::new(0., b.y, 0.),
                    Vec3::new(0., 0., b.z),
                );
                Periodicity::Triclinic(Triclinic::new(m))
            },
            _ => Periodicity::Vacuum(Vacuum),
        }
    };

    // Open input
    let mut reader =
        TrajectoryReader::new(traj_file).map_err(|e| format!("Error opening trajectory: {e}"))?;

    // Single-file .trc writer
    let mut trc_writer: Option<TrajectoryWriter> = if single && outfmt == "trc" {
        let p = format!("{name_pfx}.trc");
        Some(
            TrajectoryWriter::new(&p, "frameout sub-trajectory", false, false)
                .map_err(|e| format!("Error creating {p}: {e}"))?,
        )
    } else {
        None
    };

    let (t_start, t_end, t_dt) = time_range;
    let mut frame_idx = 0usize;
    let mut written = 0usize;

    loop {
        match reader.read_frame() {
            Ok(None) => break,
            Err(e) => {
                eprintln!("Read error: {e}");
                break;
            },
            Ok(Some(frame)) => {
                let t = frame.time;
                let time_ok = t >= t_start - 1e-9
                    && t <= t_end + 1e-9
                    && (t_dt <= 0.0 || ((t - t_start) / t_dt).round() * t_dt + t_start - t < 1e-6);

                if frame_spec.matches(frame_idx) && time_ok {
                    // 1. Gather
                    let mut pos = frame.positions.clone();
                    if do_gather && frame.box_dims.x > 0.0 {
                        if let Some(ref t) = topo {
                            let per = make_periodicity(frame.box_dims);
                            gather_molecules(&mut pos, &t.molecules, &per, Vec3::ZERO);
                        }
                    }

                    // 2. Rotational fit
                    if let (Some(ref refpos), Some(ref fidx)) = (&reference, &fit_indices) {
                        superimpose(&mut pos, refpos, fidx, None);
                    }

                    // 3. Filter atoms
                    let out_pos: Vec<Vec3> = if let Some(ref idx) = write_indices {
                        idx.iter().map(|&i| pos[i]).collect()
                    } else {
                        pos
                    };

                    let box_opt = if frame.box_dims.x > 0.0 {
                        Some(frame.box_dims)
                    } else {
                        None
                    };

                    // Topology for labelling (only when writing all atoms)
                    let topo_for_write = if write_indices
                        .as_ref()
                        .is_some_and(|idx| idx.len() < topo_ref.map_or(0, |t| t.num_atoms()))
                    {
                        None
                    } else {
                        topo_ref
                    };

                    // 4. Write via newen-io
                    let tag = format!("{frame_idx:06}_{:.6}", frame.time).replace('.', "_");

                    if single {
                        if let Some(ref mut tw) = trc_writer {
                            tw.write_trc_frame(frame.step, frame.time, &out_pos, box_opt)
                                .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                        } else {
                            // cnf single: last frame wins
                            let title = if notimeblock {
                                "frameout".to_string()
                            } else {
                                format!("t={:.9} ps step={}", frame.time, frame.step)
                            };
                            write_g96(
                                format!("{name_pfx}.cnf"),
                                &title,
                                &out_pos,
                                None,
                                box_opt,
                                topo_for_write,
                            )
                            .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                        }
                    } else {
                        let ext = match outfmt {
                            "pdb" => "pdb",
                            "trc" => "trc",
                            _ => "cnf",
                        };
                        let path = format!("{name_pfx}_{tag}.{ext}");
                        match outfmt {
                            "pdb" => {
                                let title = format!("t={:.6} ps step={}", frame.time, frame.step);
                                write_pdb_positions(
                                    &path,
                                    &title,
                                    &out_pos,
                                    box_opt,
                                    topo_for_write,
                                )
                                .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                            },
                            "trc" => {
                                let mut tw = TrajectoryWriter::new(&path, "frameout", false, false)
                                    .map_err(|e| format!("Error: {e}"))?;
                                tw.write_trc_frame(frame.step, frame.time, &out_pos, box_opt)
                                    .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                                tw.flush().ok();
                            },
                            _ => {
                                let title = if notimeblock {
                                    "frameout".to_string()
                                } else {
                                    format!("t={:.9} ps step={}", frame.time, frame.step)
                                };
                                write_g96(&path, &title, &out_pos, None, box_opt, topo_for_write)
                                    .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                            },
                        }
                        eprintln!("# frame {frame_idx} t={:.4} ps → {path}", frame.time);
                    }
                    written += 1;
                }
                frame_idx += 1;
            },
        }
    }

    if let Some(ref mut tw) = trc_writer {
        tw.flush().ok();
    }
    eprintln!("# frameout: {written} frames written ({frame_idx} total read)");
    Ok(())
}
