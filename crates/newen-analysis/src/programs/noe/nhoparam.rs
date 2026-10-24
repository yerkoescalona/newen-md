//! nhoparam — N-H bond order parameters for NMR analysis.
//!
//! Calculates the Lipari-Szabo order parameter S² for N-H bonds:
//!
//!   S² = ½ [ 3·Σᵢⱼ ⟨μᵢ(t)μⱼ(t)⟩² − 1 ]
//!
//! where μ is the unit vector along the N-H bond after rotational fitting.
//!
//! Reference: Lipari & Szabo (1982) J. Am. Chem. Soc. 104:4546.
//! Algorithm: faithful port of newen nhoparam (GROMOS analysis tools).
//!
//! Usage:
//!   nhoparam @topo <top> @traj <trc> @atoms <spec>
//!            [@ref <cnf>] [@atomsfit <spec>] [@winframe <n>] [@pbc]

use crate::fit::superimpose;
use crate::programs::io_err;
use newen_core::{
    gather::gather_molecules,
    math::{Periodicity, Rectangular, Vacuum, Vec3},
    selection::AtomSelection,
    stat::Stat,
};
use newen_io::args::{Arguments, Program, Status};
use newen_io::{
    coordinate::read_coordinates,
    topology::{build_topology, read_topology_file},
    trajectory::TrajectoryReader,
};
use std::io::Write;

pub const USAGE: &str = "# nhoparam
\t@topo          <topology>
\t@traj          <trajectory>
\t@atoms         <nitrogen atoms to analyse (e.g. 1:N)>
\t[@ref          <reference coordinates for the fit (default: the first frame)>]
\t[@atomsfit     <atoms used for the rotational fit (default: @atoms)>]
\t[@winframe     <window size for the running average (default: 10)>]
\t[@pbc          <gather the molecules before computing the vectors>]";

pub const KNOWN: &[&str] = &[
    "topo", "traj", "ref", "atoms", "atomsfit", "winframe", "pbc",
];

pub const PROGRAM: Program = Program {
    name: "nhoparam",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

const H_MASS: f64 = 1.10; // upper bound for hydrogen mass (Da)

/// S² from the accumulated tensor sum and number of frames.
fn s2_from_tensor(tensor: &[[f64; 3]; 3], n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for a in 0..3 {
        for b in 0..3 {
            let ave = tensor[a][b] / n as f64;
            sum += ave * ave;
        }
    }
    0.5 * (3.0 * sum - 1.0)
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    if argv.len() < 2 {
        return Err(USAGE.to_string());
    }
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let topo_file = args.value("topo")?;
    let traj_file = args.value("traj")?;
    let atom_spec = args.value("atoms")?;
    let ref_file = args.values("ref").first().cloned();
    let winframe: usize = args.get("winframe", 10)?;
    let do_pbc = args.has("pbc");

    // Load topology
    let topo_data =
        read_topology_file(topo_file).map_err(|e| format!("Error reading topology: {e}"))?;
    let topo = build_topology(topo_data);

    // Resolve N atoms selection
    let n_sel = AtomSelection::from_string(atom_spec, &topo)
        .map_err(|e| format!("Error in @atoms '{atom_spec}': {e}"))?;
    let n_indices: Vec<usize> = n_sel.indices().to_vec();
    eprintln!("# N atoms selected: {}", n_indices.len());

    // Verify all selected atoms are nitrogen (mass > H_MASS)
    for &ni in &n_indices {
        let mass = topo.mass.get(ni).copied().unwrap_or(0.0);
        if mass < 12.0 {
            eprintln!("Warning: atom {ni} has mass {mass:.3} — may not be nitrogen");
        }
    }

    // For each N, find bonded H atoms (mass < H_MASS)
    let mut nh_pairs: Vec<(usize, Vec<usize>)> = Vec::new(); // (N_idx, [H_idx, ...])
    for &ni in &n_indices {
        let mut h_atoms: Vec<usize> = Vec::new();
        for bond in &topo.moltypes[0].bonds {
            let (a, b) = (bond.i, bond.j);
            let partner = if a == ni {
                b
            } else if b == ni {
                a
            } else {
                continue;
            };
            let h_mass = topo.mass.get(partner).copied().unwrap_or(99.0);
            if h_mass < H_MASS {
                h_atoms.push(partner);
            }
        }
        if h_atoms.is_empty() {
            eprintln!("Warning: N atom {ni} has no bonded H — skipping");
        } else {
            nh_pairs.push((ni, h_atoms));
        }
    }
    eprintln!("# N-H pairs found: {}", nh_pairs.len());

    // Resolve fit atoms
    let fit_spec_str = args
        .values("atomsfit")
        .first()
        .map_or(atom_spec, |s| s.as_str());
    let fit_sel = AtomSelection::from_string(fit_spec_str, &topo)
        .map_err(|e| format!("Error in @atomsfit '{fit_spec_str}': {e}"))?;
    let fit_indices: Vec<usize> = fit_sel.indices().to_vec();

    // Read all trajectory frames
    let mut reader =
        TrajectoryReader::new(traj_file).map_err(|e| format!("Error opening trajectory: {e}"))?;
    let frames = reader
        .read_all_frames()
        .map_err(|e| format!("Error reading trajectory: {e}"))?;
    if frames.is_empty() {
        return Err("no frames in trajectory".into());
    }
    eprintln!("# Frames: {}", frames.len());

    // Reference positions: from @ref file or first frame
    let reference: Vec<Vec3> = if let Some(ref path) = ref_file {
        read_coordinates(path)
            .map_err(|e| format!("Error reading @ref: {e}"))?
            .positions
    } else {
        frames[0].positions.clone()
    };

    // Accumulate order parameter tensors (one 3×3 per N-H pair)
    let np = nh_pairs.len();
    let mut tensors: Vec<[[f64; 3]; 3]> = vec![[[0.0; 3]; 3]; np];
    let mut win_tensors: Vec<[[f64; 3]; 3]> = vec![[[0.0; 3]; 3]; np];
    let mut win_stats: Vec<Stat> = (0..np).map(|_| Stat::new()).collect();
    let mut n_frames = 0usize;
    let mut win_count = 0usize;

    for frame in &frames {
        let mut pos = frame.positions.clone();

        // PBC gathering
        if do_pbc {
            let periodicity = if frame.box_dims.x > 0.0 {
                Periodicity::Rectangular(Rectangular::new(frame.box_dims))
            } else {
                Periodicity::Vacuum(Vacuum)
            };
            gather_molecules(&mut pos, &topo.molecules, &periodicity, Vec3::ZERO);
        }

        // Rotational fit
        superimpose(&mut pos, &reference, &fit_indices, None);

        // Accumulate N-H tensors
        for (k, (ni, h_list)) in nh_pairs.iter().enumerate() {
            // Average unit vector over bonded H atoms
            let mut mu = Vec3::ZERO;
            for &hi in h_list {
                let v = pos[hi] - pos[*ni];
                let len = v.length();
                if len > 1e-10 {
                    mu += v / len;
                }
            }
            // Normalize average (handles both single and multiple H)
            let mu_len = mu.length();
            if mu_len < 1e-10 {
                continue;
            }
            let mu = mu / mu_len;

            let mv = [mu.x, mu.y, mu.z];
            for a in 0..3 {
                for b in 0..3 {
                    tensors[k][a][b] += mv[a] * mv[b];
                    win_tensors[k][a][b] += mv[a] * mv[b];
                }
            }
        }

        n_frames += 1;
        win_count += 1;

        // Window averaging
        if win_count >= winframe {
            for (k, ws) in win_stats.iter_mut().enumerate() {
                ws.add(s2_from_tensor(&win_tensors[k], win_count));
                win_tensors[k] = [[0.0; 3]; 3];
            }
            win_count = 0;
        }
    }

    // Flush incomplete window
    if win_count > 0 {
        for (k, ws) in win_stats.iter_mut().enumerate() {
            ws.add(s2_from_tensor(&win_tensors[k], win_count));
        }
    }

    // Output
    writeln!(out, "# nhoparam — N-H order parameters (S²)").map_err(io_err)?;
    writeln!(
        out,
        "# topo={topo_file}  traj={traj_file}  frames={n_frames}  winframe={winframe}"
    )
    .map_err(io_err)?;
    writeln!(
        out,
        "# {:>10} {:>6} {:>8} {:>10} {:>10} {:>10} {:>10}",
        "ATOM", "RESNUM", "RESNAME", "S2", "WINAV", "STDDEV", "EE"
    )
    .map_err(io_err)?;

    for (k, (ni, _h_list)) in nh_pairs.iter().enumerate() {
        let s2 = s2_from_tensor(&tensors[k], n_frames);
        let res_nr = topo.residue_nr(*ni).unwrap_or(0);
        let res_name = topo.residue_name(*ni).unwrap_or("UNK");
        let atom_name = topo.atom_name(*ni).unwrap_or("N");
        let label = format!("{}:{}", 1, ni + 1); // mol:atom (1-based)

        let (winav, stddev, ee) = if win_stats[k].n() >= 2 {
            (win_stats[k].ave(), win_stats[k].rmsd(), win_stats[k].ee())
        } else if win_stats[k].n() == 1 {
            (win_stats[k].ave(), 0.0, 0.0)
        } else {
            (s2, 0.0, 0.0)
        };

        writeln!(
            out,
            "  {:>10} {:>6} {:>8} {:>10.6} {:>10.6} {:>10.6} {:>10.6}",
            label,
            res_nr,
            format!("{res_name}/{atom_name}"),
            s2,
            winav,
            stddev,
            ee
        )
        .map_err(io_err)?;
    }
    Ok(())
}
