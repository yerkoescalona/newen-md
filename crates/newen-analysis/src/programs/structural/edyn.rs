//! edyn — essential dynamics: covariance matrix of the fitted atom positions, its eigenvalues
//! and eigenvectors, and the projections of the trajectory on selected eigenvectors
//! (gromos++ `edyn`). Writes gromos++'s files into the working directory: `AVE.pdb`,
//! `COVAR.out`, `COVATOM.out`, `EIVAL.out`, `EIFLUC.out`, `EIVEC.out`, and per selected
//! eigenvector `EVCOMP_n.out`, `EVPRJ_n.out`, `DXPRJ_n`, `PRJMAX_n.pdb`, `PRJMIN_n.pdb`, plus
//! `ESSDYN.out`.
//!
//! Usage: edyn @topo <top> @pbc <r|v> @atoms <spec> [@ref <cnf>] [@eigenvalues <n…>] [@skip] @traj <trc…>

use crate::args::Arguments;
use crate::distribution::trim_float;
use crate::fit::superimpose;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::ordered_atoms;
use nalgebra::DMatrix;
use newen_core::Vec3;
use newen_io::args::{Program, Status};
use newen_io::coordinate::read_coordinates;
use newen_io::pdb::write_pdb_positions;
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::io::Write;

pub const USAGE: &str = "# edyn
\t@topo         <molecular topology file>
\t@pbc          <boundary type>
\t@atoms        <atoms to be considered>
\t@ref          <reference coordinates>
\t[@eigenvalues <list of eigenvalues for which data is written>]
\t@traj         <trajectory files>
\t[@skip         (skip the (time-consuming) projections)]";

pub const KNOWN: &[&str] = &["topo", "traj", "atoms", "pbc", "ref", "eigenvalues", "skip"];

pub const PROGRAM: Program = Program {
    name: "edyn",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

fn write_file(name: &str, text: &str) -> Result<(), String> {
    std::fs::write(name, text).map_err(|e| format!("{name}: {e}"))
}

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
    let mut atoms = ordered_atoms(&args.values("atoms").join(";"), &topo)?;
    if atoms.is_empty() {
        return Err("No atoms specified!".into());
    }
    atoms.sort();
    atoms.dedup();
    let skip = args.count("skip") >= 0;
    let size = atoms.len();
    let ndim = 3 * size;
    let mut sel: Vec<usize> = Vec::new();
    for v in args.values("eigenvalues") {
        let t: usize = v
            .parse::<usize>()
            .map_err(|_| format!("@eigenvalues: bad number {v}"))?;
        if t >= 1 && t - 1 < ndim {
            sel.push(t - 1);
        } else {
            writeln!(out, "You specified {} atoms.\nThis leads to a covariance matrix with dimensions {ndim}x{ndim}\nSo we will have at most {ndim} eigenvalues.\nIgnoring request for eigenvalue {t}\n", size).map_err(io_err)?;
        }
    }
    writeln!(out, "Selected {} eigenvalues", sel.len()).map_err(io_err)?;
    let pbc = Pbc::from_args(&args)?;
    // reference: @ref or the first frame of the trajectory, gathered and shifted to the COG of
    // the selected atoms (gromos++ Reference + RotationalFit put both there)
    let (mut reference, ref_box) = match args.values("ref").first() {
        Some(r) => {
            let c = read_coordinates(r).map_err(|e| format!("@ref: {e}"))?;
            (c.positions, c.box_dims)
        },
        None => {
            let f = TrajectoryReader::new(args.value("traj")?)
                .map_err(|e| e.to_string())?
                .read_frame()
                .map_err(|e| e.to_string())?
                .ok_or("empty trajectory")?;
            (f.positions, f.box_dims)
        },
    };
    {
        let periodicity = pbc.periodicity(ref_box);
        pbc.gather(&topo, &mut reference, &periodicity);
        let cog = atoms
            .iter()
            .map(|&a| reference[a])
            .fold(Vec3::ZERO, |s, p| s + p)
            / size as f64;
        for p in reference.iter_mut() {
            *p -= cog;
        }
    }
    let mut num_frames = 0usize;
    let mut avpos = vec![Vec3::ZERO; size];
    let mut cov = DMatrix::<f64>::zeros(ndim, ndim);
    writeln!(out, "reading trajectory...").map_err(io_err)?;
    let mut fitted_frames: Vec<Vec<Vec3>> = Vec::new();
    for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            num_frames += 1;
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            superimpose(&mut pos, &reference, &atoms, None);
            let pvec: Vec<Vec3> = atoms.iter().map(|&a| pos[a]).collect();
            for (i, p) in pvec.iter().enumerate() {
                avpos[i] += *p;
            }
            for ii in 0..size {
                for jj in 0..=ii {
                    for iii in 0..3 {
                        for jjj in 0..3 {
                            cov[(3 * ii + iii, 3 * jj + jjj)] += pvec[ii][iii] * pvec[jj][jjj];
                        }
                    }
                }
            }
            if !skip {
                fitted_frames.push(pvec);
            }
        }
    }
    if num_frames == 0 {
        return Err("no frames".into());
    }
    for p in avpos.iter_mut() {
        *p /= num_frames as f64;
    }
    write_pdb_positions("AVE.pdb", "AVE", &avpos, None, None).map_err(|e| e.to_string())?;
    writeln!(out, "building matrix...").map_err(io_err)?;
    if num_frames < ndim {
        writeln!(out, "The number of dimensions ({ndim}) is larger than the total number of frames ({num_frames})!\nThis may lead to poor results.\n").map_err(io_err)?;
    }
    for ii in 0..size {
        for jj in 0..=ii {
            for iii in 0..3 {
                for jjj in 0..3 {
                    let v = cov[(3 * ii + iii, 3 * jj + jjj)] / num_frames as f64
                        - avpos[ii][iii] * avpos[jj][jjj];
                    cov[(3 * ii + iii, 3 * jj + jjj)] = v;
                }
            }
        }
    }
    let tcov: f64 = (0..ndim).map(|z| cov[(z, z)]).sum();
    if tcov < 0.0 {
        return Err(" trace of the covariance matrix is negative. Cannot go on. Something might be wrong with your trajectory.".into());
    }
    for ii in 0..ndim {
        for jj in 0..=ii {
            let v = cov[(ii, jj)];
            cov[(jj, ii)] = v;
        }
    }
    let mut ocov = format!("# Covariance matrix of {ndim} x {ndim}\n");
    let mut ocov2 = String::from("# Covariance matrix reduced to atomic correlations\n");
    for i in 0..ndim {
        for j in 0..ndim {
            ocov.push_str(&format!("{:>15}\n", trim_float(cov[(i, j)])));
        }
    }
    for i in (0..ndim).step_by(3) {
        for j in (0..ndim).step_by(3) {
            let fac: f64 = (0..3).map(|k| cov[(i + k, j + k)]).sum();
            ocov2.push_str(&format!("{:>15}\n", trim_float(fac)));
        }
    }
    write_file("COVAR.out", &ocov)?;
    write_file("COVATOM.out", &ocov2)?;
    writeln!(out, "diagonalizing matrix...").map_err(io_err)?;
    let eig = cov.clone().symmetric_eigen();
    // gromos++: eigenvalues sorted descending, eigenvectors as columns
    let mut order: Vec<usize> = (0..ndim).collect();
    order.sort_by(|&a, &b| eig.eigenvalues[b].partial_cmp(&eig.eigenvalues[a]).unwrap());
    let eigen: Vec<f64> = order.iter().map(|&i| eig.eigenvalues[i]).collect();
    let evec = |j: usize, i: usize| -> f64 { eig.eigenvectors[(j, order[i])] };
    let tdcov: f64 = eigen.iter().sum();
    if tdcov < 0.0 {
        return Err(" trace of the diagonalized matrix is negative. Cannot go on. Something might be wrong with your trajectory.".into());
    }
    if (tdcov - tcov).abs() > 0.01 * (tdcov + tcov) {
        return Err(" trace of the covariance and the diagonalized matrix deviate too much. Cannot go on. Something went wrong during the diagonalization. Check your trajectory.".into());
    }
    let mut oeig = String::from("Eigenvalues\n");
    for (i, e) in eigen.iter().enumerate() {
        oeig.push_str(&format!("{} {}\n", i + 1, trim_float(*e)));
    }
    write_file("EIVAL.out", &oeig)?;
    let mut orel = String::from("Relative Fluctuations of the Eigenvalues\n");
    let mut refl = 0.0;
    for (i, e) in eigen.iter().enumerate() {
        refl += e / tdcov;
        orel.push_str(&format!("{} {}\n", i + 1, trim_float(refl)));
    }
    write_file("EIFLUC.out", &orel)?;
    let mut oeiv = String::from("Eigenvectors\n");
    let mut x = 0;
    // gromos++ loops ii ≤ NDIM — one column past the last eigenvector, which reads as zeros here
    for ii in 0..=ndim {
        for jj in 0..ndim {
            let v = if ii < ndim { evec(jj, ii) } else { 0.0 };
            oeiv.push_str(&format!("{:>13} ", trim_float(v)));
            x += 1;
            if x == 6 {
                oeiv.push('\n');
                x = 0;
            }
        }
    }
    write_file("EIVEC.out", &oeiv)?;
    if skip {
        return Ok(());
    }
    for &s in &sel {
        let mut outf = String::new();
        for (j, &atom) in atoms.iter().enumerate() {
            let tmp: f64 = (0..3)
                .map(|k| evec(3 * j + k, s).powi(2))
                .sum::<f64>()
                .sqrt();
            outf.push_str(&format!("{}  {}\n", atom + 1, trim_float(tmp)));
        }
        write_file(&format!("EVCOMP_{}.out", s + 1), &outf)?;
    }
    let mut minprj = vec![100000.0f64; sel.len()];
    let mut maxprj = vec![-100000.0f64; sel.len()];
    let mut avs = vec![0.0f64; sel.len()];
    let mut avsq = vec![0.0f64; sel.len()];
    writeln!(out, "putting EV's in EIG").map_err(io_err)?;
    let mut prj_files: Vec<String> = vec![String::new(); sel.len()];
    for (f, pvec) in fitted_frames.iter().enumerate() {
        for (i, &s) in sel.iter().enumerate() {
            let mut proj = 0.0;
            for j in 0..size {
                for k in 0..3 {
                    proj += (pvec[j][k] - avpos[j][k]) * evec(3 * j + k, s);
                }
            }
            let proj = -proj;
            prj_files[i].push_str(&format!("{} {}\n", f + 1, trim_float(proj)));
            maxprj[i] = maxprj[i].max(proj);
            minprj[i] = minprj[i].min(proj);
            avs[i] += proj;
            avsq[i] += proj * proj;
        }
    }
    for (i, &s) in sel.iter().enumerate() {
        write_file(&format!("EVPRJ_{}.out", s + 1), &prj_files[i])?;
        let mut dis = String::new();
        for z in 0..60 {
            let dx = (z as f64 - 4.0) * ((maxprj[i] - minprj[i]) / 50.0) + minprj[i];
            dis.push_str(&format!("{}\n", trim_float(dx)));
        }
        write_file(&format!("DXPRJ_{}", s + 1), &dis)?;
    }
    let mut output = String::from("# Projection of trajectory along the eigenvectors\n#\n");
    output.push_str(&format!(
        "#  EV{:>15}{:>15}{:>15}{:>15}\n",
        "average", "fluctuation", "min. proj.", "max. proj."
    ));
    let nf = num_frames as f64;
    for (i, &s) in sel.iter().enumerate() {
        let av = avs[i] / nf;
        let sig = (avsq[i] / nf - av * av).sqrt();
        output.push_str(&format!(
            "{:>4}{:>15}{:>15}{:>15}{:>15}\n",
            s + 1,
            trim_float(av),
            trim_float(sig),
            trim_float(minprj[i]),
            trim_float(maxprj[i])
        ));
        for (which, ext) in [(maxprj[i], "PRJMAX"), (minprj[i], "PRJMIN")] {
            let pos: Vec<Vec3> = (0..size)
                .map(|j| {
                    Vec3::new(
                        avpos[j].x + evec(3 * j, s) * which,
                        avpos[j].y + evec(3 * j + 1, s) * which,
                        avpos[j].z + evec(3 * j + 2, s) * which,
                    )
                })
                .collect();
            write_pdb_positions(format!("{ext}_{}.pdb", s + 1), ext, &pos, None, None)
                .map_err(|e| e.to_string())?;
        }
    }
    write_file("ESSDYN.out", &output)?;
    Ok(())
}
