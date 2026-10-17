//! dfgrid — the distance-field grid around a target atom (gromos++ `dfgrid`).
//!
//! Usage: dfgrid @topo <top> @pbc <r> [cog] @atom <atom> [@distatoms <atoms>] @gridspacing <nm>
//!        @proteinoffset <penalty> @proteincutoff <nm> @proteinatoms <atoms> [@max <n>]
//!        [@nogrid] [@smooth <rounds>] [@protect <nm>] [@notimeblock] [@time <t0 dt>]
//!        [@stride <n> | @frames <n…>] @traj <trc…>
//!
//! The distance field is the length of the shortest grid path from the target atom, where grid
//! points within `@proteincutoff` of a protein atom (farther than `@protect` from the target)
//! cost `@proteinoffset` extra to enter — the distance-field restraint of gromosXX. gromos++'s
//! scan-the-open-set Dijkstra is kept so that equal-length paths resolve identically. Per
//! written frame `grid%05d.cnf` holds the solute, the target (VA), the `@distatoms` (DA), their
//! minimal paths (PA) and every grid point closer than `@max` (AR); with `@distatoms` the field
//! distances go to stdout, flagging atoms inside the protein region. Only the `cnf` output format
//! and real (non-virtual) atoms are supported; the box must be rectangular.

use crate::args::Arguments;
use crate::distribution::cpp_g;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::{atoms_to_string, ordered_atoms};
use crate::time::Time;
use newen_core::math::Vec3;
use newen_core::Topology;
use newen_io::args::{Program, Status};
use newen_io::coordinate::{format_g96, G96Atom};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::collections::BTreeSet;
use std::io::Write;

pub const USAGE: &str = "# dfgrid
\t@topo     <molecular topology file>
\t@pbc      <boundary type> [<gather method>]
\t@atom     <atom specifier for the target (zero-distance) point>
\t[@distatoms     <atom specifier for atoms for which to output the df distance>]
\t@gridspacing   <grid spacing>
\t@proteinoffset <penalty for being in the protein>
\t@proteincutoff <cutoff to determine gridpoints within the protein>
\t@proteinatoms  <atomspecifier for everything to be penalized as protein>
\t[@max     <maximum distance: do not write out grid points with higher distances> (default: 1)]
\t[@nogrid  <do not write out grid coordinate file>]
\t[@smooth  <number of rounds to smoothen the forces at the edge of the protein>]
\t[@protect <radius around the target atom that will not be flagged as protein>]
\t[@outformat   <output coordinates format (cnf)>]
\t[@notimeblock <do not write timestep block>]
\t[@time    <time and dt>]
\t[@stride  <write every nth frame> (default: 1)]
\t[@frames  <select frames to write out, starts at 0> (default: 0)]
\t@traj     <input trajectory files>";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "traj",
    "atom",
    "proteincutoff",
    "proteinoffset",
    "proteinatoms",
    "gridspacing",
    "smooth",
    "outformat",
    "stride",
    "notimeblock",
    "time",
    "max",
    "protect",
    "frames",
    "distatoms",
    "nogrid",
];

pub const PROGRAM: Program = Program {
    name: "dfgrid",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

struct Grid {
    n: [i64; 3],
    spacing: f64,
    half: Vec3,
}

impl Grid {
    fn ndim(&self) -> i64 {
        self.n[0] * self.n[1] * self.n[2]
    }

    fn coords_of(&self, ix: i64, iy: i64, iz: i64) -> Vec3 {
        Vec3::new(
            ix as f64 * self.spacing - self.half.x,
            iy as f64 * self.spacing - self.half.y,
            iz as f64 * self.spacing - self.half.z,
        )
    }

    fn coords(&self, idx: i64) -> Vec3 {
        let nxy = self.n[0] * self.n[1];
        let nz = idx / nxy;
        let ny = (idx % nxy) / self.n[0];
        let nx = (idx % nxy) % self.n[0];
        self.coords_of(nx, ny, nz)
    }

    fn index(&self, ix: i64, iy: i64, iz: i64) -> i64 {
        ix + self.n[0] * iy + self.n[0] * self.n[1] * iz
    }

    /// the grid point nearest to `p`
    fn nearest(&self, p: Vec3) -> (i64, i64, i64) {
        (
            (-(-p.x - self.half.x) / self.spacing + 0.5) as i64,
            (-(-p.y - self.half.y) / self.spacing + 0.5) as i64,
            (-(-p.z - self.half.z) / self.spacing + 0.5) as i64,
        )
    }

    /// the periodic neighbour `j` (−x, +x, −y, +y, −z, +z) of point `i`
    fn neighbour(&self, i: i64, j: usize) -> i64 {
        let (n0, nxy, n) = (self.n[0], self.n[0] * self.n[1], self.ndim());
        match j {
            0 => {
                if i % n0 == 0 {
                    i - 1 + n0
                } else {
                    i - 1
                }
            },
            1 => {
                if (i + 1) % n0 == 0 {
                    i + 1 - n0
                } else {
                    i + 1
                }
            },
            2 => {
                if (i - (i % n0)) % nxy == 0 {
                    i - n0 + nxy
                } else {
                    i - n0
                }
            },
            3 => {
                if (i + n0 - (i % n0)) % nxy == 0 {
                    i + n0 - nxy
                } else {
                    i + n0
                }
            },
            4 => {
                if i - (i % nxy) == 0 {
                    i - nxy + n
                } else {
                    i - nxy
                }
            },
            _ => {
                if i + nxy - (i % nxy) == n {
                    i + nxy - n
                } else {
                    i + nxy
                }
            },
        }
    }
}

/// gromos++ `AtomSpecifier`: the union of the specifiers in the order given, without repeats
fn selection(specs: &[String], topo: &Topology) -> Result<Vec<usize>, String> {
    let mut v: Vec<usize> = Vec::new();
    for s in specs {
        for a in ordered_atoms(s, topo)? {
            if !v.contains(&a) {
                v.push(a);
            }
        }
    }
    Ok(v)
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let parsed = read_topology_file(args.value("topo")?).map_err(|e| e.to_string())?;
    let mut topo = build_topology(parsed);
    let trajs = args.values("traj");
    let first_traj = trajs.first().ok_or("no @traj given")?;
    let first_frame = TrajectoryReader::new(first_traj)
        .map_err(|e| format!("@traj {first_traj}: {e}"))?
        .read_frame()
        .map_err(|e| format!("{first_traj}: {e}"))?
        .ok_or(format!("@traj {first_traj}: no frame"))?;
    solvate_to_atoms(&mut topo, first_frame.positions.len()).map_err(|e| e.to_string())?;
    let mut time = Time::from_args(&args)?;
    let notimeblock = args.count("notimeblock") >= 0;
    if args.count("stride") >= 0 && args.count("frames") >= 0 {
        return Err("choose either @stride or @frames".into());
    }
    let frames: Vec<usize> = args
        .values("frames")
        .iter()
        .map(|s| {
            s.parse::<usize>()
                .map_err(|_| format!("@frames: '{s}' is not a number"))
        })
        .collect::<Result<_, _>>()?;
    let pbc = Pbc::from_args(&args)?;
    let stride: usize = args.get("stride", 1)?;
    if stride == 0 {
        return Err("@stride must be positive".into());
    }
    let smooth: usize = args.get("smooth", 1)?;
    let protect: f64 = args.get("protect", -1.0)?;
    let max: i64 = args.get("max", 1)?;
    let spacing: f64 = args.require("gridspacing")?;
    let offset: f64 = args.require("proteinoffset")?;
    let cutoff: f64 = args.require("proteincutoff")?;
    let write_grid = args.count("nogrid") < 0;
    if args.has("outformat") && args.value("outformat")? != "cnf" {
        return Err("only the cnf output format is supported".into());
    }
    let atom_spec = args
        .values("atom")
        .first()
        .ok_or("@atom: no atom selected")?;
    if atom_spec.trim_start().starts_with("va(") {
        return Err("@atom: virtual atoms are not supported".into());
    }
    let target = ordered_atoms(atom_spec, &topo)?;
    if target.is_empty() {
        return Err("@atom: no atom selected".into());
    }
    if target.len() > 1 {
        return Err("more than 1 atom in selection, need exactly one".into());
    }
    let target = target[0];
    let distatoms = selection(args.values("distatoms"), &topo)?;
    let proteinatoms = selection(args.values("proteinatoms"), &topo)?;
    let b0 = first_frame.box_dims;
    if b0.x <= 0.0 || b0.y <= 0.0 || b0.z <= 0.0 {
        return Err("only rectangular boxshape supported".into());
    }
    let n = [
        (b0.x / spacing) as i64 + 1,
        (b0.y / spacing) as i64 + 1,
        (b0.z / spacing) as i64 + 1,
    ];
    let n_solute = topo.num_solute_atoms();
    let res_off = 1
        + (0..n_solute)
            .filter_map(|i| topo.residue_nr(i))
            .max()
            .unwrap_or(0);
    let mut text = String::new();
    if !distatoms.is_empty() {
        text.push_str("# ");
        for &a in &distatoms {
            text.push_str(&format!("{} ", atoms_to_string(&[a], &topo)));
        }
        text.push_str("# in protein\n");
    }
    let mut num_frames = 0usize;
    let mut skip = 0usize;
    for traj in trajs {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let t = if notimeblock {
                frame.time
            } else {
                time.next(frame.time)
            };
            let write_frame = frames.contains(&num_frames);
            if (skip == 0 && frames.is_empty()) || write_frame {
                let b = frame.box_dims;
                if b.x <= 0.0 || b.y <= 0.0 || b.z <= 0.0 {
                    return Err("only rectangular boxshape supported".into());
                }
                let grid = Grid {
                    n,
                    spacing,
                    half: Vec3::new(b.x / 2.0, b.y / 2.0, b.z / 2.0),
                };
                let ndim = grid.ndim();
                let periodicity = pbc.periodicity(b);
                let mut pos = frame.positions.clone();
                pbc.gather(&topo, &mut pos, &periodicity);
                // the image nearest to the origin
                let image = |p: Vec3| -periodicity.nearest_image(Vec3::ZERO, p);
                let ppos = image(pos[target]);
                let (nx, ny, nz) = grid.nearest(ppos);
                let startpos = grid.coords_of(nx, ny, nz);
                let start = grid.index(nx, ny, nz);
                if start < 0 || start >= ndim {
                    return Err("the target atom lies outside the grid".into());
                }
                let dist_to_i = (startpos - ppos).length();
                let far = 4.0 * (b.x + b.y + b.z);
                let mut distance = vec![far; ndim as usize];
                let mut from = vec![0i64; ndim as usize];
                let mut protein: BTreeSet<i64> = BTreeSet::new();
                for &a in &proteinatoms {
                    let pp = image(pos[a]);
                    let lo = |c: f64, h: f64| (-(cutoff - c - h) / spacing) as i64;
                    let hi = |c: f64, h: f64| (-(-cutoff - c - h) / spacing) as i64 + 1;
                    for ix in lo(pp.x, grid.half.x)..hi(pp.x, grid.half.x) {
                        for iy in lo(pp.y, grid.half.y)..hi(pp.y, grid.half.y) {
                            for iz in lo(pp.z, grid.half.z)..hi(pp.z, grid.half.z) {
                                let gpos = grid.coords_of(ix, iy, iz);
                                let d = gpos - pp;
                                let e = startpos - gpos;
                                if d.length() < cutoff && e.length() > protect {
                                    let j = grid.index(ix, iy, iz);
                                    if (0..ndim).contains(&j) {
                                        protein.insert(j);
                                    }
                                }
                            }
                        }
                    }
                }
                // gromos++'s Dijkstra: the open set is scanned for its minimum each step
                let mut current = start;
                distance[start as usize] = dist_to_i;
                if protein.contains(&start) {
                    distance[start as usize] += offset;
                }
                let mut visited = vec![false; ndim as usize];
                let mut n_visited = 0i64;
                let mut todo: BTreeSet<i64> = BTreeSet::new();
                while n_visited != ndim {
                    if !visited[current as usize] {
                        visited[current as usize] = true;
                        n_visited += 1;
                    }
                    todo.remove(&current);
                    let currdist = distance[current as usize];
                    for k in 0..6 {
                        let mut newdist = currdist + spacing;
                        let j = grid.neighbour(current, k);
                        if !visited[j as usize] {
                            if protein.contains(&j) && !protein.contains(&current) {
                                newdist += offset;
                            }
                            if newdist < distance[j as usize] {
                                distance[j as usize] = newdist;
                                from[j as usize] = current;
                            }
                            todo.insert(j);
                        }
                    }
                    let mut next = current;
                    let mut min = far + 1.0;
                    for &c in &todo {
                        if distance[c as usize] < min {
                            next = c;
                            min = distance[c as usize];
                        }
                    }
                    if next == current {
                        break;
                    }
                    current = next;
                }
                for _ in 0..smooth {
                    let mut remove = BTreeSet::new();
                    for &p in &protein {
                        for k in 0..6 {
                            let j = grid.neighbour(p, k);
                            if !protein.contains(&j)
                                && distance[p as usize] > distance[j as usize] + spacing
                            {
                                distance[p as usize] = distance[j as usize] + spacing;
                                remove.insert(p);
                            }
                        }
                    }
                    for p in remove {
                        protein.remove(&p);
                    }
                }
                let mut s_points = Vec::new();
                let mut minpaths: Vec<Vec<i64>> = Vec::new();
                for j in 0..ndim {
                    if distance[j as usize] < max as f64 {
                        s_points.push(grid.coords(j));
                    }
                }
                if !distatoms.is_empty() {
                    text.push_str(&format!("{num_frames} "));
                    let mut in_protein = String::new();
                    for &a in &distatoms {
                        let p = image(pos[a]);
                        let (x, y, z) = grid.nearest(p);
                        let gridpos = grid.coords_of(x, y, z);
                        let idx = grid.index(x, y, z);
                        if idx < 0 || idx >= ndim {
                            return Err("a @distatoms atom lies outside the grid".into());
                        }
                        let dist = distance[idx as usize] + (gridpos - p).length();
                        text.push_str(&format!("{} ", cpp_g(dist, 6)));
                        if protein.contains(&idx) {
                            in_protein.push_str(&format!("{} ", atoms_to_string(&[a], &topo)));
                        }
                        let mut path = vec![idx];
                        let mut prev = from[idx as usize];
                        while prev != start && (path.len() as i64) < ndim {
                            path.push(prev);
                            prev = from[prev as usize];
                        }
                        path.push(prev);
                        minpaths.push(path);
                    }
                    text.push_str(&format!("# {in_protein}\n"));
                }
                if write_grid {
                    let mut atoms: Vec<G96Atom> = (0..n_solute)
                        .map(|i| G96Atom {
                            res_num: topo.residue_nr(i).unwrap_or(1),
                            res_name: topo.residue_name(i).unwrap_or("").to_string(),
                            atom_name: topo.atom_name(i).unwrap_or("").to_string(),
                            atom_num: i + 1,
                            pos: pos[i],
                        })
                        .collect();
                    let mut res = 1;
                    let push = |name: &str, res: usize, p: Vec3, atoms: &mut Vec<G96Atom>| {
                        atoms.push(G96Atom {
                            res_num: res + res_off,
                            res_name: name.to_string(),
                            atom_name: name.to_string(),
                            atom_num: atoms.len() + 1,
                            pos: p,
                        });
                    };
                    push("VA", res, pos[target], &mut atoms);
                    res += 1;
                    for &a in &distatoms {
                        push("DA", res, pos[a], &mut atoms);
                        res += 1;
                    }
                    for path in &minpaths {
                        for &idx in path {
                            push("PA", res, grid.coords(idx), &mut atoms);
                        }
                        res += 1;
                    }
                    for &p in &s_points {
                        push("AR", res, p, &mut atoms);
                    }
                    let title = format!(
                        "{traj}, Frame{num_frames}\nAR: gridpoints, VA: center atom, DA: distatoms\n\n"
                    );
                    let mut text = format_g96(&title, &atoms, Some(b));
                    if !notimeblock {
                        let ts = format!(
                            "TIMESTEP\n{:>18}{:>20.9}\n#if @time flag is used the value for step refers to the\n#step-th configuration in the original trajectory file\nEND\n",
                            frame.step, t
                        );
                        let at = text.find("END\n").map(|i| i + 4).unwrap_or(0);
                        text.insert_str(at, &ts);
                    }
                    let name = format!("grid{num_frames:05}.cnf");
                    std::fs::write(&name, text).map_err(|e| format!("{name}: {e}"))?;
                }
            }
            num_frames += 1;
            skip += 1;
            skip %= stride;
        }
    }
    write!(out, "{text}").map_err(io_err)?;
    Ok(())
}
