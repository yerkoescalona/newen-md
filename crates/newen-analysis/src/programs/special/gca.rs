//! gca — generate conformations: set a bond length, bond angle or dihedral of a molecule to
//! given values by moving the atoms on one side (gromos++ `gca`).
//!
//! Usage: gca @topo <top> @pbc <r|v> [cog|nog] @prop <properties with values> [@mobile first|last] @traj <cnf>
//!
//! A property is `d%…%value`, `a%…%value`, `t%…%value` or `type%atoms%step%min%max`; every
//! combination of values is written as one frame (cnf; other `@outformat`s are not
//! supported). The property must be a bond, angle or (improper) dihedral of the topology.

use crate::args::Arguments;
use crate::distribution::trim_float;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::{Kind, Property, PropertyContainer};
use newen_core::{Mat3, Vec3};
use newen_io::args::{Program, Status};
use newen_io::coordinate::{format_g96, read_g96_labeled, G96Atom};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use std::collections::BTreeSet;
use std::io::Write;

pub const USAGE: &str = "# gca
\t@topo       <molecular topology file>
\t@pbc        <boundary type> [<gather method>]
\t@prop       <properties to change>
\t[@mobile       <which part of the molecule should be mobile: first or last (default)>]
\t[@outformat <output coordinates format>]
\t@traj       <input coordinate trajectory>";

pub const KNOWN: &[&str] = &["topo", "pbc", "prop", "outformat", "traj", "mobile"];

pub const PROGRAM: Program = Program {
    name: "gca",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// gromos++ `PositionUtils::rotateAround(v, a)`: rot2 · m3 · rot1, matrices built column-wise.
fn rotate_around(v: Vec3, a_deg: f64) -> Mat3 {
    let a = a_deg.to_radians();
    let r = v.length();
    let r_yz = (v.y * v.y + v.z * v.z).sqrt();
    // gromos++ `Matrix(u, v, w)` takes its three vectors as *columns*
    let rot1 = Mat3::from_cols(
        Vec3::new(r_yz / r, 0.0, v.x / r),
        Vec3::new(-v.x * v.y / r / r_yz, v.z / r_yz, v.y / r),
        Vec3::new(-v.x * v.z / r / r_yz, -v.y / r_yz, v.z / r),
    );
    let rot2 = rot1.transpose();
    let m3 = Mat3::from_cols(
        Vec3::new(a.cos(), a.sin(), 0.0),
        Vec3::new(-a.sin(), a.cos(), 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    rot2 * m3 * rot1
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    if let Some(f) = args.values("outformat").first() {
        if f != "cnf" {
            return Err(format!("@outformat {f}: only cnf is supported"));
        }
    }
    let mut topo =
        build_topology(read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?);
    let traj = args.value("traj")?.to_string();
    let data = read_g96_labeled(&traj).map_err(|e| format!("@traj: {e}"))?;
    solvate_to_atoms(&mut topo, data.atoms.len()).map_err(|e| e.to_string())?;
    let pbc = Pbc::from_args(&args)?;
    let mobile = args
        .values("mobile")
        .first()
        .cloned()
        .unwrap_or_else(|| "last".into());
    if mobile != "first" && mobile != "last" {
        return Err(format!("@mobile argument '{mobile}' unknown.\n"));
    }
    let mut props = PropertyContainer::parse(args.values("prop"), &topo)?;
    let mut title = format!("gca has modified coordinates in {traj} such that");
    let (mut min, mut max, mut step) = (Vec::new(), Vec::new(), Vec::new());
    for p in &props.props {
        match p.args.len() {
            1 => {
                step.push(1.0);
                min.push(p.args[0]);
                max.push(p.args[0]);
                title.push_str(&format!(
                    "\n{}\tis set to {}",
                    p.title(&topo),
                    trim_float(p.args[0])
                ));
            },
            3 => {
                step.push(p.args[0]);
                min.push(p.args[1]);
                max.push(p.args[2]);
                title.push_str(&format!(
                    "\n{}\tgoes from {:>8} to {:>8} with a step size of {:>8}",
                    p.title(&topo),
                    trim_float(p.args[1]),
                    trim_float(p.args[2]),
                    trim_float(p.args[0])
                ));
            },
            _ => return Err("properties: specify single value or step%min%max values".into()),
        }
    }
    // every combination, first property varying fastest (gromos++'s odometer)
    let n = props.props.len();
    let mut combination: Vec<Vec<f64>> = Vec::new();
    let mut single = min.clone();
    let mut index = 0usize;
    let mut flag = true;
    while flag && single[n - 1] <= max[n - 1] {
        combination.push(single.clone());
        single[index] += step[index];
        while single[index] > max[index] {
            single[index] = min[index];
            index += 1;
            if index >= n {
                flag = false;
                break;
            }
            single[index] += step[index];
            if single[index] <= max[index] {
                index = 0;
            }
        }
    }
    let bonds: BTreeSet<(usize, usize)> = topo.all_bonds_global().map(|b| (b.i, b.j)).collect();
    let angles: BTreeSet<(usize, usize, usize)> =
        topo.all_angles_global().map(|a| (a.i, a.j, a.k)).collect();
    let dihedrals: BTreeSet<(usize, usize, usize, usize)> = topo
        .all_proper_dihedrals_global()
        .map(|d| (d.i, d.j, d.k, d.l))
        .chain(
            topo.all_improper_dihedrals_global()
                .map(|d| (d.i, d.j, d.k, d.l)),
        )
        .collect();
    let mol_of = |a: usize| topo.molecules.iter().position(|r| r.contains(&a));
    for p in &props.props {
        let a = &p.atoms;
        if a.iter().map(|&x| mol_of(x)).collect::<BTreeSet<_>>().len() != 1 {
            return Err(format!(
                " Covalent interactions are always within one molecule: {}",
                p.title(&topo)
            ));
        }
        let ok = match p.kind {
            Kind::Distance => bonds.contains(&(a[0].min(a[1]), a[0].max(a[1]))),
            Kind::Angle => {
                angles.contains(&(a[0], a[1], a[2])) || angles.contains(&(a[2], a[1], a[0]))
            },
            Kind::Torsion => {
                dihedrals.contains(&(a[0], a[1], a[2], a[3]))
                    || dihedrals.contains(&(a[3], a[2], a[1], a[0]))
            },
            _ => return Err(format!("cannot modify property of type {}", p.type_name())),
        };
        if !ok {
            let what = match p.kind {
                Kind::Distance => "Bond",
                Kind::Angle => "Angle",
                _ => "(improper) Dihedral",
            };
            return Err(format!(" {what} not found in topology: {}", p.title(&topo)));
        }
    }
    let neighbours = |a: usize| -> Vec<usize> {
        bonds
            .iter()
            .filter_map(|&(i, j)| {
                if i == a {
                    Some(j)
                } else if j == a {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    };
    let mut out_frames = String::new();
    let periodicity = pbc.periodicity(data.box_dims);
    let mut pos: Vec<Vec3> = data.atoms.iter().map(|a| a.pos).collect();
    pbc.gather(&topo, &mut pos, &periodicity);
    for (stepnum, c) in combination.iter().enumerate() {
        for (pi, p) in props.props.iter_mut().enumerate() {
            let target = c[pi];
            let value = p.calc(&pos, &periodicity);
            let atoms_moved = atoms_to_change(p, &mobile, &neighbours);
            match p.kind {
                Kind::Distance => {
                    let v = (pos[p.atoms[1]] - pos[p.atoms[0]]).normalize() * (target - value);
                    for &a in &atoms_moved {
                        pos[a] += v;
                    }
                },
                Kind::Angle => {
                    let v1 = pos[p.atoms[0]] - pos[p.atoms[1]];
                    let v2 = pos[p.atoms[2]] - pos[p.atoms[1]];
                    let rot = rotate_around(v1.cross(v2), target - value);
                    let centre = pos[p.atoms[1]];
                    for &a in &atoms_moved {
                        pos[a] = rot * (pos[a] - centre) + centre;
                    }
                },
                Kind::Torsion => {
                    let v = pos[p.atoms[2]] - pos[p.atoms[1]];
                    let rot = rotate_around(v, target - value);
                    let centre = pos[p.atoms[2]];
                    for &a in &atoms_moved {
                        pos[a] = rot * (pos[a] - centre) + centre;
                    }
                },
                _ => unreachable!(),
            }
        }
        let atoms: Vec<G96Atom> = data
            .atoms
            .iter()
            .zip(&pos)
            .map(|(a, p)| {
                let mut a = a.clone();
                a.pos = *p;
                a
            })
            .collect();
        let frame = format_g96(
            "",
            &atoms,
            (data.box_dims != Vec3::ZERO).then_some(data.box_dims),
        );
        // OutG96S frame: TIMESTEP then POSITION + GENBOX (the TITLE is written once)
        let body = frame
            .split_once(
                "END
",
            )
            .map(|x| x.1)
            .unwrap_or("")
            .to_string();
        out_frames.push_str(&format!(
            "TIMESTEP\n{:>18}{:>20.9}\n#if @time flag is used the value for step refers to the\n#step-th configuration in the original trajectory file\nEND\n",
            stepnum, stepnum as f64
        ));
        out_frames.push_str(&body);
    }
    write!(out, "TITLE\n{title}\nEND\n{out_frames}").map_err(io_err)?;
    Ok(())
}

/// gromos++ `atoms_to_change`: breadth-first over the bonds from the mobile end of the property,
/// never crossing the property's own atoms.
fn atoms_to_change(
    p: &Property,
    mobile: &str,
    neighbours: &dyn Fn(usize) -> Vec<usize>,
) -> Vec<usize> {
    let atoms = &p.atoms;
    let mut frontier = BTreeSet::new();
    if mobile == "last" {
        frontier.insert(*atoms.last().unwrap());
        if matches!(p.kind, Kind::Torsion) {
            frontier.insert(atoms[2]);
        }
    } else {
        frontier.insert(atoms[0]);
        if matches!(p.kind, Kind::Torsion) {
            frontier.insert(atoms[1]);
        }
    }
    let mut result: Vec<usize> = Vec::new();
    while !frontier.is_empty() {
        let mut next = BTreeSet::new();
        for &b in &frontier {
            result.push(b);
            for n in neighbours(b) {
                if !atoms.contains(&n) && !result.contains(&n) {
                    next.insert(n);
                }
            }
        }
        frontier = next;
    }
    result
}
