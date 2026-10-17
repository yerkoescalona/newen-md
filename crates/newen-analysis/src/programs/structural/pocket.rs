//! pocket — binding-site vectors from a pocket centre out to the protein surface (gromos++
//! `pocket`).
//!
//! Usage: pocket @topo <top> @pbc <r> @center <cart(x,y,z)|atom(<atom>)> @protein <atoms>
//!        [@reject <atoms>] @radius <nm> @vec_number_factor <n> [@radH <nm>] [@hemisphere]
//!        [@volume_and_area] [@final_vector_coords] @traj <trc…>
//!
//! Vectors from the centre to the vertices of an icosahedral sphere mesh of radius `@radius`
//! (`@vec_number_factor` subdivisions per face) are truncated where they first enter a protein
//! atom's van der Waals sphere (radius ½(2C12/C6)^⅙ from the topology, `@radH` for solute atoms
//! without C6). Per frame `lengths.txt` and `charges.txt` (the truncating atom's charge) are
//! written; `@volume_and_area` adds `volume.txt` and `area.txt` of the enclosed mesh and
//! `@final_vector_coords` the truncated vectors (`vector_coords.txt`). Frames are not gathered:
//! gromos++ parses `@pbc` but never applies it. An `atom(…)` centre is evaluated on the first
//! frame — gromos++ evaluates it before any frame is read, which gives the origin.

use crate::args::Arguments;
use crate::distribution::cpp_g;
use crate::pbc::Pbc;
use crate::programs::io_err;
use crate::property::{ordered_atoms, parse_vector, VectorSpec};
use newen_core::math::Vec3;
use newen_core::Topology;
use newen_io::args::{Program, Status};
use newen_io::topology::{build_topology, read_topology_file, solvate_to_atoms};
use newen_io::trajectory::TrajectoryReader;
use std::collections::HashMap;
use std::io::Write;

pub const USAGE: &str = "# pocket
\t@topo                 <molecular topology file>
\t@pbc                    <boundary type> [<gathermethod>]
\t@center                 <pocket center>
\t@protein                <protein atoms to consider for vector generation>
\t@[reject                <protein atoms to discard for vector generation>]
\t@radius                 <max. length of the generated binding site vectors (nm)>
\t@vec_number_factor      <factor that decides the number of binding site vectors>
\t[@radH                  <radius to be used for hydrogen atoms>]
\t[@hemisphere            <keep only the z>0 hemisphere initial vectors>]
\t[@volume_and_area       <compute enclosed volume and surface area>]
\t[@final_vector_coords    <coordinates of the truncated vectors>]
\t@traj                   <trajectory files>";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "center",
    "protein",
    "reject",
    "radius",
    "vec_number_factor",
    "radH",
    "hemisphere",
    "volume_and_area",
    "final_vector_coords",
    "traj",
];

pub const PROGRAM: Program = Program {
    name: "pocket",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

struct Mesh {
    vertices: Vec<Vec3>,
    triangles: Vec<[usize; 3]>,
}

fn normalize(v: Vec3) -> Vec3 {
    let n = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
    Vec3::new(v.x / n, v.y / n, v.z / n)
}

fn round8(x: f64) -> f64 {
    let f = 1e8;
    (x * f).round() / f + 0.0
}

/// gromos++ `generate_sphere_mesh`: an icosahedron with every face subdivided `n` times and
/// projected onto the sphere; vertices are shared on their coordinates rounded to 8 decimals.
fn sphere_mesh(radius: f64, subdivisions: usize) -> Mesh {
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    let ico: Vec<Vec3> = [
        (-1.0, phi, 0.0),
        (1.0, phi, 0.0),
        (-1.0, -phi, 0.0),
        (1.0, -phi, 0.0),
        (0.0, -1.0, phi),
        (0.0, 1.0, phi),
        (0.0, -1.0, -phi),
        (0.0, 1.0, -phi),
        (phi, 0.0, -1.0),
        (phi, 0.0, 1.0),
        (-phi, 0.0, -1.0),
        (-phi, 0.0, 1.0),
    ]
    .iter()
    .map(|&(x, y, z)| normalize(Vec3::new(x, y, z)))
    .collect();
    let faces: [[usize; 3]; 20] = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    let n = subdivisions.max(1);
    let mut map: HashMap<(u64, u64, u64), usize> = HashMap::new();
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for face in &faces {
        let (v1, v2, v3) = (ico[face[0]], ico[face[1]], ico[face[2]]);
        let mut grid: Vec<Vec<usize>> = (0..=n).map(|i| vec![0; n + 1 - i]).collect();
        for i in 0..=n {
            for j in 0..=(n - i) {
                let u = i as f64 / n as f64;
                let v = j as f64 / n as f64;
                let w = 1.0 - u - v;
                let p = normalize(Vec3::new(
                    u * v1.x + v * v2.x + w * v3.x,
                    u * v1.y + v * v2.y + w * v3.y,
                    u * v1.z + v * v2.z + w * v3.z,
                ));
                let rp = Vec3::new(round8(p.x), round8(p.y), round8(p.z));
                let key = (rp.x.to_bits(), rp.y.to_bits(), rp.z.to_bits());
                let idx = *map.entry(key).or_insert_with(|| {
                    vertices.push(Vec3::new(rp.x * radius, rp.y * radius, rp.z * radius));
                    vertices.len() - 1
                });
                grid[i][j] = idx;
            }
        }
        for i in 0..n {
            for j in 0..(n - i) {
                let (a, b, c) = (grid[i][j], grid[i + 1][j], grid[i][j + 1]);
                triangles.push([a, b, c]);
                if j + i + 1 < n {
                    let d = grid[i + 1][j + 1];
                    triangles.push([b, d, c]);
                }
            }
        }
    }
    Mesh {
        vertices,
        triangles,
    }
}

fn triangle_area(v1: Vec3, v2: Vec3, v3: Vec3) -> f64 {
    let a = Vec3::new(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
    let b = Vec3::new(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);
    let cr = Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    );
    0.5 * (cr.x * cr.x + cr.y * cr.y + cr.z * cr.z).sqrt()
}

fn pyramid_volume(v1: Vec3, v2: Vec3, v3: Vec3) -> f64 {
    ((v1.x * (v2.y * v3.z - v2.z * v3.y) - v1.y * (v2.x * v3.z - v2.z * v3.x)
        + v1.z * (v2.x * v3.y - v2.y * v3.x))
        / 6.0)
        .abs()
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
    if let Some(first) = args.values("traj").first() {
        if let Ok(Some(frame)) = TrajectoryReader::new(first)
            .map_err(|e| format!("@traj {first}: {e}"))?
            .read_frame()
        {
            solvate_to_atoms(&mut topo, frame.positions.len()).map_err(|e| e.to_string())?;
        }
    }
    let pbc = Pbc::from_args(&args)?;
    let use_hemisphere = args.count("hemisphere") >= 0;
    let volume_and_area = args.count("volume_and_area") >= 0;
    let final_vector_coords = args.count("final_vector_coords") >= 0;
    let center_spec = parse_vector(args.value("center")?, &topo)?;
    let mut protein = selection(args.values("protein"), &topo)?;
    let reject = selection(args.values("reject"), &topo)?;
    protein.retain(|a| !reject.contains(a));
    let radius: f64 = args.require("radius")?;
    let vec_number_factor: usize = args.get("vec_number_factor", 0)?;
    let rad_h: f64 = args.get("radH", 0.11)?;
    let n_solute = topo.num_solute_atoms();
    // gromos++ compute_atomic_radii_vdw, then radH for solute atoms left without a radius
    let radii: Vec<f64> = protein
        .iter()
        .map(|&i| {
            let lj = topo.lj_parameter(topo.iac[i], topo.iac[i]);
            let r = if lj.c6 >= 1.0e-20 {
                0.5 * (((2.0 * lj.c12) / lj.c6).ln() / 6.0).exp()
            } else {
                0.0
            };
            if r == 0.0 && i < n_solute {
                rad_h
            } else {
                r
            }
        })
        .collect();
    let charges: Vec<f64> = protein.iter().map(|&i| topo.charge[i]).collect();
    let mesh = sphere_mesh(radius, vec_number_factor);
    let mut points = mesh.vertices;
    let mut triangles = mesh.triangles;
    let mut stdout = format!(
        "Generated {} points and {} triangles.\n",
        points.len(),
        triangles.len()
    );
    if use_hemisphere {
        let mut old_to_new = vec![None; points.len()];
        let mut kept = Vec::new();
        for (i, p) in points.iter().enumerate() {
            if p.z >= 0.0 {
                old_to_new[i] = Some(kept.len());
                kept.push(*p);
            }
        }
        triangles = triangles
            .iter()
            .filter_map(|t| Some([old_to_new[t[0]]?, old_to_new[t[1]]?, old_to_new[t[2]]?]))
            .collect();
        points = kept;
        stdout.push_str(&format!(
            "Hemisphere flag active: kept {} points and {} triangles above xy-plane.\n",
            points.len(),
            triangles.len()
        ));
    }
    let mut charges_out = String::new();
    let mut lengths_out = String::new();
    let mut volume_out = String::new();
    let mut area_out = String::new();
    let mut coords_out = String::new();
    let mut center: Option<Vec3> = None;
    let mut snapshot = 0usize;
    for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let c = *center.get_or_insert_with(|| match &center_spec {
                VectorSpec::Cart(v) => *v,
                other => other.eval(&frame.positions, &pbc.periodicity(frame.box_dims)),
            });
            let pos: Vec<Vec3> = frame.positions.iter().map(|p| *p + (-1.0 * c)).collect();
            let mut dists = Vec::with_capacity(points.len());
            let mut qs = Vec::with_capacity(points.len());
            let mut finals = Vec::with_capacity(points.len());
            for &v in &points {
                let v_abs = v.length();
                let mut min_dist = f64::INFINITY;
                let mut t1 = f64::INFINITY;
                let mut final_vector = v;
                let mut selected_charge = 0.0;
                for (k, &a) in protein.iter().enumerate() {
                    let s = pos[a];
                    let r = radii[k];
                    let d = s.length();
                    if d < 1e-6 {
                        continue;
                    }
                    let cos_alpha = s.dot(v) / (d * v_abs);
                    let mut touching = false;
                    let mut temp = Vec3::ZERO;
                    if s == v {
                        t1 = d - r;
                        touching = true;
                        temp = v / v_abs * t1;
                    } else if (0.0..=1.0).contains(&cos_alpha) {
                        let y = d * (1.0 - cos_alpha * cos_alpha).sqrt();
                        if y < r {
                            let x = (r * r - y * y).sqrt();
                            let t = s.dot(v).abs() / v_abs;
                            t1 = t - x;
                            let used = if t1 < radius { t1 } else { radius };
                            temp = v / v_abs * used;
                            touching = temp.z >= 0.0;
                        } else if (y - r).abs() < 1e-6 {
                            t1 = s.dot(v) / v_abs;
                            let used = if t1 < radius { t1 } else { radius };
                            temp = v / v_abs * used;
                            touching = true;
                        }
                    }
                    if touching && t1 < min_dist {
                        min_dist = t1;
                        final_vector = temp;
                        selected_charge = charges[k];
                    }
                }
                if !min_dist.is_finite() || min_dist > radius {
                    min_dist = radius;
                    final_vector = v / v_abs * radius;
                }
                dists.push((min_dist * 1000.0).round() / 1000.0);
                qs.push(selected_charge);
                finals.push(final_vector);
            }
            if volume_and_area {
                let (mut area, mut volume) = (0.0, 0.0);
                for t in &triangles {
                    area += triangle_area(finals[t[0]], finals[t[1]], finals[t[2]]);
                    volume += pyramid_volume(finals[t[0]], finals[t[1]], finals[t[2]]);
                }
                area_out.push_str(&format!("SNAPSHOT_{snapshot:05} {area:.6}\n"));
                volume_out.push_str(&format!("SNAPSHOT_{snapshot:05} {volume:.6}\n"));
            }
            charges_out.push_str(&format!("SNAPSHOT_{snapshot}"));
            lengths_out.push_str(&format!("SNAPSHOT_{snapshot}"));
            if final_vector_coords {
                coords_out.push_str(&format!("SNAPSHOT_{snapshot}\n"));
            }
            for i in 0..qs.len() {
                charges_out.push_str(&format!(" {}", cpp_g(qs[i], 6)));
                lengths_out.push_str(&format!(" {}", cpp_g(dists[i], 6)));
                if final_vector_coords {
                    coords_out.push_str(&format!(
                        "{:.6};{:.6};{:.6}\n",
                        finals[i].x, finals[i].y, finals[i].z
                    ));
                }
            }
            charges_out.push('\n');
            lengths_out.push('\n');
            snapshot += 1;
        }
    }
    let write =
        |name: &str, text: &str| std::fs::write(name, text).map_err(|e| format!("{name}: {e}"));
    write("charges.txt", &charges_out)?;
    write("lengths.txt", &lengths_out)?;
    if volume_and_area {
        write("volume.txt", &volume_out)?;
        write("area.txt", &area_out)?;
    }
    if final_vector_coords {
        write("vector_coords.txt", &coords_out)?;
    }
    write!(out, "{stdout}").map_err(io_err)?;
    Ok(())
}
