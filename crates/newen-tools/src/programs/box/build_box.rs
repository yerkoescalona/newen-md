//! build_box - Generate a condensed-phase system on a regular grid
//!
//! Usage: build_box @topo <topology> @pos <molecule.g96> @nsm <N> @dens <density>
//!
//! Replicates a single molecule on an N x N x N grid into a cubic box whose
//! size is chosen so that the resulting system matches the target density.

use crate::programs::{io_err, parse_clap};
use clap::Parser;
use newen::math::Vec3;
use newen_io::args::{Program, Status};
use newen_io::topology::read_topology_file;
use newen_io::{read_g96_labeled, G96Atom};
use std::io::Write;

pub const USAGE: &str = "build_box - Generate a condensed-phase system on a grid\n\nUsage: build_box @topo <topology> @pos <molecule.g96> @nsm <N> @dens <density>\n\nArguments:\n  @topo   Molecular topology file for a single molecule\n  @pos    Input coordinate file for a single molecule\n  @nsm    Number of molecules per dimension (total = nsm^3)\n  @dens   Density of the liquid (kg/m^3)";

pub const KNOWN: &[&str] = &["topo", "pos", "nsm", "dens"];

pub const PROGRAM: Program = Program {
    name: "build_box",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

/// Conversion factor from atomic mass units to kg per nm^3 at unit density
/// (newen `build_box.cc`: `vtot = nsm * (weight * 1.66056) / density`).
const AMU_TO_KG_PER_NM3: f64 = 1.66056;

#[derive(Parser)]
#[command(
    name = "build_box",
    version,
    about = "Generate a condensed-phase system on a grid"
)]
struct Args {
    /// Molecular topology file for a single molecule
    #[arg(long)]
    topo: String,

    /// Input coordinate file for a single molecule
    #[arg(long)]
    pos: String,

    /// Number of molecules per dimension (total = nsm^3)
    #[arg(long)]
    nsm: usize,

    /// Density of the liquid (kg/m^3)
    #[arg(long)]
    dens: f64,
}

/// Computes (box_side, cell_side) for `nsm` copies per dimension of a molecule
/// of the given `weight` (u) so that the system matches the target `density`
/// (kg/m^3). Mirrors newen `build_box.cc`: `vtot = nsm^3 * (weight * 1.66056) / density`.
fn box_dimensions(nsm: usize, weight: f64, density: f64) -> (f64, f64) {
    let nsm_total = nsm * nsm * nsm;
    let vtot = nsm_total as f64 * (weight * AMU_TO_KG_PER_NM3) / density;
    let box_side = vtot.cbrt();
    let cell = box_side / nsm as f64;
    (box_side, cell)
}

/// Mass-weighted center of geometry (center of mass)
fn center_of_mass(atoms: &[G96Atom], masses: &[f64]) -> Vec3 {
    let mut com = Vec3::ZERO;
    let mut total_mass = 0.0;
    for (atom, &mass) in atoms.iter().zip(masses) {
        com += atom.pos * mass;
        total_mass += mass;
    }
    com / total_mass
}

/// A molecule replicated on a grid, with the box that gives the requested density.
#[derive(Debug, Clone)]
pub struct BuildBoxResult {
    pub atoms: Vec<G96Atom>,
    /// Cubic box edge (nm).
    pub box_side: f64,
    /// Total number of copies placed (`nsm³`).
    pub n_molecules: usize,
    /// Mass of one molecule (u), as summed from the topology.
    pub weight: f64,
}

/// Replicate `molecule` on an `nsm × nsm × nsm` grid at the target `density` (kg/m³).
///
/// `masses` is the molecule's own topology masses — their sum is the molecular weight the box
/// volume follows from. Each copy is placed with its centre of mass at the centre of its cell.
pub fn build_box(
    molecule: &[G96Atom],
    masses: &[f64],
    nsm: usize,
    density: f64,
) -> Result<BuildBoxResult, String> {
    if nsm == 0 {
        return Err("Error: @nsm must be a positive integer".to_string());
    }
    if molecule.len() != masses.len() {
        return Err(format!(
            "Error: topology has {} atoms but the coordinate file has {}",
            masses.len(),
            molecule.len()
        ));
    }
    let weight: f64 = masses.iter().sum();
    let n_molecules = nsm * nsm * nsm;
    let (box_side, cell) = box_dimensions(nsm, weight, density);

    // Center the molecule's center of mass in the first grid cell
    let cell_center = Vec3::new(cell, cell, cell) * 0.5;
    let com = center_of_mass(molecule, masses);
    let mut centered_atoms = molecule.to_vec();
    for atom in &mut centered_atoms {
        atom.pos = atom.pos - com + cell_center;
    }

    // Replicate on an nsm x nsm x nsm grid
    let mut atoms: Vec<G96Atom> = Vec::with_capacity(centered_atoms.len() * n_molecules);
    let mut atom_num = 1usize;
    let mut res_num = 1usize;
    for i in 0..nsm {
        for j in 0..nsm {
            for k in 0..nsm {
                let shift = Vec3::new(i as f64 * cell, j as f64 * cell, k as f64 * cell);
                for atom in &centered_atoms {
                    atoms.push(G96Atom {
                        res_num,
                        res_name: atom.res_name.clone(),
                        atom_name: atom.atom_name.clone(),
                        atom_num,
                        pos: atom.pos + shift,
                    });
                    atom_num += 1;
                }
                res_num += 1;
            }
        }
    }

    Ok(BuildBoxResult {
        atoms,
        box_side,
        n_molecules,
        weight,
    })
}

/// build_box's own GROMOS96 layout (a `BOX` block, and its own POSITION column widths — *not*
/// `format_g96`'s; changing that is a separate, byte-visible decision).
pub fn format(result: &BuildBoxResult, pos_file: &str, density: f64) -> String {
    let mut s = String::new();
    s.push_str("TITLE\n");
    s.push_str(&format!(
        "build_box: {} copies of {}\nDensity : {} kg/m^3\tMolecular weight : {:.4} u\n",
        result.n_molecules, pos_file, density, result.weight
    ));
    s.push_str("END\nPOSITION\n");
    for atom in &result.atoms {
        s.push_str(&format!(
            "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}\n",
            atom.res_num,
            atom.res_name,
            atom.atom_name,
            atom.atom_num,
            atom.pos.x,
            atom.pos.y,
            atom.pos.z
        ));
    }
    s.push_str("END\nBOX\n");
    s.push_str(&format!(
        "{:15.9}{:15.9}{:15.9}\n",
        result.box_side, result.box_side, result.box_side
    ));
    s.push_str("END\n");
    s
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let Some(args) = parse_clap::<Args>(argv, out)? else {
        return Ok(());
    };

    // Read topology to get the molecule's total mass
    let parsed_topo = read_topology_file(&args.topo)
        .map_err(|e| format!("Error reading topology '{}': {e:?}", args.topo))?;

    // Read the single-molecule coordinates
    let mol_atoms = read_g96_labeled(&args.pos)
        .map_err(|e| format!("Error reading '{}': {e}", args.pos))?
        .atoms;

    let result = build_box(&mol_atoms, &parsed_topo.masses, args.nsm, args.dens)?;
    let cell = result.box_side / args.nsm as f64;

    eprintln!("build_box - Building condensed-phase system");
    eprintln!("  Topology:        {}", args.topo);
    eprintln!("  Molecule:        {}", args.pos);
    eprintln!("  Copies per dim:  {}", args.nsm);
    eprintln!("  Total copies:    {}", result.n_molecules);
    eprintln!("  Molecular mass:  {:.4} u", result.weight);
    eprintln!("  Target density:  {} kg/m^3", args.dens);
    eprintln!("  Box side:        {:.6} nm", result.box_side);
    eprintln!("  Cell side:       {:.6} nm", cell);
    eprintln!("  Total atoms:     {}", result.atoms.len());
    eprintln!();

    write!(out, "{}", format(&result, &args.pos, args.dens)).map_err(io_err)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_dimensions_match_target_density() {
        // 27 SPC waters (M = 18.0154 u) at 1000 kg/m^3
        let (box_side, cell) = box_dimensions(3, 18.0154, 1000.0);
        let nsm_total = 27.0;
        let vtot = nsm_total * (18.0154 * AMU_TO_KG_PER_NM3) / 1000.0;
        assert!((box_side - vtot.cbrt()).abs() < 1e-12);
        assert!((cell * 3.0 - box_side).abs() < 1e-12);
    }

    #[test]
    fn center_of_mass_is_mass_weighted() {
        let masses = vec![16.0, 1.0, 1.0];
        let atoms = vec![
            G96Atom {
                res_num: 1,
                res_name: "SOL".into(),
                atom_name: "OW".into(),
                atom_num: 1,
                pos: Vec3::new(0.0, 0.0, 0.0),
            },
            G96Atom {
                res_num: 1,
                res_name: "SOL".into(),
                atom_name: "HW1".into(),
                atom_num: 2,
                pos: Vec3::new(1.0, 0.0, 0.0),
            },
            G96Atom {
                res_num: 1,
                res_name: "SOL".into(),
                atom_name: "HW2".into(),
                atom_num: 3,
                pos: Vec3::new(0.0, 1.0, 0.0),
            },
        ];
        let com = center_of_mass(&atoms, &masses);
        let total_mass: f64 = masses.iter().sum();
        assert!((com.x - 1.0 / total_mass).abs() < 1e-12);
        assert!((com.y - 1.0 / total_mass).abs() < 1e-12);
        assert!((com.z - 0.0).abs() < 1e-12);
    }
}
