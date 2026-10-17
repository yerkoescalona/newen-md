//! sim_box - Solvate a solute in a box of pre-equilibrated solvent
//!
//! Usage: sim_box @topo <topology> @pbc <r|t> @pos <solute.cnf> @solvent <solvent.g96>
//!               [@minwall <distance>] [@thresh <distance>] [@boxsize]
//!
//! Solvates a molecular system the way gromos++ `sim_box` does (`programs/sim_box.cc`):
//! 1. box from the longest solute atom–atom distance plus twice `@minwall` (cubic with one
//!    value; with three values the x/xy/xyz maximum distances, which assumes a solute rotated
//!    with its longest axis along z), or the solute file's box with `@boxsize`;
//! 2. solute shifted to its centre of geometry, the solvent template centred likewise and
//!    replicated `int(box / template) + 1` times per dimension around the origin;
//! 3. a solvent molecule is kept when its centre of geometry lies inside the box and is farther
//!    than `@thresh` from every non-hydrogen solute atom (hydrogens: mass 1.008 in `@topo`);
//! 4. the solvated system is written as GROMOS96 with a GENBOX block.
//!
//! Reference: `tests/sim_box_reference.rs` reproduces gromos++'s output for methanol in SPC.

use crate::programs::io_err;
use newen::math::Vec3;
use newen_io::args::{Arguments, Program, Status};
use newen_io::coordinate::format_g96;
use newen_io::{read_g96_labeled, G96Atom};
use std::io::Write;

pub const USAGE: &str = "sim_box - Solvate a solute in a box of pre-equilibrated solvent\n\nUsage: sim_box @topo <topology> @pbc <r|t> @pos <solute.cnf> @solvent <solvent.g96> \\\n               [@minwall <distance>] [@thresh <distance>] [@boxsize]\n\nArguments:\n  @topo      Molecular topology file\n  @pbc       Periodic boundary conditions:\n             r = rectangular box\n             t = truncated octahedron (not yet implemented)\n  @pos       Input coordinate file for the solute\n  @solvent   Input coordinate file for pre-equilibrated solvent\n  @minwall   Minimum solute-to-wall distance (nm)\n             - 1 value: cubic box\n             - 3 values: rectangular box (x y z)\n  @thresh    Minimum solvent-to-solute distance (default: 0.23 nm)\n  @boxsize   Use box dimensions from solute coordinate file\n\nDescription:\n  Solvates a solute molecule in a pre-equilibrated solvent box.\n  The solvent configuration should contain a BOX block.\n\n  Box sizing options:\n  1. @boxsize: Use box from solute file\n  2. @minwall: Calculate box from solute extent + minwall distance\n\n  Solvent removal criteria:\n  - Solvent molecule center of geometry must be:\n    a) Inside the target box\n    b) At least @thresh distance from any solute atom\n\nExamples:\n  # Cubic box with 1.4 nm minimum wall distance\n  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @minwall 1.4\n\n  # Use box size from solute file\n  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @boxsize\n\n  # Rectangular box with different wall distances\n  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @minwall 1.0 1.5 2.0\n\nOutput:\n  GROMOS96 format to stdout containing solvated system";

pub const KNOWN: &[&str] = &[
    "topo", "pbc", "pos", "solvent", "minwall", "thresh", "boxsize",
];

pub const PROGRAM: Program = Program {
    name: "sim_box",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// How to size the box and how close solvent may come to the solute.
#[derive(Debug, Clone)]
pub struct SimBoxOptions {
    /// Minimum solute-to-wall distance (nm): one value cubic, three rectangular, none with
    /// `use_boxsize`.
    pub minwall: Vec<f64>,
    /// Minimum solvent-to-solute distance (nm).
    pub thresh: f64,
    /// Take the box from the solute's own coordinate file instead of `minwall`.
    pub use_boxsize: bool,
}

impl Default for SimBoxOptions {
    fn default() -> Self {
        Self {
            minwall: Vec::new(),
            thresh: 0.23, // Default from GROMOS++
            use_boxsize: false,
        }
    }
}

/// A solvated system: solute as read, then the solvent that survived the filter.
#[derive(Debug, Clone)]
pub struct SolvateResult {
    /// Solute atoms followed by the kept solvent atoms, solvent renumbered per molecule as `SOLV`.
    pub atoms: Vec<G96Atom>,
    pub box_dims: Vec3,
    pub n_solvent_molecules: usize,
    pub atoms_per_molecule: usize,
}

#[derive(Debug)]
struct SimBoxArgs {
    topo: String,
    pbc: char,
    pos: String,
    solvent: String,
    opts: SimBoxOptions,
}

fn parse_args(args: &Arguments) -> Result<SimBoxArgs, String> {
    let mut sb_args = SimBoxArgs {
        topo: args.value("topo")?.to_string(),
        pbc: 'r',
        pos: args.value("pos")?.to_string(),
        solvent: args.value("solvent")?.to_string(),
        opts: SimBoxOptions {
            thresh: args.get("thresh", 0.23)?,
            use_boxsize: args.has("boxsize"),
            ..Default::default()
        },
    };

    if args.has("pbc") {
        let pbc = args.value("pbc")?;
        let mut chars = pbc.chars();
        match (chars.next(), chars.next()) {
            (Some(c), None) => sb_args.pbc = c,
            _ => return Err(format!("Invalid @pbc: {pbc}")),
        }
        if sb_args.pbc != 'r' && sb_args.pbc != 't' {
            return Err(format!("Invalid @pbc: {pbc} (must be 'r' or 't')"));
        }
    }

    for v in args.values("minwall") {
        sb_args.opts.minwall.push(
            v.parse()
                .map_err(|_| format!("Invalid @minwall value: {v}"))?,
        );
    }

    Ok(sb_args)
}

/// Calculate center of geometry
fn calc_cog(atoms: &[G96Atom]) -> Vec3 {
    if atoms.is_empty() {
        return Vec3::ZERO;
    }
    let sum = atoms.iter().fold(Vec3::ZERO, |acc, a| acc + a.pos);
    sum / atoms.len() as f64
}

/// Shift all atoms by a vector
fn shift_atoms(atoms: &mut [G96Atom], shift: Vec3) {
    for atom in atoms.iter_mut() {
        atom.pos += shift;
    }
}

/// Largest atom–atom distance over the first `dim` dimensions, with the two atoms.
fn calc_max_size(atoms: &[G96Atom], dim: usize) -> (f64, usize, usize) {
    let mut max2 = 0.0;
    let mut pair = (0, 0);
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let d = atoms[i].pos - atoms[j].pos;
            let d2 = match dim {
                1 => d.x * d.x,
                2 => d.x * d.x + d.y * d.y,
                _ => d.length_squared(),
            };
            if d2 > max2 {
                max2 = d2;
                pair = (i, j);
            }
        }
    }
    (max2.sqrt(), pair.0, pair.1)
}

/// gromos++ treats a mass of 1.008 as hydrogen.
fn is_hydrogen(mass: f64) -> bool {
    (mass - 1.008).abs() < 1e-3
}

/// Nearest image of `r` in a rectangular box centred on the origin.
fn nearest_image(r: Vec3, box_dims: Vec3) -> Vec3 {
    Vec3::new(
        r.x - box_dims.x * (r.x / box_dims.x).round(),
        r.y - box_dims.y * (r.y / box_dims.y).round(),
        r.z - box_dims.z * (r.z / box_dims.z).round(),
    )
}

/// Solvate `solute` in a box of `solvent`, gromos++ `sim_box`'s algorithm.
///
/// `solute_masses` only decides which solute atoms the clash filter ignores (hydrogens), so it
/// may be the topology's masses truncated to the solute. `solvent_box` is the template's own box;
/// `solute_box` is needed only for [`SimBoxOptions::use_boxsize`].
pub fn solvate(
    solute: &[G96Atom],
    solute_masses: &[f64],
    solute_box: Option<Vec3>,
    solvent: &[G96Atom],
    solvent_box: Vec3,
    opts: &SimBoxOptions,
) -> Result<SolvateResult, String> {
    if !opts.minwall.is_empty() && opts.minwall.len() != 1 && opts.minwall.len() != 3 {
        return Err(format!(
            "@minwall requires 1 or 3 values, got {}",
            opts.minwall.len()
        ));
    }
    if !opts.use_boxsize && opts.minwall.is_empty() {
        return Err("Either @boxsize or @minwall must be specified".to_string());
    }
    if opts.use_boxsize && !opts.minwall.is_empty() {
        return Err("Cannot specify both @boxsize and @minwall".to_string());
    }
    if solute_masses.len() < solute.len() {
        return Err(format!(
            "Error: topology has {} solute atoms, coordinate file {}",
            solute_masses.len(),
            solute.len()
        ));
    }
    let solute_is_h: Vec<bool> = solute_masses[..solute.len()]
        .iter()
        .map(|&m| is_hydrogen(m))
        .collect();

    // Determine target box size
    let target_box: Vec3 = if opts.use_boxsize {
        solute_box.ok_or("Error: @boxsize specified but solute file has no BOX block")?
    } else if opts.minwall.len() == 1 {
        // gromos++: box = longest solute atom-atom distance + 2 * minwall
        let (max_dist, _, _) = calc_max_size(solute, 3);
        let size = max_dist + 2.0 * opts.minwall[0];
        Vec3::new(size, size, size)
    } else {
        // three minwall values: gromos++ expects the solute rotated (longest axis along z)
        // and takes the x, xy and xyz maximum distances for K, L and M respectively
        let (dx, _, _) = calc_max_size(solute, 1);
        let (dxy, _, _) = calc_max_size(solute, 2);
        let (dxyz, _, _) = calc_max_size(solute, 3);
        Vec3::new(
            dx + 2.0 * opts.minwall[0],
            dxy + 2.0 * opts.minwall[1],
            dxyz + 2.0 * opts.minwall[2],
        )
    };

    // Move solute to center
    let mut solute_atoms = solute.to_vec();
    let solute_cog = calc_cog(&solute_atoms);
    shift_atoms(&mut solute_atoms, Vec3::ZERO - solute_cog);

    // Centre the solvent template on its centre of geometry (gromos++ does the same).
    let mut template = solvent.to_vec();
    let solvent_cog = calc_cog(&template);
    shift_atoms(&mut template, Vec3::ZERO - solvent_cog);

    // Number of template copies per dimension: gromos++ `int(box / template) + 1`, laid out
    // symmetrically around the origin.
    let needed = [
        (target_box.x / solvent_box.x) as i32 + 1,
        (target_box.y / solvent_box.y) as i32 + 1,
        (target_box.z / solvent_box.z) as i32 + 1,
    ];
    let start = Vec3::new(
        -0.5 * (needed[0] - 1) as f64 * solvent_box.x,
        -0.5 * (needed[1] - 1) as f64 * solvent_box.y,
        -0.5 * (needed[2] - 1) as f64 * solvent_box.z,
    );
    let mut all_solvent: Vec<G96Atom> = Vec::new();
    for ix in 0..needed[0] {
        for iy in 0..needed[1] {
            for iz in 0..needed[2] {
                let shift = start
                    + Vec3::new(
                        ix as f64 * solvent_box.x,
                        iy as f64 * solvent_box.y,
                        iz as f64 * solvent_box.z,
                    );
                for atom in &template {
                    let mut shifted = atom.clone();
                    shifted.pos = atom.pos + shift;
                    all_solvent.push(shifted);
                }
            }
        }
    }

    // Atoms per solvent molecule from the template's residue numbering (3 for water).
    let atoms_per_molecule = {
        let first = template.first().map(|a| a.res_num);
        let n = template
            .iter()
            .take_while(|a| Some(a.res_num) == first)
            .count();
        n.max(1)
    };
    let num_molecules = all_solvent.len() / atoms_per_molecule;

    // Keep a molecule when its centre of geometry is inside the box (nearest image to the
    // origin is itself) and farther than @thresh from every non-hydrogen solute atom.
    let mut kept_solvent: Vec<G96Atom> = Vec::new();
    let thresh2 = opts.thresh * opts.thresh;
    let min_init = target_box.length_squared();
    for mol_idx in 0..num_molecules {
        let mol_atoms =
            &all_solvent[mol_idx * atoms_per_molecule..(mol_idx + 1) * atoms_per_molecule];
        let cog = calc_cog(mol_atoms);
        let check = nearest_image(cog, target_box);
        if check != cog {
            continue; // outside the box
        }
        let mut min2 = min_init;
        for (atom, &is_h) in solute_atoms.iter().zip(&solute_is_h) {
            if is_h {
                continue;
            }
            let d2 = (check - atom.pos).length_squared();
            if d2 < min2 {
                min2 = d2;
            }
        }
        if min2 > thresh2 {
            kept_solvent.extend(mol_atoms.iter().cloned());
        }
    }
    let n_solvent_molecules = kept_solvent.len() / atoms_per_molecule;

    // gromos++ OutG96S layout: solute as read, solvent residues SOLV numbered per molecule from 1.
    let mut atoms = solute_atoms;
    for (i, atom) in kept_solvent.into_iter().enumerate() {
        atoms.push(G96Atom {
            res_num: i / atoms_per_molecule + 1,
            res_name: "SOLV".to_string(),
            ..atom
        });
    }

    Ok(SolvateResult {
        atoms,
        box_dims: target_box,
        n_solvent_molecules,
        atoms_per_molecule,
    })
}

/// The TITLE gromos++ writes for a solvated system.
pub fn title(pos_file: &str, solvent_file: &str, n_solvent_molecules: usize) -> String {
    format!("Solvating {pos_file} in {solvent_file}\nAdded {n_solvent_molecules} solvent molecules")
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let sb_args = parse_args(&args)?;

    if sb_args.pbc == 't' {
        return Err("Error: Truncated octahedron (@pbc t) not yet implemented\n       Use rectangular boxes (@pbc r) for now".to_string());
    }

    eprintln!("sim_box - Solvating system");
    eprintln!("  Topology:   {}", sb_args.topo);
    eprintln!("  PBC:        {}", sb_args.pbc);
    eprintln!("  Solute:     {}", sb_args.pos);
    eprintln!("  Solvent:    {}", sb_args.solvent);
    eprintln!("  Threshold:  {} nm", sb_args.opts.thresh);
    eprintln!();

    eprintln!("Reading solute coordinates...");
    let solute =
        read_g96_labeled(&sb_args.pos).map_err(|e| format!("Error reading solute: {e}"))?;
    let solute_box = (solute.box_dims != Vec3::ZERO).then_some(solute.box_dims);
    eprintln!("  Solute: {} atoms", solute.atoms.len());

    // Solute masses: the clash filter ignores hydrogens, as gromos++ does.
    let parsed = newen_io::topology::read_topology_file(&sb_args.topo)
        .map_err(|e| format!("Error reading topology: {e}"))?;

    eprintln!("Reading solvent coordinates...");
    let solvent =
        read_g96_labeled(&sb_args.solvent).map_err(|e| format!("Error reading solvent: {e}"))?;
    if solvent.box_dims == Vec3::ZERO {
        return Err("Error: Solvent file must contain a BOX block".to_string());
    }
    eprintln!(
        "  Solvent: {} atoms in box ({:.3}, {:.3}, {:.3}) nm",
        solvent.atoms.len(),
        solvent.box_dims.x,
        solvent.box_dims.y,
        solvent.box_dims.z
    );

    let result = solvate(
        &solute.atoms,
        &parsed.masses,
        solute_box,
        &solvent.atoms,
        solvent.box_dims,
        &sb_args.opts,
    )?;

    let n_solute = solute.atoms.len();
    eprintln!();
    eprintln!("Final system:");
    eprintln!(
        "  {} solute atoms + {} solvent atoms = {}",
        n_solute,
        result.atoms.len() - n_solute,
        result.atoms.len()
    );
    eprintln!(
        "  Box: ({:.6}, {:.6}, {:.6}) nm",
        result.box_dims.x, result.box_dims.y, result.box_dims.z
    );
    eprintln!();

    let title = title(&sb_args.pos, &sb_args.solvent, result.n_solvent_molecules);
    write!(
        out,
        "{}",
        format_g96(&title, &result.atoms, Some(result.box_dims))
    )
    .map_err(io_err)?;
    Ok(())
}
