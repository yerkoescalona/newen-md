//! ion - Replace solvent molecules by ions
//!
//! Usage: ion @topo <top> @pbc <type> @pos <cnf> @positive <N> <name> [@negative <N> <name>]
//!            @potential <cutoff> [@mindist <d>]
//!    or: ion @f <argfile>
//!
//! Places counter-ions by replacing solvent molecules at positions with the
//! lowest (for positive ions) or highest (for negative ions) electrostatic
//! potential from solute charges.

use crate::programs::{io_err, parse_clap};
use clap::Parser;
use newen_io::args::{Program, Status};
use newen_io::topology::{read_topology_file, ParsedTopology};
use newen_io::{read_g96_labeled, G96Atom, LabeledCoordinateData};
use std::collections::HashSet;
use std::io::Write;

pub const USAGE: &str = "ion - Replace solvent molecules by ions\n\nUsage: ion @topo <top> @pbc <type> @pos <cnf> @positive <N> <name> [@negative <N> <name>]\n           @potential <cutoff> [@mindist <d>]\n   or: ion @f <argfile>\n\nArguments:\n  @topo       Molecular topology file\n  @pbc        Periodic boundary conditions (r = rectangular)\n  @pos        Input coordinate file (solvated system)\n  @positive   Place positive ions: <number> <ion_name>\n  @negative   Place negative ions: <number> <ion_name>\n  @potential  Cutoff for the electrostatic potential (nm)\n  @mindist    Minimum distance between ions (nm)";

pub const KNOWN: &[&str] = &[
    "topo",
    "pbc",
    "pos",
    "positive",
    "negative",
    "potential",
    "mindist",
];

pub const PROGRAM: Program = Program {
    name: "ion",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

#[derive(Parser)]
#[command(name = "ion", version, about = "Replace solvent molecules by ions")]
struct Args {
    /// Molecular topology file
    #[arg(long)]
    topo: String,

    /// Periodic boundary conditions (r = rectangular)
    #[arg(long)]
    pbc: String,

    /// Input coordinate file (solvated system)
    #[arg(long)]
    pos: String,

    /// Place positive ions: <number> <ion_name>
    #[arg(long, num_args = 2)]
    positive: Option<Vec<String>>,

    /// Place negative ions: <number> <ion_name>
    #[arg(long, num_args = 2)]
    negative: Option<Vec<String>>,

    /// Cutoff for electrostatic potential calculation (nm)
    #[arg(long)]
    potential: Option<f64>,

    /// Minimum distance between ions (nm)
    #[arg(long, default_value = "0.0")]
    mindist: f64,
}

/// Minimum image distance squared under rectangular PBC.
fn dist_sq_pbc(a: (f64, f64, f64), b: (f64, f64, f64), box_dim: &[f64; 3]) -> f64 {
    let mut dx = a.0 - b.0;
    let mut dy = a.1 - b.1;
    let mut dz = a.2 - b.2;
    dx -= (dx / box_dim[0]).round() * box_dim[0];
    dy -= (dy / box_dim[1]).round() * box_dim[1];
    dz -= (dz / box_dim[2]).round() * box_dim[2];
    dx * dx + dy * dy + dz * dz
}

/// Which ions to place and how.
#[derive(Debug, Clone, Default)]
pub struct IonOptions {
    /// `(count, residue/atom name)` of the positive ions, placed at the lowest potential.
    pub positive: Option<(usize, String)>,
    /// `(count, residue/atom name)` of the negative ions, placed at the highest potential.
    pub negative: Option<(usize, String)>,
    /// Cutoff for the electrostatic potential (nm). Required.
    pub potential: Option<f64>,
    /// Minimum distance between ions (nm).
    pub mindist: f64,
}

/// A system with some solvent molecules replaced by ions.
#[derive(Debug, Clone)]
pub struct IonResult {
    /// Solute, then the ions, then the solvent that was not replaced.
    pub atoms: Vec<G96Atom>,
    pub box_dims: newen::math::Vec3,
    pub n_positive: usize,
    pub n_negative: usize,
    /// Solvent molecules left after the replacement.
    pub n_solvent_molecules: usize,
    /// The potential (kJ/(mol·e)) at each replaced position, in placement order.
    pub potentials: Vec<f64>,
}

/// Replace solvent molecules by ions, gromos++ `ion`'s greedy rule.
///
/// Positive ions go where the electrostatic potential is lowest, negative ions where it is
/// highest; `mindist` keeps later ions away from the ones already placed. The number of atoms per
/// solvent molecule comes from the topology's SOLVENTATOM block.
pub fn add_ions(
    topo: &ParsedTopology,
    system: &LabeledCoordinateData,
    opts: &IonOptions,
) -> Result<IonResult, String> {
    let (n_positive, pos_name) = opts.positive.clone().unwrap_or((0, String::new()));
    let (n_negative, neg_name) = opts.negative.clone().unwrap_or((0, String::new()));
    if n_positive == 0 && n_negative == 0 {
        return Err("Error: specify @positive and/or @negative".to_string());
    }
    let cutoff = opts.potential.ok_or("Error: @potential cutoff required")?;

    let n_solute = topo.n_atoms;
    let atoms = &system.atoms;
    let box_dim = [system.box_dims.x, system.box_dims.y, system.box_dims.z];
    if atoms.len() < n_solute {
        return Err(format!(
            "Error: topology has {n_solute} atoms, the coordinate file {}",
            atoms.len()
        ));
    }
    let n_solvent_atoms = atoms.len() - n_solute;

    // Solvent molecule size from the topology (gromos++ reads SOLVENTATOM; 3 for SPC water).
    let atoms_per_solvent = topo.solvent_atoms.len().max(1);
    if !n_solvent_atoms.is_multiple_of(atoms_per_solvent) {
        eprintln!(
            "Warning: solvent atoms ({}) not divisible by {}",
            n_solvent_atoms, atoms_per_solvent
        );
    }
    let n_solvent_mols = n_solvent_atoms / atoms_per_solvent;

    if n_solvent_mols < n_positive + n_negative {
        return Err(format!(
            "Error: not enough solvent molecules ({}) for {} ions",
            n_solvent_mols,
            n_positive + n_negative
        ));
    }

    let parsed_topo = topo;
    let mindist = opts.mindist;

    // Calculate electrostatic potential at each solvent molecule's first atom
    // from all solute charges within the cutoff
    let cutoff_sq = cutoff * cutoff;
    let mut potentials: Vec<f64> = Vec::with_capacity(n_solvent_mols);

    // Coulomb constant: 1/(4*pi*eps0) in GROMOS units (kJ/mol * nm / e^2)
    let coulomb_k = newen_core::units::four_pi_eps_i;

    for mol_i in 0..n_solvent_mols {
        let ow_idx = n_solute + mol_i * atoms_per_solvent; // first atom (oxygen)
        let ow_pos = (
            atoms[ow_idx].pos.x,
            atoms[ow_idx].pos.y,
            atoms[ow_idx].pos.z,
        );

        // gromos++ `ion` (utils::Energy, ε = 1 so no reaction field): the electrostatic
        // energy of the molecule's first atom with every atom of every *other* molecule
        // within the cutoff, each pair shifted by −1/R_c. Solvent charges come from the
        // topology's SOLVENTATOM block.
        let mut phi = 0.0;
        let charge_of = |k: usize| -> f64 {
            if k < n_solute {
                parsed_topo.charges[k]
            } else {
                parsed_topo
                    .solvent_atoms
                    .get((k - n_solute) % atoms_per_solvent)
                    .map_or(0.0, |s| s.charge)
            }
        };
        let own = ow_idx..ow_idx + atoms_per_solvent;
        for (k, atom) in atoms.iter().enumerate() {
            if own.contains(&k) {
                continue;
            }
            let q_k = charge_of(k);
            if q_k.abs() < 1e-10 {
                continue;
            }
            let atom_pos = (atom.pos.x, atom.pos.y, atom.pos.z);
            let d2 = dist_sq_pbc(ow_pos, atom_pos, &box_dim);
            if d2 < cutoff_sq && d2 > 1e-10 {
                phi += coulomb_k * q_k * (1.0 / d2.sqrt() - 1.0 / cutoff);
            }
        }
        potentials.push(phi);
    }

    // Select solvent molecules to replace
    let mut selected: Vec<usize> = Vec::new(); // indices into solvent molecule list
    let mut excluded: HashSet<usize> = HashSet::new();
    let mindist_sq = mindist * mindist;

    // For positive ions: pick molecules with LOWEST potential (most negative)
    for _ in 0..n_positive {
        let mut best_idx = None;
        let mut best_pot = f64::MAX;

        #[allow(clippy::needless_range_loop)] // one index over several arrays
        for j in 0..n_solvent_mols {
            if excluded.contains(&j) {
                continue;
            }
            if potentials[j] < best_pot {
                best_pot = potentials[j];
                best_idx = Some(j);
            }
        }

        let idx =
            best_idx.ok_or("Error: not enough suitable solvent positions for positive ions")?;

        selected.push(idx);
        let ow_i = n_solute + idx * atoms_per_solvent;
        let pos_i = (atoms[ow_i].pos.x, atoms[ow_i].pos.y, atoms[ow_i].pos.z);

        // Exclude nearby solvent molecules
        for j in 0..n_solvent_mols {
            let ow_j = n_solute + j * atoms_per_solvent;
            let pos_j = (atoms[ow_j].pos.x, atoms[ow_j].pos.y, atoms[ow_j].pos.z);
            if dist_sq_pbc(pos_i, pos_j, &box_dim) < mindist_sq {
                excluded.insert(j);
            }
        }
    }

    // For negative ions: pick molecules with HIGHEST potential (most positive)
    for _ in 0..n_negative {
        let mut best_idx = None;
        let mut best_pot = f64::MIN;

        #[allow(clippy::needless_range_loop)] // one index over several arrays
        for j in 0..n_solvent_mols {
            if excluded.contains(&j) {
                continue;
            }
            if potentials[j] > best_pot {
                best_pot = potentials[j];
                best_idx = Some(j);
            }
        }

        let idx =
            best_idx.ok_or("Error: not enough suitable solvent positions for negative ions")?;

        selected.push(idx);
        let ow_i = n_solute + idx * atoms_per_solvent;
        let pos_i = (atoms[ow_i].pos.x, atoms[ow_i].pos.y, atoms[ow_i].pos.z);

        for j in 0..n_solvent_mols {
            let ow_j = n_solute + j * atoms_per_solvent;
            let pos_j = (atoms[ow_j].pos.x, atoms[ow_j].pos.y, atoms[ow_j].pos.z);
            if dist_sq_pbc(pos_i, pos_j, &box_dim) < mindist_sq {
                excluded.insert(j);
            }
        }
    }

    // Report selected positions
    let replaced_set: HashSet<usize> = selected.iter().copied().collect();
    let mut chosen_potentials = Vec::with_capacity(selected.len());
    for (i, &mol_idx) in selected.iter().enumerate() {
        let ow = n_solute + mol_idx * atoms_per_solvent;
        let ion_type = if i < n_positive { &pos_name } else { &neg_name };
        chosen_potentials.push(potentials[mol_idx]);
        eprintln!(
            "Ion {}: {} at solvent mol {} ({:.3}, {:.3}, {:.3}), potential = {:.3} kJ/(mol*e)",
            i + 1,
            ion_type,
            mol_idx + 1,
            atoms[ow].pos.x,
            atoms[ow].pos.y,
            atoms[ow].pos.z,
            potentials[mol_idx]
        );
    }

    // Solute as read, then the ions, then the solvent that was not replaced.
    let mut result_atoms: Vec<G96Atom> = Vec::with_capacity(atoms.len());
    let mut atom_counter = 1;
    for a in atoms.iter().take(n_solute) {
        result_atoms.push(G96Atom {
            atom_num: atom_counter,
            ..a.clone()
        });
        atom_counter += 1;
    }

    // Ions sit at the first atom of the water they replace; their residue numbers continue the
    // solute's. Positive ions first, then negative.
    let first_ion_res = if n_solute > 0 {
        atoms[n_solute - 1].res_num
    } else {
        0
    } + 1;
    for (i, &mol_idx) in selected.iter().enumerate() {
        let ow = n_solute + mol_idx * atoms_per_solvent;
        let ion_name = if i < n_positive { &pos_name } else { &neg_name };
        result_atoms.push(G96Atom {
            res_num: first_ion_res + i,
            res_name: ion_name.clone(),
            atom_name: ion_name.clone(),
            atom_num: atom_counter,
            pos: atoms[ow].pos,
        });
        atom_counter += 1;
    }

    for mol_i in 0..n_solvent_mols {
        if replaced_set.contains(&mol_i) {
            continue;
        }
        for k in 0..atoms_per_solvent {
            let a = &atoms[n_solute + mol_i * atoms_per_solvent + k];
            result_atoms.push(G96Atom {
                atom_num: atom_counter,
                ..a.clone()
            });
            atom_counter += 1;
        }
    }

    Ok(IonResult {
        atoms: result_atoms,
        box_dims: system.box_dims,
        n_positive,
        n_negative,
        n_solvent_molecules: n_solvent_mols - n_positive - n_negative,
        potentials: chosen_potentials,
    })
}

/// ion's own GROMOS96 layout (a `BOX` block and its own POSITION column widths).
pub fn format(result: &IonResult, pos_file: &str, opts: &IonOptions) -> String {
    let pos_name = opts
        .positive
        .as_ref()
        .map(|(_, n)| n.as_str())
        .unwrap_or("");
    let neg_name = opts
        .negative
        .as_ref()
        .map(|(_, n)| n.as_str())
        .unwrap_or("");

    let mut s = String::new();
    s.push_str("TITLE\n");
    s.push_str(&format!(
        "    ion: replaced {} solvent molecules in {}\n",
        result.n_positive + result.n_negative,
        pos_file
    ));
    if result.n_positive > 0 {
        s.push_str(&format!(
            "    {} positive ions ({})\n",
            result.n_positive, pos_name
        ));
    }
    if result.n_negative > 0 {
        s.push_str(&format!(
            "    {} negative ions ({})\n",
            result.n_negative, neg_name
        ));
    }
    s.push_str("END\nPOSITION\n");
    for a in &result.atoms {
        s.push_str(&format!(
            "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}\n",
            a.res_num, a.res_name, a.atom_name, a.atom_num, a.pos.x, a.pos.y, a.pos.z
        ));
    }
    s.push_str("END\nBOX\n");
    s.push_str(&format!(
        "{:15.9}{:15.9}{:15.9}\n",
        result.box_dims.x, result.box_dims.y, result.box_dims.z
    ));
    s.push_str("END\n");
    s
}

/// `<count> <name>` as clap collected it.
fn ion_spec(v: &Option<Vec<String>>, what: &str) -> Result<Option<(usize, String)>, String> {
    match v {
        None => Ok(None),
        Some(v) => {
            let n: usize = v[0]
                .parse()
                .map_err(|_| format!("Error: invalid {what} ion count '{}'", v[0]))?;
            Ok(Some((n, v[1].clone())))
        },
    }
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let Some(args) = parse_clap::<Args>(argv, out)? else {
        return Ok(());
    };

    let opts = IonOptions {
        positive: ion_spec(&args.positive, "positive")?,
        negative: ion_spec(&args.negative, "negative")?,
        potential: args.potential,
        mindist: args.mindist,
    };

    let parsed_topo = read_topology_file(&args.topo)
        .map_err(|e| format!("Error reading topology '{}': {e:?}", args.topo))?;
    eprintln!("Topology: {} solute atoms", parsed_topo.n_atoms);

    let data =
        read_g96_labeled(&args.pos).map_err(|e| format!("Error reading '{}': {e}", args.pos))?;
    eprintln!(
        "Coordinates: {} total atoms ({} solute, {} solvent)",
        data.atoms.len(),
        parsed_topo.n_atoms,
        data.atoms.len() - parsed_topo.n_atoms.min(data.atoms.len())
    );
    eprintln!(
        "Box: ({:.3}, {:.3}, {:.3}) nm",
        data.box_dims.x, data.box_dims.y, data.box_dims.z
    );

    let result = add_ions(&parsed_topo, &data, &opts)?;

    write!(out, "{}", format(&result, &args.pos, &opts)).map_err(io_err)?;

    eprintln!(
        "Output: {} solute + {} ions + {} solvent molecules ({} atoms)",
        parsed_topo.n_atoms,
        result.n_positive + result.n_negative,
        result.n_solvent_molecules,
        result.atoms.len()
    );
    Ok(())
}
