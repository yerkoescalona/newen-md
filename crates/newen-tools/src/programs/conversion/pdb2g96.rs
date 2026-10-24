//! pdb2g96 - Convert PDB files to GROMOS96 format (topology-guided)
//!
//! Usage: pdb2g96 @topo <top> @pdb <pdb> @lib <lib> [@gch]
//!    or: pdb2g96 @f <argfile>
//!
//! Reads a topology and a PDB file, then writes GROMOS96 coordinates
//! with atoms in topology order. Uses a library file for name mapping.
//! The @gch flag requests hydrogen coordinate generation (missing H atoms
//! get (0,0,0) and a warning; future: geometry-based generation).

use crate::programs::{io_err, parse_clap};
use clap::Parser;
use newen_io::args::{Program, Status};
use newen_io::pdb::PDBStructure;
use newen_io::topology::{read_topology_file, ParsedTopology};
use newen_io::G96Atom;
use std::collections::HashMap;
use std::io::Write;

pub const USAGE: &str = "pdb2g96 - Convert a PDB file to GROMOS96 coordinates\n\nUsage: pdb2g96 @topo <topology> @pdb <file.pdb> [@lib <residue library>] [@gch]\n\nArguments:\n  @topo  Molecular topology file (gives the atom order)\n  @pdb   Input PDB file\n  @lib   Residue/atom name library\n  @gch   Suppress warnings about missing hydrogens";

pub const KNOWN: &[&str] = &["topo", "pdb", "lib", "gch"];

pub const PROGRAM: Program = Program {
    name: "pdb2g96",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

/// Coordinates in the topology's atom order, plus what could not be matched.
#[derive(Debug, Clone)]
pub struct Pdb2g96Result {
    pub atoms: Vec<G96Atom>,
    /// Atoms with no PDB counterpart; their positions are (0, 0, 0).
    pub missing: usize,
    pub matched: usize,
}

#[derive(Parser)]
#[command(name = "pdb2g96", version, about)]
struct Args {
    /// Molecular topology file
    #[arg(long)]
    topo: String,

    /// Input PDB file
    #[arg(long)]
    pdb: String,

    /// Library for atom and residue name mapping
    #[arg(long)]
    lib: Option<String>,

    /// (Re)generate hydrogen coordinates
    #[arg(long, action = clap::ArgAction::SetTrue)]
    gch: bool,
}

/// Library mapping: PDB residue/atom names → topology names.
pub struct Library {
    residue_map: HashMap<String, String>, // pdb_res → topo_res
    atom_map: HashMap<String, HashMap<String, String>>, // topo_res → {pdb_atom → topo_atom}
}

impl Default for Library {
    fn default() -> Self {
        Self::new()
    }
}

impl Library {
    pub fn new() -> Self {
        Self {
            residue_map: HashMap::new(),
            atom_map: HashMap::new(),
        }
    }

    pub fn read(path: &str) -> Result<Self, String> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| format!("Cannot read library file '{}': {}", path, e))?;
        let mut lib = Self::new();
        let mut in_block = "";

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if trimmed == "END" {
                in_block = "";
                continue;
            }
            if trimmed == "TITLE" {
                in_block = "TITLE";
                continue;
            }
            if trimmed == "RESIDUENAMELIB" || trimmed == "RESIDUES" {
                in_block = "RESIDUES";
                continue;
            }
            if trimmed == "ATOMNAMELIB" || trimmed == "ATOMS" {
                in_block = "ATOMS";
                continue;
            }
            match in_block {
                "RESIDUES" => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 2 {
                        lib.residue_map
                            .insert(parts[0].to_string(), parts[1].to_string());
                    }
                },
                "ATOMS" => {
                    // Format: topo_residue pdb_atom topo_atom
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 3 {
                        lib.atom_map
                            .entry(parts[0].to_string())
                            .or_default()
                            .insert(parts[1].to_string(), parts[2].to_string());
                    }
                },
                _ => {},
            }
        }
        Ok(lib)
    }

    /// Map a PDB residue name to topology residue name.
    fn map_residue<'a>(&'a self, pdb_name: &'a str) -> &'a str {
        self.residue_map
            .get(pdb_name)
            .map(|s| s.as_str())
            .unwrap_or(pdb_name)
    }

    /// Map a PDB atom name to topology atom name for a given topology residue.
    fn map_atom<'a>(&'a self, topo_res: &str, pdb_atom: &'a str) -> &'a str {
        self.atom_map
            .get(topo_res)
            .and_then(|m| m.get(pdb_atom))
            .map(|s| s.as_str())
            .unwrap_or(pdb_atom)
    }
}

/// A PDB atom with its coordinates (already in nm from the PDB parser).
struct PdbAtomEntry {
    name: String,
    x: f64,
    y: f64,
    z: f64,
}

/// A PDB residue: name + list of atoms.
struct PdbResidue {
    name: String,
    atoms: Vec<PdbAtomEntry>,
}

/// Match topology atoms against a PDB residue's atoms, assigning coordinates.
/// Returns (matched, missing) counts.
fn match_atoms(
    topo_atoms: &[(usize, String)],
    topo_res_name: &str,
    pdb_res: &PdbResidue,
    used: &mut [bool],
    lib: &Library,
    coords: &mut [(f64, f64, f64)],
    masses: &[f64],
    gch: bool,
) -> (usize, usize) {
    let mut matched = 0;
    let mut missing = 0;
    for (global_idx, topo_atom_name) in topo_atoms {
        let mut found = false;
        for (j, pdb_atom) in pdb_res.atoms.iter().enumerate() {
            if used[j] {
                continue;
            }
            let mapped_name = lib.map_atom(topo_res_name, &pdb_atom.name);
            if mapped_name == topo_atom_name || pdb_atom.name == *topo_atom_name {
                coords[*global_idx] = (pdb_atom.x, pdb_atom.y, pdb_atom.z);
                used[j] = true;
                found = true;
                matched += 1;
                break;
            }
        }
        if !found {
            let is_hydrogen = masses[*global_idx] < 2.0;
            if !gch || !is_hydrogen {
                eprintln!(
                    "Warning: Could not find atom '{}' in residue '{}' — set to (0,0,0)",
                    topo_atom_name, topo_res_name
                );
            }
            missing += 1;
        }
    }
    (matched, missing)
}

/// Put the PDB's coordinates into the topology's atom order.
///
/// `lib` maps PDB residue/atom names onto the topology's; without one the names must match.
/// Atoms the PDB does not have keep (0, 0, 0) and are counted in
/// [`Pdb2g96Result::missing`].
pub fn pdb2g96(
    parsed_topo: &ParsedTopology,
    pdb: &PDBStructure,
    lib: &Library,
    gch: bool,
) -> Result<Pdb2g96Result, String> {
    // Group PDB atoms by residue
    let mut pdb_residues: Vec<PdbResidue> = Vec::new();
    for res in &pdb.residues {
        let mapped_name = lib.map_residue(res.name.trim());
        pdb_residues.push(PdbResidue {
            name: mapped_name.to_string(),
            atoms: res
                .atoms
                .iter()
                .map(|a| PdbAtomEntry {
                    name: a.name.trim().to_string(),
                    x: a.coord.x,
                    y: a.coord.y,
                    z: a.coord.z,
                })
                .collect(),
        });
    }

    // Build topology atom list: for each atom, (residue_number, residue_name, atom_name)
    // Group topology atoms by residue for matching
    let n_atoms = parsed_topo.n_atoms;
    let mut topo_residues: Vec<(usize, String, Vec<(usize, String)>)> = Vec::new();
    // (res_number, res_name, [(global_atom_idx, atom_name)])
    {
        let mut current_res_nr = 0;
        for i in 0..n_atoms {
            let res_nr = parsed_topo.residue_numbers[i];
            let atom_name = &parsed_topo.atom_names[i];
            if res_nr != current_res_nr {
                let res_name = if res_nr > 0 && res_nr <= parsed_topo.residue_names.len() {
                    parsed_topo.residue_names[res_nr - 1].clone()
                } else {
                    "UNK".to_string()
                };
                topo_residues.push((res_nr, res_name, Vec::new()));
                current_res_nr = res_nr;
            }
            if let Some(last) = topo_residues.last_mut() {
                last.2.push((i, atom_name.clone()));
            }
        }
    }

    // Match: for each topology residue, find matching PDB residue and assign coordinates
    // GROMOS topologies can have end-group residues (NH3+, COO-) that are not
    // separate residues in the PDB. These atoms belong to the adjacent PDB residue.
    let mut coords: Vec<(f64, f64, f64)> = vec![(0.0, 0.0, 0.0); n_atoms];
    let mut pdb_idx = 0; // current PDB residue index
    let mut missing_count = 0;
    let mut matched_count = 0;
    // Track which PDB atoms have been used (by pdb_idx)
    let mut used_map: HashMap<usize, Vec<bool>> = HashMap::new();

    for (_, topo_res_name, topo_atoms) in &topo_residues {
        if pdb_idx >= pdb_residues.len() {
            // No more PDB residues — try the last one for C-terminal end-groups
            if pdb_idx > 0 {
                let prev_idx = pdb_idx - 1;
                let pdb_res = &pdb_residues[prev_idx];
                let used = used_map
                    .entry(prev_idx)
                    .or_insert_with(|| vec![false; pdb_res.atoms.len()]);
                eprintln!("Note: Topology residue '{}' — looking for atoms in previous PDB residue '{}' (C-terminal end-group)",
                    topo_res_name, pdb_res.name);
                let (m, mi) = match_atoms(
                    topo_atoms,
                    topo_res_name,
                    pdb_res,
                    used,
                    lib,
                    &mut coords,
                    &parsed_topo.masses,
                    gch,
                );
                matched_count += m;
                missing_count += mi;
            } else {
                eprintln!(
                    "Warning: No more PDB residues for topology residue '{}'",
                    topo_res_name
                );
                missing_count += topo_atoms.len();
            }
            continue;
        }

        let pdb_res = &pdb_residues[pdb_idx];

        // Check if this topology residue matches the current PDB residue
        let names_match = pdb_res.name == *topo_res_name;

        // If names don't match, this might be an end-group residue (NH3+, COO-)
        // whose atoms are part of the current PDB residue. Don't advance pdb_idx.
        // If names match, we'll advance pdb_idx after matching.
        let advance_pdb = names_match;

        if !names_match {
            eprintln!(
                "Note: Topology residue '{}' — looking for atoms in PDB residue '{}' (end-group)",
                topo_res_name, pdb_res.name
            );
        }

        // Get or create the used-atoms tracker for this PDB residue
        let used = used_map
            .entry(pdb_idx)
            .or_insert_with(|| vec![false; pdb_res.atoms.len()]);

        for (global_idx, topo_atom_name) in topo_atoms {
            let mut found = false;

            for (j, pdb_atom) in pdb_res.atoms.iter().enumerate() {
                if used[j] {
                    continue;
                }

                let mapped_name = lib.map_atom(topo_res_name, &pdb_atom.name);
                if mapped_name == topo_atom_name || pdb_atom.name == *topo_atom_name {
                    coords[*global_idx] = (pdb_atom.x, pdb_atom.y, pdb_atom.z);
                    used[j] = true;
                    found = true;
                    matched_count += 1;
                    break;
                }
            }

            if !found {
                let is_hydrogen = parsed_topo.masses[*global_idx] < 2.0;
                if !gch || !is_hydrogen {
                    eprintln!(
                        "Warning: Could not find atom '{}' in residue '{}' — set to (0,0,0)",
                        topo_atom_name, topo_res_name
                    );
                }
                missing_count += 1;
            }
        }

        if advance_pdb {
            // Warn about unused PDB atoms in this residue
            let used = used_map.get(&pdb_idx).unwrap();
            for (j, pdb_atom) in pdb_res.atoms.iter().enumerate() {
                if !used[j] {
                    eprintln!(
                        "Warning: Ignored PDB atom '{}' in residue '{}'",
                        pdb_atom.name, pdb_res.name
                    );
                }
            }
            pdb_idx += 1;
        }
    }

    // The topology's atom order, with the coordinates found for each.
    let mut atoms = Vec::with_capacity(n_atoms);
    #[allow(clippy::needless_range_loop)] // one index over several arrays
    for i in 0..n_atoms {
        let res_nr = parsed_topo.residue_numbers[i];
        let res_name = if res_nr > 0 && res_nr <= parsed_topo.residue_names.len() {
            parsed_topo.residue_names[res_nr - 1].clone()
        } else {
            "UNK".to_string()
        };
        let (x, y, z) = coords[i];
        atoms.push(G96Atom {
            res_num: res_nr,
            res_name,
            atom_name: parsed_topo.atom_names[i].clone(),
            atom_num: i + 1,
            pos: newen::math::Vec3::new(x, y, z),
        });
    }

    Ok(Pdb2g96Result {
        atoms,
        missing: missing_count,
        matched: matched_count,
    })
}

/// pdb2g96's own GROMOS96 layout: TITLE and POSITION, no box (a PDB carries none).
pub fn format(result: &Pdb2g96Result, pdb_file: &str) -> String {
    let mut s = String::new();
    s.push_str("TITLE\n");
    s.push_str(&format!("    pdb2g96: Reordered atoms from {pdb_file}\n"));
    if result.missing > 0 {
        s.push_str(&format!(
            "    {} atoms could not be matched\n",
            result.missing
        ));
    }
    s.push_str("END\nPOSITION\n");
    for a in &result.atoms {
        s.push_str(&format!(
            "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}\n",
            a.res_num, a.res_name, a.atom_name, a.atom_num, a.pos.x, a.pos.y, a.pos.z
        ));
    }
    s.push_str("END\n");
    s
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let Some(args) = parse_clap::<Args>(argv, out)? else {
        return Ok(());
    };

    let parsed_topo = read_topology_file(&args.topo)
        .map_err(|e| format!("Error reading topology '{}': {e:?}", args.topo))?;
    eprintln!(
        "Topology: {} atoms, {} residues",
        parsed_topo.n_atoms,
        parsed_topo.residue_names.len()
    );

    let pdb = PDBStructure::read_pdb(&args.pdb)
        .map_err(|e| format!("Error reading PDB '{}': {e}", args.pdb))?;
    eprintln!("PDB: {} residues", pdb.residues.len());

    let lib = match args.lib {
        Some(ref lib_path) => Library::read(lib_path).map_err(|e| format!("Error: {e}"))?,
        None => Library::new(),
    };

    let result = pdb2g96(&parsed_topo, &pdb, &lib, args.gch)?;
    eprintln!(
        "Matched {} atoms, {} missing (set to 0,0,0)",
        result.matched, result.missing
    );

    write!(out, "{}", format(&result, &args.pdb)).map_err(io_err)?;
    Ok(())
}
