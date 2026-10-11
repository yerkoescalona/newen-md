//! com_top - Combine multiple GROMOS topology files
//!
//! Usage: com_top @topo <file1> [N:file2] ... @param <n> @solv <n>
//!    or: com_top @f <argfile>
//!
//! Combines multiple topology files into a single topology.
//! Use prefix 'N:' before a file to repeat that topology N times.
//! Parameters and solvent are taken from the topology indicated by @param/@solv.

use crate::programs::{io_err, parse_clap};
use clap::Parser;
use gromos_io::args::{Program, Status};
use gromos_io::topology::{read_topology_file, write_parsed_topology, ParsedTopology};
use std::collections::HashMap;
use std::io::Write;

pub const USAGE: &str = "com_top - Combine multiple GROMOS topology files\n\nUsage: com_top @topo <file1> [N:file2] ... @param <n> @solv <n>\n   or: com_top @f <argfile>\n\nArguments:\n  @topo   Topology files (prefix N: to repeat, e.g. 2:Na.top)\n  @param  Index (1-based) of the topology to take force-field parameters from\n  @solv   Index (1-based) of the topology to take the solvent from";

pub const KNOWN: &[&str] = &["topo", "param", "solv"];

pub const PROGRAM: Program = Program {
    name: "com_top",
    known: KNOWN,
    usage: USAGE,
    status: Status::Verified,
    run,
};

/// Combine multiple GROMOS topology files.
///
/// Supports GROMOS @-prefix args: com_top @topo GB3.top 2:Na.top @param 1 @solv 1
#[derive(Parser)]
#[command(name = "com_top", version, about)]
struct Args {
    /// Topology files (prefix N: to repeat, e.g. 2:Na.top)
    #[arg(long, num_args = 1..)]
    topo: Vec<String>,

    /// Index (1-based) of topology to take force field parameters from
    #[arg(long, default_value = "1")]
    param: usize,

    /// Index (1-based) of topology to take solvent from
    #[arg(long, default_value = "1")]
    solv: usize,
}

/// Parse a topology spec: "file.top" or "N:file.top"
fn parse_topo_spec(spec: &str) -> Result<(usize, String), String> {
    if let Some(pos) = spec.find(':') {
        let count: usize = spec[..pos]
            .parse()
            .map_err(|_| format!("Error: invalid repeat count in '{spec}'"))?;
        if count == 0 {
            return Err(format!("Error: repeat count 0 not allowed in '{spec}'"));
        }
        Ok((count, spec[pos + 1..].to_string()))
    } else {
        Ok((1, spec.to_string()))
    }
}

/// Combine `parts` into one topology, gromos++ `com_top`'s merge.
///
/// Each part is `(repeat, topology)`; the parts are laid out in order with atom, residue and
/// molecule numbering offset as they go. `param` and `solv` are 1-based indices *into the
/// repeat-expanded sequence* and say which topology contributes the force-field parameter blocks
/// and the solvent.
pub fn com_top(
    parts: &[(usize, &ParsedTopology)],
    param: usize,
    solv: usize,
) -> Result<ParsedTopology, String> {
    let total: usize = parts.iter().map(|(n, _)| n).sum();
    if total == 0 {
        return Err("Error: @topo requires at least one topology file".to_string());
    }
    if param < 1 || param > total {
        return Err(format!("Error: @param {param} out of range 1..{total}"));
    }
    if solv < 1 || solv > total {
        return Err(format!("Error: @solv {solv} out of range 1..{total}"));
    }

    // Which part provides the parameters and the solvent, in the expanded sequence.
    let mut count = 0;
    let mut param_topo = parts[0].1;
    let mut solv_topo = parts[0].1;
    for (repeat, topo) in parts {
        for _ in 0..*repeat {
            count += 1;
            if count == param {
                param_topo = topo;
            }
            if count == solv {
                solv_topo = topo;
            }
        }
    }

    Ok(merge(parts, param_topo, solv_topo))
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    let Some(args) = parse_clap::<Args>(argv, out)? else {
        return Ok(());
    };

    if args.topo.is_empty() {
        return Err("Error: @topo requires at least one topology file".to_string());
    }

    // Parse topology specs
    let specs: Vec<(usize, String)> = args
        .topo
        .iter()
        .map(|s| parse_topo_spec(s))
        .collect::<Result<_, _>>()?;

    // Read topology files (read each unique file once)
    let mut file_cache: HashMap<String, ParsedTopology> = HashMap::new();
    for (_, filename) in &specs {
        if !file_cache.contains_key(filename) {
            let parsed = read_topology_file(filename)
                .map_err(|e| format!("Error reading '{filename}': {e:?}"))?;
            eprintln!(
                "  Read '{}': {} atoms, {} residues",
                filename,
                parsed.n_atoms,
                parsed.residue_names.len()
            );
            file_cache.insert(filename.clone(), parsed);
        }
    }

    // Build title
    let mut title = String::from("COM_TOP: Combined topology using:\n");
    let mut c = 0;
    for (repeat, filename) in &specs {
        let start = c + 1;
        c += repeat;
        if *repeat > 1 {
            title.push_str(&format!("  {} .. {} : {}\n", start, c, filename));
        } else {
            title.push_str(&format!("  {} : {}\n", start, filename));
        }
    }

    let parts: Vec<(usize, &ParsedTopology)> =
        specs.iter().map(|(n, f)| (*n, &file_cache[f])).collect();
    let combined = com_top(&parts, args.param, args.solv)?;

    eprintln!(
        "Combined topology: {} atoms, {} residues, {} molecules",
        combined.n_atoms,
        combined.residue_names.len(),
        combined.solute_molecules.len()
    );

    write!(out, "{}", write_parsed_topology(&combined, &title)).map_err(io_err)?;
    Ok(())
}

fn merge(
    specs: &[(usize, &ParsedTopology)],
    param_topo: &ParsedTopology,
    solv_topo: &ParsedTopology,
) -> ParsedTopology {
    let mut combined = ParsedTopology {
        physical_constants: param_topo.physical_constants,
        n_atoms: 0,
        atom_names: Vec::new(),
        residue_numbers: Vec::new(),
        residue_names: Vec::new(),
        atom_type_names: param_topo.atom_type_names.clone(),
        masses: Vec::new(),
        charges: Vec::new(),
        iac: Vec::new(),
        chargegroup_codes: Vec::new(),
        exclusions: Vec::new(),
        one_four_pairs: Vec::new(),
        bonds: Vec::new(),
        bondh: Vec::new(),
        bond_parameters: param_topo.bond_parameters.clone(),
        angles: Vec::new(),
        angle_parameters: param_topo.angle_parameters.clone(),
        proper_dihedrals: Vec::new(),
        dihedral_parameters: param_topo.dihedral_parameters.clone(),
        improper_dihedrals: Vec::new(),
        improper_dihedral_parameters: param_topo.improper_dihedral_parameters.clone(),
        lj_parameters: param_topo.lj_parameters.clone(),
        temperature_groups: Vec::new(),
        pressure_groups: Vec::new(),
        solvent_atoms: solv_topo.solvent_atoms.clone(),
        solvent_constraints: solv_topo.solvent_constraints.clone(),
        solute_molecules: Vec::new(),
    };

    let mut atom_offset: usize = 0;
    let mut residue_offset: usize = 0;

    for (repeat, topo) in specs {
        for _ in 0..*repeat {
            // Merge residue names
            combined
                .residue_names
                .extend(topo.residue_names.iter().cloned());

            // Merge atoms
            for i in 0..topo.n_atoms {
                combined.atom_names.push(topo.atom_names[i].clone());
                combined
                    .residue_numbers
                    .push(topo.residue_numbers[i] + residue_offset);
                combined.masses.push(topo.masses[i]);
                combined.charges.push(topo.charges[i]);
                combined.iac.push(topo.iac[i]);
                combined.chargegroup_codes.push(topo.chargegroup_codes[i]);

                // Renumber exclusions
                let excl = topo.exclusions.get(i).cloned().unwrap_or_default();
                combined
                    .exclusions
                    .push(excl.iter().map(|&e| e + atom_offset).collect());

                // Renumber 1-4 pairs
                let pairs = topo.one_four_pairs.get(i).cloned().unwrap_or_default();
                combined
                    .one_four_pairs
                    .push(pairs.iter().map(|&p| p + atom_offset).collect());
            }

            // Merge bonds with renumbered atoms
            for &(i, j, t) in &topo.bonds {
                combined.bonds.push((i + atom_offset, j + atom_offset, t));
            }
            for &(i, j) in &topo.bondh {
                combined.bondh.push((i + atom_offset, j + atom_offset));
            }

            // Merge angles
            for &(i, j, k, t) in &topo.angles {
                combined
                    .angles
                    .push((i + atom_offset, j + atom_offset, k + atom_offset, t));
            }

            // Merge proper dihedrals
            for &(i, j, k, l, t) in &topo.proper_dihedrals {
                combined.proper_dihedrals.push((
                    i + atom_offset,
                    j + atom_offset,
                    k + atom_offset,
                    l + atom_offset,
                    t,
                ));
            }

            // Merge improper dihedrals
            for &(i, j, k, l, t) in &topo.improper_dihedrals {
                combined.improper_dihedrals.push((
                    i + atom_offset,
                    j + atom_offset,
                    k + atom_offset,
                    l + atom_offset,
                    t,
                ));
            }

            // Temperature / pressure groups: renumber
            for &tg in &topo.temperature_groups {
                combined.temperature_groups.push(tg + atom_offset);
            }
            for &pg in &topo.pressure_groups {
                combined.pressure_groups.push(pg + atom_offset);
            }

            // Solute molecules: renumber
            for &mol_last in &topo.solute_molecules {
                combined.solute_molecules.push(mol_last + atom_offset);
            }

            residue_offset += topo.residue_names.len();
            atom_offset += topo.n_atoms;
        }
    }

    combined.n_atoms = combined.masses.len();
    combined
}
