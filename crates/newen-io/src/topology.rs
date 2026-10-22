//! Parser for GROMOS topology files (.topo/.top)
//!
//! # IAC (Integer Atom Code) indexing convention
//!
//! **All IAC values stored in memory are 0-indexed** (0..N_types-1).
//! GROMOS file formats (SOLUTEATOM, CGPARAMETERS/LJPARAMETERS, SOLVENTATOM)
//! write IAC as 1-indexed.  The subtraction happens exactly once, at each
//! file-read site, so callers never need to worry about it:
//!
//! ```text
//!  file value → subtract 1 → topo.iac[i]  (0-indexed, matches LJ matrix rows/cols)
//!                          → PerturbedAtom.a_iac  (ptp.rs also subtracts 1)
//!                          → lj.get(iac_i, iac_j)
//! ```
//!
//! This matches GROMOS, which also stores IAC 0-indexed internally
//! (`in_topology.cc` line 1695: `topo.add_solute_atom(…, t-1, …)` and
//! lines 1147-1148: `--i; --j;` before indexing the LJ matrix).
//!
//! Format blocks:
//! - SOLUTEATOM: atom definitions (mass, charge, IAC)
//! - BOND/BONDSTRETCHTYPE: covalent bonds
//! - BONDANGLE/BONDANGLEBENDTYPE: bond angles
//! - DIHEDRAL/TORSDIHEDRALTYPE: torsional angles
//! - CGPARAMETERS: LJ interaction parameters
//! - TEMPERATUREGROUPS: temperature coupling groups

use crate::IoError;
use newen_core::topology::{
    Angle, AngleParameters, Bond, BondParameters, Dihedral, DihedralParameters,
    ImproperDihedralParameters, LJParameters, Topology,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::Path;

/// Parsed solvent atom template
#[derive(Debug, Clone)]
pub struct ParsedSolventAtom {
    pub iac: usize,
    pub name: String,
    pub mass: f64,
    pub charge: f64,
}

/// Parsed solvent constraint
#[derive(Debug, Clone)]
pub struct ParsedSolventConstraint {
    pub i: usize, // 0-indexed within solvent molecule
    pub j: usize,
    pub length: f64,
}

/// Parsed topology data from GROMOS .topo file
#[derive(Debug, Clone)]
pub struct ParsedTopology {
    /// Physical constants parsed from the PHYSICALCONSTANTS block.
    /// If the block is absent, falls back to the newen_core defaults.
    pub physical_constants: newen_core::units::PhysicalConstants,
    pub n_atoms: usize,
    pub atom_names: Vec<String>,      // PANM from SOLUTEATOM
    pub residue_numbers: Vec<usize>,  // MRES from SOLUTEATOM
    pub residue_names: Vec<String>,   // from RESNAME block
    pub atom_type_names: Vec<String>, // from ATOMTYPENAME block
    pub masses: Vec<f64>,
    pub charges: Vec<f64>,
    pub iac: Vec<usize>,
    pub chargegroup_codes: Vec<usize>, // CGC values: 1 = end of chargegroup
    pub exclusions: Vec<Vec<usize>>,
    pub one_four_pairs: Vec<Vec<usize>>, // 1-4 pairs from INE14 (0-based)
    pub bonds: Vec<(usize, usize, usize)>, // (i, j, type)
    /// The `(i, j)` of the bonds read from BONDH (bonds to hydrogen, what NTC = 2 constrains).
    pub bondh: Vec<(usize, usize)>,
    pub bond_parameters: Vec<BondParameters>,
    pub angles: Vec<(usize, usize, usize, usize)>, // (i, j, k, type)
    pub angle_parameters: Vec<AngleParameters>,
    pub proper_dihedrals: Vec<(usize, usize, usize, usize, usize)>, // (i, j, k, l, type)
    pub dihedral_parameters: Vec<DihedralParameters>,
    pub improper_dihedrals: Vec<(usize, usize, usize, usize, usize)>, // (i, j, k, l, type)
    pub improper_dihedral_parameters: Vec<ImproperDihedralParameters>,
    pub lj_parameters: HashMap<(usize, usize), LJParameters>,
    pub temperature_groups: Vec<usize>, // Last atom index of each group
    pub pressure_groups: Vec<usize>, // Last atom index of each pressure group (GROMOS: PRESSUREGROUPS)
    pub solvent_atoms: Vec<ParsedSolventAtom>,
    pub solvent_constraints: Vec<ParsedSolventConstraint>,
    pub solute_molecules: Vec<usize>, // last atom index of each molecule
}

/// Read GROMOS topology file
pub fn read_topology_file<P: AsRef<Path>>(path: P) -> Result<ParsedTopology, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    parse_topology(reader)
}

/// [`read_topology_file`] over an already-open reader — for callers holding the text, not a path
/// (a program that built a topology in memory and reads it back, the Python binding).
pub fn parse_topology(reader: impl BufRead) -> Result<ParsedTopology, IoError> {
    let mut lines = reader.lines();

    let mut n_atoms = 0;
    let mut atom_names = Vec::new();
    let mut residue_numbers = Vec::new();
    let mut residue_names = Vec::new();
    let mut atom_type_names = Vec::new();
    let mut masses = Vec::new();
    let mut charges = Vec::new();
    let mut iac = Vec::new();
    let mut chargegroup_codes = Vec::new();
    let mut exclusions = Vec::new();
    let mut one_four_pairs = Vec::new();
    let mut bonds = Vec::new();
    let mut bondh = Vec::new();
    let mut bond_parameters = Vec::new();
    let mut angles = Vec::new();
    let mut angle_parameters = Vec::new();
    let mut lj_parameters = HashMap::new();
    let mut temperature_groups = Vec::new();
    let mut pressure_groups = Vec::new();
    let mut proper_dihedrals = Vec::new();
    let mut dihedral_parameters = Vec::new();
    let mut improper_dihedrals = Vec::new();
    let mut improper_dihedral_parameters = Vec::new();
    let mut solvent_atoms = Vec::new();
    let mut solvent_constraints = Vec::new();
    let mut solute_molecules = Vec::new();
    let mut physical_constants = newen_core::units::PhysicalConstants {
        four_pi_eps_i: newen_core::units::four_pi_eps_i,
        kB: newen_core::units::kB,
        hBar: newen_core::units::hBar,
        spd_l: newen_core::units::spd_l,
    };

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match trimmed {
            "ATOMTYPENAME" => {
                parse_name_block(&mut lines, &mut atom_type_names)?;
            },
            "RESNAME" => {
                parse_name_block(&mut lines, &mut residue_names)?;
            },
            "SOLUTEATOM" => {
                parse_soluteatom(
                    &mut lines,
                    &mut n_atoms,
                    &mut atom_names,
                    &mut residue_numbers,
                    &mut masses,
                    &mut charges,
                    &mut iac,
                    &mut chargegroup_codes,
                    &mut exclusions,
                    &mut one_four_pairs,
                )?;
            },
            "BONDSTRETCHTYPE" => {
                parse_bond_types(&mut lines, &mut bond_parameters)?;
            },
            "BOND" => {
                parse_bonds(&mut lines, &mut bonds)?;
            },
            "BONDH" => {
                let start = bonds.len();
                parse_bonds(&mut lines, &mut bonds)?;
                bondh.extend(bonds[start..].iter().map(|&(i, j, _)| (i, j)));
            },
            "BONDANGLEBENDTYPE" => {
                parse_angle_types(&mut lines, &mut angle_parameters)?;
            },
            "BONDANGLE" | "BONDANGLEH" => {
                parse_angles(&mut lines, &mut angles)?;
            },
            "TORSDIHEDRALTYPE" | "TORSDIHEDRALTYPECODE" => {
                parse_dihedral_types(&mut lines, &mut dihedral_parameters)?;
            },
            "DIHEDRAL" | "DIHEDRALH" => {
                parse_dihedrals(&mut lines, &mut proper_dihedrals)?;
            },
            "IMPDIHEDRALTYPE" | "IMPDIHEDRALTYPECODE" => {
                parse_improper_dihedral_types(&mut lines, &mut improper_dihedral_parameters)?;
            },
            "IMPDIHEDRAL" | "IMPDIHEDRALH" => {
                parse_dihedrals(&mut lines, &mut improper_dihedrals)?;
            },
            "CGPARAMETERS" | "LJPARAMETERS" => {
                parse_lj_parameters(&mut lines, &mut lj_parameters)?;
            },
            "TEMPERATUREGROUPS" => {
                parse_temperature_groups(&mut lines, &mut temperature_groups)?;
            },
            "PRESSUREGROUPS" => {
                parse_pressure_groups(&mut lines, &mut pressure_groups)?;
            },
            "SOLVENTATOM" => {
                parse_solventatom(&mut lines, &mut solvent_atoms)?;
            },
            "SOLVENTCONSTR" => {
                parse_solventconstr(&mut lines, &mut solvent_constraints)?;
            },
            "SOLUTEMOLECULES" => {
                parse_solutemolecules(&mut lines, &mut solute_molecules)?;
            },
            "PHYSICALCONSTANTS" => {
                parse_physical_constants(&mut lines, &mut physical_constants)?;
            },
            _ => {
                // Skip unknown blocks
                continue;
            },
        }
    }

    Ok(ParsedTopology {
        physical_constants,
        n_atoms,
        atom_names,
        residue_numbers,
        residue_names,
        atom_type_names,
        masses,
        charges,
        iac,
        chargegroup_codes,
        exclusions,
        one_four_pairs,
        bonds,
        bondh,
        bond_parameters,
        angles,
        angle_parameters,
        proper_dihedrals,
        dihedral_parameters,
        improper_dihedrals,
        improper_dihedral_parameters,
        lj_parameters,
        temperature_groups,
        pressure_groups,
        solvent_atoms,
        solvent_constraints,
        solute_molecules,
    })
}

/// Parse ATOMTYPENAME or RESNAME block: count line, then one name per line.
fn parse_name_block<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    names: &mut Vec<String>,
) -> Result<(), IoError> {
    let mut count_read = false;
    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }
        if !count_read {
            // First non-comment line is the count — skip it
            count_read = true;
            continue;
        }
        names.push(trimmed.to_string());
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn parse_soluteatom<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    n_atoms: &mut usize,
    atom_names: &mut Vec<String>,
    residue_numbers: &mut Vec<usize>,
    masses: &mut Vec<f64>,
    charges: &mut Vec<f64>,
    iac: &mut Vec<usize>,
    chargegroup_codes: &mut Vec<usize>,
    exclusions: &mut Vec<Vec<usize>>,
    one_four_pairs: &mut Vec<Vec<usize>>,
) -> Result<(), IoError> {
    // First line after SOLUTEATOM is atom count
    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First non-comment line is atom count
        if *n_atoms == 0 {
            *n_atoms = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom count: {}", trimmed)))?;
            continue;
        }

        // Parse atom line (format: ATNM MRES PANM IAC MASS CG CGC INE ...)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 8 {
            let mres: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid MRES: {}", parts[1])))?;
            // IAC in GROMOS .topo is 1-indexed (1..=N_types).
            // Subtract 1 here so topo.iac is always 0-indexed (0..N_types) throughout the
            // codebase — matching array offsets, the LJ matrix rows/cols, and pttopo
            // PerturbedAtom.a_iac / b_iac (which the pttopo reader also converts to 0-indexed).
            let atom_iac: usize = parts[3]
                .parse::<usize>()
                .map_err(|_| IoError::ParseError(format!("Invalid IAC: {}", parts[3])))?
                .saturating_sub(1);
            let mass: f64 = parts[4]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid mass: {}", parts[4])))?;
            let charge: f64 = parts[5]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid charge: {}", parts[5])))?;
            let cgc: usize = parts[6]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CGC: {}", parts[6])))?;
            let n_exclusions: usize = parts[7].parse().map_err(|_| {
                IoError::ParseError(format!("Invalid exclusion count: {}", parts[7]))
            })?;

            atom_names.push(parts[2].to_string());
            residue_numbers.push(mres);
            iac.push(atom_iac);
            masses.push(mass);
            charges.push(charge);
            chargegroup_codes.push(cgc);

            // Exclusion indices are on the SAME line, starting at parts[8]
            // Format: ATNM MRES PANM IAC MASS CG CGC INE [excl1 excl2 ...]
            // gromosXX reads the rest of the atom as a token stream: INE exclusions, then
            // INE14 and its partners, wherever the line breaks fall (gromos++ wraps six per
            // line, other writers put the counts on their own lines).
            let mut tokens: Vec<String> = parts[8..].iter().map(|t| t.to_string()).collect();
            let mut need = n_exclusions + 1; // exclusions, then the INE14 count
            let mut n14: Option<usize> = None;
            loop {
                if let Some(count_pos) = (tokens.len() > n_exclusions).then_some(n_exclusions) {
                    if n14.is_none() {
                        n14 = Some(tokens[count_pos].parse::<usize>().map_err(|_| {
                            IoError::ParseError(format!(
                                "Invalid INE14 count: {}",
                                tokens[count_pos]
                            ))
                        })?);
                        need = n_exclusions + 1 + n14.unwrap();
                    }
                }
                if tokens.len() >= need && n14.is_some() {
                    break;
                }
                match lines.next() {
                    Some(Ok(next_line)) => {
                        let t = next_line.trim();
                        if t.is_empty() || t.starts_with('#') {
                            continue;
                        }
                        if t == "END" {
                            return Err(IoError::ParseError(format!(
                                "SOLUTEATOM: atom {} ends before its exclusion lists are complete",
                                atom_names.len()
                            )));
                        }
                        tokens.extend(t.split_whitespace().map(|x| x.to_string()));
                    },
                    _ => {
                        return Err(IoError::ParseError(
                            "SOLUTEATOM: unexpected end of file".to_string(),
                        ))
                    },
                }
            }
            let excl: Result<Vec<usize>, IoError> = tokens[..n_exclusions]
                .iter()
                .map(|t| {
                    t.parse::<usize>()
                        .map(|x| x - 1)
                        .map_err(|_| IoError::ParseError(format!("Invalid exclusion: {t}")))
                })
                .collect();
            exclusions.push(excl?);
            let n14 = n14.unwrap_or(0);
            let pairs: Result<Vec<usize>, IoError> = tokens
                [n_exclusions + 1..n_exclusions + 1 + n14]
                .iter()
                .map(|t| {
                    t.parse::<usize>()
                        .map(|x| x - 1)
                        .map_err(|_| IoError::ParseError(format!("Invalid 1-4 partner: {t}")))
                })
                .collect();
            one_four_pairs.push(pairs?);
        }
    }

    Ok(())
}

fn parse_bond_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    bond_parameters: &mut Vec<BondParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is type count
        if n_types == 0 {
            n_types = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid bond type count: {}", trimmed))
            })?;
            continue;
        }

        // Parse bond type (format: CB CHB B0)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let k_quartic: f64 = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CB: {}", parts[0])))?;
            let k_harmonic: f64 = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CHB: {}", parts[1])))?;
            let r0: f64 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid B0: {}", parts[2])))?;

            bond_parameters.push(BondParameters {
                k_quartic,
                k_harmonic,
                r0,
            });
        }
    }

    Ok(())
}

fn parse_bonds<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    bonds: &mut Vec<(usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_bonds = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is bond count
        if n_bonds == 0 {
            n_bonds = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond count: {}", trimmed)))?;
            continue;
        }

        // Parse bond (format: IB JB ICB)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let i: usize = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond atom i: {}", parts[0])))?;
            let j: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond atom j: {}", parts[1])))?;
            let bond_type: usize = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond type: {}", parts[2])))?;

            // Convert from 1-indexed to 0-indexed
            bonds.push((i - 1, j - 1, bond_type - 1));
        }
    }

    Ok(())
}

fn parse_angle_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    angle_parameters: &mut Vec<AngleParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is type count
        if n_types == 0 {
            n_types = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid angle type count: {}", trimmed))
            })?;
            continue;
        }

        // Parse angle type (format: CT CHT T0)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let k_cosine: f64 = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CT: {}", parts[0])))?;
            let k_harmonic: f64 = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CHT: {}", parts[1])))?;
            let theta0_deg: f64 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid T0: {}", parts[2])))?;

            angle_parameters.push(AngleParameters {
                k_cosine,
                k_harmonic: k_harmonic * 180.0 * 180.0
                    / (std::f64::consts::PI * std::f64::consts::PI), // Convert from kJ/(mol·deg²) to kJ/(mol·rad²)
                theta0: theta0_deg.to_radians(),
            });
        }
    }

    Ok(())
}

fn parse_angles<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    angles: &mut Vec<(usize, usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_angles = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is angle count
        if n_angles == 0 {
            n_angles = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle count: {}", trimmed)))?;
            continue;
        }

        // Parse angle (format: IT JT KT ICT)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 4 {
            let i: usize = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom i: {}", parts[0])))?;
            let j: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom j: {}", parts[1])))?;
            let k: usize = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom k: {}", parts[2])))?;
            let angle_type: usize = parts[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle type: {}", parts[3])))?;

            // Skip angles with duplicate atoms (invalid topology entries)
            if i == j || i == k || j == k {
                log::warn!(
                    "Skipping invalid angle with duplicate atoms: {}-{}-{} (type {})",
                    i,
                    j,
                    k,
                    angle_type
                );
                continue;
            }

            // Convert from 1-indexed to 0-indexed
            angles.push((i - 1, j - 1, k - 1, angle_type - 1));
        }
    }

    Ok(())
}

fn parse_dihedrals<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    dihedrals: &mut Vec<(usize, usize, usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_dihedrals = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }

        if n_dihedrals == 0 {
            n_dihedrals = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid dihedral count: {}", trimmed)))?;
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 5 {
            let i: usize = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[0])))?;
            let j: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[1])))?;
            let k: usize = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[2])))?;
            let l: usize = parts[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[3])))?;
            let dtype: usize = parts[4]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid type: {}", parts[4])))?;
            dihedrals.push((i - 1, j - 1, k - 1, l - 1, dtype - 1));
        }
    }
    Ok(())
}

fn parse_dihedral_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    params: &mut Vec<DihedralParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }

        if n_types == 0 {
            n_types = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid dihedral type count: {}", trimmed))
            })?;
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 3 {
            let cp: f64 = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CP: {}", parts[0])))?;
            let pd: f64 = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid PD: {}", parts[1])))?;
            let np: i32 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid NP: {}", parts[2])))?;
            params.push(DihedralParameters {
                k: cp,
                pd: pd * std::f64::consts::PI / 180.0,
                cospd: (pd * std::f64::consts::PI / 180.0).cos(),
                m: np,
            });
        }
    }
    Ok(())
}

fn parse_improper_dihedral_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    params: &mut Vec<ImproperDihedralParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }

        if n_types == 0 {
            n_types = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid improper type count: {}", trimmed))
            })?;
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 2 {
            let cq: f64 = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CQ: {}", parts[0])))?;
            let q0: f64 = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid Q0: {}", parts[1])))?;
            params.push(ImproperDihedralParameters {
                q0: q0 * std::f64::consts::PI / 180.0, // Convert to radians
                k: cq * 180.0 * 180.0 / (std::f64::consts::PI * std::f64::consts::PI), // Convert from kJ/(mol·deg²) to kJ/(mol·rad²)
            });
        }
    }
    Ok(())
}

fn parse_lj_parameters<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    lj_params: &mut HashMap<(usize, usize), LJParameters>,
) -> Result<(), IoError> {
    let mut n_pairs = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is pair count
        if n_pairs == 0 {
            n_pairs = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid LJ pair count: {}", trimmed)))?;
            continue;
        }

        // Parse LJ parameters (format: IAC JAC C12 C6 [CS12 CS6])
        // IAC/JAC are 1-indexed in the file; subtract 1 so the LJ matrix is 0-indexed,
        // consistent with topo.iac and pttopo PerturbedAtom.a_iac / b_iac.
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 4 {
            let iac: usize = parts[0]
                .parse::<usize>()
                .map_err(|_| IoError::ParseError(format!("Invalid IAC: {}", parts[0])))?
                .saturating_sub(1);
            let jac: usize = parts[1]
                .parse::<usize>()
                .map_err(|_| IoError::ParseError(format!("Invalid JAC: {}", parts[1])))?
                .saturating_sub(1);
            let c12: f64 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid C12: {}", parts[2])))?;
            let c6: f64 = parts[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid C6: {}", parts[3])))?;

            // Parse optional 1-4 parameters (CS12, CS6)
            let (cs12, cs6) = if parts.len() >= 6 {
                let cs12: f64 = parts[4]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid CS12: {}", parts[4])))?;
                let cs6: f64 = parts[5]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid CS6: {}", parts[5])))?;
                (cs12, cs6)
            } else {
                (c12, c6) // Default: same as normal
            };

            let params = LJParameters::new_with_14(c6, c12, cs6, cs12);
            lj_params.insert((iac, jac), params);

            // LJ matrix is symmetric
            if iac != jac {
                lj_params.insert((jac, iac), params);
            }
        }
    }

    Ok(())
}

fn parse_temperature_groups<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    temp_groups: &mut Vec<usize>,
) -> Result<(), IoError> {
    let mut n_groups = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is group count
        if n_groups == 0 {
            n_groups = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid temp group count: {}", trimmed))
            })?;
            continue;
        }

        // Next line(s) contain last atom indices for each group
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        for part in parts {
            if let Ok(last_atom) = part.parse::<usize>() {
                temp_groups.push(last_atom - 1); // Convert to 0-indexed
            }
        }
    }

    Ok(())
}

/// Parse PRESSUREGROUPS block — cumulative last-atom indices (1-indexed).
/// GROMOS format: first line = count, then last atom of each group.
fn parse_pressure_groups<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    pressure_groups: &mut Vec<usize>,
) -> Result<(), IoError> {
    let mut n_groups = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First data line is group count
        if n_groups == 0 {
            n_groups = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid pressure group count: {}", trimmed))
            })?;
            continue;
        }

        // Subsequent lines contain last atom indices (1-indexed, cumulative)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        for part in parts {
            if let Ok(last_atom) = part.parse::<usize>() {
                pressure_groups.push(last_atom); // Keep 1-indexed; convert in build_topology
            }
        }
    }

    Ok(())
}

fn parse_solventatom<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    solvent_atoms: &mut Vec<ParsedSolventAtom>,
) -> Result<(), IoError> {
    let mut nram: usize = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }

        if nram == 0 {
            nram = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid NRAM: {}", trimmed)))?;
            if nram == 0 {
                break;
            }
            continue;
        }

        // Format: I ANMS IACS MASS CGS
        // IAC is 1-indexed in file; subtract 1 for 0-indexed convention.
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 5 {
            let atom_iac: usize = parts[2]
                .parse::<usize>()
                .map_err(|_| IoError::ParseError(format!("Invalid solvent IAC: {}", parts[2])))?
                .saturating_sub(1);
            let mass: f64 = parts[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid solvent mass: {}", parts[3])))?;
            let charge: f64 = parts[4].parse().map_err(|_| {
                IoError::ParseError(format!("Invalid solvent charge: {}", parts[4]))
            })?;
            solvent_atoms.push(ParsedSolventAtom {
                iac: atom_iac,
                name: parts[1].to_string(),
                mass,
                charge,
            });
        }
    }
    Ok(())
}

fn parse_solventconstr<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    constraints: &mut Vec<ParsedSolventConstraint>,
) -> Result<(), IoError> {
    let mut ncons: usize = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }

        if ncons == 0 {
            ncons = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid NCONS: {}", trimmed)))?;
            if ncons == 0 {
                break;
            }
            continue;
        }

        // Format: ICONS JCONS CONS
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 3 {
            let i: usize = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid ICONS: {}", parts[0])))?;
            let j: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid JCONS: {}", parts[1])))?;
            let length: f64 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CONS: {}", parts[2])))?;
            constraints.push(ParsedSolventConstraint {
                i: i - 1, // 0-indexed within molecule
                j: j - 1,
                length,
            });
        }
    }
    Ok(())
}

/// Parse SOLUTEMOLECULES block: NSPM followed by last-atom indices for each molecule.
fn parse_solutemolecules<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    molecules: &mut Vec<usize>,
) -> Result<(), IoError> {
    let mut nspm: usize = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }

        if nspm == 0 {
            nspm = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid NSPM: {}", trimmed)))?;
            continue;
        }

        // The last atom index of each molecule (1-based); gromos++ writes several per line.
        for tok in trimmed.split_whitespace() {
            let last_atom: usize = tok
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid molecule end: {}", tok)))?;
            molecules.push(last_atom);
        }
    }
    Ok(())
}

/// Parse the PHYSICALCONSTANTS block.
///
/// Reads FPEPSI, HBAR, SPDL, BOLTZ in that order (gromosXX convention).
/// Any value that fails to parse is left at its default.
fn parse_physical_constants<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    constants: &mut newen_core::units::PhysicalConstants,
) -> Result<(), IoError> {
    let mut values: Vec<f64> = Vec::with_capacity(4);

    for line in lines.by_ref() {
        let line = line.map_err(|e| IoError::ParseError(e.to_string()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END" {
            break;
        }
        if let Ok(v) = trimmed.parse::<f64>() {
            values.push(v);
        }
    }

    // gromosXX order: FPEPSI, HBAR, SPDL, BOLTZ
    if let Some(&fpepsi) = values.first() {
        constants.four_pi_eps_i = fpepsi;
    }
    if let Some(&hbar) = values.get(1) {
        constants.hBar = hbar;
    }
    if let Some(&spdl) = values.get(2) {
        constants.spd_l = spdl;
    }
    if let Some(&boltz) = values.get(3) {
        constants.kB = boltz;
    }

    Ok(())
}

/// Convert ParsedTopology to Topology.
///
/// Builds solute topology and stores the solvent template.
/// To expand solvent molecules, call `topo.solvate(nsm)` afterwards
/// (GROMOS convention: NSM comes from the IMD SYSTEM block).
pub fn build_topology(parsed: ParsedTopology) -> Topology {
    use newen_core::topology::{Atom, ChargeGroup, SolventAtomTemplate, SolventConstraintTemplate};

    let mut topo = Topology::new();
    let n_solute = parsed.n_atoms;

    // --- Populate solute atoms ---
    topo.mass = parsed.masses.clone();
    topo.charge = parsed.charges.clone();
    topo.iac = parsed.iac.clone();

    for i in 0..n_solute {
        let atom_name = parsed
            .atom_names
            .get(i)
            .cloned()
            .unwrap_or_else(|| format!("ATOM{}", i + 1));
        let res_nr = parsed.residue_numbers.get(i).copied().unwrap_or(1);
        let res_name = if res_nr > 0 && res_nr <= parsed.residue_names.len() {
            parsed.residue_names[res_nr - 1].clone()
        } else {
            "UNK".to_string()
        };
        topo.moltypes[0].atoms.push(Atom {
            name: atom_name,
            residue_nr: res_nr,
            residue_name: res_name,
            iac: parsed.iac[i],
            mass: parsed.masses[i],
            charge: parsed.charges[i],
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: true,
        });
    }

    // --- Store solvent template (expanded later by topo.solvate(nsm)) ---
    for sa in &parsed.solvent_atoms {
        topo.solvent_atom_template.push(SolventAtomTemplate {
            iac: sa.iac,
            name: sa.name.clone(),
            mass: sa.mass,
            charge: sa.charge,
        });
    }
    for sc in &parsed.solvent_constraints {
        topo.solvent_constraint_template
            .push(SolventConstraintTemplate {
                i: sc.i,
                j: sc.j,
                length: sc.length,
            });
    }

    topo.compute_inverse_masses();

    // --- Store chargegroup codes and build solute chargegroups ---
    topo.chargegroup_codes = parsed.chargegroup_codes.clone();
    {
        let mut current_cg_atoms = Vec::new();
        for i in 0..n_solute {
            current_cg_atoms.push(i);
            let cgc = if i < parsed.chargegroup_codes.len() {
                parsed.chargegroup_codes[i]
            } else {
                1 // default: each atom is its own chargegroup
            };
            if cgc == 1 {
                topo.chargegroups.push(ChargeGroup {
                    atoms: current_cg_atoms.clone(),
                });
                current_cg_atoms.clear();
            }
        }
        if !current_cg_atoms.is_empty() {
            topo.chargegroups.push(ChargeGroup {
                atoms: current_cg_atoms,
            });
        }

        // Build atom_to_chargegroup mapping (solute only at this point)
        topo.atom_to_chargegroup = vec![0; n_solute];
        for (cg_idx, cg) in topo.chargegroups.iter().enumerate() {
            for &atom in &cg.atoms {
                if atom < n_solute {
                    topo.atom_to_chargegroup[atom] = cg_idx;
                }
            }
        }
        // Record how many CGs are solute (before solvate adds solvent CGs)
        topo.num_solute_chargegroups = topo.chargegroups.len();
    }

    // --- Populate solute molecules from SOLUTEMOLECULES block ---
    if !parsed.solute_molecules.is_empty() {
        let mut start = 0usize;
        for &last_atom in &parsed.solute_molecules {
            // last_atom is 1-indexed from the file, but no subtraction is needed:
            // a 1-indexed last atom N is the correct *exclusive* upper bound for a
            // 0-indexed Rust range (the -1 for 0-indexing and the +1 for exclusivity cancel).
            topo.molecules.push(start..last_atom);
            start = last_atom;
        }
    }

    // --- Initialize exclusions for solute atoms ---
    topo.exclusions = vec![Vec::new(); n_solute];

    for (i, excl_list) in parsed.exclusions.iter().enumerate() {
        for &j in excl_list {
            topo.exclusions[i].push(j);
            if j < n_solute {
                topo.exclusions[j].push(i);
            }
        }
    }

    // --- Store 1-4 pairs from INE14 and merge into exclusions ---
    // GROMOS convention: 1-4 pairs are excluded from the regular pairlist
    // (merged into all_exclusion) but computed separately with cs6/cs12.
    // They are NOT included in topo.exclusions (used for RF excluded interactions).
    // Stored symmetrically, like `exclusions`, so `Topology::is_excluded_or_14`
    // can query one side. INE14 lists each pair once (from the lower atom); the
    // list is only ever used for membership tests, never to enumerate pairs.
    topo.one_four_pairs = vec![Vec::new(); n_solute];
    for (i, pairs_14) in parsed.one_four_pairs.iter().enumerate() {
        for &j in pairs_14 {
            topo.one_four_pairs[i].push(j);
            if j < n_solute {
                topo.one_four_pairs[j].push(i);
            }
        }
    }

    // Build bonds
    let bondh: std::collections::HashSet<(usize, usize)> = parsed.bondh.iter().copied().collect();
    for (i, j, bond_type) in parsed.bonds {
        topo.moltypes[0].bonds.push(Bond {
            i,
            j,
            bond_type,
            hydrogen: bondh.contains(&(i, j)),
        });
        topo.exclusions[i].push(j);
        topo.exclusions[j].push(i);
    }

    topo.bond_parameters = parsed.bond_parameters;

    // Build angles
    for (i, j, k, angle_type) in parsed.angles {
        topo.moltypes[0].angles.push(Angle {
            i,
            j,
            k,
            angle_type,
        });
        topo.exclusions[i].push(k);
        topo.exclusions[k].push(i);
    }

    topo.angle_parameters = parsed.angle_parameters;

    // Build proper dihedrals
    for (i, j, k, l, dihedral_type) in parsed.proper_dihedrals {
        topo.moltypes[0].proper_dihedrals.push(Dihedral {
            i,
            j,
            k,
            l,
            dihedral_type,
        });
    }
    topo.dihedral_parameters = parsed.dihedral_parameters;

    // Build improper dihedrals
    for (i, j, k, l, dihedral_type) in parsed.improper_dihedrals {
        topo.moltypes[0].improper_dihedrals.push(Dihedral {
            i,
            j,
            k,
            l,
            dihedral_type,
        });
    }
    topo.improper_dihedral_parameters = parsed.improper_dihedral_parameters;

    // Build LJ parameter matrix.
    //
    // IMPORTANT: size the matrix from the LJ parameter *table* (which covers all
    // force-field atom types, e.g. 45 types in GROMOS 54A7), NOT from the atom
    // types actually present in this topology (which may be a small subset).
    //
    // Rationale: perturbed B-state IAC values (from pttopo) may refer to atom types
    // that don't appear in the regular topology but ARE listed in the LJ table.
    // Using topo.iac.max() would give a matrix too small to index those types,
    // causing out-of-bounds memory access via the unsafe get_unchecked path.
    let max_iac_table = parsed
        .lj_parameters
        .keys()
        .flat_map(|&(i, j)| [i, j])
        .max()
        .unwrap_or(0);
    // Also include atom IAC values (solute + solvent) to handle topologies without
    // a full LJPARAMETERS block.
    let max_iac = max_iac_table
        .max(topo.iac.iter().copied().max().unwrap_or(0))
        .max(
            topo.solvent_atom_template
                .iter()
                .map(|sa| sa.iac)
                .max()
                .unwrap_or(0),
        );
    let n_types = max_iac + 1;
    topo.lj_parameters = vec![vec![LJParameters::default(); n_types]; n_types];

    for ((iac, jac), params) in parsed.lj_parameters {
        if iac < n_types && jac < n_types {
            topo.lj_parameters[iac][jac] = params;
            topo.lj_parameters[jac][iac] = params;
            log::debug!(
                "LJ params[{}][{}]: c6={:.6e}, c12={:.6e}, cs6={:.6e}, cs12={:.6e}",
                iac,
                jac,
                params.c6,
                params.c12,
                params.cs6,
                params.cs12
            );
        }
    }

    log::debug!(
        "Built topology: {} solute atoms, {} chargegroups, LJ {}x{}, solvent template: {} atoms",
        n_solute,
        topo.chargegroups.len(),
        n_types,
        n_types,
        topo.solvent_atom_template.len()
    );

    // Populate the solute MoleculeType + MoleculeInstance now that solute.atoms
    // and the flat arrays are both fully built.
    topo.init_solute_moltype();

    // --- Populate solute pressure groups from PRESSUREGROUPS block ---
    // GROMOS: pressure_groups is a boundary vector [last_atom_1, last_atom_2, ...]
    // Each value is 1-indexed in the file, but no subtraction is needed: a 1-indexed
    // last atom N is the correct *exclusive* upper bound for a 0-indexed Rust range
    // (the -1 for 0-indexing and the +1 for exclusivity cancel).
    if !parsed.pressure_groups.is_empty() {
        let mut start = 0usize;
        for &last_atom in &parsed.pressure_groups {
            topo.pressure_groups.push(start..last_atom);
            start = last_atom;
        }
        log::debug!(
            "Parsed {} solute pressure groups from topology",
            topo.pressure_groups.len()
        );
    } else {
        // Default: all solute atoms as one pressure group (GROMOS fallback)
        if n_solute > 0 {
            topo.pressure_groups.push(0..n_solute);
            log::debug!("Default pressure group: all {} solute atoms", n_solute);
        }
    }

    // --- Solute temperature groups from TEMPERATUREGROUPS (0-based last atoms; that the last
    // one is the last solute atom is an MD requirement, checked where a run is prepared);
    // without the block the solute is one group (gromosXX in_topology.cc:271-276). `solvate`
    // adds one group per solvent molecule.
    {
        let mut start = 0usize;
        for &last in &parsed.temperature_groups {
            if last + 1 > start && last < n_solute {
                topo.temperature_groups.push((start..=last).collect());
                start = last + 1;
            }
        }
        if parsed.temperature_groups.is_empty() && n_solute > 0 {
            topo.temperature_groups.push((0..n_solute).collect());
        }
    }

    // Sorted and deduplicated: `Topology::pair_kind_of` relies on the order
    for list in topo.exclusions.iter_mut().chain(&mut topo.one_four_pairs) {
        list.sort_unstable();
        list.dedup();
    }

    topo
}

/// Write GROMOS topology file
pub fn write_topology_file<P: AsRef<Path>>(
    path: P,
    topo: &Topology,
    title: &str,
) -> Result<(), IoError> {
    use std::io::Write;

    let file = File::create(path.as_ref())
        .map_err(|e| IoError::WriteError(format!("Cannot create topology file: {}", e)))?;
    let mut writer = std::io::BufWriter::new(file);

    // TITLE block
    writeln!(writer, "TITLE").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "{}", title).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // PHYSICALCONSTANTS block (optional but good practice)
    writeln!(writer, "PHYSICALCONSTANTS").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(
        writer,
        "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)"
    )
    .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  {:.7}", newen_core::units::four_pi_eps_i) // gromosXX legacy value
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# HBAR: Planck's constant HBAR = H/(2* PI)")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  {}", newen_core::units::hBar)
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# SPDL: Speed of light (nm/ps)")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  {}", newen_core::units::spd_l)
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# BOLTZ: Boltzmann's constant kB")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  {:.8}", newen_core::units::kB) // gromosXX legacy value
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // TOPVERSION block
    writeln!(writer, "TOPVERSION").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "2.0").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // ATOMTYPENAME block
    writeln!(writer, "ATOMTYPENAME").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# NTYP: number of atom types")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    let n_types = topo.num_atom_types();
    writeln!(writer, "{}", n_types).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# TYPE: atom type names").map_err(|e| IoError::WriteError(e.to_string()))?;
    for i in 1..=n_types {
        writeln!(writer, "TYPE{}", i).map_err(|e| IoError::WriteError(e.to_string()))?;
    }
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // RESNAME block
    writeln!(writer, "RESNAME").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# NRES: number of residues")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "1").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# residue name").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "MOL").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // SOLUTEATOM block
    let n_atoms = topo.num_solute_atoms();
    writeln!(writer, "SOLUTEATOM").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "#   NRP: number of solute atoms")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "{}", n_atoms).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# ATNM: atom number").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# MRES: residue number").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# PANM: atom name").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# IAC: integer atom code").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# MASS: mass").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# CG: charge").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# CGC: charge group code").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# INE: number of exclusions")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "#ATNM MRES PANM IAC  MASS      CG     CGC INE")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    for (i, atom) in topo.moltypes[0].atoms.iter().enumerate() {
        let n_exclusions = topo.exclusions.get(i).map_or(0, |e| e.len());
        let cg_code = topo.atom_to_chargegroup.get(i).map_or(1, |&c| c + 1);

        writeln!(
            writer,
            "{:5} {:4} {:>4} {:3} {:8.4} {:8.5} {:4} {:3}",
            i + 1,
            atom.residue_nr,
            atom.name,
            atom.iac,
            atom.mass,
            atom.charge,
            cg_code,
            n_exclusions
        )
        .map_err(|e| IoError::WriteError(e.to_string()))?;

        // Write exclusions if any
        if let Some(exclusions) = topo.exclusions.get(i) {
            if !exclusions.is_empty() {
                let mut excl_vec: Vec<_> = exclusions.iter().collect();
                excl_vec.sort();
                for (j, excl) in excl_vec.iter().enumerate() {
                    write!(writer, "{:5}", *excl + 1)
                        .map_err(|e| IoError::WriteError(e.to_string()))?;
                    if (j + 1) % 10 == 0 {
                        writeln!(writer).map_err(|e| IoError::WriteError(e.to_string()))?;
                    }
                }
                writeln!(writer).map_err(|e| IoError::WriteError(e.to_string()))?;
            }
        }
    }
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // BONDSTRETCHTYPE block
    if !topo.bond_parameters.is_empty() {
        writeln!(writer, "BONDSTRETCHTYPE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NBTY: number of bond types")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.bond_parameters.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  CB     CHB       B0")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        for params in &topo.bond_parameters {
            writeln!(
                writer,
                "{:10.5e} {:10.5e} {:10.7}",
                params.k_quartic, params.k_harmonic, params.r0
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BOND block
    if !&topo.moltypes[0].bonds.is_empty() {
        writeln!(writer, "BOND").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NBH: number of bonds")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.moltypes[0].bonds.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  IB   JB  ICB").map_err(|e| IoError::WriteError(e.to_string()))?;
        for bond in &topo.moltypes[0].bonds {
            writeln!(
                writer,
                "{:5}{:5}{:5}",
                bond.i + 1,
                bond.j + 1,
                bond.bond_type + 1
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BONDANGLEBENDTYPE block
    if !topo.angle_parameters.is_empty() {
        writeln!(writer, "BONDANGLEBENDTYPE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NTTY: number of angle types")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.angle_parameters.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#   CT      CHT     T0[deg]")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        for params in &topo.angle_parameters {
            writeln!(
                writer,
                "{:10.5e} {:10.5e} {:10.4}",
                params.k_cosine,
                params.k_harmonic,
                params.theta0.to_degrees()
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BONDANGLE block
    if !&topo.moltypes[0].angles.is_empty() {
        writeln!(writer, "BONDANGLE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NTHEH: number of angles")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.moltypes[0].angles.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  IT   JT   KT  ICT").map_err(|e| IoError::WriteError(e.to_string()))?;
        for angle in &topo.moltypes[0].angles {
            writeln!(
                writer,
                "{:5}{:5}{:5}{:5}",
                angle.i + 1,
                angle.j + 1,
                angle.k + 1,
                angle.angle_type + 1
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    writer
        .flush()
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    Ok(())
}

/// Give `topo` the solvent molecules a coordinate set of `n_atoms` atoms implies (the analysis
/// programs read every atom of a frame, like gromos++ `InG96::select("ALL")`): whole solvent
/// molecules beyond the solute, or an error when the count does not fit.
pub fn solvate_to_atoms(topo: &mut newen_core::Topology, n_atoms: usize) -> Result<(), IoError> {
    let n_solute = topo.num_solute_atoms();
    if n_atoms < n_solute {
        return Err(IoError::FormatError(format!(
            "{n_atoms} atoms in the coordinates, {n_solute} solute atoms in the topology"
        )));
    }
    let per = topo.solvent_atom_template.len();
    let extra = n_atoms - n_solute;
    if extra == 0 {
        return Ok(());
    }
    if per == 0 || !extra.is_multiple_of(per) {
        return Err(IoError::FormatError(format!(
            "{extra} atoms beyond the solute do not form whole solvent molecules of {per} atoms"
        )));
    }
    let want = extra / per;
    if topo.num_solvent_molecules() != want {
        topo.solvate(want);
    }
    Ok(())
}

/// Write a parsed topology in gromos++ `OutTopology` layout — every block gromosXX reads,
/// 1-based indices in the file, exclusion and 1-4 lists six per line. The inverse of
/// [`read_topology_file`]: what this writes reads back identical.
pub fn write_parsed_topology(topo: &ParsedTopology, title: &str) -> String {
    use crate::table::{cpp_e, cpp_g};
    use std::f64::consts::PI;
    use std::fmt::Write as _;
    let mut o = String::new();
    // every tenth entry of a list gets a "# n" marker, as gromos++ writes them
    let marker = |o: &mut String, i: usize| {
        if i > 0 && i.is_multiple_of(10) {
            let _ = writeln!(o, "# {i}");
        }
    };
    let _ = write!(o, "TITLE\n{title}\nEND\n");
    let pc = &topo.physical_constants;
    let _ = write!(
        o,
        "PHYSICALCONSTANTS\n# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n{}\n# HBAR: Planck's constant HBAR = H/(2* PI)\n{}\n# SPDL: Speed of light (nm/ps)\n{}\n# BOLTZ: Boltzmann's constant kB\n{}\nEND\n",
        cpp_g(pc.four_pi_eps_i, 10),
        cpp_g(pc.hBar, 10),
        cpp_g(pc.spd_l, 10),
        cpp_g(pc.kB, 10)
    );
    o.push_str("TOPVERSION\n2.0\nEND\n");
    let _ = write!(
        o,
        "ATOMTYPENAME\n# NRATT: number of van der Waals atom types\n{}\n# TYPE: atom type names\n",
        topo.atom_type_names.len()
    );
    for (i, t) in topo.atom_type_names.iter().enumerate() {
        marker(&mut o, i);
        let _ = writeln!(o, "{t}");
    }
    o.push_str("END\n");
    let _ = write!(
        o,
        "RESNAME\n# NRAA2: number of residues in a solute molecule\n{}\n# AANM: residue names\n",
        topo.residue_names.len()
    );
    for (i, r) in topo.residue_names.iter().enumerate() {
        marker(&mut o, i);
        let _ = writeln!(o, "{r}");
    }
    o.push_str("END\n");
    let _ = write!(
        o,
        "SOLUTEATOM\n#   NRP: number of solute atoms\n{:5}\n",
        topo.n_atoms
    );
    o.push_str("#  ATNM: atom number\n#  MRES: residue number\n#  PANM: atom name of solute atom\n#   IAC: integer (van der Waals) atom type code\n#  MASS: mass of solute atom\n#    CG: charge of solute atom\n#   CGC: charge group code (0 or 1)\n#   INE: number of excluded atoms\n# INE14: number of 1-4 interactions\n# ATNM MRES PANM IAC     MASS       CG  CGC INE\n#                                           INE14\n");
    for i in 0..topo.n_atoms {
        let excl: Vec<usize> = topo
            .exclusions
            .get(i)
            .map(|e| e.iter().map(|x| x + 1).collect())
            .unwrap_or_default();
        let one_four: Vec<usize> = topo
            .one_four_pairs
            .get(i)
            .map(|e| e.iter().map(|x| x + 1).collect())
            .unwrap_or_default();
        let _ = write!(
            o,
            "{:6} {:>4} {:>4} {:>3} {:8.5} {:8.5} {:>2} {:>5}",
            i + 1,
            topo.residue_numbers[i],
            topo.atom_names[i],
            topo.iac[i] + 1,
            topo.masses[i],
            topo.charges[i],
            topo.chargegroup_codes[i],
            excl.len()
        );
        for (k, x) in excl.iter().enumerate() {
            if k.is_multiple_of(6) && k != 0 {
                let _ = write!(o, "\n{:47}", "");
            }
            let _ = write!(o, " {:>5}", x);
        }
        let _ = write!(o, "\n{:46}{}", "", one_four.len());
        for (k, x) in one_four.iter().enumerate() {
            if k.is_multiple_of(6) && k != 0 {
                let _ = write!(o, "\n{:45}", "");
            }
            let _ = write!(o, " {:>5}", x);
        }
        o.push('\n');
    }
    o.push_str("END\n");
    let _ = write!(o, "BONDSTRETCHTYPE\n#  NBTY: number of covalent bond types\n{}\n#  CB:  quartic force constant\n#  CHB: harmonic force constant\n#  B0:  bond length at minimum energy\n#         CB         CHB         B0\n", topo.bond_parameters.len());
    for (i, b) in topo.bond_parameters.iter().enumerate() {
        marker(&mut o, i);
        let _ = writeln!(
            o,
            "{:>16} {:>15} {:>15}",
            cpp_e(b.k_quartic, 5),
            cpp_e(b.k_harmonic, 5),
            cpp_e(b.r0, 5)
        );
    }
    o.push_str("END\n");
    let is_h = |i: usize| (topo.masses[i] - 1.008).abs() < 1e-9;
    let (bondh, bond): (Vec<&(usize, usize, usize)>, Vec<&(usize, usize, usize)>) =
        topo.bonds.iter().partition(|b| is_h(b.0) || is_h(b.1));
    for (header, v) in [
        ("BONDH\n#  NBONH: number of bonds involving H atoms in solute\n{n}\n#  IBH, JBH: atom sequence numbers of atoms forming a bond\n#  ICBH: bond type code\n#   IBH    JBH ICBH\n", &bondh),
        ("BOND\n#  NBON: number of bonds NOT involving H atoms in solute\n{n}\n#  IB, JB: atom sequence numbers of atoms forming a bond\n#  ICB: bond type code\n#    IB     JB  ICB\n", &bond),
    ] {
        o.push_str(&header.replace("{n}", &v.len().to_string()));
        for (k, b) in v.iter().enumerate() {
            marker(&mut o, k);
            let _ = writeln!(o, "{:>7} {:>6} {:>4}", b.0 + 1, b.1 + 1, b.2 + 1);
        }
        o.push_str("END\n");
    }
    let _ = write!(o, "BONDANGLEBENDTYPE\n#  NTTY: number of bond angle types\n{}\n#  CT:  force constant (based on potential\n#       harmonic in the angle cosine)\n#  CHT: force constant (based on potential\n#       harmonic in the angle)\n#  T0:  bond angle at minimum energy in degrees\n#         CT         CHT          T0\n", topo.angle_parameters.len());
    for (i, a) in topo.angle_parameters.iter().enumerate() {
        marker(&mut o, i);
        // CHT is stored per rad², the file carries it per degree²
        let _ = writeln!(
            o,
            "{:>16} {:>15} {:>15}",
            cpp_e(a.k_cosine, 5),
            cpp_e(a.k_harmonic * PI * PI / (180.0 * 180.0), 5),
            cpp_e(a.theta0.to_degrees(), 5)
        );
    }
    o.push_str("END\n");
    let (angleh, angle): (
        Vec<&(usize, usize, usize, usize)>,
        Vec<&(usize, usize, usize, usize)>,
    ) = topo
        .angles
        .iter()
        .partition(|a| is_h(a.0) || is_h(a.1) || is_h(a.2));
    for (header, v) in [
        ("BONDANGLEH\n#  NTHEH: number of bond angles involving H atoms in solute\n{n}\n#  ITH, JTH, KTH: atom sequence numbers\n#    of atoms forming a bond angle in solute\n#  ICTH: bond angle type code\n#   ITH    JTH    KTH ICTH\n", &angleh),
        ("BONDANGLE\n#  NTHE: number of bond angles NOT\n#        involving H atoms in solute\n{n}\n#  IT, JT, KT: atom sequence numbers of atoms\n#     forming a bond angle\n#  ICT: bond angle type code\n#    IT     JT     KT  ICT\n", &angle),
    ] {
        o.push_str(&header.replace("{n}", &v.len().to_string()));
        for (k, a) in v.iter().enumerate() {
            marker(&mut o, k);
            let _ = writeln!(
                o,
                "{:>7} {:>6} {:>6} {:>4}",
                a.0 + 1,
                a.1 + 1,
                a.2 + 1,
                a.3 + 1
            );
        }
        o.push_str("END\n");
    }
    let _ = write!(o, "IMPDIHEDRALTYPE\n#  NQTY: number of improper dihedrals\n{}\n#  CQ: force constant of improper dihedral per degrees square\n#  Q0: improper dihedral angle at minimum energy in degrees\n#            CQ             Q0\n", topo.improper_dihedral_parameters.len());
    for (i, p) in topo.improper_dihedral_parameters.iter().enumerate() {
        marker(&mut o, i);
        // CQ is stored per rad², the file carries it per degree²
        let _ = writeln!(
            o,
            "{:>15} {:>14}",
            cpp_e(p.k * PI * PI / (180.0 * 180.0), 5),
            cpp_e(p.q0.to_degrees(), 5)
        );
    }
    o.push_str("END\n");
    type Dih<'a> = Vec<&'a (usize, usize, usize, usize, usize)>;
    let (imph, imp): (Dih, Dih) = topo
        .improper_dihedrals
        .iter()
        .partition(|d| is_h(d.0) || is_h(d.1) || is_h(d.2) || is_h(d.3));
    let (dihh, dih): (Dih, Dih) = topo
        .proper_dihedrals
        .iter()
        .partition(|d| is_h(d.0) || is_h(d.1) || is_h(d.2) || is_h(d.3));
    let write_dihedrals = |o: &mut String, header: &str, v: &Dih| {
        o.push_str(&header.replace("{n}", &v.len().to_string()));
        for (k, d) in v.iter().enumerate() {
            marker(o, k);
            let _ = writeln!(
                o,
                "{:>7} {:>6} {:>6} {:>6} {:>4}",
                d.0 + 1,
                d.1 + 1,
                d.2 + 1,
                d.3 + 1,
                d.4 + 1
            );
        }
        o.push_str("END\n");
    };
    write_dihedrals(&mut o, "IMPDIHEDRALH\n#  NQHIH: number of improper dihedrals\n#         involving H atoms in the solute\n{n}\n#  IQH,JQH,KQH,LQH: atom sequence numbers\n#     of atoms forming an improper dihedral\n#  ICQH: improper dihedral type code\n#   IQH    JQH    KQH    LQH ICQH\n", &imph);
    write_dihedrals(&mut o, "IMPDIHEDRAL\n#  NQHI: number of improper dihedrals NOT\n#    involving H atoms in solute\n{n}\n#  IQ,JQ,KQ,LQ: atom sequence numbers of atoms\n#    forming an improper dihedral\n#  ICQ: improper dihedral type code\n#    IQ     JQ     KQ     LQ  ICQ\n", &imp);
    let _ = write!(o, "TORSDIHEDRALTYPE\n#  NPTY: number of dihedral types\n{}\n#  CP: force constant\n#  PD: phase-shift angle\n#  NP: multiplicity\n#       CP         PD  NP\n", topo.dihedral_parameters.len());
    for (i, d) in topo.dihedral_parameters.iter().enumerate() {
        marker(&mut o, i);
        let _ = writeln!(o, "{:>10.5} {:>10.5} {:>3}", d.k, d.pd.to_degrees(), d.m);
    }
    o.push_str("END\n");
    write_dihedrals(&mut o, "DIHEDRALH\n#  NPHIH: number of dihedrals involving H atoms in solute\n{n}\n#  IPH, JPH, KPH, LPH: atom sequence numbers\n#    of atoms forming a dihedral\n#  ICPH: dihedral type code\n#   IPH    JPH    KPH    LPH ICPH\n", &dihh);
    write_dihedrals(&mut o, "DIHEDRAL\n#  NPHI: number of dihedrals NOT involving H atoms in solute\n{n}\n#  IP, JP, KP, LP: atom sequence numbers\n#     of atoms forming a dihedral\n#  ICP: dihedral type code\n#    IP     JP     KP     LP  ICP\n", &dih);
    o.push_str("CROSSDIHEDRALH\n#  NPHIH: number of cross dihedrals involving H atoms in solute\n0\n#  APH, BPH, CPH, DPH, EPH, FPH, GPH, HPH: atom sequence numbers\n#    of atoms forming a dihedral\n#  ICCH: dihedral type code\n#   APH    BPH    CPH    DPH    EPH    FPH    GPH    HPH ICCH\nEND\n");
    o.push_str("CROSSDIHEDRAL\n#  NPPC: number of cross dihedrals NOT involving H atoms in solute\n0\n#  AP, BP, CP, DP, EP, FP, GP, HP: atom sequence numbers\n#     of atoms forming a dihedral\n#  ICC: dihedral type code\n#    AP     BP     CP     DP     EP     FP     GP     HP  ICC\nEND\n");
    let n_types = topo.atom_type_names.len();
    let _ = write!(o, "LJPARAMETERS\n#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2\n{}\n#  IAC,JAC: integer (van der Waals) atom type code\n#  C12: r**(-12) term in nonbonded interactions\n#   C6: r**(-6) term in nonbonded interactions\n# CS12: r**(-12) term in 1-4 nonbonded interactions\n#  CS6: r**(-6) term in 1-4 nonbonded interactions\n# IAC  JAC           C12            C6          CS12           CS6\n", n_types * (n_types + 1) / 2);
    for i in 0..n_types {
        for j in 0..=i {
            let lj = topo
                .lj_parameters
                .get(&(j, i))
                .or_else(|| topo.lj_parameters.get(&(i, j)));
            let (c6, c12, cs6, cs12) = lj
                .map(|p| (p.c6, p.c12, p.cs6, p.cs12))
                .unwrap_or((0.0, 0.0, 0.0, 0.0));
            let _ = writeln!(
                o,
                "{:>5}{:>5}{:>14}{:>14}{:>14}{:>14}",
                j + 1,
                i + 1,
                cpp_e(c12, 6),
                cpp_e(c6, 6),
                cpp_e(cs12, 6),
                cpp_e(cs6, 6)
            );
        }
        o.push_str("#\n");
    }
    o.push_str("END\n");
    let write_groups = |o: &mut String, header: &str, groups: &[usize]| {
        let _ = writeln!(o, "{header}{:>10}", groups.len());
        for (i, g) in groups.iter().enumerate() {
            let _ = write!(o, " {:>5}", g);
            if (i + 1).is_multiple_of(10) {
                o.push('\n');
            }
        }
        if !groups.len().is_multiple_of(10) {
            o.push('\n');
        }
        o.push_str("END\n");
    };
    write_groups(&mut o, "SOLUTEMOLECULES\n# NSPM: number of separate molecules in solute block\n# NSP[1...NSPM]: atom sequence number of last atom\n#                of the successive submolecules\n#      NSPM  NSP[1...NSPM]\n", &topo.solute_molecules);
    let temp: Vec<usize> = topo.temperature_groups.iter().map(|g| g + 1).collect();
    write_groups(&mut o, "TEMPERATUREGROUPS\n# NSTM: number of temperature atom groups\n# NST[1...NSTM]: atom sequence number of last atom\n#                of the successive temperature atom groups\n#      NSTM  NST[1...NSTM]\n", &temp);
    write_groups(&mut o, "PRESSUREGROUPS\n# NSVM: number of pressure atom groups\n# NSV[1...NSVM]: atom sequence number of last atom\n#                of the successive pressure atom groups\n#      NSVM  NSV[1...NSVM]\n", &topo.pressure_groups);
    o.push_str("LJEXCEPTIONS\n# This block defines special LJ-interactions based on atom numbers\n# This overrules the normal LJ-parameters (including 1-4 interactions)\n# NEX: number of exceptions\n0\n# AT1  AT2           C12            C6\nEND\n");
    let _ = write!(o, "SOLVENTATOM\n#  NRAM: number of atoms per solvent molecule\n{}\n#     I: solvent atom sequence number\n#  IACS: integer (van der Waals) atom type code\n#  ANMS: atom name of solvent atom\n#  MASS: mass of solvent atom\n#   CGS: charge of solvent atom\n#  I  ANMS IACS      MASS        CGS\n", topo.solvent_atoms.len());
    for (idx, a) in topo.solvent_atoms.iter().enumerate() {
        let _ = writeln!(
            o,
            "{:>4}  {:>4} {:>3} {:>10.5} {:>10.5}",
            idx + 1,
            a.name,
            a.iac + 1,
            a.mass,
            a.charge
        );
    }
    o.push_str("END\n");
    let _ = write!(o, "SOLVENTCONSTR\n#  NCONS: number of constraints\n{}\n#  ICONS, JCONS: atom sequence numbers forming constraint\n#   CONS constraint length\n#ICONS JCONS         CONS\n", topo.solvent_constraints.len());
    for c in &topo.solvent_constraints {
        let _ = writeln!(o, "{:>5} {:>4} {:>14.7}", c.i + 1, c.j + 1, c.length);
    }
    o.push_str("END\n# end of topology file\n");
    o
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Minimal 4-atom GROMOS topology (cg16-style).
    /// SOLUTEATOM format: ATNM MRES PANM IAC MASS CG CGC INE IEXCL...
    const CG16_TOPO: &str = "\
TITLE
  4-atom test topology (cg16)
END
SOLUTEATOM
# ATNM MRES PANM IAC   MASS      CG   CGC INE IEXCL
4
    1    1   C     6  62.000000  0.000000  0  2  2  3
                                              0
    2    1   H     1  12.000000  0.000000  0  1  1
                                              0
    3    1   H     1  12.000000  0.000000  0  1  1
                                              0
    4    1   H     1  12.000000  0.000000  0  0
                                              0
END
BONDSTRETCHTYPE
# CB        CHB       B0
1
  1000.000  500.000   0.150
END
BOND
# IB  JB  ICB
3
  1   2   1
  2   3   1
  3   4   1
END
BONDANGLEBENDTYPE
# CT        CHT       T0
1
   50.000  100.000  109.500
END
BONDANGLE
# IT  JT  KT  ICT
2
  1   2   3   1
  2   3   4   1
END
";

    fn write_tmp(content: &str, suffix: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("newen_test_{suffix}.tmp"));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_parse_cg16_topology() {
        let path = write_tmp(CG16_TOPO, "cg16_topo");
        let parsed = read_topology_file(&path).expect("Failed to parse inline topology");

        assert_eq!(parsed.n_atoms, 4);
        assert_eq!(parsed.bonds.len(), 3);
        assert_eq!(parsed.angles.len(), 2);

        // First atom: IAC=6 in file (1-indexed) → 5 in memory (0-indexed)
        assert_eq!(parsed.iac[0], 5);
        assert!((parsed.masses[0] - 62.0).abs() < 1e-6);
        assert!((parsed.charges[0] - 0.0).abs() < 1e-6);

        // First bond: atoms 1-2 in file → 0-1 (0-indexed), type 1 → 0
        assert_eq!(parsed.bonds[0], (0, 1, 0));

        // Bond parameters
        assert!((parsed.bond_parameters[0].k_harmonic - 500.0).abs() < 1e-6);
        assert!((parsed.bond_parameters[0].r0 - 0.15).abs() < 1e-6);

        // BONDANGLEBENDTYPE conversion (GROMOS in_topology.cc:854-855):
        // CT untouched, CHT *= (180/pi)^2, T0 deg -> rad
        let deg_to_rad_sq = 180.0 * 180.0 / (std::f64::consts::PI * std::f64::consts::PI);
        let angle = &parsed.angle_parameters[0];
        assert!((angle.k_cosine - 50.0).abs() < 1e-6);
        assert!((angle.k_harmonic - 100.0 * deg_to_rad_sq).abs() < 1e-6);
        assert!((angle.theta0 - 109.5_f64.to_radians()).abs() < 1e-12);

        std::fs::remove_file(path).ok();
    }

    /// Minimal topology exercising TORSDIHEDRALTYPE and IMPDIHEDRALTYPE unit conversions.
    const DIHEDRAL_TYPES_TOPO: &str = "\
TITLE
  dihedral/improper type conversion test
END
SOLUTEATOM
# ATNM MRES PANM IAC   MASS      CG   CGC INE IEXCL
2
    1    1   C     6  12.000000  0.000000  0  1  2
                                              0
    2    1   H     1   1.000000  0.000000  0  0
                                              0
END
TORSDIHEDRALTYPE
# CP        PD        NP
1
    5.000   180.000    2
END
IMPDIHEDRALTYPE
# CQ        Q0
1
    0.051     0.000
END
";

    #[test]
    fn test_parse_dihedral_and_improper_type_conversions() {
        let path = write_tmp(DIHEDRAL_TYPES_TOPO, "dihedral_types_topo");
        let parsed = read_topology_file(&path).expect("Failed to parse inline topology");

        let deg_to_rad_sq = 180.0 * 180.0 / (std::f64::consts::PI * std::f64::consts::PI);

        // TORSDIHEDRALTYPE conversion (GROMOS read_block_TORSDIHEDRALTYPE):
        // CP (k) untouched; PD deg -> rad, plus its cosine
        assert_eq!(parsed.dihedral_parameters.len(), 1);
        let dihedral = &parsed.dihedral_parameters[0];
        assert!((dihedral.k - 5.0).abs() < 1e-6);
        assert!((dihedral.pd - 180.0_f64.to_radians()).abs() < 1e-12);
        assert!((dihedral.cospd - 180.0_f64.to_radians().cos()).abs() < 1e-12);
        assert!((dihedral.cospd - (-1.0)).abs() < 1e-12);
        assert_eq!(dihedral.m, 2);

        // IMPDIHEDRALTYPE conversion (GROMOS in_topology.cc:1055-1056):
        // CQ (k) *= (180/pi)^2, Q0 deg -> rad
        assert_eq!(parsed.improper_dihedral_parameters.len(), 1);
        let improper = &parsed.improper_dihedral_parameters[0];
        assert!((improper.k - 0.051 * deg_to_rad_sq).abs() < 1e-6);
        assert!((improper.q0 - 0.0_f64.to_radians()).abs() < 1e-12);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn temperature_groups_are_read_as_gromos_pp_reads_them() {
        // gromos++ InTopology.cc:1405-1429 accepts a block that stops before the last solute
        // atom; only the MD engine refuses it
        let short = format!("{CG16_TOPO}\nTEMPERATUREGROUPS\n1\n3\nEND\n");
        let topo = build_topology(parse_topology(short.as_bytes()).unwrap());
        assert_eq!(topo.temperature_groups, vec![vec![0, 1, 2]]);

        let full = format!("{CG16_TOPO}\nTEMPERATUREGROUPS\n2\n1 4\nEND\n");
        let topo = build_topology(parse_topology(full.as_bytes()).unwrap());
        assert_eq!(topo.temperature_groups, vec![vec![0], vec![1, 2, 3]]);
    }

    #[test]
    fn test_parse_missing_topology_returns_error() {
        let result = read_topology_file("/nonexistent/path/file.topo");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_topology_bonds_are_zero_indexed() {
        let path = write_tmp(CG16_TOPO, "cg16_topo_idx");
        let parsed = read_topology_file(&path).unwrap();

        // All bond indices must be 0-indexed (< n_atoms)
        for (i, j, _t) in &parsed.bonds {
            assert!(*i < parsed.n_atoms);
            assert!(*j < parsed.n_atoms);
        }

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_parse_angles_filters_duplicate_atoms() {
        // Topology with angles that have duplicate atoms (e.g., from bad make_top C-terminus)
        let topo_str = "\
TITLE
  test duplicate-atom angles
END
SOLUTEATOM
4
    1    1   C     6  12.000000  0.000000  0  0
                                              0
    2    1   O     8  16.000000  0.000000  0  0
                                              0
    3    1   O2    8  16.000000  0.000000  0  0
                                              0
    4    1   CA    6  12.000000  0.000000  0  0
                                              0
END
BONDSTRETCHTYPE
1
  1000.000  500.000   0.150
END
BOND
3
  4   1   1
  1   2   1
  1   3   1
END
BONDANGLEBENDTYPE
1
   50.000  100.000  109.500
END
BONDANGLE
# includes two invalid angles (duplicate atoms) that should be filtered
5
  4   1   2   1
  4   1   3   1
  2   1   3   1
  1   1   2   1
  1   1   3   1
END
";
        let path = write_tmp(topo_str, "dup_angle_topo");
        let parsed = read_topology_file(&path).unwrap();

        // Only 3 valid angles should remain (the two with atom 1 repeated are filtered)
        assert_eq!(parsed.angles.len(), 3);
        // Verify none of the parsed angles have duplicate atoms
        for (i, j, k, _) in &parsed.angles {
            assert_ne!(i, j, "angle has i==j");
            assert_ne!(i, k, "angle has i==k");
            assert_ne!(j, k, "angle has j==k");
        }

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_parse_angles_internal_function() {
        // Test parse_angles directly with duplicate-atom entries
        let input = "5\n  1  2  3  1\n  2  2  3  1\n  1  3  3  1\n  4  1  2  1\n  1  1  3  1\n";
        let mut lines = input.lines().map(|l| Ok(l.to_string()));
        let mut angles = Vec::new();

        parse_angles(&mut lines, &mut angles).unwrap();

        // 5 entries declared, but 3 have duplicate atoms → only 2 valid
        assert_eq!(angles.len(), 2);
        assert_eq!(angles[0], (0, 1, 2, 0)); // 1-2-3 → 0-1-2
        assert_eq!(angles[1], (3, 0, 1, 0)); // 4-1-2 → 3-0-1
    }
}
