//! MTB (Molecular Topology Building Block) file parser
//!
//! Parses GROMOS MTB files that define building blocks for residues,
//! terminal groups, and solvents.

use crate::IoError;
use std::fs;
use std::path::Path;

/// A single atom in a building block.
#[derive(Debug, Clone)]
pub struct BbAtom {
    /// 1-based atom number within the building block
    pub number: i32,
    /// Atom name (e.g., "N", "CA", "CB")
    pub name: String,
    /// Integer atom type code (IAC) — index into force field atom types
    pub iac: i32,
    /// Mass code — index into MASSATOMTYPECODE table
    pub mass_code: usize,
    /// Partial charge
    pub charge: f64,
    /// Charge group membership (1 = end of charge group)
    pub chargegroup: usize,
    /// Exclusion list: atom numbers excluded from non-bonded interactions
    pub exclusions: Vec<i32>,
    /// Expected number of exclusions (MAE field, used for continuation line parsing)
    pub expected_exclusions: usize,
}

/// Preceding exclusion entry (atoms from previous residue).
#[derive(Debug, Clone)]
pub struct PrecedingExclusion {
    /// Atom index (-NLIN to 0, referring to atoms in previous residue or overlap)
    pub atom: i32,
    /// Number of forward exclusions
    pub num_exclusions: usize,
    /// Forward-excluded atom numbers (within current building block)
    pub excluded_atoms: Vec<i32>,
}

/// A solute building block (amino acid, nucleotide, ion, etc.).
#[derive(Debug, Clone)]
pub struct BbSolute {
    /// Residue name (e.g., "ALA", "GLY", "NA+")
    pub name: String,
    /// All atoms (NMAT total, 1-indexed in the file)
    pub atoms: Vec<BbAtom>,
    /// Number of preceding exclusions (NLIN)
    pub num_preceding_exclusions: usize,
    /// Preceding exclusion specifications
    pub preceding_exclusions: Vec<PrecedingExclusion>,
    /// Bonds: (atom_i, atom_j, bond_type_code). Atom indices can be negative (linking).
    pub bonds: Vec<(i32, i32, usize)>,
    /// Angles: (atom_i, atom_j, atom_k, angle_type_code)
    pub angles: Vec<(i32, i32, i32, usize)>,
    /// Improper dihedrals: (i, j, k, l, type_code)
    pub improper_dihedrals: Vec<(i32, i32, i32, i32, usize)>,
    /// Proper dihedrals: (i, j, k, l, type_code)
    pub proper_dihedrals: Vec<(i32, i32, i32, i32, i32)>,
    /// LJ exceptions: (atom_i, atom_j, type_code)
    pub lj_exceptions: Vec<(i32, i32, usize)>,
}

/// An end-group building block (N-terminal or C-terminal patch).
#[derive(Debug, Clone)]
pub struct BbEnd {
    /// End group name (e.g., "NH3+", "COO-")
    pub name: String,
    /// All atoms (NMAT total)
    pub atoms: Vec<BbAtom>,
    /// Number of replacing atoms (NREP):
    /// - Positive: last NREP atoms replace first NREP of next residue (N-terminal)
    /// - Negative: replaces last |NREP| atoms of previous residue (C-terminal)
    pub n_replace: i32,
    /// Bonds
    pub bonds: Vec<(i32, i32, usize)>,
    /// Angles
    pub angles: Vec<(i32, i32, i32, usize)>,
    /// Improper dihedrals
    pub improper_dihedrals: Vec<(i32, i32, i32, i32, usize)>,
    /// Proper dihedrals
    pub proper_dihedrals: Vec<(i32, i32, i32, i32, i32)>,
    /// LJ exceptions
    pub lj_exceptions: Vec<(i32, i32, usize)>,
}

/// A solvent building block (e.g., SPC water).
#[derive(Debug, Clone)]
pub struct BbSolvent {
    /// Solvent name (e.g., "H2O")
    pub name: String,
    /// Atoms: (name, iac, mass_code, charge)
    pub atoms: Vec<(String, usize, usize, f64)>,
    /// Constraints: (atom_i, atom_j, distance) — 1-indexed
    pub constraints: Vec<(usize, usize, f64)>,
}

/// Complete parsed MTB file.
#[derive(Debug, Clone)]
pub struct BuildingBlocks {
    /// Force field code (e.g., "54A7")
    pub force_field: String,
    /// Physical constants
    pub fpepsi: f64,
    pub hbar: f64,
    pub spdl: f64,
    pub boltz: f64,
    /// Number of nearest-neighbour exclusions when linking
    pub link_exclusions: usize,
    /// Solute building blocks
    pub solute_blocks: Vec<BbSolute>,
    /// End-group building blocks
    pub end_blocks: Vec<BbEnd>,
    /// Solvent building blocks
    pub solvent_blocks: Vec<BbSolvent>,
}

impl BuildingBlocks {
    /// Find a solute building block by name.
    pub fn find_solute(&self, name: &str) -> Option<&BbSolute> {
        self.solute_blocks.iter().find(|bb| bb.name == name)
    }

    /// Find an end-group building block by name.
    pub fn find_end(&self, name: &str) -> Option<&BbEnd> {
        self.end_blocks.iter().find(|bb| bb.name == name)
    }

    /// Find a solvent building block by name.
    pub fn find_solvent(&self, name: &str) -> Option<&BbSolvent> {
        self.solvent_blocks.iter().find(|bb| bb.name == name)
    }

    /// Find a building block by name. Returns (index, is_end_group).
    /// Searches end groups first, then solute blocks.
    pub fn find_block(&self, name: &str) -> Option<BlockRef<'_>> {
        if let Some(bb) = self.end_blocks.iter().find(|bb| bb.name == name) {
            return Some(BlockRef::End(bb));
        }
        if let Some(bb) = self.solute_blocks.iter().find(|bb| bb.name == name) {
            return Some(BlockRef::Solute(bb));
        }
        None
    }
}

/// Reference to a found building block.
#[derive(Debug)]
pub enum BlockRef<'a> {
    Solute(&'a BbSolute),
    End(&'a BbEnd),
}

/// Parse an MTB file from a path.
pub fn read_mtb_file<P: AsRef<Path>>(path: P) -> Result<BuildingBlocks, IoError> {
    let content = fs::read_to_string(path.as_ref()).map_err(|e| {
        IoError::FileNotFound(format!(
            "Cannot read MTB file '{}': {}",
            path.as_ref().display(),
            e
        ))
    })?;
    parse_mtb(&content)
}

/// Parse MTB content from a string.
pub fn parse_mtb(content: &str) -> Result<BuildingBlocks, IoError> {
    let mut result = BuildingBlocks {
        force_field: String::new(),
        fpepsi: newen_core::units::four_pi_eps_i,
        hbar: newen_core::units::hBar,
        spdl: newen_core::units::spd_l,
        boltz: newen_core::units::kB,
        link_exclusions: 2,
        solute_blocks: Vec::new(),
        end_blocks: Vec::new(),
        solvent_blocks: Vec::new(),
    };

    // Split into blocks
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();

        match line {
            "FORCEFIELD" => {
                i += 1;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    if !l.is_empty() && !l.starts_with('#') {
                        result.force_field = l.to_string();
                    }
                    i += 1;
                }
            },
            "PHYSICALCONSTANTS" | "TOPPHYSCON" => {
                i += 1;
                let mut vals = Vec::new();
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    if !l.is_empty() && !l.starts_with('#') {
                        if let Ok(v) = l.parse::<f64>() {
                            vals.push(v);
                        }
                    }
                    i += 1;
                }
                if vals.len() >= 4 {
                    result.fpepsi = vals[0];
                    result.hbar = vals[1];
                    result.spdl = vals[2];
                    result.boltz = vals[3];
                }
            },
            "LINKEXCLUSIONS" => {
                i += 1;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    if !l.is_empty() && !l.starts_with('#') {
                        if let Ok(v) = l.parse::<usize>() {
                            result.link_exclusions = v;
                        }
                    }
                    i += 1;
                }
            },
            "MTBUILDBLSOLUTE" => {
                let bb = parse_solute_block(&lines, &mut i)?;
                result.solute_blocks.push(bb);
            },
            "MTBUILDBLEND" => {
                let bb = parse_end_block(&lines, &mut i)?;
                result.end_blocks.push(bb);
            },
            "MTBUILDBLSOLVENT" => {
                let bb = parse_solvent_block(&lines, &mut i)?;
                result.solvent_blocks.push(bb);
            },
            _ => {},
        }
        i += 1;
    }

    Ok(result)
}

/// Parse data lines: skip comments and empty lines, collect non-comment tokens.
fn next_data_line<'a>(lines: &[&'a str], i: &mut usize) -> Option<&'a str> {
    while *i < lines.len() {
        let line = lines[*i].trim();
        if line == "END" {
            return None;
        }
        *i += 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        return Some(line);
    }
    None
}

fn parse_solute_block(lines: &[&str], i: &mut usize) -> Result<BbSolute, IoError> {
    *i += 1; // skip MTBUILDBLSOLUTE line

    // Skip comment lines starting with #@ (block type annotations)
    while *i < lines.len() {
        let l = lines[*i].trim();
        if l.starts_with("#@") || l.starts_with("# building") || l.starts_with("# RNME") {
            *i += 1;
            continue;
        }
        break;
    }

    // Read name
    let name = match next_data_line(lines, i) {
        Some(l) => l.to_string(),
        None => return Err(IoError::ParseError("Expected residue name".to_string())),
    };

    // Skip "# number of atoms" comment
    // Read NMAT NLIN
    let (nmat, nlin) = match next_data_line(lines, i) {
        Some(l) => {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() < 2 {
                return Err(IoError::ParseError(format!(
                    "Expected NMAT NLIN in block {}, got: {}",
                    name, l
                )));
            }
            (
                tokens[0].parse::<usize>().map_err(|_| {
                    IoError::ParseError(format!("Invalid NMAT in {}: {}", name, tokens[0]))
                })?,
                tokens[1].parse::<usize>().map_err(|_| {
                    IoError::ParseError(format!("Invalid NLIN in {}: {}", name, tokens[1]))
                })?,
            )
        },
        None => {
            return Err(IoError::ParseError(format!(
                "Expected NMAT NLIN in {}",
                name
            )))
        },
    };

    // Read preceding exclusions
    let mut preceding_exclusions = Vec::new();
    for _ in 0..nlin {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 2 {
                let atom: i32 = tokens[0].parse().map_err(|_| {
                    IoError::ParseError(format!("Invalid pexcl atom in {}: {}", name, tokens[0]))
                })?;
                let mae: usize = tokens[1].parse().map_err(|_| {
                    IoError::ParseError(format!("Invalid MAE in {}: {}", name, tokens[1]))
                })?;
                let mut excluded = Vec::new();
                for tok in tokens.iter().take(2 + mae).skip(2) {
                    if let Ok(v) = tok.parse::<i32>() {
                        excluded.push(v);
                    }
                }
                preceding_exclusions.push(PrecedingExclusion {
                    atom,
                    num_exclusions: mae,
                    excluded_atoms: excluded,
                });
            }
        }
    }

    // Read atoms
    let mut atoms = Vec::with_capacity(nmat);
    for _ in 0..nmat {
        if let Some(l) = next_data_line(lines, i) {
            let mut atom = parse_bb_atom(l, &name)?;
            // Handle continuation lines for exclusions
            while atom.exclusions.len() < atom.expected_exclusions {
                if let Some(cont) = next_data_line(lines, i) {
                    for token in cont.split_whitespace() {
                        if let Ok(v) = token.parse::<i32>() {
                            atom.exclusions.push(v);
                        }
                    }
                } else {
                    break;
                }
            }
            atoms.push(atom);
        }
    }

    // Read bonds
    let bonds = parse_bonds(lines, i, &name)?;

    // Read angles
    let angles = parse_angles(lines, i, &name)?;

    // Read improper dihedrals
    let impropers = parse_improper_dihedrals(lines, i, &name)?;

    // Read proper dihedrals
    let dihedrals = parse_proper_dihedrals(lines, i, &name)?;

    // Skip #@FREELINE if present
    while *i < lines.len() {
        let l = lines[*i].trim();
        if l.starts_with("#@") {
            *i += 1;
            continue;
        }
        break;
    }

    // Read LJ exceptions
    let lj_exceptions = parse_lj_exceptions(lines, i, &name)?;

    // Skip to END
    while *i < lines.len() && lines[*i].trim() != "END" {
        *i += 1;
    }

    Ok(BbSolute {
        name,
        atoms,
        num_preceding_exclusions: nlin,
        preceding_exclusions,
        bonds,
        angles,
        improper_dihedrals: impropers,
        proper_dihedrals: dihedrals,
        lj_exceptions,
    })
}

fn parse_end_block(lines: &[&str], i: &mut usize) -> Result<BbEnd, IoError> {
    *i += 1; // skip MTBUILDBLEND line

    // Skip comment lines
    while *i < lines.len() {
        let l = lines[*i].trim();
        if l.starts_with("#@") || l.starts_with("# building") || l.starts_with("# RNME") {
            *i += 1;
            continue;
        }
        break;
    }

    // Read name
    let name = match next_data_line(lines, i) {
        Some(l) => l.to_string(),
        None => return Err(IoError::ParseError("Expected end group name".to_string())),
    };

    // Read NMAT NREP
    let (nmat, nrep) = match next_data_line(lines, i) {
        Some(l) => {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() < 2 {
                return Err(IoError::ParseError(format!(
                    "Expected NMAT NREP in end block {}, got: {}",
                    name, l
                )));
            }
            (
                tokens[0].parse::<usize>().map_err(|_| {
                    IoError::ParseError(format!("Invalid NMAT in {}: {}", name, tokens[0]))
                })?,
                tokens[1].parse::<i32>().map_err(|_| {
                    IoError::ParseError(format!("Invalid NREP in {}: {}", name, tokens[1]))
                })?,
            )
        },
        None => {
            return Err(IoError::ParseError(format!(
                "Expected NMAT NREP in {}",
                name
            )))
        },
    };

    // Read atoms (all NMAT atoms, some are replacing)
    let mut atoms = Vec::with_capacity(nmat);
    for _ in 0..nmat {
        if let Some(l) = next_data_line(lines, i) {
            let atom = parse_bb_atom(l, &name)?;
            atoms.push(atom);
        }
    }

    // Read bonds
    let bonds = parse_bonds(lines, i, &name)?;

    // Read angles
    let angles = parse_angles(lines, i, &name)?;

    // Read improper dihedrals
    let impropers = parse_improper_dihedrals(lines, i, &name)?;

    // Read proper dihedrals
    let dihedrals = parse_proper_dihedrals(lines, i, &name)?;

    // Skip #@FREELINE
    while *i < lines.len() {
        let l = lines[*i].trim();
        if l.starts_with("#@") {
            *i += 1;
            continue;
        }
        break;
    }

    // Read LJ exceptions
    let lj_exceptions = parse_lj_exceptions(lines, i, &name)?;

    // Skip to END
    while *i < lines.len() && lines[*i].trim() != "END" {
        *i += 1;
    }

    Ok(BbEnd {
        name,
        atoms,
        n_replace: nrep,
        bonds,
        angles,
        improper_dihedrals: impropers,
        proper_dihedrals: dihedrals,
        lj_exceptions,
    })
}

fn parse_solvent_block(lines: &[&str], i: &mut usize) -> Result<BbSolvent, IoError> {
    *i += 1; // skip MTBUILDBLSOLVENT line

    // Skip comments
    while *i < lines.len() {
        let l = lines[*i].trim();
        if l.starts_with("#@") || l.starts_with("#solvent") || l.starts_with("#RNMES") {
            *i += 1;
            continue;
        }
        break;
    }

    // Read name
    let name = match next_data_line(lines, i) {
        Some(l) => l.to_string(),
        None => return Err(IoError::ParseError("Expected solvent name".to_string())),
    };

    // Read number of atoms
    let natoms = match next_data_line(lines, i) {
        Some(l) => l.trim().parse::<usize>().map_err(|_| {
            IoError::ParseError(format!("Invalid natoms in solvent {}: {}", name, l))
        })?,
        None => {
            return Err(IoError::ParseError(format!(
                "Expected natoms in solvent {}",
                name
            )))
        },
    };

    // Read atoms
    let mut atoms = Vec::with_capacity(natoms);
    for _ in 0..natoms {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 5 {
                let _num: usize = tokens[0].parse().unwrap_or(0);
                let aname = tokens[1].to_string();
                let iac: usize = tokens[2].parse().unwrap_or(0);
                let mass: usize = tokens[3].parse().unwrap_or(0);
                let charge: f64 = tokens[4].parse().unwrap_or(0.0);
                atoms.push((aname, iac, mass, charge));
            }
        }
    }

    // Read constraints
    let nconstr = match next_data_line(lines, i) {
        Some(l) => l.trim().parse::<usize>().unwrap_or(0),
        None => 0,
    };

    let mut constraints = Vec::with_capacity(nconstr);
    for _ in 0..nconstr {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 3 {
                let ai: usize = tokens[0].parse().unwrap_or(0);
                let aj: usize = tokens[1].parse().unwrap_or(0);
                let dist: f64 = tokens[2].parse().unwrap_or(0.0);
                constraints.push((ai, aj, dist));
            }
        }
    }

    // Skip to END
    while *i < lines.len() && lines[*i].trim() != "END" {
        *i += 1;
    }

    Ok(BbSolvent {
        name,
        atoms,
        constraints,
    })
}

/// Parse a building block atom line.
/// Regular: ATOM_NUM ANM IACM MASS CHARGE CGMICGM MAE MSAE...
/// Trailing: ATOM_NUM ANM IACM MASS CHARGE CGMICGM (no MAE)
fn parse_bb_atom(line: &str, block_name: &str) -> Result<BbAtom, IoError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 6 {
        return Err(IoError::ParseError(format!(
            "Not enough fields for atom in {}: {}",
            block_name, line
        )));
    }

    let number: i32 = tokens[0].parse().map_err(|_| {
        IoError::ParseError(format!(
            "Invalid atom number in {}: {}",
            block_name, tokens[0]
        ))
    })?;
    let name = tokens[1].to_string();
    let iac: i32 = tokens[2].parse().map_err(|_| {
        IoError::ParseError(format!("Invalid IAC in {}: {}", block_name, tokens[2]))
    })?;
    let mass_code: usize = tokens[3].parse().map_err(|_| {
        IoError::ParseError(format!(
            "Invalid mass code in {}: {}",
            block_name, tokens[3]
        ))
    })?;
    let charge: f64 = tokens[4].parse().map_err(|_| {
        IoError::ParseError(format!("Invalid charge in {}: {}", block_name, tokens[4]))
    })?;
    let chargegroup: usize = tokens[5].parse().map_err(|_| {
        IoError::ParseError(format!(
            "Invalid chargegroup in {}: {}",
            block_name, tokens[5]
        ))
    })?;

    // Read exclusions if present (trailing atoms have only 6 fields)
    let mut exclusions = Vec::new();
    let mut expected_exclusions = 0;
    if tokens.len() > 6 {
        let mae: usize = tokens[6].parse().unwrap_or(0);
        expected_exclusions = mae;
        for tok in tokens.iter().skip(7) {
            if let Ok(v) = tok.parse::<i32>() {
                exclusions.push(v);
            }
        }
    }

    Ok(BbAtom {
        number,
        name,
        iac,
        mass_code,
        charge,
        chargegroup,
        exclusions,
        expected_exclusions,
    })
}

fn parse_bonds(
    lines: &[&str],
    i: &mut usize,
    name: &str,
) -> Result<Vec<(i32, i32, usize)>, IoError> {
    // Read NB (number of bonds)
    let nb = match next_data_line(lines, i) {
        Some(l) => l
            .trim()
            .parse::<usize>()
            .map_err(|_| IoError::ParseError(format!("Invalid NB in {}: {}", name, l)))?,
        None => return Ok(Vec::new()),
    };

    let mut bonds = Vec::with_capacity(nb);
    for _ in 0..nb {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 3 {
                let ai: i32 = tokens[0].parse().unwrap_or(0);
                let aj: i32 = tokens[1].parse().unwrap_or(0);
                let btype: usize = tokens[2].parse().unwrap_or(0);
                bonds.push((ai, aj, btype));
            }
        }
    }
    Ok(bonds)
}

fn parse_angles(
    lines: &[&str],
    i: &mut usize,
    name: &str,
) -> Result<Vec<(i32, i32, i32, usize)>, IoError> {
    let na = match next_data_line(lines, i) {
        Some(l) => l
            .trim()
            .parse::<usize>()
            .map_err(|_| IoError::ParseError(format!("Invalid NBA in {}: {}", name, l)))?,
        None => return Ok(Vec::new()),
    };

    let mut angles = Vec::with_capacity(na);
    for _ in 0..na {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 4 {
                let ai: i32 = tokens[0].parse().unwrap_or(0);
                let aj: i32 = tokens[1].parse().unwrap_or(0);
                let ak: i32 = tokens[2].parse().unwrap_or(0);
                let atype: usize = tokens[3].parse().unwrap_or(0);
                angles.push((ai, aj, ak, atype));
            }
        }
    }
    Ok(angles)
}

fn parse_improper_dihedrals(
    lines: &[&str],
    i: &mut usize,
    name: &str,
) -> Result<Vec<(i32, i32, i32, i32, usize)>, IoError> {
    let n = match next_data_line(lines, i) {
        Some(l) => l
            .trim()
            .parse::<usize>()
            .map_err(|_| IoError::ParseError(format!("Invalid NIDA in {}: {}", name, l)))?,
        None => return Ok(Vec::new()),
    };

    let mut result = Vec::with_capacity(n);
    for _ in 0..n {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 5 {
                let ai: i32 = tokens[0].parse().unwrap_or(0);
                let aj: i32 = tokens[1].parse().unwrap_or(0);
                let ak: i32 = tokens[2].parse().unwrap_or(0);
                let al: i32 = tokens[3].parse().unwrap_or(0);
                let dtype: usize = tokens[4].parse().unwrap_or(0);
                result.push((ai, aj, ak, al, dtype));
            }
        }
    }
    Ok(result)
}

fn parse_proper_dihedrals(
    lines: &[&str],
    i: &mut usize,
    name: &str,
) -> Result<Vec<(i32, i32, i32, i32, i32)>, IoError> {
    let n = match next_data_line(lines, i) {
        Some(l) => l
            .trim()
            .parse::<usize>()
            .map_err(|_| IoError::ParseError(format!("Invalid NDA in {}: {}", name, l)))?,
        None => return Ok(Vec::new()),
    };

    let mut result = Vec::with_capacity(n);
    for _ in 0..n {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 5 {
                let ai: i32 = tokens[0].parse().unwrap_or(0);
                let aj: i32 = tokens[1].parse().unwrap_or(0);
                let ak: i32 = tokens[2].parse().unwrap_or(0);
                let al: i32 = tokens[3].parse().unwrap_or(0);
                // Type can be negative (for deletion)
                let dtype: i32 = tokens[4].parse().unwrap_or(0);
                result.push((ai, aj, ak, al, dtype));
            }
        }
    }
    Ok(result)
}

fn parse_lj_exceptions(
    lines: &[&str],
    i: &mut usize,
    _name: &str,
) -> Result<Vec<(i32, i32, usize)>, IoError> {
    let n = match next_data_line(lines, i) {
        Some(l) => l.trim().parse::<usize>().unwrap_or(0),
        None => return Ok(Vec::new()),
    };

    let mut result = Vec::with_capacity(n);
    for _ in 0..n {
        if let Some(l) = next_data_line(lines, i) {
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() >= 3 {
                let ai: i32 = tokens[0].parse().unwrap_or(0);
                let aj: i32 = tokens[1].parse().unwrap_or(0);
                let ltype: usize = tokens[2].parse().unwrap_or(0);
                result.push((ai, aj, ltype));
            }
        }
    }
    Ok(result)
}
