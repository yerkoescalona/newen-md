//! IFP (Interaction Function Parameter) file parser
//!
//! Parses GROMOS IFP files that define force field parameters:
//! mass codes, bond types, angle types, dihedral types,
//! improper dihedral types, and Lennard-Jones parameters.

use crate::IoError;
use std::fs;
use std::path::Path;

/// Mass atom type: maps a mass code to an actual mass and name.
#[derive(Debug, Clone)]
pub struct MassAtomType {
    pub code: usize,
    pub mass: f64,
    pub name: String,
}

/// Bond stretch type parameters.
#[derive(Debug, Clone)]
pub struct BondStretchType {
    pub code: usize,
    /// Quartic force constant
    pub k_quartic: f64,
    /// Harmonic force constant
    pub k_harmonic: f64,
    /// Ideal bond length
    pub r0: f64,
}

/// Bond angle bend type parameters.
#[derive(Debug, Clone)]
pub struct BondAngleBendType {
    pub code: usize,
    /// Non-harmonic force constant
    pub k_non_harmonic: f64,
    /// Harmonic force constant
    pub k_harmonic: f64,
    /// Ideal bond angle (degrees)
    pub theta0: f64,
}

/// Improper dihedral type parameters.
#[derive(Debug, Clone)]
pub struct ImproperDihedralType {
    pub code: usize,
    /// Force constant
    pub k: f64,
    /// Ideal improper dihedral angle (degrees)
    pub xi0: f64,
}

/// Torsional dihedral type parameters.
#[derive(Debug, Clone)]
pub struct TorsionalDihedralType {
    pub code: usize,
    /// Force constant
    pub k: f64,
    /// Phase shift (degrees)
    pub delta: f64,
    /// Multiplicity
    pub multiplicity: usize,
}

/// Single atom LJ parameters (per atom type).
#[derive(Debug, Clone)]
pub struct SingleAtomLJ {
    /// 1-based atom type code
    pub iac: usize,
    /// Atom type name (e.g., "O", "CH2", "N")
    pub name: String,
    /// sqrt(C6) normal
    pub sqrt_c6: f64,
    /// sqrt(C12) for interaction set 1
    pub sqrt_c12_1: f64,
    /// sqrt(C12) for interaction set 2
    pub sqrt_c12_2: f64,
    /// sqrt(C12) for interaction set 3
    pub sqrt_c12_3: f64,
    /// CS6 for 1-4 interactions
    pub cs6: f64,
    /// CS12 for 1-4 interactions
    pub cs12: f64,
    /// Selection matrix row: which C12 set to use for each pair (1-based indices)
    pub selection: Vec<usize>,
}

/// Mixed atom LJ pair (overrides for specific atom type pairs).
#[derive(Debug, Clone)]
pub struct MixedAtomLJPair {
    pub iac_i: usize,
    pub iac_j: usize,
    pub c6: f64,
    pub c12: f64,
    pub cs6: f64,
    pub cs12: f64,
}

/// Complete parsed IFP file.
#[derive(Debug, Clone)]
pub struct ForceFieldParameters {
    /// Force field name
    pub force_field: String,
    /// Mass atom type table
    pub mass_types: Vec<MassAtomType>,
    /// Bond stretch type parameters
    pub bond_types: Vec<BondStretchType>,
    /// Bond angle bend type parameters
    pub angle_types: Vec<BondAngleBendType>,
    /// Improper dihedral type parameters
    pub improper_types: Vec<ImproperDihedralType>,
    /// Torsional dihedral type parameters
    pub torsion_types: Vec<TorsionalDihedralType>,
    /// Single atom LJ parameters (one per atom type)
    pub atom_lj: Vec<SingleAtomLJ>,
    /// Mixed atom LJ pair overrides
    pub mixed_lj: Vec<MixedAtomLJPair>,
    /// Number of atom types
    pub num_atom_types: usize,
}

impl ForceFieldParameters {
    /// Get mass for a mass code.
    pub fn get_mass(&self, code: usize) -> Option<f64> {
        self.mass_types
            .iter()
            .find(|m| m.code == code)
            .map(|m| m.mass)
    }

    /// Get bond stretch parameters by type code (1-based).
    pub fn get_bond_type(&self, code: usize) -> Option<&BondStretchType> {
        self.bond_types.iter().find(|b| b.code == code)
    }

    /// Get angle bend parameters by type code (1-based).
    pub fn get_angle_type(&self, code: usize) -> Option<&BondAngleBendType> {
        self.angle_types.iter().find(|a| a.code == code)
    }

    /// Get improper dihedral parameters by type code (1-based).
    pub fn get_improper_type(&self, code: usize) -> Option<&ImproperDihedralType> {
        self.improper_types.iter().find(|d| d.code == code)
    }

    /// Get torsional dihedral parameters by type code (1-based).
    pub fn get_torsion_type(&self, code: usize) -> Option<&TorsionalDihedralType> {
        self.torsion_types.iter().find(|d| d.code == code)
    }

    /// Get atom LJ parameters by IAC (1-based).
    pub fn get_atom_lj(&self, iac: usize) -> Option<&SingleAtomLJ> {
        self.atom_lj.iter().find(|a| a.iac == iac)
    }

    /// Compute C6 for a pair of atom types.
    pub fn compute_c6(&self, iac_i: usize, iac_j: usize) -> f64 {
        // Check mixed pairs first
        for mp in &self.mixed_lj {
            if (mp.iac_i == iac_i && mp.iac_j == iac_j) || (mp.iac_i == iac_j && mp.iac_j == iac_i)
            {
                return mp.c6;
            }
        }
        let ai = self.get_atom_lj(iac_i);
        let aj = self.get_atom_lj(iac_j);
        match (ai, aj) {
            (Some(a), Some(b)) => a.sqrt_c6 * b.sqrt_c6,
            _ => 0.0,
        }
    }

    /// Compute C12 for a pair of atom types (using selection matrix).
    pub fn compute_c12(&self, iac_i: usize, iac_j: usize) -> f64 {
        // Check mixed pairs first
        for mp in &self.mixed_lj {
            if (mp.iac_i == iac_i && mp.iac_j == iac_j) || (mp.iac_i == iac_j && mp.iac_j == iac_i)
            {
                return mp.c12;
            }
        }
        let ai = self.get_atom_lj(iac_i);
        let aj = self.get_atom_lj(iac_j);
        match (ai, aj) {
            (Some(a), Some(b)) => {
                // Use selection matrix from atom i to determine which C12 set for j
                let j_idx = iac_j - 1;
                let set = if j_idx < a.selection.len() {
                    a.selection[j_idx]
                } else {
                    1
                };
                let c12_i = match set {
                    1 => a.sqrt_c12_1,
                    2 => a.sqrt_c12_2,
                    3 => a.sqrt_c12_3,
                    _ => a.sqrt_c12_1,
                };
                let i_idx = iac_i - 1;
                let set_j = if i_idx < b.selection.len() {
                    b.selection[i_idx]
                } else {
                    1
                };
                let c12_j = match set_j {
                    1 => b.sqrt_c12_1,
                    2 => b.sqrt_c12_2,
                    3 => b.sqrt_c12_3,
                    _ => b.sqrt_c12_1,
                };
                c12_i * c12_j
            },
            _ => 0.0,
        }
    }

    /// Compute C6/C12 for 1-4 interactions.
    pub fn compute_c6_14(&self, iac_i: usize, iac_j: usize) -> f64 {
        for mp in &self.mixed_lj {
            if (mp.iac_i == iac_i && mp.iac_j == iac_j) || (mp.iac_i == iac_j && mp.iac_j == iac_i)
            {
                return mp.cs6;
            }
        }
        let ai = self.get_atom_lj(iac_i);
        let aj = self.get_atom_lj(iac_j);
        match (ai, aj) {
            (Some(a), Some(b)) => a.cs6 * b.cs6,
            _ => 0.0,
        }
    }

    pub fn compute_c12_14(&self, iac_i: usize, iac_j: usize) -> f64 {
        for mp in &self.mixed_lj {
            if (mp.iac_i == iac_i && mp.iac_j == iac_j) || (mp.iac_i == iac_j && mp.iac_j == iac_i)
            {
                return mp.cs12;
            }
        }
        let ai = self.get_atom_lj(iac_i);
        let aj = self.get_atom_lj(iac_j);
        match (ai, aj) {
            (Some(a), Some(b)) => a.cs12 * b.cs12,
            _ => 0.0,
        }
    }
}

/// Parse an IFP file from a path.
pub fn read_ifp_file<P: AsRef<Path>>(path: P) -> Result<ForceFieldParameters, IoError> {
    let content = fs::read_to_string(path.as_ref()).map_err(|e| {
        IoError::FileNotFound(format!(
            "Cannot read IFP file '{}': {}",
            path.as_ref().display(),
            e
        ))
    })?;
    parse_ifp(&content)
}

/// Parse IFP content from a string.
pub fn parse_ifp(content: &str) -> Result<ForceFieldParameters, IoError> {
    let mut result = ForceFieldParameters {
        force_field: String::new(),
        mass_types: Vec::new(),
        bond_types: Vec::new(),
        angle_types: Vec::new(),
        improper_types: Vec::new(),
        torsion_types: Vec::new(),
        atom_lj: Vec::new(),
        mixed_lj: Vec::new(),
        num_atom_types: 0,
    };

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
            "MASSATOMTYPECODE" => {
                i += 1;
                // Read NRMATY NMATY header
                let mut header_read = false;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    i += 1;
                    if l.is_empty() || l.starts_with('#') {
                        continue;
                    }
                    if !header_read {
                        // First data line: NRMATY NMATY
                        header_read = true;
                        continue;
                    }
                    // Data lines: N ATMAS ATMASN
                    let tokens: Vec<&str> = l.split_whitespace().collect();
                    if tokens.len() >= 3 {
                        let code: usize = tokens[0].parse().unwrap_or(0);
                        let mass: f64 = tokens[1].parse().unwrap_or(0.0);
                        let name = tokens[2].to_string();
                        result.mass_types.push(MassAtomType { code, mass, name });
                    }
                }
            },
            "BONDSTRETCHTYPECODE" => {
                i += 1;
                let mut header_read = false;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    i += 1;
                    if l.is_empty() || l.starts_with('#') {
                        continue;
                    }
                    if !header_read {
                        header_read = true;
                        continue;
                    }
                    let tokens: Vec<&str> = l.split_whitespace().collect();
                    if tokens.len() >= 4 {
                        let code: usize = tokens[0].parse().unwrap_or(0);
                        let k_quartic: f64 = tokens[1].parse().unwrap_or(0.0);
                        let k_harmonic: f64 = tokens[2].parse().unwrap_or(0.0);
                        let r0: f64 = tokens[3].parse().unwrap_or(0.0);
                        result.bond_types.push(BondStretchType {
                            code,
                            k_quartic,
                            k_harmonic,
                            r0,
                        });
                    }
                }
            },
            "BONDANGLEBENDTYPECODE" => {
                i += 1;
                let mut header_read = false;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    i += 1;
                    if l.is_empty() || l.starts_with('#') {
                        continue;
                    }
                    if !header_read {
                        header_read = true;
                        continue;
                    }
                    let tokens: Vec<&str> = l.split_whitespace().collect();
                    if tokens.len() >= 4 {
                        let code: usize = tokens[0].parse().unwrap_or(0);
                        let k_non_harmonic: f64 = tokens[1].parse().unwrap_or(0.0);
                        let k_harmonic: f64 = tokens[2].parse().unwrap_or(0.0);
                        let theta0: f64 = tokens[3].parse().unwrap_or(0.0);
                        result.angle_types.push(BondAngleBendType {
                            code,
                            k_non_harmonic,
                            k_harmonic,
                            theta0,
                        });
                    }
                }
            },
            "IMPDIHEDRALTYPECODE" => {
                i += 1;
                let mut header_read = false;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    i += 1;
                    if l.is_empty() || l.starts_with('#') {
                        continue;
                    }
                    if !header_read {
                        header_read = true;
                        continue;
                    }
                    let tokens: Vec<&str> = l.split_whitespace().collect();
                    if tokens.len() >= 3 {
                        let code: usize = tokens[0].parse().unwrap_or(0);
                        let k: f64 = tokens[1].parse().unwrap_or(0.0);
                        let xi0: f64 = tokens[2].parse().unwrap_or(0.0);
                        result
                            .improper_types
                            .push(ImproperDihedralType { code, k, xi0 });
                    }
                }
            },
            "TORSDIHEDRALTYPECODE" => {
                i += 1;
                let mut header_read = false;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    i += 1;
                    if l.is_empty() || l.starts_with('#') {
                        continue;
                    }
                    if !header_read {
                        header_read = true;
                        continue;
                    }
                    let tokens: Vec<&str> = l.split_whitespace().collect();
                    if tokens.len() >= 4 {
                        let code: usize = tokens[0].parse().unwrap_or(0);
                        let k: f64 = tokens[1].parse().unwrap_or(0.0);
                        let delta: f64 = tokens[2].parse().unwrap_or(0.0);
                        let multiplicity: usize = tokens[3].parse().unwrap_or(0);
                        result.torsion_types.push(TorsionalDihedralType {
                            code,
                            k,
                            delta,
                            multiplicity,
                        });
                    }
                }
            },
            "SINGLEATOMLJPAIR" => {
                i += 1;
                let mut num_types = 0usize;
                let mut num_types_read = false;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    i += 1;
                    if l.is_empty() || l.starts_with('#') {
                        continue;
                    }
                    if !num_types_read {
                        num_types = l.trim().parse().unwrap_or(0);
                        result.num_atom_types = num_types;
                        num_types_read = true;
                        continue;
                    }

                    // Each atom type spans 4 lines:
                    // Line 1: IAC TYPE SQRT(C6) SQRT(C12(1)) SQRT(C12(2)) SQRT(C12(3))
                    // Line 2: #CS6 CS12 parameters LJ14PAIR  (comment - already skipped)
                    // Line 3: CS6 CS12
                    // Lines 4+: selection matrix rows (may span multiple lines)

                    let tokens: Vec<&str> = l.split_whitespace().collect();
                    if tokens.len() >= 6 {
                        let iac: usize = tokens[0].parse().unwrap_or(0);
                        let name = tokens[1].to_string();
                        let sqrt_c6: f64 = tokens[2].parse().unwrap_or(0.0);
                        let sqrt_c12_1: f64 = tokens[3].parse().unwrap_or(0.0);
                        let sqrt_c12_2: f64 = tokens[4].parse().unwrap_or(0.0);
                        let sqrt_c12_3: f64 = tokens[5].parse().unwrap_or(0.0);

                        // Next data line: CS6 CS12
                        let mut cs6 = 0.0;
                        let mut cs12 = 0.0;
                        while i < lines.len() && lines[i].trim() != "END" {
                            let l2 = lines[i].trim();
                            i += 1;
                            if l2.is_empty() || l2.starts_with('#') {
                                continue;
                            }
                            let t2: Vec<&str> = l2.split_whitespace().collect();
                            if t2.len() >= 2 {
                                cs6 = t2[0].parse().unwrap_or(0.0);
                                cs12 = t2[1].parse().unwrap_or(0.0);
                            }
                            break;
                        }

                        // Read selection matrix values (num_types values, may span multiple lines)
                        let mut selection = Vec::with_capacity(num_types);
                        while selection.len() < num_types && i < lines.len() {
                            let l3 = lines[i].trim();
                            if l3 == "END" {
                                break;
                            }
                            i += 1;
                            if l3.is_empty() || l3.starts_with('#') {
                                // After selection matrix, we may hit #--- separator
                                if l3.starts_with("#---") {
                                    break;
                                }
                                continue;
                            }
                            for token in l3.split_whitespace() {
                                if let Ok(v) = token.parse::<usize>() {
                                    selection.push(v);
                                }
                            }
                        }

                        // Skip any remaining comment lines until #--- or next atom
                        while i < lines.len() {
                            let l4 = lines[i].trim();
                            if l4 == "END" || l4.starts_with("#---") {
                                if l4.starts_with("#---") {
                                    i += 1;
                                }
                                break;
                            }
                            if !l4.is_empty() && !l4.starts_with('#') {
                                break;
                            }
                            i += 1;
                        }

                        result.atom_lj.push(SingleAtomLJ {
                            iac,
                            name,
                            sqrt_c6,
                            sqrt_c12_1,
                            sqrt_c12_2,
                            sqrt_c12_3,
                            cs6,
                            cs12,
                            selection,
                        });
                    }
                }
            },
            "MIXEDATOMLJPAIR" => {
                i += 1;
                while i < lines.len() && lines[i].trim() != "END" {
                    let l = lines[i].trim();
                    i += 1;
                    if l.is_empty() || l.starts_with('#') {
                        continue;
                    }
                    let tokens: Vec<&str> = l.split_whitespace().collect();
                    if tokens.len() >= 6 {
                        let iac_i: usize = tokens[0].parse().unwrap_or(0);
                        let iac_j: usize = tokens[1].parse().unwrap_or(0);
                        let c6: f64 = tokens[2].parse().unwrap_or(0.0);
                        let c12: f64 = tokens[3].parse().unwrap_or(0.0);
                        let cs6: f64 = tokens[4].parse().unwrap_or(0.0);
                        let cs12: f64 = tokens[5].parse().unwrap_or(0.0);
                        result.mixed_lj.push(MixedAtomLJPair {
                            iac_i,
                            iac_j,
                            c6,
                            c12,
                            cs6,
                            cs12,
                        });
                    }
                }
            },
            _ => {},
        }
        i += 1;
    }

    Ok(result)
}
