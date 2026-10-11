//! Parser for GROMOS distance restraint files (.distrest).
//!
//! Reads DISTANCERESSPEC and PERTDISRESSPEC blocks.
//! Only virtual atom type 0 (explicit atom) is supported; other types are
//! logged as warnings and skipped in DISTANCERESSPEC, and are an error in PERTDISRESSPEC
//! (a perturbed restraint left out would silently change the energy and dH/dλ).

use std::io::BufRead;
use std::path::Path;

use gromos_core::topology::{DistanceRestraintSpec, PerturbedDistanceRestraintSpec};

use crate::IoError;

/// Parse a GROMOS distance restraint file (.distrest).
///
/// Returns `(unperturbed, perturbed)` lists.  Atoms are converted from
/// 1-indexed (GROMOS) to 0-indexed (Rust).
pub fn read_distanceres<P: AsRef<Path>>(
    path: P,
) -> Result<
    (
        Vec<DistanceRestraintSpec>,
        Vec<PerturbedDistanceRestraintSpec>,
    ),
    IoError,
> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;

    let mut unperturbed = Vec::new();
    let mut perturbed = Vec::new();

    #[derive(PartialEq)]
    enum Block {
        None,
        Dist,
        Pert,
    }
    let mut block = Block::None;
    let mut header_lines = 0usize; // skip DISH/DISC line

    for line in reader.lines() {
        let line = line.map_err(|e| IoError::ParseError(e.to_string()))?;
        let trimmed = line.trim();

        if trimmed == "DISTANCERESSPEC" {
            block = Block::Dist;
            header_lines = 1;
            continue;
        }
        if trimmed == "PERTDISRESSPEC" {
            block = Block::Pert;
            header_lines = 1;
            continue;
        }
        if trimmed == "END" {
            block = Block::None;
            continue;
        }
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if block == Block::None {
            continue;
        }

        // Skip the DISH/DISC header data line
        if header_lines > 0 {
            header_lines -= 1;
            continue;
        }

        let tokens: Vec<&str> = trimmed.split_whitespace().collect();

        match block {
            Block::Dist => {
                // Format: i j k l type  i j k l type  r0 w0 rah
                //          0 1 2 3  4   5 6 7 8  9    10 11 12
                if tokens.len() < 13 {
                    log::warn!("distanceres: short DISTANCERESSPEC line, skipping: {trimmed}");
                    continue;
                }
                let type1: i32 = tokens[4].parse().unwrap_or(0);
                let type2: i32 = tokens[9].parse().unwrap_or(0);
                if type1 != 0 || type2 != 0 {
                    log::warn!(
                        "distanceres: virtual atom type {}/{} not supported (only type 0), skipping",
                        type1, type2
                    );
                    continue;
                }
                let atom1_1: usize = tokens[0].parse::<usize>().map_err(|_| {
                    IoError::ParseError(format!("Invalid atom index: {}", tokens[0]))
                })?;
                let atom1 = atom1_1 - 1;
                let atom2_1: usize = tokens[5].parse::<usize>().map_err(|_| {
                    IoError::ParseError(format!("Invalid atom index: {}", tokens[5]))
                })?;
                let atom2 = atom2_1 - 1;
                let r0: f64 = tokens[10]
                    .parse::<f64>()
                    .map_err(|_| IoError::ParseError(format!("Invalid r0: {}", tokens[10])))?;
                let w0: f64 = tokens[11]
                    .parse::<f64>()
                    .map_err(|_| IoError::ParseError(format!("Invalid w0: {}", tokens[11])))?;
                let rah: i32 = tokens[12]
                    .parse::<i32>()
                    .map_err(|_| IoError::ParseError(format!("Invalid rah: {}", tokens[12])))?;
                unperturbed.push(DistanceRestraintSpec {
                    atom1,
                    atom2,
                    r0,
                    w0,
                    rah,
                });
            },
            Block::Pert => {
                // Format: i j k l type  i j k l type  n m  A_r0 A_w0 B_r0 B_w0 rah
                //          0 1 2 3  4   5 6 7 8  9   10 11  12   13   14   15   16
                if tokens.len() < 17 {
                    log::warn!("distanceres: short PERTDISRESSPEC line, skipping: {trimmed}");
                    continue;
                }
                let type1: i32 = tokens[4].parse().unwrap_or(0);
                let type2: i32 = tokens[9].parse().unwrap_or(0);
                if type1 != 0 || type2 != 0 {
                    return Err(IoError::ParseError(format!(
                        "PERTDISRESSPEC: virtual atom types {type1}/{type2} are not supported \
                         (only type 0, an explicit atom): {trimmed}"
                    )));
                }
                let atom1_1: usize = tokens[0].parse::<usize>().map_err(|_| {
                    IoError::ParseError(format!("Invalid atom index: {}", tokens[0]))
                })?;
                let atom1 = atom1_1 - 1;
                let atom2_1: usize = tokens[5].parse::<usize>().map_err(|_| {
                    IoError::ParseError(format!("Invalid atom index: {}", tokens[5]))
                })?;
                let atom2 = atom2_1 - 1;
                let n: i32 = tokens[10]
                    .parse::<i32>()
                    .map_err(|_| IoError::ParseError(format!("Invalid n: {}", tokens[10])))?;
                let m: i32 = tokens[11]
                    .parse::<i32>()
                    .map_err(|_| IoError::ParseError(format!("Invalid m: {}", tokens[11])))?;
                let a_r0: f64 = tokens[12]
                    .parse::<f64>()
                    .map_err(|_| IoError::ParseError(format!("Invalid A_r0: {}", tokens[12])))?;
                let a_w0: f64 = tokens[13]
                    .parse::<f64>()
                    .map_err(|_| IoError::ParseError(format!("Invalid A_w0: {}", tokens[13])))?;
                let b_r0: f64 = tokens[14]
                    .parse::<f64>()
                    .map_err(|_| IoError::ParseError(format!("Invalid B_r0: {}", tokens[14])))?;
                let b_w0: f64 = tokens[15]
                    .parse::<f64>()
                    .map_err(|_| IoError::ParseError(format!("Invalid B_w0: {}", tokens[15])))?;
                let rah: i32 = tokens[16]
                    .parse::<i32>()
                    .map_err(|_| IoError::ParseError(format!("Invalid rah: {}", tokens[16])))?;
                perturbed.push(PerturbedDistanceRestraintSpec {
                    atom1,
                    atom2,
                    n,
                    m,
                    a_r0,
                    b_r0,
                    a_w0,
                    b_w0,
                    rah,
                });
            },
            Block::None => {},
        }
    }

    Ok((unperturbed, perturbed))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_the_nacl_restraint_of_the_reference_suite() {
        let p = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../gromos-md/tests/gromosXX_references/nacl_1water_distres/nacl_1water_distres.distrest"
        );
        let (unpert, pert) = read_distanceres(p).unwrap();
        assert_eq!(unpert.len(), 1);
        assert!(pert.is_empty());
        let r = &unpert[0];
        assert!((r.r0 - 0.3).abs() < 1e-12);
        assert!((r.w0 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_missing_file_is_reported_as_such() {
        assert!(matches!(
            read_distanceres("/nonexistent/x.distrest"),
            Err(IoError::FileNotFound(_))
        ));
    }

    #[test]
    fn a_perturbed_restraint_on_a_virtual_atom_is_an_error() {
        let text = "PERTDISRESSPEC\n 0.1 0.153\n    5    0    0    0    5         9    0    0    0    0   \
                    1 1 0.19   2.0    0.21   1.0 1\nEND\n";
        let path =
            std::env::temp_dir().join(format!("gromos_pertdisres_{}.tmp", std::process::id()));
        std::fs::write(&path, text).unwrap();
        let err = read_distanceres(&path).unwrap_err().to_string();
        assert!(err.contains("virtual atom types 5/0"), "{err}");
        std::fs::remove_file(path).ok();
    }
}
