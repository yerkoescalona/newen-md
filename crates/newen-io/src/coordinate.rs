//! Parser for GROMOS coordinate files (.conf/.cnf/.g96)
//!
//! Supports both legacy BOX and GROMOS GENBOX block formats,
//! as well as POSITION/POSITIONRED and VELOCITY/VELOCITYRED blocks.
//!
//! # Format Example
//! ```text
//! TITLE
//!   System description
//! END
//! POSITION
//!     1 RES    ATOM      1    x        y        z
//! END
//! VELOCITY (optional)
//!     1 RES    ATOM      1    vx       vy       vz
//! END
//! GENBOX
//! # box type (0=vacuum, 1=rectangular)
//! 1
//! # box dimensions
//! 1.8652 1.8652 1.8652
//! 90.0 90.0 90.0
//! 0.0 0.0 0.0
//! 0.0 0.0 0.0
//! END
//! ```

use crate::IoError;
use newen_core::configuration::{Box as SimBox, Configuration};
use newen_core::math::Vec3;
use std::io::BufRead;
use std::path::Path;

/// Raw coordinate data read from a .conf/.cnf/.g96 file.
///
/// This is a lightweight struct without Configuration dependencies,
/// suitable for use from the CLI binary.
#[derive(Debug, Clone)]
pub struct CoordinateData {
    pub positions: Vec<Vec3>,
    pub velocities: Vec<Vec3>,
    /// Box dimensions (0,0,0 for vacuum)
    pub box_dims: Vec3,
    /// Box type from GENBOX (0=vacuum, 1=rectangular, 2=triclinic, 3=truncated octahedron)
    pub box_type: i32,
    /// λ of a `PERTDATA` block (gromosXX writes it into the final configuration of a perturbed
    /// run; a run with NRDGL=1 starts from it).
    pub pertdata: Option<f64>,
    /// The data lines of an `AEDSSEARCH` block (the A-EDS search state of the run that wrote
    /// the file); parsed by [`parse_aedssearch`] when a run continues the search.
    pub aedssearch: Option<Vec<String>>,
    /// The data lines of an `NHCVARIABLES` block (the Nosé–Hoover variables of the run that
    /// wrote the file, one line per bath); parsed by [`parse_nhcvariables`] when a run reads them.
    pub nhcvariables: Option<Vec<String>>,
}

/// Read raw coordinate data from a GROMOS coordinate file (.conf/.cnf/.g96).
///
/// Handles POSITION, POSITIONRED, VELOCITY, VELOCITYRED, BOX, GENBOX, PERTDATA, AEDSSEARCH and
/// NHCVARIABLES blocks.
/// Returns raw data without constructing a Configuration object.
pub fn read_coordinates<P: AsRef<Path>>(path: P) -> Result<CoordinateData, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    parse_coordinates(reader)
}

/// [`read_coordinates`] over an already-open reader — for callers holding the text, not a path.
pub fn parse_coordinates(reader: impl BufRead) -> Result<CoordinateData, IoError> {
    let mut positions = Vec::new();
    let mut velocities = Vec::new();
    let mut box_dims = Vec3::ZERO;
    let mut box_type: i32 = -1; // unknown until parsed
    let mut pertdata: Option<f64> = None;
    let mut aedssearch: Option<Vec<String>> = None;
    let mut nhcvariables: Option<Vec<String>> = None;

    #[derive(PartialEq)]
    enum Section {
        None,
        Position,
        PositionRed,
        Velocity,
        VelocityRed,
        Box,
        GenBox,
        PertData,
        AedsSearch,
        NhcVariables,
    }
    let mut section = Section::None;
    let mut genbox_line_num: usize = 0;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        match trimmed {
            "POSITION" => {
                section = Section::Position;
                continue;
            },
            "POSITIONRED" => {
                section = Section::PositionRed;
                continue;
            },
            "VELOCITY" => {
                section = Section::Velocity;
                continue;
            },
            "VELOCITYRED" => {
                section = Section::VelocityRed;
                continue;
            },
            "BOX" => {
                section = Section::Box;
                continue;
            },
            "GENBOX" => {
                section = Section::GenBox;
                genbox_line_num = 0;
                continue;
            },
            "PERTDATA" => {
                section = Section::PertData;
                continue;
            },
            "AEDSSEARCH" => {
                section = Section::AedsSearch;
                aedssearch = Some(Vec::new());
                continue;
            },
            "NHCVARIABLES" => {
                section = Section::NhcVariables;
                nhcvariables = Some(Vec::new());
                continue;
            },
            "END" => {
                section = Section::None;
                continue;
            },
            _ => {},
        }

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match section {
            Section::Position => {
                // Full format: first 24 chars are residue/atom metadata, then x y z
                if line.len() < 24 {
                    continue;
                }
                let coords = &line[24..];
                let parts: Vec<&str> = coords.split_whitespace().collect();
                if parts.len() >= 3 {
                    positions.push(parse_vec3(&parts[..3])?);
                }
            },
            Section::PositionRed => {
                // Reduced format: just x y z (last 3 whitespace-separated values)
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 3 {
                    let len = parts.len();
                    positions.push(parse_vec3(&parts[len - 3..len])?);
                }
            },
            Section::Velocity => {
                // `G96Writer::write_velocities` (the `md` binary's final configuration) writes the
                // three components without the 24-character label
                let all: Vec<&str> = trimmed.split_whitespace().collect();
                if all.len() == 3 {
                    velocities.push(parse_vec3(&all)?);
                    continue;
                }
                if line.len() < 24 {
                    continue;
                }
                let coords = &line[24..];
                let parts: Vec<&str> = coords.split_whitespace().collect();
                if parts.len() >= 3 {
                    velocities.push(parse_vec3(&parts[..3])?);
                }
            },
            Section::VelocityRed => {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 3 {
                    let len = parts.len();
                    velocities.push(parse_vec3(&parts[len - 3..len])?);
                }
            },
            Section::Box => {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 3 {
                    box_dims = parse_vec3(&parts[..3])?;
                    if box_type < 0 {
                        box_type = if box_dims == Vec3::ZERO { 0 } else { 1 };
                    }
                }
            },
            Section::GenBox => {
                // GENBOX format:
                //   line 0: box_type (int)
                //   line 1: box dimensions (3 floats)
                //   line 2: angles (3 floats) - ignored for rectangular
                //   line 3: origin (3 floats) - ignored
                //   line 4: origin2 (3 floats) - ignored
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                match genbox_line_num {
                    0 => {
                        box_type = parts[0].parse::<i32>().map_err(|_| {
                            IoError::ParseError(format!("Invalid GENBOX type: {}", parts[0]))
                        })?;
                    },
                    1 if parts.len() >= 3 => {
                        box_dims = parse_vec3(&parts[..3])?;
                    },
                    _ => {}, // angles, origin lines - skip for now
                }
                genbox_line_num += 1;
            },
            // gromosXX `_read_pertdata` (in_configuration.cc:2351-2371): one value, λ. Only a
            // run with NRDGL=1 reads it, so an unreadable value is left as `None` here.
            Section::PertData => {
                if pertdata.is_none() {
                    pertdata = trimmed
                        .split_whitespace()
                        .next()
                        .and_then(|x| x.parse::<f64>().ok());
                }
            },
            Section::AedsSearch => {
                if let Some(lines) = aedssearch.as_mut() {
                    lines.push(trimmed.to_string());
                }
            },
            Section::NhcVariables => {
                if let Some(lines) = nhcvariables.as_mut() {
                    lines.push(trimmed.to_string());
                }
            },
            Section::None => {},
        }
    }

    if positions.is_empty() {
        return Err(IoError::FormatError(
            "No positions found in coordinate file".to_string(),
        ));
    }

    Ok(CoordinateData {
        positions,
        velocities,
        box_dims,
        box_type,
        pertdata,
        aedssearch,
        nhcvariables,
    })
}

/// The Nosé–Hoover variables of an `NHCVARIABLES` block's data lines: one row per bath, every
/// value of the line (gromosXX `_print_nose_hoover_chain_variables`,
/// `out_configuration.cc:3426-3443`, writes the whole chain of each bath on one line).
pub fn parse_nhcvariables(lines: &[String]) -> Result<Vec<Vec<f64>>, IoError> {
    lines
        .iter()
        .map(|l| {
            l.split_whitespace()
                .map(|x| {
                    x.parse::<f64>().map_err(|_| {
                        IoError::ParseError(format!("bad value '{x}' in NHCVARIABLES block"))
                    })
                })
                .collect()
        })
        .collect()
}

/// The A-EDS search state of an `AEDSSEARCH` block's data lines for `num_states` states
/// (gromosXX `_read_aedssearch`, `in_configuration.cc:3687-3752`): E_max, E_min, the E_max
/// candidate, the round-trip count, the previous state and the RESTREMIN flag one per line, then
/// one line per state with E_i^R, the log average, the free energy, the visited flag, the visit
/// count, the two mean energies, the Welford sum, the standard deviation and the frame count.
/// Unlike gromosXX (which checks the stream once, after the last line) every line and the
/// number of state lines are checked.
pub fn parse_aedssearch(
    lines: &[String],
    num_states: usize,
) -> Result<newen_core::configuration::AedsSearch, IoError> {
    let bad = |what: &str| IoError::ParseError(format!("bad line in AEDSSEARCH block: {what}"));
    if lines.len() != 6 + num_states {
        return Err(IoError::ParseError(format!(
            "AEDSSEARCH block has {} lines, NUMSTATES = {num_states} needs {}",
            lines.len(),
            6 + num_states
        )));
    }
    fn first<T: std::str::FromStr>(line: &str) -> Option<T> {
        line.split_whitespace().next()?.parse().ok()
    }
    let mut s = newen_core::configuration::AedsSearch::new(
        first(&lines[0]).ok_or_else(|| bad("EMAX"))?,
        first(&lines[1]).ok_or_else(|| bad("EMIN"))?,
        vec![0.0; num_states],
        first::<i64>(&lines[5]).ok_or_else(|| bad("FULLEMIN"))? == 1,
    );
    s.searchemax = first(&lines[2]).ok_or_else(|| bad("SEARCHEMAX"))?;
    s.emaxcounts = first(&lines[3]).ok_or_else(|| bad("EMAXCOUNTS"))?;
    s.oldstate = first(&lines[4]).ok_or_else(|| bad("OLDSTATE"))?;
    s.initialise = false;
    for (i, line) in lines[6..].iter().enumerate() {
        let v: Vec<&str> = line.split_whitespace().collect();
        let what = format!("state {}", i + 1);
        if v.len() != 10 {
            return Err(bad(&what));
        }
        let real = |k: usize| v[k].parse::<f64>().map_err(|_| bad(&what));
        let count = |k: usize| v[k].parse::<u32>().map_err(|_| bad(&what));
        s.eir[i] = real(0)?;
        s.lnexpde[i] = real(1)?;
        s.statefren[i] = real(2)?;
        s.visited[i] = v[3].parse::<i64>().map_err(|_| bad(&what))? == 1;
        s.visitcounts[i] = count(4)?;
        s.avgenergy[i] = real(5)?;
        s.eiravgenergy[i] = real(6)?;
        s.bigs[i] = real(7)?;
        s.stdevenergy[i] = real(8)?;
        s.framecounts[i] = count(9)?;
    }
    Ok(s)
}

/// An atom with its label information and coordinates, as read from a G96 file.
#[derive(Debug, Clone)]
pub struct G96Atom {
    pub res_num: usize,
    pub res_name: String,
    pub atom_name: String,
    pub atom_num: usize,
    pub pos: Vec3,
}

/// Coordinate data with full atom labels, as read from a G96 file.
#[derive(Debug, Clone)]
pub struct LabeledCoordinateData {
    pub atoms: Vec<G96Atom>,
    pub box_dims: Vec3,
    pub box_type: i32,
}

/// Read a GROMOS coordinate file, preserving atom labels (residue number, name, atom name).
///
/// This is the labeled counterpart to [`read_coordinates`]. Use this when downstream
/// code needs to know which atom is which (e.g. sim_box, ion, pdb2g96).
/// A labelled coordinate file in gromos++ `OutG96S` layout: TITLE, POSITION (residue number,
/// residue name, atom name, serial, x y z) and GENBOX (rectangular when `box_dims` is given,
/// vacuum otherwise). Serial numbers are renumbered from 1.
pub fn format_g96(title: &str, atoms: &[G96Atom], box_dims: Option<Vec3>) -> String {
    let mut o = String::new();
    o.push_str("TITLE\n");
    o.push_str(title);
    if !title.ends_with('\n') {
        o.push('\n');
    }
    o.push_str("END\nPOSITION\n");
    for (k, a) in atoms.iter().enumerate() {
        o.push_str(&format!(
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}\n",
            a.res_num,
            a.res_name,
            a.atom_name,
            k + 1,
            a.pos.x,
            a.pos.y,
            a.pos.z
        ));
    }
    o.push_str("END\nGENBOX\n");
    let b = box_dims.unwrap_or(Vec3::ZERO);
    o.push_str(&format!("{:>8}\n", if box_dims.is_some() { 1 } else { 0 }));
    o.push_str(&format!("{:15.9}{:15.9}{:15.9}\n", b.x, b.y, b.z));
    let angles = if box_dims.is_some() { 90.0 } else { 0.0 };
    o.push_str(&format!("{:15.9}{:15.9}{:15.9}\n", angles, angles, angles));
    o.push_str(&format!("{:15.9}{:15.9}{:15.9}\n", 0.0, 0.0, 0.0));
    o.push_str(&format!("{:15.9}{:15.9}{:15.9}\nEND\n", 0.0, 0.0, 0.0));
    o
}

/// A VELOCITY block in the layout of [`format_g96`]'s POSITION block, labelled by `atoms`
/// (velocities in nm/ps). Readers find blocks by name, so it may follow the GENBOX block.
pub fn format_g96_velocity_block(atoms: &[G96Atom], velocities: &[Vec3]) -> String {
    let mut o = String::from("VELOCITY\n");
    for (k, (a, v)) in atoms.iter().zip(velocities).enumerate() {
        o.push_str(&format!(
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}\n",
            a.res_num,
            a.res_name,
            a.atom_name,
            k + 1,
            v.x,
            v.y,
            v.z
        ));
    }
    o.push_str("END\n");
    o
}

pub fn read_g96_labeled<P: AsRef<Path>>(path: P) -> Result<LabeledCoordinateData, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    parse_g96_labeled(reader)
}

/// [`read_g96_labeled`] over an already-open reader — for callers holding the text, not a path.
pub fn parse_g96_labeled(reader: impl BufRead) -> Result<LabeledCoordinateData, IoError> {
    let mut atoms = Vec::new();
    let mut box_dims = Vec3::ZERO;
    let mut box_type: i32 = -1;

    #[derive(PartialEq)]
    enum Section {
        None,
        Position,
        Box,
        GenBox,
        Skip,
    }
    let mut section = Section::None;
    let mut genbox_line_num: usize = 0;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        match trimmed {
            "POSITION" | "POSITIONRED" => {
                section = Section::Position;
                continue;
            },
            "BOX" => {
                section = Section::Box;
                continue;
            },
            "GENBOX" => {
                section = Section::GenBox;
                genbox_line_num = 0;
                continue;
            },
            "END" => {
                section = Section::None;
                continue;
            },
            "TITLE" | "VELOCITY" | "VELOCITYRED" | "REFPOSITION" | "POSRESSPEC"
            | "LATTICESHIFTS" | "STOCHINT" => {
                section = Section::Skip;
                continue;
            },
            _ => {},
        }

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match section {
            Section::Position => {
                // Parse using whitespace splitting (robust for varying column widths)
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 7 {
                    let n = parts.len();
                    let res_num: usize = parts[0].parse().unwrap_or(1);
                    let res_name = parts[1].to_string();
                    let atom_name = parts[2].to_string();
                    let atom_num: usize = parts[3].parse().unwrap_or(0);
                    let x: f64 = parts[n - 3]
                        .parse()
                        .map_err(|_| IoError::ParseError(format!("Invalid x: {}", parts[n - 3])))?;
                    let y: f64 = parts[n - 2]
                        .parse()
                        .map_err(|_| IoError::ParseError(format!("Invalid y: {}", parts[n - 2])))?;
                    let z: f64 = parts[n - 1]
                        .parse()
                        .map_err(|_| IoError::ParseError(format!("Invalid z: {}", parts[n - 1])))?;
                    atoms.push(G96Atom {
                        res_num,
                        res_name,
                        atom_name,
                        atom_num,
                        pos: Vec3::new(x, y, z),
                    });
                }
            },
            Section::Box => {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 3 {
                    box_dims = parse_vec3(&parts[..3])?;
                    if box_type < 0 {
                        box_type = if box_dims == Vec3::ZERO { 0 } else { 1 };
                    }
                }
            },
            Section::GenBox => {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                match genbox_line_num {
                    0 => {
                        box_type = parts[0].parse::<i32>().map_err(|_| {
                            IoError::ParseError(format!("Invalid GENBOX type: {}", parts[0]))
                        })?;
                    },
                    1 if parts.len() >= 3 => {
                        box_dims = parse_vec3(&parts[..3])?;
                    },
                    _ => {},
                }
                genbox_line_num += 1;
            },
            _ => {},
        }
    }

    if atoms.is_empty() {
        return Err(IoError::FormatError(
            "No positions found in coordinate file".to_string(),
        ));
    }

    Ok(LabeledCoordinateData {
        atoms,
        box_dims,
        box_type,
    })
}

fn parse_vec3(parts: &[&str]) -> Result<Vec3, IoError> {
    let x: f64 = parts[0]
        .parse()
        .map_err(|_| IoError::ParseError(format!("Invalid coordinate: {}", parts[0])))?;
    let y: f64 = parts[1]
        .parse()
        .map_err(|_| IoError::ParseError(format!("Invalid coordinate: {}", parts[1])))?;
    let z: f64 = parts[2]
        .parse()
        .map_err(|_| IoError::ParseError(format!("Invalid coordinate: {}", parts[2])))?;
    Ok(Vec3::new(x, y, z))
}

/// Read GROMOS coordinate file and return a Configuration object.
///
/// This is a convenience wrapper around [`read_coordinates`] that constructs
/// a full Configuration with the specified number of temperature/energy groups.
pub fn read_coordinate_file<P: AsRef<Path>>(
    path: P,
    num_temp_groups: usize,
    num_energy_groups: usize,
) -> Result<Configuration, IoError> {
    let data = read_coordinates(path)?;
    let n_atoms = data.positions.len();

    let mut conf = Configuration::new(n_atoms, num_temp_groups, num_energy_groups);
    conf.current_mut().pos = data.positions;

    if !data.velocities.is_empty() {
        if data.velocities.len() != n_atoms {
            return Err(IoError::FormatError(format!(
                "Velocity count ({}) doesn't match atom count ({})",
                data.velocities.len(),
                n_atoms
            )));
        }
        conf.current_mut().vel = data.velocities;
    }

    if data.box_dims != Vec3::ZERO {
        conf.current_mut().box_config =
            SimBox::rectangular(data.box_dims.x, data.box_dims.y, data.box_dims.z);
    }

    Ok(conf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Minimal 4-atom GROMOS .conf content (cg16-style).
    /// Position lines: first 24 chars are residue/atom metadata, then x y z.
    const CG16_CONF: &str = "\
TITLE
  4-atom test system (cg16)
END
POSITION
    1 RES    C         1   4.197156491  0.505921049  2.679733124
    1 RES    H         2   4.300000000  0.600000000  2.700000000
    1 RES    H         3   4.100000000  0.400000000  2.650000000
    1 RES    H         4   4.200000000  0.500000000  2.800000000
END
BOX
  10.0  10.0  10.0
END
";

    fn write_tmp(content: &str, suffix: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("newen_test_{suffix}.tmp"));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_parse_cg16() {
        let path = write_tmp(CG16_CONF, "cg16_conf");
        let conf = read_coordinate_file(&path, 1, 1).expect("Failed to parse inline conf");

        assert_eq!(conf.current().pos.len(), 4);

        let pos0 = conf.current().pos[0];
        assert!((pos0.x - 4.197156491).abs() < 1e-5);
        assert!((pos0.y - 0.505921049).abs() < 1e-5);
        assert!((pos0.z - 2.679733124).abs() < 1e-5);

        let box_dims = conf.current().box_config.dimensions();
        assert!((box_dims.x - 10.0).abs() < 1e-5);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_parse_missing_file_returns_error() {
        let result = read_coordinate_file("/nonexistent/path/file.conf", 1, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_conf_with_velocities() {
        let content = "\
TITLE
  test
END
POSITION
    1 RES    C         1   1.000000000  2.000000000  3.000000000
END
VELOCITY
    1 RES    C         1   0.010000000  0.020000000  0.030000000
END
BOX
  5.0  5.0  5.0
END
";
        let path = write_tmp(content, "conf_vel");
        let conf = read_coordinate_file(&path, 1, 1).expect("Failed to parse conf with velocities");

        assert_eq!(conf.current().pos.len(), 1);
        assert_eq!(conf.current().vel.len(), 1);
        assert!((conf.current().vel[0].x - 0.01).abs() < 1e-5);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_read_coordinates_genbox_rectangular() {
        let content = "\
TITLE
  GENBOX test
END
POSITION
    1 WAT    OW        1   0.100000000  0.200000000  0.300000000
    1 WAT    HW1       2   0.150000000  0.250000000  0.350000000
END
VELOCITY
    1 WAT    OW        1   0.001000000  0.002000000  0.003000000
    1 WAT    HW1       2   0.004000000  0.005000000  0.006000000
END
GENBOX
# box type (1=rectangular)
1
# box dimensions
1.8652 1.8652 1.8652
90.0 90.0 90.0
0.0 0.0 0.0
0.0 0.0 0.0
END
";
        let path = write_tmp(content, "genbox_rect");
        let data = read_coordinates(&path).expect("Failed to parse GENBOX file");

        assert_eq!(data.positions.len(), 2);
        assert_eq!(data.velocities.len(), 2);
        assert_eq!(data.box_type, 1);
        assert!((data.box_dims.x - 1.8652).abs() < 1e-9);
        assert!((data.box_dims.y - 1.8652).abs() < 1e-9);
        assert!((data.positions[0].x - 0.1).abs() < 1e-9);
        assert!((data.velocities[1].z - 0.006).abs() < 1e-9);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn pertdata_is_read_and_absent_is_none() {
        let with = "\
TITLE
  pertdata
END
POSITION
    1 AR     AR        1   0.000000000  0.000000000  0.000000000
END
PERTDATA
      0.2500000
END
";
        let path = write_tmp(with, "pertdata_with");
        assert_eq!(read_coordinates(&path).unwrap().pertdata, Some(0.25));
        std::fs::remove_file(path).ok();
        let path = write_tmp(CG16_CONF, "pertdata_without");
        assert_eq!(read_coordinates(&path).unwrap().pertdata, None);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_read_coordinates_genbox_vacuum() {
        let content = "\
TITLE
  vacuum test
END
POSITION
    1 AR     AR        1   0.000000000  0.000000000  0.000000000
    1 AR     AR        2   0.350000000  0.000000000  0.000000000
END
GENBOX
# box type (0=vacuum)
0
# box dimensions
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
END
";
        let path = write_tmp(content, "genbox_vac");
        let data = read_coordinates(&path).expect("Failed to parse vacuum GENBOX");

        assert_eq!(data.positions.len(), 2);
        assert_eq!(data.box_type, 0);
        assert_eq!(data.box_dims, Vec3::ZERO);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_read_coordinates_positionred() {
        let content = "\
TITLE
  reduced format test
END
POSITIONRED
   0.100000000   0.200000000   0.300000000
   0.400000000   0.500000000   0.600000000
END
VELOCITYRED
   0.001000000   0.002000000   0.003000000
   0.004000000   0.005000000   0.006000000
END
GENBOX
0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
END
";
        let path = write_tmp(content, "posred");
        let data = read_coordinates(&path).expect("Failed to parse POSITIONRED");

        assert_eq!(data.positions.len(), 2);
        assert_eq!(data.velocities.len(), 2);
        assert!((data.positions[0].x - 0.1).abs() < 1e-9);
        assert!((data.velocities[1].y - 0.005).abs() < 1e-9);

        std::fs::remove_file(path).ok();
    }
}
