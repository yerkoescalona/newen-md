//! GROMOS-RS I/O: File format readers and writers
//!
//! This crate handles all file I/O for GROMOS-RS:
//! - .topo/.top - Topology files
//! - .conf/.cnf - Coordinate files  
//! - .imd - Input parameter files (simulation settings)
//! - .trc/.trj - Trajectory files (coordinates over time)
//! - .tre - Energy files (ENE/ENA format)
//! - .trf - Force files
//! - .ptp - Perturbation topology files (FEP)

pub mod args;
pub mod coordinate;
pub mod distanceres;
pub mod dlg;
pub mod energy;
pub mod energy_binary;
pub mod energy_traj;
pub mod force;
pub mod free_energy;
pub mod g96;
pub mod gz;
pub mod ifp;
pub mod imd;
pub mod imd_write;
pub mod input;
pub mod jobs;
pub mod jvalue;
pub mod mk_script;
pub mod mtb;
pub mod output;
pub mod pbc;
pub mod pdb;
pub mod posres;
pub mod ptp;
pub mod script_template;
pub mod table;
pub mod topology;
pub mod trajectory;
pub mod trajectory_binary;

// Re-export commonly used types
pub use args::{run_main, Arguments, Program, Status};
pub use coordinate::{
    parse_g96_labeled, read_coordinates, read_g96_labeled, CoordinateData, G96Atom,
    LabeledCoordinateData,
};
pub use dlg::{DlgWriter, LambdaDerivativeFrame};
pub use energy::read_energy_trajectory_native;
pub use energy::{EnergyBlock, EnergyFrame, EnergyReader, EnergyWriter, VolumePressure};
pub use energy_binary::{BinaryEnergyReader, BinaryEnergyWriter};
pub use force::ForceWriter;
pub use free_energy::{
    read_free_energy_trajectory, FreeEnergyFrame, FreeEnergyTrajectory, FreeEnergyWriter,
};
pub use g96::{write_por, write_rpr};
pub use imd::{
    parse_imd_str, EdsBlock, EdsStatePair, ImdParameters, LambdasBlock, LambdasRow, NonbondedExtra,
    PrecalcLamBlock, PressureParameters, TempBathParameters,
};
pub use imd_write::{write_imd, write_imd_file};
pub use input::{GamdBlock, ReplicaBlock};
pub use output::{GamdBoostWriter, GamdStatsWriter};
pub use pdb::write_pdb_positions;
pub use ptp::{read_pttopo, PerturbedTopology, PtpWriter};
pub use trajectory::TrajectoryWriter;
pub use trajectory_binary::{
    BinaryFrame, BinaryTrajectoryReader, BinaryTrajectoryWriter, DcdReader, DcdWriter,
};

use std::io;

/// Common error type for I/O operations
#[derive(Debug)]
pub enum IoError {
    FileNotFound(String),
    ParseError(String),
    FormatError(String),
    WriteError(String),
    Io(io::Error),
}

impl From<io::Error> for IoError {
    fn from(err: io::Error) -> Self {
        IoError::Io(err)
    }
}

impl std::fmt::Display for IoError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            IoError::FileNotFound(path) => write!(f, "File not found: {}", path),
            IoError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            IoError::FormatError(msg) => write!(f, "Format error: {}", msg),
            IoError::WriteError(msg) => write!(f, "Write error: {}", msg),
            IoError::Io(err) => write!(f, "I/O error: {}", err),
        }
    }
}

impl std::error::Error for IoError {}

/// Pre-process argv to support GROMOS `@key` conventions.
///
/// Converts `@key` → `--key` and expands `@f argfile` contents inline.
/// Blank lines and `#` comments in argfiles are skipped.
/// This allows using clap with the standard GROMOS argument syntax.
///
/// The process-free half of [`gromos_args`]: the caller supplies argv and gets an error back
/// instead of an exit, so an in-process caller (a `Program`'s `run`, the Python binding) and the
/// command line expand `@f` identically.
pub fn normalize(raw: Vec<String>) -> Result<Vec<String>, String> {
    let mut expanded = vec![raw.first().cloned().unwrap_or_default()];

    let mut i = 1;
    while i < raw.len() {
        if raw[i] == "@f" {
            i += 1;
            if i >= raw.len() {
                return Err("Error: @f requires a filename".to_string());
            }
            let content = std::fs::read_to_string(&raw[i])
                .map_err(|e| format!("Error: cannot read argfile '{}': {}", raw[i], e))?;
            for line in content.lines() {
                // gromos++ erases everything from the first '#' on (Arguments.cc:170)
                let line = line.split('#').next().unwrap_or("").trim();
                if line.is_empty() {
                    continue;
                }
                for token in line.split_whitespace() {
                    if let Some(key) = token.strip_prefix('@') {
                        expanded.push(format!("--{}", key));
                    } else {
                        expanded.push(token.to_string());
                    }
                }
            }
        } else if let Some(key) = raw[i].strip_prefix('@') {
            expanded.push(format!("--{}", key));
        } else {
            expanded.push(raw[i].clone());
        }
        i += 1;
    }
    Ok(expanded)
}

/// [`normalize`] over the process arguments; a bad `@f` prints and exits 1, as gromos++ does.
pub fn gromos_args() -> Vec<String> {
    normalize(std::env::args().collect()).unwrap_or_else(|e| args::fail(&e))
}

#[cfg(test)]
mod normalize_tests {
    use super::normalize;

    fn argv(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn at_keys_become_long_flags_and_values_pass_through() {
        let out = normalize(argv(&[
            "prog", "@topo", "a.top", "@pbc", "r", "--traj", "a.trc",
        ]));
        assert_eq!(
            out.unwrap(),
            argv(&["prog", "--topo", "a.top", "--pbc", "r", "--traj", "a.trc"])
        );
    }

    #[test]
    fn argfile_is_expanded_with_comments_stripped() {
        let dir = std::env::temp_dir().join("gromos_normalize_argfile");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("args.arg");
        std::fs::write(
            &path,
            "@topo a.top   # the topology\n\n# only a comment\n@pbc r\n",
        )
        .unwrap();

        let out = normalize(argv(&[
            "prog",
            "@f",
            path.to_str().unwrap(),
            "@traj",
            "a.trc",
        ]));
        assert_eq!(
            out.unwrap(),
            argv(&["prog", "--topo", "a.top", "--pbc", "r", "--traj", "a.trc"])
        );
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn a_missing_argfile_is_an_error_not_an_exit() {
        let err = normalize(argv(&["prog", "@f", "/nonexistent/args.arg"])).unwrap_err();
        assert!(err.contains("cannot read argfile"), "{err}");

        let err = normalize(argv(&["prog", "@f"])).unwrap_err();
        assert!(err.contains("@f requires a filename"), "{err}");
    }

    #[test]
    fn an_empty_argv_does_not_panic() {
        assert_eq!(normalize(Vec::new()).unwrap(), vec![String::new()]);
    }
}
