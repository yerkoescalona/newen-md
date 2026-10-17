//! IMD file text-level manipulation for mk_script.
//!
//! mk_script reads a base .imd file and overrides specific parameter values
//! per job. Parameters are identified by name (matching the GROMOS convention)
//! and located within their IMD block.
//!
//! The override mechanism works at the text level: each IMD block has a known
//! layout mapping parameter names to (data_line_index, field_index) positions.
//! This avoids needing a full IMD serializer while supporting the flexible
//! .jobs column scheme.

use crate::jobs::JobSpec;
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;

/// Mapping of IMD parameter names to their location in the IMD text.
///
/// Each entry: (block_name, data_line_index (1-based), field_index (0-based))
/// This covers the most commonly overridden parameters. Parameters that use
/// array notation like TEMP0[1] are handled specially.
fn param_location(name: &str) -> Option<(&'static str, usize, usize)> {
    // (block, data_line, field_index)
    Some(match name {
        // INITIALISE block: 1 data line with 10 fields
        "NTIVEL" => ("INITIALISE", 1, 0),
        "NTISHK" => ("INITIALISE", 1, 1),
        "NTINHT" => ("INITIALISE", 1, 2),
        "NTINHB" => ("INITIALISE", 1, 3),
        "NTISHI" => ("INITIALISE", 1, 4),
        "NTIRTC" => ("INITIALISE", 1, 5),
        "NTICOM" => ("INITIALISE", 1, 6),
        "NTISTI" => ("INITIALISE", 1, 7),
        "IG" => ("INITIALISE", 1, 8),
        "TEMPI" => ("INITIALISE", 1, 9),

        // STEP block: 1 data line
        "NSTLIM" => ("STEP", 1, 0),
        "T" => ("STEP", 1, 1),
        "DT" => ("STEP", 1, 2),

        // BOUNDCOND block: 1 data line
        "NTB" => ("BOUNDCOND", 1, 0),
        "NDFMIN" => ("BOUNDCOND", 1, 1),

        // POSITIONRES block: 1 data line
        "NTPOR" => ("POSITIONRES", 1, 0),
        "NTPORB" => ("POSITIONRES", 1, 1),
        "NTPORS" => ("POSITIONRES", 1, 2),
        "CPOR" => ("POSITIONRES", 1, 3),

        // ENERGYMIN block: 1 data line
        "NTEM" => ("ENERGYMIN", 1, 0),
        "NCYC" => ("ENERGYMIN", 1, 1),
        "DELE" => ("ENERGYMIN", 1, 2),
        "DX0" => ("ENERGYMIN", 1, 3),
        "DXM" => ("ENERGYMIN", 1, 4),
        "NMIN" => ("ENERGYMIN", 1, 5),
        "FLIM" => ("ENERGYMIN", 1, 6),

        // COMTRANSROT block: 1 data line
        "NSCM" => ("COMTRANSROT", 1, 0),

        // CONSTRAINT block
        "NTC" => ("CONSTRAINT", 1, 0),

        // NONBONDED block: 1 data line
        "NLRELE" => ("NONBONDED", 1, 0),

        // WRITETRAJ block: 1 data line
        "NTWX" => ("WRITETRAJ", 1, 0),
        "NTWSE" => ("WRITETRAJ", 1, 1),
        "NTWV" => ("WRITETRAJ", 1, 2),
        "NTWF" => ("WRITETRAJ", 1, 3),
        "NTWE" => ("WRITETRAJ", 1, 4),
        "NTWG" => ("WRITETRAJ", 1, 5),
        "NTWB" => ("WRITETRAJ", 1, 6),

        // PRINTOUT block: 1 data line
        "NTPR" => ("PRINTOUT", 1, 0),

        // PRESSURESCALE/BAROSTAT
        "NTP" => ("PRESSURESCALE", 1, 0),
        "COUPLE" => ("PRESSURESCALE", 1, 0), // alias

        _ => return None,
    })
}

/// Apply job-specific overrides to a base IMD file (as text lines).
///
/// Parameters are matched by name to their block/field location.
/// Array parameters like TEMP0[N] and PRSBTH[N] are handled specially.
///
/// Returns the modified IMD content as a string.
pub fn apply_job_overrides(base_imd: &str, job: &JobSpec) -> String {
    // Build a per-block, per-data-line override index
    // Key: (block_name, data_line_idx) → Vec<(field_idx, new_value)>
    let mut overrides: HashMap<(String, usize), Vec<(usize, String)>> = HashMap::new();

    // Collect TEMP0[N] / TAU[N] overrides for special MULTIBATH handling
    let mut temp0_overrides: HashMap<usize, String> = HashMap::new();
    let mut tau_overrides: HashMap<usize, String> = HashMap::new();

    for (name, value) in &job.params {
        // Handle TEMP0[N] / TAU[N] array parameters (gromos++ mk_script.cc:4797)
        let indexed = |prefix: &str| {
            name.strip_prefix(prefix)
                .and_then(|rest| rest.strip_suffix(']'))
                .and_then(|idx| idx.parse::<usize>().ok())
        };
        if let Some(idx) = indexed("TEMP0[") {
            temp0_overrides.insert(idx, value.clone());
            continue;
        }
        if let Some(idx) = indexed("TAU[") {
            tau_overrides.insert(idx, value.clone());
            continue;
        }

        // Handle standard parameters via the location table
        if let Some((block, line_idx, field_idx)) = param_location(name) {
            overrides
                .entry((block.to_string(), line_idx))
                .or_default()
                .push((field_idx, value.clone()));
        }
    }

    let mut output = String::with_capacity(base_imd.len());
    let mut current_block = String::new();
    let mut data_line_idx: usize = 0;

    for line in base_imd.lines() {
        let trimmed = line.trim();

        if trimmed == "END" {
            current_block.clear();
            data_line_idx = 0;
            writeln!(output, "{}", line).unwrap();
            continue;
        }

        // Detect block start. Only outside a block: IMD blocks do not nest, and a data line
        // that is a bare integer (MULTIBATH's ALGORITHM/NBATHS/DOFSET) otherwise reads as a
        // block name and resets the data-line counter.
        if current_block.is_empty()
            && !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && trimmed
                .chars()
                .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_')
        {
            current_block = trimmed.to_string();
            data_line_idx = 0;
            writeln!(output, "{}", line).unwrap();
            continue;
        }

        if trimmed.starts_with('#') || trimmed.is_empty() {
            writeln!(output, "{}", line).unwrap();
            continue;
        }

        data_line_idx += 1;

        // Special handling for MULTIBATH TEMP0/TAU line: TEMP0(1..NBATHS) TAU(1..NBATHS)
        if current_block == "MULTIBATH"
            && data_line_idx == 3
            && !(temp0_overrides.is_empty() && tau_overrides.is_empty())
        {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            let n_baths = parts.len() / 2;
            let mut new_line = String::new();
            for i in 0..n_baths {
                let pick = |ov: &HashMap<usize, String>, field: usize| {
                    ov.get(&(i + 1))
                        .cloned()
                        .unwrap_or_else(|| parts[field].to_string())
                };
                let temp = pick(&temp0_overrides, i * 2);
                let tau = pick(&tau_overrides, i * 2 + 1);
                write!(new_line, "{:>10}{:>8}", temp, tau).unwrap();
            }
            writeln!(output, "{}", new_line).unwrap();
            continue;
        }

        // Check if this data line has overrides
        let key = (current_block.clone(), data_line_idx);
        if let Some(field_overrides) = overrides.get(&key) {
            let mut parts: Vec<String> =
                trimmed.split_whitespace().map(|s| s.to_string()).collect();
            for (field_idx, value) in field_overrides {
                if *field_idx < parts.len() {
                    parts[*field_idx] = value.clone();
                }
            }
            // Reconstruct with reasonable spacing
            let new_line: String = parts.iter().map(|p| format!("{:>12}", p)).collect();
            writeln!(output, "{}", new_line).unwrap();
            continue;
        }

        writeln!(output, "{}", line).unwrap();
    }

    output
}

/// The `WRITETRAJ` switches of a job's IMD file: which trajectories the run produces.
///
/// gromos++ decides from these both which `@tr*` flags the run script passes and which files it
/// gzips afterwards (`mk_script.cc:3853`), which is why every analysis input in the GROMOS
/// tutorials is a `.gz`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct WriteTraj {
    pub ntwx: i64,
    pub ntwse: i64,
    pub ntwv: i64,
    pub ntwf: i64,
    pub ntwe: i64,
    pub ntwg: i64,
    pub ntwb: i64,
}

impl WriteTraj {
    /// Read the switches from an IMD file's `WRITETRAJ` block; all zero when the block is absent.
    pub fn from_imd(imd: &str) -> Self {
        let mut in_block = false;
        for line in imd.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if trimmed == "END" {
                if in_block {
                    break;
                }
                continue;
            }
            if !in_block {
                in_block = trimmed == "WRITETRAJ";
                continue;
            }
            let v: Vec<i64> = trimmed
                .split_whitespace()
                .map(|t| t.parse().unwrap_or(0))
                .collect();
            let at = |i: usize| v.get(i).copied().unwrap_or(0);
            return WriteTraj {
                ntwx: at(0),
                ntwse: at(1),
                ntwv: at(2),
                ntwf: at(3),
                ntwe: at(4),
                ntwg: at(5),
                ntwb: at(6),
            };
        }
        WriteTraj::default()
    }
}

/// Configuration for mk_script generation.
#[derive(Debug, Clone)]
pub struct MkScriptConfig {
    /// System name (used in filename templates)
    pub system: String,
    /// Path to MD binary
    pub bin: String,
    /// Working directory
    pub dir: String,
    /// Input files: key → path
    pub files: HashMap<String, String>,
}

/// Generate a shell run script for a single job.
pub fn generate_run_script(
    config: &MkScriptConfig,
    template: &crate::script_template::ScriptTemplate,
    job: &JobSpec,
    prev_job: Option<&JobSpec>,
    writetraj: &WriteTraj,
) -> String {
    let num = job.job_id;
    let sys = &config.system;

    let input_name = template.expanded_filename("input", sys, num);
    let output_name = template.expanded_filename("output", sys, num);
    let fin_name = template.expanded_filename("coord", sys, num);

    // Initial coordinates: first job uses @files coord, rest use previous job's output coords
    let coord = if let Some(prev) = prev_job {
        template.expanded_filename("coord", sys, prev.job_id)
    } else {
        config
            .files
            .get("coord")
            .cloned()
            .unwrap_or_else(|| "input.cnf".to_string())
    };

    let topo = config
        .files
        .get("topo")
        .cloned()
        .unwrap_or_else(|| "system.top".to_string());

    let mut script = String::new();
    writeln!(script, "#!/bin/sh").unwrap();
    writeln!(script, "# mk_script generated run script").unwrap();
    writeln!(script, "# Job {} of system {}", num, sys).unwrap();
    writeln!(script).unwrap();

    // MD command
    writeln!(script, "{} \\", config.bin).unwrap();
    writeln!(script, "  @topo {} \\", topo).unwrap();
    writeln!(script, "  @conf {} \\", coord).unwrap();
    writeln!(script, "  @fin {} \\", fin_name).unwrap();
    writeln!(script, "  @input {} \\", input_name).unwrap();

    // One output flag per WRITETRAJ switch that is on, as gromos++ does.
    let outputs: Vec<(&str, &str, i64)> = vec![
        ("@trc", "outtrx", writetraj.ntwx),
        ("@trv", "outtrv", writetraj.ntwv),
        ("@trf", "outtrf", writetraj.ntwf),
        ("@tre", "outtre", writetraj.ntwe),
        ("@trg", "outtrg", writetraj.ntwg),
        ("@bae", "outbae", writetraj.ntwb),
    ];
    let produced: Vec<String> = outputs
        .iter()
        .filter(|(_, _, on)| *on > 0)
        .map(|(flag, kind, _)| {
            let name = template.expanded_filename(kind, sys, num);
            writeln!(script, "  {} {} \\", flag, name).unwrap();
            name
        })
        .collect();

    // Position restraint files (check if NTPOR is set and > 0)
    let has_posres = job.get_i32("NTPOR").map(|v| v > 0).unwrap_or(false);
    if has_posres {
        if let Some(por) = config.files.get("posresspec") {
            writeln!(script, "  @posresspec {} \\", por).unwrap();
        }
        if let Some(rpr) = config.files.get("refpos") {
            writeln!(script, "  @refpos {} \\", rpr).unwrap();
        }
    }

    writeln!(script, "  > {}", output_name).unwrap();

    // gromos++ compresses every trajectory it produced (`mk_script.cc:3853`)
    for name in &produced {
        writeln!(script, "gzip {}", name).unwrap();
    }

    // Chain to next job
    if let Some(lastcmd) = template.misc.get("lastcommand") {
        let next_num = num + 1;
        let next_script = crate::script_template::ScriptTemplate::expand(lastcmd, sys, next_num);
        writeln!(script).unwrap();
        writeln!(script, "# chain to next job").unwrap();
        writeln!(script, "{}", next_script).unwrap();
    }

    script
}

#[cfg(test)]
mod tests {
    use super::*;

    fn job(params: &[(&str, &str)]) -> JobSpec {
        JobSpec {
            job_id: 1,
            params: params
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect(),
            subdir: ".".to_string(),
            run_after: 0,
        }
    }

    /// The thermalisation joblist of the LiveCoMS GB3 tutorial: MULTIBATH's ALGORITHM, NBATHS
    /// and DOFSET are bare integers on their own lines, which must not read as block names.
    const MULTIBATH_IMD: &str = "\
MULTIBATH
#          ALGORITHM
                   0
#  NBATHS
         2
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
        60     0.1      60     0.1
#   DOFSET
         2
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
       567         1         1       19091        2         2
END
";

    #[test]
    fn temp0_overrides_reach_the_multibath_line() {
        let out = apply_job_overrides(
            MULTIBATH_IMD,
            &job(&[("TEMP0[1]", "300.0"), ("TEMP0[2]", "300.0")]),
        );
        let temps: Vec<&str> = out
            .lines()
            .nth(6)
            .expect("TEMP0/TAU line")
            .split_whitespace()
            .collect();
        assert_eq!(temps, ["300.0", "0.1", "300.0", "0.1"]);
        // the DOFSET line below it is untouched
        assert!(out.contains("       567         1         1       19091        2         2"));
    }

    #[test]
    fn tau_overrides_reach_the_multibath_line() {
        let out = apply_job_overrides(MULTIBATH_IMD, &job(&[("TAU[2]", "0.4")]));
        let temps: Vec<&str> = out
            .lines()
            .nth(6)
            .expect("TEMP0/TAU line")
            .split_whitespace()
            .collect();
        assert_eq!(temps, ["60", "0.1", "60", "0.4"]);
    }

    #[test]
    fn writetraj_drives_the_output_flags_and_the_gzip_lines() {
        let imd = "TITLE\nt\nEND\nWRITETRAJ\n#  NTWX NTWSE NTWV NTWF NTWE NTWG NTWB\n  100 0 50 0 100 0 0\nEND\n";
        let wt = WriteTraj::from_imd(imd);
        assert_eq!(
            (wt.ntwx, wt.ntwv, wt.ntwf, wt.ntwe),
            (100, 50, 0, 100),
            "{wt:?}"
        );

        let template = crate::script_template::ScriptTemplate {
            filenames: [
                ("script", "%system%_%number%.run"),
                ("input", "%system%_%number%.imd"),
                ("output", "%system%_%number%.omd"),
                ("coord", "%system%_%number%.cnf"),
                ("outtrx", "%system%_%number%.trc"),
                ("outtrv", "%system%_%number%.trv"),
                ("outtre", "%system%_%number%.tre"),
            ]
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
            misc: HashMap::new(),
        };
        let config = MkScriptConfig {
            system: "s".to_string(),
            bin: "md".to_string(),
            dir: ".".to_string(),
            files: [("topo", "s.top"), ("coord", "in.cnf")]
                .into_iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect(),
        };
        let script = generate_run_script(&config, &template, &job(&[]), None, &wt);
        for expected in ["@trc s_1.trc", "@trv s_1.trv", "@tre s_1.tre"] {
            assert!(script.contains(expected), "{expected} missing:\n{script}");
        }
        assert!(!script.contains("@trf"), "NTWF=0 must not write forces");
        // gromos++ compresses each trajectory it produced — which is why the tutorials'
        // analysis inputs are all .gz
        assert!(script.contains("gzip s_1.trc"));
        assert!(script.contains("gzip s_1.trv"));
        assert!(script.contains("gzip s_1.tre"));
        assert!(!script.contains("gzip s_1.trf"));
    }

    #[test]
    fn a_bare_integer_data_line_does_not_shift_later_blocks() {
        let imd =
            format!("{MULTIBATH_IMD}POSITIONRES\n#  NTPOR NTPORB NTPORS CPOR\n1 1 0 2.5E4\nEND\n");
        let out = apply_job_overrides(&imd, &job(&[("NTPOR", "0"), ("CPOR", "0.0")]));
        let posres = out.lines().rev().nth(1).expect("POSITIONRES data line");
        assert_eq!(
            posres.split_whitespace().collect::<Vec<_>>(),
            ["0", "1", "0", "0.0"]
        );
    }
}
