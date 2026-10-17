//! mk_script - Generate GROMOS simulation scripts from job lists and templates
//!
//! Two modes of operation:
//!   @joblist mode:  Read a .jobs file with per-job parameter overrides
//!                   (equilibration, custom protocols)
//!   @script mode:   Generate N identical sequential scripts
//!                   (production runs)
//!
//! Usage:
//!   mk_script @sys eq_GB3 @bin ./md @dir ./eq @files topo t.top input eq.imd
//!             coord c.cnf posresspec c.por refpos c.rpr
//!             @template mk_script.lib @joblist equilibration.jobs
//!
//!   mk_script @sys md_GB3 @bin ./md @dir ./md @files topo t.top input md.imd
//!             coord c.cnf @template mk_script.lib @script 1 21

use newen_io::args::{Program, Status};
use newen_io::jobs::read_joblist;
use newen_io::mk_script::{apply_job_overrides, generate_run_script, MkScriptConfig, WriteTraj};
use newen_io::script_template::read_script_template;
use std::collections::HashMap;
use std::fs;
use std::io::Write;

pub const USAGE: &str = "mk_script - Generate simulation run scripts from a template\n\nUsage: mk_script @sys <name> @bin <md binary> @dir <workdir> @files <key value ...>\n                 @template <file> [@joblist <file>] [@script <first> <n>]\n\nArguments:\n  @sys       System name\n  @bin       Path to the md binary\n  @dir       Working directory\n  @files     Input files as key/value pairs (input, topo, coord, ...)\n  @template  Script template library (required)\n  @joblist   Job list for an equilibration series\n  @script    First script number and how many to write\n\nThe scripts are written to disk; stdout carries the progress report.";

pub const KNOWN: &[&str] = &[
    "sys", "bin", "dir", "files", "template", "joblist", "version", "script", "queue", "cmd",
];

pub const PROGRAM: Program = Program {
    name: "mk_script",
    known: KNOWN,
    usage: USAGE,
    status: Status::StandIn,
    run,
};

#[allow(dead_code)]
fn print_usage() {
    eprintln!("mk_script - Generate GROMOS simulation scripts");
    eprintln!();
    eprintln!("Usage:");
    eprintln!("  mk_script @sys <name> @bin <md_binary> @dir <workdir>");
    eprintln!("            @files topo <file> input <file> coord <file>");
    eprintln!("                   [posresspec <file>] [refpos <file>]");
    eprintln!("            @template <mk_script.lib>");
    eprintln!("            @joblist <jobs_file>      (equilibration mode)");
    eprintln!("       or:  @script <first> <last>    (production mode)");
}

/// Parse the @-style arguments for mk_script.
/// This tool uses its own parser because @files takes sub-key/value pairs.
fn parse_mk_args(argv: &[String]) -> Result<MkScriptArgs, String> {
    let mut args = MkScriptArgs::default();
    let mut i = 1;

    while i < argv.len() {
        match argv[i].as_str() {
            "--sys" => {
                i += 1;
                args.system = argv.get(i).cloned().ok_or("Missing @sys value")?;
            },
            "--bin" => {
                i += 1;
                args.bin = argv.get(i).cloned().ok_or("Missing @bin value")?;
            },
            "--dir" => {
                i += 1;
                args.dir = argv.get(i).cloned().ok_or("Missing @dir value")?;
            },
            "--files" => {
                // Read key-value pairs until we hit another --flag
                i += 1;
                while i < argv.len() && !argv[i].starts_with("--") {
                    let key = argv[i].clone();
                    i += 1;
                    let val = argv
                        .get(i)
                        .cloned()
                        .ok_or(format!("Missing value for @files key '{}'", key))?;
                    args.files.insert(key, val);
                    i += 1;
                }
                continue; // Don't increment i again
            },
            "--template" => {
                i += 1;
                args.template = argv.get(i).cloned().ok_or("Missing @template value")?;
            },
            "--version" => {
                i += 1;
                args.version = argv.get(i).cloned().ok_or("Missing @version value")?;
            },
            "--joblist" => {
                i += 1;
                args.joblist = argv.get(i).cloned();
            },
            "--script" => {
                i += 1;
                let first: usize = argv
                    .get(i)
                    .ok_or("Missing @script first")?
                    .parse()
                    .map_err(|_| "Invalid @script first")?;
                i += 1;
                let last: usize = argv
                    .get(i)
                    .ok_or("Missing @script last")?
                    .parse()
                    .map_err(|_| "Invalid @script last")?;
                args.script_range = Some((first, last));
            },
            // gromos++ prints the usage and stops; `run` reports it the same way every
            // Arguments program does
            "-h" | "--help" => return Err(USAGE.to_string()),
            other => {
                return Err(format!("Unknown argument: {}", other));
            },
        }
        i += 1;
    }

    if args.system.is_empty() {
        return Err("Missing required argument: @sys".to_string());
    }
    if args.template.is_empty() {
        return Err("Missing required argument: @template".to_string());
    }
    if args.joblist.is_none() && args.script_range.is_none() {
        return Err("Either @joblist or @script must be specified".to_string());
    }

    Ok(args)
}

#[derive(Default)]
struct MkScriptArgs {
    system: String,
    bin: String,
    dir: String,
    files: HashMap<String, String>,
    template: String,
    version: String,
    joblist: Option<String>,
    script_range: Option<(usize, usize)>,
}

// Everything mk_script reports goes to stderr and its scripts go to disk; stdout
// carries nothing, so the sink is unused.
pub fn run(argv: &[String], _out: &mut dyn Write) -> Result<(), String> {
    if argv.len() < 2 || argv.iter().any(|a| a == "-h" || a == "--help") {
        return Err(USAGE.to_string());
    }

    let mk_args = parse_mk_args(argv).map_err(|e| format!("Error: {e}\n\n{USAGE}"))?;

    // Read template library
    let template = read_script_template(&mk_args.template)
        .map_err(|e| format!("Error reading template: {e}"))?;

    // Read base IMD
    let imd_path = mk_args
        .files
        .get("input")
        .cloned()
        .unwrap_or_else(|| "input.imd".to_string());
    let base_imd = fs::read_to_string(&imd_path)
        .map_err(|e| format!("Error reading IMD file '{imd_path}': {e}"))?;

    let config = MkScriptConfig {
        system: mk_args.system.clone(),
        bin: mk_args.bin.clone(),
        dir: mk_args.dir.clone(),
        files: mk_args.files.clone(),
    };

    eprintln!("mk_script: system={}", mk_args.system);

    if let Some(joblist_path) = &mk_args.joblist {
        // === Joblist mode (equilibration) ===
        let joblist = read_joblist(joblist_path)
            .map_err(|e| format!("Error reading joblist '{joblist_path}': {e}"))?;

        eprintln!(
            "  {} jobs from {}, headers: {:?}",
            joblist.jobs.len(),
            joblist_path,
            joblist.headers
        );

        for (idx, job) in joblist.jobs.iter().enumerate() {
            let prev = if idx > 0 {
                Some(&joblist.jobs[idx - 1])
            } else {
                None
            };

            // Generate modified IMD
            let modified_imd = apply_job_overrides(&base_imd, job);
            let imd_name = template.expanded_filename("input", &mk_args.system, job.job_id);
            fs::write(&imd_name, &modified_imd)
                .map_err(|e| format!("Error writing IMD '{imd_name}': {e}"))?;

            // Generate run script
            let script = generate_run_script(
                &config,
                &template,
                job,
                prev,
                &WriteTraj::from_imd(&modified_imd),
            );
            let script_name = template.expanded_filename("script", &mk_args.system, job.job_id);
            fs::write(&script_name, &script)
                .map_err(|e| format!("Error writing script '{script_name}': {e}"))?;

            // Make script executable
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let _ = fs::set_permissions(&script_name, fs::Permissions::from_mode(0o755));
            }

            eprintln!("  Job {}: {} + {}", job.job_id, imd_name, script_name);
        }
    } else if let Some((first, last)) = mk_args.script_range {
        // === Script mode (production) ===
        eprintln!("  Scripts {} to {}", first, last);

        for num in first..=last {
            // For production: no IMD overrides, just generate sequential scripts
            let dummy_job = newen_io::jobs::JobSpec {
                job_id: num,
                params: HashMap::new(),
                subdir: ".".to_string(),
                run_after: if num > first { num - 1 } else { 0 },
            };
            let prev = if num > first {
                Some(newen_io::jobs::JobSpec {
                    job_id: num - 1,
                    params: HashMap::new(),
                    subdir: ".".to_string(),
                    run_after: 0,
                })
            } else {
                None
            };

            // Copy base IMD as-is for each job
            let imd_name = template.expanded_filename("input", &mk_args.system, num);
            fs::write(&imd_name, &base_imd)
                .map_err(|e| format!("Error writing IMD '{imd_name}': {e}"))?;

            let script = generate_run_script(
                &config,
                &template,
                &dummy_job,
                prev.as_ref(),
                &WriteTraj::from_imd(&base_imd),
            );
            let script_name = template.expanded_filename("script", &mk_args.system, num);
            fs::write(&script_name, &script)
                .map_err(|e| format!("Error writing script '{script_name}': {e}"))?;

            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let _ = fs::set_permissions(&script_name, fs::Permissions::from_mode(0o755));
            }
        }

        eprintln!("  Generated {} scripts", last - first + 1);
    }

    eprintln!("Done!");
    Ok(())
}
