//! Parser for GROMOS mk_script job list files (.jobs)
//!
//! Job list files define a series of simulation jobs with per-job parameter
//! overrides. The JOBSCRIPTS block header line defines which IMD parameters
//! are overridden — any valid IMD parameter name can appear as a column.
//!
//! Fixed columns: `job_id` (first), `subdir` (second-to-last), `run_after` (last).
//! All columns between `job_id` and `subdir` are parameter overrides.

use crate::IoError;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// A single job specification from a JOBSCRIPTS block.
///
/// Parameters are stored as string key-value pairs because the set of
/// overridable parameters is open-ended (any IMD parameter name is valid).
#[derive(Debug, Clone)]
pub struct JobSpec {
    pub job_id: usize,
    /// Parameter overrides: column_name → value (as string)
    pub params: HashMap<String, String>,
    /// Subdirectory for this job ("." = current)
    pub subdir: String,
    /// Job to run after (0 = first/independent)
    pub run_after: usize,
}

impl JobSpec {
    /// Get a parameter as f64, or None if not present.
    pub fn get_f64(&self, key: &str) -> Option<f64> {
        self.params.get(key).and_then(|v| v.parse().ok())
    }

    /// Get a parameter as i32, or None if not present.
    pub fn get_i32(&self, key: &str) -> Option<i32> {
        self.params.get(key).and_then(|v| v.parse().ok())
    }

    /// Get a parameter as string, or None if not present.
    pub fn get_str(&self, key: &str) -> Option<&str> {
        self.params.get(key).map(|s| s.as_str())
    }
}

/// Parsed job list file.
#[derive(Debug, Clone)]
pub struct JobList {
    pub title: String,
    /// Column headers (parameter names), excluding job_id/subdir/run_after
    pub headers: Vec<String>,
    pub jobs: Vec<JobSpec>,
}

/// Parse a GROMOS mk_script job list file.
///
/// Format:
/// ```text
/// TITLE
/// description
/// END
/// JOBSCRIPTS
/// job_id PARAM1 PARAM2 ... subdir run_after
/// 1      val1   val2   ... .      0
/// ...
/// END
/// ```
///
/// The first column must be `job_id`, the last two must be `subdir` and `run_after`.
/// All columns in between are parameter names that will be used to override
/// values in the base IMD file.
pub fn read_joblist<P: AsRef<Path>>(path: P) -> Result<JobList, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;

    let mut title = String::new();
    let mut headers: Vec<String> = Vec::new();
    let mut jobs = Vec::new();

    enum Section {
        None,
        Title,
        JobScripts,
    }
    let mut section = Section::None;
    let mut header_parsed = false;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match trimmed {
            "TITLE" => {
                section = Section::Title;
                continue;
            },
            "JOBSCRIPTS" => {
                section = Section::JobScripts;
                header_parsed = false;
                continue;
            },
            "END" => {
                section = Section::None;
                continue;
            },
            _ => {},
        }

        match section {
            Section::Title => {
                if !title.is_empty() {
                    title.push('\n');
                }
                title.push_str(trimmed);
            },
            Section::JobScripts => {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();

                if !header_parsed {
                    if parts.len() < 3 {
                        return Err(IoError::ParseError(
                            "JOBSCRIPTS header needs at least: job_id ... subdir run_after"
                                .to_string(),
                        ));
                    }
                    if parts[0] != "job_id" {
                        return Err(IoError::ParseError(format!(
                            "JOBSCRIPTS header must start with 'job_id', got '{}'",
                            parts[0]
                        )));
                    }
                    if parts[parts.len() - 1] != "run_after" || parts[parts.len() - 2] != "subdir" {
                        return Err(IoError::ParseError(
                            "JOBSCRIPTS header must end with 'subdir run_after'".to_string(),
                        ));
                    }
                    headers = parts[1..parts.len() - 2]
                        .iter()
                        .map(|s| s.to_string())
                        .collect();
                    header_parsed = true;
                    continue;
                }

                if parts.len() != headers.len() + 3 {
                    return Err(IoError::ParseError(format!(
                        "JOBSCRIPTS data line has {} columns, expected {}",
                        parts.len(),
                        headers.len() + 3,
                    )));
                }

                let job_id: usize = parts[0]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid job_id: {}", parts[0])))?;

                let mut params = HashMap::new();
                for (i, header) in headers.iter().enumerate() {
                    params.insert(header.clone(), parts[i + 1].to_string());
                }

                let subdir = parts[parts.len() - 2].to_string();
                let run_after: usize = parts[parts.len() - 1].parse().map_err(|_| {
                    IoError::ParseError(format!("Invalid run_after: {}", parts[parts.len() - 1]))
                })?;

                jobs.push(JobSpec {
                    job_id,
                    params,
                    subdir,
                    run_after,
                });
            },
            _ => {},
        }
    }

    if jobs.is_empty() {
        return Err(IoError::FormatError(
            "No jobs found in JOBSCRIPTS block".to_string(),
        ));
    }

    Ok(JobList {
        title,
        headers,
        jobs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_joblist() {
        let content = "\
TITLE
test jobs
END
JOBSCRIPTS
job_id NTIVEL TEMPI TEMP0[1] TEMP0[2] COUPLE NTPOR CPOR subdir run_after
1      1      60.0      60.0     60.0       1     1  2.5E4    .      0
2      0      0.0      120.0    120.0       1     1  2.5E3    .      1
END
";
        let path = "/tmp/test_jobs.dat";
        std::fs::write(path, content).unwrap();

        let jl = read_joblist(path).unwrap();
        assert_eq!(jl.jobs.len(), 2);
        assert_eq!(
            jl.headers,
            vec!["NTIVEL", "TEMPI", "TEMP0[1]", "TEMP0[2]", "COUPLE", "NTPOR", "CPOR"]
        );

        assert_eq!(jl.jobs[0].job_id, 1);
        assert_eq!(jl.jobs[0].get_i32("NTIVEL"), Some(1));
        assert_eq!(jl.jobs[0].get_f64("TEMPI"), Some(60.0));
        assert_eq!(jl.jobs[0].get_f64("CPOR"), Some(25000.0));
        assert_eq!(jl.jobs[0].subdir, ".");

        assert_eq!(jl.jobs[1].job_id, 2);
        assert_eq!(jl.jobs[1].get_i32("NTIVEL"), Some(0));
        assert_eq!(jl.jobs[1].run_after, 1);
    }

    #[test]
    fn test_custom_columns() {
        let content = "\
TITLE
custom test
END
JOBSCRIPTS
job_id NSTLIM DT NLRELE subdir run_after
1      50000  0.001 2   .  0
2      100000 0.002 1   .  1
END
";
        let path = "/tmp/test_jobs_custom.dat";
        std::fs::write(path, content).unwrap();

        let jl = read_joblist(path).unwrap();
        assert_eq!(jl.headers, vec!["NSTLIM", "DT", "NLRELE"]);
        assert_eq!(jl.jobs[0].get_i32("NSTLIM"), Some(50000));
        assert_eq!(jl.jobs[0].get_f64("DT"), Some(0.001));
        assert_eq!(jl.jobs[1].get_i32("NLRELE"), Some(1));
    }
}
