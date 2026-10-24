//! Parser for GROMOS mk_script library/template files (.lib)
//!
//! Template files define filename patterns and miscellaneous settings
//! for generating simulation run scripts. Variables like `%system%` and
//! `%number%` are expanded when generating actual filenames.

use crate::IoError;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// Parsed mk_script template library.
#[derive(Debug, Clone)]
pub struct ScriptTemplate {
    /// Filename templates: type → pattern (e.g. "script" → "%system%_%number%.run")
    pub filenames: HashMap<String, String>,
    /// Miscellaneous settings: type → value (e.g. "workdir" → "...")
    pub misc: HashMap<String, String>,
}

impl ScriptTemplate {
    /// Get a filename template, falling back to a default.
    pub fn filename(&self, kind: &str) -> String {
        self.filenames
            .get(kind)
            .cloned()
            .unwrap_or_else(|| default_filename(kind))
    }

    /// Expand a template string, replacing %system% and %number%.
    pub fn expand(template: &str, system: &str, number: usize) -> String {
        template
            .replace("%system%", system)
            .replace("%number%", &format!("{}", number))
    }

    /// Get the expanded filename for a given kind, system, and job number.
    pub fn expanded_filename(&self, kind: &str, system: &str, number: usize) -> String {
        Self::expand(&self.filename(kind), system, number)
    }
}

fn default_filename(kind: &str) -> String {
    match kind {
        "script" => "jmd%system%_%number%.sh".to_string(),
        "input" => "imd%system%_%number%.dat".to_string(),
        "output" => "omd%system%_%number%.out".to_string(),
        "coord" => "o%system%sxmd_%number%.dat".to_string(),
        "outtrx" => "o%system%trmd_%number%.dat".to_string(),
        "outtrv" => "o%system%tvmd_%number%.dat".to_string(),
        "outtre" => "o%system%temd_%number%.dat".to_string(),
        "outtrg" => "o%system%tgmd_%number%.dat".to_string(),
        other => format!("{}_%system%_%number%.dat", other),
    }
}

/// Parse a GROMOS mk_script template library file.
pub fn read_script_template<P: AsRef<Path>>(path: P) -> Result<ScriptTemplate, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;

    let mut filenames = HashMap::new();
    let mut misc = HashMap::new();

    enum Section {
        None,
        Filenames,
        Miscellaneous,
        Skip,
    }
    let mut section = Section::None;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match trimmed {
            "TITLE" | "LINKADDITION" => {
                section = Section::Skip;
                continue;
            },
            "FILENAMES" => {
                section = Section::Filenames;
                continue;
            },
            "MISCELLANEOUS" => {
                section = Section::Miscellaneous;
                continue;
            },
            "END" => {
                section = Section::None;
                continue;
            },
            _ => {},
        }

        match section {
            Section::Filenames => {
                let parts: Vec<&str> = trimmed.splitn(2, char::is_whitespace).collect();
                if parts.len() == 2 {
                    filenames.insert(parts[0].to_string(), parts[1].trim().to_string());
                }
            },
            Section::Miscellaneous => {
                let parts: Vec<&str> = trimmed.splitn(2, char::is_whitespace).collect();
                if parts.len() == 2 {
                    misc.insert(parts[0].to_string(), parts[1].trim().to_string());
                }
            },
            _ => {},
        }
    }

    Ok(ScriptTemplate { filenames, misc })
}
