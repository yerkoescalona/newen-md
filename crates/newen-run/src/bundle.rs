//! Run bundles — a run on disk (3.9 A5).
//!
//! An `.imd` describes a run only together with its topology, coordinates and auxiliary files.
//! The reference systems already carry that as `input.toml`:
//!
//! ```toml
//! [input]
//! topology      = "system.topo"
//! configuration = "start.cnf"
//! parameters    = "run.imd"        # GROMOS front-end …
//! recipe        = "run.recipe.toml" # … or the recipe itself (wins when both are present)
//! pttopo        = "system.ptp"     # optional
//! posresspec    = "posres.por"     # optional
//! refpos        = "refpos.rpr"     # optional
//! distrest      = "restraints.dsr" # optional
//! ```
//!
//! Paths are relative to the bundle file. Other sections (`[system]`, `[expected]`, …) are
//! ignored, so the reference directories are valid bundles as they are.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use newen_io::imd::{read_imd_file, ImdParameters};

use crate::recipe::{Diagnostics, PassthroughPolicy};
use crate::{RunError, RunInputs, RunRecipe};

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
struct InputSection {
    topology: Option<String>,
    configuration: Option<String>,
    parameters: Option<String>,
    recipe: Option<String>,
    pttopo: Option<String>,
    posresspec: Option<String>,
    refpos: Option<String>,
    distrest: Option<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
struct BundleFile {
    input: InputSection,
}

/// A run's files, with every path resolved against the bundle's directory.
#[derive(Debug, Clone, PartialEq)]
pub struct RunBundle {
    pub dir: PathBuf,
    pub topology: PathBuf,
    pub configuration: PathBuf,
    /// The GROMOS parameter file, if the bundle names one.
    pub parameters: Option<PathBuf>,
    /// The recipe file, if the bundle names one (takes precedence over `parameters`).
    pub recipe: Option<PathBuf>,
    /// Auxiliary inputs, resolved.
    pub inputs: RunInputs,
}

fn resolve(dir: &Path, p: &Option<String>) -> Option<PathBuf> {
    p.as_ref().map(|s| dir.join(s))
}

/// Read a bundle description (`input.toml` shape) and resolve its paths.
pub fn read_bundle<P: AsRef<Path>>(path: P) -> Result<RunBundle, RunError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|e| RunError::Io {
        what: format!("bundle '{}'", path.display()),
        source: newen_io::IoError::Io(e),
    })?;
    let file: BundleFile = toml::from_str(&text).map_err(|e| RunError::Serde(e.to_string()))?;
    let dir = path.parent().unwrap_or(Path::new(".")).to_path_buf();
    let i = &file.input;
    let missing =
        |what: &str| RunError::Inconsistent(format!("bundle '{}' names no {what}", path.display()));
    Ok(RunBundle {
        topology: resolve(&dir, &i.topology).ok_or_else(|| missing("topology"))?,
        configuration: resolve(&dir, &i.configuration).ok_or_else(|| missing("configuration"))?,
        parameters: resolve(&dir, &i.parameters),
        recipe: resolve(&dir, &i.recipe),
        inputs: RunInputs {
            pttopo: resolve(&dir, &i.pttopo),
            posresspec: resolve(&dir, &i.posresspec),
            refpos: resolve(&dir, &i.refpos),
            distrest: resolve(&dir, &i.distrest),
        },
        dir,
    })
}

/// The recipe a bundle describes: its recipe file if present, else its parameter file through
/// `RunRecipe::from_imd_with`. The bundle's auxiliary paths are written into `recipe.inputs`.
/// Read a GROMOS `.imd` file — the one place a parameter file is opened outside `newen-io`
/// (G6): the `md` binary, `Recipe.from_imd` and the bundle loader all come through here.
pub fn read_imd<P: AsRef<Path>>(path: P) -> Result<ImdParameters, RunError> {
    let path = path.as_ref();
    read_imd_file(path).map_err(|e| RunError::Io {
        what: format!("parameters '{}'", path.display()),
        source: e,
    })
}

pub fn load_bundle<P: AsRef<Path>>(
    path: P,
    policy: &PassthroughPolicy,
) -> Result<(RunBundle, RunRecipe, Diagnostics), RunError> {
    let bundle = read_bundle(path)?;
    let (mut recipe, diagnostics) = if let Some(recipe_path) = &bundle.recipe {
        let text = std::fs::read_to_string(recipe_path).map_err(|e| RunError::Io {
            what: format!("recipe '{}'", recipe_path.display()),
            source: newen_io::IoError::Io(e),
        })?;
        (RunRecipe::from_toml(&text)?, Diagnostics::default())
    } else if let Some(params) = &bundle.parameters {
        RunRecipe::from_imd_with(&read_imd(params)?, policy)?
    } else {
        return Err(RunError::Inconsistent(format!(
            "bundle '{}' names neither a recipe nor a parameter file",
            bundle.dir.display()
        )));
    };
    recipe.inputs = bundle.inputs.clone();
    Ok((bundle, recipe, diagnostics))
}

/// Write a run as a bundle: `input.toml`, `run.recipe.toml` and `run.imd` (the GROMOS form, so
/// gromosXX can run the same bundle) into `dir`. Topology/coordinate/auxiliary paths are
/// recorded as given (relative paths are kept relative to `dir`).
pub fn write_bundle(
    dir: &Path,
    recipe: &RunRecipe,
    topology: &Path,
    configuration: &Path,
    n_atoms: Option<usize>,
) -> Result<PathBuf, RunError> {
    let io = |what: &str, e: std::io::Error| RunError::Io {
        what: format!("{what} in '{}'", dir.display()),
        source: newen_io::IoError::Io(e),
    };
    std::fs::create_dir_all(dir).map_err(|e| io("bundle directory", e))?;
    std::fs::write(dir.join("run.recipe.toml"), recipe.to_toml()?)
        .map_err(|e| io("run.recipe.toml", e))?;
    std::fs::write(dir.join("run.imd"), recipe.to_imd_string(n_atoms))
        .map_err(|e| io("run.imd", e))?;
    let rel = |p: &Path| -> String {
        p.strip_prefix(dir)
            .unwrap_or(p)
            .to_string_lossy()
            .into_owned()
    };
    let file = BundleFile {
        input: InputSection {
            topology: Some(rel(topology)),
            configuration: Some(rel(configuration)),
            parameters: Some("run.imd".into()),
            recipe: Some("run.recipe.toml".into()),
            pttopo: recipe.inputs.pttopo.as_deref().map(rel),
            posresspec: recipe.inputs.posresspec.as_deref().map(rel),
            refpos: recipe.inputs.refpos.as_deref().map(rel),
            distrest: recipe.inputs.distrest.as_deref().map(rel),
        },
    };
    let text = toml::to_string_pretty(&file).map_err(|e| RunError::Serde(e.to_string()))?;
    let path = dir.join("input.toml");
    std::fs::write(&path, text).map_err(|e| io("input.toml", e))?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_written_bundle_reads_back_with_resolved_paths() {
        let dir = std::env::temp_dir().join(format!("newen_bundle_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let topo = dir.join("sys.top");
        let conf = dir.join("sys.cnf");
        std::fs::write(&topo, "").unwrap();
        std::fs::write(&conf, "").unwrap();
        let recipe = RunRecipe::default();
        let manifest = write_bundle(&dir, &recipe, &topo, &conf, Some(10)).unwrap();
        assert_eq!(manifest, dir.join("input.toml"));
        assert!(dir.join("run.recipe.toml").exists() && dir.join("run.imd").exists());
        let bundle = read_bundle(&manifest).unwrap();
        assert_eq!(bundle.topology, topo);
        assert_eq!(bundle.configuration, conf);
        assert_eq!(
            bundle.recipe.as_deref(),
            Some(dir.join("run.recipe.toml").as_path())
        );
        let (_, loaded, _) = load_bundle(&manifest, &PassthroughPolicy::default()).unwrap();
        assert_eq!(loaded.forcefield, recipe.forcefield);
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn a_missing_manifest_is_an_io_error() {
        assert!(matches!(
            read_bundle("/nonexistent/input.toml"),
            Err(RunError::Io { .. })
        ));
    }
}
