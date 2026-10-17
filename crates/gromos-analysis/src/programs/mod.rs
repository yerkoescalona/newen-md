//! One module per program: the body that both the binary and the Python binding run.
//!
//! Each module exposes `USAGE`, `KNOWN`, `pub fn run(argv, out)` and `pub const PROGRAM`. A
//! program never reads the process arguments, never writes to the process stdout and never exits:
//! argv arrives already normalised (`@key` → `--key`, `@f` expanded), output goes to `out` and a
//! failure is an `Err`. `gromos_io::run_main` is what turns that back into a command-line binary.
//!
//! The modules are grouped the way `src/bin/` and the `[[bin]]` sections of `Cargo.toml` group the
//! binaries, and every program is re-exported flat, so a caller writes `programs::rmsd`. Binaries
//! whose body is still in `main()` have no module here yet.

use gromos_io::args::Program;

pub mod distribution;
pub mod energy;
pub mod free_energy;
pub mod noe;
pub mod special;
pub mod structural;
pub mod trajectory;

// Every program is reachable by its own name as well as through its category.
pub use distribution::{bilayer_dist, bilayer_oparam, rdf};
pub use energy::{ene_ana, ener};
pub use free_energy::{
    bar, dfmult, dg_ener, dgslv_pbsolv, eds_update_1, eds_update_2, ext_ti_ana, ext_ti_merge, mbar,
};
pub use noe::{jepot, jval, nhoparam};
pub use special::{dfgrid, gca, matrix_overlap};
pub use structural::{atominfo, edyn, pocket, rmsd, rmsf};
pub use trajectory::{frameout, mdf, tser};

/// Every program this crate provides.
pub static PROGRAMS: &[&Program] = &[
    // Structural
    &atominfo::PROGRAM,
    &edyn::PROGRAM,
    &pocket::PROGRAM,
    &rmsd::PROGRAM,
    &rmsf::PROGRAM,
    // Energy
    &ene_ana::PROGRAM,
    &ener::PROGRAM,
    // Distribution
    &bilayer_dist::PROGRAM,
    &bilayer_oparam::PROGRAM,
    &rdf::PROGRAM,
    // Free energy
    &bar::PROGRAM,
    &dfmult::PROGRAM,
    &dg_ener::PROGRAM,
    &dgslv_pbsolv::PROGRAM,
    &eds_update_1::PROGRAM,
    &eds_update_2::PROGRAM,
    &ext_ti_ana::PROGRAM,
    &ext_ti_merge::PROGRAM,
    &mbar::PROGRAM,
    // Trajectory manipulation
    &frameout::PROGRAM,
    &mdf::PROGRAM,
    &tser::PROGRAM,
    // NOE / NMR
    &jepot::PROGRAM,
    &jval::PROGRAM,
    &nhoparam::PROGRAM,
    // Special
    &dfgrid::PROGRAM,
    &gca::PROGRAM,
    &matrix_overlap::PROGRAM,
];

/// The program called `name`, if this build has it.
pub fn find(name: &str) -> Option<&'static Program> {
    PROGRAMS.iter().copied().find(|p| p.name == name)
}

/// A failed write to `out`. Programs print into a caller-owned sink, so this is a real error, not
/// the broken-pipe-is-fine case `println!` papers over.
pub(crate) fn io_err(e: std::io::Error) -> String {
    format!("cannot write output: {e}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_are_unique_and_findable() {
        let mut names: Vec<&str> = PROGRAMS.iter().map(|p| p.name).collect();
        let n = names.len();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), n, "duplicate program name in PROGRAMS");
        for p in PROGRAMS {
            assert_eq!(find(p.name).map(|f| f.name), Some(p.name));
        }
        assert!(find("no_such_program").is_none());
    }

    /// The `[[bin]]` path of every registered program, from the manifest.
    fn bin_paths() -> std::collections::HashMap<&'static str, &'static str> {
        let manifest = include_str!("../../Cargo.toml");
        let mut paths = std::collections::HashMap::new();
        let mut name = "";
        for line in manifest.lines().map(str::trim) {
            if let Some(n) = line
                .strip_prefix("name = \"")
                .and_then(|l| l.strip_suffix('"'))
            {
                name = n;
            } else if let Some(p) = line
                .strip_prefix("path = \"")
                .and_then(|l| l.strip_suffix('"'))
            {
                paths.insert(name, p);
            }
        }
        paths
    }

    /// A registered program is also a `[[bin]]` — a registry entry without a binary breaks the
    /// command line, and a program added without a registry entry is invisible to Python.
    #[test]
    fn every_program_is_a_bin() {
        let paths = bin_paths();
        for p in PROGRAMS {
            assert!(
                paths.contains_key(p.name),
                "{} is registered but has no [[bin]] in Cargo.toml",
                p.name
            );
        }
    }

    /// A program's body and its shim must be filed under the same category —
    /// `src/programs/<cat>/<name>.rs` beside `src/bin/<cat>/<name>.rs`. Only registered programs
    /// are checked: a binary whose body is still in `main()` has no module to be filed anywhere.
    #[test]
    fn program_and_bin_agree_on_category() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        for (name, path) in bin_paths() {
            if find(name).is_none() {
                continue;
            }
            let category = path
                .strip_prefix("src/bin/")
                .and_then(|p| p.split('/').next())
                .unwrap_or_else(|| {
                    panic!("{name}'s [[bin]] path is not under src/bin/<category>/")
                });
            let body = root.join(format!("src/programs/{category}/{name}.rs"));
            assert!(
                body.exists(),
                "{name}'s shim is in bin/{category}/ but its body is not in programs/{category}/"
            );
        }
    }

    /// `run` must pass `PROGRAM.known` to `Arguments`, or the registry advertises keys the program
    /// rejects. Cheap proxy: every program accepts its own advertised keys on an empty value list
    /// far enough to fail on something other than "unknown argument".
    #[test]
    fn advertised_keys_are_accepted() {
        for p in PROGRAMS {
            for key in p.known {
                let argv: Vec<String> = vec![p.name.to_string(), format!("--{key}")];
                let mut out = Vec::new();
                if let Err(e) = (p.run)(&argv, &mut out) {
                    assert!(
                        !e.contains("unknown argument"),
                        "{} rejects its own advertised @{key}: {e}",
                        p.name
                    );
                }
            }
        }
    }
}
