//! One module per program: the body that both the binary and the Python binding run.
//!
//! Each module exposes `USAGE`, `KNOWN`, `pub fn run(argv, out)` and `pub const PROGRAM`. A
//! program never reads the process arguments, never writes to the process stdout and never exits:
//! argv arrives already normalised (`@key` → `--key`, `@f` expanded), output goes to `out` and a
//! failure is an `Err`. `newen_io::run_main` is what turns that back into a command-line binary.

use newen_io::args::Program;

pub mod r#box;
pub mod conversion;
pub mod topology;
pub mod utilities;

// The categories group the bodies the way `src/bin/` and the `[[bin]]` sections of `Cargo.toml`
// group the binaries. A program is still reached by its own name alone — `programs::make_top`.
pub use conversion::{amber2gromos, pdb2g96, traj2pdb};
pub use r#box::{
    bin_box, build_box, check_box, copy_box, duplicate, explode, inbox, ion, ran_box,
    ran_solvation, sim_box, unify_box,
};
pub use topology::{
    addvirt_top, check_top, com_top, con_top, link_top, make_pt_top, make_sasa_top, make_top,
    pert_top, pt_top, red_top, top2build,
};
pub use utilities::{gch, mk_script, prep_eds, prep_noe, prep_posres, prep_xray};

/// Every program this crate provides.
pub static PROGRAMS: &[&Program] = &[
    // Topology construction
    &addvirt_top::PROGRAM,
    &check_top::PROGRAM,
    &com_top::PROGRAM,
    &con_top::PROGRAM,
    &link_top::PROGRAM,
    &make_pt_top::PROGRAM,
    &make_sasa_top::PROGRAM,
    &make_top::PROGRAM,
    &pert_top::PROGRAM,
    &pt_top::PROGRAM,
    &red_top::PROGRAM,
    &top2build::PROGRAM,
    // Box / coordinate preparation
    &bin_box::PROGRAM,
    &build_box::PROGRAM,
    &check_box::PROGRAM,
    &copy_box::PROGRAM,
    &duplicate::PROGRAM,
    &explode::PROGRAM,
    &inbox::PROGRAM,
    &ion::PROGRAM,
    &ran_box::PROGRAM,
    &ran_solvation::PROGRAM,
    &sim_box::PROGRAM,
    &unify_box::PROGRAM,
    // Coordinate conversion
    &amber2gromos::PROGRAM,
    &pdb2g96::PROGRAM,
    &traj2pdb::PROGRAM,
    // Utilities
    &gch::PROGRAM,
    &mk_script::PROGRAM,
    &prep_eds::PROGRAM,
    &prep_noe::PROGRAM,
    &prep_posres::PROGRAM,
    &prep_xray::PROGRAM,
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

/// Parse `argv` with clap, keeping the command line's behaviour for `--help` / `--version`.
///
/// clap reports those two as an `Error` it would print to stdout and exit 0 on. `Ok(None)` means
/// exactly that happened — the text is already in `out` and the program is done; the caller
/// returns `Ok(())`. Anything else is a real parse failure and becomes the `Err` `run` returns.
pub(crate) fn parse_clap<A: clap::Parser>(
    argv: &[String],
    out: &mut dyn std::io::Write,
) -> Result<Option<A>, String> {
    match A::try_parse_from(argv) {
        Ok(a) => Ok(Some(a)),
        Err(e) if e.use_stderr() => Err(e.render().to_string().trim_end().to_string()),
        Err(e) => {
            write!(out, "{}", e.render()).map_err(io_err)?;
            Ok(None)
        },
    }
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

    /// Every registered program is also a `[[bin]]`, and every `[[bin]]` that has a module here is
    /// registered — a program added to `programs/` without a registry entry is invisible to
    /// Python, and a registry entry without a binary breaks the command line.
    #[test]
    fn every_program_is_a_bin() {
        let manifest = include_str!("../../Cargo.toml");
        let bins: Vec<&str> = manifest
            .lines()
            .map(str::trim)
            .filter_map(|l| l.strip_prefix("name = \""))
            .filter_map(|l| l.strip_suffix('"'))
            .collect();
        for p in PROGRAMS {
            assert!(
                bins.contains(&p.name),
                "{} is registered but has no [[bin]] in Cargo.toml",
                p.name
            );
        }
    }

    /// A program's body and its shim must be filed under the same category —
    /// `src/programs/<cat>/<name>.rs` beside `src/bin/<cat>/<name>.rs`. The two trees drifted once
    /// before, `programs/` flat while `bin/` was grouped.
    #[test]
    fn program_and_bin_agree_on_category() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let manifest = include_str!("../../Cargo.toml");
        let mut name = "";
        for line in manifest.lines().map(str::trim) {
            if let Some(n) = line
                .strip_prefix("name = \"")
                .and_then(|l| l.strip_suffix('"'))
            {
                name = n;
            } else if let Some(path) = line
                .strip_prefix("path = \"")
                .and_then(|l| l.strip_suffix('"'))
            {
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
