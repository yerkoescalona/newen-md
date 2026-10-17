//! 9a-1: Empirical margin measurement for CellList vs Standard pairlist.
//!
//! Runs `water_216_box` (rectangular, 648 atoms, 100 steps) twice — once
//! with `ALGORITHM standard` and once with `ALGORITHM grid_cell` — and
//! measures the per-step relative energy deviation |ΔE_total|/|E_total|.
//!
//! The pair *set* is proven identical (see newen-core pairlist unit tests).
//! The open question is whether the reordered summation stays < 1e-8 over
//! a full 100-step trajectory under MD's positive Lyapunov exponent.
//!
//! Decision gate (recorded in the test output):
//!   margin ≪ 1e-8 → CellList is safe to use; lower the auto-select
//!             threshold from 5000 to a useful value in from_imd().
//!   margin ≥ 1e-8 → add canonical-order sort to CellList output so the
//!             pair iteration order matches Standard's exactly (bit-identical).
//!
//! This test always reports the observed margin in eprintln! output
//! (visible with `cargo test -- --nocapture`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn ref_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

fn md_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_md"))
}

// ─── Total energies from the native .tre layout ──────────────────────────────

fn parse_enertrj_etotal(path: &Path) -> Vec<f64> {
    newen_io::read_energy_trajectory_native(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .into_iter()
        .map(|f| f.total)
        .collect()
}

// ─── Run md binary with a patched PAIRLIST algorithm keyword ────────────────

fn run_with_algorithm(system: &str, algorithm_keyword: &str) -> Vec<f64> {
    let sys_dir = ref_root().join(system);
    assert!(
        sys_dir.exists(),
        "reference system not found: {}",
        sys_dir.display()
    );

    let toml = fs::read_to_string(sys_dir.join("input.toml")).expect("input.toml missing");

    fn toml_str(content: &str, key: &str) -> String {
        let prefix = format!("{key} = ");
        for line in content.lines() {
            let t = line.trim();
            if t.starts_with(&prefix) {
                if let (Some(a), Some(b)) = (t.find('"'), t.rfind('"')) {
                    if a < b {
                        return t[a + 1..b].to_string();
                    }
                }
            }
        }
        panic!("{key} not found in input.toml");
    }

    let topo = sys_dir.join(toml_str(&toml, "topology"));
    let conf = sys_dir.join(toml_str(&toml, "configuration"));
    let params_path = sys_dir.join(toml_str(&toml, "parameters"));

    let out = std::env::temp_dir().join(format!(
        "newen_margin_{}_{}_{}",
        system,
        algorithm_keyword,
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(&out).unwrap();

    // Patch the PAIRLIST ALGORITHM in a temp copy of the input file
    let original_input = fs::read_to_string(&params_path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", params_path.display()));

    let patched = patch_pairlist_algorithm(&original_input, algorithm_keyword);
    let patched_input = out.join("patched.in");
    fs::write(&patched_input, patched).unwrap();

    let tre = out.join("energies.tre");

    let result = Command::new(md_bin())
        .arg("@topo")
        .arg(&topo)
        .arg("@conf")
        .arg(&conf)
        .arg("@input")
        .arg(&patched_input)
        .arg("@fin")
        .arg(out.join("final.conf"))
        .arg("@tre")
        .arg(&tre)
        .arg("@trf")
        .arg(out.join("forces.trf"))
        .arg("@trc")
        .arg(out.join("trajectory.trc"))
        .output()
        .expect("failed to execute md binary");

    assert!(
        result.status.success(),
        "md failed for {system} algorithm={algorithm_keyword}\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr),
    );

    let totals = parse_enertrj_etotal(&tre);
    let _ = fs::remove_dir_all(&out);
    totals
}

/// Replace the ALGORITHM token on the PAIRLIST data line.
///
/// The GROMOS PAIRLIST block has one data line:
///   ALGORITHM  NSNB  RCUTP  RCUTL  SIZE  TYPE
/// We replace only the first whitespace-delimited token on that line.
fn patch_pairlist_algorithm(input: &str, new_algorithm: &str) -> String {
    let mut in_pairlist = false;
    let mut data_line_done = false;
    let mut out = String::with_capacity(input.len());

    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed == "PAIRLIST" {
            in_pairlist = true;
            data_line_done = false;
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if in_pairlist && trimmed == "END" {
            in_pairlist = false;
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if in_pairlist && !data_line_done && !trimmed.starts_with('#') && !trimmed.is_empty() {
            // This is the data line; replace the first token (ALGORITHM keyword)
            let mut tokens: Vec<&str> = trimmed.split_whitespace().collect();
            if !tokens.is_empty() {
                tokens[0] = new_algorithm;
            }
            // Preserve indentation: detect leading whitespace from original line
            let leading: String = line.chars().take_while(|c| c.is_whitespace()).collect();
            out.push_str(&leading);
            out.push_str(&tokens.join("        "));
            out.push('\n');
            data_line_done = true;
            continue;
        }
        out.push_str(line);
        out.push('\n');
    }
    out
}

// ─── Test ────────────────────────────────────────────────────────────────────

/// 9a-1: Measure per-step |ΔE_total|/|E_total| between Standard and CellList
/// over the full 100-step water_216_box trajectory.
///
/// Asserts the margin is below 1e-8 (the reference-test tolerance).
/// If this fails, apply canonical-order sort to CellList output (see §9a-1).
#[test]
fn celllist_vs_standard_margin_water_216() {
    let system = "water_216_box";

    let standard = run_with_algorithm(system, "standard");
    let cell_list = run_with_algorithm(system, "grid_cell");

    assert_eq!(
        standard.len(),
        cell_list.len(),
        "step count mismatch: standard={} cell_list={}",
        standard.len(),
        cell_list.len()
    );
    assert!(!standard.is_empty(), "no energy frames produced");

    let mut max_rel = 0.0f64;
    let mut max_abs = 0.0f64;
    let mut max_step = 0usize;

    for (step, (e_std, e_cl)) in standard.iter().zip(&cell_list).enumerate() {
        let abs_diff = (e_cl - e_std).abs();
        let rel = if e_std.abs() > 1e-30 {
            abs_diff / e_std.abs()
        } else {
            abs_diff
        };
        if rel > max_rel {
            max_rel = rel;
            max_abs = abs_diff;
            max_step = step;
        }
    }

    eprintln!(
        "\n9a-1 margin ({} steps): max |ΔE|/|E| = {:.3e}  (abs={:.3e} kJ/mol at step {})",
        standard.len(),
        max_rel,
        max_abs,
        max_step
    );
    eprintln!(
        "Reference tolerance: 1e-8 | Margin is {} the tolerance",
        if max_rel < 1e-8 {
            "BELOW — CellList is safe"
        } else {
            "ABOVE — need canonical sort"
        }
    );

    assert!(
        max_rel < 1e-8,
        "CellList/Standard margin {:.3e} exceeds 1e-8 at step {} \
         (ΔE_abs={:.3e} kJ/mol). \
         Apply canonical-order sort to CellList output per item 9a-1.",
        max_rel,
        max_step,
        max_abs
    );
}
