//! End-to-end integration tests against GROMOS double-precision references.
//!
//! Runs the `md` binary on each reference system (10 steps) and compares its output with the
//! native gromosXX's: every value of every `ENERGY03` and `VOLUMEPRESSURE03` block of the `.tre`
//! and, where the system has a reference `.trg`, of every `FREEENERDERIVS03` block; forces
//! (free and constraint), trajectory and final configuration.
//!
//! Reference data: `tests/gromosXX_references/{system}/expected/`
//! Systems matching GROMOS are active; known mismatches are `#[ignore]`.
//!
//! Run passing:  cargo test -p newen-md --test test_gromosXX_references
//! Run all:      cargo test -p newen-md --test test_gromosXX_references -- --include-ignored

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─── Tolerances ─────────────────────────────────────────────────────────────

const ENERGY_REL_TOL: f64 = 1e-8;
const ENERGY_ABS_TOL: f64 = 1e-10; // for near-zero energies
const FORCE_ABS_TOL: f64 = 1e-6; // kJ/(mol*nm)
const DHDL_REL_TOL: f64 = 1e-6; // dH/dλ relative tolerance
/// The final configuration after the whole run: the frame-by-frame agreement is 1e-8 relative in
/// the energies, and the positions accumulate that over every step, so this is looser than the
/// per-step position tolerance.
const FINAL_POSITION_ABS_TOL: f64 = 1e-6; // nm
/// Trajectory frames are written with nine decimals, so agreement is exact up to that rounding.
const TRAJECTORY_POSITION_ABS_TOL: f64 = 1e-8; // nm
const FINAL_VELOCITY_ABS_TOL: f64 = 1e-4; // nm/ps

// ─── Paths ──────────────────────────────────────────────────────────────────

fn ref_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

fn md_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_md"))
}

// ─── Minimal TOML field extraction ──────────────────────────────────────────

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

/// The `"key" = "value"` lines of the optional `[recipe_overrides]` table: recipe fields no
/// `.imd` block carries, passed as `md @recipe_set key=value` (`RunRecipe::set`).
fn recipe_overrides(content: &str) -> Vec<(String, String)> {
    let mut lines = content.lines().map(str::trim);
    if !lines.any(|l| l == "[recipe_overrides]") {
        return Vec::new();
    }
    lines
        .take_while(|l| !l.starts_with('['))
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            let (k, v) = l.split_once('=').expect("recipe_overrides: key = value");
            let unquote = |s: &str| s.trim().trim_matches('"').to_string();
            (unquote(k), unquote(v))
        })
        .collect()
}

fn toml_str_opt(content: &str, key: &str) -> Option<String> {
    let prefix = format!("{key} = ");
    for line in content.lines() {
        let t = line.trim();
        if t.starts_with(&prefix) {
            if let (Some(a), Some(b)) = (t.find('"'), t.rfind('"')) {
                if a < b {
                    return Some(t[a + 1..b].to_string());
                }
            }
        }
    }
    None
}

// ─── Parser: ENERGY03 / FREEENERDERIVS03 blocks, section by section ────────
//
// Both blocks are gromosXX's `_print_energyred_helper` body (FREEENERDERIVS03 first carries
// `# lambda`): `# totals` (52 values), `# baths` (count, then kinetic / COM / internal-rotational
// per bath), `# bonded` (count, then one row per energy group), `# nonbonded` (one row per group
// pair), `# special` (one row per group), then `# eds`, `# gamd`, `# precalclam` (each a count
// and rows) and `# ABdih`. Every value is compared; count lines exactly.

/// Sections whose first line is an integer count (`# temperature`: the number of baths, then
/// per bath total, COM and internal-rotational temperature and the scaling factor).
const COUNTED_SECTIONS: [&str; 6] = [
    "baths",
    "bonded",
    "eds",
    "gamd",
    "precalclam",
    "temperature",
];
/// Comment words that open a section; any other comment line (column headers, `# numstates`)
/// stays in the current one. The last four are `VOLUMEPRESSURE03`'s (gromosXX
/// `_print_volumepressurered_helper`): `# pressure` holds the pressure, virial and kinetic
/// energy (trace / 3) and then the pressure, virial and kinetic-energy tensors.
const SECTIONS: [&str; 14] = [
    "lambda",
    "totals",
    "baths",
    "bonded",
    "nonbonded",
    "special",
    "eds",
    "gamd",
    "precalclam",
    "ABdih",
    "mass",
    "temperature",
    "volume",
    "pressure",
];

/// One energy block: its time (from the preceding `TIMESTEP`) and its sections in file order,
/// each a list of rows.
#[derive(Debug, Default)]
struct EnergyBlock {
    time: f64,
    sections: Vec<(String, Vec<Vec<f64>>)>,
}

fn parse_energy_blocks(path: &Path, block: &str) -> Vec<EnergyBlock> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut frames = Vec::new();
    let mut current: Option<EnergyBlock> = None;
    let mut in_timestep = false;
    let mut time = 0.0;
    for line in content.lines() {
        let t = line.trim();
        if t == "TIMESTEP" {
            in_timestep = true;
            continue;
        }
        if t == block {
            current = Some(EnergyBlock {
                time,
                ..Default::default()
            });
            continue;
        }
        if t == "END" {
            in_timestep = false;
            if let Some(f) = current.take() {
                frames.push(f);
            }
            continue;
        }
        if in_timestep {
            if let Some(v) = t.split_whitespace().nth(1).and_then(|x| x.parse().ok()) {
                time = v;
            }
            continue;
        }
        let Some(f) = current.as_mut() else { continue };
        if t.is_empty() {
            continue;
        }
        if let Some(comment) = t.strip_prefix('#') {
            let word = comment.split_whitespace().next().unwrap_or("");
            if SECTIONS.contains(&word) {
                f.sections.push((word.to_string(), Vec::new()));
            }
            continue;
        }
        let row: Vec<f64> = t
            .split_whitespace()
            .map(|x| {
                x.parse()
                    .unwrap_or_else(|_| panic!("{}: not a number {x:?} in {block}", path.display()))
            })
            .collect();
        match f.sections.last_mut() {
            Some((_, rows)) => rows.push(row),
            None => panic!(
                "{}: values before the first section of {block}",
                path.display()
            ),
        }
    }
    frames
}

/// One value that differs from gromosXX.
struct Mismatch {
    frame: usize,
    block: &'static str,
    section: String,
    row: usize,
    col: usize,
    expected: f64,
    actual: f64,
}

/// Compare two blocks section by section, value by value; `rel_tol` relative (with
/// `ENERGY_ABS_TOL` as the floor), count lines exactly. Structural differences (a missing
/// section, a different number of rows or columns) are reported as mismatches with NaN.
fn compare_energy_block(
    frame: usize,
    block: &'static str,
    expected: &EnergyBlock,
    actual: &EnergyBlock,
    rel_tol: f64,
    out: &mut Vec<Mismatch>,
) {
    let mut push = |section: &str, row: usize, col: usize, e: f64, a: f64| {
        out.push(Mismatch {
            frame,
            block,
            section: section.to_string(),
            row,
            col,
            expected: e,
            actual: a,
        })
    };
    for (k, (name, exp_rows)) in expected.sections.iter().enumerate() {
        let Some((act_name, act_rows)) = actual.sections.get(k) else {
            push(&format!("{name} (missing)"), 0, 0, f64::NAN, f64::NAN);
            continue;
        };
        if act_name != name {
            push(
                &format!("{name} (found {act_name})"),
                0,
                0,
                f64::NAN,
                f64::NAN,
            );
            continue;
        }
        if exp_rows.len() != act_rows.len() {
            push(
                &format!("{name} (row count)"),
                0,
                0,
                exp_rows.len() as f64,
                act_rows.len() as f64,
            );
        }
        let counted = COUNTED_SECTIONS.contains(&name.as_str());
        for (r, (er, ar)) in exp_rows.iter().zip(act_rows).enumerate() {
            if er.len() != ar.len() {
                push(name, r, 0, er.len() as f64, ar.len() as f64);
                continue;
            }
            for (c, (&e, &a)) in er.iter().zip(ar).enumerate() {
                let ok = if counted && r == 0 {
                    e == a
                } else {
                    (e - a).abs() <= (e.abs() * rel_tol).max(ENERGY_ABS_TOL)
                };
                if !ok {
                    push(name, r, c, e, a);
                }
            }
        }
    }
    if actual.sections.len() > expected.sections.len() {
        for (name, _) in &actual.sections[expected.sections.len()..] {
            push(&format!("{name} (extra)"), 0, 0, f64::NAN, f64::NAN);
        }
    }
}

fn report_mismatches(system: &str, mismatches: &[Mismatch]) {
    if mismatches.is_empty() {
        return;
    }
    let mut msg = format!(
        "{system}: {} value(s) differ from gromosXX (frame, block, section, row, column: expected vs actual):\n",
        mismatches.len()
    );
    for m in mismatches.iter().take(40) {
        msg.push_str(&format!(
            "  [{}] {} {} row {} col {}: {:.10e} vs {:.10e}\n",
            m.frame, m.block, m.section, m.row, m.col, m.expected, m.actual
        ));
    }
    if mismatches.len() > 40 {
        msg.push_str(&format!("  … and {} more\n", mismatches.len() - 40));
    }
    let mut by_section: Vec<(String, usize)> = Vec::new();
    for m in mismatches {
        let key = format!("{} {}", m.block, m.section);
        match by_section.iter_mut().find(|(k, _)| *k == key) {
            Some((_, n)) => *n += 1,
            None => by_section.push((key, 1)),
        }
    }
    msg.push_str("  by section:");
    for (k, n) in &by_section {
        msg.push_str(&format!(" [{k}: {n}]"));
    }
    panic!("{msg}");
}

// ─── Force frame ────────────────────────────────────────────────────────────

/// One frame of per-atom forces: Vec of (fx, fy, fz) for each atom.
type ForceFrame = Vec<[f64; 3]>;

// ─── Parser: FREEFORCERED / CONSFORCERED blocks from .trf files ─────────────
//
// Both GROMOS and newen write these blocks with 3 floats per line.
// Lines starting with '#' are comments (atom count markers).

fn parse_trf(path: &Path, block: &str) -> Vec<ForceFrame> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut frames = Vec::new();
    let mut in_block = false;
    let mut current_frame: ForceFrame = Vec::new();

    for line in content.lines() {
        let t = line.trim();
        if t == block {
            in_block = true;
            current_frame.clear();
            continue;
        }
        if t == "END" && in_block {
            frames.push(current_frame.clone());
            in_block = false;
            continue;
        }
        if in_block && !t.starts_with('#') && !t.is_empty() {
            let vals: Vec<f64> = t
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if vals.len() == 3 {
                current_frame.push([vals[0], vals[1], vals[2]]);
            }
        }
    }
    frames
}

fn assert_force_close(actual: f64, expected: f64, label: &str) {
    let diff = (actual - expected).abs();
    assert!(
        diff <= FORCE_ABS_TOL,
        "{label}: expected {expected:.10e}, got {actual:.10e}, diff={diff:.2e}, tol={FORCE_ABS_TOL:.2e}"
    );
}

// ─── Parser: GROMOS .conf position file ─────────────────────────────────────
//
// Reads the POSITION block. Each data line: <label x4> <seq> x y z  (nm).
// The first 24 chars are a label field that is ignored by gromosXX convention.

// ─── Test driver ────────────────────────────────────────────────────────────

fn run_reference(system: &str) {
    let sys_dir = ref_root().join(system);
    if !sys_dir.exists() {
        eprintln!(
            "SKIP {system}: reference data not found at {}",
            sys_dir.display()
        );
        return;
    }

    // Read input file paths from input.toml
    let toml = fs::read_to_string(sys_dir.join("input.toml")).expect("input.toml missing");
    let topo = sys_dir.join(toml_str(&toml, "topology"));
    let conf = sys_dir.join(toml_str(&toml, "configuration"));
    let params = sys_dir.join(toml_str(&toml, "parameters"));

    // Output directory
    let out = std::env::temp_dir().join(format!("newen_reftest_{}_{}", system, std::process::id()));
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(&out).unwrap();

    let tre = out.join("energies.tre");
    let trf = out.join("forces.trf");
    let trg = out.join("free_energy.trg");

    // Run md binary
    let mut cmd = Command::new(md_bin());
    cmd.arg("@topo")
        .arg(&topo)
        .arg("@conf")
        .arg(&conf)
        .arg("@input")
        .arg(&params)
        .arg("@fin")
        .arg(out.join("final.conf"))
        .arg("@tre")
        .arg(&tre)
        .arg("@trf")
        .arg(&trf)
        .arg("@trc")
        .arg(out.join("trajectory.trc"))
        .arg("@trv")
        .arg(out.join("velocities.trv"))
        .arg("@trg")
        .arg(&trg);

    // Optional position restraints
    if let Some(por) = toml_str_opt(&toml, "posresspec") {
        cmd.arg("@posresspec").arg(sys_dir.join(por));
    }
    if let Some(rpr) = toml_str_opt(&toml, "refpos") {
        cmd.arg("@refpos").arg(sys_dir.join(rpr));
    }

    // Optional distance restraints
    if let Some(dr) = toml_str_opt(&toml, "distrest") {
        cmd.arg("@distrest").arg(sys_dir.join(dr));
    }

    // Optional GaMD specification
    if let Some(g) = toml_str_opt(&toml, "gamd") {
        cmd.arg("@gamd").arg(sys_dir.join(g));
    }

    // Optional perturbation topology
    if let Some(pt) = toml_str_opt(&toml, "pttopo") {
        cmd.arg("@pttopo").arg(sys_dir.join(pt));
    }

    for (key, value) in recipe_overrides(&toml) {
        cmd.arg("@recipe_set").arg(format!("{key}={value}"));
    }

    let t0 = std::time::Instant::now();
    let result = cmd.output().expect("failed to execute md");
    let elapsed = t0.elapsed();

    eprintln!("  {system}: {:.3}ms", elapsed.as_secs_f64() * 1000.0);

    assert!(
        result.status.success(),
        "{system}: md exited with {}\nstdout: {}\nstderr: {}",
        result.status,
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );

    // Energies: every value of every ENERGY03 block, and of every FREEENERDERIVS03 block where
    // the system has a reference .trg.
    let expected = parse_energy_blocks(&sys_dir.join("expected/energies.tre"), "ENERGY03");
    let actual = parse_energy_blocks(&tre, "ENERGY03");

    // A reference whose forces come from a defective gromosXX branch (EDS FORM 2/3 weights,
    // patches/D2_eds_weights.patch): only what does not depend on those forces is compared —
    // the end-state energies of frame 0 and, where gromosXX's V_R is right, slot 35.
    if let Some(scope) = toml_str_opt(&toml, "compare") {
        let with_vr = match scope.as_str() {
            "frame0_eds_states" => false,
            "frame0_eds_states_vr" => true,
            other => panic!("{system}: unknown compare scope {other:?}"),
        };
        let section = |b: &EnergyBlock, name: &str| -> Vec<Vec<f64>> {
            b.sections
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, rows)| rows.clone())
                .unwrap_or_else(|| panic!("{system}: no `# {name}` section"))
        };
        let (e0, a0) = (&expected[0], &actual[0]);
        let mut mismatches = Vec::new();
        let (exp_eds, act_eds) = (section(e0, "eds"), section(a0, "eds"));
        assert_eq!(exp_eds[0], act_eds[0], "{system}: number of end states");
        assert_eq!(exp_eds.len(), act_eds.len(), "{system}: end-state rows");
        for (r, (er, ar)) in exp_eds.iter().zip(&act_eds).enumerate().skip(1) {
            for (c, (&e, &a)) in er.iter().zip(ar).enumerate() {
                if (e - a).abs() > (e.abs() * ENERGY_REL_TOL).max(ENERGY_ABS_TOL) {
                    mismatches.push(Mismatch {
                        frame: 0,
                        block: "ENERGY03",
                        section: "eds".to_string(),
                        row: r,
                        col: c,
                        expected: e,
                        actual: a,
                    });
                }
            }
        }
        if with_vr {
            let (e, a) = (section(e0, "totals")[34][0], section(a0, "totals")[34][0]);
            if (e - a).abs() > (e.abs() * ENERGY_REL_TOL).max(ENERGY_ABS_TOL) {
                mismatches.push(Mismatch {
                    frame: 0,
                    block: "ENERGY03",
                    section: "totals (slot 35, eds_vr)".to_string(),
                    row: 34,
                    col: 0,
                    expected: e,
                    actual: a,
                });
            }
        }
        report_mismatches(system, &mismatches);
        let _ = fs::remove_dir_all(&out);
        return;
    }

    // Both engines run exactly NSTLIM steps and write one frame per step, so the counts must
    // match. (They did not until 0.0.42: our loop ran NSTLIM+1 steps, and this assertion used to
    // allow the extra frame — which is precisely why the step count went unnoticed.)
    assert_eq!(
        actual.len(),
        expected.len(),
        "{system}: frame count differs from gromosXX (expected {}, got {})",
        expected.len(),
        actual.len()
    );

    let mut mismatches = Vec::new();
    for (i, (exp, act)) in expected.iter().zip(&actual).enumerate() {
        compare_energy_block(i, "ENERGY03", exp, act, ENERGY_REL_TOL, &mut mismatches);
    }
    let expected_vp =
        parse_energy_blocks(&sys_dir.join("expected/energies.tre"), "VOLUMEPRESSURE03");
    let actual_vp = parse_energy_blocks(&tre, "VOLUMEPRESSURE03");
    assert_eq!(
        actual_vp.len(),
        expected_vp.len(),
        "{system}: VOLUMEPRESSURE03 frame count differs from gromosXX"
    );
    for (i, (exp, act)) in expected_vp.iter().zip(&actual_vp).enumerate() {
        compare_energy_block(
            i,
            "VOLUMEPRESSURE03",
            exp,
            act,
            ENERGY_REL_TOL,
            &mut mismatches,
        );
    }
    if let Some(fe_path) = toml_str_opt(&toml, "free_energy") {
        let expected_trg = sys_dir.join(fe_path);
        if expected_trg.exists() {
            let exp_fe = parse_energy_blocks(&expected_trg, "FREEENERDERIVS03");
            let act_fe = if trg.exists() {
                parse_energy_blocks(&trg, "FREEENERDERIVS03")
            } else {
                Vec::new()
            };
            // Both engines write steps 0 … NSTLIM−1 every NTWG steps, so the frames pair one to
            // one; a reference written with NTWG=0 (empty) needs an empty or absent .trg.
            assert_eq!(
                act_fe.len(),
                exp_fe.len(),
                "{system}: FREEENERDERIVS03 frame count differs from gromosXX (expected {}, got {})",
                exp_fe.len(),
                act_fe.len()
            );
            for (i, (ef, af)) in exp_fe.iter().zip(&act_fe).enumerate() {
                if (af.time - ef.time).abs() > 1e-9 {
                    mismatches.push(Mismatch {
                        frame: i,
                        block: "FREEENERDERIVS03",
                        section: "(TIMESTEP time)".to_string(),
                        row: 0,
                        col: 0,
                        expected: ef.time,
                        actual: af.time,
                    });
                    continue;
                }
                compare_energy_block(i, "FREEENERDERIVS03", ef, af, DHDL_REL_TOL, &mut mismatches);
            }
            // the engine's file also reads back through the crate's own reader
            if trg.exists() {
                let read = newen_io::read_free_energy_trajectory(&trg)
                    .unwrap_or_else(|e| panic!("{system}: {}: {e}", trg.display()));
                assert_eq!(
                    read.frames.len(),
                    act_fe.len(),
                    "{system}: newen-io reads a different number of .trg frames"
                );
            }
        }
    }
    report_mismatches(system, &mismatches);

    // Every frame's TIMESTEP (step and time, t = T + step·DT) in each file gromosXX wrote: the
    // energies, dH/dλ, forces, positions and the final configuration — the frame count of each
    // with it.
    for (file, ours) in [
        ("expected/energies.tre".to_string(), tre.clone()),
        ("expected/forces.trf".to_string(), trf.clone()),
        (
            "expected/trajectory.trc".to_string(),
            out.join("trajectory.trc"),
        ),
        ("expected/final.conf".to_string(), out.join("final.conf")),
    ]
    .into_iter()
    .chain(toml_str_opt(&toml, "free_energy").map(|fe| (fe, trg.clone())))
    {
        let theirs = sys_dir.join(&file);
        if !theirs.exists() {
            continue;
        }
        let expected = parse_timesteps(&theirs);
        let actual = if ours.exists() {
            parse_timesteps(&ours)
        } else {
            Vec::new()
        };
        assert_eq!(
            actual.len(),
            expected.len(),
            "{system}: {file}: {} TIMESTEP blocks, gromosXX wrote {}",
            actual.len(),
            expected.len()
        );
        for (i, (&(es, et), &(a_s, at))) in expected.iter().zip(&actual).enumerate() {
            assert!(
                es == a_s && (et - at).abs() <= 1e-9,
                "{system}: {file} frame {i}: TIMESTEP {a_s} {at} vs gromosXX {es} {et}"
            );
        }
    }

    // Compare forces if reference .trf exists: the free forces, and the constraint forces
    // gromosXX writes when SHAKE runs
    let expected_trf = sys_dir.join("expected/forces.trf");
    for block in ["FREEFORCERED", "CONSFORCERED"] {
        if !(expected_trf.exists() && trf.exists()) {
            break;
        }
        let exp_forces = parse_trf(&expected_trf, block);
        let act_forces = parse_trf(&trf, block);
        if exp_forces.is_empty() != act_forces.is_empty() {
            panic!(
                "{system}: {block} written by gromosXX: {}, by md: {}",
                !exp_forces.is_empty(),
                !act_forces.is_empty()
            );
        }

        // a converged minimisation repeats its last frame (`out_configuration.cc:624-630`)
        assert_eq!(
            act_forces.len(),
            exp_forces.len(),
            "{system}: {block} frame count differs from gromosXX (expected {}, got {})",
            exp_forces.len(),
            act_forces.len()
        );
        for (frame_idx, (ef, af)) in exp_forces.iter().zip(&act_forces).enumerate() {
            assert_eq!(
                ef.len(),
                af.len(),
                "{system}[frame {frame_idx}]: atom count mismatch (expected {}, got {})",
                ef.len(),
                af.len()
            );
            for (atom_idx, (e, a)) in ef.iter().zip(af.iter()).enumerate() {
                assert_force_close(
                    a[0],
                    e[0],
                    &format!("{system}[{frame_idx}] {block} atom {atom_idx} fx"),
                );
                assert_force_close(
                    a[1],
                    e[1],
                    &format!("{system}[{frame_idx}] {block} atom {atom_idx} fy"),
                );
                assert_force_close(
                    a[2],
                    e[2],
                    &format!("{system}[{frame_idx}] {block} atom {atom_idx} fz"),
                );
            }
        }
    }

    // ── Trajectory ───────────────────────────────────────────────────────────
    // A frame's structure must belong with the energies of the same frame: gromosXX writes both at
    // the same point of the step, so its frame 0 is the input configuration. Ours wrote the state
    // *after* the step until 0.0.43, which left the two outputs one step out of step.
    let expected_trc = sys_dir.join("expected/trajectory.trc");
    if expected_trc.exists() {
        let (ours, _) = parse_configuration(&out.join("trajectory.trc"));
        let (theirs, _) = parse_configuration(&expected_trc);

        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: trajectory position count differs from gromosXX"
        );
        for (i, (x, y)) in ours.iter().zip(&theirs).enumerate() {
            let d = fold_into_box(&expected_trc, [x[0] - y[0], x[1] - y[1], x[2] - y[2]]);
            for (c, dc) in d.iter().enumerate() {
                assert!(
                    dc.abs() <= TRAJECTORY_POSITION_ABS_TOL,
                    "{system}: trajectory position {} component {c}: {:.10e} vs gromosXX {:.10e}",
                    i + 1,
                    x[c],
                    y[c]
                );
            }
        }
    }

    // ── Velocity trajectory (@trv / NTWV) ────────────────────────────────────
    let expected_trv = sys_dir.join("expected/velocities.trv");
    if expected_trv.exists() {
        let (_, ours) = parse_configuration(&out.join("velocities.trv"));
        let (_, theirs) = parse_configuration(&expected_trv);
        // `@trv` was accepted and ignored until 0.0.44; the count is what proves it is written.
        // The *values* are not compared yet: velocities are half-step quantities in leap-frog and
        // gromosXX's frame 0 is the input velocities, while ours is one half-step further on — the
        // residual from frame 1 onwards is ~1e-4 nm/ps, the size of one thermostat scaling. Which
        // half-step each frame carries is an open item (1.5b).
        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: velocity trajectory count differs from gromosXX (is @trv written at all?)"
        );
        assert!(
            !ours.is_empty(),
            "{system}: @trv produced no velocity frames"
        );
    }

    // ── Final configuration ──────────────────────────────────────────────────
    // The end state is what a continuation run starts from, and it is the only place a wrong
    // step count shows up — the per-step frames look right either way.
    let expected_fin = sys_dir.join("expected/final.conf");
    if expected_fin.exists() {
        let ours = parse_configuration(&out.join("final.conf"));
        let theirs = parse_configuration(&expected_fin);
        // Which periodic image a configuration file records is not physics: gromosXX puts charge
        // groups back into the box at the start of each step, so an atom that the last position
        // update pushed just outside is written on one side by one engine and on the other by the
        // other. Differences are compared modulo a box vector; anything else is a real deviation.

        for (what, a, b) in [
            ("position", &ours.0, &theirs.0),
            ("velocity", &ours.1, &theirs.1),
        ] {
            if b.is_empty() {
                continue;
            }
            assert_eq!(
                a.len(),
                b.len(),
                "{system}: final {what} count differs from gromosXX"
            );
            for (i, (x, y)) in a.iter().zip(b).enumerate() {
                for (c, (xc, yc)) in x.iter().zip(y).enumerate() {
                    let tol = if what == "position" {
                        FINAL_POSITION_ABS_TOL
                    } else {
                        FINAL_VELOCITY_ABS_TOL
                    };
                    let mut diff = xc - yc;
                    if what == "position" {
                        diff =
                            fold_into_box(&expected_fin, [x[0] - y[0], x[1] - y[1], x[2] - y[2]])
                                [c];
                    }
                    assert!(
                        diff.abs() <= tol,
                        "{system}: final {what} of atom {} component {c}: {xc:.10e} vs gromosXX {yc:.10e} (tol {tol:.1e})",
                        i + 1
                    );
                }
            }
        }
    }

    // A perturbed run ends its final configuration with PERTDATA (λ, `%15.7f`), the start of a
    // run with NRDGL=1; an unperturbed one has none.
    if expected_fin.exists() {
        let ours = parse_pertdata(&out.join("final.conf"));
        let theirs = parse_pertdata(&expected_fin);
        assert_eq!(
            ours.is_some(),
            theirs.is_some(),
            "{system}: final configuration PERTDATA present in ours: {}, in gromosXX's: {}",
            ours.is_some(),
            theirs.is_some()
        );
        if let (Some(a), Some(b)) = (ours, theirs) {
            assert!(
                (a - b).abs() <= 5e-8,
                "{system}: final configuration PERTDATA {a} vs gromosXX {b}"
            );
        }
        // An A-EDS search (FORMS 2-6): the AEDSSEARCH state a continuation starts from, printed
        // with 9 decimals (`out_configuration.cc:3985-4004`).
        let ours = parse_block_values(&out.join("final.conf"), "AEDSSEARCH");
        let theirs = parse_block_values(&expected_fin, "AEDSSEARCH");
        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: final configuration AEDSSEARCH has {} values, gromosXX's {}",
            ours.len(),
            theirs.len()
        );
        for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            assert!(
                (a - b).abs() <= 1e-8 * b.abs().max(1.0),
                "{system}: AEDSSEARCH value {}: {a} vs gromosXX {b}",
                i + 1
            );
        }
        // A Nosé–Hoover run: the ζ of every bath a continuation starts from (NHCVARIABLES, 9
        // decimals, `out_configuration.cc:3426-3443`); gromosXX writes it for chains only (D41).
        let ours = parse_block_values(&out.join("final.conf"), "NHCVARIABLES");
        let theirs = parse_block_values(&expected_fin, "NHCVARIABLES");
        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: final configuration NHCVARIABLES has {} values, gromosXX's {}",
            ours.len(),
            theirs.len()
        );
        for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            assert!(
                (a - b).abs() <= 1e-8 * b.abs().max(1.0),
                "{system}: NHCVARIABLES value {}: {a} vs gromosXX {b}",
                i + 1
            );
        }
        // A GaMD search (SEARCH 1/2): the statistics a continuation starts from, printed with 9
        // decimals (`out_configuration.cc:4013-4034`).
        let ours = parse_block_values(&out.join("final.conf"), "GAMDSTAT");
        let theirs = parse_block_values(&expected_fin, "GAMDSTAT");
        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: final configuration GAMDSTAT has {} values, gromosXX's {}",
            ours.len(),
            theirs.len()
        );
        for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            assert!(
                (a - b).abs() <= 1e-8 * b.abs().max(1.0),
                "{system}: GAMDSTAT value {}: {a} vs gromosXX {b}",
                i + 1
            );
        }
        // Time-averaged distance restraints (NTDIR < 0): the ⟨r⁻³⟩ a continuation starts from.
        let ours = parse_block_values(&out.join("final.conf"), "DISRESEXPAVE");
        let theirs = parse_block_values(&expected_fin, "DISRESEXPAVE");
        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: final configuration DISRESEXPAVE has {} values, gromosXX's {}",
            ours.len(),
            theirs.len()
        );
        for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            assert!(
                (a - b).abs() <= 1e-5 * b.abs().max(1.0),
                "{system}: DISRESEXPAVE value {}: {a} vs gromosXX {b}",
                i + 1
            );
        }
    }

    let _ = fs::remove_dir_all(&out);
}

/// The (step, time) of every `TIMESTEP` block of a file, in file order.
fn parse_timesteps(path: &Path) -> Vec<(u64, f64)> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.lines().map(str::trim);
    let mut out = Vec::new();
    while let Some(l) = lines.next() {
        if l != "TIMESTEP" {
            continue;
        }
        let row = lines
            .find(|l| !l.is_empty() && !l.starts_with('#'))
            .unwrap_or_else(|| panic!("{}: empty TIMESTEP block", path.display()));
        let mut v = row.split_whitespace();
        let bad = format!("{}: bad TIMESTEP line {row:?}", path.display());
        let step: u64 = v.next().and_then(|x| x.parse().ok()).expect(&bad);
        let time: f64 = v.next().and_then(|x| x.parse().ok()).expect(&bad);
        out.push((step, time));
    }
    out
}

/// Every number of a `.cnf` block (empty if the block is absent).
fn parse_block_values(path: &Path, block: &str) -> Vec<f64> {
    let Ok(text) = fs::read_to_string(path) else {
        return Vec::new();
    };
    text.lines()
        .skip_while(|l| l.trim() != block)
        .skip(1)
        .map(str::trim)
        .take_while(|l| *l != "END")
        .filter(|l| !l.starts_with('#'))
        .flat_map(|l| l.split_whitespace().map(|x| x.parse::<f64>().unwrap()))
        .collect()
}

/// The λ of a `.cnf`'s PERTDATA block, if it has one.
fn parse_pertdata(path: &Path) -> Option<f64> {
    let text = fs::read_to_string(path).ok()?;
    text.lines()
        .skip_while(|l| l.trim() != "PERTDATA")
        .skip(1)
        .map(str::trim)
        .find(|l| !l.is_empty() && !l.starts_with('#'))
        .and_then(|l| l.split_whitespace().next()?.parse().ok())
}

/// Fold a coordinate difference into the periodic images of a reference `.cnf`'s box, so that a
/// comparison does not depend on which image an engine happened to write. Configurations are
/// written in the frame of the input file, so NTB = −1 is folded with the truncated-octahedron
/// rule *in the cube frame* (the cubic minimum image, then the body-centred (±L/2)³ shift), not
/// with the triclinic lattice the engine works in.
fn fold_into_box(path: &Path, d: [f64; 3]) -> [f64; 3] {
    let Some((ntb, edges, angles)) = parse_box(path) else {
        return d;
    };
    if edges[0] <= 0.0 || ntb == 0 {
        return d;
    }
    if ntb == 2 {
        // a triclinic cell: fold the fractional coordinates
        use newen_core::math::{genbox_box, Vec3};
        let h = genbox_box(Vec3::from(edges), Vec3::from(angles));
        let f = h.inverse() * Vec3::from(d);
        return (h * (f - f.round())).to_array();
    }
    let mut r = [0.0; 3];
    for c in 0..3 {
        let l = edges[c];
        r[c] = if l > 0.0 {
            d[c] - l * (d[c] / l).round()
        } else {
            d[c]
        };
    }
    if ntb == -1 {
        let l = edges[0];
        if r[0].abs() + r[1].abs() + r[2].abs() > 0.75 * l {
            for (c, rc) in r.iter_mut().enumerate() {
                let _ = c;
                *rc -= 0.5 * l * rc.signum();
            }
        }
    }
    r
}

/// NTB, the box edge lengths and the angles of a `.cnf`'s GENBOX block, if it has one.
fn parse_box(path: &Path) -> Option<(i32, [f64; 3], [f64; 3])> {
    let text = fs::read_to_string(path).ok()?;
    let mut lines = text
        .lines()
        .skip_while(|l| l.trim() != "GENBOX")
        .skip(1)
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let ntb: i32 = lines.next()?.split_whitespace().next()?.parse().ok()?;
    let mut row = || -> Option<[f64; 3]> {
        let v: Vec<f64> = lines
            .next()?
            .split_whitespace()
            .filter_map(|x| x.parse().ok())
            .collect();
        (v.len() >= 3).then_some([v[0], v[1], v[2]])
    };
    let edges = row()?;
    let angles = row().unwrap_or([90.0; 3]);
    Some((ntb, edges, angles))
}

/// POSITION and VELOCITY blocks of a `.cnf`, as `[x, y, z]` per atom.
fn parse_configuration(path: &Path) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let (mut pos, mut vel) = (Vec::new(), Vec::new());
    let mut block: Option<&str> = None;
    for line in text.lines() {
        let t = line.trim();
        if t == "POSITION" || t == "POSITIONRED" {
            block = Some("pos");
            continue;
        }
        if t == "VELOCITY" || t == "VELOCITYRED" {
            block = Some("vel");
            continue;
        }
        if t == "END" {
            block = None;
            continue;
        }
        let Some(kind) = block else { continue };
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let v: Vec<f64> = t
            .split_whitespace()
            .filter_map(|x| x.parse::<f64>().ok())
            .collect();
        if v.len() >= 3 {
            let xyz = [v[v.len() - 3], v[v.len() - 2], v[v.len() - 1]];
            if kind == "pos" {
                pos.push(xyz);
            } else {
                vel.push(xyz);
            }
        }
    }
    (pos, vel)
}

// ─── Test declarations ──────────────────────────────────────────────────────

macro_rules! ref_test {
    ($name:ident, $sys:literal) => {
        #[test]
        fn $name() {
            run_reference($sys);
        }
    };
    (ignore: $name:ident, $sys:literal) => {
        #[test]
        #[ignore = "known mismatch — see "]
        fn $name() {
            run_reference($sys);
        }
    };
    (ignore: $name:ident, $sys:literal, $reason:literal) => {
        #[test]
        #[ignore = $reason]
        fn $name() {
            run_reference($sys);
        }
    };
}

// ── Level 0: pair interactions, vacuum ──────────────────────────────────────
ref_test!(pair_lj, "pair_lj");
ref_test!(pair_lj_mixed, "pair_lj_mixed");
ref_test!(nacl_pair, "nacl_pair");

// ── Level 1: single molecule, bonded terms, PBC+RF ──────────────────────────
ref_test!(water_single, "water_single");
ref_test!(water_single_genvel, "water_single_genvel");
ref_test!(benzene_vacuum, "benzene_vacuum");
ref_test!(nacl_pair_box, "nacl_pair_box");
// The direct Ewald sum (NLRELE 2, atomic cutoff, Gaussian a = 0.4 nm, cubic A2); expected
// output from the released binary run with @develop (see its input.toml).
ref_test!(nacl_pair_box_ewald, "nacl_pair_box_ewald");
// P3M (NLRELE 3): NSHAPE 3, ASHAPE 1.4 nm, grid 32³, NASORD 3, NFDORD 2, NALIAS 3, P3M Ã2
// (NA2CLC 3); expected output from the released binary. The FFTs differ (PhastFT vs FFTW), so
// the agreement is within the tolerances, not byte for byte.
ref_test!(nacl_pair_box_p3m_fixed, "nacl_pair_box_p3m_fixed");
ref_test!(butane_vacuum, "butane_vacuum");
// COVALENTFORM 0 0 1: the Chebyshev-polynomial dihedral (NTBDN = 1).
ref_test!(butane_vacuum_ntbdn1, "butane_vacuum_ntbdn1");

// ── Level 2: PBC, pairlist, SHAKE, solvent, twin-range ──────────────────────
ref_test!(water_3_box, "water_3_box");
ref_test!(nacl_1water_box, "nacl_1water_box");
ref_test!(nacl_3water_box, "nacl_3water_box");
ref_test!(water_3_box_twinrange, "water_3_box_twinrange");
ref_test!(water_10_box, "water_10_box");
ref_test!(nacl_3water_cutoff, "nacl_3water_cutoff");
ref_test!(nacl_water_box, "nacl_water_box");
ref_test!(nacl_water_box_shifted, "nacl_water_box_shifted");
// GENBOX NTB 2 with angles 75/80/65: the box vectors from lengths and angles
ref_test!(nacl_water_box_triclinic, "nacl_water_box_triclinic");
// P3M in the same cell: the wave vectors Σ l_i b_i (patched oracle, LS3)
ref_test!(
    nacl_water_box_triclinic_p3m_fixed,
    "nacl_water_box_triclinic_p3m_fixed"
);

// ── Constraint algorithms ────────────────────────────────────────────────────
ref_test!(nacl_1water_settle, "nacl_1water_settle");
ref_test!(nacl_1water_lincs, "nacl_1water_lincs");
ref_test!(aladip_vacuum_lincs, "aladip_vacuum_lincs");

// ── Level 2: grid_cell pairlist (1000 SPC molecules) ────────────────────────
ref_test!(water_1000_spc_gridcell, "water_1000_spc_gridcell");

// ── Level 3: bulk water ──────────────────────────────────────────────────────
ref_test!(water_216_box, "water_216_box");
// Same system with the MULTIBATH block *absent*: gromosXX runs no temperature coupling. Until
// 3.9 step 2 the parser turned the absent block into a Berendsen bath (A18); the
// defaults now mean "absent", and this guards that they keep meaning it.
ref_test!(water_216_nve_nobath, "water_216_nve_nobath");
ref_test!(water_216_box_com, "water_216_box_com");
ref_test!(water_216_box_com_rot, "water_216_box_com_rot");
// The direct Ewald sum in bulk water: excluded pairs, the twin-range LJ-only shell; and with the
// surface term of a dielectric boundary (EPSLS = 61). Expected output from @develop.
ref_test!(water_216_box_ewald, "water_216_box_ewald");
ref_test!(water_216_box_ewald_epsls, "water_216_box_ewald_epsls");
// The same under Berendsen pressure coupling: the surface-term virial stored as −2Ξ with the whole
// tensor (patched oracle, LS2).
ref_test!(
    water_216_npt_ewald_epsls_fixed,
    "water_216_npt_ewald_epsls_fixed"
);
// P3M in bulk water: the finite-difference field (NFDORD 2) with NSHAPE 3, and the Gaussian
// (a = 0.4 nm) with ik differentiation (NFDORD 0).
ref_test!(water_216_box_p3m_fixed, "water_216_box_p3m_fixed");
ref_test!(
    water_216_box_p3m_gauss_fixed,
    "water_216_box_p3m_gauss_fixed"
);

// ── Level 3: thermostats / barostats ────────────────────────────────────────
ref_test!(water_216_nvt, "water_216_nvt");
// The final configuration's NHCVARIABLES block, plain Nosé–Hoover included (D41).
ref_test!(water_216_nvt_nosehoover, "water_216_nvt_nosehoover");
ref_test!(water_216_nvt_nhc_chain, "water_216_nvt_nhc_chain");
// Steps 100-199 of a 100 + 100 run, continued (NTINHT 0) against gromosXX with D41 patched
// (`patches/D41_nh_restart.patch`): every variable of every bath restored. The continued run is
// checked against the single run by `tests/nh_restart.rs`.
ref_test!(
    water_216_nvt_nosehoover_cont_fixed,
    "water_216_nvt_nosehoover_cont_fixed"
);
ref_test!(water_216_nvt_nhc_cont_fixed, "water_216_nvt_nhc_cont_fixed");
ref_test!(water_216_npt, "water_216_npt");
// Berendsen barostat (VIRIAL 2) with all-bond and solvent SHAKE: the constraint virial enters
// only the atomic virial (gromosXX `shake.h:380`).
ref_test!(aladip_solvated_npt_shake, "aladip_solvated_npt_shake");

// ── Alanine dipeptide ────────────────────────────────────────────────────────
ref_test!(aladip_vacuum, "aladip_vacuum");
ref_test!(aladip_solvated, "aladip_solvated");
ref_test!(aladip_trunc_oct, "aladip_trunc_oct");

// The reaction-field constant (1.8): EPSRF = 0 is a conducting boundary, C_RF = −1; and the RF
// constants come from RCRF (1.2), not from RCUTL (1.4).
ref_test!(aladip_solvated_epsrf0, "aladip_solvated_epsrf0");
ref_test!(aladip_solvated_rcrf, "aladip_solvated_rcrf");

// ── Restraints ───────────────────────────────────────────────────────────────
ref_test!(nacl_1water_distres, "nacl_1water_distres");

// ── FEP / TI ────────────────────────────────────────────────────────────────
ref_test!(ch4_water_fep, "ch4_water_fep");
ref_test!(ch4_water_fep_l000, "ch4_water_fep_l000");
ref_test!(ch4_water_fep_l025, "ch4_water_fep_l025");
ref_test!(ch4_water_fep_l075, "ch4_water_fep_l075");
ref_test!(ch4_water_fep_l100, "ch4_water_fep_l100");
// Perturbed bonded (regular and soft-core), atom types, charges, λ-mixed masses and a
// PERTATOMPAIR entry, all at once — passes since the 2026-08-29/30 FEP fixes.
ref_test!(aladip_vacuum_fep, "aladip_vacuum_fep");
// Charged perturbed atoms (54a7 CH3OH -> dummies): the soft-core reaction-field terms that
// `ch4_water_fep` (zero charge) never exercised. Found and fixed the 1/cutoff³ error in the
// perturbed pair kernel's softened RF term (0.0.34).
ref_test!(meoh_water_fep, "meoh_water_fep");
// Two energy groups (atoms 1-6, 7-12) and a `.trg` every step: every bonded, nonbonded and dH/dλ
// term booked into the energy group (pair) of its first atom(s), per term.
ref_test!(aladip_vacuum_fep_groups, "aladip_vacuum_fep_groups");
// One energy group and a `.trg` every step: the group-0 derivative cells carry the totals.
ref_test!(aladip_vacuum_fep_trg, "aladip_vacuum_fep_trg");
// Individual λ (LAMBDAS, 1.8): vdw/vdw_soft polynomials on the solute-solvent pair (its crf rows
// are inert: CH4 is uncharged in both states); the same block switched off (must equal
// `ch4_water_fep_l025`); and NLAM = 2 on the soft-core kernel.
ref_test!(ch4_water_fep_lambdas, "ch4_water_fep_lambdas");
ref_test!(ch4_water_fep_lambdas_off, "ch4_water_fep_lambdas_off");
ref_test!(ch4_water_fep_nlam2, "ch4_water_fep_nlam2");
// NRDGL = 1: λ comes from the configuration's PERTDATA (0.25, `ch4_water_fep_l025`'s final
// configuration), not from RLAM (0).
ref_test!(ch4_water_fep_nrdgl, "ch4_water_fep_nrdgl");
// Two energy groups: bonded, mass and nonbonded individual λ per group (pair); NLAM = 2.
ref_test!(aladip_vacuum_fep_lambdas, "aladip_vacuum_fep_lambdas");
ref_test!(aladip_vacuum_fep_nlam2, "aladip_vacuum_fep_nlam2");
// NLAM = 2 with the LAMBDAS rows of `aladip_vacuum_fep_lambdas`: the per-kind (1−L)^(n−1)·L′
// factors with L_lj ≠ L_crf.
ref_test!(
    aladip_vacuum_fep_lambdas_nlam2,
    "aladip_vacuum_fep_lambdas_nlam2"
);
// Perturbed SHAKE (NTC = 3): λ-mixed constraint lengths and their dH/dλ in the constraints
// column, with and without individual λ.
ref_test!(aladip_vacuum_fep_shake, "aladip_vacuum_fep_shake");
ref_test!(
    aladip_vacuum_fep_shake_lambdas,
    "aladip_vacuum_fep_shake_lambdas"
);
// PERTBONDSTRETCHH and PERTBONDSTRETCH sharing atom 4: gromosXX shakes the H block first; the
// group-2 bond 8-10 changes length, so its bond λ and constraint dH/dλ are those of group 2.
ref_test!(
    aladip_vacuum_fep_shake_hblock,
    "aladip_vacuum_fep_shake_hblock"
);
// Perturbed distance restraints: energy and dH/dλ per group, with and without individual λ.
ref_test!(aladip_vacuum_fep_distres, "aladip_vacuum_fep_distres");
ref_test!(
    aladip_vacuum_fep_distres_lambdas,
    "aladip_vacuum_fep_distres_lambdas"
);
// NTDIR = −2: the r⁻³ running average (TAUDIR 0.1) in energy, force and dH/dλ, perturbed and
// unperturbed restraints alike, and the DISRESEXPAVE block of the final configuration.
ref_test!(
    aladip_vacuum_fep_distres_tavg,
    "aladip_vacuum_fep_distres_tavg"
);

// Extended TI (PRECALCLAM, 1.8): A/B energies and derivatives at NRLAM extra λ points in the
// `# precalclam` / `# ABdih` sections of the .tre and .trg. CH4: nonbonded LJ only (11 points);
// methanol: the CRF pair, RF self and RF exclusion ext terms; aladip: bonded AB_*, AB_kinetic and
// ABdih; with LAMBDAS (6 points on [0.2, 0.7]): the grid ignores the polynomials except the mass λ
// in AB_kinetic; with perturbed SHAKE: the AB_bond derivatives; with perturbed distance restraints:
// AB_disres.
ref_test!(ch4_water_fep_precalc, "ch4_water_fep_precalc");
ref_test!(meoh_water_fep_precalc, "meoh_water_fep_precalc");
ref_test!(aladip_vacuum_fep_precalc, "aladip_vacuum_fep_precalc");
ref_test!(
    aladip_vacuum_fep_lambdas_precalc,
    "aladip_vacuum_fep_lambdas_precalc"
);
ref_test!(
    aladip_vacuum_fep_shake_precalc,
    "aladip_vacuum_fep_shake_precalc"
);
ref_test!(
    aladip_vacuum_fep_distres_precalc,
    "aladip_vacuum_fep_distres_precalc"
);
// PERTATOMPAIR rows with a perturbed atom, one per state type (excluded / normal / 1-4) in A and
// in B: the ext energies of the filled state only, the other state's parameters set to 0. Against
// gromosXX with the argument shift patched (D1, `patches/D1_pair_kernel_args.patch`).
ref_test!(
    aladip_vacuum_fep_pairpert_precalc_fixed,
    "aladip_vacuum_fep_pairpert_precalc_fixed"
);
// Twin range with the pairlist every 5 steps and output every 2: gromosXX computes the long-range
// perturbed pairs (and their ext energies) only on pairlist steps, the ext part only when that
// step is also an output step, and reuses both until the next update.
ref_test!(
    ch4_water_fep_precalc_twinrange,
    "ch4_water_fep_precalc_twinrange"
);
// A charge change (ΔQ = −1): Na⁺ → neutral dummy in `nacl_water_box` under reaction field. The
// perturbed RF self term and the RF excluded Na–Cl pair move only when the charge does; plain
// TI, PRECALCLAM 11 0 1, and LAMBDAS with crf ≠ vdw rows on (1,1) and (1,2).
ref_test!(nacl_water_fep, "nacl_water_fep");
ref_test!(nacl_water_fep_precalc, "nacl_water_fep_precalc");
ref_test!(nacl_water_fep_lambdas, "nacl_water_fep_lambdas");
// The same with the Na–Cl pair a 1-4 pair (`shared/nacl_water_box_14.topo`) and PRECALCLAM: the
// perturbed 1-4 term, its dH/dλ and ext rows.
ref_test!(nacl_water_fep_14, "nacl_water_fep_14");

// Slow growth (DLAMT ≠ 0): λ += DLAMT·Δt after every step, the masses and every perturbed term
// following it; the .trg carries the λ the step ran at, the .tre mass the masses at the next λ,
// the final PERTDATA the λ after the last step. With perturbed SHAKE; and continued from the
// first run's final configuration (NRDGL = 1), λ running past 1 unbounded.
ref_test!(aladip_vacuum_fep_slowgrowth, "aladip_vacuum_fep_slowgrowth");
ref_test!(
    aladip_vacuum_fep_slowgrowth_shake,
    "aladip_vacuum_fep_slowgrowth_shake"
);
ref_test!(
    aladip_vacuum_fep_slowgrowth_cont,
    "aladip_vacuum_fep_slowgrowth_cont"
);
// Started at STEP T = 10 ps: every TIMESTEP time is T + step·DT, the physics unchanged.
ref_test!(
    aladip_vacuum_fep_slowgrowth_t10,
    "aladip_vacuum_fep_slowgrowth_t10"
);
// In solvent: soft-core pairs at a moving λ in a periodic box with solvent SHAKE.
ref_test!(ch4_water_fep_slowgrowth, "ch4_water_fep_slowgrowth");
// Slow growth under steepest descent: λ moves on every step but the one the minimiser stops at.
ref_test!(aladip_vacuum_em_slowgrowth, "aladip_vacuum_em_slowgrowth");

// ── EDS (1.8) ────────────────────────────────────────────────────────────────
// Enveloping distribution sampling through the recipe: per-state nonbonded terms, the reference
// Hamiltonian V_R and the mixed forces, whole-block. FORM 1 against the released gromosXX.
ref_test!(aladip_vacuum_eds, "aladip_vacuum_eds");
ref_test!(aladip_vacuum_eds3, "aladip_vacuum_eds3");
ref_test!(aladip_eds_solv, "aladip_eds_solv");
ref_test!(two_waters_eds, "two_waters_eds");
// TI-EDS: perturbed atoms next to the EDS atoms, the per-state soft-core pairs and ∂V_i/∂λ.
// newen's EDS atom owns the excluded and 1-4 EDS–perturbed pairs, gromosXX's perturbed atom:
// the same terms summed in another order.
ref_test!(aladip_eds_solv_ti, "aladip_eds_solv_ti");
ref_test!(aladip_eds_solv_ti_groups, "aladip_eds_solv_ti_groups");
// FORMS 2 and 3 against gromosXX with its force weights corrected
// (`tests/gromosXX_references/patches/D2_eds_weights.patch`): the documented formulas.
ref_test!(
    aladip_vacuum_eds_multis_fixed,
    "aladip_vacuum_eds_multis_fixed"
);
ref_test!(
    aladip_vacuum_eds_pairs_fixed,
    "aladip_vacuum_eds_pairs_fixed"
);
ref_test!(
    aladip_eds_solv_ti_multis_fixed,
    "aladip_eds_solv_ti_multis_fixed"
);
// The same FORMS against the released gromosXX, whose force weights are defective: only the
// end-state energies of frame 0 (and, for FORM 3, V_R) — see the `compare` key of their
// input.toml.
ref_test!(aladip_vacuum_eds_multis, "aladip_vacuum_eds_multis");
ref_test!(aladip_vacuum_eds_pairs, "aladip_vacuum_eds_pairs");

// ── Accelerated EDS (1.8) ────────────────────────────────────────────────────
// FORM 1 with V_mix below, inside and above the [E_min, E_max] window: V_R*, the scaled EDS forces.
ref_test!(aladip_vacuum_aeds_below, "aladip_vacuum_aeds_below");
ref_test!(aladip_vacuum_aeds_inside, "aladip_vacuum_aeds_inside");
ref_test!(aladip_vacuum_aeds_above, "aladip_vacuum_aeds_above");
// The parameter searches, whole-block over 300 steps, and their final AEDSSEARCH block, against
// gromosXX with D16, D25 and D32 patched (`tests/gromosXX_references/patches/`): the columns a
// frame's V_R* used, log-sum-exp prevalences, a sampled state as the E_min target. `_sigma0`:
// BMAXTYPE 2 with a target state without variance, held by RESTREMIN 1. `_cont1/2`: a search in
// two jobs, the second continuing the first's AEDSSEARCH block (NTIAEDSS 0). `_f6_ovf`: FORM 6
// with E_2^R = 6000, where e^{−β(V_2 − E_2^R)} overflows.
ref_test!(aladip_vacuum_aeds_f2_fixed, "aladip_vacuum_aeds_f2_fixed");
ref_test!(aladip_vacuum_aeds_f4_fixed, "aladip_vacuum_aeds_f4_fixed");
ref_test!(aladip_vacuum_aeds_f5_fixed, "aladip_vacuum_aeds_f5_fixed");
ref_test!(aladip_vacuum_aeds_f6_fixed, "aladip_vacuum_aeds_f6_fixed");
ref_test!(
    aladip_vacuum_aeds_f6_ovf_fixed,
    "aladip_vacuum_aeds_f6_ovf_fixed"
);
ref_test!(
    aladip_vacuum_aeds_sigma0_fixed,
    "aladip_vacuum_aeds_sigma0_fixed"
);
ref_test!(
    aladip_vacuum_aeds_cont1_fixed,
    "aladip_vacuum_aeds_cont1_fixed"
);
ref_test!(
    aladip_vacuum_aeds_cont2_fixed,
    "aladip_vacuum_aeds_cont2_fixed"
);

// ── Energy minimization ──────────────────────────────────────────────────────
// gromosXX writes the converged energies once more at convergence; `md` does the same since 0.0.34,
// the forces too since 0.0.68. With output every 5 steps: the repeat only on the schedule of the
// step before, the stopping step's positions only on their own.
ref_test!(aladip_vacuum_em, "aladip_vacuum_em");
ref_test!(aladip_vacuum_em_ntw5, "aladip_vacuum_em_ntw5");
ref_test!(aladip_vacuum_em_shake, "aladip_vacuum_em_shake");
ref_test!(aladip_solvated_em_noshake, "aladip_solvated_em_noshake");
ref_test!(aladip_solvated_em_shake, "aladip_solvated_em_shake");
ref_test!(aladip_solvated_em_posres, "aladip_solvated_em_posres");
ref_test!(aladip_solvated_em, "aladip_solvated_em");
// Three energy groups (1-4, 5-12, water): position restraints split over two groups.
ref_test!(
    aladip_solvated_em_posres_groups,
    "aladip_solvated_em_posres_groups"
);

// ── Level 5: the input styles real GROMOS files come in ─────────────────────
// Added 2026-08-30 after the LiveCoMS tutorial suite exposed a class of defect the levels above
// cannot see: every system there uses one temperature bath, one value per line, and molecules that
// are never wrapped across the box. See README-livecoms.md in the reference directory.
ref_test!(aladip_multibath, "aladip_multibath");
ref_test!(aladip_multibath_collapsed, "aladip_multibath_collapsed");
ref_test!(aladip_wrapped, "aladip_wrapped");
ref_test!(aladip_harmonic_covalent, "aladip_harmonic_covalent");
// DOFSET `12 2 2  72 1 1`: the bath numbers do not ascend with the atoms (gromosXX `in_bath`).
ref_test!(aladip_multibath_swapped, "aladip_multibath_swapped");
// ROTTRANS on the solute of the two-bath system: 3 COM and 3 internal-rotational degrees of
// freedom leave the solute's baths (gromosXX `topology.cc:1370-1402`).
ref_test!(aladip_multibath_rottrans, "aladip_multibath_rottrans");
// Two Nosé–Hoover baths (plain and chain) with different TEMP0 and TAU: one ζ or chain per bath
// (gromosXX `nosehoover_thermostat.cc:120-230`).
ref_test!(aladip_multibath_nosehoover, "aladip_multibath_nosehoover");
ref_test!(aladip_multibath_nhc, "aladip_multibath_nhc");
// NDFMIN 3 over two baths whose degrees of freedom are integers (12, 13): gromosXX spreads it
// through an `int` total (`multibath.cc:145-152`), exact here. The half-integer case (14.5, 10.5)
// is `_half_fixed`, generated with D40 patched (the total summed exactly).
ref_test!(
    aladip_vacuum_multibath_ndfmin3,
    "aladip_vacuum_multibath_ndfmin3"
);
ref_test!(
    aladip_vacuum_multibath_ndfmin3_half_fixed,
    "aladip_vacuum_multibath_ndfmin3_half_fixed"
);
// Three energy groups: harmonic bonds/angles and dihedrals booked by first atom across groups.
ref_test!(
    aladip_harmonic_covalent_groups,
    "aladip_harmonic_covalent_groups"
);
// GaMD, fixed parameters (SEARCH 0), one interaction group over all atoms. FORM 3: the dihedral
// boost, ΔV = kD(V_D − ED)²/2 ≈ 1.23 kJ/mol in the special energy and ENER[51], the proper and
// improper dihedral forces scaled by 1 + kD(V_D − ED).
ref_test!(aladip_gamd, "aladip_gamd");
// FORM 1 with ET above V_T: both boosts act on the forces, and ΔV is the dihedral one only —
// gromosXX assigns `gamd_DV` twice (`gamd.cc:364, 381`), reproduced.
ref_test!(aladip_gamd_dual, "aladip_gamd_dual");
// SEARCH 1, NTIGAMDS 1: a cMD search from scratch, no boost; the statistics end in the final
// configuration's GAMDSTAT block.
ref_test!(aladip_gamd_cmd_search, "aladip_gamd_cmd_search");
