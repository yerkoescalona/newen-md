//! End-to-end MD benchmarks using reference systems.
//!
//! Measures wall-clock time for complete MD steps (forces + integration + constraints).
//! Uses the same systems as the reference tests to enable meaningful performance tracking.
//!
//! Run: cargo bench -p newen-md
//! Compare: cargo bench -p newen-md -- --baseline main

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::Duration;

fn ref_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

fn md_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_md"))
}

/// Run md binary on a reference system for N steps, return elapsed time.
fn run_md_steps(system: &str, nstlim: u32) -> Duration {
    let sys_dir = ref_root().join(system);
    let toml_content = fs::read_to_string(sys_dir.join("input.toml")).unwrap();

    let topo = sys_dir.join(toml_str(&toml_content, "topology"));
    let conf = sys_dir.join(toml_str(&toml_content, "configuration"));
    let params = sys_dir.join(toml_str(&toml_content, "parameters"));

    let out = std::env::temp_dir().join(format!("newen_bench_{}_{}", system, std::process::id()));
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(&out).unwrap();

    // Create modified input with desired NSTLIM and no trajectory output
    let input_text = fs::read_to_string(&params).unwrap();
    let modified_input = patch_nstlim(&input_text, nstlim);
    let bench_input = out.join("bench.in");
    fs::write(&bench_input, &modified_input).unwrap();

    let start = std::time::Instant::now();

    let result = Command::new(md_bin())
        .arg("@topo")
        .arg(&topo)
        .arg("@conf")
        .arg(&conf)
        .arg("@input")
        .arg(&bench_input)
        .arg("@fin")
        .arg(out.join("final.conf"))
        .output()
        .expect("failed to run md");

    let elapsed = start.elapsed();

    assert!(
        result.status.success(),
        "md failed on {system}: {}",
        String::from_utf8_lossy(&result.stdout)
    );

    let _ = fs::remove_dir_all(&out);
    elapsed
}

/// Extract a string value from TOML content
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
    panic!("{key} not found");
}

/// Patch NSTLIM in the input file and disable trajectory output for pure speed measurement.
fn patch_nstlim(input: &str, nstlim: u32) -> String {
    let mut result = String::new();
    let mut in_step = false;
    let mut in_writetraj = false;
    let mut step_comment_seen = false;
    let mut traj_comment_seen = false;

    for line in input.lines() {
        let trimmed = line.trim();

        if trimmed == "STEP" {
            in_step = true;
            step_comment_seen = false;
            result.push_str(line);
            result.push('\n');
            continue;
        }
        if trimmed == "WRITETRAJ" {
            in_writetraj = true;
            traj_comment_seen = false;
            result.push_str(line);
            result.push('\n');
            continue;
        }
        if trimmed == "END" {
            in_step = false;
            in_writetraj = false;
            result.push_str(line);
            result.push('\n');
            continue;
        }

        if in_step && !trimmed.starts_with('#') && !trimmed.is_empty() {
            if !step_comment_seen {
                step_comment_seen = true;
                result.push_str(line);
                result.push('\n');
                continue;
            }
            // This is the data line: NSTLIM T DT
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 3 {
                result.push_str(&format!(
                    "       {}       {}     {}\n",
                    nstlim, parts[1], parts[2]
                ));
            } else {
                result.push_str(line);
                result.push('\n');
            }
            continue;
        }

        if in_writetraj && !trimmed.starts_with('#') && !trimmed.is_empty() {
            if !traj_comment_seen {
                traj_comment_seen = true;
                result.push_str(line);
                result.push('\n');
                continue;
            }
            // Disable all trajectory output: NTWX NTWSE NTWV NTWF NTWE NTWG NTWB = 0
            result.push_str(
                "         0         0         0         0         0         0         0\n",
            );
            continue;
        }

        result.push_str(line);
        result.push('\n');
    }
    result
}

// ─── Benchmarks ─────────────────────────────────────────────────────────────

fn bench_water_216_nve(c: &mut Criterion) {
    let mut group = c.benchmark_group("md_end_to_end");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    // water_216_box: 648 atoms, NVE, SHAKE, twin-range pairlist
    group.bench_function("water_216_nve/100steps", |b| {
        b.iter(|| run_md_steps(black_box("water_216_box"), 100))
    });

    group.bench_function("water_216_nve/1000steps", |b| {
        b.iter(|| run_md_steps(black_box("water_216_box"), 1000))
    });

    group.finish();
}

fn bench_water_216_nvt(c: &mut Criterion) {
    let mut group = c.benchmark_group("md_end_to_end");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    // water_216_nvt: 648 atoms, NVT (Berendsen thermostat)
    group.bench_function("water_216_nvt/100steps", |b| {
        b.iter(|| run_md_steps(black_box("water_216_nvt"), 100))
    });

    group.finish();
}

fn bench_water_216_npt(c: &mut Criterion) {
    let mut group = c.benchmark_group("md_end_to_end");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    // water_216_npt: 648 atoms, NPT (thermostat + barostat + virial)
    group.bench_function("water_216_npt/100steps", |b| {
        b.iter(|| run_md_steps(black_box("water_216_npt"), 100))
    });

    group.finish();
}

fn bench_aladip_solvated(c: &mut Criterion) {
    let mut group = c.benchmark_group("md_end_to_end");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    // aladip_solvated: 72 atoms, SHAKE + solute-solvent + all bonded terms
    group.bench_function("aladip_solvated/100steps", |b| {
        b.iter(|| run_md_steps(black_box("aladip_solvated"), 100))
    });

    group.finish();
}

fn bench_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("md_scaling");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(20));

    // Compare systems of different sizes to measure scaling behavior
    for (system, label) in [
        ("nacl_water_box", "62atoms"),  // 62 atoms
        ("aladip_solvated", "72atoms"), // 72 atoms
        ("water_216_box", "648atoms"),  // 648 atoms
    ] {
        group.bench_with_input(BenchmarkId::new("100steps", label), &system, |b, &sys| {
            b.iter(|| run_md_steps(black_box(sys), 100))
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_water_216_nve,
    bench_water_216_nvt,
    bench_water_216_npt,
    bench_aladip_solvated,
    bench_scaling,
);
criterion_main!(benches);
