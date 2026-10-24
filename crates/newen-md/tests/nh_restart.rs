//! Nosé–Hoover runs continued through the final configuration's NHCVARIABLES block (D41):
//! 100 + 100 steps follow the 200-step run frame by frame, for plain Nosé–Hoover and for a 3-link
//! chain, as closely as weak coupling (no thermostat state) does: what is left is the 9-decimal
//! round-off of the restart, grown over 100 steps.

use std::path::{Path, PathBuf};
use std::process::Command;

fn refs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("newen_nh_restart_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Run `md` on the 216 waters with `imd` (text) from `conf`; panics with the log on failure.
fn md(dir: &Path, imd: &str, conf: &Path, extra: &[&str]) {
    let input = dir.join("run.imd");
    std::fs::write(&input, imd).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_md"))
        .arg("@topo")
        .arg(refs().join("water_216_box/water_216_box.topo"))
        .arg("@conf")
        .arg(conf)
        .arg("@input")
        .arg(&input)
        .arg("@fin")
        .arg(dir.join("final.cnf"))
        .arg("@tre")
        .arg(dir.join("ener.tre"))
        .args(extra)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

/// Every number of every `block` of a file, one vector per occurrence.
fn blocks(path: &Path, block: &str) -> Vec<Vec<f64>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = Vec::new();
    let mut cur: Option<Vec<f64>> = None;
    for l in text.lines().map(str::trim) {
        if l == block {
            cur = Some(Vec::new());
        } else if l == "END" {
            out.extend(cur.take());
        } else if let Some(v) = cur.as_mut().filter(|_| !l.starts_with('#')) {
            v.extend(l.split_whitespace().filter_map(|x| x.parse::<f64>().ok()));
        }
    }
    out
}

/// The largest |a − b| / max(|b|, 1) over two equally long lists.
fn max_rel(a: &[f64], b: &[f64], what: &str) -> f64 {
    assert_eq!(a.len(), b.len(), "{what}: lengths differ");
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs() / y.abs().max(1.0))
        .fold(0.0, f64::max)
}

/// The first half (steps 0-99), the continuation (100-199) and the single run (0-199) of the
/// continuation reference `name`, whose `.in` is the second half; `weak` swaps the thermostat
/// for weak coupling.
fn three_inputs(name: &str, chain: bool, weak: bool) -> (String, String, String) {
    let mut second = std::fs::read_to_string(refs().join(name).join(format!("{name}.in"))).unwrap();
    if weak {
        let algorithm = if chain {
            "    nose-hoover-chains\t3\n"
        } else {
            "    nose-hoover\n"
        };
        assert!(second.contains(algorithm));
        second = second.replace(algorithm, "    weak-coupling\n");
    }
    let ntinht = if chain {
        "        0       0       1       0"
    } else {
        "        0       0       0       0"
    };
    let first = second
        .replace(
            "        0       0       0       0\n",
            &format!("{ntinht}\n"),
        )
        .replace(
            "\n        0       0       0\n",
            "\n        1       0       0\n",
        )
        .replace(
            "       100       0.2     0.002",
            "       100       0.0     0.002",
        );
    let full = first.replace(
        "       100       0.0     0.002",
        "       200       0.0     0.002",
    );
    assert_ne!(first, second);
    assert_ne!(full, first);
    (first, second, full)
}

/// The largest relative difference, per output, between the continued run and the single run
/// (frames 100-199 of the energy blocks; the final configuration).
fn continued_vs_single(name: &str, chain: bool, weak: bool) -> Vec<(&'static str, f64)> {
    let (first, second, full) = three_inputs(name, chain, weak);
    let start = refs().join("water_216_box/water_216_box.conf");
    let tag = format!("{name}_{weak}");
    let (j1, j2, single) = (
        scratch(&format!("{tag}_job1")),
        scratch(&format!("{tag}_job2")),
        scratch(&format!("{tag}_single")),
    );
    md(&j1, &first, &start, &[]);
    md(&j2, &second, &j1.join("final.cnf"), &[]);
    md(&single, &full, &start, &[]);

    let zeta = |dir: &Path| blocks(&dir.join("final.cnf"), "NHCVARIABLES");
    let links = if weak {
        0
    } else if chain {
        3
    } else {
        1
    };
    assert_eq!(zeta(&j1).concat().len(), links, "{tag}: NHCVARIABLES");
    let mut worst = Vec::new();
    for block in ["ENERGY03", "VOLUMEPRESSURE03"] {
        let cont = blocks(&j2.join("ener.tre"), block);
        let one = blocks(&single.join("ener.tre"), block);
        assert_eq!(
            cont.len(),
            100,
            "{name}: {block} frames of the continuation"
        );
        assert_eq!(one.len(), 200, "{name}: {block} frames of the single run");
        let d = cont
            .iter()
            .zip(&one[100..])
            .enumerate()
            .map(|(k, (a, b))| max_rel(a, b, &format!("{name} {block} frame {k}")))
            .fold(0.0, f64::max);
        worst.push((block, d));
    }
    for block in ["POSITION", "VELOCITY", "NHCVARIABLES"] {
        let a: Vec<f64> = blocks(&j2.join("final.cnf"), block).concat();
        let b: Vec<f64> = blocks(&single.join("final.cnf"), block).concat();
        worst.push((block, max_rel(&a, &b, &format!("{name} final {block}"))));
    }
    eprintln!("{tag}: continued vs single, largest relative difference {worst:?}");
    for d in [j1, j2, single] {
        let _ = std::fs::remove_dir_all(d);
    }
    worst
}

fn continued_run_follows_the_single_run(name: &str, chain: bool) {
    let nh = continued_vs_single(name, chain, false);
    let control = continued_vs_single(name, chain, true);
    for ((what, d), (_, c)) in nh.iter().zip(&control) {
        // NHCVARIABLES has no control: ζ itself is printed with 9 decimals
        let bound = if *what == "NHCVARIABLES" {
            1e-6
        } else {
            10.0 * c.max(1e-9)
        };
        assert!(
            *d <= bound,
            "{name}: {what} of the continued run differs from the single run by {d:e} \
             (weak coupling: {c:e})"
        );
    }
}

#[test]
fn plain_nose_hoover_continued_follows_the_single_run() {
    continued_run_follows_the_single_run("water_216_nvt_nosehoover_cont_fixed", false);
}

#[test]
fn nose_hoover_chain_continued_follows_the_single_run() {
    continued_run_follows_the_single_run("water_216_nvt_nhc_cont_fixed", true);
}
