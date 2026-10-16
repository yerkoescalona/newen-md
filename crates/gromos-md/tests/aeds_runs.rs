//! Accelerated EDS through the `md` binary (1.8): a search continued over two jobs against the
//! same search in one job, and the AEDS inputs refused instead of reproduced (gromosXX crashes or
//! computes nonsense on them).

use std::path::{Path, PathBuf};
use std::process::Command;

fn refs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("gromos_aeds_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Run `md` on the aladip toy with `imd` (text) from `conf`; returns (success, stdout + stderr).
fn md(dir: &Path, imd: &str, conf: &Path, mpt: &str) -> (bool, String) {
    let input = dir.join("run.imd");
    std::fs::write(&input, imd).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_md"))
        .arg("@topo")
        .arg(refs().join("shared/aladip.topo"))
        .arg("@conf")
        .arg(conf)
        .arg("@input")
        .arg(&input)
        .arg("@pttopo")
        .arg(refs().join("shared").join(mpt))
        .arg("@fin")
        .arg(dir.join("final.cnf"))
        .arg("@tre")
        .arg(dir.join("ener.tre"))
        .output()
        .unwrap();
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    (out.status.success(), text)
}

fn reference_imd(name: &str) -> String {
    std::fs::read_to_string(refs().join(name).join(format!("{name}.in"))).unwrap()
}

fn start_conf() -> PathBuf {
    refs().join("aladip_vacuum/aladip_vacuum.conf")
}

/// The numbers of a configuration's AEDSSEARCH block.
fn aedssearch(path: &Path) -> Vec<f64> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .skip_while(|l| l.trim() != "AEDSSEARCH")
        .skip(1)
        .take_while(|l| l.trim() != "END")
        .flat_map(|l| l.split_whitespace().map(|x| x.parse::<f64>().unwrap()))
        .collect()
}

/// Replace the data of the AEDS block.
fn with_aeds(imd: &str, body: &str) -> String {
    let start = imd.find("\nAEDS\n").unwrap();
    let end = start + imd[start..].find("\nEND").unwrap() + 4;
    format!("{}\nAEDS\n{body}\nEND{}", &imd[..start], &imd[end..])
}

/// 100 + 100 steps through the AEDSSEARCH block equal 200 steps in one job up to the round-off
/// of the 9-decimal restart (memory time ASTEPS = BSTEPS, so its per-run restart changes
/// nothing). gromosXX's own chained and single runs differ by the same amount (5·10⁻⁸ relative
/// on this toy); the chained run is compared with gromosXX's exactly by the
/// `aladip_vacuum_aeds_cont1/2` references.
#[test]
fn a_search_continued_over_two_jobs_follows_the_single_job() {
    let job1 = scratch("job1");
    let (ok, log) = md(
        &job1,
        &reference_imd("aladip_vacuum_aeds_cont1"),
        &start_conf(),
        "aladip_aeds3.mpt",
    );
    assert!(ok, "{log}");
    let job2 = scratch("job2");
    let (ok, log) = md(
        &job2,
        &reference_imd("aladip_vacuum_aeds_cont2"),
        &job1.join("final.cnf"),
        "aladip_aeds3.mpt",
    );
    assert!(ok, "{log}");
    let single = scratch("single");
    let imd = reference_imd("aladip_vacuum_aeds_cont1").replace(
        "        100       0.0     0.002",
        "        200       0.0     0.002",
    );
    let (ok, log) = md(&single, &imd, &start_conf(), "aladip_aeds3.mpt");
    assert!(ok, "{log}");
    let (a, b) = (
        aedssearch(&job2.join("final.cnf")),
        aedssearch(&single.join("final.cnf")),
    );
    assert_eq!(a.len(), 6 + 3 * 10);
    assert_eq!(a.len(), b.len());
    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        assert!(
            (x - y).abs() <= 1e-7 * y.abs().max(1.0),
            "AEDSSEARCH value {}: chained {x} vs single {y}",
            i + 1
        );
    }
    for d in [job1, job2, single] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// What gromosXX accepts and then crashes on, or computes nonsense with, is refused with a
/// message naming the setting.
#[test]
fn defective_aeds_inputs_are_refused() {
    let f2 = reference_imd("aladip_vacuum_aeds_f2");
    let eds_block = "EDS\n 1\n 1 3\n 0.5\n 0 0 0\nEND\n";
    let cases: Vec<(&str, String, &str)> = vec![
        (
            "emin_above_emax",
            with_aeds(&f2, "1\n2 3\n0 10\n0 0 0\n1 1 1 0.1 20 50"),
            "EMIN = 10 is larger than EMAX = 0",
        ),
        (
            "bsteps_below_asteps",
            with_aeds(&f2, "1\n2 3\n5 -5\n0 0 0\n1 1 1 0.1 50 20"),
            "BSTEPS = 20 < ASTEPS = 50",
        ),
        (
            "with_eds",
            format!("{f2}{eds_block}"),
            "AEDS cannot be used in combination with EDS",
        ),
        (
            "no_aedssearch_block",
            with_aeds(&f2, "1\n2 3\n5 -5\n0 0 0\n0 1 1 0.1 20 50"),
            "no AEDSSEARCH block in configuration",
        ),
        (
            "zero_variance",
            with_aeds(
                &reference_imd("aladip_vacuum_aeds_sigma0"),
                "1\n3 2\n0 0\n-20 0\n1 0 2 2 100 100",
            ),
            "use BMAXTYPE 1 or RESTREMIN 1",
        ),
    ];
    for (name, imd, what) in cases {
        let dir = scratch(name);
        let mpt = if name == "zero_variance" {
            "aladip_aeds2.mpt"
        } else {
            "aladip_aeds3.mpt"
        };
        let (ok, log) = md(&dir, &imd, &start_conf(), mpt);
        assert!(!ok, "{name}: md accepted it\n{log}");
        assert!(log.contains(what), "{name}: expected {what:?} in\n{log}");
        let _ = std::fs::remove_dir_all(dir);
    }
}
