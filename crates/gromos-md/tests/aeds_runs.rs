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

/// Run `md` on the aladip toy with `imd` (text) from `conf`, `extra` appended to the command
/// line; returns (success, stdout + stderr).
fn md(dir: &Path, imd: &str, conf: &Path, mpt: &str, extra: &[&str]) -> (bool, String) {
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
        .args(extra)
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
        &[],
    );
    assert!(ok, "{log}");
    let job2 = scratch("job2");
    let (ok, log) = md(
        &job2,
        &reference_imd("aladip_vacuum_aeds_cont2"),
        &job1.join("final.cnf"),
        "aladip_aeds3.mpt",
        &[],
    );
    assert!(ok, "{log}");
    let single = scratch("single");
    let imd = reference_imd("aladip_vacuum_aeds_cont1").replace(
        "        100       0.0     0.002",
        "        200       0.0     0.002",
    );
    let (ok, log) = md(&single, &imd, &start_conf(), "aladip_aeds3.mpt", &[]);
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
        // D15 is reached through gromosXX's never-visited target state 1 (D32)
        let extra: &[&str] = if name == "zero_variance" {
            &["@recipe_set", "aeds.target_state=gromosxx-state1-default"]
        } else {
            &[]
        };
        let (ok, log) = md(&dir, &imd, &start_conf(), mpt, extra);
        assert!(!ok, "{name}: md accepted it\n{log}");
        assert!(log.contains(what), "{name}: expected {what:?} in\n{log}");
        let _ = std::fs::remove_dir_all(dir);
    }
}

/// One energy frame of an A-EDS `.tre`: ENER (1-based slots) and the `# eds` rows.
struct Frame {
    ener: Vec<f64>,
    eds: Vec<Vec<f64>>,
}

fn energy_frames(path: &Path) -> Vec<Frame> {
    let text = std::fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = text.lines().map(str::trim).collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        if lines[i] == "ENERGY03" {
            let mut f = Frame {
                ener: Vec::new(),
                eds: Vec::new(),
            };
            i += 1;
            while lines[i] != "END" {
                if lines[i] == "# totals" {
                    i += 1;
                    while !lines[i].starts_with('#') {
                        f.ener.push(lines[i].parse().unwrap());
                        i += 1;
                    }
                    continue;
                }
                if lines[i] == "# eds" {
                    let n: usize = lines[i + 2].parse().unwrap();
                    for k in 0..n {
                        let row = lines[i + 4 + k].split_whitespace();
                        f.eds.push(row.map(|x| x.parse().unwrap()).collect());
                    }
                    i += 4 + n;
                    continue;
                }
                i += 1;
            }
            out.push(f);
        }
        i += 1;
    }
    out
}

/// V_mix = −1/β ln Σ_i e^{−β(V_i − E_i^R)} and V_R* (the quadratic boost between E_min and
/// E_max), written out from the A-EDS definition, from the printed e_i of frame `k` and the
/// printed parameters of frame `p`; returns the largest absolute deviation from frame k's
/// printed V_mix (slot 34) and V_R* (slot 35).
fn rebuild_error(frames: &[Frame], k: usize, p: usize) -> f64 {
    let beta = 1.0 / (0.00831441 * 300.0);
    let x: Vec<f64> = (0..frames[k].eds.len())
        .map(|i| -beta * (frames[k].eds[i][0] - frames[p].eds[i][3]))
        .collect();
    let m = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let v_mix = -(m + x.iter().map(|xi| (xi - m).exp()).sum::<f64>().ln()) / beta;
    let (emax, emin) = (frames[p].ener[35], frames[p].ener[36]);
    let v_r = if v_mix <= emin {
        v_mix
    } else if v_mix >= emax {
        v_mix - 0.5 * (emax - emin)
    } else {
        v_mix - (v_mix - emin).powi(2) / (2.0 * (emax - emin))
    };
    let printed = &frames[k].ener;
    (v_mix - printed[33]).abs().max((v_r - printed[34]).abs())
}

/// D16: by default a frame's V_mix and V_R* follow from that frame's own e_i, offsets, E_max
/// and E_min (to the printed digits), what the LiveCoMS t_05 `search.py` assumes; gromosXX's
/// columns (`aeds.energy_columns = gromosxx-post-update`) need the previous frame's
/// parameters instead.
#[test]
fn a_frame_carries_the_parameters_its_reference_energy_used() {
    let imd = reference_imd("aladip_vacuum_aeds_f4");
    for (columns, own_frame) in [("intended", true), ("gromosxx-post-update", false)] {
        let dir = scratch(&format!("columns_{columns}"));
        let set = format!("aeds.energy_columns={columns}");
        let (ok, log) = md(
            &dir,
            &imd,
            &start_conf(),
            "aladip_aeds3.mpt",
            &["@recipe_set", &set],
        );
        assert!(ok, "{log}");
        let frames = energy_frames(&dir.join("ener.tre"));
        assert_eq!(frames.len(), 300);
        let own = (0..frames.len())
            .map(|k| rebuild_error(&frames, k, k))
            .fold(0.0, f64::max);
        let previous = (1..frames.len())
            .map(|k| rebuild_error(&frames, k, k - 1))
            .fold(0.0, f64::max);
        // ten significant digits printed, |V| ≲ 30 kJ/mol
        let printed = 1e-7;
        if own_frame {
            assert!(own < printed, "{columns}: own frame off by {own}");
            assert!(
                previous > 1e-3,
                "{columns}: previous frame fits too ({previous})"
            );
        } else {
            assert!(
                previous < printed,
                "{columns}: previous frame off by {previous}"
            );
            assert!(own > 1e-3, "{columns}: own frame fits too ({own})");
        }
        let _ = std::fs::remove_dir_all(dir);
    }
}

/// D32: gromosXX's never-visited target state 1 (σ = 0) makes the FORM 3, BMAXTYPE 2,
/// RESTREMIN 0 search divide by b_max = 0 (refused above); with the sampled state as the
/// target the same input runs, and its G is the sampled state's mean, not 0.
#[test]
fn the_e_min_target_is_a_sampled_state() {
    let imd = with_aeds(
        &reference_imd("aladip_vacuum_aeds_sigma0"),
        "1\n3 2\n0 0\n-20 0\n1 0 2 2 100 100",
    );
    let dir = scratch("target_intended");
    let (ok, log) = md(&dir, &imd, &start_conf(), "aladip_aeds2.mpt", &[]);
    assert!(ok, "{log}");
    let frames = energy_frames(&dir.join("ener.tre"));
    for (k, f) in frames.iter().enumerate().skip(1) {
        // only state 2 is ever visited (E_1^R = −20): G is the mean of its V_2 − E_2^R
        let (emin, g, sigma) = (f.ener[36], f.ener[37], f.ener[38]);
        assert!(
            emin.is_finite() && g > 1.0 && sigma > 0.0,
            "frame {k}: {emin} {g} {sigma}"
        );
        let mean: f64 = frames[..=k].iter().map(|f| f.eds[1][0]).sum::<f64>() / (k + 1) as f64;
        assert!(
            (g - mean).abs() < 1e-6,
            "frame {k}: G {g} vs mean V_2 {mean}"
        );
    }
    let _ = std::fs::remove_dir_all(dir);
}
