//! Integration tests for the `frameout` binary.
//!
//! Each test builds a minimal .trc input in a temp dir, runs the binary,
//! and verifies the output.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn frameout_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_frameout"))
}

fn temp_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("frameout_test_{}_{}", name, std::process::id()));
    fs::create_dir_all(&d).unwrap();
    d
}

/// Write a minimal 3-column .trc with `n_frames` frames and `n_atoms` atoms.
fn write_test_trc(path: &PathBuf, n_frames: usize, n_atoms: usize) {
    let mut content = "TITLE\n  test trajectory\nEND\n".to_string();
    for f in 0..n_frames {
        let step = f;
        let time = f as f64 * 0.002;
        content.push_str(&format!(
            "TIMESTEP\n{:15}{:20.9}\nEND\nPOSITIONRED\n",
            step, time
        ));
        for a in 0..n_atoms {
            let x = (f as f64) * 0.1 + (a as f64) * 0.01;
            content.push_str(&format!("{:15.9}{:15.9}{:15.9}\n", x, 0.5, 1.0));
        }
        content.push_str("END\nGENBOX\n   3.000000000   3.000000000   3.000000000\nEND\n");
    }
    fs::write(path, content).unwrap();
}

// ─── Format tests ────────────────────────────────────────────────────────────

#[test]
fn all_frames_cnf_output() {
    let dir = temp_dir("all_cnf");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 3, 2);

    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@spec",
            "ALL",
            "@out",
            "cnf",
            "@name",
            dir.join("frame").to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "frameout exited with {status}");

    // Should produce 3 separate .cnf files
    let cnf_files: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "cnf"))
        .collect();
    assert_eq!(
        cnf_files.len(),
        3,
        "expected 3 cnf files, got {}",
        cnf_files.len()
    );

    // Each cnf must contain POSITION block
    for f in &cnf_files {
        let content = fs::read_to_string(f.path()).unwrap();
        assert!(
            content.contains("POSITIONRED") || content.contains("POSITION"),
            "{:?} missing POSITION block",
            f.path()
        );
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn single_trc_output() {
    let dir = temp_dir("single_trc");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 5, 3);
    let out = dir.join("sub");

    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@spec",
            "ALL",
            "@out",
            "trc",
            "@single",
            "@name",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let out_trc = dir.join("sub.trc");
    assert!(out_trc.exists(), "sub.trc not created");

    let content = fs::read_to_string(&out_trc).unwrap();
    let n_timesteps = content.matches("TIMESTEP").count();
    assert_eq!(
        n_timesteps, 5,
        "expected 5 TIMESTEP blocks, got {n_timesteps}"
    );
    fs::remove_dir_all(dir).ok();
}

#[test]
fn every_n_frames() {
    let dir = temp_dir("every_n");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 10, 2);

    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@spec",
            "EVERY",
            "3",
            "@out",
            "cnf",
            "@name",
            dir.join("f").to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    // EVERY 3 from 10 frames: frames 0,3,6,9 → 4 files
    let cnf_count = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "cnf"))
        .count();
    assert_eq!(
        cnf_count, 4,
        "EVERY 3 of 10 should give 4 frames, got {cnf_count}"
    );
    fs::remove_dir_all(dir).ok();
}

#[test]
fn spec_frame_list() {
    let dir = temp_dir("spec_list");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 8, 2);

    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@frames",
            "0,2,4",
            "@out",
            "cnf",
            "@name",
            dir.join("s").to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let cnf_count = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "cnf"))
        .count();
    assert_eq!(cnf_count, 3, "frames 0,2,4 → 3 files, got {cnf_count}");
    fs::remove_dir_all(dir).ok();
}

#[test]
fn time_range_filter() {
    let dir = temp_dir("time_range");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 10, 2); // t = 0.000, 0.002, 0.004 ... 0.018 ps

    // Select t in [0.004, 0.010] → frames 2,3,4,5 → 4 files
    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@time",
            "0.004",
            "0.010",
            "@out",
            "cnf",
            "@name",
            dir.join("t").to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let cnf_count = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "cnf"))
        .count();
    assert_eq!(
        cnf_count, 4,
        "t in [0.004,0.010] → 4 frames, got {cnf_count}"
    );
    fs::remove_dir_all(dir).ok();
}

#[test]
fn pdb_output_format() {
    let dir = temp_dir("pdb_out");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 1, 3);

    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@out",
            "pdb",
            "@name",
            dir.join("p").to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let pdb_files: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "pdb"))
        .collect();
    assert_eq!(pdb_files.len(), 1);

    let content = fs::read_to_string(pdb_files[0].path()).unwrap();
    let atom_count = content.lines().filter(|l| l.starts_with("ATOM")).count();
    assert_eq!(atom_count, 3, "expected 3 ATOM records, got {atom_count}");
    assert!(content.contains("END"));
    fs::remove_dir_all(dir).ok();
}

#[test]
fn single_cnf_overwrites_with_last_frame() {
    let dir = temp_dir("single_cnf");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 4, 2);
    let out_name = dir.join("out");

    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@single",
            "@out",
            "cnf",
            "@name",
            out_name.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let cnf = dir.join("out.cnf");
    assert!(cnf.exists(), "out.cnf not created");
    // Only one file, not four
    let cnf_count = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "cnf"))
        .count();
    assert_eq!(cnf_count, 1);
    fs::remove_dir_all(dir).ok();
}

#[test]
fn trc_positions_are_3column() {
    // Verify the output .trc uses standard 3-column POSITIONRED
    let dir = temp_dir("trc_3col");
    let trc = dir.join("in.trc");
    write_test_trc(&trc, 1, 2);

    let status = Command::new(frameout_bin())
        .args([
            "@traj",
            trc.to_str().unwrap(),
            "@out",
            "trc",
            "@single",
            "@name",
            dir.join("out").to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let content = fs::read_to_string(dir.join("out.trc")).unwrap();
    let pos_line = content
        .lines()
        .skip_while(|l| !l.starts_with("POSITIONRED"))
        .nth(1)
        .unwrap()
        .trim()
        .to_string();
    let cols: Vec<&str> = pos_line.split_whitespace().collect();
    assert_eq!(
        cols.len(),
        3,
        "POSITIONRED must be 3-column, got: {pos_line}"
    );
    fs::remove_dir_all(dir).ok();
}
