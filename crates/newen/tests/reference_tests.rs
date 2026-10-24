#![allow(clippy::needless_range_loop)] // reference comparisons index several arrays in lockstep
//! Integration tests validating newen against GROMOS reference data
//!
//! These tests parse the reference topology/coordinates/parameters,
//! compute forces using the newen force kernels, and compare
//! against the expected forces/energies from GROMOS (double precision).
//!
//! Reference data is in GROMOS_references/ and must NOT be modified.
//!
//! Test levels (simplest first):
//!   Level 0: pair_lj, nacl_pair, pair_lj_mixed (2 atoms, vacuum, no PBC)
//!   Level 1: water_single, aladip_vacuum, benzene_vacuum
//!   Level 2: water_3_box, nacl_water_box
//!   Level 3: water_216_box, water_216_nvt, water_216_npt
//!   Level 4: aladip_solvated

use newen::forces::nonbonded::{
    lj_crf_innerloop, CRFParameters, EnergyGroups, ForceStorage, LJParamMatrix,
};
use newen::io::coordinate::read_coordinate_file;
use newen::io::topology::{build_topology, read_topology_file};
use newen::math::Vec3;
use newen::topology::LJParameters;

use std::path::{Path, PathBuf};

// Initialize env_logger once for all tests (safe to call multiple times)
fn init_logging() {
    let _ = env_logger::try_init();
}

// ─── Helpers ────────────────────────────────────────────────────────────────

/// Path to gromosXX_references (in newen-md/tests/)
fn ref_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("newen-md/tests/gromosXX_references")
}

/// Parse the FREEFORCERED block from forces.trf at a given step
fn parse_forces_trf(path: &Path, step: usize) -> Vec<Vec3> {
    let content = std::fs::read_to_string(path).expect("Cannot read forces.trf");
    let mut forces = Vec::new();
    let mut current_step = None;
    let mut in_freeforce = false;

    for line in content.lines() {
        let trimmed = line.trim();

        if trimmed == "TIMESTEP" {
            // Next non-comment/non-empty line has step number
            current_step = None;
            continue;
        }

        if current_step.is_none()
            && trimmed != "END"
            && !trimmed.is_empty()
            && !trimmed.starts_with('#')
        {
            // Parse step number (first field on the line)
            if let Some(step_str) = trimmed.split_whitespace().next() {
                if let Ok(s) = step_str.parse::<usize>() {
                    current_step = Some(s);
                }
            }
            continue;
        }

        if trimmed == "FREEFORCERED" {
            if current_step == Some(step) {
                in_freeforce = true;
                forces.clear();
            }
            continue;
        }

        if trimmed == "END" {
            if in_freeforce {
                return forces;
            }
            continue;
        }

        if in_freeforce {
            let parts: Vec<f64> = trimmed
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if parts.len() == 3 {
                forces.push(Vec3::new(parts[0], parts[1], parts[2]));
            }
        }
    }

    forces
}

/// Parse ENERGY03 block from energies.tre at a given step.
///
/// Reads the `# totals` section (one float per line) and the `# nonbonded`
/// section (6 floats per line: e_lj, e_crf, ...) summing across molecule-type
/// pair rows to produce system totals.
fn parse_energies_tre(path: &Path, step: usize) -> EnergyValues {
    let content = std::fs::read_to_string(path).expect("Cannot read energies.tre");
    let mut current_step = None;
    let mut in_energy03 = false;
    let mut in_nonbonded = false;
    let mut totals: Vec<f64> = Vec::new();
    let mut e_lj = 0.0_f64;
    let mut e_crf = 0.0_f64;

    for line in content.lines() {
        let trimmed = line.trim();

        if trimmed == "TIMESTEP" {
            current_step = None;
            in_energy03 = false;
            in_nonbonded = false;
            continue;
        }

        if current_step.is_none()
            && trimmed != "END"
            && !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && !trimmed.starts_with("ENEVERSION")
            && !trimmed.starts_with("2023")
        {
            if let Some(step_str) = trimmed.split_whitespace().next() {
                if let Ok(s) = step_str.parse::<usize>() {
                    current_step = Some(s);
                }
            }
            continue;
        }

        if trimmed == "ENERGY03" && current_step == Some(step) {
            in_energy03 = true;
            totals.clear();
            e_lj = 0.0;
            e_crf = 0.0;
            continue;
        }

        if trimmed == "END" {
            if in_energy03 {
                break;
            }
            continue;
        }

        if !in_energy03 || trimmed.is_empty() {
            continue;
        }

        if trimmed.starts_with('#') {
            in_nonbonded = trimmed == "# nonbonded";
            continue;
        }

        if in_nonbonded {
            // Each row: e_lj  e_crf  ... (6 values, one row per molecule-type pair).
            let vals: Vec<f64> = trimmed
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if vals.len() >= 2 {
                e_lj += vals[0];
                e_crf += vals[1];
            }
        } else if totals.len() < 3 {
            // # totals section has one float per line.
            if let Ok(v) = trimmed.parse::<f64>() {
                totals.push(v);
            }
        }
    }

    EnergyValues {
        e_total: totals.first().copied().unwrap_or(0.0),
        e_kinetic: totals.get(1).copied().unwrap_or(0.0),
        e_potential: totals.get(2).copied().unwrap_or(0.0),
        e_lj,
        e_crf,
    }
}

#[derive(Debug)]
struct EnergyValues {
    #[allow(dead_code)]
    e_total: f64,
    #[allow(dead_code)]
    e_kinetic: f64,
    #[allow(dead_code)]
    e_potential: f64,
    e_lj: f64,
    #[allow(dead_code)]
    e_crf: f64,
}

/// Convert topology LJ parameters to the nonbonded format
fn convert_lj_params(topo_lj: &[Vec<LJParameters>]) -> LJParamMatrix {
    let nested: Vec<Vec<newen::forces::nonbonded::LJParameters>> = topo_lj
        .iter()
        .map(|row| {
            row.iter()
                .map(|p| newen::forces::nonbonded::LJParameters {
                    c6: p.c6,
                    c12: p.c12,
                    cs6: p.cs6,
                    cs12: p.cs12,
                })
                .collect()
        })
        .collect();
    LJParamMatrix::from_nested(&nested)
}

/// Assert two Vec3 are approximately equal (absolute tolerance)
fn assert_vec3_approx(actual: Vec3, expected: Vec3, tol: f64, msg: &str) {
    let dx = (actual.x - expected.x).abs();
    let dy = (actual.y - expected.y).abs();
    let dz = (actual.z - expected.z).abs();
    assert!(
        dx < tol && dy < tol && dz < tol,
        "{}: expected ({:.9e}, {:.9e}, {:.9e}), got ({:.9e}, {:.9e}, {:.9e}), max_diff={:.2e}",
        msg,
        expected.x,
        expected.y,
        expected.z,
        actual.x,
        actual.y,
        actual.z,
        dx.max(dy).max(dz)
    );
}

/// Assert f64 values are approximately equal (relative tolerance)
fn assert_f64_rel(actual: f64, expected: f64, tol: f64, msg: &str) {
    let diff = (actual - expected).abs();
    let scale = expected.abs().max(1e-15);
    let rel = diff / scale;
    assert!(
        rel < tol,
        "{}: expected {:.10e}, got {:.10e}, rel_diff={:.2e}",
        msg,
        expected,
        actual,
        rel
    );
}

// ─── Level 0: pair_lj ───────────────────────────────────────────────────────

#[test]
fn test_pair_lj_forces_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj");

    // Read topology
    let parsed_topo = read_topology_file(base.join("pair_lj.topo")).unwrap();
    let topo = build_topology(parsed_topo);

    // Read coordinates
    let conf = read_coordinate_file(base.join("pair_lj.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    // Set up force calculation
    let n_atoms = positions.len();
    assert_eq!(n_atoms, 2);

    // Build pairlist: all pairs (no exclusions for this system)
    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];

    // IAC and charges
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();

    // LJ parameters
    let lj_params = convert_lj_params(&topo.lj_parameters);

    // CRF parameters (no electrostatics needed, charges are 0)
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    // Compute forces
    let bc = newen::core::math::Vacuum;
    let mut storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        newen_core::units::four_pi_eps_i,
        &mut storage,
        EnergyGroups::single(),
    );

    // Load expected forces
    let expected_forces = parse_forces_trf(&base.join("expected/forces.trf"), 0);
    assert_eq!(expected_forces.len(), n_atoms);

    // Compare
    let tol = 1e-6; // kJ/(mol*nm)
    for i in 0..n_atoms {
        assert_vec3_approx(
            storage.forces[i],
            expected_forces[i],
            tol,
            &format!("pair_lj force atom {}", i),
        );
    }
}

#[test]
fn test_pair_lj_energy_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj");

    let parsed_topo = read_topology_file(base.join("pair_lj.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("pair_lj.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = newen::core::math::Vacuum;
    let mut storage = ForceStorage::new(positions.len());
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        newen_core::units::four_pi_eps_i,
        &mut storage,
        EnergyGroups::single(),
    );

    let energies = parse_energies_tre(&base.join("expected/energies.tre"), 0);
    let tol = 1e-8;
    assert_f64_rel(storage.e_lj, energies.e_lj, tol, "pair_lj E_Vdw");
}

// ─── Level 0: nacl_pair ─────────────────────────────────────────────────────

#[test]
fn test_nacl_pair_forces_step0() {
    init_logging();
    let base = ref_dir().join("nacl_pair");

    let parsed_topo = read_topology_file(base.join("nacl_pair.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("nacl_pair.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let n_atoms = positions.len();
    assert_eq!(n_atoms, 2);

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);

    // nacl_pair.in: EPSRF=1.0 (vacuum), RCRF=RCUTL=1.4 nm.
    // gromosXX with EPSRF=1 uses plain Coulomb: crf=0, so all correction
    // terms vanish (crf_cut=0, not 1/rc — confirmed by oracle e_crf=-463.118
    // which equals four_pi_eps_i*q1*q2/r with no shift).
    let crf = CRFParameters {
        crf_cut: 0.0,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = newen::core::math::Vacuum;
    let mut storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        newen_core::units::four_pi_eps_i,
        &mut storage,
        EnergyGroups::single(),
    );

    // Load expected forces
    let expected_forces = parse_forces_trf(&base.join("expected/forces.trf"), 0);
    assert_eq!(expected_forces.len(), n_atoms);

    // Oracle uses gromosXX legacy constants (138.9354) from the topo PHYSICALCONSTANTS block.
    let tol = 1e-9;
    for i in 0..n_atoms {
        assert_vec3_approx(
            storage.forces[i],
            expected_forces[i],
            tol,
            &format!("nacl_pair force atom {}", i),
        );
    }
}

#[test]
fn test_nacl_pair_energy_step0() {
    init_logging();
    let base = ref_dir().join("nacl_pair");

    let parsed_topo = read_topology_file(base.join("nacl_pair.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("nacl_pair.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    // EPSRF=1.0 → plain Coulomb, all CRF correction terms zero.
    let crf = CRFParameters {
        crf_cut: 0.0,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = newen::core::math::Vacuum;
    let mut storage = ForceStorage::new(positions.len());
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        newen_core::units::four_pi_eps_i,
        &mut storage,
        EnergyGroups::single(),
    );

    let energies = parse_energies_tre(&base.join("expected/energies.tre"), 0);
    assert_f64_rel(storage.e_lj, energies.e_lj, 1e-8, "nacl_pair E_Vdw");

    // Oracle uses gromosXX legacy constants from the topo — constants agree so tolerance is tight.
    assert_f64_rel(storage.e_crf, energies.e_crf, 1e-7, "nacl_pair E_crf");
}

// ─── Level 0: pair_lj_mixed ─────────────────────────────────────────────────

#[test]
fn test_pair_lj_mixed_forces_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj_mixed");

    let parsed_topo = read_topology_file(base.join("pair_lj_mixed.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("pair_lj_mixed.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let n_atoms = positions.len();
    assert_eq!(n_atoms, 2);

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = newen::core::math::Vacuum;
    let mut storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        newen_core::units::four_pi_eps_i,
        &mut storage,
        EnergyGroups::single(),
    );

    // Load expected forces
    let expected_forces = parse_forces_trf(&base.join("expected/forces.trf"), 0);
    assert_eq!(expected_forces.len(), n_atoms);

    let tol = 1e-6;
    for i in 0..n_atoms {
        assert_vec3_approx(
            storage.forces[i],
            expected_forces[i],
            tol,
            &format!("pair_lj_mixed force atom {}", i),
        );
    }
}

#[test]
fn test_pair_lj_mixed_energy_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj_mixed");

    let parsed_topo = read_topology_file(base.join("pair_lj_mixed.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("pair_lj_mixed.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = newen::core::math::Vacuum;
    let mut storage = ForceStorage::new(positions.len());
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        newen_core::units::four_pi_eps_i,
        &mut storage,
        EnergyGroups::single(),
    );

    let energies = parse_energies_tre(&base.join("expected/energies.tre"), 0);
    let tol = 1e-8;
    assert_f64_rel(storage.e_lj, energies.e_lj, tol, "pair_lj_mixed E_Vdw");
}
