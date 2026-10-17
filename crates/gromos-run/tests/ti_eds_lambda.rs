//! TI-EDS (1.8): the λ derivative the force field books per end state, ∂V_i/∂λ, and the one the
//! EDS step books for the reference state, ∂V_R/∂λ = Σ p_i ∂V_i/∂λ (`.trg` slot 35), against
//! central differences of V_i and V_R between runs at λ ± δ, at step 0 (force field and EDS
//! step, no integration), for the three forms of the reference Hamiltonian.
//!
//! The system is `aladip_eds_solv_ti` (aladip in 20 SPC, EDS atoms N4 and H11, perturbed atoms
//! CB1, C8, N10, CB12: EDS–perturbed pairlist, 1-4 and excluded pairs) with three end states, so
//! that FORMS 2 and 3 differ from FORM 1; the offsets give every state a weight.

use std::path::{Path, PathBuf};

use gromos_core::algorithm::{AlgorithmSequence, SimulationState};
use gromos_integrators::algorithms::Forcefield;
use gromos_integrators::eds::{reference_hamiltonian, Eds};
use gromos_io::coordinate::read_coordinates;
use gromos_io::imd::{read_imd_file, EdsBlock, EdsStatePair, ImdParameters};
use gromos_io::topology::{build_topology, read_topology_file};
use gromos_run::{build_sequence_from_imd, prepare_system, Coordinates, RunInputs, RunOptions};

fn refs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
}

/// `aladip_eds_solv_ti.ptp` with a third end state on N4 and H11.
const PTP: &str = "TITLE
TI-EDS finite differences: three end states
END
PERTATOMPARAM
    4
    1     1 CB      14  15.0350    0.000   12  15.0350    0.150    0.3       0.6
    8     2 C       11  12.0110    0.380   11  12.0110    0.250    0.5       0.4
   10     3 N        5  14.0067   -0.280    1  14.0067   -0.400    0.2       0.7
   12     3 CB      14  15.0350    0.000   18  15.0350    0.100    0.4       0.5
END
MPERTATOM
    2    3
  stateA  stateB  stateC
   4  N        5   -0.280    19    0.200     1   -0.500
  11  H       18    0.300    14    0.000    18    0.450
END
";

/// What the force field and the EDS step book at step 0 of a run.
struct Booked {
    vr: f64,
    vi: Vec<f64>,
    d_vr: f64,
    d_vi: Vec<f64>,
    weights: Vec<f64>,
}

/// The step-0 bookings of the run at `lambda`; `eir` replaces the offsets when given.
fn at_lambda(imd: &ImdParameters, ptp: &Path, lambda: f64, eir: Option<&[f64]>) -> Booked {
    let mut imd = imd.clone();
    imd.rlam = lambda;
    let parsed = read_topology_file(refs().join("shared/aladip.topo")).unwrap();
    let constants = parsed.physical_constants;
    let coords: Coordinates = read_coordinates(refs().join("shared/aladip.conf"))
        .unwrap()
        .into();
    let inputs = RunInputs {
        pttopo: Some(ptp.to_path_buf()),
        ..Default::default()
    };
    let prepared =
        prepare_system(&imd, build_topology(parsed), constants, coords, &inputs).unwrap();
    let built = build_sequence_from_imd(&imd, &prepared, &inputs, &RunOptions::default()).unwrap();
    let mut seq: AlgorithmSequence = built.sequence;
    let topo = prepared.topology;
    let mut conf = prepared.configuration;
    let sim = SimulationState::new(0.002, 1);
    seq.init(&topo, &mut conf, &sim).unwrap();
    use gromos_core::algorithm::Algorithm;
    seq.find_mut::<Forcefield>()
        .expect("force field")
        .apply(&topo, &mut conf, &sim)
        .unwrap();
    let eds = seq.find_mut::<Eds>().expect("EDS step");
    if let Some(eir) = eir {
        eds.eir = eir.to_vec();
    }
    eds.apply(&topo, &mut conf, &sim).unwrap();
    let (e, d) = (
        &conf.current().energies,
        &conf.current().perturbed_energy_derivatives,
    );
    let mut weights = vec![0.0; e.eds_vi.len()];
    reference_hamiltonian(&eds.form, eds.beta, &eds.eir, &e.eds_vi, &mut weights);
    Booked {
        vr: e.eds_vr,
        vi: e.eds_vi.clone(),
        d_vr: d.eds_vr,
        d_vi: d.eds_vi.clone(),
        weights,
    }
}

fn check_form(form: i32, s: Vec<f64>, pairs: Vec<EdsStatePair>, offsets: [f64; 3]) {
    let dir = refs().join("aladip_eds_solv_ti");
    let mut imd = read_imd_file(dir.join("aladip_eds_solv_ti.in")).unwrap();
    imd.eds = Some(EdsBlock {
        enabled: true,
        form,
        num_states: 3,
        s,
        pairs,
        eir: vec![0.0; 3],
    });
    let ptp = std::env::temp_dir().join(format!(
        "gromos_ti_eds_lambda_{form}_{}.ptp",
        std::process::id()
    ));
    std::fs::write(&ptp, PTP).unwrap();

    let lambda = 0.3;
    let vi = at_lambda(&imd, &ptp, lambda, None).vi;
    let eir: Vec<f64> = (0..3).map(|k| vi[k] + offsets[k]).collect();
    let Booked {
        d_vr,
        d_vi,
        weights: w,
        ..
    } = at_lambda(&imd, &ptp, lambda, Some(&eir));
    assert!(w.iter().all(|&p| p > 0.05), "FORM {form}: weights {w:?}");
    assert!(
        d_vi.iter().all(|d| d.abs() > 1e-2),
        "FORM {form}: ∂V_i/∂λ {d_vi:?}"
    );

    let central = |delta: f64| {
        let up = at_lambda(&imd, &ptp, lambda + delta, Some(&eir));
        let down = at_lambda(&imd, &ptp, lambda - delta, Some(&eir));
        let d_vr = (up.vr - down.vr) / (2.0 * delta);
        let d_vi: Vec<f64> = (0..3)
            .map(|k| (up.vi[k] - down.vi[k]) / (2.0 * delta))
            .collect();
        (d_vr, d_vi)
    };
    let (coarse, fine) = (central(1e-3), central(1e-4));
    // Richardson over δ and δ/10: error O(δ⁴)
    let rich = |c: f64, f: f64| (100.0 * f - c) / 99.0;
    let fd_vr = rich(coarse.0, fine.0);
    let tol = |x: f64| 1e-6 * x.abs().max(1.0);
    assert!(
        (fd_vr - d_vr).abs() < tol(d_vr),
        "FORM {form}: dV_R/dλ: finite difference {fd_vr:.10} vs booked {d_vr:.10}"
    );
    for (k, &booked) in d_vi.iter().enumerate() {
        let fd = rich(coarse.1[k], fine.1[k]);
        assert!(
            (fd - booked).abs() < tol(booked),
            "FORM {form}: dV_{k}/dλ: finite difference {fd:.10} vs booked {booked:.10}"
        );
    }
    let mixed: f64 = w.iter().zip(&d_vi).map(|(p, d)| p * d).sum();
    assert!((mixed - d_vr).abs() < 1e-12 * d_vr.abs().max(1.0));
    println!(
        "FORM {form}: dV_R/dλ booked {d_vr:.10}, finite difference {fd_vr:.10} (Δ {:.2e})",
        (fd_vr - d_vr).abs()
    );
    std::fs::remove_file(ptp).ok();
}

#[test]
fn single_s_books_the_lambda_derivative_of_v_r() {
    check_form(1, vec![0.3], Vec::new(), [3.0, 1.0, -1.5]);
}

#[test]
fn multi_s_books_the_lambda_derivative_of_v_r() {
    check_form(2, vec![0.3, 0.2, 0.1], Vec::new(), [9.0, 1.0, -1.5]);
}

#[test]
fn pair_s_books_the_lambda_derivative_of_v_r() {
    let pair = |i, j, s| EdsStatePair { i, j, s };
    check_form(
        3,
        Vec::new(),
        vec![pair(1, 2, 0.3), pair(2, 3, 0.1)],
        [12.0, 1.0, -1.5],
    );
}
