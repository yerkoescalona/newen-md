//! The EDS forces are −∇V_R (1.8): central differences of the potential plus V_R (the special
//! energy) against the mixed forces the force field and the EDS step produce, for the three
//! forms of the reference Hamiltonian, on the vacuum alanine dipeptide with three end states.
//!
//! 5 atoms × 3 axes, central differences at δ = 1e-4 and δ/10 combined by Richardson
//! extrapolation (error O(δ⁴)); the tolerance is the one of the native check
//! (`gromos-divergences` sims/e2_eds_fd, 1e-3 kJ mol⁻¹ nm⁻¹ through 10-digit output, where δ/10
//! was below the output precision and 1e-3 nm was used instead — with every state weighted, the
//! mixing curvature leaves up to 5e-3 at that pair here). The end-state offsets are set so that
//! all three states carry weight, else one state's force would pass for the mixture.

use std::path::{Path, PathBuf};

use gromos_core::algorithm::SimulationState;
use gromos_core::math::Vec3;
use gromos_integrators::algorithms::Forcefield;
use gromos_integrators::eds::{reference_hamiltonian, Eds};
use gromos_io::coordinate::read_coordinates;
use gromos_io::imd::{read_imd_file, EdsBlock, EdsStatePair};
use gromos_io::topology::{build_topology, read_topology_file};
use gromos_run::{build_sequence_from_imd, prepare_system, Coordinates, RunInputs, RunOptions};

fn refs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
}

/// Potential + special energy and the forces of `positions`, from the force field and the EDS
/// step of `seq` (no integration).
fn evaluate(
    seq: &mut gromos_core::algorithm::AlgorithmSequence,
    topo: &gromos_core::Topology,
    conf: &mut gromos_core::Configuration,
    positions: &[Vec3],
) -> (f64, Vec<Vec3>, Vec<f64>) {
    use gromos_core::algorithm::Algorithm;
    let sim = SimulationState::new(0.002, 1);
    conf.current_mut().pos = positions.to_vec();
    let ff = seq.find_mut::<Forcefield>().expect("force field");
    ff.apply(topo, conf, &sim).unwrap();
    let eds = seq.find_mut::<Eds>().expect("EDS step");
    eds.apply(topo, conf, &sim).unwrap();
    let e = &conf.current().energies;
    (
        e.potential_total + e.special_total,
        conf.current().force.clone(),
        e.eds_vi.clone(),
    )
}

fn check_form(form: i32, s: Vec<f64>, pairs: Vec<EdsStatePair>, offsets: [f64; 3]) {
    let dir = refs().join("aladip_vacuum_eds3");
    let mut imd = read_imd_file(dir.join("aladip_vacuum_eds3.in")).unwrap();
    imd.eds = Some(EdsBlock {
        enabled: true,
        form,
        num_states: 3,
        s,
        pairs,
        eir: vec![0.0; 3],
    });
    let parsed = read_topology_file(refs().join("shared/aladip.topo")).unwrap();
    let constants = parsed.physical_constants;
    let coords: Coordinates = read_coordinates(refs().join("aladip_vacuum/aladip_vacuum.conf"))
        .unwrap()
        .into();
    let inputs = RunInputs {
        pttopo: Some(refs().join("shared/aladip_eds3.mpt")),
        ..Default::default()
    };
    let prepared =
        prepare_system(&imd, build_topology(parsed), constants, coords, &inputs).unwrap();
    let built = build_sequence_from_imd(&imd, &prepared, &inputs, &RunOptions::default()).unwrap();
    let mut seq = built.sequence;
    let topo = prepared.topology;
    let mut conf = prepared.configuration;
    seq.init(&topo, &mut conf, &SimulationState::new(0.002, 1))
        .unwrap();
    let x0 = conf.current().pos.clone();

    // offsets that give every state a weight of at least 5 % at x0
    let (_, _, vi) = evaluate(&mut seq, &topo, &mut conf, &x0);
    let eds = seq.find_mut::<Eds>().unwrap();
    eds.eir = (0..3).map(|k| vi[k] + offsets[k]).collect();
    let mut w = [0.0; 3];
    reference_hamiltonian(&eds.form, eds.beta, &eds.eir, &vi, &mut w);
    assert!(w.iter().all(|&p| p > 0.05), "FORM {form}: weights {w:?}");
    let (_, forces, _) = evaluate(&mut seq, &topo, &mut conf, &x0);

    let mut worst: f64 = 0.0;
    for atom in [0, 3, 6, 8, 10] {
        for axis in 0..3 {
            let mut central = |delta: f64| {
                let mut e = [0.0; 2];
                for (k, sign) in [1.0, -1.0].into_iter().enumerate() {
                    let mut x = x0.clone();
                    x[atom][axis] += sign * delta;
                    e[k] = evaluate(&mut seq, &topo, &mut conf, &x).0;
                }
                -(e[0] - e[1]) / (2.0 * delta)
            };
            let richardson = (100.0 * central(1e-5) - central(1e-4)) / 99.0;
            let diff = (richardson - forces[atom][axis]).abs();
            worst = worst.max(diff);
            assert!(
                diff < 1e-3,
                "FORM {form}: atom {} axis {axis}: −dE/dx = {richardson:.9} vs F = {:.9}",
                atom + 1,
                forces[atom][axis]
            );
        }
    }
    println!("FORM {form}: max |F_FD − F| = {worst:.3e} kJ/mol/nm");
}

#[test]
fn single_s_forces_are_the_gradient_of_v_r() {
    check_form(1, vec![0.3], Vec::new(), [3.0, 1.0, -1.5]);
}

#[test]
fn multi_s_forces_are_the_gradient_of_v_r() {
    check_form(2, vec![0.3, 0.2, 0.1], Vec::new(), [9.0, 1.0, -1.5]);
}

#[test]
fn pair_s_forces_are_the_gradient_of_v_r() {
    let pair = |i, j, s| EdsStatePair { i, j, s };
    check_form(
        3,
        Vec::new(),
        vec![pair(1, 2, 0.3), pair(2, 3, 0.1)],
        [12.0, 1.0, -1.5],
    );
}
