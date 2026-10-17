//! Behaviours that the lift (3.9 step 1) made shared — each one used to exist in only
//! one of the two builders, or in neither. Fixtures are the gromosXX reference systems.

use std::path::{Path, PathBuf};

use gromos_integrators::constraints::{NtcMode, ShakeBuffers};
use gromos_io::{
    coordinate::read_coordinates,
    imd::{read_imd_file, ImdParameters},
    ptp::read_pttopo,
    topology::{build_topology, read_topology_file},
};
use gromos_run::{
    build_sequence_from_imd, prepare_system, total_dof, ConstraintSelection, Coordinates,
    ParallelPolicy, PrepareNote, RunError, RunInputs, RunOptions,
};

fn refs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
}

fn load(system: &str, topo: &str, conf: &str, input: &str) -> Option<Fixture> {
    let dir = refs().join(system);
    if !dir.exists() {
        eprintln!("SKIP: {} not found", dir.display());
        return None;
    }
    let parsed = read_topology_file(dir.join(topo)).expect("topology");
    let physical_constants = parsed.physical_constants;
    let topology = build_topology(parsed);
    let coords: Coordinates = read_coordinates(dir.join(conf))
        .expect("coordinates")
        .into();
    let imd = read_imd_file(dir.join(input)).expect("imd");
    Some(Fixture {
        topology,
        physical_constants,
        coords,
        imd,
    })
}

struct Fixture {
    topology: gromos_core::Topology,
    physical_constants: gromos_core::units::PhysicalConstants,
    coords: Coordinates,
    imd: ImdParameters,
}

#[test]
fn parallel_policy_resolves_like_the_binary() {
    assert!(!ParallelPolicy::Auto.resolve(100));
    assert!(ParallelPolicy::Auto.resolve(101));
    assert!(!ParallelPolicy::Serial.resolve(100_000));
    assert!(ParallelPolicy::Parallel.resolve(2));
    assert_eq!(ParallelPolicy::default(), ParallelPolicy::Auto);
}

/// The coordinate file decides the solvent count; the parameter file's NSM is only a hint
/// (what `md.rs` always did and the Python binding did not).
#[test]
fn nsm_comes_from_the_coordinate_file() {
    let Some(mut fx) = load(
        "nacl_water_box",
        "nacl_water_box.topo",
        "nacl_water_box.conf",
        "nacl_water_box.in",
    ) else {
        return;
    };
    assert!(fx.imd.nsm > 0, "fixture must be a solvated system");
    let true_nsm = fx.imd.nsm;
    let n_coords = fx.coords.positions.len();
    fx.imd.nsm = 1; // wrong on purpose

    let prepared = prepare_system(
        &fx.imd,
        fx.topology,
        fx.physical_constants,
        fx.coords,
        &RunInputs::default(),
    )
    .expect("prepare_system");

    assert_eq!(prepared.topology.num_atoms(), n_coords);
    assert_eq!(prepared.topology.num_solvent_molecules(), true_nsm);
    assert!(prepared.notes.iter().any(|n| matches!(
        n,
        PrepareNote::NsmAdjusted { imd: 1, coordinates } if *coordinates == true_nsm
    )));
}

/// An already-solvated topology (the Python object path) is left alone.
#[test]
fn presolvated_topology_is_not_solvated_again() {
    let Some(mut fx) = load(
        "nacl_water_box",
        "nacl_water_box.topo",
        "nacl_water_box.conf",
        "nacl_water_box.in",
    ) else {
        return;
    };
    fx.topology.solvate(fx.imd.nsm);
    let n_before = fx.topology.num_atoms();
    let prepared = prepare_system(
        &fx.imd,
        fx.topology,
        fx.physical_constants,
        fx.coords,
        &RunInputs::default(),
    )
    .expect("prepare_system");
    assert_eq!(prepared.topology.num_atoms(), n_before);
    assert!(!prepared
        .notes
        .iter()
        .any(|n| matches!(n, PrepareNote::NsmAdjusted { .. })));
}

#[test]
fn atom_count_mismatch_is_an_error_not_a_panic() {
    let Some(mut fx) = load(
        "aladip_vacuum",
        "../shared/aladip.topo",
        "aladip_vacuum.conf",
        "aladip_vacuum.in",
    ) else {
        return;
    };
    fx.coords.positions.pop();
    let err = prepare_system(
        &fx.imd,
        fx.topology,
        fx.physical_constants,
        fx.coords,
        &RunInputs::default(),
    )
    .err()
    .expect("must fail");
    assert!(matches!(err, RunError::AtomCountMismatch { .. }), "{err}");
}

/// One DOF formula: solvent constraints only when the solvent is constrained, and — new in
/// the lift — the solute's own constraints when the solute is constrained.
#[test]
fn total_dof_counts_solute_constraints() {
    let Some(fx) = load(
        "aladip_vacuum_lincs",
        "../shared/aladip.topo",
        "aladip_vacuum_lincs.conf",
        "aladip_vacuum_lincs.in",
    ) else {
        return;
    };
    let topo = &fx.topology;
    let n = topo.num_atoms();

    // NTC=1: nothing constrained.
    let mut imd = fx.imd.clone();
    imd.ntc = 1;
    let sel = ConstraintSelection::from_imd(&imd, false);
    assert_eq!(
        total_dof(topo, &sel, ConstraintSelection::ntc_mode(&imd), imd.ndfmin),
        (3 * n) as f64 - imd.ndfmin as f64
    );

    // NTC=3 (all solute bonds): 3N − every solute bond − NDFMIN.
    let mut imd = fx.imd.clone();
    imd.ntc = 3;
    let sel = ConstraintSelection::from_imd(&imd, false);
    assert!(sel.solute_constrained());
    let n_bonds = ShakeBuffers::new(topo, NtcMode::AllBonds, false)
        .solute_constraints
        .len();
    assert!(n_bonds > 0);
    assert_eq!(
        total_dof(topo, &sel, ConstraintSelection::ntc_mode(&imd), imd.ndfmin),
        (3 * n - n_bonds) as f64 - imd.ndfmin as f64
    );

    // NTC=2 (H-bonds only) subtracts fewer than NTC=3.
    let mut imd = fx.imd.clone();
    imd.ntc = 2;
    let sel = ConstraintSelection::from_imd(&imd, false);
    let dof_h = total_dof(topo, &sel, ConstraintSelection::ntc_mode(&imd), imd.ndfmin);
    assert!(dof_h > (3 * n - n_bonds) as f64 - imd.ndfmin as f64);
    assert!(dof_h < (3 * n) as f64 - imd.ndfmin as f64);
}

/// A perturbation topology given with NTG=0 is ignored (the binary's behaviour); a topology
/// that already carries one but runs with NTG=0 is rejected (only reachable from Python).
#[test]
fn perturbation_topology_follows_ntg() {
    let Some(fx) = load(
        "ch4_water_fep",
        "../shared/ch4_spc.top",
        "../shared/ch4_spc.cnf",
        "ch4_water_fep.in",
    ) else {
        return;
    };
    let ptp = refs().join("shared/ch4_spc_dummy.ptp");
    assert_ne!(fx.imd.ntg, 0, "fixture must be an FEP run");

    // NTG != 0 with @pttopo: merged, reported.
    let inputs = RunInputs {
        pttopo: Some(ptp.clone()),
        ..Default::default()
    };
    let prepared = prepare_system(
        &fx.imd,
        fx.topology.clone(),
        fx.physical_constants,
        fx.coords.clone(),
        &inputs,
    )
    .expect("prepare_system");
    assert!(!prepared.topology.perturbed_solute.atoms.is_empty());
    assert!(prepared
        .notes
        .iter()
        .any(|n| matches!(n, PrepareNote::PerturbationLoaded { .. })));
    // …and the sequence builds with the FEP fields set (no panic, no error).
    build_sequence_from_imd(&fx.imd, &prepared, &inputs, &RunOptions::default())
        .expect("build_sequence_from_imd");

    // NTG = 0 with @pttopo: ignored, noted, topology untouched.
    let mut imd0 = fx.imd.clone();
    imd0.ntg = 0;
    let prepared = prepare_system(
        &imd0,
        fx.topology.clone(),
        fx.physical_constants,
        fx.coords.clone(),
        &inputs,
    )
    .expect("prepare_system");
    assert!(prepared.topology.perturbed_solute.atoms.is_empty());
    assert!(prepared
        .notes
        .iter()
        .any(|n| matches!(n, PrepareNote::PerturbationIgnored { .. })));

    // A pre-perturbed topology with NTG = 0: inconsistent.
    let pt = read_pttopo(&ptp).expect("ptp");
    let mut topo = fx.topology.clone();
    topo.apply_perturbation(pt.perturbed_solute, pt.is_perturbed)
        .unwrap();
    let err = prepare_system(
        &imd0,
        topo,
        fx.physical_constants,
        fx.coords.clone(),
        &RunInputs::default(),
    )
    .err()
    .expect("must fail");
    assert!(matches!(err, RunError::Inconsistent(_)), "{err}");
}

/// The parameter file asks for a restraint file that was not given: a named error.
#[test]
fn missing_restraint_input_is_named() {
    let Some(fx) = load(
        "nacl_1water_distres",
        "../nacl_1water_box/nacl_1water_box.topo",
        "../nacl_1water_box/nacl_1water_box.conf",
        "nacl_1water_distres.in",
    ) else {
        return;
    };
    assert_ne!(fx.imd.ntdir, 0);
    let prepared = prepare_system(
        &fx.imd,
        fx.topology,
        fx.physical_constants,
        fx.coords,
        &RunInputs::default(),
    )
    .expect("prepare_system");
    let err = build_sequence_from_imd(
        &fx.imd,
        &prepared,
        &RunInputs::default(),
        &RunOptions::default(),
    )
    .err()
    .expect("must fail without @distrest");
    assert!(
        matches!(
            err,
            RunError::MissingInput {
                flag: "distrest",
                ..
            }
        ),
        "{err}"
    );
}
