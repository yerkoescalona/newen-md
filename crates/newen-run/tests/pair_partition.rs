//! The MPI pair decomposition without MPI: the shares of `size` ranks, each built the way `md`
//! builds them under `mpirun` (`RunOptions::pair_partition`), sum to the unpartitioned pair
//! terms — forces, energies, virial. A reducer that only records its storage stands in for
//! the cross-rank sum. (`md` under `mpirun -np 2/4` reproduces the serial energies to all
//! printed digits — BENCHMARKING.md §6.)

use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use newen_core::algorithm::SimulationState;
use newen_forces::nonbonded::ForceStorage;
use newen_integrators::algorithms::Forcefield;
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};
use newen_run::{
    build_sequence_from_imd, prepare_system, read_imd, start, Coordinates, RunInputs, RunOptions,
};

fn shared() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../newen-md/tests/gromosXX_references")
}

/// The pair-term storage a rank would hand to MPI, for `partition`.
fn pair_terms(partition: Option<(usize, usize)>) -> ForceStorage {
    let dir = shared().join("aladip_solvated");
    let imd = read_imd(dir.join("aladip_solvated.in")).unwrap();
    let parsed = read_topology_file(shared().join("shared/aladip.topo")).unwrap();
    let constants = parsed.physical_constants;
    let topo = build_topology(parsed);
    let coords: Coordinates = read_coordinates(shared().join("shared/aladip.conf"))
        .unwrap()
        .into();
    let inputs = RunInputs::default();
    let mut prepared = prepare_system(&imd, topo, constants, coords, &inputs).unwrap();
    let options = RunOptions {
        pair_partition: partition,
        ..Default::default()
    };
    let mut built = build_sequence_from_imd(&imd, &prepared, &inputs, &options).unwrap();
    let captured: Arc<Mutex<Option<ForceStorage>>> = Arc::new(Mutex::new(None));
    let slot = Arc::clone(&captured);
    built
        .sequence
        .find_mut::<Forcefield>()
        .unwrap()
        .set_nonbonded_reducer(Box::new(move |s: &mut ForceStorage| {
            *slot.lock().unwrap() = Some(s.clone());
        }));
    let state = SimulationState::new(imd.dt, imd.nstlim);
    start(
        &mut built.sequence,
        &prepared.topology,
        &mut prepared.configuration,
        &state,
    )
    .unwrap();
    let got = captured.lock().unwrap().take().expect("the reducer ran");
    got
}

#[test]
fn rank_shares_sum_to_the_full_pair_terms() {
    let full = pair_terms(None);
    assert!(
        full.e_lj != 0.0 && full.e_crf != 0.0,
        "the system has pair terms"
    );
    for size in [2usize, 3, 5] {
        let mut sum = ForceStorage::with_groups(full.forces.len(), full.groups.n_groups);
        let mut nonempty = 0;
        for rank in 0..size {
            let share = pair_terms(Some((rank, size)));
            if share.e_lj != 0.0 {
                nonempty += 1;
            }
            sum.merge(&share);
        }
        assert!(nonempty >= 2, "size {size}: the work is actually split");
        let tol = 1e-9;
        assert!(
            (sum.e_lj - full.e_lj).abs() <= tol * full.e_lj.abs(),
            "size {size}: e_lj {} vs {}",
            sum.e_lj,
            full.e_lj
        );
        assert!(
            (sum.e_crf - full.e_crf).abs() <= tol * full.e_crf.abs(),
            "size {size}: e_crf {} vs {}",
            sum.e_crf,
            full.e_crf
        );
        for (i, (a, b)) in sum.forces.iter().zip(&full.forces).enumerate() {
            assert!(
                (*a - *b).length() <= 1e-8 * (1.0 + b.length()),
                "size {size}: force on atom {i}: {a:?} vs {b:?}"
            );
        }
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (sum.virial[a][b] - full.virial[a][b]).abs()
                        <= 1e-8 * (1.0 + full.virial[a][b].abs()),
                    "size {size}: virial"
                );
            }
        }
        // the energy-group-pair matrices travel with the shares and add up the same way
        assert!(
            full.groups.is_grouped(),
            "the system has several energy groups"
        );
        for (m, (a, b)) in [
            (&sum.groups.lj, &full.groups.lj),
            (&sum.groups.crf, &full.groups.crf),
        ]
        .into_iter()
        .enumerate()
        {
            for (k, (x, y)) in a.iter().zip(b).enumerate() {
                assert!(
                    (x - y).abs() <= 1e-9 * (1.0 + y.abs()),
                    "size {size}: group matrix {m} cell {k}: {x} vs {y}"
                );
            }
        }
    }
}
