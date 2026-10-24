//! The long-range solvent–solvent list may be stored as one sentinel pair per
//! molecule pair and evaluated with a single shared periodic shift (the way the
//! short-range solvent list always is) — but only when
//! `sentinel_long_range_is_exact` says the shared shift reproduces the per-atom
//! minimum image for every atom pair. This test checks that claim on a real
//! solvated box: forces, energies and virial from the sentinel path must equal
//! those from the fully expanded per-atom-image path.

use newen_core::configuration::{Box as SimBox, Configuration};
use newen_core::math::{Rectangular, Vec3};
use newen_core::pairlist::{
    sentinel_long_range_is_exact, PairlistContainer, StandardPairlistAlgorithm,
};
use newen_forces::nonbonded::{
    lj_crf_innerloop, solvent_innerloop, CRFParameters, EnergyGroups, ForceStorage, LJParamMatrix,
    LJParameters as KernelLJ,
};
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};
use std::path::PathBuf;

fn shared(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../newen-md/tests/gromosXX_references/shared")
        .join(name)
}

#[test]
fn sentinel_long_range_matches_expanded_per_atom_images() {
    let mut topo = build_topology(read_topology_file(shared("ch4_spc.top")).unwrap());
    topo.solvate(999);
    let coords = read_coordinates(shared("ch4_spc.cnf")).unwrap();
    assert_eq!(coords.positions.len(), topo.num_atoms());

    let n = topo.num_atoms();
    let mut conf = Configuration::new(n, 1, 1);
    conf.current_mut().pos = coords.positions.clone();
    let d = coords.box_dims;
    conf.current_mut().box_config = SimBox::rectangular(d.x, d.y, d.z);
    let periodicity = Rectangular::new(Vec3::new(d.x, d.y, d.z));

    // 3.1 nm box: 0.9 nm cutoff is inside the exactness bound (half-box 1.55 >
    // 0.9 + 0.6), 1.4 nm is not. Both facts are part of what is being tested.
    assert!(sentinel_long_range_is_exact(
        &conf.current().box_config,
        0.9
    ));
    assert!(!sentinel_long_range_is_exact(
        &conf.current().box_config,
        1.4
    ));

    let mut pairlist = PairlistContainer::new(0.6, 0.9, 0.0);
    StandardPairlistAlgorithm::new(true).update(&topo, &conf, &mut pairlist, &periodicity);
    assert!(pairlist.solvent_long_is_sentinel);
    assert!(
        pairlist.solvent_long.len() > 10_000,
        "expected many long-range molecule pairs"
    );

    let apm = topo.atoms_per_solvent();
    assert_eq!(apm, 3);
    let charges = topo.charge.clone();
    let iac: Vec<u32> = topo.iac.iter().map(|&t| t as u32).collect();
    // Topology and kernel use distinct `LJParameters` types; convert as the forcefield does.
    let lj_nested: Vec<Vec<KernelLJ>> = topo
        .lj_parameters
        .iter()
        .map(|row| row.iter().map(KernelLJ::from).collect())
        .collect();
    let lj = LJParamMatrix::from_nested(&lj_nested);
    let crf = CRFParameters {
        crf_cut: 0.9,
        crf_2cut3i: 0.5,
        crf_cut3i: 1.0,
        cutoff_sq: 0.81,
    };
    let four_pi_eps_i = 138.9354;

    // Sentinel path: one shared shift per molecule pair, parameter tables.
    let mut sentinel = ForceStorage::new(n);
    solvent_innerloop(
        &conf.current().pos,
        &charges,
        &iac,
        &pairlist.solvent_long,
        &lj,
        &crf,
        &periodicity,
        apm,
        four_pi_eps_i,
        &mut sentinel,
        EnergyGroups::single(),
    );

    // Reference path: every atom pair explicitly, per-atom minimum image.
    let expanded: Vec<(u32, u32)> = pairlist
        .solvent_long
        .iter()
        .flat_map(|&(i0, j0)| {
            (0..apm as u32).flat_map(move |a| (0..apm as u32).map(move |b| (i0 + a, j0 + b)))
        })
        .collect();
    let mut reference = ForceStorage::new(n);
    lj_crf_innerloop(
        &conf.current().pos,
        &charges,
        &iac,
        &expanded,
        &lj,
        &crf,
        &periodicity,
        four_pi_eps_i,
        &mut reference,
        EnergyGroups::single(),
    );

    let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
    assert!(
        rel(sentinel.e_lj, reference.e_lj) < 1e-11,
        "e_lj {} vs {}",
        sentinel.e_lj,
        reference.e_lj
    );
    assert!(
        rel(sentinel.e_crf, reference.e_crf) < 1e-11,
        "e_crf {} vs {}",
        sentinel.e_crf,
        reference.e_crf
    );
    let mut max_df = 0.0f64;
    for i in 0..n {
        max_df = max_df.max(
            (sentinel.forces[i] - reference.forces[i])
                .abs()
                .max_element(),
        );
    }
    assert!(max_df < 1e-8, "max force deviation {max_df}");
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                rel(sentinel.virial[a][b], reference.virial[a][b]) < 1e-9,
                "virial[{a}][{b}]"
            );
        }
    }
}
