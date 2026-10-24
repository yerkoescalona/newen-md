//! Invariants that `Topology::is_excluded` / `is_excluded_or_14` rely on.
//!
//! Both query a single side of the per-atom lists and stop early once the last
//! (largest) entry is below the partner, which is only correct if every list is
//! sorted and symmetric. These tests load real reference
//! topologies — one with dihedrals (so 1-4 pairs are present) and one that is
//! then solvated — and check the invariants directly, plus agreement with the
//! two-sided formulation the single-sided lookups replaced.

use newen_core::topology::Topology;
use newen_io::topology::{build_topology, read_topology_file};
use std::path::PathBuf;

fn shared(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../newen-md/tests/gromosXX_references/shared")
        .join(name)
}

fn load(name: &str) -> Topology {
    build_topology(read_topology_file(shared(name)).expect("reference topology reads"))
}

fn assert_sorted_and_symmetric(lists: &[Vec<usize>], what: &str) {
    for (i, list) in lists.iter().enumerate() {
        assert!(
            list.windows(2).all(|w| w[0] < w[1]),
            "{what}[{i}] is not strictly sorted: {list:?}"
        );
        for &j in list {
            assert!(
                lists[j].contains(&i),
                "{what} not symmetric: {j} in {what}[{i}] but {i} not in {what}[{j}]"
            );
        }
    }
}

/// The two-sided test the single-sided lookup replaced; kept as the oracle.
fn two_sided_excluded_or_14(topo: &Topology, i: usize, j: usize) -> bool {
    topo.exclusions[i].binary_search(&j).is_ok()
        || topo.exclusions[j].binary_search(&i).is_ok()
        || (!topo.one_four_pairs.is_empty()
            && (topo.one_four_pairs[i].contains(&j) || topo.one_four_pairs[j].contains(&i)))
}

fn check(topo: &Topology) {
    assert_sorted_and_symmetric(&topo.exclusions, "exclusions");
    assert_sorted_and_symmetric(&topo.one_four_pairs, "one_four_pairs");

    let n = topo.num_atoms();
    let mut excluded_pairs = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let expected = two_sided_excluded_or_14(topo, i, j);
            assert_eq!(topo.is_excluded_or_14(i, j), expected, "pair ({i}, {j})");
            assert_eq!(topo.is_excluded_or_14(j, i), expected, "pair ({j}, {i})");
            excluded_pairs += expected as usize;
        }
    }
    assert!(excluded_pairs > 0, "test system has no exclusions at all");
}

#[test]
fn aladip_solute_with_one_four_pairs() {
    let topo = load("aladip.topo");
    assert!(
        topo.one_four_pairs.iter().any(|l| !l.is_empty()),
        "aladip must have 1-4 pairs for this test to mean anything"
    );
    check(&topo);
}

#[test]
fn ch4_spc_solvated_lists_stay_symmetric() {
    let mut topo = load("ch4_spc.top");
    topo.solvate(50);
    assert!(topo.num_atoms() > 150);
    // Solvation must keep every per-atom list addressable for every atom;
    // the atomic-cutoff pairlist queries solvent atoms too.
    assert_eq!(topo.exclusions.len(), topo.num_atoms());
    assert!(
        topo.one_four_pairs.is_empty() || topo.one_four_pairs.len() == topo.num_atoms(),
        "one_four_pairs must be empty or cover all {} atoms, has {}",
        topo.num_atoms(),
        topo.one_four_pairs.len()
    );
    check(&topo);
}
