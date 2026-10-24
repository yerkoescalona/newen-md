//! ** P2.8-1 exit criterion, updated for P2.8-4.** `zones.rs` (P2.7 Step 4) states the
//! double-counting contract; `zone_partition_reference.rs` (P2.7 Step 5) proves the ownership
//! table itself partitions real pairs with no gaps or overlaps. Neither, until P2.8-1, proved
//! anything about `LjCrfInteraction` — the double-counting contract was a table nothing consulted.
//!
//! # Why this isn't `partitioned + excluded == full` anymore
//!
//! P2.8-1's original version of this test partitioned every real pair into exactly two classes —
//! classical (`owner() == Classical`, `LjCrfInteraction` evaluates it in full) and excluded
//! (everything else, assumed zero contribution from `LjCrfInteraction`) — and checked
//! `partitioned + excluded(full LJ+CRF) == full`. **P2.8-4 broke that second assumption on
//! purpose**: gromosXX's real `QMMM_Interaction::modify_exclusions` computes classical LJ for
//! inner-outer pairs *unconditionally* (only their electrostatics moves to embedding), so
//! `LjCrfInteraction` now adds an LJ-only supplement for exactly those pairs
//! (`ZonePartition::lj_only_should_evaluate`, `interaction.rs`). Naively summing the *full*
//! LJ+CRF of the "excluded" set on top of that double-counts the LJ portion.
//!
//! The correct three-way decomposition, since every real pair is now in exactly one of three
//! disjoint sets:
//! - **classical** (`owner() == Classical`): full LJ+CRF, in `partitioned.energy`.
//! - **lj-only-supplement** (`lj_only_should_evaluate`): LJ only (CRF zeroed), also already in
//!   `partitioned.energy` — but their CRF is still real physics that has to appear *somewhere*
//!   to reconstruct the unpartitioned baseline.
//! - **fully excluded** (neither of the above — provider-only, e.g. inner-inner/inner-buffer
//!   pairs when `QMLJ` is off): zero contribution from `LjCrfInteraction`, entirely the QM/ML
//!   provider's business.
//!
//! So: `full == partitioned.energy + Σ(lj-only-supplement pairs' CRF) + Σ(fully-excluded pairs'
//! full LJ+CRF)`. That's what this test now checks — still a pure pairwise-sum identity over
//! disjoint sets, still a statement about the *filter and the two-pass wiring*, not a claim about
//! any QM/ML provider's physics.

use newen_core::configuration::{Box as SimBox, Configuration};
use newen_core::math::{Periodicity, Rectangular};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::{ConfigurationSpatialIndex, SpatialIndex};
use newen_forces::nonbonded::{
    lj_crf_innerloop_novirial, CRFParameters, EnergyGroups, ForceStorage, LjCrfInteraction,
};
use newen_forces::provider::PotentialProvider;
use newen_forces::zones::ZonePartition;
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};
use std::collections::HashSet;
use std::path::PathBuf;

/// Real `QMZONE` from `md_burnn/meoh.qmm`: the 6 methanol atoms.
const QM_ZONE_SIZE: usize = 6;
/// Real `BUFFERZONE` radius from `md_burnn/meoh.qmm`'s header line (`0 1 0.5`), in nm.
const BUFFER_RADIUS: f64 = 0.5;
/// Cutoff for this test's classical pairlist — deliberately the same `RCUTQM` from
/// `md_burnn/md.imd` used elsewhere in this suite, so the pair set matches
/// `zone_partition_reference.rs` exactly.
const CUTOFF: f64 = 1.4;

fn tutorial_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../.local/gromos_tutorial_livecoms/tutorial_files/t_06")
}

#[test]
fn partitioned_plus_excluded_reconstructs_full_system_energy() {
    let root = tutorial_root();
    let (topo_path, cnf_path) = (
        root.join("topo/meoh_54a7.top"),
        root.join("eq/eq_meoh_5.cnf"),
    );
    if !topo_path.exists() || !cnf_path.exists() {
        eprintln!(
            "skipping: real BuRNN tutorial not cloned at {}",
            root.display()
        );
        return;
    }

    let mut topo = build_topology(read_topology_file(&topo_path).expect("read t_06 topology"));
    let coords = read_coordinates(&cnf_path).expect("read t_06 eq_meoh_5.cnf");
    let n_atoms = coords.positions.len();
    let n_solute = topo.num_solute_atoms();
    let atoms_per_solvent = topo.solvent_atom_template.len();
    topo.solvate((n_atoms - n_solute) / atoms_per_solvent);

    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions.clone();
    // `LjCrfInteraction::contribute` derives its own periodicity from `conf.current().box_config`
    // (`periodicity_of`, interaction.rs) rather than taking one as a parameter.
    // `Configuration::new` otherwise defaults to `Box::vacuum()`, which would silently make the
    // classical provider use unwrapped distances while `index` below correctly wraps at the box
    // edge — corrupting every pair whose minimum image crosses a boundary. Must match
    // `coords.box_dims` exactly, or this test fails for reasons that have nothing to do with the
    // zone-partition filter it's meant to check (caught empirically: an early draft without this
    // line reconstructed full-system energy off by ~3.7e3 kJ/mol out of ~8e5).
    conf.current_mut().box_config =
        SimBox::rectangular(coords.box_dims.x, coords.box_dims.y, coords.box_dims.z);
    let periodicity = Periodicity::Rectangular(Rectangular::new(coords.box_dims));
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

    // Same inner/buffer construction as zone_partition_reference.rs: whole solvent molecules
    // within BUFFER_RADIUS of the QM zone become buffer, everything else stays outer.
    let inner: Vec<usize> = (0..QM_ZONE_SIZE).collect();
    let inner_sel = AtomSelection::from_indices(inner.clone(), n_atoms).unwrap();

    let mut buffer_molecule = vec![false; (n_atoms - n_solute) / atoms_per_solvent];
    for (i, j, _) in index.neighbor_pairs(&inner_sel, BUFFER_RADIUS) {
        for atom in [i, j] {
            if atom >= n_solute {
                buffer_molecule[(atom - n_solute) / atoms_per_solvent] = true;
            }
        }
    }
    let buffer: Vec<usize> = buffer_molecule
        .iter()
        .enumerate()
        .filter(|(_, &b)| b)
        .flat_map(|(m, _)| {
            let start = n_solute + m * atoms_per_solvent;
            start..start + atoms_per_solvent
        })
        .collect();
    assert!(!buffer.is_empty());
    let partition = ZonePartition::new(n_atoms, &inner, &buffer);

    let region = AtomSelection::all(n_atoms);

    // ── Full, unpartitioned classical energy: today's default (no ZonePartition set) ─────
    let mut unpartitioned = LjCrfInteraction::new(CUTOFF, 1.0, 66.6, 0.0);
    let full = unpartitioned
        .contribute(&region, &topo, &conf, &index)
        .unwrap();

    // ── Partitioned classical energy: Provider/Embedding-owned pairs skipped ──────────────
    let mut partitioned_interaction =
        LjCrfInteraction::new(CUTOFF, 1.0, 66.6, 0.0).with_zone_partition(partition.clone());
    let partitioned = partitioned_interaction
        .contribute(&region, &topo, &conf, &index)
        .unwrap();

    // ── The three-way split (P2.8-4) ───────────────────────────────────────────────────────
    let all_pairs: Vec<(usize, usize, newen_core::math::Vec3)> =
        index.neighbor_pairs(&region, CUTOFF);
    let classical_pairs: Vec<(u32, u32)> = all_pairs
        .iter()
        .filter(|(i, j, _)| partition.classical_should_evaluate(*i, *j))
        .map(|(i, j, _)| (*i as u32, *j as u32))
        .collect();
    let lj_only_pairs: Vec<(u32, u32)> = all_pairs
        .iter()
        .filter(|(i, j, _)| partition.lj_only_should_evaluate(*i, *j))
        .map(|(i, j, _)| (*i as u32, *j as u32))
        .collect();
    let fully_excluded_pairs: Vec<(u32, u32)> = all_pairs
        .iter()
        .filter(|(i, j, _)| {
            !partition.classical_should_evaluate(*i, *j)
                && !partition.lj_only_should_evaluate(*i, *j)
        })
        .map(|(i, j, _)| (*i as u32, *j as u32))
        .collect();
    assert!(
        !lj_only_pairs.is_empty() && !fully_excluded_pairs.is_empty(),
        "a real three-zone system must exercise both the LJ-only supplement (inner-outer pairs) \
         and full exclusion (inner-inner/inner-buffer pairs, QMLJ off), or this test proves \
         nothing about either"
    );

    // Sanity: the three sets are a genuine partition of every pair LjCrfInteraction itself would
    // find — no pair lost, none claimed twice.
    let classical_set: HashSet<_> = classical_pairs.iter().collect();
    let lj_only_set: HashSet<_> = lj_only_pairs.iter().collect();
    let excluded_set: HashSet<_> = fully_excluded_pairs.iter().collect();
    assert_eq!(
        classical_pairs.len() + lj_only_pairs.len() + fully_excluded_pairs.len(),
        all_pairs.len()
    );
    assert!(classical_set.is_disjoint(&lj_only_set));
    assert!(classical_set.is_disjoint(&excluded_set));
    assert!(lj_only_set.is_disjoint(&excluded_set));

    let charges = topo.charge.clone();
    let iac_u32: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let lj_nested: Vec<Vec<newen_forces::nonbonded::LJParameters>> = topo
        .lj_parameters
        .iter()
        .map(|row| {
            row.iter()
                .map(newen_forces::nonbonded::LJParameters::from)
                .collect()
        })
        .collect();
    let lj_mat = newen_forces::nonbonded::LJParamMatrix::from_nested(&lj_nested);
    let crf = CRFParameters::new(CUTOFF, 1.0, 66.6, 0.0);

    // Σ(lj-only-supplement pairs' full LJ+CRF, real charges) — this is what those pairs would
    // contribute if they were ordinary classical pairs, i.e. what "full" already counts them as.
    let mut lj_only_full_storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop_novirial(
        &conf.current().pos,
        &charges,
        &iac_u32,
        &lj_only_pairs,
        &lj_mat,
        &crf,
        &periodicity,
        newen_core::units::four_pi_eps_i,
        &mut lj_only_full_storage,
        EnergyGroups::single(),
    );
    let lj_only_full_energy = lj_only_full_storage.e_lj + lj_only_full_storage.e_crf;

    // Σ(lj-only-supplement pairs' LJ only) — exactly what `partitioned.energy` already includes
    // for these pairs (mirrors `contribute()`'s own zeroed-inner-charge trick).
    let mut lj_only_charges = charges.clone();
    for atom in partition.atoms_in(newen_forces::zones::Zone::Inner) {
        lj_only_charges[atom] = 0.0;
    }
    let mut lj_only_supplement_storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop_novirial(
        &conf.current().pos,
        &lj_only_charges,
        &iac_u32,
        &lj_only_pairs,
        &lj_mat,
        &crf,
        &periodicity,
        newen_core::units::four_pi_eps_i,
        &mut lj_only_supplement_storage,
        EnergyGroups::single(),
    );
    let lj_only_supplement_energy =
        lj_only_supplement_storage.e_lj + lj_only_supplement_storage.e_crf;

    // The CRF remainder for lj-only-supplement pairs: real physics not yet accounted for by
    // `partitioned.energy` (which only has their LJ).
    let lj_only_crf_remainder = lj_only_full_energy - lj_only_supplement_energy;

    // Σ(fully-excluded pairs' full LJ+CRF) — zero contribution from `LjCrfInteraction`, entirely
    // provider-owned; this is the leftover physics the QM/ML provider would supply in a real run.
    let mut excluded_storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop_novirial(
        &conf.current().pos,
        &charges,
        &iac_u32,
        &fully_excluded_pairs,
        &lj_mat,
        &crf,
        &periodicity,
        newen_core::units::four_pi_eps_i,
        &mut excluded_storage,
        EnergyGroups::single(),
    );
    let fully_excluded_energy = excluded_storage.e_lj + excluded_storage.e_crf;

    // ── The exit criterion ─────────────────────────────────────────────────────────────
    let reconstructed = partitioned.energy + lj_only_crf_remainder + fully_excluded_energy;
    assert!(
        (reconstructed - full.energy).abs() < 1e-8 * full.energy.abs().max(1.0),
        "partitioned ({}) + lj-only CRF remainder ({}) + fully-excluded ({}) = {} != full ({}); \
         the filter or the two-pass wiring dropped or double-counted a pair",
        partitioned.energy,
        lj_only_crf_remainder,
        fully_excluded_energy,
        reconstructed,
        full.energy
    );
    assert!(
        partitioned.energy != full.energy,
        "sanity: partitioning must actually remove some energy on a real three-zone system"
    );

    eprintln!(
        "P2.8-1/P2.8-4 hold on real data: full={:.6}, partitioned={:.6} (classical {} pairs + \
         lj-only-supplement {} pairs), fully-excluded={:.6} kJ/mol over {} pairs, of {} total \
         within {CUTOFF} nm.",
        full.energy,
        partitioned.energy,
        classical_pairs.len(),
        lj_only_pairs.len(),
        fully_excluded_energy,
        fully_excluded_pairs.len(),
        all_pairs.len(),
    );
}
