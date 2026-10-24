//! ** P2.7 Step 5 (revised)** — the double-counting contract, checked on the real system.
//!
//! # Why this replaces the originally-planned oracle
//!
//! Step 5 was originally scoped as: recompute the classical buffer energy with newen and
//! compare against "the MM baseline the BuRNN training pipeline used". **That premise is false**,
//! and the tutorial's own source disproves it:
//!
//! - `train_dataset_tutorial/mopac.py::get_burnn_energy()` subtracts two **MOPAC** (QM) energies —
//!   `compl_energy_kj − buffer_energy_kj`. No classical energy appears (assumption **A2**).
//! - Grepping the whole training pipeline for `lj|coulomb|classical|forcefield` returns nothing:
//!   no MM term enters training at any point.
//! - `get_reference_energies()` — the function that would have produced the vacuum references —
//!   is marked `# not finished`, is **never called**, and never returns a value; the references
//!   are hand-supplied through the constructor instead.
//!
//! So there is no MM baseline to reproduce, and inventing one would have manufactured a
//! false-confidence test. P2.7 Step 5 explicitly said to drop it in this case.
//!
//! # What is actually checkable, and is checked here
//!
//! The real question Step 5 was reaching for — *"how do we detect if the two sides' logics
//! differ?"* — has a genuine answer at the level the contract lives: **the ownership table must
//! partition the interactions**. If the three owners overlap, terms get double-counted; if they
//! leave a gap, physics silently vanishes. Both are checkable on real data without any QM
//! reference, and both are exactly the failure modes assumption **A5** says nothing detects at
//! runtime.
//!
//! Fixture is the real BuRNN tutorial system with its real zone parameters (gitignored `.local/`
//! clone; skips gracefully when absent).

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Rectangular};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::{ConfigurationSpatialIndex, SpatialIndex};
use newen_forces::provider::Embedding;
use newen_forces::zones::{PairOwner, Zone, ZonePartition};
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};
use std::path::PathBuf;

/// Real `QMZONE` from `md_burnn/meoh.qmm`: the 6 methanol atoms.
const QM_ZONE_SIZE: usize = 6;
/// Real `BUFFERZONE` radius from `md_burnn/meoh.qmm`'s header line (`0 1 0.5`), in nm.
/// The block lists all 3507 solvent atoms as the *eligible pool*; this radius picks the actual
/// buffer each step, which is why the buffer must be computed here rather than read off.
const BUFFER_RADIUS: f64 = 0.5;
/// Real `RCUTQM` from `md_burnn/md.imd`.
const RCUTQM: f64 = 1.4;

fn tutorial_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../.local/gromos_tutorial_livecoms/tutorial_files/t_06")
}

#[test]
fn ownership_table_partitions_real_burnn_system_without_gaps_or_overlaps() {
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
    let periodicity = Periodicity::Rectangular(Rectangular::new(coords.box_dims));
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

    // ── Build the real zone assignment ───────────────────────────────────────────────────
    // Buffer = solvent molecules with any atom within BUFFER_RADIUS of the QM zone. Whole
    // molecules, not individual atoms: splitting a water across the boundary would leave a
    // fractional charge, the same reason GROMOS uses charge groups as its cutoff unit.
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

    assert!(
        !buffer.is_empty(),
        "a solvated system must have solvent within {BUFFER_RADIUS} nm of the solute"
    );
    let partition = ZonePartition::new(n_atoms, &inner, &buffer);

    // All three zones must be populated, or the test is not exercising what it claims.
    let n_inner = partition.atoms_in(Zone::Inner).len();
    let n_buffer = partition.atoms_in(Zone::Buffer).len();
    let n_outer = partition.atoms_in(Zone::Outer).len();
    assert_eq!(n_inner, QM_ZONE_SIZE);
    assert!(n_buffer > 0 && n_outer > 0, "need a real three-zone system");
    assert_eq!(
        n_inner + n_buffer + n_outer,
        n_atoms,
        "zones must cover every atom exactly once"
    );

    // ── The contract: every interacting pair has exactly one owner ───────────────────────
    let pairs = index.neighbor_pairs(&inner_sel, RCUTQM);
    assert!(!pairs.is_empty());

    let (mut n_provider, mut n_classical, mut n_embedding) = (0usize, 0usize, 0usize);
    for (i, j, _) in &pairs {
        let owner = partition.owner(*i, *j);

        // No overlap: a pair the classical field evaluates must never also be embedded, and
        // neither may touch a pair the QM/ML provider already accounted for. This is the
        // double-counting guard A5 says cannot be detected at runtime — so it is asserted
        // statically here instead.
        let classical = partition.classical_should_evaluate(*i, *j);
        let embedded = partition.embedding_should_evaluate(Embedding::Electrostatic, *i, *j);
        assert!(
            !(classical && embedded),
            "pair ({i},{j}) claimed by both the classical field and embedding"
        );
        if owner == PairOwner::Provider {
            assert!(
                !classical && !embedded,
                "pair ({i},{j}) is provider-owned but also claimed classically/by embedding"
            );
        }

        // No gap: every pair is claimed by exactly one of the three.
        let claims = [owner == PairOwner::Provider, classical, embedded]
            .iter()
            .filter(|&&c| c)
            .count();
        assert_eq!(
            claims, 1,
            "pair ({i},{j}) has {claims} owners, expected exactly 1"
        );

        match owner {
            PairOwner::Provider => n_provider += 1,
            PairOwner::Classical => n_classical += 1,
            PairOwner::Embedding => n_embedding += 1,
        }
    }

    assert_eq!(n_provider + n_classical + n_embedding, pairs.len());
    // A real three-zone system must exercise all three branches, otherwise the partition is
    // trivially consistent and proves nothing.
    assert!(n_provider > 0, "expected inner-inner/inner-buffer pairs");
    assert!(
        n_embedding > 0,
        "expected inner-outer pairs needing embedding"
    );

    eprintln!(
        "Step 5 (revised) holds on real data: {n_atoms} atoms → inner {n_inner} / buffer \
         {n_buffer} / outer {n_outer}; {} QM-cutoff pairs partitioned as provider {n_provider}, \
         classical {n_classical}, embedding {n_embedding} — no gaps, no overlaps.",
        pairs.len()
    );
}
