//! `ProviderOrchestrator` running two real providers together on the real `t_06` system —
//! `LjCrfInteraction` (zone-partitioned) and `ElectrostaticEmbedding` — the classical and
//! embedding halves of a QM/MM step. Earlier orchestrator tests proved transparency with exactly
//! one provider on a plain two-atom fixture; this is the first with two providers on a genuinely
//! partitioned, real system.
//!
//! # Scope, precisely
//!
//! A full QM/MM step needs three energy sources: the classical field, the QM↔MM embedding, and
//! the QM/ML provider's own energy for the inner zone. This environment has no QM engine (no
//! `xtb`/`mopac` binary) and no `libtorch` (`SchNetInteraction` needs `--features ml`, confirmed
//! unbuildable here — `cargo build -p newen-forces --features ml` fails with "Cannot find a
//! libtorch install"), so the third source is unavailable and not stood in for. This test
//! registers the two sources that *are* real and available, and checks their combination is
//! correct — not that the result is a complete QM/MM energy.

use newen_core::configuration::{Box as SimBox, Configuration};
use newen_core::math::{Periodicity, Rectangular};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::{ConfigurationSpatialIndex, SpatialIndex};
use newen_forces::nonbonded::{ElectrostaticEmbedding, LjCrfInteraction};
use newen_forces::orchestrator::ProviderOrchestrator;
use newen_forces::provider::PotentialProvider;
use newen_forces::zones::ZonePartition;
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};
use std::path::PathBuf;

const QM_ZONE_SIZE: usize = 6;
const BUFFER_RADIUS: f64 = 0.5;
const RCUTQM: f64 = 1.4;

fn tutorial_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../.local/gromos_tutorial_livecoms/tutorial_files/t_06")
}

#[test]
fn orchestrator_combines_classical_and_embedding_on_the_real_system() {
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
    conf.current_mut().box_config =
        SimBox::rectangular(coords.box_dims.x, coords.box_dims.y, coords.box_dims.z);
    let periodicity = Periodicity::Rectangular(Rectangular::new(coords.box_dims));
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

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

    // Stand-in QM charges (same values `embedding_gather_reference.rs` uses) — a real run would
    // refresh these from a QM engine each step; the orchestrator wiring under test is unaffected
    // by where the charges came from.
    let mut qm_charges = vec![0.0; n_atoms];
    for (i, q) in [0.27, -0.64, 0.04, 0.04, 0.04, 0.41].iter().enumerate() {
        qm_charges[i] = *q;
    }

    // ── Two independent, direct calls — the oracle ─────────────────────────────────────────
    let mut classical_alone =
        LjCrfInteraction::new(RCUTQM, 1.0, 66.6, 0.0).with_zone_partition(partition.clone());
    let classical_region = AtomSelection::all(n_atoms);
    let classical_direct = classical_alone
        .contribute(&classical_region, &topo, &conf, &index)
        .unwrap();

    let mut embedding_alone = ElectrostaticEmbedding::new(RCUTQM, qm_charges.clone());
    let embedding_direct = embedding_alone
        .contribute(&inner_sel, &topo, &conf, &index)
        .unwrap();

    // ── The same two providers, through the orchestrator ───────────────────────────────────
    let mut orchestrator = ProviderOrchestrator::new();
    orchestrator.register(
        AtomSelection::all(n_atoms),
        Box::new(LjCrfInteraction::new(RCUTQM, 1.0, 66.6, 0.0).with_zone_partition(partition)),
    );
    orchestrator.register(
        inner_sel,
        Box::new(ElectrostaticEmbedding::new(RCUTQM, qm_charges)),
    );
    let combined = orchestrator
        .evaluate(&topo, &conf, &index)
        .expect("orchestrator should combine both providers without an index-contract violation");

    // Transparency, extended to two providers: the orchestrator's sum must equal the two direct
    // calls added by hand — no interaction between providers beyond addition.
    let expected_energy = classical_direct.energy + embedding_direct.energy;
    assert!(
        (combined.energy - expected_energy).abs() < 1e-8 * expected_energy.abs().max(1.0),
        "combined {} vs classical {} + embedding {} = {}",
        combined.energy,
        classical_direct.energy,
        embedding_direct.energy,
        expected_energy
    );

    // Forces from both sources must all be present: LjCrfInteraction's classical/lj-only-
    // supplement atoms, and ElectrostaticEmbedding's forces reaching MM atoms outside its own
    // `region` — exactly the case the orchestrator's index check must let through rather than
    // reject.
    assert!(
        combined.forces.len() >= embedding_direct.forces.len(),
        "combined forces must be at least as many atoms as embedding touched on its own"
    );
    let mm_embedding_atoms = embedding_direct
        .forces
        .iter()
        .filter(|&&(i, _)| i >= QM_ZONE_SIZE)
        .count();
    assert!(
        mm_embedding_atoms > 0,
        "sanity: embedding must be reaching real MM atoms for this test to mean anything"
    );

    eprintln!(
        "orchestrator combines classical+embedding correctly on the real system: combined={:.6} \
         (classical={:.6} + embedding={:.6}) kJ/mol, {} total force-bearing atoms; the QM/ML \
         provider's own energy is not part of this (needs a real QM engine, not available here).",
        combined.energy,
        classical_direct.energy,
        embedding_direct.energy,
        combined.forces.len(),
    );
}
