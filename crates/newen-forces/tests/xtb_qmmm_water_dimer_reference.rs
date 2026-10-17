//! A real, combined QM/MM evaluation — `XtbInteraction` + `ElectrostaticEmbedding` + a
//! zone-partitioned `LjCrfInteraction`, together, on a real system, with real xtb-derived
//! embedding charges. `qmmm_orchestrator_reference.rs` already proved the classical+embedding
//! half combines correctly, but explicitly could not include the QM zone's own energy ("this
//! environment has no QM engine... not that the result is a complete QM/MM energy") and used
//! hard-coded stand-in charges for the embedding term. `xtb` closes both: this is the first test
//! in this repo where all three real pieces run together — a real QM engine, real QM-derived
//! charges, and real classical MM — so the QM<->MM coupling can actually be observed, not just
//! asserted correct at the provider level.
//!
//! **System: `water_dimer` (Poliak dataset), not `t_06`.** `t_06`'s 54a7 force field is
//! united-atom; its QM zone's true per-atom element identity lives in a `.qmm` file's
//! `QMZONE`/`QMZ` column, which nothing here parses yet. `water_dimer` (two explicit-atom SPC
//! waters, checked into `gromosXX_qmmm_references/water_dimer_mechst/`, already used by
//! `water_dimer_qmmm_mechst_reference.rs`) needs no such mapping — every atom is real, atomic
//! numbers are trivially `[8,1,1]` per water.
//!
//! **Charges: read from xtb's own `charges` file, not a new `Contribution` channel.**
//! `provider.rs` explicitly defers a charge-output channel on `Contribution` (it would make
//! providers non-additive — one provider's output feeding another's input within a step). Rather
//! than build that now, this calls `XtbInteraction::contribute()` once (its real return value is
//! the QM zone's own energy and forces), then separately reads the `charges` file xtb already
//! writes into the same `work_dir` as a side effect — a documented stopgap for this one demo, not
//! a new architectural pattern.
//!
//! **What this does *not* claim:** it won't match gromosXX's real published number for this
//! system (their reference run used AM1, not xtb — a known, already-accepted method mismatch,
//! same as elsewhere in this repo); no link atoms (neither water crosses the QM/MM boundary with
//! a covalent bond, so none are needed here). This is single-point, not a dynamical run — the
//! waters are SHAKE-constrained rigid bodies in the real `.imd`, and a real multi-step run would
//! need the constraint solver wired into the sequence too; that's a separate, larger follow-up.
//!
//! **Two embedding paths, both demonstrated here.** The first test below still uses the static-
//! charge path (c) (`ElectrostaticEmbedding`, one xtb call for real Mulliken charges, held fixed)
//! — kept as-is since it's a real, useful demonstration of the orchestrator combining three
//! independent providers. The second test uses `XtbInteraction`'s own `Embedding::Electrostatic`
//! (path (a) — real SCF polarization, real MM forces from xtb's `pcgrad`, no static-charge
//! approximation and no separate embedding provider at all): the more physically correct
//! approach whenever the engine supports it, per `nonbonded/embedding.rs`'s module docs on why
//! the two are alternatives, not additions.
//!
//! Skips gracefully if `xtb` isn't on `PATH`, same as every other xtb test in this repo.

use std::path::PathBuf;

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_forces::nonbonded::{ElectrostaticEmbedding, LjCrfInteraction, XtbInteraction};
use newen_forces::orchestrator::ProviderOrchestrator;
use newen_forces::provider::{Embedding, PotentialProvider};
use newen_forces::zones::ZonePartition;
use newen_io::coordinate::read_coordinates;
use newen_io::imd::read_imd_file;
use newen_io::topology::{build_topology, read_topology_file};

const QM_ZONE_SIZE: usize = 3;

fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/gromosXX_qmmm_references/water_dimer_mechst")
}

fn xtb_available() -> bool {
    std::process::Command::new("xtb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

#[test]
fn combined_qmmm_shows_real_bidirectional_coupling_on_a_real_system() {
    if !xtb_available() {
        eprintln!("skipping: xtb not found on PATH");
        return;
    }

    let root = fixture_root();
    let mut topo = build_topology(
        read_topology_file(root.join("water_dimer.top")).expect("read water_dimer.top"),
    );
    let coords = read_coordinates(root.join("water_dimer.cnf")).expect("read water_dimer.cnf");
    let imd = read_imd_file(root.join("water_dimer.imd")).expect("read water_dimer.imd");

    let n_solute = topo.num_solute_atoms();
    let atoms_per_solvent = topo.solvent_atom_template.len();
    topo.solvate((coords.positions.len() - n_solute) / atoms_per_solvent);

    let n_atoms = topo.num_atoms();
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions.clone();
    let periodicity = Periodicity::Vacuum(Vacuum);
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

    // Solute water (0,1,2) is the QM zone; the one NSM=1 solvent water (3,4,5) is MM.
    let inner: Vec<usize> = (0..QM_ZONE_SIZE).collect();
    let inner_sel = AtomSelection::from_indices(inner.clone(), n_atoms).unwrap();
    let partition = ZonePartition::new(n_atoms, &inner, &[]);

    // ── Real xtb call #1: QM zone's own energy/forces, and (as a side effect) real Mulliken
    // charges for the embedding term ────────────────────────────────────────────────────────
    let charge_work_dir = std::env::temp_dir().join("newen_xtb_qmmm_water_dimer_charges");
    let mut qm_for_charges = XtbInteraction::new(&charge_work_dir, 2, 0, 1, vec![8, 1, 1, 8, 1, 1])
        .expect("work_dir creatable");
    let qm_direct = qm_for_charges
        .contribute(&inner_sel, &topo, &conf, &index)
        .expect("xtb QM-zone calculation should succeed");

    let charges_content = std::fs::read_to_string(charge_work_dir.join("charges"))
        .expect("xtb should have written a charges file as a side effect");
    let mulliken_charges: Vec<f64> = charges_content
        .lines()
        .map(|l| l.trim().parse().expect("charges file: one float per line"))
        .collect();
    assert_eq!(
        mulliken_charges.len(),
        QM_ZONE_SIZE,
        "expected one Mulliken charge per QM-zone atom"
    );
    let mut region_charges = vec![0.0; n_atoms];
    for (local, &global) in inner.iter().enumerate() {
        region_charges[global] = mulliken_charges[local];
    }
    eprintln!("real xtb Mulliken charges for the QM zone: {mulliken_charges:?}");

    // ── The three real providers, evaluated directly (the oracle for the orchestrator sum) ──
    let mut classical = LjCrfInteraction::new(imd.rcrf, 1.0, imd.epsrf, imd.appak)
        .with_zone_partition(partition.clone());
    let classical_region = AtomSelection::all(n_atoms);
    let classical_direct = classical
        .contribute(&classical_region, &topo, &conf, &index)
        .unwrap();

    let mut embedding = ElectrostaticEmbedding::new(imd.rcrf, region_charges.clone());
    let embedding_direct = embedding
        .contribute(&inner_sel, &topo, &conf, &index)
        .unwrap();

    // qm_direct (above) is the third.

    // ── Through the orchestrator ─────────────────────────────────────────────────────────────
    let qm_work_dir = std::env::temp_dir().join("newen_xtb_qmmm_water_dimer_orchestrated");
    let mut orchestrator = ProviderOrchestrator::new();
    orchestrator.register(
        classical_region,
        Box::new(
            LjCrfInteraction::new(imd.rcrf, 1.0, imd.epsrf, imd.appak)
                .with_zone_partition(partition),
        ),
    );
    orchestrator.register(
        inner_sel.clone(),
        Box::new(ElectrostaticEmbedding::new(imd.rcrf, region_charges)),
    );
    orchestrator.register(
        inner_sel,
        Box::new(
            XtbInteraction::new(&qm_work_dir, 2, 0, 1, vec![8, 1, 1, 8, 1, 1])
                .expect("work_dir creatable"),
        ),
    );
    let combined = orchestrator.evaluate(&topo, &conf, &index).expect(
        "orchestrator should combine all three providers without an index-contract violation",
    );

    // Transparency: the orchestrator's sum must equal the three direct calls added by hand.
    let expected_energy = classical_direct.energy + embedding_direct.energy + qm_direct.energy;
    assert!(
        (combined.energy - expected_energy).abs() < 1e-6 * expected_energy.abs().max(1.0),
        "combined {} vs classical {} + embedding {} + xtb {} = {}",
        combined.energy,
        classical_direct.energy,
        embedding_direct.energy,
        qm_direct.energy,
        expected_energy
    );

    // The concrete, checkable version of "QM and MM influence each other": forces present on
    // both sides, not just the QM zone.
    let qm_atoms_with_force = combined
        .forces
        .iter()
        .filter(|&&(i, _)| i < QM_ZONE_SIZE)
        .count();
    let mm_atoms_with_force = combined
        .forces
        .iter()
        .filter(|&&(i, _)| i >= QM_ZONE_SIZE)
        .count();
    assert!(
        qm_atoms_with_force > 0,
        "QM zone atoms should feel forces (xtb's own + embedding's)"
    );
    assert!(
        mm_atoms_with_force > 0,
        "MM atoms should feel a real reaction force from the QM zone's charges via \
         ElectrostaticEmbedding — this is the actual QM->MM coupling this test exists to show"
    );

    // Sanity bounds, not an external oracle (none exists for this xtb-based combination —
    // P8, same tier as every other real-QM-engine test here): the QM zone's own energy
    // for two waters' worth of electrons should be a real, large negative number, and the
    // embedding term should be nonzero given real nonzero charges on both sides.
    assert!(
        qm_direct.energy < -1000.0,
        "QM zone energy {} kJ/mol doesn't look like a real xtb single-water electronic energy",
        qm_direct.energy
    );
    assert!(
        embedding_direct.energy.abs() > 1e-6,
        "embedding energy is suspiciously exactly zero"
    );

    eprintln!(
        "Combined QM/MM (real xtb + real Mulliken-charge embedding + classical MM):\n\
         \x20 E_QM        = {:>14.6} kJ/mol\n\
         \x20 E_classical = {:>14.6} kJ/mol\n\
         \x20 E_embedding = {:>14.6} kJ/mol\n\
         \x20 E_total     = {:>14.6} kJ/mol\n\
         \x20 {} QM-zone atoms with force, {} MM atoms with force (real bidirectional coupling)",
        qm_direct.energy,
        classical_direct.energy,
        embedding_direct.energy,
        combined.energy,
        qm_atoms_with_force,
        mm_atoms_with_force,
    );
}

/// Path (a) — real electrostatic embedding, no static-charge approximation and no separate
/// embedding provider. One `XtbInteraction` call carries the QM zone's own energy *and* the real
/// QM<->MM coupling (its SCF is polarized by the real MM charges via `pcharge`, and MM forces
/// come back via `pcgrad`). Zone-partitioned `LjCrfInteraction` supplies the classical part;
/// `ZonePartition::owner()` unconditionally routes inner-outer electrostatics to `Embedding`
/// regardless of which path handles it, so the classical CRF pass already skips those pairs —
/// no `ElectrostaticEmbedding` needed or registered here.
#[test]
fn path_a_real_electrostatic_embedding_on_the_real_water_dimer() {
    if !xtb_available() {
        eprintln!("skipping: xtb not found on PATH");
        return;
    }

    let root = fixture_root();
    let mut topo = build_topology(
        read_topology_file(root.join("water_dimer.top")).expect("read water_dimer.top"),
    );
    let coords = read_coordinates(root.join("water_dimer.cnf")).expect("read water_dimer.cnf");
    let imd = read_imd_file(root.join("water_dimer.imd")).expect("read water_dimer.imd");

    let n_solute = topo.num_solute_atoms();
    let atoms_per_solvent = topo.solvent_atom_template.len();
    topo.solvate((coords.positions.len() - n_solute) / atoms_per_solvent);

    let n_atoms = topo.num_atoms();
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions.clone();
    let periodicity = Periodicity::Vacuum(Vacuum);
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

    let inner: Vec<usize> = (0..QM_ZONE_SIZE).collect();
    let inner_sel = AtomSelection::from_indices(inner.clone(), n_atoms).unwrap();
    let partition = ZonePartition::new(n_atoms, &inner, &[]);

    let mut classical =
        LjCrfInteraction::new(imd.rcrf, 1.0, imd.epsrf, imd.appak).with_zone_partition(partition);
    let classical_direct = classical
        .contribute(&AtomSelection::all(n_atoms), &topo, &conf, &index)
        .unwrap();

    let work_dir = std::env::temp_dir().join("newen_xtb_qmmm_water_dimer_path_a");
    let mut qm_embedded = XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1, 8, 1, 1])
        .expect("work_dir creatable")
        .with_embedding(Embedding::Electrostatic)
        .with_cutoff(imd.rcrf);
    let qm_direct = qm_embedded
        .contribute(&inner_sel, &topo, &conf, &index)
        .expect("xtb electrostatic-embedding calculation should succeed");

    let total = classical_direct.energy + qm_direct.energy;

    let mm_atoms_with_force = qm_direct
        .forces
        .iter()
        .filter(|&&(i, _)| i >= QM_ZONE_SIZE)
        .count();
    assert!(
        mm_atoms_with_force > 0,
        "MM atoms should carry a real force straight from xtb's own pcgrad — no separate \
         embedding provider involved in this path"
    );
    assert!(
        qm_direct.energy < -1000.0,
        "QM zone energy {} kJ/mol doesn't look like a real xtb electronic energy",
        qm_direct.energy
    );

    eprintln!(
        "Path (a) combined QM/MM (real xtb electrostatic embedding + classical MM):\n\
         \x20 E_QM+embedding = {:>14.6} kJ/mol (one xtb call, real SCF polarization)\n\
         \x20 E_classical    = {:>14.6} kJ/mol\n\
         \x20 E_total        = {:>14.6} kJ/mol\n\
         \x20 {} MM atoms with force straight from xtb's pcgrad",
        qm_direct.energy, classical_direct.energy, total, mm_atoms_with_force,
    );
}
