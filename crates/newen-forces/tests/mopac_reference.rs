//! `MopacInteraction` on a real multi-element molecule (methanol: C, O, H), not just the
//! water-only fixtures in `nonbonded/mopac.rs`'s own unit tests.
//!
//! **Geometry provenance.** The BuRNN tutorial (`.local/gromos_tutorial_livecoms/`, gitignored)
//! ships 1722 real archived MOPAC `.aux` files with real methanol+water cluster geometries — but
//! that repository has no explicit license file for its data (unlike the Poliak/Zenodo QM/MM
//! dataset, which is CC BY 4.0 and already used by `gromosXX_qmmm_references/water_dimer_mechst/`).
//! Rather than check in geometry lifted from an ambiguously-licensed source, this fixture is a
//! **standard, self-constructed** methanol geometry (typical bond lengths: C-O 1.43 Å, C-H
//! 1.09 Å, O-H 0.96 Å, tetrahedral angles around C) — real chemistry, just not extracted from
//! any external dataset.
//!
//! **Oracle**: pinned from a real run of this machine's own installed MOPAC 23.1.2 (same pattern
//! `nonbonded/mopac.rs`'s water test and `nonbonded/xtb.rs`'s pinned tests use — see those
//! modules' docs for why cross-version MOPAC comparison isn't reliable at a fixed geometry).
//!
//! Skips gracefully if `mopac` isn't on `PATH`.

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::Topology;
use newen_core::units::KCAL_TO_KJ;
use newen_forces::nonbonded::{MopacInteraction, MopacMethod};
use newen_forces::provider::PotentialProvider;

fn mopac_available() -> bool {
    std::process::Command::new("mopac")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Standard-geometry methanol (CH3OH) — see module docs for why this is self-constructed rather
/// than lifted from the BuRNN tutorial archive.
fn methanol_positions_nm() -> Vec<Vec3> {
    vec![
        Vec3::new(0.0000000000, 0.0000000000, 0.0000000000), // C
        Vec3::new(0.1430000000, 0.0000000000, 0.0000000000), // O
        Vec3::new(-0.0363849477, 0.1027479225, 0.0000000000), // H
        Vec3::new(-0.0363849477, -0.0513739613, 0.0889823111), // H
        Vec3::new(-0.0363849477, -0.0513739613, -0.0889823111), // H
        Vec3::new(0.1726656315, 0.0874509053, 0.0262352716), // H
    ]
}

#[test]
fn methanol_energy_and_charges_match_pinned_mopac_oracle() {
    if !mopac_available() {
        eprintln!("skipping: mopac not found on PATH");
        return;
    }

    let work_dir = std::env::temp_dir().join("newen_mopac_methanol_reference");
    let mut interaction =
        MopacInteraction::new(work_dir, MopacMethod::Pm7, 0, 1, vec![6, 8, 1, 1, 1, 1])
            .expect("work_dir creatable");

    let topo = Topology::new();
    let mut conf = Configuration::new(6, 1, 1);
    conf.current_mut().pos = methanol_positions_nm();
    let region = AtomSelection::all(6);
    let periodicity = Periodicity::Vacuum(Vacuum);
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

    let c = interaction
        .contribute(&region, &topo, &conf, &index)
        .expect("mopac calculation should succeed");

    // Pinned against a real run of this machine's own installed MOPAC 23.1.2.
    let expected_energy = -47.4973590171721 * KCAL_TO_KJ;
    assert!(
        (c.energy - expected_energy).abs() < 1e-4,
        "energy {} vs pinned oracle {}",
        c.energy,
        expected_energy
    );
    assert_eq!(c.forces.len(), 6, "all 6 atoms should carry a force");

    // `MopacInteraction` doesn't expose charges through `Contribution` (same deferred-channel
    // reasoning as `XtbInteraction` — see `provider.rs`), so read them the same way
    // `xtb_qmmm_water_dimer_reference.rs` reads xtb's `charges` side-effect file: directly from
    // the `.aux` this same call just wrote.
    let aux_content = std::fs::read_to_string(
        std::env::temp_dir()
            .join("newen_mopac_methanol_reference")
            .join("region.aux"),
    )
    .expect("mopac should have written region.aux as a side effect");
    let charges: Vec<f64> = aux_content
        .lines()
        .skip_while(|l| !l.trim().starts_with("ATOM_CHARGES"))
        .skip(1)
        .take_while(|l| !l.trim().is_empty() && !l.contains('['))
        .flat_map(|l| l.split_whitespace())
        .map(|tok| tok.parse::<f64>().expect("charges: one float per token"))
        .collect();
    assert_eq!(charges.len(), 6, "expected one Mulliken charge per atom");
    let total_charge: f64 = charges.iter().sum();
    assert!(
        total_charge.abs() < 1e-3,
        "methanol is neutral, but Mulliken charges summed to {total_charge}: {charges:?}"
    );
    eprintln!(
        "methanol PM7: E={:.6} kJ/mol, charges={charges:?} (sum={total_charge:.2e}), \
         {} atoms with force",
        c.energy,
        c.forces.len()
    );
}
