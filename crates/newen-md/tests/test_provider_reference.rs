//! `PotentialProvider` against a real gromosXX reference — not a Rust-internal cross-check.
//!
//! `test_gromosXX_references.rs` validates the `md` binary end to end by running it as a
//! subprocess. This test validates `LjCrfInteraction` (`newen-forces/src/nonbonded/
//! interaction.rs`) directly and in-process, against the *same* real gromosXX-produced
//! `expected/energies.tre` oracle data, reusing the `pair_lj` reference system already in
//! `tests/gromosXX_references/`.
//!
//! `pair_lj` (two argon atoms, vacuum, zero charge, no exclusions) is the ideal fixture for
//! this provider's current scope: `LjCrfInteraction` only handles solute-solute atom-level
//! pairs (see its module docs) and `ConfigurationSpatialIndex` doesn't filter GROMOS
//! exclusions — `pair_lj` has none to get wrong, and zero charge makes the CRF term trivial,
//! so the comparison is a clean, real oracle check on the LJ math specifically.
//!
//! Run: `cargo test -p newen-md --test test_provider_reference`

use std::fs;
use std::path::PathBuf;

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Rectangular, Vacuum};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_forces::nonbonded::LjCrfInteraction;
use newen_forces::provider::PotentialProvider;
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};

const ENERGY_REL_TOL: f64 = 1e-8;

fn ref_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

/// Extract the potential-energy total (`ENERGY03` field index 2) from the first frame of a
/// real gromosXX `.tre` file. Deliberately independent of `test_gromosXX_references.rs`'s
/// `parse_energy03` (small enough to duplicate; this file should stay self-contained and not
/// take on a cross-test-file dependency for one field).
fn first_frame_potential(path: &std::path::Path) -> f64 {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut in_energy = false;
    let mut vals: Vec<f64> = Vec::new();
    for line in content.lines() {
        let t = line.trim();
        match t {
            "ENERGY03" => {
                in_energy = true;
                vals.clear();
            },
            "END" if in_energy => break,
            _ if in_energy && !t.starts_with('#') && !t.is_empty() => {
                if let Ok(v) = t.parse::<f64>() {
                    vals.push(v);
                }
            },
            _ => {},
        }
    }
    assert!(
        vals.len() >= 3,
        "expected at least 3 ENERGY03 values (total, kinetic, potential) in {}",
        path.display()
    );
    vals[2] // e_potential
}

#[test]
fn lj_crf_interaction_matches_real_gromosxx_energy() {
    let sys_dir = ref_root().join("pair_lj");

    let topo_data = read_topology_file(sys_dir.join("pair_lj.topo")).expect("read topology");
    let topo = build_topology(topo_data);

    let coords = read_coordinates(sys_dir.join("pair_lj.conf")).expect("read coordinates");

    let n_atoms = topo.num_atoms();
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions.clone();

    let periodicity = if coords.box_type == 0 {
        Periodicity::Vacuum(Vacuum)
    } else {
        Periodicity::Rectangular(Rectangular::new(coords.box_dims))
    };
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
    let region = AtomSelection::all(n_atoms);

    // RCUTP from pair_lj.in's PAIRLIST block; RF params from NONBONDED (irrelevant here —
    // both atoms carry zero charge, so CRF is trivially zero regardless of RF settings).
    let mut interaction = LjCrfInteraction::new(0.8, 1.0, 1.0, 0.0);
    let contribution = interaction
        .contribute(&region, &topo, &conf, &index)
        .expect("LjCrfInteraction::contribute should succeed on pair_lj");

    let expected_potential = first_frame_potential(&sys_dir.join("expected/energies.tre"));

    // pair_lj has zero charges and no other terms (no bonds, no exclusions, no restraints),
    // so gromosXX's reported potential energy is exactly this provider's LJ+CRF contribution.
    let rel_diff = (contribution.energy - expected_potential).abs() / expected_potential.abs();
    assert!(
        rel_diff < ENERGY_REL_TOL,
        "LjCrfInteraction energy {} vs real gromosXX energy {} (rel diff {})",
        contribution.energy,
        expected_potential,
        rel_diff
    );
}
