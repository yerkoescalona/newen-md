//! `SchNetInteraction` on a real molecule from GROMOS's own NN(QM)/MM tutorial.
//!
//! Not a GROMOS oracle test — there is no gromosXX-produced reference to match against for
//! an ML potential (P8), and this pass's model is architecturally real (SchNetPack
//! 2) but untrained (see `nonbonded/schnet.rs` module docs). What this *does* check, that the
//! synthetic 3-atom unit test in `schnet.rs` doesn't: the seam holds up on a real molecule, a
//! real periodic box, and the real QM-zone/cutoff convention from GROMOS's own tutorial —
//! `#![cfg(feature = "ml")]`.
//!
//! Fixture: `.local/gromos_tutorial_livecoms/tutorial_files/t_06` (Tutorial 6, "BuRNN"),
//! cloned into `.local/` — gitignored, not checked into this repo, so this test skips
//! gracefully if it isn't present locally (matches how the `ml`-feature model path is
//! already handled). The real `md_burnn/meoh.qmm` `QMZONE` block defines the QM region as
//! atoms 1-6 of the system — CA(C), OB(O), HA1/HA2/HA3(H), HB1(H), `QMZ` = `[6,8,1,1,1,1]` —
//! and `eq/eq_meoh_5.cnf`'s `POSITION` block confirms those are literally the file's first 6
//! atoms (methanol, followed by SPC solvent). The cutoff (1.4 nm) is `md_burnn/md.imd`'s
//! `QMMM` block's `RCUTQM`. None of this is invented for the test — it's read directly off
//! the tutorial's own files.

#![cfg(feature = "ml")]

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Rectangular, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::Topology;
use newen_forces::nonbonded::SchNetInteraction;
use newen_forces::provider::PotentialProvider;
use newen_io::coordinate::read_coordinates;
use std::path::PathBuf;

fn tutorial_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../.local/gromos_tutorial_livecoms/tutorial_files/t_06")
}

fn model_path() -> String {
    std::env::var("TOY_SCHNET_MODEL").unwrap_or_else(|_| "/tmp/toy_schnet.pt".to_string())
}

#[test]
fn schnet_provider_on_real_burnn_methanol_qm_zone() {
    let cnf_path = tutorial_root().join("eq/eq_meoh_5.cnf");
    let model = model_path();

    if !cnf_path.exists() {
        eprintln!(
            "skipping: real BuRNN tutorial not cloned at {} \
             (git clone github.com/biomos/gromos_tutorial_livecoms into .local/)",
            cnf_path.display()
        );
        return;
    }
    if !std::path::Path::new(&model).exists() {
        eprintln!("skipping: no model at {model} (run scripts/export_toy_schnet.py first)");
        return;
    }

    let coords = read_coordinates(&cnf_path).expect("read real BuRNN eq_meoh_5.cnf");
    let n_atoms = coords.positions.len();
    assert!(
        n_atoms > 6,
        "expected a solvated system (methanol + SPC water), got {n_atoms} atoms"
    );

    // The real QMZONE from meoh.qmm: atoms 1-6 (CA, OB, HA1, HA2, HA3, HB1), Z = [6,8,1,1,1,1].
    let qm_zone: Vec<usize> = (0..6).collect();
    let atomic_numbers = vec![6i64, 8, 1, 1, 1, 1];

    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions.clone();
    let periodicity = Periodicity::Rectangular(Rectangular::new(coords.box_dims));
    let region = AtomSelection::from_indices(qm_zone, n_atoms).unwrap();
    let topo = Topology::new(); // SchNetInteraction::contribute doesn't read topology

    // Real RCUTQM from md_burnn/md.imd's QMMM block.
    let mut interaction =
        SchNetInteraction::load(&model, 1.4, atomic_numbers).expect("model should load");

    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
    let contribution = interaction
        .contribute(&region, &topo, &conf, &index)
        .expect("SchNetInteraction should run on the real BuRNN QM zone");

    assert_eq!(contribution.forces.len(), 6);
    assert!(contribution.energy.is_finite());
    for (_, f) in &contribution.forces {
        assert!(f.x.is_finite() && f.y.is_finite() && f.z.is_finite());
    }

    // Same validation tier as the synthetic-fixture test in schnet.rs (P8: no
    // GROMOS oracle for an ML potential), just on the real geometry instead of a toy one.
    let mut energy_at = |positions: &[Vec3]| -> f64 {
        let mut conf = Configuration::new(n_atoms, 1, 1);
        conf.current_mut().pos = positions.to_vec();
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        interaction
            .contribute(&region, &topo, &conf, &index)
            .unwrap()
            .energy
    };

    let h = 1e-4;
    let mut max_abs_diff = 0.0f64;
    for atom in 0..6 {
        for axis in 0..3 {
            let mut plus = coords.positions.clone();
            let mut minus = coords.positions.clone();
            let delta = match axis {
                0 => Vec3::new(h, 0.0, 0.0),
                1 => Vec3::new(0.0, h, 0.0),
                _ => Vec3::new(0.0, 0.0, h),
            };
            plus[atom] += delta;
            minus[atom] -= delta;

            let finite_diff_force = -(energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
            let (_, model_force) = contribution.forces[atom];
            let model_component = match axis {
                0 => model_force.x,
                1 => model_force.y,
                _ => model_force.z,
            };
            max_abs_diff = max_abs_diff.max((finite_diff_force - model_component).abs());
        }
    }

    assert!(
        max_abs_diff < 5e-3,
        "max |finite-diff - model force| = {max_abs_diff} on the real BuRNN geometry (tolerance 5e-3)"
    );
}
