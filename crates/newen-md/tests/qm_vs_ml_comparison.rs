//! QM/MM vs ML/MM — the comparison step of the generate/train/compare pipeline.
//!
//! Runs real `XtbInteraction` (`Embedding::None`) and a *trained* `SchNetInteraction` (loaded
//! from `/tmp/trained_water_schnet.pt`, produced by `scripts/train_qmmm_schnet.py` on data from
//! `cargo run -p newen-md --example generate_qm_training_data`) on the same held-out water-
//! monomer configurations, and compares energy/force agreement.
//!
//! **What this does and doesn't claim.** This proves the *pipeline* works end-to-end — real QM
//! data generated, a real model trained on it, the trained model loaded and run as an ML
//! potential, and its predictions actually correlating with the QM engine it was trained to
//! approximate. It is **not** a chemical-accuracy claim: ~1000 frames, a small (32-feature,
//! 2-interaction) model, and a few hundred epochs of CPU training on a single water molecule is a
//! toy-scale fit, not a production BuRNN-quality potential. Tolerances below are generous and
//! documented as such.
//!
//! **Held-out, not a reserved split.** The comparison trajectory starts from a geometry the
//! generator script doesn't produce (a deliberately different perturbation), rather than holding
//! back part of the original dataset — simpler, and unambiguously "not seen in training."
//!
//! Skips gracefully if `xtb` isn't on `PATH` or the trained model file doesn't exist, printing
//! the commands to run first.

#![cfg(feature = "ml")]

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::{MolTypeAtom, Topology};
use newen_forces::nonbonded::{SchNetInteraction, XtbInteraction};
use newen_forces::provider::PotentialProvider;
use newen_integrators::{Integrator, LeapFrog};

const MODEL_PATH: &str = "/tmp/trained_water_schnet.pt";

fn xtb_available() -> bool {
    std::process::Command::new("xtb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

fn water_topology() -> Topology {
    let masses = [15.9994, 1.008, 1.008];
    let mut topo = Topology::new();
    topo.charge = vec![0.0; 3];
    topo.iac = vec![0; 3];
    topo.mass = masses.to_vec();
    topo.inverse_mass = masses.iter().map(|m| 1.0 / m).collect();
    topo.exclusions = vec![Vec::new(); 3];
    topo.moltypes[0].atoms = (0..3)
        .map(|i| MolTypeAtom {
            name: format!("A{i}"),
            residue_nr: 1,
            residue_name: "WAT".to_string(),
            iac: 0,
            mass: masses[i],
            charge: 0.0,
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: false,
        })
        .collect();
    topo
}

#[test]
fn trained_schnet_tracks_real_xtb_on_held_out_water_configurations() {
    if !xtb_available() {
        eprintln!("skipping: xtb not found on PATH");
        return;
    }
    if !std::path::Path::new(MODEL_PATH).exists() {
        eprintln!(
            "skipping: no trained model at {MODEL_PATH} — run first:\n\
             \x20 cargo run -p newen-md --example generate_qm_training_data\n\
             \x20 (cd py-newen && uv sync --group ml) && source py-newen/.venv/bin/activate && \\\n\
             \x20 python3 scripts/train_qmmm_schnet.py"
        );
        return;
    }

    let topo = water_topology();
    let n_atoms = 3;
    let atomic_numbers = vec![8i64, 1, 1];
    let periodicity = Periodicity::Vacuum(Vacuum);
    let region = AtomSelection::all(n_atoms);

    let xtb_work_dir = std::env::temp_dir().join("newen_qm_vs_ml_comparison");
    let mut xtb = XtbInteraction::new(xtb_work_dir, 2, 0, 1, atomic_numbers.clone())
        .expect("work_dir creatable");
    let mut ml = SchNetInteraction::load(MODEL_PATH, 1.0, atomic_numbers)
        .expect("trained model should load");

    // A held-out starting geometry: a deliberately different, larger perturbation than the
    // generator's own (which draws from N(0, 0.01)) — a point the training trajectories are
    // unlikely to have sampled closely.
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = vec![
        Vec3::new(0.012, -0.008, 0.005),
        Vec3::new(0.0758602 + 0.015, 0.010, 0.0504284 - 0.006),
        Vec3::new(0.0758602 - 0.009, -0.004, -0.0504284 + 0.011),
    ];
    conf.current_mut().vel = vec![Vec3::ZERO; n_atoms];
    conf.copy_current_to_old();

    let xtb_energy_forces = |conf: &Configuration, xtb: &mut XtbInteraction| -> (f64, Vec<Vec3>) {
        let index = ConfigurationSpatialIndex::new(conf, &periodicity);
        let c = xtb
            .contribute(&region, &topo, conf, &index)
            .expect("xtb calculation should succeed");
        let mut forces = vec![Vec3::ZERO; n_atoms];
        for (i, f) in c.forces {
            forces[i] = f;
        }
        (c.energy, forces)
    };

    let dt = 2e-4;
    let n_steps = 40;
    let mut integrator = LeapFrog::new();

    let (e0, f0) = xtb_energy_forces(&conf, &mut xtb);
    conf.current_mut().force = f0;
    conf.copy_current_to_old();

    let mut energy_sq_err = 0.0;
    let mut force_sq_err = 0.0;
    let mut force_count = 0usize;
    let mut n_frames = 0usize;

    let mut compare =
        |conf: &Configuration, xtb: &mut XtbInteraction, ml: &mut SchNetInteraction| {
            let (e_qm, f_qm) = xtb_energy_forces(conf, xtb);
            let index = ConfigurationSpatialIndex::new(conf, &periodicity);
            let ml_contribution = ml
                .contribute(&region, &topo, conf, &index)
                .expect("trained model forward pass should succeed");
            let mut f_ml = vec![Vec3::ZERO; n_atoms];
            for (i, f) in ml_contribution.forces {
                f_ml[i] = f;
            }

            energy_sq_err += (e_qm - ml_contribution.energy).powi(2);
            for i in 0..n_atoms {
                force_sq_err += (f_qm[i] - f_ml[i]).length_squared();
                force_count += 3;
            }
            n_frames += 1;
            (e_qm, f_qm)
        };

    let (_, f0_qm) = compare(&conf, &mut xtb, &mut ml);
    conf.current_mut().force = f0_qm;

    for _ in 0..n_steps {
        integrator.step(dt, &topo, &mut conf);
        let (_, f_qm) = compare(&conf, &mut xtb, &mut ml);
        conf.current_mut().force = f_qm;
    }

    let energy_rmse = (energy_sq_err / n_frames as f64).sqrt();
    let force_rmse = (force_sq_err / force_count as f64).sqrt();

    eprintln!(
        "QM/MM vs ML/MM over {n_frames} held-out frames: energy RMSE = {energy_rmse:.3} kJ/mol, \
         force component RMSE = {force_rmse:.3} kJ/mol/nm (E0_qm={e0:.3} kJ/mol)"
    );

    // Generous, documented tolerances (see module docs) — this checks the trained model tracks
    // xtb's PES at all, not chemical accuracy. A completely untrained/random model would show
    // energy RMSE orders of magnitude larger than the ~tens-of-kJ/mol scale energy variation
    // actually present in this trajectory.
    assert!(
        energy_rmse < 200.0,
        "energy RMSE {energy_rmse:.3} kJ/mol too large — trained model isn't tracking xtb at all"
    );
    assert!(
        force_rmse < 5000.0,
        "force RMSE {force_rmse:.3} kJ/mol/nm too large — trained model isn't tracking xtb at all"
    );
}
