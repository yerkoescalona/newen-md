//! `XtbInteraction` driving an actual MD step loop, not just single-point force evaluation.
//!
//! Mirrors `schnet_nve_loop.rs`'s rationale exactly, but with a real semiempirical QM engine
//! instead of an untrained neural net: leapfrog is symplectic, so total energy should stay
//! bounded over a trajectory for any smooth potential, which makes NVE conservation a genuine
//! test of the loop (force eval -> integrate -> force eval, repeated) rather than of chemical
//! accuracy. Skips gracefully if `xtb` isn't on `PATH`.
//!
//! Each step spawns a real `xtb` subprocess, so this uses far fewer steps than the SchNet
//! version's 500 to keep wall time reasonable.

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::{MolTypeAtom, Topology};
use newen_forces::nonbonded::XtbInteraction;
use newen_forces::provider::PotentialProvider;
use newen_integrators::{Integrator, LeapFrog};

fn xtb_available() -> bool {
    std::process::Command::new("xtb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// A real water monomer: O, H, H with real atomic masses — a plain vacuum molecule, no
/// embedding (`XtbInteraction` doesn't support it — see `nonbonded/xtb.rs`).
fn water_topology() -> Topology {
    let mut topo = Topology::new();
    let masses = [15.9994, 1.008, 1.008];

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
fn xtb_driven_nve_conserves_energy_over_a_trajectory() {
    if !xtb_available() {
        eprintln!("skipping: xtb not found on PATH");
        return;
    }

    let topo = water_topology();
    let n_atoms = 3;
    let work_dir = std::env::temp_dir().join("newen_xtb_nve_loop");
    let mut interaction =
        XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1]).expect("work_dir creatable");

    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0758602, 0.0, 0.0504284),
        Vec3::new(0.0758602, 0.0, -0.0504284),
    ];
    conf.current_mut().vel = vec![Vec3::ZERO; n_atoms];
    conf.copy_current_to_old();

    let periodicity = Periodicity::Vacuum(Vacuum);
    let region = AtomSelection::all(n_atoms);

    let potential_energy =
        |conf: &Configuration, interaction: &mut XtbInteraction| -> (f64, Vec<Vec3>) {
            let index = ConfigurationSpatialIndex::new(conf, &periodicity);
            let c = interaction
                .contribute(&region, &topo, conf, &index)
                .unwrap();
            let mut forces = vec![Vec3::ZERO; n_atoms];
            for (i, f) in c.forces {
                forces[i] = f;
            }
            (c.energy, forces)
        };

    let kinetic_energy = |conf: &Configuration| -> f64 {
        (0..n_atoms)
            .map(|i| 0.5 * topo.mass[i] * conf.current().vel[i].length_squared())
            .sum()
    };

    // dt chosen empirically for a real GFN2-xTB water PES (much stiffer than the untrained
    // SchNet toy model `schnet_nve_loop.rs` tunes for): scanned down from 1e-3 ps until the O-H
    // stretch stopped blowing up leapfrog's bounded-oscillation guarantee.
    let dt = 2e-4;
    let n_steps = 120;
    let mut integrator = LeapFrog::new();

    let (e0, f0) = potential_energy(&conf, &mut interaction);
    conf.current_mut().force = f0;
    conf.copy_current_to_old();

    let mut total_energies = Vec::with_capacity(n_steps + 1);
    total_energies.push(e0 + kinetic_energy(&conf));

    for _ in 0..n_steps {
        integrator.step(dt, &topo, &mut conf);
        let (pe, f) = potential_energy(&conf, &mut interaction);
        conf.current_mut().force = f;
        total_energies.push(pe + kinetic_energy(&conf));
    }

    let mean: f64 = total_energies.iter().sum::<f64>() / total_energies.len() as f64;
    let max_abs_dev = total_energies
        .iter()
        .map(|&e| (e - mean).abs())
        .fold(0.0, f64::max);
    let mean_abs = mean.abs().max(1e-6);
    let relative_drift = max_abs_dev / mean_abs;

    eprintln!(
        "xtb NVE loop: {n_steps} steps, E0={e0:.6}, E_kin_final={:.6}, mean_total={mean:.6}, \
         max |deviation|={max_abs_dev:.6} ({:.4}% of mean)",
        kinetic_energy(&conf),
        relative_drift * 100.0
    );

    assert!(
        relative_drift < 0.005,
        "energy fluctuation {:.4}% of mean exceeds 0.5% over {n_steps} steps — \
         the step loop (force eval -> integrate -> force eval) is not behaving like NVE dynamics",
        relative_drift * 100.0
    );

    let half = total_energies.len() / 2;
    let first_half_mean: f64 = total_energies[..half].iter().sum::<f64>() / half as f64;
    let second_half_mean: f64 =
        total_energies[half..].iter().sum::<f64>() / (total_energies.len() - half) as f64;
    let half_drift = (second_half_mean - first_half_mean).abs() / mean_abs;
    assert!(
        half_drift < 0.002,
        "first-half vs second-half mean energy differs by {:.4}% of mean — looks like drift, not oscillation",
        half_drift * 100.0
    );
}
