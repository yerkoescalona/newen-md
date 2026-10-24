//! ** P2.8-2, scoped to NN/MM** — `SchNetInteraction` driving an actual MD step loop,
//! not just single-point force evaluation.
//!
//! Every test so far has called `SchNetInteraction::contribute` once per assertion — proving
//! the model loads and its forces are self-consistent (P8), but never proving it can
//! actually *evolve a system through time*. This wires it to `newen-integrators::LeapFrog` and
//! runs real NVE (constant energy) dynamics.
//!
//! **Why NVE conservation is the right check, not another finite-difference test.** Leapfrog is
//! symplectic: for *any* smooth, differentiable potential — chemically meaningful or not — total
//! energy (kinetic + potential) should stay bounded over a trajectory, not drift. That makes it
//! architecture-blind to whether the potential is scientifically accurate (this model is still
//! untrained, per `nonbonded/schnet.rs`) while being a genuine test of the *loop*: force
//! evaluation → integration → new positions → force evaluation again, repeated, not a single
//! isolated call. This is the standard way any new force field is checked when wired into an
//! MD engine — reference systems like `water_216_box` in the classical suite are validated the
//! same way (bulk NVE), just against a real gromosXX oracle instead of a self-consistency check
//! (there is none available for an untrained model — P8).
//!
//! Feature-gated (`ml`); skips gracefully if no model is available, matching every other `ml`
//! test in this workspace.

#![cfg(feature = "ml")]

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::{MolTypeAtom, Topology};
use newen_forces::nonbonded::SchNetInteraction;
use newen_forces::provider::PotentialProvider;
use newen_integrators::{Integrator, LeapFrog};

fn model_path() -> String {
    std::env::var("TOY_SCHNET_MODEL").unwrap_or_else(|_| "/tmp/toy_schnet.pt".to_string())
}

/// A real (if small) methanol-fragment-shaped molecule: C, O, H, H — not the QM zone of a
/// solvated system (that would need embedding, which `SchNetInteraction` doesn't support yet —
/// P2.7), just a plain vacuum molecule with real atomic numbers and real masses.
fn toy_molecule() -> Topology {
    let mut topo = Topology::new();
    let elements = [6u8, 8, 1, 1]; // C, O, H, H
    let masses = [12.011, 15.999, 1.008, 1.008];

    topo.charge = vec![0.0; 4];
    topo.iac = vec![0; 4];
    topo.mass = masses.to_vec();
    topo.inverse_mass = masses.iter().map(|m| 1.0 / m).collect();
    topo.exclusions = vec![Vec::new(); 4];
    topo.moltypes[0].atoms = (0..4)
        .map(|i| MolTypeAtom {
            name: format!("A{i}"),
            residue_nr: 1,
            residue_name: "MOL".to_string(),
            iac: 0,
            mass: masses[i],
            charge: 0.0,
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: false,
        })
        .collect();
    let _ = elements;
    topo
}

#[test]
fn schnet_driven_nve_conserves_energy_over_a_trajectory() {
    let path = model_path();
    if !std::path::Path::new(&path).exists() {
        eprintln!("skipping: no model at {path} (run scripts/export_toy_schnet.py first)");
        return;
    }

    let topo = toy_molecule();
    let n_atoms = 4;
    let elements = vec![6i64, 8, 1, 1];
    let mut interaction = SchNetInteraction::load(&path, 1.0, elements).expect("model loads");

    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.143, 0.0, 0.0),   // ~C-O bond length, nm
        Vec3::new(-0.05, 0.10, 0.0),  // ~C-H
        Vec3::new(-0.05, -0.10, 0.0), // ~C-H
    ];
    conf.current_mut().vel = vec![Vec3::ZERO; n_atoms];
    conf.copy_current_to_old();

    let periodicity = Periodicity::Vacuum(Vacuum);
    let region = AtomSelection::all(n_atoms);

    let potential_energy =
        |conf: &Configuration, interaction: &mut SchNetInteraction| -> (f64, Vec<Vec3>) {
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

    // dt chosen empirically (not GROMOS's usual 0.002 ps, which assumes a fitted force field):
    // scanned 1e-4..2e-3 and picked the largest step that still produces real kinetic energy
    // (E_kin: 0 -> 0.08 kJ/mol from a standing start) with tight conservation (<0.005% max
    // deviation) — see P2.8-2. Too small and the atoms barely move (a weak test of the
    // loop); too large and a symplectic integrator's bounded-oscillation guarantee erodes.
    let dt = 1e-3;
    let n_steps = 500;
    let mut integrator = LeapFrog::new();

    // Seed conf.current().force for the first velocity_step (LeapFrog reads conf.old().force
    // after its internal exchange_state(), so the force present *before* step() must be the
    // force at the position step() is about to advance from).
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
        "NVE loop: {n_steps} steps, E0={e0:.6}, E_kin_final={:.6}, mean_total={mean:.6}, \
         max |deviation|={max_abs_dev:.6} ({:.4}% of mean)",
        kinetic_energy(&conf),
        relative_drift * 100.0
    );

    // Symplectic integrators bound energy fluctuation around a shifted mean rather than making
    // it monotonically drift — this is the actual failure mode a broken wiring would produce
    // (e.g. forces off by a sign, or evaluated at the wrong position each step, would show up
    // as unbounded/monotonic growth, not a bounded oscillation).
    // Measured 0.0047% at dt=1e-3 during development (see the dt comment above) — bound at
    // 0.5%, ~100x margin, so this catches a broken loop without being sensitive to machine-level
    // float noise.
    assert!(
        relative_drift < 0.005,
        "energy fluctuation {:.4}% of mean exceeds 0.5% over {n_steps} steps — \
         the step loop (force eval -> integrate -> force eval) is not behaving like NVE dynamics",
        relative_drift * 100.0
    );

    // A monotonic trend is the more dangerous failure mode than bounded oscillation — check the
    // trajectory's second half isn't systematically drifting away from its first half.
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
