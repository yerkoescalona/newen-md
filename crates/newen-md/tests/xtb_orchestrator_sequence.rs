//! `ProviderOrchestrator` driven by the real `AlgorithmSequence` step loop — the P2.8-6 exit
//! criterion.
//!
//! Every existing `ProviderOrchestrator`/provider test (`test_orchestrator_reference.rs`,
//! `xtb_nve_loop.rs`, `schnet_nve_loop.rs`) either calls `orchestrator.evaluate()` directly or
//! drives a bespoke `LeapFrog` loop — none of them go through the actual `AlgorithmSequence` +
//! `Forcefield` machinery both the `md` CLI binary and `pyo3-newen` build and run via
//! `run_step()`. This does: a real sequence (`Forcefield` -> `ProviderOrchestratorAlgorithm` ->
//! `LeapFrogVelocity` -> `LeapFrogPosition` -> `TemperatureCalculation` -> `EnergyCalculation`),
//! proving the orchestrator's energy/force/virial actually flow through that pipeline, not a
//! shortcut around it.
//!
//! The classical `Forcefield` here sees a zero-charge, zero-LJ, bond-free water monomer, so its
//! own contribution is legitimately zero — `XtbInteraction` (via the orchestrator) carries all
//! the potential energy. Deliberately narrower than "both nonzero simultaneously" (same
//! incremental-proof discipline `schnet_nve_loop.rs` used); what's new here is the *pipeline*,
//! not a second nonzero force source. Validated by NVE energy conservation (no gromosXX oracle
//! exists for this combination — P8), reading `conf.old.energies.total` after each
//! `run_step()` — the GROMOS convention this codebase's own `EnergyCalculation`/
//! `TemperatureCalculation` already use (results land in `old()` after the leapfrog position
//! step's internal state exchange).
//!
//! Skips gracefully if `xtb` isn't on `PATH`, same as P2.8-2b's tests.

use newen_core::algorithm::{AlgorithmSequence, SimulationState};
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::pairlist::{PairlistAlgorithm, PairlistContainer, StandardPairlistAlgorithm};
use newen_core::selection::AtomSelection;
use newen_core::topology::{LJParameters, MolTypeAtom, Topology};
use newen_forces::nonbonded::{CRFParameters, XtbInteraction};
use newen_forces::orchestrator::ProviderOrchestrator;
use newen_forces::orchestrator_algorithm::ProviderOrchestratorAlgorithm;
use newen_integrators::algorithms::{
    EnergyCalculation, Forcefield, LeapFrogPosition, LeapFrogVelocity, TemperatureCalculation,
};

fn xtb_available() -> bool {
    std::process::Command::new("xtb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// A water monomer with zero charge/LJ/bonds, so the classical `Forcefield` in the sequence
/// contributes exactly zero — `XtbInteraction` (via the orchestrator) carries all the potential
/// energy, isolating what this test actually needs to prove (see module docs).
fn water_topology() -> Topology {
    let masses = [15.9994, 1.008, 1.008];
    let mut topo = Topology::new();
    topo.charge = vec![0.0; 3];
    topo.iac = vec![0; 3];
    topo.mass = masses.to_vec();
    topo.inverse_mass = masses.iter().map(|m| 1.0 / m).collect();
    topo.exclusions = vec![Vec::new(); 3];
    topo.lj_parameters = vec![vec![LJParameters {
        c6: 0.0,
        c12: 0.0,
        cs6: 0.0,
        cs12: 0.0,
    }]];
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
fn orchestrator_driven_nve_conserves_energy_through_the_real_algorithm_sequence() {
    if !xtb_available() {
        eprintln!("skipping: xtb not found on PATH");
        return;
    }

    let topo = water_topology();
    let n_atoms = 3;

    let crf_params = CRFParameters::new(1.4, 1.0, 1.0, 0.0);
    let lj_params = Forcefield::convert_lj_parameters(&topo);
    let pairlist = PairlistContainer::new(1.4, 1.4, 0.0);
    let pairlist_algorithm = PairlistAlgorithm::Standard(StandardPairlistAlgorithm::new(false));
    let forcefield = Forcefield::new(
        lj_params,
        crf_params,
        Periodicity::Vacuum(Vacuum),
        pairlist,
        pairlist_algorithm,
    );

    let work_dir = std::env::temp_dir().join("newen_xtb_orchestrator_sequence");
    let mut orchestrator = ProviderOrchestrator::new();
    orchestrator.register(
        AtomSelection::all(n_atoms),
        Box::new(
            XtbInteraction::new(work_dir, 2, 0, 1, vec![8, 1, 1]).expect("work_dir creatable"),
        ),
    );
    let orchestrator_algorithm =
        ProviderOrchestratorAlgorithm::new(orchestrator, Periodicity::Vacuum(Vacuum));

    let mut sequence = AlgorithmSequence::new();
    sequence.push(Box::new(forcefield));
    sequence.push(Box::new(orchestrator_algorithm));
    sequence.push(Box::new(LeapFrogVelocity::new()));
    sequence.push(Box::new(LeapFrogPosition::new()));
    sequence.push(Box::new(TemperatureCalculation::new()));
    sequence.push(Box::new(EnergyCalculation::new()));

    let mut conf = newen_core::configuration::Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0758602, 0.0, 0.0504284),
        Vec3::new(0.0758602, 0.0, -0.0504284),
    ];
    conf.current_mut().vel = vec![Vec3::ZERO; n_atoms];
    conf.copy_current_to_old();

    // Same dt as `xtb_nve_loop.rs`, tuned for the same real GFN2-xTB water PES; only the driving
    // mechanism differs here (a real `AlgorithmSequence` instead of a bespoke `LeapFrog` loop).
    let dt = 2e-4;
    let n_steps = 120;
    let mut sim_state = SimulationState::new(dt, n_steps);

    sequence
        .init(&topo, &mut conf, &sim_state)
        .expect("sequence init should succeed");

    let mut total_energies = Vec::with_capacity(n_steps);
    for _ in 0..n_steps {
        sequence
            .run_step(&topo, &mut conf, &sim_state)
            .expect("run_step should succeed");
        sim_state.advance();
        total_energies.push(conf.old().energies.total());
    }

    let mean: f64 = total_energies.iter().sum::<f64>() / total_energies.len() as f64;
    let max_abs_dev = total_energies
        .iter()
        .map(|&e| (e - mean).abs())
        .fold(0.0, f64::max);
    let mean_abs = mean.abs().max(1e-6);
    let relative_drift = max_abs_dev / mean_abs;

    eprintln!(
        "orchestrator-via-AlgorithmSequence NVE: {n_steps} steps, mean_total={mean:.6}, \
         max |deviation|={max_abs_dev:.6} ({:.4}% of mean)",
        relative_drift * 100.0
    );

    assert!(
        relative_drift < 0.005,
        "energy fluctuation {:.4}% of mean exceeds 0.5% over {n_steps} steps — the orchestrator's \
         contribution is not flowing correctly through the real AlgorithmSequence step loop",
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
