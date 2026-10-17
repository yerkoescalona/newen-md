//! Generates real QM training data for an ML potential — the "QM/MM for training" half of the
//! QM/MM → ML/MM pipeline.
//!
//! Drives a water monomer (O, H, H) through several independent short real `xtb`-driven NVE
//! trajectories (same proven-stable setup `crates/newen-md/tests/xtb_nve_loop.rs` already
//! validates: `dt=2e-4`, `XtbInteraction` with `Embedding::None` — xtb's own real intramolecular
//! chemistry keeps bond lengths physically bounded with no classical bond potential or SHAKE
//! needed), and records every step's energy/forces/positions to a training dataset file.
//!
//! **Training target is deliberately the *isolated* QM-zone energy** (`Embedding::None`), not an
//! electrostatically-embedded one. `SchNetInteraction`'s architecture
//! has no environment/charge input channel, so training on embedded energies would bake in one
//! specific environment's contribution into a model that structurally can't represent it (same
//! reasoning as the BuRNN training target in `zones.rs`'s docs).
//!
//! Multiple independent trajectories, not one long one — a single initial condition only explores
//! one energy shell, poor coverage for training.
//!
//! Output format (no `serde_json` dependency needed — this file is both written and read by code
//! we fully control): header line `n_atoms Z1 Z2 ... Zn`, then one line per frame:
//! `energy px1 py1 pz1 fx1 fy1 fz1 ...` — GROMOS-native units throughout (nm, kJ/mol, kJ/mol/nm),
//! matching what `nonbonded/schnet.rs`'s Rust side expects with no conversion.
//!
//! Usage: `cargo run -p newen-md --example generate_qm_training_data [output_path]`

use std::io::Write;

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Vacuum, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::{MolTypeAtom, Topology};
use newen_forces::nonbonded::XtbInteraction;
use newen_forces::provider::PotentialProvider;
use newen_integrators::{Integrator, LeapFrog};
use rand_distr::{Distribution, Normal};

const N_TRAJECTORIES: usize = 10;
const STEPS_PER_TRAJECTORY: usize = 90;
const DT: f64 = 2e-4;

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

fn equilibrium_positions_nm() -> Vec<Vec3> {
    vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0758602, 0.0, 0.0504284),
        Vec3::new(0.0758602, 0.0, -0.0504284),
    ]
}

fn main() {
    let out_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "/tmp/qm_water_training_data.txt".to_string());

    let topo = water_topology();
    let n_atoms = 3;
    let atomic_numbers = vec![8i64, 1, 1];
    let periodicity = Periodicity::Vacuum(Vacuum);
    let region = AtomSelection::all(n_atoms);

    let work_dir = std::env::temp_dir().join("newen_qm_training_data_gen");
    let mut interaction =
        XtbInteraction::new(work_dir, 2, 0, 1, atomic_numbers.clone()).expect("work_dir creatable");

    let potential_energy =
        |conf: &Configuration, interaction: &mut XtbInteraction| -> (f64, Vec<Vec3>) {
            let index = ConfigurationSpatialIndex::new(conf, &periodicity);
            let c = interaction
                .contribute(&region, &topo, conf, &index)
                .expect("xtb calculation should succeed");
            let mut forces = vec![Vec3::ZERO; n_atoms];
            for (i, f) in c.forces {
                forces[i] = f;
            }
            (c.energy, forces)
        };

    let mut out = std::fs::File::create(&out_path)
        .unwrap_or_else(|e| panic!("cannot create {out_path}: {e}"));
    writeln!(
        out,
        "{n_atoms} {} {} {}",
        atomic_numbers[0], atomic_numbers[1], atomic_numbers[2]
    )
    .unwrap();

    let mut rng = rand::thread_rng();
    // Small perturbation around equilibrium (nm) and a modest random initial velocity (nm/ps) —
    // enough to explore a spread of bond lengths/angles across trajectories without starting so
    // far from equilibrium that the geometry is chemically nonsensical.
    let pos_perturb = Normal::new(0.0, 0.01).unwrap();
    let vel_scale = Normal::new(0.0, 0.15).unwrap();

    let mut total_frames = 0usize;
    for traj in 0..N_TRAJECTORIES {
        let base = equilibrium_positions_nm();
        let mut conf = Configuration::new(n_atoms, 1, 1);
        conf.current_mut().pos = base
            .iter()
            .map(|&p| {
                p + Vec3::new(
                    pos_perturb.sample(&mut rng),
                    pos_perturb.sample(&mut rng),
                    pos_perturb.sample(&mut rng),
                )
            })
            .collect();
        conf.current_mut().vel = (0..n_atoms)
            .map(|_| {
                Vec3::new(
                    vel_scale.sample(&mut rng),
                    vel_scale.sample(&mut rng),
                    vel_scale.sample(&mut rng),
                )
            })
            .collect();
        conf.copy_current_to_old();

        let (e0, f0) = potential_energy(&conf, &mut interaction);
        conf.current_mut().force = f0;
        conf.copy_current_to_old();
        write_frame(&mut out, e0, &conf.current().pos, &conf.current().force);
        total_frames += 1;

        let mut integrator = LeapFrog::new();
        for _ in 0..STEPS_PER_TRAJECTORY {
            integrator.step(DT, &topo, &mut conf);
            let (pe, f) = potential_energy(&conf, &mut interaction);
            conf.current_mut().force = f;
            write_frame(&mut out, pe, &conf.current().pos, &conf.current().force);
            total_frames += 1;
        }

        eprintln!(
            "trajectory {}/{N_TRAJECTORIES}: {} frames, E0={e0:.3} kJ/mol",
            traj + 1,
            STEPS_PER_TRAJECTORY + 1
        );
    }

    eprintln!("wrote {total_frames} frames to {out_path}");
}

fn write_frame(out: &mut std::fs::File, energy: f64, pos: &[Vec3], force: &[Vec3]) {
    write!(out, "{energy:.10}").unwrap();
    for i in 0..pos.len() {
        write!(
            out,
            " {:.10} {:.10} {:.10} {:.10} {:.10} {:.10}",
            pos[i].x, pos[i].y, pos[i].z, force[i].x, force[i].y, force[i].z
        )
        .unwrap();
    }
    writeln!(out).unwrap();
}
