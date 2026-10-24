//! `ProviderOrchestratorAlgorithm` — a [`ProviderOrchestrator`] driven by the real
//! [`AlgorithmSequence`](newen_core::algorithm::AlgorithmSequence) step loop.
//!
//! Every existing `ProviderOrchestrator` caller evaluates it directly against a loaded
//! `Topology`/`Configuration` (test harnesses, `xtb_nve_loop.rs`'s bespoke `LeapFrog` loop);
//! nothing runs it as part of the same `Algorithm` pipeline `Forcefield` uses. This wraps it as
//! an `Algorithm` so it can be pushed into a real sequence right after `Forcefield`, the same way
//! `md.rs`/`pyo3-newen` already build one.
//!
//! Unlike `Forcefield`, which *overwrites* `conf.current_mut().force` (it computes the whole
//! classical force), this **adds** to whatever's already there — `Forcefield` (or nothing) must
//! run first in the sequence.
//!
//! **Energy accounting:** the summed contribution goes into `state.energies.special_total`, not
//! `potential_total`. `update_potential_total()` (called both by `Forcefield` and, again, by the
//! final `EnergyCalculation` algorithm) only sums bond/angle/dihedral/improper/cross_dihedral/lj/
//! crf/ls/sasa — anything added to `potential_total` directly here would be silently wiped out by
//! that later recompute. `Energy::total()` adds `special_total` on top separately, so this is
//! where an additive term that doesn't fit the standard breakdown actually survives to the
//! reported total (same slot position restraints already use). A dedicated per-term field is
//! real follow-up work, not done here.
//!
//! **Overwrites, not accumulates, `special_total`** — unlike `force`, which `Forcefield`
//! unconditionally recomputes from scratch every step (so `+=` after it means "this step's
//! classical force plus this step's orchestrator force"), nothing resets `special_total` to a
//! fresh per-step baseline when position restraints aren't configured — `Forcefield` only *sets*
//! it (`forcefield.rs`'s posres block), conditionally. `+=` here would silently compound across
//! steps whenever that block doesn't run. Overwriting means this algorithm currently claims
//! exclusive ownership of `special_total`; running it alongside active position restraints isn't
//! supported yet (a real follow-up, not attempted here).

use newen_core::algorithm::{Algorithm, SimulationState};
use newen_core::configuration::{BoxType, Configuration};
use newen_core::math::{Periodicity, Rectangular, Triclinic};
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_core::topology::Topology;

use crate::orchestrator::ProviderOrchestrator;

pub struct ProviderOrchestratorAlgorithm {
    orchestrator: ProviderOrchestrator,
    periodicity: Periodicity,
}

impl ProviderOrchestratorAlgorithm {
    pub fn new(orchestrator: ProviderOrchestrator, periodicity: Periodicity) -> Self {
        Self {
            orchestrator,
            periodicity,
        }
    }
}

impl Algorithm for ProviderOrchestratorAlgorithm {
    fn apply(
        &mut self,
        topo: &Topology,
        conf: &mut Configuration,
        _sim: &SimulationState,
    ) -> Result<(), String> {
        // Refresh periodicity from the current box, mirroring `Forcefield::apply()` exactly (NPT
        // scales the box at the end of the previous step).
        if !matches!(self.periodicity, Periodicity::Vacuum(_)) {
            let box_cfg = &conf.current().box_config;
            match box_cfg.box_type {
                BoxType::Rectangular => {
                    let dims = box_cfg.dimensions();
                    if dims.x > 0.0 && dims.y > 0.0 && dims.z > 0.0 {
                        self.periodicity = Periodicity::Rectangular(Rectangular::new(dims));
                    }
                },
                BoxType::Triclinic | BoxType::TruncatedOctahedral => {
                    self.periodicity = Periodicity::Triclinic(Triclinic::new(box_cfg.vectors));
                },
                _ => {},
            }
        }

        let index = ConfigurationSpatialIndex::new(conf, &self.periodicity);
        let (contribution, per_term) = self
            .orchestrator
            .evaluate_with_terms(topo, conf, &index)
            .map_err(|e| e.to_string())?;

        let state = conf.current_mut();
        for (i, f) in contribution.forces {
            state.force[i] += f;
        }
        state.energies.special_total = contribution.energy;
        state.energies.term_energies = per_term;
        state.virial_tensor += contribution.virial;

        Ok(())
    }

    fn name(&self) -> &str {
        "ProviderOrchestrator"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::provider::{Contribution, Embedding, PotentialProvider, ProviderError};
    use newen_core::configuration::Configuration;
    use newen_core::math::{Vacuum, Vec3};
    use newen_core::selection::AtomSelection;
    use newen_core::spatial_index::SpatialIndex;

    /// A provider that always returns a fixed contribution, to check the `Algorithm` wrapper's
    /// own bookkeeping (`+=` not `=`, `special_total` not `potential_total`) in isolation from
    /// any real physics.
    struct FixedProvider {
        energy: f64,
        force: Vec3,
    }

    impl PotentialProvider for FixedProvider {
        fn contribute(
            &mut self,
            region: &AtomSelection,
            _topo: &Topology,
            _conf: &Configuration,
            _neigh: &dyn SpatialIndex,
        ) -> Result<Contribution, ProviderError> {
            Ok(Contribution {
                energy: self.energy,
                forces: region.iter().map(|i| (i, self.force)).collect(),
                virial: newen_core::math::Mat3::ZERO,
                extra: Default::default(),
            })
        }

        fn name(&self) -> &str {
            "fixed"
        }

        fn embedding(&self) -> Embedding {
            Embedding::None
        }
    }

    /// A minimal real 2-atom topology — `Topology::new()` alone has zero atoms
    /// (`num_atoms()` is derived from `moltypes[0].atoms`, not just allocated Vecs).
    fn two_atom_topology() -> Topology {
        use newen_core::topology::MolTypeAtom;
        let mut topo = Topology::new();
        topo.charge = vec![0.0, 0.0];
        topo.iac = vec![0, 0];
        topo.mass = vec![1.0, 1.0];
        topo.inverse_mass = vec![1.0, 1.0];
        topo.exclusions = vec![Vec::new(), Vec::new()];
        topo.moltypes[0].atoms = (0..2)
            .map(|i| MolTypeAtom {
                name: format!("A{i}"),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 0,
                mass: 1.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            })
            .collect();
        topo
    }

    #[test]
    fn adds_to_existing_force_but_overwrites_special_total() {
        let topo = two_atom_topology();
        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos = vec![Vec3::ZERO, Vec3::new(0.3, 0.0, 0.0)];
        // Simulate a preceding `Forcefield` having already written this step's classical force
        // (always a fresh overwrite, per `forcefield.rs`) — this algorithm must add on top.
        conf.current_mut().force = vec![Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO];
        // ...and a *stale* `special_total` left over from an earlier step (nothing resets it
        // when position restraints aren't configured) — this algorithm must overwrite it, not
        // compound on top, or `Energy::total()` would drift upward every step regardless of
        // physics. See module docs.
        conf.current_mut().energies.special_total = 999.0;

        let mut orchestrator = ProviderOrchestrator::new();
        orchestrator.register(
            AtomSelection::all(2),
            Box::new(FixedProvider {
                energy: 10.0,
                force: Vec3::new(0.0, 2.0, 0.0),
            }),
        );
        let mut algorithm =
            ProviderOrchestratorAlgorithm::new(orchestrator, Periodicity::Vacuum(Vacuum));
        let sim = SimulationState::new(0.001, 1);

        algorithm.apply(&topo, &mut conf, &sim).unwrap();

        assert_eq!(conf.current().force[0], Vec3::new(1.0, 2.0, 0.0));
        assert_eq!(conf.current().force[1], Vec3::new(0.0, 2.0, 0.0));
        assert_eq!(conf.current().energies.special_total, 10.0);
    }
}
