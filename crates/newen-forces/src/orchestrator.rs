//! `ProviderOrchestrator` — sums `Vec<Box<dyn PotentialProvider>>` into one [`Contribution`].
//!
//! ** P2.8-2.** Every piece up to here (`provider.rs`, `zones.rs`, `LjCrfInteraction`,
//! `ElectrostaticEmbedding`, `SchNetInteraction`) was validated in isolation; nothing iterated
//! more than one provider. This is that missing piece: register `(region, provider)` pairs,
//! evaluate all of them against the same shared `Topology`/`Configuration`/`SpatialIndex`, and
//! sum the results — the step that turns "a QM/MM force term" into "a QM/MM-shaped force
//! evaluation."
//!
//! # The index-validation this exists to do (P2.6 review finding)
//!
//! `Contribution` is a scattered accumulator rather than `&mut Configuration` specifically so a
//! caller *can* check which indices a provider touched before trusting them (`provider.rs`
//! module docs, P5). This module is that caller. A provider that declares
//! [`Embedding::None`] (`provider.rs`) is asserting "I only touch atoms in my own region" — the
//! orchestrator enforces that assertion instead of trusting it, exactly the failure mode the
//! scattered-`Contribution` design was chosen to make checkable. A provider with
//! [`Embedding::Mechanical`] or [`Embedding::Electrostatic`] is exempted from that check by
//! definition: touching atoms outside `region` is the whole point (e.g.
//! [`crate::nonbonded::ElectrostaticEmbedding`] puts forces on MM atoms it was never handed as
//! its `region`).

use crate::provider::{Contribution, Embedding, PotentialProvider, ProviderError, ProviderExtra};
use newen_core::configuration::Configuration;
use newen_core::math::{Mat3, Vec3};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::SpatialIndex;
use newen_core::topology::Topology;

/// Registers providers over their regions and sums their [`Contribution`]s.
///
/// Deliberately dumb: no scheduling, no caching, no parallelism — evaluation order is
/// registration order, and correctness (this PLAN step's whole point) comes from the per-provider
/// index check, not from anything clever. Async/cancellable evaluation for slow external
/// providers is explicitly out of scope (`provider.rs` module docs already flag this as
/// deliberately unsolved).
#[derive(Default)]
pub struct ProviderOrchestrator {
    entries: Vec<(AtomSelection, Box<dyn PotentialProvider>)>,
    /// One label per entry — the name each provider's energy is reported under
    /// (`evaluate_with_terms`); defaults to the provider's own `name()`.
    labels: Vec<String>,
}

impl ProviderOrchestrator {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a provider scoped to `region`. Evaluated in registration order.
    pub fn register(&mut self, region: AtomSelection, provider: Box<dyn PotentialProvider>) {
        let label = provider.name().to_string();
        self.register_labelled(label, region, provider);
    }

    /// `register`, reporting the provider's energy under `label` (a run's registry name for the
    /// term — 3.9 G10) instead of the provider's own `name`.
    pub fn register_labelled(
        &mut self,
        label: impl Into<String>,
        region: AtomSelection,
        provider: Box<dyn PotentialProvider>,
    ) {
        self.labels.push(label.into());
        self.entries.push((region, provider));
    }

    /// The labels, in evaluation order.
    pub fn labels(&self) -> &[String] {
        &self.labels
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Evaluate every registered provider and sum the results.
    ///
    /// Fails fast on the first provider error or index-contract violation — a QM/MM step with a
    /// broken term shouldn't silently produce a partial, wrong energy.
    pub fn evaluate(
        &mut self,
        topo: &Topology,
        conf: &Configuration,
        neigh: &dyn SpatialIndex,
    ) -> Result<Contribution, ProviderError> {
        self.evaluate_with_terms(topo, conf, neigh).map(|(c, _)| c)
    }

    /// `evaluate`, also returning every provider's own energy under its label, in evaluation
    /// order — the per-term view the summed `Contribution` erases (3.9 G10).
    pub fn evaluate_with_terms(
        &mut self,
        topo: &Topology,
        conf: &Configuration,
        neigh: &dyn SpatialIndex,
    ) -> Result<(Contribution, Vec<(String, f64)>), ProviderError> {
        let n_atoms = topo.num_atoms();
        let mut energy = 0.0;
        let mut force_acc = vec![Vec3::ZERO; n_atoms];
        let mut touched = vec![false; n_atoms];
        let mut virial = Mat3::ZERO;
        let mut per_term = Vec::with_capacity(self.entries.len());

        for ((region, provider), label) in self.entries.iter_mut().zip(&self.labels) {
            let contribution = provider.contribute(region, topo, conf, neigh)?;
            per_term.push((label.clone(), contribution.energy));

            if provider.embedding() == Embedding::None {
                let mut in_region = vec![false; n_atoms];
                for i in region.iter() {
                    in_region[i] = true;
                }
                for &(i, _) in &contribution.forces {
                    if !in_region[i] {
                        return Err(ProviderError::InvalidRegion(format!(
                            "provider \"{}\" declares Embedding::None but placed a force on atom \
                             {i}, outside the region it was given — this is exactly the \
                             unchecked-mutable-access failure mode Contribution exists to catch",
                            provider.name()
                        )));
                    }
                }
            }

            energy += contribution.energy;
            for &(i, f) in &contribution.forces {
                force_acc[i] += f;
                touched[i] = true;
            }
            virial += contribution.virial;
        }

        let forces = (0..n_atoms)
            .filter(|&i| touched[i])
            .map(|i| (i, force_acc[i]))
            .collect();

        Ok((
            Contribution {
                energy,
                forces,
                virial,
                extra: ProviderExtra::default(),
            },
            per_term,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nonbonded::LjCrfInteraction;
    use newen_core::configuration::Configuration;
    use newen_core::math::{Periodicity, Vacuum};
    use newen_core::spatial_index::ConfigurationSpatialIndex;
    use newen_core::topology::{LJParameters as TopoLJParameters, MolTypeAtom};

    fn two_atom_fixture() -> (Topology, Configuration) {
        let mut topo = Topology::new();
        topo.charge = vec![0.4, -0.4];
        topo.iac = vec![0, 0];
        topo.mass = vec![12.0, 12.0];
        topo.inverse_mass = vec![1.0 / 12.0; 2];
        topo.exclusions = vec![newen_core::topology::Exclusions::default(); 2];
        topo.lj_parameters = vec![vec![TopoLJParameters {
            c6: 0.0026,
            c12: 0.0000023,
            cs6: 0.0026,
            cs12: 0.0000023,
        }]];
        topo.moltypes[0].atoms = (0..2)
            .map(|i| MolTypeAtom {
                name: format!("A{i}"),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 0,
                mass: 12.0,
                charge: topo.charge[i],
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            })
            .collect();

        let mut conf = Configuration::new(2, 1, 1);
        conf.current_mut().pos = vec![Vec3::ZERO, Vec3::new(0.35, 0.0, 0.0)];
        (topo, conf)
    }

    /// Transparency: one provider through the orchestrator must exactly match calling that same
    /// provider directly — the orchestrator's own bookkeeping (accumulation, index check) must
    /// not perturb the result. This is the Rust-internal half of the P2.8-2 exit criterion; the
    /// real-gromosXX-oracle half lives in `test_orchestrator_reference.rs`.
    #[test]
    fn single_provider_matches_direct_call() {
        let (topo, conf) = two_atom_fixture();
        let region = AtomSelection::all(2);
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let direct = LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0)
            .contribute(&region, &topo, &conf, &index)
            .unwrap();

        let mut orchestrator = ProviderOrchestrator::new();
        orchestrator.register(region, Box::new(LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0)));
        let orchestrated = orchestrator.evaluate(&topo, &conf, &index).unwrap();

        assert_eq!(orchestrated.energy, direct.energy);
        assert_eq!(orchestrated.forces.len(), direct.forces.len());
        for (i, f) in &direct.forces {
            let of = orchestrated
                .forces
                .iter()
                .find(|&&(j, _)| j == *i)
                .unwrap()
                .1;
            assert_eq!(of, *f);
        }
    }

    /// Two providers on disjoint regions: the orchestrator must sum both, and each individually
    /// obeys its own `Embedding::None` contract (all forces stay within its own region), so
    /// nothing should be flagged.
    #[test]
    fn two_disjoint_providers_sum_correctly() {
        let (topo, conf) = two_atom_fixture();
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        // Two providers both scoped to the whole system (LjCrfInteraction has no way to scope
        // narrower than "the pairs it finds"), so this checks summation, not disjoint regions.
        let mut orchestrator = ProviderOrchestrator::new();
        orchestrator.register(
            AtomSelection::all(2),
            Box::new(LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0)),
        );
        orchestrator.register(
            AtomSelection::all(2),
            Box::new(LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0)),
        );
        let doubled = orchestrator.evaluate(&topo, &conf, &index).unwrap();

        let single = LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0)
            .contribute(&AtomSelection::all(2), &topo, &conf, &index)
            .unwrap();

        assert_eq!(doubled.energy, 2.0 * single.energy);
    }

    /// A provider that lies about its embedding — claims `Embedding::None` but puts a force
    /// outside its region — must be caught, not silently accumulated. This is the P2.6 review
    /// finding this whole module exists to enforce.
    #[test]
    fn embedding_none_provider_touching_outside_region_is_rejected() {
        struct Lying;
        impl PotentialProvider for Lying {
            fn contribute(
                &mut self,
                _region: &AtomSelection,
                _topo: &Topology,
                _conf: &Configuration,
                _neigh: &dyn SpatialIndex,
            ) -> Result<Contribution, ProviderError> {
                // Claims region {0} but puts a force on atom 1.
                Ok(Contribution {
                    energy: 1.0,
                    forces: vec![(1, Vec3::new(1.0, 0.0, 0.0))],
                    virial: Mat3::ZERO,
                    extra: ProviderExtra::default(),
                })
            }
            fn name(&self) -> &str {
                "Lying"
            }
            // Embedding::None is the trait default — deliberately not overridden.
        }

        let (topo, conf) = two_atom_fixture();
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let region = AtomSelection::from_indices(vec![0], 2).unwrap();

        let mut orchestrator = ProviderOrchestrator::new();
        orchestrator.register(region, Box::new(Lying));

        let result = orchestrator.evaluate(&topo, &conf, &index);
        assert!(matches!(result, Err(ProviderError::InvalidRegion(_))));
    }

    /// The same lying provider is fine once it honestly declares `Embedding::Electrostatic` —
    /// touching atoms outside `region` is then expected, not a contract violation.
    #[test]
    fn embedding_electrostatic_provider_may_touch_outside_region() {
        struct Honest;
        impl PotentialProvider for Honest {
            fn contribute(
                &mut self,
                _region: &AtomSelection,
                _topo: &Topology,
                _conf: &Configuration,
                _neigh: &dyn SpatialIndex,
            ) -> Result<Contribution, ProviderError> {
                Ok(Contribution {
                    energy: 1.0,
                    forces: vec![(1, Vec3::new(1.0, 0.0, 0.0))],
                    virial: Mat3::ZERO,
                    extra: ProviderExtra::default(),
                })
            }
            fn name(&self) -> &str {
                "Honest"
            }
            fn embedding(&self) -> Embedding {
                Embedding::Electrostatic
            }
        }

        let (topo, conf) = two_atom_fixture();
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let region = AtomSelection::from_indices(vec![0], 2).unwrap();

        let mut orchestrator = ProviderOrchestrator::new();
        orchestrator.register(region, Box::new(Honest));

        let result = orchestrator.evaluate(&topo, &conf, &index).unwrap();
        assert_eq!(result.energy, 1.0);
        assert_eq!(result.forces.len(), 1);
    }

    #[test]
    fn empty_orchestrator_is_zero() {
        let (topo, conf) = two_atom_fixture();
        let periodicity = Periodicity::Vacuum(Vacuum);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let mut orchestrator = ProviderOrchestrator::new();
        assert!(orchestrator.is_empty());
        let result = orchestrator.evaluate(&topo, &conf, &index).unwrap();
        assert_eq!(result.energy, 0.0);
        assert!(result.forces.is_empty());
    }
}
