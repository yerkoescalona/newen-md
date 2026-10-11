//! `LjCrfInteraction` — the classical LJ+CRF nonbonded term as a [`PotentialProvider`].
//!
//! Named `*Interaction`, not `*Provider`: gromosXX's own C++ base class for force terms is
//! `Interaction` (`Nonbonded_Interaction`, `QMMM_Interaction` — gromosXX shapes QM/MM as one too).
//! `PotentialProvider` stays the trait name (already fixed in `architecture.md`'s
//! five-contract taxonomy); concrete implementations read as GROMOS terms.
//!
//! This wraps the *exact* existing [`lj_crf_innerloop_novirial`] — the same function
//! [`crate::energy::single_point_energy`] already calls — to prove the provider shape is
//! transparent, not a rewrite. Solvent's charge-group "sentinel pair per molecule pair"
//! pairlist shape (see `gromos_core::pairlist` docs) doesn't fit the generic
//! [`SpatialIndex`] query, which returns literal atom-index pairs — so **this first pass
//! covers solute-solute (atom-level) LJ+CRF only**. Solvent needs its own adapter/expansion
//! and is deferred to when providers are actually wired into the MD loop.

use gromos_core::configuration::{BoxType, Configuration};
use gromos_core::math::{Periodicity, Rectangular, Vacuum, Vec3};
use gromos_core::selection::AtomSelection;
use gromos_core::spatial_index::SpatialIndex;
use gromos_core::topology::Topology;

use crate::nonbonded::{
    lj_crf_innerloop_novirial, CRFParameters, EnergyGroups, ForceStorage, LJParamMatrix,
    LJParameters,
};
use crate::provider::{Contribution, PotentialProvider, ProviderError, ProviderExtra};
use crate::zones::{Zone, ZonePartition};
use gromos_core::math::Mat3;

/// Classical LJ + reaction-field Coulomb, evaluated over solute-solute atom pairs.
#[derive(Debug, Clone)]
pub struct LjCrfInteraction {
    /// Nonbonded cutoff (nm).
    pub cutoff: f64,
    /// Interior dielectric constant (always 1.0 in GROMOS).
    pub epsilon: f64,
    /// Reaction-field dielectric constant.
    pub rf_epsilon: f64,
    /// Reaction-field inverse Debye screening length (nm^-1).
    pub rf_kappa: f64,
    /// ** P2.8-1.** When set, pairs whose [`ZonePartition::classical_should_evaluate`]
    /// is `false` (owned by a QM/ML provider or by electrostatic embedding, see
    /// [`crate::zones`]) are skipped here — the double-counting contract `zones.rs` states but
    /// nothing previously consulted. `None` means "no partitioning" (every pair is classical),
    /// which is what every existing reference test and call site still gets.
    pub zone_partition: Option<ZonePartition>,
}

impl LjCrfInteraction {
    /// Construct from the reaction-field parameters (mirrors [`crate::energy::EnergyParams`]).
    /// No zone partition — use [`Self::with_zone_partition`] to add one.
    pub fn new(cutoff: f64, epsilon: f64, rf_epsilon: f64, rf_kappa: f64) -> Self {
        Self {
            cutoff,
            epsilon,
            rf_epsilon,
            rf_kappa,
            zone_partition: None,
        }
    }

    /// Skip pairs owned by a QM/ML provider or by electrostatic embedding (P2.8-1).
    pub fn with_zone_partition(mut self, partition: ZonePartition) -> Self {
        self.zone_partition = Some(partition);
        self
    }
}

fn periodicity_of(conf: &Configuration) -> Periodicity {
    let box_config = &conf.current().box_config;
    match box_config.box_type {
        BoxType::Vacuum => Periodicity::Vacuum(Vacuum),
        _ => Periodicity::Rectangular(Rectangular::new(box_config.dimensions())),
    }
}

impl PotentialProvider for LjCrfInteraction {
    fn contribute(
        &mut self,
        region: &AtomSelection,
        topo: &Topology,
        conf: &Configuration,
        neigh: &dyn SpatialIndex,
    ) -> Result<Contribution, ProviderError> {
        let n_atoms = topo.num_atoms();
        if n_atoms == 0 {
            return Err(ProviderError::InvalidRegion(
                "topology has no atoms".to_string(),
            ));
        }

        let charges = topo.charge.clone();
        let iac_u32: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
        let lj_nested: Vec<Vec<LJParameters>> = topo
            .lj_parameters
            .iter()
            .map(|row| row.iter().map(LJParameters::from).collect())
            .collect();
        let lj_mat = LJParamMatrix::from_nested(&lj_nested);
        let crf = CRFParameters::new(self.cutoff, self.epsilon, self.rf_epsilon, self.rf_kappa);
        let periodicity = periodicity_of(conf);

        // Shift vectors aren't needed here: lj_crf_innerloop_novirial recomputes the
        // minimum image itself from `periodicity` for each pair.
        let raw_pairs = neigh.neighbor_pairs(region, self.cutoff);

        let pairlist_u32: Vec<(u32, u32)> = raw_pairs
            .iter()
            .filter(|(i, j, _)| {
                self.zone_partition
                    .as_ref()
                    .is_none_or(|p| p.classical_should_evaluate(*i, *j))
            })
            .map(|(i, j, _)| (*i as u32, *j as u32))
            .collect();

        let mut storage = ForceStorage::new(n_atoms);
        if !pairlist_u32.is_empty() {
            lj_crf_innerloop_novirial(
                &conf.current().pos,
                &charges,
                &iac_u32,
                &pairlist_u32,
                &lj_mat,
                &crf,
                &periodicity,
                gromos_core::units::four_pi_eps_i,
                &mut storage,
                EnergyGroups::single(),
            );
        }

        // ** P2.8-4.** Pairs whose combined LJ+CRF pass was skipped above
        // (`classical_should_evaluate` false — Provider- or Embedding-owned) but whose LJ term
        // gromosXX still computes classically (`ZonePartition::lj_owner`, sourced from the real
        // `QMMM_Interaction::modify_exclusions`): inner-outer pairs always qualify; inner-inner /
        // inner-buffer pairs qualify only when `QMLJ` is on. Evaluated with CRF zeroed by running
        // the same innerloop against a charge array with every `Zone::Inner` atom's charge set to
        // 0 — safe because every `lj_only_should_evaluate` pair has an `Inner` endpoint, and an
        // `Inner` atom is never part of a `classical_should_evaluate` pair (see `zones.rs`), so
        // this can't perturb the pass above or double-count anything.
        if let Some(partition) = &self.zone_partition {
            let lj_only_pairs: Vec<(u32, u32)> = raw_pairs
                .iter()
                .filter(|(i, j, _)| partition.lj_only_should_evaluate(*i, *j))
                .map(|(i, j, _)| (*i as u32, *j as u32))
                .collect();
            if !lj_only_pairs.is_empty() {
                let mut lj_only_charges = charges.clone();
                for atom in partition.atoms_in(Zone::Inner) {
                    lj_only_charges[atom] = 0.0;
                }
                lj_crf_innerloop_novirial(
                    &conf.current().pos,
                    &lj_only_charges,
                    &iac_u32,
                    &lj_only_pairs,
                    &lj_mat,
                    &crf,
                    &periodicity,
                    gromos_core::units::four_pi_eps_i,
                    &mut storage,
                    EnergyGroups::single(),
                );
            }
        }

        let forces = (0..n_atoms)
            .map(|i| (i, storage.forces[i]))
            .filter(|(_, f)| *f != Vec3::ZERO)
            .collect();

        Ok(Contribution {
            energy: storage.e_lj + storage.e_crf,
            forces,
            virial: Mat3::ZERO, // "novirial" path never populates storage.virial
            extra: ProviderExtra::default(),
        })
    }

    fn name(&self) -> &str {
        "LJ+CRF"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::energy::{single_point_energy, EnergyParams};
    use gromos_core::spatial_index::ConfigurationSpatialIndex;
    use gromos_core::topology::{LJParameters as TopoLJParameters, MolTypeAtom};

    /// Four solute atoms, no exclusions/1-4/chargegroups, in vacuum: `single_point_energy`'s
    /// `.lj`/`.crf` totals then equal the pure pairlist LJ+CRF sum with no self/excluded/1-4
    /// contamination, making it a genuine second, independent computation path (its own
    /// chargegroup-aware `StandardPairlistAlgorithm` construction, vs. this module's
    /// brute-force `ConfigurationSpatialIndex`) to cross-check against.
    fn four_atom_fixture() -> (Topology, Vec<Vec3>) {
        let mut topo = Topology::new();
        let charges = [0.4, -0.4, 0.3, -0.3];
        let iac = [0usize, 0, 0, 0];

        topo.charge = charges.to_vec();
        topo.iac = iac.to_vec();
        topo.mass = vec![12.0; 4];
        topo.inverse_mass = vec![1.0 / 12.0; 4];
        topo.exclusions = vec![gromos_core::topology::Exclusions::default(); 4];
        topo.lj_parameters = vec![vec![TopoLJParameters {
            c6: 0.0026,
            c12: 0.0000023,
            cs6: 0.0026,
            cs12: 0.0000023,
        }]];

        topo.moltypes[0].atoms = (0..4)
            .map(|i| MolTypeAtom {
                name: format!("A{i}"),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: iac[i],
                mass: 12.0,
                charge: charges[i],
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            })
            .collect();

        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(0.0, 0.35, 0.0),
            Vec3::new(0.3, 0.35, 0.0),
        ];
        (topo, positions)
    }

    #[test]
    fn wrapper_transparency_matches_direct_call() {
        let (topo, positions) = four_atom_fixture();
        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos = positions.clone();

        let region = AtomSelection::all(4);
        let periodicity = periodicity_of(&conf);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let mut interaction = LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0);
        let contribution = interaction
            .contribute(&region, &topo, &conf, &index)
            .unwrap();

        // Direct call: same function, same inputs, computed by hand outside the wrapper.
        let charges = topo.charge.clone();
        let iac_u32: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
        let lj_mat = LJParamMatrix::from_nested(&[topo.lj_parameters[0]
            .iter()
            .map(LJParameters::from)
            .collect()]);
        let crf = CRFParameters::new(1.4, 1.0, 78.5, 0.0);
        let pairlist_u32: Vec<(u32, u32)> = index
            .neighbor_pairs(&region, 1.4)
            .into_iter()
            .map(|(i, j, _)| (i as u32, j as u32))
            .collect();
        let mut storage = ForceStorage::new(4);
        lj_crf_innerloop_novirial(
            &positions,
            &charges,
            &iac_u32,
            &pairlist_u32,
            &lj_mat,
            &crf,
            &periodicity,
            gromos_core::units::four_pi_eps_i,
            &mut storage,
            EnergyGroups::single(),
        );

        assert_eq!(contribution.energy, storage.e_lj + storage.e_crf);
        for (i, f) in &contribution.forces {
            assert_eq!(*f, storage.forces[*i]);
        }
    }

    #[test]
    fn cross_checks_against_single_point_energy_oracle() {
        let (topo, positions) = four_atom_fixture();
        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos = positions.clone();

        let region = AtomSelection::all(4);
        let periodicity = periodicity_of(&conf);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let mut interaction = LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0);
        let contribution = interaction
            .contribute(&region, &topo, &conf, &index)
            .unwrap();

        // Independent path: single_point_energy builds its own StandardPairlistAlgorithm
        // pairlist (chargegroup/atom-based) and calls the same innerloop separately.
        let params = EnergyParams {
            cutoff: 1.4,
            epsilon: 1.0,
            rf_epsilon: 78.5,
            rf_kappa: 0.0,
            pairlist_freq: 1,
            ntf: [true, true, true, true],
            atoms_per_solvent: 3,
            quartic_bonds: true,
        };
        let oracle = single_point_energy(&topo, &positions, Vec3::ZERO, &params);

        // single_point_energy's `.crf` also includes gromosXX's unconditional per-atom RF
        // self-energy (added by `rf_excluded_interactions` regardless of whether there are
        // any exclusions — it's not an "excluded pair" correction, it's the self-term every
        // solute atom gets). `LjCrfInteraction` deliberately doesn't compute that (it's a
        // separate GROMOS term, not part of the pairlist LJ+CRF loop), so back it out here
        // to compare like with like: pure pairlist LJ+CRF vs. pure pairlist LJ+CRF.
        let crf = CRFParameters::new(1.4, 1.0, 78.5, 0.0);
        let self_energy: f64 = topo
            .charge
            .iter()
            .map(|q| -0.5 * q * q * gromos_core::units::four_pi_eps_i * crf.crf_cut)
            .sum();

        assert!((contribution.energy - (oracle.lj + oracle.crf - self_energy)).abs() < 1e-10);
    }

    /// ** P2.8-1 exit criterion, isolated.** Partition the 4-atom fixture into inner
    /// {0,1} / outer {2,3}: pair (2,3) is the only [`crate::zones::PairOwner::Classical`] pair,
    /// the rest are Provider- or Embedding-owned and must be skipped. The partitioned energy
    /// must equal exactly the direct-call energy of the single surviving pair — proving the
    /// filter drops precisely what the ownership table says it should, no more, no less.
    #[test]
    fn zone_partition_restricts_to_classically_owned_pairs_only() {
        use crate::zones::ZonePartition;

        let (topo, positions) = four_atom_fixture();
        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos = positions.clone();

        let region = AtomSelection::all(4);
        let periodicity = periodicity_of(&conf);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        let partition = ZonePartition::new(4, &[0, 1], &[]); // inner {0,1}, outer {2,3}
        let mut interaction =
            LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0).with_zone_partition(partition.clone());
        let partitioned = interaction
            .contribute(&region, &topo, &conf, &index)
            .unwrap();

        // Hand-built oracle matching P2.8-4's two-pass rule: (2,3) is outer-outer, full
        // classical LJ+CRF; (0,2)/(0,3)/(1,2)/(1,3) are inner-outer, LJ-only (CRF always
        // classical-owned but zeroed here per gromosXX's real modify_exclusions rule — this
        // provider computes it as a supplement, see the module's contribute() comment); (0,1) is
        // inner-inner, fully provider-owned (QMLJ defaults off), absent from both passes.
        let charges = topo.charge.clone();
        let iac_u32: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
        let lj_mat = LJParamMatrix::from_nested(&[topo.lj_parameters[0]
            .iter()
            .map(LJParameters::from)
            .collect()]);
        let crf = CRFParameters::new(1.4, 1.0, 78.5, 0.0);
        let mut storage = ForceStorage::new(4);
        lj_crf_innerloop_novirial(
            &positions,
            &charges,
            &iac_u32,
            &[(2u32, 3u32)],
            &lj_mat,
            &crf,
            &periodicity,
            gromos_core::units::four_pi_eps_i,
            &mut storage,
            EnergyGroups::single(),
        );
        let mut lj_only_charges = charges.clone();
        lj_only_charges[0] = 0.0;
        lj_only_charges[1] = 0.0; // zero both inner atoms' charges: LJ-only for the supplement
        lj_crf_innerloop_novirial(
            &positions,
            &lj_only_charges,
            &iac_u32,
            &[(0u32, 2u32), (0u32, 3u32), (1u32, 2u32), (1u32, 3u32)],
            &lj_mat,
            &crf,
            &periodicity,
            gromos_core::units::four_pi_eps_i,
            &mut storage,
            EnergyGroups::single(),
        );

        assert_eq!(partitioned.energy, storage.e_lj + storage.e_crf);
        // Pair (0,1) — inner-inner, fully provider-owned — must be the only one absent from both
        // the combined classical pass and the LJ-only supplement.
        assert!(!partition.classical_should_evaluate(0, 1));
        assert!(!partition.lj_only_should_evaluate(0, 1));
    }

    /// Isolates the P2.8-4 supplement itself: with no buffer, an inner-outer pair must get
    /// exactly its LJ term (matching a direct `lj_crf_interaction` call with `q_prod` forced to
    /// 0) and exactly zero CRF — proving the zeroed-charge trick used inside `contribute()`
    /// really does suppress CRF rather than merely approximating it.
    #[test]
    fn inner_outer_supplement_is_lj_only_with_zero_crf() {
        use crate::zones::ZonePartition;

        let (topo, positions) = four_atom_fixture();
        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos = positions.clone();

        let region = AtomSelection::all(4);
        let periodicity = periodicity_of(&conf);
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

        // Inner = {0}, everything else outer: isolates a single inner-outer pair per atom 0.
        let partition = ZonePartition::new(4, &[0], &[]);
        let mut interaction =
            LjCrfInteraction::new(1.4, 1.0, 78.5, 0.0).with_zone_partition(partition);
        let partitioned = interaction
            .contribute(&region, &topo, &conf, &index)
            .unwrap();

        // Oracle: every pair evaluated with CRF forced to zero (q_prod = 0), since (2,3) is the
        // only pair not touching atom 0, and it's outer-outer (fully classical) — so to isolate
        // "is the (0,*) contribution LJ-only", zero just atom 0's charge and compare the total
        // against a full LJ+CRF run with the same charge zeroed for every pair.
        let charges = topo.charge.clone();
        let iac_u32: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
        let lj_mat = LJParamMatrix::from_nested(&[topo.lj_parameters[0]
            .iter()
            .map(LJParameters::from)
            .collect()]);
        let crf = CRFParameters::new(1.4, 1.0, 78.5, 0.0);

        // Full brute-force pairlist (all 6 pairs of 4 atoms) with atom 0's charge zeroed: this
        // makes every (0,*) pair LJ-only and leaves (1,2),(1,3),(2,3) as full classical LJ+CRF —
        // exactly what the partitioned provider should also produce, since inner={0} only ever
        // appears in inner-outer (LJ-only-supplement) pairs here (no inner-inner, no buffer).
        let mut charges_zero_0 = charges.clone();
        charges_zero_0[0] = 0.0;
        let all_pairs: Vec<(u32, u32)> = vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        let mut oracle = ForceStorage::new(4);
        lj_crf_innerloop_novirial(
            &positions,
            &charges_zero_0,
            &iac_u32,
            &all_pairs,
            &lj_mat,
            &crf,
            &periodicity,
            gromos_core::units::four_pi_eps_i,
            &mut oracle,
            EnergyGroups::single(),
        );

        assert_eq!(partitioned.energy, oracle.e_lj + oracle.e_crf);
    }
}
