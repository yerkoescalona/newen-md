//! Zone partitioning and **interaction ownership** — who computes which atom pair.
//!
//! ** P2.7 Step 4.** In a partitioned QM/MM or BuRNN system, more than one provider can
//! plausibly claim the same atom pair, and nothing detects that at runtime: two independently
//! written providers share no contract (assumption **A5**). Double-counting is therefore
//! prevented by an *explicit, static* ownership table applied by the orchestrator — not by
//! hidden logic inside either provider, and not by hoping their conventions happen to agree.
//!
//! # Why this exact table (derived, not assumed)
//!
//! The BuRNN training target is read directly from the tutorial's own data-generation code,
//! `train_dataset_tutorial/mopac.py::get_burnn_energy()`:
//!
//! ```text
//! E_burnn = [E_QM(inner+buffer) − N_h2o·E_h2o_vac − E_solute_vac]
//!         − [E_QM(buffer)       − N_h2o·E_h2o_vac]
//! ```
//!
//! Both terms come from MOPAC — this is **QM minus QM**, not QM minus MM (assumption **A2**;
//! the vacuum terms are per-species self-energy normalisation). With the same waters in both
//! calls this reduces to `E_QM(inner+buffer) − E_QM(buffer) − E_solute_vac`, so the network
//! carries exactly: the inner region's internal energy, the full inner↔buffer interaction, and
//! the buffer's polarisation response — and explicitly **not** the buffer's own internal energy,
//! which the subtraction removes. `get_burnn_forces()` agrees per-atom: inner atoms take the
//! complex force unchanged, buffer atoms take `f_complex − f_buffer`.
//!
//! Everything the network does not carry must therefore still be computed classically, which
//! fixes the ownership of all six pair classes ([`PairOwner`]).
//!
//! # Van der Waals across the QM boundary (P2.8-4, resolved)
//!
//! Originally left undecided here — the tutorial's `md.imd` `QMMM` block lists seven field names
//! against six values, so `QMLJ`'s position was ambiguous from the tutorial alone. Resolved by
//! cloning `gromosXX` (`.local/gromosXX/md++/src`, gitignored) and reading the real source:
//! `io/parameter/in_parameter.cc::read_QMMM` (field order `NTQMMM NTQMSW RCUTQM NTWQMMM QMLJ
//! QMCON MMSCALE`) and `interaction/qmmm/qmmm_interaction.cc::modify_exclusions`, whose own
//! comment states the rule plainly:
//!
//! ```text
//!      QM  Buf  MM
//! QM    R    R  R/K
//! Buf   R    K   K
//! MM    R    K   K
//! R/K - keep in mechanical embedding, remove otherwise
//! Copy QM-QM only if requested to do QM-QM LJ
//! Copy QM-MM only if electrostatic or polarizable embedding
//! ```
//!
//! `QMLJ` (`qm_lj` on [`ZonePartition`]) gates the `copy` flag **only** for QM-QM (inner-inner)
//! and QM-buffer (inner-buffer) pairs (`case QM-QM`/the buffer branch of `case QM-MM` in
//! `modify_exclusions`) — i.e. it decides whether *provider-owned* pairs also get a classical LJ
//! supplement. It does **not** gate inner-outer pairs: those go through the unconditional
//! `else if (qmmm > mechanical) { copy = erase = true; }` branch, so **inner-outer LJ is always
//! classical**, regardless of `QMLJ` or the embedding scheme — only inner-outer *electrostatics*
//! is replaced (by [`Embedding::Electrostatic`]). `t_06` sets `QMLJ=0`, matching this crate's
//! previous (undocumented) assumption for inner-inner/inner-buffer, but that assumption was never
//! correct for inner-outer — see [`ZonePartition::lj_owner`], which is now the gromosXX-faithful
//! answer instead of a guess.
//!
//! # What this module still deliberately does not decide
//!
//! - **Excluded/1-4 pairs**, which follow the usual GROMOS rules independently of zoning.

use crate::provider::Embedding;
use newen_core::selection::{AtomSelection, SelectionError};
use newen_core::topology::Topology;

/// Which region of a partitioned system an atom belongs to.
///
/// Named after the BuRNN scheme (Gómez-Flores et al. 2022, *J. Chem. Theory Comput.* 18:1213).
/// A plain QM/MM system uses only [`Zone::Inner`] and [`Zone::Outer`] — no buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Zone {
    /// Treated fully by the QM/ML provider.
    Inner,
    /// Dual-described: the provider sees it, *and* it keeps its classical description.
    Buffer,
    /// Classical MM only.
    #[default]
    Outer,
}

/// Which computation is responsible for a given atom pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairOwner {
    /// The QM/ML provider. Classical providers **must skip** these pairs.
    Provider,
    /// The classical force field, using topology charges as usual.
    Classical,
    /// Classical Coulomb, but with QM-derived charges — [`Embedding::Electrostatic`]
    /// (see [`crate::nonbonded::ElectrostaticEmbedding`]).
    Embedding,
}

/// Per-atom zone assignment plus the ownership rule derived from it.
#[derive(Debug, Clone)]
pub struct ZonePartition {
    zones: Vec<Zone>,
    /// GROMOS `QMLJ` (`io/parameter/in_parameter.cc::read_QMMM`). Gates whether inner-inner and
    /// inner-buffer pairs get a classical LJ supplement on top of the provider's own energy.
    /// Defaults `false` — matches every reference system in this codebase (`t_06` sets `QMLJ=0`).
    /// Does **not** affect inner-outer pairs; see [`Self::lj_owner`].
    qm_lj: bool,
}

impl ZonePartition {
    /// All atoms classical; assign zones with [`Self::set_zone`].
    pub fn all_outer(n_atoms: usize) -> Self {
        Self {
            zones: vec![Zone::Outer; n_atoms],
            qm_lj: false,
        }
    }

    /// Build from explicit inner/buffer index lists; everything else is [`Zone::Outer`].
    pub fn new(n_atoms: usize, inner: &[usize], buffer: &[usize]) -> Self {
        let mut p = Self::all_outer(n_atoms);
        for &i in inner {
            p.zones[i] = Zone::Inner;
        }
        for &i in buffer {
            p.zones[i] = Zone::Buffer;
        }
        p
    }

    /// Build from selector strings (the `AtomSelection::from_string` grammar — `1:res(LIG:a)`,
    /// name lists, `not()`/`minus()`) instead of hand-counted index lists — the "renumber the
    /// ligand every time the topology changes" ergonomics problem `AtomSelection::from_string`
    /// already solved for atom selection generally, applied here to zones (P3.7).
    /// `buffer = None` gives a plain inner/outer QM/MM split, no buffer zone.
    pub fn from_selections(
        topo: &Topology,
        inner: &str,
        buffer: Option<&str>,
    ) -> Result<Self, SelectionError> {
        let n_atoms = topo.num_atoms();
        let inner_indices = AtomSelection::from_string(inner, topo)?.indices().to_vec();
        let buffer_indices = match buffer {
            Some(spec) => AtomSelection::from_string(spec, topo)?.indices().to_vec(),
            None => Vec::new(),
        };
        Ok(Self::new(n_atoms, &inner_indices, &buffer_indices))
    }

    /// Set GROMOS `QMLJ` — `true` adds a classical LJ supplement for inner-inner/inner-buffer
    /// pairs on top of the provider's own energy (see module docs). Does not affect inner-outer
    /// pairs, which are always classical for LJ regardless of this flag.
    pub fn with_qm_lj(mut self, qm_lj: bool) -> Self {
        self.qm_lj = qm_lj;
        self
    }

    pub fn set_zone(&mut self, atom: usize, zone: Zone) {
        self.zones[atom] = zone;
    }

    pub fn zone(&self, atom: usize) -> Zone {
        self.zones[atom]
    }

    pub fn len(&self) -> usize {
        self.zones.len()
    }

    pub fn is_empty(&self) -> bool {
        self.zones.is_empty()
    }

    /// Atoms in a given zone, ascending.
    pub fn atoms_in(&self, zone: Zone) -> Vec<usize> {
        (0..self.zones.len())
            .filter(|&i| self.zones[i] == zone)
            .collect()
    }

    /// **The CRF/electrostatics and provider-energy ownership table.** Symmetric in `i`/`j`.
    /// This is *not* LJ ownership — gromosXX splits the two (see module docs and
    /// [`Self::lj_owner`]); this table governs the provider's own energy and, for inner-outer
    /// pairs, which side owns the *electrostatics* specifically.
    ///
    /// | pair | owner | why |
    /// |---|---|---|
    /// | inner–inner | [`PairOwner::Provider`] | in `E_QM(inner+buffer)`, survives the subtraction |
    /// | inner–buffer | [`PairOwner::Provider`] | the inner↔buffer interaction the NN is trained on |
    /// | buffer–buffer | [`PairOwner::Classical`] | removed by `− E_QM(buffer)`, so MM must supply it |
    /// | inner–outer | [`PairOwner::Embedding`] | outside the NN's input; classical with QM charges |
    /// | buffer–outer | [`PairOwner::Classical`] | never seen by the NN |
    /// | outer–outer | [`PairOwner::Classical`] | pure MM |
    pub fn owner(&self, i: usize, j: usize) -> PairOwner {
        use Zone::*;
        match (self.zones[i], self.zones[j]) {
            (Inner, Inner) | (Inner, Buffer) | (Buffer, Inner) => PairOwner::Provider,
            (Inner, Outer) | (Outer, Inner) => PairOwner::Embedding,
            (Buffer, Buffer) | (Buffer, Outer) | (Outer, Buffer) | (Outer, Outer) => {
                PairOwner::Classical
            },
        }
    }

    /// Should a classical provider evaluate this pair?
    ///
    /// The one call an orchestrator needs when filtering a pairlist before handing it to
    /// [`crate::nonbonded::LjCrfInteraction`]: `true` for anything the classical force field
    /// still owns, `false` for pairs the QM/ML provider has already accounted for.
    pub fn classical_should_evaluate(&self, i: usize, j: usize) -> bool {
        self.owner(i, j) == PairOwner::Classical
    }

    /// Should the given embedding scheme evaluate this pair?
    ///
    /// [`Embedding::Electrostatic`] owns exactly the inner↔outer pairs; anything else owns none
    /// (an [`Embedding::None`] provider by definition ignores its environment).
    pub fn embedding_should_evaluate(&self, embedding: Embedding, i: usize, j: usize) -> bool {
        embedding == Embedding::Electrostatic && self.owner(i, j) == PairOwner::Embedding
    }

    /// **The LJ ownership table (P2.8-4).** Distinct from [`Self::owner`] — sourced from
    /// gromosXX's real `QMMM_Interaction::modify_exclusions` (see module docs), not derived from
    /// the BuRNN energy decomposition the way [`Self::owner`] is (the NN's *own* training target
    /// says nothing about the classical field's LJ bookkeeping).
    ///
    /// | pair | `owner()` | `lj_owner()` | why they differ |
    /// |---|---|---|---|
    /// | inner–outer | `Embedding` | **`Classical`** | LJ is always classical here; only the electrostatics moves to embedding |
    /// | inner–inner, inner–buffer | `Provider` | `Classical` iff `qm_lj`, else `Provider` | `QMLJ` gates a classical LJ *supplement* on top of the provider's own energy |
    /// | buffer–buffer, buffer–outer, outer–outer | `Classical` | `Classical` | unchanged, gromosXX never touches these pairs' exclusions |
    pub fn lj_owner(&self, i: usize, j: usize) -> PairOwner {
        use Zone::*;
        match (self.zones[i], self.zones[j]) {
            (Inner, Outer) | (Outer, Inner) => PairOwner::Classical,
            (Inner, Inner) | (Inner, Buffer) | (Buffer, Inner) => {
                if self.qm_lj {
                    PairOwner::Classical
                } else {
                    PairOwner::Provider
                }
            },
            (Buffer, Buffer) | (Buffer, Outer) | (Outer, Buffer) | (Outer, Outer) => {
                PairOwner::Classical
            },
        }
    }

    /// Pairs the classical field must give an LJ-only supplement to: `lj_owner() == Classical`
    /// (LJ is classical) but `owner() != Classical` (the pair was skipped by
    /// [`Self::classical_should_evaluate`], so its combined LJ+CRF pass never ran). This is
    /// exactly the gap [`crate::nonbonded::LjCrfInteraction`]'s zone-partitioned pass must close
    /// with CRF zeroed — inner-outer pairs always qualify; inner-inner/inner-buffer pairs qualify
    /// only when `qm_lj` is set.
    pub fn lj_only_should_evaluate(&self, i: usize, j: usize) -> bool {
        self.lj_owner(i, j) == PairOwner::Classical && self.owner(i, j) != PairOwner::Classical
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use newen_core::topology::Atom;

    /// 6 atoms: 0-1 inner, 2-3 buffer, 4-5 outer.
    fn partition() -> ZonePartition {
        ZonePartition::new(6, &[0, 1], &[2, 3])
    }

    /// 4 atoms, 2 residues: LIG (0-1), SOL (2-3) — enough to exercise
    /// `AtomSelection::from_string`'s residue-name grammar.
    fn ligand_solvent_topo() -> Topology {
        let atoms: &[(&str, usize, &str)] = &[
            ("C1", 1, "LIG"),
            ("C2", 1, "LIG"),
            ("OW", 2, "SOL"),
            ("HW", 2, "SOL"),
        ];
        let mut topo = Topology::new();
        for &(name, res_nr, res_name) in atoms {
            topo.moltypes[0].atoms.push(Atom {
                name: name.into(),
                residue_nr: res_nr,
                residue_name: res_name.into(),
                iac: 0,
                mass: 12.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
            topo.iac.push(0);
            topo.mass.push(12.0);
            topo.charge.push(0.0);
        }
        // Populate molecules[0] so molecule_nr()/residue_name() work — same as
        // newen_core::selection's own `aladip_topo()` test fixture.
        topo.init_solute_moltype();
        topo
    }

    #[test]
    fn from_selections_matches_hand_built_partition() {
        let topo = ligand_solvent_topo();
        let by_name = ZonePartition::from_selections(&topo, "1:res(LIG:a)", None)
            .expect("selector should resolve");
        let by_index = ZonePartition::new(4, &[0, 1], &[]);
        for atom in 0..4 {
            assert_eq!(
                by_name.zone(atom),
                by_index.zone(atom),
                "atom {atom}: selector-built partition disagrees with the hand-built one"
            );
        }
    }

    #[test]
    fn from_selections_with_buffer() {
        let topo = ligand_solvent_topo();
        let p = ZonePartition::from_selections(&topo, "1:res(LIG:a)", Some("1:res(SOL:a)"))
            .expect("selectors should resolve");
        assert_eq!(p.zone(0), Zone::Inner);
        assert_eq!(p.zone(1), Zone::Inner);
        assert_eq!(p.zone(2), Zone::Buffer);
        assert_eq!(p.zone(3), Zone::Buffer);
    }

    #[test]
    fn from_selections_propagates_a_bad_selector_as_an_error() {
        let topo = ligand_solvent_topo();
        let err = ZonePartition::from_selections(&topo, "not a valid selector!!", None)
            .expect_err("a nonsense selector must not silently resolve to an empty/wrong zone");
        assert!(matches!(err, SelectionError::ParseError(_)));
    }

    /// Step 4's exit criterion: every one of the six pair classes matches the decomposition
    /// derived from `get_burnn_energy()` in A2, term by term.
    #[test]
    fn ownership_matches_burnn_energy_decomposition() {
        let p = partition();

        // In the NN: survives `E_QM(inner+buffer) − E_QM(buffer)`.
        assert_eq!(p.owner(0, 1), PairOwner::Provider, "inner-inner");
        assert_eq!(p.owner(0, 2), PairOwner::Provider, "inner-buffer");

        // Removed by the `− E_QM(buffer)` term, so the classical field must still supply it.
        // This is the single most double-counting-prone class: it is inside the NN's input
        // geometry yet deliberately *not* in its energy.
        assert_eq!(p.owner(2, 3), PairOwner::Classical, "buffer-buffer");

        // Never in the NN's input at all.
        assert_eq!(p.owner(2, 4), PairOwner::Classical, "buffer-outer");
        assert_eq!(p.owner(4, 5), PairOwner::Classical, "outer-outer");

        // Outside the NN, but the inner zone carries QM-derived charges.
        assert_eq!(p.owner(0, 4), PairOwner::Embedding, "inner-outer");
    }

    #[test]
    fn ownership_is_symmetric() {
        let p = partition();
        for i in 0..6 {
            for j in 0..6 {
                assert_eq!(p.owner(i, j), p.owner(j, i), "pair ({i},{j}) asymmetric");
            }
        }
    }

    /// The orchestrator's filter must exclude precisely what the provider already counted —
    /// no more (would drop real physics), no less (would double-count).
    #[test]
    fn classical_filter_excludes_exactly_the_provider_owned_pairs() {
        let p = partition();
        let excluded: Vec<(usize, usize)> = (0..6)
            .flat_map(|i| (i + 1..6).map(move |j| (i, j)))
            .filter(|&(i, j)| !p.classical_should_evaluate(i, j))
            .collect();

        // inner-inner (0,1); inner-buffer (0,2),(0,3),(1,2),(1,3); inner-outer (0,4),(0,5),(1,4),(1,5)
        assert_eq!(
            excluded,
            vec![
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (0, 5),
                (1, 2),
                (1, 3),
                (1, 4),
                (1, 5)
            ]
        );
        // Nothing involving only buffer/outer may ever be excluded from the classical field.
        assert!(excluded.iter().all(|&(i, _)| p.zone(i) == Zone::Inner));
    }

    #[test]
    fn embedding_owns_only_inner_outer_and_only_when_electrostatic() {
        let p = partition();
        assert!(p.embedding_should_evaluate(Embedding::Electrostatic, 0, 4));
        assert!(!p.embedding_should_evaluate(Embedding::Electrostatic, 0, 2)); // inner-buffer: NN
        assert!(!p.embedding_should_evaluate(Embedding::Electrostatic, 2, 4)); // buffer-outer: MM
        assert!(!p.embedding_should_evaluate(Embedding::None, 0, 4));
        assert!(!p.embedding_should_evaluate(Embedding::Mechanical, 0, 4));
    }

    /// A plain QM/MM system (no buffer) must degrade to the classic two-zone picture with no
    /// special cases: QM internal to the provider, QM↔MM embedded, MM internal classical.
    #[test]
    fn degrades_to_plain_qm_mm_without_a_buffer() {
        let p = ZonePartition::new(4, &[0, 1], &[]);
        assert_eq!(p.owner(0, 1), PairOwner::Provider);
        assert_eq!(p.owner(0, 2), PairOwner::Embedding);
        assert_eq!(p.owner(2, 3), PairOwner::Classical);
    }

    #[test]
    fn all_outer_is_fully_classical() {
        let p = ZonePartition::all_outer(4);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(p.owner(i, j), PairOwner::Classical);
            }
        }
        assert_eq!(p.atoms_in(Zone::Outer).len(), 4);
    }

    /// ** P2.8-4's central finding, as a test.** `owner` and `lj_owner` must diverge
    /// on inner-outer pairs specifically: `owner()` says `Embedding` (electrostatics moves to the
    /// embedding term) but `lj_owner()` says `Classical` (LJ stays classical unconditionally,
    /// per gromosXX's real `modify_exclusions` — see module docs). Getting this backwards is
    /// exactly the "wrong LJ ownership silently changes energies" failure mode the plan warned
    /// against before the source was read.
    #[test]
    fn inner_outer_lj_is_classical_even_though_electrostatics_is_embedded() {
        let p = partition();
        assert_eq!(
            p.owner(0, 4),
            PairOwner::Embedding,
            "electrostatics: embedded"
        );
        assert_eq!(
            p.lj_owner(0, 4),
            PairOwner::Classical,
            "LJ: always classical here"
        );
        assert!(
            p.lj_only_should_evaluate(0, 4),
            "inner-outer needs the LJ-only supplement, since classical_should_evaluate skips it"
        );
    }

    /// `QMLJ=0` (this crate's default, matching every reference system) leaves inner-inner and
    /// inner-buffer LJ fully owned by the provider — no classical supplement, no double-counting.
    #[test]
    fn qm_lj_off_leaves_inner_pairs_fully_provider_owned_for_lj() {
        let p = partition();
        assert_eq!(p.lj_owner(0, 1), PairOwner::Provider, "inner-inner");
        assert_eq!(p.lj_owner(0, 2), PairOwner::Provider, "inner-buffer");
        assert!(!p.lj_only_should_evaluate(0, 1));
        assert!(!p.lj_only_should_evaluate(0, 2));
    }

    /// `QMLJ=1` adds a classical LJ supplement on top of the provider's own energy for exactly
    /// the two cases gromosXX's `copy = sim.param().qmmm.qm_lj` covers — inner-inner and
    /// inner-buffer — without touching inner-outer (already unconditionally classical) or
    /// anything buffer/outer (never gated by `QMLJ` at all).
    #[test]
    fn qm_lj_on_adds_classical_lj_supplement_only_where_gromosxx_does() {
        let p = ZonePartition::new(6, &[0, 1], &[2, 3]).with_qm_lj(true);
        assert_eq!(p.lj_owner(0, 1), PairOwner::Classical, "inner-inner");
        assert_eq!(p.lj_owner(0, 2), PairOwner::Classical, "inner-buffer");
        assert!(p.lj_only_should_evaluate(0, 1));
        assert!(p.lj_only_should_evaluate(0, 2));

        // Unaffected by qm_lj either way.
        assert_eq!(
            p.lj_owner(0, 4),
            PairOwner::Classical,
            "inner-outer: always classical"
        );
        assert!(
            p.lj_only_should_evaluate(0, 4),
            "inner-outer still needs the supplement regardless of qm_lj \u{2014} owner() stays Embedding"
        );
        assert_eq!(
            p.lj_owner(2, 3),
            PairOwner::Classical,
            "buffer-buffer: untouched by qm_lj"
        );
        assert!(
            !p.lj_only_should_evaluate(2, 3),
            "already classical via owner(), no supplement needed"
        );
    }

    /// `lj_only_should_evaluate` must never overlap with `classical_should_evaluate` — each pair
    /// gets exactly one classical LJ+CRF path (combined) or LJ-only supplement, never both,
    /// mirroring the no-double-counting guarantee `classical_filter_excludes_exactly_the_provider_owned_pairs`
    /// already proves for the combined case.
    #[test]
    fn lj_only_supplement_never_overlaps_the_combined_classical_pass() {
        for qm_lj in [false, true] {
            let p = ZonePartition::new(6, &[0, 1], &[2, 3]).with_qm_lj(qm_lj);
            for i in 0..6 {
                for j in 0..6 {
                    if i == j {
                        continue;
                    }
                    assert!(
                        !(p.classical_should_evaluate(i, j) && p.lj_only_should_evaluate(i, j)),
                        "pair ({i},{j}) claimed by both the combined pass and the LJ-only supplement (qm_lj={qm_lj})"
                    );
                }
            }
        }
    }
}
