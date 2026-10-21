# newen-forces — stage contract

## Job
L1 pure compute. Every force and energy calculation for classical terms.
Single source of physics — both MD engine and analysis facade call this layer.

## Inputs (consumes from)
`newen-core`: Topology, Configuration, BoundaryCondition, pairlist algorithms.

## Outputs (public API)
Bonded + nonbonded force/energy; pairlist; restraint forces; `single_point_energy()`;
`PotentialProvider` trait + implementors (QM/MM/ML seam, 2.6-P2.8).

## Status
- **Pair kernels verified (2026-10-22):** `lj_crf_pair` (one LJ + CRF pair at `r²`; `lj_crf_interaction`
  wraps it), `rf_excluded_pair`, `rf_self_term` and `CRFParameters::new` carry `/// spec:` keys.
  Lean model in `spec/NewenSpec/Nonbonded.lean` (ℝ theorems in
  `spec/NewenTheory/Interactions/`); `tests/nonbonded_spec.rs` compares the Rust with the model's
  `Float` evaluation bit for bit on `spec/vectors/nonbonded.json`, and a unit test in
  `innerloops.rs` checks every lane of `lj_crf_interaction_x4` against `lj_crf_pair`. Keep the
  operation order of these kernels as it is: a reassociation fails those tests.
- LJ + CRF nonbonded ✓; 1-4 interactions ✓; twin-range (RCUTP/RCUTL) ✓
- All bonded types ✓: quartic/harmonic bonds, cos-harmonic/harmonic angles, dihedrals, impropers, cross-dihedrals
- Perturbed nonbonded ✓: soft-core LJ+CRF, pairlist correction, self/excluded/1-4/PERTATOMPAIR corrections
- EDS end states ✓ (`nonbonded/eds.rs`, 1.8): per-state pairs, 1-4, RF self/excluded; with
  perturbed partners (TI-EDS, 0.0.58) the soft-core A/B pair per state and its ∂V_s/∂λ
- Position restraints ✓; distance restraints ✓
- `energy.rs` ✓ — `single_point_energy(topo, positions, box_dims, params)` — used by `ener` analysis binary
- ⚠️ Perturbed RF self-term needs second-sourcing from GROMOS book (`perturbed_nonbonded_term.cc:596,749,1444`)
- Dihedral restraints: **not yet** — next after CellList wiring (P1.6)
- **QM/MM + ML seam (2.6-P2.8), pieces built and reference-tested, not yet wired into an
  MD binary:**
  - `PotentialProvider` trait (`provider.rs`) — `contribute(region, topo, conf, neigh) ->
    Contribution` (scattered energy+forces+virial), plus `Embedding {None, Mechanical,
    Electrostatic}` declaring how a provider treats atoms outside its region.
  - `nonbonded::LjCrfInteraction` — classical LJ+CRF as a provider, wrapping
    `lj_crf_innerloop_novirial` with zero math changes. Solute-solute atom-level pairs only
    (solvent's charge-group pairlist shape isn't ported to this seam). Honours an optional
    `zone_partition: Option<ZonePartition>` (P2.8-1) — skips pairs a QM/ML provider or
    electrostatic embedding already owns; `None` (default) is today's unpartitioned behaviour.
    When set, also runs a second LJ-only pass (CRF zeroed via a per-call charge-array trick) for
    `ZonePartition::lj_only_should_evaluate` pairs — inner-outer pairs get classical LJ
    unconditionally in real gromosXX, only their electrostatics moves to embedding (P2.8-4).
  - `nonbonded::ElectrostaticEmbedding` — QM→MM Coulomb coupling, forces on **both** sides
    (Poliak et al. 2025 path (c): host computes pairwise Coulomb from QM-derived charges). Virial
    computed too (P2.8-3) — `trace(virial) == energy` exactly for a `1/r` pair (Euler's theorem),
    verified closed-form, FD, and on the real periodic `t_06` system (1363 pairs).
  - `nonbonded::SchNetInteraction` (feature `ml`, needs libtorch 2.11.0 — see its module docs for
    exact version pinning and known environment traps) — first ML provider, TorchScript via `tch`.
  - `zones.rs` — `Zone {Inner, Buffer, Outer}` + `ZonePartition` derives `PairOwner {Provider,
    Classical, Embedding}` per atom pair from the BuRNN training-data decomposition (not
    invented) — the double-counting contract an orchestrator must apply. `owner()` is the
    CRF/provider-energy table; `lj_owner()` is a *separate* LJ-specific table (sourced from real
    gromosXX, P2.8-4) — they diverge on inner-outer pairs (LJ always classical, CRF embedded).
    `qm_lj: bool` (default `false`) mirrors GROMOS `QMLJ`.
  - `newen_core::spatial_index::SpatialIndex` — query-based neighbor service independent of the
    MD pairlist's charge-group/twin-range shape, used by every provider above.
  - `orchestrator.rs` — `ProviderOrchestrator::register(region, provider)` +
    `evaluate(topo, conf, neigh)` sums `Vec<Box<dyn PotentialProvider>>` into one `Contribution`,
    enforcing that an `Embedding::None` provider never places a force outside its own `region`
    (hard error, not silently accumulated — the P2.6 review finding this exists to check).
  - `qmmm_orchestrator_reference.rs` — `LjCrfInteraction` (zone-partitioned) and
    `ElectrostaticEmbedding` registered together in one `ProviderOrchestrator`, on the real
    3513-atom `t_06` system; combined energy matches the two direct calls added by hand.
  - **Not done:** no *binary* drives `ProviderOrchestrator` yet. The QM/ML provider's own energy
    for the inner zone is still missing from every real-system test — confirmed, not assumed: this
    environment has no `xtb`/`mopac` binary, and `cargo build -p newen-forces --features ml`
    fails outright (no `libtorch` install), so `SchNetInteraction` cannot even be exercised here.
    QM charges refreshed per step needs a real QM engine too, not built speculatively.
  - `.local/gromosXX` (gitignored, `https://github.com/biomos/gromosXX`) is now cloned alongside
    `.local/gromos_tutorial_livecoms` — the real gromosXX C++ source is available locally for
    second-sourcing any future uncertain physics, not just the tutorial.
  - `tests/gromosXX_qmmm_references/water_dimer_mechst/` — a real (not gitignored, small,
    committed with attribution) gromosXX QM/MM run's input + reported energy, from the dataset
    accompanying Poliak et al. 2025 (Zenodo 10.5281/zenodo.14549978, CC BY 4.0).
    `water_dimer_qmmm_mechst_reference.rs` reproduces it to `1e-8`. Only the `NTQMMM=-1`
    (mechanical, constant-charge) slice is usable without an actual QM program (`mndo` — checked,
    not apt-installable here); everything else in that dataset needing dynamic QM charges stays
    out of reach for the same reason `xtb`/`libtorch` do.

- **3.9 step 5 (2026-08-29).** `ProviderOrchestrator::register_labelled` /
  `evaluate_with_terms` report every provider's energy under a label (the run's registry name for
  the term) into `Energy.term_energies` — the per-term view the summed `Contribution` erases (G10).
  `XtbInteraction::with_timeout` + `qm_subprocess::run_subprocess_with_timeout` bound every xtb
  call (child killed, clear error). The `newen-run` `xtb` term is built on these.

## Key files
```
src/bonded/          — bonded force calculations (split from bonded.rs)
src/nonbonded/       — LJ+CRF, rf_excluded, pairlist loops (split from nonbonded.rs)
  innerloops.rs        — `lj_crf_pair` (scalar pair kernel), `lj_crf_interaction_x4` (4 lanes), the loops
  excluded.rs          — `rf_excluded_pair`, `rf_self_term`, RF excluded and 1-4 loops
  params.rs            — `CRFParameters::new`
  interaction.rs      — LjCrfInteraction (PotentialProvider impl, zone-partition-aware)
  embedding.rs         — ElectrostaticEmbedding (QM/MM Coulomb, forces on MM atoms)
  schnet.rs            — SchNetInteraction (feature `ml`)
src/provider.rs       — PotentialProvider trait, Contribution, Embedding
src/orchestrator.rs    — ProviderOrchestrator: sums Vec<Box<dyn PotentialProvider>>, enforces the
                         Embedding::None index-scoping contract
src/zones.rs          — Zone, ZonePartition, PairOwner (QM/MM/BuRNN double-counting contract)
src/energy.rs         — single_point_energy() entry point for analysis
src/restraints.rs     — position + distance restraints
tests/nonbonded_spec.rs — pair kernels = Lean `Float` model, bit for bit
tests/                — QM/MM seam reference tests against the real BuRNN tutorial (t_06),
                         gitignored under .local/, skip gracefully if not cloned
```

## Crate-specific rules
- **Single source of physics.** The moment `ener.rs` or any analysis program re-implements LJ, you've recreated gromos++.
- **Bonded force vectors:** `v = pos(i) - pos(j)` (GROMOS convention).
- **Do NOT port the Martina solute/solvent misclassification bug** (`extended_grid_pairlist_algorithm.cc:1309`).
- **Perturbed pairlist approach:** newen uses one combined pairlist + correction (subtract state-A, add perturbed) — mathematically equivalent to separate perturbed pairlists when `alpha_lj` is correct.
- **`perturbed_lj_crf_interaction` LJ is cutoff-independent.** `crf.cutoff_sq` only affects CRF soft-core.
- Grep `interaction/` for `bug|fixme|wrong|hack` before porting any new subsystem.
- **`PotentialProvider` impls apply no topology exclusions.** `LjCrfInteraction` computes raw
  pairwise LJ+CRF for whatever pairs `SpatialIndex::neighbor_pairs` returns — it does not consult
  `Topology::exclusions`. Fine for its documented solute-solute scope; do not point it at a full
  solvated system's `AtomSelection::all` (bonded intramolecular pairs would be double-counted
  against the real bonded-force loop, though not against each other — see the P2.8-1 exit test
  for a correct restricted use).
- **A bare `Configuration::new()` defaults to `Box::vacuum()`.** Any test that builds a
  `Configuration` by hand and feeds it to a `PotentialProvider` must also set
  `conf.current_mut().box_config` — `LjCrfInteraction::contribute` derives its own periodicity
  from `conf`, not from whatever `Periodicity`/`SpatialIndex` the caller built separately, so a
  mismatch here silently uses unwrapped distances for cross-boundary pairs (caught building the
  P2.8-1 exit test: ~0.5% energy discrepancy from exactly this).

## The provider constraints P1–P9

The nine problems a QM/MM/ML provider has to solve, written before any of it was built and used as
the specification the seam was designed against. **They are cited by number from the code**
(`provider.rs`, `orchestrator.rs`, `spatial_index.rs`, the QM/ML reference tests) — a comment saying
"P5" means the constraint below, so keep the numbering stable. The design they describe is built,
so this is the record of *why* the shapes are what they are, not a plan.

- **(P1) The abstraction is data-in/data-out, not file-in/file-out.** The inverse of gromosXX, whose
  `QM_Worker` base class is shaped around `write_input`/`run_calculation`/`process_output` — the
  abstraction itself assumes "format a text file → fork/exec → parse text". Our trait takes arrays
  (positions, charges, region) and returns energy + gradients; a file-based external engine is *one
  implementation* that formats and parses internally. The file dance must never leak into the trait
  shape. This inversion is the whole ballgame.
- **(P2) The model-export boundary.** Most trained potentials live in PyTorch. The on-ramp is train
  in Python → export → load in Rust, but not every architecture exports cleanly: equivariant
  message-passing models (NequIP/Allegro/MACE, on `e3nn`) use custom ops and dynamic neighbour
  graphs that may not lower to ONNX. TorchScript via `tch` is more permissive and heavier. Which
  models are first-class follows from this choice.
- **(P3) Multiple, differently-parameterised neighbour structures.** The MD pairlist, an ML model's
  radial graph (different cutoff, edge lists, sometimes triplets) and the QM-zone gather (MM charges
  within cutoff of the QM region, rebuilt every step) are *three different spatial queries*. Hence
  `SpatialIndex` is a general query service, not an MD pairlist — see `newen-core/src/spatial_index.rs`.
- **(P4) Units and conventions, typed.** QM works in Hartree/Bohr, ML models in eV/Å or kcal/mol,
  GROMOS in kJ/mol·nm. gromosXX manages this with `print_unit_factors` and discipline. A unit
  mismatch must be a compile- or return-time error, not a silent 23×/27.2× force bug — a classic and
  brutal ML-MD failure mode.
- **(P5) Forces flow onto MM atoms too.** Electrostatic embedding means the QM density polarises and
  the MM point charges feel a force, so a provider returns gradients on **both** the QM region and
  the embedding MM atoms — contributions on an *arbitrary atom subset*, not "the QM region only".
  The accumulator interface accepts partial, scattered contributions, which is also what the
  classical providers want.
- **(P6) Stateful, iterative coupling.** Polarisable embedding is an SCF loop (mutual QM↔MM
  polarisation) and QM workers reuse the previous wavefunction as an initial guess, so the trait
  takes `&mut self` and cannot assume pure functions.
- **(P7) Uncertainty / committee / active learning.** gromosXX's `nn_worker` already loads validation
  models. Ensemble disagreement → uncertainty → trigger a DFT recompute is where ML-MD is going, so a
  provider may return uncertainty alongside energy and forces (`ProviderExtra`), leaving room for a
  scheduler to escalate ML → DFT at runtime.
- **(P8) Determinism boundary.** QM and ML paths are not bit-deterministic, so they sit in a
  *separate test tier* from the classical bit-for-bit suite: validated against reference energies and
  forces within a stated tolerance. Cited by every QM/ML reference test to say which tier it is in.
- **(P9) Process lifecycle and failure** for external engines: crashes, timeouts, NaNs. `Result` plus
  a supervised subprocess and a fallback.
