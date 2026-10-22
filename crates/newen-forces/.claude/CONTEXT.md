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
- **Bonded kernels verified (2026-10-22):** `quartic_bond`, `harmonic_bond`, `cosine_angle`,
  `harmonic_angle`, `dihedral_chebyshev`, `dihedral_new`, `improper_dihedral` (one pure kernel per
  form in `bonded/`, called by the `calculate_*_grouped` loops) carry `/// spec:` keys. Lean model
  in `spec/NewenSpec/Bonded.lean`, theory in `spec/NewenTheory/Interactions/Bonded.lean`;
  `tests/bonded_spec.rs` compares them bit for bit on `spec/vectors/bonded.json` (acos/cos/sin from
  the same C library as Lean's `Float`) and checks every force against the energy by finite
  differences. The guards are `ZERO_DISTANCE`, `ZERO_DISTANCE2`, `ZERO_SINE`, `NEAR_LINEAR` in
  `bonded/mod.rs`. Same rule: keep the operation order.
- **Perturbed kernels verified (2026-10-22):** the free-energy loops in `nonbonded/perturbed.rs`
  and `bonded/perturbed.rs` call pure scalar kernels with `/// spec:` keys:
  `PerturbedLambdaParams::from_lambda` (the λ table, `powf`), `perturbed_lj_crf_pair` (soft-core
  LJ + CRF at `r²`; `perturbed_lj_crf_interaction[_scaled]` wrap it), `rf_soft_interaction`,
  `perturbed_rf_excluded_pair`, `perturbed_rf_self_term`, `lj_crf_soft_interaction_ext`,
  `rf_soft_interaction_ext`; `chebyshev`, `perturbed_quartic_bond_term`, `perturbed_angle_term`,
  `perturbed_improper_term`, `perturbed_dihedral_term`, `soft_perturbed_term`. Lean models in
  `spec/NewenSpec/SoftCore.lean` and `PerturbedBonded.lean`; over ℝ each kernel's dV/dλ is proven
  to be the derivative of its energy (`softcore_dvdl_is_deriv`,
  `perturbedLjCrfPairModel_dvdl_is_deriv`, `pert_*_dvdl_is_deriv`), and the end states λ = 0, 1 are
  the plain terms. `tests/softcore_spec.rs` and `tests/perturbed_bonded_spec.rs` compare the Rust
  with the `Float` models bit for bit and check dV/dλ and the forces against central differences.
  Keep the operation order of these kernels as it is. The extraction moved no bit: all 121
  gromosXX reference systems give byte-identical output files before and after. Since
  2026-10-22 `perturbed_lj_crf_pair` has no distance guard and takes `C12/C6` for any `C6 ≠ 0`,
  as gromosXX and the ext kernel (overlapping soft-core atoms are finite; a state with `C6 = 0`
  has no soft-core length there and diverges at `r = 0`, as in gromosXX).
- **Lattice-sum kernels verified (2026-10-22):** `nonbonded/lattice_sum.rs` holds the
  pieces Ewald (NLRELE 2) and P3M (NLRELE 3) share, ported from gromosXX in its operation order:
  `ChargeShape` (NSHAPE −1 Gaussian, 0..=10 truncated polynomials), `shape_switch` (`η`, `η'`),
  `shape_fourier` (`γ̂`, `γ̂'`), `lj_ls_pair` (LJ + real-space pair, gromosXX's LJ force form),
  `ls_excluded_pair`, `self_constants` (`A1`, `A3`), `reciprocal_box`, `lattice_vector`,
  `a2_cubic`, `a2_shells` (bounded by `max_shells`), `charge_sums`, `self_and_a_terms`,
  `self_and_a_virials`, `surface_term`, with `/// spec:` keys. C `erfc` is declared `extern "C"`
  there (Rust's std has none). Lean model `spec/NewenSpec/LatticeSum.lean`, theory
  `spec/NewenTheory/Interactions/LatticeSum.lean`; `tests/lattice_sum_spec.rs` compares the phase-1
  kernels bit for bit on `spec/vectors/lattice_sum.json`, checks `η'`, `γ̂'` and the pair forces by
  finite differences and `γ̂`, `A1`, `A3` against Vol. 2's shape table by quadrature. Keep the
  operation order.
- **Direct Ewald sum (NLRELE 2), wired (2026-10-22):** `nonbonded/ewald.rs`: `k_vectors` (`|l_j| ≤
  NK_j`, `k ≠ 0`, `|k|² ≤ KCUT` as gromosXX compares), `ewald_kspace` (energy, forces, `Ã2`,
  virial `Ξ`), and `EwaldSum` (the per-run state: rebuilds the wave vectors and `A1`, `A2`, `A3`
  when the box changes, then adds the k-space, self, A and surface terms to a `ForceStorage`, the
  virial as `−2 Ξ`). `nonbonded/ls_real_space.rs`: the real-space loops (pairlist pairs, excluded
  pairs, 1-4 pairs, LJ-only long-range shell), energies to the storage's `e_crf`. `ewald_spec.rs`
  compares the kernels with the Lean model bit for bit on `spec/vectors/ewald.json` and checks
  the forces and virials by finite differences, translation invariance, and the triclinic
  periodicity. Differences from gromosXX (see `docs/src/verification/index.md`, phase 2): the
  triclinic wave vectors (gromosXX multiplies by the transpose of the reciprocal basis), the
  surface virial (gromosXX adds `+Ξ`'s diagonal, not `−2 Ξ`). `1/ε₀` is a parameter: a run passes
  `FPEPSI · 4 · π` (gromosXX `in_topology.cc:1408`); `newen_core::units::EPS0_INV` is gromosXX's
  literal `eps0_i = 1745.9137`, used only without a PHYSICALCONSTANTS block.
- **P3M (NLRELE 3), wired (2026-10-22):** `nonbonded/p3m.rs`: the per-point kernels with
  `/// spec:` keys (`assignment_1d`, `nearest_grid_point`, `assignment_weights`, `Mesh::new`,
  `assign_charges`, `selfterm_fp`, `assign_squared_charges`, `fd_coefficients`, `sinc_power`
  (C `pow`, `extern "C"`), `influence_at` with `influence_function`, `quality_at`,
  `influence_box_correction`, `force_error_estimate`, `kspace_energy`, `potential_hat`,
  `field_ik`, `field_fd`, `interpolate_forces`, `a2_tilde_p3m`), `Fft3` (PhastFT along z, y, x,
  both directions unnormalised; powers of two only, `Fft3::new` panics otherwise, the run checks
  refuse first), and `P3m` (per-run state: the influence function at the first evaluation,
  refused above FACCUR, re-evaluated every NQEVAL steps, corrected to first order in a changed
  box). The self, A and surface terms are shared with Ewald (`ewald.rs::add_constant_terms`,
  `self_constants_of_box`), and `ewald.rs::LatticeSum` is the enum the force field drives.
  `tests/p3m_spec.rs`: the kernels bit for bit against `spec/vectors/p3m.json`, `Fft3` within
  1e-12 of the Lean direct transform (measured 3e-15 on 16³), the k-space energy within 1e-11 of
  the Lean pipeline, P3M against the direct Ewald sum within the force-error estimate, forces
  against central differences of the energy. Differences from gromosXX: see
  `docs/src/verification/index.md`, phase 3 (the NVT A term, the triclinic finite-difference
  field, the triclinic wave vectors). The FFT is PhastFT only: the unwired SPME code
  (`electrostatics.rs`, `pme.rs`, `pme_mpi.rs`, not GROMOS's P3M), its `use-fftw` feature and the
  `rustfft`/`fftw` dependencies were removed (2026-10-22).
- **Lattice sums against the patched oracle (2026-10-22):** the corrected gromosXX defects are
  pinned by `_fixed` references: LS1 (NVT A term), LS2 (surface virial,
  `water_216_npt_ewald_epsls_fixed`), LS3 (triclinic wave vectors,
  `nacl_water_box_triclinic_p3m_fixed`, ik differentiation). The triclinic finite-difference
  field (NFDORD > 0) is corrected in newen but has no patched oracle yet.
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
  lattice_sum.rs       — charge shapes, `lj_ls_pair`, `ls_excluded_pair`, `self_constants`, `A2`, self/A/surface terms
  ewald.rs             — `k_vectors`, `ewald_kspace`, `EwaldSum` (NLRELE 2), the shared constant terms, `LatticeSum`
  p3m.rs               — P3M (NLRELE 3): assignment, influence function, `Fft3` (PhastFT), fields, `P3m`
  ls_real_space.rs     — the real-space loops of a lattice-sum run
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
tests/lattice_sum_spec.rs — lattice-sum kernels = Lean `Float` model, bit for bit; FD, quadrature
tests/ewald_spec.rs   — Ewald kernels = Lean `Float` model, bit for bit; forces and virials by FD
tests/p3m_spec.rs     — P3M kernels = Lean `Float` model, bit for bit; FFT within a tolerance; vs Ewald
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
