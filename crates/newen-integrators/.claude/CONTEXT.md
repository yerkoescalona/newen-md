# newen-integrators — stage contract

## Job
L2 steppers. Mutating algorithms that advance state: integrators, constraints, thermostats,
barostats, energy minimisation.

## Inputs (consumes from)
`newen-core` + `newen-forces`.

## Outputs (public API)
`LeapFrog`, `SteepestDescent`, `SHAKE`/`SETTLE`/`LINCS`, Berendsen + Nosé-Hoover thermostats,
Berendsen barostat, COM removal.

## Status
- **Thermostat kernels verified (2026-10-22):** `berendsen_scale` and `nose_hoover_scale` are the
  pure per-bath steps `apply` calls, with one threshold `ZERO_TEMPERATURE = 1e-12` (gromosXX
  `math::epsilon`; was `1e-10` and `f64::EPSILON`). Lean model in `spec/NewenSpec/Thermostat.lean`
  (ℝ theorems in `spec/NewenTheory/Dynamics/Thermostats.lean`); `tests/thermostat_spec.rs`
  compares the Rust with the model's `Float` evaluation bit for bit on `spec/vectors/thermostat.json`.
  Keep the operation order of both kernels as it is: a reassociation fails that test.
- **Nosé–Hoover restart (0.0.69, D41):** `NoseHooverThermostat::init` restores ζ from
  `Configuration::nose_hoover` when `read_variables` (NTINHT = 0), every link of every bath,
  plain NH included (`restored_zeta`);
  `apply` mirrors ζ back into the configuration for the final NHCVARIABLES block.
- **MPI seam (2026-08-30):** `Forcefield::pair_partition = Some((rank, size))` keeps the pairs
  whose first atom index ≡ rank (mod size) after every pairlist update (and for the initial list);
  `set_nonbonded_reducer(f)` is called once per step with the pair-term `ForceStorage` (forces,
  E_LJ, E_CRF, virial, dH/dλ) — the cross-rank sum lives in `md`, not here. Terms every rank
  computes in full (excluded-pair/self RF, 1-4, non-pairlist perturbed corrections) accumulate in
  a separate `local_storage` and are added after the reduction. `Algorithm::as_any_mut` +
  `AlgorithmSequence::find_mut::<T>()` reach the built force field.
- SHAKE ✓ (solute NTC>1 + solvent NTCS>0, virial, skip optimisation, NTISHK)
- SETTLE ✓ (analytical 3-site rigid water, Miyamoto & Kollman 1992)
- LINCS ✓ (wired, reference-tested for solvent + solute)
- COM removal ✓: translational (NTICOM + NSCM) + rotational; uses min-image wrapped positions (gromosXX convention); periodic COM rotation suppressed in PBC
- Berendsen thermostat ✓ (MULTIBATH)
- Nosé-Hoover ✓: single NHC (algorithm=1) + chain NHC (algorithm=N≥2); both reference-tested
- Berendsen barostat ✓ (PRESSURESCALE), pressure/virial ✓
- Steepest Descent EM ✓ (adaptive step, FLIM, DELE, NMIN; works with SHAKE + posres)
- **Perturbed nonbonded** ✓: soft-core LJ+CRF correction (`perturbed_pairlist_correction`), self-energy, excluded-pair, and 1-4 corrections; wired via `Forcefield`; `ch4_water_fep` reference test passes.
- EDS ✓ (1.8): `eds::Eds`, the reference Hamiltonian after the force field (FORMS 1-3,
  reference-tested), with perturbed atoms (TI-EDS, 0.0.58) Σ p_i ∂V_i/∂λ into the derivative
  record's `eds_vr`
- A-EDS ✓ (1.8, 0.0.60): `Eds::accelerated`, V_R* and the fkfac-scaled EDS forces; FORMS 2-6
  search through the pure functions `memory_time`, `initialise_search`, `offset_update`,
  `track_state`, `emax_emin_update`, `advanced_offset_update` (`search_step`), on the
  configuration's `AedsSearch` (the `AEDSSEARCH` block); reference-tested and replayed against
  gromosXX with the D16, D25 and D32 patches (`newen-run/tests/aeds_replay.rs`). Since 0.0.63
  three gromosXX defects are corrected: the energy columns before the update (D16), log-sum-exp
  prevalences (D25), a sampled state as the E_min target (D32)
- GaMD ✓ (gromosXX's algorithm, `gamd::Gamd`, run at the end of `Forcefield::apply` when
  `Forcefield::gamd` is set by `md`): FORMS 1-3, SEARCH 0-2, the statistics and `GAMDSTAT`;
  one accelerated interaction group (plus the group gromosXX books the RF excluded solute pairs
  to), reference-tested. The older `GamdParameters`/`GamdRunner` serve only the `gamd` binary.
- REMD: **stubs only, untested**

## Key files
```
src/algorithms/       — Forcefield, LeapFrog*, Temperature, Energy
src/algorithms/berendsen_thermostat.rs, nosehoover_thermostat.rs — the thermostat kernels
                        (`berendsen_scale`, `nose_hoover_scale`) and their `Algorithm`s
src/algorithms/multibath.rs — atom → bath mapping, `ZERO_TEMPERATURE`
tests/thermostat_spec.rs — kernels = Lean `Float` model, bit for bit
src/constraints.rs    — SHAKE, SETTLE, LINCS
src/thermostats.rs    — Berendsen, Nosé-Hoover (+ chain NHC)
src/barostats.rs      — Berendsen, Parrinello-Rahman
```

## Decisions logged
- **Nosé-Hoover "small flexible constraints hack"** (`nosehoover_thermostat.cc:129,185`, also `berendsen_thermostat.cc:105`) **not ported** — only applies with flex-shake, which we don't support. Deliberate GROMOS quirk.
- **NTINHT convention:** 0 = read from file, 1 = generate from scratch (opposite of intuition — matches gromosXX, don't "fix" it).
- **IMD algorithm codes:** 0 = Berendsen, 1 = NHC, N = chain length (matches gromosXX after IMD parser bug fix).
- **FEP/TI:** second-source the perturbed RF self-term before porting (see newen-forces rules + porting.md).
- **Effective soft-core alpha = per-atom × global** (`in_perturbation.cc:1308`). The `.ptp` PERTATOMPARAM column ALJ/ACRF is a **per-atom scaling factor**, not the final alpha. gromosXX multiplies it by the global ALPHLJ/ALPHC from the `.imd` PERTURBATION block before storing: `lj_soft *= param.perturbation.soft_vdw`. `Forcefield.build_pert_info` must replicate this: `alpha_lj = pa.lj_soft * self.global_alphlj`. Failing to do so gives the wrong soft-core distance and off-by-2× errors in repulsive-regime pairs. Both `global_alphlj` and `global_alphc` must be set from `imd.alphlj`/`imd.alphc` before calling `init()`.
- **gromosXX echoes effective (post-multiplication) alphas** in its startup output under `PERTURBATION TOPOLOGY`. If the echo says `LJ(soft) 0.25` but the `.ptp` file has `ALJ 0.5`, that is normal — it shows `0.5 × ALPHLJ_global = 0.5 × 0.5 = 0.25`. Trust the echo for debugging, not the raw `.ptp` value.
- **Debugging FEP discrepancies: per-pair comparison is the right tool.** If total energies differ but the formula looks identical, add per-pair `eprintln!` for 5 representative pairs (one close/repulsive, one mid-range, one long-range) in both codes. If per-pair values match but totals differ, the pairlist counts or pairlist membership differ — check `GXX_OUTERLOOP call=N pairs=M` vs Rust diagnostic. If per-pair values also differ, the parameter (usually alpha) is wrong. In this case: start from gromosXX's echoed topology output, work backward from a single pair's energy to derive what parameter value gromosXX must be using, then trace where that value originates in the C++ source.
