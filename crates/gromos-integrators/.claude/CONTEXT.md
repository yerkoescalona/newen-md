# gromos-integrators — stage contract

## Job
L2 steppers. Mutating algorithms that advance state: integrators, constraints, thermostats,
barostats, energy minimisation.

## Inputs (consumes from)
`gromos-core` + `gromos-forces`.

## Outputs (public API)
`LeapFrog`, `SteepestDescent`, `SHAKE`/`SETTLE`/`LINCS`, Berendsen + Nosé-Hoover thermostats,
Berendsen barostat, COM removal.

## Status
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
  reference-tested); GaMD / REMD: **stubs only, untested**

## Key files
```
src/algorithms/       — Forcefield, LeapFrog*, Temperature, Energy
src/constraints.rs    — SHAKE, SETTLE, LINCS
src/thermostats.rs    — Berendsen, Nosé-Hoover (+ chain NHC)
src/barostats.rs      — Berendsen, Parrinello-Rahman
```

## Decisions logged
- **Nosé-Hoover "small flexible constraints hack"** (`nosehoover_thermostat.cc:129,185`, also `berendsen_thermostat.cc:105`) **not ported** — only applies with flex-shake, which we don't support. Deliberate GROMOS quirk.
- **NTINHT convention:** 0 = read from file, 1 = generate from scratch (opposite of intuition — matches gromosXX, don't "fix" it).
- **IMD algorithm codes:** 0 = Berendsen, 1 = NHC, N = chain length (matches gromosXX after IMD parser bug fix).
- **FEP/TI:** second-source the perturbed RF self-term before porting (see gromos-forces rules + porting.md).
- **Effective soft-core alpha = per-atom × global** (`in_perturbation.cc:1308`). The `.ptp` PERTATOMPARAM column ALJ/ACRF is a **per-atom scaling factor**, not the final alpha. gromosXX multiplies it by the global ALPHLJ/ALPHC from the `.imd` PERTURBATION block before storing: `lj_soft *= param.perturbation.soft_vdw`. `Forcefield.build_pert_info` must replicate this: `alpha_lj = pa.lj_soft * self.global_alphlj`. Failing to do so gives the wrong soft-core distance and off-by-2× errors in repulsive-regime pairs. Both `global_alphlj` and `global_alphc` must be set from `imd.alphlj`/`imd.alphc` before calling `init()`.
- **gromosXX echoes effective (post-multiplication) alphas** in its startup output under `PERTURBATION TOPOLOGY`. If the echo says `LJ(soft) 0.25` but the `.ptp` file has `ALJ 0.5`, that is normal — it shows `0.5 × ALPHLJ_global = 0.5 × 0.5 = 0.25`. Trust the echo for debugging, not the raw `.ptp` value.
- **Debugging FEP discrepancies: per-pair comparison is the right tool.** If total energies differ but the formula looks identical, add per-pair `eprintln!` for 5 representative pairs (one close/repulsive, one mid-range, one long-range) in both codes. If per-pair values match but totals differ, the pairlist counts or pairlist membership differ — check `GXX_OUTERLOOP call=N pairs=M` vs Rust diagnostic. If per-pair values also differ, the parameter (usually alpha) is wrong. In this case: start from gromosXX's echoed topology output, work backward from a single pair's energy to derive what parameter value gromosXX must be using, then trace where that value originates in the C++ source.
