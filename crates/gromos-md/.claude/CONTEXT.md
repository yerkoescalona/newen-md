# gromos-md — stage contract

## Job
L3 orchestration. The `md` binary: composes gromos-core + gromos-forces + gromos-integrators +
gromos-io into a full MD simulation loop matching gromosXX output bit-for-bit.

## Inputs (consumes from)
All lower crates (gromos-core, gromos-forces, gromos-integrators, gromos-io).

## Outputs
`md` binary: runs GROMOS MD simulations; writes trajectories, energies, force trajectories
matching gromosXX to reference tolerances.

## Status
- **Net charge (1.8, 0.0.59):** `md` prints `Perturbation net charge: A … B … (ΔQ = …)` on a
  perturbed run and warns when ΔQ ≠ 0 under a reaction field (`fep_null_perturbation` pins
  ΔQ = 0 and no warning for the null perturbation).
- **MPI (2026-08-30):** `md` built with `--features use-mpi` (`LIBCLANG_PATH=/usr/lib/llvm-19/lib`)
  runs under `mpirun -np N`: every rank runs the whole driver on the same inputs, the pair terms are
  split by `RunOptions::pair_partition` and summed each step (`mod mpi_support` at the bottom of
  `md.rs`: reduce to rank 0 + broadcast, 3N+13 doubles), rank 0 alone writes files and stdout.
  `np=1` is bit-identical to the serial binary. The old `md_mpi`/`md_mpi_cuda`/`repex_mpi` bins are
  gone. Timing: `scripts/bench_mpi.py`, BENCHMARKING.md §6.
- Full MD loop ✓: NVE / NVT (Berendsen + NHC) / NPT
- Force trajectory output (@trf) ✓: FREEFORCERED + CONSFORCERED, atom-by-atom tolerance 1e-6
- **37/40 reference tests passing** (3 ignored: `aladip_vacuum_fep` known FEP mismatch,
  `aladip_vacuum_em` EM frame-count off-by-one, `water_216_nve_nobath` — a *correct* reference the
  engine fails because the IMD parser turns an absent MULTIBATH block into a Berendsen bath, 3.9 A18 — the full matrix is below)
- `PairlistAlgorithm::from_imd(imd.algorithm, …)` dispatches Standard/CellList here (9a-0/
  9a-1); CellList only activates on explicit `ALGORITHM grid_cell`, never via a size heuristic — all
  37 active reference systems still use `standard`.
- `test_provider_reference.rs` / `test_orchestrator_reference.rs` — gromosXX-reference tests for the
  `PotentialProvider`/`ProviderOrchestrator` seam (2.6/P2.8-2), not the main `md` binary path.
  `schnet_nve_loop.rs` (feature `ml`) — real leapfrog loop driving `SchNetInteraction`;
  `xtb_orchestrator_sequence.rs` — `ProviderOrchestratorAlgorithm` inside a real `AlgorithmSequence`
  (2.8-6 ✓). Neither is wired into the `md` CLI yet.
- **3.9 step 1 ✓ (2026-08-29): the run assembly left this crate.** `md.rs` (2060 → 1359
  lines) now calls `gromos_run::{prepare_system, build_sequence_from_imd, start}` — the same functions
  `pyo3-gromos` calls — and only keeps CLI parsing, printing (from `PrepareNote`/`BuildSummary`),
  writers, and the GaMD/EDS in-loop application (still applied *after* `run_step`, still parsed
  out-of-band; 3.9 A6). Do not add algorithm construction back here, and do not add a
  second IMD reader.
- **3.9 step 2 ✓ (2026-08-29).** The binary now runs from a `RunRecipe`
  (`gromos_run::build_sequence_from_imd` = `from_imd` → `build_plan` → `validate_plan` →
  `instantiate`); it allows `GAMD`/`EDS` blocks to pass through the recipe because it applies them
  itself, and rejects any other unmodelled block. `md @dump` prints recipe + plan as JSON and exits;
  every run writes `<tre>.recipe.toml` (the effective run, diagnostics as comments — GROMACS's
  `mdout.mdp`). The absent-MULTIBATH hazard is fixed at the parser defaults: `water_216_nve_nobath`
  is a regular passing reference (38/40, 2 ignored).

- **3.9 steps 2–4 ✓ (2026-08-29): the recipe is the only entry point.** `md.rs` reads the
  `.imd` through `gromos_run::read_imd`, builds a `RunRecipe` (`GAMD`/`EDS` allowed through as
  passthrough blocks, applied out-of-band), and runs `prepare_system → build_plan → validate_plan →
  instantiate → start` — the same calls `py-gromos`'s `Simulation` makes. `md @dump` prints recipe
  + plan; `<tre>.recipe.toml` is written next to the energies. `just lint` (G6) fails on any second
  builder or reader.
- **EDS through the recipe (1.8, 0.0.57).** The `EDS` block is modelled (`RunRecipe::eds`, the
  `eds` step after the force field); the out-of-band EDS path (`EdsBlock`, `eds_stats.dat`,
  `eds_vr.dat`, the `eds` binary) is gone and only `GAMD` still passes through the recipe.

## Key files
```
src/bin/md.rs                               — main MD driver, CLI, simulation setup
tests/test_gromosXX_references.rs           — integration tests vs gromosXX
tests/gromos_check_suite.rs                 — gromosXX's check/*.t.cc per-term energies, ported
tests/run_references.py                     — generate gromosXX reference data
tests/gromosXX_references/                  — reference input + expected output
tests/test_provider_reference.rs            — LjCrfInteraction vs gromosXX (provider seam)
tests/test_orchestrator_reference.rs        — ProviderOrchestrator vs gromosXX (multi-provider seam)
tests/schnet_nve_loop.rs                    — SchNetInteraction + LeapFrog NVE loop (feature `ml`)
```

## Crate-specific rules
- **Reference test oracle = gromosXX.** A change that alters wired, tested output is a regression until proven a deliberate, documented decision.
- **DON'T edit `gromosXX_references/*/expected/`** — those are ground truth.
- **All simulation parameters from `@input` .imd/.in file** (not CLI flags).
- **CLI:** clap `#[derive(Parser)]` with `gromos_args()` pre-processor (`@key → --key`, `@f argfile` expansion). No custom arg parsers.
- **Tolerances:** force_abs=1e-6, energy_rel=1e-8, position_abs=1e-9.

## Reference systems

The per-system status matrix — which of the 45 levels pass and what each isolates — lives with the
systems themselves: `tests/gromosXX_references/README.md`.

## How to add a reference test
1. Create `tests/gromosXX_references/<name>/` with the `.imd` and an `input.toml` (copy a sibling's;
   inputs may live in `shared/`). Set `NTWG=1` in WRITETRAJ when the system perturbs (dH/dλ frames).
2. `python3 scripts/regen_gromosXX_references.py <name>` — runs the native gromosXX
   (`.local/gromosXX/md++/BUILD/program/md`) and writes `expected/`. Never re-run it without names:
   the other systems' `expected/` files are the record.
3. Add `ref_test!(name, "dir")` in `tests/test_gromosXX_references.rs` (or `ref_test!(ignore: …,
   "reason")` while it fails) and the name to `REFERENCE_SYSTEMS` in `py-gromos/tests/test_gromosXX_references.py`.