# newen-md — stage contract

## Job
L3 orchestration. The `md` binary: composes newen-core + newen-forces + newen-integrators +
newen-io into a full MD simulation loop matching gromosXX output bit-for-bit.

## Inputs (consumes from)
newen-run (the run assembly) and the lower crates it drives (newen-core, newen-forces,
newen-integrators, newen-io, newen-logic for the output schedule), plus the `newen` facade.

## Outputs
`md` binary: runs GROMOS MD simulations; writes trajectories, energies, force trajectories
matching gromosXX to reference tolerances.

## Status
- **Nosé–Hoover restart (0.0.69, D41):** the final configuration carries NHCVARIABLES (the ζ of
  every bath, `G96Writer::write_nhcvariables`, from `Configuration::nose_hoover`) for a chain and
  for plain Nosé–Hoover. References `water_216_nvt_{nosehoover,nhc}_cont_fixed` (D41 patch; the
  released-binary reproducer lives in gromos-divergences); `tests/nh_restart.rs`: 100 + 100 steps
  follow 200 as closely as weak coupling does.
- **Output times and schedules (0.0.68):** every written time is T + step·DT
  (`SimulationState::time_of`, STEP T from the recipe); forces follow NTWF; a converged
  minimisation repeats the previous step's force, energy and dH/dλ frame on that step's
  NTWF/NTWE/NTWG schedule and writes the stopping step's positions only on the NTWX schedule; the
  final configuration has gromosXX's TIMESTEP block and title layout. The comparator checks every
  TIMESTEP block and the force frame count. References `aladip_vacuum_fep_slowgrowth_t10`,
  `aladip_vacuum_em_ntw5`. Since 2026-10-22 the decisions are `newen_logic::schedule::output`
  (`due`, `minimisation_stop`, `final_step`; proven to write gromosXX's top-of-loop and final
  frames), and NTWV ≠ 0 in a minimisation is refused as gromosXX refuses it
  (`check_parameter.cc:62-65`), so a minimisation writes no velocity frames.
- **Slow growth (1.8, 0.0.61):** `md` calls `newen_run::advance_lambda` after each step (not on
  the step the minimiser stops at), writes the λ the step ran at into the `.trg` (`old_lambda`)
  and, on a converged minimisation, the repeated `.trg` frame next to the repeated `.tre` one.
  References `aladip_vacuum_fep_slowgrowth{,_shake,_cont}`, `ch4_water_fep_slowgrowth`,
  `aladip_vacuum_em_slowgrowth`.
- **Net charge (1.8, 0.0.59):** `md` prints `Perturbation net charge: A … B … (ΔQ = …)` on a
  perturbed run and warns when ΔQ ≠ 0 under a reaction field (`fep_null_perturbation` pins
  ΔQ = 0 and no warning for the null perturbation). Since 0.0.65 a charge change has native
  references: `nacl_water_fep{,_precalc,_lambdas,_14}` (Na⁺ → dummy, ΔQ = −1, under RF: the
  perturbed RF self term, the RF excluded and the 1-4 Na–Cl pair).
- **MPI (2026-08-30):** `md` built with `--features use-mpi` (`LIBCLANG_PATH=/usr/lib/llvm-19/lib`)
  runs under `mpirun -np N`: every rank runs the whole driver on the same inputs, the pair terms are
  split by `RunOptions::pair_partition` and summed each step (`mod mpi_support` at the bottom of
  `md.rs`: reduce to rank 0 + broadcast, 3N+13 doubles), rank 0 alone writes files and stdout.
  `np=1` is bit-identical to the serial binary. The old `md_mpi`/`md_mpi_cuda`/`repex_mpi` bins are
  gone. Timing: `scripts/bench_mpi.py`, BENCHMARKING.md §6.
- Full MD loop ✓: NVE / NVT (Berendsen + NHC) / NPT
- Force trajectory output (@trf) ✓: FREEFORCERED + CONSFORCERED, atom-by-atom tolerance 1e-6
- **135 reference tests registered, none ignored, all passing** in
  `tests/test_gromosXX_references.rs` (the matrix is below)
- **Lattice sums and boxes (2026-10-22, 2026-10-23):** `md` writes GENBOX (NTB, lengths, angles)
  into the trajectory and the final configuration, and prints "SPME is experimental: no gromosXX
  oracle; …" with the header when the plan's lattice sum is experimental (NLRELE 5).
- `PairlistAlgorithm::from_imd(imd.algorithm, …)` dispatches Standard/CellList here (9a-0/
  9a-1); CellList only activates on explicit `ALGORITHM grid_cell`, never via a size heuristic —
  every active reference system uses `standard` except `water_1000_spc_gridcell`.
- `test_provider_reference.rs` / `test_orchestrator_reference.rs` — gromosXX-reference tests for the
  `PotentialProvider`/`ProviderOrchestrator` seam (2.6/P2.8-2), not the main `md` binary path.
  `schnet_nve_loop.rs` (feature `ml`) — real leapfrog loop driving `SchNetInteraction`;
  `xtb_orchestrator_sequence.rs` — `ProviderOrchestratorAlgorithm` inside a real `AlgorithmSequence`
  (2.8-6 ✓). Neither is wired into the `md` CLI yet.
- **3.9 step 1 ✓ (2026-08-29): the run assembly left this crate.** `md.rs` (2060 → 1359
  lines) now calls `newen_run::{prepare_system, build_sequence_from_imd, start}` — the same functions
  `pyo3-newen` calls — and only keeps CLI parsing, printing (from `PrepareNote`/`BuildSummary`),
  writers, and the GaMD/EDS in-loop application (then still applied *after* `run_step`; since
  moved: EDS into the recipe in 0.0.57, GaMD into the force field in 0.0.69). Do not add
  algorithm construction back here, and do not add a second IMD reader.
- **3.9 step 2 ✓ (2026-08-29).** The binary now runs from a `RunRecipe`
  (`newen_run::build_sequence_from_imd` = `from_imd` → `build_plan` → `validate_plan` →
  `instantiate`); it allows `GAMD`/`EDS` blocks to pass through the recipe because it applies them
  itself, and rejects any other unmodelled block. `md @dump` prints recipe + plan as JSON and exits;
  every run writes `<tre>.recipe.toml` (the effective run, diagnostics as comments — GROMACS's
  `mdout.mdp`). The absent-MULTIBATH hazard is fixed at the parser defaults: `water_216_nve_nobath`
  is a regular passing reference (38/40, 2 ignored).

- **3.9 steps 2–4 ✓ (2026-08-29): the recipe is the only entry point.** `md.rs` reads the
  `.imd` through `newen_run::read_imd`, builds a `RunRecipe` (`GAMD`/`EDS` allowed through as
  passthrough blocks, applied out-of-band), and runs `prepare_system → build_plan → validate_plan →
  instantiate → start` — the same calls `py-newen`'s `Simulation` makes. `md @dump` prints recipe
  + plan; `<tre>.recipe.toml` is written next to the energies. `just lint` (G6) fails on any second
  builder or reader.
- **EDS through the recipe (1.8, 0.0.57).** The `EDS` block is modelled (`RunRecipe::eds`, the
  `eds` step after the force field); the out-of-band EDS path (`EdsBlock`, `eds_stats.dat`,
  `eds_vr.dat`, the `eds` binary) is gone and only `GAMD` still passes through the recipe.
- **A-EDS through the recipe (1.8, 0.0.60).** The `AEDS` block runs (FORMS 1-6); `md` writes the
  search state as the final configuration's `AEDSSEARCH` block (FORMS 2-6) and a run with
  NTIAEDSS = 0 continues it. References `aladip_vacuum_aeds_*`; `tests/aeds_runs.rs` checks a
  two-job search against one job and the refused inputs. Since 0.0.63 D16, D25, D32 are
  corrected; the searches are pinned by `aladip_vacuum_aeds_*_fixed` (the patched oracle; the
  released-binary reproducers live in gromos-divergences); `aeds_runs.rs` also checks that a
  frame's V_mix/V_R* follow from its own columns and that the E_min target is a sampled state.

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
- **CLI:** `md` still parses its `@key`s in a hand-rolled loop over `env::args()` (`parse_args`,
  `src/bin/md.rs`), without `newen_io::normalize`, so no `@f` expansion: an exception to the
  workspace CLI rule (overview, "CLI arg parsing"), not a pattern to copy. New binaries follow
  the rule.
- **Tolerances:** force_abs=1e-6, energy_rel=1e-8, position_abs=1e-9.

## Reference systems

The per-system status matrix — which of the 135 registered systems pass and what each isolates —
lives with the systems themselves: `tests/gromosXX_references/README.md`.

## How to add a reference test
1. Create `tests/gromosXX_references/<name>/` with the `.imd` and an `input.toml` (copy a sibling's;
   inputs may live in `shared/`). Set `NTWG=1` in WRITETRAJ when the system perturbs (dH/dλ frames).
2. `python3 scripts/regen_gromosXX_references.py <name>` — runs the native gromosXX
   (`.local/gromosXX/md++/BUILD/program/md`, released 247c95c4) and writes `expected/`. Never
   re-run it without names: the other systems' `expected/` files are the record. Where newen
   corrects a gromosXX defect, the system is `<name>_fixed` and is generated with
   `--binary .local/gromosXX-patched/md++/BUILD/program/md` (D1, D2, D16, D25, D32, D40, D41,
   LS1–LS5; `patches/BUILD.md`), its `input.toml` naming the patch. An Ewald input needs
   gromosXX's `@develop` flag, which the script does not pass: run `md` by hand and record the
   command in the `input.toml` header (see `nacl_pair_box_ewald/input.toml`).
3. Add `ref_test!(name, "dir")` in `tests/test_gromosXX_references.rs` (or `ref_test!(ignore: …,
   "reason")` while it fails) and the name to `REFERENCE_SYSTEMS` in `py-newen/tests/test_gromosXX_references.py`.