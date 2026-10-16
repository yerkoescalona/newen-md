# pyo3-gromos — stage contract

## Job
L4 Python bindings (PyO3). A thin wrapper over the Rust core — no physics, no data duplication.

## Inputs (consumes from)
All Rust crates.

## Outputs
Python-callable API for running simulations and analysing trajectories.

## Status
- **Charge-change corrections (1.8, 0.0.65–0.0.66):** `_charge_correction` and `_rf_self_energy`
  (`programs.rs`) bind `gromos_analysis::free_energy::{charge_correction, rf_self_energy}` for
  `gromos.analysis.charge_correction` / `rf_self_energy`.
- **Perturbations and parallel windows (1.8, 0.0.64):** `_make_pt_top_mutate`, `_ptp_net_charge`,
  `_pt_top` (`programs.rs`) bind gromos-tools' `make_pt_top::{mutate, make_pt_top}` and
  `pt_top::pt_top` for `gromos.free_energy.Perturbation`. `System` pickles when made by
  `System.from_files` (as `from_files` on the absolute paths; `module = "gromos"`), which
  `run_windows(parallel=N)` needs to hand it to worker processes; one built in memory raises
  `TypeError`. `Simulation` raises a `UserWarning` when slow growth takes λ out of [0, 1]
  (`advance_lambda`'s return, in the constructor, `step` and `run`) or the run starts outside it
  (`PrepareNote::LambdaOutOfRange`).
- **Extended TI in process (1.8, 0.0.62):** `Simulation.precalclam` (the grid, the `# precalclam`
  energy and derivative rows, `# ABdih`) of the last step; `_precalc_bar_data`, `_precalc_bar`,
  `_precalc_mbar` (`programs.rs`) take per-window frame stacks and call
  `ext_ti_ana::bar_file` → `bar::bar` / `mbar::{reduced_potentials, mbar}`
  (`gromos.free_energy.FreeEnergyWindows.bar/.mbar`).
- **Slow growth (1.8, 0.0.61):** `Simulation` calls `gromos_run::advance_lambda` after the
  constructor's step 0 and after every `step`/`run` step, as `md` does (no "converged" stop, so
  under steepest descent λ also moves on the step `md` stops at).
- `System`, `Recipe`/`Term`/`Algorithm`/`Plan`, `Simulation(system, recipe[, plan=])` ✓ (3.9,
  2026-08-29); `InputParameters` and the `AlgorithmSequence.*` presets are deprecation shims
- `sim.run(steps, ene_freq)` → numpy array; energy decomposition (bond/angle/dihedral/improper/lj/coulomb) ✓ (P3.3)
- Python tests: 121 (118 passed, 2 documented position-mismatch skips) ✓ (`py-gromos/tests/`)
- `.pyi` stubs ✓
- **3.6 ✓ done** — `build_simulation()` (`simulation.rs`) now dispatches every gromosXX
  feature that was already implemented lower in the stack but silently unwired from Python: NTIVEL=1
  velocity generation, SETTLE/LINCS (`ConstraintSelection::from_imd`), Nosé-Hoover (+chain)
  thermostat (`push_thermostat()`), distance/position restraints (`distrest=`/`posresspec=`/
  `refpos=` kwargs), triclinic/truncated-octahedron box (`NTB=-1`, plus `sim.forces` frame
  rotation). `REFERENCE_SYSTEMS` (`py-gromos/tests/test_gromosXX_references.py`) now covers all 37
  active Rust reference systems except the 2 explicitly-deferred FEP ones. **Known, tracked
  divergence:** `gromos-rs`'s own `.trc` position output disagrees with gromosXX for
  triclinic/EM-shift-frame systems (pre-existing, Rust-side, not a Python binding bug) — see
  `POSITION_MISMATCH_SYSTEMS` in the test file.
- **3.5 ✓ done (M1, M2, SD)**:
  - **M1**: `InputParameters.nve/nvt/npt` take `constraints="none"|"hbonds"|"allbonds"` (default
    `"none"`) → `ntc=1|2|3` (`ntc_from_constraints()`, `gromos-io/imd.rs`); readable `.constraints`
    getter. No builder change needed — the existing `shake_enabled = imd.ntc > 1 || …` guard picks it
    up. Verified empirically: factory `nvt(..., constraints="hbonds")` on `aladip_solvated` stays
    stable at `dt=0.002`; `constraints="none"` on the same system still diverges (documented
    contrast, not a regression). Test: `test_constrained_system_stable_with_factory_params`.
  - **M2**: `temperature` getter now uses `PySimulation.total_dof` (constraint-aware, computed once
    via the shared `compute_total_dof()` in `algorithm_sequence.rs` — same value both builders and
    the thermostat use) instead of bare `n_atoms*3`.
  - **SD**: `build_simulation()` branches on `imd.ntem > 0` — pushes `SteepestDescentAlgorithm`
    instead of leap-frog, skips thermostat/barostat/`TemperatureCalculation`/COM removal (mirrors
    the `md` binary's own EM sequence). New `SteepestDescent` algorithm-sequence building block +
    `AlgorithmSequence.minimize()` preset; `from_parameters()` dispatches to it when `ntem > 0` — the
    composable path and `Simulation(system, params)` now agree. `shake_algorithm_from_imd()`
    extracted so EM and standard-MD SHAKE construction can't drift apart. Verified: `aladip_vacuum`
    under `steepest_descent(steps=50)` actually drops -13.5→-47.98 kJ/mol (previously a silent no-op
    at `dt=0`). Tests: `test_steepest_descent_via_simulation`, `test_steepest_descent_via_algorithm_sequence`.
- `sim.volume`/`sim.pressure` getters (mirror `temperature`'s shape). **`pressure` is only
  physically meaningful under NPT** — the virial term is only populated by `PressureCalculation`,
  which only NPT's sequence includes; under NVE/NVT it returns the misleading kinetic-only term, not
  zero. Documented in the getter, `.pyi`, API reference, and `02_short_md.ipynb`. Test:
  `test_volume_and_pressure_getters`.
- Remaining P3 items:
  - [ ] Expose ForceField evaluation (single-point energy/force)
  - [ ] Rich `__repr__` / `_repr_html_` for Jupyter (Topology, Configuration, Energy)
  - [x] 3.7 — ML potential binding: `PySchNetPotential` (`--features ml`, libtorch via the
    venv's torch — INSTALL.md §5); the ML term is `Term("schnet", ...)` on the recipe since 3.9 step 3.
- The July 2026 "two builders / scattered defaults" audit is done by 3.9 (steps 0–4): one
  builder in `gromos-run`, defaults derived from `ImdParameters::default()`, the recipe as the one
  entry point.

- **Composition audit + target model (2026-08-29): 3.8/3.9.** Three IMD→sequence builders
  existed (`md.rs`, `build_simulation`, `resolve_algorithm_sequence`), the descriptor enum is closed,
  and force terms arrive as `Simulation` kwargs. Target: one `RunRecipe` in `gromos-run`, one
  builder, IMD and Python objects as front-ends, providers as `terms`. Do not add new kwargs to
  `Simulation` or new `AlgorithmDescriptor` variants without reading 3.9 first.
- **3.9 step 1 ✓ (2026-08-29).** `simulation.rs::build_simulation` contains no algorithm
  construction: it calls `gromos_run::{prepare_system, build_sequence_from_imd, start}` (the `md`
  binary's code), inserts the optional ML term after `Forcefield` via `AlgorithmSequence::insert`,
  and maps `RunError` onto the builtin exceptions the binding always raised (`run_err`). Python now
  honours the topology's `four_pi_eps_i`, takes NSM from the coordinate file, and uses the binary's
  parallel-kernel policy. FEP: `Topology.apply_perturbation(path)` (no `Simulation` kwarg),
  which returns the solute's net charge in state A and B (0.0.59).
  `resolve_algorithm_sequence` (the descriptor path) survives until steps 3–4 with the same
  parallel policy/`four_pi_eps_i` so `test_front_end_parity.py` stays exact.
- **3.9 step 2 ✓ (2026-08-29).** `build_simulation` goes through the recipe path
  (`gromos_run::build_sequence_from_imd`, strict passthrough policy: no unmodelled block reaches
  Python silently); `Simulation.recipe_toml`, `.plan_json`, `.diagnostics` expose the effective
  recipe, the resolved plan and the parser notes (`None`/empty on the descriptor path). Step 3 turns
  the recipe into `gromos.Recipe`/`Term`/`Algorithm` and retires `InputParameters` + the kwargs.
- **3.9 step 3 ✓ (2026-08-29).** `src/recipe.rs` binds the `gromos-run` types as
  `Recipe`/`Term`/`Algorithm`/`Plan` through `pythonize` (serde ↔ Python; dict/TOML/JSON/pickle
  round trips, `update(**groups)` deep-merges through `serde_json::Value` and re-deserialises with
  `deny_unknown_fields`, so a typo is a `RecipeError`). `simulation.rs` has one entry
  (`build_from_recipe`: `prepare_system` → `build_sequence_from_recipe|from_plan` → `start`);
  every other constructor form (`InputParameters`, the restraint/ML kwargs, `from_sequence`) is a
  deprecation shim in Rust (`parameters::deprecated`) that *translates into a recipe* — not a
  second path (`test_front_end_parity.py` proves it with `array_equal`). Exceptions are
  `create_exception!` classes re-exported by `gromos.exceptions`. `System` accepts a solute
  topology with a solvated coordinate file (`system::atom_count_ok`). Deviations from the plan
  text: shims live in Rust, not `_deprecation.py`; the plan is a `Plan` class from
  `recipe.plan(system)` rather than `AlgorithmSequence.from_recipe`; provenance (model checksum)
  and the per-term energy slot (G10) are not done. The descriptor path
  (`algorithm_sequence.rs`) is now dead weight — step 4 deletes it.
- **3.9 step 4 ✓ (2026-08-29).** The descriptor path is gone: `algorithm_sequence.rs`
  (1473 → ~100 lines) keeps only the deprecated `AlgorithmSequence.nve/nvt/npt/minimize/
  from_parameters` names, each returning the `Plan` of the parameters (`gromos_run::build_plan` on
  the topology); `Simulation.from_sequence` takes that `Plan`. `Simulation.recipe`/`.plan` are
  always present. `Vec3`/`Frame`/`rmsd`/`rdf` are f64 (`gromos_core::math::Vec3`). Recipe- and
  atom-count errors are `gromos.exceptions.RecipeError`. The IMD is read through
  `gromos_run::read_imd` only. **The recipe is the only entry point**: `simulation.rs` has one
  constructor path (`build_from_recipe`), and `just lint` (G6) fails on any `AlgorithmSequence::new()`
  / `.push(Box::new(` / `read_imd_file(` outside `gromos-run`.
- **3.9 step 5 ✓ (2026-08-29).** `Term("xtb", …)` needed **no change in this crate** — the
  proof of the door. `Simulation.term_energies` (dict, registry name → kJ/mol, from
  `Energy.term_energies`) is the G10 slot; `potential_energy` stays the classical force field.

## Crate-specific rules
- **Thin wrapper only.** Zero physics, zero data structures that duplicate the Rust core.
- **API design reference:** study Polars pyo3 patterns (`PyDataFrame` / `PyExpr` / `PyLazyFrame`) at `.local/polars/py-polars/src/`.
- `py-gromos` (the Python package, separate maturin build) sits alongside this crate; `maturin develop` must pass after any change here.
