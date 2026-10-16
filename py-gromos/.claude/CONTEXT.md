# py-gromos — stage contract

## Job
L4+ user-facing Python package. The importable `gromos` package: Python wrappers, high-level
simulation runners, analysis helpers, notebooks, and examples. Built with `maturin` on top of
`pyo3-gromos`. Contains no Rust code — only Python that calls the compiled extension.

## Inputs (consumes from)
`pyo3-gromos` (compiled `.so` extension, via `maturin develop` or installed wheel).

## Outputs
- `python/gromos/` — importable Python package
  - `__init__.py` — top-level API surface
  - `md_runners.py` — high-level simulation runner wrappers
  - `tools/` — the gromos++ system-preparation programs, in process (3.10)
  - `analysis.py` — the old subprocess stub; `gromos.analysis` replaces it in 3.10 step 6
  - `gromos.pyi` — type stubs for IDE support
- `notebooks/` — Jupyter notebooks (education + demonstration)
- `tests/` — Python-level integration tests

## Status
- Basic wrappers ✓; `maturin develop` builds ✓
- P3.3 ✓ — `EnergyTimeseries` (`python/gromos/timeseries.py`): `sim.run()` → numpy array wrapper,
  `block_average()`, `to_dataframe()` (polars/pandas/dict), `plot()` (plotly/matplotlib); backends
  configurable via `gromos.timeseries.config` (default polars + plotly).
- P3.4 ✓ — notebooks rewritten: `01_load_and_inspect.ipynb` (System/Topology/Config, RDF, mass hist,
  atom-count guard, manual `solvate()`+`System()` composition), `02_short_md.ipynb` (native params,
  `run()`→`EnergyTimeseries`, component breakdown, block-average, NVE-vs-NVT). Legacy notebooks deleted.
  `md_runners.py` marked deprecated. `notebooks` dep group in `pyproject.toml`; CI/Makefile use it.
- P3.5 ✓ (M1/M2/SD) — `02_short_md.ipynb` §2a demonstrates the constraint fix live:
  `InputParameters.nvt(..., constraints="none")` on `aladip_solvated` diverges to ~95,000 K within 60
  steps (baked-in output, kept as documented contrast); `constraints="hbonds"` stays stable and tracks
  300 K. `SteepestDescent` + `AlgorithmSequence.minimize()` now importable from `gromos`; `.pyi` updated.
- P3.6 ✓ — `tests/test_gromosXX_references.py`'s `REFERENCE_SYSTEMS` now covers all 37 active Rust
  reference systems except the 2 explicitly-deferred FEP ones (was missing 18, all real Python-API
  defects in `pyo3-gromos`, not test gaps — see that crate's CONTEXT.md). 121 tests (118 passed, 2
  documented `POSITION_MISMATCH_SYSTEMS` skips for a pre-existing Rust-side `.trc` divergence).
- **3.9 step 0 ✓ (2026-08-29)** — `tests/test_front_end_parity.py`: every reference system
  through `Simulation(...)` and through `AlgorithmSequence.from_parameters` → `np.array_equal`
  (never `allclose`), plus a `same_path_twice` determinism baseline. 27/37 systems are bit-identical
  between the two builders; the 10 divergences are `xfail(strict=True)` entries in
  `EXPECTED_DIVERGENCE`, each naming the feature path C lacks — they must be *deleted* when the shared
  builder lands, never left to rot. `test_gromosXX_references.py` gained `EXPECTED_ENGINE_FAILURES`
  (`water_216_nve_nobath`, the absent-MULTIBATH parser bug) and `REFERENCE_PARAMS`.
- **3.9 step 1 ✓ (2026-08-29)** — `Simulation` is built by the `md` binary's own code
  (`gromos-run`); `Topology.apply_perturbation(path)` brings FEP to Python (`ch4_water_fep` in
  `REFERENCE_SYSTEMS`, passing; `aladip_vacuum_fep` strict xfail). Full suite: 208 passed, 16 skipped,
  18 xfailed (10 descriptor-path divergences + 2 FEP-on-path-C + nobath ×3 + aladip_vacuum_fep ×3).
- **3.9 step 2 ✓ (2026-08-29)** — `water_216_nve_nobath` passes (parser defaults fixed; xfail
  removed); `tests/test_recipe.py` pins `Simulation.recipe_toml` / `.plan_json` / `.diagnostics`.
  Suite: 211 + 4 passed, 16 skipped, 15 xfailed.
- **3.9 step 3 ✓ (2026-08-29)** — `Recipe`/`Term`/`Algorithm`/`Plan`, `gromos.terms()`/
  `algorithms()`/`build_info()`, `gromos.exceptions`. `test_gromosXX_references.py` builds every
  reference through `Simulation(system, Recipe.from_imd(...).with_inputs(...))`;
  `test_front_end_parity.py` adds path B (Recipe) and path D (Recipe + its own Plan) against path A,
  no divergences, plus the registry guards `test_every_kind_has_a_parity_case` (exemptions in
  `UNCOVERED_KINDS` must name the step that removes them) and `test_pyi_lists_every_kind`
  (`_AlgorithmKind`/`_TermKind` literals in `gromos.pyi`). `gromos.pyi` is `stubtest`-clean:
  `MYPYPATH=python uv run python -m mypy.stubtest gromos.gromos --allowlist
  stubtest_allowlist_no_ml.txt` (pyo3 classes are `@final` and construct through `__new__`; slot
  methods take positional-only arguments). Deprecated forms warn (`pytest.warns`; the pytest config
  ignores `DeprecationWarning` otherwise). Suite: 308 passed, 16 skipped, 15 xfailed.
- **3.9 step 4 ✓ (2026-08-29)** — the descriptor path is gone from the package surface:
  no `Forcefield`/`LeapFrog…`/`BerendsenThermostat`… classes, no `md_runners`, `analysis` only by
  explicit import; `AlgorithmSequence.*` presets return a `Plan` (deprecated), `Simulation.from_sequence`
  takes one. `test_front_end_parity.py`: paths A/B/D + the `from_sequence` shim, no xfail table.
  `Vec3`/`Frame`/`rmsd`/`rdf` are f64. Suite: 300 passed, 8 skipped, 3 xfailed.
- **3.9 step 5 ✓ (2026-08-29)** — `tests/test_xtb_term.py` is the physics oracle for
  `Term("xtb", …)` (water-dimer fixture from `gromos-forces/tests/gromosXX_qmmm_references`; skips
  without `xtb`): exact additivity against `XtbPotential`, per-term visibility
  (`sim.term_energies`), NVE drift, and the named rejections. Pattern for the next term: no Python
  change beyond the `_TermKind` literal and a test. Suite: 306 passed, 8 skipped, 3 xfailed.
- **3.10 steps 0–5 ✓ (0.0.49)** — `gromos.tools`: `run(name, **kwargs)` over all 33
  gromos-tools programs in process (the same Rust code the binaries run), `programs()` with each
  program's keys/usage/status, typed wrappers for the preparation chain (`make_top`, `com_top`,
  `red_top`, `pdb2g96`, `sim_box`, `build_box`, `ion`, `check_top`) taking and returning
  `Topology`/`Configuration`/`CheckReport`, and `System.solvate()`/`.add_ions()`/`.check()` over
  the same verbs. `gromos.ProgramError` carries the program's partial stdout. `tests/test_tools.py`
  (26) pins text ↔ object agreement across the boundary; the Rust side pins library ↔ binary.
- **Free-energy windows (1.8, 0.0.59)** — `gromos.free_energy`: `LambdaSchedule(points,
  interaction_lambdas)` → one `Recipe` per λ, `run_windows(system, recipes, workdir, dry_run)` in
  process (kept windows `<workdir>/l0.250.npz` are read, not rerun), `collect(workdir)`,
  `FreeEnergyWindows.ti()/.report()` (the `gromos_analysis::ti` estimator).
  `tests/test_free_energy.py` reproduces every dH/dλ sample of the four `ch4_water_fep_l*`
  references. `system_builder.py` no longer has a second `EnergyTimeseries`. User guide:
  `docs/user-guide/free-energy.md`.
- **BAR/MBAR of the windows (1.8, 0.0.62)** — a window run with `PRECALCLAM` keeps
  `Simulation.precalclam` of every sampled frame (`Window.precalclam`, also in its `.npz`);
  `FreeEnergyWindows.bar(skip)`/`.mbar(skip)` → `Estimate` (total and per neighbouring pair),
  `.bar_data(skip)`: `ext_ti_ana`'s `@bar_data` energies, `bar` and `mbar` in Rust through
  `_precalc_bar_data`/`_precalc_bar`/`_precalc_mbar`, no estimator arithmetic in Python. The grid
  must hold every window's λ; `LAMBDAS` windows are refused. Tests: the nine `*precalc*`
  references' rows frame by frame (rel 1e-9), `bar()`/`mbar()` = `analysis.bar`/`mbar` on the same
  frames (1e-10), TI/BAR/MBAR agreement on nine CH4 windows.
- **Charge changes end to end (1.8, 0.0.64)** — `Perturbation.mutate(system, atoms, charge_b=…,
  iac_b=…, mass_b=…)` / `.from_ptp(path, topology)` (the `.ptp` text from gromos-tools'
  `make_pt_top`, net charges as `Topology.apply_perturbation` gives them, `.end_state()` by
  `pt_top`); `run_windows(..., perturbation=, parallel=N)`: N spawn worker processes (a
  `Simulation` is unsendable and holds the GIL, so threads would not overlap), the same numbers
  bit for bit; `FreeEnergyWindows.window_statistics(skip)`. `scripts/ti_ch4.py` runs gromos-rs
  through the module and gromosXX on the same recipes' `.imd`. `tests/test_free_energy_charge.py`:
  CH4 → CH4⁺ in 999 SPC, seven windows, TI/BAR/MBAR within 2σ, the `dGslv_pbsolv` correction on
  the end states against gromos++ (0.01 kJ/mol).
- Charge-change correction chain (1.8, 0.0.66): `gromos.analysis.charge_correction`,
  `rf_self_energy` (the reaction-field self term the raw ΔG holds), `corrected_charging_dg`
  ([raw − self] + PB(B) − PB(A) + C1 + C2, or + the five analytic terms); user guide "Correcting
  a charge change". `scripts/na_hydration.py` runs Na⁺ → dummy in SPC through the whole chain
  (`docs/free_energy_validation.md`).
- Remaining P3 items:
  - [ ] `gromos.analysis` over the 82 analysis programs — 3.10 steps 6–8. The stub
        `analysis.py` (subprocess, mostly `NotImplementedError`) is deleted in step 6.
  - [ ] Rich `__repr__` / `_repr_html_` for Jupyter (Topology, Configuration, Energy)

## Key files
```
python/gromos/__init__.py        — top-level API (the module contract; System/Recipe/Plan/Term/Simulation)
python/gromos/gromos.pyi         — type stubs (stubtest-clean; `_AlgorithmKind`/`_TermKind` literals)
python/gromos/exceptions.py      — RecipeError / PlanError / MissingFeatureError / RunError
python/gromos/timeseries.py      — EnergyTimeseries over run()'s array
python/gromos/free_energy.py     — λ windows: Perturbation, LambdaSchedule, run_windows (parallel), TI, BAR, MBAR
python/gromos/tools/             — the gromos++ programs: run()/programs() + the typed prep chain
python/gromos/analysis.py        — the old subprocess stub; replaced by gromos/analysis/ in 3.10 step 6
python/gromos/system_builder.py  — P3 design sketch: ForceField→BuildingBlock→Topology algebra
stubtest_allowlist_no_ml.txt     — the ml-only names, for stubtest on a default build
tests/test_gromosXX_references.py, test_front_end_parity.py, test_recipe.py, test_xtb_term.py,
                                 test_tools.py (gromos.tools: text ↔ object, the prep chain)
docs/                            — mkdocs (user-guide/recipe.md is the concept page)
notebooks/                       — 03 two argon atoms (the first-run tutorial; asserts its own numbers),
                                   01 load & inspect, 02 short MD (on Recipe). The superseded API
                                   mockup is docs/development/api-design-mockup.ipynb (not executed)
pyproject.toml                   — package metadata (maturin)
```

## Design target (P3 algebra)

```python
ff = ForceField("54a7")
mol = ff.molecule(["STA", "ALA", "GLY", "END"])
system = mol * 10 + ff.solvent("SPC") * 216
solvated = system.neutralize("Na+").solvate("SPC", density=900)
result = solvated.minimize().equilibrate().run(steps=1_000_000)
df = result.energies.to_dataframe()
```

Three documented design shapes in system_builder.py:
- Shape A (OpenMM-style): `ff.build_topology(seq)` — implemented
- Shape B (builder/fluent): `* / +` algebra — implemented
- Shape C (declarative dict): `build_system(forcefield="54a7", molecules=[...])` — implemented

Motivation: vsomm_modeler required manual tracking of atom counts, binary
paths as plain strings, and `input_xxx` dicts per step — all encapsulated now.

- **Composition model (2026-08-29): 3.8/3.9.** The Python API should build a `Recipe`
  (run control ⟂ forcefield+terms ⟂ constraints ⟂ ensemble) that the engine consumes; IMD files are
  the serialisation of the same recipe. `test_basic.py`/`test_advanced_features.py` are placeholders
  for a past API (eleven skipped classes reference non-existent types) — replace, don't extend.

## Crate-specific rules
- **No physics, no data structures** — this layer only wraps what pyo3-gromos exposes.
- **API design reference:** model method chaining and DataFrame interop on Polars patterns (`.local/polars`).
- Changes here require `maturin develop` in this directory to rebuild the extension before testing.
- `gromos.abi3.so` is a compiled artifact — never edit it directly; rebuild via maturin.
