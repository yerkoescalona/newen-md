# Learning Guide

A guide to understanding the py-gromos codebase and where things are headed.

## What exists and what is stub

The cleanest way to understand the current state is to look at `__init__.py`:
it imports only what actually works from the Rust extension. Everything that
raises `NotImplementedError` or shells out to a binary lives in sub-modules
(`analysis`, `system_builder`), importable explicitly (`import gromos.analysis`) but not
in the default namespace. The legacy `md_runners` subprocess wrappers were removed in
3.9 step 4 — `Simulation` is the way to run MD from Python.

| Name | Backed by | Status |
|------|-----------|--------|
| `System`, `Topology`, `Configuration` | Rust | ✅ working |
| `Recipe`, `Term`, `Algorithm`, `Plan` + `terms()`/`algorithms()`/`build_info()` | Rust (`gromos-run` types via pythonize) | ✅ working |
| `Simulation` | Rust | ✅ working |
| `InputParameters` + factories, `AlgorithmSequence` | Rust | ⚠️ deprecated (one release): translated into a `Recipe`, warn |
| `Vec3`, `Energy`, `Frame`, `rmsd`, `rdf` | Rust | ✅ working |
| `analysis.*` | Python subprocess → gromos++ programs | 🔜 programs not yet ported |
| `system_builder.*` (ForceField, molecule, …) | Python stub | 🔜 design sketch, raises NotImplementedError |

## Source layout

```
py-gromos/
├── python/gromos/
│   ├── __init__.py        ← re-exports working names; the module contract
│   ├── gromos.pyi         ← type stubs for the Rust extension (stubtest-clean)
│   ├── exceptions.py      ← RecipeError / PlanError / MissingFeatureError / RunError
│   ├── timeseries.py      ← EnergyTimeseries over run()'s array
│   ├── analysis.py        ← future analysis wrappers (mostly stub; explicit import)
│   └── system_builder.py  ← future system-builder design sketch
├── tests/
│   ├── test_gromosXX_references.py    ← 40-system reference validation suite (the Recipe front-end)
│   ├── test_front_end_parity.py       ← every front-end against every other, exact equality
│   ├── test_recipe.py                 ← Recipe / Term / Plan / bundles / deprecation shims
│   └── test_xtb_term.py               ← the QM term's physics oracle (skips without xtb)
└── docs/                  ← this documentation
```

The Rust side lives in `crates/pyo3-gromos/src/`:

| File | Exposes |
|------|---------|
| `system.rs` | `System` |
| `topology.rs` | `Topology` |
| `py_conf.rs` | `Configuration` |
| `recipe.rs` | `Recipe`, `Term`, `Algorithm`, `Plan`, `terms`/`algorithms`/`build_info`, the exceptions |
| `parameters.rs` | `InputParameters` + factories (deprecated shims) |
| `simulation.rs` | `Simulation` |
| `algorithm_sequence.rs` | `AlgorithmSequence` — deprecated shim whose presets return the `Plan` of the parameters |
| `lib.rs` | `Vec3`, `Energy`, `Frame`, `rmsd`, `rdf`, module assembly |

## Running the test suite

```bash
make test-python           # build + full suite

# Just the reference energy/force/position tests
.venv/bin/pytest py-gromos/tests/test_gromosXX_references.py -v

# Just the workflow tests
.venv/bin/pytest py-gromos/tests/test_gromosXX_references.py \
    -k "system_constructor or factory_workflow or factory_nvt" -v
```

The reference suite validates against double-precision gromosXX output for
every system in its `REFERENCE_SYSTEMS` list: vacuum pairs, single molecules, small solvated boxes, bulk water
(NVE/NVT/NPT), and solvated alanine dipeptide.

## Roadmap

### Energy reporters ✓ done

`sim.run(steps, ene_freq)` streams energies into a NumPy array without writing
a `.tre` file; `EnergyTimeseries` (`gromos.timeseries`) wraps it with
`block_average()`, `to_dataframe()` (polars/pandas/dict), and `plot()`
(plotly/matplotlib). The `angle`/`dihedral`/`improper` energy components are
wired end to end — including a fix in the underlying Rust bonded-force
combiner (`gromos-forces/src/bonded/mod.rs`), which lumped all four bonded
terms into `bond_total` even outside the Python bindings.

### Working notebooks ✓ done

`notebooks/01_load_and_inspect.ipynb` and `02_short_md.ipynb` use the real
`from_files → nvt → run` path against reference systems, including the
`Topology.solvate()` + `System(topo, conf)` path for systems `from_files()`
can't load directly (unsolvated topology vs. solvated configuration).

### Constraints and energy minimization ✓ done

Two gaps found while writing those notebooks, both closed:

- **Constraints on the factories.** `InputParameters.nve/nvt/npt` take
  `constraints="none"|"hbonds"|"allbonds"` (default `"none"`, matching GROMOS's
  own `NTC=1` default). Previously factory-built params had no way to turn SHAKE
  on, so a constrained system (e.g. one with solute H-bonds) would silently run
  unconstrained and diverge. `02_short_md.ipynb` §2a demonstrates the failure
  mode and the fix side by side, live.
- **Energy minimization via `Simulation`.** `InputParameters.steepest_descent(steps)`
  now actually minimizes when run through `Simulation` — previously it silently
  fell back to plain leap-frog with `dt=0` (a no-op). Today: `Recipe.minimize(steps)`;
  the plan it builds is `recipe.plan(system)`.

### FUTURE — System builder algebra

The full design is in `system_builder.py`'s module docs. Key idea: replace the traditional 8-binary
gromos++ pipeline with a composable Python API:

```python
ff     = ForceField.load("54A7")
system = molecule("ALA", ff) * 10 + solvent("SPC", n=2000)
system.neutralize(ion="CL")
```

Open decisions (recorded in `system_builder.py` comments) include:
native topology building vs. subprocess `make_top`, the `.solvate()` / `.pack()`
coordinate pipeline, and the `+` / `*` assembly algebra.

### FUTURE — SoA memory layout

The Rust core stores atoms as `Vec<Vec3>` (array-of-structs). A migration to
structure-of-arrays (`struct Soa { x: Vec<f64>, y: Vec<f64>, z: Vec<f64> }`)
would allow true zero-copy NumPy sharing and cleaner SIMD in the nonbonded loop.

## References

- GROMOS force fields and theory: [www.gromos.net](https://www.gromos.net)
- PyO3 (Rust-Python bindings): [pyo3.rs](https://pyo3.rs)
- Maturin (build tool): [github.com/PyO3/maturin](https://github.com/PyO3/maturin)
