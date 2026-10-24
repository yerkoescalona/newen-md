# py-newen

Python bindings for the newen molecular dynamics engine.

⚠️ **Educational / research project, early alpha.** Bit-for-bit validated against double-precision
gromosXX on 40 reference systems (`tests/test_gromosXX_references.py`), not a production MD
package. See the repository README for the licensing and provenance disclaimer.

## Five lines to a running simulation

```python
from newen import System, Recipe, Simulation

system = System.from_files("water_216.topo", "equilibrated.cnf")
recipe = Recipe.nvt(dt=0.002, steps=5000, temperature=300.0)   # or Recipe.from_imd("run.imd")
sim    = Simulation(system, recipe)
sim.step(1000)
print(sim.total_energy, sim.positions.shape)                     # kJ/mol, (648, 3)
```

No `.imd` files to author, no temporary directories, no subprocess: the engine runs in-process
and every quantity is a NumPy array or a float.

## The model

| Object | What it is |
|---|---|
| `System` | topology + coordinates (`.topo` + `.cnf`), atom counts checked |
| `Recipe` | the run as data — the content of a GROMOS `.imd` grouped by concern (`.control`, `.ensemble`, `.constraints`, …), plus additive `terms` and auxiliary `inputs`. Immutable; `update(**groups)` / `with_*` return copies; a typo is a `RecipeError`. `to_imd()` is what gromosXX would run |
| `Plan` | `recipe.plan(system)`: the MD step as an ordered, validated list of `Algorithm`s; edit it and run it with `Simulation(system, recipe, plan=plan)` |
| `Term` | an additive energy/force provider over a region: `Term("xtb", ...)` (real GFN-xTB, needs `xtb` on PATH), `Term("schnet", ...)` (TorchScript SchNetPack, needs a `--features ml` build) |
| `Simulation` | `step(n)`, `run(steps, ene_freq)` → `(n_frames, 12)` energies; `positions`, `velocities`, `forces`, `energies`, `term_energies`, `temperature`, `pressure`, `recipe`, `plan`, `diagnostics` |
| `EnergyTimeseries` | wraps `run()`'s array: named columns, block averages, dataframes (polars/pandas), plots (plotly/matplotlib) |

`newen.terms()` / `newen.algorithms()` / `newen.build_info()` say what this build knows.
`newen.exceptions` has one class per failure kind (`RecipeError`, `PlanError`,
`MissingFeatureError`, `RunError`), each a subclass of the builtin raised before.

The `md` binary and `Simulation` are the same code path: both hand a recipe to the `newen-run`
crate (`prepare_system → build_plan → validate_plan → instantiate → start`). The Rust reference
suite drives the binary, this package's suite drives the binding, and
`tests/test_front_end_parity.py` checks every front-end against every other with exact equality.

Deprecated for one release (they warn and are translated into a recipe): `InputParameters`, the
`distrest=`/`posresspec=`/`refpos=`/`ml_*=` keywords, `AlgorithmSequence.*` presets and
`Simulation.from_sequence`. Migration table: `docs/user-guide/quick-start.md`.

## Installation

Prerequisites: Python ≥ 3.13 (see `pyproject.toml`), a Rust toolchain (https://rustup.rs), `uv`.
The full dependency list (QM engines, the ML stack) is in the repository's `INSTALL.md`.

```bash
cd py-newen
uv sync                          # creates .venv with the dev dependencies
uv run maturin develop --release # builds crates/pyo3-newen into the venv
uv run pytest tests/ -q          # 306 passed / 8 skipped / 3 xfailed at 0.0.32
```

`maturin` is in the `dev` group, so `uv run maturin` uses the venv's own copy and a later
`uv sync --all-groups` (what `make sync` runs) keeps it; `uv sync --extra education` adds the
viewer and plotting dependencies of `newen.education`.

Optional: `uv sync --group ml` + `uv run maturin develop --release --features ml` for
`Term("schnet", ...)` (libtorch via the venv's torch — `INSTALL.md` §5).

## Documentation

- [User guide](docs/user-guide/): [installation](docs/user-guide/installation.md),
  [quick start](docs/user-guide/quick-start.md) (with the migration table),
  [the recipe model](docs/user-guide/recipe.md), [learning guide](docs/user-guide/learning-guide.md)
- [API reference](docs/api/reference.md)
- [Development](docs/development/): building, contributing, a PyO3 bindings tutorial
- Notebooks: `notebooks/03_two_argon_atoms.ipynb` (start here: a first MD run, about 30 minutes),
  `notebooks/01_load_and_inspect.ipynb`, `notebooks/02_short_md.ipynb`

`uv run mkdocs serve` in this directory renders the site at http://localhost:8000.

## Testing

```bash
uv run pytest tests/ -q                                   # Python suite (needs the built extension)
MYPYPATH=python uv run python -m mypy.stubtest newen.newen \
    --allowlist stubtest_allowlist_no_ml.txt              # the .pyi against the extension
cargo test -p newen-run -p pyo3-newen                   # the crates behind the binding
```

`tests/test_xtb_term.py` (the QM term's physics oracle) and `tests/test_qm_vs_ml_comparison.py`
skip without `xtb` on PATH; `tests/test_ml_potential.py` skips without a `--features ml` build.

## Layout

```
python/newen/          the package: __init__.py (the module contract), newen.pyi (stubs),
                        exceptions.py, timeseries.py, analysis.py (explicit import only),
                        system_builder.py (design sketch)
tests/                  test_gromosXX_references.py, test_front_end_parity.py, test_recipe.py,
                        test_xtb_term.py, test_simulation*.py, ...
docs/                   mkdocs site        notebooks/   Jupyter notebooks
../crates/pyo3-newen   the Rust binding   ../crates/newen-run   the run builder both entry points use
```

## License

GPL-2.0 — see the repository LICENSE.
