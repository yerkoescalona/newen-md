# py-newen

Python bindings for the newen molecular dynamics engine.

## Five lines to a running simulation

```python
from newen import System, Recipe, Simulation

system = System.from_files("water_216.topo", "equilibrated.cnf")
recipe = Recipe.nvt(dt=0.002, steps=5000, temperature=300.0)   # or Recipe.from_imd("run.imd")
sim    = Simulation(system, recipe)
sim.step(1000)

print(sim.total_energy)     # kJ/mol
print(sim.positions.shape)  # (648, 3)
```

No `.imd` files to author. No temporary directories. All in memory.

## What works today

| Class / function | Description |
|-----------------|-------------|
| `System` | Paired topology + coordinates, atom-count validated |
| `Topology` | Load `.topo` — atoms, masses, charges, bonds; `.solvate(nsm)` |
| `Configuration` | Load `.cnf` / `.g96` — positions, velocities, box |
| `Recipe` | The run as data: load `.imd` (`from_imd`) or build via `nve / nvt / npt / minimize` factories (with a `constraints=` knob for SHAKE); immutable `update(**groups)`; `to_imd` / `to_toml` / `to_bundle`; `terms` (`Term("schnet", ...)`) and `inputs` (restraint files) |
| `Simulation` | `Simulation(system, recipe)` + `step(n)` / `run(steps, ene_freq)` + full energy/position/force access; `sim.recipe`, `sim.plan`, `sim.diagnostics` |
| `Plan` | `recipe.plan(system)`: the MD step as an editable, validated list of `Algorithm`s (`remove("remove_com")`, `insert_after(...)`); run it with `Simulation(system, recipe, plan=plan)` |
| `newen.terms()`, `newen.algorithms()`, `newen.build_info()` | What this build knows: term and algorithm kinds, their parameters and ordering rules, cargo features |
| `EnergyTimeseries` | Wraps `sim.run()`'s output — `block_average()`, `to_dataframe()` (polars/pandas/dict), `plot()` (plotly/matplotlib) |
| `Vec3`, `Energy`, `Frame` | Utility types |
| `rmsd`, `rdf` | Analysis functions |

The engine is validated bit-for-bit against double-precision gromosXX output
(energy rel. tolerance 1 × 10⁻⁸, force abs. tolerance 1 × 10⁻⁶ kJ mol⁻¹ nm⁻¹)
across 21 reference systems.

## What is coming next

| Milestone | Description |
|-----------|-------------|
| **P3 remaining** | ForceField single-point evaluation; rich `_repr_html_` for Jupyter |
| **FUTURE** | System builder algebra: `molecule("ALA", ff) * 10 + solvent("SPC", n=2000)` |

## Architecture

```
Python user code
      │
      ▼
newen/__init__.py   — clean re-export of working names only
      │
      ▼
newen/newen.abi3.so   — PyO3 extension compiled from pyo3-newen (Rust)
      │
      ▼
newen workspace
  newen-core · newen-io · newen-forces · newen-integrators
```

Positions, velocities, and forces are returned as `(N, 3)` float64 NumPy arrays.

## Documentation

- **[Quick Start](user-guide/quick-start.md)** — first simulation step by step
- **[API Reference](api/reference.md)** — every class and method
- **[Installation](user-guide/installation.md)** — build from source

## License

GPL-2.0 — same as GROMOS.
