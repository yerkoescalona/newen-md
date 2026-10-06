# API Reference

All classes and functions exported by `import gromos`.

---

## System

Pairs a `Topology` with a `Configuration`. Validates that both describe the same
number of atoms at construction time.

```python
from gromos import System
```

### Constructors

```python
System(topology: Topology, configuration: Configuration) -> System
```
Raises `ValueError` if `topology.n_atoms != configuration.n_atoms`.

```python
System.from_files(topo_file: str, conf_file: str) -> System  # staticmethod
```
Reads both files and validates in one call.

**Example**
```python
system = System.from_files("water_216.topo", "equilibrated.cnf")

topo   = Topology("water_216.topo")
conf   = Configuration("equilibrated.cnf")
system = System(topo, conf)   # identical result
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_atoms` | `int` | Total atom count |
| `charge` | `int` | Net integer charge (e) |
| `positions` | `ndarray (N,3) f64` | Positions, nm |
| `velocities` | `ndarray (N,3) f64` | Velocities, nm/ps |
| `box` | `tuple[float,float,float]` | Box dimensions (Lx, Ly, Lz), nm |
| `topology` | `Topology` | Underlying topology |
| `configuration` | `Configuration` | Underlying configuration |

### Methods

**`write(path: str)`**  
Write current coordinates to a GROMOS `.cnf` file.

---

## Topology

Reads a GROMOS topology file (`.topo`).

```python
from gromos import Topology
topo = Topology("system.topo")
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_atoms` | `int` | Total atom count (solute + solvent) |
| `n_solute_atoms` | `int` | Solute-only count |
| `n_solvent_atoms` | `int` | Solvent-only count |
| `masses` | `ndarray (N,) f64` | Atom masses, g/mol |
| `charges` | `ndarray (N,) f64` | Partial charges, e |

### Methods

**`solvate(nsm: int)`**  
Append `nsm` copies of the solvent template to the topology in-place.

---

## Configuration

Reads a GROMOS coordinate file (`.cnf` / `.g96`).

```python
from gromos import Configuration
conf = Configuration("initial.cnf")
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_atoms` | `int` | Atom count in this file |
| `positions` | `ndarray (N,3) f64` | Positions, nm |
| `velocities` | `ndarray (N,3) f64` | Velocities, nm/ps |
| `box_dimensions` | `tuple[float,float,float]` | (Lx, Ly, Lz), nm |

---

## Recipe

The one description of a run — the data of a GROMOS `.imd` grouped by concern, plus
additive `terms` and auxiliary `inputs`. Immutable: `update`/`with_*` return a new recipe.
`Recipe`, `Term`, `Algorithm` and `Plan` are the engine's own (`gromos-run`) types, so a
typo'd field or an unknown kind is a `RecipeError`, never a silent default.

```python
from gromos import Recipe, Term, Plan, Algorithm, terms, algorithms, build_info
```

### Constructors

```python
Recipe.from_imd(path: str, allow_passthrough: list[str] | None = None)   # a GROMOS .imd
Recipe.from_toml(text: str) / Recipe.from_json(text: str) / Recipe.from_dict(d: dict)
Recipe.from_bundle(path: str, allow_passthrough=None)                    # a run bundle's input.toml
Recipe(**groups)                                                         # defaults + groups
```

An `.imd` block the engine does not model is an error unless named in `allow_passthrough`;
an absent optional block is reported in `recipe.diagnostics` together with what gromosXX does
without it.

### Factory methods (staticmethod)

```python
Recipe.nve(dt: float, steps: int, constraints: str = "none") -> Recipe
Recipe.nvt(dt: float, steps: int, temperature: float, constraints: str = "none") -> Recipe
Recipe.npt(dt: float, steps: int, temperature: float, pressure: float, constraints: str = "none") -> Recipe
Recipe.minimize(steps: int) -> Recipe
```
`nve`: a bath with coupling time −1 (gromosXX's "no coupling"). `nvt`: Berendsen thermostat,
τ = 0.1 ps. `npt`: Berendsen thermostat + barostat, `pressure` in bar, compressibility
4.575 × 10⁻⁴ nm² kJ⁻¹ mol (water). `constraints`: `"none"` | `"hbonds"` | `"allbonds"`
(solute SHAKE, GROMOS `NTC`).

### Groups (read-only dicts) and derived views

`recipe.version`, `.title`, `.system`, `.control`, `.boundary`, `.forcefield`, `.constraints`,
`.ensemble`, `.minimisation`, `.perturbation`, `.outputs`, `.execution`, `.inputs`,
`.passthrough`, `.terms` (`list[Term]`), `.diagnostics` (`list[str]`), `is_minimization()`.

### Building a new recipe

```python
recipe.update(**groups) -> Recipe            # deep merge; unknown group/field -> RecipeError
recipe.with_term(term: Term) -> Recipe       # add an additive term (e.g. Term("schnet", ...))
recipe.without_terms() -> Recipe
recipe.with_inputs(pttopo=None, posresspec=None, refpos=None, distrest=None) -> Recipe
recipe.with_execution(parallel: str) -> Recipe   # "auto" | "serial" | "parallel"
```

### Serialisation

```python
recipe.to_dict() / .to_toml() / .to_json()
recipe.to_imd(n_atoms: int | None = None) -> str          # what gromosXX would run
recipe.save_imd(path: str, n_atoms=None)
recipe.to_bundle(directory: str, system: System, topology_path: str, configuration_path: str) -> str
```

### The plan

```python
recipe.plan(system: System) -> Plan
```
Stage 1 of the builder: the MD step as an ordered, validated list of `Algorithm`s. Address
entries by index or kind: `plan["forcefield"]`, `plan[-1]`, `"remove_com" in plan`, `len(plan)`,
`plan.kinds`; edit with `insert(i, alg)`, `insert_after(target, alg)`, `insert_before(target,
alg)`, `remove(target)`, `replace(target, alg)`; `plan.validate()` raises `PlanError` on a broken
GROMOS step order; `to_json()` / `Plan.from_json(text)` / `to_dicts()`. Run it with
`Simulation(system, recipe, plan=plan)` (re-validated).

### Terms

```python
Term("xtb", region: str, elements: list[int], gfn: int = 2, charge: int = 0, multiplicity: int = 1,
     work_dir: str | None = None, timeout_s: int = 600, coupling: str = "delta")   # real xtb, no feature
Term("schnet", model: str, cutoff: float, elements: list[int], region: str,
     buffer: str | None = None, coupling: str = "delta")                          # --features ml
```
`elements` is one atomic number per global atom index (it must cover the region). A term is
additive over the classical force field; its energy is reported under its kind in
`Simulation.term_energies`.

### Registries

```python
terms()      -> [{"kind", "params", "feature", "available"}, ...]
algorithms() -> [{"kind", "params", "rules"}, ...]      # rules: unique/required/first/last/after/before/...
build_info() -> {"version", "recipe_version", "features"}
```

### Exceptions (`gromos.exceptions`)

`RecipeError(ValueError)`, `PlanError(ValueError)`, `MissingFeatureError(RuntimeError)` — a term
needs a cargo feature this build lacks (`Term("schnet", ...)` on a non-`ml` build) —
`RunError(RuntimeError)`.

### Deprecated: `InputParameters`

`InputParameters(path)` / `.from_file` / `.nve` / `.nvt` / `.npt` / `.steepest_descent` still
work for one release, warn, and are translated into a recipe when passed to `Simulation`.
See the migration table in the [Quick Start](../user-guide/quick-start.md).

---

## Simulation

Prepares the system, builds (or takes) the plan of the recipe, instantiates it and runs the
MD loop — through the same `gromos-run` code the `md` binary uses.

```python
from gromos import Simulation
```

### Constructors

```python
# Recommended
Simulation(system: System, recipe: Recipe)
Simulation(system: System, recipe: Recipe, plan=plan)      # an edited recipe.plan(system)

# Three-argument forms
Simulation(topo: Topology, conf: Configuration, recipe: Recipe)
Simulation(topo_file: str, conf_file: str, input_file: str)
Simulation.from_files(topo_file: str, conf_file: str, input_file: str)
Simulation.from_bundle(path: str, allow_passthrough=None)   # everything from a run bundle

# Deprecated (one release; warn and are translated into a recipe)
Simulation(system, params: InputParameters)
Simulation(..., distrest=, posresspec=, refpos=)            # -> recipe.with_inputs(...)
Simulation(..., ml_potential=, ml_region=, ml_buffer=)      # -> recipe.with_term(Term("schnet", ...))
Simulation.from_sequence(topo, conf, params, sequence: Plan)
```

### Running

**`step(n_steps: int)`**  
Advance the simulation by `n_steps` steps. All state properties are updated after
this call.

**`run(n_steps: int, ene_freq: int = 100) -> ndarray (n_frames, 13) f64`**  
Advance `n_steps` steps in Rust, sampling energies every `ene_freq`-th step, and
return one `(n_frames, 13)` NumPy array — no Python-side per-step loop, no `.tre`
file. Columns: `[time, kinetic, potential, total, volume, pressure, bond, angle,
improper, dihedral, lj, coulomb, dhdl]` (`dhdl` is dH/dλ, 0 without a perturbation). Frame 0 is the state before any of these steps
ran. Wrap the result in `EnergyTimeseries` (see below) for named-column access,
block-averaging, dataframes, and plots.

```python
energies = sim.run(5000, ene_freq=100)   # (51, 13) array
from gromos import EnergyTimeseries
ts = EnergyTimeseries(energies)
ts.total          # ndarray — total energy column
ts.plot("bond", "angle")
```

### State properties (read-only)

| Property | Type | Description |
|----------|------|-------------|
| `total_energy` | `float` | Total energy, kJ/mol |
| `kinetic_energy` | `float` | Kinetic energy, kJ/mol |
| `potential_energy` | `float` | Potential energy, kJ/mol |
| `temperature` | `float` | Instantaneous temperature, K (uses constraint-aware degrees of freedom, not bare `3*n_atoms`) |
| `volume` | `float` | Current box volume, nm³ |
| `pressure` | `float` | Instantaneous pressure, bar. **Only physically meaningful under NPT** — `P = (2*KE - virial)/(3*V)`, and the virial is only populated by `PressureCalculation`, which only NPT's sequence includes. Under NVE/NVT this returns a kinetic-only number, not zero — don't mistake it for a real value. |
| `positions` | `ndarray (N,3) f64` | Positions, nm |
| `velocities` | `ndarray (N,3) f64` | Velocities, nm/ps |
| `forces` | `ndarray (N,3) f64` | Forces, kJ/(mol·nm) |
| `time` | `float` | Current time, ps |
| `current_step` | `int` | Step counter |
| `dt` | `float` | Timestep, ps |
| `n_atoms` | `int` | Atom count |
| `n_solute_atoms` | `int` | Solute atom count |
| `n_solvent_atoms` | `int` | Solvent atom count |
| `algorithm_names` | `list[str]` | Names of algorithms in the sequence |
| `term_energies` | `dict[str, float]` | Energy of each additive term at the last step, keyed by registry name (`xtb`, `xtb:1` for a repeated kind); their sum is what the terms add to `total_energy` |
| `dhdl` | `float` | ∂H/∂λ at the last step (kJ/mol) — the TI integrand; 0 without a perturbation |
| `dhdl_terms` | `dict[str, float]` | ∂H/∂λ split: `lj`, `crf`, `bonded` |
| `recipe` | `Recipe` | The effective recipe the engine ran (`recipe_toml`: as TOML) |
| `plan` | `Plan` | The plan the engine instantiated, a frozen snapshot (`plan_json`: as JSON) |
| `diagnostics` | `list[str]` | Absent optional blocks (and what that means), passed-through blocks |
| `energies` | `Energy` | Full energy object, all components (bond/angle/dihedral/improper/lj/coulomb) |

### Energy minimization

`Recipe.minimize(steps)` runs real steepest-descent minimization through `Simulation` —
the plan swaps in `steepest_descent` for the leap-frog integrator and skips the
thermostat/barostat/temperature calculation (GROMOS convention: no velocities during EM,
`total_energy == potential_energy`). The minimizer converges (and becomes a no-op) once
the energy change between steps drops below its tolerance — a flat tail in `sim.run()`'s
output is expected, not a bug.

```python
recipe = Recipe.minimize(steps=500)
sim = Simulation(system, recipe)
energies = sim.run(500, ene_freq=10)
print(sim.plan.kinds)
# ['forcefield', 'steepest_descent', 'energy_calculation']
```

---

## AlgorithmSequence (deprecated)

The descriptor path (`Forcefield`, `LeapFrogIntegrator`, `BerendsenThermostat`, … building
blocks and a second builder) was removed in 3.9 step 4. The preset names survive one
release: each returns the `Plan` of `params` (what `Recipe.from_imd(...).plan(system)`
builds) and warns. `AlgorithmSequence` itself cannot be instantiated.

```python
AlgorithmSequence.nve / nvt / npt / minimize / from_parameters(topo: Topology, params: InputParameters) -> Plan
Simulation.from_sequence(topo, conf, params, sequence: Plan)   # deprecated
```

Today:

```python
plan = recipe.plan(system)
plan.remove("remove_com")
sim = Simulation(system, recipe, plan=plan)
```

---

## Vec3

A 3D vector. Useful for single-atom inspection; use NumPy arrays from `Simulation`
or `System` for bulk work.

```python
from gromos import Vec3
v = Vec3(1.0, 2.0, 3.0)
```

### Properties

`x`, `y`, `z` — `float` (read-only)

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `length()` | `float` | Euclidean norm |
| `normalize()` | `Vec3` | Unit vector |
| `dot(other: Vec3)` | `float` | Dot product |
| `cross(other: Vec3)` | `Vec3` | Cross product |

### Operators

`+`, `-`, `* scalar` between `Vec3` objects.

---

## Energy

Snapshot of energy components (kJ/mol).

```python
from gromos import Energy
e = Energy()
```

### Properties

`total`, `kinetic`, `potential`, `bond`, `angle`, `dihedral`, `improper`, `lj`,
`coulomb` — `float`

---

## EnergyTimeseries

Wraps the `(n_frames, 13)` NumPy array returned by `Simulation.run()`.

```python
from gromos import EnergyTimeseries

energies = sim.run(steps=5000, ene_freq=100)
ts = EnergyTimeseries(energies)
```

### Column access (by name)

`time`, `kinetic`, `potential`, `total`, `volume`, `pressure`, `bond`, `angle`,
`improper`, `dihedral`, `lj`, `coulomb` — each returns the `(n_frames,)` NumPy
column, e.g. `ts.total`, `ts.bond`.

### Methods

**`to_dataframe(backend: str | None = None)`**  
Returns a dataframe. `backend`: `"polars"` (default) | `"pandas"` | `"dict"`.
Default is set by `gromos.timeseries.config.dataframe_backend`.

**`plot(*components: str, backend: str | None = None)`**  
Plots one or more energy components against time. `backend`: `"plotly"`
(default) | `"matplotlib"`. Default is set by `gromos.timeseries.config.plot_backend`.
Returns the figure/axes object.

**`block_average(component: str, block_size: int) -> tuple[float, float]`**  
Block-averages `component`, returning `(mean, standard_error)`. Splits the
trajectory into non-overlapping blocks; the standard error is the standard
deviation of block means over `sqrt(n_blocks)`. Use to sanity-check
autocorrelation: the error should grow with block size until it plateaus.

**Example**
```python
ts.plot("bond", "angle")                       # plotly figure
ts.plot("total", backend="matplotlib")          # matplotlib axes
mean, err = ts.block_average("total", block_size=10)
df = ts.to_dataframe()                          # polars.DataFrame
```

### Global defaults

```python
import gromos.timeseries
gromos.timeseries.config.dataframe_backend = "pandas"   # or "polars" / "dict"
gromos.timeseries.config.plot_backend = "matplotlib"     # or "plotly"
```

---

## Frame

A trajectory frame (time + step metadata + optional positions).

```python
from gromos import Frame
f = Frame(time=1.5, step=100)
print(f.time, f.step, f.n_atoms)
```

---

## Free functions

### `rmsd`

Root-mean-square deviation between two coordinate arrays, nm.

```python
from gromos import rmsd
import numpy as np

r = rmsd(
    positions,   # (N, 3) float64
    reference,   # (N, 3) float64
)
```

### `rdf`

Radial distribution function between two atom index groups.

```python
from gromos import rdf

r_vals, g_vals = rdf(
    positions,                     # (N, 3) float64
    group1=[0, 3, 6],              # list[int] — indices in group 1
    group2=[1, 4, 7],              # list[int] — indices in group 2
    n_bins=100,
    r_max=1.5,                     # nm
)
# r_vals: (n_bins,) float64 — bin centres, nm
# g_vals: (n_bins,) float64 — g(r) values
```

---

## Planned API (not yet implemented)

### FUTURE — System builder algebra

Design rationale in `system_builder.py`. Nothing below is implemented.

```python
from gromos import ForceField

ff     = ForceField.load("54A7")
system = molecule("ALA", ff) * 10 + solvent("SPC", n=2000)
system.neutralize(ion="CL")
system.write("prepared.topo", "prepared.cnf")
```

---
