# Quick Start

## Installation

```bash
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromos-rs
make build-python        # creates .venv, compiles Rust, installs gromos
source .venv/bin/activate
```

See [Installation](installation.md) for details and troubleshooting.

---

## Loading a system

A `System` is a topology + coordinate file pair.  
The atom count is validated at construction — mismatched files raise `ValueError` immediately.

```python
from gromos import System

system = System.from_files("water_216.topo", "equilibrated.cnf")
print(system)
# System(n_atoms=648, charge=0, box=(3.105, 3.105, 3.105))

print(system.n_atoms)           # 648
print(system.charge)            # 0  — integer, units of e
print(system.positions.shape)   # (648, 3) — float64, nm
print(system.velocities.shape)  # (648, 3) — float64, nm/ps
print(system.box)               # (3.105, 3.105, 3.105) nm
```

You can also build a `System` from pre-loaded objects:

```python
from gromos import Topology, Configuration, System

topo   = Topology("water_216.topo")
conf   = Configuration("equilibrated.cnf")
system = System(topo, conf)   # same validation
```

### Inspecting topology and configuration directly

```python
topo = Topology("water_216.topo")
print(topo.n_atoms)          # total atoms (solute + solvent)
print(topo.n_solute_atoms)   # solute-only
print(topo.n_solvent_atoms)  # solvent-only
print(topo.masses.shape)     # (N,) float64, g/mol
print(topo.charges.shape)    # (N,) float64, e

conf = Configuration("equilibrated.cnf")
print(conf.positions.shape)   # (N, 3) float64, nm
print(conf.velocities.shape)  # (N, 3) float64, nm/ps
print(conf.box_dimensions)    # (Lx, Ly, Lz) nm
```

---

## Creating a recipe

A `Recipe` is the one description of a run — the data of a GROMOS `.imd` grouped by
concern (`recipe.control`, `.boundary`, `.forcefield`, `.constraints`, `.ensemble`,
`.outputs`, …), plus additive `terms` and auxiliary `inputs`. Use a factory — no `.imd`
file authoring required.

```python
from gromos import Recipe

# NVT — Berendsen thermostat, τ = 0.1 ps
recipe = Recipe.nvt(dt=0.002, steps=5000, temperature=300.0)

# NVE — microcanonical, no thermostat
recipe = Recipe.nve(dt=0.002, steps=5000)

# NPT — Berendsen thermostat + barostat, water compressibility by default
recipe = Recipe.npt(dt=0.002, steps=5000, temperature=300.0, pressure=1.0)

# Energy minimisation — steepest descent
recipe = Recipe.minimize(steps=500)

# Or load an existing GROMOS input file. A block the engine does not model is an
# error (pass allow_passthrough=[...] to accept it); an absent optional block is
# reported in recipe.diagnostics with what gromosXX does without it.
recipe = Recipe.from_imd("run.imd")

print(recipe.control["dt"], recipe.control["steps"])
print(recipe.ensemble)           # {"thermostat": {...} | None, "barostat": {...} | None}
print(recipe.diagnostics)
```

A recipe is immutable: `update(**groups)` deep-merges and returns a new one, and a
typo is a `RecipeError`, never a silent default.

```python
longer = recipe.update(control={"steps": 20000})     # everything else unchanged
recipe.update(control={"stepz": 1})                  # RecipeError: unknown field `stepz`

recipe.to_toml()                                     # the run, as text
recipe.to_imd(n_atoms=system.n_atoms)                # what gromosXX would run
recipe.save_imd("run.imd", n_atoms=system.n_atoms)
Recipe.from_toml(recipe.to_toml()) == recipe         # True — lossless
```

### Constraints (SHAKE)

`nve`/`nvt`/`npt` take `constraints="none"|"hbonds"|"allbonds"` (default `"none"`,
matching GROMOS's own `NTC=1` default). A system with flexible solute H-bonds (a
peptide, for example) run with `constraints="none"` at a normal MD timestep will
diverge — this is a real trap, not a hypothetical one:

```python
# Diverges within ~50-100 steps: solute H-bonds are unconstrained.
recipe = Recipe.nvt(dt=0.002, steps=5000, temperature=300.0)

# Stable: SHAKE constrains solute H-bonds.
recipe = Recipe.nvt(dt=0.002, steps=5000, temperature=300.0, constraints="hbonds")
print(recipe.constraints)   # {"solute": "hbonds", "solute_algorithm": "shake", ...}
```

This only sets solute constraints (`NTC`). Solvent rigidity (SETTLE) isn't yet
exposed on the factories — `recipe.update(constraints={...})` or load a validated
`.in` file via `Recipe.from_imd()` if you need it.

---

## Running a simulation

```python
from gromos import Simulation

# Recommended: two-argument form
sim = Simulation(system, params)

# Legacy three-argument forms also work
sim = Simulation(topo, conf, params)
sim = Simulation("water_216.topo", "equilibrated.cnf", "run.imd")
```

### Advancing and reading state

```python
sim.step(100)   # advance 100 MD steps

# Energies (kJ/mol)
print(sim.total_energy)
print(sim.kinetic_energy)
print(sim.potential_energy)

# Positions and velocities ((N, 3) float64 NumPy arrays)
print(sim.positions.shape)   # nm
print(sim.velocities.shape)  # nm/ps
print(sim.forces.shape)      # kJ/(mol·nm)

# Clock
print(sim.time)          # ps
print(sim.current_step)  # step number

# Thermostat temperature
print(sim.temperature)   # K
```

### Collecting energies over a run

`sim.run(steps, ene_freq)` streams energies as one `(n_frames, 13)` NumPy array —
no Python-side per-step loop, no `.tre` file. Wrap it in `EnergyTimeseries` for
named-column access, block-averaging, dataframes, and plots.

```python
from gromos import EnergyTimeseries

energies = sim.run(5000, ene_freq=100)   # (51, 13) array
ts = EnergyTimeseries(energies)

print(f"Mean: {ts.total.mean():.1f} kJ/mol")
ts.plot("total", "kinetic", "potential")       # plotly figure (default backend)
ts.plot("bond", "angle", backend="matplotlib")  # or matplotlib

mean, err = ts.block_average("total", block_size=10)
print(f"Block-averaged: {mean:.1f} ± {err:.1f} kJ/mol")

df = ts.to_dataframe()   # polars.DataFrame by default
```

Defaults (`"polars"` for dataframes, `"plotly"` for plots) are configurable
globally via `gromos.timeseries.config`.

### Energy minimization

`Recipe.minimize(steps)` runs real minimization through `Simulation` — the plan
swaps in `steepest_descent` for the leap-frog integrator and has no
thermostat/barostat/kinetic energy (GROMOS convention: `total_energy ==
potential_energy` during EM).

```python
recipe = Recipe.minimize(steps=500)
sim = Simulation(system, recipe)
energies = sim.run(500, ene_freq=10)
print(sim.plan.kinds)
# ['forcefield', 'steepest_descent', 'energy_calculation']
```

The minimizer converges once the energy change between steps drops below its
tolerance and becomes a no-op after that — a flat tail at the end of the energy
trace is expected, not a bug.

---

## Inspecting and editing the plan

Every `Simulation` runs a fixed, ordered list of algorithms each step — the
*plan*. `recipe.plan(system)` returns it as data (one fully resolved `Algorithm`
per entry, already validated against the GROMOS ordering rules); edit it and hand
it back with `plan=`.

```python
# What is running
print(sim.plan.kinds)
# ['remove_com', 'forcefield', 'leap_frog_velocity', 'thermostat',
#  'leap_frog_position', 'temperature_calculation', 'energy_calculation']
print(sim.plan["forcefield"])          # every parameter, resolved

# Edit before constructing the simulation
plan = recipe.plan(system)
plan.remove("remove_com")              # disable COM motion removal
print("remove_com" in plan, len(plan))

sim = Simulation(system, recipe, plan=plan)   # re-validated: a broken order is a PlanError
sim.step(100)
```

`gromos.algorithms()` lists every kind this build knows with its parameters and
ordering rules; `gromos.terms()` lists the additive terms.

## Additive terms: QM and ML potentials

A `Term` is an energy/force provider added on top of the classical force field
(`coupling="delta"`) over a region of the system. `Term("xtb", ...)` runs a real
GFN-xTB subprocess (needs `xtb` on PATH); `Term("schnet", ...)` a TorchScript
SchNetPack model (needs a `--features ml` build — `Term(...).available` says
whether this build has it).

```python
from gromos import Term

qm = Term(
    "xtb",
    region="1:a",                 # molecule 1, all atoms (gromos-rs atom-specifier syntax)
    elements=[8, 1, 1, 8, 1, 1],  # atomic number per global atom index
    gfn=2, charge=0, multiplicity=1,
    timeout_s=600,                # every xtb call is bounded; expiry is a RunError, not a hang
)
sim = Simulation(system, recipe.with_term(qm))
sim.step(10)
print(sim.term_energies)          # {"xtb": -13307.5}  — each term on its own, kJ/mol
print(sim.total_energy)           # classical potential + terms + kinetic
```

Rules the plan enforces before anything runs: `coupling="replace"` is not available
yet (`PlanError` naming 2.8), a term without a virial cannot meet a barostat,
and two terms of one kind report as `xtb:0` / `xtb:1` (each with its own work directory).

## Free-energy perturbation from Python

A recipe with a perturbation topology in its `inputs` and `perturbation.enabled = 1` runs FEP;
`sim.dhdl` is the thermodynamic-integration integrand ∂H/∂λ at the current step — the whole
derivative: kinetic (λ-dependent masses), potential, perturbed constraints and restraints —
`sim.dhdl_terms` its `lj` / `crf` / `bonded` / `kinetic` / `constraints` / `distanceres` split,
and `run()`'s last column is `dhdl`:

```python
recipe = Recipe.from_imd("run.imd").with_inputs(pttopo="solute_dummy.ptp")
means = {}
for lam in (0.0, 0.25, 0.5, 0.75, 1.0):
    sim = Simulation(system, recipe.update(perturbation={"lambda": lam}))
    frames = sim.run(5000, ene_freq=10)                  # (501, 13)
    means[lam] = EnergyTimeseries(frames).dhdl[100:].mean()   # drop the first 1 ps
# ΔG = ∫ <∂H/∂λ> dλ — numpy.trapezoid(list(means.values()), list(means))
```

The `md` binary writes the same quantity, as the first value of each frame, to a `.trg` file
(`@trg`, every `NTWG` steps); `gromos.analysis.ti_trapezoid(trg_files, lambdas)` integrates a set
of windows. With a `PRECALCLAM` block the `.tre`/`.trg` also carry the energies on a λ grid, which
`gromos.analysis.ext_ti_ana` (gromos++'s program) reweights and writes for `gromos.analysis.bar`
and `gromos.analysis.mbar`. `scripts/ti_ch4.py` is a complete two-engine study with all three
estimators (`--estimator ti|bar|mbar`).

## Migrating from `InputParameters`

The pre-recipe forms still work for one release and emit a `DeprecationWarning`
naming the replacement. They are translations into a recipe, not a second code
path — `tests/test_front_end_parity.py` checks the two are bit-identical.

| Before (deprecated) | Now |
|---|---|
| `InputParameters("run.imd")`, `InputParameters.from_file(...)` | `Recipe.from_imd("run.imd")` |
| `InputParameters.nve/nvt/npt(...)` | `Recipe.nve/nvt/npt(...)` — same arguments |
| `InputParameters.steepest_descent(steps)` | `Recipe.minimize(steps)` |
| `params.dt`, `params.nstlim` | `recipe.control["dt"]`, `recipe.control["steps"]` |
| `params.temperature`, `params.constraints` | `recipe.ensemble["thermostat"]`, `recipe.constraints["solute"]` |
| `Simulation(system, params)` | `Simulation(system, recipe)` |
| `Simulation(..., distrest=, posresspec=, refpos=)` | `Simulation(system, recipe.with_inputs(distrest=..., posresspec=..., refpos=...))` |
| `Simulation(..., ml_potential=p, ml_region=r, ml_buffer=b)` | `Simulation(system, recipe.with_term(Term("schnet", model=..., cutoff=..., elements=[...], region=r, buffer=b)))` |
| `AlgorithmSequence.nvt(topo, params)` … `Simulation.from_sequence(topo, conf, params, seq)` | `plan = recipe.plan(system)`; edit it; `Simulation(system, recipe, plan=plan)` |
| `sim.recipe_toml`, `sim.plan_json` | still there; `sim.recipe` / `sim.plan` are the objects |

---

## Writing output

```python
# Write current coordinates to a GROMOS .cnf file
system.write("output.cnf")
```

---

## Units

| Quantity | Unit |
|----------|------|
| Length | nm |
| Time | ps |
| Energy | kJ/mol |
| Force | kJ/(mol·nm) |
| Velocity | nm/ps |
| Temperature | K |
| Pressure | bar |
| Mass | g/mol |
| Charge | elementary charge (e) |

---

## Next steps

- **[API Reference](../api/reference.md)** — every class and method documented
- **[Installation](installation.md)** — build options, troubleshooting
