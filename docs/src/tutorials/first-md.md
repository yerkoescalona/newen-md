# Running Your First MD Simulation

Two argon atoms in vacuum, with no solvent and no heat bath: the smallest system that still
shows how an MD step works. In about 30 minutes you run it from Python, read one MD step as
data, break that step on purpose, measure how the energy error depends on the time step, and
read the Lennard-Jones well depth off a run.

The runnable version is the notebook
[`py-newen/notebooks/03_two_argon_atoms.ipynb`](../../../py-newen/notebooks/03_two_argon_atoms.ipynb).
This page follows its four sections and gives the numbers you should see. The notebook also
asks you to predict each result before you run it.

> newen is an educational reimplementation of GROMOS. It is not intended for research
> work; see the repository README.

## Before you start

From a checkout of the repository:

```sh
cd py-newen
uv sync --group notebooks      # Python dependencies and Jupyter; builds the engine
uv run --group notebooks jupyter lab notebooks/03_two_argon_atoms.ipynb
```

The first `uv sync` compiles the Rust extension, which takes a minute or more. The input files
come from `crates/newen-md/tests/gromosXX_references/pair_lj/`, a reference system that the
test suite checks against gromosXX.

## 1. The system and the recipe

A **system** is what is simulated: a *topology* (what the atoms are and how they interact)
plus a *configuration* (where they are, how fast they move, and the box). A **recipe** says
how to run it: time step, number of steps, ensemble, cut-offs.

```python
from newen import Recipe, Simulation, System

system = System.from_files("pair_lj.topo", "pair_lj.conf")  # two argon atoms, 0.35 nm apart
recipe = Recipe.from_imd("pair_lj.in")                      # NVE, dt = 0.002 ps
sim = Simulation(system, recipe)
```

Argon atoms interact through the Lennard-Jones potential, `V(r) = C12/r^12 - C6/r^6`, with
`C6 = 6.2647e-3 kJ mol⁻¹ nm⁶` and `C12 = 9.847e-6 kJ mol⁻¹ nm¹²` from the topology. At
0.35 nm that gives −0.493955 kJ/mol, and `sim.potential_energy` agrees to every printed digit.
The force is repulsive: the atoms start closer than the bottom of the well.

Creating a `Simulation` already takes step 0. The engine evaluates energies and forces at the
starting positions, then moves the atoms once, so `sim.potential_energy` and `sim.forces`
belong to the starting positions.

## 2. The plan: one MD step as data

The **plan** is the ordered list of algorithms that make up one MD step:

```python
plan = recipe.plan(system)
plan.kinds
# ['remove_com', 'forcefield', 'leap_frog_velocity', 'leap_frog_position',
#  'temperature_calculation', 'energy_calculation']
```

`forcefield` computes the force `F(t)` at the current positions `x(t)`. `leap_frog_velocity`
is the **kick**, `v(t + dt/2) = v(t - dt/2) + F(t)/m · dt`, and `leap_frog_position` is the
**drift**, `x(t + dt) = x(t) + v(t + dt/2) · dt`. Force, kick, drift is the **leap-frog
step**: velocities live at half steps and positions at whole steps.

The library refuses a plan that breaks this order, and says why:

| Edit | `plan.validate()` |
|---|---|
| drop `leap_frog_velocity` | `leap_frog_position requires leap_frog_velocity` |
| move `forcefield` after `leap_frog_position` | `leap_frog_velocity must come after forcefield` |
| put `leap_frog_position` before `leap_frog_velocity` | `leap_frog_velocity must come before leap_frog_position` |

The kick needs the force at positions the drift has not moved yet, and the drift needs the new
velocity. `Simulation(system, recipe, plan=...)` applies the same checks to an edited plan.

## 3. The time step and the energy error

With no heat bath the total energy should stay constant. Leap-frog keeps it constant only
approximately. Measure the error as the range of the total energy, largest minus smallest,
over 8 ps of simulated time:

```python
import numpy as np
from newen import EnergyTimeseries

def energy_range(system, recipe, dt, duration):
    n_steps = round(duration / dt)
    sim = Simulation(system, recipe.update(control={"dt": dt, "steps": n_steps}))
    return np.ptp(EnergyTimeseries(sim.run(n_steps, ene_freq=1)).total)
```

| dt (ps) | 0.0005 | 0.001 | 0.002 | 0.004 | 0.008 | 0.016 | 0.032 | 0.064 | 0.128 |
|---|---|---|---|---|---|---|---|---|---|
| range (kJ/mol) | 4.33e-6 | 1.73e-5 | 6.93e-5 | 2.77e-4 | 1.11e-3 | 4.47e-3 | 1.83e-2 | 8.20e-2 | 8.03e-1 |

Each doubling of `dt` multiplies the error by four. The log–log slope up to 0.032 ps is 2.01:
leap-frog is a **second-order** method, with an energy error that grows as `dt²`. At 0.128 ps
the points leave that line. The atoms swing apart and back once every 1.98 ps, so a step of
0.128 ps leaves only about 15 steps per swing, and the step is no longer small.

## 4. The well depth

The Lennard-Jones minimum lies at `r_min = (2·C12/C6)^(1/6)` = 0.3827 nm. Its depth, the
**well depth** `ε`, is `C6² / (4·C12)` = 0.99641 kJ/mol.

The lowest potential energy in an 8 ps run at `dt = 0.0005 ps` is −0.99641 kJ/mol, which is
`−ε`. The total energy, −0.494 kJ/mol, lies above the bottom of the well, so the atoms swing
between turning points at 0.350 nm and 0.470 nm. On every swing they pass `r_min`, where the
potential energy is `−ε` and the kinetic energy is largest.

## Where next

- `py-newen/notebooks/01_load_and_inspect.ipynb` loads larger systems: bulk water and a
  solvated alanine dipeptide.
- `py-newen/notebooks/02_short_md.ipynb` adds solvent, heat baths and constraints.
- [Pairlist & Neighbor Search](../algorithms/pairlist.md) explains how the engine finds the
  pairs that `forcefield` evaluates.
