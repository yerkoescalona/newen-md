# Experiments with a check: `newen.education.experiments`

Each experiment runs a small simulation, draws what it measured over the theory it should
match, and returns the claim as a number with the band it must fall in. A notebook or a CI job
asserts it, so a change to the engine that breaks the physics fails loudly. The notebook
`notebooks/05_experiments.ipynb` runs all of them in a few seconds.

```python
from newen.education import experiments as ex

r = ex.energy_flow()
r.plot()                  # a matplotlib Figure: measurement over theory
r.check()                 # Check(total-energy range / kinetic-energy gain: measured 0.0226, ...: PASS)
assert all(r.checks)
r.data["total"]           # every measured series, as numpy arrays
```

`Result.check()` never raises; `Check.assert_passes()` raises an `AssertionError` whose message
is one line: what was measured, the value, what theory expects and the band.

## One core, no second codebase

The experiments compute no physics. Every number comes from the engine: the `Simulation`'s
positions, velocities, energies and temperature, `EnergyTimeseries.block_average`, and
`newen.rdf` (the gromos++ `rdf` port). The module arranges these and compares them with
`newen.education.theory`, a handful of textbook formulas with their own unit tests. A test
reads the module's source and fails if it gains a histogram or minimum-image arithmetic of its
own.

## The toy systems

They are built in memory as GROMOS topology and coordinate text, read by
`Topology.from_string` and `Configuration.from_string`; each returns a `System` and a `Recipe`
you can change with `recipe.update(...)`.

| Builder | What it is |
|---|---|
| `lj_gas(n, density, temperature, dim, seed, speeds)` | argon on a lattice in a periodic box, reduced units (σ 0.3405 nm, ε 0.9961 kJ/mol); `dim=2` keeps every atom in one plane and counts 2n − 2 degrees of freedom; `speeds="equal"` starts every atom at one speed with zero net momentum |
| `argon_pair(distance)` | two argon atoms at rest in vacuum |
| `diatomic(m1, m2, period_fs, r0, stretch)` | an O–H harmonic bond whose force constant gives the requested period, started stretched |

`with_bath(recipe, temperature_K, tau=0.5, algorithm="berendsen")` couples any of these recipes
to one heat bath: `"berendsen"` or `"nose_hoover"` (a chain of three).

## The experiments

| Function | What it shows | Check |
|---|---|---|
| `energy_flow()` | a 2D lattice melts: kinetic and potential energy trade, the total stays put | total-energy range under 5% of the kinetic-energy gain (the residual, about 2%, comes from pairs crossing the LJ cut-off and does not shrink with dt) |
| `time_step_scan()` | the energy error of leap-frog grows as dt² (argon pair, 8 ps per step size) | log–log slope in [1.9, 2.1] |
| `energy_landscape()` | the pair swings between turning points on V(r), drawn from the topology's C6/C12 | lowest potential energy within 0.1% of −ε |
| `thermostat_fingerprint()` | Berendsen and a Nosé–Hoover chain hold 120 K with different kinetic-energy fluctuations | Berendsen σ/mean under half of √(2/N_df); Nosé–Hoover within 15% of it |
| `block_error()` | correlated samples: the error of a mean grows with the block length | error at the longest block over twice the naive one; the curve rises, with at most one dip of at most 10% |
| `maxwell_boltzmann()` | equal starting speeds spread into the 2D Maxwell–Boltzmann law | late mean speed within 3% of √(πk_BT/2m) |
| `liquid_structure()` | a Lennard-Jones liquid has a shell of neighbours near σ and g(r) → 1 beyond | first peak in [1.0σ, 1.2σ] and above 1.5; mean g beyond 2.5σ in [0.9, 1.1] |
| `water_structure(topology, configuration, recipe)` | the same on a water box given as GROMOS files | O–O peak in [0.26, 0.29] nm; mean g beyond 0.6 nm in [0.9, 1.1] |
| `heat_bath()` | a hot 2D gas meets a colder bath: Berendsen relaxes without overshoot, a Nosé–Hoover chain overshoots | see below |

### A hot gas meets a cold bath

`heat_bath()` equilibrates a 64-atom 2D gas at 240 K, then quenches the same state three times
to a 120 K bath: Berendsen with τ = 0.2 ps and τ = 1 ps, and a Nosé–Hoover chain with τ = 1 ps,
sampling the energies every step. It claims three things:

- **The bath temperature is reached.** The mean temperature of the second half is within 3% of
  120 K for both Berendsen runs and within 5% for Nosé–Hoover.
- **What the gas loses, the bath takes.** Berendsen removes kinetic energy at the rate
  (N_df k_B / 2τ)(T − T₀), so the heat given to the bath is Q = (N_df k_B / 2τ) ∫ (T − T₀) dt;
  the drop of the gas's total energy matches Q within 5%.
- **Berendsen relaxes on a time scale set by τ.** The fitted time constant over τ lies in
  [0.9, 2.0], and the ratio of the two fitted constants is within [0.5, 1.6] of the ratio of
  the τ.

The fitted constant is longer than τ because the potential energy cools too: it is
τ · 2C_V/(N_df k_B) (Berendsen et al. 1984), about 1.2–1.4 τ here, and τ itself only for an
ideal gas. The left panel shows the three temperatures with the kinetic-only exponentials, the
right one the total-energy drop and Q on top of each other.

```python
heat_bath(temperature_K=120.0, start_K=240.0, taus=(0.2, 1.0), steps=8000,
          equilibration=1500, n=64, seed=7)
```

Every experiment takes its sizes as parameters (atoms, steps, frames), with seeds fixed, so the
same call gives the same numbers on every machine.

## Install

The plots need matplotlib, part of the `education` extra:

```bash
pip install "newen[education]"
```

Importing `newen.education` works without the extra; only the live viewer
([Watching a Run](education-view.md)) needs `anywidget`.
