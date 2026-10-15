# Free-energy windows

A free-energy difference by thermodynamic integration (TI) needs the same perturbed run at
several values of λ, the coupling parameter between state A (λ = 0) and state B (λ = 1):

$$
\Delta G_{A \to B} = \int_0^1 \left\langle \frac{\partial H}{\partial \lambda} \right\rangle_\lambda d\lambda
$$

`gromos.free_energy` builds those runs from one recipe, runs them in process and integrates
the result. It adds no physics: each window is an ordinary `Simulation`, and dH/dλ is the
value the engine writes first in each `.trg` frame (`Simulation.dhdl`, the `dhdl` column of
`Simulation.run()`).

## The perturbation

The perturbation topology (`.ptp`) comes in as a recipe input, as `md @pttopo` takes it, and
the recipe's `perturbation` group must be on (`PERTURBATION NTG = 1`):

```python
from gromos import Recipe, System, Topology

system = System.from_files("ch4_spc.top", "ch4_spc.cnf")
base = Recipe.from_imd("ch4_water_fep.in").with_inputs(pttopo="ch4_spc_dummy.ptp")
```

`Topology.apply_perturbation(path)` merges a `.ptp` into a topology object instead and returns
the solute's net charge in state A and state B:

```python
net_a, net_b = Topology("ch4_spc.top").apply_perturbation("ch4_spc_dummy.ptp")  # (0.0, 0.0)
```

A perturbation that changes the net charge under a reaction field runs, but the `md` binary
warns: the ΔG of a charge change carries a finite-size artefact and needs a correction
(`gromos.analysis.dgslv_pbsolv`). `pt_top` and `make_pt_top` print each perturbed molecule's
net charge in both states as well.

## Windows

```python
from gromos.free_energy import LambdaSchedule, run_windows, collect

recipes = LambdaSchedule([0.0, 0.25, 0.5, 0.75, 1.0]).recipes(base)
windows = run_windows(system, recipes, workdir="ti_ch4")
dg, error = windows.ti(skip=2)       # kJ/mol
print(windows.report(skip=2))
```

- `LambdaSchedule(points, interaction_lambdas=None)` sets `perturbation.lambda` of each window;
  `interaction_lambdas` are `LAMBDAS` rows (`{"kind": "vdw", "g1": 1, "g2": 2, "a": …, "e": …}`),
  the same individual λ_i(λ) in every window.
- `run_windows(system, recipes, workdir=None, dry_run=False)` runs the windows one after the
  other. A window samples dH/dλ every `outputs["ntwg"]` steps at the steps a `.trg` of the same
  run holds (0, NTWG, … below `control["steps"]`), so its series is the `.trg`'s first column.
  With a `workdir`, each window's energies are kept as `<workdir>/l0.250.npz`, and a window
  whose file exists is read instead of run again; `dry_run=True` returns the job list without
  running anything. `collect(workdir)` reads the kept windows back.
- `FreeEnergyWindows.ti(skip)` integrates ⟨∂H/∂λ⟩ with the trapezoid rule; its error comes
  from the block-averaging error of each window. It is the estimator of
  `gromos.analysis.ti_trapezoid`, which reads `.trg` files instead.

BAR and MBAR need the energies at the other λ values (`PRECALCLAM`), which the `md` binary
writes into its `.tre`/`.trg`: `gromos.analysis.ext_ti_ana`, `bar` and `mbar` read those files
(`scripts/ti_ch4.py --estimator bar|mbar`).

## Checked against gromosXX

`tests/test_free_energy.py` runs the CH4 → dummy windows at λ = 0, 0.25, 0.5, 0.75 and 1
(CH4 in 999 SPC, 10 steps) from Python: every dH/dλ sample agrees with gromosXX's `.trg` for the
four windows that have a reference (relative 1e-6), and `ti()` equals `analysis.ti_trapezoid`
on those `.trg` files.
