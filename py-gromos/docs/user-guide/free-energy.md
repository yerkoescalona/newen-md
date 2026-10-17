# Free-energy windows

A free-energy difference by thermodynamic integration (TI), BAR or MBAR needs the same perturbed run at
several values of λ, the coupling parameter between state A (λ = 0) and state B (λ = 1):

$$
\Delta G_{A \to B} = \int_0^1 \left\langle \frac{\partial H}{\partial \lambda} \right\rangle_\lambda d\lambda
$$

`gromos.free_energy` builds those runs from one recipe, runs them in process and estimates
ΔG from them. It adds no physics: each window is an ordinary `Simulation`, and dH/dλ is the
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
(`gromos.analysis.dgslv_pbsolv`, below). `pt_top` and `make_pt_top` print each perturbed
molecule's net charge in both states as well.

A `.ptp` that perturbs a bond the topology does not have is refused (`ValueError`, gromosXX's
"Perturbation of a non-existing bond"), and so are the pre-2005 `…03` block names (`PERTBOND03`,
`PERTATOM03`, …) and any block the reader does not apply (`PERTCONSTRAINT03`, `PERTPOLPARAM`),
as gromosXX refuses the blocks it does not read. In the soft-core blocks type 0 is a term absent
in that state: `PERTBONDSOFT` with state-A type 0 makes a bond.

Slow growth (`perturbation={"dlamt": …}`, PERTURBATION DLAMT ≠ 0) runs: λ moves by DLAMT·Δt
after every step, and `Simulation` warns (`UserWarning`) when that takes λ out of [0, 1].
NSCALE ≠ 0 (the `.ptp` SCALEDINTERACTIONS block) is refused.

### Making a perturbation: `Perturbation`

`gromos.free_energy.Perturbation` holds a `.ptp` and the solute's net charge in its two states.
`Perturbation.mutate` writes one with gromos-tools' `make_pt_top`, from a topology (state A) and
a copy of it with new charges, types or masses on chosen atoms (state B); `from_ptp` reads a
file:

```python
from gromos.free_energy import Perturbation

plus = Perturbation.mutate(system, "1:1", charge_b=1.0, alpha_lj=0.0, alpha_crf=0.0)
plus.net_charge_a, plus.net_charge_b, plus.delta_charge     # (0.0, 1.0, 1.0)
dummy = Perturbation.mutate(system, "1:CM", iac_b=22)       # IAC as in the topology file
kept = Perturbation.from_ptp("ch4_spc_dummy.ptp", system)   # charges need the topology
```

- `atoms` is a gromos++ atom specifier or a sequence of 0-based indices (solute atoms only);
  `charge_b`, `iac_b` and `mass_b` take one value for every atom or one per atom, and keep state
  A's when left out. `alpha_lj`/`alpha_crf` are the soft-core parameters (`@softpar`).
- `plus.end_state(system, "B")` is the topology in state B (`pt_top`), e.g. for
  `dgslv_pbsolv`; `plus.write(path)` writes the file; `str(plus)` is its text.
- `run_windows(..., perturbation=plus)` gives it to every window as `inputs["pttopo"]`
  (written to `<workdir>/perturbation.ptp`), so the recipes need no file.

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
- `run_windows(system, recipes, perturbation=None, workdir=None, dry_run=False, parallel=1)`
  runs the windows one after the other, or `parallel` at a time in worker processes. A window samples dH/dλ every `outputs["ntwg"]` steps at the steps a `.trg` of the same
  run holds (0, NTWG, … below `control["steps"]`), so its series is the `.trg`'s first column.
  With a `workdir`, each window's energies are kept as `<workdir>/l0.250.npz`, and a window
  whose file exists is read instead of run again; `dry_run=True` returns the job list without
  running anything. `collect(workdir)` reads the kept windows back.
- `FreeEnergyWindows.ti(skip)` integrates ⟨∂H/∂λ⟩ with the trapezoid rule; its error comes
  from the block-averaging error of each window. It is the estimator of
  `gromos.analysis.ti_trapezoid`, which reads `.trg` files instead.

### Parallel windows

`parallel=N` runs N windows at a time, each in its own process (`multiprocessing` with *spawn*;
a script needs the `if __name__ == "__main__":` guard). Not threads: a `Simulation` holds the
interpreter while it runs and belongs to the thread that made it, so threads would take turns.
A worker reads the system again with `System.from_files` on the same files (the only `System`
that pickles; one built in memory is refused) and runs its window exactly as the sequential loop
does, so the numbers are the same, bit for bit. Each worker uses the whole thread pool, as a
sequential run does (the parallel kernels' sums depend on the number of threads), so processes
pay off when one window does not keep the machine busy: small systems, or recipes with
`execution.parallel = "serial"`.

## BAR and MBAR

Bennett's acceptance ratio and MBAR need each frame's energy at the other windows' λ as well.
The engine computes them with extended TI (`PRECALCLAM`): on every NTWG step, the energies of
every perturbed term at `nr_lambdas` λ points from `min_lambda` to `max_lambda`. Choose the grid
so that it holds every window's λ, and sample energies and dH/dλ together (NTWE = NTWG):

```python
base = base.update(
    outputs={"ntwg": 5, "energy_every": 5},
    perturbation={"precalclam": {"nr_lambdas": 5, "min_lambda": 0.0, "max_lambda": 1.0}},
)
windows = run_windows(system, LambdaSchedule([0.0, 0.25, 0.5, 0.75, 1.0]).recipes(base),
                      workdir="fep_ch4")
bar = windows.bar(skip=2)            # Estimate: .dg, .error, .steps
dg, error = windows.mbar(skip=2)
for (lam_i, lam_j), (step, step_error) in zip(zip(bar.lambdas, bar.lambdas[1:]), bar.steps):
    print(f"{lam_i:.2f} → {lam_j:.2f}: {step:8.3f} ± {step_error:.3f} kJ/mol")
```

- Each window keeps `Simulation.precalclam` of every sampled frame (`Window.precalclam`, and in
  its `.npz`, so `collect()` gives windows BAR and MBAR still work on).
- `bar(skip, temperature=None)` and `mbar(skip, temperature=None)` drop the first `skip` frames
  of each window, as `ti(skip)` does. The temperature is the windows' first bath unless given.
- The energy of a frame at λ is what `ext_ti_ana @bar_data` computes from the same rows
  (`windows.bar_data(skip)` returns them), and the estimators are `gromos.analysis.bar` (the
  corrected analytical error, with the statistical inefficiency of each side) and
  `gromos.analysis.mbar` (the asymptotic error, which assumes uncorrelated frames and so comes
  out smaller than BAR's on frames a few steps apart). No estimator arithmetic is done in Python.
- Windows run with `LAMBDAS` rows (individual λ per interaction) are refused for now.

The `md` binary writes the same data into its `.tre`/`.trg`; `gromos.analysis.ext_ti_ana`,
`bar` and `mbar` read those files (`scripts/ti_ch4.py --estimator bar|mbar`).

## A charge change and its finite-size correction

Charging an atom under a reaction field gives a ΔG that depends on the box and the cut-off;
`dGslv_pbsolv` computes the correction from Poisson–Boltzmann solves on the end states
(CB/NPBC − LS/PBC, plus LS/PBC[FFT] − RF/PBC[FFT] for a reaction field):

```python
from gromos import analysis

def pb(state):
    path = f"state_{state}.top"
    plus.end_state(system, state).write(path)
    es = base.forcefield["electrostatics"]
    return analysis.dgslv_pbsolv(
        topo=path, coord="ch4_spc.cnf", pbc="r", atoms="1:a", atoms_to_charge="1:a",
        scheme="RF", eps_solv=66.6, eps_rf=es["epsilon_rf"], rcut=es["rf_cutoff"],
        gridspacing=0.1, npbc_size=5.0, probe_iac=5,
    )["documented_correction"]

dg, error = windows.bar(skip=20)
corrected = dg + pb("B") - pb("A")
```

`tests/test_free_energy_charge.py` does exactly this for CH4 → CH4⁺ in 999 SPC: seven windows
(600 steps each, four worker processes), TI, BAR and MBAR within twice their combined errors
(about −242 kJ/mol), and a correction of −77.98 kJ/mol whose parts agree with gromos++'s
`dGslv_pbsolv` on the same state to 0.01 kJ/mol.

## Correcting a charge change

Kastenholz and Hünenberger (J. Chem. Phys. 124, 224501 (2006)) split what a periodic, cut-off
simulation gets wrong about charging an ion into five terms, each added to the raw ΔG:

| term | what it corrects | where it comes from |
|------|------------------|---------------------|
| A1 | no ion–solvent interaction beyond the cut-off | `dgslv_pbsolv` (and `charge_correction`, for a sphere) |
| A2 | wrong solvent polarisation inside the cut-off | `dgslv_pbsolv` (and `charge_correction`, for a sphere) |
| B | finite box, periodic copies of the ion | `dgslv_pbsolv` (and `charge_correction`, for a sphere) |
| C1 | the potential summed atom by atom: the solvent's quadrupole-moment trace | `charge_correction` only |
| C2 | the interface potential (lattice sums; 0 with a cut-off) | `charge_correction` only |

`dGslv_pbsolv` is continuum electrostatics on the real solute: it gives the polarisation terms,
A1 + A2 + B under a reaction field (B alone for a lattice sum), and no continuum sees C1 or C2,
which come from the discrete solvent molecules. `gromos.analysis.charge_correction` returns all
five for a spherical ion, their sum `total`, and `beyond_pb` = C1 + C2, the part to add on top of
`dGslv_pbsolv`.

Under a reaction field one more term comes off first: the **reaction-field self energy**.
GROMOS (gromosXX and gromos-rs alike) books E_self = −½ q_i² FPEPSI (1 − C_rf/2)/R_C for every
charged atom, and the perturbed loop books its λ-derivative, so the raw TI/BAR/MBAR ΔG of a
charge change contains

ΔG_self = −½ FPEPSI (1 − C_rf/2)/R_C · Σ_i (q_B,i² − q_A,i²)

exactly (C_rf from RCRF, EPSRF, RKAPPA and EPSCS as in the engine). The Kastenholz–Hünenberger
scheme is built for a raw ΔG without it (Reif, Hünenberger & Oostenbrink, J. Chem. Theory
Comput. 8, 3705 (2012); Reif & Oostenbrink, J. Comput. Chem. 35, 227 (2014), which apply the
scheme to GROMOS's reaction field). It is not part of A1 or of `dGslv_pbsolv`'s correction; it
is a constant of the cut-off and the charges. The corrected value is therefore

ΔG_cor = [ΔG_raw − ΔG_self] + A1 + A2 + B + C1 + C2, or, with `dGslv_pbsolv`,
ΔG_cor = [ΔG_raw − ΔG_self] + PB(B) − PB(A) + C1 + C2 (C2 = 0 under a reaction field).

`gromos.analysis.rf_self_energy` gives ΔG_self and `gromos.analysis.corrected_charging_dg` the
whole chain. For q 0 → 1 at R_C = 1.4 nm (κ = 0, EPSCS 1) ΔG_self is −73.83 (EPSRF 61), −73.88
(66.6) and −73.96 kJ/mol (78.5), next to A1 = −48.81, −48.88 and −48.99: forgetting to remove
it leaves the corrected hydration free energy of a monovalent ion about 74 kJ/mol too negative.
The full chain for the CH4 → CH4⁺ windows and `pb()` above:

```python
dg, error = windows.bar(skip=20)          # raw ΔG: windows.ti(...) or windows.mbar(...) alike
es = base.forcefield["electrostatics"]
rf_self = analysis.rf_self_energy(        # what the raw ΔG holds of the RF self term
    plus.net_charge_a,                    # one perturbed atom: its charges in A and B
    plus.net_charge_b,
    es["rf_cutoff"],
    es["epsilon_rf"],
)
cc = analysis.charge_correction(
    plus.delta_charge,                    # Q_B − Q_A (e)
    3.106,                                # box edge (nm): ch4_spc.cnf's cubic box
    es["rf_cutoff"],
    es["epsilon_rf"],
    "SPC",
    charge_a=plus.net_charge_a,
    ion_radius=0.2,                       # nm; the terms barely move within ±0.05 nm
)
corrected = analysis.corrected_charging_dg(dg, rf_self, cc, pb_a=pb("A"), pb_b=pb("B"))
# = (dg - rf_self) + (pb("B") - pb("A")) + cc["beyond_pb"]
```

For an ion small enough to be a sphere, `corrected_charging_dg(dg, rf_self, cc)` =
`(dg - rf_self) + cc["total"]` is the same correction without the Poisson–Boltzmann solves.
With several perturbed atoms, pass their charges atom by atom: the self term sums q² per atom,
not the net charge. The terms, with the paper's equation numbers (its Sec. II):

- **A1** = −(8πε₀)⁻¹ q² (1 − ε_S⁻¹)/R_C (Eq. 25), ε_S = `eps_solvent` (78, water).
- **A2** solves the paper's radial integral equation for the polarisation inside the cut-off
  (Eqs. 26–29, kernels A3 and A6) with ε_S = `eps_model`, the model's permittivity (66.6 for
  SPC; give it for the other models).
- **B** is Eq. 32 for a lattice sum and the paper's fit to continuum electrostatics for a cut-off
  (Eq. 33, fitted for ε_S = 78 and L/(2R_C) up to 1.7).
- **C1** = −q ΔΦ_ODL times a scheme factor (Eqs. 37, 41, 42), ΔΦ_ODL = ρ γ/(6ε₀) (Eq. 34), with
  γ = Σ qᵢ rᵢ² over the solvent's charges around its oxygen: 0.0082 e nm² for SPC (ΔΦ_ODL =
  79.76 kJ mol⁻¹ e⁻¹ at 1 g cm⁻³), 0.008476 for SPC/E, 0.007641 for TIP3P, 0.009295 for TIP4P.
- **C2** is 0 with a cut-off (Eq. 46) and the paper's SPC fit for lattice sums (Eq. 47; other
  solvents are refused).

A1, A2 and B are quadratic in the charge and C1, C2 linear, so a change q_A → q_B scales them
with q_B² − q_A² and q_B − q_A (`charge_a` is q_A). `scheme="RF"` (or `"BM"`) is GROMOS's
reaction field with charge groups; `"BA"` an atomic cut-off, `"CM"` a plain Coulomb cut-off,
`"LS"` a lattice sum with tinfoil boundaries. `eps_rf=0` is a conducting boundary, as EPSRF 0 is
in the engine. The corrected value is intrinsic: it excludes the surface term q Δχ° (the paper's
Eq. 49, +70.3 q kJ/mol for SPC), the cavity and the standard-state terms.

`tests/test_charge_correction.py` and the Rust tests check every term against the paper's
Tables III and IV (Na⁺ in SPC, R_I = 0.2 nm, ten lattice-sum boxes, CM/BM/BA at R_C = 0.8–1.6
nm): A1, C1 and the lattice-sum terms to 0.01 kJ/mol, B with a cut-off to 0.02, A2 under the
reaction field to 0.06. A2 for the plain cut-off is the exception: our converged solve lies
0.07–0.25 kJ/mol below the table, the gap shrinking as R_C grows, most likely the table's own
discretisation.

## Checked against gromosXX

`tests/test_free_energy.py` runs the CH4 → dummy windows at λ = 0, 0.25, 0.5, 0.75 and 1
(CH4 in 999 SPC, 10 steps) from Python: every dH/dλ sample agrees with gromosXX's `.trg` for the
four windows that have a reference (relative 1e-6), and `ti()` equals `analysis.ti_trapezoid`
on those `.trg` files. Every `PRECALCLAM` reference (nine systems: CH4 and methanol in water,
alanine dipeptide with bonded, SHAKE, distance-restraint, `LAMBDAS` and atom-pair
perturbations) run as a window gives gromosXX's `# precalclam` and `# ABdih` rows of every
frame, energies and derivatives (relative 1e-9, the printed precision). Nine CH4 windows on a
nine-point grid (300 steps each) give `bar()` and `mbar()` equal to `analysis.bar` and
`analysis.mbar` on the same frames (1e-10), and TI, BAR and MBAR within twice their combined
errors of each other.

!!! note "The reaction-field self term"
    The raw ΔG of a charge change under a reaction field contains the GROMOS self term
    −½ FPEPSI (1 − C_rf/2)/R_C Σ (q_B² − q_A²) (gromosXX `nonbonded_term.cc:79-99` and the
    perturbed loop's derivative of it). The self term is not a part of the continuum response
    that A1 describes: it is a constant of the cut-off that the Kastenholz–Hünenberger raw ΔG
    does not contain. Subtract it (`rf_self_energy`) before adding A1 + A2 + B + C1 + C2 or the
    `dGslv_pbsolv` correction (Reif, Hünenberger & Oostenbrink, JCTC 2012; Reif & Oostenbrink,
    JCC 2014). Leaving it in shifts a monovalent ion's corrected value by about −74 kJ/mol at
    R_C = 1.4 nm.
