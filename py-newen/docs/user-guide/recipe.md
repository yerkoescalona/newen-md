# The recipe model

A run in newen is described **once**, as data, and that description is shared by every
entry point: the `md` command-line binary, a GROMOS `.imd` file, and the Python objects. This
page explains the pieces and why they are shaped this way.

## One description, two front-ends

```
.imd file  ──read losslessly──▶  Recipe  ◀──build / edit──  Python (Recipe.nvt(...), .update(...))
                                    │
                                    │ recipe.plan(system)          stage 1 — the MD step as data
                                    ▼
                                  Plan  = [forcefield, orchestrator?, leap_frog_velocity, thermostat?, ...]
                                    │ validate                     GROMOS step-order invariants
                                    │ instantiate                  the real Rust algorithms
                                    ▼
                               Simulation                          step / run
```

All of this lives in one Rust crate, `newen-run`, which the `md` binary and the Python binding
both call. The Rust reference suite drives the binary against double-precision gromosXX output
on 40 systems; the Python suite drives the binding; `tests/test_front_end_parity.py` checks the
Python front-ends against each other with exact equality. A change to how a run is built is
therefore seen from both sides, and there is no second builder that could drift.

## `Recipe` — the run as data

A `Recipe` holds exactly what a GROMOS `.imd` holds, grouped by concern:

| Group | GROMOS blocks | Examples of fields |
|---|---|---|
| `system` | SYSTEM | `npm`, `nsm` (a hint — the coordinate file decides) |
| `control` | STEP, INITIALISE, COMTRANSROT | `steps`, `dt`, `seed`, `generate_velocities`, `com_motion` |
| `boundary` | BOUNDCOND | `kind` (`vacuum` / `rectangular` / `triclinic` / `truncated_octahedron`), `ndfmin` |
| `forcefield` | FORCE, PAIRLIST, NONBONDED, restraint blocks | `bonds`…`nonbonded` switches, `pairlist`, `electrostatics`, `restraints`, `terms` |
| `constraints` | CONSTRAINT | `solute` (`none`/`hbonds`/`allbonds`), algorithms, tolerances |
| `ensemble` | MULTIBATH, PRESSURESCALE | `thermostat` (baths, τ), `barostat` |
| `minimisation` | ENERGYMIN | steepest-descent controls |
| `perturbation` | PERTURBATION | λ, soft-core, `pttopo` handling |
| `outputs` | WRITETRAJ, PRINTOUT | write frequencies |
| `execution` | — | `parallel`: `auto` / `serial` / `parallel` (not physics) |
| `inputs` | the `md @pttopo @posresspec @refpos @distrest` files | auxiliary file paths |
| `passthrough` | anything unmodelled the caller allowed | raw block text |

Three properties matter:

- **Lossless.** `Recipe.from_imd(path).to_imd()` reproduces the file, and gromosXX runs the
  rewritten file identically (checked on every reference `.imd` with the real gromosXX,
  `scripts/roundtrip_imd_gromosXX.py`). What the engine *does not* model is an error unless you
  name it in `allow_passthrough`; an optional block that is *absent* is reported in
  `recipe.diagnostics` together with what gromosXX does without it (no bath, no COM-motion removal).
- **Immutable, strict.** `recipe.update(control={"steps": 20000})` returns a new recipe with that
  one field changed; `recipe.update(control={"stepz": 1})` is a `RecipeError`. There is no
  silent default anywhere between you and the engine.
- **Serialisable.** `to_toml()` / `to_json()` / `to_dict()` / pickling round-trip exactly;
  `to_bundle(dir, system, topo_path, conf_path)` writes `input.toml` + `run.recipe.toml` +
  `run.imd`, and `Recipe.from_bundle` / `Simulation.from_bundle` read it back. The `md` binary
  writes the effective recipe next to its `.tre` for the same reason.

Factories give you a starting point without authoring a file: `Recipe.nve(dt, steps)`,
`Recipe.nvt(dt, steps, temperature)`, `Recipe.npt(dt, steps, temperature, pressure)`,
`Recipe.minimize(steps)` — each with `constraints="none"|"hbonds"|"allbonds"`. Their defaults
are derived from the engine's own `.imd` defaults, and a GROMOS convention is kept where it
exists: `Recipe.nve` carries a bath with coupling time −1, which is how gromosXX itself says
"no coupling" (the `repr` still says NVE).

## `Plan` — the MD step as data

`recipe.plan(system)` is stage one of the builder: the ordered algorithms of a step, every
parameter resolved (cut-offs, degrees of freedom, `four_pi_eps_i` from the topology, …), already
validated against the GROMOS ordering rules:

```python
plan = recipe.plan(system)
plan.kinds
# ['remove_com', 'forcefield', 'leap_frog_velocity', 'thermostat', 'leap_frog_position',
#  'shake', 'temperature_calculation', 'energy_calculation']
plan["forcefield"]                    # Algorithm(kind='forcefield', cutoff_long=1.4, ...)
plan.remove("remove_com")
plan.insert_after("forcefield", Algorithm("orchestrator", terms=[...]))
plan.validate()                        # PlanError if an invariant is broken
sim = Simulation(system, recipe, plan=plan)   # re-validated, then instantiated
```

`newen.algorithms()` lists every kind with its parameters and rules (which kinds are unique,
required, first/last, what must follow what). A plan that would not run is rejected with a
named reason before anything is built — the second builder that used to accept anything is gone.

## `Term` — additive potentials

A term is an energy/force provider over a region, added on top of the classical force field
(`coupling="delta"`). The kinds this build knows are in `newen.terms()`:

```python
Term("xtb", region="1:a", elements=[8, 1, 1, 8, 1, 1], gfn=2, charge=0, multiplicity=1,
     work_dir=None, timeout_s=600)                  # real GFN-xTB subprocess, no feature needed
Term("schnet", model="model.pt", cutoff=0.5, elements=[...], region="1:a", buffer=None)
                                                    # TorchScript SchNetPack, --features ml
```

`region` uses the newen atom-specifier syntax (`1:a` = molecule 1, all atoms; `s:a` = all
solvent; `1-3,7` = global atom numbers). `elements` is one atomic number per global atom index.
The plan enforces what the physics needs: `coupling="replace"` waits for the zone-aware force
field (2.8) and is rejected by name; a term without a virial cannot meet a barostat; a
term whose cargo feature this build lacks constructs (`Term(...).available == False`) but
`Simulation` raises `MissingFeatureError`. Each term's energy is visible on its own,
`sim.term_energies`, so two terms cannot hide a compensating error.

Adding a term kind to the engine is three Rust files (`TermSpec` variant, registry arms, one
`instantiate` arm) — measured with `xtb` in 3.9 step 5. Nothing in the binding changes.

## Reading a simulation back

`sim.recipe`, `sim.plan` and `sim.diagnostics` are what the engine actually used (`sim.recipe_toml`
and `sim.plan_json` as text). `sim.potential_energy` is the classical force field,
`sim.term_energies` the terms, `sim.total_energy` everything plus kinetic.

## Migrating

The pre-recipe forms (`InputParameters`, the `distrest=`/`posresspec=`/`refpos=`/`ml_*=`
keywords, `AlgorithmSequence.*`, `Simulation.from_sequence`) keep working for one release: each
emits a `DeprecationWarning` naming the replacement and is translated into a recipe — not a
second path. The table is in the [Quick Start](quick-start.md#migrating-from-inputparameters).
