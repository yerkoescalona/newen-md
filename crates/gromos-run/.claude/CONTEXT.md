# gromos-run — stage contract

## Job
L3 orchestration, as a library: the **one** place that turns "what to run" (today: an
`ImdParameters` + topology + coordinates + auxiliary inputs) into a running `AlgorithmSequence`.
The `md` binary (`gromos-md`) and the Python binding (`pyo3-gromos`) both call it — that is the
structural fix for the three drifting builders 3.8 found (drift guard G1 of 3.9).

## Inputs (consumes from)
gromos-core, gromos-forces, gromos-integrators, gromos-io. No clap/log/env_logger, no MPI/CUDA.

## Outputs
```
RunRecipe::from_imd(&ImdParameters)          -> (RunRecipe, Diagnostics)   the .imd front-end
RunRecipe::to_imd() / to_imd_string(n_atoms) -> ImdParameters / .imd text  (gromosXX-readable)
RunRecipe::{to_toml,from_toml,to_json,from_json}, ::{nve,nvt,npt,minimize}
prepare_system(imd, topology, physical_constants, coords, &RunInputs)  -> Prepared (+ notes)
build_plan(&RunRecipe, &Topology, &RunInputs, four_pi_eps_i)          -> Vec<AlgorithmSpec>  stage 1
validate_plan(&[AlgorithmSpec])                                         -> ordering invariants  (G9)
instantiate(&[AlgorithmSpec], &Topology, &Configuration, &Periodicity) -> Instantiated         stage 2
                                                                           (reads ONLY the plan, G8)
build_sequence_from_recipe / build_sequence_from_imd                    -> Built { sequence, recipe, plan, diagnostics, summary }
start(&mut sequence, &topo, &mut conf, &state)                          -> init + step 0        stage 3
total_dof(topo, &ConstraintSelection, ntc, ndfmin)                      -> the one DOF formula
AlgorithmSpec::{KINDS, name, examples, rules}, TermSpec::{KINDS, name, examples, feature}
RunError                                                                 -> every failure, as data
```

## Status
- **A-EDS (1.8, 0.0.60):** `RunRecipe::aeds` (the `AEDS` block), `check_aeds` (refusals, called
  from `check_eds`), `EdsPlan::aeds` (FORM 1, s = 1, the search settings). `prepare_system` takes
  the end states from an `AEDS` block as from `EDS`, and with NTIAEDSS = 0 in FORMS 2-6 parses the
  configuration's `AEDSSEARCH` block into `Configuration::aeds` (required; with NTIAEDSS = 1
  ignored, `PrepareNote::AedsSearchIgnored`). `Coordinates::aedssearch` carries the block's lines.
- **Net charge of a perturbation (1.8, 0.0.59):** `prepare_system` notes the solute's net
  charge in state A and B (`PrepareNote::PerturbationNetCharge`) on every perturbed run and
  `NetChargeChanged` when ΔQ ≠ 0 under a reaction field (NLRELE −1/1/4): a warning, no refusal.
- **MPI (2026-08-30):** `RunOptions::pair_partition: Option<(usize, usize)>` sets
  `Forcefield::pair_partition` on the built sequence; the caller installs the reducer
  (`sequence.find_mut::<Forcefield>()`). `tests/pair_partition.rs` checks that rank shares sum to
  the full pair terms without MPI.
- **3.9 step 1 ✓ (2026-08-29)** — assembly lifted verbatim from `md.rs:414-572` and
  `md.rs:987-1298`; the binary's `println!`/`process::exit` sites became `RunError` variants and
  `PrepareNote`s / `BuildSummary` that the binary prints. Divergences resolved *toward the binary*
  (3.9 A10/A11): `four_pi_eps_i` from the topology's PHYSICALCONSTANTS, NSM from the
  coordinate file, `ParallelPolicy::Auto` (= `n_atoms > 100`) for both callers, one DOF formula
  that now also counts solute constraints (no reference combines NTC>1 with a live thermostat, so
  no reference moved). FEP reaches Python through `Topology::apply_perturbation` (gromos-core) +
  `RunInputs.pttopo` (binary).
- **3.9 step 2 ✓ (2026-08-29)** — the recipe is in front of the builder:
  `recipe.rs` (`RunRecipe`, lossless for every modelled `.imd` block, serde with
  `deny_unknown_fields`, `version`, `passthrough` allowlist, `Diagnostics`), `plan.rs`
  (`AlgorithmSpec` fully resolved, `build_plan`, `validate_plan` with per-kind `KindRules`,
  `TermSpec` with `coupling`), `build.rs` (`instantiate` takes no recipe; `start`). The
  `.imd` parser (gromos-io) became lossless (`ntirtc`, `ntisti`, `ntwse/ntwg/ntwb`, `ntpp`,
  `ncyc`, `ntcs0`, DOF sets, PRESSURESCALE couple/SEMI, NONBONDED lattice-sum lines) and strict
  (garbage numbers are errors naming block and token), `write_imd` exists, and the defaults
  mean "absent" (no MULTIBATH → no bath; no COMTRANSROT → `nscm = 0`, gromosXX `parameter.h`).
  `legacy_builder.rs` (`cfg(test)`) is the frozen step-1 builder: `build::tests::
  plan_matches_legacy_builder_bit_for_bit` proves the recipe path identical on all 41 references.
- Step 3 (next): Python `Recipe`/`Term`/`Algorithm` over this data (pythonize), shims for the
  factories/kwargs/presets, `gromos.exceptions`. See 3.9.

## Key files
```
src/lib.rs          — pipeline doc + re-exports
src/prepare.rs      — prepare_system: pttopo merge, truncated-octahedron transform, NSM from
                      coordinates, validation, initial velocities
src/build.rs        — build_sequence_from_imd (the GROMOS step order), start, periodicity_of
src/constraints.rs  — ConstraintSelection::from_imd (NTC/NTCP/NTCS dispatch)
src/dof.rs          — total_dof
src/inputs.rs       — RunInputs (aux file paths), RunOptions / ParallelPolicy (not physics)
src/error.rs        — RunError
```

## Crate-specific rules
- **No `println!`, no `process::exit`, no logging.** Everything the binary prints comes back as
  data or as a `RunError`.
- **No second builder, no second IMD reader.** Callers that need a different sequence edit the
  plan (`build_plan` → `Plan` → `build_sequence_from_plan`); they never push algorithms
  themselves, and they read parameter files only through `read_imd` / `RunRecipe::from_imd`.
  `just lint` (G6, 3.9 step 4) greps for violations.
- **Behaviour is defined by the reference suite from both sides:** `cargo test -p gromos-md
  --test test_gromosXX_references` (drives the binary) and `py-gromos/tests/` (drives the binding)
  must both stay green; `py-gromos/tests/test_front_end_parity.py` compares the Python front-ends
  with `np.array_equal`.
- Defaults are derived from `ImdParameters::default()` (3.9 G7, step 2) — no second table.
- **Adding a term** (3.9 step 5, measured with `xtb`): a `TermSpec` variant in `recipe.rs`,
  its registry arms in `plan.rs` (`KINDS`, `name`, `examples`, `feature`, `provides_virial`,
  `coupling`) and one arm in `build.rs::instantiate_orchestrator`, registered with
  `register_labelled(<registry name>, …)` so `Energy.term_energies` reports it by name (G10).
  Nothing else — no binding change, no `md.rs` change.
- **This crate is the only entry point** for building a run: the `md` binary and `py-gromos`'s
  `Simulation` both call `prepare_system` → `build_plan` → `validate_plan` → `instantiate` → `start`
  (steps 1–4). Every `.imd`/Python front-end is a translation into a `RunRecipe`.
