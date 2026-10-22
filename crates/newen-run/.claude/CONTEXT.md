# newen-run — stage contract

## Job
L3 orchestration, as a library: the **one** place that turns "what to run" (today: an
`ImdParameters` + topology + coordinates + auxiliary inputs) into a running `AlgorithmSequence`.
The `md` binary (`newen-md`) and the Python binding (`pyo3-newen`) both call it — that is the
structural fix for the three drifting builders 3.8 found (drift guard G1 of 3.9).

## Inputs (consumes from)
newen-core, newen-forces, newen-integrators, newen-io. No clap/log/env_logger, no MPI/CUDA.

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
- **Nosé–Hoover restart (0.0.69, D41):** `prepare_system` parses the configuration's NHCVARIABLES
  into `Configuration::nose_hoover` for a Nosé–Hoover run with NTINHT = 0 (with NTINHT = 1:
  `PrepareNote::NoseHooverVariablesIgnored`; without a block,
  `PrepareNote::NoseHooverVariablesMissing`, a warning); the plan's `Thermostat` carries
  `read_variables`. Every link of every bath is restored, plain Nosé–Hoover included. The
  legacy-builder parity test skips a run that restores NHCVARIABLES (the frozen builder never
  did).
- **gromosXX's refusals (0.0.68):** `check_start_time` (STEP T outside ">=0", and T = −1, the
  time from the TIMESTEP block, which newen does not read), `check_pairlist` (grid pairlist with
  NTB = −1, gromosXX's FeatureChecker message) and `check_dihedral_phases` (NTBDN = 1 with a
  TORSDIHEDRALTYPE phase other than 0°/180°, `in_topology.cc:960-962`) in `plan.rs`, run by
  `build_plan` and `build_sequence_from_plan`. No reference input trips them.
- **More refusals (2026-10-22):** `check_com_motion` (NTICOM outside 0..=2,
  `in_parameter.cc:1166`), `check_outputs` (NTWSE ≠ 0, parsed but never read),
  `check_perturbation` (NLAM ≤ 0 whatever NTG says, `in_parameter.cc:882`; RCRF ≤ 0 in a perturbed
  run, where the soft-core reaction field divides by the softened cutoff). `thermostat_baths`
  refuses a repeated DOFSET LAST (the ranges must increase strictly) and ROTTRANS atoms
  `1..=RTCLAST` spread over more than one bath pair (`topology.cc:1370-1394`); one bath is
  counted by `bath_dofs` like several. newen-io's topology reader refuses LJEXCEPTIONS with
  NEX > 0. No reference input trips them.
- **Ewald and P3M runs (2026-10-22):** `check_electrostatics` (`plan.rs`, called by
  `build_plan` and `build_sequence_from_plan`) checks the NONBONDED lattice-sum ranges gromosXX
  reads for every NLRELE (NSHAPE −1..10, ASHAPE > 0, NA2CLC 0..4, TOLA2 > 0, EPSLS 0 or ≥ 1,
  NK* > 0, KCUT > 0, NG* > 0, NASORD 1..5, NFDORD 0..5, NALIAS > 0, NQEVAL ≥ 0, FACCUR > 0,
  NRDGRD/NWRGRD 0 or 1; NA2CLC 1 only with Ewald, 4 only with P3M; an odd grid and influence-
  function I/O are errors whatever the NLRELE), refuses a P3M grid that is not a power of two in
  every dimension (PhastFT; gromosXX accepts any even size), and for NLRELE 2 and 3 the
  features gromosXX's FeatureChecker keeps locked against `ewald` and `p3m` (the same list): the
  charge-group cutoff (PAIRLIST TYPE atomic is required), the grid-cell pairlist, vacuum, the
  truncated octahedron, several energy groups, perturbation, EDS, A-EDS, GaMD, QM/ML terms, the
  conjugate-gradient minimiser, LINCS/SETTLE for the solvent. `ForcefieldPlan::lattice_sum`
  (`LatticeSumPlan`) carries NSHAPE, ASHAPE, NK, KCUT, NA2CLC (4 stored as 3, as gromosXX), TOLA2,
  EPSLS and, for P3M, `p3m: Option<P3mPlan>` (grid, NASORD, NFDORD, NALIAS, NQEVAL, FACCUR);
  `instantiate_forcefield` builds a `newen_forces::nonbonded::ewald::LatticeSum::{Ewald, P3m}`
  from it with `1/ε₀ = FPEPSI · 4 · π`; the RF excluded/self terms are off for NLRELE 2/3. A plan
  whose force field disagrees with the recipe's NLRELE (no lattice sum, Ewald or P3M) is refused
  (`check_plan_against_recipe`). References: `nacl_pair_box_ewald`, `water_216_box_ewald`,
  `water_216_box_ewald_epsls`; `nacl_pair_box_p3m_fixed`, `water_216_box_p3m_fixed`,
  `water_216_box_p3m_gauss_fixed` (gromosXX + `patches/LS1_p3m_nvt_a_term.patch`);
  `water_216_npt_ewald_epsls_fixed` (LS2, surface virial), `nacl_water_box_triclinic_p3m_fixed`
  (LS3, triclinic wave vectors, ik differentiation).
- **GENBOX angles and the cell height (2026-10-22):** `Coordinates` carries the GENBOX angles and
  Euler angles; `prepare_system` builds NTB = 2 with an angle ≠ 90° as `SimBox::triclinic`
  (`newen_core::math::genbox_box`, gromosXX `in_configuration.cc:1977-2005`;
  `Prepared::triclinic_box`, used by `periodicity_of`), keeps NTB = 2 with right angles
  rectangular, refuses NTB = 1 with an angle ≠ 90° (gromosXX's message) and any non-zero Euler
  angle (gromosXX stores them unapplied on input but rotates by them on output).
  `check_cutoffs_fit_box` refuses a periodic run whose largest of RCUTP, RCUTL, RCRF is not below
  half the smallest cell height (`newen_core::math::box_heights`; gromosXX
  `boundary_check_cutoff` checks RCUTL alike, `boundary_checks.cc:35-101`). `md` writes GENBOX
  (NTB, lengths, angles) into the trajectory and the final configuration, a truncated octahedron as
  its cube. Reference: `nacl_water_box_triclinic`.
- **Slow growth (1.8, 0.0.61):** `advance_lambda(&recipe, &mut topo)` (`build.rs`) is gromosXX's
  `Slow_Growth`: λ += DLAMT·Δt and `Topology::update_for_lambda`, a no-op unless NTG = 1 and
  DLAMT ≠ 0. It needs `&mut Topology`, so it is not an `Algorithm`: the callers (`md`, `Simulation`)
  call it after every step, `start`'s included. `check_slow_growth` (`plan.rs`) refuses gromosXX's
  locks; DLAMT < 0 is refused whatever NTG says. Since 0.0.64 `advance_lambda` returns the new
  λ on the step that takes it out of [0, 1] and `prepare_system` notes
  `PrepareNote::LambdaOutOfRange { lambda }` when a run starts outside it (a `PERTDATA` λ); `md`
  warns and `Simulation` raises a Python `UserWarning`, both with `lambda_warning`'s text. λ stays
  unbounded, as in gromosXX (never a refusal).
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
  no reference moved). FEP reaches Python through `Topology::apply_perturbation` (newen-core) +
  `RunInputs.pttopo` (binary).
- **3.9 step 2 ✓ (2026-08-29)** — the recipe is in front of the builder:
  `recipe.rs` (`RunRecipe`, lossless for every modelled `.imd` block, serde with
  `deny_unknown_fields`, `version`, `passthrough` allowlist, `Diagnostics`), `plan.rs`
  (`AlgorithmSpec` fully resolved, `build_plan`, `validate_plan` with per-kind `KindRules`,
  `TermSpec` with `coupling`), `build.rs` (`instantiate` takes no recipe; `start`). The
  `.imd` parser (newen-io) became lossless (`ntirtc`, `ntisti`, `ntwse/ntwg/ntwb`, `ntpp`,
  `ncyc`, `ntcs0`, DOF sets, PRESSURESCALE couple/SEMI, NONBONDED lattice-sum lines) and strict
  (garbage numbers are errors naming block and token), `write_imd` exists, and the defaults
  mean "absent" (no MULTIBATH → no bath; no COMTRANSROT → `nscm = 0`, gromosXX `parameter.h`).
  `legacy_builder.rs` (`cfg(test)`) is the frozen step-1 builder: `build::tests::
  plan_matches_legacy_builder_bit_for_bit` proves the recipe path identical on all 41 references.
- Step 3 (next): Python `Recipe`/`Term`/`Algorithm` over this data (pythonize), shims for the
  factories/kwargs/presets, `newen.exceptions`. See 3.9.

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
- **Behaviour is defined by the reference suite from both sides:** `cargo test -p newen-md
  --test test_gromosXX_references` (drives the binary) and `py-newen/tests/` (drives the binding)
  must both stay green; `py-newen/tests/test_front_end_parity.py` compares the Python front-ends
  with `np.array_equal`.
- Defaults are derived from `ImdParameters::default()` (3.9 G7, step 2) — no second table.
- **Adding a term** (3.9 step 5, measured with `xtb`): a `TermSpec` variant in `recipe.rs`,
  its registry arms in `plan.rs` (`KINDS`, `name`, `examples`, `feature`, `provides_virial`,
  `coupling`) and one arm in `build.rs::instantiate_orchestrator`, registered with
  `register_labelled(<registry name>, …)` so `Energy.term_energies` reports it by name (G10).
  Nothing else — no binding change, no `md.rs` change.
- **This crate is the only entry point** for building a run: the `md` binary and `py-newen`'s
  `Simulation` both call `prepare_system` → `build_plan` → `validate_plan` → `instantiate` → `start`
  (steps 1–4). Every `.imd`/Python front-end is a translation into a `RunRecipe`.
