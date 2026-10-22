# newen-logic — stage contract

## Job
The **decision logic** of newen, written so that it can be proven equal to a Lean model:
which steps rebuild the pairlist or write output, how many degrees of freedom a bath has, which atom
belongs to which bath, which feature combinations are refused, how indices convert. Anything that
computes an energy, a force or a velocity is physics and stays in newen-forces / newen-integrators.
Every function here has a Lean model in `spec/NewenSpec/`, a proof that the Rust equals the model
(Charon + Aeneas), and a proof that the model means what the physics says (`spec/NewenTheory/`).
Explained for newcomers in `docs/src/verification/index.md`.

## The scheme: three layers, one direction of dependency

```
spec/NewenTheory/   physics, organised like the physics (GROMOS Vol. 2, textbooks)
                     definitions over ℝ / ℕ and theorems; knows nothing about newen
        ▲  bridge theorems: "the code's model is an instance of the theory"
spec/NewenSpec/     what newen computes, organised like this crate (one Lean file per module)
                     models written to mirror the Rust; bridge theorems to NewenTheory
        ▲  Aeneas refinement (integer logic) / differential testing (float kernels, other crates)
crates/newen-logic  the Rust; NewenSpec/Extracted/ is its Aeneas translation (generated)
```

- `NewenTheory` never imports `NewenSpec`. The divergence dossier (gromos-divergences repository)
  imports `NewenTheory` to state gromosXX's transcriptions against the same definitions.
- The bridge theorems are where discretisation and bookkeeping enter (output steps, pairlist steps,
  indices); that is where the 39 gromosXX divergences lived.

## Inputs (consumes from)
Nothing. No workspace crate, no external crate in the default build (the one Charon translates; `serde` is an optional feature for the run plan's JSON), `#![forbid(unsafe_code)]`. It sits at the bottom of
the dependency graph; newen-core, newen-run, newen-forces and newen-integrators call it.

## Outputs
```
schedule::output_misses_longrange(nsnb, ntwg) -> bool     spec: schedule.output_misses_longrange
thermo::dof::total_dof(n_atoms: u64, solvent_constraints: u64, solute_constraints: u64,
                       ndfmin: i64) -> i64              spec: thermo.dof.total_dof
thermo::bath::baths_of(&[DofRange], atom) -> Option<(usize, usize)>
                                                        spec: thermo.bath.baths_of
thermo::bath::groups_starting_in(&[usize], first, last) -> u64
                                                        spec: thermo.bath.groups_starting_in
thermo::bath::bath_half_dof(bath, &[DofRange], &[usize], &[Constraint], rottrans) -> HalfDof
                                                        spec: thermo.bath.bath_half_dof
thermo::bath::bath_half_dofs(n_baths, &[DofRange], &[usize], &[Constraint], rottrans)
                       -> Vec<HalfDof>                  spec: thermo.bath.bath_half_dofs
thermo::bath::ndfmin_spread(&[i64], ndfmin: i64) -> Spread
                                                        spec: thermo.bath.ndfmin_spread
topology::pairs::contains_sorted(&[usize], j) -> bool   spec: topology.pairs.contains_sorted
topology::pairs::pair_kind(exclusions_i, one_four_i, atom_pair_i: &[usize], j) -> PairKind
                                                        spec: topology.pairs.pair_kind
topology::groups::group_of(last_atoms: &[usize], atom) -> Option<usize>
                                                        spec: topology.groups.group_of
topology::groups::check_last_atoms(last_atoms: &[usize], n_atoms) -> LastAtomsCheck
                                                        spec: topology.groups.check_last_atoms
topology::groups::energy_cell(eg_i, eg_j, n_groups) -> Option<usize>
                                                        spec: topology.groups.energy_cell
```
`DofRange` is newen-integrators' `BathRange` (re-exported there, with `in_bath` a wrapper of
`baths_of`); newen-run's `bath_dofs` adapts the topology to `bath_half_dofs` + `ndfmin_spread`.
newen-core's `Topology::is_excluded` / `pair_kind_of` / `is_excluded_or_14` (every pairlist
builder) call `contains_sorted` / `pair_kind` (re-exporting `PairKind`); newen-run's
`assign_energy_groups` calls `check_last_atoms` and `group_of`; newen-forces'
`EnergyGroups::cell` calls `energy_cell`.

## Crate rules
1. **Decisions only.** No floats, no energies or forces. A function that needs a float threshold
   takes the already-decided boolean from its caller (see `newen-run::precalclam_longrange_gaps`).
2. **Plain data in, plain data out.** Integers, `bool`, slices, `Vec`, small structs defined here.
   Never `Topology`/`Configuration`: they are too rich for Aeneas; callers adapt in a few lines.
3. **Aeneas-friendly shape.** Safe Rust, no traits or generics unless needed, simple loops (no
   nested loops for now), no iterator adaptor chains, explicit guards before `%` and `/`.
4. **One path everywhere.** Rust module `x/y.rs` ↔ `spec/NewenSpec/X/Y.lean` ↔ `SPEC.toml` id
   `x.y.<fn>`; the theory it bridges to is named in the module doc comment.
5. **A function joins only with all of:** a theory definition (and theorem where reachable), a
   model, a bridge theorem, the Aeneas refinement, an exhaustive or random Rust test, a `SPEC.toml`
   entry, and byte-identical gromosXX references for the code that now calls it.
6. **Verification gates promotion, not entry.** Experimental decision logic may live in its caller
   behind an `experimental` switch with only the self-consistency tests; it moves here when it
   becomes a default or goes into a paper or tutorial.

## Layout (now, and the agreed target)

```
src/
├── schedule/      L3 timing: pairlist steps, output steps, continuation counters    (D34; now: schedule.rs)
├── thermo/        L2: degrees of freedom per bath, atom → bath mapping   (now: dof.rs, bath.rs)
├── topology/      L0 bookkeeping: index conversions, exclusions and 1-4 ownership,
│                  groups, perturbed atoms and bond matching  (now: pairs.rs, groups.rs)
└── plan/          L3 composition: step-order rules, the compatibility/refusal matrix (target)
```

```
spec/NewenTheory/                       physics-shaped (target unless noted)
├── Foundations/   units, constants
├── StatMech/      canonical ensemble, free energy (Zwanzig, TI identity), BAR, MBAR
├── Interactions/  bonded, Lennard-Jones, Coulomb, reaction field, soft-core, λ-mixing
├── Dynamics/      leap-frog, constraints, thermostats (incl. equipartition and dof), barostats;
│                  pairlist and output steps (now: Schedule.lean, DegreesOfFreedom.lean)
├── EnhancedSampling/  EDS, A-EDS
├── Corrections/   Kastenholz–Hünenberger charging terms
└── Topology/      pair kinds, groups of consecutive atoms (now: Pairs.lean, Groups.lean)
```

## Status
- **Pilot (2026-10-17):** `schedule::output_misses_longrange`, the integer core of
  `newen-run::precalclam_longrange_gaps` (D34). Meaning (`pairlistStep`, `RowsComplete`,
  `Misses`, `rows_complete_iff`) in `spec/NewenTheory/Dynamics/Schedule.lean`; model and bridge
  `missesModel_iff` in `spec/NewenSpec/Schedule.lean`; Aeneas refinement
  `output_misses_longrange_iff` in `spec/NewenSpec/Refinement.lean`; exhaustive test
  `tests/schedule_spec.rs` (1 ≤ nsnb, ntwg ≤ 200, against a step-by-step simulation and
  Lean-generated vectors). 119/120 gromosXX references unchanged.
- **Degrees of freedom (2026-10-18):** `thermo::dof::total_dof`, the integer core of
  `newen-run::total_dof` (which still counts the constraints and converts to `f64`). Signed
  `u64`/`i64` arithmetic: a non-physical input gives a negative count instead of a wrapped
  `usize`. Meaning (`Molecule`, `internalDof`, `dof`; additivity `internalDof_append`,
  identical molecules `internalDof_replicate`, physicality `dof_nonneg_iff`) in
  `spec/NewenTheory/Dynamics/DegreesOfFreedom.lean`; model and bridge `totalDofModel_eq_dof`
  (solute plus `k` solvent copies) in `spec/NewenSpec/Thermo/Dof.lean`; Aeneas refinement
  `total_dof_eq` / `total_dof_eq_dof` in `spec/NewenSpec/Refinement.lean`, for
  `3·n_atoms + solvent + solute + |ndfmin| ≤ i64::MAX`, and `total_dof_overflow` (a `3·n_atoms`
  past `i64::MAX` fails, never returns); test `tests/dof_spec.rs` (molecule-by-molecule count,
  128-bit arithmetic on large counts, Lean vectors `spec/vectors/dof.json`). Per-bath counting
  followed on 2026-10-21 (below).
- **Per-bath degrees of freedom (2026-10-21):** `thermo::bath` (`baths_of`, `groups_starting_in`,
  `bath_half_dof`, `bath_half_dofs`, `ndfmin_spread`), the integer core of
  `newen-run::bath_dofs` and of newen-integrators' `in_bath` (whose `BathRange` is now
  `DofRange`). Counts in half degrees of freedom, so a constraint shared between two baths stays
  an integer; the NDFMIN spread as exact numerators over one denominator (gromosXX truncates the
  total, D40). Meaning (`BathSystem`, `comHalf`, `irHalf`, atom by atom; partition
  `bathHalf_partition`, one bath `bathHalf_single` tied to `dof`; `spread_sum`, `spread_zero`,
  `spread_total_nonpos`) in `spec/NewenTheory/Dynamics/DegreesOfFreedom.lean`; models and
  bridges (`bathHalfDofModel_eq` under `Covers`/`RangesWellFormed`, `ndfminSpreadModel_eq`) in
  `spec/NewenSpec/Thermo/Bath.lean`; Aeneas refinements `baths_of_eq`, `groups_starting_in_eq`,
  `bath_half_dof_eq`, `bath_half_dofs_eq`, `ndfmin_spread_eq` (and `*_eq_theory`) in
  `spec/NewenSpec/Refinement.lean`, loops by `loop.spec_decr_nat` with an "accumulated + model
  of the rest" invariant, for `BathHalfDofInRange` / `SpreadInRange` (no `i64` overflow); test
  `tests/bath_spec.rs` (atom-by-atom count on every system of ≤ 6 atoms with 1–3 ranges and
  baths, exact spread by 128-bit cross-multiplication, Lean vectors `spec/vectors/bath.json`,
  the two DOFSET splits of `aladip_vacuum_multibath_ndfmin3`). A private helper `end_in_bath`
  (one constraint end) keeps the constraint loop free of nested matches. Not modelled: virtual
  atoms, position/angle/dihedral constraints. gromosXX references unchanged.
- **Topology bookkeeping (2026-10-22):** `topology::pairs` (`contains_sorted`, `pair_kind`) and
  `topology::groups` (`group_of`, `check_last_atoms`, `energy_cell`), the integer cores of
  newen-core's pair tests (every pairlist builder), newen-run's `assign_energy_groups` and
  newen-forces' `EnergyGroups::cell`. Meaning (`pairKind` by membership, `pairKind_symm` for
  symmetric lists; `groupOf` as the least index, `LastAtomsOk`, `groupOf_some_iff`,
  `groupOf_contiguous`, `groupOf_mono`, `energyCell_injective`) in
  `spec/NewenTheory/Topology/{Pairs,Groups}.lean`; models and bridges
  (`containsSortedModel_iff`, `pairKindModel_eq` for ascending lists; `groupOfModel_eq`,
  `checkLastAtomsModel_iff`, `energyCellModel_eq`, and `bathsOf_eq_groupOf` tying
  `thermo.bath.baths_of` to the same rule) in `spec/NewenSpec/Topology/{Pairs,Groups}.lean`;
  Aeneas refinements `contains_sorted_eq`, `pair_kind_eq`, `group_of_eq`, `check_last_atoms_eq`,
  `energy_cell_eq` (and `*_iff` / `*_eq_theory`) in `spec/NewenSpec/Refinement.lean`,
  `energy_cell` for `EnergyCellInRange`; test `tests/topology_spec.rs` (membership and least
  index by brute force, Lean vectors `spec/vectors/topology.json`, butane's lists and charge
  groups, three FORCE `NRE` lines). `one_four_pairs` is now sorted and deduplicated like
  `exclusions` (asserted in newen-io's `exclusion_invariants`), and the exclusion search is
  linear after the early-out, as in gromosXX (it was a binary search). Not modelled: the list
  construction itself (symmetry and order are unit-tested, not proven).
- **NewenTheory (2026-10-17):** created as a second default library of the `spec/` Lake project,
  core Lean only (imports neither `NewenSpec` nor Aeneas); the pilot's meaning moved there with
  no change to any statement beyond its namespace. `SPEC.toml` names it in the `theory` field.
- **Toolchain:** Lean v4.35.0 and Mathlib v4.35.0, fixed by the Aeneas Lean library; Aeneas +
  Charon from release `nightly-2026.10.17-e62dd83` (`~/.local/aeneas`), Charon on Rust
  `nightly-2026-09-17`. `NewenTheory` must use the same toolchain.

## Key files
`src/schedule.rs`, `tests/schedule_spec.rs`, `spec/NewenTheory/Dynamics/Schedule.lean`,
`spec/NewenSpec/{Schedule,Refinement}.lean`, `src/thermo/dof.rs`, `tests/dof_spec.rs`,
`spec/NewenTheory/Dynamics/DegreesOfFreedom.lean`, `spec/NewenSpec/Thermo/Dof.lean`, `src/thermo/bath.rs`,
`tests/bath_spec.rs`, `spec/NewenSpec/Thermo/Bath.lean`, `src/topology/{pairs,groups}.rs`,
`tests/topology_spec.rs`, `spec/NewenTheory/Topology/{Pairs,Groups}.lean`,
`spec/NewenSpec/Topology/{Pairs,Groups}.lean`, `spec/NewenSpec/Extracted/NewenLogic.lean`
(generated), `spec/scripts/extract.sh`, `spec/SPEC.toml`, `docs/src/verification/index.md`.
