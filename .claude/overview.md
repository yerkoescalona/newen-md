# newen — workspace guide

newen ports **gromosXX/md++** (MD engine) + **gromosPlsPls/gromos++** (analysis + tools) into
one Rust workspace with a single shared core (renamed from gromos-rs in 0.0.67; repository
`yerkoescalona/newen-md`). The architectural invariants live in `architecture.md`.
The porting discipline lives in `porting.md`. What is proven or tested against a Lean model lives in
`spec/` and is explained in `docs/src/verification/index.md`.

**Work is tracked privately**, as draft items on a private GitHub Project, never as issues on this
public repository: an item is a deliverable, its checklist is the steps, and no item is ready
without a **Gate** — the objective, checkable condition for calling it done. Un-gated ideas,
including the long-horizon architectural bets, stay in the icebox. **History lives in
`CHANGELOG.md`**, whose item index resolves the stable ids (`1.5`, `2.8`, `3.9 A16`) that code
comments cite.

## When working in a crate, read its stage contract

| Crate | Layer | Stage contract |
|-------|-------|----------------|
| newen-logic | decision logic, proven (Aeneas) | `crates/newen-logic/.claude/CONTEXT.md` |
| newen-core | L0 data core | `crates/newen-core/.claude/CONTEXT.md` |
| newen-io | L0 I/O service | `crates/newen-io/.claude/CONTEXT.md` |
| newen-forces | L1 pure compute | `crates/newen-forces/.claude/CONTEXT.md` |
| newen-integrators | L2 steppers | `crates/newen-integrators/.claude/CONTEXT.md` |
| newen-md | L3 orchestration | `crates/newen-md/.claude/CONTEXT.md` |
| newen-run | L3 orchestration (library) | `crates/newen-run/.claude/CONTEXT.md` |
| newen-tools | L3 system builder | `crates/newen-tools/.claude/CONTEXT.md` |
| newen-analysis | L4 analysis facade | `crates/newen-analysis/.claude/CONTEXT.md` |
| pyo3-newen | L4 Python bindings | `crates/pyo3-newen/.claude/CONTEXT.md` |
| py-newen | L4+ Python package | `py-newen/.claude/CONTEXT.md` |

Each contract covers: the crate's job, what it consumes/produces, current status, key files,
and crate-specific rules that extend the workspace-wide guides. `crates/newen` is the facade that
re-exports the Rust crates; it has no contract.

newen-logic sits below every other crate (it depends on nothing; newen-core, newen-forces,
newen-integrators, newen-run and newen-md call it). Its contract: decisions only (integers and
booleans, no floats, plain data in and out, the Rust subset Aeneas translates), and a function
joins only with its Lean theory, model, bridge theorem, Aeneas refinement, test and `SPEC.toml`
entry.

## Prime directives

1. **The reference suite is the oracle.** A change that alters wired, tested output is a regression
   until proven a deliberate, documented decision. See `porting.md`.
2. **One core, no second codebase.** One `(Topology, State)` in newen-core; the MD engine and the
   analysis facade both *call* it. The moment analysis re-implements a physics term, gromos++ has been
   recreated. See `architecture.md`.
3. **Shape force/state/neighbor work as composable providers over shared data** — so QM/MM/ML and
   GPU become additive, not forks. See `architecture.md`.
4. **Don't port the bugs, don't diverge by accident.** Tests before wiring; derivation before
   transcription; grep before porting; document every deliberate quirk. See `porting.md`.
5. **One behaviour per feature.** newen computes the intended physics. Where gromosXX is wrong, the
   oracle is patched, not newen: no option, flag or recipe value reproduces a known defect, so the
   code stays self-consistent and its intent unambiguous. A corrected defect has a dossier id
   (D1–D59, never renumbered), a `_fixed` reference from the patched oracle, and its
   released-binary reproducer in `~/Development/gromos-divergences`. See `porting.md`, "One
   behaviour per feature".

## Reference sources (read-only)

| Source | Path |
|--------|------|
| gromosXX / md++, released (247c95c4) | `.local/gromosXX/md++/src` |
| gromosXX, patched oracle (247c95c4 + D1, D2, D16, D25, D32, D40, D41, LS1–LS5) | `.local/gromosXX-patched/md++` |
| gromos++ | `.local/gromosPlsPls/gromos++/src` (patched with D4, D5, D8: `.local/gromosPlsPls-patched`) |
| How the patched binaries were built | `crates/newen-md/tests/gromosXX_references/patches/BUILD.md` |
| GROMOS theory book | `.local/doc/gromos_book` |
| Tutorials | `.local/gromos_tutorial_livecoms/tutorial_files` |
| Force fields | `.local/gromosXX/forcefields` |

## How to verify

**Against gromosXX.** `cargo build --release --bin md`, then `cargo test -p newen-md --test
test_gromosXX_references` (135 references, none ignored; the per-system table is
`crates/newen-md/tests/gromosXX_references/README.md`), and `cargo test -p newen-md --test
gromos_check_suite` (gromosXX's own per-term oracle). A reference named `*_fixed` holds the output
of the patched oracle (the patch is named in its `input.toml`); every other one, of the released
binary. Ewald references were generated with `@develop` (gromosXX marks Ewald as under
development). See `crates/newen-md/.claude/CONTEXT.md` for how to add a reference test.
MPI: `LIBCLANG_PATH=/usr/lib/llvm-19/lib cargo build --release --features use-mpi --bin md`, then
`mpirun -np 2 target/release/md @topo … ` must reproduce the serial energies (`scripts/bench_mpi.py`, BENCHMARKING.md §6).

**Against the Lean spec** (`spec/`: `NewenTheory` the physics, `NewenSpec` the models of the code,
`NewenSpec/Extracted` the Aeneas translation, `vectors/` the Lean-generated test inputs;
`spec/SPEC.toml` lists 101 functions: 18 integer functions of newen-logic, proven by Aeneas, 80
float kernels tested bit for bit, 3 FFT paths within a tolerance):

```
cd spec && lake build              # every proof; must end without a `sorry` warning
cd spec && lake exe spec-vectors   # rewrite spec/vectors/*.json from the Lean models
spec/scripts/extract.sh            # after changing newen-logic: regenerate NewenSpec/Extracted/, then lake build
cargo test -p newen-logic          # schedule_spec, dof_spec, bath_spec, topology_spec
cargo test -p newen-core --test periodicity_spec
cargo test -p newen-integrators --test thermostat_spec --test dynamics_spec
cargo test -p newen-forces --test nonbonded_spec --test bonded_spec --test softcore_spec \
  --test perturbed_bonded_spec --test lattice_sum_spec --test ewald_spec --test p3m_spec --test spme_spec
```

The Lean package cache (Mathlib, Aeneas, about 7.8 GB) lives in `spec/.lake/packages` in this
checkout; `lake build` reuses it. The toolchain and the Aeneas/Charon install are in
`docs/src/verification/index.md`, "Running it".

The same behaviour from the Python side (both suites must stay green — they drive one builder):
`cd py-newen && uv run maturin develop --release && uv run pytest tests/ -q`, and the stubs
against the extension: `MYPYPATH=python uv run python -m mypy.stubtest newen.newen --allowlist
stubtest_allowlist_no_ml.txt`. `just lint` = clippy + the G6 drift greps (run the recipe body by
hand if `just` is not installed).

## Global decisions

- **f64 everywhere** (not f32)
- **gromosXX `@` CLI convention:** `@topo @conf @input @fin @trc @tre @trf @trv @verb ...`
- **Tolerances:** force_abs=1e-6, energy_rel=1e-8, position_abs=1e-9
- **CLI arg parsing:** `newen_io::args::Arguments` (the gromos++ reader: multi-valued `@key`s,
  absent/flag/valued distinguished) or clap `#[derive(Parser)]`, both over argv already put through
  `newen_io::normalize` (`@key → --key`, `@f argfile` expansion). No third parser, and no
  hand-rolled `env::args()` loop.
- **A program's body is a library function (3.10):** `<crate>::programs::<name>::run(argv,
  out)` with a `Program` record in the crate's registry; `src/bin/**` is `fn main() {
  newen_io::run_main(&PROGRAM) }` and nothing else. A program never reads the process arguments,
  never writes to the process stdout and never exits — `just lint` greps for all three (G11).
- **One run description (3.9, done 2026-08-29):** `newen_run::RunRecipe` is the only
  description of a run; `.imd` files and `py-newen` objects are front-ends to it, both built by
  `newen-run` (`prepare_system → build_plan → validate_plan → instantiate → start`). No second
  IMD→algorithm builder, no second IMD reader — `just lint` greps for both (G6). Adding a force term
  = a `TermSpec` variant + its registry arms + one `instantiate` arm, nothing in the binding.
- **Electrostatics (NONBONDED NLRELE):** −1 and 1 reaction field (−1 is read as 1, as gromosXX
  does), 0 none, 2 Ewald, 3 P3M (PhastFT, grids must be powers of two), 5 SPME (newen only,
  experimental). 4, gromosXX's shifted reaction field, is not ported and refused for every run.
  The lattice sums carry
  gromosXX's feature locks (atomic cutoff required; no vacuum, truncated octahedron, several
  energy groups, perturbation, EDS, A-EDS, GaMD, QM/ML terms, conjugate gradients, LINCS/SETTLE
  for the solvent); `newen-run`'s contract lists every refusal.
- **Boxes:** NTB 2 reads GENBOX's angles (a triclinic cell) and writes them back; NTB 1 with an
  angle ≠ 90° and any non-zero Euler angle are refused, and so is a cutoff that reaches half the
  smallest cell height.
- **Verification gates promotion, not entry.** New methods may run in an experimental zone with
  self-consistency tests only, marked as experimental in the run's output and plan (SPME is the
  first); code becomes a default, or goes into a paper or tutorial, only with its Lean model,
  theorem and test. Where verified code goes: `architecture.md`, "Where verified code goes".
- **Doc style:** Rust → KaTeX + `[^label]` footnotes; Python → NumPy docstrings + `.. math::`
- **On commit:** update CHANGELOG.md and Cargo.toml version.

## Working in this checkout

- **No worktrees** (neither `git worktree add` nor isolated agents): one checkout, one
  `target/`, one `spec/.lake`. Each extra copy carries its own build caches (several GB) and
  hides uncommitted work.
- **Agents run one after another** in this checkout, on the current branch; when two must run at
  once, they get disjoint files and each is told about the other.
- **Commit each step as it passes, by pathspec** (`git commit -- <paths>`), so a step's commit holds
  only its own files.

## Before committing

Always run these two checks before any `git commit`. If either fails, fix it first — do not commit broken format or type errors.

```
cargo fmt          # fix formatting in-place
cargo check --tests  # type check workspace including test code (matches CI coverage)
```

CI runs both but with `continue-on-error: true` (advisory only). The commit is the real gate.
CI must never auto-modify code — only check and report.
