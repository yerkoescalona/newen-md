# newen — workspace guide

newen ports **gromosXX/md++** (MD engine) + **gromosPlsPls/gromos++** (analysis + tools) into
one Rust workspace with a single shared core. The architectural invariants live in `architecture.md`.
The porting discipline lives in `porting.md`.

**Work is tracked on the project board** (`Icebox` → `Ready` → `Doing` → `Done`): a card is a
deliverable, its checklist is the steps, and no card reaches `Ready` without a **Gate** — the
objective, checkable condition for calling it done. `Icebox` holds the un-gated ideas, including the
long-horizon architectural bets. **History lives in `CHANGELOG.md`**, whose item index resolves the
stable ids (`1.5`, `2.8`, `3.9 A16`) that code comments cite.

## When working in a crate, read its stage contract

| Crate | Layer | Stage contract |
|-------|-------|----------------|
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
and crate-specific rules that extend the workspace-wide guides.

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
   code stays self-consistent and its intent unambiguous. See `porting.md`, "One behaviour per
   feature".

## Reference sources (read-only)

| Source | Path |
|--------|------|
| gromosXX / md++ | `.local/gromosXX/md++/src` |
| gromos++ | `.local/gromosPlsPls/gromos++/src` |
| GROMOS theory book | `.local/doc/gromos_book` |
| Tutorials | `.local/gromos_tutorial_livecoms/tutorial_files` |
| Force fields | `.local/gromosXX/forcefields` |

## How to verify

`cargo build --release --bin md`, then `cargo test -p newen-md --test test_gromosXX_references`,
and `cargo test -p newen-md --test gromos_check_suite` (gromosXX's own per-term oracle).
MPI: `LIBCLANG_PATH=/usr/lib/llvm-19/lib cargo build --release --features use-mpi --bin md`, then
`mpirun -np 2 target/release/md @topo … ` must reproduce the serial energies (`scripts/bench_mpi.py`, BENCHMARKING.md §6).
See `crates/newen-md/.claude/CONTEXT.md` for how to add a reference test.

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
- **Doc style:** Rust → KaTeX + `[^label]` footnotes; Python → NumPy docstrings + `.. math::`
- **On commit:** update CHANGELOG.md and Cargo.toml version.

## Before committing

Always run these two checks before any `git commit`. If either fails, fix it first — do not commit broken format or type errors.

```
cargo fmt          # fix formatting in-place
cargo check --tests  # type check workspace including test code (matches CI coverage)
```

CI runs both but with `continue-on-error: true` (advisory only). The commit is the real gate.
CI must never auto-modify code — only check and report.
