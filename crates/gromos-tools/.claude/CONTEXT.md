# gromos-tools — stage contract

## Job
L3 system builder. System construction binaries: the gromos++-style tools for preparing simulation inputs.

## Inputs (consumes from)
`gromos-core` + `gromos-io` (+ `gromos-forces` for single-point energy — not yet).

## Outputs
Standalone binaries: make_top, com_top, check_top, pdb2g96, sim_box, ion, mk_script,
make_pt_top, prep_posres, build_box.

## Status (audit 2026-08-30, every tool run on real inputs against gromos++ where a counterpart exists)

**Verified against gromos++ output** (reference tests under `tests/`, fixtures + commands in
`tests/data/<tool>/README.md`; gromos++ = `.local/gromosPlsPls`):

| tool | what is compared | test |
|------|------------------|------|
| `make_top` | topology reads back identical: atoms, exclusions, 1-4, every bonded term and parameter, solvent — peptide `NH3+ ALA GLY COO-` and a two-MTB methanol; gromos++'s own reader accepts ours | `topology_references.rs` |
| `com_top` | same, aladip + CH4 | `topology_references.rs` |
| `red_top` | same, atoms 1–6 of aladip | `topology_references.rs` |
| `make_pt_top` | `PERTATOMPARAM` reads back identical; prints each perturbed molecule's net charge A/B (1.8, 0.0.59; `pt_top` puts it in its title). Typed verbs (0.0.64): `make_pt_top::make_pt_top` (two topologies → `.ptp` text + net charges), `make_pt_top::mutate` (state B = state A with new charges/IAC/masses on chosen atoms), `pt_top::pt_top` (state B of a topology); `gromos.free_energy.Perturbation` calls them | `topology_references.rs` |
| `sim_box` | box, solvent count, every position (352 SPC around methanol) | `sim_box_reference.rs` |
| `copy_box` | every position and label, `@dir x` | `box_tool_references.rs` |
| `inbox` | every position | `box_tool_references.rs` |
| `ion` | the water gromos++ replaces (`@potential`, its energy rule: all atoms of other molecules, shifted Coulomb) | `box_tool_references.rs` |
| `build_box` | positions and box (methanol, `@nsm 3 @dens 790`) — checked by hand, identical | — |
| `check_top` | runs; gromos++ additionally cross-checks dihedral types against the building blocks (it flags one in the tutorial's methanol topology that we accept) | — |
| `atominfo` | gromos++'s columns, 1-based IAC | — |

**Run, but not gromos++'s program** — same name, simplified or different semantics; results were
not cross-checked. Treat as utilities, not ports: `gch` (checks constraint distances; gromos++
regenerates hydrogen coordinates), `check_box` (box statistics; gromos++ reports close atom pairs),
`bin_box` (bins a trajectory; gromos++ builds a binary mixture), `unify_box` (rescales to a box;
gromos++ converts between box types), `pert_top` (compares two topologies; gromos++ modifies
selected atoms), `con_top` (concatenates; gromos++ converts force-field parameters), `link_top`,
`addvirt_top`, `prep_eds`, `prep_noe`, `amber2gromos`, `ran_box`, `ran_solvation` (random
placement with other argument models), `mk_script` (needs `@template`).

**No gromos++ counterpart:** `traj2pdb`, `top2build`, `prep_posres`, `prep_xray` (gromos++'s
needs clipper), `mpi_scaling`.

Conventions: `make_top`, `com_top`, `check_top`, `pdb2g96`, `ion`, `build_box`, `prep_posres` parse
with clap, the rest with `gromos_io::args::Arguments`; both see argv already normalised by
`gromos_io::normalize` (`@flag` → `--flag`, `@f` expanded), so every program takes either form and
every program accepts `@f`. An unknown `@key` is an error (0.0.49; the hand-rolled parsers used to
ignore them).

Adding a reference test: put gromos++'s output under `tests/data/<tool>/` with a README that records
the exact command, and compare in a `tests/*_references.rs` (`read_topology_file` for topologies —
formatting-independent — or the `parse` helper of `box_tool_references.rs` for coordinates).

## Shape (3.10, 0.0.49)
Every program's body is `src/programs/<category>/<name>.rs`: `pub const PROGRAM` (name, accepted
`@key`s, usage, status, `run`) and `pub fn run(argv, out) -> Result<(), String>`. `src/bin/**` is a
one-line `gromos_io::run_main(&PROGRAM)`. The four categories — `topology`, `box`, `conversion`,
`utilities` — are the same on both sides and in the `[[bin]]` sections of `Cargo.toml`; the module
is `r#box`, `box` being a keyword. `programs/mod.rs` re-exports every program flat, so a caller
writes `programs::make_top`, not `programs::topology::make_top`, and `super::` inside a body is the
category, so shared helpers are `crate::programs::{io_err, parse_clap}`. The preparation chain adds
a typed verb beside `run` — `make_top`, `com_top::com_top`, `red_top::reduce`, `pdb2g96`,
`sim_box::solvate`, `build_box::build_box`, `ion::add_ions`, `check_top::check_top` — returning
`ParsedTopology` / `LabeledCoordinateData` / a report, with `run` = parse args + read files + verb
+ serialise.
`PROGRAMS` is the registry `pyo3-gromos` binds as `gromos.tools`.

## Key files
```
src/lib.rs          — the registry (`PROGRAMS`, `find`)
src/programs/<cat>/ — one module per program: PROGRAM, run(argv, out), and the typed verb
src/programs/mod.rs — the four categories, the flat re-exports, PROGRAMS, io_err, parse_clap
src/bin/<cat>/      — one-line `run_main(&PROGRAM)` shims; nothing else belongs here
tests/library_parity.rs — library run == binary stdout, byte for byte, for every program
tests/typed_verbs.rs    — the typed verb's result, serialised, == what run writes
tests/              — reference tests against gromos++ output (`tests/data/<tool>/README.md`)
```

## Crate-specific rules
- **Every program is a library function** (not a binary-only `main`): `src/programs/<category>/<name>.rs`
  holds `pub const PROGRAM` + `pub fn run(argv, out)`, the prep chain adds a typed verb, and
  `src/bin/<category>/<name>.rs` is a one-line `run_main(&PROGRAM)` — 3.10. Do not add a new binary
  with its compute in `main()`.
- **A new program lands in four places, all in the same category:** the body under
  `src/programs/<category>/`, a `pub mod` in that category's `mod.rs` and the flat re-export plus a
  `PROGRAMS` entry in `src/programs/mod.rs`, the shim under `src/bin/<category>/`, and a `[[bin]]`
  in `Cargo.toml`. `every_program_is_a_bin` catches a missing `[[bin]]` and
  `program_and_bin_agree_on_category` catches a body and its shim filed under different categories
  — the drift that put a flat `programs/` next to a grouped `bin/`.
- **File parsing always delegates to gromos-io.** Never duplicate format code in tool binaries.
- **AtomSpecifier:** build the gromos++ grammar as a facade over gromos-core primitives (P2.1), not as a standalone reimplementation.
