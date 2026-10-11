# gromos-analysis — stage contract

## Job
L4 facade. GROMOS-style analysis programs as a thin layer over the shared engine core.
Zero physics re-implementation — all physics comes from `gromos-forces`.

## Inputs (consumes from)
`gromos-core` + `gromos-io` + `gromos-forces` (for `single_point_energy`).

## Outputs
Analysis binaries producing results comparable to GROMOS analysis tools.

## Status — implemented
- `atominfo` ✓ — reads topology + AtomSpecifier, prints TITLE+ATOMS; verified against GROMOS atominfo
- `rmsd` ✓ — Kabsch fit (`fit.rs`), @atomspec, @ref, @pbc, @nofit
- `nhoparam` ✓ — N-H order parameters S², rotational fit, window averaging, `ee()`
- `frameout` ✓ — full GROMOS feature parity: @pbc gather, @include SOLUTE/SOLVENT/ALL,
  @ref/@atomsfit rotational fit, @spec ALL/EVERY/SPEC, @time range, cnf/pdb/trc output, @single;
  all I/O through gromos-io; 8 integration tests
- `ener` ✓ — energy recalculation from trajectory via `gromos-forces::energy::single_point_energy`
- `ext_ti_ana` ✓ — trapezoidal ΔG ± ee() from `.trg` files
- `bar` ✓ — BAR iteration (numerically stable log-sum-exp), bootstrap error
- `ext_ti_merge` ✓ — λ interpolation, trapezoidal ΔG
- `fit.rs` ✓ — Kabsch rotational fit (Horn 1987 quaternion); `superimpose()`, `rmsd()`, 7 unit tests

- **Ports of gromos++ programs with reference tests (2026-08-30, `tests/gromospp_references.rs`,
  fixtures + exact gromos++ commands in `tests/data/README.md`):** `tser` (rewritten: property
  specifiers, distributions), `mdf`, `dg_ener`, `dfmult`, `matrix_overlap`. In gromos-tools:
  `explode`, `duplicate`, `pt_top`.

## Shared machinery (use it; do not re-implement per binary)
- `gromos_io::args::Arguments` — gromos++ `Arguments` (`@key` multi-valued, `count`, `@f` files).
- `gromos_io::pbc::Pbc` — `@pbc <r|v> [cog|nog]`, `periodicity(box)`, `gather()` = gromos++ `coggather`.
- `gromos_analysis::time::Time` — `@time t0 dt` or the frame's time; `Time::fmt` = gromos++ layout.
- `gromos_analysis::property` — `d/a/t/tp/o/op/pr/pa` specifiers, `PropertyContainer`,
  `atoms_to_string` (gromos++ `AtomSpecifier::toString`), `ordered_atoms` (order-preserving).
- `gromos_analysis::distribution` — gromos++ `Distribution`; `cpp_g` = C++ default float format.
- `gromos_analysis::lnexp` — log-exponential averages / covariance / statistical inefficiency.
- `gromos_io::table::read_columns` — numeric column files; `gromos_io::coordinate::format_g96`.
- `gromos_core::Stat::ee_strict` — gromos++'s error estimate (`nan` below 200 values).
- gromos-core `nearest_image(ri, rj)` returns the minimum-image **vector** ri − rj; gromos++'s
  `nearestImage(ref, pos)` is the image **position** `ref − that vector`. Every port must mind this.

## Status — stubs / parked
- `visco`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy` — stubs; parked
- Ported with reference tests (batch 2): `bilayer_dist`, `bilayer_oparam`, `jval`, `edyn`
  (eigenvector projections — gromos++'s are a bug, see CHANGELOG 0.0.40), `gca`.
- Ported with reference tests (batch 3): `eds_update_1`, `eds_update_2` (`src/eds.rs` shared),
  `jepot` (`gromos_io::jvalue` readers, `time::TimeSpec` shared with `jval`), `pocket` (no
  gathering — gromos++ never applies `@pbc` there), `dfgrid` (gromos++'s open-set Dijkstra kept
  for identical tie-breaking; `cnf` output, real atoms only); `make_sasa_top` lives in
  gromos-tools.
- Not ported, with reasons (CHANGELOG 0.0.41): `cos_dipole`, `cos_epsilon` (COS special
  trajectories), `prep_bb` (interactive), `prep_xray_le` (X-ray).
- `dGslv_pbsolv` (1.8): the finite-difference side (`@schemeELEC LS`, coordinate and PQR input)
  over `src/pb/` (gromos++ `PB_Parameters`, `FDPoissonBoltzmann`, the two ICCG solvers, the
  atomic radii); stdout and the `@verbose` log byte-identical to a patched gromos++ build on four
  cases (`tests/data/dgslv_pbsolv/README.md`). `@schemeELEC RF` (the FFT solver) is refused until
  it is ported. `pb::dgslv_fd` / `pb::solve_fd` are the library API (potentials and ΔG as numbers).

## Shape — 26 of the 83 programs are library functions
The programs converted so far follow gromos-tools: the body is
`src/programs/<category>/<name>.rs` (`pub const PROGRAM` + `pub fn run(argv, out)`), the binary is
a one-line `gromos_io::run_main(&PROGRAM)`, and the two trees carry the same categories as the
`[[bin]]` sections of `Cargo.toml`. `programs/mod.rs` re-exports every program flat, so a caller
writes `programs::rmsd`; `super::` inside a body is the category, so the shared helper is
`crate::programs::io_err`. `PROGRAMS` is the registry a Python binding can bind.

Converted (3.10 step 6): the 18 that already parsed with `Arguments` — `bilayer_dist`,
`bilayer_oparam`, `dfgrid`, `dfmult`, `dg_ener`, `edyn`, `eds_update_1`, `eds_update_2`, `ene_ana`,
`gca`, `jepot`, `jval`, `matrix_overlap`, `mdf`, `pocket`, `rmsd`, `rmsf`, `tser` — and the 7 that
read `gromos_args()` and walked the vector by hand, whose parsers became `Arguments` here:
`atominfo`, `bar`, `ener`, `ext_ti_ana`, `ext_ti_merge`, `frameout`, `nhoparam`. Written as a
library function from the start: `dGslv_pbsolv` (module `dgslv_pbsolv`, the file keeps the
gromos++ name because the registry checks `programs/<category>/<bin name>.rs`).

**The other 57 still have their compute in `main()`** (3.10 step 8, by directory, in batches). They
are ordinary binaries with hand-rolled `env::args()` parsers; nothing here applies to them until
they are converted.

## Key files
```
src/fit.rs                          — Kabsch fit, superimpose(), rmsd()
src/pb/                             — gromos++ pb library: PB_Parameters, FD Poisson–Boltzmann,
                                      ICCG solvers, PB radii (dGslv_pbsolv)
src/programs/mod.rs                 — the categories, the flat re-exports, PROGRAMS, io_err
src/programs/<cat>/<name>.rs        — one module per converted program: PROGRAM, run(argv, out)
src/bin/<cat>/<name>.rs             — a converted program's one-line run_main shim; the other 57
                                      binaries still hold their own compute
src/bin/dynamics/, clustering/, xray/ — no program converted yet
tests/library_parity.rs             — library run == binary stdout, for every registered program
tests/gromospp_references.rs        — the 19 against gromos++ output (tests/data/README.md)
tests/test_frameout.rs              — 8 frameout integration tests
```

## Crate-specific rules
- **A converted program is a library function**, never a binary-only `main`: body in
  `src/programs/<category>/`, `pub mod` in that category's `mod.rs`, the flat re-export and a
  `PROGRAMS` entry in `src/programs/mod.rs`, the shim in `src/bin/<category>/`, and a
  `library_parity.rs` case. `every_program_is_a_bin`, `program_and_bin_agree_on_category` and
  `every_registered_program_has_a_case` fail if one of those is missing.
- **A program never owns the process:** no `gromos_args()`, no `println!`, no `process::exit`.
  argv arrives normalised, output goes to `out`, a failure is an `Err`. A program that writes its
  results to files gets `_out` and a `library_parity` case on the failure path.
- **Zero physics re-implementation.** All force/energy calls go through `gromos-forces`.
- **All I/O through gromos-io.** No hand-rolled `BufWriter`/`File` in binary entry points;
  use `write_g96`, `write_pdb_positions`, `TrajectoryWriter::write_trc_frame`.
- **Kabsch fit lives here** (`src/fit.rs`), not in gromos-core — it is analysis-only.
- **PBC gathering lives in gromos-core** (`src/gather.rs`) — shared with MD engine.
