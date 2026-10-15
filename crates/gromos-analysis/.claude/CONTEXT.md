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
- `ext_ti_ana`, `bar`, `ext_ti_merge` ✓ verified (1.8) — the gromos++ programs, byte for byte
  on the oracle cases (`tests/data/README.md`: A1–A5 and the five-window CH4 chain): PRECALCLAM
  reweighting, BAR and its error, the merge with its error integral. Three gromos++ 9cd9f6b3
  defects (BAR error ÷ (ln⟨f⟩)² twice, merge error without point N−2, dihedral ∂H/∂λ without
  dΛ/dλ) are corrected by default and reproduced under `@parity gromospp`
  (`programs::free_energy::Parity`); the corrected output is checked against patched gromos++
  (`A2_fixed`, `chain_fixed`). Deliberate differences in each module's header. Library API:
  `programs::ext_ti_ana::{Options, analyse, write_files}`, `programs::bar::{read_bar_file, bar,
  do_bar}`, `programs::ext_ti_merge::{merge, integrate, integrate_error}`; Python
  `gromos.analysis.ext_ti_ana/bar/ext_ti_merge`.
- `mbar` (no gromos++ counterpart, 1.8) — MBAR on `ext_ti_ana @bar_data` files, against pymbar
  (`tests/mbar_reference.rs`); `programs::mbar::mbar`, Python `gromos.analysis.mbar`.
- `ti::trapezoid_from_trg` — the plain trapezoid TI over `.trg` totals the old gromos-rs
  `ext_ti_ana` was (`scripts/ti_ch4.py --estimator ti`, `gromos.analysis.ti_trapezoid`);
  `TiWindow::from_samples` is the same window over in-memory dH/dλ
  (`gromos.free_energy.FreeEnergyWindows.ti`, 0.0.59).
- `energy_sequence::EnergySequence` — the `@en_files`/`@fr_files` frame stream `ene_ana` and
  `ext_ti_ana` share.
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
- `dGslv_pbsolv` (1.8) ✓ verified: over `src/pb/` (gromos++ `PB_Parameters`,
  `FDPoissonBoltzmann`, the two ICCG solvers, the atomic radii; for `@schemeELEC RF` the FFT
  Poisson solver `FFTPoisson`/`FFTPoissonIterator` with its grid, occupancy, hat charges, LS and
  RF vacuum fields and tensors, on rustfft). Against a patched gromos++ build on six cases
  (`tests/data/dgslv_pbsolv/README.md`): everything from the FD solves byte-identical, the FFT
  numbers within the last printed digit (byte-identical where they were made). `@epsRF 1` (SC) is
  refused: gromos++ aborts there. Library API: `pb::dgslv_fd` / `pb::solve_fd`, `pb::dgslv_fft` /
  `pb::solve_fft`, and the whole program as `programs::dgslv_pbsolv::evaluate` (bound in Python
  as `gromos.analysis.dgslv_pbsolv`).

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
