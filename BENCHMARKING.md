# BENCHMARKING.md — gromos-rs vs gromosXX, one core → threads → MPI

A checklist-driven plan and, below Phase 2, the measured results so far (2026-08-28/29).
Every `- [ ]` is a step still to do, and every **Assumption** is something the numbers silently
depend on. Do not report a number whose assumptions are unchecked.

The plan has four phases with a gate between each. Do them in order — a thread-scaling number is
meaningless if single-core parity is unproven, and an MPI number is meaningless if the MPI binary
does not exist (today it doesn't, on either side).

```
Phase 0  prerequisites & fixes  ──gate──▶  Phase 1  parity (same physics?)
                                                    │
                                                  gate
                                                    ▼
Phase 4  MPI (measured §6.3)  ◀──gate──  Phase 3  threads  ◀──gate──  Phase 2  one core vs one core
```

Phase X (spin-off, future): the integrated GPU — a probe with numbers is at the end of this file;
not on the critical path.

---

## At a glance: what was changed, what is only suggested

**Applied** (all in the tree, every one verified against the 37-system gromosXX reference
suite; details and measurements in "Optimisations applied so far" under Phase 2):

| # | Change | Where | Effect |
|---|--------|-------|--------|
| 1 | Long-range solute uses the charge-group-grouped kernel and a reused buffer | `forcefield.rs` | water_1000 long-range 1.59 → 0.60 s |
| 2 | Cell-list neighbour scan without a hash set | `pairlist.rs` | none measurable (kept, simpler) |
| 3 | `NTWX/NTWE/NTPR = 0` no longer panics (divide by zero) | `md.rs` | made benchmarking possible |
| 4 | Timing is opt-in; `md @verb 1` prints a gromosXX-style `TIMING` block | `algorithm.rs`, `forcefield.rs`, `md.rs` | removed 2 clock reads/algorithm/step |
| 5 | Cell list sized from IMD `SIZE` (was parsed, never used) + true inter-cell distance pruning | `pairlist.rs`, `md.rs` | water_1000 pairlist 3.15 → 2.15 s |
| 6 | Minimum image rounds ties-to-even (`rint`, as gromosXX) | `math.rs` | water_216 kernel 1.65 → 1.15 s |
| 7 | Parallel kernels: one chunk per thread, not per 1024 pairs; `par_chunks(0)` panic fixed | `innerloops.rs` | ~7 % |
| 8 | Four pairs per SIMD register; one shared pair-physics routine for all kernels | `innerloops.rs` | water_216 kernel 1.15 → 0.90 s |
| 9 | Solvent kernel batched the same way | `innerloops.rs` | ch4 short-range solvent 4.0 vs XX 5.8 s |
| 10 | Pairlist stored as `u32`; per-step copy deleted | `pairlist.rs` + 5 consumers | ch4_water_fep −2.5 s |
| 11 | Exclusion test: one sorted lookup + `last < j` early-out; `solvate` bug fixed | `topology.rs`, io reader | water_216 pairlist 0.28 → 0.19 s |
| 12 | Long-range solvent as sentinel molecule pairs when provably exact | `pairlist.rs`, `forcefield.rs` | 24k box 1.01× → 1.08× |
| 13 | Long-range kernels skip the virial when nothing reads it | `forcefield.rs` | 24k box 1.08× → 1.19× |
| 14 | Pairlist construction in parallel, both algorithms (per-range private lists, appended in order) | `pairlist.rs` | water_1000 8-thread pairlist 1.89 → 0.47 s |
| 15 | Long-range block uses the parallel kernels the short-range block already used | `forcefield.rs` | water_1000 8-thread long-range 0.35 → 0.09 s; S(8) 1.8× → 3.9× with 14 |
| 16 | Charge-group pair-group build in parallel (per-chunk runs, stitched at boundaries) | `forcefield.rs` | water_1000 8-thread 0.20 → 0.07 s; step 1.74 → 1.35 s |
| — | `PAIRLIST TYPE atomic` honoured; `ALGORITHM grid` → cell list (was O(N²)) | `pairlist.rs`, callers | fidelity |
| — | Stale harness paths fixed; `bench_engines.py` (interleaved, n/mean/std); `make_solvated_box.py` | `scripts/`, `benches/` | infrastructure |

**Tried and reverted** (documented under Phase 2 so they are not retried blind): in-batch
vectorised minimum image (slower), CSR per-atom pairlist layout (no gain), an invented
`CELLS_PER_CUTOFF` constant (replaced by the real IMD parameter).

**Suggested, not done:**

- Phase 4 MPI: done 2026-08-30 — `md --features use-mpi` under `mpirun`, np=1 bit-identical to serial, 1.6× at np=4 on 8k–24k atoms (§6.3)
  re-export). Phase 3 is done through 8 threads plus the SMT point; the 81 000-atom box
  (`make_solvated_box.py --replicate 3`) has not been run.
- Threads: the remaining serial phases are SHAKE, integration and the thermostat (both engines);
  `ParallelPairlistAlgorithm` in `pairlist.rs` is now redundant with fix 14 and could be removed.
- A thread-count flag for `md` (today only `RAYON_NUM_THREADS`); a Rust
  build without `target-cpu=native` as the A1 control; `performance` governor (needs sudo).
- The remaining 10 % on `ch4_water_fep` (NSNB=1 regime) and the Rust-side repeat variance there.
- Porting gromosXX's `Extended_Grid` plane sweep — judged unnecessary (same pair set as the cell list).
- Phase X (spin-off, future): an integrated-GPU force provider — feasible in f64 via wgpu/Vulkan,
  ~2–4× the 8 cores on arithmetic, capped by serial SHAKE. Parallel SHAKE/SETTLE on the CPU
  should come first regardless. See the last section.

---

## 0. Ground truth (state found on 2026-08-28)

| Item | Found | Consequence |
|------|-------|-------------|
| Machine | AMD Ryzen 7 PRO 8700GE, 8 physical cores / 16 threads (SMT on), 1 socket | Thread scaling tops out at 8 real cores; CPUs `n` and `n+8` are siblings |
| CPU governor | `powersave` | Timings will be noisy/low until set to `performance` |
| gromosXX binary | `.local/gromosXX/md++/BUILD/program/md`, v1.6.0, built 2026-08-22 via `../configure` **with no options** | Serial only: no OpenMP, no MPI. `-O3 -DNDEBUG`, **no `-march=native`** |
| gromosXX `md_mpi` | Exists in the same BUILD dir but links **no libmpi** and has **no `MPI_Init` symbol** | It is a stub build (`XXMPI` undefined). Cannot be used for MPI runs |
| gromos-rs `md` | `cargo build --release --bin md` → `target/release/md`; `opt-level=3`, `lto="fat"`, `codegen-units=1`, `panic="abort"`, **`-C target-cpu=native`** (`.cargo/config.toml`) | Rust gets AVX2/Zen4 codegen that gromosXX currently does not — an unfair default |
| gromos-rs threads | No `@nt`-style flag. Rayon global pool, never configured → defaults to **all 16 logical CPUs** | Control with `RAYON_NUM_THREADS`; an unset run is a 16-thread run |
| gromos-rs `md_mpi` | `crates/gromos-md/src/bin/md_mpi.rs`, feature `use-mpi`. **Does not build**: (a) `mpi-sys` needs `libclang` (bindgen) — missing; (b) it uses `MpiControl`, which only exists in `gromos-integrators/src/mpi.rs` behind `#[cfg(feature = "use-mpi")]` — a feature that crate never declares, so the module is always empty; (c) the `gromos` facade has no `pub mod mpi` re-export | Rust MPI is blocked by build plumbing before any physics question |
| gromos-rs MPI code | `broadcast_positions`/`reduce_forces` cast `&[Vec3]` (f64) to `&[f32]` | Data corruption the moment the gate is fixed. Must be fixed *before* any MPI run |
| MPI toolchain | Open MPI 5.0.7 (`mpirun`, `mpicc` present) | OK for both sides once binaries exist |
| Profiling tools | `taskset` present; `hyperfine`, `perf`, `numactl` **absent** | Pinning works; use a shell loop for repeats, install `perf` only if a breakdown is wanted |
| Reference systems | `crates/gromos-md/tests/gromosXX_references/` — 37 passing parity tests at NSTLIM=100 | Already-validated inputs; but 100 steps is too short to time |
| Existing bench harness | `crates/gromos-md/benches/md_bench.rs` read `tests/GROMOS_references` (dir is `gromosXX_references`); `scripts/benchmark.sh` did `cd gromos-rs` from inside the repo; `scripts/regen_gromosXX_references.py` used `build/` (dir is `BUILD/`) — **all three fixed 2026-08-28**; `scripts/bench_engines.py` is the engine-vs-engine harness | Benchmarks are automated now |
| Pairlists | gromos-rs has `StandardPairlistAlgorithm` (O(N²)), `CellListPairlistAlgorithm`, `ParallelPairlistAlgorithm`; gromosXX has `standard`, `grid_cell`, … | Both sides must be forced to the *same* `PAIRLIST ALGORITHM` per run |

---

## 1. Assumptions — check each before trusting any number

Each assumption says what it silently affects and how to verify it. Tick the box only after verifying on this machine.

### A. Fairness of the comparison

- [x] **A1 — Same compiler aggressiveness.** Rust builds with `target-cpu=native`; gromosXX is built with plain `-O3`. Either rebuild gromosXX with `-march=native` (preferred — see Phase 0) or build Rust once without `target-cpu=native` as a control. *Verify:* `grep MY_CXXFLAGS <build>/Makefile` shows `-march=native`; and/or a Rust control build with `RUSTFLAGS=""`.
- [x] **A2 — Same precision.** Both are f64 (`overview.md`: "f64 everywhere"; gromosXX `input.toml`: "double precision"). *Verify:* `md @version` on gromosXX reports double; no f32 in Rust physics paths (audit 2026-08-28 confirmed f32 only in vector aliases, GPU, MPI, FFI).
- [ ] **A3 — Same input file, same algorithm choices.** Both binaries read the same `.imd`. In particular `PAIRLIST ALGORITHM` (`standard` vs `grid_cell`), `NSNB` (pairlist update frequency), cutoffs, constraint algorithm, and `WRITETRAJ` frequencies must be identical *and* honoured by both. *Verify:* diff the two engines' echoed parameters in their log output; the 37-system parity suite proves parsing agreement at 100 steps.
- [x] **A4 — Same amount of I/O.** Trajectory/energy writing is I/O-bound and would dominate short runs. Set `NTWX=0 NTWV=0 NTWF=0 NTWE=0 NTWB=0` (or identical non-zero values on both). *Verify:* output directories after the run contain only the final configuration.
- [ ] **A5 — Same work per step.** Pairlist rebuild frequency (`NSNB`) is the biggest per-step cost lever; if one engine rebuilds every step and the other every 5, the comparison is void. Also SHAKE/SETTLE iteration counts, and whether virial/pressure is computed (`PRESSURESCALE`, `NTB`). *Verify:* same `.imd`, and compare the detailed timing breakdowns (gromosXX `DETAILED_TIMING` block in the `.omd`; Rust `@verb 2` `Timer` lines).
- [ ] **A6 — Startup cost is excluded or amortised.** Both engines parse topology/configuration at startup (Rust logs "init" time separately at `@verb 1`). Either subtract it or run enough steps (≥ 2000) that it is < 2 % of wall time. *Verify:* run NSTLIM=100 and NSTLIM=2000; per-step cost should be constant.

### B. The machine

- [ ] **A7 — Frequency is stable.** *(Partly mitigated: `bench_engines.py` now interleaves engines per repeat, A, B, A, B, so drift hits both equally; a block-ordered run on 2026-08-29 showed 9–11 % spread on whichever engine ran first after a compile.)* Governor is `powersave` today; Zen boost clocks vary with thermal state and with how many cores are active (a 1-thread run boosts higher than an 8-thread run, which *deflates* measured parallel speedup). *Verify:* set `performance`, run the same 1-core benchmark 5× and require < 3 % spread; optionally disable boost for scaling runs (`/sys/devices/system/cpu/cpufreq/boost`).
- [ ] **A8 — Pinning is real.** Without `taskset` the scheduler may place two threads on SMT siblings (`n` and `n+8`) and halve their throughput. *Verify:* `taskset -cp <pid>` during the run; for k threads pin to `0-(k-1)` (physical cores) — never mix in 8–15 unless deliberately measuring SMT.
- [ ] **A9 — The machine is otherwise idle.** IDE indexers, `cargo` daemons, browser — **and other Claude Code sessions**: a run on 2026-08-28 was invalidated by an unrelated session's Python job taking 648 % CPU, inflating both engines by ~45 % and blowing the spread out to 28–65 %. *Verify:* `uptime` load average < 0.5 **and** `ps -eo pcpu,args --sort=-pcpu | head` shows nothing large before each batch. `scripts/bench_engines.py` reports spread; treat > 5 % as "redo".
- [x] **A10 — Threads ≠ cores.** *(Measured: 16 threads is no faster than 8 for gromos-rs and 60–70 % slower for gromosXX-OpenMP.)* `RAYON_NUM_THREADS=16` uses SMT siblings; the honest ceiling for a scaling curve is 8. Report 16 as a separate "SMT" point, not as part of the curve.

### C. What is actually parallel (so scaling numbers are interpreted correctly)

- [ ] **A11 — Rust thread parallelism covers only part of the step.** Confirmed in code: the LJ/CRF innerloop (`innerloops.rs`) splits charge-group pair groups across Rayon (`par_chunks`) **only when there are ≥ 64 groups** — tiny systems (`aladip_solvated`, 72 atoms) run serially regardless of thread count. `ParallelPairlistAlgorithm` exists but whether `md` selects it is unverified. Bonded terms, constraints, integration, thermostats: serial. *Verify:* `grep -rn "par_iter\|par_chunks" crates/*/src` and check which are reachable from `md.rs`; expect Amdahl-limited scaling.
- [ ] **A12 — gromosXX OpenMP parallelism is also nonbonded-centric.** `omp_nonbonded_interaction.cc` plus a few `#pragma omp` in constraints/integration (54 pragmas total). Comparable scope to A11 — good for a fair threads-vs-threads comparison, but only after an `--enable-openmp` rebuild.
- [ ] **A13 — MPI in both codebases is a *nonbonded master/worker* design, not domain decomposition.** gromosXX `mpi_nonbonded_master.cc` broadcasts positions and reduces forces every step; the Rust port mirrors it. Expect communication to dominate on one node for small systems; MPI is only worth measuring on ≥ 3000-atom systems. MPI is deliberately ranked last (board card V8).

### D. The systems

- [x] **A14 — The reference systems are too small to show scaling.** *(24 000-atom box built 2026-08-29 with `sim_box` via `scripts/make_solvated_box.py --replicate 2` → `bench/systems/ch4_water_7975`.)* 648 atoms (`water_216_box`) fits in L2; 3000 atoms (`water_1000_spc_gridcell`) is the largest validated system. For thread/MPI curves a ≥ 8 000-atom and a ≥ 27 000-atom water box are needed. *Verify:* generate them (Phase 0) and confirm they run on both engines.
- [x] **A15 — Synthetic lattice water is not equilibrated.** *(Avoided: `sim_box` tiles the equilibrated 1000-water reference box instead of a lattice; 500 steps run stably on both engines.)* `generate_water_box.py` places SPC on a cubic lattice. Timing is still valid (pairlist/nonbonded cost depends on density and cutoff, not on equilibration) **provided the run does not blow up**. *Verify:* total energy stays finite and the constraint algorithm converges for the full NSTLIM; if not, pre-equilibrate 10 ps with one engine and use that `.cnf` for both.
- [x] **A16 — `grid_cell` has a box-size constraint.** *(6.2 nm box, 1.4 nm cutoff: satisfied; gromosXX accepts the input.)* gromosXX requires `2·RCUTL ≤ box − cell_size` (noted in `water_1000_spc_gridcell/input.toml`). Larger boxes satisfy it trivially; the 648-atom box may not. *Verify:* gromosXX accepts the `.imd` without a PAIRLIST error.

---

## 2. Phase 0 — Prerequisites and fixes

### 0.1 Repair the stale harness (all three are one-line path fixes)

- [x] `crates/gromos-md/benches/md_bench.rs:16` — `tests/GROMOS_references` → `tests/gromosXX_references`; run `cargo bench -p gromos-md -- --test` to confirm it executes.
- [x] `scripts/benchmark.sh` — remove `cd gromos-rs` / `cd ..` (script lives inside the repo); scope to `cargo bench -p gromos-md` rather than `--workspace`.
- [x] `scripts/regen_gromosXX_references.py:23` — `build/program/md` → `BUILD/program/md`.
- [x] Add a `bench` target to `Justfile` that runs the Phase 2 protocol below (`just bench`, extra args pass through).

### 0.2 Build the Rust side

- [x] `cargo build --release --bin md` → `target/release/md`. Record `git rev-parse HEAD` and `rustc --version`.
- [ ] Control build without native codegen (for A1): `RUSTFLAGS="" cargo build --release --bin md --target-dir target-generic`.

### 0.3 Build gromosXX three ways (autotools; `configure` is present)

- [x] **Native serial** (the fair single-core baseline):
  ```sh
  cd .local/gromosXX/md++ && mkdir -p BUILD_native && cd BUILD_native
  ../configure CXXFLAGS="-O3 -march=native" && make -j8
  grep -n "march" Makefile   # A1: confirm the flag survived configure
  ```
- [ ] **OpenMP** (for Phase 3): `mkdir BUILD_omp && cd BUILD_omp && ../configure --enable-openmp CXXFLAGS="-O3 -march=native" && make -j8`; confirm `ldd program/md | grep gomp`.
- [ ] **MPI** (for Phase 4): `mkdir BUILD_mpi && cd BUILD_mpi && ../configure --enable-mpi CXX=mpicxx CXXFLAGS="-O3 -march=native" && make -j8`; confirm `nm -D program/md_mpi | grep MPI_Init`.
- [ ] Record for each: `program/md @version` output and the exact configure line.

### 0.4 Machine setup (repeat before every session)

- [ ] `sudo cpupower frequency-set -g performance` (or `echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor`); verify with `cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor`.
- [ ] Decide SMT policy: leave on and pin to `0-7` (default here), **or** `echo off | sudo tee /sys/devices/system/cpu/smt/control`. Record which.
- [ ] `uptime` load < 0.5; close IDE indexers.
- [ ] Optional: `sudo apt install hyperfine linux-perf` for repeat management and a hot-spot breakdown.

### 0.5 Benchmark inputs

- [x] Create `bench/` (git-ignored outputs, committed inputs) with one directory per system. For each, copy the reference `.imd` and patch: `NSTLIM` → 2000 (small) / 500 (large), all `WRITETRAJ` frequencies → 0, `NTPR`-style printing → 0 or NSTLIM. The `patch_nstlim` helper in `md_bench.rs` already does part of this — reuse it.
- [x] Generate larger water boxes — with the project's `sim_box` (gromos++ solvation tool): `scripts/make_solvated_box.py --replicate n` builds a CH4 solute file whose GENBOX is n× the *equilibrated* `water_1000_spc` box and runs `sim_box @boxsize @thresh 0.33` (n=2 → 7 975 waters / 23 926 atoms; n=3 → ~81 000 atoms), reusing `ch4_spc.top`. Run with `bench_engines.py --ref-dir bench/systems`. Requires `cargo build --release -p gromos-tools --bin sim_box`.
- [ ] Write the final matrix of systems:

| System | Atoms | Pairlist | Purpose |
|--------|------:|----------|---------|
| `aladip_solvated` | 72 | standard | Tiny — measures per-step overhead, exposes A11 (no parallel path) |
| `water_216_box` | 648 | standard | Reference-suite workhorse; NVE |
| `water_1000_spc_gridcell` | 3 000 | grid_cell | Largest validated system; first meaningful parallel point |
| `ch4_water_7975` (sim_box) | 23 926 | grid_cell, NSNB=5, twin-range | Production regime; thread scaling |
| `ch4_water_27k` (sim_box, `--replicate 3`) | ~81 000 | grid_cell | Thread + MPI scaling; not yet run |

---

## 3. Phase 1 — Parity gate: are we timing the same physics?

Speed numbers are only comparable between runs that produce the same trajectory. This phase costs minutes and prevents the classic mistake of "faster because it's doing less".

- [x] `cargo test -p gromos-md --test test_gromosXX_references` passes (37 tests, 2026-08-28 baseline).
- [ ] For each benchmark system and each benchmark `.imd` (the *patched* one, not the 100-step reference), run **both** engines once with energy output enabled (`NTWE=NSTLIM` — one frame at the end) and compare final-step energies within `energy_rel = 1e-8` and final positions within `position_abs = 1e-9` (the project tolerances). This extends parity from 100 steps to the benchmark length.
- [x] Confirm both engines echo the same `PAIRLIST ALGORITHM` and `NSNB` in their logs (A3, A5) — `SIZE` and `TYPE` were parsed but ignored by gromos-rs until 2026-08-29; both now honoured.
- [ ] For the synthetic boxes: energies finite for the full run on both engines (A15).
- [ ] **Gate:** any mismatch → stop. It is a correctness bug or a parameter ignored by one engine, and no timing below is valid until it is explained.

---

## 4. Phase 2 — One core vs one core

### Protocol

```sh
# Rust
RAYON_NUM_THREADS=1 taskset -c 2 target/release/md \
    @topo $T @conf $C @input bench.imd @fin out.cnf @verb 1 2>&1 | tee rust.log
# → "Simulation wall time: X s" (info level; excludes init, which is logged separately)

# gromosXX (native serial build)
OMP_NUM_THREADS=1 taskset -c 2 .local/gromosXX/md++/BUILD_native/program/md \
    @topo $T @conf $C @input bench.imd @fin out.cnf > gromosXX.omd
# → TIMING / DETAILED_TIMING block at the end of the .omd
```

- [x] Pin both to the *same* physical core (core 2 avoids core 0's interrupt load).
- [x] 1 warm-up run, then **5 timed repeats** per (engine, system), interleaved across engines; report **median** as the headline and **min** as the "ideal" figure. Spread > 5 % → machine not stable (A7/A9), redo.
- [x] Derive `µs/step = wall / NSTLIM` and `ns/day = NSTLIM·dt[ps]·86400 / (1000·wall[s])`.
- [ ] Also time the Rust *generic* build (no `target-cpu=native`) on `water_1000` once — this isolates how much of any Rust advantage is codegen (A1).
- [x] Pull the per-phase breakdown on `water_1000`: gromosXX `DETAILED_TIMING`; Rust `@verb 2` Timer lines (pairlist / nonbonded / bonded / constraints / integration). This is the table that tells *where* the difference is.

### History of the pass (medians unless noted; gromosXX `-march=native` throughout)

| Date | State | water_216 | water_1000 | ch4_water_fep |
|------|-------|----------:|-----------:|--------------:|
| 08-28 | first measurement | 2.49 s (0.48×) | 9.24 s (0.77×) | — |
| 08-28 | + long-range CG kernel, timing opt-in | 2.34 s | 7.63 s (0.97×) | — |
| 08-29 | + cell list from `SIZE`, distance pruning | — | 7.14 s (0.94×) | — |
| 08-29 | + `rint` minimum image, per-thread chunks | 1.60 s (0.74×) | 7.14 s | — |
| 08-29 | + SIMD batch kernel, solvent batching | 1.30 s (0.92×) | 6.85 s (1.02×) | 61.3 s (0.82×) |
| 08-29 | + `u32` pairlist, single-sided exclusions | **1.18 s (0.98×)** | **5.84 s (1.12×)** | **45.0 s (0.94×)** |

Ratios are gromosXX / Rust (> 1 means gromos-rs is faster). Earlier rows were taken under
varying machine load and block-ordered repeats; only the last row is the quiet, interleaved
protocol, which is why it — and only it — is quoted in "Where it stands".

A note on an early misreading: the first `water_1000` breakdown showed the LJ/CRF kernel at
parity and the whole deficit in the pairlist, and this document briefly concluded that
The early SIMD priority was misplaced (board card V3). That comparison was flattered by the charge-group-
grouped kernel (one `nearest_image` per CG pair) that `water_1000`'s large box enables and
`water_216`'s does not; on `water_216` the same kernel was 2× slower than gromosXX's. Both
readings were "right" for their system, and both the pairlist and the kernel needed work.

### Where it stands (2026-08-29, quiet machine, interleaved repeats)

`taskset -c 2`, `RAYON_NUM_THREADS=1`, both engines `-march=native`, load < 2, engines
interleaved per repeat. **n = timed repeats; ± is the sample standard deviation (n−1).** Raw
values are in `bench/results/*.json`; the harness prints n, mean, std, median and min.

| System | n | Rust mean ± std | gromosXX mean ± std | gromosXX / Rust |
|--------|--:|----------------:|--------------------:|----------------:|
| water_216_box (648 atoms, standard pairlist, 2000 steps) | 5 | 1.188 ± 0.019 s | 1.153 ± 0.002 s | **0.97×** |
| water_1000_spc_gridcell (3 000 atoms, grid_cell, 2000 steps) | 5 | 5.852 ± 0.033 s | 6.524 ± 0.021 s | **1.11×** |
| ch4_water_fep (2 998 atoms, 999 SPC, NSNB=1, twin-range 0.8/1.4, 2000 steps) | 3 | 47.37 ± 4.40 s | 42.46 ± 0.11 s | **0.90×** (median 0.94×) |
| ch4_water_7975 (sim_box; 23 926 atoms, production: NSNB=5, twin-range 0.8/1.4, grid_cell, 500 steps) | 3 | 24.87 ± 0.09 s | 28.29 ± 0.61 s | **1.14×** |

gromos-rs is 14 % faster than gromosXX on the production-regime box (the earlier script-built
box gave 1.19×; the `sim_box` box is the one to quote) and 11 % faster on the 3 000-atom cell-list
system; level (3 % behind) on the smallest box; and 10 % behind on the
`NSNB=1` twin-range stress case, where the Rust standard deviation (± 4.4 s: one of three repeats
ran ~8 s slow) is itself the open question — the other three systems have std ≤ 1.6 %.

**With threads (Phase 3, after fixes 14–16):** at 8 cores gromos-rs is 1.38× faster than
gromosXX-OpenMP on the production box and 1.06× on `water_1000`; before those fixes it was
0.53× and 0.43×. Full tables under Phase 3.

The production box is built by the project's own **`sim_box`** (the gromos++ solvation program,
ported in gromos-tools): `scripts/make_solvated_box.py` only prepares its inputs (a CH4 at the
origin with a target box of exactly 2× the equilibrated 1000-water box, so the tiling is
seamless) and writes the MD input. An earlier version of the script replicated the box itself;
that was replaced so the benchmark exercises the real tool. Both engines agree on E_pot at step 0
to the printed precision on the `sim_box` system.

### Optimisations applied so far

1. **Long-range twin-range interactions: 1.59 s → 0.60 s (2.6× on that phase).**
   The long-range path called the per-atom-pair kernel `lj_crf_innerloop` while the short-range
   path used `lj_crf_innerloop_cg_grouped`, which computes one `nearest_image` per charge-group
   pair instead of one per atom pair — for SPC water that is one per nine. Long-range now builds
   the same CG-group metadata (`build_cg_groups`, extracted from the inlined short-range version)
   and uses the grouped kernel under the identical safety condition, which is keyed on the
   *long-range* cutoff and so already covered these pairs. It also reuses one `ForceStorage`
   instead of allocating a fresh `n_atoms` buffer on every pairlist update.
   Verified: all 37 reference tests still pass.

2. **`neighbor_cells` allocated a `HashSet` + `Vec` per cell per pairlist update**, purely to
   deduplicate wrapped cell indices. Replaced with an allocation-free
   `neighbor_cells_into` writing into a `[usize; 27]`, using the fact that the neighbour set is
   the Cartesian product of the per-axis distinct indices. Correct for short axes, where `±1`
   aliases. **Measured effect: none** (~10 800 hash sets was not the bottleneck) — kept because it
   removes allocation from a hot path and is strictly simpler, but recorded honestly as not a win.

3. **Bug fixed, not an optimisation:** `md` panicked with "attempt to calculate the remainder with
   a divisor of zero" whenever `NTWX`/`NTWE`/`NTPR` was 0. GROMOS treats 0 as "never write", and a
   benchmark input needs exactly that (A4). Five `step % freq` sites now go through a `due()`
   helper. This bug made it impossible to benchmark with output disabled.

4. **Removed unconditional timing overhead:** `AlgorithmSequence::run_step` called
   `Instant::now()` twice per algorithm per step even with trace logging off, and `Forcefield`
   did the same for three phases. Timing is now opt-in (`enable_timing`), off in production runs.

5. **Cell list rewritten to actually prune: pairlist 3.15 s → 2.15 s.** Two defects, both
   divergences from gromosXX rather than Rust problems:
   - *Cell size was hardcoded.* gromos-rs sized cells at one full cutoff, giving a 3×3×3 grid
     on this box — where the ±1 neighbourhood covers **every** cell, so the cell list pruned
     nothing and ran O(N_cg²) with grid overhead on top. gromosXX takes the spacing from the IMD
     `PAIRLIST` block's `SIZE` field (`grid_cell_pairlist.cc:136-141`), with `auto` resolving to
     half the *short* cutoff (`in_parameter.cc:1426`). **gromos-rs parsed `SIZE` into
     `ImdParameters::size` and then never used it** — a dead input. It is now plumbed through
     `PairlistContainer::grid_size` and honoured, so `SIZE auto` gives 0.4 nm cells (7×7×7) here.
   - *Neighbour selection was a cubic shell.* Scanning ±r cells tests the enclosing cube of the
     cutoff sphere. gromosXX prunes on the true cell-to-cell distance (`minIm`,
     `grid_cell_pairlist.cc:541-546`). `cell_pair_offsets` now precomputes, once per update, the
     offsets whose minimum inter-cell distance is within the cutoff.

   Verified by the `test_cell_list_matches_standard_*` set-equality tests against the O(N²)
   oracle, plus the gromosXX bit-comparison on `water_1000_spc_gridcell`.

   > An earlier version of this fix used an invented `CELLS_PER_CUTOFF` constant tuned by
   > measurement (2.04 s). It was removed: it performed the same as the faithful port within
   > noise, and inventing a tuning knob GROMOS already exposes as a user parameter is exactly the
   > "diverge by accident" failure `porting.md` warns about.

6. **`nearest_image` used the wrong rounding mode: nonbonded kernel 1.65 s → 1.15 s (30 %).**
   `Rectangular::nearest_image` reduced the minimum image with `f64::round()`, which is
   round-half-*away-from-zero*. x86 has no native instruction for that mode, so LLVM expands it
   into several instructions per component, three components per pair. gromosXX uses `rint`
   (`boundary_implementation.cc:267`) — round-half-to-*even*, which is the one mode hardware
   implements natively as a single `roundsd`. Switching to `round_ties_even()` is therefore both
   faster **and** the faithful port; the two differ only when a separation is exactly half a box
   length. This was worth ~47 % of the entire nonbonded kernel.

7. **Parallel kernels allocated a force buffer per 1024 pairs: ~7 %.** The flat parallel path used
   `par_chunks(1024).fold(|| ForceStorage::new(n_atoms), …)`, so a 62 k-pair system built and
   merged 61 full per-atom force buffers *per step* regardless of thread count. Chunks are now
   sized `pairs / num_threads`, one accumulator per thread — which is what the CG-grouped
   parallel path already did. Also fixed a latent `par_chunks(0)` panic there: its chunk size was
   `groups.len() / num_threads`, which truncates to zero when there are fewer groups than threads.

8. **Four pairs per SIMD register: kernel 1.15 s → 0.90 s on `water_216`, 4.23 s → 3.6 s on
   `water_1000`.** The divide and square root — the expensive part of the pair arithmetic — now
   run four-wide through `wide::f64x4` (`lj_crf_interaction_x4`). Every expression keeps the
   scalar kernel's association order and `vdivpd`/`vsqrtpd` are correctly rounded, so each lane
   is bit-identical to the scalar result; the batched and remainder paths scatter in pair order,
   so no reference output changed. Along the way the three copies of the pair physics (flat,
   charge-group-grouped, solvent) collapsed into one `process_pair_slice`, parameterised by how
   the separation and the pair parameters are obtained — the first time the file's "single
   source of truth" promise has actually held.

9. **Solvent kernel batched through the same routine.** On `ch4_water_fep` short-range
   solvent–solvent is 4.0 s against gromosXX's 5.8 s.

10. **Pairlist stored as `u32`: the per-step copy is gone.** `PairlistContainer` held
    `(usize, usize)` and the forcefield narrowed the whole list to `(u32, u32)` on every update
    — 2.5 s on `ch4_water_fep`, 0.3 s on `water_1000`, with no gromosXX counterpart — and the
    16-byte pushes doubled the pairlist's write traffic. Storing `u32` natively deleted the copy
    and took `water_1000`'s pairlist from 2.15 s to 1.97 s. Mechanical change across six files;
    every consumer already wanted `u32`.

11. **Exclusion test: one sorted lookup instead of two, with gromosXX's early-out —
    `water_216` pairlist 0.28 s → 0.19 s (gromosXX 0.16 s).** `is_excluded_or_14` is the
    pairlist's innermost test (~36 M calls per 2000 steps on `water_216`). It searched
    `exclusions[i]` *and* `exclusions[j]` and probed `one_four_pairs` from both sides, although
    every construction path builds those lists symmetrically, so one side answers for both. The
    single-sided lookup now also stops on `last < j` before searching, as gromosXX's
    `Exclusions::is_excluded` does (`topology/exclusions.h:79`) — for a partner in another
    molecule that settles it without a search. The reader now stores 1-4 pairs symmetrically
    too. Sortedness and symmetry are pinned by `crates/gromos-io/tests/exclusion_invariants.rs`,
    which also checks the new lookup against the old two-sided one on every atom pair of
    `aladip` (has 1-4 pairs) and a solvated `ch4_spc`.

    Writing that test exposed a latent bug: `Topology::solvate` extended `exclusions` to the
    solvated size but left `one_four_pairs` at the solute length. The charge-group pairlist never
    asks about solvent atoms so it never showed; the atomic-cutoff pairlist (`TYPE atomic`, wired
    up in this pass) asks about every atom and would have indexed past the end on any solvated
    system. Fixed in `solvate`; the test asserts the invariant.

12. **Long-range solvent–solvent as sentinel molecule pairs when provably exact: production box
    1.01× → 1.08×.** The short-range solvent list has always stored one first-atom pair per
    molecule pair and let the solvent kernel expand it with a single shared periodic shift; the
    long-range list stored all nine atom pairs with a per-atom minimum image, because a shared
    shift is only exact when `half_box > cutoff + 2·radius`. `sentinel_long_range_is_exact`
    now decides that per update (both pairlist algorithms, same rule, so the CellList-vs-Standard
    equality tests still hold) and records the representation on the container
    (`solvent_long_is_sentinel`); the forcefield dispatches on it. `ch4_water_fep`'s 3.1 nm box
    with a 1.4 nm cutoff fails the bound and keeps the expanded path — gromosXX's behaviour —
    so it is unchanged. Verified by `gromos-forces/tests/solvent_longrange_sentinel.rs`: on the
    real 999-water box, sentinel vs expanded per-atom images agree on energies, forces and virial.
    Deliberate deviation from gromosXX (which always expands), documented per `porting.md`.

13. **Long-range kernels skip the virial when nothing reads it: 1.08× → 1.19×.** The
    long-range block always called the virial kernel variants, although the short-range path
    already selects `_novirial` when there is no pressure coupling — nine extra multiply-adds
    per pair for a tensor nobody consumed. Long-range 9.0 s → 7.3 s on the production box.

14. **Cell-list pairlist construction in parallel.** Once charge groups are binned, cells are
    independent, so the traversal runs over contiguous cell ranges on the Rayon pool, each range
    filling private lists that are appended in cell order afterwards. Pair order is therefore
    identical to the serial traversal, force accumulation order is unchanged, and every reference
    output is bit-for-bit what it was — only the wall time changes. Several ranges per thread
    keep the load balanced. The O(N²) `standard` algorithm got the same treatment over
    charge-group ranges (its per-pair `log::debug!` calls in the hot loop went too). Measured
    in Phase 3.

15. **Long-range evaluation uses the parallel kernels.** The long-range block called the serial
    kernels unconditionally while the short-range block dispatched to the Rayon versions when
    `parallel_nonbonded` was set — so the long-range phase stayed flat with threads. It now
    mirrors the short-range selection (grouped/flat, parallel/serial, virial/not). Measured in
    Phase 3.

16. **Charge-group pair-group build in parallel.** The last serial phase at 8 threads on
    `water_1000` (0.20 s, 11 % of the step): runs of identical charge-group pairs are now found per
    contiguous chunk of the pairlist in parallel and stitched where a run straddles a chunk
    boundary, so the group boundaries are exactly the serial ones. 0.20 → 0.07 s; the 8-thread
    step 1.74 → 1.35 s.

**Tried and reverted, recorded so it is not retried blind:**

- *Vectorising the minimum image inside the batch.* A stub experiment on `ch4_water_fep`
  showed the per-pair `nearest_image` at ~27 % of the long-range phase, so the reduction was
  moved into the batch as `f64x4` ops. Measured **slower**: with the batch kept as
  `[Vec3; 4]` (AoS) the shuffles to and from per-axis registers cost more than the twelve scalar
  ops they replaced (`water_216` kernel 0.90 → 1.09 s, `ch4` long-range 31 → 37 s); a full
  structure-of-arrays pipeline (image, r², force products all vectorised) recovered only to
  29.7 s on `ch4` and still lost on `water_216`. The scalar `round_ties_even` was checked to
  compile to `vroundsd` with no libm call, so that was not the reason either. Reverted; the
  batch keeps per-pair `Vec3`s and vectorises only the arithmetic.
- *Two stub experiments produced no output* because constant parameters / zero forces made
  the dynamics blow up before the timing block printed. A stub must keep the run stable.

### Final per-phase breakdown (2026-08-29, quiet machine, one run each, seconds)

Rust from `md @verb 1`'s `TIMING` block, gromosXX from its own. Single runs taken back to
back, so totals sit a little above the harness medians (thermal state), but the two engines
are affected alike and the phase ratios are what matter.

**water_216_box** (648 atoms, standard O(N²) pairlist) — Rust 1.23 s, gromosXX 1.18 s

| Phase | Rust | gromosXX |
|-------|-----:|---------:|
| pairlist update | 0.19 | 0.15 |
| shortrange solute | 0.88 | 0.88 |
| longrange | 0.12 | 0.09 |
| bonded + RF excluded | 0.03 | 0.03 |

**water_1000_spc_gridcell** (3 000 atoms, grid_cell) — Rust 5.90 s, gromosXX 6.85 s

| Phase | Rust | gromosXX |
|-------|-----:|---------:|
| pairlist update | 1.60 | 1.84 |
| CG pair groups | 0.17 | — |
| shortrange solute | 3.53 | 4.29 |
| longrange | 0.42 | 0.41 |
| bonded + RF excluded | 0.12 | 0.12 |

**ch4_water_7996** (23 989 atoms, 7 996 SPC, grid_cell, NSNB=5, twin-range 0.8/1.4, 500 steps) —
Rust 25.6 s (was 26.7 s before fixes 12–13), gromosXX 28.1 s

| Phase | Rust before 12–13 | Rust now | gromosXX |
|-------|-----:|-----:|---------:|
| pairlist update | 7.6 | 7.5 | 7.5 |
| shortrange solvent–solvent | 7.8 | 8.6 | 10.0 |
| longrange solvent–solvent | 9.2 | **7.3** | 8.1 |
| SHAKE | 1.4 | 1.5 | 2.1 |

**ch4_water_fep** (2 998 atoms, 999 SPC, standard pairlist, NSNB=1, twin-range 0.8/1.4) —
Rust 50.6 s, gromosXX 48.8 s

| Phase | Rust | gromosXX |
|-------|-----:|---------:|
| pairlist update | 17.1 | 18.4 |
| shortrange solvent–solvent | 3.9 | 5.5 |
| longrange solvent–solvent | **27.1** | **23.4** |
| SHAKE | 0.7 | 1.0 |

What the tables say:

- **The kernels are ahead or level everywhere.** Short-range solute is 18 % faster than gromosXX
  on `water_1000`, short-range solvent 29 % faster on `ch4`, and level on `water_216`.
- **`water_216`'s remaining 4 %** is the O(N²) standard pairlist (0.19 vs 0.15) and the
  long-range list (0.12 vs 0.09) — both tiny-system, cache-resident effects.
- **`ch4`'s remaining 6 % is one phase: long-range solvent–solvent, 27.1 vs 23.4 s.** That
  path evaluates ~1.4 M expanded atom pairs per step through the flat kernel with a per-pair
  minimum image, because at this box size (3.1 nm) a 1.4 nm cutoff is too close to half the box
  for the shared-shift trick to be safe. gromosXX takes the same per-atom images and is simply
  faster at it. Note the regime: `NSNB=1` rebuilds the pairlist and the whole 1.4 nm shell every
  step; with the usual `NSNB=5` this phase is amortised five-fold.
- **In the production regime every phase is now ahead or level.** Long-range solvent–solvent
  was the deficit (9.2 vs 8.1 s) until fixes 12–13 took it to 7.3 s; short-range solvent (−14 %),
  and SHAKE (−29 %) were already ahead, the pairlist is level. The stored pairlist shrank from
  11.5 M to 1.5 M entries.

### Why the pairlist was slow (superseded by fix 5, kept for the reasoning)

- With box 3.1057 nm and cutoff 0.9 nm, `grid_dim_for_axis` yields a **3×3×3 grid**, and with
  `±1` neighbours every cell is adjacent to every other. **The cell list prunes nothing on this
  system** — it degenerates to O(N_cg²) with grid bookkeeping on top. Getting real pruning needs
  cells of ~cutoff/2 with a ±2 scan (6×6×6 → ~58 % of pairs tested).
- An experiment that skipped exclusion checks made the pairlist phase *slower* (3.14 → 3.96 s)
  because far more pairs were pushed — which identifies the dominant cost as **memory traffic
  pushing pairs**, not the exclusion lookups.
- `PairlistContainer` stores `(usize, usize)` = 16 bytes/pair, then `Forcefield` copies the whole
  list into `(u32, u32)` = 8 bytes/pair (the 0.34 s "pairlist cache" line). Storing `u32` natively
  would halve the write traffic *and* delete the copy. This is a cross-crate change and is the
  single highest-value next step.
- `is_excluded_or_14` costs up to four pointer-chasing lookups per atom pair
  (`exclusions[i]`, `exclusions[j]`, `one_four_pairs[i]`, `one_four_pairs[j]`). For pure water
  `one_four_pairs` is `n_atoms` *empty* vectors, so `!is_empty()` passes and the last two lookups
  are pure waste. A cached "any 1-4 pairs at all" flag would skip them.

> **Caveat on the re-measurement:** the head-to-head totals above are the clean 5-repeat run taken
> before the machine became busy. A later confirmation run was invalidated by an unrelated process
> at 648 % CPU (A9) — both engines inflated ~45 %, spread 28–65 %. The per-phase breakdown is
> unaffected in its *proportions*, but **the post-optimisation wall-clock totals still need one
> clean re-run on an idle machine.**

---

## 5. Phase 3 — Thread scaling

### Protocol

- [x] Threads `k ∈ {1, 2, 4, 8}` pinned to `taskset -c 0-$((k-1))` (harness does this); `k = 16` **SMT** point measured (A10): no gain for Rust, a loss for gromosXX.
- [x] Rust: `RAYON_NUM_THREADS=$k`. gromosXX: `OMP_NUM_THREADS=$k` with the `BUILD_omp` binary — its `k=1` time is 2.5 % above `BUILD_native` (OpenMP overhead), recorded.
- [x] Systems: `water_1000` and `ch4_water_7975` (24k); the 81k box is not yet run.
- [x] 3 repeats per point; mean ± std reported, medians in the JSON.
- [x] Speedup and efficiency are in the tables (`scripts/scaling_table.py`). Amdahl serial fraction from S(8), `f = (8/S − 1)/7`: Rust 0.12 on `water_1000` (S = 4.4) and 0.13 on the 24k box (S = 4.1); gromosXX-omp 0.20 and 0.19. Rust's serial share matches the phases that are still serial (SHAKE, integration, thermostat, COM removal ≈ 10–13 % of a 1-thread step on the 24k box), so nothing hidden is limiting; the next thread-scaling gain is SHAKE/SETTLE over molecules.

### Results (2026-08-29, n = 3 interleaved repeats, mean ± std, pinned to cores 0..k−1)

gromosXX from the `--enable-openmp` build (`BUILD_omp`); its 1-thread time is 2.5 % above the
serial build's, the OpenMP overhead. Rust: `RAYON_NUM_THREADS=k`.

**water_1000_spc_gridcell** (3 000 atoms, 2000 steps)

| k | Rust s | S(k) | E(k) | gromosXX-omp s | S(k) | E(k) | XX / Rust |
|--:|-------:|-----:|-----:|---------------:|-----:|-----:|----------:|
| 1 | 5.94 ± 0.04 | 1.00 | 1.00 | 6.69 ± 0.01 | 1.00 | 1.00 | 1.13× |
| 2 | 4.45 ± 0.01 | 1.33 | 0.67 | 3.71 ± 0.01 | 1.80 | 0.90 | 0.83× |
| 4 | 3.62 ± 0.03 | 1.64 | 0.41 | 2.13 ± 0.01 | 3.14 | 0.78 | 0.59× |
| 8 | 3.29 ± 0.03 | 1.80 | 0.23 | 1.41 ± 0.02 | 4.73 | 0.59 | 0.43× |

**The single-core lead does not survive threads.** gromosXX's OpenMP build scales to 4.7× on
8 cores; gromos-rs reaches 1.8× and is 2.3× *slower* than gromosXX at 8 threads. Assumption A11
predicted an Amdahl cap from the serial pairlist (~27 % of a step → ≤ 3.7×), but 1.8× is well
below even that, so the parallel part itself is scaling badly too. The per-phase breakdown at
1 vs 8 threads (below) is what says where.

**Where the threads go — per-phase, water_1000, 2000 steps, 1 vs 8 threads (seconds)**

| Phase | Rust 1 | Rust 8 | gromosXX-omp 1 | gromosXX-omp 8 |
|-------|-------:|-------:|---------------:|---------------:|
| pairlist update | 1.73 | **1.89** | 1.97 | **0.33** |
| CG pair groups | 0.19 | 0.19 | — | — |
| shortrange solute | 3.67 | 0.84 | 4.54 | 0.74 |
| longrange | 0.33 | **0.35** | 0.45 | **0.085** |
| total step | 6.12 | 3.48 | 7.32 | 1.65 |

Two phases do not move at all with threads in gromos-rs: **the pairlist** (55 % of the 8-thread
step) and **the long-range evaluation**. gromosXX parallelises both with OpenMP. The short-range
kernel scales 4.4× (gromosXX 6.1×). Fixes 14–15 below address the two serial phases.

**water_1000_spc_gridcell after fixes 14–15** (2026-08-29, n = 3, mean ± std; both engines'
1-thread times are ~5 % above the earlier run — same warm machine for both)

| k | rust s | S(k) | E(k) | gromosXX-omp s | S(k) | E(k) | gromosXX-omp / rust (time; >1 = rust faster) |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 6.30 ± 0.04 (n=3) | 1.00 | 1.00 | 7.20 ± 0.04 (n=3) | 1.00 | 1.00 | 1.14× |
| 2 | 3.91 ± 0.09 (n=3) | 1.61 | 0.81 | 3.96 ± 0.08 (n=3) | 1.82 | 0.91 | 1.01× |
| 4 | 2.25 ± 0.07 (n=3) | 2.80 | 0.70 | 2.49 ± 0.15 (n=3) | 2.90 | 0.72 | 1.11× |
| 8 | 1.60 ± 0.07 (n=3) | 3.93 | 0.49 | 1.78 ± 0.06 (n=3) | 4.06 | 0.51 | 1.11× |

Rust now scales 3.9× on 8 cores (was 1.8×) and is ahead of gromosXX-OpenMP at every thread
count. Per-phase at 8 threads: pairlist 1.89 → **0.47 s**, long-range 0.35 → **0.09 s**,
short-range 0.81 s; the remaining serial piece is the charge-group pair-group build (0.20 s,
11 % of the 8-thread step), which is the next candidate.

**water_1000_spc_gridcell after fix 16** (2026-08-29, n = 3, mean ± std; a warmer machine than
the run above — the k = 4 Rust point has one slow repeat, hence its std)

| k | rust s | S(k) | E(k) | gromosXX-omp s | S(k) | E(k) | gromosXX-omp / rust (time; >1 = rust faster) |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 6.39 ± 0.15 (n=3) | 1.00 | 1.00 | 7.16 ± 0.33 (n=3) | 1.00 | 1.00 | 1.12× |
| 2 | 3.73 ± 0.06 (n=3) | 1.71 | 0.86 | 3.81 ± 0.02 (n=3) | 1.88 | 0.94 | 1.02× |
| 4 | 2.62 ± 0.22 (n=3) | 2.44 | 0.61 | 2.41 ± 0.14 (n=3) | 2.96 | 0.74 | 0.92× |
| 8 | 1.44 ± 0.07 (n=3) | 4.44 | 0.55 | 1.53 ± 0.09 (n=3) | 4.68 | 0.58 | 1.06× |

**ch4_water_7975** (23 926 atoms, 500 steps) — *before* fixes 14–15 (k=1 from the 1-thread
run above: Rust 24.87 ± 0.09, gromosXX-native 28.29 ± 0.61; the OpenMP build at k=1 is ~2.5 %
slower than native)

| k | Rust s | S(k) | gromosXX-omp s | S(k) | XX / Rust |
|--:|-------:|-----:|---------------:|-----:|----------:|
| 1 | 24.87 ± 0.09 | 1.00 | ≈ 29.0 (native 28.29) | 1.00 | 1.14× |
| 2 | 18.9 ± 0.5 | 1.32 | 16.3 ± 0.3 | 1.78 | 0.86× |
| 4 | 19.40 ± 0.35 | 1.28 | 11.67 ± 0.41 | 2.49 | 0.60× |
| 8 | 19.34 ± 0.41 | 1.29 | 10.29 ± 0.93 | 2.82 | 0.53× |

Rust is flat from 2 threads on: with only the short-range kernel parallel, the serial pairlist
(and the serial long-range block) is the whole step on a big box. Re-measured after fixes 14–15
below.

**ch4_water_7975 after fixes 14–15** (2026-08-29, n = 3, mean ± std; k=1 from a separate run
on the same binary)

| k | rust s | S(k) | E(k) | gromosXX-omp s | S(k) | E(k) | gromosXX-omp / rust (time; >1 = rust faster) |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 25.46 ± 0.68 (n=3) | 1.00 | 1.00 | 28.78 ± 1.07 (n=3) | 1.00 | 1.00 | 1.13× |
| 2 | 14.20 ± 0.07 (n=3) | 1.79 | 0.90 | 17.69 ± 0.22 (n=3) | 1.63 | 0.81 | 1.25× |
| 4 | 8.64 ± 0.04 (n=3) | 2.95 | 0.74 | 11.15 ± 0.12 (n=3) | 2.58 | 0.65 | 1.29× |
| 8 | 6.15 ± 0.04 (n=3) | 4.14 | 0.52 | 8.52 ± 0.29 (n=3) | 3.38 | 0.42 | 1.38× |

**SMT point (A10), 16 threads on the 8 physical cores + siblings, n = 3:** Rust
7.31 ± 0.55 s (vs 6.15 s at 8), gromosXX-omp 13.51 ± 1.26 s (vs 8.52 s at 8); on `water_1000`
Rust 1.46 ± 0.06 s (vs 1.44 s at 8), gromosXX-omp 2.62 ± 0.72 s (vs 1.53 s at 8). SMT
buys gromos-rs nothing and costs gromosXX's OpenMP build 60–70 % with large variance — the
honest ceiling is 8 threads, as A10 said; leave `RAYON_NUM_THREADS` at the physical core count.

- [x] **Gate / decision (Phase 3 → 4):** `E(8)` is 0.52 on the 24k box and 0.55 on `water_1000`
  — threads are well short of saturation, so MPI on this single node (Phase 4) could only prove
  *that MPI works*, not that it is faster. Phase 4 is measured in §6.3 (2026-08-30) and is not
  scheduled; the next thread-scaling gain is SHAKE/SETTLE over molecules. **Phase 3 closed
  2026-08-29.**

gromos-rs is ahead at every thread count on the production box and the lead grows with
threads: 1.13× → 1.25× → 1.29× → 1.38×. 8-core speedup 4.1× (E = 0.52) for Rust against
3.4× (E = 0.42) for gromosXX; the serial remainder in both is SHAKE, integration and the
thermostat, plus per-step memory traffic on 24 000 atoms.

### Kernel determinism (3.9 step 0, 2026-08-29)

Measured with `scripts/kernel_determinism.py` (release `md`, `RAYON_NUM_THREADS` = 1 and 8,
**n = 3 runs per thread count**, every ENERTRJ column over every step compared exactly). Question
(3.9 A2/A3): the parallel nonbonded kernels reduce through rayon `fold/reduce`
(`innerloops.rs`) — is the result reproducible run-to-run, and is it thread-count independent?

| system | comparison | n pairs | bit-identical | worst column | max abs Δ (kJ/mol) | max rel Δ |
|---|---|---|---|---|---|---|
| water_216_box (100 steps) | 1 thr vs 1 thr, run-to-run | 3 | yes | — | 0 | 0 |
| water_216_box (100 steps) | 8 thr vs 8 thr, run-to-run | 3 | yes | — | 0 | 0 |
| water_216_box (100 steps) | 1 thr vs 8 thr | 9 | no | angle | 1.0e-10 | 1.6e-13 |
| ch4_water_fep (50 steps) | 1 thr vs 1 thr, run-to-run | 3 | yes | — | 0 | 0 |
| ch4_water_fep (50 steps) | 8 thr vs 8 thr, run-to-run | 3 | yes | — | 0 | 0 |
| ch4_water_fep (50 steps) | 1 thr vs 8 thr | 9 | no | kinetic | 1.0e-9 | 2.7e-13 |

Reading: at a fixed thread count the reduction tree is stable in practice (`par_chunks` with
`len / threads` chunks merged by `reduce`), so results are reproducible; across thread counts the
chunking — hence the summation order — changes, and the difference sits in the last printed digit
(~1e-13 relative, five orders below the 1e-8 reference tolerance). Consequence for the front-end
parity tests (3.9 G2): pinning the thread count is sufficient for exact (`==`) same-build
comparisons; `execution.parallel = Serial` remains the documented, stronger setting. Since 3.9 step 1 the Python path uses the binary's own policy (`execution.parallel = auto`: parallel
above 100 atoms); `same_path_twice` in `py-gromos/tests/test_front_end_parity.py` is bit-identical
on all 40 reference systems, and the front-ends (`InputParameters` shim, `Recipe`, `Recipe` + `Plan`)
agree exactly under that same policy.


## 6. Phase 4 — MPI: first make it exist, then make it correct, then time it

### 6.1 gromosXX side (straightforward — a rebuild)

- [ ] `BUILD_mpi/program/md_mpi` exists and `nm -D … | grep MPI_Init` is non-empty (Phase 0.3).
- [ ] `mpirun -np 1 md_mpi …` on `water_1000` reproduces the serial `md` energies to `1e-8` (the master/worker split with one worker must be a no-op).
- [ ] `mpirun -np 2` and `-np 4` reproduce the same energies (forces are summed in a different order → expect agreement to ~1e-10 relative, not bit-exact; record the actual deviation).
- [ ] Time `-np ∈ {1, 2, 4, 8}` on `water_8000` and `water_27000`, with `--bind-to core` and `OMP_NUM_THREADS=1`.

### 6.2 gromos-rs side — done 2026-08-30 (0.0.38): `md` carries MPI

What existed was not salvageable: `md_mpi.rs` spoke an API that no longer exists, the
`gromos-integrators` MPI module reinterpreted `&[Vec3]` as `&[f32]`, and nothing ever distributed
work. All of it was removed (`md_mpi`, `md_mpi_cuda`, `repex_mpi`, `mpi_scaling`, `integrators::
{mpi, remd_mpi}`) and MPI was put where it belongs — on the algorithm sequence every front-end
already uses:

- `Forcefield::pair_partition = Some((rank, size))` keeps, after every pairlist update (and for
  the initial list), only the pairs whose first atom index ≡ rank (mod size). Membership does not
  depend on list order, so every rank builds the full pairlist itself and no pairlist is ever
  communicated. The terms every rank computes in full (excluded-pair and self reaction field,
  1-4, the non-pairlist perturbed corrections) go to a separate storage that is not reduced.
- `Forcefield::set_nonbonded_reducer(f)`: `f` is called once per step with the pair-term storage
  (forces, E_LJ, E_CRF, virial, dH/dλ). Under `mpirun` it reduces to rank 0 and broadcasts the
  sum, so every rank continues with bit-identical numbers and runs the whole integrator itself —
  the only traffic per step is that one reduce + broadcast of 3N + 13 doubles.
- `RunOptions::pair_partition` carries the decomposition through `gromos-run`; `md` built with
  `--features use-mpi` (needs `LIBCLANG_PATH=/usr/lib/llvm-19/lib` for bindgen) initialises MPI,
  sets the partition, installs the reducer, and lets only rank 0 write files and stdout. The same
  `@topo @conf @input` command line, the same `bench.imd` (A3); no separate binary.
- `crates/gromos-run/tests/pair_partition.rs` checks the decomposition without MPI: the shares of
  2, 3 and 5 ranks sum to the unpartitioned pair terms (forces, energies, virial) to 1e-9.

### 6.3 gromos-rs correctness and timing — measured 2026-08-30

- [x] `mpirun -np 1 md` == serial `md` **bit-for-bit** on `water_1000` (`.tre`, `.trf`, `final.cnf`
      byte-identical: with one rank the reducer is the identity).
- [x] `-np 2, 4`: last-frame energies identical to the serial run at all ten printed digits;
      forces and final positions differ by at most 1e-9 (the printed resolution) — the cross-rank
      sum adds the per-atom force contributions in a different order.
- [x] Timed with `--bind-to core`, `RAYON_NUM_THREADS=1`, `scripts/bench_mpi.py` (n = 3 repeats,
      median and std). Results below.
- [ ] Hybrid point (`-np 2` × `RAYON_NUM_THREADS=4`): not measured.
- [ ] gromosXX's `md_mpi` (6.1): still the stub build; the gromosXX columns stay empty.

### Results (2026-08-30; 8-core desktop, `mpirun --bind-to core`, `RAYON_NUM_THREADS=1` per rank, `scripts/bench_mpi.py`)

| System | steps | run | np | wall s (median) | std | n | speedup vs serial | max rel. energy deviation vs serial |
|--------|--:|-----|--:|--:|--:|--:|--:|--:|
| water_1000 (3 000 atoms) | 500 | serial | 1 | 1.55 | 0.00 | 3 | 1.00 | 0.0e+00 |
| water_1000 (3 000 atoms) | 500 | np1 | 1 | 1.49 | 0.00 | 3 | 1.04 | 0.0e+00 |
| water_1000 (3 000 atoms) | 500 | np2 | 2 | 1.27 | 0.00 | 3 | 1.22 | 2.4e-09 |
| water_1000 (3 000 atoms) | 500 | np4 | 4 | 1.13 | 0.00 | 3 | 1.37 | 1.3e-09 |
| water_1000 (3 000 atoms) | 500 | np8 | 8 | 1.18 | 0.01 | 3 | 1.31 | 9.4e-10 |
| ch4_water_7975 (7 975 atoms) | 200 | serial | 1 | 9.60 | 0.25 | 3 | 1.00 | 0.0e+00 |
| ch4_water_7975 (7 975 atoms) | 200 | np1 | 1 | 9.56 | 0.01 | 3 | 1.00 | 0.0e+00 |
| ch4_water_7975 (7 975 atoms) | 200 | np2 | 2 | 7.09 | 0.02 | 3 | 1.35 | 0.0e+00 |
| ch4_water_7975 (7 975 atoms) | 200 | np4 | 4 | 6.00 | 0.01 | 3 | 1.60 | 0.0e+00 |
| ch4_water_7975 (7 975 atoms) | 200 | np8 | 8 | 6.44 | 0.01 | 3 | 1.49 | 0.0e+00 |
| ch4_water_7995 (23 986 atoms, `make_solvated_box --replicate 2`) | 100 | serial | 1 | 4.85 | 0.02 | 3 | 1.00 | 0.0e+00 |
| ch4_water_7995 (23 986 atoms, `make_solvated_box --replicate 2`) | 100 | np1 | 1 | 4.86 | 0.01 | 3 | 1.00 | 0.0e+00 |
| ch4_water_7995 (23 986 atoms, `make_solvated_box --replicate 2`) | 100 | np2 | 2 | 3.58 | 0.01 | 3 | 1.35 | 0.0e+00 |
| ch4_water_7995 (23 986 atoms, `make_solvated_box --replicate 2`) | 100 | np4 | 4 | 3.01 | 0.00 | 3 | 1.61 | 0.0e+00 |
| ch4_water_7995 (23 986 atoms, `make_solvated_box --replicate 2`) | 100 | np8 | 8 | 3.21 | 0.01 | 3 | 1.51 | 0.0e+00 |

Reading: with one rank the MPI build is the serial build (bit-identical output); the pair kernels
are the only distributed part, so the speedup is bounded by their share of the step — every rank
still builds the pairlist, evaluates the bonded terms, applies constraints and integrates. The
reduce + broadcast of 3N+13 doubles per step is cheap next to that redundancy. Shared-memory
threading (`RAYON_NUM_THREADS`, §4) remains the better single-node choice; the MPI seam is what a
multi-node run needs, and the redundant per-rank work is the next thing to cut (pairlist and
bonded terms partitioned the same way). gromosXX's `md_mpi` columns stay empty until §6.1 is done.

## 7. Recording a run (do this for every table row)

Each result must carry enough context to be reproduced. Keep a `bench/results/<date>/` directory with one `meta.txt` per session:

- [ ] `git rev-parse HEAD` of gromos-rs; `rustc --version`; contents of `.cargo/config.toml` (or the `RUSTFLAGS` override used)
- [ ] gromosXX: `md @version` (version + build date), the configure line, `MY_CXXFLAGS` from the Makefile
- [ ] `lscpu | grep -E 'Model name|MHz'`, governor, SMT state, kernel (`uname -r`), Open MPI version
- [ ] The exact `bench.imd` used (copy it in), the `taskset`/env line, the five raw wall times

---

## 8. What "done" looks like

- [~] Phase 1 parity: the 37-system reference suite passes after every change and both engines agree on E_pot at step 0 on the generated boxes; the explicit end-of-run energy comparison at benchmark length (Phase 1, second item) has not been done.
- [x] A one-core table (Section 4) with n / mean ± std, breakdowns on all systems, and the ratio attributed per phase.
- [x] A thread-scaling curve (Section 5) for both engines on 2 systems (3 000 and 24 000 atoms), with fitted serial fraction — **Phase 3 closed 2026-08-29**.
- [ ] Either an MPI table (Section 6) or a precise statement of which 6.2 item is still blocking, with the error text.
- [ ] `CHANGELOG.md` entry and a `docs/` link to the results directory.

### Quiet-machine re-measurement — done 2026-08-29

Taken with load 0.8, no competing process, engines interleaved per repeat; results are the
"Where it stands" table. The Rust-side 17 % spread on `ch4_water_fep` (one slow repeat of
three) is the one open measurement question.

---

## Phase X — spin-off, for a future round: the integrated GPU

*Not part of the current plan.* Phases 0–3 are the single-node CPU story and Phase 4 is MPI;
this section records a feasibility probe (2026-08-29) so the question "should the APU's GPU be
used?" has data behind it when it is picked up. Nothing here is implemented.

Question: can the APU's integrated GPU speed up the force calculation? Facts measured on this
machine (`bench/gpu_probe`, wgpu over Vulkan/RADV, nothing installed beyond the stock Mesa driver):

| Item | Found |
|------|-------|
| Device | AMD Radeon **760M** (Phoenix, RDNA3, 8 CUs, ≤ 2.7 GHz), shares the 32 GB system RAM |
| Compute stacks | Vulkan (RADV, Mesa 25.0) only. No ROCm/HIP (gfx1103 is not officially supported anyway), no OpenCL |
| 64-bit floats in shaders | **yes** — `SHADER_F64` and `SHADER_INT64` exposed to wgpu |
| Existing `gromos-forces/src/gpu` | CUDA-only (`cudarc` + `.cu`), never compiled, f32-centric — not reusable on AMD |
| Dispatch overhead | ~0.2 ms per submit+dispatch (matters at 0.7 ms/step on `water_1000`, not at 12 ms/step on the 24k box) |

**Pair-arithmetic throughput** (the LJ + CRF expression: one divide, one sqrt, ~15 FMAs;
1 M independent pairs × 64 dependent evaluations, so the CPU numbers are latency-bound and
understate its SIMD-batched real rate by ~2×):

| | M pair-evaluations / s |
|--|--:|
| CPU f64, 1 core (dependent chain) | 68 |
| CPU f64, 8 cores | 511 |
| **GPU f64** | **2 134** |
| **GPU f32** | **45 951** |

So on raw arithmetic the iGPU is ~2–4× the eight CPU cores in **f64** (RDNA3 runs f64 at ~1/16
of f32, and f64 divide/sqrt are slow, yet 8 CUs × 64 lanes still win) and ~40–90× in **f32**.

**What that does and does not mean.** The real kernel is not arithmetic-bound: it gathers
positions and parameters per pair and scatters forces to two atoms. On a GPU the scatter must be
avoided (no f64 atomics; f32 atomics are an extension), which means the standard design — one
thread per atom over a *per-atom* neighbour list, each pair computed from both sides (2× the
arithmetic, no atomics; the CSR layout that was tried and dropped on the CPU is exactly the GPU
layout). With 46 G f32 or 2 G f64 evaluations/s to spend, the 2× is affordable. Data movement is
cheap on an APU (positions in, forces out: ~0.6 MB per step for 24 000 atoms over shared memory).

The ceiling is Amdahl, again. On the 24k box at 8 threads the step is 12.3 ms, of which the
nonbonded kernels are roughly 4–5 ms; the rest is the (now parallel) pairlist, **SHAKE (~3 ms,
serial)**, integration and thermostat. Moving the nonbonded work to the GPU therefore caps at
~1.6× on that box unless SHAKE/SETTLE and the pairlist move too — and parallelising SHAKE over
molecules on the CPU (or using SETTLE for water) is a cheaper win that helps every configuration.

**Precision is the real decision.** An f64 GPU path keeps the project's f64-everywhere rule and
the reference tolerances (energy 1e-8, force 1e-6) — it would be validated by the existing
reference suite, just like every CPU change. An f32 path is 20× faster still but *cannot* meet
those tolerances; it would be a GROMACS-style mixed-precision mode (f32 pair forces, f64
accumulation and integration) validated by energy-conservation drift rather than bit-comparison,
and must live behind its own feature flag (e.g. `gpu-f32`) — never the default.

**If this phase is ever started, the recommended shape:** a `PotentialProvider` (`provider.rs` — the additive
provider pattern was designed for exactly this) behind a `gpu` feature using `wgpu`, f64 first
(reference-suite-verified), per-atom neighbour lists uploaded every `NSNB` steps, forces read back
per step; f32 as a second flag later. Order of value on this machine: parallel SHAKE/SETTLE on the
CPU first, then the GPU provider. Not started; the probe is kept in `bench/gpu_probe`.

## 7. Thermodynamic integration, engine vs engine (2026-08-29)

`scripts/ti_ch4.py --steps 5000`: CH4 → dummy in 999 SPC (54a7, the vol. 7 tutorial system, the
inputs proven bit-for-bit at λ = 0, 0.25, 0.5, 0.75, 1), eleven λ windows of 10 ps each, NVT 300 K,
run by the native gromosXX and by `md` from the same `.imd` files, integrated with the project's
`ext_ti_ana` (trapezoid over ⟨∂H/∂λ⟩, first 20 % of each window dropped). One run per window — the
per-window means differ by up to 25 kJ/mol at λ = 0.9 (independent 10 ps trajectories in the
soft-core end-point region, ee 16–20 kJ/mol), so the ΔG agreement is statistical, not bit-for-bit;
the bit-for-bit statement is the reference suite. Binary: `md` 0.0.33 (the fixes of 0.0.33/0.0.34 do
not touch a zero-charge, atom-1, equal-mass perturbation).

| λ | ⟨∂H/∂λ⟩ gromosxx ± ee (n) | ⟨∂H/∂λ⟩ gromos-rs ± ee (n) | Δ (rs − XX) |
|---|---|---|---|
| 0.00 |    13.549 ± 0.609 (400) |    13.918 ± 0.328 (401) | +0.368 |
| 0.10 |    11.941 ± 0.769 (400) |    13.336 ± 0.698 (401) | +1.395 |
| 0.20 |     9.445 ± 0.795 (400) |    10.795 ± 0.924 (401) | +1.350 |
| 0.30 |     9.233 ± 0.329 (400) |     9.384 ± 0.931 (401) | +0.151 |
| 0.40 |     6.129 ± 1.633 (400) |     6.828 ± 0.995 (401) | +0.699 |
| 0.50 |     3.704 ± 1.259 (400) |     1.895 ± 2.193 (401) | -1.808 |
| 0.60 |    -0.771 ± 1.767 (400) |     1.967 ± 2.813 (401) | +2.738 |
| 0.70 |    -2.957 ± 2.711 (400) |   -11.261 ± 3.833 (401) | -8.304 |
| 0.80 |   -13.256 ± 6.245 (400) |   -35.797 ± 1.931 (401) | -22.541 |
| 0.90 |  -107.265 ± 16.436 (400) |   -81.536 ± 20.108 (401) | +25.729 |
| 1.00 |   -65.329 ± 5.188 (400) |   -62.112 ± 9.043 (401) | +3.217 |

- **gromosxx**: ΔG = -10.97 kJ/mol; wall time 0 s for 11 windows (all cached)
- **gromos-rs**: ΔG = -10.85 kJ/mol; wall time 167 s for 11 windows
- ΔG(gromos-rs) − ΔG(gromosXX) = +0.120 kJ/mol
- Reference: ΔG_hyd(CH4) experiment ≈ +8.4 kJ/mol, i.e. ΔG(CH4→dummy) ≈ −8 kJ/mol; the tutorial's 54a7 result is of that size. Windows are independent short runs, not a converged study.

Wall time (this machine, 1 thread for gromosXX, `execution.parallel = auto` for `md`): gromosXX
≈ 40 s per window (25.8 s / 1000 steps), `md` ≈ 12 s per window (5.8 s / 1000 steps).

## 9. Free-energy features, cost (1.8)

`scripts/bench_engines.py`, one thread (`taskset -c 0`), 500 steps, engines interleaved per
repeat, n = 5, mean ± sample std (n−1), internal wall time (initialisation excluded). gromosXX
is the plain serial build (`--engines rust,plain`). Machine: governor `powersave` (the script
warns), load ≈ 1.

**1.8 step 2 (0.0.54, PRECALCLAM), 2026-10-12.** Step-1 binary = 0.0.53 (`--rust-md`).

| system | engine | n | mean ± std | gromosXX / gromos-rs |
|---|---|--:|--:|--:|
| `ch4_water_fep_l025`, outputs off | gromos-rs 0.0.53 | 5 | 11.332 ± 0.106 s | — |
| `ch4_water_fep_l025`, outputs off | gromos-rs 0.0.54 | 5 | 11.302 ± 0.081 s | — |
| `ch4_water_fep_l025`, NTWE = NTWG = 1 | gromos-rs 0.0.54 | 5 | 11.308 ± 0.062 s | |
| `ch4_water_fep_l025`, NTWE = NTWG = 1 | gromosXX | 5 | 11.929 ± 0.021 s | **1.055 ± 0.006** |
| `ch4_water_fep_precalc`, `PRECALCLAM 81 0 1`, NTWE = NTWG = 1 | gromos-rs 0.0.54 | 5 | 12.150 ± 0.054 s | |
| `ch4_water_fep_precalc`, `PRECALCLAM 81 0 1`, NTWE = NTWG = 1 | gromosXX | 5 | 13.136 ± 0.024 s | **1.081 ± 0.005** |

Without the block, 0.0.54 / 0.0.53 = 0.997 ± 0.012: no change beyond noise. The worst case
(81 λ points, energies every step) costs gromos-rs 7.4 % ± 0.7 % and gromosXX 10.1 % ± 0.3 %
over the same run without the block, so the ratio moves from 1.055 to 1.081. The precalc system
is a scratch copy of the reference with NRLAM 11 → 81; both `.tre` files carry 500
`# precalclam` frames.
