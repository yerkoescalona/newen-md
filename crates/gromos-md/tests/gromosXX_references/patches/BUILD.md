# Patched gromosXX / gromos++ oracles (item 1.8): how they were built

## Sources

Worktrees (new branches, no commits; the patches live only in the working trees and in the files here):

    git -C /home/yerko/Development/gromos-rs/.local/gromosXX worktree add \
        /home/yerko/Development/gromos-rs/.local/gromosXX-patched -b divergences-fix 247c95c4
    git -C /home/yerko/Development/gromos-rs/.local/gromosPlsPls worktree add \
        /home/yerko/Development/gromos-rs/.local/gromosPlsPls-patched -b divergences-fix 9cd9f6b3

The autotools outputs are untracked in the originals (`configure`, `aclocal.m4`, `config.h.in`, `config/`,
`m4/`, `autom4te.cache/`, every `Makefile.in`). They were copied from the originals with `cp -a` (the list is
`git -C <orig> status --short --ignored`, minus the BUILD* directories), so the same generated configure
script is used rather than a fresh autoreconf.

Patches applied to the working trees (each file applies with `git apply` from the worktree root):

    md++:     D1_pair_kernel_args.patch, D2_eds_weights.patch
    gromos++: D4_bar_error.patch, D5_ext_ti_ana_dlam.patch, D8_ext_ti_merge_loop.patch

## Configuration of the originals

- `.local/gromosXX/md++/BUILD/config.log`: `$ ../configure` (no options), g++, `LIBS='-lz -lm -lfftw3'`,
  `CXXFLAGS=' -I/usr/include'`, `AM_CXXFLAGS = -std=c++11 -O3 -DNDEBUG`; non-MPI (`program/md_mpi` is the stub),
  no OpenMP, no CUDA. (`BUILD_native` and `BUILD_omp` exist but are not the oracle.)
- `.local/gromosPlsPls/gromos++/BUILD/config.log`: `$ ../configure` (no options), g++,
  `CXXFLAGS = -std=c++11`, `AM_CXXFLAGS = -g -O3 -funroll-loops -fno-gcse -DNDEBUG`, `LIBS = -lfftw3 -lz -lm`,
  no OpenMP.

## Build commands (2026-10-13)

    cd /home/yerko/Development/gromos-rs/.local/gromosXX-patched/md++
    mkdir BUILD && cd BUILD && ../configure > configure.out 2>&1 && make -j8 > make.out 2>&1

    cd /home/yerko/Development/gromos-rs/.local/gromosPlsPls-patched/gromos++
    mkdir BUILD && cd BUILD && ../configure > configure.out 2>&1 && make -j8 > make.out 2>&1

Both exit 0. Check of sameness: the compiler/flag/library lines (`CXX CXXFLAGS CPPFLAGS LIBS LDFLAGS
AM_CXXFLAGS MY_CXXFLAGS DEFS`) of the top-level, `src/` and `program(s)/` Makefiles are identical to the
originals'; `config.h` differs only in `MD_DATE` / `GROMOS_DATE` (build time stamps).

Binaries:

    /home/yerko/Development/gromos-rs/.local/gromosXX-patched/md++/BUILD/program/md
    /home/yerko/Development/gromos-rs/.local/gromosPlsPls-patched/gromos++/BUILD/programs/{bar,ext_ti_ana,ext_ti_merge}

md prints `build date : Tue Oct 13 12:10:47 CEST 2026` in its stdout header (the original: Sat Aug 22 11:35:53).

## Verification scripts (this directory)

- `compare_md.py OUT SYSTEM...`: both md binaries, invocation as `scripts/regen_gromosXX_references.py`,
  file-by-file byte comparison (`byte_identity_unaffected.txt`).
- `d1_frame0_check.py OUT MD`: D1 gate against `kernel.py` intended mode (`d1_frame0_check.csv`).
- `e2_driver.py MD OUT` + `e2_summarize.py`: the unchanged E2 run.py against another binary (`e2_patched.csv`,
  full CSVs in `e2_patched_run/`). With the unpatched binary the driver reproduces the published
  `gromos-divergences/data/e2_*.csv` byte for byte (e2_vr.csv except its command column).
- `gpp_check.py WORK`: gromos++ orig vs patched on the chain, the A2 dihedral variants and the synthetic merges
  (`gpp_check.txt`, `gpp_value_check.csv`).
- `eds_regen.sh`, `chain_regen.sh`, `a2_regen.sh`: the commands behind `eds_refs_fixed/`, `ti_oracle/chain_fixed/`
  and `ti_oracle/A2_fixed/`.
Python: `uv run --no-project --with pandas --with numpy --with scipy python <script>` (system python has no pandas).
