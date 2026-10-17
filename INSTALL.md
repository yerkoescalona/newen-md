# Installation guide

This lists every dependency you might need for newen, from "just build and test the
Rust/Python workspace" up to "compile the real gromosXX/gromos++ reference binaries and run
the real QM/ML engines." Install only the layers you need — layers 1-2 cover `make ci` and
`make test-python`; everything after that is optional and additive.

Verified on Debian 13 (trixie) / Ubuntu; package names should match on any recent
Debian-derived distro. If you're on something else, the packages are the same libraries under
different names (see comments).

## 1. Core toolchain (required)

```bash
# Rust — rustup, not apt: the workspace pins rust-version = "1.91" in Cargo.toml,
# which is newer than most distro-packaged rustc/cargo.
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup component add rustfmt clippy

# Python tooling — py-newen is managed by uv (see py-newen/pyproject.toml), not raw pip.
sudo apt install -y python3-venv python3-dev
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Verify: `cargo build --workspace` and `make check` (repo root) should both succeed with
nothing else installed.

## 2. Python extension build (required for `make test-python`)

py-newen compiles a Rust extension via maturin, orchestrated by `uv`:

```bash
cd py-newen
uv sync --all-groups      # installs deps + compiles the extension (make build-python)
uv run pytest tests/ -v   # make test-python
```

## 3. C/C++ toolchain — only if you need to build the gromosXX/gromos++ reference sources

Needed to compile `.local/gromosXX` (md++) and `.local/gromosPlsPls` (gromos++) from source
(see `.claude/overview.md` "Reference sources"). Not needed to run `cargo test` — the
reference tests compare against fixtures already checked into
`crates/newen-md/tests/gromosXX_references/`, generated ahead of time by
`scripts/regen_gromosXX_references.py`.

```bash
sudo apt install -y build-essential cmake autoconf automake libtool pkg-config
sudo apt install -y libgsl-dev libfftw3-dev
```

Then, for each of `.local/gromosXX/md++` and `.local/gromosPlsPls/gromos++`:

```bash
./Config.sh
mkdir BUILD && cd BUILD
../configure
make -j"$(nproc)"
```

`.local/` is gitignored and not part of this repo. Populate it with:

```bash
mkdir -p .local
git clone https://github.com/biomos/gromosXX .local/gromosXX
git clone https://github.com/biomos/gromosPlsPls .local/gromosPlsPls
git clone https://github.com/biomos/gromos_tutorial_livecoms .local/gromos_tutorial_livecoms
```

Two reference sources are **not** available via a public clone and require registering at
[gromos.net](https://www.gromos.net):
- `.local/gromosXX/forcefields` (force-field parameter sets, e.g. 54a7/53a6)
- `.local/doc/gromos_book` (the GROMOS theory manual)

## 4. Real QM engines — for `--features qmmm` / the xtb and MOPAC reference tests

Both are packaged (no license step, no build-from-source) on Debian trixie:

```bash
sudo apt install -y xtb mopac
```

Without these, `crates/newen-forces/tests/xtb_qmmm_water_dimer_reference.rs` and
`mopac_reference.rs` self-skip at runtime (they check `which xtb` / `which mopac` first) —
`cargo test` still reports success, just without exercising the real engine.

## 5. Real ML potential — for `--features ml` (SchNetPack via TorchScript/`tch`)

Version-pinned exactly — `tch = "0.24.0"` links only against **libtorch 2.11.0**; any other
version fails to build rather than silently misbehaving. Full rationale in
`crates/newen-forces/src/nonbonded/schnet.rs` module docs.

`torch==2.11.0` (CPU wheel) and `schnetpack` are declared as the `ml` dependency-group in
`py-newen/pyproject.toml`, pulled from the official PyTorch CPU index via
`[tool.uv.sources]`/`[[tool.uv.index]]` — reproducible via `uv sync`, no hand-built throwaway
venv:

```bash
cd py-newen
uv sync --group ml            # resolves torch==2.11.0+cpu, schnetpack, into py-newen/.venv

cd ..
source py-newen/.venv/bin/activate
export LIBTORCH_USE_PYTORCH=1
export LD_LIBRARY_PATH="$(python3 -c 'import torch,os;print(os.path.dirname(torch.__file__))')/lib:$LD_LIBRARY_PATH"

cargo test -p newen-forces --features ml schnet::
```

Note `uv sync --group ml` only syncs the default + `ml` groups (dropping `dev`/`notebooks` from
the venv) — use `uv sync --all-groups --group ml` to keep everything installed together.

This is CPU-only unless you build a CUDA-enabled libtorch yourself — the PyPI wheel pinned
above is the `cpu` variant. There's no GPU-specific setup documented here because `use-cuda`
(the `cudarc` feature) needs an NVIDIA CUDA toolkit, which is a separate, larger lift not
covered by this stack.

## 6. Optional distributed/GPU features

Only needed if you're building with `--features use-mpi` or `--features use-cuda`:

```bash
# use-mpi (crate `mpi = "0.8"`)
sudo apt install -y libopenmpi-dev openmpi-bin
sudo apt install -y libclang-19-dev          # bindgen, for mpi-sys
export LIBCLANG_PATH=/usr/lib/llvm-19/lib     # then:
cargo build --release --features use-mpi --bin md
mpirun -np 4 target/release/md @topo sys.top @conf sys.cnf @input md.imd @fin out.cnf @tre out.tre
# every rank runs the whole driver; the nonbonded pair terms are split by rank and summed each
# step; rank 0 writes the files. np=1 is bit-identical to the serial build (BENCHMARKING.md §6).

# use-cuda (crate `cudarc`) needs the NVIDIA CUDA toolkit — see
# https://developer.nvidia.com/cuda-downloads. Requires an NVIDIA GPU; not covered further
# here since it's hardware-gated.
```

## 7. Dev/CI tooling (optional, for contributing)

```bash
pip install pre-commit maturin black ruff pytest
pre-commit install
cargo install cargo-edit cargo-audit cargo-outdated cargo-llvm-cov
```

## Verifying your install

| You installed | Verify with |
|---|---|
| Layer 1 | `cargo build --workspace` |
| Layer 2 | `make test-python` |
| Layer 3 | `.local/gromosXX/md++/BUILD/program/md` runs and prints usage |
| Layer 4 | `cargo test -p newen-forces xtb`/`mopac` shows real runs, not "skipping: ... not found" |
| Layer 5 | `cargo test -p newen-forces --features ml schnet::` passes |
| Layer 6 | `cargo build --features use-mpi` (or `use-cuda`) succeeds |
