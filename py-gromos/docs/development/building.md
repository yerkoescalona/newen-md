# Building py-gromos

## Quick build

From the repository root:

```bash
make build-python     # debug build — fast compile, slower runtime
make build-release    # release build — slower compile, faster runtime
```

Both targets create `.venv` if it does not exist, install `maturin`, and call
`maturin develop` inside `py-gromos/`. From inside `py-gromos/` the same targets
are available via its own `Makefile`.

## What maturin does

[Maturin](https://github.com/PyO3/maturin) is the build tool for PyO3-based packages. It:

1. Runs `cargo build` on the `py-gromos` Rust crate (which depends on `pyo3-gromos`).
2. Copies the resulting `.so` into `python/gromos/gromos.abi3.so`.
3. Installs the package as editable (`pip install -e`).

The extension is built as a stable-ABI wheel (`abi3`, Python ≥ 3.12) so one
binary works across Python minor versions.

## Rust crate structure

```
py-gromos/src/lib.rs          ← thin entry point, re-exports pyo3-gromos
crates/pyo3-gromos/src/
    lib.rs                    ← module assembly, Vec3/Energy/Frame/rmsd/rdf
    system.rs                 ← System
    topology.rs               ← Topology
    py_conf.rs                ← Configuration
    recipe.rs                 ← Recipe / Term / Algorithm / Plan, registries, exceptions (3.9)
    parameters.rs             ← InputParameters + factories (deprecated shims)
    simulation.rs             ← Simulation (one entry: prepare → plan → instantiate → start)
    algorithm_sequence.rs     ← AlgorithmSequence + building blocks (deprecated; step 4 deletes it)
```

Changes to any file under `crates/pyo3-gromos/` require a rebuild.

## Build configurations

| Command | Optimisation | Compile time | Use for |
|---------|-------------|-------------|---------|
| `make build-python` | `-O0` (debug) | ~40 s | development, quick iteration |
| `make build-release` | `-O3` + LTO | ~2–3 min | benchmarks, production use |

## Incremental builds

Cargo is incremental by default — only changed crates recompile. A full cold
build from scratch takes ~40 s (debug) or ~3 min (release) on a modern machine.

## Running tests after a build

```bash
make test-python              # build-python + pytest
.venv/bin/pytest py-gromos/tests/ -v   # pytest only, no rebuild
```

## Building a wheel for distribution

```bash
cd py-gromos
../.venv/bin/maturin build --release
# wheel written to: ../target/wheels/gromos-<version>-cp312-abi3-linux_x86_64.whl
# (<version> is the workspace version in the root Cargo.toml)
pip install ../target/wheels/gromos-*.whl
```

## Troubleshooting

**`ModuleNotFoundError: No module named 'gromos'`**  
The extension has not been built yet. Run `make build-python`.

**`ImportError: ... undefined symbol`**  
Version mismatch between the `.so` and the Python interpreter. Rebuild with
`make build-python` to match the active Python.

**Linker errors on Linux**  
Install system dependencies:
```bash
sudo apt install build-essential libssl-dev pkg-config
```

**`maturin: command not found`**  
Run `pip install maturin` inside your virtual environment, or let
`make build-python` handle it.
