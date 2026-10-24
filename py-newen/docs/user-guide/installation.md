# Installation

## Prerequisites

- Python 3.12 or later
- Rust stable toolchain
- A C compiler (gcc / clang / MSVC)

### Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env   # or restart your shell
rustc --version           # should print stable version
```

## From source (recommended)

The root `Makefile` handles the virtual environment and Rust compilation in one step.

```bash
git clone https://github.com/yerkoescalona/newen-md.git
cd newen-md

make build-python         # creates .venv, compiles Rust, installs newen
source .venv/bin/activate
python -c "import newen; print(newen.__version__)"
```

### What `make build-python` does

1. Creates `.venv` with Python 3.12+ if it does not exist.
2. Installs `maturin` and the `dev` extras (`pytest`, `ruff`, etc.) into `.venv`.
3. Runs `maturin develop --extras dev` inside `py-newen/`, which compiles the
   Rust extension and installs it as an editable package.

### Rebuild after Rust changes

```bash
make build-python         # re-compiles and re-installs
```

For an optimised build (slower to compile, faster to run):

```bash
make build-release
```

## Running tests

```bash
make test-python          # builds then runs py-newen/tests/
```

Or without rebuilding (after an initial `make build-python`):

```bash
.venv/bin/pytest py-newen/tests/ -v
```

## Manual setup (without the Makefile)

If you prefer to manage your own environment:

```bash
python -m venv my-env
source my-env/bin/activate
pip install maturin

cd py-newen
maturin develop --extras dev   # debug build
# or
maturin develop --release --extras dev  # optimised build
```

## Troubleshooting

**`cargo: command not found`** — Rust is not on `PATH`. Run `source $HOME/.cargo/env`.

**`maturin: command not found`** — run `pip install maturin` inside your environment.

**Linker errors on Linux** — install `gcc` and `libssl-dev`:
```bash
sudo apt install build-essential libssl-dev pkg-config
```

**Python version mismatch** — `pyproject.toml` requires Python ≥ 3.12. Check with
`python --version`.

## Platform support

| Platform | Status |
|----------|--------|
| Linux x86_64 | ✅ tested in CI |
| macOS x86_64 | ✅ |
| macOS Apple Silicon | ✅ |
| Windows x86_64 | should work, not tested |

## Next steps

- **[Quick Start](quick-start.md)** — run your first simulation
- **[API Reference](../api/reference.md)** — full class documentation
