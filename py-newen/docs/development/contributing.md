# Contributing

## Setup

```bash
git clone https://github.com/yerkoescalona/newen-md.git
cd newen-md
make build-python
source .venv/bin/activate
```

## Workflow

Before committing, run the pre-commit checks defined in the root `Makefile`:

```bash
make check       # cargo fmt + cargo check --tests (fast, run every commit)
make test-python # build + full Python test suite
```

Before pushing:

```bash
make ci          # fmt-check + clippy + cargo test --workspace
```

## Where things live

| What you want to change | Where |
|------------------------|-------|
| Python API exposed to users | `crates/pyo3-newen/src/` |
| Python module packaging | `py-newen/python/newen/__init__.py` |
| Type stubs | `py-newen/python/newen/newen.pyi` |
| Python tests | `py-newen/tests/` |
| Documentation | `py-newen/docs/` |
| Rust MD engine | `crates/newen-{core,forces,integrators,io}/` |
| Reference test data | `crates/newen-md/tests/gromosXX_references/` |

## Adding a new Python-exposed feature

1. Implement in Rust under `crates/pyo3-newen/src/`.
2. Register the class/function in `crates/pyo3-newen/src/lib.rs`.
3. Re-export from `py-newen/python/newen/__init__.py`.
4. Add the type signature to `py-newen/python/newen/newen.pyi`.
5. Write a test in `py-newen/tests/`.
6. Document in `py-newen/docs/api/reference.md`.

## Adding a reference test system

Reference data lives in `crates/newen-md/tests/gromosXX_references/`. Each
system directory contains:

```
<system_name>/
    input.toml          ← paths + tolerances + reference values
    <system>.topo
    <system>.cnf
    <system>.in
    expected/
        energies.tre
        trajectory.trc
        forces.trf
        final.conf
```

See `crates/newen-md/.claude/CONTEXT.md` for the full procedure.

## Code style

- **Rust**: `cargo fmt` (enforced in CI), `cargo clippy` with workspace-level lint flags.
- **Python**: `ruff check` + `ruff format` (`make fmt` inside `py-newen/`).
- **Comments**: only when the *why* is non-obvious. No docstrings restating the signature.

## Documentation

```bash
cd py-newen
make docs        # build static site → site/
make docs-serve  # live-reload at http://localhost:8000
```

Docs are written in Markdown under `py-newen/docs/` and rendered with
[MkDocs Material](https://squidfunk.github.io/mkdocs-material/).
