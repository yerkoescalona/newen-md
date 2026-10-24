# newen — local dev convenience wrapper (not used by GitHub CI, which
# runs its own native steps per .github/workflows/*.yml)
# Usage: make check        (fast, before commit)
#        make ci           (full local suite: fmt/lint/test)
#        make build-python (compile + install Python extension into .venv)
#        make test-python  (run py-newen test suite)

VENV := .venv
VENV_BIN := $(VENV)/bin

.PHONY: check ci fmt fmt-check lint-rust lint test \
        .venv build-python build-release test-python \
        docs docs-serve

# ── Rust ─────────────────────────────────────────────────────────────────────

# Run before every commit (matches .claude/overview.md rule)
check:
	cargo fmt
	cargo check --tests

# Full CI equivalent — run before pushing
ci: fmt-check lint-rust test

fmt:
	cargo fmt

fmt-check:
	cargo fmt --check

lint-rust:
	cargo clippy --workspace --all-targets -- \
		-W warnings \
		-A clippy::upper_case_acronyms \
		-A non_snake_case

test:
	cargo test --workspace

# ── Python ───────────────────────────────────────────────────────────────────
# Managed via `uv` — see py-newen/pyproject.toml for pinned deps
# (runtime deps + the `dev` dependency-group). Creates/uses py-newen/.venv.

# Sync deps and compile the Rust extension (dev/debug build)
build-python:
	cd py-newen && uv sync --all-groups

# Compile with optimisations
build-release: build-python
	cd py-newen && uv run maturin develop --release

# Run the Python test suite (build first)
test-python: build-python
	cd py-newen && uv run pytest tests/ -v

# ── Docs ─────────────────────────────────────────────────────────────────────

# Create a root virtualenv for tooling that isn't part of the py-newen
# package itself (mkdocs is not a project dependency).
.venv:
	python3 -m venv $(VENV)
	$(VENV_BIN)/pip install --quiet --upgrade pip

# Build static documentation site into py-newen/site/
docs: .venv
	$(VENV_BIN)/pip install --quiet mkdocs-material
	cd py-newen && $(CURDIR)/$(VENV_BIN)/mkdocs build --strict

# Serve documentation locally with live-reload (http://localhost:8000)
docs-serve: .venv
	$(VENV_BIN)/pip install --quiet mkdocs-material
	cd py-newen && $(CURDIR)/$(VENV_BIN)/mkdocs serve
