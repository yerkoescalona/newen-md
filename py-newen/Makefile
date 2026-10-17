.DEFAULT_GOAL := help

SHELL=bash

##@ General

.PHONY: help
help: ## Display this help
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Setup

.PHONY: sync
sync: ## Sync deps and compile the Rust extension (dev/debug build) into .venv
	uv sync --all-groups

##@ Building

.PHONY: build
build: sync ## Compile and install Python extension for development

.PHONY: build-release
build-release: sync ## Compile and install with optimizations
	uv run maturin develop --release

##@ Testing

.PHONY: test
test: build ## Run Python tests
	uv run pytest tests/ -v

.PHONY: test-verbose
test-verbose: build ## Run tests with verbose output
	uv run pytest tests/ -vv -s

.PHONY: test-fast
test-fast: ## Run tests without rebuilding
	uv run pytest tests/ -v

.PHONY: test-coverage
test-coverage: build ## Run tests with coverage
	uv run pytest tests/ --cov=newen --cov-report=html --cov-report=term

##@ Code Quality

.PHONY: fmt
fmt: sync ## Format Python code
	uv run ruff check --fix python/ tests/ examples/
	uv run ruff format python/ tests/ examples/

.PHONY: fmt-check
fmt-check: sync ## Check Python formatting without modifying files
	uv run ruff check python/ tests/ examples/
	uv run ruff format --check python/ tests/ examples/

.PHONY: lint
lint: sync ## Lint Python code
	uv run ruff check python/ tests/ examples/

.PHONY: type-check
type-check: sync ## Run type checking with mypy
	uv run mypy python/ --ignore-missing-imports

.PHONY: pre-commit
pre-commit: fmt lint test ## Run pre-commit checks

##@ Notebooks

.PHONY: notebooks
notebooks: build ## Run Jupyter notebooks
	uv run jupyter notebook notebooks/

.PHONY: notebooks-test
notebooks-test: build ## Execute all notebooks headlessly (same as CI)
	cd notebooks && for nb in *.ipynb; do \
		echo "Testing $$nb..."; \
		uv run jupyter nbconvert --to notebook --execute "$$nb" --output temp.ipynb || exit 1; \
		rm temp.ipynb; \
	done

##@ Documentation

.PHONY: docs
docs: sync ## Build static documentation site into site/
	uv run --with mkdocs-material mkdocs build --strict

.PHONY: docs-serve
docs-serve: sync ## Serve documentation locally with live-reload (http://localhost:8000)
	uv run --with mkdocs-material mkdocs serve

##@ Cleaning

.PHONY: clean
clean: ## Clean build artifacts (preserves .venv)
	rm -rf .pytest_cache
	rm -rf .coverage htmlcov
	rm -rf .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.so" -delete

.PHONY: clean-cache
clean-cache: ## Clean Python cache files only
	rm -rf .pytest_cache
	rm -rf .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete

##@ Advanced

.PHONY: profile
profile: build-release ## Profile Python code
	uv run python -m cProfile -o profile.stats examples/basic_usage.py
	uv run python -m pstats profile.stats

.PHONY: publish-test
publish-test: ## Publish to TestPyPI
	uv run --with twine twine upload --repository testpypi dist/*

.PHONY: publish
publish: ## Publish to PyPI
	uv run --with twine twine upload dist/*

.PHONY: all
all: fmt lint test build-release ## Run formatter, linter, tests and build release
