# Contributing to gromos-rs

Thank you for your interest in contributing to gromos-rs! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Code Standards](#code-standards)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Pre-commit Hooks](#pre-commit-hooks)

## Getting Started

### Prerequisites

- Rust 1.91+ via [rustup](https://rustup.rs/) — the pinned distro `rustc`/`cargo` is too old
- Python 3.12+ and [uv](https://docs.astral.sh/uv/) (py-gromos is managed by `uv`, not raw pip)
- Git
- (Optional) C/C++ toolchain + GSL/FFTW to build the gromosXX/gromos++ reference sources
- (Optional) `xtb`/`mopac` for real QM/MM tests, libtorch 2.11.0 for the real ML potential
- (Optional) MPI runtime for parallel features
- (Optional) CUDA toolkit for GPU features

See **[INSTALL.md](INSTALL.md)** for the full dependency list, including exact package names
and version pins for every optional layer above.

### Setup Development Environment

```bash
# Clone the repository
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromos-rs

# Run the setup script
./scripts/setup-dev.sh

# Or manually:
cargo build --workspace

cd py-gromos
uv sync --all-groups
cd ..
```

## Development Workflow

### 1. Create a Branch

```bash
git checkout -b feature/my-feature
# or
git checkout -b fix/issue-123
```

### 2. Make Changes

Edit code, add tests, update documentation as needed.

### 3. Run Quality Checks

```bash
# Format code
make fmt

# Run linters
make lint

# Run tests
make test

# Or run everything at once
make pre-commit
```

### 4. Commit Changes

We follow [Conventional Commits](https://www.conventionalcommits.org/):

```bash
git commit -m "feat: add new feature"
git commit -m "fix: resolve issue #123"
git commit -m "docs: update README"
git commit -m "test: add tests for feature X"
git commit -m "refactor: improve performance"
git commit -m "chore: update dependencies"
```

**Commit Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `test`: Adding or updating tests
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `style`: Code style changes (formatting)
- `chore`: Maintenance tasks
- `ci`: CI/CD changes

### 5. Push and Create Pull Request

```bash
git push origin feature/my-feature
```

Then create a Pull Request on GitHub.

## Code Standards

**Read these first:**
- [Architecture](.claude/architecture.md) - how new and refactored code is shaped
- [Porting](.claude/porting.md) - how gromosXX is translated without inheriting its defects

### Rust Code

- **Format:** Use `rustfmt` (run `cargo fmt`)
- **Lint:** Fix all `clippy` warnings (run `cargo clippy`)
- **Style:**
  - `snake_case` for functions and regular variables
  - `PascalCase` for types and structs
  - `SCREAMING_SNAKE_CASE` for constants
  - **Scientific notation allowed**: `kT`, `pH`, `deltaG`

**Example:**
```rust
pub struct SystemState {
    temperature: f64,
    pressure: f64,
}

impl SystemState {
    pub fn new(temperature: f64, pressure: f64) -> Self {
        Self {
            temperature,
            pressure,
        }
    }

    pub fn total_energy(&self) -> f64 {
        // Implementation
        0.0
    }
}
```

### Python Code

- **Format:** Use `black` with line length 100
- **Lint:** Fix all `ruff` warnings
- **Type hints:** Add type hints for public APIs
- **Docstrings:** Use NumPy-style docstrings, with `.. math::` for formulas
- **Naming:** Standard Python conventions with scientific exceptions
  - Regular variables: `snake_case`
  - Scientific variables: `kT`, `pH`, `deltaG` allowed

**Example:**
```python
def calculate_energy(
    positions: np.ndarray,
    charges: np.ndarray,
    cutoff: float = 1.2,
) -> float:
    """
    Calculate total electrostatic energy.

    Parameters
    ----------
    positions : np.ndarray
        Atomic positions (N, 3)
    charges : np.ndarray
        Atomic charges (N,)
    cutoff : float, optional
        Cutoff distance in nm (default: 1.2)

    Returns
    -------
    float
        Total energy in kJ/mol
    """
    # Implementation
    return 0.0
```

### Documentation

- **Rust:** Add doc comments with `///` for public APIs
- **Python:** Add docstrings for all public functions/classes
- **Examples:** Include usage examples in documentation
- **Markdown:** Follow markdown best practices

## Testing

### Writing Tests

**Rust tests:**
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feature() {
        let result = my_function(42);
        assert_eq!(result, 84);
    }
}
```

**Python tests:**
```python
def test_feature():
    result = gromos.my_function(42)
    assert result == 84
```

### Running Tests

```bash
# All tests
make test

# Rust only
make test-rust

# Python only
make test-python

# With coverage
make test-coverage

# Specific test
cd gromos-rs && cargo test test_name
cd py-gromos && pytest tests/test_file.py::test_function
```

### Test Requirements

- All new features must have tests
- Bug fixes should include regression tests
- Tests should be fast (< 1 second per test)
- Mark slow tests with `#[ignore]` (Rust) or `@pytest.mark.slow` (Python)

## Submitting Changes

### Pull Request Checklist

Before submitting a PR, ensure:

- [ ] Code is formatted (`make fmt`)
- [ ] All linters pass (`make lint`)
- [ ] All tests pass (`make test`)
- [ ] New tests are added for new features
- [ ] Documentation is updated
- [ ] Commit messages follow conventional commits
- [ ] PR description explains the changes

### Pull Request Process

1. **Create PR** with clear title and description
2. **CI checks** will run automatically
3. **Address review comments** if any
4. **Squash commits** if requested
5. **Maintainer will merge** once approved

### PR Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
How has this been tested?

## Checklist
- [ ] Code follows project style
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] All CI checks pass
```

## Pre-commit Hooks

We use pre-commit hooks to ensure code quality. Install them:

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run manually on all files
pre-commit run --all-files
```

The hooks will:
- Format Rust code with `rustfmt`
- Lint Rust code with `clippy`
- Format Python code with `black`
- Lint Python code with `ruff`
- Check for typos
- Check YAML/TOML syntax
- Check for large files
- Check for merge conflicts

## Code Review Guidelines

### For Contributors

- Be responsive to feedback
- Keep PRs focused and small
- Write clear commit messages
- Update documentation

### For Reviewers

- Be constructive and respectful
- Check for correctness
- Check for test coverage
- Check for documentation
- Consider performance implications

## Questions?

- Open an issue for bugs or feature requests
- Join discussions on GitHub Discussions
- Check the [CI/CD Guide](gromos-rs/docs/development/ci-cd.md) for CI/CD info
- See [Developer Tools](gromos-rs/docs/development/developer-tools.md) for tooling setup
- Read [Pre-commit Hooks](gromos-rs/docs/development/pre-commit.md) for hook configuration

## License

By contributing, you agree that your contributions will be licensed under the GPL-2.0 license.
