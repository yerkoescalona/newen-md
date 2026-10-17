#!/bin/bash
# Development environment setup script

set -e

echo "==> Setting up gromosXX development environment..."

# Check for Rust
if ! command -v cargo &> /dev/null; then
    echo "Rust not found. Installing Rust..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
else
    echo "✓ Rust found: $(rustc --version)"
fi

# Check for Python
if ! command -v python3 &> /dev/null; then
    echo "Python 3 not found. Please install Python 3.9 or later."
    exit 1
else
    echo "✓ Python found: $(python3 --version)"
fi

# Install Rust components
echo "==> Installing Rust components..."
rustup component add rustfmt clippy

# Install Python tools
echo "==> Installing Python tools..."
python3 -m pip install --upgrade pip
pip install maturin black ruff pytest

# Optional: Install cargo tools
echo "==> Installing optional Cargo tools..."
echo "This may take a while..."
cargo install cargo-edit cargo-audit cargo-outdated cargo-llvm-cov 2>/dev/null || echo "Some tools already installed"

# Build newen
echo "==> Building newen..."
cd newen
cargo build
cd ..

# Setup py-newen
echo "==> Setting up py-newen..."
cd py-newen
make venv
make install-dev
make build
cd ..

echo ""
echo "==> Development environment setup complete! ✓"
echo ""
echo "Quick commands:"
echo "  make help           - Show all available commands"
echo "  make build          - Build both projects"
echo "  make test           - Run all tests"
echo "  make pre-commit     - Run pre-commit checks"
echo ""
echo "See CI_CD_GUIDE.md for detailed documentation."
