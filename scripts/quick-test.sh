#!/bin/bash
# Quick test script - runs fast subset of tests

set -e

echo "==> Running quick test suite..."

# Test Rust library (library only, faster)
echo "==> Testing newen (lib only)..."
cd newen
cargo test --lib
cd ..

# Test Python bindings (basic tests only)
echo "==> Testing py-newen (basic tests)..."
cd py-newen
if [ ! -d ".venv" ]; then
    echo "Virtual environment not found. Run 'make venv' first."
    exit 1
fi
.venv/bin/pytest tests/test_basic.py -v
cd ..

echo ""
echo "==> Quick tests passed! ✓"
echo ""
echo "For comprehensive testing, run: make test-all"
