#!/bin/bash
# Run the newen-md criterion benchmarks and compare with a saved baseline.
#
#   scripts/benchmark.sh            # save baseline "main"
#   scripts/benchmark.sh <name>     # compare against baseline <name>
#
# For the engine-vs-engine (gromosXX) protocol see BENCHMARKING.md.

set -e

cd "$(dirname "$0")/.."

BASELINE=${1:-main}

echo "==> Running newen-md benchmarks..."
if [ "$BASELINE" != "main" ]; then
    echo "Comparing against baseline: $BASELINE"
    cargo bench -p newen-md -- --baseline "$BASELINE"
else
    echo "Saving baseline as: $BASELINE"
    cargo bench -p newen-md -- --save-baseline "$BASELINE"
fi

echo ""
echo "==> Benchmarks complete!"
echo "Results saved in target/criterion/"
