#!/usr/bin/env python3
"""
Turn `bench_engines.py --threads ...` JSON output into the Phase 3 tables of
BENCHMARKING.md: per engine, mean ± std at each thread count, speedup
S(k) = t(1)/t(k), parallel efficiency E(k) = S(k)/k, and the cross-engine ratio.

    scripts/scaling_table.py bench/results/phase3_water_1000.json [more.json ...]

Several files may be given for the same system (e.g. k=1 in one run, k=2,4,8 in
another); runs are merged by (system, engine, threads).
"""

from __future__ import annotations

import json
import statistics as st
import sys
from collections import defaultdict


def load(paths: list[str]):
    runs: dict[tuple[str, str, int], list[float]] = defaultdict(list)
    for path in paths:
        for r in json.load(open(path)):
            vals = [v for v in r["internals"] if v is not None] or r["externals"]
            runs[(r["system"], r["engine"], r["threads"])].extend(vals)
    return runs


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    runs = load(sys.argv[1:])
    systems = sorted({k[0] for k in runs})
    for system in systems:
        # Rust first so the cross-engine column reads gromosXX / Rust (> 1 = newen faster).
        engines = sorted({k[1] for k in runs if k[0] == system}, key=lambda e: (e != "rust", e))
        threads = sorted({k[2] for k in runs if k[0] == system})
        print(f"\n**{system}** (n per cell = timed repeats; mean ± std s)\n")
        header = "| k |"
        for e in engines:
            header += f" {e} s | S(k) | E(k) |"
        if len(engines) == 2:
            header += f" {engines[1]} / {engines[0]} (time; >1 = {engines[0]} faster) |"
        print(header)
        print("|" + "--:|" * (header.count("|") - 1))
        base = {e: st.fmean(runs[(system, e, 1)]) for e in engines if (system, e, 1) in runs}
        for k in threads:
            row = f"| {k} |"
            means = {}
            for e in engines:
                key = (system, e, k)
                if key not in runs:
                    row += " — | — | — |"
                    continue
                v = runs[key]
                m = st.fmean(v)
                sd = st.stdev(v) if len(v) > 1 else 0.0
                means[e] = m
                if e in base:
                    s = base[e] / m
                    row += f" {m:.2f} ± {sd:.2f} (n={len(v)}) | {s:.2f} | {s / k:.2f} |"
                else:
                    row += f" {m:.2f} ± {sd:.2f} (n={len(v)}) | — | — |"
            if len(engines) == 2 and all(e in means for e in engines):
                row += f" {means[engines[1]] / means[engines[0]]:.2f}× |"
            print(row)


if __name__ == "__main__":
    main()
