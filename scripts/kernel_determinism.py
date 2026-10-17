#!/usr/bin/env python3
"""Measure run-to-run and thread-count determinism of the newen `md` binary.

3.9 step 0 / assumption A2–A3: the parallel nonbonded kernels reduce through rayon
`fold/reduce`, whose merge tree depends on work-stealing, so results may differ (a) between
thread counts and (b) between two runs at the *same* thread count. This script measures both
on reference systems and prints a Markdown table (n, max |Δ| and max relative Δ per energy
column) for BENCHMARKING.md.

Usage (from the repo root, after `cargo build --release --bin md`):
    python3 scripts/kernel_determinism.py --system water_216_box --system ch4_water_fep \
        --threads 1,8 --repeats 3 --steps 100
"""

from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
REF_DIR = REPO_ROOT / "crates" / "newen-md" / "tests" / "gromosXX_references"
MD_BIN = REPO_ROOT / "target" / "release" / "md"

# ENERTRJ column layout written by crates/newen-md/src/bin/md.rs (see the header it emits).
COLUMNS = [
    "time", "kinetic", "potential", "total", "temperature", "volume", "pressure",
    "bond", "angle", "improper", "dihedral", "lj", "coulomb_real", "coulomb_recip",
    "coulomb_self", "shake", "distres",
]


def toml_str(text: str, key: str) -> str | None:
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(f"{key} ="):
            value = line.split("=", 1)[1].strip()
            return value.strip('"')
    return None


def patch_steps(imd_text: str, steps: int) -> str:
    """Set NSTLIM in the STEP block (first value of the data line after the comment)."""
    out, in_step = [], False
    for line in imd_text.splitlines():
        if line.strip() == "STEP":
            in_step = True
        elif in_step and not line.startswith("#") and line.strip() and line.strip() != "END":
            parts = line.split()
            parts[0] = str(steps)
            line = "    " + "    ".join(parts)
            in_step = False
        out.append(line)
    return "\n".join(out) + "\n"


def parse_enertrj(path: Path) -> np.ndarray:
    rows = []
    in_block = False
    for line in path.read_text().splitlines():
        s = line.strip()
        if s == "ENERTRJ":
            in_block = True
            continue
        if s == "END":
            in_block = False
            continue
        if in_block and s and not s.startswith("#"):
            rows.append([float(x) for x in s.split()])
    return np.array(rows, dtype=np.float64)


def run_md(system: str, threads: int, steps: int | None, workdir: Path) -> np.ndarray:
    sys_dir = REF_DIR / system
    toml = (sys_dir / "input.toml").read_text()
    topo = sys_dir / toml_str(toml, "topology")
    conf = sys_dir / toml_str(toml, "configuration")
    params = sys_dir / toml_str(toml, "parameters")
    if steps is not None:
        patched = workdir / f"{system}.in"
        patched.write_text(patch_steps(params.read_text(), steps))
        params = patched
    tre = workdir / f"{system}_t{threads}.tre"
    cmd = [str(MD_BIN), "@topo", str(topo), "@conf", str(conf), "@input", str(params),
           "@tre", str(tre), "@trc", str(workdir / "traj.trc")]
    for key, flag in (("pttopo", "@pttopo"), ("posresspec", "@posresspec"),
                      ("refpos", "@refpos"), ("distrest", "@distrest")):
        aux = toml_str(toml, key)
        if aux:
            cmd += [flag, str(sys_dir / aux)]
    env = dict(os.environ, RAYON_NUM_THREADS=str(threads))
    result = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=workdir)
    if result.returncode != 0:
        sys.exit(f"{system} @ {threads} threads failed:\n{result.stderr[-2000:]}")
    return parse_enertrj(tre)


def compare(a: np.ndarray, b: np.ndarray) -> tuple[str, float, float]:
    """(column, max |Δ|, max relative Δ) over all steps, worst column by relative Δ."""
    n = min(len(a), len(b))
    a, b = a[:n, 1:], b[:n, 1:]
    diff = np.abs(a - b)
    scale = np.maximum(np.abs(a), np.abs(b))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(scale > 0, diff / scale, 0.0)
    col = int(np.argmax(rel.max(axis=0)))
    return COLUMNS[col + 1], float(diff[:, col].max()), float(rel[:, col].max())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", action="append", required=True)
    ap.add_argument("--threads", default="1,8")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--steps", type=int, default=None, help="override NSTLIM")
    args = ap.parse_args()
    thread_counts = [int(t) for t in args.threads.split(",")]

    if not MD_BIN.exists():
        sys.exit(f"missing {MD_BIN}; run `cargo build --release --bin md`")

    print("| system | comparison | n pairs | bit-identical | worst column | max abs Δ | max rel Δ |")
    print("|---|---|---|---|---|---|---|")
    with tempfile.TemporaryDirectory(prefix="newen_determinism_") as tmp:
        tmp = Path(tmp)
        for system in args.system:
            runs: dict[int, list[np.ndarray]] = {}
            for t in thread_counts:
                runs[t] = [run_md(system, t, args.steps, tmp) for _ in range(args.repeats)]
            # run-to-run at fixed thread count
            for t in thread_counts:
                pairs = list(itertools.combinations(runs[t], 2))
                results = [compare(a, b) for a, b in pairs]
                identical = all(np.array_equal(a, b) for a, b in pairs)
                worst = max(results, key=lambda r: r[2])
                print(f"| {system} | {t} thr vs {t} thr (run-to-run) | {len(pairs)} | "
                      f"{'yes' if identical else 'NO'} | {worst[0]} | {worst[1]:.3e} | {worst[2]:.3e} |")
            # across thread counts
            for t1, t2 in itertools.combinations(thread_counts, 2):
                pairs = [(a, b) for a in runs[t1] for b in runs[t2]]
                results = [compare(a, b) for a, b in pairs]
                identical = all(np.array_equal(a, b) for a, b in pairs)
                worst = max(results, key=lambda r: r[2])
                print(f"| {system} | {t1} thr vs {t2} thr | {len(pairs)} | "
                      f"{'yes' if identical else 'NO'} | {worst[0]} | {worst[1]:.3e} | {worst[2]:.3e} |")


if __name__ == "__main__":
    main()
