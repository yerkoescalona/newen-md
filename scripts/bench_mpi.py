#!/usr/bin/env python3
"""MPI scaling of newen `md` (feature `use-mpi`): wall time and energy agreement across
`mpirun -np N` for one system (BENCHMARKING.md §6).

    scripts/bench_mpi.py --system water_1000_spc_gridcell --np 1,2,4 --steps 200 --repeats 3
    scripts/bench_mpi.py --ref-dir bench/systems --system ch4_water_7975 --np 1,2,4,8

Needs `LIBCLANG_PATH=/usr/lib/llvm-19/lib cargo build --release --features use-mpi --bin md`
(the binary is target/release/md; a serial build of the same binary is the np=1 baseline —
pass `--serial <path>` to compare against a build without the feature).
Every run reports the median wall time over the repeats (n and std printed), and the maximum
relative deviation of the last energy frame (.tre totals) from the np=1 run.
"""
from __future__ import annotations

import argparse
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REF_DIR = REPO / "crates/newen-md/tests/gromosXX_references"
AUX_INPUTS = [("pttopo", "@pttopo"), ("posresspec", "@posresspec"), ("refpos", "@refpos"), ("distrest", "@distrest")]


def toml_get(text: str, key: str) -> str | None:
    m = re.search(rf'^\s*{key}\s*=\s*"([^"]+)"', text, re.M)
    return m.group(1) if m else None


def patch_imd(text: str, nstlim: int) -> str:
    """NSTLIM ← nstlim, energies every step, no trajectories, as bench_engines.py does."""
    out, block = [], None
    for line in text.split("\n"):
        t = line.strip()
        if t and not t.startswith("#") and t.isupper() and block is None:
            block = t
        elif t == "END":
            block = None
        elif block == "STEP" and not t.startswith("#"):
            parts = t.split()
            if len(parts) >= 3 and parts[0].isdigit():
                line = f"{nstlim:>10} {parts[1]:>15} {parts[2]:>15}"
        elif block == "WRITETRAJ" and not t.startswith("#"):
            parts = t.split()
            if len(parts) >= 7:
                line = "         0        0        0        0        1        0        0"
        out.append(line)
    return "\n".join(out)


def last_totals(tre: Path) -> list[float]:
    """The `# totals` numbers of the last ENERGY03 block."""
    text = tre.read_text().split("ENERGY03")[-1]
    section = text.split("# totals", 1)[1] if "# totals" in text else text
    vals = []
    for line in section.split("\n"):
        t = line.strip()
        if t.startswith("#") and vals:
            break
        if t.startswith("#") or not t:
            continue
        if t == "END":
            break
        vals += [float(x) for x in t.split()]
    return vals


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", required=True)
    ap.add_argument("--ref-dir", type=Path, default=REF_DIR)
    ap.add_argument("--np", default="1,2,4")
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--binary", type=Path, default=REPO / "target/release/md")
    ap.add_argument("--serial", type=Path, help="md built without use-mpi, as the reference")
    ap.add_argument("--threads", type=int, default=1, help="RAYON threads per rank")
    a = ap.parse_args()

    sys_dir = a.ref_dir / a.system
    toml = (sys_dir / "input.toml").read_text()
    topo = (sys_dir / toml_get(toml, "topology")).resolve()
    conf = (sys_dir / toml_get(toml, "configuration")).resolve()
    params = sys_dir / toml_get(toml, "parameters")
    work = REPO / "bench" / "work" / f"mpi_{a.system}"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    imd = work / "bench.imd"
    imd.write_text(patch_imd(params.read_text(), a.steps))

    def cmd_for(binary: Path, np: int, tag: str) -> list[str]:
        out = work / tag
        out.mkdir(exist_ok=True)
        c = [str(binary), "@topo", str(topo), "@conf", str(conf), "@input", str(imd),
             "@fin", str(out / "final.cnf"), "@tre", str(out / "energies.tre"), "@verb", "1"]
        for key, flag in AUX_INPUTS:
            aux = toml_get(toml, key)
            if aux:
                c += [flag, str((sys_dir / aux).resolve())]
        if np > 1 or binary == a.binary:
            c = ["mpirun", "--bind-to", "core", "-np", str(np)] + c
        return c

    env = {"RAYON_NUM_THREADS": str(a.threads), "OMP_NUM_THREADS": str(a.threads), "RUST_LOG": "info",
           "PATH": "/usr/bin:/bin:/usr/local/bin"}
    import os
    env = {**os.environ, **env}
    rows = []
    baseline = None
    configs = []
    if a.serial:
        configs.append(("serial", a.serial, 1))
    configs += [(f"np{n}", a.binary, int(n)) for n in a.np.split(",")]
    for tag, binary, np in configs:
        walls = []
        for _ in range(a.repeats):
            r = subprocess.run(cmd_for(binary, np, tag), capture_output=True, text=True, env=env, cwd=work)
            if r.returncode != 0:
                print(f"{tag}: exit {r.returncode}\n{r.stderr[-2000:]}", file=sys.stderr)
                return 1
            m = re.search(r"Simulation wall time:\s*([0-9.]+)\s*s", r.stderr)
            if not m:
                print(f"{tag}: no wall time in stderr", file=sys.stderr)
                return 1
            walls.append(float(m.group(1)))
        totals = last_totals(work / tag / "energies.tre")
        if baseline is None:
            baseline = totals
            dev = 0.0
        else:
            dev = max(abs(x - y) / max(abs(y), 1e-12) for x, y in zip(totals, baseline) if y != 0.0)
        med = statistics.median(walls)
        std = statistics.pstdev(walls) if len(walls) > 1 else 0.0
        rows.append((tag, np, med, std, len(walls), dev, totals[0]))
    base_t = rows[0][2]
    print(f"| run | np | wall s (median) | std | n | speedup vs first | max rel. dev. of last-frame energies vs first | E_total |")
    print("|---|--:|--:|--:|--:|--:|--:|--:|")
    for tag, np, med, std, n, dev, etot in rows:
        print(f"| {tag} | {np} | {med:.2f} | {std:.2f} | {n} | {base_t / med:.2f} | {dev:.1e} | {etot:.6e} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
