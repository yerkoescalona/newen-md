#!/usr/bin/env python3
"""
Engine-vs-engine benchmark driver: newen `md` vs gromosXX `md`.

Implements the protocol in BENCHMARKING.md. Patches a reference `.imd` to the
requested step count with all trajectory output disabled, runs each engine
pinned to a fixed core set, and reports median/min wall time per engine.

Repeats are **interleaved across engines** (A, B, A, B, ...) rather than run
as blocks, so slow thermal or boost-clock drift hits both engines equally
instead of penalising whichever ran first (BENCHMARKING.md assumption A7).

Usage
-----
    scripts/bench_engines.py --system water_216_box --steps 2000 --repeats 5
    scripts/bench_engines.py --system water_1000_spc_gridcell --threads 1,2,4,8
    scripts/bench_engines.py --list
    scripts/bench_engines.py --engines rust --rust-md bin/md_old --system ch4_water_fep

`--rust-md PATH` (repeatable) adds another newen build as an engine named
`rust-<file name>`, run on the same patched input and interleaved with the
others, so two Rust builds can be compared directly.

Wall time is taken from each engine's own report (newen "Simulation wall
time", gromosXX "Wall time simulation") when available, so start-up parsing is
excluded; the externally measured total is recorded alongside it.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REF_DIR = REPO / "crates" / "newen-md" / "tests" / "gromosXX_references"
RUST_MD = REPO / "target" / "release" / "md"

# gromosXX build variants, in preference order for the serial baseline.
GROMOSXX_BUILDS = {
    "native": REPO / ".local/gromosXX/md++/BUILD_native/program/md",
    "plain": REPO / ".local/gromosXX/md++/BUILD/program/md",
    "omp": REPO / ".local/gromosXX/md++/BUILD_omp/program/md",
    "mpi": REPO / ".local/gromosXX/md++/BUILD_mpi/program/md_mpi",
}

# Auxiliary inputs some systems need (FEP perturbation topology, restraint
# specs); same input.toml keys and flags the reference generator uses.
AUX_INPUTS = (
    ("pttopo", "@pttopo"),
    ("posresspec", "@posresspec"),
    ("refpos", "@refpos"),
    ("distrest", "@distrest"),
)


# ── input patching ──────────────────────────────────────────────────────────


def patch_imd(text: str, nstlim: int, energy_every: int = 0) -> str:
    """Set NSTLIM and zero every WRITETRAJ / PRINTOUT frequency.

    Trajectory writing is I/O bound and would dominate a short run
    (BENCHMARKING.md assumption A4), so the benchmark input writes nothing —
    except NTWE = NTWG = `energy_every` when it is set (a PRECALCLAM run needs
    them, and the energies at the extra λ points are computed on those steps).
    """
    out: list[str] = []
    block = None
    data_seen = 0

    for line in text.splitlines():
        stripped = line.strip()

        if stripped in ("STEP", "WRITETRAJ", "PRINTOUT"):
            block, data_seen = stripped, 0
            out.append(line)
            continue
        if stripped == "END":
            block = None
            out.append(line)
            continue

        if block and stripped and not stripped.startswith("#"):
            data_seen += 1
            if block == "STEP" and data_seen == 1:
                parts = stripped.split()
                if len(parts) >= 3:
                    out.append(f"      {nstlim}       {parts[1]}     {parts[2]}")
                    continue
            elif block == "WRITETRAJ" and data_seen == 1:
                e = energy_every
                out.append(f"         0         0         0         0         {e}         {e}         0")
                continue
            elif block == "PRINTOUT" and data_seen == 1:
                out.append(f"     {nstlim}         0")
                continue

        out.append(line)

    return "\n".join(out) + "\n"


def toml_get(text: str, key: str) -> str | None:
    m = re.search(rf'^{key}\s*=\s*"([^"]+)"', text, re.MULTILINE)
    return m.group(1) if m else None


# ── timing extraction ───────────────────────────────────────────────────────


def rust_wall(stderr: str) -> float | None:
    m = re.search(r"Simulation wall time:\s*([0-9.]+)\s*s", stderr)
    return float(m.group(1)) if m else None


def gromosxx_wall(stdout: str) -> float | None:
    # md++ epilogue: "Wall time simulation (s):   1.190" — excludes initialisation,
    # matching what newen reports as "Simulation wall time".
    m = re.search(r"Wall time simulation \(s\):\s*([0-9.]+)", stdout)
    return float(m.group(1)) if m else None


# ── running ─────────────────────────────────────────────────────────────────


@dataclass
class Prepared:
    """One (engine, system, threads) configuration ready to run repeatedly."""

    engine: str
    kind: str  # "rust" | "gromosXX"
    cmd: list[str]
    env: dict[str, str]
    work: Path


@dataclass
class Run:
    engine: str
    system: str
    steps: int
    threads: int
    externals: list[float] = field(default_factory=list)
    internals: list[float | None] = field(default_factory=list)

    def _vals(self) -> list[float]:
        vals = [v for v in self.internals if v is not None]
        return vals or self.externals

    def n(self) -> int:
        return len(self._vals())

    def mean(self) -> float:
        return statistics.fmean(self._vals())

    def std(self) -> float:
        """Sample standard deviation (n-1); 0 for a single value."""
        vals = self._vals()
        return statistics.stdev(vals) if len(vals) > 1 else 0.0

    def best(self) -> float:
        return min(self._vals())

    def median(self) -> float:
        return statistics.median(self._vals())

    def spread(self) -> float:
        vals = self._vals()
        return (max(vals) - min(vals)) / statistics.median(vals) * 100 if vals else 0.0


def prepare(engine: str, kind: str, binary: Path, system: str, steps: int, threads: int,
            ref_dir: Path = REF_DIR, energy_every: int = 0) -> Prepared:
    sys_dir = ref_dir / system
    toml = (sys_dir / "input.toml").read_text()
    topo = (sys_dir / toml_get(toml, "topology")).resolve()
    conf = (sys_dir / toml_get(toml, "configuration")).resolve()
    params = sys_dir / toml_get(toml, "parameters")

    work = REPO / "bench" / "work" / f"{engine}_{system}_{threads}t"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)

    imd = work / "bench.imd"
    imd.write_text(patch_imd(params.read_text(), steps, energy_every))

    env = dict(os.environ)
    env["RAYON_NUM_THREADS"] = str(threads)
    env["OMP_NUM_THREADS"] = str(threads)
    if kind == "rust":
        env["RUST_LOG"] = "info"

    cmd = [
        str(binary),
        "@topo", str(topo),
        "@conf", str(conf),
        "@input", str(imd),
        "@fin", str(work / "final.cnf"),
    ]
    for key, flag in AUX_INPUTS:
        aux = toml_get(toml, key)
        if aux:
            cmd += [flag, str((sys_dir / aux).resolve())]
    if energy_every:
        cmd += ["@tre", str(work / "bench.tre"), "@trg", str(work / "bench.trg")]
    if kind == "rust":
        cmd += ["@verb", "1"]

    # Pin to physical cores 0..threads-1 (CPUs 8-15 are SMT siblings here).
    cpus = ",".join(str(c) for c in range(threads))
    cmd = ["taskset", "-c", cpus] + cmd

    return Prepared(engine, kind, cmd, env, work)


def run_once(p: Prepared, system: str) -> tuple[float, float | None]:
    t0 = time.perf_counter()
    proc = subprocess.run(p.cmd, capture_output=True, text=True, cwd=p.work, env=p.env)
    elapsed = time.perf_counter() - t0

    if proc.returncode != 0:
        sys.stderr.write(
            f"\n{p.engine} failed on {system} (exit {proc.returncode}):\n"
            f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}\n"
        )
        raise SystemExit(1)

    internal = rust_wall(proc.stderr) if p.kind == "rust" else gromosxx_wall(proc.stdout)
    return elapsed, internal


def run_interleaved(preps: list[Prepared], system: str, steps: int, threads: int, repeats: int) -> list[Run]:
    runs = {p.engine: Run(p.engine, system, steps, threads) for p in preps}

    # One warm-up each (page cache, allocator), untimed.
    for p in preps:
        run_once(p, system)

    for rep in range(repeats):
        for p in preps:
            print(f"  rep {rep + 1}/{repeats}  {p.engine:<18} {system:<26} {threads} thread(s)",
                  file=sys.stderr, flush=True)
            external, internal = run_once(p, system)
            runs[p.engine].externals.append(external)
            runs[p.engine].internals.append(internal)

    return [runs[p.engine] for p in preps]


# ── reporting ───────────────────────────────────────────────────────────────


def report(runs: list[Run]) -> None:
    print()
    print(f"{'engine':<16} {'system':<26} {'thr':>4} {'n':>3} {'mean s':>9} {'± std':>8} "
          f"{'median s':>9} {'min s':>8} {'us/step':>9}")
    print("-" * 100)
    for r in runs:
        print(f"{r.engine:<16} {r.system:<26} {r.threads:>4} {r.n():>3} {r.mean():>9.3f} "
              f"{r.std():>8.3f} {r.median():>9.3f} {r.best():>8.3f} {r.mean()/r.steps*1e6:>9.1f}")

    by_key: dict[tuple[str, int], dict[str, Run]] = {}
    for r in runs:
        by_key.setdefault((r.system, r.threads), {})[r.engine] = r

    pairs = [(k, v) for k, v in by_key.items() if len(v) > 1]
    if pairs:
        print()
        print(f"{'system':<26} {'thr':>4} {'other':<24} {'other/rust (mean)':>18} "
              f"{'± (propagated)':>15} {'(median)':>9}")
        print("-" * 101)
        for (system, threads), engines in pairs:
            rust = engines.get("rust")
            if not rust:
                continue
            for name, other in engines.items():
                if name == "rust":
                    continue
                ratio = other.mean() / rust.mean()
                # First-order error propagation on the ratio of two means.
                rel = ((other.std() / other.mean()) ** 2 + (rust.std() / rust.mean()) ** 2) ** 0.5
                print(f"{system:<26} {threads:>4} {name:<24} {ratio:>17.3f}x {ratio * rel:>14.3f} "
                      f"{other.median()/rust.median():>8.2f}x")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", "-s", action="append",
                    help="reference system name (repeatable)")
    ap.add_argument("--steps", type=int, default=2000, help="NSTLIM (default 2000)")
    ap.add_argument("--repeats", "-r", type=int, default=5, help="timed repeats (default 5)")
    ap.add_argument("--threads", "-t", default="1",
                    help="comma-separated thread counts (default 1)")
    ap.add_argument("--engines", "-e", default="rust,native",
                    help="comma-separated: rust + any of " + ",".join(GROMOSXX_BUILDS))
    ap.add_argument("--rust-md", type=Path, action="append", default=[],
                    help="another newen md binary to time as engine rust-<file name> "
                         "(repeatable)")
    ap.add_argument("--ref-dir", type=Path, default=REF_DIR,
                    help="directory of system subdirectories with input.toml "
                         "(default: the gromosXX reference suite; use bench/systems for generated boxes)")
    ap.add_argument("--energy-every", type=int, default=0,
                    help="keep NTWE = NTWG = N and write the .tre/.trg (both engines); "
                         "default 0 writes nothing")
    ap.add_argument("--list", action="store_true", help="list available systems")
    ap.add_argument("--json", type=Path, help="write raw results to this JSON file")
    args = ap.parse_args()

    if args.list:
        for d in sorted(args.ref_dir.iterdir()):
            if (d / "input.toml").exists():
                print(d.name)
        return

    systems = args.system or ["water_216_box"]
    threads = [int(t) for t in args.threads.split(",")]
    engines = args.engines.split(",")

    binaries: dict[str, tuple[str, Path]] = {}
    for e in engines:
        if e == "rust":
            binaries["rust"] = ("rust", RUST_MD)
        elif e in GROMOSXX_BUILDS:
            binaries[f"gromosXX-{e}"] = ("gromosXX", GROMOSXX_BUILDS[e])
        else:
            raise SystemExit(f"unknown engine: {e}")
    for path in args.rust_md:
        binaries[f"rust-{path.name}"] = ("rust", path.resolve())

    for name, (_, path) in binaries.items():
        if not path.exists():
            raise SystemExit(f"{name} binary not found: {path}")

    gov = Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
    if gov.exists() and gov.read_text().strip() != "performance":
        print(f"WARNING: cpu governor is '{gov.read_text().strip()}', not 'performance' "
              f"— timings will be noisier (BENCHMARKING.md A7)", file=sys.stderr)

    runs: list[Run] = []
    for system in systems:
        for t in threads:
            preps = [prepare(name, kind, path, system, args.steps, t, args.ref_dir,
                             args.energy_every)
                     for name, (kind, path) in binaries.items()]
            runs.extend(run_interleaved(preps, system, args.steps, t, args.repeats))

    report(runs)

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps([asdict(r) for r in runs], indent=2))
        print(f"\nraw results → {args.json}")


if __name__ == "__main__":
    main()
