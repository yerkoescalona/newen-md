#!/usr/bin/env python3
"""Does gromosXX accept what newen writes? (3.9 assumption A7, guard G3)

For every reference system: read its `.in` with newen, write it back with `write_imd`
(`cargo run -p newen-io --example imd_roundtrip`), run the *real* gromosXX on the rewritten
file with the system's own topology/coordinates/auxiliary inputs, and compare the energies it
produces with a fresh gromosXX run of the **original** `.in` made in the same session — the
numbers must be identical, because the rewritten file describes the same run.

Why not the committed `expected/`: for 9 of the 40 systems gromosXX's output on the *original*
input differs from `expected/` today in ~1e-32 near-zero terms (centre-of-mass kinetic
energies), while two fresh runs agree exactly — noise between the environment that generated
the references and this one, invisible at the 1e-8 reference tolerance. The comparison to
`expected/` is reported as an informational column.

Usage (from the repo root):
    python3 scripts/roundtrip_imd_gromosXX.py [--md-binary PATH] [--system NAME ...]

The default binary is the CODATA-patched `.local/gromosXX/md++/BUILD/program/md`, the one the
committed references were generated with (BUILD_native/BUILD_omp differ in the last digit).
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
REF_DIR = REPO_ROOT / "crates" / "newen-md" / "tests" / "gromosXX_references"
DEFAULT_MD = REPO_ROOT / ".local" / "gromosXX" / "md++" / "BUILD" / "program" / "md"
ROUNDTRIP = REPO_ROOT / "target" / "release" / "examples" / "imd_roundtrip"


def toml_str(text: str, key: str) -> str | None:
    m = re.search(rf'^{key}\s*=\s*"(.+)"', text, re.M)
    return m.group(1) if m else None


def numeric_lines(path: Path) -> list[str]:
    """Energy-file lines that carry numbers, with the title dropped (it is free text)."""
    out = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s in ("TITLE", "END", "ENERGY03", "VOLUMEPRESSURE03"):
            continue
        if re.match(r"^[-+0-9.eE ]+$", s):
            out.append(s)
    return out


def build_tool() -> None:
    subprocess.run(
        ["cargo", "build", "--release", "-p", "newen-io", "--example", "imd_roundtrip"],
        check=True, cwd=REPO_ROOT, capture_output=True,
    )


def run_gromosxx(md_bin: Path, sys_dir: Path, toml: str, params: Path, workdir: Path) -> Path | None:
    """Run gromosXX with `params` in `workdir`; returns the energy file or None on failure."""
    topo = sys_dir / toml_str(toml, "topology")
    conf = sys_dir / toml_str(toml, "configuration")
    workdir.mkdir(parents=True, exist_ok=True)
    tre = workdir / "e.tre"
    cmd = [str(md_bin), "@topo", str(topo), "@conf", str(conf), "@input", str(params),
           "@fin", str(workdir / "fin.cnf"), "@tre", str(tre), "@trf", str(workdir / "f.trf"),
           "@trc", str(workdir / "t.trc"), "@trg", str(workdir / "g.trg")]
    for key in ("posresspec", "refpos", "distrest", "pttopo"):
        aux = toml_str(toml, key)
        if aux:
            cmd += [f"@{key}", str(sys_dir / aux)]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    if res.returncode != 0 or not tre.exists():
        tail = (res.stderr or res.stdout).strip().splitlines()[-1:] or ["?"]
        (workdir / "FAILED").write_text(tail[0])
        return None
    return tre


def differ(a: list[str], b: list[str]) -> str:
    if a == b:
        return "yes"
    n = sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))
    return f"NO ({n} lines differ, {len(a)} vs {len(b)})"


def run_system(sys_dir: Path, md_bin: Path) -> tuple[str, str, str, str]:
    """Returns (system, accepted, identical to fresh original run, identical to expected/)."""
    name = sys_dir.name
    toml = (sys_dir / "input.toml").read_text()
    params = sys_dir / toml_str(toml, "parameters")
    expected = sys_dir / "expected" / "energies.tre"
    with tempfile.TemporaryDirectory(prefix=f"imd_rt_{name}_") as tmp:
        tmp = Path(tmp)
        rewritten = tmp / f"{name}.rewritten.in"
        rt = subprocess.run([str(ROUNDTRIP), str(params), str(rewritten)], capture_output=True, text=True)
        if rt.returncode != 0:
            return name, f"write failed: {rt.stderr.strip()[:80]}", "-", "-"
        got_tre = run_gromosxx(md_bin, sys_dir, toml, rewritten, tmp / "rewritten")
        if got_tre is None:
            return name, f"REJECTED: {(tmp / 'rewritten' / 'FAILED').read_text()[:90]}", "-", "-"
        ref_tre = run_gromosxx(md_bin, sys_dir, toml, params, tmp / "original")
        if ref_tre is None:
            return name, "yes", "original .imd itself failed?!", "-"
        got = numeric_lines(got_tre)
        vs_fresh = differ(got, numeric_lines(ref_tre))
        vs_expected = differ(got, numeric_lines(expected)) if expected.exists() else "no expected/"
        return name, "yes", vs_fresh, vs_expected


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md-binary", type=Path, default=DEFAULT_MD)
    ap.add_argument("--system", action="append")
    args = ap.parse_args()
    if not args.md_binary.exists():
        sys.exit(f"gromosXX binary not found: {args.md_binary}")
    build_tool()
    dirs = sorted(p for p in REF_DIR.iterdir() if p.is_dir() and (p / "input.toml").exists())
    if args.system:
        dirs = [d for d in dirs if d.name in set(args.system)]
    print(f"gromosXX: {args.md_binary}")
    print("| system | gromosXX accepts rewritten .imd | energies identical to a fresh run of the original .imd | identical to committed expected/ (informational) |")
    print("|---|---|---|---|")
    bad = 0
    for d in dirs:
        name, accepted, vs_fresh, vs_expected = run_system(d, args.md_binary)
        print(f"| {name} | {accepted} | {vs_fresh} | {vs_expected} |")
        if accepted != "yes" or vs_fresh != "yes":
            bad += 1
    print(f"\n{len(dirs) - bad}/{len(dirs)} systems: accepted, and identical to a fresh run of the original")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
