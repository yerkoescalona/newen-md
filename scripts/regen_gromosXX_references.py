#!/usr/bin/env python3
"""Regenerate gromosXX reference data using the patched CODATA binary.

Run from the repo root:
    python3 scripts/regen_gromosXX_references.py             # every system
    python3 scripts/regen_gromosXX_references.py meoh_water_fep ch4_water_fep_l000   # only these
    python3 scripts/regen_gromosXX_references.py --binary PATH name   # another md binary

The script runs the gromosXX md binary (patched to use CODATA 2018 constants)
on every system in crates/gromos-md/tests/gromosXX_references/ and overwrites
the expected/ output files in place.

Requires the patched binary at:
    .local/gromosXX/md++/BUILD/program/md
unless --binary PATH names another one.
"""

import argparse
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
REF_DIR = REPO_ROOT / "crates" / "gromos-md" / "tests" / "gromosXX_references"
MD_BIN = REPO_ROOT / ".local" / "gromosXX" / "md++" / "BUILD" / "program" / "md"


def parse_toml_str(text: str, key: str) -> str | None:
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(f"{key} ="):
            parts = line.split("=", 1)[1].strip()
            if parts.startswith('"') and parts.endswith('"'):
                return parts[1:-1]
    return None


def run_system(sys_dir: Path, md_bin: Path = MD_BIN) -> bool:
    name = sys_dir.name
    toml_path = sys_dir / "input.toml"
    if not toml_path.exists():
        print(f"  SKIP {name}: no input.toml")
        return True

    toml = toml_path.read_text()

    topo = sys_dir / parse_toml_str(toml, "topology")
    conf = sys_dir / parse_toml_str(toml, "configuration")
    params = sys_dir / parse_toml_str(toml, "parameters")

    with tempfile.TemporaryDirectory(prefix=f"regen_{name}_") as tmp:
        tmp = Path(tmp)
        tre = tmp / "energies.tre"
        trf = tmp / "forces.trf"
        trc = tmp / "trajectory.trc"
        trv = tmp / "velocities.trv"
        trg = tmp / "free_energy.trg"
        fin = tmp / "final.conf"

        cmd = [
            str(md_bin),
            "@topo", str(topo),
            "@conf", str(conf),
            "@input", str(params),
            "@fin", str(fin),
            "@tre", str(tre),
            "@trf", str(trf),
            "@trc", str(trc),
            "@trv", str(trv),
            "@trg", str(trg),
        ]

        posresspec = parse_toml_str(toml, "posresspec")
        if posresspec:
            cmd += ["@posresspec", str(sys_dir / posresspec)]

        refpos = parse_toml_str(toml, "refpos")
        if refpos:
            cmd += ["@refpos", str(sys_dir / refpos)]

        distrest = parse_toml_str(toml, "distrest")
        if distrest:
            cmd += ["@distrest", str(sys_dir / distrest)]

        gamd = parse_toml_str(toml, "gamd")
        if gamd:
            cmd += ["@gamd", str(sys_dir / gamd)]

        pttopo = parse_toml_str(toml, "pttopo")
        if pttopo:
            cmd += ["@pttopo", str(sys_dir / pttopo)]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"  FAIL {name}: md exited {result.returncode}")
            print(result.stderr[-2000:])
            return False

        # Copy outputs into expected/
        expected = sys_dir / "expected"
        expected.mkdir(exist_ok=True)

        copied = []
        # gromosXX writes a title-only .trv when NTWV = 0; that is not a reference.
        if trv.exists() and "VELOCITYRED" not in trv.read_text():
            trv.unlink()

        for src, dst_name in [
            (tre, "energies.tre"),
            (trf, "forces.trf"),
            (trc, "trajectory.trc"),
            (trv, "velocities.trv"),
            (fin, "final.conf"),
        ]:
            if src.exists():
                shutil.copy2(src, expected / dst_name)
                copied.append(dst_name)

        fe_path = parse_toml_str(toml, "free_energy")
        if fe_path and trg.exists():
            dst = sys_dir / fe_path
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(trg, dst)
            copied.append(fe_path)

        print(f"  OK  {name}: {', '.join(copied)}")
        return True


def main():
    parser = argparse.ArgumentParser(description="Regenerate gromosXX reference outputs.")
    parser.add_argument("--binary", type=Path, default=MD_BIN, help=f"md binary (default {MD_BIN})")
    parser.add_argument("systems", nargs="*", help="systems to regenerate (default: all)")
    args = parser.parse_args()
    md_bin = args.binary.resolve()
    if not md_bin.exists():
        print(f"ERROR: gromosXX binary not found at {md_bin}")
        sys.exit(1)

    print(f"Using binary: {md_bin}")
    print(f"Reference dir: {REF_DIR}")
    print()

    systems = sorted(p for p in REF_DIR.iterdir() if p.is_dir() and p.name != "shared")
    # Optional: regenerate only the named systems (new references) — the committed expected/
    # files of the others are never rewritten by accident.
    if args.systems:
        wanted = set(args.systems)
        unknown = wanted - {p.name for p in systems}
        if unknown:
            print(f"ERROR: unknown system(s): {', '.join(sorted(unknown))}")
            sys.exit(1)
        systems = [p for p in systems if p.name in wanted]

    failed = []
    for sys_dir in systems:
        ok = run_system(sys_dir, md_bin)
        if not ok:
            failed.append(sys_dir.name)

    print()
    print(f"Done: {len(systems) - len(failed)}/{len(systems)} systems regenerated.")
    if failed:
        print(f"FAILED: {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
