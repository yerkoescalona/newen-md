#!/usr/bin/env python3
"""
Build a large solvated benchmark system for BENCHMARKING.md Phase 0.5, using
the project's own `sim_box` (the gromos++ solvation program, ported in
newen-tools) — this script only prepares its inputs and writes the MD input.

`sim_box` tiles the *equilibrated* 1000-SPC reference box to fill a target
box, drops solvent overlapping the solute, and writes the coordinates. The
target box is an exact multiple of the solvent box so the tiling is seamless
(assumption A15: equilibrated water, no lattice).

Writes `bench/systems/<name>/{<name>.cnf,<name>.imd,input.toml}` in the layout
`scripts/bench_engines.py --ref-dir bench/systems` expects.

    scripts/make_solvated_box.py --replicate 2   # ~8 000 waters, ~24 000 atoms, 6.2 nm
    scripts/make_solvated_box.py --replicate 3   # ~27 000 waters, ~81 000 atoms, 9.3 nm

Requires `cargo build --release -p newen-tools --bin sim_box`.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REF = REPO / "crates/newen-md/tests/gromosXX_references"
WATER_CNF = REF / "water_1000_spc_gridcell/water_1000_spc.cnf"
CH4_TOP = REF / "shared/ch4_spc.top"
CH4_CNF = REF / "shared/ch4_spc.cnf"
FEP_IMD = REF / "ch4_water_fep/ch4_water_fep.in"
SIM_BOX = REPO / "target/release/sim_box"

# United-atom CH4 must not sit inside a water's LJ core: solvent molecules with
# any atom closer than this (nm) are removed by sim_box's @thresh.
CH4_CLEARANCE = 0.33


def genbox_lengths(cnf: Path) -> tuple[list[str], list[float]]:
    """GENBOX block lines of a .cnf and the box lengths it declares."""
    text = cnf.read_text()
    lines = re.search(r"^GENBOX\n(.*?)^END", text, re.M | re.S).group(1).splitlines()
    data = [i for i, l in enumerate(lines) if not l.startswith("#")]
    return lines, [float(v) for v in lines[data[1]].split()]


def block(text: str, name: str) -> str:
    m = re.search(rf"^{name}\n.*?^END\n", text, re.M | re.S)
    if not m:
        raise SystemExit(f"block {name} not found in {FEP_IMD}")
    return m.group(0).rstrip("\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replicate", "-n", type=int, default=2,
                    help="target box = n × the 1000-water box per axis (default 2)")
    ap.add_argument("--steps", type=int, default=500, help="NSTLIM for the benchmark input")
    ap.add_argument("--out", type=Path, default=REPO / "bench/systems")
    args = ap.parse_args()
    if not SIM_BOX.exists():
        raise SystemExit(f"{SIM_BOX} not found: cargo build --release -p newen-tools --bin sim_box")

    n = args.replicate
    genbox, (lx, ly, lz) = genbox_lengths(WATER_CNF)
    target = [lx * n, ly * n, lz * n]

    # Solute coordinate file: the CH4 at the origin with the target GENBOX, so
    # sim_box @boxsize tiles the solvent box exactly n times per axis.
    work = args.out / "_solute"
    work.mkdir(parents=True, exist_ok=True)
    data = [i for i, l in enumerate(genbox) if not l.startswith("#")]
    genbox = list(genbox)
    genbox[data[1]] = "".join(f"{v:15.9f}" for v in target)
    solute = work / "ch4_solute.cnf"
    solute.write_text("\n".join([
        "TITLE", f"united-atom CH4 at the origin; GENBOX = {n} x water_1000_spc box", "END",
        "POSITION", "# first 24 chars ignored",
        f"{1:5d} {'CH4':<5} {'CM':<6}{1:6d}{0.0:15.9f}{0.0:15.9f}{0.0:15.9f}", "END",
        "GENBOX", *genbox, "END",
    ]) + "\n")

    # ── the project's solvation tool does the actual work ────────────────────
    proc = subprocess.run(
        [str(SIM_BOX), "@topo", str(CH4_TOP), "@pbc", "r", "@pos", str(solute),
         "@solvent", str(WATER_CNF), "@boxsize", "@thresh", str(CH4_CLEARANCE)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise SystemExit(f"sim_box failed:\n{proc.stderr}")
    m = re.search(r"Kept (\d+) solvent molecules \((\d+) atoms\)", proc.stderr)
    if not m:
        raise SystemExit(f"could not read sim_box summary:\n{proc.stderr}")
    nsm, solvent_atoms = int(m.group(1)), int(m.group(2))
    natoms = 1 + solvent_atoms

    name = f"ch4_water_{nsm}"
    out = args.out / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.cnf").write_text(proc.stdout)

    # ── MD input: the FEP reference minus perturbation, production pairlist settings ──
    fep = FEP_IMD.read_text()
    imd = "\n".join([
        "TITLE", f"{name}: production-regime benchmark (NSNB=5, twin-range 0.8/1.4, grid_cell)", "END",
        "SYSTEM", "#      NPM      NSM", f"         1     {nsm}", "END",
        "INITIALISE", "#  NTIVEL  NTISHK  NTINHT  NTINHB", "        1       0       0       0",
        "#  NTISHI  NTIRTC  NTICOM", "        1       0       0", "#  NTISTI", "        0",
        "#      IG   TEMPI", "   210185   300.0", "END",
        "STEP", "#   NSTLIM         T        DT", f"      {args.steps}       0.0     0.002", "END",
        block(fep, "BOUNDCOND"),
        "MULTIBATH", "    weak-coupling", "#   NBATHS", "    1", "#   TEMP0  TAU", "    300.0  0.1",
        "#   DOFSET", "    1", "#   LAST   COM-BATH  IR-BATH", f"    {natoms}     1         1", "END",
        block(fep, "COMTRANSROT"),
        "PRINTOUT", "#     NTPR      NTPP", f"     {args.steps}         0", "END",
        "WRITETRAJ", "#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB",
        "         0         0         0         0         0         0         0", "END",
        block(fep, "CONSTRAINT"),
        "FORCE", "#      NTF array", "# bonds    angles    imp.     dihe     charge nonbonded",
        "# H        H         H        H", "     0        1         1        1     1  1",
        "# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)", f"     2       1     {natoms}", "END",
        "PAIRLIST", "#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE",
        "        grid_cell       5       0.8     1.4     auto    chargegroup", "END",
        block(fep, "NONBONDED"),
        block(fep, "POSITIONRES"),
    ]) + "\n"
    (out / f"{name}.imd").write_text(imd)

    topo_rel = Path("..") / ".." / ".." / CH4_TOP.relative_to(REPO)
    (out / "input.toml").write_text(
        f"# Generated by scripts/make_solvated_box.py --replicate {n} via sim_box\n"
        f"[system]\nname = \"{name}\"\n"
        f"description = \"CH4 in {nsm} SPC; sim_box tiling of water_1000_spc {n}x{n}x{n}; production regime\"\n"
        f"natoms = {natoms}\n\n[input]\ntopology = \"{topo_rel}\"\nconfiguration = \"{name}.cnf\"\n"
        f"parameters = \"{name}.imd\"\n"
    )
    print(f"{out}: {nsm} waters, {natoms} atoms, box {target[0]:.3f} nm (sim_box)")


if __name__ == "__main__":
    main()
