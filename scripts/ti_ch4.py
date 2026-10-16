#!/usr/bin/env python3
"""
Free energy of CH4 -> dummy in SPC water (the vol. 7 tutorial system), through gromosXX and
gromos-rs, from the reference inputs that are proven bit-for-bit at five λ values
(`crates/gromos-md/tests/gromosXX_references/ch4_water_fep*`).

The windows are one `gromos.free_energy.LambdaSchedule` over the reference recipe (NVT 300 K,
velocities generated with seed `--seed` + window index, dH/dλ written every `--ntwg` steps).
gromos-rs runs them in process (`gromos.free_energy.run_windows`, `--parallel` worker
processes); gromosXX runs the same recipes written out as `.imd` files (`Recipe.to_imd`). The
windows are analysed with one of three estimators (`--estimator`):

- `ti` (default): the trapezoid over ⟨∂H/∂λ⟩_λ, error from the per-window block-averaging
  errors, the first `--skip` frames of each window dropped (`FreeEnergyWindows.ti`; gromosXX's
  `.trg` files through `gromos.analysis.ti_trapezoid`, the same estimator);
- `bar`: every window also precalculates its energies on a λ grid (PRECALCLAM) and BAR combines
  neighbouring windows (`FreeEnergyWindows.bar`; gromosXX's files through
  `gromos.analysis.ext_ti_ana @bar_data` + `gromos.analysis.bar`);
- `mbar`: the same energies, all windows at once (`FreeEnergyWindows.mbar`;
  `gromos.analysis.mbar`).

gromos++'s ext_ti_ana has no frame skip, so `--skip` applies to `ti` only, for both engines.
ΔG(CH4 -> dummy) is the negative hydration free energy: experiment ΔG_hyd(CH4) ≈ +8.4 kJ/mol,
so ΔG here is expected around −8 kJ/mol for 54a7 (the tutorial's own result).

    cd py-gromos && uv run python ../scripts/ti_ch4.py --steps 5000      # 11 windows, both engines
    cd py-gromos && uv run python ../scripts/ti_ch4.py --engine gromosxx --estimator mbar

Needs the py-gromos extension (`uv run maturin develop --release`) and (for gromosXX) the native
build at `.local/gromosXX/md++/BUILD/program/md`. Re-running reads the windows already in
`--workdir` instead of running them again. Writes `<workdir>/README.md` (the table) and
`<workdir>/results.json`.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np

from gromos import Recipe, System, analysis
from gromos.free_energy import LambdaSchedule, Perturbation, run_windows

REPO = Path(__file__).resolve().parent.parent
REF = REPO / "crates/gromos-md/tests/gromosXX_references"
TEMPLATE = REF / "ch4_water_fep/ch4_water_fep.in"
TOPO = REF / "shared/ch4_spc.top"
CONF = REF / "shared/ch4_spc.cnf"
PTTOPO = REF / "shared/ch4_spc_dummy.ptp"
GROMOSXX = REPO / ".local/gromosXX/md++/BUILD/program/md"
ENGINES = ["gromosxx", "gromos-rs"]


def window_recipes(
    lambdas: list[float], steps: int, ntwe: int, ntwg: int, seed: int, nrlam: int
) -> list[Recipe]:
    """The reference recipe as a production window: NVT 300 K (τ = 0.1 ps), velocities at
    300 K, energies every `ntwe` and dH/dλ every `ntwg` steps, pairlist every 5 steps, and a
    PRECALCLAM grid of `nrlam` points on [0, 1] (none for 0); one per λ, seeds seed + index."""
    base = Recipe.from_imd(str(TEMPLATE)).update(
        control={"steps": steps, "generate_velocities": True, "initial_temperature": 300.0},
        ensemble={"thermostat": {"baths": [{"temperature": 300.0, "tau": 0.1}]}},
        outputs={
            "trajectory_every": 0,
            "forces_every": 0,
            "energy_every": ntwe,
            "ntwg": ntwg,
            "print_every": 1000,
        },
        forcefield={"pairlist": {"every": 5}},
    )
    if nrlam:
        grid = {"nr_lambdas": nrlam, "min_lambda": 0.0, "max_lambda": 1.0}
        base = base.update(perturbation={"precalclam": grid})
    recipes = LambdaSchedule(lambdas).recipes(base)
    return [r.update(control={"seed": seed + i}) for i, r in enumerate(recipes)]


def grid_size(lambdas: list[float]) -> int:
    """The PRECALCLAM NRLAM of the coarsest grid on [0, 1] that holds every window λ."""
    for n in range(2, 1002):
        step = 1.0 / (n - 1)
        if all(abs(lam / step - round(lam / step)) < 1e-6 for lam in lambdas):
            return n
    raise SystemExit("the windows are on no PRECALCLAM grid of up to 1001 points: give --nrlam")


def run_gromos_rs(
    recipes: list[Recipe], wdir: Path, estimator: str, skip: int, parallel: int
) -> dict:
    system = System.from_files(str(TOPO), str(CONF))
    t0 = time.perf_counter()
    windows = run_windows(
        system,
        recipes,
        perturbation=Perturbation.from_ptp(PTTOPO, system),
        workdir=wdir,
        parallel=parallel,
    )
    elapsed = time.perf_counter() - t0
    stats = windows.window_statistics(skip if estimator == "ti" else 0)
    if estimator == "ti":
        dg, err = windows.ti(skip)
    elif estimator == "bar":
        dg, err = windows.bar()
    else:
        dg, err = windows.mbar()
    return {"windows": stats, "dG": dg, "dG_err": err, "wall_s": elapsed}


def run_gromosxx(recipes: list[Recipe], wdir: Path, estimator: str, skip: int, nrlam: int) -> dict:
    """The recipes as `.imd` files through the native binary, analysed from its files."""
    trgs, elapsed = [], 0.0
    lambdas = [float(r.perturbation["lambda"]) for r in recipes]
    for lam, recipe in zip(lambdas, recipes, strict=True):
        out = wdir / f"L_{lam:.3f}"
        out.mkdir(parents=True, exist_ok=True)
        trg = out / "window.trg"
        trgs.append(trg)
        if trg.exists() and trg.stat().st_size > 1000:
            continue
        imd = out / "window.imd"
        recipe.save_imd(str(imd))
        cmd = [
            str(GROMOSXX),
            "@topo",
            str(TOPO),
            "@pttopo",
            str(PTTOPO),
            "@conf",
            str(CONF),
            "@input",
            str(imd),
            "@fin",
            str(out / "final.cnf"),
            "@tre",
            str(out / "window.tre"),
            "@trg",
            str(trg),
        ]
        t0 = time.perf_counter()
        with open(out / "run.log", "w") as log:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
        elapsed += time.perf_counter() - t0
        if rc != 0 or not trg.exists():
            raise SystemExit(f"gromosxx failed on λ={lam:.3f}: see {out / 'run.log'}")
        print(f"gromosxx   λ={lam:.3f}  {time.perf_counter() - t0:7.1f} s", flush=True)
    if estimator == "ti":
        r = analysis.ti_trapezoid(trgs, lambdas, skip)
        stats = [
            {"lambda": lam, "mean": m, "ee": e, "n": n}
            for lam, m, e, n in zip(r["lambda"], r["mean"], r["ee"], r["n"], strict=True)
        ]
        return {"windows": stats, "dG": r["dg"], "dG_err": r["dg_error"], "wall_s": elapsed}
    # ext_ti_ana @bar_data on every window, then bar on the neighbours or mbar on all
    preds, files = [], []
    for trg, lam in zip(trgs, lambdas, strict=True):
        r = analysis.ext_ti_ana(
            en_files=[trg.with_suffix(".tre")],
            fr_files=[trg],
            temp=300.0,
            nrlambdas=nrlam,
            minlam=0.0,
            maxlam=1.0,
            slam=lam,
            nlams=1,
            bar_data=lambdas,
            lam_precision=3,
            outdir=trg.parent,
        )
        preds.append(r)
        files.append(next(f for f in r["files"] if Path(f).name.startswith("bar_data")))
    # ⟨∂H/∂λ⟩ per window from the .trg totals, as for gromos-rs
    t = analysis.ti_trapezoid(trgs, lambdas, 0)
    stats = [
        {"lambda": lam, "mean": m, "ee": e, "n": n}
        for lam, m, e, n in zip(t["lambda"], t["mean"], t["ee"], t["n"], strict=True)
    ]
    if estimator == "bar":
        b = analysis.bar(files, 300.0, convergence=1e-10)
        return {
            "windows": stats,
            "dG": b["dg_total"],
            "dG_err": b["error_total"],
            "wall_s": elapsed,
        }
    order = sorted(range(len(lambdas)), key=lambda i: lambdas[i])
    bars = [preds[i]["predictions"][0]["bar"] for i in order]
    u_kn = []
    for k, lk in enumerate(lambdas[i] for i in order):
        row = [
            b["e_s"] if m == k else b["e_p"][int(np.argmin(np.abs(b["plam"] - lk)))]
            for m, b in enumerate(bars)
        ]
        u_kn.append(np.concatenate(row))
    m = analysis.mbar(np.array(u_kn), [len(b["e_s"]) for b in bars], 300.0)
    return {
        "windows": stats,
        "dG": float(m["dg"][0, -1]),
        "dG_err": float(m["ddg"][0, -1]),
        "wall_s": elapsed,
    }


def table(args: argparse.Namespace, results: dict, engines: list[str], lambdas: list[float]) -> str:
    energies = args.estimator != "ti"
    skip, nrlam = results["skip_frames"], results["nrlam"]
    dropped = (
        f"first {skip} frames dropped" if not energies else f"PRECALCLAM {nrlam}, no frames dropped"
    )
    two = len(engines) == 2
    lines = [
        f"# {args.estimator.upper()}: CH4 -> dummy in 999 SPC (54a7), {args.steps} steps/window (dt 2 fs), "
        f"{len(lambdas)} windows, dH/dλ every {args.ntwg} steps, {dropped}",
        "",
        "| λ | "
        + " | ".join(f"⟨∂H/∂λ⟩ {e} ± ee (n)" for e in engines)
        + (" | Δ (rs − XX) |" if two else " |"),
        "|---|" + "---|" * len(engines) + ("---|" if two else ""),
    ]
    per = {e: {round(w["lambda"], 6): w for w in results["engines"][e]["windows"]} for e in engines}
    for lam in lambdas:
        cells = []
        for e in engines:
            w = per[e].get(round(lam, 6))
            if not w:
                cells.append("—")
            elif w["ee"] is None:
                cells.append(f"{w['mean']:9.3f} ({w['n']})")
            else:
                cells.append(f"{w['mean']:9.3f} ± {w['ee']:.3f} ({w['n']})")
        row = f"| {lam:.2f} | " + " | ".join(cells)
        if two:
            a, b = per[engines[0]].get(round(lam, 6)), per[engines[1]].get(round(lam, 6))
            row += f" | {b['mean'] - a['mean']:+.3f} |" if a and b else " | — |"
        else:
            row += " |"
        lines.append(row)
    lines.append("")
    for e in engines:
        r = results["engines"][e]
        if r["dG"] is None:
            lines.append(f"- **{e}**: no ΔG (fewer than two windows with frames)")
            continue
        err = f" ± {r['dG_err']:.2f}" if r.get("dG_err") is not None else ""
        lines.append(
            f"- **{e}**: ΔG = {r['dG']:.2f}{err} kJ/mol; wall time {r['wall_s']:.0f} s for "
            f"{len(lambdas)} windows" + ("" if r["wall_s"] else " (all kept)")
        )
    if two and all(results["engines"][e]["dG"] is not None for e in engines):
        a, b = (results["engines"][e]["dG"] for e in engines)
        lines.append(f"- ΔG(gromos-rs) − ΔG(gromosXX) = {b - a:+.3f} kJ/mol")
    lines.append(
        "- Reference: ΔG_hyd(CH4) experiment ≈ +8.4 kJ/mol, i.e. ΔG(CH4→dummy) ≈ −8 kJ/mol; the "
        "tutorial's 54a7 result is of that size. Windows are independent short runs, not a converged study."
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--engine", choices=["both", *ENGINES], default="both")
    ap.add_argument("--steps", type=int, default=5000, help="MD steps per window (dt = 2 fs)")
    ap.add_argument("--lambdas", default=",".join(f"{i / 10:.1f}" for i in range(11)))
    ap.add_argument("--ntwg", type=int, default=10, help="dH/dlambda output frequency (steps)")
    ap.add_argument("--ntwe", type=int, default=100, help="energy output frequency (steps)")
    ap.add_argument(
        "--skip",
        type=int,
        default=None,
        help="dH/dlambda frames to drop per window (default: 20 %%)",
    )
    ap.add_argument("--seed", type=int, default=210185)
    ap.add_argument("--workdir", type=Path, default=REPO / "bench/work/ti_ch4")
    ap.add_argument("--estimator", choices=["ti", "bar", "mbar"], default="ti")
    ap.add_argument(
        "--nrlam",
        type=int,
        default=None,
        help="PRECALCLAM grid size for bar/mbar (default: the coarsest grid holding every window)",
    )
    ap.add_argument(
        "--parallel", type=int, default=1, help="gromos-rs windows run at once (processes)"
    )
    args = ap.parse_args()

    lambdas = [float(x) for x in args.lambdas.split(",")]
    engines = ENGINES if args.engine == "both" else [args.engine]
    if "gromosxx" in engines and not GROMOSXX.exists():
        raise SystemExit(f"gromosxx: binary not found at {GROMOSXX}")
    energies = args.estimator != "ti"
    nrlam = (args.nrlam or grid_size(lambdas)) if energies else 0
    # ext_ti_ana reads the .tre and the .trg frame by frame together
    ntwe = args.ntwg if energies else args.ntwe
    n_frames = args.steps // args.ntwg
    skip = args.skip if args.skip is not None else n_frames // 5
    recipes = window_recipes(lambdas, args.steps, ntwe, args.ntwg, args.seed, nrlam)
    kind = f"precalclam{nrlam}" if energies else "plain"

    results: dict = {
        "steps": args.steps,
        "ntwg": args.ntwg,
        "skip_frames": skip,
        "lambdas": lambdas,
        "estimator": args.estimator,
        "nrlam": nrlam,
        "engines": {},
    }
    for engine in engines:
        wdir = args.workdir / engine / kind
        if engine == "gromos-rs":
            r = run_gromos_rs(recipes, wdir, args.estimator, skip, args.parallel)
        else:
            r = run_gromosxx(recipes, wdir, args.estimator, skip, nrlam)
        print(f"{engine:10s} {len(lambdas)} windows  {r['wall_s']:7.1f} s", flush=True)
        results["engines"][engine] = r

    text = table(args, results, engines, lambdas)
    args.workdir.mkdir(parents=True, exist_ok=True)
    (args.workdir / "README.md").write_text(text)
    (args.workdir / "results.json").write_text(json.dumps(results, indent=2))
    print("\n" + text)


if __name__ == "__main__":
    main()
