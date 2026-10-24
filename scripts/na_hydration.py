#!/usr/bin/env python3
"""
Na+ -> neutral dummy in SPC water under a reaction field: the raw ΔG by TI, BAR and MBAR, and
the dGslv_pbsolv finite-size terms of the two end states (`docs/free_energy_validation.md`).

Inputs in `docs/free_energy_validation/`: `na_em.imd`, `na_eq.imd`, `na_window.imd` (written
for the 3.0 nm box; the atom and molecule counts are set here for other boxes). Every step
writes into `--workdir` and reads what an earlier step left there:

    build     the topology of one Na+ in SPC (red_top: atom 1 of `nacl_water_box.topo`, the
              FE-9 reference system), the ion solvated in a cubic box of edge --box (sim_box,
              minwall = box/2, the solvent box --solvent), and the FE-9 perturbation
              `shared/nacl_na_dummy.ptp` (Na+ IAC 34, +1 -> DUM IAC 19, 0; mass kept; α 0.5/0.5)
    prepare   steepest descent and the state-A equilibration (`md`, --md)
    windows   the λ windows of --replicas, --workers at a time, each one through
              newen.free_energy.run_windows in its own process (one engine thread each)
    analyse   TI, BAR and MBAR per replica after --skip-ps, mean ± std over replicas
    pb        newen.analysis.dgslv_pbsolv on both end states of the equilibrated box
    correct   the chain per estimator: [raw − rf_self_energy] + PB(B) − PB(A) + C1 + C2, and
              ΔG°_slv(Na⁺) = −that + ΔG_cmp (`correction.json`; needs `analysis` and `pb`)

    cd py-newen && uv run python ../scripts/na_hydration.py build --box 3.0 --workdir W
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np

from newen import Recipe, System, Topology, analysis, tools
from newen.free_energy import LambdaSchedule, Perturbation, collect, run_windows

REPO = Path(__file__).resolve().parent.parent
INPUTS = REPO / "docs/free_energy_validation"
REF = REPO / "crates/newen-md/tests/gromosXX_references"
TOPO = REF / "nacl_water_box/nacl_water_box.topo"
PTP = REF / "shared/nacl_na_dummy.ptp"
SOLVENT = REPO / ".local/gromos_tutorial_livecoms/tutorial_files/t_01/box/spc.cnf"
LAMBDAS = [i / 10 for i in range(11)]
SEED = 210185


def counts(workdir: Path) -> tuple[int, int]:
    """(atoms, SPC molecules) of the solvated box."""
    meta = json.loads((workdir / "build.json").read_text())
    return int(meta["n_atoms"]), int(meta["n_solvent"])


def for_box(recipe: Recipe, n_atoms: int, n_solvent: int) -> Recipe:
    """A template recipe with this box's molecule and atom counts (one bath, groups Na / all)."""
    r = recipe.update(
        system={"nsm": n_solvent},
        forcefield={"energy_groups": {"count": 2, "last_atoms": [1, n_atoms]}},
    )
    if r.ensemble.get("thermostat"):
        r = r.update(ensemble={"thermostat": {"dof_sets": [[n_atoms, 1, 1]]}})
    return r


def cmd_build(a: argparse.Namespace) -> None:
    w = a.workdir
    w.mkdir(parents=True, exist_ok=True)
    top = tools.red_top(TOPO, last=1)
    top.write(str(w / "na_spc.top"))
    (w / "na.cnf").write_text(
        "TITLE\none Na+ at the origin\nEND\nPOSITION\n"
        "    1 NA+   NA         1    0.000000000    0.000000000    0.000000000\nEND\n"
    )
    cnf = tools.sim_box(top, w / "na.cnf", solvent=a.solvent, minwall=a.box / 2)
    (w / "na_spc.cnf").write_text(cnf.to_g96())
    pert = Perturbation.from_ptp(PTP, top)
    pert.write(w / "na_dummy.ptp")
    meta = {
        "box": cnf.box_dimensions,
        "n_atoms": cnf.n_atoms,
        "n_solvent": (cnf.n_atoms - 1) // 3,
        "net_charge_a": pert.net_charge_a,
        "net_charge_b": pert.net_charge_b,
        "solvent_box": str(a.solvent),
    }
    (w / "build.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta))


def cmd_prepare(a: argparse.Namespace) -> None:
    w = a.workdir
    n_atoms, n_solvent = counts(w)
    top = w / "na_spc.top"
    for stage, conf in (("em", w / "na_spc.cnf"), ("eq", w / "em.cnf")):
        recipe = for_box(Recipe.from_imd(str(INPUTS / f"na_{stage}.imd")), n_atoms, n_solvent)
        imd = w / f"{stage}.imd"
        recipe.save_imd(str(imd))
        cmd = [str(a.md), "@topo", str(top), "@conf", str(conf), "@input", str(imd)]
        cmd += ["@fin", str(w / f"{stage}.cnf"), "@tre", str(w / f"{stage}.tre")]
        if stage == "eq":
            cmd += ["@trc", str(w / "eq.trc")]
        t0 = time.perf_counter()
        with open(w / f"{stage}.log", "w") as log:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=w).returncode
        print(f"{stage}: rc {rc}, {time.perf_counter() - t0:.1f} s", flush=True)
        if rc != 0:
            raise SystemExit(f"{stage} failed: {w / f'{stage}.log'}")


def window_recipes(a: argparse.Namespace, rep: int) -> list[Recipe]:
    n_atoms, n_solvent = counts(a.workdir)
    base = for_box(Recipe.from_imd(str(INPUTS / "na_window.imd")), n_atoms, n_solvent)
    base = base.update(control={"steps": a.steps}).with_inputs(
        pttopo=str(a.workdir / "na_dummy.ptp")
    )
    recipes = LambdaSchedule(LAMBDAS).recipes(base)
    return [r.update(control={"seed": SEED + 1000 * rep + i}) for i, r in enumerate(recipes)]


def _window(top: str, cnf: str, recipe: dict[str, Any], out: str) -> tuple[float, float]:
    """One window in this worker process: (λ, wall seconds)."""
    r = Recipe.from_dict(recipe)
    t0 = time.perf_counter()
    run_windows(System.from_files(top, cnf), [r], workdir=out)
    return float(r.perturbation["lambda"]), time.perf_counter() - t0


def cmd_windows(a: argparse.Namespace) -> None:
    os.environ.setdefault("RAYON_NUM_THREADS", "1")
    w = a.workdir
    jobs = []
    for rep in a.replicas:
        out = w / f"rep{rep}"
        for r in window_recipes(a, rep):
            lam = float(r.perturbation["lambda"])
            if not (out / f"l{lam:.3f}.npz").exists():
                jobs.append((rep, lam, str(w / "na_spc.top"), str(w / "eq.cnf"), r.to_dict(), out))
    # longest first is not known; replicas interleaved so that each finishes early
    jobs.sort(key=lambda j: (j[1], j[0]))
    ctx = multiprocessing.get_context("spawn")
    t0 = time.perf_counter()
    with (
        ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx) as pool,
        open(w / "windows.jsonl", "a") as log,
    ):
        futures = {pool.submit(_window, *j[2:5], str(j[5])): j[:2] for j in jobs}
        for f in as_completed(futures):
            rep, lam = futures[f]
            _, wall = f.result()
            rec = {
                "rep": rep,
                "lambda": lam,
                "steps": a.steps,
                "wall_s": wall,
                "workers": a.workers,
            }
            log.write(json.dumps(rec) + "\n")
            log.flush()
            print(f"rep{rep} λ={lam:.1f} {wall:8.0f} s", flush=True)
    print(f"{len(jobs)} windows, {time.perf_counter() - t0:.0f} s", flush=True)


def cmd_analyse(a: argparse.Namespace) -> None:
    w = a.workdir
    reps = sorted(p for p in w.glob("rep*") if p.is_dir())
    out: dict[str, Any] = {"replicas": {}}
    for rep in reps:
        fw = collect(rep)
        energies = [x.energies for x in fw if x.energies is not None]
        if len(energies) != len(fw):
            raise SystemExit(f"{rep}: a window has no energy time series")
        ntwg = int(np.diff(energies[0].time[:2])[0] / 0.002 + 0.5)
        skip = int(round(a.skip_ps / (0.002 * ntwg)))
        ti = fw.ti(skip)
        bar, mbar = fw.bar(skip), fw.mbar(skip)
        stats = fw.window_statistics(skip)
        out["replicas"][rep.name] = {
            "skip_frames": skip,
            "frames": [len(x.dhdl) for x in fw],
            "TI": list(ti),
            "BAR": [bar.dg, bar.error],
            "MBAR": [mbar.dg, mbar.error],
            "windows": stats,
            "volume_nm3": [float(np.mean(e.volume[skip:])) for e in energies],
        }
    summary: dict[str, Any] = {}
    for est in ("TI", "BAR", "MBAR"):
        v = np.array([r[est] for r in out["replicas"].values()])
        summary[est] = {
            "n": len(v),
            "mean": float(v[:, 0].mean()),
            "std": float(v[:, 0].std(ddof=1)) if len(v) > 1 else None,
            "replicas": v[:, 0].tolist(),
            "mean_error": float(v[:, 1].mean()),
        }
    per = np.array([[s["mean"] for s in r["windows"]] for r in out["replicas"].values()])
    vol = np.array([r["volume_nm3"] for r in out["replicas"].values()])
    summary["dhdl"] = {
        "lambda": LAMBDAS,
        "mean": per.mean(axis=0).tolist(),
        "std": per.std(axis=0, ddof=1).tolist() if len(per) > 1 else None,
    }
    summary["edge_nm"] = {
        "per_window_mean": np.cbrt(vol).mean(axis=0).tolist(),
        "mean": float(np.cbrt(vol).mean()),
        "std": float(np.cbrt(vol).std(ddof=1)),
    }
    out["summary"] = summary
    (w / "analysis.json").write_text(json.dumps(out, indent=2))
    for est in ("TI", "BAR", "MBAR"):
        s = summary[est]
        reps_txt = ", ".join(f"{x:.2f}" for x in s["replicas"])
        print(
            f"{est:5s} n={s['n']} {s['mean']:.2f} ± {s['std'] or 0:.2f} [{reps_txt}] err {s['mean_error']:.2f}"
        )
    print(f"edge {summary['edge_nm']['mean']:.4f} ± {summary['edge_nm']['std']:.4f} nm")


def cmd_pb(a: argparse.Namespace) -> None:
    w = a.workdir
    top = Topology.from_file(str(w / "na_spc.top"))
    pert = Perturbation.from_ptp(w / "na_dummy.ptp", top)
    states = {"A": w / "na_spc.top", "B": w / "na_spc_B.top"}
    pert.end_state(top, "B").write(str(states["B"]))
    res: dict[str, Any] = {"settings": {}, "states": {}}
    common = {
        "coord": str(w / a.coord),
        "pbc": "r",
        "atoms": "1:a",
        "atoms_to_charge": "1:a",
        "eps_solv": a.eps_solv,
        "gridspacing": 0.1,
        "npbc_size": 5.0,
        "probe_iac": a.probe_iac,
    }
    schemes = {"LS": {"scheme": "LS"}, "RF": {"scheme": "RF", "eps_rf": a.eps_rf, "rcut": a.rcut}}
    res["settings"] = {"common": common, "schemes": schemes}
    for state, topo in states.items():
        for name, scheme in schemes.items():
            r = analysis.dgslv_pbsolv(
                topo=str(topo), **common, **scheme, verbose=str(w / f"pb_{state}_{name}.log")
            )
            (w / f"pb_{state}_{name}.stdout").write_text(r["stdout"])
            keep = ("dg_npbc", "dg_pbc", "dg_fft_ls_pbc", "dg_fft_rf_pbc", "documented_correction")
            res["states"][f"{state}_{name}"] = {k: r[k] for k in keep}
    for name in schemes:
        res[f"term_{name}"] = (
            res["states"][f"B_{name}"]["documented_correction"]
            - res["states"][f"A_{name}"]["documented_correction"]
        )
    (w / "pb.json").write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != "settings"}, indent=2, default=str))


def cmd_correct(a: argparse.Namespace) -> None:
    w = a.workdir
    an = json.loads((w / "analysis.json").read_text())
    pb = json.loads((w / "pb.json").read_text())
    meta = json.loads((w / "build.json").read_text())
    qa, qb = float(meta["net_charge_a"]), float(meta["net_charge_b"])
    edge = an["summary"]["edge_nm"]["mean"]
    self_term = analysis.rf_self_energy(qa, qb, a.rcut, a.eps_rf)
    cc = analysis.charge_correction(qb - qa, edge, a.rcut, a.eps_rf, "SPC", charge_a=qa)
    pb_a = pb["states"]["A_RF"]["documented_correction"]
    pb_b = pb["states"]["B_RF"]["documented_correction"]
    out: dict[str, Any] = {
        "rf_self": self_term,
        "charge_correction": cc,
        "box_edge_nm": edge,
        "pb_a": pb_a,
        "pb_b": pb_b,
        "dg_cmp": a.dg_cmp,
        "estimators": {},
    }
    for est in ("TI", "BAR", "MBAR"):
        raws = an["summary"][est]["replicas"]
        cor = [analysis.corrected_charging_dg(r, self_term, cc, pb_a=pb_a, pb_b=pb_b) for r in raws]
        sph = [analysis.corrected_charging_dg(r, self_term, cc) for r in raws]
        slv = [-c + a.dg_cmp for c in cor]
        slv_sph = [-c + a.dg_cmp for c in sph]

        def ms(v: list[float]) -> list[float | None]:
            x = np.array(v)
            return [float(x.mean()), float(x.std(ddof=1)) if len(x) > 1 else None]

        out["estimators"][est] = {
            "corrected_pb": cor,
            "corrected_sphere": sph,
            "dg_slv_pb": slv,
            "dg_slv_sphere": slv_sph,
            "corrected_pb_mean_std": ms(cor),
            "dg_slv_pb_mean_std": ms(slv),
            "dg_slv_sphere_mean_std": ms(slv_sph),
            "mean_error": an["summary"][est]["mean_error"],
        }
    (w / "correction.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("build", "prepare", "windows", "analyse", "pb", "correct"):
        p = sub.add_parser(name)
        p.add_argument("--workdir", type=Path, required=True)
        if name == "build":
            p.add_argument("--box", type=float, default=3.0, help="cubic box edge (nm)")
            p.add_argument("--solvent", type=Path, default=SOLVENT)
        if name == "prepare":
            p.add_argument("--md", type=Path, default=REPO / "target/release/md")
        if name == "windows":
            p.add_argument("--replicas", type=lambda s: [int(x) for x in s.split(",")])
            p.add_argument("--steps", type=int, default=10_000)
            p.add_argument("--workers", type=int, default=33)
        if name == "analyse":
            p.add_argument("--skip-ps", type=float, default=4.0)
        if name == "pb":
            p.add_argument("--coord", default="eq.cnf")
            p.add_argument("--eps-solv", type=float, default=66.6)
            p.add_argument("--probe-iac", type=int, default=5)
        if name in ("pb", "correct"):
            p.add_argument("--eps-rf", type=float, default=61.0)
            p.add_argument("--rcut", type=float, default=1.4)
        if name == "correct":
            p.add_argument("--dg-cmp", type=float, default=7.95, help="standard state (kJ/mol)")
    a = ap.parse_args()
    {
        "build": cmd_build,
        "prepare": cmd_prepare,
        "windows": cmd_windows,
        "analyse": cmd_analyse,
        "pb": cmd_pb,
        "correct": cmd_correct,
    }[a.cmd](a)


if __name__ == "__main__":
    main()
