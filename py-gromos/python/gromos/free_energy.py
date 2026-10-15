"""
Free-energy windows from Python: one perturbed :class:`~gromos.Recipe` becomes one recipe per λ
(:class:`LambdaSchedule`), :func:`run_windows` runs them in process and keeps each window's
energies with their dH/dλ column, and :class:`FreeEnergyWindows` integrates ⟨∂H/∂λ⟩ over λ
(thermodynamic integration, the estimator of ``gromos.analysis.ti_trapezoid``).

>>> from gromos import Recipe, System
>>> from gromos.free_energy import LambdaSchedule, run_windows
>>> system = System.from_files("ch4_spc.top", "ch4_spc.cnf")
>>> base = Recipe.from_imd("ch4_water_fep.in").with_inputs(pttopo="ch4_spc_dummy.ptp")
>>> recipes = LambdaSchedule([0.0, 0.25, 0.5, 0.75, 1.0]).recipes(base)
>>> windows = run_windows(system, recipes, workdir="ti")
>>> dg, error = windows.ti()

A window samples dH/dλ every ``outputs["ntwg"]`` steps (every step when 0) at the steps a
``.trg`` of the same run holds: 0, NTWG, 2·NTWG, … below ``control["steps"]``. A window whose
file already exists in ``workdir`` is read, not run again.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from .exceptions import RecipeError
from .gromos import Recipe, Simulation, System, _ti_trapezoid_series
from .timeseries import EnergyTimeseries

if TYPE_CHECKING:
    import os
    from collections.abc import Iterator, Sequence

__all__ = ["LambdaSchedule", "Window", "FreeEnergyWindows", "run_windows", "collect"]


def window_name(lam: float) -> str:
    """``l0.250``: the name of the window at λ (its file in a ``workdir`` is ``<name>.npz``)."""
    return f"l{lam:.3f}"


@dataclass(frozen=True)
class LambdaSchedule:
    """The λ points of a set of windows and, optionally, the individual λ of each interaction.

    Parameters
    ----------
    points : sequence of float
        λ of each window, in [0, 1].
    interaction_lambdas : sequence of dict, optional
        ``LAMBDAS`` rows, as ``Recipe.perturbation["lambdas"]["rows"]`` holds them
        (``{"kind": "vdw", "g1": 1, "g2": 2, "a": 0, "b": 0, "c": 0.5, "d": 0.5, "e": 0}``):
        the same polynomials λ_i(λ) in every window.
    """

    points: Sequence[float]
    interaction_lambdas: Sequence[dict[str, Any]] | None = None

    def __post_init__(self) -> None:
        pts = [float(x) for x in self.points]
        if not pts:
            raise ValueError("a LambdaSchedule needs at least one λ point")
        bad = [x for x in pts if not 0.0 <= x <= 1.0]
        if bad:
            raise ValueError(f"λ points outside [0, 1]: {bad}")
        names = [window_name(x) for x in pts]
        if len(set(names)) != len(names):
            raise ValueError(f"λ points closer than the window names resolve (0.001): {pts}")
        object.__setattr__(self, "points", tuple(pts))

    @property
    def names(self) -> list[str]:
        """The window names, one per point (:func:`window_name`)."""
        return [window_name(x) for x in self.points]

    def recipes(self, base: Recipe) -> list[Recipe]:
        """One recipe per point: ``base`` with ``perturbation.lambda`` set (and the
        ``LAMBDAS`` rows, if given). ``base`` must be a perturbed recipe (``NTG = 1``)."""
        if not base.perturbation.get("enabled"):
            raise RecipeError("the base recipe has no perturbation (PERTURBATION NTG = 0)")
        out = []
        for lam in self.points:
            change: dict[str, Any] = {"lambda": lam}
            if self.interaction_lambdas is not None:
                change["lambdas"] = {
                    "enabled": True,
                    "rows": [dict(r) for r in self.interaction_lambdas],
                }
            out.append(base.update(perturbation=change))
        return out


@dataclass
class Window:
    """One λ window: its recipe and, once run or read, its energies (``dhdl`` column)."""

    lam: float
    recipe: Recipe | None = None
    energies: EnergyTimeseries | None = None
    path: Path | None = None

    @property
    def name(self) -> str:
        return window_name(self.lam)

    @property
    def dhdl(self) -> np.ndarray:
        """dH/dλ (kJ/mol) of every sampled step — the first value of each ``.trg`` frame."""
        if self.energies is None:
            raise ValueError(f"window {self.name} has not been run")
        return self.energies.dhdl


@dataclass
class FreeEnergyWindows:
    """The windows of one transformation, in ascending λ."""

    windows: list[Window] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.windows.sort(key=lambda w: w.lam)

    def __len__(self) -> int:
        return len(self.windows)

    def __iter__(self) -> Iterator[Window]:
        return iter(self.windows)

    def __getitem__(self, i: int) -> Window:
        return self.windows[i]

    @property
    def lambdas(self) -> list[float]:
        return [w.lam for w in self.windows]

    def mean_dhdl(self, skip: int = 0) -> list[float]:
        """⟨∂H/∂λ⟩ of each window after dropping its first ``skip`` samples."""
        return [float(np.mean(w.dhdl[skip:])) for w in self.windows]

    def ti(self, skip: int = 0) -> tuple[float, float]:
        """ΔG = ∫⟨∂H/∂λ⟩dλ by the trapezoid rule and its error (kJ/mol), from the
        block-averaging error of each window (the rmsd below 4 samples), after dropping the
        first ``skip`` samples of each window."""
        r = self._ti(skip)
        if r["dg"] is None:
            raise ValueError("TI needs at least two windows with samples")
        return float(r["dg"]), float(r["dg_error"])

    def report(self, skip: int = 0) -> str:
        """A table of ⟨∂H/∂λ⟩ ± error (n) per window and the TI result."""
        r = self._ti(skip)
        lines = ["    λ    ⟨∂H/∂λ⟩ (kJ/mol)      error      n"]
        for lam, m, e, n in zip(r["lambda"], r["mean"], r["ee"], r["n"], strict=True):
            lines.append(f"{lam:6.3f} {m:18.6f} {e:10.6f} {n:6d}")
        if r["dg"] is not None:
            lines.append(f"ΔG (TI, trapezoid) = {r['dg']:.6f} ± {r['dg_error']:.6f} kJ/mol")
        return "\n".join(lines)

    def _ti(self, skip: int) -> dict[str, Any]:
        return _ti_trapezoid_series(
            self.lambdas, [[float(x) for x in w.dhdl] for w in self.windows], skip
        )


def _sample_every(recipe: Recipe) -> int:
    return int(recipe.outputs["ntwg"]) or 1


def _save(window: Window) -> None:
    assert window.path is not None and window.energies is not None
    window.path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(window.path, energies=window.energies.array, lam=window.lam)


def _load(path: Path) -> Window:
    with np.load(path) as f:
        return Window(lam=float(f["lam"]), energies=EnergyTimeseries(f["energies"]), path=path)


def run_windows(
    system: System,
    recipes: Sequence[Recipe],
    *,
    workdir: str | os.PathLike[str] | None = None,
    dry_run: bool = False,
) -> FreeEnergyWindows:
    """Run one window per recipe in process, one after the other.

    Parameters
    ----------
    system : System
        The unperturbed system; each recipe's ``inputs["pttopo"]`` brings the perturbation (or
        ``system.topology.apply_perturbation`` did).
    recipes : sequence of Recipe
        One per window (:meth:`LambdaSchedule.recipes`); λ is ``perturbation["lambda"]``.
    workdir : path, optional
        Where each window's energies are kept (``<workdir>/<name>.npz``); a window whose file
        exists is read instead of run.
    dry_run : bool
        Run nothing: return the job list, each window with its recipe and file (and its
        energies when the file already exists).

    Returns
    -------
    FreeEnergyWindows
    """
    lams = [float(r.perturbation["lambda"]) for r in recipes]
    if len({window_name(x) for x in lams}) != len(lams):
        raise ValueError(f"two windows at the same λ: {lams}")
    out = []
    for lam, recipe in zip(lams, recipes, strict=True):
        if not recipe.perturbation.get("enabled"):
            raise RecipeError(f"window λ={lam}: the recipe has no perturbation (NTG = 0)")
        path = None if workdir is None else Path(workdir) / f"{window_name(lam)}.npz"
        if path is not None and path.exists():
            w = _load(path)
            w.recipe = recipe
            out.append(w)
            continue
        w = Window(lam=lam, recipe=recipe, path=path)
        if not dry_run:
            steps = int(recipe.control["steps"])
            sim = Simulation(system, recipe)
            # frames at 0, every, … below `steps`: the frames of the run's .trg
            w.energies = EnergyTimeseries(
                sim.run(max(steps - 1, 0), ene_freq=_sample_every(recipe))
            )
            if path is not None:
                _save(w)
        out.append(w)
    return FreeEnergyWindows(out)


def collect(workdir: str | os.PathLike[str]) -> FreeEnergyWindows:
    """The windows :func:`run_windows` kept in ``workdir``."""
    paths = sorted(Path(workdir).glob("l*.npz"))
    if not paths:
        raise FileNotFoundError(f"no windows in {workdir}")
    return FreeEnergyWindows([_load(p) for p in paths])
