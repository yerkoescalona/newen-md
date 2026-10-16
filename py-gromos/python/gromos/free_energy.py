"""
Free-energy windows from Python: one perturbed :class:`~gromos.Recipe` becomes one recipe per λ
(:class:`LambdaSchedule`), :func:`run_windows` runs them in process and keeps each window's
energies with their dH/dλ column (and, with ``PRECALCLAM``, the energies at the other λ), and
:class:`FreeEnergyWindows` gives ΔG by thermodynamic integration (the estimator of
``gromos.analysis.ti_trapezoid``), BAR and MBAR (those of ``gromos.analysis.bar`` and ``mbar``
on ``ext_ti_ana @bar_data``).

>>> from gromos import Recipe, System
>>> from gromos.free_energy import LambdaSchedule, run_windows
>>> system = System.from_files("ch4_spc.top", "ch4_spc.cnf")
>>> base = Recipe.from_imd("ch4_water_fep.in").with_inputs(pttopo="ch4_spc_dummy.ptp")
>>> recipes = LambdaSchedule([0.0, 0.25, 0.5, 0.75, 1.0]).recipes(base)
>>> windows = run_windows(system, recipes, workdir="ti")
>>> dg, error = windows.ti()

BAR and MBAR need each window's energies at the neighbouring λ: a ``PRECALCLAM`` grid that
holds every window's λ, here five points on [0, 1]:

>>> recipes = LambdaSchedule([0.0, 0.25, 0.5, 0.75, 1.0]).recipes(
...     base.update(
...         perturbation={"precalclam": {"nr_lambdas": 5, "min_lambda": 0.0, "max_lambda": 1.0}}
...     )
... )
>>> windows = run_windows(system, recipes, workdir="fep")
>>> dg, error = windows.bar()
>>> dg, error = windows.mbar()

A window samples dH/dλ every ``outputs["ntwg"]`` steps (every step when 0) at the steps a
``.trg`` of the same run holds: 0, NTWG, 2·NTWG, … below ``control["steps"]``. A window whose
file already exists in ``workdir`` is read, not run again.

The perturbation can also be made here instead of read from a ``.ptp`` file: a
:class:`Perturbation` (``make_pt_top`` on the topology and a copy with new charges, types or
masses on chosen atoms; ``pt_top`` for its end states) goes to :func:`run_windows`, which can
run the windows in several processes:

>>> plus = Perturbation.mutate(system, "1:1", charge_b=1.0)
>>> plus.net_charge_a, plus.net_charge_b
(0.0, 1.0)
>>> recipes = LambdaSchedule([0.0, 0.5, 1.0]).recipes(Recipe.from_imd("ch4_water_fep.in"))
>>> windows = run_windows(system, recipes, perturbation=plus, workdir="charge", parallel=3)
"""

from __future__ import annotations

import multiprocessing
import pickle
import re
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from .exceptions import RecipeError
from .gromos import (
    Recipe,
    Simulation,
    System,
    Topology,
    _make_pt_top_mutate,
    _precalc_bar,
    _precalc_bar_data,
    _precalc_mbar,
    _precalclam_longrange_warning,
    _pt_top,
    _ptp_net_charge,
    _ti_trapezoid_series,
)
from .timeseries import EnergyTimeseries

if TYPE_CHECKING:
    import os
    from collections.abc import Iterator, Sequence

__all__ = [
    "Perturbation",
    "LambdaSchedule",
    "Window",
    "Precalclam",
    "Estimate",
    "FreeEnergyWindows",
    "run_windows",
    "collect",
]


def window_name(lam: float) -> str:
    """``l0.250``: the name of the window at λ (its file in a ``workdir`` is ``<name>.npz``)."""
    return f"l{lam:.3f}"


def _topology_of(system: System | Topology) -> Topology:
    return system.topology if isinstance(system, System) else system


@dataclass(frozen=True)
class Perturbation:
    """A perturbation topology (``.ptp``) and the solute's net charge in its two states.

    Made by :meth:`from_ptp` (a file) or :meth:`mutate` (gromos-tools' ``make_pt_top`` on a
    topology and a copy of it with new charges, types or masses); every number here comes from
    the Rust programs, the text is what they wrote. :func:`run_windows` takes it in place of
    ``recipe.inputs["pttopo"]``.

    Attributes
    ----------
    text : str
        The ``.ptp`` file's content.
    net_charge_a, net_charge_b : float or None
        The solute's net charge (e) in state A and B, as ``Topology.apply_perturbation``
        returns them; ``None`` when no topology was given to :meth:`from_ptp`.
    path : Path or None
        The file it was read from, or last written to.
    """

    text: str
    net_charge_a: float | None = None
    net_charge_b: float | None = None
    path: Path | None = None

    @classmethod
    def from_ptp(
        cls, path: str | os.PathLike[str], topology: System | Topology | None = None
    ) -> Perturbation:
        """Read a ``.ptp`` file; with ``topology`` (or a :class:`~gromos.System`), also the net
        charges of its two states."""
        p = Path(path)
        text = p.read_text()
        qa = qb = None
        if topology is not None:
            qa, qb = _ptp_net_charge(_topology_of(topology), text)
        return cls(text, qa, qb, p)

    @classmethod
    def mutate(
        cls,
        system: System | Topology,
        atoms: str | Sequence[int],
        *,
        charge_b: float | Sequence[float] | None = None,
        iac_b: int | Sequence[int] | None = None,
        mass_b: float | Sequence[float] | None = None,
        alpha_lj: float = 1.51,
        alpha_crf: float = 0.5,
    ) -> Perturbation:
        """The perturbation from ``system``'s topology (state A) to the same topology with new
        parameters on ``atoms`` (state B), written by ``make_pt_top``.

        Parameters
        ----------
        system : System or Topology
            State A.
        atoms : str or sequence of int
            A gromos++ atom specifier (``"1:1"``, ``"1:CM"``) or 0-based atom indices; solute
            atoms only.
        charge_b, iac_b, mass_b : number or sequence, optional
            State B's charge (e), IAC (the topology file's 1-based number) and mass (u): one
            value for every atom or one per atom. Left out, state A's is kept; at least one is
            needed.
        alpha_lj, alpha_crf : float
            The soft-core parameters of every perturbed atom (``make_pt_top @softpar``).
        """

        def values(v: Any) -> list[Any] | None:
            if v is None:
                return None
            return [v] if np.ndim(v) == 0 else list(v)

        topology = _topology_of(system)
        text, _ = _make_pt_top_mutate(
            topology,
            atoms if isinstance(atoms, str) else [int(i) for i in atoms],
            charge_b=values(charge_b),
            iac_b=values(iac_b),
            mass_b=values(mass_b),
            alpha_lj=alpha_lj,
            alpha_crf=alpha_crf,
        )
        qa, qb = _ptp_net_charge(topology, text)
        return cls(text, qa, qb)

    @property
    def delta_charge(self) -> float | None:
        """ΔQ = Q_B − Q_A (e); a change under a reaction field needs a finite-size correction
        (``gromos.analysis.dgslv_pbsolv``)."""
        if self.net_charge_a is None or self.net_charge_b is None:
            return None
        return self.net_charge_b - self.net_charge_a

    def write(self, path: str | os.PathLike[str]) -> Perturbation:
        """Write the ``.ptp`` file; returns the perturbation with its :attr:`path`."""
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(self.text)
        return Perturbation(self.text, self.net_charge_a, self.net_charge_b, p)

    def end_state(self, system: System | Topology, state: str = "B") -> Topology:
        """``system``'s topology in end state ``"A"`` (unchanged) or ``"B"`` (``pt_top``: the
        perturbed atoms' state-B types, masses and charges, bonded types and pairs), e.g. for
        ``gromos.analysis.dgslv_pbsolv`` on that state."""
        topology = _topology_of(system)
        if state.upper() == "A":
            return topology
        if state.upper() != "B":
            raise ValueError(f"state must be 'A' or 'B', not {state!r}")
        b, _ = _pt_top(topology, self.text)
        return b

    def __str__(self) -> str:
        return self.text


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
class Precalclam:
    """A window's extended-TI data: :attr:`gromos.Simulation.precalclam` of every sampled step,
    stacked, and what BAR and MBAR need to know of the run.

    Attributes
    ----------
    grid : tuple of (int, float, float)
        ``PRECALCLAM``'s (NRLAM, MINLAM, MAXLAM).
    nlam : float
        The run's NLAM (the λ exponent).
    temperature : float or None
        The first bath's temperature (K); ``None`` without a thermostat.
    individual_lambdas : bool
        Whether the run had ``LAMBDAS`` rows, which BAR and MBAR here do not take.
    energies, derivatives : ndarray
        (frames, NRLAM, 12): the ``.tre``'s and the ``.trg``'s ``# precalclam`` rows.
    abdih, free_abdih : ndarray
        (frames, 2): the ``# ABdih`` rows (A, B).
    """

    grid: tuple[int, float, float]
    nlam: float
    temperature: float | None
    individual_lambdas: bool
    energies: np.ndarray
    derivatives: np.ndarray
    abdih: np.ndarray
    free_abdih: np.ndarray


@dataclass
class Window:
    """One λ window: its recipe and, once run or read, its energies (``dhdl`` column) and,
    with ``PRECALCLAM``, its :class:`Precalclam` data."""

    lam: float
    recipe: Recipe | None = None
    energies: EnergyTimeseries | None = None
    path: Path | None = None
    precalclam: Precalclam | None = None

    @property
    def name(self) -> str:
        return window_name(self.lam)

    @property
    def dhdl(self) -> np.ndarray:
        """dH/dλ (kJ/mol) of every sampled step — the first value of each ``.trg`` frame."""
        if self.energies is None:
            raise ValueError(f"window {self.name} has not been run")
        return self.energies.dhdl


@dataclass(frozen=True)
class Estimate:
    """ΔG over the windows by one estimator (kJ/mol): the total, from the first window's λ to
    the last, and each neighbouring pair's. Unpacks as ``dg, error = estimate``.

    Attributes
    ----------
    method : str
        ``"BAR"`` or ``"MBAR"``.
    lambdas : tuple of float
        The windows' λ, ascending.
    dg, error : float
        G(λ_last) − G(λ_first) and its error.
    steps : tuple of (float, float)
        (ΔG, error) from ``lambdas[k]`` to ``lambdas[k + 1]``.
    """

    method: str
    lambdas: tuple[float, ...]
    dg: float
    error: float
    steps: tuple[tuple[float, float], ...]

    def __iter__(self) -> Iterator[float]:
        return iter((self.dg, self.error))


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

    def window_statistics(self, skip: int = 0) -> list[dict[str, Any]]:
        """Per window, after dropping its first ``skip`` samples: ``lambda``, ``mean``
        (⟨∂H/∂λ⟩, kJ/mol), ``ee`` (its block-averaging error, as :meth:`ti` uses it) and ``n``
        (the samples)."""
        r = self._ti(skip)
        return [
            {"lambda": lam, "mean": m, "ee": e, "n": n}
            for lam, m, e, n in zip(r["lambda"], r["mean"], r["ee"], r["n"], strict=True)
        ]

    def report(self, skip: int = 0) -> str:
        """A table of ⟨∂H/∂λ⟩ ± error (n) per window and the TI result."""
        r = self._ti(skip)
        lines = ["    λ    ⟨∂H/∂λ⟩ (kJ/mol)      error      n"]
        for w in self.window_statistics(skip):
            lines.append(f"{w['lambda']:6.3f} {w['mean']:18.6f} {w['ee']:10.6f} {w['n']:6d}")
        if r["dg"] is not None:
            lines.append(f"ΔG (TI, trapezoid) = {r['dg']:.6f} ± {r['dg_error']:.6f} kJ/mol")
        return "\n".join(lines)

    def bar_data(self, skip: int = 0, temperature: float | None = None) -> list[dict[str, Any]]:
        """The per-frame energies BAR and MBAR use, after dropping the first ``skip`` frames of
        each window: per window ``slam`` (its λ), ``plam`` (every window's λ), ``e_s`` (the
        energy at ``slam``) and ``e_p`` (at each ``plam``), kJ/mol — what ``ext_ti_ana
        @bar_data`` writes for the same frames."""
        return _precalc_bar_data(*self._precalc(temperature), skip)

    def bar(
        self,
        skip: int = 0,
        temperature: float | None = None,
        *,
        max_iterations: int = 500,
        convergence: float = 1e-10,
    ) -> Estimate:
        """ΔG between neighbouring windows by Bennett's acceptance ratio and their sum
        (``gromos.analysis.bar``: the self-consistent iteration, the analytical error, errors
        summed in quadrature), after dropping the first ``skip`` frames of each window.

        Every window must have been run with ``PRECALCLAM`` on a grid holding every window's λ.
        A ``UserWarning`` (once) when a window's recipe has the twin-range schedule that leaves
        the long-range terms out of some rows (see :func:`run_windows`).

        Parameters
        ----------
        skip : int
            Frames to drop at the start of each window (as :meth:`ti`).
        temperature : float, optional
            K; default the windows' first bath.
        max_iterations, convergence : int, float
            The iteration's limit and the relative change of ΔG that ends it.
        """
        r = _precalc_bar(
            *self._precalc(temperature),
            skip,
            max_iterations=max_iterations,
            convergence=convergence,
        )
        return Estimate(
            "BAR",
            tuple(sorted(self.lambdas)),
            float(r["dg_total"]),
            float(r["error_total"]),
            tuple((float(g), float(e)) for g, e in zip(r["dg"], r["error"], strict=True)),
        )

    def mbar(
        self,
        skip: int = 0,
        temperature: float | None = None,
        *,
        max_iterations: int = 10000,
        tolerance: float = 1e-10,
    ) -> Estimate:
        """ΔG from every window's frames at every window's λ by MBAR
        (``gromos.analysis.mbar``: asymptotic errors), after dropping the first ``skip`` frames
        of each window. The same requirements and parameters as :meth:`bar`; ``tolerance``
        ends the self-consistent iteration (largest change of a reduced free energy)."""
        r = _precalc_mbar(
            *self._precalc(temperature),
            skip,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )
        dg, ddg = r["dg"], r["ddg"]
        k = len(r["lambdas"])
        return Estimate(
            "MBAR",
            tuple(float(x) for x in r["lambdas"]),
            float(dg[0, k - 1]),
            float(ddg[0, k - 1]),
            tuple((float(dg[i, i + 1]), float(ddg[i, i + 1])) for i in range(k - 1)),
        )

    def _precalc(
        self, temperature: float | None
    ) -> tuple[list[dict[str, Any]], tuple[int, float, float], float, float]:
        """The arguments of the ``_precalc_*`` estimators, checked."""
        for text in _longrange_warnings([w.recipe for w in self.windows]):
            warnings.warn(text, UserWarning, stacklevel=3)
        if len(self.windows) < 2:
            raise ValueError("BAR and MBAR need at least two windows")
        missing = [w.name for w in self.windows if w.precalclam is None]
        if missing:
            raise ValueError(
                f"windows without PRECALCLAM energies: {missing}; run them with "
                "perturbation['precalclam'] on a grid holding every window's λ"
            )
        pcs = [w.precalclam for w in self.windows if w.precalclam is not None]
        first = pcs[0]
        if any(p.grid != first.grid or p.nlam != first.nlam for p in pcs):
            raise ValueError("the windows' PRECALCLAM grids or NLAM differ")
        if any(p.individual_lambdas for p in pcs):
            raise ValueError("BAR and MBAR of windows with LAMBDAS rows are not supported")
        if temperature is None:
            temps = {p.temperature for p in pcs}
            if len(temps) != 1 or first.temperature is None:
                raise ValueError(
                    f"the windows' temperatures are {sorted(temps, key=str)}: give temperature="
                )
            temperature = first.temperature
        windows = [
            {
                "lam": w.lam,
                "energies": p.energies,
                "derivatives": p.derivatives,
                "abdih": p.abdih,
                "free_abdih": p.free_abdih,
            }
            for w, p in zip(self.windows, pcs, strict=True)
        ]
        return windows, first.grid, first.nlam, float(temperature)

    def _ti(self, skip: int) -> dict[str, Any]:
        return _ti_trapezoid_series(
            self.lambdas, [[float(x) for x in w.dhdl] for w in self.windows], skip
        )


def _longrange_warnings(recipes: Sequence[Recipe | None]) -> list[str]:
    """The distinct PRECALCLAM long-range warnings of these recipes (``None``: not known)."""
    out: list[str] = []
    for r in recipes:
        text = None if r is None else _precalclam_longrange_warning(r)
        if text is not None and text not in out:
            out.append(text)
    return out


def _sample_every(recipe: Recipe) -> int:
    return int(recipe.outputs["ntwg"]) or 1


def _precalclam(recipe: Recipe, frames: list[dict[str, np.ndarray]]) -> Precalclam:
    """The ``Simulation.precalclam`` of every sampled step, stacked, with the run's grid, NLAM,
    temperature and whether it had ``LAMBDAS`` rows."""
    p = recipe.perturbation
    pc = p["precalclam"]
    thermostat = recipe.ensemble.get("thermostat")
    lambdas = p.get("lambdas")
    return Precalclam(
        grid=(int(pc["nr_lambdas"]), float(pc["min_lambda"]), float(pc["max_lambda"])),
        nlam=float(p["nlam"]),
        temperature=float(thermostat["baths"][0]["temperature"]) if thermostat else None,
        individual_lambdas=bool(lambdas and lambdas.get("enabled")),
        **{
            key: np.array([f[key] for f in frames])
            for key in ("energies", "derivatives", "abdih", "free_abdih")
        },
    )


def _save(window: Window) -> None:
    assert window.path is not None and window.energies is not None
    window.path.parent.mkdir(parents=True, exist_ok=True)
    extra: dict[str, Any] = {}
    pc = window.precalclam
    if pc is not None:
        extra = {
            "pc_grid": np.array(pc.grid[1:]),
            "pc_nr_lambdas": pc.grid[0],
            "pc_nlam": pc.nlam,
            "pc_temperature": np.nan if pc.temperature is None else pc.temperature,
            "pc_individual_lambdas": pc.individual_lambdas,
            "pc_energies": pc.energies,
            "pc_derivatives": pc.derivatives,
            "pc_abdih": pc.abdih,
            "pc_free_abdih": pc.free_abdih,
        }
    np.savez(window.path, energies=window.energies.array, lam=window.lam, **extra)


def _load(path: Path) -> Window:
    with np.load(path) as f:
        pc = None
        if "pc_grid" in f:
            t = float(f["pc_temperature"])
            pc = Precalclam(
                grid=(int(f["pc_nr_lambdas"]), float(f["pc_grid"][0]), float(f["pc_grid"][1])),
                nlam=float(f["pc_nlam"]),
                temperature=None if np.isnan(t) else t,
                individual_lambdas=bool(f["pc_individual_lambdas"]),
                energies=f["pc_energies"],
                derivatives=f["pc_derivatives"],
                abdih=f["pc_abdih"],
                free_abdih=f["pc_free_abdih"],
            )
        return Window(
            lam=float(f["lam"]),
            energies=EnergyTimeseries(f["energies"]),
            path=path,
            precalclam=pc,
        )


def _run_window(system: System, recipe: Recipe, lam: float, path: Path | None) -> Window:
    """Run one window in this process (also the worker of ``run_windows(parallel=N)``)."""
    w = Window(lam=lam, recipe=recipe, path=path)
    steps = int(recipe.control["steps"])
    every = _sample_every(recipe)
    with warnings.catch_warnings():
        # run_windows has warned once for all windows
        for text in _longrange_warnings([recipe]):
            warnings.filterwarnings("ignore", message=re.escape(text), category=UserWarning)
        sim = Simulation(system, recipe)
    # frames at 0, every, … below `steps`: the frames of the run's .trg
    rows = [sim.run(0)[0]]
    pre = [sim.precalclam]
    for _ in range(max(steps - 1, 0) // every):
        rows.append(sim.run(every, ene_freq=every)[-1])
        pre.append(sim.precalclam)
    w.energies = EnergyTimeseries(np.array(rows))
    if pre[0] is not None:
        w.precalclam = _precalclam(recipe, [p for p in pre if p is not None])
    if path is not None:
        _save(w)
    return w


def run_windows(
    system: System,
    recipes: Sequence[Recipe],
    *,
    perturbation: Perturbation | None = None,
    workdir: str | os.PathLike[str] | None = None,
    dry_run: bool = False,
    parallel: int = 1,
) -> FreeEnergyWindows:
    """Run one window per recipe in process, one after the other or ``parallel`` at a time.

    Parameters
    ----------
    system : System
        The unperturbed system; the perturbation comes from ``perturbation``, from each
        recipe's ``inputs["pttopo"]`` or from the system's topology
        (``Topology.apply_perturbation`` before ``System(topology, configuration)``).
    recipes : sequence of Recipe
        One per window (:meth:`LambdaSchedule.recipes`); λ is ``perturbation["lambda"]``.
    perturbation : Perturbation, optional
        Written to ``<workdir>/perturbation.ptp`` and given to every window as its
        ``inputs["pttopo"]``; without ``workdir`` the file is temporary and gone once the
        windows have run. A recipe that already names a perturbation topology is refused.
    workdir : path, optional
        Where each window's energies (and :class:`Precalclam` data) are kept
        (``<workdir>/<name>.npz``); a window whose file exists is read instead of run.
    dry_run : bool
        Run nothing: return the job list, each window with its recipe and file (and its
        energies when the file already exists).
    parallel : int
        Windows run at the same time, each in its own process (``multiprocessing``, *spawn*):
        a ``Simulation`` holds the interpreter and is bound to the thread that made it, so
        threads would run one after the other. A worker builds the system again with
        ``System.from_files`` on the same files (the only form a ``System`` pickles in), runs
        its window exactly as the sequential loop does and returns it: the numbers are the
        same. Each worker runs the engine with the whole thread pool, as a sequential run does
        (the parallel kernels' sums depend on the number of threads); more than one process
        pays off when one window does not keep the machine busy (small systems, or recipes
        with ``execution.parallel = "serial"``). A script that uses it needs the
        ``if __name__ == "__main__":`` guard *spawn* requires.

    Returns
    -------
    FreeEnergyWindows

    Warns
    -----
    UserWarning
        Once, when a recipe has ``PRECALCLAM`` with a twin-range pairlist whose update (NSNB)
        and output (NTWG) steps interleave: as in gromosXX, some λ_k rows miss the long-range
        terms and BAR/MBAR from them is biased. NTWG a multiple of NSNB avoids it.
    """
    if parallel < 1:
        raise ValueError(f"parallel must be at least 1, not {parallel}")
    lams = [float(r.perturbation["lambda"]) for r in recipes]
    if len({window_name(x) for x in lams}) != len(lams):
        raise ValueError(f"two windows at the same λ: {lams}")
    for lam, recipe in zip(lams, recipes, strict=True):
        if not recipe.perturbation.get("enabled"):
            raise RecipeError(f"window λ={lam}: the recipe has no perturbation (NTG = 0)")
        if perturbation is not None and recipe.inputs.get("pttopo"):
            raise RecipeError(
                f"window λ={lam}: the recipe names a perturbation topology "
                f"({recipe.inputs['pttopo']}) and perturbation= gives another"
            )
    for text in _longrange_warnings(recipes):
        warnings.warn(text, UserWarning, stacklevel=2)
    with tempfile.TemporaryDirectory(prefix="gromos_fe_") as tmp:
        if perturbation is not None:
            where = Path(workdir) if workdir is not None else Path(tmp)
            ptp = perturbation.write(where / "perturbation.ptp").path
            recipes = [r.with_inputs(pttopo=str(ptp)) for r in recipes]
        out: list[Window] = []
        todo: list[tuple[int, Recipe, float, Path | None]] = []
        for lam, recipe in zip(lams, recipes, strict=True):
            path = None if workdir is None else Path(workdir) / f"{window_name(lam)}.npz"
            if path is not None and path.exists():
                w = _load(path)
                w.recipe = recipe
                out.append(w)
                continue
            out.append(Window(lam=lam, recipe=recipe, path=path))
            if not dry_run:
                todo.append((len(out) - 1, recipe, lam, path))
        if parallel == 1 or len(todo) < 2:
            for i, recipe, lam, path in todo:
                out[i] = _run_window(system, recipe, lam, path)
        else:
            try:
                pickle.dumps(system)
            except (TypeError, pickle.PicklingError) as e:
                raise TypeError(
                    f"run_windows(parallel={parallel}) sends the system to worker processes: {e}"
                ) from e
            ctx = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(max_workers=min(parallel, len(todo)), mp_context=ctx) as pool:
                jobs = {
                    i: pool.submit(_run_window, system, recipe, lam, path)
                    for i, recipe, lam, path in todo
                }
                for i, job in jobs.items():
                    out[i] = job.result()
    return FreeEnergyWindows(out)


def collect(workdir: str | os.PathLike[str]) -> FreeEnergyWindows:
    """The windows :func:`run_windows` kept in ``workdir``."""
    paths = sorted(Path(workdir).glob("l*.npz"))
    if not paths:
        raise FileNotFoundError(f"no windows in {workdir}")
    return FreeEnergyWindows([_load(p) for p in paths])
