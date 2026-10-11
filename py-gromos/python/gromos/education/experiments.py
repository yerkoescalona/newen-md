"""Nine small experiments, each a measurement from the engine set beside the theory it should match.

Every function runs a short simulation, keeps what it measured in :attr:`Result.data`, draws it
over the textbook curve with :meth:`Result.plot` and states the claim as one or more
:class:`Check` objects that a notebook or a test asserts. The numbers come from the engine
(``Simulation`` properties, :class:`gromos.EnergyTimeseries`, :func:`gromos.rdf`); this module
only arranges and compares them.

Examples
--------
>>> from gromos.education import experiments as ex
>>> r = ex.energy_flow()  # doctest: +SKIP
>>> r.plot()  # doctest: +SKIP
>>> r.check().assert_passes()  # doctest: +SKIP
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from gromos import EnergyTimeseries, Recipe, Simulation, System, rdf

from . import theory
from .systems import (
    ARGON_MASS,
    EPSILON,
    K_B,
    SIGMA,
    argon_pair,
    lj_gas,
    system_from_frame,
    with_bath,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from matplotlib.figure import Figure
    from matplotlib.typing import RcKeyType

__all__ = [
    "Check",
    "Result",
    "block_error",
    "energy_flow",
    "energy_landscape",
    "heat_bath",
    "liquid_structure",
    "maxwell_boltzmann",
    "thermostat_fingerprint",
    "time_step_scan",
    "water_structure",
]

_STYLE: dict[RcKeyType, Any] = {
    "figure.dpi": 110,
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
}
_KINETIC, _POTENTIAL, _TOTAL, _THEORY, _MEASURED = (
    "#e8590c",
    "#1c7ed6",
    "#212529",
    "#868e96",
    "#7048e8",
)
_SECOND = "#f08c00"


@dataclass
class Check:
    """One claim of an experiment: a measured number against what theory expects.

    The check passes when ``measured`` lies inside ``band`` (when given) and within
    ``rel_tol``/``abs_tol`` of ``expected`` (when either is given). At least one criterion
    must be set. A NaN measurement never passes.

    Parameters
    ----------
    name : str
        What was measured, with its unit.
    measured : float
        The value the simulation produced.
    expected : float or None
        The value theory gives, or a reference value for a ``band`` check.
    rel_tol, abs_tol : float, optional
        Allowed deviation from ``expected``: relative to ``|expected|``, absolute, or the
        larger of the two when both are set.
    band : tuple of float, optional
        Inclusive ``(low, high)`` interval for ``measured``; use ``math.inf`` for an open side.
    """

    name: str
    measured: float
    expected: float | None
    rel_tol: float | None = None
    abs_tol: float | None = None
    band: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        if self.band is None and self.rel_tol is None and self.abs_tol is None:
            raise ValueError(f"check {self.name!r} needs a band, rel_tol or abs_tol")
        if (self.rel_tol is not None or self.abs_tol is not None) and self.expected is None:
            raise ValueError(f"check {self.name!r} has a tolerance but no expected value")
        self.measured = float(self.measured)
        if self.expected is not None:
            self.expected = float(self.expected)

    @property
    def passed(self) -> bool:
        """Whether the measurement satisfies every criterion that is set."""
        if not math.isfinite(self.measured):
            return False
        if self.band is not None and not self.band[0] <= self.measured <= self.band[1]:
            return False
        if self.rel_tol is not None or self.abs_tol is not None:
            assert self.expected is not None
            allowed = max(
                (self.rel_tol or 0.0) * abs(self.expected),
                self.abs_tol or 0.0,
            )
            if abs(self.measured - self.expected) > allowed:
                return False
        return True

    def __bool__(self) -> bool:
        return self.passed

    def _criterion(self) -> str:
        parts = []
        if self.band is not None:
            parts.append(f"in [{self.band[0]:.4g}, {self.band[1]:.4g}]")
        if self.rel_tol is not None:
            parts.append(f"within {100 * self.rel_tol:.3g}%")
        if self.abs_tol is not None:
            parts.append(f"within ±{self.abs_tol:.3g}")
        return " and ".join(parts)

    def _summary(self) -> str:
        expected = "—" if self.expected is None else f"{self.expected:.4g}"
        verdict = "PASS" if self.passed else "FAIL"
        return (
            f"{self.name}: measured {self.measured:.4g}, expected {expected}, "
            f"required {self._criterion()}: {verdict}"
        )

    def assert_passes(self) -> None:
        """Raise :class:`AssertionError` with a one-line summary unless the check passes."""
        if not self.passed:
            raise AssertionError(self._summary())

    def __repr__(self) -> str:
        return f"Check({self._summary()})"


@dataclass
class Result:
    """What an experiment measured, how to draw it and what it claims.

    Attributes
    ----------
    title : str
        The claim, in words.
    checks : list of Check
        The claims, as numbers.
    data : dict of str to numpy.ndarray
        Every measured series and the theory curves drawn beside them.
    """

    title: str
    checks: list[Check]
    data: dict[str, np.ndarray]
    _draw: Callable[[Result, Any], None] = field(repr=False)
    _size: tuple[float, float] = field(default=(6.4, 3.6), repr=False)
    _panels: int = field(default=1, repr=False)
    _sharey: bool = field(default=True, repr=False)

    def plot(self) -> Figure:
        """Draw the measurement over the theory.

        Returns
        -------
        matplotlib.figure.Figure
            The figure; it is not shown or saved.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        with plt.rc_context(_STYLE):
            fig, axes = plt.subplots(
                1,
                self._panels,
                figsize=self._size,
                sharey=self._panels > 1 and self._sharey,
                squeeze=False,
            )
            # A notebook's inline backend turns interactive mode on as the first figure is
            # made; leaving the style context would turn it off and hide later figures.
            interactive = matplotlib.is_interactive()
            self._draw(self, axes[0] if self._panels > 1 else axes[0, 0])
            fig.tight_layout()
        matplotlib.interactive(interactive)
        return fig

    def check(self) -> Check | list[Check]:
        """The experiment's check, or its list of checks when there are several.

        Never raises; call :meth:`Check.assert_passes` or ``assert all(result.checks)``.
        """
        return self.checks[0] if len(self.checks) == 1 else list(self.checks)


def _arrays(**values: Any) -> dict[str, np.ndarray]:
    return {key: np.asarray(value) for key, value in values.items()}


# -- 1. where the energy goes ----------------------------------------------------------------


def energy_flow(n: int = 64, steps: int = 3000, seed: int = 7) -> Result:
    """A 2D lattice of argon melts; kinetic and potential energy trade, the total stays put.

    The total still wanders by about 2% of what the kinetic energy gains. That residual comes
    from pairs crossing the truncated Lennard-Jones cut-off, where the energy jumps by
    :math:`V(r_c)`; it does not shrink when the time step is halved.

    Parameters
    ----------
    n : int
        Number of atoms of the 2D gas (:func:`~gromos.education.lj_gas`).
    steps : int
        Number of constant-energy steps; energies are sampled every 10.
    seed : int
        Seed of the initial velocities.

    Returns
    -------
    Result
        ``data``: ``time``, ``kinetic``, ``potential``, ``total``. Check: the range of the
        total energy is under 5% of what the kinetic energy gains.
    """
    system, recipe = lj_gas(n=n, seed=seed)
    ts = EnergyTimeseries(Simulation(system, recipe).run(steps, ene_freq=10))
    gain = float(ts.kinetic[-1] - ts.kinetic[0])
    spread = float(np.ptp(ts.total))
    check = Check("total-energy range / kinetic-energy gain", spread / gain, 0.0, band=(0, 0.05))

    def draw(result: Result, ax: Any) -> None:
        d = result.data
        ax.plot(d["time"], d["kinetic"], color=_KINETIC, lw=1.2, label="kinetic")
        ax.plot(d["time"], d["potential"], color=_POTENTIAL, lw=1.2, label="potential")
        ax.plot(d["time"], d["total"], color=_TOTAL, lw=2, label="total")
        ax.set(xlabel="time (ps)", ylabel="energy (kJ/mol)", title=result.title)
        ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1, 0.36), ncol=3, fontsize=9)
        ax.text(
            0.98,
            0.57,
            f"total stays within {spread:.2f} kJ/mol while kinetic gains {gain:+.0f}\n"
            "the residual: pairs crossing the LJ cut-off; halving dt does not shrink it",
            transform=ax.transAxes,
            ha="right",
            fontsize=9,
            color=_TOTAL,
        )

    return Result(
        "Where the energy goes: a 2D lattice melts",
        [check],
        _arrays(time=ts.time, kinetic=ts.kinetic, potential=ts.potential, total=ts.total),
        draw,
    )


# -- 2. time step and accuracy ---------------------------------------------------------------


def time_step_scan(
    dts: Sequence[float] = (0.0005, 0.001, 0.002, 0.004, 0.008, 0.016, 0.032),
    duration: float = 8.0,
    distance: float = 0.35,
) -> Result:
    """Halve the time step and the energy error falls fourfold: leap-frog is second order.

    Parameters
    ----------
    dts : sequence of float
        Time steps to try (ps).
    duration : float
        Length of each run (ps).
    distance : float
        Initial separation of the argon pair (nm).

    Returns
    -------
    Result
        ``data``: ``dt``, ``energy_range`` (range of the total energy over the run, kJ/mol),
        ``slope``. Check: the log–log slope lies in [1.9, 2.1].
    """
    dt = np.asarray(dts, dtype=float)
    system, recipe = argon_pair(distance)
    spreads = []
    for step in dt:
        n = round(duration / step)
        run = Simulation(system, recipe.update(control={"dt": float(step), "steps": n}))
        spreads.append(np.ptp(EnergyTimeseries(run.run(n, ene_freq=1)).total))
    spread = np.asarray(spreads)
    slope = float(np.polyfit(np.log(dt), np.log(spread), 1)[0])
    check = Check("log–log slope of energy range against dt", slope, 2.0, band=(1.9, 2.1))

    def draw(result: Result, ax: Any) -> None:
        d = result.data
        ax.loglog(d["dt"], d["energy_range"], "o-", color=_MEASURED, label="measured")
        ax.loglog(
            d["dt"],
            d["energy_range"][0] * (d["dt"] / d["dt"][0]) ** 2,
            "--",
            color=_THEORY,
            label="theory: slope 2",
        )
        ax.set(
            xlabel="time step dt (ps)",
            ylabel=f"total-energy range, {duration:g} ps (kJ/mol)",
            title=result.title,
        )
        ax.text(0.05, 0.85, f"fitted slope {float(d['slope']):.2f}", transform=ax.transAxes)
        ax.legend(frameon=False, loc="lower right")

    return Result(
        "Time step and accuracy: error grows as dt²",
        [check],
        _arrays(dt=dt, energy_range=spread, slope=slope),
        draw,
        (5.2, 3.6),
    )


# -- 3. the energy landscape of two atoms ----------------------------------------------------


def energy_landscape(distance: float = 0.35, dt: float = 0.001, steps: int = 4000) -> Result:
    """Two argon atoms swing between turning points on the Lennard-Jones curve.

    The curve is drawn from the topology's own C6 and C12 (:meth:`gromos.Topology.lj_parameters`).

    Parameters
    ----------
    distance : float
        Initial separation (nm); the pair starts at rest.
    dt : float
        Time step (ps).
    steps : int
        Number of steps.

    Returns
    -------
    Result
        ``data``: ``r`` and ``potential`` along the run, ``r_curve``/``v_curve`` (the
        potential), ``total`` (the energy of the start) and ``turning_points`` (the shortest
        and longest distance reached). Check: the lowest potential energy reached is within
        0.1% of :math:`-\\varepsilon`.
    """
    system, recipe = argon_pair(distance)
    topology = system.topology
    c6, c12, _, _ = topology.lj_parameters(int(topology.iac[0]), int(topology.iac[1]))
    r_min, epsilon = theory.lj_minimum(c6, c12)
    sim = Simulation(system, recipe.update(control={"dt": dt, "steps": steps}))
    distances, energies = [], []
    for _ in range(steps):
        p = sim.positions
        distances.append(np.linalg.norm(p[1] - p[0]))
        sim.step(1)
        energies.append(sim.potential_energy)  # evaluated at the positions the step began from
    r = np.asarray(distances)
    potential = np.asarray(energies)
    r_curve = np.linspace(0.95 * min(r.min(), r_min), 1.3 * max(r.max(), r_min), 400)
    total = float(theory.lj_energy(distance, c6, c12))
    check = Check("lowest potential energy (kJ/mol)", potential.min(), -epsilon, rel_tol=1e-3)

    def draw(result: Result, ax: Any) -> None:
        d = result.data
        ax.plot(d["r_curve"], d["v_curve"], color=_POTENTIAL, lw=2, label="V(r) from the topology")
        ax.axhline(float(d["total"]), ls="--", color=_THEORY, label="total energy")
        every = max(1, len(d["r"]) // 100)
        ax.plot(
            d["r"][::every],
            d["potential"][::every],
            "o",
            ms=3,
            color=_MEASURED,
            label=f"the run, every {every} steps",
        )
        ax.axvspan(*d["turning_points"], color=_MEASURED, alpha=0.08)
        ax.annotate(
            f"−ε = {-epsilon:.3f} kJ/mol (lowest seen {d['potential'].min():.3f})",
            xy=(r_min, -epsilon),
            xytext=(r_min * 1.12, -0.95 * epsilon),
            arrowprops={"arrowstyle": "->"},
            fontsize=8,
        )
        ax.set(
            ylim=(-1.1 * epsilon, 0.4 * epsilon),
            xlabel="distance r (nm)",
            ylabel="energy (kJ/mol)",
            title=result.title,
        )
        ax.legend(frameon=False, fontsize=8, loc="upper right")

    return Result(
        "Two atoms swing between turning points",
        [check],
        _arrays(
            r=r,
            potential=potential,
            r_curve=r_curve,
            v_curve=theory.lj_energy(r_curve, c6, c12),
            total=total,
            turning_points=[r.min(), r.max()],
        ),
        draw,
        (5.6, 3.6),
    )


# -- 4. thermostats: same temperature, different physics -------------------------------------


def thermostat_fingerprint(
    temperature_K: float = 120.0, steps: int = 12000, n: int = 64, tau: float = 0.5
) -> Result:
    """Berendsen and a Nosé–Hoover chain hold the same temperature with different fluctuations.

    Only the Nosé–Hoover chain samples the canonical ensemble, whose kinetic energy is gamma
    distributed with relative width :math:`\\sqrt{2/N_{df}}`; Berendsen suppresses the
    fluctuations. The first quarter of each run is discarded.

    Parameters
    ----------
    temperature_K : float
        Bath temperature (K).
    steps : int
        Steps per thermostat; the kinetic energy is sampled every 5.
    n : int
        Number of atoms of the 2D gas.
    tau : float
        Coupling time of both thermostats (ps).

    Returns
    -------
    Result
        ``data``: ``berendsen`` and ``nose_hoover`` (kinetic energies), ``ke_axis`` and
        ``canonical_pdf``. Checks: the Berendsen width is under half the canonical one; the
        Nosé–Hoover width is within 15% of it.
    """
    system, recipe = lj_gas(n=n)
    recipe = recipe.update(control={"steps": steps})
    dof = 2 * n - 2
    width = theory.canonical_kinetic_energy_relative_width(dof)
    kinetic = {}
    algorithms: tuple[Literal["berendsen", "nose_hoover"], ...] = ("berendsen", "nose_hoover")
    for name in algorithms:
        coupled = with_bath(recipe, temperature_K, tau, name)
        ke = EnergyTimeseries(Simulation(system, coupled).run(steps, ene_freq=5)).kinetic
        kinetic[name] = ke[len(ke) // 4 :]
    ratios = {name: float(ke.std() / ke.mean()) for name, ke in kinetic.items()}
    checks = [
        Check(
            "Berendsen kinetic-energy σ/mean",
            ratios["berendsen"],
            width,
            band=(0.0, 0.5 * width),
        ),
        Check("Nosé–Hoover kinetic-energy σ/mean", ratios["nose_hoover"], width, rel_tol=0.15),
    ]
    mean_ke = dof / 2 * K_B * temperature_K
    ke_axis = np.linspace(0.6, 1.4, 300) * mean_ke

    def draw(result: Result, ax: Any) -> None:
        d = result.data
        for key, label, colour in (
            ("berendsen", f"Berendsen (τ {tau:g} ps)", _SECOND),
            ("nose_hoover", "Nosé–Hoover chain", _MEASURED),
        ):
            ke = d[key]
            ax.hist(
                ke,
                bins=40,
                density=True,
                alpha=0.5,
                color=colour,
                label=f"{label}: σ/mean {ke.std() / ke.mean():.3f}",
            )
        ax.plot(
            d["ke_axis"],
            d["canonical_pdf"],
            color=_TOTAL,
            lw=2,
            label=f"canonical: σ/mean {width:.3f}",
        )
        ax.set(xlabel="kinetic energy (kJ/mol)", ylabel="probability density", title=result.title)
        ax.legend(frameon=False, fontsize=8)

    return Result(
        f"Same temperature ({temperature_K:.0f} K), different physics",
        checks,
        _arrays(
            berendsen=kinetic["berendsen"],
            nose_hoover=kinetic["nose_hoover"],
            ke_axis=ke_axis,
            canonical_pdf=theory.canonical_kinetic_energy_pdf(ke_axis, temperature_K, dof),
        ),
        draw,
    )


# -- 5. how sure are we? block averaging -----------------------------------------------------


def block_error(
    duration_ps: float = 80.0,
    discard_ps: float = 20.0,
    n: int = 108,
    density: float = 0.8,
    temperature: float = 1.0,
) -> Result:
    """Successive samples are correlated, so the naive error of a mean is too small.

    A 3D Lennard-Jones liquid under a Nosé–Hoover chain; the potential energy is sampled every
    5 steps and its mean's standard error estimated with
    :meth:`gromos.EnergyTimeseries.block_average` for blocks of 1, 2, 4, … samples, up to a
    sixteenth of the run. Block size 1 is the naive error.

    Parameters
    ----------
    duration_ps : float
        Length of the run (ps).
    discard_ps : float
        Equilibration discarded from the start (ps); the liquid starts as a lattice.
    n : int
        Number of atoms.
    density, temperature : float
        Reduced density and temperature.

    Returns
    -------
    Result
        ``data``: ``block_ps`` (block length), ``error`` and ``naive``. Checks: the error at
        the longest block is more than twice the naive one; the curve rises up to it, with at
        most one dip, of at most 10%.
    """
    system, recipe = lj_gas(n=n, density=density, temperature=temperature, dim=3)
    dt = recipe.control["dt"]
    steps = round(duration_ps / dt)
    every = 5
    recipe = with_bath(
        recipe.update(control={"steps": steps}), temperature * EPSILON / K_B, 0.5, "nose_hoover"
    )
    rows = Simulation(system, recipe).run(steps, ene_freq=every)
    ts = EnergyTimeseries(rows[round(discard_ps / (dt * every)) :])
    sizes = [2**k for k in range(int(math.log2(max(1, len(ts) // 16))) + 1)]
    error = np.asarray([ts.block_average("potential", size)[1] for size in sizes])
    naive = error[0]
    drops = [(a - b) / a for a, b in zip(error[:-1], error[1:], strict=True) if b < a]
    checks = [
        Check("longest-block error / naive error", error[-1] / naive, 1.0, band=(2.0, math.inf)),
        Check("number of dips in the block-error curve", len(drops), 0, band=(0, 1)),
        Check("largest dip in the block-error curve", max(drops, default=0.0), 0.0, band=(0, 0.1)),
    ]

    def draw(result: Result, ax: Any) -> None:
        d = result.data
        ax.semilogx(d["block_ps"], d["error"], "o-", color=_MEASURED, label="block standard error")
        ax.axhline(
            float(d["naive"]), ls="--", color=_THEORY, label="naive error (samples independent)"
        )
        ax.set(
            xlabel="block length (ps)",
            ylabel="error of mean potential (kJ/mol)",
            title=result.title,
        )
        ax.text(
            0.03,
            0.93,
            f"longest block: {d['error'][-1] / float(d['naive']):.1f} × the naive error",
            transform=ax.transAxes,
            fontsize=9,
        )
        ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, 0.88))

    return Result(
        "How sure are we? Block averaging",
        checks,
        _arrays(block_ps=np.asarray(sizes) * every * dt, error=error, naive=naive),
        draw,
        (5.6, 3.6),
    )


# -- 6. speeds become Maxwell–Boltzmann ------------------------------------------------------


def maxwell_boltzmann(
    n: int = 64, temperature: float = 1.0, steps: int = 7000, density: float = 0.3
) -> Result:
    """Every atom starts at the same speed; collisions produce the Maxwell–Boltzmann law.

    Speeds are taken at the start, over a short window around 0.6 ps, and every 20 steps over
    the second half of the run, whose mean temperature sets the theory curve.

    That temperature counts :math:`2n - 2` degrees of freedom (the momentum is conserved) while
    the speeds are those of :math:`n` atoms, so :math:`\\langle v^2 \\rangle` is
    :math:`(n-1)/n` of :math:`2 k_B T / m` and the mean speed sits about 0.8% below theory for
    64 atoms, before sampling noise.

    Parameters
    ----------
    n : int
        Number of atoms of the 2D gas.
    temperature : float
        Reduced temperature of the initial speeds.
    steps : int
        Number of constant-energy steps.
    density : float
        Reduced density.

    Returns
    -------
    Result
        ``data``: ``start``, ``early``, ``late`` (speeds, nm/ps), ``temperature`` (late mean,
        K), ``v`` and ``pdf`` (the 2D Maxwell–Boltzmann curve). Check: the late mean speed is
        within 3% of :math:`\\sqrt{\\pi k_B T / 2m}`.
    """
    system, recipe = lj_gas(n=n, density=density, temperature=temperature, speeds="equal")
    dt = recipe.control["dt"]
    sim = Simulation(system, recipe.update(control={"steps": steps}))

    def speeds() -> np.ndarray:
        return np.linalg.norm(sim.velocities[:, :2], axis=1)

    early_step = round(0.6 / dt)
    start, early, late, temperatures = speeds(), [], [], []
    for k in range(1, steps + 1):
        sim.step(1)
        if abs(k - early_step) <= 5:
            early.append(speeds())
        if k >= steps // 2 and k % 20 == 0:
            late.append(speeds())
            temperatures.append(sim.temperature)
    late_speeds = np.concatenate(late)
    t_late = float(np.mean(temperatures))
    expected = theory.maxwell_boltzmann_mean_speed(t_late, ARGON_MASS, 2)
    check = Check("late mean speed (nm/ps)", late_speeds.mean(), expected, rel_tol=0.03)
    v_max = 2.5 * expected
    v = np.linspace(0.0, v_max, 300)

    def draw(result: Result, axes: Any) -> None:
        d = result.data
        t_early = early_step * dt
        panels = (
            ("t = 0", d["start"]),
            (f"t ≈ {t_early:.1f} ps", d["early"]),
            (f"second half (from {steps // 2 * dt:.0f} ps)", d["late"]),
        )
        for ax, (label, sample) in zip(axes, panels, strict=True):
            ax.hist(
                sample,
                bins=30,
                range=(0, v_max),
                density=True,
                color=_MEASURED,
                alpha=0.55,
                label="simulation",
            )
            ax.plot(
                d["v"],
                d["pdf"],
                color=_TOTAL,
                lw=1.8,
                label=f"Maxwell–Boltzmann, {float(d['temperature']):.0f} K",
            )
            ax.set(title=label, xlabel="speed (nm/ps)")
        axes[0].set_ylabel("probability density")
        axes[0].set_ylim(0, 2.5 * d["pdf"].max())  # the t = 0 spike runs off the top
        axes[2].text(
            0.97,
            0.6,
            f"mean {d['late'].mean():.4f}\ntheory {expected:.4f} nm/ps",
            transform=axes[2].transAxes,
            ha="right",
            fontsize=8,
        )
        axes[2].legend(frameon=False, fontsize=8)
        axes[0].figure.suptitle(result.title, fontsize=11, fontweight="bold")

    return Result(
        "Every atom starts at the same speed; collisions produce Maxwell–Boltzmann",
        [check],
        _arrays(
            start=start,
            early=np.concatenate(early),
            late=late_speeds,
            temperature=t_late,
            v=v,
            pdf=theory.maxwell_boltzmann_speed_pdf(v, t_late, ARGON_MASS, 2),
        ),
        draw,
        (9.6, 3.4),
        3,
    )


# -- 7. structure of a liquid ----------------------------------------------------------------


def _average_rdf(
    sim: Simulation,
    atoms: list[int],
    box: tuple[float, float, float],
    frames: int,
    every: int,
    n_bins: int,
    r_max: float,
) -> tuple[np.ndarray, np.ndarray]:
    curves = []
    for _ in range(frames):
        sim.step(every)
        r, g = rdf(sim.positions, atoms, atoms, n_bins, r_max, box)
        curves.append(g)
    return r, np.mean(curves, axis=0)


def _draw_rdf(xlabel: str, scale: float, tail_from: float) -> Callable[[Result, Any], None]:
    def draw(result: Result, ax: Any) -> None:
        d = result.data
        ax.plot(d["r"] / scale, d["g"], color=_MEASURED, lw=2, label="g(r) from gromos.rdf")
        ax.axhline(1.0, color=_THEORY, lw=0.8, label="ideal gas: g = 1")
        ax.axvspan(tail_from / scale, d["r"][-1] / scale, color=_THEORY, alpha=0.08)
        peak = int(np.argmax(d["g"]))
        ax.annotate(
            f"first peak {d['r'][peak] / scale:.2f}, height {d['g'][peak]:.2f}",
            xy=(d["r"][peak] / scale, d["g"][peak]),
            xytext=(d["r"][peak] / scale * 1.3, d["g"][peak] * 0.9),
            arrowprops={"arrowstyle": "->"},
            fontsize=8,
        )
        ax.text(
            0.97,
            0.25,
            f"tail mean {float(d['tail']):.3f}",
            transform=ax.transAxes,
            ha="right",
            fontsize=8,
        )
        ax.set(xlabel=xlabel, ylabel="g(r)", title=result.title)
        ax.legend(frameon=False, fontsize=8)

    return draw


def liquid_structure(
    n: int = 216,
    density: float = 0.8,
    temperature: float = 1.0,
    frames: int = 20,
    equilibration_steps: int = 2000,
    every: int = 50,
    n_bins: int = 100,
) -> Result:
    """A liquid has shells of neighbours: g(r) peaks near σ and levels off at 1.

    A 3D Lennard-Jones liquid under a Nosé–Hoover chain; :func:`gromos.rdf` over all atoms,
    averaged over ``frames`` frames, out to half the box.

    Parameters
    ----------
    n : int
        Number of atoms; at the default density, 216 gives a box wide enough to see g(r)
        beyond 2.5σ.
    density, temperature : float
        Reduced density and temperature.
    frames : int
        Number of frames averaged.
    equilibration_steps : int
        Steps run before the first frame; the liquid starts as a lattice.
    every : int
        Steps between frames.
    n_bins : int
        Number of g(r) bins.

    Returns
    -------
    Result
        ``data``: ``r`` (nm), ``g``, ``tail`` (mean g beyond 2.5σ). Checks: the first peak lies
        between 1.0σ and 1.2σ and is higher than 1.5; the tail is within [0.9, 1.1].
    """
    system, recipe = lj_gas(n=n, density=density, temperature=temperature, dim=3)
    recipe = with_bath(
        recipe.update(control={"steps": equilibration_steps + frames * every}),
        temperature * EPSILON / K_B,
        0.5,
        "nose_hoover",
    )
    sim = Simulation(system, recipe)
    sim.step(equilibration_steps)
    box = system.box
    r, g = _average_rdf(sim, list(range(n)), box, frames, every, n_bins, min(box) / 2)
    peak = int(np.argmax(g))
    tail = float(g[r > 2.5 * SIGMA].mean())
    r_min = 2 ** (1 / 6)
    checks = [
        Check("first peak position (σ)", r[peak] / SIGMA, r_min, band=(1.0, 1.2)),
        Check("first peak height", g[peak], None, band=(1.5, math.inf)),
        Check("mean g(r) beyond 2.5σ", tail, 1.0, band=(0.9, 1.1)),
    ]
    return Result(
        "A liquid has shells of neighbours",
        checks,
        _arrays(r=r, g=g, tail=tail),
        _draw_rdf("r / σ", SIGMA, 2.5 * SIGMA),
        (5.6, 3.6),
    )


def water_structure(
    topology_path: str,
    configuration_path: str,
    recipe_path: str,
    equilibration_steps: int = 2500,
    frames: int = 20,
    every: int = 25,
    n_bins: int = 45,
    r_max: float = 0.9,
    tail_from: float = 0.6,
    oxygen: str = "OW",
) -> Result:
    """The oxygen–oxygen g(r) of a water box given as GROMOS files.

    Parameters
    ----------
    topology_path, configuration_path, recipe_path : str
        GROMOS topology, coordinates and ``.imd`` input; the box must be rectangular.
    equilibration_steps : int
        Steps run before the first frame.
    frames : int
        Number of frames averaged.
    every : int
        Steps between frames.
    n_bins : int
        Number of g(r) bins.
    r_max : float
        Largest distance (nm); at most half the shortest box edge.
    tail_from : float
        Distance (nm) beyond which g(r) is averaged for the tail check.
    oxygen : str
        Atom name of the water oxygen.

    Returns
    -------
    Result
        ``data``: ``r`` (nm), ``g``, ``tail``. Checks: the first peak lies in
        [0.26, 0.29] nm; the mean g(r) beyond ``tail_from`` is within [0.9, 1.1].
    """
    system = System.from_files(topology_path, configuration_path)
    recipe = Recipe.from_imd(recipe_path)
    atoms = [i for i, name in enumerate(system.topology.atom_names) if name == oxygen]
    sim = Simulation(system, recipe.update(control={"steps": equilibration_steps + frames * every}))
    sim.step(equilibration_steps)
    r, g = _average_rdf(sim, atoms, system.box, frames, every, n_bins, r_max)
    peak = int(np.argmax(g))
    tail = float(g[r > tail_from].mean())
    checks = [
        Check("first O–O peak position (nm)", r[peak], 0.28, band=(0.26, 0.29)),
        Check(f"mean g(r) beyond {tail_from:g} nm", tail, 1.0, band=(0.9, 1.1)),
    ]
    return Result(
        "Water structure: g(r) of oxygen",
        checks,
        _arrays(r=r, g=g, tail=tail),
        _draw_rdf("O–O distance r (nm)", 1.0, tail_from),
        (5.6, 3.6),
    )


# -- 8. a hot gas meets a cold bath ----------------------------------------------------------


def _heat_to_bath(
    time: np.ndarray, temperature: np.ndarray, bath_K: float, tau: float, dof: int
) -> np.ndarray:
    """Energy (kJ/mol) a Berendsen bath has taken from the gas up to each sample.

    Berendsen removes kinetic energy at the rate :math:`(N_{df} k_B / 2\\tau)(T - T_0)`; the
    rate is integrated over the samples with the trapezoid rule.
    """
    rate = dof * K_B / (2.0 * tau) * (np.asarray(temperature, dtype=float) - bath_K)
    pieces = 0.5 * (rate[1:] + rate[:-1]) * np.diff(time)
    return np.concatenate([[0.0], np.cumsum(pieces)])


def _relaxation_time(time: np.ndarray, temperature: np.ndarray, bath_K: float) -> float:
    """Time constant of :math:`T - T_0`, fitted from the start until it has fallen to 20%."""
    deviation = (temperature - bath_K) / (temperature[:5].mean() - bath_K)
    below = np.flatnonzero(deviation < 0.2)
    end = int(below[0]) if below.size else len(deviation)
    slope = np.polyfit(time[:end], np.log(deviation[:end]), 1)[0]
    return float(-1.0 / slope)


def heat_bath(
    temperature_K: float = 120.0,
    start_K: float = 240.0,
    taus: Sequence[float] = (0.2, 1.0),
    steps: int = 8000,
    equilibration: int = 1500,
    n: int = 64,
    seed: int = 7,
) -> Result:
    """A hot 2D gas is put in contact with a colder bath; the heat it loses, the bath takes.

    The gas is equilibrated at ``start_K`` (Berendsen, τ 0.1 ps), and from that one state it is
    quenched three times to ``temperature_K``: Berendsen with each τ in ``taus``, and a
    Nosé–Hoover chain with the last τ. Energies are sampled every step.

    Berendsen scales the velocities each step by :math:`\\lambda^2 = 1 + (\\Delta t/\\tau)
    (T_0/T - 1)`, so it takes kinetic energy out at the rate :math:`(N_{df} k_B / 2\\tau)
    (T - T_0)` and the heat given to the bath is

    .. math::

        Q(t) = \\frac{N_{df} k_B}{2\\tau} \\int_0^t \\bigl(T(t') - T_0\\bigr)\\,dt',

    which by the first law equals the drop of the gas's total energy. The temperature relaxes
    without overshoot, but more slowly than τ: the potential energy re-equilibrates too, so the
    observed time constant is :math:`\\tau \\cdot 2 C_V / (N_{df} k_B)` (Berendsen et al.
    1984), at least τ and τ itself only for an ideal gas; the same factor is the ratio of the
    total-energy drop to the kinetic-energy drop over the quench. The fitted τ is therefore
    checked against a band above τ, not for equality. The Nosé–Hoover chain reaches the same
    temperature but overshoots and oscillates.

    Parameters
    ----------
    temperature_K : float
        Bath temperature (K).
    start_K : float
        Temperature of the hot gas (K).
    taus : sequence of float
        Berendsen coupling times (ps); the last is also the Nosé–Hoover one.
    steps : int
        Steps of each quench.
    equilibration : int
        Steps at ``start_K`` before the quenches.
    n : int
        Number of atoms of the 2D gas.
    seed : int
        Seed of the initial velocities.

    Returns
    -------
    Result
        ``data``: ``time`` (ps); ``taus``; ``berendsen`` (temperatures, one row per τ, K);
        ``nose_hoover`` (K); ``total_energy`` and ``heat_to_bath`` (kJ/mol, one row per
        Berendsen τ); ``kinetic_only`` (:math:`T_0 + (T_{start} - T_0) e^{-t/\\tau}` per τ) and
        ``tau_fit`` (fitted time constants, ps). Checks: the mean temperature of the second half
        is within 3% of the bath for each Berendsen run and within 5% for Nosé–Hoover; for each
        Berendsen run the drop of the total energy to its mean over the last quarter is within
        5% of the mean :math:`Q` there; the fitted τ over τ lies in [0.9, 2.0]; the last fitted τ
        over the first lies in [0.5, 1.6] times the ratio of the τ.
    """
    taus = tuple(float(tau) for tau in taus)
    system, recipe = lj_gas(n=n, dim=2, temperature=start_K * K_B / EPSILON, seed=seed)
    dof = 2 * n - 2
    hot = Simulation(
        system, with_bath(recipe.update(control={"steps": equilibration}), start_K, 0.1)
    )
    hot.step(equilibration)
    start = system_from_frame(system, hot.positions, hot.velocities)
    quench = recipe.update(control={"steps": steps})

    def run(tau: float, algorithm: Literal["berendsen", "nose_hoover"]) -> EnergyTimeseries:
        sim = Simulation(start, with_bath(quench, temperature_K, tau, algorithm))
        return EnergyTimeseries(sim.run(steps, ene_freq=1))

    runs = [run(tau, "berendsen") for tau in taus]
    chain = run(taus[-1], "nose_hoover")
    time = runs[0].time
    berendsen = np.asarray([theory.temperature_from_kinetic(ts.kinetic, dof) for ts in runs])
    nose_hoover = np.asarray(theory.temperature_from_kinetic(chain.kinetic, dof))
    total = np.asarray([ts.total for ts in runs])
    heat = np.asarray(
        [
            _heat_to_bath(time, t, temperature_K, tau, dof)
            for t, tau in zip(berendsen, taus, strict=True)
        ]
    )
    t_start = berendsen[:, :5].mean(axis=1)
    kinetic_only = temperature_K + (t_start[:, None] - temperature_K) * np.exp(
        -time[None, :] / np.asarray(taus)[:, None]
    )
    tau_fit = np.asarray([_relaxation_time(time, t, temperature_K) for t in berendsen])
    half, last_quarter = len(time) // 2, len(time) - len(time) // 4

    checks = [
        Check(
            f"Berendsen τ {tau:g} ps: mean T, second half (K)",
            t[half:].mean(),
            temperature_K,
            rel_tol=0.03,
        )
        for t, tau in zip(berendsen, taus, strict=True)
    ]
    checks.append(
        Check(
            "Nosé–Hoover: mean T, second half (K)",
            nose_hoover[half:].mean(),
            temperature_K,
            rel_tol=0.05,
        )
    )
    checks += [
        Check(
            f"Berendsen τ {tau:g} ps: total-energy drop against heat to the bath (kJ/mol)",
            e[0] - e[last_quarter:].mean(),
            q[last_quarter:].mean(),
            rel_tol=0.05,
        )
        for e, q, tau in zip(total, heat, taus, strict=True)
    ]
    checks += [
        Check(f"Berendsen τ {tau:g} ps: fitted τ / τ", fit / tau, 1.0, band=(0.9, 2.0))
        for fit, tau in zip(tau_fit, taus, strict=True)
    ]
    ratio = taus[-1] / taus[0]
    checks.append(
        Check(
            "fitted τ, last over first",
            tau_fit[-1] / tau_fit[0],
            ratio,
            band=(0.5 * ratio, 1.6 * ratio),
        )
    )

    def draw(result: Result, axes: Any) -> None:
        d = result.data
        left, right = axes
        colours = (_SECOND, _TOTAL, _KINETIC, _POTENTIAL)
        left.axhline(temperature_K, ls="--", color=_THEORY, lw=1, label=f"bath {temperature_K:g} K")
        for k, (tau, fit) in enumerate(zip(d["taus"], d["tau_fit"], strict=True)):
            left.plot(
                d["time"],
                d["kinetic_only"][k],
                ls=":",
                color=_THEORY,
                lw=1.2,
                label="kinetic only: exp(−t/τ)" if k == 0 else None,
            )
            left.plot(
                d["time"],
                d["berendsen"][k],
                color=colours[k % len(colours)],
                lw=1.2,
                label=f"Berendsen τ {tau:g} ps (fit {fit:.2f} ps)",
            )
        left.plot(
            d["time"],
            d["nose_hoover"],
            color=_MEASURED,
            lw=1.2,
            label=f"Nosé–Hoover chain τ {d['taus'][-1]:g} ps",
        )
        left.set(
            xlim=(0, min(d["time"][-1], 10 * d["taus"][-1])),
            xlabel="time (ps)",
            ylabel="temperature (K)",
            title=f"The gas cools (first {min(d['time'][-1], 10 * d['taus'][-1]):.0f} ps)",
        )
        left.legend(frameon=False, fontsize=8)
        drop = d["total_energy"][-1, 0] - d["total_energy"][-1]
        right.plot(d["time"], drop, color=_TOTAL, lw=2, label="E_total(0) − E_total(t): the gas")
        right.plot(
            d["time"],
            d["heat_to_bath"][-1],
            ls="--",
            color=_KINETIC,
            lw=1.6,
            label="Q(t): heat taken by the bath",
        )
        right.set(
            xlabel="time (ps)",
            ylabel="energy (kJ/mol)",
            title=f"First law, Berendsen τ {d['taus'][-1]:g} ps",
        )
        right.legend(frameon=False, fontsize=8, loc="lower right")
        left.figure.suptitle(result.title, fontsize=11, fontweight="bold")

    return Result(
        "A hot gas meets a cold bath",
        checks,
        _arrays(
            time=time,
            taus=taus,
            berendsen=berendsen,
            nose_hoover=nose_hoover,
            total_energy=total,
            heat_to_bath=heat,
            kinetic_only=kinetic_only,
            tau_fit=tau_fit,
        ),
        draw,
        (10.4, 3.8),
        2,
        False,
    )
