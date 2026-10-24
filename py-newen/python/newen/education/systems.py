"""Toy systems for teaching, built in memory as GROMOS text and read by the engine.

Every builder writes the topology and configuration a GROMOS user would write by hand, reads
them through :meth:`newen.Topology.from_string` and :meth:`newen.Configuration.from_string`,
and returns a :class:`newen.System` with a :class:`newen.Recipe` to run it. Nothing here
computes a force or an energy: the engine does.

The Lennard-Jones systems use reduced units mapped onto argon: lengths in :data:`SIGMA`,
energies in :data:`EPSILON`, masses in :data:`ARGON_MASS`, so a reduced temperature
:math:`T^*` is :math:`T = T^* \\varepsilon / k_B` (about 120 K for :math:`T^* = 1`) and the
reduced time unit is :data:`TAU`.

Examples
--------
>>> from newen import Simulation
>>> from newen.education import lj_gas
>>> system, recipe = lj_gas(n=64, density=0.45, temperature=1.0)
>>> sim = Simulation(system, recipe)
>>> sim.step(100)
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from newen import Configuration, Recipe, System, Topology

__all__ = [
    "ARGON_MASS",
    "EPSILON",
    "HYDROGEN_MASS",
    "K_B",
    "OXYGEN_MASS",
    "SIGMA",
    "TAU",
    "argon_pair",
    "diatomic",
    "lj_gas",
    "system_from_frame",
    "with_bath",
]

#: Boltzmann constant (kJ/mol/K): the engine's own value (newen-core ``units::kB``, gromosXX
#: ``math::k_Boltzmann``), so temperatures computed here match ``Simulation.temperature``.
K_B = 0.00831441
#: Lennard-Jones diameter of argon (nm).
SIGMA = 0.3405
#: Lennard-Jones well depth of argon (kJ/mol).
EPSILON = 0.9961
#: Mass of argon (u).
ARGON_MASS = 39.948
#: Mass of oxygen (u).
OXYGEN_MASS = 15.9994
#: Mass of hydrogen (u).
HYDROGEN_MASS = 1.008
#: Reduced time unit :math:`\sigma \sqrt{m / \varepsilon}` (ps), about 2.16 ps.
TAU = SIGMA * math.sqrt(ARGON_MASS / EPSILON)

_C6 = 4.0 * EPSILON * SIGMA**6
_C12 = 4.0 * EPSILON * SIGMA**12
_EMPTY_BLOCKS = (
    "BONDSTRETCHTYPE",
    "BONDH",
    "BOND",
    "BONDANGLEBENDTYPE",
    "BONDANGLEH",
    "BONDANGLE",
    "IMPDIHEDRALTYPE",
    "IMPDIHEDRALH",
    "IMPDIHEDRAL",
    "TORSDIHEDRALTYPE",
    "DIHEDRALH",
    "DIHEDRAL",
    "CROSSDIHEDRALH",
    "CROSSDIHEDRAL",
)


def _topology_text(
    title: str,
    atoms: list[tuple[str, int, float, list[int]]],
    type_names: list[str],
    lj: dict[tuple[int, int], tuple[float, float]],
    residue: str,
    bonds: list[tuple[int, int, int]] | None = None,
    bond_types: list[tuple[float, float, float]] | None = None,
    molecule_ends: list[int] | None = None,
) -> str:
    """A GROMOS topology: one residue, every atom its own charge group, no charges.

    ``atoms`` holds ``(name, 1-based IAC, mass, 1-based excluded atoms)``; ``lj`` maps a
    1-based type pair ``(i, j)``, ``i <= j``, to ``(c6, c12)``; ``bonds`` are
    ``(i, j, type)`` 1-based and go to the BONDH block when either atom is a hydrogen;
    ``molecule_ends`` is the last atom of each molecule (default: every atom its own).
    """
    n = len(atoms)
    lines = [
        "TITLE",
        title,
        "END",
        "PHYSICALCONSTANTS",
        "138.9354",
        "0.0635078",
        "299792.458",
        f"{K_B}",
        "END",
        "TOPVERSION",
        "2.0",
        "END",
        "ATOMTYPENAME",
        str(len(type_names)),
        *type_names,
        "END",
        "RESNAME",
        "1",
        residue,
        "END",
        "SOLUTEATOM",
        str(n),
    ]
    for i, (name, iac, mass, excluded) in enumerate(atoms, start=1):
        exclusions = "".join(f"{j:6d}" for j in excluded)
        lines.append(
            f"{i:5d}{1:5d}{name:>6s}{iac:4d}{mass:12.5f}{0.0:12.5f}{1:3d}{len(excluded):3d}"
            + exclusions
        )
        lines.append(f"{0:47d}")
    lines.append("END")
    bonds = bonds or []
    bond_types = bond_types or []
    with_h = [b for b in bonds if "H" in (atoms[b[0] - 1][0][0], atoms[b[1] - 1][0][0])]
    without_h = [b for b in bonds if b not in with_h]
    for block in _EMPTY_BLOCKS:
        if block == "BONDSTRETCHTYPE":
            rows = [f"{cb:.6e} {chb:.6e} {b0:.6f}" for cb, chb, b0 in bond_types]
        elif block in ("BONDH", "BOND"):
            chosen = with_h if block == "BONDH" else without_h
            rows = [f"{i:6d}{j:6d}{t:6d}" for i, j, t in chosen]
        else:
            rows = []
        lines += [block, str(len(rows)), *rows, "END"]
    lines += ["LJPARAMETERS", str(len(lj))]
    for (i, j), (c6, c12) in sorted(lj.items()):
        lines.append(f"{i:5d}{j:5d}  {c12:.6e}  {c6:.6e}  {c12:.6e}  {c6:.6e}")
    lines.append("END")
    ends = molecule_ends or list(range(1, n + 1))
    molecules = [str(len(ends))] + [
        " ".join(str(k) for k in ends[s : s + 10]) for s in range(0, len(ends), 10)
    ]
    lines += [
        "SOLUTEMOLECULES",
        *molecules,
        "END",
        "TEMPERATUREGROUPS",
        "1",
        str(n),
        "END",
        "PRESSUREGROUPS",
        *molecules,
        "END",
        "LJEXCEPTIONS",
        "0",
        "END",
        "SOLVENTATOM",
        "0",
        "END",
        "SOLVENTCONSTR",
        "0",
        "END",
    ]
    return "\n".join(lines) + "\n"


def _configuration_text(
    title: str,
    names: list[str],
    residue: str,
    positions: np.ndarray,
    velocities: np.ndarray,
    box: float | tuple[float, float, float] | None,
) -> str:
    """A GROMOS96 configuration; ``box`` is a cube edge or three edges (nm), ``None`` for vacuum."""

    def row(i: int, values: np.ndarray) -> str:
        return f"{1:5d} {residue:<5s} {names[i]:<5s}{i + 1:7d}" + "".join(
            f"{v:15.9f}" for v in values
        )

    n = len(names)
    edges = (0.0, 0.0, 0.0) if box is None else (box,) * 3 if isinstance(box, float) else box
    return (
        "\n".join(
            [
                "TITLE",
                title,
                "END",
                "POSITION",
                "# first 24 chars ignored",
                *(row(i, positions[i]) for i in range(n)),
                "END",
                "VELOCITY",
                *(row(i, velocities[i]) for i in range(n)),
                "END",
                "GENBOX",
                "0" if box is None else "1",
                " ".join(f"{e:.9f}" for e in edges),
                "90.0 90.0 90.0",
                "0.0 0.0 0.0",
                "0.0 0.0 0.0",
                "END",
            ]
        )
        + "\n"
    )


def _system(topology: str, configuration: str) -> System:
    return System(Topology.from_string(topology), Configuration.from_string(configuration))


def _balanced_directions(n: int, dim: int, rng: np.random.Generator) -> np.ndarray:
    """``n`` random unit vectors that sum to zero, in random order.

    In 2D, angles evenly spaced round the circle from a random offset. In 3D, antipodal pairs
    of random directions, plus three at 120 degrees in a random plane when ``n`` is odd.
    """
    if n < 2:
        raise ValueError("equal speeds with zero net momentum need at least two atoms")
    if dim == 2:
        angles = rng.uniform(0.0, 2.0 * math.pi) + 2.0 * math.pi * np.arange(n) / n
        unit = np.column_stack([np.cos(angles), np.sin(angles)])
    else:
        triple = n % 2 * 3
        random = rng.normal(size=((n - triple) // 2, 3))
        random /= np.linalg.norm(random, axis=1, keepdims=True)
        unit = np.concatenate([random, -random])
        if triple:
            a, b = np.linalg.qr(rng.normal(size=(3, 2)))[0].T
            angles = 2.0 * math.pi * np.arange(3) / 3
            in_plane = np.cos(angles)[:, None] * a + np.sin(angles)[:, None] * b
            unit = np.concatenate([unit, in_plane])
    return rng.permutation(unit)


def lj_gas(
    n: int = 64,
    density: float = 0.45,
    temperature: float = 1.0,
    dim: Literal[2, 3] = 2,
    seed: int = 7,
    speeds: Literal["maxwell", "equal"] = "maxwell",
) -> tuple[System, Recipe]:
    """Argon atoms on a lattice in a periodic box, in two or three dimensions.

    The atoms start on a square (``dim=2``) or simple-cubic (``dim=3``) lattice filling a
    cubic box whose edge gives the requested reduced density. In two dimensions every atom
    shares one :math:`z` and has :math:`v_z = 0`; nothing pushes along :math:`z`, so the gas
    stays planar, and ``boundary.ndfmin = n + 2`` removes the :math:`n` frozen :math:`z`
    degrees of freedom and the two of the conserved momentum, leaving :math:`2n - 2`. In three
    dimensions ``ndfmin = 3`` (the conserved momentum), leaving :math:`3n - 3`. The initial
    velocities carry no net momentum, and the centre-of-mass motion is not touched again.

    Parameters
    ----------
    n : int
        Number of atoms.
    density : float
        Reduced density :math:`\\rho^* = n \\sigma^d / L^d`.
    temperature : float
        Reduced temperature :math:`T^* = k_B T / \\varepsilon` of the initial velocities.
    dim : {2, 3}
        Number of dimensions the atoms move in.
    seed : int
        Seed of the initial velocities.
    speeds : {"maxwell", "equal"}
        ``"maxwell"`` draws each velocity component from a Gaussian of variance
        :math:`k_B T / m`; ``"equal"`` gives every atom the speed
        :math:`\\sqrt{d \\, k_B T / m}` in a random direction, the directions chosen so
        that they sum to zero.

    Returns
    -------
    System
        The gas.
    Recipe
        Constant energy, time step :math:`0.003\\,\\tau`, 1000 steps, cut-off
        :math:`2.5\\sigma` with the pair list rebuilt every step. Change it with
        ``recipe.update(...)``.

    Raises
    ------
    ValueError
        For a ``dim`` other than 2 or 3, an unknown ``speeds``, or a box shorter than twice the
        cut-off.
    """
    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim!r}")
    if speeds not in ("maxwell", "equal"):
        raise ValueError(f"speeds must be 'maxwell' or 'equal', got {speeds!r}")
    edge = SIGMA * (n / density) ** (1.0 / dim)
    cutoff = 2.5 * SIGMA
    if edge < 2.0 * cutoff:
        raise ValueError(
            f"box edge {edge:.3f} nm is shorter than twice the cut-off ({2 * cutoff:.3f} nm); "
            "use more atoms or a lower density"
        )

    side = math.ceil(n ** (1.0 / dim) - 1e-9)
    index = np.arange(n)
    spacing = edge / side
    positions = np.full((n, 3), edge / 2.0)
    positions[:, 0] = (index % side + 0.5) * spacing
    positions[:, 1] = (index // side % side + 0.5) * spacing
    if dim == 3:
        positions[:, 2] = (index // side**2 + 0.5) * spacing

    rng = np.random.default_rng(seed)
    kt_over_m = temperature * EPSILON / ARGON_MASS
    velocities = np.zeros((n, 3))
    if speeds == "maxwell":
        velocities[:, :dim] = rng.normal(0.0, math.sqrt(kt_over_m), (n, dim))
        velocities[:, :dim] -= velocities[:, :dim].mean(axis=0)
    else:
        velocities[:, :dim] = math.sqrt(dim * kt_over_m) * _balanced_directions(n, dim, rng)

    names = ["AR"] * n
    topology = _topology_text(
        f"{n} argon atoms, {dim}D",
        [("AR", 1, ARGON_MASS, [])] * n,
        ["AR"],
        {(1, 1): (_C6, _C12)},
        "AR",
    )
    configuration = _configuration_text(
        f"{n} argon atoms, {dim}D, reduced density {density}, temperature {temperature}",
        names,
        "AR",
        positions,
        velocities,
        edge,
    )
    recipe = Recipe.nve(dt=0.003 * TAU, steps=1000).update(
        control={"com_motion": {"initial": 0, "every": 0}},
        boundary={"ndfmin": n + 2 if dim == 2 else 3},
        forcefield={
            "pairlist": {"cutoff_short": cutoff, "cutoff_long": cutoff, "every": 1},
            "electrostatics": {"rf_cutoff": cutoff},
        },
    )
    return _system(topology, configuration), recipe


def argon_pair(distance: float = 0.35) -> tuple[System, Recipe]:
    """Two argon atoms at rest in vacuum, ``distance`` apart along :math:`x`.

    Parameters
    ----------
    distance : float
        Initial separation (nm). Below the minimum :math:`2^{1/6}\\sigma \\approx 0.382` nm the
        pair is pushed apart and oscillates about the minimum.

    Returns
    -------
    System
        The pair, in a vacuum box.
    Recipe
        Constant energy, time step 0.002 ps, 4000 steps, the pair list rebuilt every step with
        a 1.4 nm cut-off.
    """
    atoms: list[tuple[str, int, float, list[int]]] = [
        ("AR1", 1, ARGON_MASS, []),
        ("AR2", 1, ARGON_MASS, []),
    ]
    topology = _topology_text(
        "two argon atoms", atoms, ["AR"], {(1, 1): (_C6, _C12)}, "AR2", molecule_ends=[2]
    )
    positions = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
    configuration = _configuration_text(
        f"two argon atoms {distance} nm apart",
        ["AR1", "AR2"],
        "AR2",
        positions,
        np.zeros((2, 3)),
        None,
    )
    recipe = Recipe.nve(dt=0.002, steps=4000).update(
        control={"com_motion": {"initial": 0, "every": 0}},
        boundary={"kind": "vacuum"},
        forcefield={
            "pairlist": {"cutoff_short": 0.8, "cutoff_long": 1.4, "every": 1},
            "electrostatics": {"rf_cutoff": 1.4},
        },
    )
    return _system(topology, configuration), recipe


def diatomic(
    m1: float = OXYGEN_MASS,
    m2: float = HYDROGEN_MASS,
    period_fs: float = 9.8,
    r0: float = 0.1,
    stretch: float = 0.01,
) -> tuple[System, Recipe]:
    """Two atoms joined by one harmonic bond, in vacuum, started stretched and at rest.

    The force constant is chosen so that the bond vibrates with the requested period:

    .. math::

        \\mu = \\frac{m_1 m_2}{m_1 + m_2}, \\qquad
        \\omega = \\frac{2\\pi}{T}, \\qquad
        K_H = \\mu \\omega^2,

    with :math:`V(b) = \\tfrac{1}{2} K_H (b - b_0)^2` (``covalent_form.harmonic_bonds``).
    The quartic constant of the same bond type is :math:`K_Q = K_H / (2 b_0^2)`, which has the
    same curvature at :math:`b_0`. The defaults are an O–H bond with a period of 9.8 fs, so
    leap-frog is stable up to :math:`\\Delta t = T/\\pi \\approx 3.12` fs.

    Parameters
    ----------
    m1, m2 : float
        Masses (u). The first atom is named ``O``, the second ``H``.
    period_fs : float
        Vibration period (fs).
    r0 : float
        Bond length at rest :math:`b_0` (nm).
    stretch : float
        Initial stretch :math:`b - b_0` (nm).

    Returns
    -------
    System
        The two atoms, in a vacuum box.
    Recipe
        Constant energy, time step 0.0005 ps, 2000 steps, harmonic bonds, no constraints.
    """
    mu = m1 * m2 / (m1 + m2)
    omega = 2.0 * math.pi / (period_fs * 1e-3)
    k_harmonic = mu * omega**2
    k_quartic = k_harmonic / (2.0 * r0**2)
    atoms = [("O", 1, m1, [2]), ("H", 2, m2, [])]
    zero = (0.0, 0.0)
    topology = _topology_text(
        "a harmonic diatomic",
        atoms,
        ["O", "H"],
        {(1, 1): zero, (1, 2): zero, (2, 2): zero},
        "OH",
        bonds=[(1, 2, 1)],
        bond_types=[(k_quartic, k_harmonic, r0)],
        molecule_ends=[2],
    )
    positions = np.array([[0.0, 0.0, 0.0], [r0 + stretch, 0.0, 0.0]])
    configuration = _configuration_text(
        f"a diatomic stretched by {stretch} nm", ["O", "H"], "OH", positions, np.zeros((2, 3)), None
    )
    recipe = Recipe.nve(dt=0.0005, steps=2000).update(
        control={"com_motion": {"initial": 0, "every": 0}},
        boundary={"kind": "vacuum"},
        forcefield={
            "covalent_form": {"harmonic_bonds": True},
            "pairlist": {"cutoff_short": 0.8, "cutoff_long": 1.4, "every": 1},
            "electrostatics": {"rf_cutoff": 1.4},
        },
    )
    return _system(topology, configuration), recipe


def with_bath(
    recipe: Recipe,
    temperature_K: float,
    tau: float = 0.5,
    algorithm: Literal["berendsen", "nose_hoover"] = "berendsen",
) -> Recipe:
    """``recipe`` coupled to one heat bath, with the thermostat block :meth:`Recipe.nvt` writes.

    Parameters
    ----------
    recipe : Recipe
        The recipe to couple; its time step and number of steps are kept.
    temperature_K : float
        Bath temperature (K).
    tau : float
        Coupling time (ps).
    algorithm : {"berendsen", "nose_hoover"}
        Berendsen weak coupling, or a Nosé–Hoover chain of three thermostats.

    Returns
    -------
    Recipe
        A new recipe; ``recipe`` is not changed.

    Raises
    ------
    ValueError
        For an unknown ``algorithm``.
    """
    engine_algorithm: str | dict[str, int]
    if algorithm == "berendsen":
        engine_algorithm = "berendsen"
    elif algorithm == "nose_hoover":
        engine_algorithm = {"nose_hoover_chain": 3}
    else:
        raise ValueError(f"algorithm must be 'berendsen' or 'nose_hoover', got {algorithm!r}")
    thermostat = (
        Recipe.nvt(
            dt=recipe.control["dt"], steps=recipe.control["steps"], temperature=temperature_K
        )
        .update(
            ensemble={
                "thermostat": {
                    "algorithm": engine_algorithm,
                    "baths": [{"tau": tau, "temperature": temperature_K}],
                }
            }
        )
        .ensemble["thermostat"]
    )
    if algorithm == "nose_hoover":
        # a fresh start: the thermostat variables are initialised to zero (NTINHT = 1), not read
        # from a configuration that has none
        return recipe.update(ensemble={"thermostat": thermostat}, control={"ntinht": 1})
    return recipe.update(ensemble={"thermostat": thermostat})


def system_from_frame(system: System, positions: np.ndarray, velocities: np.ndarray) -> System:
    """``system`` with new positions and velocities, such as a running simulation's.

    Parameters
    ----------
    system : System
        Supplies the topology, the atom names and the box (rectangular or vacuum).
    positions, velocities : numpy.ndarray
        ``(n_atoms, 3)`` arrays (nm, nm/ps), e.g. ``sim.positions`` and ``sim.velocities``.

    Returns
    -------
    System
        A new system; ``Simulation(system_from_frame(system, sim.positions, sim.velocities),
        recipe)`` continues from where ``sim`` stands.
    """
    box = system.box
    configuration = _configuration_text(
        "newen.education frame",
        system.topology.atom_names,
        "MOL",
        np.asarray(positions, dtype=float),
        np.asarray(velocities, dtype=float),
        None if not any(box) else box,
    )
    return System(system.topology, Configuration.from_string(configuration))
