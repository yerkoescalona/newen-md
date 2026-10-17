"""The preparation chain, taking and returning objects instead of file paths.

Each function calls the same Rust code its program's text path calls (a Rust test pins the two
together), so ``str(tools.sim_box(...))`` is exactly ``tools.run("sim_box", ...)``. Every input
that names a file may be given as a path *or* as the object it would have been read from.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

from ..newen import (
    Configuration,
    Topology,
    _build_box,
    _check_top,
    _com_top,
    _ion,
    _make_top,
    _pdb2g96,
    _red_top,
    _sim_box,
)

__all__ = [
    "make_top",
    "com_top",
    "red_top",
    "pdb2g96",
    "sim_box",
    "build_box",
    "ion",
    "check_top",
    "CheckReport",
]

_PathLike = str | os.PathLike


def _topology(value: Topology | _PathLike) -> Topology:
    """A `Topology`, whether given as one or as the file to read it from."""
    return value if isinstance(value, Topology) else Topology.from_file(os.fspath(value))


def _configuration(value: Configuration | _PathLike) -> Configuration:
    """A `Configuration`, whether given as one or as the file to read it from."""
    return value if isinstance(value, Configuration) else Configuration.from_file(os.fspath(value))


@dataclass(frozen=True)
class CheckReport:
    """What ``check_top`` found. Truthy when the topology passed.

    ``text`` is the report as the program prints it; ``errors`` and ``warnings`` are the same
    findings as lists.
    """

    errors: list[str]
    warnings: list[str]
    text: str

    def __bool__(self) -> bool:
        return not self.errors

    def __str__(self) -> str:
        return self.text

    def __repr__(self) -> str:
        return f"CheckReport(errors={len(self.errors)}, warnings={len(self.warnings)})"


def make_top(
    *,
    build: _PathLike | Sequence[_PathLike],
    param: _PathLike,
    seq: Sequence[str],
    solv: str | None = None,
) -> Topology:
    """Build a topology from building blocks (gromos++ ``make_top``).

    Args:
        build: Building-block (MTB) file, or several — later files extend the first.
        param: Force-field parameter (IFP) file.
        seq: The residue sequence, e.g. ``["NH3+", "ALA", "GLY", "COO-"]``.
        solv: Solvent building block to append, e.g. ``"H2O"``.

    >>> top = make_top(
    ...     build="54a7.mtb", param="54a7.ifp", seq=["NH3+", "ALA", "COO-"], solv="H2O"
    ... )  # doctest: +SKIP
    """
    files = [build] if isinstance(build, (str, os.PathLike)) else list(build)
    return _make_top(
        build=[os.fspath(f) for f in files],
        param=os.fspath(param),
        seq=list(seq),
        solv=solv,
    )


def com_top(
    *parts: Topology | _PathLike | tuple[int, Topology | _PathLike],
    param: int = 1,
    solv: int = 1,
) -> Topology:
    """Combine topologies (gromos++ ``com_top``).

    A part is a topology, or ``(repeat, topology)`` to include it several times. ``param`` and
    ``solv`` are 1-based indices into the repeat-expanded sequence, saying which topology
    contributes the force-field parameter blocks and the solvent.

    >>> combined = com_top(protein_top, (8, "na.top"), param=1, solv=1)  # doctest: +SKIP
    """
    expanded: list[tuple[int, Topology]] = []
    for part in parts:
        if isinstance(part, tuple):
            repeat, topo = part
            expanded.append((int(repeat), _topology(topo)))
        else:
            expanded.append((1, _topology(part)))
    return _com_top(expanded, param=param, solv=solv)


def red_top(topology: Topology | _PathLike, *, first: int = 1, last: int) -> Topology:
    """Keep only atoms ``first..last`` of a topology, 1-based and inclusive (gromos++ ``red_top``).

    >>> head = red_top(top, last=6)  # doctest: +SKIP
    """
    return _red_top(_topology(topology), first=first, last=last)


def pdb2g96(
    topology: Topology | _PathLike,
    *,
    pdb: _PathLike,
    lib: _PathLike | None = None,
    gch: bool = False,
    strict: bool = False,
) -> Configuration:
    """Put a PDB's coordinates into the topology's atom order (gromos++ ``pdb2g96``).

    Atoms the PDB does not have are set to (0, 0, 0). With ``strict=True`` any such atom is an
    error instead; otherwise the count is reported on stderr, as the program does.

    Args:
        lib: Residue/atom name library, when the PDB's names differ from the topology's.
        gch: Suppress the warnings about missing hydrogens.
    """
    conf, missing = _pdb2g96(
        _topology(topology),
        os.fspath(pdb),
        lib=os.fspath(lib) if lib is not None else None,
        gch=gch,
    )
    if strict and missing:
        raise ValueError(f"pdb2g96: {missing} atoms had no match in {os.fspath(pdb)}")
    return conf


def sim_box(
    topology: Topology | _PathLike,
    solute: Configuration | _PathLike,
    *,
    solvent: Configuration | _PathLike,
    minwall: float | Sequence[float] | None = None,
    thresh: float = 0.23,
    boxsize: bool = False,
) -> Configuration:
    """Solvate a solute in a box of pre-equilibrated solvent (gromos++ ``sim_box``).

    Args:
        minwall: Minimum solute-to-wall distance (nm): one value for a cubic box, three for a
            rectangular one. Mutually exclusive with ``boxsize``.
        thresh: Minimum solvent-to-solute distance (nm).
        boxsize: Take the box from the solute's own coordinates instead of ``minwall``.

    >>> solvated = sim_box(top, "solute.cnf", solvent="spc.cnf", minwall=1.0)  # doctest: +SKIP
    """
    if minwall is None:
        walls = None
    elif isinstance(minwall, (int, float)):
        walls = [float(minwall)]
    else:
        walls = [float(v) for v in minwall]
    return _sim_box(
        _topology(topology),
        _configuration(solute),
        _configuration(solvent),
        minwall=walls,
        thresh=thresh,
        boxsize=boxsize,
    )


def build_box(
    topology: Topology | _PathLike,
    molecule: Configuration | _PathLike,
    *,
    nsm: int,
    dens: float,
) -> Configuration:
    """Replicate one molecule on an ``nsm × nsm × nsm`` grid at a target density (gromos++
    ``build_box``).

    Args:
        nsm: Copies per dimension — the box holds ``nsm ** 3`` molecules.
        dens: Target density (kg/m³); it sets the box size.
    """
    return _build_box(_topology(topology), _configuration(molecule), nsm=nsm, dens=dens)


def ion(
    topology: Topology | _PathLike,
    system: Configuration | _PathLike,
    *,
    positive: tuple[int, str] | None = None,
    negative: tuple[int, str] | None = None,
    potential: float | None = None,
    mindist: float = 0.0,
) -> Configuration:
    """Replace solvent molecules by counter-ions (gromos++ ``ion``).

    Positive ions go where the electrostatic potential is lowest, negative ions where it is
    highest. **The topology is not updated** — as in gromos++, adding the ions' topology is the
    caller's job (``com_top``); ``System.add_ions`` does both.

    Args:
        positive: ``(count, name)``, e.g. ``(2, "NA+")``.
        negative: ``(count, name)``, e.g. ``(2, "CL-")``.
        potential: Cutoff for the electrostatic potential (nm). Required.
        mindist: Minimum distance between the ions placed (nm).

    >>> ionised = ion(top, solvated, positive=(2, "NA+"), potential=0.8)  # doctest: +SKIP
    """
    return _ion(
        _topology(topology),
        _configuration(system),
        positive=positive,
        negative=negative,
        potential=potential,
        mindist=mindist,
    )


def check_top(topology: Topology | _PathLike) -> CheckReport:
    """Run every consistency check over a topology (gromos++ ``check_top``).

    Returns a :class:`CheckReport`, which is falsy when the topology has errors — so
    ``if not check_top(top): ...`` reads the obvious way. Unlike the command-line program this
    does not fail; inspect ``.errors`` or ``print(report)``.

    >>> report = check_top("system.top")  # doctest: +SKIP
    >>> if not report:
    ...     print(report)  # doctest: +SKIP
    """
    errors, warnings, text = _check_top(_topology(topology))
    return CheckReport(errors=errors, warnings=warnings, text=text)
