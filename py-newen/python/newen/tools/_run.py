"""Calling a GROMOS program by name: keywords in, its own text out."""

from __future__ import annotations

import os
from typing import Any

from ..newen import ProgramError, _programs, _run_program

__all__ = ["run", "programs", "argv_from", "ProgramError"]

# What a keyword value may be. A path-like is accepted anywhere a filename is, so
# `topo=Path("a.top")` works.
_Scalar = str | int | float | os.PathLike


def _tokens(key: str, value: Any) -> list[str]:
    """One keyword as gromos++ argv tokens.

    ``True`` is a flag (``@nots``), a sequence is a multi-valued key (``@minwall 1.0 1.5 2.0``,
    ``@traj a.trc b.trc``), anything else is a single value. ``False`` and ``None`` mean "not
    given", so a caller can pass an optional argument through unconditionally.
    """
    flag = f"@{key}"
    if value is None or value is False:
        return []
    if value is True:
        return [flag]
    if isinstance(value, (str, os.PathLike)):
        return [flag, os.fspath(value)]
    if isinstance(value, (list, tuple)):
        return [flag, *(os.fspath(v) if isinstance(v, os.PathLike) else str(v) for v in value)]
    return [flag, str(value)]


def argv_from(**kwargs: Any) -> list[str]:
    """Keywords as a gromos++ argument list.

    >>> argv_from(topo="a.top", pbc="r", minwall=[1.0, 1.5], nots=True)
    ['@topo', 'a.top', '@pbc', 'r', '@minwall', '1.0', '1.5', '@nots']
    """
    argv: list[str] = []
    for key, value in kwargs.items():
        argv.extend(_tokens(key, value))
    return argv


def run(name: str, /, *, _argv: list[str] | None = None, **kwargs: Any) -> str:
    """Run the program `name` and return what it wrote to standard output.

    Keywords are the program's own ``@key`` names (see ``programs()`` for the accepted keys and
    the usage text). Pass ``_argv`` instead to hand the argument list over verbatim — useful for
    reproducing a command line exactly, including ``@f`` argument files.

    Raises :class:`ProgramError` if the program fails; the exception's ``stdout`` attribute holds
    whatever it had written by then (``check_top`` prints its whole report before failing).

    >>> run("red_top", topo="a.top", atoms="1-6")  # doctest: +SKIP
    'TITLE\\n...'
    """
    argv = _argv if _argv is not None else argv_from(**kwargs)
    return _run_program(name, argv)


def programs() -> list[dict[str, Any]]:
    """Every program this build provides.

    Each row has ``name``, ``known`` (the accepted ``@key``s), ``usage`` and ``status``:

    - ``verified`` — a reference test compares our output with gromos++'s
    - ``ported`` — ported from gromos++, no reference test yet
    - ``stand_in`` — same name, different semantics; a utility, not a port
    - ``stub`` — recognised and parsed, does not do the work yet
    - ``no_counterpart`` — no gromos++ program of this name

    >>> {p["name"] for p in programs()} >= {"make_top", "sim_box", "ion"}
    True
    """
    return list(_programs())
