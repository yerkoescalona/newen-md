"""
newen.tools — the gromos++ system-preparation programs, called in process.

Every program runs the same Rust code its command-line binary runs (3.10), so nothing
here shells out, nothing has to find a binary on disk, and the text you get back is byte for byte
what ``sim_box @topo …`` would print.

Two ways to call them:

``run(name, **kwargs)`` — any of the programs, returning its output as text. Keywords are the
gromos++ ``@key`` names:

>>> from newen import tools
>>> text = tools.run(
...     "sim_box", topo="a.top", pbc="r", pos="a.cnf", solvent="spc.cnf", minwall=[1.0]
... )  # doctest: +SKIP

Typed wrappers for the preparation chain — they take and return objects, so a workflow needs no
files in between:

>>> top = tools.make_top(
...     build="54a7.mtb", param="54a7.ifp", seq=["NH3+", "ALA", "GLY", "COO-"], solv="H2O"
... )  # doctest: +SKIP
>>> cnf = tools.pdb2g96(top, pdb="peptide.pdb")  # doctest: +SKIP
>>> cnf = tools.sim_box(top, cnf, solvent="spc.cnf", minwall=[1.0])  # doctest: +SKIP
>>> cnf = tools.ion(top, cnf, positive=(2, "NA+"), potential=0.8)  # doctest: +SKIP
>>> system = newen.System(top, cnf)  # doctest: +SKIP

``programs()`` lists what this build has, including how faithful each is to gromos++'s program of
the same name (``verified`` / ``ported`` / ``stand_in`` / ``stub`` / ``no_counterpart``).
"""

from ._run import ProgramError, argv_from, programs, run
from ._typed import (
    build_box,
    check_top,
    com_top,
    ion,
    make_top,
    pdb2g96,
    red_top,
    sim_box,
)

__all__ = [
    # the general entry point
    "run",
    "programs",
    "argv_from",
    "ProgramError",
    # the preparation chain, typed
    "make_top",
    "com_top",
    "red_top",
    "pdb2g96",
    "sim_box",
    "build_box",
    "ion",
    "check_top",
]
