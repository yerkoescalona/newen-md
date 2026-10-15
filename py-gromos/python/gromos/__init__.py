"""
gromos — Python bindings for the GROMOS-RS molecular dynamics engine.

Quick start
-----------
>>> from gromos import System, Recipe, Simulation
>>>
>>> system = System.from_files("water_216.topo", "equilibrated.cnf")
>>> recipe = Recipe.nvt(dt=0.002, steps=1000, temperature=300.0)  # or Recipe.from_imd("run.imd")
>>> sim = Simulation(system, recipe)
>>> sim.step(1000)
>>> print(sim.total_energy)  # kJ/mol

The recipe is the one description of a run — the same data a GROMOS ``.imd`` file carries,
grouped by concern (``recipe.control``, ``.ensemble``, ``.constraints`` …), plus additive
``terms`` (``Term("xtb", ...)``, ``Term("schnet", ...)``) and auxiliary ``inputs``. ``recipe.plan(system)`` is the MD
step as editable data; ``recipe.to_imd(n_atoms=...)`` is what gromosXX would run.
``gromos.terms()`` / ``gromos.algorithms()`` list what this build knows.

All units follow the GROMOS convention: nm, ps, kJ/mol, K.
"""

# ── Rust extension (always available) ────────────────────────────────────────
from .gromos import (
    Algorithm,
    # Deprecated: its presets return the Plan of the parameters; use recipe.plan(system)
    AlgorithmSequence,
    Configuration,
    Energy,
    # .tre/.trg read through the energy library, shape kept
    EnergyTrajectory,
    Frame,
    InputParameters,  # deprecated: Recipe.from_imd / Recipe.nve|nvt|npt|minimize
    Plan,
    # The run as data (3.9)
    Recipe,
    Simulation,
    # Simulation objects
    System,
    Term,
    # I/O objects
    Topology,
    # Math / utility
    Vec3,
    # QM potential (real xtb — always available, no --features ml needed)
    XtbPotential,
    algorithms,
    build_info,
    rdf,
    rmsd,
    terms,
)

# ── ML potential API (only present if the extension was built --features ml) ─
# `Term("schnet", ...)` always constructs (`.available` says whether this build can run it);
# `Simulation` raises `MissingFeatureError` on a build without the feature.
try:
    from .gromos import SchNetPotential, resolve_zone_partition

    _HAS_ML = True
except ImportError:
    _HAS_ML = False

# ── Exceptions (subclasses of the builtins the binding raised before) ─────────
# ── The gromos++ system-preparation programs, in process ───────
# `tools.run(name, **kwargs)` for any of the 33; typed wrappers for the preparation chain
# (`make_top`, `pdb2g96`, `sim_box`, `ion`, `check_top`, …) that take and return objects.
# ── λ windows of a perturbed run, in process: `free_energy.LambdaSchedule`, `run_windows` ─
from . import exceptions, free_energy, tools
from .exceptions import MissingFeatureError, PlanError, ProgramError, RecipeError, RunError

# `gromos.analysis` (subprocess wrappers around gromos++ programs, mostly stubs) and
# `gromos.system_builder` (design sketch) stay importable explicitly — `import gromos.analysis` —
# but are not part of the default namespace.
# ── Energy timeseries wrapper for Simulation.run() ────────────────────────────
from .timeseries import EnergyTimeseries

try:
    from importlib.metadata import version as _dist_version

    __version__ = _dist_version("gromos")
except Exception:  # not installed as a distribution (e.g. imported from a checkout)
    __version__ = "0.0.0+unknown"

__all__ = [
    # Math / utility
    "Vec3",
    "Energy",
    "Frame",
    "rmsd",
    "rdf",
    # I/O
    "Topology",
    "Configuration",
    "InputParameters",
    # The run as data
    "Recipe",
    "Term",
    "Algorithm",
    "Plan",
    "terms",
    "algorithms",
    "build_info",
    # Exceptions
    "exceptions",
    "RecipeError",
    "PlanError",
    "MissingFeatureError",
    "RunError",
    "ProgramError",
    # The gromos++ programs
    "tools",
    # Simulation
    "System",
    "Simulation",
    "AlgorithmSequence",  # deprecated
    # QM potential
    "XtbPotential",
    # Timeseries
    "EnergyTimeseries",
    "EnergyTrajectory",
    # Free-energy windows
    "free_energy",
    # Meta
    "__version__",
]

if _HAS_ML:
    __all__ += ["SchNetPotential", "resolve_zone_partition"]
