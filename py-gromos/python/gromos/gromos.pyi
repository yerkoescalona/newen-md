"""Type stubs for the gromos native extension module (gromos.abi3.so)."""

from __future__ import annotations

from typing import Any, Literal, Self, final, overload

import numpy as np
import numpy.typing as npt

__all__ = [
    "Algorithm",
    "AlgorithmSequence",
    "Configuration",
    "Energy",
    "EnergyTrajectory",
    "Frame",
    "InputParameters",
    "MissingFeatureError",
    "Plan",
    "PlanError",
    "Recipe",
    "RecipeError",
    "RunError",
    "SchNetPotential",
    "Simulation",
    "System",
    "Term",
    "Topology",
    "Vec3",
    "XtbPotential",
    "__author__",
    "__version__",
    "algorithms",
    "build_info",
    "rdf",
    "resolve_zone_partition",
    "rmsd",
    "terms",
]

# =============================================================================
# Vec3
# =============================================================================

@final
class Vec3:
    def __new__(cls, x: float, y: float, z: float) -> Self: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def z(self) -> float: ...
    def length(self) -> float: ...
    def dot(self, other: Vec3) -> float: ...
    def cross(self, other: Vec3) -> Vec3: ...
    def normalize(self) -> Vec3: ...
    def __repr__(self) -> str: ...
    def __add__(self, value: Vec3, /) -> Vec3: ...
    def __sub__(self, value: Vec3, /) -> Vec3: ...
    def __mul__(self, value: float, /) -> Vec3: ...

# =============================================================================
# Energy
# =============================================================================

@final
class Energy:
    def __new__(cls) -> Self: ...
    @property
    def total(self) -> float: ...
    @property
    def kinetic(self) -> float: ...
    @property
    def potential(self) -> float: ...
    @property
    def bond(self) -> float: ...
    @property
    def angle(self) -> float: ...
    @property
    def dihedral(self) -> float: ...
    @property
    def improper(self) -> float: ...
    @property
    def lj(self) -> float: ...
    @property
    def coulomb(self) -> float: ...
    def __repr__(self) -> str: ...

# =============================================================================
# EnergyTrajectory
# =============================================================================

@final
class EnergyTrajectory:
    """A `.tre`/`.trg` read through the energy library, with the shape it declares.

    The layout is established before any value is believed: the file's own
    self-description if it has one, else the built-in layout for its `ENEVERSION`,
    else the library on trust — `warnings` says which
    (`docs/src/reference/energy-library.md`).
    """

    def __new__(
        cls,
        path: str,
        library: str | None = None,
        properties: list[str] | None = None,
        free_energy: bool = False,
    ) -> Self: ...
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
    def __getitem__(self, name: str, /) -> npt.NDArray[np.float64]:
        """One requested property over every frame."""
        ...

    @property
    def times(self) -> npt.NDArray[np.float64]: ...
    @property
    def steps(self) -> npt.NDArray[np.float64]: ...
    @property
    def warnings(self) -> list[str]: ...
    @property
    def version(self) -> str | None: ...
    @property
    def profile_version(self) -> int: ...
    @property
    def n_energy_groups(self) -> int: ...
    @property
    def n_baths(self) -> int: ...
    @property
    def subblocks(self) -> list[str]: ...
    @property
    def group_pairs(self) -> list[tuple[int, int]]:
        """The 0-based energy-group pair of each row of a `matrix_`-sized subblock."""
        ...

    def table(self, name: str) -> npt.NDArray[np.float64]:
        """One subblock over every frame: `(n_frames, rows, cols)`."""
        ...

    def to_long(self, name: str) -> dict[str, npt.NDArray[Any]]:
        """Tidy columns of one subblock: frame, time, row, col, group_i, group_j, value."""
        ...

    def to_dict(self) -> dict[str, npt.NDArray[np.float64]]: ...
    def last(self, name: str) -> npt.NDArray[np.float64]:
        """The last frame of a subblock, `(rows, cols)`."""
        ...

# =============================================================================
# Frame
# =============================================================================

@final
class Frame:
    def __new__(cls, time: float, step: int) -> Self: ...
    @property
    def time(self) -> float: ...
    @property
    def step(self) -> int: ...
    @property
    def n_atoms(self) -> int: ...
    def positions_array(self) -> npt.NDArray[np.float64]: ...
    def __repr__(self) -> str: ...

# =============================================================================
# Topology
# =============================================================================

@final
class Topology:
    def __new__(cls, topo_file: str) -> Self: ...
    @staticmethod
    def from_file(topo_file: str) -> Topology: ...
    @staticmethod
    def from_string(text: str) -> Topology:
        """Read a topology from the text of a GROMOS topology file — what ``str(topology)``
        writes, this reads."""
        ...
    def write(self, path: str, title: str | None = None) -> None:
        """Write the topology to a GROMOS topology file."""
        ...
    def __str__(self) -> str:
        """The topology as a GROMOS topology file (gromos++ ``OutTopology`` layout)."""
        ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def charge(self) -> int:
        """Total integer charge (e), the sum of the atomic charges rounded."""
        ...
    @property
    def n_solute_atoms(self) -> int: ...
    @property
    def n_solvent_atoms(self) -> int: ...
    @property
    def masses(self) -> npt.NDArray[np.float64]: ...
    @property
    def charges(self) -> npt.NDArray[np.float64]: ...
    def solvate(self, nsm: int) -> None: ...
    def apply_perturbation(self, pttopo_file: str) -> tuple[float, float]:
        """Merge a GROMOS perturbation topology (.ptp) for FEP/TI runs and its MPERTATOM EDS
        end states (what `md @pttopo` does with NTG != 0 / EDS = 1). Call before solvate() and
        before building a Simulation. Returns the solute's net charge (e) in state A and B.
        Raises ValueError when a perturbed bond is not a bond of the topology."""
        ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...

# =============================================================================
# Configuration
# =============================================================================

@final
class Configuration:
    def __new__(cls, conf_file: str) -> Self: ...
    @staticmethod
    def from_file(conf_file: str) -> Configuration: ...
    @staticmethod
    def from_string(text: str) -> Configuration:
        """Read coordinates from the text of a GROMOS96 file. Velocities are not part of the
        POSITION block and come back empty."""
        ...
    def to_g96(self, title: str | None = None) -> str:
        """The configuration as GROMOS96 text, with an optional TITLE."""
        ...
    def __str__(self) -> str:
        """The configuration as a GROMOS96 coordinate file."""
        ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def positions(self) -> npt.NDArray[np.float64]: ...
    @property
    def velocities(self) -> npt.NDArray[np.float64]: ...
    @property
    def box_dimensions(self) -> tuple[float, float, float]: ...
    def __repr__(self) -> str: ...

# =============================================================================
# InputParameters
# =============================================================================

@final
class InputParameters:
    def __new__(cls, input_file: str) -> Self: ...
    @staticmethod
    def from_file(input_file: str) -> InputParameters: ...
    @staticmethod
    def nve(dt: float, steps: int, constraints: str = "none") -> InputParameters: ...
    @staticmethod
    def nvt(
        dt: float, steps: int, temperature: float, constraints: str = "none"
    ) -> InputParameters: ...
    @staticmethod
    def npt(
        dt: float,
        steps: int,
        temperature: float,
        pressure: float,
        constraints: str = "none",
    ) -> InputParameters: ...
    @staticmethod
    def steepest_descent(steps: int) -> InputParameters: ...
    @property
    def dt(self) -> float: ...
    @property
    def nstlim(self) -> int: ...
    @property
    def nsm(self) -> int: ...
    @property
    def cutoff(self) -> float: ...
    @property
    def rcutp(self) -> float: ...
    @property
    def temperature(self) -> float: ...
    @property
    def ntc(self) -> int: ...
    @property
    def constraints(self) -> str: ...
    @property
    def ntb(self) -> int: ...
    @property
    def nsnb(self) -> int: ...
    @property
    def ntwx(self) -> int: ...
    @property
    def ntwe(self) -> int: ...
    def __repr__(self) -> str: ...

_TermKind = Literal["schnet", "xtb"]
"""The term kinds this build knows (``gromos.terms()``)."""

_AlgorithmKind = Literal[
    "remove_com",
    "forcefield",
    "orchestrator",
    "eds",
    "leap_frog_velocity",
    "thermostat",
    "leap_frog_position",
    "shake",
    "settle",
    "lincs",
    "rottrans",
    "steepest_descent",
    "temperature_calculation",
    "pressure_calculation",
    "barostat",
    "energy_calculation",
]
"""The algorithm kinds (``gromos.algorithms()``), in GROMOS step order."""

# =============================================================================
# The run as data (PLAN.md 3.9): Recipe / Term / Algorithm / Plan
# =============================================================================

class RecipeError(ValueError):
    """A recipe, term or algorithm is invalid, or an input cannot be represented as a run."""

class PlanError(ValueError):
    """An algorithm plan violates the GROMOS step-order invariants."""

class MissingFeatureError(RuntimeError):
    """A term needs a cargo feature this build lacks (e.g. ``ml``)."""

class RunError(RuntimeError):
    """The engine failed to initialise or to take a step."""

class ProgramError(RuntimeError):
    """A ``gromos.tools`` program failed. ``stdout`` carries whatever it had written before
    failing — for ``check_top`` that is the report itself."""

    stdout: str

@final
class Term:
    """An additive energy term (``Term("schnet", model=..., cutoff=..., ...)``)."""

    def __new__(cls, kind: _TermKind, **params: Any) -> Self: ...
    @property
    def kind(self) -> str: ...
    @property
    def feature(self) -> str | None:
        """The cargo feature this term needs, or None."""
        ...
    @property
    def available(self) -> bool:
        """Whether this build can run the term."""
        ...
    def to_dict(self) -> dict[str, Any]: ...
    @staticmethod
    def from_dict(d: dict[str, Any]) -> Term: ...
    def __eq__(self, value: object, /) -> bool: ...
    def __repr__(self) -> str: ...

@final
class Algorithm:
    """One fully-resolved entry of a :class:`Plan` (``Algorithm("remove_com", every=100)``)."""

    def __new__(cls, kind: _AlgorithmKind, **params: Any) -> Self: ...
    @property
    def kind(self) -> str: ...
    def to_dict(self) -> dict[str, Any]: ...
    @staticmethod
    def from_dict(d: dict[str, Any]) -> Algorithm: ...
    def __eq__(self, value: object, /) -> bool: ...
    def __repr__(self) -> str: ...

@final
class Plan:
    """The MD step as an ordered, editable list of :class:`Algorithm`; ``recipe.plan(system)``
    returns it validated, ``Simulation(system, recipe, plan=plan)`` runs it. Targets are an
    index or a kind name."""

    def __new__(cls, algorithms: list[Algorithm] | None = None) -> Self: ...
    @property
    def kinds(self) -> list[str]: ...
    def __len__(self) -> int: ...
    def __getitem__(self, key: int | str, /) -> Algorithm: ...
    def __contains__(self, key: object, /) -> bool: ...
    def to_list(self) -> list[Algorithm]: ...
    def insert(self, index: int, algorithm: Algorithm) -> None: ...
    def insert_after(self, target: int | str, algorithm: Algorithm) -> None: ...
    def insert_before(self, target: int | str, algorithm: Algorithm) -> None: ...
    def remove(self, target: int | str) -> Algorithm: ...
    def replace(self, target: int | str, algorithm: Algorithm) -> None: ...
    def validate(self) -> None:
        """Raise :class:`PlanError` if the step order breaks a GROMOS invariant."""
        ...
    def to_json(self) -> str: ...
    @staticmethod
    def from_json(text: str) -> Plan: ...
    def to_dicts(self) -> list[dict[str, Any]]: ...
    def __repr__(self) -> str: ...

@final
class Recipe:
    """The one description of a run — the data of a GROMOS ``.imd`` grouped by concern, plus
    additive ``terms`` and auxiliary ``inputs``. Immutable: every ``with_*``/``update`` call
    returns a new recipe."""

    def __new__(cls, **groups: Any) -> Self: ...
    @staticmethod
    def from_imd(path: str, allow_passthrough: list[str] | None = None) -> Recipe: ...
    @staticmethod
    def from_dict(d: dict[str, Any]) -> Recipe: ...
    @staticmethod
    def from_toml(text: str) -> Recipe: ...
    @staticmethod
    def from_json(text: str) -> Recipe: ...
    @staticmethod
    def from_bundle(path: str, allow_passthrough: list[str] | None = None) -> Recipe: ...
    @staticmethod
    def nve(dt: float, steps: int, constraints: str = "none") -> Recipe: ...
    @staticmethod
    def nvt(dt: float, steps: int, temperature: float, constraints: str = "none") -> Recipe: ...
    @staticmethod
    def npt(
        dt: float,
        steps: int,
        temperature: float,
        pressure: float,
        constraints: str = "none",
    ) -> Recipe: ...
    @staticmethod
    def minimize(steps: int) -> Recipe: ...
    @property
    def version(self) -> int: ...
    @property
    def title(self) -> str: ...
    @property
    def diagnostics(self) -> list[str]:
        """Absent optional blocks (and what that means) and passed-through blocks."""
        ...
    @property
    def terms(self) -> list[Term]: ...
    @property
    def system(self) -> dict[str, Any]: ...
    @property
    def control(self) -> dict[str, Any]: ...
    @property
    def boundary(self) -> dict[str, Any]: ...
    @property
    def forcefield(self) -> dict[str, Any]: ...
    @property
    def constraints(self) -> dict[str, Any]: ...
    @property
    def ensemble(self) -> dict[str, Any]: ...
    @property
    def minimisation(self) -> dict[str, Any]: ...
    @property
    def perturbation(self) -> dict[str, Any]: ...
    @property
    def eds(self) -> dict[str, Any] | None:
        """The ``EDS`` block (``enabled``, ``form``, ``num_states``, ``s``, ``pairs``, ``eir``),
        or ``None`` without one; the end states come from ``inputs["pttopo"]``'s ``MPERTATOM``."""
        ...
    @property
    def aeds(self) -> dict[str, Any] | None:
        """The ``AEDS`` block (``enabled``, ``form``, ``num_states``, ``emax``, ``emin``,
        ``eir``, ``ntiaedss``, ``restremin``, ``bmaxtype``, ``bmax``, ``asteps``, ``bsteps``), or
        ``None`` without one; the end states come from ``inputs["pttopo"]``'s ``MPERTATOM``.
        ``energy_columns``, ``prevalence`` and ``target_state`` (no ``.imd`` block; present only
        when not ``"intended"``) reproduce gromosXX's A-EDS defects (``"gromosxx-post-update"``,
        ``"gromosxx-plain-exp"``, ``"gromosxx-state1-default"``)."""
        ...
    @property
    def outputs(self) -> dict[str, Any]: ...
    @property
    def execution(self) -> dict[str, Any]: ...
    @property
    def inputs(self) -> dict[str, Any]: ...
    @property
    def passthrough(self) -> dict[str, Any]: ...
    def is_minimization(self) -> bool: ...
    def update(self, **groups: Any) -> Recipe:
        """A copy with the given groups deep-merged (``recipe.update(control={"steps": 10})``);
        an unknown group or field raises :class:`RecipeError`."""
        ...
    def with_term(self, term: Term) -> Recipe: ...
    def without_terms(self) -> Recipe: ...
    def with_inputs(
        self,
        pttopo: str | None = None,
        posresspec: str | None = None,
        refpos: str | None = None,
        distrest: str | None = None,
    ) -> Recipe: ...
    def with_execution(self, parallel: str) -> Recipe:
        """``parallel``: ``"auto"`` | ``"serial"`` | ``"parallel"``."""
        ...
    def to_dict(self) -> dict[str, Any]: ...
    def to_toml(self) -> str: ...
    def to_json(self) -> str: ...
    def to_imd(self, n_atoms: int | None = None) -> str:
        """The GROMOS ``.imd`` text gromosXX would run for this recipe."""
        ...
    def save_imd(self, path: str, n_atoms: int | None = None) -> None: ...
    def to_bundle(
        self,
        directory: str,
        system: System,
        topology_path: str,
        configuration_path: str,
    ) -> str:
        """Write ``input.toml`` + ``run.recipe.toml`` + ``run.imd`` into ``directory``;
        returns the ``input.toml`` path."""
        ...
    def plan(self, system: System) -> Plan:
        """Stage 1 of the builder on ``system``: the validated MD step as editable data."""
        ...
    def __eq__(self, value: object, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __reduce__(self) -> tuple[Any, tuple[str]]: ...

def terms() -> list[dict[str, Any]]:
    """The term kinds this build knows: ``[{"kind", "params", "feature", "available"}, ...]``."""
    ...

def algorithms() -> list[dict[str, Any]]:
    """The algorithm kinds and their ordering rules: ``[{"kind", "params", "rules"}, ...]``."""
    ...

def build_info() -> dict[str, Any]:
    """``{"version", "recipe_version", "features"}`` of this extension build."""
    ...

# =============================================================================
# Simulation
# =============================================================================

@final
class Simulation:
    """A running system: ``Simulation(system, recipe)``; ``plan=`` runs an edited
    ``recipe.plan(system)``. The ``InputParameters`` form, the ``distrest``/``posresspec``/
    ``refpos`` and ``ml_*`` keywords, and ``from_sequence`` are deprecated (they warn and are
    translated into a recipe)."""

    @overload
    def __new__(
        cls,
        arg1: System,
        arg2: Recipe | InputParameters,
        arg3: None = None,
        *,
        plan: Plan | None = None,
        distrest: str | None = None,
        posresspec: str | None = None,
        refpos: str | None = None,
        ml_potential: SchNetPotential | None = None,
        ml_region: str | None = None,
        ml_buffer: str | None = None,
    ) -> Self: ...
    @overload
    def __new__(
        cls,
        arg1: str | Topology,
        arg2: str | Configuration,
        arg3: str | Recipe | InputParameters,
        *,
        plan: Plan | None = None,
        distrest: str | None = None,
        posresspec: str | None = None,
        refpos: str | None = None,
        ml_potential: SchNetPotential | None = None,
        ml_region: str | None = None,
        ml_buffer: str | None = None,
    ) -> Self: ...
    @staticmethod
    def from_files(topo_file: str, conf_file: str, input_file: str) -> Simulation: ...
    @staticmethod
    def from_bundle(path: str, allow_passthrough: list[str] | None = None) -> Simulation:
        """Everything from a run bundle (``input.toml``): topology, coordinates, recipe."""
        ...
    @staticmethod
    def from_sequence(
        topo: Topology,
        conf: Configuration,
        params: InputParameters,
        sequence: Plan,
    ) -> Simulation:
        """Deprecated: ``Simulation(System(topo, conf), recipe, plan=sequence)``."""
        ...
    def step(self, n_steps: int) -> None: ...
    def run(self, n_steps: int, ene_freq: int = 100) -> npt.NDArray[np.float64]: ...
    @property
    def time(self) -> float: ...
    @property
    def current_step(self) -> int: ...
    @property
    def dt(self) -> float: ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def algorithm_names(self) -> list[str]: ...
    @property
    def recipe(self) -> Recipe:
        """The effective recipe — what the engine actually ran."""
        ...
    @property
    def plan(self) -> Plan:
        """The algorithm plan the engine instantiated (a frozen snapshot)."""
        ...
    @property
    def recipe_toml(self) -> str:
        """The effective recipe as TOML — the text the `md` binary writes next to its `.tre`."""
        ...
    @property
    def plan_json(self) -> str:
        """The plan as JSON, one fully-resolved entry per algorithm."""
        ...
    @property
    def diagnostics(self) -> list[str]:
        """Which optional parameter blocks were absent (and what that means) and which
        blocks passed through; empty when nothing was defaulted."""
        ...
    @property
    def positions(self) -> npt.NDArray[np.float64]: ...
    @property
    def velocities(self) -> npt.NDArray[np.float64]: ...
    @property
    def forces(self) -> npt.NDArray[np.float64]: ...
    @property
    def energies(self) -> Energy: ...
    @property
    def total_energy(self) -> float: ...
    @property
    def potential_energy(self) -> float: ...
    @property
    def kinetic_energy(self) -> float: ...
    @property
    def precalclam(self) -> dict[str, npt.NDArray[np.float64]] | None:
        """The extended-TI (``PRECALCLAM``) data of the last step, ``None`` without it:
        ``lambdas`` (the grid, N), ``energies`` and ``derivatives`` (N × 12, the ``.tre``'s and
        the ``.trg``'s ``# precalclam`` rows), ``abdih`` and ``free_abdih`` (A, B). Computed
        on the steps NTWG divides only."""
        ...
    @property
    def eds_vr(self) -> float:
        """EDS reference-state energy V_R (kJ/mol) at the last step (``.tre`` slot 35, part of
        ``total_energy`` through the special energy); 0 without EDS."""
        ...
    @property
    def eds_vmix(self) -> float:
        """A-EDS mixed energy V_mix (kJ/mol) at the last step (``.tre`` slot 34, in no total);
        0 without A-EDS."""
        ...
    @property
    def eds_emax(self) -> float:
        """A-EDS E_max (kJ/mol) after the last step's search update (``.tre`` slot 36)."""
        ...
    @property
    def eds_emin(self) -> float:
        """A-EDS E_min (kJ/mol) after the last step's search update (``.tre`` slot 37)."""
        ...
    @property
    def eds_eir(self) -> list[float]:
        """A-EDS offsets E_i^R (kJ/mol) after the last step's search update, one per state
        (the ``.tre``'s ``# eds`` offset column); zeros in plain EDS, empty without EDS."""
        ...
    @property
    def eds_vi(self) -> list[float]:
        """EDS end-state energies V_i (kJ/mol) at the last step, one per state; empty
        without EDS."""
        ...
    @property
    def eds_dvi(self) -> list[float]:
        """∂V_i/∂λ of the EDS end states (kJ/mol) at the last step, one per state: the
        EDS–perturbed pairs of a run with perturbed atoms; empty without EDS."""
        ...
    @property
    def dhdl(self) -> float:
        """dH/dλ (kJ/mol) at the last step — the TI integrand: the whole derivative (kinetic,
        potential, constraints and restraints), the first value of a ``.trg`` frame; 0 unless
        the recipe's ``perturbation`` is enabled (a ``pttopo`` in ``inputs`` and ``NTG=1``)."""
        ...
    @property
    def dhdl_terms(self) -> dict[str, float]:
        """dH/dλ split by term at the last step (kJ/mol): ``lj``, ``crf``, ``bonded``,
        ``kinetic`` (λ-dependent masses), ``constraints`` (perturbed SHAKE),
        ``distanceres`` (perturbed distance restraints) and ``eds`` (∂V_R/∂λ of an EDS run
        with perturbed atoms, ``.trg`` slot 35)."""
        ...
    @property
    def term_energies(self) -> dict[str, float]:
        """Energy of every additive term (``Term("xtb", ...)``) at the last step, keyed by the
        term's registry name (``"xtb"``, ``"xtb:1"`` for a repeated kind), in plan order. Their
        sum is what the terms add to ``total_energy``; ``potential_energy`` is the classical
        force field alone."""
        ...
    @property
    def temperature(self) -> float: ...
    @property
    def volume(self) -> float: ...
    @property
    def pressure(self) -> float: ...
    @property
    def n_solute_atoms(self) -> int: ...
    @property
    def n_solvent_atoms(self) -> int: ...
    def __repr__(self) -> str: ...

# =============================================================================
# AlgorithmSequence — deprecated shim (PLAN.md 3.9 step 4): use recipe.plan(system)
# =============================================================================

@final
class AlgorithmSequence:
    """Deprecated (one release). The descriptor path is gone: each preset returns the
    :class:`Plan` of ``params`` (what ``Recipe.from_imd(...).plan(system)`` builds) and
    warns. Cannot be instantiated."""

    @staticmethod
    def nve(topo: Topology, params: InputParameters) -> Plan: ...
    @staticmethod
    def nvt(topo: Topology, params: InputParameters) -> Plan: ...
    @staticmethod
    def npt(topo: Topology, params: InputParameters) -> Plan: ...
    @staticmethod
    def minimize(topo: Topology, params: InputParameters) -> Plan: ...
    @staticmethod
    def from_parameters(topo: Topology, params: InputParameters) -> Plan: ...

# =============================================================================
# System
# =============================================================================

@final
class System:
    def __new__(cls, topology: Topology, configuration: Configuration) -> Self: ...
    @staticmethod
    def from_files(topo_file: str, conf_file: str) -> System: ...
    def solvate(
        self,
        solvent: Configuration,
        *,
        minwall: list[float] | None = None,
        thresh: float = 0.23,
        boxsize: bool = False,
    ) -> System:
        """Solvate this system in a box of pre-equilibrated solvent (`gromos.tools.sim_box`).

        Returns a new System — the topology follows the coordinates, so both halves agree."""
        ...
    def add_ions(
        self,
        *,
        forcefield: tuple[str, str],
        positive: tuple[int, str] | None = None,
        negative: tuple[int, str] | None = None,
        potential: float | None = None,
        mindist: float = 0.0,
    ) -> System:
        """Replace solvent molecules by counter-ions (`gromos.tools.ion`), topology included.

        `forcefield` is the `(mtb, ifp)` pair `make_top` takes; each ion's topology is built from
        it and combined with the solute's, which gromos++'s `ion` leaves to the caller."""
        ...
    def check(self) -> tuple[list[str], list[str], str]:
        """Run every consistency check over this system's topology; `gromos.tools.check_top`
        wraps the result as a `CheckReport`."""
        ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def charge(self) -> int: ...
    @property
    def positions(self) -> npt.NDArray[np.float64]: ...
    @property
    def velocities(self) -> npt.NDArray[np.float64]: ...
    @property
    def box(self) -> tuple[float, float, float]: ...
    @property
    def topology(self) -> Topology: ...
    @property
    def configuration(self) -> Configuration: ...
    def write(self, path: str) -> None: ...
    def __reduce__(self) -> tuple[Any, tuple[str, str]]:
        """A `System.from_files` system pickles as `from_files` on its two (absolute) paths; one
        built in memory raises `TypeError`."""
        ...
    def __repr__(self) -> str: ...

# Free functions
# =============================================================================

# =============================================================================
# The gromos++-style programs (PLAN.md 3.10) — wrapped by `gromos.tools`
# =============================================================================

def _run_program(name: str, argv: list[str]) -> str:
    """Run the program `name` with an already-`@key`-form argument list; returns its stdout."""
    ...

def _programs() -> list[dict[str, Any]]:
    """Every program this build provides: name, known keys, usage, status."""
    ...

def _sim_box(
    topology: Topology,
    solute: Configuration,
    solvent: Configuration,
    *,
    minwall: list[float] | None = None,
    thresh: float = 0.23,
    boxsize: bool = False,
) -> Configuration: ...
def _ion(
    topology: Topology,
    system: Configuration,
    *,
    positive: tuple[int, str] | None = None,
    negative: tuple[int, str] | None = None,
    potential: float | None = None,
    mindist: float = 0.0,
) -> Configuration: ...
def _build_box(
    topology: Topology, molecule: Configuration, *, nsm: int, dens: float
) -> Configuration: ...
def _com_top(parts: list[tuple[int, Topology]], *, param: int = 1, solv: int = 1) -> Topology: ...
def _red_top(topology: Topology, *, first: int, last: int) -> Topology: ...
def _pdb2g96(
    topology: Topology, pdb: str, *, lib: str | None = None, gch: bool = False
) -> tuple[Configuration, int]: ...
def _make_top(
    *, build: list[str], param: str, seq: list[str], solv: str | None = None
) -> Topology: ...
def _check_top(topology: Topology) -> tuple[list[str], list[str], str]: ...
def _make_pt_top_mutate(
    topology: Topology,
    atoms: str | list[int],
    *,
    charge_b: list[float] | None = None,
    iac_b: list[int] | None = None,
    mass_b: list[float] | None = None,
    alpha_lj: float = 1.51,
    alpha_crf: float = 0.5,
    title: str | None = None,
) -> tuple[str, list[tuple[int, float, float]]]:
    """`make_pt_top` from `topology` to a copy with new state-B parameters on `atoms` (a gromos++
    specifier or 0-based indices; `iac_b` 1-based): the `.ptp` text and each perturbed
    molecule's (0-based index, Q_A, Q_B)."""
    ...

def _ptp_net_charge(topology: Topology, text: str) -> tuple[float, float]:
    """The solute's net charge in A and B with the `.ptp` text merged into a copy of
    `topology` (what `Topology.apply_perturbation` returns)."""
    ...

def _precalclam_longrange_warning(recipe: Recipe) -> str | None:
    """The `PrecalclamLongRangeGaps` warning a run of `recipe` raises, or `None`
    (`gromos.free_energy` warns once for a set of windows with it)."""
    ...

def _pt_top(
    topology: Topology, text: str, *, first_atom: int = 0
) -> tuple[Topology, list[tuple[int, float, float]]]:
    """`pt_top`: `topology` in state B of the `.ptp` text, and each perturbed molecule's
    (0-based index, Q_A, Q_B)."""
    ...

def _ext_ti_ana(argv: list[str], *, seed: int | None = None, write: bool = True) -> dict[str, Any]:
    """`ext_ti_ana` in process: the predictions as arrays and the files written;
    `gromos.analysis.ext_ti_ana` builds `argv` from keywords."""
    ...

def _bar(
    files: list[str],
    temperature: float,
    *,
    max_iterations: int = 500,
    convergence: float = 1e-5,
    bootstrap: int = 0,
    seed: int | None = None,
    printdist: str | None = None,
    parity: str | None = None,
) -> dict[str, Any]:
    """`bar` in process on `ext_ti_ana @bar_data` files: per-interval and total ΔG, kJ/mol."""
    ...

def _ext_ti_merge(argv: list[str]) -> dict[str, Any]:
    """`ext_ti_merge` in process: the merged ⟨∂H/∂λ⟩ curve and its integral."""
    ...

def _mbar(
    u_kn: npt.NDArray[np.float64],
    n_k: list[int],
    temperature: float,
    *,
    max_iterations: int = 10000,
    tolerance: float = 1e-10,
) -> dict[str, Any]:
    """MBAR on potentials in kJ/mol: free energies and their errors, kJ/mol."""
    ...

def _ti_trapezoid(
    files: list[str], lambdas: list[float] | None = None, skip: int = 0
) -> dict[str, Any]:
    """Plain multi-window TI over `.trg` files (trapezoid rule)."""
    ...

def _ti_trapezoid_series(
    lambdas: list[float], series: list[list[float]], skip: int = 0
) -> dict[str, Any]:
    """The same TI over in-memory dH/dλ series, one per λ (`gromos.free_energy`)."""
    ...

def _precalc_bar_data(
    windows: list[dict[str, Any]],
    grid: tuple[int, float, float],
    nlam: float,
    temperature: float,
    skip: int = 0,
) -> list[dict[str, Any]]:
    """The per-frame energies BAR/MBAR use, from in-memory PRECALCLAM frames
    (`gromos.free_energy`): per window `slam`, `plam`, `e_s`, `e_p`, kJ/mol."""
    ...

def _precalc_bar(
    windows: list[dict[str, Any]],
    grid: tuple[int, float, float],
    nlam: float,
    temperature: float,
    skip: int = 0,
    *,
    max_iterations: int = 500,
    convergence: float = 1e-5,
) -> dict[str, Any]:
    """`bar` over in-memory PRECALCLAM frames: what `_bar` returns."""
    ...

def _precalc_mbar(
    windows: list[dict[str, Any]],
    grid: tuple[int, float, float],
    nlam: float,
    temperature: float,
    skip: int = 0,
    *,
    max_iterations: int = 10000,
    tolerance: float = 1e-10,
) -> dict[str, Any]:
    """MBAR over in-memory PRECALCLAM frames: what `_mbar` returns, plus `lambdas`."""
    ...

def _dgslv_pbsolv(argv: list[str]) -> dict[str, Any]:
    """`dGslv_pbsolv` in process: potentials and free energies of every solve, the table and
    the log; `gromos.analysis.dgslv_pbsolv` builds `argv` from keywords."""
    ...

def _charge_correction(
    scheme: str,
    charge_a: float,
    charge_b: float,
    box_edge: float,
    cutoff: float,
    eps_rf: float,
    solvent: str,
    *,
    eps_solvent: float,
    eps_model: float | None,
    ion_radius: float,
    density: float,
    quadrupole_trace: float | None,
    four_pi_eps_i: float,
) -> dict[str, float]:
    """The analytic charging-free-energy correction terms (kJ/mol);
    `gromos.analysis.charge_correction` fills in the defaults."""
    ...

def rmsd(
    positions: npt.NDArray[np.float64],
    reference: npt.NDArray[np.float64],
) -> float: ...
def rdf(
    positions: npt.NDArray[np.float64],
    group1: list[int],
    group2: list[int],
    n_bins: int,
    r_max: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]: ...

# =============================================================================
# ML potential API (only present in a --features ml build — see __init__.py)
# =============================================================================

@final
class SchNetPotential:
    """A trained SchNetPack 2 TorchScript model, to be attached to a
    `Simulation`'s QM/ML zone via `ml_potential=`/`ml_region=`/`ml_buffer=`."""

    def __new__(cls, model_path: str, cutoff: float, elements: list[int]) -> Self: ...
    def evaluate(self, positions: npt.NDArray[np.float64]) -> tuple[float, npt.NDArray[np.float64]]:
        """Real model energy (kJ/mol) and forces (kJ/mol/nm) for an Nx3
        positions array (nm) — a standalone call, not via a Simulation."""
        ...
    def __repr__(self) -> str: ...

@final
class XtbPotential:
    """A real `xtb` (GFN-xTB) QM engine, callable directly for comparison
    against a trained `SchNetPotential`. Always available (no --features ml
    needed — this wraps a subprocess call to the real `xtb` binary)."""

    def __new__(
        cls,
        work_dir: str,
        elements: list[int],
        gfn: int = 2,
        charge: int = 0,
        multiplicity: int = 1,
    ) -> Self: ...
    def evaluate(self, positions: npt.NDArray[np.float64]) -> tuple[float, npt.NDArray[np.float64]]:
        """Real xtb energy (kJ/mol) and forces (kJ/mol/nm) for an Nx3
        positions array (nm) — isolated cluster, vacuum, Embedding::None."""
        ...
    def __repr__(self) -> str: ...

def resolve_zone_partition(
    topology: Topology, inner: str, buffer: str | None = None
) -> tuple[list[int], list[int], list[int]]:
    """Resolve an inner/buffer/outer zone split by selector string against a
    real topology. Returns (inner_indices, buffer_indices, outer_indices)."""
    ...

# Module metadata
__version__: str
__author__: str
