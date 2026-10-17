# PROTOTYPE — design sketch only. Nothing here is implemented or stable.
# Open decisions tracked in py-gromos/docs/development/api-design-mockup.ipynb
# and ("Compositional topology construction in py-gromos").
# Do not use in production code; do not import in tests.

"""
gromos.system_builder — System building algebra for GROMOS.

Design target (3.x, the py-gromos API thread):
    ForceField → BuildingBlock → Topology as an algebra;
    Python expresses verbs; Rust core owns every invariant.

This module is a **design sketch** — the classes below show the intended
public API with full NumPy-style docstrings and runnable stub bodies.
Stubs marked *subprocess* delegate to gromospp binaries today; stubs
marked *native* will eventually call into pyo3-gromos.

Motivation (from vsomm_modeler experience)
-------------------------------------------
The old pattern looked like this in every method::

    input_make_top = {}
    input_make_top['make_top'] = gromospp_bin + "make_top"
    input_make_top['mtb']  = vsomm_dir + forcefield + "/*.mtb"
    input_make_top['seq']  = " ".join(mol)
    input_make_top['output'] = workdir + "/mol_0.top"
    command = "{make_top} @build {mtb} @param {ifp} @seq {seq} > {output}"
    subprocess.check_output(command, shell=True, ...)

Resulting pain:
- Atom counts, charges, water counts tracked manually → drift.
- Binary paths threaded through every method as plain strings.
- Topology and coordinates always out of step (built separately, matched by file name).
- No type boundary — a `.top` path and a `.cnf` path are both just `str`.
- Workdir cleanup forgotten → temp files everywhere.

The new algebra solves each of these:

    ff  = ForceField("54a7")                          # binary paths: looked up once
    mol = ff.molecule(["STA", "ALA", "GLY", "END"])   # topology + atoms owned here
    sys = mol * 10 + ff.solvent("SPC") * 216          # algebra
    conf = sys.build_coordinates(seed=42)             # geometry follows topology
    result = sim.minimize(conf).equilibrate().run()   # pipeline, not dicts

Alternatives considered
-----------------------
Three API shapes are documented in Notes sections throughout this file:

- **Shape A — OpenMM-style composition**: ``forcefield.build_topology(sequence)``.
  Clean for single-protein MD; awkward when topology must be built from many
  independently minimised molecules.

- **Shape B — Builder / fluent chaining**: ``SystemBuilder().add(...).solvate(...).build()``.
  Explicit pipeline; good for Jupyter; makes it obvious that each step is lazy.

- **Shape C — Declarative dict / dataclass**: ``System(molecules=[...], solvent=...)``.
  Readable for config files; loses the algebra, makes ``mol * 10`` impossible.

The module implements **Shape A** with the ``*`` / ``+`` algebra from
**Shape B**, and every public class can be serialised to a dict for **Shape C**
round-trip. See each class's ``Notes`` for per-decision discussion.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

    # `MDResult.energies`: the package's one energy timeseries (`Simulation.run()`'s columns)
    from .timeseries import EnergyTimeseries

# ---------------------------------------------------------------------------
# Internal helpers — binary discovery
# ---------------------------------------------------------------------------

_GROMOS_ENV_XX = "GROMOSXX_BIN"
_GROMOS_ENV_PP = "GROMOSPP_BIN"


def _find_bin(name: str, env_var: str) -> Path:
    """Return the path to a GROMOS binary.

    Search order: environment variable *env_var* → ``target/release`` of the
    gromos-rs workspace → ``PATH``.

    Parameters
    ----------
    name : str
        Binary name, e.g. ``"make_top"``, ``"md"``.
    env_var : str
        Environment variable that overrides the search, e.g. ``"GROMOSPP_BIN"``.

    Returns
    -------
    Path
        Resolved path to the binary.

    Raises
    ------
    FileNotFoundError
        If the binary cannot be found.
    """
    env_dir = os.environ.get(env_var)
    if env_dir:
        candidate = Path(env_dir) / name
        if candidate.exists():
            return candidate

    # gromos-rs workspace layout: py-gromos is one level inside the workspace
    workspace_root = Path(__file__).parent.parent.parent.parent
    for subdir in ("target/release", "target/debug"):
        candidate = workspace_root / subdir / name
        if candidate.exists():
            return candidate

    found = shutil.which(name)
    if found:
        return Path(found)

    raise FileNotFoundError(
        f"Cannot find '{name}'. Set {env_var} or build gromos-rs with 'cargo build --release'."
    )


def _run(cmd: str, workdir: Path) -> str:
    """Run a shell command, log it, and return stdout.

    Parameters
    ----------
    cmd : str
        Full shell command string.
    workdir : Path
        Working directory; ``log.txt`` is appended here.

    Returns
    -------
    str
        Combined stdout output.

    Raises
    ------
    RuntimeError
        On non-zero exit with the captured stderr/stdout.
    """
    log = workdir / "log.txt"
    with log.open("a") as f:
        f.write(f"\n[+] {cmd}\n")
    try:
        out = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT, universal_newlines=True
        )
    except subprocess.CalledProcessError as exc:
        with log.open("a") as f:
            f.write(exc.output)
        raise RuntimeError(f"Command failed:\n{cmd}\n\n{exc.output}") from exc
    with log.open("a") as f:
        f.write(out)
    return out


# ---------------------------------------------------------------------------
# ForceField
# ---------------------------------------------------------------------------


class ForceField:
    """GROMOS force field — the root of the building algebra.

    Knows the ``.ifp`` (interaction function parameters) and one or more
    ``.mtb`` (molecular topology building blocks) files.  Every
    :class:`BuildingBlock` and every :class:`MoleculeTopology` is created
    through this object so that the parameter set is always consistent.

    Parameters
    ----------
    name_or_ifp : str or Path
        Either a named force field (e.g. ``"54a7"``, ``"2016h66"``) looked up
        from ``GROMOS_FF_DIR`` or the path to an existing ``.ifp`` file.
    mtb : str or Path or list of (str or Path), optional
        Path(s) to ``.mtb`` file(s).  If *name_or_ifp* is a named force
        field the mtb is resolved automatically from the same directory.
    ff_dir : str or Path, optional
        Base directory for named force fields.  Falls back to the
        ``GROMOS_FF_DIR`` environment variable, then to
        ``gromos-rs/.local/forcefields/``.

    Attributes
    ----------
    name : str
        Force-field identifier string (e.g. ``"54a7"``).
    ifp : Path
        Resolved path to the ``.ifp`` file.
    mtb_files : list of Path
        Resolved paths to all ``.mtb`` files.

    Notes
    -----
    **Shape A vs Shape C**: In Shape C you would write
    ``{"forcefield": "54a7"}`` and let the system resolve it later.  Here we
    resolve eagerly so that a typo in the force-field name fails at
    construction, not buried inside a pipeline step.

    **Binary path management**: The :class:`ForceField` object stores the
    resolved paths to ``make_top``, ``com_top``, ``sim_box``, ``ran_box``,
    ``ion``, ``gch`` once at construction time.  Every downstream object
    receives a reference to this :class:`ForceField` rather than path strings.

    Examples
    --------
    Load by name (requires ``GROMOS_FF_DIR`` or the standard workspace layout):

    >>> ff = ForceField("54a7")
    >>> ff.name
    '54a7'

    Load from explicit paths:

    >>> ff = ForceField("54a7.ifp", mtb="54a7.mtb")

    Inspect available building blocks:

    >>> ff.building_blocks()  # doctest: +SKIP
    ['ALA', 'GLY', 'PRO', 'SER', ..., 'SPC', 'NA+', 'CL-']
    """

    def __init__(
        self,
        name_or_ifp: str | Path,
        mtb: None | str | Path | list = None,
        ff_dir: None | str | Path = None,
    ) -> None:
        p = Path(name_or_ifp)
        if p.suffix == ".ifp" and p.exists():
            self.name = p.stem
            self.ifp = p
            self.mtb_files = (
                [Path(m) for m in (mtb if isinstance(mtb, list) else [mtb])] if mtb else []
            )
        else:
            self.name = str(name_or_ifp)
            self.ifp, self.mtb_files = self._resolve_named(self.name, ff_dir, mtb)

        # Resolve gromospp binaries once
        self._make_top = _find_bin("make_top", _GROMOS_ENV_PP)
        self._com_top = _find_bin("com_top", _GROMOS_ENV_PP)
        self._sim_box = _find_bin("sim_box", _GROMOS_ENV_PP)
        self._ran_box = _find_bin("ran_box", _GROMOS_ENV_PP)
        self._ion_bin = _find_bin("ion", _GROMOS_ENV_PP)
        self._gch = _find_bin("gch", _GROMOS_ENV_PP)
        self._pdb2g96 = _find_bin("pdb2g96", _GROMOS_ENV_PP)
        self._md = _find_bin("md", _GROMOS_ENV_XX)

    # ------------------------------------------------------------------
    # Building-block algebra
    # ------------------------------------------------------------------

    def molecule(self, sequence: Sequence[str]) -> MoleculeTopology:
        """Create a molecule topology from a building-block sequence.

        Calls ``make_top`` internally.  The result owns its own ``.top`` file
        inside a managed temporary directory.

        Parameters
        ----------
        sequence : sequence of str
            Ordered building-block codes as they appear in the ``.mtb``.
            Must include start and end groups (e.g. ``["STA", "ALA", "END"]``).

        Returns
        -------
        MoleculeTopology
            Topology for a single molecule, including atom count and charge.

        Notes
        -----
        **Alternative**: ``ff["ALA"]`` (``__getitem__``) returns a
        :class:`BuildingBlock` object.  Molecules can then be composed with
        ``+``::

            mol = ff["STA"] + ff["ALA"] + ff["GLY"] + ff["END"]

        Both approaches call the same ``make_top`` backend.  The sequence API
        is simpler for static sequences; the ``+`` algebra is cleaner when
        building programmatically from lists.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> ala_gly = ff.molecule(["STA", "ALA", "GLY", "END"])
        >>> ala_gly.n_atoms
        21
        >>> ala_gly.charge
        0
        """
        return MoleculeTopology._from_ff(self, sequence)

    def solvent(self, name: str = "SPC") -> MoleculeTopology:
        """Return the solvent building block as a MoleculeTopology.

        Parameters
        ----------
        name : str, default ``"SPC"``
            Solvent code in the ``.mtb`` (``"SPC"``, ``"TIP3P"``, …).

        Returns
        -------
        MoleculeTopology
            Single-residue solvent topology.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> spc = ff.solvent("SPC")
        >>> spc.n_atoms
        3
        """
        return self.molecule([name])

    def building_blocks(self) -> list[str]:
        """Return the list of building-block codes in all loaded .mtb files.

        Returns
        -------
        list of str
            All residue codes (e.g. ``['ALA', 'GLY', 'SPC', 'NA+', …]``).
        """
        # subprocess: grep MTBUILDBLSOLUTE / MTBUILDBLSOLVENT from mtb
        codes: list[str] = []
        for mtb in self.mtb_files:
            with mtb.open() as f:
                for line in f:
                    stripped = line.strip()
                    if stripped.startswith("MTBUILDBLSOLUTE") or stripped.startswith(
                        "MTBUILDBLSOLVENT"
                    ):
                        pass  # next non-comment line is the code
        return codes  # TODO: implement parser

    def __getitem__(self, code: str) -> BuildingBlock:
        """Return a single building block by code.

        Parameters
        ----------
        code : str
            Building-block code, e.g. ``"ALA"``, ``"SPC"``, ``"NA+"``.

        Returns
        -------
        BuildingBlock
            The building block descriptor.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> ff["ALA"].n_atoms
        10
        >>> ff["ALA"].charge
        0
        """
        return BuildingBlock(code=code, ff=self)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _mtb_arg(self) -> str:
        return " ".join(str(m) for m in self.mtb_files)

    def _resolve_named(
        self,
        name: str,
        ff_dir: str | Path | None,
        extra_mtb,
    ) -> tuple[Path, list[Path]]:
        candidates = []
        if ff_dir:
            candidates.append(Path(ff_dir))
        env = os.environ.get("GROMOS_FF_DIR")
        if env:
            candidates.append(Path(env))
        workspace = Path(__file__).parent.parent.parent.parent / ".local" / "forcefields"
        candidates += [
            workspace / name,
            workspace / "official",
            workspace / "gromos96",
        ]

        ifp: Path | None = None
        mtb_files: list[Path] = []

        for base in candidates:
            candidate_ifp = base / f"{name}.ifp"
            if candidate_ifp.exists():
                ifp = candidate_ifp
                mtb_files = sorted(base.glob(f"{name}*.mtb"))
                break

        if ifp is None:
            raise FileNotFoundError(
                f"Force field '{name}' not found. "
                f"Set GROMOS_FF_DIR or provide explicit ifp/mtb paths."
            )

        if extra_mtb:
            extras = [extra_mtb] if not isinstance(extra_mtb, list) else extra_mtb
            mtb_files += [Path(m) for m in extras]

        return ifp, mtb_files

    def __repr__(self) -> str:
        return f"ForceField(name={self.name!r}, n_mtb={len(self.mtb_files)})"


# ---------------------------------------------------------------------------
# BuildingBlock
# ---------------------------------------------------------------------------


@dataclass
class BuildingBlock:
    """A single building-block entry from a GROMOS .mtb file.

    Not constructed directly — use :meth:`ForceField.__getitem__` or
    :meth:`ForceField.building_blocks`.

    Parameters
    ----------
    code : str
        Residue code (e.g. ``"ALA"``).
    ff : ForceField
        The force field this building block belongs to.

    Attributes
    ----------
    code : str
    n_atoms : int or None
        Atom count (populated lazily by parsing the .mtb).
    charge : int or None
        Integer charge (populated lazily).
    atom_count : int or None
        Alias for ``n_atoms``.

    Notes
    -----
    **Design decision — lazy vs eager parsing**: The .mtb files can be large
    (>10 MB for 2016H66).  Parsing every block at ``ForceField.__init__`` time
    would slow down import.  Instead, per-block attributes are populated on
    first access.  A full parse is only triggered by ``ff.building_blocks()``.

    Examples
    --------
    >>> ff = ForceField("54a7")  # doctest: +SKIP
    >>> bb = ff["GLY"]
    >>> bb.code
    'GLY'
    >>> bb.n_atoms  # triggers lazy mtb parse   # doctest: +SKIP
    7
    """

    code: str
    ff: ForceField
    n_atoms: int | None = field(default=None, repr=False)
    charge: int | None = field(default=None, repr=False)

    @property
    def atom_count(self) -> int | None:
        return self.n_atoms

    def __add__(self, other: BuildingBlock) -> MoleculeTopology:
        """Concatenate two building blocks into the start of a molecule.

        Parameters
        ----------
        other : BuildingBlock
            The next block in sequence.

        Returns
        -------
        MoleculeTopology
            A partial sequence that can be extended further with ``+``.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> mol = ff["STA"] + ff["ALA"] + ff["END"]
        >>> mol.sequence
        ['STA', 'ALA', 'END']
        """
        return MoleculeTopology._from_sequence(self.ff, [self.code, other.code])

    def __repr__(self) -> str:
        return f"BuildingBlock(code={self.code!r})"


# ---------------------------------------------------------------------------
# MoleculeTopology
# ---------------------------------------------------------------------------


class MoleculeTopology:
    """Topology of a single molecule built from a building-block sequence.

    Wraps a GROMOS ``.top`` file produced by ``make_top``.  Knows its own
    atom count, charge, and sequence.  Supports the ``*`` / ``+`` algebra
    for composing :class:`SystemTopology` objects.

    Parameters
    ----------
    top_file : Path
        Path to the ``.top`` file produced by ``make_top``.
    sequence : list of str
        Building-block codes used to create this topology.
    ff : ForceField
        Parent force field.
    _workdir : Path, optional
        Temporary directory that owns *top_file*.  If provided, this object
        takes ownership and may clean it up.

    Attributes
    ----------
    top_file : Path
    sequence : list of str
    ff : ForceField
    n_atoms : int
        Total atom count (read from the ``.top`` header).
    charge : int
        Net charge in elementary charges.

    Notes
    -----
    **Immutability**: :class:`MoleculeTopology` is intentionally immutable
    after construction.  ``mol * 10`` does not copy the ``.top`` file — it
    just records that the topology should be replicated 10 times when
    :class:`SystemTopology` is built.  This avoids the vsomm_modeler pattern
    where a dozen ``.top`` files had to be cleaned up manually.

    **Shape A vs Shape B**: In Shape B you would write
    ``builder.add_molecule(["ALA", "GLY"], n=10)``.  Here the algebra is
    ``ff.molecule(["STA", "ALA", "GLY", "END"]) * 10``.  Both express the
    same intent; the ``*`` form composes with ``+`` more naturally.

    Examples
    --------
    >>> ff = ForceField("54a7")  # doctest: +SKIP
    >>> ala = ff.molecule(["STA", "ALA", "END"])
    >>> gly = ff.molecule(["STA", "GLY", "END"])

    Algebra:

    >>> system_top = ala * 5 + gly * 5 + ff.solvent("SPC") * 216
    >>> isinstance(system_top, SystemTopology)
    True
    """

    def __init__(
        self,
        top_file: Path,
        sequence: list[str],
        ff: ForceField,
        _workdir: Path | None = None,
    ) -> None:
        self.top_file = top_file
        self.sequence = sequence
        self.ff = ff
        self._workdir = _workdir
        self.n_atoms: int = self._read_n_atoms()
        self.charge: int = self._read_charge()

    # ------------------------------------------------------------------
    # Algebra
    # ------------------------------------------------------------------

    def __mul__(self, n: int) -> SystemTopology:
        """Replicate this molecule *n* times.

        Parameters
        ----------
        n : int
            Number of copies.

        Returns
        -------
        SystemTopology
            A system consisting of *n* copies of this molecule.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> spc = ff.solvent("SPC")
        >>> water_box = spc * 216
        """
        return SystemTopology(molecules=[(self, n)], ff=self.ff)

    def __add__(self, other: MoleculeTopology | SystemTopology) -> SystemTopology:
        """Add another molecule or system topology.

        Parameters
        ----------
        other : MoleculeTopology or SystemTopology

        Returns
        -------
        SystemTopology

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> protein = ff.molecule(["STA", "ALA", "GLY", "END"])
        >>> solvent = ff.solvent("SPC")
        >>> system = protein * 10 + solvent * 500
        """
        left = SystemTopology(molecules=[(self, 1)], ff=self.ff)
        return left + other

    # ------------------------------------------------------------------
    # Coordinate building
    # ------------------------------------------------------------------

    def build_coordinates(self, seed: int = 42) -> Configuration:
        """Generate initial coordinates for this molecule.

        Calls ``pdb2g96`` → ``gca`` → ``gch`` in a managed temp directory.

        Parameters
        ----------
        seed : int, default 42
            Random seed for the random-rotation step.

        Returns
        -------
        Configuration
            The coordinates for this molecule.

        Notes
        -----
        **Alternative** (:meth:`SystemTopology.build_coordinates`): For
        multi-molecule systems, coordinates are built per molecule then
        assembled by :meth:`SystemTopology.assemble`.  Only call this method
        directly for a single isolated molecule.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> mol = ff.molecule(["STA", "ALA", "END"])
        >>> conf = mol.build_coordinates(seed=123)
        >>> conf.n_atoms
        12
        """
        raise NotImplementedError("pdb2g96 → gca → gch pipeline — to be implemented")

    # ------------------------------------------------------------------
    # Minimisation
    # ------------------------------------------------------------------

    def minimize(
        self,
        steps: int = 50_000,
        conf: Configuration | None = None,
        seed: int = 42,
    ) -> Configuration:
        """Energy-minimise this molecule.

        Builds coordinates if *conf* is ``None``, then runs steepest-descent
        minimisation using the gromos-rs ``md`` binary.

        Parameters
        ----------
        steps : int, default 50 000
            Maximum steepest-descent steps.
        conf : Configuration, optional
            Starting configuration.  If ``None``, built from scratch.
        seed : int, default 42
            Random seed for coordinate generation.

        Returns
        -------
        Configuration
            Minimised configuration.

        Notes
        -----
        **Pipeline design**: Rather than exposing raw .imd strings as
        vsomm_modeler did (``imds.minimize_molecule.format(...)``), the
        minimisation input is constructed from the topology's atom count so
        the user never sees the IMD format.  To access the raw IMD, call
        ``result.imd_file``.

        The three-stage minimisation sequence from vsomm_modeler::

            (["", "min_", 50000, 0, 0],)  # no charges, no SHAKE
            (["min_", "min2_", 50000, 1, 0],)  # charges, no SHAKE
            ["min2_", "min3_", 50000, 1, 3]  # charges + SHAKE

        is reproduced by calling ``minimize(steps=50000)`` — it internally
        runs all three stages.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> mol = ff.molecule(["STA", "ALA", "END"])
        >>> min_conf = mol.minimize(steps=50_000)
        >>> min_conf.potential_energy  # kJ/mol   # doctest: +SKIP
        -312.7
        """
        if conf is None:
            conf = self.build_coordinates(seed=seed)
        raise NotImplementedError("Three-stage minimisation — to be implemented")

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @classmethod
    def _from_ff(cls, ff: ForceField, sequence: Sequence[str]) -> MoleculeTopology:
        workdir = Path(tempfile.mkdtemp(prefix="gromos_mol_"))
        top_out = workdir / "mol.top"
        mtb_arg = ff._mtb_arg()
        seq_arg = " ".join(sequence)
        cmd = f"{ff._make_top} @build {mtb_arg} @param {ff.ifp} @seq {seq_arg} > {top_out}"
        _run(cmd, workdir)
        return cls(top_file=top_out, sequence=list(sequence), ff=ff, _workdir=workdir)

    @classmethod
    def _from_sequence(cls, ff: ForceField, codes: list[str]) -> MoleculeTopology:
        return cls._from_ff(ff, codes)

    def _read_n_atoms(self) -> int:
        if not self.top_file.exists():
            return 0
        # parse SOLUTECOMP / SOLVENTCOMP block atom counts from .top
        # TODO: call gromos-io parser via pyo3
        return 0

    def _read_charge(self) -> int:
        if not self.top_file.exists():
            return 0
        # TODO: call gromos-io parser via pyo3
        return 0

    def __repr__(self) -> str:
        return (
            f"MoleculeTopology(sequence={self.sequence!r}, "
            f"n_atoms={self.n_atoms}, charge={self.charge:+d})"
        )


# ---------------------------------------------------------------------------
# SystemTopology
# ---------------------------------------------------------------------------


class SystemTopology:
    """Combined topology for a multi-molecule, optionally solvated system.

    Built by the ``+`` / ``*`` algebra on :class:`MoleculeTopology` objects.
    Wraps ``com_top`` for the combined ``.top`` file.

    Parameters
    ----------
    molecules : list of (MoleculeTopology, int)
        Each entry is ``(mol_topology, count)``.
    ff : ForceField
        Parent force field.
    solvent : MoleculeTopology, optional
        Solvent molecule topology.
    n_solvent : int, default 0
        Number of solvent molecules.
    ions : list of (str, int), optional
        Counter-ion specifications as ``[(ion_code, count), …]``.

    Attributes
    ----------
    n_atoms : int
        Total atom count including solvent and ions.
    charge : int
        Net system charge.
    molecules : list of (MoleculeTopology, int)
    ff : ForceField

    Notes
    -----
    **No workdir argument**: Unlike vsomm_modeler, ``workdir`` is not a
    constructor parameter.  It is created lazily when
    :meth:`build_coordinates` or :meth:`assemble` is called, and cleaned up
    automatically unless ``debug=True`` is passed to those methods.

    **Solvent / ion algebra alternatives**:

    Shape A (explicit solvation call)::

        system = (protein * 10).__add__(solvent * 500)
        solvated = system.add_solvent("SPC", n=2000, density=900)

    Shape B (constructor keyword)::

        SystemTopology(molecules=[(protein, 10)], solvent=spc, n_solvent=2000)

    This class supports both:  the ``+`` operator and the
    :meth:`solvate` method produce the same result.

    Examples
    --------
    >>> ff = ForceField("54a7")  # doctest: +SKIP
    >>> protein = ff.molecule(["STA", "ALA", "GLY", "END"])
    >>> spc = ff.solvent("SPC")

    Algebra:

    >>> system = protein * 5 + spc * 216
    >>> system.n_atoms
    0

    Explicit solvation:

    >>> system = (protein * 5).solvate("SPC", n=216, density=900, seed=42)
    """

    def __init__(
        self,
        molecules: list[tuple[MoleculeTopology, int]],
        ff: ForceField,
        solvent: MoleculeTopology | None = None,
        n_solvent: int = 0,
        ions: list[tuple[str, int]] | None = None,
    ) -> None:
        self.molecules = molecules
        self.ff = ff
        self._solvent = solvent
        self._n_solvent = n_solvent
        self._ions: list[tuple[str, int]] = ions or []
        self._top_file: Path | None = None
        self._workdir: Path | None = None

    # ------------------------------------------------------------------
    # Algebra
    # ------------------------------------------------------------------

    def __add__(self, other: MoleculeTopology | SystemTopology) -> SystemTopology:
        """Combine two topologies.

        Parameters
        ----------
        other : MoleculeTopology or SystemTopology

        Returns
        -------
        SystemTopology
        """
        extra = [(other, 1)] if isinstance(other, MoleculeTopology) else other.molecules
        return SystemTopology(molecules=self.molecules + extra, ff=self.ff)

    def __mul__(self, n: int) -> SystemTopology:
        """Scale all molecule counts by *n*.

        Parameters
        ----------
        n : int

        Returns
        -------
        SystemTopology
        """
        return SystemTopology(molecules=[(m, k * n) for m, k in self.molecules], ff=self.ff)

    # ------------------------------------------------------------------
    # Solvation and ionisation
    # ------------------------------------------------------------------

    def solvate(
        self,
        solvent: str = "SPC",
        n: int | None = None,
        density: float = 900.0,
        boxsize: tuple[float, float, float] | None = None,
        seed: int = 42,
        debug: bool = False,
    ) -> SolvatedSystem:
        """Add solvent molecules to the system.

        Calls ``ran_box`` or ``sim_box`` internally (depending on whether
        *boxsize* or *density* is given).

        Parameters
        ----------
        solvent : str, default ``"SPC"``
            Solvent building-block code.
        n : int, optional
            Number of solvent molecules.  If ``None`` and *density* is given,
            the count is determined by ``ran_box``.
        density : float, default 900.0
            Initial density in kg m⁻³.  Ignored when *boxsize* is given.
        boxsize : tuple of (float, float, float), optional
            Box edge lengths in nm ``(a, b, c)``.  If given, *density* is
            ignored and ``sim_box`` is used.
        seed : int, default 42
            Random seed for ``ran_box``.
        debug : bool, default False
            If ``True``, preserve all intermediate files.

        Returns
        -------
        SolvatedSystem
            Solvated topology + coordinates ready for minimisation.

        Notes
        -----
        **ran_box vs sim_box**:  In the original GROMOS pipeline, ``ran_box``
        places molecules at random positions at a target density, while
        ``sim_box`` fills a pre-defined box with solvent.  This method selects
        the appropriate binary automatically based on the arguments.

        The vsomm_modeler code required the user to track *water_molecules*,
        *counterions*, and *atom counts* manually to build the correct .imd
        file.  Here, :class:`SolvatedSystem` owns those counts.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> protein = ff.molecule(["STA", "ALA", "END"])
        >>> solvated = (protein * 5).solvate("SPC", density=900, seed=42)
        >>> solvated.n_solvent
        0
        """
        spc_mol = self.ff.solvent(solvent)
        return SolvatedSystem(
            topology=self,
            solvent=spc_mol,
            n_solvent=n or 0,
            density=density,
            boxsize=boxsize,
            seed=seed,
            debug=debug,
        )

    def neutralize(
        self,
        ion: str = "Na+",
        seed: int = 42,
    ) -> SystemTopology:
        """Add counter-ions to neutralise the system charge.

        Parameters
        ----------
        ion : str, default ``"Na+"``
            Counter-ion building-block code (``"Na+"``, ``"Cl-"``,
            ``"Ca2+"``, …).
        seed : int, default 42
            Random seed for ion placement.

        Returns
        -------
        SystemTopology
            New topology with ions added.

        Notes
        -----
        The valence-aware charge neutralisation logic that vsomm_modeler
        embedded in ``set_counterions()`` (including the special-case for
        Ca²⁺ with odd system charges) is encapsulated here.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> charged_system = ff.molecule(["STA", "ASP", "END"]) * 5
        >>> neutral = charged_system.neutralize("Na+")
        >>> neutral.charge
        0
        """
        valence = {"Na+": 1, "Ca2+": 2, "Cl-": -1, "Mg2+": 2}
        v = valence.get(ion, 1)
        n_ions = abs(self.charge) // v
        new_ions = self._ions + [(ion, n_ions)]
        return SystemTopology(molecules=self.molecules, ff=self.ff, ions=new_ions)

    # ------------------------------------------------------------------
    # System-level topology file
    # ------------------------------------------------------------------

    def top_file(self, workdir: Path | None = None) -> Path:
        """Build and return the combined ``.top`` file (calls ``com_top``).

        Parameters
        ----------
        workdir : Path, optional
            Directory in which to write the file.  Created as a temp dir if
            ``None``.

        Returns
        -------
        Path
            Path to ``system.top``.
        """
        if self._top_file and self._top_file.exists():
            return self._top_file
        if workdir is None:
            workdir = Path(tempfile.mkdtemp(prefix="gromos_sys_"))
            self._workdir = workdir
        out = workdir / "system.top"
        topo_arg = " ".join(" ".join([str(mol.top_file)] * n) for mol, n in self.molecules)
        cmd = f"{self.ff._com_top} @topo {topo_arg} @param 1 @solv 1 > {out}"
        _run(cmd, workdir)
        self._top_file = out
        return out

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def n_atoms(self) -> int:
        """Total atom count (solute only; excludes solvent unless solvated)."""
        return sum(mol.n_atoms * n for mol, n in self.molecules)

    @property
    def charge(self) -> int:
        """Net system charge."""
        return sum(mol.charge * n for mol, n in self.molecules)

    def __repr__(self) -> str:
        mol_summary = ", ".join(
            f"{mol.sequence[0] if mol.sequence else '?'}×{n}" for mol, n in self.molecules
        )
        return f"SystemTopology([{mol_summary}], n_atoms={self.n_atoms}, charge={self.charge:+d})"


# ---------------------------------------------------------------------------
# SolvatedSystem
# ---------------------------------------------------------------------------


class SolvatedSystem:
    """A fully solvated (and optionally ionised) system ready for MD.

    Produced by :meth:`SystemTopology.solvate`.  Knows the combined
    topology, initial coordinates, solvent count, and ion count.  All values
    needed to write a valid ``.imd`` file are computed automatically from the
    topology.

    Parameters
    ----------
    topology : SystemTopology
    solvent : MoleculeTopology
    n_solvent : int
    density : float
        Target density in kg m⁻³.
    boxsize : tuple of float, optional
        Box edge lengths in nm; overrides *density*.
    seed : int
    debug : bool

    Attributes
    ----------
    n_solute_atoms : int
        Atom count excluding solvent and ions.
    n_solvent : int
        Number of solvent molecules.
    n_ions : int
        Total number of counter-ions.
    n_atoms : int
        Total atom count.
    top_file : Path
    conf_file : Path

    Notes
    -----
    **IMD generation**: The vsomm_modeler code required the user to know
    ``mol``, ``watermol``, ``bbsatom``, ``cationatom``, ``lastatom`` — five
    interrelated integers that had to be updated whenever the system changed.
    Here they are computed as properties::

        self.n_solute_atoms   → bbsatom
        self.n_ions           → counted from self.topology._ions
        self.n_solute_atoms + self.n_ions  → cationatom
        self.n_atoms          → lastatom

    and passed automatically to any MD runner.

    Examples
    --------
    >>> ff = ForceField("54a7")  # doctest: +SKIP
    >>> protein = ff.molecule(["STA", "ALA", "END"])
    >>> solvated = (protein * 5).solvate("SPC", density=900, seed=42)
    >>> solvated.n_atoms
    0
    >>> runner = solvated.minimize()
    >>> runner.potential_energy  # doctest: +SKIP
    -12345.6
    """

    def __init__(
        self,
        topology: SystemTopology,
        solvent: MoleculeTopology,
        n_solvent: int,
        density: float,
        boxsize: tuple[float, float, float] | None,
        seed: int,
        debug: bool,
    ) -> None:
        self.topology = topology
        self._solvent = solvent
        self.n_solvent = n_solvent
        self.density = density
        self.boxsize = boxsize
        self.seed = seed
        self.debug = debug
        self._workdir: Path | None = None
        self._top_file: Path | None = None
        self._conf_file: Path | None = None

    @property
    def n_solute_atoms(self) -> int:
        return self.topology.n_atoms

    @property
    def n_ions(self) -> int:
        return sum(n for _, n in self.topology._ions)

    @property
    def n_atoms(self) -> int:
        return self.n_solute_atoms + self.n_ions + self.n_solvent * self._solvent.n_atoms

    def ionize(
        self,
        ion: str = "Na+",
        seed: int | None = None,
        potential_cutoff: float = 1.4,
        min_dist: float = 0.25,
        random_placement: bool = False,
    ) -> SolvatedSystem:
        """Place counter-ions into the solvated box.

        Parameters
        ----------
        ion : str, default ``"Na+"``
            Ion code.
        seed : int, optional
            Random seed for placement.  Inherits from system seed if ``None``.
        potential_cutoff : float, default 1.4
            Electrostatic potential cutoff for ion placement (nm).
        min_dist : float, default 0.25
            Minimum distance from solute for ion placement (nm).
        random_placement : bool, default False
            If ``True``, place ions randomly (``@random`` flag).

        Returns
        -------
        SolvatedSystem
            Copy with ions placed.
        """
        raise NotImplementedError("ionize — calls gromos-rs 'ion' binary, to be implemented")

    def minimize(
        self,
        steps: int = 50_000,
        n_threads: int = 1,
    ) -> MDResult:
        """Minimise the solvated system.

        Parameters
        ----------
        steps : int, default 50 000
            Maximum steepest-descent steps.
        n_threads : int, default 1
            OMP thread count.

        Returns
        -------
        MDResult
        """
        raise NotImplementedError("System minimisation — to be implemented")

    def equilibrate(
        self,
        stages: list[dict] | None = None,
        n_threads: int = 1,
    ) -> MDResult:
        """Run the multi-stage equilibration protocol.

        Parameters
        ----------
        stages : list of dict, optional
            Each dict may contain ``steps``, ``temperature``, ``ntivel``.
            Default reproduces the vsomm_modeler three-stage protocol::

                [
                    {"steps": 2_000_000, "temperature": 420, "ntivel": 1},
                    {"steps": 500_000, "temperature": 360, "ntivel": 0},
                    {"steps": 500_000, "temperature": 300, "ntivel": 0},
                ]

        n_threads : int, default 1
            OMP thread count.

        Returns
        -------
        MDResult
        """
        default_stages = [
            {"steps": 2_000_000, "temperature": 420, "ntivel": 1},
            {"steps": 500_000, "temperature": 360, "ntivel": 0},
            {"steps": 500_000, "temperature": 300, "ntivel": 0},
        ]
        stages = stages or default_stages
        raise NotImplementedError("System equilibration — to be implemented")

    def run(
        self,
        steps: int,
        dt: float = 0.002,
        temperature: float = 300.0,
        pressure: float | None = None,
        n_threads: int = 1,
        traj_freq: int = 100,
        ene_freq: int = 100,
    ) -> MDResult:
        """Run production MD.

        Parameters
        ----------
        steps : int
            Number of MD steps.
        dt : float, default 0.002
            Time step in ps.
        temperature : float, default 300.0
            Target temperature in K.
        pressure : float, optional
            Target pressure in bar.  If ``None``, runs NVT.
        n_threads : int, default 1
            OMP thread count.
        traj_freq : int, default 100
            Steps between trajectory frames.
        ene_freq : int, default 100
            Steps between energy frames.

        Returns
        -------
        MDResult
            Contains trajectory and energy timeseries.

        Examples
        --------
        >>> ff = ForceField("54a7")  # doctest: +SKIP
        >>> solvated = ...  # doctest: +SKIP
        >>> result = solvated.run(steps=1_000_000, temperature=300)
        >>> df = result.energies.to_dataframe()
        >>> df["total"].mean()
        -48231.4
        """
        raise NotImplementedError("Production MD — to be implemented")

    def __repr__(self) -> str:
        return (
            f"SolvatedSystem(n_solute={self.n_solute_atoms}, "
            f"n_solvent={self.n_solvent}, n_ions={self.n_ions})"
        )


# ---------------------------------------------------------------------------
# Configuration (enhanced, beyond the pyo3 binding)
# ---------------------------------------------------------------------------


class Configuration:
    """Atomic positions, velocities, and box for a GROMOS system.

    Wraps a ``.cnf`` file.  Extends the thin pyo3 binding with coordinate
    manipulation and format conversion.

    Parameters
    ----------
    cnf_file : str or Path
        Path to ``.cnf`` file.

    Attributes
    ----------
    n_atoms : int
    positions : numpy.ndarray, shape (n_atoms, 3)
        Positions in nm.
    velocities : numpy.ndarray, shape (n_atoms, 3)
        Velocities in nm ps⁻¹.
    box : numpy.ndarray, shape (3, 3)
        Box matrix in nm.
    potential_energy : float or None
        Potential energy in kJ mol⁻¹ if available (from minimisation output).

    Notes
    -----
    **Relationship to pyo3 ``Configuration``**: The pyo3 binding exposes a
    thin ``Configuration`` that reads ``.cnf`` files.  This class adds:

    - ``to_pdb()`` / ``to_gro()`` conversion (calls ``frameout``).
    - ``gather(pbc="r")`` for PBC gathering.
    - Arithmetic: ``conf + displacement`` shifts all positions.
    - ``conf.select("1:CA")`` returns a view for a subset of atoms.

    Examples
    --------
    >>> conf = Configuration("min_system.cnf")  # doctest: +SKIP
    >>> conf.n_atoms
    648
    >>> conf.positions.shape
    (648, 3)
    >>> conf.box
    array([[3.1, 0. , 0. ],
           [0. , 3.1, 0. ],
           [0. , 0. , 3.1]])
    """

    def __init__(self, cnf_file: str | Path) -> None:
        self.cnf_file = Path(cnf_file)
        if not self.cnf_file.exists():
            raise FileNotFoundError(self.cnf_file)
        self.potential_energy: float | None = None
        self._positions: np.ndarray | None = None
        self._velocities: np.ndarray | None = None
        self._box: np.ndarray | None = None

    @property
    def n_atoms(self) -> int:
        # TODO: delegate to pyo3 Configuration
        return 0

    @property
    def positions(self) -> np.ndarray:
        if self._positions is None:
            self._positions = np.empty((0, 3))  # TODO: parse via pyo3
        return self._positions

    @property
    def velocities(self) -> np.ndarray:
        if self._velocities is None:
            self._velocities = np.empty((0, 3))
        return self._velocities

    @property
    def box(self) -> np.ndarray:
        if self._box is None:
            self._box = np.zeros((3, 3))
        return self._box

    def to_pdb(self, output: str | Path, pbc: str = "r") -> Path:
        """Write a PDB file using ``frameout``.

        Parameters
        ----------
        output : str or Path
            Output ``.pdb`` path.
        pbc : str, default ``"r"``
            Periodic boundary condition handling (``"r"`` rectangular,
            ``"v"`` vacuum, …).

        Returns
        -------
        Path
            The written PDB file.
        """
        raise NotImplementedError("frameout → pdb — to be implemented")

    def __repr__(self) -> str:
        return f"Configuration(cnf_file={self.cnf_file!r}, n_atoms={self.n_atoms})"


# ---------------------------------------------------------------------------
# MDResult
# ---------------------------------------------------------------------------


class MDResult:
    """Output of a single MD run segment.

    Returned by :meth:`SolvatedSystem.minimize`,
    :meth:`SolvatedSystem.equilibrate`, and :meth:`SolvatedSystem.run`.

    Attributes
    ----------
    conf : Configuration
        Final configuration (the ``.cnf`` / ``.fin`` file).
    energies : EnergyTimeseries or None
        Energy timeseries if the run produced a ``.tre`` file.
    traj_file : Path or None
        Path to the trajectory ``.trc`` file (if written).
    tre_file : Path or None
        Path to the energy trajectory ``.tre`` file.
    imd_file : Path
        Path to the ``.imd`` file used for this run.
    n_steps : int
        Actual steps completed.

    Notes
    -----
    **Method chaining**: :class:`MDResult` intentionally mirrors the
    :class:`SolvatedSystem` API so that runs can be chained::

        result = solvated.minimize(steps=50_000).equilibrate().run(steps=1_000_000, temperature=300)
        df = result.energies.to_dataframe()

    Each step returns an :class:`MDResult`; calling ``.equilibrate()`` on an
    :class:`MDResult` starts from its ``.conf``.

    Examples
    --------
    >>> result = solvated.run(steps=10_000)  # doctest: +SKIP
    >>> result.energies.total.mean()
    -48231.4
    >>> result.conf.n_atoms
    648
    >>> df = result.energies.to_dataframe()
    >>> df.write_csv("energies.csv")  # polars   # doctest: +SKIP
    """

    def __init__(
        self,
        conf: Configuration,
        energies: EnergyTimeseries | None,
        traj_file: Path | None,
        tre_file: Path | None,
        imd_file: Path,
        n_steps: int,
    ) -> None:
        self.conf = conf
        self.energies = energies
        self.traj_file = traj_file
        self.tre_file = tre_file
        self.imd_file = imd_file
        self.n_steps = n_steps

    def equilibrate(self, **kwargs) -> MDResult:
        """Continue with an equilibration stage from the current configuration.

        Returns
        -------
        MDResult
        """
        raise NotImplementedError

    def run(self, steps: int, **kwargs) -> MDResult:
        """Continue with a production run from the current configuration.

        Returns
        -------
        MDResult
        """
        raise NotImplementedError

    def __repr__(self) -> str:
        return (
            f"MDResult(n_steps={self.n_steps}, "
            f"n_atoms={self.conf.n_atoms}, "
            f"has_traj={self.traj_file is not None})"
        )


# ---------------------------------------------------------------------------
# Convenience alias for the full pipeline (Shape B sketch)
# ---------------------------------------------------------------------------


def build_system(
    forcefield: str,
    molecules: list[dict],
    solvent: str = "SPC",
    density: float = 900.0,
    ion: str = "Na+",
    seed: int = 42,
    debug: bool = False,
) -> SolvatedSystem:
    """High-level factory: build a solvated system from a declarative spec.

    This is the **Shape C / declarative** entry point.  It composes
    :class:`ForceField`, :class:`MoleculeTopology`, :class:`SystemTopology`,
    and :class:`SolvatedSystem` automatically.

    Parameters
    ----------
    forcefield : str
        Named force field (e.g. ``"54a7"``).
    molecules : list of dict
        Each dict has ``"sequence"`` (list of str) and ``"n"`` (int)::

            [
                {"sequence": ["STA", "ALA", "GLY", "END"], "n": 10},
                {"sequence": ["STA", "PRO", "END"], "n": 5},
            ]

    solvent : str, default ``"SPC"``
        Solvent code.
    density : float, default 900.0
        Initial density in kg m⁻³.
    ion : str, default ``"Na+"``
        Counter-ion code for charge neutralisation.
    seed : int, default 42
        Random seed.
    debug : bool, default False
        Preserve intermediate files.

    Returns
    -------
    SolvatedSystem
        Ready for :meth:`~SolvatedSystem.minimize` → :meth:`~SolvatedSystem.equilibrate`
        → :meth:`~SolvatedSystem.run`.

    Notes
    -----
    This function is equivalent to::

        ff = ForceField(forcefield)
        sys = sum(
            (ff.molecule(m["sequence"]) * m["n"] for m in molecules), start=SystemTopology([], ff)
        )
        return sys.neutralize(ion, seed=seed).solvate(solvent, density=density, seed=seed)

    For programmatic use (e.g. GA-optimised sequences, dynamic n values) prefer
    the algebra API directly.

    Examples
    --------
    Reproduce the vsomm_modeler use case in ~5 lines:

    >>> solvated = build_system(  # doctest: +SKIP
    ...     forcefield="54a7",
    ...     molecules=[{"sequence": ["STA", "ALA", "GLY", "END"], "n": 200}],
    ...     solvent="SPC",
    ...     density=900,
    ...     ion="Na+",
    ...     seed=42,
    ... )
    >>> result = solvated.minimize().equilibrate().run(steps=1_000_000)
    >>> result.energies.to_dataframe().write_csv("energies.csv")  # doctest: +SKIP
    """
    ff = ForceField(forcefield)
    mols: list[tuple[MoleculeTopology, int]] = [
        (ff.molecule(m["sequence"]), m["n"]) for m in molecules
    ]
    system = SystemTopology(molecules=mols, ff=ff)
    system = system.neutralize(ion, seed=seed)
    return system.solvate(solvent, density=density, seed=seed, debug=debug)
