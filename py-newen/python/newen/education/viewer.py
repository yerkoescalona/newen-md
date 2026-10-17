"""Watch a simulation run: lit spheres beside energy and temperature strip charts.

``view(sim)`` steps a :class:`newen.Simulation`, shows each sampled frame as it arrives, and
then replays the whole run in the browser. ``view(recording)`` replays a :class:`Recording`
made earlier by :func:`record`. Two renderers draw the spheres: three.js (the default; glossy
spheres with soft shadows, a 2D system on a table, a 3D box by its edges) and 3Dmol.js
(ambient occlusion and outlines). Both ship with the package, so the widget works without an
internet connection.

Examples
--------
>>> from newen import Recipe, Simulation, System
>>> from newen.education import record, view
>>> system = System.from_files("water.topo", "water.conf")  # doctest: +SKIP
>>> sim = Simulation(system, Recipe.from_imd("water.imd"))  # doctest: +SKIP
>>> view(sim, frames=100, every=5, system=system)  # doctest: +SKIP
"""

from __future__ import annotations

import time
import warnings
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, Any

import anywidget
import numpy as np
import traitlets

if TYPE_CHECKING:
    from collections.abc import Iterator

    from newen import Simulation, System, Topology

__all__ = [
    "ELEMENTS",
    "SPEED_RAMP",
    "Recording",
    "SphereView",
    "atom_names",
    "element_of",
    "record",
    "speed_colours",
    "view",
]

_HERE = Path(__file__).parent

#: Element -> (RGB colour, sphere radius in nm). Atoms whose element is not listed use ``"X"``.
ELEMENTS: dict[str, tuple[tuple[int, int, int], float]] = {
    "O": ((224, 49, 49), 0.115),
    "H": ((255, 255, 255), 0.075),
    "AR": ((143, 91, 214), 0.17),
    "C": ((110, 110, 110), 0.17),
    "N": ((48, 80, 248), 0.155),
    "X": ((176, 176, 176), 0.15),
}

#: Colour ramp of :func:`speed_colours`, slow to fast: deep blue, violet, magenta, orange, gold.
SPEED_RAMP = ("#1d3fa8", "#6a3fc4", "#c2418f", "#f07d3a", "#ffd56b")
_RAMP_RGB = np.array([[int(h[i : i + 2], 16) for i in (1, 3, 5)] for h in SPEED_RAMP], float)

_ENERGY_SEND_INTERVAL_S = 0.05
_BACKGROUND = {"three": "#0c121c", "3dmol": "#ffffff"}
_SOURCES = {
    "dmol": "_vendor/3Dmol-min.js",
    "three": "_vendor/three.module.min.js",
    "renderer:three": "_renderer_three.js",
    "renderer:3dmol": "_renderer_3dmol.js",
}


@cache
def _source(key: str) -> str:
    return (_HERE / _SOURCES[key]).read_text(encoding="utf-8")


def speed_colours(speed: Any, vmax: float) -> np.ndarray:
    """RGB colours of speeds on :data:`SPEED_RAMP`: slow is deep blue, fast is pale gold.

    Parameters
    ----------
    speed : array_like
        Speeds, any shape, e.g. ``(n_atoms,)`` in nm/ps.
    vmax : float
        The speed drawn at the gold end; faster atoms are drawn gold too.

    Returns
    -------
    numpy.ndarray
        uint8 colours, shape ``speed.shape + (3,)``.
    """
    speed = np.asarray(speed, dtype=np.float64)
    x = speed / vmax if vmax > 0 else np.zeros_like(speed)
    s = np.clip(x, 0.0, 0.9999) * (len(_RAMP_RGB) - 1)
    i = s.astype(int)
    f = (s - i)[..., None]
    rgb: np.ndarray = (_RAMP_RGB[i] * (1 - f) + _RAMP_RGB[i + 1] * f).astype(np.uint8)
    return rgb


def _vmax(speeds: np.ndarray) -> float:
    """The legend's top speed: the 99th percentile, so a rare fast atom does not wash it out."""
    return float(np.percentile(speeds, 99)) if np.size(speeds) else 0.0


def _planar(z: Any, vz: Any = None) -> bool:
    """True when every ``z`` equals the first within 1e-9 nm and every ``vz`` is 0."""
    z = np.asarray(z, dtype=np.float64)
    if z.size == 0 or np.any(np.abs(z - z.flat[0]) > 1e-9):
        return False
    return vz is None or bool(np.all(np.abs(np.asarray(vz)) <= 1e-9))


def element_of(name: str) -> str:
    """The element key of :data:`ELEMENTS` for a GROMOS atom name.

    Parameters
    ----------
    name : str
        Atom name from the topology, e.g. ``"OW"``, ``"HW1"``, ``"AR1"``.

    Returns
    -------
    str
        ``"AR"`` for names starting with AR, else the first letter when it is O, H, C or N,
        else ``"X"``. Ion names carrying a charge sign (``"NA+"``, ``"CL-"``) are ``"X"``.
    """
    upper = name.strip().upper()
    if upper.startswith("AR"):
        return "AR"
    if "+" in upper or "-" in upper or not upper or upper[0] not in "OHCN":
        return "X"
    return upper[0]


@dataclass
class Recording:
    """What a run produced, for viewing or saving.

    Attributes
    ----------
    frames : numpy.ndarray
        Positions, shape ``(n_frames, n_atoms, 3)``, nm.
    energy : numpy.ndarray
        One row per frame, shape ``(n_frames, 5)``: time (ps), kinetic, potential and total
        energy (kJ/mol), temperature (K). The energies are those of the step just taken, so
        the positions of a frame are one time step later than its time.
    box : tuple of float
        Rectangular box edges in nm; ``(0, 0, 0)`` in vacuum.
    names : list of str
        Atom names.
    elements : list of str
        Element key per atom (see :func:`element_of`).
    radii : numpy.ndarray
        Sphere radius per atom, nm.
    speeds : numpy.ndarray or None
        Speed of every atom in every frame, shape ``(n_frames, n_atoms)``, nm/ps, when known
        (:func:`record` keeps them); ``view(recording, colors="speed")`` needs them.
    """

    frames: np.ndarray
    energy: np.ndarray
    box: tuple[float, float, float]
    names: list[str]
    elements: list[str]
    radii: np.ndarray
    speeds: np.ndarray | None = None

    def save_cnf(self, path: str | Path, frame: int = 0) -> None:
        """Write one frame as a GROMOS96 ``.cnf`` (TITLE, POSITION, GENBOX), in nm.

        :class:`newen.Configuration` reads it back; positions keep 9 decimals.

        Parameters
        ----------
        path : str or pathlib.Path
            File to write.
        frame : int, default 0
            Index of the frame to write.
        """
        lines = ["TITLE", f"newen.education frame {frame} at t = {self.energy[frame, 0]:.4f} ps"]
        lines += ["END", "POSITION"]
        for i, (name, xyz) in enumerate(zip(self.names, self.frames[frame], strict=True)):
            x, y, z = xyz
            lines.append(f"{1:>5} {'MOL':<5} {name[:6]:<6}{i + 1:>6}{x:15.9f}{y:15.9f}{z:15.9f}")
        rectangular = all(b > 0 for b in self.box)
        angle = 90.0 if rectangular else 0.0
        lines += ["END", "GENBOX", f"{int(rectangular):>8}"]
        lines.append("".join(f"{b:15.9f}" for b in self.box))
        lines.append(f"{angle:15.9f}" * 3)
        lines += [f"{0.0:15.9f}" * 3] * 2 + ["END", ""]
        Path(path).write_text("\n".join(lines))

    def save_pdb(self, path: str | Path) -> None:
        """Write every frame as a ``MODEL``/``ENDMDL`` record of a PDB file, in Angstrom.

        Parameters
        ----------
        path : str or pathlib.Path
            File to write.
        """
        out = [_cryst1(self.box)] if all(b > 0 for b in self.box) else []
        for k, xyz in enumerate(self.frames):
            out.append(f"MODEL     {k + 1:>4}")
            out += _atom_records(self.names, self.elements, xyz)
            out.append("ENDMDL")
        out += ["END", ""]
        Path(path).write_text("\n".join(out))


def _cryst1(box: tuple[float, float, float]) -> str:
    a, b, c = (10 * x for x in box)
    return f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1"


def _atom_records(names: list[str], elements: list[str], xyz_nm: np.ndarray) -> list[str]:
    records = []
    for i, (name, elem, (x, y, z)) in enumerate(zip(names, elements, 10 * xyz_nm, strict=True)):
        symbol = "" if elem == "X" else elem
        records.append(
            f"ATOM  {(i + 1) % 100000:>5} {name[:4]:<4} MOL A{1:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {symbol:>2}"
        )
    return records


def _blocks(text: str) -> dict[str, list[str]]:
    """Blocks of GROMOS text by name, as their non-comment lines."""
    blocks: dict[str, list[str]] = {}
    current: list[str] | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if current is None:
            if stripped:
                current = blocks.setdefault(stripped, [])
        elif stripped == "END":
            current = None
        elif stripped and not stripped.startswith("#"):
            current.append(stripped)
    return blocks


def atom_names(topology: Topology) -> list[str]:
    """Atom names of a topology, solute first, then each solvent molecule.

    Parameters
    ----------
    topology : newen.Topology
        The topology. Its ``atom_names`` are used when the extension has them; otherwise
        its text form is read for the SOLUTEATOM and SOLVENTATOM blocks.

    Returns
    -------
    list of str
        One name per atom.
    """
    names = getattr(topology, "atom_names", None)
    return list(names) if names is not None else _parsed_atom_names(topology)


def _parsed_atom_names(topology: Topology) -> list[str]:
    """Atom names read from the SOLUTEATOM and SOLVENTATOM blocks of ``str(topology)``."""
    blocks = _blocks(str(topology))
    tokens = " ".join(blocks.get("SOLUTEATOM", [])).split()
    names: list[str] = []
    if tokens:
        i = 1
        for _ in range(int(tokens[0])):
            # ATNM MRES PANM IAC MASS CG CGC INE <INE exclusions> INE14 <INE14 partners>
            names.append(tokens[i + 2])
            i += 8 + int(tokens[i + 7])
            i += 1 + int(tokens[i])
    solvent = [line.split()[1] for line in blocks.get("SOLVENTATOM", [])[1:]]
    if solvent:
        names += solvent * (topology.n_solvent_atoms // len(solvent))
    return names


def _describe(
    n_atoms: int, system: System | Topology | None, volume: float = 0.0
) -> tuple[list[str], tuple[float, float, float]]:
    """Atom names and box for ``n_atoms`` atoms, from ``system`` when given."""
    if system is None:
        edge = volume ** (1 / 3) if volume > 0 else 0.0
        return ["X"] * n_atoms, (edge, edge, edge)
    topology = system.topology if hasattr(system, "topology") else system
    a, b, c = system.box if hasattr(system, "box") else (0.0, 0.0, 0.0)
    names = atom_names(topology)
    if len(names) != n_atoms:
        warnings.warn(
            f"the topology names {len(names)} atoms but the simulation has {n_atoms}; "
            "drawing every atom as an unknown element",
            stacklevel=3,
        )
        names = ["X"] * n_atoms
    return names, (a, b, c)


def _speeds(sim: Simulation) -> np.ndarray:
    return np.linalg.norm(np.asarray(sim.velocities, dtype=np.float64), axis=1)


def _samples(
    sim: Simulation, frames: int, every: int
) -> Iterator[tuple[np.ndarray, list[float], np.ndarray]]:
    """Step ``sim`` ``frames`` times by ``every`` steps; yield positions, energy row, speeds."""
    for _ in range(frames):
        sim.step(every)
        e = sim.energies
        row = [sim.time - sim.dt, e.kinetic, e.potential, e.total, sim.temperature]
        yield np.array(sim.positions, dtype=np.float64), row, _speeds(sim)


def record(
    sim: Simulation, frames: int = 300, every: int = 10, *, system: System | None = None
) -> Recording:
    """Run ``sim`` and keep its positions, energies and temperature.

    Parameters
    ----------
    sim : newen.Simulation
        The simulation; it is advanced by ``frames * every`` steps.
    frames : int, default 300
        Number of samples.
    every : int, default 10
        MD steps between samples.
    system : newen.System or newen.Topology, optional
        Where the atom names (hence colours and radii) and the box come from: a
        ``Simulation`` does not expose them. Without it every atom is drawn as an unknown
        element and the box is taken as a cube of the simulation's volume.

    Returns
    -------
    Recording
        ``frames`` positions, energy rows and atom speeds.
    """
    names, box = _describe(sim.n_atoms, system, sim.volume)
    positions, rows, speeds = [], [], []
    for xyz, row, speed in _samples(sim, frames, every):
        positions.append(xyz)
        rows.append(row)
        speeds.append(speed)
    return _recording(names, box, positions, rows, speeds)


def _recording(names, box, positions, rows, speeds) -> Recording:
    elements = [element_of(n) for n in names]
    n = len(names)
    return Recording(
        frames=np.stack(positions) if positions else np.empty((0, n, 3)),
        energy=np.array(rows, dtype=np.float64).reshape(-1, 5),
        box=box,
        names=names,
        elements=elements,
        radii=np.array([ELEMENTS[e][1] for e in elements]),
        speeds=np.stack(speeds) if speeds else np.empty((0, n)),
    )


def _angstrom_bytes(xyz_nm: Any) -> bytes:
    return (np.asarray(xyz_nm, dtype=np.float64) * 10).astype(np.float32).tobytes()


def _rgb_bytes(rgb: Any) -> bytes:
    return np.ascontiguousarray(rgb, dtype=np.uint8).tobytes()


class SphereView(anywidget.AnyWidget):
    """Atoms as lit spheres beside energy and temperature charts.

    ``renderer`` picks the drawing: ``"three"`` (three.js: glossy spheres with soft shadows, a
    table under a ``flat`` system, the box by its edges) or ``"3dmol"`` (3Dmol.js: ambient
    occlusion and outlines). Both are shipped with the package and reach the browser as string
    traits, so nothing is fetched from the network.

    Coordinates travel as float32 bytes in Angstrom. ``positions`` is shown at the next repaint
    (the newest wins when updates arrive faster than the screen); ``trajectory`` is replayed in
    the browser at ``fps`` frames per second while ``playing`` (``fps <= 0``: every repaint).
    ``colors`` (uint8 RGB per atom) and ``color_trajectory`` (one such block per frame of
    ``trajectory``) set the sphere colours; ``vmax > 0`` draws the speed legend.
    ``energy`` holds rows ``[t, kinetic, potential, total, T]`` for the charts and readouts.

    ``status`` is set by the browser: ``"ready"`` once the first frame is drawn, or a message
    naming the step that failed (loading three.js, loading the renderer module, loading
    3Dmol.js, creating the scene, drawing the first frame) and the error. The same message is
    shown in the stage, and the charts still draw. It stays ``""`` until a front end has run.
    """

    _esm = _HERE / "_view.js"
    renderer = traitlets.Enum(("three", "3dmol"), "three").tag(sync=True)
    renderer_src = traitlets.Unicode().tag(sync=True)
    three_src = traitlets.Unicode().tag(sync=True)
    dmol_src = traitlets.Unicode().tag(sync=True)
    pdb = traitlets.Unicode("").tag(sync=True)
    styles: traitlets.List[Any] = traitlets.List().tag(sync=True)
    view_style = traitlets.Dict(
        {"style": "outline ambientOcclusion", "color": "black", "width": 0.03, "strength": 1.5}
    ).tag(sync=True)
    background = traitlets.Unicode().tag(sync=True)
    box_angstrom = traitlets.List(traitlets.Float()).tag(sync=True)
    flat = traitlets.Bool(False).tag(sync=True)
    width = traitlets.Int(900).tag(sync=True)
    height = traitlets.Int(600).tag(sync=True)
    radii = traitlets.Bytes(b"").tag(sync=True)
    colors = traitlets.Bytes(b"").tag(sync=True)
    color_trajectory = traitlets.Bytes(b"").tag(sync=True)
    vmax = traitlets.Float(0.0).tag(sync=True)
    ramp = traitlets.List(traitlets.Unicode(), list(SPEED_RAMP)).tag(sync=True)
    positions = traitlets.Bytes(b"").tag(sync=True)
    trajectory = traitlets.Bytes(b"").tag(sync=True)
    playing = traitlets.Bool(False).tag(sync=True)
    fps = traitlets.Float(30.0).tag(sync=True)
    energy: traitlets.List[Any] = traitlets.List().tag(sync=True)
    status = traitlets.Unicode("").tag(sync=True)

    #: The run behind this view, when :func:`view` made it.
    recording: Recording | None = None
    _shown_by_view = False

    @traitlets.default("renderer_src")
    def _default_renderer_src(self) -> str:
        return _source(f"renderer:{self.renderer}")

    @traitlets.default("three_src")
    def _default_three_src(self) -> str:
        return _source("three") if self.renderer == "three" else ""

    @traitlets.default("dmol_src")
    def _default_dmol_src(self) -> str:
        return _source("dmol") if self.renderer == "3dmol" else ""

    @traitlets.default("background")
    def _default_background(self) -> str:
        return _BACKGROUND[self.renderer]

    @traitlets.observe("renderer")
    def _renderer_changed(self, _change: dict[str, Any]) -> None:
        self.renderer_src = self._default_renderer_src()
        self.three_src = self._default_three_src()
        self.dmol_src = self._default_dmol_src()

    def set_positions(self, xyz_nm: Any) -> None:
        """Show new coordinates, shape ``(n_atoms, 3)`` in nm."""
        self.positions = _angstrom_bytes(xyz_nm)

    def play(self, frames_nm: Any, fps: float = 30.0) -> None:
        """Replay ``frames_nm``, shape ``(n_frames, n_atoms, 3)`` in nm, in the browser."""
        self.trajectory = _angstrom_bytes(frames_nm)
        self.fps = fps
        self.playing = True

    def _repr_mimebundle_(self, **kwargs: Any) -> Any:
        # view() already displayed this widget; its own cell's result must not show it again.
        if self._shown_by_view:
            return {"text/plain": ""}
        return super()._repr_mimebundle_(**kwargs)


def _styles(radii_nm: np.ndarray, rgb: np.ndarray) -> list[list[dict[str, Any]]]:
    """3Dmol sphere styles, one per distinct (radius, colour), each selecting its atoms."""
    groups: dict[tuple[float, str], list[int]] = {}
    for i, (r, c) in enumerate(zip(radii_nm, rgb, strict=True)):
        key = (round(float(r) * 10, 4), "#{:02x}{:02x}{:02x}".format(*(int(v) for v in c)))
        groups.setdefault(key, []).append(i)
    return [
        [{"index": idx}, {"sphere": {"radius": r, "color": colour}}]
        for (r, colour), idx in groups.items()
    ]


def _widget(
    xyz_nm: np.ndarray,
    names: list[str],
    box: tuple[float, float, float],
    radii: Any,
    colors: Any,
    **traits: Any,
) -> SphereView:
    n = len(names)
    elements = [element_of(name) for name in names]
    radii = np.array([ELEMENTS[e][1] for e in elements]) if radii is None else np.asarray(radii)
    rgb = np.array([ELEMENTS[e][0] for e in elements]) if colors is None else np.asarray(colors)
    if radii.shape != (n,):
        raise ValueError(f"radii must have shape ({n},), got {radii.shape}")
    if rgb.shape != (n, 3):
        raise ValueError(f"colors must have shape ({n}, 3), got {rgb.shape}")
    pdb = "\n".join([*_atom_records(names, elements, np.asarray(xyz_nm)), "END", ""])
    widget = SphereView(
        pdb=pdb,
        styles=_styles(radii, rgb),
        radii=(10 * radii).astype(np.float32).tobytes(),
        colors=_rgb_bytes(rgb),
        box_angstrom=[10 * b for b in box] if all(b > 0 for b in box) else [],
        **traits,
    )
    widget.set_positions(xyz_nm)
    return widget


def _show(widget: SphereView) -> None:
    """Display ``widget`` now when running in a Jupyter kernel, and not again as the cell's result."""
    try:
        from IPython import get_ipython
        from IPython.display import display
    except ImportError:
        return
    shell = get_ipython()
    if shell is None or getattr(shell, "kernel", None) is None:
        return
    display(widget)
    widget._shown_by_view = True

    def _forget(*_: Any) -> None:
        widget._shown_by_view = False
        shell.events.unregister("post_run_cell", _forget)

    shell.events.register("post_run_cell", _forget)


def view(
    target: Simulation | Recording,
    *,
    frames: int = 300,
    every: int = 10,
    renderer: str = "three",
    radii: Any = None,
    colors: Any = None,
    flat: bool | None = None,
    width: int = 900,
    height: int = 600,
    background: str | None = None,
    system: System | Topology | None = None,
) -> SphereView:
    """Show a simulation as it runs, or replay a recording.

    Parameters
    ----------
    target : newen.Simulation or Recording
        A ``Simulation`` is advanced by ``frames * every`` steps; the widget is displayed
        first and shows each sample as it arrives, then replays the whole run. A
        ``Recording`` is replayed.
    frames : int, default 300
        Number of samples (``Simulation`` only).
    every : int, default 10
        MD steps between samples (``Simulation`` only).
    renderer : {"three", "3dmol"}, default "three"
        ``"three"``: glossy spheres with soft shadows, a 2D system on a table, a 3D box by its
        edges. ``"3dmol"``: ambient occlusion and outlines, the molecular look.
    radii : array_like, optional
        Sphere radius per atom, nm, shape ``(n_atoms,)``. Default: by element.
    colors : array_like or "speed", optional
        RGB colour per atom, uint8, shape ``(n_atoms, 3)``; or ``"speed"``, which colours
        every frame by each atom's speed on :data:`SPEED_RAMP` (see :func:`speed_colours`) up
        to the 99th percentile of the run's speeds, and draws the legend. Default: by element.
    flat : bool, optional
        Draw the atoms on a table, as a 2D system. Default: inferred, true when every atom
        of the first frame shares one z (within 1e-9 nm) and none moves along z.
    width, height : int, default 900, 600
        Size in pixels of the whole widget and of the 3D view.
    background : str, optional
        Background colour of the 3D view. Default: ``"#0c121c"`` for three.js, white for
        3Dmol.js.
    system : newen.System or newen.Topology, optional
        Atom names and box for a ``Simulation``, which does not expose them (see
        :func:`record`).

    Returns
    -------
    SphereView
        The widget; ``.recording`` holds the run it shows.
    """
    by_speed = isinstance(colors, str)
    if by_speed and colors != "speed":
        raise ValueError(f'colors must be an array or "speed", got {colors!r}')
    traits = {"renderer": renderer, "width": width, "height": height}
    if background is not None:
        traits["background"] = background

    if isinstance(target, Recording):
        if by_speed and (target.speeds is None or target.speeds.shape != target.frames.shape[:2]):
            raise ValueError('colors="speed" needs a recording with speeds, as record() makes')
        vmax = _vmax(target.speeds) if by_speed and target.speeds is not None else 0.0
        rgb = (
            speed_colours(target.speeds[0], vmax)
            if by_speed and target.speeds is not None
            else colors
        )
        radii = target.radii if radii is None else radii
        flat = _planar(target.frames[..., 2]) if flat is None else flat
        widget = _widget(
            target.frames[0], target.names, target.box, radii, rgb, flat=flat, **traits
        )
        widget.energy = target.energy.tolist()
        widget.recording = target
        _replay(widget, target, vmax)
        return widget

    names, box = _describe(target.n_atoms, system, target.volume)
    xyz = np.asarray(target.positions)
    vmax = _vmax(_speeds(target)) if by_speed else 0.0
    rgb = speed_colours(_speeds(target), vmax) if by_speed else colors
    if flat is None:
        flat = _planar(xyz[:, 2], np.asarray(target.velocities)[:, 2])
    widget = _widget(xyz, names, box, radii, rgb, flat=flat, vmax=vmax, **traits)
    _show(widget)
    positions, rows, speeds, sent = [], [], [], 0.0
    for xyz, row, speed in _samples(target, frames, every):
        positions.append(xyz)
        rows.append(row)
        speeds.append(speed)
        with widget.hold_sync():
            if by_speed:
                vmax = max(vmax, _vmax(speed))
                widget.vmax = vmax
                widget.colors = _rgb_bytes(speed_colours(speed, vmax))
            widget.set_positions(xyz)
        if time.perf_counter() - sent >= _ENERGY_SEND_INTERVAL_S:
            widget.energy, sent = list(rows), time.perf_counter()
    widget.energy = list(rows)
    widget.recording = _recording(names, box, positions, rows, speeds)
    assert widget.recording.speeds is not None
    _replay(widget, widget.recording, _vmax(widget.recording.speeds) if by_speed else 0.0)
    return widget


def _replay(widget: SphereView, rec: Recording, vmax: float) -> None:
    """Hand ``rec`` to the browser for replay; ``vmax > 0`` colours every frame by speed."""
    with widget.hold_sync():
        if vmax > 0:
            widget.vmax = vmax
            widget.color_trajectory = _rgb_bytes(speed_colours(rec.speeds, vmax))
        widget.play(rec.frames)
