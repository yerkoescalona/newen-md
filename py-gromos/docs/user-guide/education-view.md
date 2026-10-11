# Watching a run: `gromos.education.viewer`

`view()` shows a simulation in a Jupyter notebook as it runs: the atoms as lit spheres
(element colours, or colours by speed) in a 3D view you can rotate and zoom, and beside it the
kinetic, potential and total energy and the temperature as strip charts that follow the frame
on screen. The notebook `notebooks/04_live_view.ipynb` is a short tour.

## Install

The viewer needs the `education` extra, which adds `anywidget` (and through it `ipywidgets`):

```bash
pip install "gromos[education]"
```

It works offline. The 3D drawing is done by [three.js](https://threejs.org/) (MIT) or, on
request, [3Dmol.js](https://3dmol.csb.pitt.edu/) (BSD-3-Clause); both ship inside the package
and are handed to the browser by the widget itself, so nothing is fetched from the internet.

## Use

```python
from gromos import Recipe, Simulation, System
from gromos.education import record, view

system = System.from_files("water_216_box.topo", "water_216_box.conf")
sim = Simulation(system, Recipe.from_imd("water_216_box.in"))

view(sim, frames=100, every=5, system=system)    # runs 500 steps, showing each sample
```

`view(sim)` displays the widget at once, advances the simulation by `frames * every` steps and
sends every sample to the widget as it is taken. When the run ends, the widget replays it in
the browser with **Play/Pause** and a frame slider; replaying needs no Python.

```python
rec = record(sim, frames=200, every=5, system=system)   # run without drawing
view(rec)                                               # replay
rec.save_pdb("run.pdb")                                 # every frame, MODEL/ENDMDL, Angstrom
rec.save_cnf("last.cnf", frame=-1)                      # one frame as a GROMOS96 .cnf, nm
```

`save_pdb` is the file to hand to a desktop viewer (VMD, PyMOL, ChimeraX): every frame as one
multi-model PDB; `save_cnf` writes a frame that `gromos.Configuration` reads back.

| Call | Does |
|---|---|
| `record(sim, frames, every, *, system=None)` | Steps `sim`; returns a `Recording`. |
| `view(target, *, frames=300, every=10, renderer="three", radii=None, colors=None, flat=None, width=900, height=600, background=None, system=None)` | Shows a `Simulation` live, or replays a `Recording`. Returns the `SphereView` widget. |
| `Recording.frames` | Positions, `(n_frames, n_atoms, 3)`, nm. |
| `Recording.energy` | `(n_frames, 5)`: time (ps), kinetic, potential, total (kJ/mol), temperature (K). |
| `Recording.box`, `.names`, `.elements`, `.radii` | Box edges (nm), atom names, element per atom, sphere radius per atom (nm). |
| `Recording.speeds` | Speed of every atom in every frame, `(n_frames, n_atoms)`, nm/ps. |
| `speed_colours(speed, vmax)` | uint8 RGB per speed on the ramp deep blue → violet → magenta → orange → pale gold. |
| `SphereView.set_positions(xyz_nm)` / `.play(frames_nm, fps)` | Show one frame / replay frames in the browser. |

**Atom names.** A `Simulation` does not expose its atom names or box, so pass the `System`
(or `Topology`) it was built from as `system=`. The names set colour and size: O red 0.115 nm,
H white 0.075 nm, AR purple 0.17 nm, C grey 0.17 nm, N blue 0.155 nm, anything else light grey
0.15 nm. Without `system=` every atom is drawn light grey and the box is a cube of the
simulation's volume.

**Your own colours and sizes.** `colors` is an `(n_atoms, 3)` array of 0–255 RGB values,
`radii` an `(n_atoms,)` array in nm. `colors="speed"` colours every frame by each atom's speed,
slow deep blue to fast pale gold, up to the 99th percentile of the run's speeds, and draws that
scale in the corner of the view.

## Two renderers

| `renderer=` | Look | Pick it for |
|---|---|---|
| `"three"` (default) | three.js: glossy spheres lit by a key light with soft shadows, on a dark stage. A 2D system sits on a table; a 3D box is drawn by its edges. | Gases, liquids and toy systems, where motion and speed are the point; the 2D experiments. |
| `"3dmol"` | 3Dmol.js: ambient occlusion and black outlines on a white background, the look of a molecular viewer. | Molecules, where shape and packing matter more than motion. |

Both draw the same atoms with the same colours and sizes, replay in the browser and work
offline. `flat` chooses the table: by default it is on when every atom of the first frame
shares one z and none moves along z (as in `lj_gas(dim=2)`); pass `flat=True` or `flat=False`
to choose. 3Dmol.js ignores it. `background` changes the stage colour of either.

**Energies and positions.** The energies the engine reports after a step are those of the step
just taken, so a frame's positions are one time step later than the time in its energy row.

## If the stage stays blank

The widget writes what went wrong into the stage: a short message naming the step that failed
(loading three.js, loading the renderer module, loading 3Dmol.js, creating the scene) and the
browser's error. The charts and readouts still draw. From Python, `w.status` holds the same
message once the front end has run, or `"ready"` when the first frame was drawn.

- **After upgrading `gromos`, restart the kernel** and re-run the cell, so the browser gets the
  new widget code.
- **In VS Code**, set `"jupyter.widgetScriptSources": ["jsdelivr.com", "unpkg.com"]` in the
  settings. VS Code fetches anywidget's own loader from a CDN; this is the only network access,
  and three.js and 3Dmol.js still come from the package.
