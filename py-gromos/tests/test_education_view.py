"""gromos.education.viewer: recording a run, the widget's synced state, and the saved files."""

import importlib
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("anywidget")

from gromos import Configuration, Recipe, Simulation, System, Topology  # noqa: E402
from gromos.education import Recording, SphereView, lj_gas, record, view  # noqa: E402

REF_DIR = (
    Path(__file__).parent.parent.parent / "crates" / "gromos-md" / "tests" / "gromosXX_references"
)
view_module = importlib.import_module("gromos.education.viewer")
EDUCATION = Path(view_module.__file__).parent


def _simulation(name: str) -> tuple[Simulation, System]:
    d = REF_DIR / name
    system = System.from_files(str(d / f"{name}.topo"), str(d / f"{name}.conf"))
    return Simulation(system, Recipe.from_imd(str(d / f"{name}.in"))), system


@pytest.fixture(scope="module")
def argon() -> Recording:
    sim, system = _simulation("pair_lj")
    return record(sim, frames=50, every=2, system=system)


def test_record_shapes_time_and_energy(argon):
    assert argon.frames.shape == (50, 2, 3)
    assert argon.energy.shape == (50, 5)
    assert np.all(np.diff(argon.energy[:, 0]) > 0)
    total = argon.energy[:, 3]
    assert np.abs(total - total[0]).max() < 1e-3
    assert argon.names == ["AR1", "AR2"]
    assert argon.elements == ["AR", "AR"]
    np.testing.assert_allclose(argon.radii, [0.17, 0.17])


def test_view_of_a_recording_syncs_bytes_pdb_and_styles(argon):
    w = view(argon)
    assert isinstance(w, SphereView)
    assert len(w.positions) == 2 * 3 * 4
    assert len(w.trajectory) == 50 * 2 * 3 * 4
    np.testing.assert_allclose(
        np.frombuffer(w.trajectory, np.float32).reshape(50, 2, 3), 10 * argon.frames, atol=1e-4
    )
    assert sum(line.startswith("ATOM") for line in w.pdb.splitlines()) == 2
    radius = {}
    for selection, style in w.styles:
        for i in selection["index"]:
            radius[i] = style["sphere"]["radius"]
    assert radius == {0: pytest.approx(1.7), 1: pytest.approx(1.7)}
    assert len(w.energy) == 50 and len(w.energy[0]) == 5
    assert w.playing


def test_view_of_a_simulation_runs_it_and_fills_the_charts():
    sim, system = _simulation("water_216_box")
    start = sim.current_step
    w = view(sim, frames=20, every=5, system=system)
    assert sim.current_step - start == 100
    assert len(w.energy) == 20
    assert w.recording.frames.shape == (20, 648, 3)
    assert len(w.trajectory) == 20 * 648 * 3 * 4
    np.testing.assert_allclose(w.box_angstrom, [18.652] * 3)
    by_element = dict(zip(w.recording.elements, w.recording.radii, strict=True))
    assert by_element == {"O": pytest.approx(0.115), "H": pytest.approx(0.075)}


def test_save_cnf_reads_back_as_a_configuration(argon, tmp_path):
    path = tmp_path / "frame.cnf"
    argon.save_cnf(path, frame=7)
    conf = Configuration(str(path))
    # GROMOS96 text in nm with 9 decimals: 1e-3 nm is far above the rounding.
    np.testing.assert_allclose(conf.positions, argon.frames[7], atol=1e-3)


def test_save_pdb_writes_one_model_per_frame(argon, tmp_path):
    path = tmp_path / "run.pdb"
    argon.save_pdb(path)
    lines = path.read_text().splitlines()
    assert sum(line.startswith("MODEL") for line in lines) == 50
    assert sum(line.startswith("ENDMDL") for line in lines) == 50
    first = next(line for line in lines if line.startswith("ATOM"))
    # Angstrom in the PDB, nm in the recording.
    assert float(first[30:38]) == pytest.approx(10 * argon.frames[0, 0, 0], abs=1e-3)


FRONT_END = ("_view.js", "_renderer_three.js", "_renderer_3dmol.js")
VENDORED = ("3Dmol-min.js", "3Dmol-LICENSE", "three.module.min.js", "three-LICENSE")
THREE_BANNER = (
    "/**\n * @license\n * Copyright 2010-2024 Three.js Authors\n * SPDX-License-Identifier: MIT"
)


def test_vendored_libraries_and_no_network_in_the_front_end():
    dmol = (EDUCATION / "_vendor" / "3Dmol-min.js").read_text()
    assert dmol.startswith("/*! For license information please see 3Dmol-min.js.LICENSE.txt */")
    assert "BSD" in (EDUCATION / "_vendor" / "3Dmol-LICENSE").read_text()
    three = (EDUCATION / "_vendor" / "three.module.min.js").read_text()
    assert three.startswith(THREE_BANNER)
    assert "The MIT License" in (EDUCATION / "_vendor" / "three-LICENSE").read_text()
    for name in FRONT_END:
        js = (EDUCATION / name).read_text()
        assert "http://" not in js and "https://" not in js, name


def test_front_end_files_contain_no_http():
    for name in FRONT_END:
        assert "http" not in (EDUCATION / name).read_text(), name


def test_status_starts_empty_and_is_synced(argon):
    w = view(argon)
    assert w.status == ""
    assert w.trait_metadata("status", "sync") is True
    assert w.get_state("status") == {"status": ""}


def test_the_shell_reports_load_failures_and_falls_back_to_data_urls():
    js = (EDUCATION / "_view.js").read_text()
    assert "data:text/javascript;base64" in js
    assert "TextEncoder" in js and "window.__gromos_view_load" in js
    load = js[js.index("(async (name, step) => {") :]
    assert "try {" in load and "} catch (e) {" in load
    assert load.index('importSource("three"') < load.index("} catch (e) {")
    assert load.index("await create(") < load.index("} catch (e) {")
    assert 'setStatus("ready")' in load
    # 3Dmol.js loads through the shell's fallback, not its own script tag.
    dmol = (EDUCATION / "_renderer_3dmol.js").read_text()
    assert "loadScript(src)" in dmol and 'createElement("script")' not in dmol


@pytest.mark.parametrize("renderer", ["three", "3dmol"])
def test_both_renderers_draw_the_same_recording(argon, renderer):
    w = view(argon, renderer=renderer)
    assert w.renderer == renderer
    assert w.renderer_src == (EDUCATION / f"_renderer_{renderer}.js").read_text()
    assert "create(" in w.renderer_src
    if renderer == "three":
        assert w.three_src.startswith(THREE_BANNER)
        assert w.dmol_src == ""
        assert w.background == "#0c121c"
    else:
        assert w.three_src == ""
        assert w.dmol_src.startswith("/*! For license information")
        assert w.background == "#ffffff"
    np.testing.assert_allclose(np.frombuffer(w.radii, np.float32), [1.7, 1.7], rtol=1e-6)
    rgb = np.frombuffer(w.colors, np.uint8).reshape(2, 3)
    assert rgb.tolist() == [list(view_module.ELEMENTS["AR"][0])] * 2
    assert len(w.trajectory) == 50 * 2 * 3 * 4


def test_flat_is_inferred_for_a_2d_gas_and_not_for_a_3d_one():
    gas_2d = record(Simulation(*lj_gas(n=16, dim=2)), 5, 2)
    gas_3d = record(Simulation(*lj_gas(n=64, dim=3)), 5, 2)
    assert view(gas_2d).flat is True
    assert view(gas_3d).flat is False
    assert view(gas_2d, flat=False).flat is False
    assert view(Simulation(*lj_gas(n=16, dim=2)), frames=2, every=1).flat is True


def test_colors_by_speed_use_the_ramp_and_send_the_legend():
    rec = record(Simulation(*lj_gas(n=16, dim=2)), 5, 2)
    assert rec.speeds.shape == (5, 16)
    w = view(rec, colors="speed")
    vmax = np.percentile(rec.speeds, 99)
    assert w.vmax == pytest.approx(vmax)
    rgb = np.frombuffer(w.colors, np.uint8).reshape(16, 3)
    np.testing.assert_array_equal(rgb, view_module.speed_colours(rec.speeds[0], vmax))
    assert len(w.color_trajectory) == 5 * 16 * 3
    assert w.ramp == list(view_module.SPEED_RAMP)
    with pytest.raises(ValueError, match="speed"):
        view(
            Recording(rec.frames, rec.energy, rec.box, rec.names, rec.elements, rec.radii),
            colors="speed",
        )


def test_speed_colours_run_from_deep_blue_to_pale_gold():
    rgb = view_module.speed_colours(np.array([0.0, 0.5, 1.0, 5.0]), 1.0)
    assert rgb.dtype == np.uint8 and rgb.shape == (4, 3)
    assert rgb[0].tolist() == [0x1D, 0x3F, 0xA8]
    assert rgb[1].tolist() == [0xC2, 0x41, 0x8F]
    assert rgb[3].tolist() == rgb[2].tolist()
    assert np.abs(rgb[2].astype(int) - [0xFF, 0xD5, 0x6B]).max() <= 1


def test_colors_by_speed_are_accepted(argon):
    speed = np.linalg.norm(np.diff(argon.frames, axis=0)[0], axis=1)
    shade = (255 * speed / speed.max()).astype(np.uint8)
    colors = np.stack([shade, np.zeros_like(shade), 255 - shade], axis=1)
    w = view(argon, colors=colors)
    hexes = {style["sphere"]["color"] for _, style in w.styles}
    assert hexes == {"#{:02x}{:02x}{:02x}".format(*c) for c in colors}
    with pytest.raises(ValueError, match="colors"):
        view(argon, colors=colors[:1])


def test_atom_names_binding_agrees_with_the_topology_text():
    if not hasattr(Topology, "atom_names"):
        pytest.skip("this build of the extension has no Topology.atom_names")
    d = REF_DIR / "water_216_box"
    topology = Topology(str(d / "water_216_box.topo"))
    assert list(topology.atom_names) == view_module._parsed_atom_names(topology)
