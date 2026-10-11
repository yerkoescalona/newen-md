"""
`Configuration` and `System` construction.

Two properties are pinned here. A configuration read from text keeps what one read from a file
keeps — velocities included — and writes them back out. A `System` of a solute-only topology
and a solvated coordinate file is solvated when it is built, so every count it reports agrees
with the coordinates, and a `Simulation` built from it has the same atoms.
"""

from pathlib import Path

import numpy as np
import pytest

import gromos
from gromos import Configuration, Recipe, Simulation, System

REF_DIR = (
    Path(__file__).parent.parent.parent / "crates" / "gromos-md" / "tests" / "gromosXX_references"
)
SHARED = REF_DIR / "shared"

CONFS_WITH_VELOCITIES = sorted(
    p for p in REF_DIR.rglob("*.conf") if "\nVELOCITY\n" in p.read_text().replace("\r", "")
)


def _id(path: Path) -> str:
    return str(path.relative_to(REF_DIR))


def test_there_are_reference_configurations_with_velocities():
    assert len(CONFS_WITH_VELOCITIES) > 10


@pytest.mark.parametrize("path", CONFS_WITH_VELOCITIES, ids=_id)
def test_from_string_keeps_what_from_file_keeps(path: Path):
    from_file = Configuration(str(path))
    from_text = Configuration.from_string(path.read_text())

    assert from_file.velocities.shape == from_file.positions.shape
    assert np.array_equal(from_text.velocities, from_file.velocities)
    assert np.array_equal(from_text.positions, from_file.positions)
    assert from_text.box_dimensions == from_file.box_dimensions


@pytest.mark.parametrize("path", CONFS_WITH_VELOCITIES, ids=_id)
def test_velocities_round_trip_through_the_text(path: Path):
    # The text is gromosXX's layout, nine decimals (out_configuration.cc). Inputs written with
    # more digits (water_216_box.conf has 16) come back rounded to those nine; from then on
    # the text is a fixed point.
    conf = Configuration.from_string(path.read_text())
    again = Configuration.from_string(str(conf))

    assert again.velocities.shape == conf.velocities.shape
    assert np.allclose(again.velocities, conf.velocities, rtol=0, atol=5e-10)
    assert np.allclose(again.positions, conf.positions, rtol=0, atol=5e-10)
    assert str(Configuration.from_string(str(again))) == str(again)
    assert np.array_equal(Configuration.from_string(str(again)).velocities, again.velocities)


def test_a_configuration_without_velocities_writes_no_velocity_block():
    # an energy minimisation writes no velocities
    path = REF_DIR / "aladip_vacuum_em" / "expected" / "final.conf"
    conf = Configuration.from_string(path.read_text())
    assert conf.velocities.size == 0
    assert "VELOCITY" not in str(conf)


# ── System: a solute topology with a solvated configuration ─────────────────────────────────


@pytest.fixture(params=["from_files", "constructor"])
def aladip_system(request) -> System:
    topo, conf = str(SHARED / "aladip.topo"), str(SHARED / "aladip.conf")
    if request.param == "from_files":
        return System.from_files(topo, conf)
    return System(gromos.Topology(topo), Configuration(conf))


def test_system_solvates_a_solute_topology_to_the_configuration(aladip_system: System):
    # aladip.topo: 12 solute atoms, 3-atom SPC solvent; aladip.conf: 72 atoms (12 + 20 x 3)
    assert aladip_system.n_atoms == 72 == len(aladip_system.positions)
    assert aladip_system.topology.n_atoms == 72
    assert aladip_system.topology.n_solute_atoms == 12
    assert aladip_system.topology.n_solvent_atoms == 60


def test_simulation_of_a_solvated_system_has_the_same_atoms(aladip_system: System):
    recipe = Recipe.from_imd(str(REF_DIR / "aladip_solvated" / "aladip_solvated.in"))
    sim = Simulation(aladip_system, recipe)
    assert sim.n_atoms == 72 == aladip_system.n_atoms


def test_system_still_refuses_a_genuine_mismatch():
    with pytest.raises(ValueError, match="must match"):
        System.from_files(
            str(REF_DIR / "water_216_box" / "water_216_box.topo"),
            str(REF_DIR / "nacl_pair" / "nacl_pair.conf"),
        )
