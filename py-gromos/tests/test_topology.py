"""
`Topology.iac` and `Topology.lj_parameters`: the Lennard-Jones parameters as the engine holds
them, read from the topology rather than typed in by hand.
"""

from pathlib import Path

import pytest

from gromos import Topology

REF_DIR = (
    Path(__file__).parent.parent.parent / "crates" / "gromos-md" / "tests" / "gromosXX_references"
)


def test_pair_lj_parameters_are_the_file_values():
    # pair_lj.topo: two argon atoms of IAC 1 (file, 1-based); LJPARAMETERS row
    # "1 1  C12=9.847000e-06  C6=6.264700e-03  CS12=9.847000e-06  CS6=6.264700e-03"
    topo = Topology(str(REF_DIR / "pair_lj" / "pair_lj.topo"))

    assert len(topo.iac) == topo.n_atoms == 2
    assert len(topo.atom_names) == topo.n_atoms
    assert topo.atom_names == ["AR1", "AR2"]
    assert topo.iac.dtype.kind == "i"
    assert topo.iac[0] == topo.iac[1] == 0  # 0-based: file IAC 1 is stored as 0

    c6, c12, cs6, cs12 = topo.lj_parameters(topo.iac[0], topo.iac[1])
    assert c6 == 6.2647e-3
    assert c12 == 9.847e-6
    assert (cs6, cs12) == (6.2647e-3, 9.847e-6)


def test_lj_parameters_rejects_an_unknown_type():
    topo = Topology(str(REF_DIR / "pair_lj" / "pair_lj.topo"))
    with pytest.raises(IndexError, match="out of range"):
        topo.lj_parameters(0, 1)


def test_iac_and_names_cover_the_solvent_once_solvated():
    topo = Topology(str(REF_DIR / "shared" / "aladip.topo"))
    topo.solvate(20)
    assert len(topo.iac) == len(topo.atom_names) == topo.n_atoms == 72
    assert topo.atom_names[12:15] == ["OW", "HW1", "HW2"]
    assert topo.atom_names[-3:] == ["OW", "HW1", "HW2"]


def test_atom_names_of_a_many_molecule_solute():
    # water_216_box.topo: 216 waters, all solute — every molecule is named, not only the first
    topo = Topology(str(REF_DIR / "water_216_box" / "water_216_box.topo"))
    assert len(topo.atom_names) == topo.n_atoms == 648
    assert topo.atom_names == ["OW", "HW1", "HW2"] * 216
