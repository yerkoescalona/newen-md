"""
`newen.tools` — the gromos++ system-preparation programs, called in process (3.10).

The property these tests pin: **there is one program, reachable two ways.** `tools.run(name, …)`
returns the text the command-line binary prints; the typed wrappers return objects built by the
same Rust code. A test that finds the two disagreeing has found a real defect, not a formatting
difference — the Rust side already proves `run`'s bytes equal the binary's
(`crates/newen-tools/tests/library_parity.rs`) and equal the typed verb's serialised result
(`typed_verbs.rs`); these tests carry that guarantee across the Python boundary.
"""

from pathlib import Path

import pytest

import newen
from newen import tools

REPO_ROOT = Path(__file__).parent.parent.parent
SHARED = REPO_ROOT / "crates" / "newen-md" / "tests" / "gromosXX_references" / "shared"
TOOLS_DATA = REPO_ROOT / "crates" / "newen-tools" / "tests" / "data"
WATER_1000 = (
    REPO_ROOT
    / "crates"
    / "newen-md"
    / "tests"
    / "gromosXX_references"
    / "water_1000_spc_gridcell"
    / "water_1000_spc.cnf"
)

ALADIP_TOP = SHARED / "aladip.topo"
ALADIP_CNF = SHARED / "aladip.conf"
CH4_TOP = SHARED / "ch4_spc.top"
MEOH_TOP = SHARED / "meoh_spc.top"
MEOH_CNF = TOOLS_DATA / "sim_box" / "meoh.cnf"
ALADIP_PDB = TOOLS_DATA / "pdb2g96" / "aladip.pdb"
MTB = TOOLS_DATA / "make_top" / "54a7.mtb"
IFP = TOOLS_DATA / "make_top" / "54a7.ifp"

pytestmark = pytest.mark.skipif(
    not ALADIP_TOP.exists(), reason="reference fixtures not in the checkout"
)


# ── the registry ─────────────────────────────────────────────────────────────


def test_every_program_is_listed_with_its_keys_and_status():
    rows = tools.programs()
    by_name = {r["name"]: r for r in rows}

    # the whole newen-tools crate, not a hand-picked subset (atominfo lives in
    # newen-analysis — one program, one crate, since 0.0.49)
    assert len(rows) == 33
    assert {"make_top", "com_top", "red_top", "pdb2g96", "sim_box", "ion", "check_top"} <= set(
        by_name
    )

    sim_box = by_name["sim_box"]
    assert sim_box["status"] == "verified"
    assert "minwall" in sim_box["known"]
    assert sim_box["usage"].startswith("sim_box")

    # a status a caller can act on: what is a real port and what only shares the name
    assert {r["status"] for r in rows} <= {
        "verified",
        "ported",
        "stand_in",
        "stub",
        "no_counterpart",
    }
    assert by_name["gch"]["status"] == "stand_in"
    assert by_name["ran_box"]["status"] == "stub"


def test_keywords_become_gromospp_arguments():
    assert tools.argv_from(topo="a.top", pbc="r") == ["@topo", "a.top", "@pbc", "r"]
    # a sequence is a multi-valued key, True is a flag
    assert tools.argv_from(minwall=[1.0, 1.5, 2.0]) == ["@minwall", "1.0", "1.5", "2.0"]
    assert tools.argv_from(write=True) == ["@write"]
    # False and None mean "not given", so an optional argument can be passed unconditionally
    assert tools.argv_from(write=False, thresh=None) == []
    assert tools.argv_from(topo=Path("a.top")) == ["@topo", "a.top"]


def test_an_unknown_program_is_an_error_naming_it():
    with pytest.raises(newen.ProgramError, match="no GROMOS program named 'nope'"):
        tools.run("nope")


def test_an_unknown_key_is_rejected_rather_than_ignored():
    with pytest.raises(newen.ProgramError, match="unknown argument"):
        tools.run("red_top", topo=str(ALADIP_TOP), natoms=6, nosuchkey="x")


# ── the text path ────────────────────────────────────────────────────────────


def test_run_returns_what_the_program_prints():
    text = tools.run("red_top", topo=str(ALADIP_TOP), atoms="1-6")
    assert text.startswith("TITLE")
    assert "SOLUTEATOM" in text
    # the reduced topology really is six atoms
    assert newen.Topology.from_string(text).n_solute_atoms == 6


def test_a_failing_program_raises_with_the_output_it_had_written():
    # check_top writes its whole report and *then* fails; the report is the answer
    with pytest.raises(newen.ProgramError) as excinfo:
        tools.run("check_top", topo=str(TOOLS_DATA / "no_such.top"))
    assert isinstance(excinfo.value.stdout, str)


def test_argv_is_accepted_verbatim_for_reproducing_a_command_line():
    by_kwargs = tools.run("red_top", topo=str(ALADIP_TOP), natoms=6)
    by_argv = tools.run("red_top", _argv=["@topo", str(ALADIP_TOP), "@natoms", "6"])
    assert by_kwargs == by_argv


# ── typed results, and that they are the same thing as the text ──────────────


def test_red_top_object_and_text_are_the_same_topology():
    top = newen.Topology.from_file(str(ALADIP_TOP))
    reduced = tools.red_top(top, last=6)
    text = tools.run("red_top", topo=str(ALADIP_TOP), atoms="1-6")

    assert reduced.n_solute_atoms == 6
    # only the TITLE differs: `run` names the input file, the object says newen
    assert str(reduced).splitlines()[3:] == text.splitlines()[3:]


def test_com_top_combines_and_accepts_paths_or_objects():
    a = newen.Topology.from_file(str(ALADIP_TOP))
    b = newen.Topology.from_file(str(CH4_TOP))
    combined = tools.com_top(a, b)
    assert combined.n_solute_atoms == a.n_solute_atoms + b.n_solute_atoms

    # a path where an object is expected, and a repeat count
    from_paths = tools.com_top(str(ALADIP_TOP), str(CH4_TOP))
    assert str(from_paths) == str(combined)

    twice = tools.com_top((2, a))
    assert twice.n_solute_atoms == 2 * a.n_solute_atoms


def test_red_top_rejects_a_range_outside_the_topology():
    top = newen.Topology.from_file(str(ALADIP_TOP))
    with pytest.raises(newen.ProgramError, match="outside"):
        tools.red_top(top, last=top.n_solute_atoms + 1)


@pytest.mark.skipif(not WATER_1000.exists(), reason="solvent box not in the checkout")
def test_sim_box_solvates_and_agrees_with_its_text():
    top = newen.Topology.from_file(str(MEOH_TOP))
    solute = newen.Configuration.from_file(str(MEOH_CNF))
    solvent = newen.Configuration.from_file(str(WATER_1000))

    solvated = tools.sim_box(top, solute, solvent=solvent, minwall=1.0)

    # 352 SPC around methanol, the gromos++ figure (crates/newen-tools/tests/data/sim_box)
    assert solvated.n_atoms == solute.n_atoms + 352 * 3
    box = solvated.box_dimensions
    assert box[0] == pytest.approx(2.207681950, abs=1e-8)

    text = tools.run(
        "sim_box",
        topo=str(MEOH_TOP),
        pbc="r",
        pos=str(MEOH_CNF),
        solvent=str(WATER_1000),
        minwall=[1.0],
    )

    # same atoms, same box; only the TITLE differs (`run` names the input files)
    def body(g96: str) -> list[str]:
        lines = g96.splitlines()
        return lines[lines.index("POSITION") :]

    assert body(solvated.to_g96("x")) == body(text)


def test_ion_replaces_a_water_and_leaves_the_topology_alone():
    top = newen.Topology.from_file(str(ALADIP_TOP))
    system = newen.Configuration.from_file(str(ALADIP_CNF))

    ionised = tools.ion(top, system, positive=(1, "NA+"), potential=1.4)

    # one 3-atom water became one ion
    assert ionised.n_atoms == system.n_atoms - 2
    assert "NA+" in str(ionised)
    # gromos++ leaves the topology to the caller, and so do we
    assert top.n_solute_atoms == newen.Topology.from_file(str(ALADIP_TOP)).n_solute_atoms


def test_ion_without_a_potential_cutoff_is_an_error():
    top = newen.Topology.from_file(str(ALADIP_TOP))
    system = newen.Configuration.from_file(str(ALADIP_CNF))
    with pytest.raises(newen.ProgramError, match="@potential"):
        tools.ion(top, system, positive=(1, "NA+"))


def test_build_box_fills_a_box_to_a_density():
    top = newen.Topology.from_file(str(MEOH_TOP))
    molecule = newen.Configuration.from_file(str(MEOH_CNF))

    box = tools.build_box(top, molecule, nsm=2, dens=790.0)
    assert box.n_atoms == molecule.n_atoms * 8  # 2^3 copies
    assert box.box_dimensions[0] > 0.0


@pytest.mark.skipif(not ALADIP_PDB.exists(), reason="PDB fixture not in the checkout")
def test_pdb2g96_reorders_a_pdb_into_the_topology_order():
    top = newen.Topology.from_file(str(ALADIP_TOP))
    conf = tools.pdb2g96(top, pdb=str(ALADIP_PDB))
    assert conf.n_atoms == top.n_solute_atoms
    # every atom found a match, so no position is the (0,0,0) placeholder
    assert not (conf.positions == 0.0).all(axis=1).any()

    # strict turns unmatched atoms into an error instead of silent zeros; this fixture has none
    tools.pdb2g96(top, pdb=str(ALADIP_PDB), strict=True)


@pytest.mark.skipif(not MTB.exists(), reason="force-field files not in the checkout")
def test_make_top_builds_a_topology_and_agrees_with_its_text():
    top = tools.make_top(
        build=str(MTB), param=str(IFP), seq=["NH3+", "ALA", "GLY", "COO-"], solv="H2O"
    )
    assert top.n_solute_atoms > 0
    assert top.masses.shape == (top.n_atoms,)

    text = tools.run(
        "make_top",
        build=str(MTB),
        param=str(IFP),
        seq=["NH3+", "ALA", "GLY", "COO-"],
        solv="H2O",
    )
    # the object is the text, parsed: same atoms and same masses
    assert newen.Topology.from_string(text).n_solute_atoms == top.n_solute_atoms


@pytest.mark.skipif(not MTB.exists(), reason="force-field files not in the checkout")
def test_make_top_rejects_an_unknown_solvent_without_writing_a_topology():
    with pytest.raises(newen.ProgramError, match="NO_SUCH_SOLVENT"):
        tools.make_top(build=str(MTB), param=str(IFP), seq=["ALA"], solv="NO_SUCH_SOLVENT")


def test_check_top_reports_rather_than_failing():
    report = tools.check_top(str(ALADIP_TOP))
    assert isinstance(report.errors, list)
    assert isinstance(report.warnings, list)
    assert report.text
    # falsy when the topology has errors, so `if not report:` reads the obvious way
    assert bool(report) == (not report.errors)
    assert str(report) == report.text


# ── the chain, with no files in between ──────────────────────────────────────


@pytest.mark.skipif(
    not (WATER_1000.exists() and MTB.exists()), reason="fixtures not in the checkout"
)
def test_the_preparation_chain_passes_objects_not_paths():
    """make_top → sim_box → System, with no file anywhere between the steps."""
    meoh = newen.Topology.from_file(str(MEOH_TOP))
    conf = newen.Configuration.from_file(str(MEOH_CNF))
    solvent = newen.Configuration.from_file(str(WATER_1000))

    solvated = tools.sim_box(meoh, conf, solvent=solvent, minwall=1.0)

    # the topology follows the coordinates, then the two make a System
    n_solvent_atoms = solvated.n_atoms - meoh.n_solute_atoms
    meoh.solvate(n_solvent_atoms // 3)
    system = newen.System(meoh, solvated)
    assert system.n_atoms == solvated.n_atoms


def test_ion_leaves_the_topology_one_atom_short_until_the_ion_is_added():
    """The gromos++ contract, made visible.

    `ion` swaps a 3-atom water for a 1-atom ion in the *coordinates* only — the ion's own
    topology is the caller's to add (`com_top` with the ion's `.top`). Until then the topology
    and the configuration disagree, and `System` says so rather than running the wrong physics.
    `System.add_ions` does both halves.
    """
    meoh = newen.Topology.from_file(str(MEOH_TOP))
    conf = newen.Configuration.from_file(str(MEOH_CNF))
    solvent = newen.Configuration.from_file(str(WATER_1000))

    solvated = tools.sim_box(meoh, conf, solvent=solvent, minwall=1.0)
    ionised = tools.ion(meoh, solvated, positive=(1, "NA+"), potential=1.4)

    assert ionised.n_atoms == solvated.n_atoms - 2  # −3 water atoms, +1 ion

    meoh.solvate((ionised.n_atoms - meoh.n_solute_atoms - 1) // 3)
    with pytest.raises(newen.RecipeError, match="must match"):
        newen.System(meoh, ionised)


def test_a_configuration_round_trips_through_its_own_text():
    conf = newen.Configuration.from_file(str(ALADIP_CNF))
    again = newen.Configuration.from_string(str(conf))
    assert again.n_atoms == conf.n_atoms
    assert again.positions == pytest.approx(conf.positions)
    assert again.box_dimensions == pytest.approx(conf.box_dimensions)


def test_a_topology_round_trips_through_its_own_text():
    top = newen.Topology.from_file(str(ALADIP_TOP))
    again = newen.Topology.from_string(str(top))
    assert again.n_solute_atoms == top.n_solute_atoms
    assert again.masses == pytest.approx(top.masses)
    assert again.charges == pytest.approx(top.charges)


# ── System verbs: the same programs, as methods ──────────────────────────────


@pytest.mark.skipif(not WATER_1000.exists(), reason="solvent box not in the checkout")
def test_system_solvate_equals_the_module_function():
    top = newen.Topology.from_file(str(MEOH_TOP))
    conf = newen.Configuration.from_file(str(MEOH_CNF))
    solvent = newen.Configuration.from_file(str(WATER_1000))

    system = newen.System(top, conf)
    solvated = system.solvate(solvent, minwall=[1.0])
    by_function = tools.sim_box(top, conf, solvent=solvent, minwall=1.0)

    # the same coordinates, to the last bit — one verb, called two ways
    assert solvated.positions == pytest.approx(by_function.positions, abs=0.0)
    assert solvated.n_atoms == by_function.n_atoms
    # and the topology now matches, which the loose configuration cannot promise
    assert solvated.topology.n_atoms == solvated.n_atoms


@pytest.mark.skipif(not WATER_1000.exists(), reason="solvent box not in the checkout")
def test_system_solvate_leaves_the_input_system_untouched():
    top = newen.Topology.from_file(str(MEOH_TOP))
    conf = newen.Configuration.from_file(str(MEOH_CNF))
    solvent = newen.Configuration.from_file(str(WATER_1000))

    system = newen.System(top, conf)
    before = system.n_atoms
    system.solvate(solvent, minwall=[1.0])
    assert system.n_atoms == before


@pytest.mark.skipif(
    not (WATER_1000.exists() and MTB.exists()), reason="fixtures not in the checkout"
)
def test_system_add_ions_updates_the_topology_too():
    """What `tools.ion` deliberately leaves to the caller, `System.add_ions` does."""
    top = newen.Topology.from_file(str(MEOH_TOP))
    conf = newen.Configuration.from_file(str(MEOH_CNF))
    solvent = newen.Configuration.from_file(str(WATER_1000))

    solvated = newen.System(top, conf).solvate(solvent, minwall=[1.0])
    ionised = solvated.add_ions(forcefield=(str(MTB), str(IFP)), positive=(1, "NA+"), potential=1.4)

    # one 3-atom water became one ion, and the topology knows about it
    assert ionised.n_atoms == solvated.n_atoms - 2
    assert ionised.topology.n_atoms == ionised.n_atoms
    # the coordinates are what tools.ion produces
    by_function = tools.ion(
        solvated.topology, solvated.configuration, positive=(1, "NA+"), potential=1.4
    )
    assert ionised.positions == pytest.approx(by_function.positions, abs=0.0)


def test_system_check_is_the_check_top_report():
    system = newen.System.from_files(str(ALADIP_TOP), str(ALADIP_CNF))
    errors, warnings, text = system.check()
    report = tools.check_top(str(ALADIP_TOP))
    assert errors == report.errors
    assert warnings == report.warnings
    assert text == report.text
