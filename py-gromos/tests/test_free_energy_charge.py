"""
A charge-changing perturbation end to end, all inside gromos-rs: ``ch4_water_fep``'s CH4 gains a
unit charge in state B (``Perturbation.mutate``, i.e. gromos-tools' ``make_pt_top``), seven λ
windows run in worker processes (``run_windows(parallel=...)``), TI, BAR and MBAR agree, and the
finite-size correction of the charge change comes from ``dGslv_pbsolv`` on the two end states
(``Perturbation.end_state``, i.e. ``pt_top``), whose numbers gromos++ computed for the same state
(``crates/gromos-analysis/tests/data/dgslv_pbsolv/README.md``, cases ``ch4q_*``).
"""

import pickle
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import gromos
from gromos import Recipe, System, Topology, analysis
from gromos.free_energy import LambdaSchedule, Perturbation, run_windows

REPO_ROOT = Path(__file__).parent.parent.parent
REF = REPO_ROOT / "crates" / "gromos-md" / "tests" / "gromosXX_references"
SHARED = REF / "shared"
PB = REPO_ROOT / "crates" / "gromos-analysis" / "tests" / "data" / "dgslv_pbsolv"

pytestmark = pytest.mark.skipif(not REF.exists(), reason="reference systems not in the checkout")

# seven windows λ = k/6 on a seven-point PRECALCLAM grid, 600 steps each at 300 K (velocities
# generated at 300 K, Berendsen τ = 0.1 ps), a frame every 5 steps; the first 20 frames dropped
LAMBDAS = [k / 6 for k in range(7)]
SKIP = 20

# gromos++ dGslv_pbsolv on the charged CH4 (ch4q_rf.log, ch4q_ls.log): the DGRESULT lines and
# the free energies that end the two FFT sections, kJ/mol, printed to 6 significant digits
GROMOSPP_RF = {
    "dg_npbc": -296.549,
    "dg_pbc": -235.623,
    "dg_fft_ls_pbc": -206.065,
    "dg_fft_rf_pbc": -189.005,
}
PB_TOLERANCE = 0.01


def ch4_system() -> System:
    return System.from_files(str(SHARED / "ch4_spc.top"), str(SHARED / "ch4_spc.cnf"))


def base_recipe(steps: int, nr_lambdas: int) -> Recipe:
    return Recipe.from_imd(str(REF / "ch4_water_fep" / "ch4_water_fep.in")).update(
        control={"steps": steps, "initial_temperature": 300.0},
        outputs={"ntwg": 5, "energy_every": 5},
        ensemble={"thermostat": {"baths": [{"tau": 0.1, "temperature": 300.0}]}},
        perturbation={
            "precalclam": {"nr_lambdas": nr_lambdas, "min_lambda": 0.0, "max_lambda": 1.0}
        },
    )


def unit_charge(system: System) -> Perturbation:
    return Perturbation.mutate(system, "1:1", charge_b=1.0, alpha_lj=0.0, alpha_crf=0.0)


@pytest.fixture(scope="module")
def charging(tmp_path_factory):
    system = ch4_system()
    plus = unit_charge(system)
    recipes = LambdaSchedule(LAMBDAS).recipes(base_recipe(600, len(LAMBDAS)))
    workdir = tmp_path_factory.mktemp("charging")
    windows = run_windows(system, recipes, perturbation=plus, workdir=workdir, parallel=4)
    return system, plus, windows, workdir


def test_mutate_is_make_pt_top_and_reports_the_net_charges(tmp_path):
    system = ch4_system()
    plus = unit_charge(system)
    assert (plus.net_charge_a, plus.net_charge_b, plus.delta_charge) == (0.0, 1.0, 1.0)
    # the fixture the gromos++ oracle ran on was written by this call
    assert plus.text == (PB / "ch4_plus.ptp").read_text()
    kept = Perturbation.from_ptp(PB / "ch4_plus.ptp", system)
    assert (kept.text, kept.net_charge_a, kept.net_charge_b) == (plus.text, 0.0, 1.0)
    assert kept.path == PB / "ch4_plus.ptp"
    # without a topology the charges are unknown; the same answer as Topology.apply_perturbation
    assert Perturbation.from_ptp(SHARED / "ch4_spc_dummy.ptp").delta_charge is None
    assert Perturbation.from_ptp(SHARED / "ch4_spc_dummy.ptp", system).delta_charge == 0.0
    written = plus.write(tmp_path / "plus.ptp")
    assert written.path == tmp_path / "plus.ptp"
    top = Topology(str(SHARED / "ch4_spc.top"))
    assert top.apply_perturbation(str(written.path)) == (0.0, 1.0)

    # the same atom by index; type and mass; one value per atom
    by_index = Perturbation.mutate(system.topology, [0], charge_b=[1.0], alpha_lj=0, alpha_crf=0)
    assert by_index.text == plus.text
    dummy = Perturbation.mutate(system, "1:CM", iac_b=22, mass_b=12.011, alpha_lj=0.5)
    row = dummy.text.splitlines()[dummy.text.splitlines().index("PERTATOMPARAM") + 4].split()
    assert row == ["1", "1", "CM", "17", "16.04300", "0.00000", "22", "12.01100", "0.00000",
                   "0.50", "0.50"]  # fmt: skip
    assert dummy.delta_charge == 0.0

    refusals: list[tuple[dict[str, Any], str]] = [
        ({"atoms": [3], "charge_b": 1.0}, "not a solute atom"),
        ({"atoms": "1:1"}, "nothing to change"),
        ({"atoms": "1:1", "iac_b": 99}, "atom types"),
        ({"atoms": "1:1", "iac_b": 0}, "start at 1"),
        ({"atoms": "1:1", "charge_b": [1.0, 2.0]}, "one per atom"),
        ({"atoms": "1:1", "mass_b": -1.0}, "positive"),
    ]
    for bad, match in refusals:
        with pytest.raises(gromos.ProgramError, match=match):
            Perturbation.mutate(system, **bad)


def test_end_states_are_pt_top():
    system = ch4_system()
    plus = unit_charge(system)
    a = plus.end_state(system, "A")
    b = plus.end_state(system)
    assert list(a.charges) == [0.0] and list(b.charges) == [1.0]
    assert list(b.masses) == list(a.masses)
    # the topology gromos++ read for the oracle (`pt_top` binary) holds the same state
    oracle = Topology(str(PB / "ch4q.top"))
    assert list(oracle.charges) == list(b.charges) and list(oracle.masses) == list(b.masses)
    with pytest.raises(ValueError, match="'A' or 'B'"):
        plus.end_state(system, "C")


def test_parallel_windows_are_the_sequential_windows(tmp_path):
    system = ch4_system()
    plus = unit_charge(system)
    recipes = LambdaSchedule([0.0, 0.5, 1.0]).recipes(base_recipe(30, 3))
    seq = run_windows(system, recipes, perturbation=plus)
    par = run_windows(system, recipes, perturbation=plus, parallel=3, workdir=tmp_path)
    assert (tmp_path / "perturbation.ptp").read_text() == plus.text
    for a, b in zip(seq, par, strict=True):
        assert a.lam == b.lam and a.energies is not None and b.energies is not None
        assert np.array_equal(a.energies.array, b.energies.array)
        assert a.precalclam is not None and b.precalclam is not None
        for key in ("energies", "derivatives", "abdih", "free_abdih"):
            assert np.array_equal(getattr(a.precalclam, key), getattr(b.precalclam, key)), key
    assert tuple(seq.bar()) == tuple(par.bar()) and seq.ti() == par.ti()
    # kept windows are read, not sent to workers again
    again = run_windows(system, recipes, perturbation=plus, parallel=3, workdir=tmp_path)
    assert all(np.array_equal(a.dhdl, b.dhdl) for a, b in zip(par, again, strict=True))


def test_parallel_needs_a_system_that_pickles(tmp_path):
    system = ch4_system()
    assert pickle.loads(pickle.dumps(system)).n_atoms == system.n_atoms
    in_memory = System(system.topology, system.configuration)
    with pytest.raises(TypeError, match="from_files"):
        pickle.dumps(in_memory)
    recipes = LambdaSchedule([0.0, 1.0]).recipes(base_recipe(10, 2))
    with pytest.raises(TypeError, match="worker processes"):
        run_windows(in_memory, recipes, perturbation=unit_charge(system), parallel=2)
    with pytest.raises(ValueError, match="at least 1"):
        run_windows(system, recipes, parallel=0)
    named = [r.with_inputs(pttopo=str(SHARED / "ch4_spc_dummy.ptp")) for r in recipes]
    with pytest.raises(gromos.RecipeError, match="names a perturbation topology"):
        run_windows(system, named, perturbation=unit_charge(system))


def test_ti_bar_and_mbar_of_the_charging_agree(charging):
    _, plus, windows, workdir = charging
    assert windows.lambdas == LAMBDAS
    assert (workdir / "perturbation.ptp").read_text() == plus.text
    ptp = str(workdir / "perturbation.ptp")
    for w in windows:
        assert w.recipe is not None and w.recipe.inputs["pttopo"] == ptp
        assert len(w.dhdl) == 120
    mean = windows.mean_dhdl(SKIP)
    # charging lowers the energy more the more charge there is: ⟨∂H/∂λ⟩ falls with λ
    assert all(b < a < 0 for a, b in zip(mean, mean[1:], strict=False))
    ti = windows.ti(SKIP)
    bar = tuple(windows.bar(SKIP))
    mbar = tuple(windows.mbar(SKIP))
    estimates = {"TI": ti, "BAR": bar, "MBAR": mbar}
    for a, (ga, ea) in estimates.items():
        assert np.isfinite(ga) and ea > 0, a
        assert -300.0 < ga < -200.0, a
        for b, (gb, eb) in estimates.items():
            assert abs(ga - gb) <= 2 * np.hypot(ea, eb), f"{a} {ga} ± {ea} vs {b} {gb} ± {eb}"


def _dgslv(topology: Path, recipe: Recipe) -> dict:
    """``dGslv_pbsolv @schemeELEC RF`` with the run's reaction field, the README's grid."""
    es = recipe.forcefield["electrostatics"]
    return analysis.dgslv_pbsolv(
        topo=topology,
        coord=SHARED / "ch4_spc.cnf",
        pbc="r",
        atoms="1:a",
        atoms_to_charge="1:a",
        scheme="RF",
        eps_solv=66.6,
        eps_rf=es["epsilon_rf"],
        rcut=es["rf_cutoff"],
        gridspacing=0.1,
        npbc_size=5.0,
        probe_iac=5,
    )


def test_the_charge_change_is_corrected_with_dgslv_pbsolv(charging, tmp_path):
    system, plus, windows, _ = charging
    recipe = windows[0].recipe
    assert recipe is not None
    corrections = {}
    for state in ("A", "B"):
        path = tmp_path / f"state_{state}.top"
        plus.end_state(system, state).write(str(path))
        corrections[state] = _dgslv(path, recipe)
    a, b = corrections["A"], corrections["B"]
    # state A carries no charge: nothing to correct
    assert a["documented_correction"] == 0.0
    # state B: gromos++'s numbers for the same state, grid and settings
    for key, value in GROMOSPP_RF.items():
        assert b[key] == pytest.approx(value, abs=PB_TOLERANCE), key
    oracle = (
        GROMOSPP_RF["dg_npbc"]
        - GROMOSPP_RF["dg_pbc"]
        + GROMOSPP_RF["dg_fft_ls_pbc"]
        - GROMOSPP_RF["dg_fft_rf_pbc"]
    )
    correction = b["documented_correction"] - a["documented_correction"]
    assert correction == pytest.approx(oracle, abs=PB_TOLERANCE)
    assert correction == pytest.approx(-77.985, abs=PB_TOLERANCE)

    # the corrected charging free energy, by each estimator
    for dg, _err in (windows.ti(SKIP), tuple(windows.bar(SKIP)), tuple(windows.mbar(SKIP))):
        corrected = dg + correction
        assert -380.0 < corrected < -270.0
