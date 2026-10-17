"""
The first new term through the new door — 3.9 step 5.

`Term("xtb", coupling="delta", ...)` runs a real GFN-xTB subprocess over a region, additively on
top of the classical force field, registered like any other term: a `TermSpec` variant, its
registry lines and one `instantiate` arm. No cargo feature. This file is the physics oracle from
the Python side, the counterpart of `crates/newen-md/tests/xtb_orchestrator_sequence.rs`:

* adding the term to a recipe adds **exactly** `XtbPotential`'s direct energy (and forces) of the
  region at step 0 — `E(recipe + term) - E(recipe) == E_xtb(region)`;
* the term's forces conserve energy through the real step loop (NVE drift, same thresholds as the
  Rust oracle: fluctuation < 0.5 % of the mean, first/second-half drift < 0.2 %);
* `coupling="replace"` is rejected with the 2.8 message; a barostat is rejected (no
  virial); a region the `elements` table does not cover is a `RecipeError`; an `xtb` call that
  exceeds `timeout_s` is a `RunError`, not a hang.

The system is the water dimer of the gromosXX QM/MM references: molecule 1 (`1:a`) is a free,
unconstrained solute water carrying the xtb term, molecule 2 a rigid SPC solvent water. Skips
when `xtb` is not on PATH.
"""

from __future__ import annotations

import shutil

import numpy as np
import pytest
from test_gromosXX_references import REPO_ROOT

from newen import Recipe, Simulation, System, Term, XtbPotential, terms
from newen.exceptions import PlanError, RecipeError, RunError

FIXTURE = (
    REPO_ROOT
    / "crates"
    / "newen-forces"
    / "tests"
    / "gromosXX_qmmm_references"
    / "water_dimer_mechst"
)
TOPO = str(FIXTURE / "water_dimer.top")
CONF = str(FIXTURE / "water_dimer.cnf")
ELEMENTS = [8, 1, 1, 8, 1, 1]  # atomic number per global atom index (XtbInteraction convention)
QM_ATOMS = slice(0, 3)  # molecule 1, the solute water

pytestmark = pytest.mark.skipif(shutil.which("xtb") is None, reason="xtb not found on PATH")


def _base_recipe(steps: int = 120) -> Recipe:
    """NVE in vacuum, dt = 0.1 fs (the water-dimer reference's own), the solute water held together
    by the xtb term alone (classical bonded terms off), one rigid SPC solvent water."""
    return Recipe.nve(dt=1e-4, steps=steps).update(
        system={"nsm": 1},
        boundary={"kind": "vacuum", "ndfmin": 6},
        forcefield={"bonds": False, "angles": False, "impropers": False, "dihedrals": False},
    )


def _non_kinetic(sim: Simulation) -> float:
    """Classical potential + terms (+ restraints): everything in the total but the kinetic part."""
    return sim.total_energy - sim.kinetic_energy


def _xtb_term(tmp_path, **overrides) -> Term:
    params = {"region": "1:a", "elements": ELEMENTS, "work_dir": str(tmp_path / "xtb-term")}
    params.update(overrides)
    return Term("xtb", **params)


def test_xtb_is_a_registered_term_without_a_feature():
    entry = next(d for d in terms() if d["kind"] == "xtb")
    assert entry["feature"] is None and entry["available"] is True
    term = Term("xtb", **entry["params"])
    assert term.available and term.feature is None
    assert term.to_dict()["gfn"] == 2 and term.to_dict()["timeout_s"] == 600
    with pytest.raises(RecipeError):
        Term("xtb", region="1:a", elements=[8, 1, 1], gfn2=True)  # a typo is an error


def test_adding_the_term_adds_exactly_its_direct_energy_and_forces(tmp_path):
    system = System.from_files(TOPO, CONF)
    base = _base_recipe()
    with_term = base.with_term(_xtb_term(tmp_path))
    plan = with_term.plan(system)
    kinds = plan.kinds
    assert kinds.index("orchestrator") == kinds.index("forcefield") + 1
    assert plan["orchestrator"].to_dict()["terms"][0]["kind"] == "xtb"

    classical = Simulation(system, base)
    coupled = Simulation(system, with_term)
    direct_energy, direct_forces = XtbPotential(
        str(tmp_path / "direct"), ELEMENTS[QM_ATOMS]
    ).evaluate(system.positions[QM_ATOMS])

    # The term is additive (coupling=delta): what it adds is its own isolated-region answer,
    # visible on its own (G10) and in the total; the classical potential is untouched.
    assert coupled.term_energies == {"xtb": pytest.approx(direct_energy, rel=1e-9)}
    assert classical.term_energies == {}
    # Step 0 has already run the velocity half-step, so compare the non-kinetic part: the
    # term's energy is what it adds, and the classical potential is untouched.
    assert _non_kinetic(coupled) - _non_kinetic(classical) == pytest.approx(direct_energy, rel=1e-9)
    assert coupled.potential_energy == classical.potential_energy
    added = coupled.forces - classical.forces
    assert np.allclose(added[QM_ATOMS], direct_forces, rtol=1e-9, atol=1e-9)
    assert np.allclose(added[3:], 0.0, atol=1e-12), (
        "the term must not touch atoms outside its region"
    )
    assert coupled.recipe.terms == [_xtb_term(tmp_path)]


def test_nve_drift_through_the_real_step_loop(tmp_path):
    """Same measure as `xtb_orchestrator_sequence.rs`: the term's forces are consistent with its
    energy, so total energy is conserved by the real `AlgorithmSequence` step loop."""
    system = System.from_files(TOPO, CONF)
    sim = Simulation(system, _base_recipe().with_term(_xtb_term(tmp_path)))
    frames = sim.run(120, ene_freq=1)
    total = frames[1:, 3]  # frame 0 is the state before the first step
    mean = total.mean()
    mean_abs = max(abs(mean), 1e-6)
    fluctuation = np.abs(total - mean).max() / mean_abs
    half = len(total) // 2
    half_drift = abs(total[half:].mean() - total[:half].mean()) / mean_abs
    print(
        f"xtb term NVE: mean_total={mean:.6f} kJ/mol, fluctuation={fluctuation * 100:.6f} %, "
        f"half drift={half_drift * 100:.6f} %, kinetic[-1]={frames[-1, 1]:.4f} kJ/mol"
    )
    assert np.all(np.isfinite(total))
    assert frames[-1, 1] > 0.0 and np.any(sim.positions != system.positions), "nothing moved"
    assert fluctuation < 0.005, (
        f"energy fluctuation {fluctuation * 100:.4f} % of the mean exceeds 0.5 %"
    )
    assert half_drift < 0.002, f"first/second-half mean energy differs by {half_drift * 100:.4f} %"


def test_two_xtb_terms_get_their_own_work_dirs_and_add_up(tmp_path):
    system = System.from_files(TOPO, CONF)
    base = _base_recipe()
    # No work_dir given: each term derives its own from its index, so they cannot collide.
    both = base.with_term(Term("xtb", region="1:a", elements=ELEMENTS)).with_term(
        Term("xtb", region="2:a", elements=ELEMENTS)
    )
    e_base = _non_kinetic(Simulation(system, base))
    sim = Simulation(system, both)
    e1, _ = XtbPotential(str(tmp_path / "d1"), ELEMENTS[0:3]).evaluate(system.positions[0:3])
    e2, _ = XtbPotential(str(tmp_path / "d2"), ELEMENTS[3:6]).evaluate(system.positions[3:6])
    # Each term is visible on its own (a compensating error between the two could not hide).
    assert sim.term_energies == {
        "xtb:0": pytest.approx(e1, rel=1e-9),
        "xtb:1": pytest.approx(e2, rel=1e-9),
    }
    assert _non_kinetic(sim) - e_base == pytest.approx(e1 + e2, rel=1e-9)


def test_replace_coupling_barostat_and_bad_elements_are_rejected(tmp_path):
    system = System.from_files(TOPO, CONF)
    replace = _base_recipe().with_term(_xtb_term(tmp_path, coupling="replace"))
    with pytest.raises(PlanError, match="2.8"):
        replace.plan(system)
    with pytest.raises(PlanError, match="2.8"):
        Simulation(system, replace)

    npt = Recipe.npt(dt=1e-4, steps=10, temperature=300.0, pressure=1.0).update(system={"nsm": 1})
    with pytest.raises(PlanError, match="virial"):
        npt.with_term(_xtb_term(tmp_path)).plan(system)

    short = _base_recipe().with_term(_xtb_term(tmp_path, elements=[8]))
    short.plan(
        system
    )  # the plan is fine: the table is checked against the region when instantiating
    with pytest.raises(RecipeError, match="elements"):
        Simulation(system, short)


def test_timeout_kills_the_subprocess_instead_of_hanging(tmp_path):
    system = System.from_files(TOPO, CONF)
    impatient = _base_recipe().with_term(_xtb_term(tmp_path, timeout_s=0))
    with pytest.raises(RunError, match="timeout"):
        Simulation(system, impatient)
