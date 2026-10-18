"""
Front-end parity — 3.9, drift guard G2.

The reference suite (`test_gromosXX_references.py`) guards the Python front-end against gromosXX;
this file guards the front-ends against *each other* — the trap named when the composition redesign was parked (see `.claude/architecture.md`), which parked the
builder unification. Every path below must be **exactly** equal (``np.array_equal``, never
``allclose``) on every reference system: per-step energy rows, final positions, final forces.
Same builder, same kernels, same order — any difference is a defect, never a tolerance question
(3.9 A1/A3).

  path A — ``Simulation(topo, conf, InputParameters, **restraint_kwargs)``   (deprecated shim)
  path B — ``Simulation(system, Recipe.from_imd(...).with_inputs(...))``      (the front-end)
  path D — ``Simulation(system, recipe, plan=recipe.plan(system))``           (the plan front-end)
  shim   — ``Simulation.from_sequence(topo, conf, params, AlgorithmSequence.from_parameters(...))``

Path C (the descriptor resolver, a second builder) was deleted in step 4 together with its
``xfail`` table: the shims that replaced it are translations into a recipe, and the tests below
prove it. ``test_same_path_twice`` is the determinism baseline (A3). All paths run with the
binary's ``ParallelPolicy::Auto`` — the same kernels at the same thread count; run-to-run
determinism at a fixed thread count is measured by ``scripts/kernel_determinism.py``
(BENCHMARKING.md).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
from test_gromosXX_references import (
    REF_DIR,
    REFERENCE_SYSTEMS,
    _get_n_steps,
    _parse_input_toml,
)

from newen import (
    AlgorithmSequence,
    Configuration,
    InputParameters,
    Recipe,
    Simulation,
    System,
    Topology,
    algorithms,
    terms,
)

# the slow-growth references that take λ past 1 warn, and so does
# ch4_water_fep_precalc_twinrange (see test_gromosXX_references.py); the plain Nose-Hoover
# references run here without their `nh_restart` override, so NTINHT 0 finds no NHCVARIABLES
pytestmark = [
    pytest.mark.filterwarnings(r"ignore:slow growth takes lambda out of \[0, 1\]"),
    pytest.mark.filterwarnings(r"ignore:PRECALCLAM with a twin-range pairlist"),
    pytest.mark.filterwarnings(r"ignore:INITIALISE NTINHT = 0 asks for the Nose-Hoover"),
]


def _load(system_name: str):
    system_dir = REF_DIR / system_name
    inputs = _parse_input_toml(system_dir)
    topo = Topology(str((system_dir / inputs["topology"]).resolve()))
    if inputs.get("pttopo"):
        topo.apply_perturbation(str((system_dir / inputs["pttopo"]).resolve()))
    conf = Configuration(str((system_dir / inputs["configuration"]).resolve()))
    params = InputParameters(str((system_dir / inputs["parameters"]).resolve()))
    kwargs = {}
    for key in ("distrest", "posresspec", "refpos"):
        if inputs.get(key):
            kwargs[key] = str((system_dir / inputs[key]).resolve())
    return topo, conf, params, kwargs


def _trace(sim: Simulation, n_steps: int):
    """Per-step energy rows (frame 0 = state after step 0), final positions and forces."""
    frames = sim.run(n_steps, ene_freq=1)
    return frames, sim.positions, sim.forces


def _run_path_a(system_name: str):
    topo, conf, params, kwargs = _load(system_name)
    sim = Simulation(topo, conf, params, **kwargs)
    return _trace(sim, _get_n_steps(REF_DIR / system_name))


def _recipe_and_system(system_name: str) -> tuple[Recipe, System]:
    """Path B's inputs: the recipe of the reference `.imd`, auxiliary files as `inputs`."""
    system_dir = REF_DIR / system_name
    inputs = _parse_input_toml(system_dir)
    system = System.from_files(
        str((system_dir / inputs["topology"]).resolve()),
        str((system_dir / inputs["configuration"]).resolve()),
    )
    aux = {
        key: str((system_dir / inputs[key]).resolve())
        for key in ("pttopo", "distrest", "posresspec", "refpos")
        if inputs.get(key)
    }
    recipe = Recipe.from_imd(str((system_dir / inputs["parameters"]).resolve())).with_inputs(**aux)
    return recipe, system


def _run_path_b(system_name: str):
    """Path B — `Simulation(system, Recipe.from_imd(...).with_inputs(...))` (step 3)."""
    recipe, system = _recipe_and_system(system_name)
    return _trace(Simulation(system, recipe), _get_n_steps(REF_DIR / system_name))


def _run_path_d(system_name: str):
    """Path D — the recipe's own (unedited) `Plan` handed back: `plan=` is stage 1 of the
    same builder, not a second one."""
    recipe, system = _recipe_and_system(system_name)
    plan = recipe.plan(system)
    return _trace(Simulation(system, recipe, plan=plan), _get_n_steps(REF_DIR / system_name))


def _run_path_shim(system_name: str):
    """The deprecated `AlgorithmSequence.from_parameters` + `Simulation.from_sequence` names,
    now a translation: `params` -> recipe, the preset -> `recipe.plan`. They carry no restraint
    inputs, so the restraint systems skip here (path B covers them)."""
    topo, conf, params, kwargs = _load(system_name)
    if kwargs:
        pytest.skip(f"{system_name}: the deprecated from_sequence shim has no restraint inputs")
    plan = AlgorithmSequence.from_parameters(topo, params)
    sim = Simulation.from_sequence(topo, conf, params, plan)
    return _trace(sim, _get_n_steps(REF_DIR / system_name))


def _first_difference(a: np.ndarray, b: np.ndarray) -> str:
    """Human-readable location of the first exact mismatch (for the divergence table)."""
    if a.shape != b.shape:
        return f"shape {a.shape} vs {b.shape}"
    idx = np.argwhere(a != b)
    if idx.size == 0:
        return "identical"
    first = tuple(int(i) for i in idx[0])
    return f"first mismatch at index {first}: {a[first]!r} vs {b[first]!r}"


def _assert_identical(label: str, a, b) -> None:
    fa, pa, xa = a
    fb, pb, xb = b
    problems = []
    if not np.array_equal(fa, fb):
        problems.append(f"energy rows — {_first_difference(fa, fb)}")
    if not np.array_equal(pa, pb):
        problems.append(f"final positions — {_first_difference(pa, pb)}")
    if not np.array_equal(xa, xb):
        problems.append(f"final forces — {_first_difference(xa, xb)}")
    assert not problems, f"{label}:\n  " + "\n  ".join(problems)


@pytest.mark.parametrize("system_name", [pytest.param(s, id=s) for s in REFERENCE_SYSTEMS])
def test_same_path_twice(system_name):
    """Determinism baseline (A3): the same construction path twice must be bit-identical."""
    _assert_identical(
        f"{system_name}: path A vs path A", _run_path_a(system_name), _run_path_a(system_name)
    )


@pytest.mark.parametrize("system_name", [pytest.param(s, id=s) for s in REFERENCE_SYSTEMS])
def test_recipe_front_end_parity(system_name):
    """Path A (deprecated `InputParameters` + restraint kwargs) vs path B (`Recipe`): the
    deprecation shim is a translation into the recipe, so the two must be bit-identical."""
    _assert_identical(
        f"{system_name}: path A (InputParameters) vs path B (Recipe)",
        _run_path_a(system_name),
        _run_path_b(system_name),
    )


@pytest.mark.parametrize("system_name", [pytest.param(s, id=s) for s in REFERENCE_SYSTEMS])
def test_from_sequence_shim_parity(system_name):
    """The deprecated `AlgorithmSequence`/`from_sequence` names vs path B."""
    _assert_identical(
        f"{system_name}: from_sequence shim vs path B (Recipe)",
        _run_path_shim(system_name),
        _run_path_b(system_name),
    )


@pytest.mark.parametrize("system_name", [pytest.param(s, id=s) for s in REFERENCE_SYSTEMS])
def test_plan_front_end_parity(system_name):
    """Path B (`Recipe`) vs path D (`Recipe` + its own `Plan` passed back through `plan=`)."""
    _assert_identical(
        f"{system_name}: path B (Recipe) vs path D (Recipe + Plan)",
        _run_path_b(system_name),
        _run_path_d(system_name),
    )


# ---------------------------------------------------------------------------------------------
# Drift guards on the registries (G5/G6): every kind the engine knows is covered by a parity
# case and named in the stubs. An exemption below must say which step removes it.
# ---------------------------------------------------------------------------------------------
UNCOVERED_KINDS: dict[str, str] = {
    "orchestrator": "needs a Term (xtb/schnet); no gromosXX reference has one — 3.9 step 5",
}

PYI = Path(__file__).resolve().parents[1] / "python" / "newen" / "newen.pyi"


def test_every_kind_has_a_parity_case():
    registry = {d["kind"] for d in algorithms()}
    covered: set[str] = set()
    for system_name in REFERENCE_SYSTEMS:
        recipe, system = _recipe_and_system(system_name)
        covered |= set(recipe.plan(system).kinds)
    stale = set(UNCOVERED_KINDS) - registry
    assert not stale, f"UNCOVERED_KINDS names kinds the engine no longer has: {sorted(stale)}"
    no_longer_needed = set(UNCOVERED_KINDS) & covered
    assert not no_longer_needed, (
        f"remove from UNCOVERED_KINDS, now covered: {sorted(no_longer_needed)}"
    )
    missing = registry - covered - set(UNCOVERED_KINDS)
    assert not missing, f"algorithm kinds no reference system exercises: {sorted(missing)}"


def _literal(name: str) -> set[str]:
    m = re.search(rf"^{name} = Literal\[(.*?)\]", PYI.read_text(), re.S | re.M)
    assert m, f"{name} not defined in {PYI}"
    return set(re.findall(r'"([^"]+)"', m.group(1)))


def test_pyi_lists_every_kind():
    assert _literal("_AlgorithmKind") == {d["kind"] for d in algorithms()}
    assert _literal("_TermKind") == {d["kind"] for d in terms()}
