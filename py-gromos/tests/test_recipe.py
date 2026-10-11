"""
The run recipe — 3.9 steps 2 and 3.

Step 2: a `Simulation` exposes the effective recipe (TOML), the algorithm plan (JSON, one
fully resolved entry per algorithm) and the parser diagnostics.

Step 3: the recipe is a first-class object. `Recipe`, `Term`, `Algorithm` and `Plan` *are* the
Rust data types (serde <-> Python through pythonize, no Python-side copy to drift); a recipe
is immutable on the Python side (`update`/`with_*` return copies); a typo is a `RecipeError`,
never a silent default; everything round-trips through TOML/JSON/dict/pickle/bundle; and the
pre-recipe forms (`InputParameters`, the restraint kwargs, `from_sequence`) warn and are
translated into a recipe — one construction path, seen from several front-ends.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pytest
from test_gromosXX_references import REF_DIR, _parse_input_toml

from gromos import (
    Algorithm,
    AlgorithmSequence,
    Configuration,
    InputParameters,
    Plan,
    Recipe,
    Simulation,
    System,
    Term,
    Topology,
    algorithms,
    build_info,
    terms,
)
from gromos.exceptions import MissingFeatureError, PlanError, RecipeError, RunError

WATER = REF_DIR / "water_216_box"
NOBATH = REF_DIR / "water_216_nve_nobath"
COM_ROT = REF_DIR / "water_216_box_com_rot"
DISTRES = REF_DIR / "nacl_1water_distres"


def _paths(system_dir: Path) -> tuple[str, str, str]:
    inputs = _parse_input_toml(system_dir)
    topo, conf, params = (
        str((system_dir / inputs[key]).resolve())
        for key in ("topology", "configuration", "parameters")
    )
    return topo, conf, params


def _load(system_dir: Path):
    """The deprecated front-end — kept for the step-2 tests, which pin what it still guarantees."""
    topo_path, conf_path, params_path = _paths(system_dir)
    return Topology(topo_path), Configuration(conf_path), InputParameters(params_path)


def _system_and_recipe(system_dir: Path) -> tuple[System, Recipe]:
    topo_path, conf_path, params_path = _paths(system_dir)
    return System.from_files(topo_path, conf_path), Recipe.from_imd(params_path)


# ---------------------------------------------------------------------------------------------
# Step 2 — what a Simulation exposes
# ---------------------------------------------------------------------------------------------


def test_recipe_and_plan_are_exposed_and_ordered():
    topo, conf, params = _load(WATER)
    sim = Simulation(topo, conf, params)

    recipe = sim.recipe_toml
    assert recipe is not None
    assert "version = 1" in recipe
    for group in ("[control]", "[boundary]", "[forcefield", "[constraints]", "[outputs]"):
        assert group in recipe, f"missing {group} in\n{recipe}"

    plan = json.loads(sim.plan_json)
    kinds = [entry["kind"] for entry in plan]
    # GROMOS step order, one entry per algorithm, and it matches what actually runs.
    # gromosXX's head is the optional COM removal, then the lattice-shift tracker for any
    # non-vacuum boundary, then the forcefield (`create_md_sequence.cc:108-131`).
    head = [k for k in kinds[:3] if k in ("remove_com", "lattice_shift")]
    assert kinds[: len(head) + 1] == [*head, "forcefield"], kinds
    assert kinds[-1] == "energy_calculation"
    assert kinds.index("leap_frog_velocity") < kinds.index("leap_frog_position")
    assert len(kinds) == len(sim.algorithm_names)
    forcefield = next(e for e in plan if e["kind"] == "forcefield")
    # Fully resolved values, not "look it up in the parameter file".
    assert forcefield["cutoff_long"] == pytest.approx(params.cutoff)
    assert forcefield["four_pi_eps_i"] == pytest.approx(138.9354)
    assert forcefield["parallel"] is True  # 648 atoms > 100: the binary's policy


def test_diagnostics_report_absent_blocks():
    topo, conf, params = _load(NOBATH)
    sim = Simulation(topo, conf, params)
    notes = sim.diagnostics
    assert any("MULTIBATH block absent" in n and "no temperature coupling" in n for n in notes), (
        notes
    )
    # …and the plan indeed has no thermostat.
    kinds = [e["kind"] for e in json.loads(sim.plan_json)]
    assert "thermostat" not in kinds

    topo, conf, params = _load(WATER)
    sim = Simulation(topo, conf, params)
    # A reference input carries every block it needs; the notes only name the optional blocks
    # it deliberately omits (no pressure coupling, no restraints, no perturbation, no EM).
    assert not any("MULTIBATH" in n or "COMTRANSROT" in n for n in sim.diagnostics), sim.diagnostics
    assert all("block absent" in n for n in sim.diagnostics), sim.diagnostics


def test_from_sequence_shim_is_a_translation():
    """The deprecated names return and consume a `Plan`: there is no descriptor path left."""
    topo, conf, params = _load(WATER)
    plan = AlgorithmSequence.from_parameters(topo, params)
    assert isinstance(plan, Plan)
    sim = Simulation.from_sequence(topo, conf, params, plan)
    assert sim.recipe == Recipe.from_imd(_paths(WATER)[2])
    assert sim.plan.kinds == plan.kinds
    assert json.loads(sim.plan_json)[-1]["kind"] == "energy_calculation"
    with pytest.raises(TypeError):
        AlgorithmSequence()


def test_factory_parameters_default_to_what_gromosxx_does_without_the_blocks():
    """No COMTRANSROT in a factory-built run means no COM-motion removal (gromosXX:
    skip_step 0), and NVE means no bath — the recipe shows both."""
    system = System.from_files(
        str((WATER / "water_216_box.topo").resolve()), str((WATER / "water_216_box.conf").resolve())
    )
    sim = Simulation(system, Recipe.nve(dt=0.002, steps=2))
    assert "remove_com" not in sim.plan.kinds
    assert "thermostat" not in sim.plan.kinds
    # gromosXX's own way of saying "no coupling" is a bath with TAU < 0; the recipe keeps it.
    bath = sim.recipe.ensemble["thermostat"]
    assert bath is None or all(b["tau"] <= 0 for b in bath["baths"])
    assert repr(sim.recipe).startswith("Recipe(NVE,")


# ---------------------------------------------------------------------------------------------
# Step 3 — the recipe as data
# ---------------------------------------------------------------------------------------------


def test_recipe_round_trips_losslessly():
    _, recipe = _system_and_recipe(WATER)
    assert recipe.version == 1 == build_info()["recipe_version"]
    assert Recipe.from_toml(recipe.to_toml()) == recipe
    assert Recipe.from_json(recipe.to_json()) == recipe
    assert Recipe.from_dict(recipe.to_dict()) == recipe
    assert pickle.loads(pickle.dumps(recipe)) == recipe
    # The groups are plain dicts of the recipe's own fields; the .imd is what gromosXX runs.
    assert set(recipe.to_dict()) >= {
        "version",
        "title",
        "system",
        "control",
        "boundary",
        "forcefield",
        "constraints",
        "ensemble",
        "outputs",
        "execution",
        "inputs",
    }
    assert recipe.control["steps"] > 0 and recipe.control["dt"] > 0
    assert recipe.control == recipe.to_dict()["control"]
    imd = recipe.to_imd(n_atoms=648)
    for block in ("STEP", "SYSTEM", "BOUNDCOND", "PAIRLIST", "END"):
        assert block in imd
    assert Recipe.from_toml(recipe.to_toml()).to_imd(n_atoms=648) == imd


def test_update_returns_a_new_recipe_and_deep_merges():
    _, recipe = _system_and_recipe(WATER)
    updated = recipe.update(control={"steps": 7})
    assert updated.control["steps"] == 7
    assert recipe.control["steps"] != 7, "the original must be untouched"
    assert updated.control["dt"] == recipe.control["dt"], "siblings survive the merge"
    assert updated != recipe
    assert recipe.update() == recipe


def test_lambdas_kind_is_a_name_or_its_number():
    """A LAMBDAS row's `kind` takes gromosXX's NTLI name or its number 1-12, as in a `.imd`."""
    _, recipe = _system_and_recipe(WATER)
    row = {"g1": 1, "g2": 1, "a": 0.0, "b": 0.0, "c": 0.1, "d": 0.9, "e": 0.0}
    by_number = recipe.update(
        perturbation={"lambdas": {"enabled": True, "rows": [{"kind": 5, **row}]}}
    )
    by_name = recipe.update(
        perturbation={"lambdas": {"enabled": True, "rows": [{"kind": "vdw", **row}]}}
    )
    assert by_number == by_name
    assert by_number.perturbation["lambdas"]["rows"][0]["kind"] == "vdw"
    with pytest.raises(RecipeError):
        recipe.update(perturbation={"lambdas": {"enabled": True, "rows": [{"kind": 13, **row}]}})


def test_a_typo_is_a_recipe_error_not_a_silent_default():
    _, recipe = _system_and_recipe(WATER)
    with pytest.raises(RecipeError):
        recipe.update(control={"stepz": 7})
    with pytest.raises(RecipeError):
        recipe.update(contorl={"steps": 7})
    with pytest.raises(RecipeError):
        Recipe.from_toml(recipe.to_toml().replace("steps =", "stepz ="))
    with pytest.raises(RecipeError):
        Term("bogus")
    with pytest.raises(RecipeError):
        Algorithm("bogus")
    # The engine's errors subclass the builtins the binding raised before they existed.
    assert issubclass(RecipeError, ValueError) and issubclass(PlanError, ValueError)
    assert issubclass(MissingFeatureError, RuntimeError) and issubclass(RunError, RuntimeError)


def test_plan_is_the_step_as_editable_validated_data():
    system, recipe = _system_and_recipe(WATER)
    plan = recipe.plan(system)
    assert isinstance(plan, Plan)
    assert len(plan) == len(plan.kinds) > 0
    assert [a.kind for a in plan] == plan.kinds  # the sequence protocol (IndexError ends it)
    assert plan.kinds[-1] == "energy_calculation"
    assert "forcefield" in plan and plan["forcefield"].kind == "forcefield"
    assert plan[-1] == plan["energy_calculation"]
    with pytest.raises(IndexError):
        plan[len(plan)]
    with pytest.raises(PlanError):
        plan["barostat"]  # not in an NVE plan
    assert Plan.from_json(plan.to_json()).kinds == plan.kinds
    assert [d["kind"] for d in plan.to_dicts()] == plan.kinds
    assert Algorithm.from_dict(plan["forcefield"].to_dict()) == plan["forcefield"]
    plan.validate()
    # Breaking a step-order invariant is a PlanError, raised before anything is instantiated.
    broken = Plan(plan.to_list())
    broken.insert(0, broken.remove("energy_calculation"))
    with pytest.raises(PlanError):
        broken.validate()
    with pytest.raises(PlanError):
        Simulation(system, recipe, plan=broken)


def test_an_edited_plan_is_what_runs():
    system, recipe = _system_and_recipe(COM_ROT)
    reference = Simulation(system, recipe)
    plan = recipe.plan(system)
    assert "remove_com" in plan and reference.plan.kinds == plan.kinds
    plan.remove("remove_com")
    sim = Simulation(system, recipe, plan=plan)
    assert "remove_com" not in sim.plan.kinds
    assert len(sim.plan) == len(sim.algorithm_names) == len(plan) == len(reference.plan) - 1
    assert sim.recipe == reference.recipe, "editing the plan does not rewrite the recipe"


def test_simulation_exposes_its_recipe_and_plan():
    system, recipe = _system_and_recipe(WATER)
    sim = Simulation(system, recipe)
    assert sim.recipe == recipe
    assert sim.recipe.terms == []
    assert sim.plan.kinds == [e["kind"] for e in json.loads(sim.plan_json)]
    assert sim.recipe.to_toml() == sim.recipe_toml
    assert sim.diagnostics == recipe.diagnostics


def test_with_inputs_and_with_execution():
    _, recipe = _system_and_recipe(WATER)
    with_restraints = recipe.with_inputs(distrest="restraints.dat")
    assert with_restraints.inputs["distrest"] == "restraints.dat"
    assert recipe.inputs["distrest"] is None
    serial = recipe.with_execution("serial")
    assert serial.execution["parallel"] == "serial"
    assert recipe.execution["parallel"] == "auto"
    with pytest.raises(RecipeError):
        recipe.with_execution("sometimes")


def test_terms_registry_and_the_schnet_term():
    registry = {d["kind"]: d for d in terms()}
    assert "schnet" in registry
    info = build_info()
    assert info["recipe_version"] == 1 and isinstance(info["features"], list)
    example = registry["schnet"]
    term = Term("schnet", **example["params"])
    assert term.kind == "schnet" and term.feature == "ml"
    assert term.available == example["available"] == ("ml" in info["features"])
    assert Term.from_dict(term.to_dict()) == term

    system, recipe = _system_and_recipe(WATER)
    with_term = recipe.with_term(term)
    assert with_term.terms == [term] and recipe.terms == []
    assert with_term.without_terms() == recipe
    assert Recipe.from_toml(with_term.to_toml()) == with_term
    assert "orchestrator" in with_term.plan(system).kinds
    if not term.available:
        with pytest.raises(MissingFeatureError):
            Simulation(system, with_term)


def test_bundle_round_trip(tmp_path):
    topo_path, conf_path, _ = _paths(WATER)
    system, recipe = _system_and_recipe(WATER)
    input_toml = recipe.to_bundle(str(tmp_path), system, topo_path, conf_path)
    assert Path(input_toml).name == "input.toml"
    assert (tmp_path / "run.recipe.toml").exists() and (tmp_path / "run.imd").exists()
    assert Recipe.from_bundle(input_toml) == recipe

    from_bundle = Simulation.from_bundle(input_toml)
    direct = Simulation(system, recipe)
    assert from_bundle.n_atoms == system.n_atoms
    assert from_bundle.recipe == direct.recipe
    from_bundle.step(2)
    direct.step(2)
    assert np.array_equal(from_bundle.positions, direct.positions)


def test_deprecated_forms_warn_and_are_translations():
    topo_path, conf_path, params_path = _paths(WATER)
    with pytest.warns(DeprecationWarning, match="Recipe.from_imd"):
        params = InputParameters(params_path)
    with pytest.warns(DeprecationWarning, match="Recipe.nvt"):
        InputParameters.nvt(0.002, 10, 300.0)
    system = System.from_files(topo_path, conf_path)
    with pytest.warns(DeprecationWarning, match="Recipe"):
        legacy = Simulation(system, params)
    modern = Simulation(system, Recipe.from_imd(params_path))
    assert legacy.recipe == modern.recipe, "the shim is a translation, not a second path"
    legacy.step(3)
    modern.step(3)
    assert np.array_equal(legacy.positions, modern.positions)

    topo, conf = Topology(topo_path), Configuration(conf_path)
    with pytest.warns(DeprecationWarning, match="recipe.plan"):
        seq = AlgorithmSequence.from_parameters(topo, params)
    with pytest.warns(DeprecationWarning, match="plan="):
        Simulation.from_sequence(topo, conf, params, seq)

    inputs = _parse_input_toml(DISTRES)
    distrest = str((DISTRES / inputs["distrest"]).resolve())
    topo_d, conf_d, params_d = _paths(DISTRES)
    with pytest.warns(DeprecationWarning, match="with_inputs"):
        shimmed = Simulation(
            System.from_files(topo_d, conf_d), Recipe.from_imd(params_d), distrest=distrest
        )
    assert shimmed.recipe.inputs["distrest"] == distrest
    assert shimmed.recipe == Recipe.from_imd(params_d).with_inputs(distrest=distrest)


def test_registries_are_the_engine_s_own():
    kinds = [d["kind"] for d in algorithms()]
    assert kinds[-1] == "energy_calculation" and "forcefield" in kinds
    rules = {d["kind"]: d["rules"] for d in algorithms()}
    assert rules["energy_calculation"]["last"] and rules["forcefield"]["required"]
    for d in algorithms():
        assert Algorithm(d["kind"], **d["params"]).kind == d["kind"]
