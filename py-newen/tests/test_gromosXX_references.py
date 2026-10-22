"""
Integration tests that exercise the newen engine through the Python
Simulation API (newen.Simulation) and compare output against gromosXX
reference data.

This validates the full Rust simulation pipeline is accessible from Python
without shelling out to the md binary.

Usage style mirrors OpenMM but with gromosXX naming conventions:

    sim = Simulation("system.topo", "initial.cnf", "run.imd")
    sim.step(100)
    print(sim.energies, sim.positions, sim.forces)
"""

import re
import tomllib
from pathlib import Path

import numpy as np
import pytest

import newen.timeseries
from newen import (
    AlgorithmSequence,
    Configuration,
    EnergyTimeseries,
    EnergyTrajectory,
    InputParameters,
    Recipe,
    Simulation,
    System,
    Topology,
)

# Two slow-growth references take λ past 1 on purpose (`aladip_vacuum_fep_slowgrowth_cont` pins
# gromosXX's unbounded λ; `aladip_vacuum_em_slowgrowth` runs all its 500 steps here, where `md`
# stops the minimiser at 13): `Simulation` warns, and these tests compare the numbers instead.
pytestmark = [
    pytest.mark.filterwarnings(r"ignore:slow growth takes lambda out of \[0, 1\]"),
    # ch4_water_fep_precalc_twinrange (NSNB 5, NTWG 2) pins the PRECALCLAM long-range quirk
    pytest.mark.filterwarnings(r"ignore:PRECALCLAM with a twin-range pairlist"),
]

# Paths
REPO_ROOT = Path(__file__).parent.parent.parent
REF_DIR = REPO_ROOT / "crates" / "newen-md" / "tests" / "gromosXX_references"

# Tolerances (same as Rust tests)
ENERGY_REL_TOL = 1e-8
FORCE_ABS_TOL = 1e-6
POSITION_ABS_TOL = 1e-9


# ============================================================================
# Helpers to parse reference (expected) data files
# ============================================================================


def _parse_input_toml(system_dir):
    """Parse input.toml to get file paths for this reference system."""
    toml_path = system_dir / "input.toml"
    text = toml_path.read_text()

    def get_val(key):
        m = re.search(rf'^{key}\s*=\s*"(.+)"', text, re.MULTILINE)
        return m.group(1) if m else None

    return {
        "topology": get_val("topology"),
        "configuration": get_val("configuration"),
        "parameters": get_val("parameters"),
        "distrest": get_val("distrest"),
        "posresspec": get_val("posresspec"),
        "refpos": get_val("refpos"),
        "pttopo": get_val("pttopo"),
    }


def _parse_expected_energies(tre_path):
    """Parse gromosXX ENERGY03 block format reference energies.

    Returns list of dicts: {total, kinetic, potential}
    """
    text = tre_path.read_text()
    frames = []
    blocks = text.split("TIMESTEP")
    for block in blocks[1:]:
        lines = block.strip().split("\n")
        energy_start = None
        for i, line in enumerate(lines):
            if "ENERGY03" in line:
                energy_start = i + 1
                break
        if energy_start is None:
            continue

        vals = []
        for i in range(energy_start, len(lines)):
            line = lines[i].strip()
            if line.startswith("#") and "totals" in line:
                continue
            if line.startswith("#") or line == "END":
                break
            try:
                vals.append(float(line))
            except ValueError:
                break

        if len(vals) >= 3:
            frames.append(
                {
                    "total": vals[0],
                    "kinetic": vals[1],
                    "potential": vals[2],
                }
            )
    return frames


def _parse_expected_positions(trc_path):
    """Parse gromosXX POSITIONRED trajectory (3 columns: x, y, z).

    Returns list of Nx3 numpy arrays, one per frame.
    """
    text = trc_path.read_text()
    frames = []
    blocks = text.split("POSITIONRED")
    for block in blocks[1:]:
        lines = block.strip().split("\n")
        positions = []
        for line in lines:
            line = line.strip()
            if line == "END" or line.startswith("TIMESTEP") or line.startswith("GENBOX"):
                break
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 3:
                try:
                    positions.append([float(x) for x in parts])
                except ValueError:
                    break
        if positions:
            frames.append(np.array(positions, dtype=np.float64))
    return frames


def _parse_expected_forces(trf_path):
    """Parse gromosXX FREEFORCERED trajectory (3 columns per atom).

    Returns list of Nx3 numpy arrays, one per frame.
    """
    text = trf_path.read_text()
    frames = []
    blocks = text.split("FREEFORCERED")
    for block in blocks[1:]:
        lines = block.strip().split("\n")
        forces = []
        for line in lines:
            line = line.strip()
            if line == "END" or line.startswith("CONSFORCERED"):
                break
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 3:
                try:
                    forces.append([float(x) for x in parts])
                except ValueError:
                    break
        if forces:
            frames.append(np.array(forces, dtype=np.float64))
    return frames


def _recipe_for(system_dir):
    """The reference run as a `Recipe`: its `.imd`, with the auxiliary files (`pttopo`,
    `distrest`, `posresspec`, `refpos`) as `inputs` — what `md @pttopo …` receives."""
    inputs = _parse_input_toml(system_dir)
    aux = {
        key: str((system_dir / inputs[key]).resolve())
        for key in ("pttopo", "distrest", "posresspec", "refpos")
        if inputs.get(key)
    }
    recipe = Recipe.from_imd(str((system_dir / inputs["parameters"]).resolve())).with_inputs(**aux)
    return _apply_recipe_overrides(recipe, system_dir)


def _apply_recipe_overrides(recipe, system_dir):
    """The `[recipe_overrides]` of input.toml (dotted keys, e.g.
    `"execution.parallel"`): recipe fields no `.imd` block carries, what the Rust
    runner passes as `md @recipe_set key=value`. Only that table is parsed: not every
    input.toml is valid TOML as a whole."""
    section = re.search(
        r"^\[recipe_overrides\]\n(.*?)(?=^\[|\Z)",
        (system_dir / "input.toml").read_text(),
        re.MULTILINE | re.DOTALL,
    )
    overrides = tomllib.loads(section.group(1)) if section else {}
    for key, value in overrides.items():
        group, *path = key.split(".")
        patch = value
        for part in reversed(path):
            patch = {part: patch}
        recipe = recipe.update(**{group: patch})
    return recipe


def _create_simulation(system_dir):
    """Create a Simulation from a reference system directory: `Simulation(system, recipe)`,
    the front-end 3.9 step 3 makes canonical (the same `newen-run` path the `md`
    binary takes)."""
    inputs = _parse_input_toml(system_dir)
    system = System.from_files(
        str((system_dir / inputs["topology"]).resolve()),
        str((system_dir / inputs["configuration"]).resolve()),
    )
    return Simulation(system, _recipe_for(system_dir))


def _get_n_steps(system_dir):
    """NSTLIM, read from the recipe."""
    return _recipe_for(system_dir).control["steps"]


def _get_write_frequencies(system_dir):
    """(ntwe, ntwx, ntwf) — energy, trajectory and force write intervals, from the recipe
    (0 means "never" in GROMOS; treated as every step here so a reference that writes nothing
    still compares every frame it does have)."""
    outputs = _recipe_for(system_dir).outputs
    return (
        outputs["energy_every"] or 1,
        outputs["trajectory_every"] or 1,
        outputs["forces_every"] or 1,
    )


# ============================================================================
# Reference test systems — same as test_gromosXX_references.rs
# ============================================================================

AEDS_SYSTEMS = [
    "aladip_vacuum_aeds_below",
    "aladip_vacuum_aeds_inside",
    "aladip_vacuum_aeds_above",
    # the searches against gromosXX with D16, D25 and D32 patched (patches/)
    "aladip_vacuum_aeds_f2_fixed",
    "aladip_vacuum_aeds_f4_fixed",
    "aladip_vacuum_aeds_f5_fixed",
    "aladip_vacuum_aeds_f6_fixed",
    "aladip_vacuum_aeds_f6_ovf_fixed",
    "aladip_vacuum_aeds_sigma0_fixed",
    "aladip_vacuum_aeds_cont1_fixed",
    "aladip_vacuum_aeds_cont2_fixed",
]

REFERENCE_SYSTEMS = [
    # Level 0: Pair interactions
    "pair_lj",
    "pair_lj_mixed",
    "nacl_pair",
    # Level 1: Single molecule
    "water_single",
    "water_single_genvel",
    "benzene_vacuum",
    "nacl_pair_box",
    "nacl_pair_box_ewald",
    "nacl_pair_box_ewald_rect",
    "nacl_pair_box_p3m_fixed",
    "butane_vacuum",
    "butane_vacuum_ntbdn1",
    "aladip_vacuum",
    # Level 2: Small solvated
    "water_3_box",
    "nacl_1water_box",
    "nacl_3water_box",
    "water_3_box_twinrange",
    "water_10_box",
    "nacl_3water_cutoff",
    "nacl_water_box",
    "nacl_water_box_shifted",
    "nacl_water_box_ewald_na2clc1",
    # Constraint algorithms
    "nacl_1water_settle",
    "nacl_1water_lincs",
    "aladip_vacuum_lincs",
    # Grid-cell pairlist
    "water_1000_spc_gridcell",
    # Level 3: Bulk
    "water_216_box",
    "water_216_nve_nobath",
    "water_216_box_com",
    "water_216_box_com_rot",
    "water_216_box_ewald",
    "water_216_box_ewald_epsls",
    "water_216_box_p3m_fixed",
    "water_216_box_p3m_gauss_fixed",
    "water_216_nvt",
    "water_216_nvt_nosehoover",
    "water_216_nvt_nhc_chain",
    # Nosé–Hoover continued from a first half (NHCVARIABLES, NTINHT 0) against the D41 patch
    "water_216_nvt_nosehoover_cont_fixed",
    "water_216_nvt_nhc_cont_fixed",
    "water_216_npt",
    "water_216_npt_p3m_fixed",
    # Level 4: Full system
    "aladip_solvated",
    "aladip_trunc_oct",
    # Restraints
    "nacl_1water_distres",
    # Steepest-descent energy minimization
    "aladip_vacuum_em",
    "aladip_vacuum_em_ntw5",
    "aladip_vacuum_em_shake",
    "aladip_solvated_em_noshake",
    "aladip_solvated_em_shake",
    "aladip_solvated_em_posres",
    "aladip_solvated_em",
    "aladip_solvated_em_posres_groups",
    # Free-energy perturbation (Topology.apply_perturbation + NTG != 0), 3.9 step 1
    "ch4_water_fep",
    "ch4_water_fep_l000",
    "ch4_water_fep_l025",
    "ch4_water_fep_l075",
    "ch4_water_fep_l100",
    "meoh_water_fep",
    "aladip_vacuum_fep",
    "aladip_vacuum_fep_groups",
    "aladip_vacuum_fep_trg",
    # Individual λ (LAMBDAS), NLAM, perturbed SHAKE and perturbed distance restraints (1.8)
    "ch4_water_fep_lambdas",
    "ch4_water_fep_lambdas_off",
    "ch4_water_fep_nlam2",
    "aladip_vacuum_fep_lambdas",
    "aladip_vacuum_fep_nlam2",
    "aladip_vacuum_fep_shake",
    "aladip_vacuum_fep_shake_lambdas",
    "aladip_vacuum_fep_distres",
    "aladip_vacuum_fep_distres_lambdas",
    # NRDGL = 1 (λ from PERTDATA), PERTBONDSTRETCHH order, time-averaged distance restraints,
    # NLAM = 2 with LAMBDAS (1.8)
    "ch4_water_fep_nrdgl",
    "aladip_vacuum_fep_shake_hblock",
    "aladip_vacuum_fep_distres_tavg",
    "aladip_vacuum_fep_lambdas_nlam2",
    # Energy groups across the solute, DOFSET ranges not in bath order (1.8)
    "aladip_harmonic_covalent_groups",
    "aladip_multibath_swapped",
    # Barostat with SHAKE (constraint virial) and ROTTRANS degrees of freedom (1.8)
    "aladip_solvated_npt_shake",
    "aladip_multibath_rottrans",
    # Two Nose-Hoover baths, plain and chain, with different TEMP0 and TAU
    "aladip_multibath_nosehoover",
    "aladip_multibath_nhc",
    # NDFMIN over two baths with integer degrees of freedom, and half-integer ones (D40 patched)
    "aladip_vacuum_multibath_ndfmin3",
    "aladip_vacuum_multibath_ndfmin3_half_fixed",
    # Extended TI (PRECALCLAM) and the reaction-field constant from RCRF / EPSRF = 0 (1.8)
    "ch4_water_fep_precalc",
    "meoh_water_fep_precalc",
    "aladip_vacuum_fep_precalc",
    "aladip_vacuum_fep_lambdas_precalc",
    "aladip_vacuum_fep_shake_precalc",
    "aladip_vacuum_fep_distres_precalc",
    "aladip_vacuum_fep_pairpert_precalc_fixed",
    "ch4_water_fep_precalc_twinrange",
    # A charge change (Na+ -> neutral dummy, dQ = -1) under reaction field: the perturbed RF self
    # term, the RF excluded Na-Cl pair (as a 1-4 pair in `_14`); plain TI, PRECALCLAM, LAMBDAS
    "nacl_water_fep",
    "nacl_water_fep_precalc",
    "nacl_water_fep_lambdas",
    "nacl_water_fep_14",
    # Slow growth (DLAMT != 0): lambda += DLAMT*dt after every step, masses and perturbed terms
    # following it; with perturbed SHAKE, continued with NRDGL = 1 past lambda = 1, in solvent,
    # and under steepest descent
    "aladip_vacuum_fep_slowgrowth",
    "aladip_vacuum_fep_slowgrowth_shake",
    "aladip_vacuum_fep_slowgrowth_cont",
    "aladip_vacuum_fep_slowgrowth_t10",
    "ch4_water_fep_slowgrowth",
    "aladip_vacuum_em_slowgrowth",
    "aladip_solvated_epsrf0",
    "aladip_solvated_rcrf",
    # EDS through the recipe (1.8): FORM 1 against the released gromosXX, FORMS 2/3 against
    # gromosXX with its force weights corrected (patches/D2_eds_weights.patch)
    "aladip_vacuum_eds",
    "aladip_vacuum_eds3",
    "aladip_eds_solv",
    "two_waters_eds",
    "aladip_vacuum_eds_multis_fixed",
    "aladip_vacuum_eds_pairs_fixed",
    # TI-EDS: EDS with perturbed atoms (FORM 2 against the D2-patched gromosXX)
    "aladip_eds_solv_ti",
    "aladip_eds_solv_ti_groups",
    "aladip_eds_solv_ti_multis_fixed",
    # Accelerated EDS (1.8): FORM 1 at three window placements, the searches (FORMS 2, 4, 5,
    # 6, and 3 with a zero-variance target held by RESTREMIN 1), a search over two jobs, FORM 6
    # with overflowing prevalences; the released gromosXX (its D16, D25, D32 reproduced by the
    # recipe overrides) and the patched one
    *AEDS_SYSTEMS,
]

# EDS references compared end-state by end-state at frame 0 (`Simulation.eds_vi`,
# `Simulation.eds_vr`): the whole-block ones above, and the released gromosXX's FORMS 2/3, whose
# defective force weights leave only frame 0's end-state energies (and FORM 3's V_R) comparable.
EDS_SYSTEMS = {
    "aladip_vacuum_eds": True,
    "aladip_vacuum_eds3": True,
    "aladip_eds_solv": True,
    "two_waters_eds": True,
    "aladip_vacuum_eds_multis_fixed": True,
    "aladip_vacuum_eds_pairs_fixed": True,
    "aladip_vacuum_eds_multis": False,
    "aladip_vacuum_eds_pairs": True,
    "aladip_eds_solv_ti": True,
    "aladip_eds_solv_ti_groups": True,
    "aladip_eds_solv_ti_multis_fixed": True,
}

# Systems with a known, pre-existing position-output mismatch that is not a
# pyo3-binding gap: `aladip_trunc_oct`'s truncated-octahedron box means
# `newen`'s own `md` binary already disagrees with the gromosXX
# `.trc` convention here (confirmed by running `md` directly and comparing —
# energies and forces match gromosXX exactly, `.trc` positions don't). The
# Rust reference suite (`test_gromosXX_references.rs`) never compares
# positions at all, so this was never caught there either. Tracked as
# follow-up, not a Python-API defect — excluded only from the position test.
#
# `aladip_vacuum_em`: gromosXX repeats the converged frame at the end of its
# trajectories (the `md` binary does the same since 0.0.34 and the Rust suite
# passes); `Simulation` has no "converged, stop" step, so its position frame 6
# is the state after one more (no-op) minimisation step, not the repeat.
# Energies and forces validate; only positions are excluded here. `aladip_vacuum_em_slowgrowth`:
# the same minimisation with slow growth, the same repeated frame.
POSITION_MISMATCH_SYSTEMS = {"aladip_trunc_oct", "aladip_vacuum_em", "aladip_vacuum_em_slowgrowth"}

# Systems whose newen result is known to be WRONG against a correct reference — an engine
# defect, not a binding gap. `strict=True` so the fix is noticed and the entry removed.
# `water_216_nve_nobath`: the IMD has no MULTIBATH block; gromosXX then runs no temperature
# coupling, but `read_imd_file` keeps `TempBathParameters::default()` (Berendsen, 300 K,
# tau 0.1) and every newen path silently thermostats (3.9 A18, fixed by step 2).
EXPECTED_ENGINE_FAILURES = {
    # `water_216_nve_nobath` (absent MULTIBATH) passed here until 3.9 step 2 fixed the
    # parser defaults; it now runs as a regular reference system.
    # Same known FEP mismatch the Rust suite `ignore`s (Reference Test Status).
}

REFERENCE_PARAMS = [
    pytest.param(
        s,
        id=s,
        marks=(
            [pytest.mark.xfail(strict=True, reason=EXPECTED_ENGINE_FAILURES[s])]
            if s in EXPECTED_ENGINE_FAILURES
            else []
        ),
    )
    for s in REFERENCE_SYSTEMS
]


# ============================================================================
# Energy tests — through Python Simulation API
# ============================================================================


@pytest.mark.parametrize("system_name", REFERENCE_PARAMS)
def test_reference_energies(system_name):
    """Compare energies from newen.Simulation against gromosXX reference."""
    system_dir = REF_DIR / system_name
    if not system_dir.exists():
        pytest.skip(f"Reference system {system_name} not found")

    expected_tre = system_dir / "expected" / "energies.tre"
    if not expected_tre.exists():
        pytest.skip(f"No expected energies for {system_name}")

    expected = _parse_expected_energies(expected_tre)
    n_expected = len(expected)
    assert n_expected > 0, f"No expected energy frames for {system_name}"

    n_steps = _get_n_steps(system_dir)
    ntwe, _, _ = _get_write_frequencies(system_dir)

    # Create simulation (constructor initializes + runs step 0)
    sim = _create_simulation(system_dir)

    # Collect energies at write intervals (matching reference output)
    actual = [
        {
            "total": sim.total_energy,
            "kinetic": sim.kinetic_energy,
            "potential": sim.potential_energy,
        }
    ]

    for step in range(1, n_steps + 1):
        sim.step(1)
        if step % ntwe == 0:
            actual.append(
                {
                    "total": sim.total_energy,
                    "kinetic": sim.kinetic_energy,
                    "potential": sim.potential_energy,
                }
            )

    # Trim to match reference frame count
    if len(actual) > n_expected:
        actual = actual[:n_expected]

    assert len(actual) == n_expected, (
        f"{system_name}: frame count mismatch: {len(actual)} vs {n_expected}"
    )

    for i, (act, exp) in enumerate(zip(actual, expected, strict=True)):
        for key in ["total", "kinetic", "potential"]:
            a, e = act[key], exp[key]
            if abs(e) > 1e-10:
                rel = abs(a - e) / abs(e)
                assert rel < ENERGY_REL_TOL, (
                    f"{system_name} frame {i} {key}: {a} vs {e} (rel={rel:.2e})"
                )
            else:
                assert abs(a - e) < 1e-10, f"{system_name} frame {i} {key}: {a} vs {e}"


# ============================================================================
# EDS end states — through Python Simulation API
# ============================================================================


def _parse_expected_eds_frame0(tre_path):
    """Frame 0 of an EDS `.tre`: the V_i of the `# eds` rows and V_R (totals slot 35)."""
    lines = tre_path.read_text().splitlines()
    totals = lines.index("# totals")
    vr = float(lines[totals + 35])
    n_at = lines.index("# numstates")
    n = int(lines[n_at + 1])
    vi = [float(lines[n_at + 3 + k].split()[0]) for k in range(n)]
    return vi, vr


@pytest.mark.parametrize("system_name", sorted(EDS_SYSTEMS))
def test_reference_eds_end_states(system_name):
    """`Simulation.eds_vi` / `eds_vr` at frame 0 are gromosXX's `# eds` rows and slot 35."""
    system_dir = REF_DIR / system_name
    if not system_dir.exists():
        pytest.skip(f"Reference system {system_name} not found")
    vi_exp, vr_exp = _parse_expected_eds_frame0(system_dir / "expected" / "energies.tre")
    sim = _create_simulation(system_dir)
    assert len(sim.eds_vi) == len(vi_exp)
    for a, e in zip(sim.eds_vi, vi_exp, strict=True):
        assert abs(a - e) <= max(abs(e) * ENERGY_REL_TOL, 1e-10), f"{system_name}: {a} vs {e}"
    if EDS_SYSTEMS[system_name]:
        assert abs(sim.eds_vr - vr_exp) <= abs(vr_exp) * ENERGY_REL_TOL
    assert _recipe_for(system_dir).eds["enabled"]


def _parse_expected_aeds(tre_path):
    """Every frame of an A-EDS `.tre`: slots 34 (V_mix), 36 (E_max), 37 (E_min) and the
    `# eds` offset column."""
    frames = []
    lines = tre_path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line == "# totals":
            vmix, emax, emin = (float(lines[i + k]) for k in (34, 36, 37))
            n_at = lines.index("# numstates", i)
            n = int(lines[n_at + 1])
            eir = [float(lines[n_at + 3 + k].split()[3]) for k in range(n)]
            frames.append((vmix, emax, emin, eir))
    return frames


@pytest.mark.parametrize("system_name", AEDS_SYSTEMS)
def test_reference_aeds_search(system_name):
    """`Simulation.eds_vmix` / `eds_emax` / `eds_emin` / `eds_eir` every frame are gromosXX's
    slots 34, 36, 37 and offset column (the released binary's written after the search update,
    the patched one's before it, as the recipe overrides say); the continuation jobs read their
    start state from the configuration's AEDSSEARCH block."""
    system_dir = REF_DIR / system_name
    expected = _parse_expected_aeds(system_dir / "expected" / "energies.tre")
    sim = _create_simulation(system_dir)
    assert _recipe_for(system_dir).aeds["enabled"]

    def close(a, e):
        return abs(a - e) <= max(abs(e) * ENERGY_REL_TOL, 1e-9)

    for k, (vmix, emax, emin, eir) in enumerate(expected):
        if k:
            sim.step(1)
        assert close(sim.eds_vmix, vmix), f"{system_name} frame {k}: {sim.eds_vmix} vs {vmix}"
        assert close(sim.eds_emax, emax), f"{system_name} frame {k}: {sim.eds_emax} vs {emax}"
        assert close(sim.eds_emin, emin), f"{system_name} frame {k}: {sim.eds_emin} vs {emin}"
        assert all(close(a, e) for a, e in zip(sim.eds_eir, eir, strict=True)), (
            f"{system_name} frame {k}: {sim.eds_eir} vs {eir}"
        )


def test_aeds_continuation_needs_the_aedssearch_block():
    """NTIAEDSS = 0 without an AEDSSEARCH block in the configuration is refused, as gromosXX
    does."""
    system_dir = REF_DIR / "aladip_vacuum_aeds_cont2_fixed"
    inputs = _parse_input_toml(system_dir)
    system = System.from_files(
        str((system_dir / inputs["topology"]).resolve()),
        str((REF_DIR / "aladip_vacuum" / "aladip_vacuum.conf").resolve()),
    )
    with pytest.raises(Exception, match="no AEDSSEARCH block"):
        Simulation(system, _recipe_for(system_dir))


# ============================================================================
# dH/dλ tests — through Python Simulation API
# ============================================================================


def _parse_expected_dhdl(trg_path, slot=1):
    """Value `slot` (1-based) of every FREEENERDERIVS03 block's totals; the first is the whole
    dH/dλ, the 35th the EDS reference state's ∂V_R/∂λ."""
    values = []
    lines = trg_path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.strip() == "# totals":
            values.append(float(lines[i + slot]))
    return values


@pytest.mark.parametrize("system_name", REFERENCE_PARAMS)
def test_reference_dhdl(system_name):
    """`Simulation.dhdl` is the whole derivative gromosXX writes first in each `.trg` frame
    (kinetic + potential + constraints + restraints), and `dhdl_terms` splits it."""
    system_dir = REF_DIR / system_name
    expected_trg = system_dir / "expected" / "free_energy.trg"
    if not expected_trg.exists():
        pytest.skip(f"No expected free-energy trajectory for {system_name}")

    expected = _parse_expected_dhdl(expected_trg)
    if not expected:
        pytest.skip(f"{system_name}: its .trg has no frames (NTWG = 0)")
    every = _recipe_for(system_dir).outputs["ntwg"] or 1

    sim = _create_simulation(system_dir)
    actual = [sim.dhdl]
    terms = [sim.dhdl_terms]
    step = 0
    while len(actual) < len(expected):
        sim.step(1)
        step += 1
        if step % every == 0:
            actual.append(sim.dhdl)
            terms.append(sim.dhdl_terms)

    for i, (a, e) in enumerate(zip(actual, expected, strict=True)):
        assert abs(a - e) <= max(abs(e) * 1e-6, 1e-8), (
            f"{system_name} frame {i} dhdl: {a} vs {e}"
        )
    for t in terms:
        assert set(t) == {"lj", "crf", "bonded", "kinetic", "constraints", "distanceres", "eds"}
    for i, (t, e) in enumerate(
        zip(terms, _parse_expected_dhdl(expected_trg, slot=35), strict=True)
    ):
        assert abs(t["eds"] - e) <= max(abs(e) * 1e-6, 1e-8), (
            f"{system_name} frame {i} dV_R/dλ: {t['eds']} vs {e}"
        )


# ============================================================================
# Position tests — through Python Simulation API
# ============================================================================


@pytest.mark.parametrize("system_name", REFERENCE_PARAMS)
def test_reference_positions(system_name):
    """Compare positions from newen.Simulation against gromosXX reference."""
    if system_name in POSITION_MISMATCH_SYSTEMS:
        pytest.skip(f"{system_name}: known pre-existing position-output mismatch")

    system_dir = REF_DIR / system_name
    if not system_dir.exists():
        pytest.skip(f"Reference system {system_name} not found")

    expected_trc = system_dir / "expected" / "trajectory.trc"
    if not expected_trc.exists():
        pytest.skip(f"No expected trajectory for {system_name}")

    expected = _parse_expected_positions(expected_trc)
    n_expected = len(expected)
    assert n_expected > 0, f"No expected position frames for {system_name}"

    n_steps = _get_n_steps(system_dir)
    _, ntwx, _ = _get_write_frequencies(system_dir)

    sim = _create_simulation(system_dir)

    # Get box dimensions for minimum image convention
    inputs = _parse_input_toml(system_dir)
    conf = Configuration(str((system_dir / inputs["configuration"]).resolve()))
    box_dims = np.array(conf.box_dimensions)

    # Collect positions at write intervals
    actual = [sim.positions.copy()]
    for step in range(1, n_steps + 1):
        sim.step(1)
        if step % ntwx == 0:
            actual.append(sim.positions.copy())

    if len(actual) > n_expected:
        actual = actual[:n_expected]

    assert len(actual) == n_expected, (
        f"{system_name}: position frame count mismatch: {len(actual)} vs {n_expected}"
    )

    for i, (act, exp) in enumerate(zip(actual, expected, strict=True)):
        diff = act - exp
        # Apply minimum image convention for periodic systems
        if box_dims.min() > 0:
            diff -= np.round(diff / box_dims) * box_dims
        max_diff = np.abs(diff).max()
        assert max_diff < POSITION_ABS_TOL, (
            f"{system_name} frame {i}: max position diff = {max_diff:.2e}"
        )


# ============================================================================
# Force tests — through Python Simulation API
# ============================================================================


@pytest.mark.parametrize("system_name", REFERENCE_PARAMS)
def test_reference_forces(system_name):
    """Compare forces from newen.Simulation against gromosXX reference."""
    system_dir = REF_DIR / system_name
    if not system_dir.exists():
        pytest.skip(f"Reference system {system_name} not found")

    expected_trf = system_dir / "expected" / "forces.trf"
    if not expected_trf.exists():
        pytest.skip(f"No expected forces for {system_name}")

    expected = _parse_expected_forces(expected_trf)
    n_expected = len(expected)
    assert n_expected > 0, f"No expected force frames for {system_name}"

    n_steps = _get_n_steps(system_dir)
    _, _, ntwf = _get_write_frequencies(system_dir)

    sim = _create_simulation(system_dir)

    # Collect forces at write intervals
    actual = [sim.forces.copy()]
    for step in range(1, n_steps + 1):
        sim.step(1)
        if step % ntwf == 0:
            actual.append(sim.forces.copy())

    if len(actual) > n_expected:
        actual = actual[:n_expected]

    assert len(actual) == n_expected, (
        f"{system_name}: force frame count mismatch: {len(actual)} vs {n_expected}"
    )

    for i, (act, exp) in enumerate(zip(actual, expected, strict=True)):
        diff = np.abs(act - exp)
        max_diff = diff.max()
        assert max_diff < FORCE_ABS_TOL, f"{system_name} frame {i}: max force diff = {max_diff:.2e}"


# ============================================================================
# 3.2 — System constructor produces identical results to topo/conf/params form
# ============================================================================


def test_system_constructor_matches_file_constructor():
    """Simulation(System, params) must be bit-identical to Simulation(topo, conf, params).

    Uses the same params object for both paths — isolates the constructor
    change from any parameter difference.
    """
    system_dir = REF_DIR / "water_216_nvt"
    if not system_dir.exists():
        pytest.skip("water_216_nvt reference not found")

    inputs = _parse_input_toml(system_dir)
    topo_path = str((system_dir / inputs["topology"]).resolve())
    conf_path = str((system_dir / inputs["configuration"]).resolve())
    params_path = str((system_dir / inputs["parameters"]).resolve())

    params = InputParameters(params_path)

    # Old three-arg path
    sim_old = Simulation(Topology(topo_path), Configuration(conf_path), params)

    # New two-arg path: Simulation(System, params)
    system = System.from_files(topo_path, conf_path)
    sim_new = Simulation(system, params)

    assert sim_new.total_energy == pytest.approx(sim_old.total_energy, rel=ENERGY_REL_TOL), (
        f"total energy: {sim_new.total_energy} vs {sim_old.total_energy}"
    )
    assert sim_new.kinetic_energy == pytest.approx(sim_old.kinetic_energy, rel=ENERGY_REL_TOL), (
        f"kinetic energy: {sim_new.kinetic_energy} vs {sim_old.kinetic_energy}"
    )
    assert sim_new.potential_energy == pytest.approx(
        sim_old.potential_energy, rel=ENERGY_REL_TOL
    ), f"potential energy: {sim_new.potential_energy} vs {sim_old.potential_energy}"


def test_factory_nvt_properties():
    """InputParameters.nvt(...) stores dt, nstlim and temperature correctly."""
    p = InputParameters.nvt(0.001, 500, 298.0)
    assert p.dt == pytest.approx(0.001)
    assert p.nstlim == 500
    assert p.temperature == pytest.approx(298.0)


def test_factory_constraints_knob():
    """3.5 M1: nve/nvt/npt accept constraints="none"/"hbonds"/"allbonds".

    Default stays "none" (NTC=1, GROMOS-faithful, backward compatible).
    """
    assert InputParameters.nve(0.001, 10).ntc == 1
    assert InputParameters.nve(0.001, 10).constraints == "none"

    assert InputParameters.nve(0.001, 10, constraints="hbonds").ntc == 2
    assert InputParameters.nve(0.001, 10, constraints="hbonds").constraints == "hbonds"

    assert InputParameters.nvt(0.001, 10, 300.0, constraints="allbonds").ntc == 3
    assert InputParameters.npt(0.001, 10, 300.0, 1.0, constraints="hbonds").ntc == 2

    with pytest.raises(ValueError):
        InputParameters.nve(0.001, 10, constraints="bogus")


def test_constrained_system_stable_with_factory_params():
    """Regression test for the bug that motivated M1/M2.

    `aladip_solvated` (solute has flexible H-bonds) diverges to NaN within
    ~50-100 steps when run with factory-built params that don't SHAKE
    (constraints="none", the old/default behavior — kept as a documented
    contrast, not a design goal). `constraints="hbonds"` must keep it stable,
    and the reported temperature (M2: constraint-aware DOF, not bare 3*n_atoms)
    must stay in a physically sane range, not explode.
    """

    def load_aladip_solvated():
        topo = Topology(str(REF_DIR / "shared" / "aladip.topo"))
        topo.solvate(20)
        conf = Configuration(str(REF_DIR / "shared" / "aladip.conf"))
        return System(topo, conf)

    params = InputParameters.nvt(dt=0.002, steps=100, temperature=300.0, constraints="hbonds")
    sim = Simulation(load_aladip_solvated(), params)

    temperatures = []
    for _ in range(100):
        sim.step(1)
        temperatures.append(sim.temperature)

    assert all(np.isfinite(t) for t in temperatures), "temperature diverged (NaN/inf)"
    assert all(50.0 < t < 600.0 for t in temperatures), (
        f"temperature outside sane range: min={min(temperatures)}, max={max(temperatures)}"
    )
    # Should track the 300 K bath target reasonably after equilibrating.
    mean_late = sum(temperatures[-30:]) / 30
    assert 250.0 < mean_late < 350.0, f"NVT not tracking 300 K target: {mean_late}"


def test_steepest_descent_via_simulation():
    """3.5: steepest-descent EM must actually run through Simulation.

    Previously `build_simulation` never checked `ntem`/`steepest_descent()` and
    silently ran plain leap-frog with dt=0 (a no-op) instead of minimizing.
    Uses `aladip_vacuum`'s own starting structure (not an `_em` reference,
    which would already be minimized and produce a flat, uninformative trace).
    """
    system_dir = REF_DIR / "aladip_vacuum"
    topo_path = str(system_dir / ".." / "shared" / "aladip.topo")
    conf_path = str(system_dir / "aladip_vacuum.conf")

    system = System.from_files(topo_path, conf_path)
    params = InputParameters.steepest_descent(steps=50)
    sim = Simulation(system, params)

    assert any("Steepest" in name for name in sim.algorithm_names), sim.algorithm_names
    assert "LeapFrogVelocity" not in sim.algorithm_names
    assert "LeapFrogPosition" not in sim.algorithm_names

    frames = sim.run(50, ene_freq=1)
    potential = frames[:, 2]

    assert np.all(np.isfinite(potential)), "EM potential energy diverged"
    assert potential[-1] < potential[0], (
        f"EM did not decrease potential energy: {potential[0]} -> {potential[-1]}"
    )
    # Once converged (dE < DELE), the algorithm is a documented no-op — the
    # trace should plateau, not keep oscillating indefinitely.
    assert np.allclose(potential[-5:], potential[-1]), "did not settle/converge"


def test_steepest_descent_via_algorithm_sequence():
    """The deprecated `AlgorithmSequence.minimize / .from_parameters` shims (now returning the
    `Plan` of the parameters) must agree with the direct Simulation(system, params) EM path."""
    system_dir = REF_DIR / "aladip_vacuum"
    topo_path = str(system_dir / ".." / "shared" / "aladip.topo")
    conf_path = str(system_dir / "aladip_vacuum.conf")

    system = System.from_files(topo_path, conf_path)
    params = InputParameters.steepest_descent(steps=30)

    seq_direct = AlgorithmSequence.minimize(system.topology, params)
    seq_auto = AlgorithmSequence.from_parameters(system.topology, params)
    assert (
        seq_direct.kinds
        == seq_auto.kinds
        == ["forcefield", "steepest_descent", "energy_calculation"]
    )

    sim = Simulation.from_sequence(system.topology, system.configuration, params, seq_auto)
    frames = sim.run(30, ene_freq=1)
    potential = frames[:, 2]
    assert np.all(np.isfinite(potential))
    assert potential[-1] < potential[0]

    # Names agreeing is not parity ('s trap, 3.9 step 0): the composable
    # path must produce the *same numbers* as the direct path, exactly.
    sim_direct = Simulation(system, params)
    frames_direct = sim_direct.run(30, ene_freq=1)
    assert np.array_equal(frames_direct, frames), (
        "AlgorithmSequence.from_parameters and Simulation(system, params) disagree on "
        "aladip_vacuum EM energies"
    )
    assert np.array_equal(sim_direct.positions, sim.positions)


def test_volume_and_pressure_getters():
    """`sim.volume`/`sim.pressure` must match the equivalent `run()` array columns
    (volume/pressure at indices 4/5), and NPT's volume must actually respond to
    the barostat while NVE/NVT's stays exactly fixed (no PressureCalculation/
    BerendsenBarostat in their sequence).
    """
    topo_path = str(REF_DIR / "water_216_box" / "water_216_box.topo")
    conf_path = str(REF_DIR / "water_216_box" / "water_216_box.conf")

    # sim.volume/sim.pressure agree with the run() array at the same frame.
    system = System.from_files(topo_path, conf_path)
    params = InputParameters.from_file(str(REF_DIR / "water_216_npt" / "water_216_npt.in"))
    sim = Simulation(system, params)
    frames = sim.run(5, ene_freq=1)
    assert sim.volume == pytest.approx(frames[-1, 4])
    assert sim.pressure == pytest.approx(frames[-1, 5])

    # NVE/NVT: fixed box. NPT: box responds to the barostat.
    def volume_trace(in_file, n_steps=100):
        s = System.from_files(topo_path, conf_path)
        p = InputParameters.from_file(str(in_file))
        sim = Simulation(s, p)
        volumes = []
        for _ in range(n_steps):
            sim.step(1)
            volumes.append(sim.volume)
        return volumes

    nve_vol = volume_trace(REF_DIR / "water_216_box" / "water_216_box.in")
    nvt_vol = volume_trace(REF_DIR / "water_216_nvt" / "water_216_nvt.in")
    npt_vol = volume_trace(REF_DIR / "water_216_npt" / "water_216_npt.in")

    assert max(nve_vol) - min(nve_vol) == 0.0
    assert max(nvt_vol) - min(nvt_vol) == 0.0
    assert max(npt_vol) - min(npt_vol) > 0.01


def test_system_factory_workflow():
    """Full 3.2 workflow: System.from_files + InputParameters.nvt + Simulation + step.

    Mirrors the intended API:
        system = System.from_files("system.topo", "initial.cnf")
        params = InputParameters.nvt(0.002, 1000, 300.0)
        sim = Simulation(system, params)
        sim.step(100)

    Uses water_3_box (9 atoms) so factory cutoff defaults are safe regardless
    of what the config was equilibrated at.
    """
    system_dir = REF_DIR / "water_3_box"
    if not system_dir.exists():
        pytest.skip("water_3_box reference not found")

    inputs = _parse_input_toml(system_dir)
    topo_path = str((system_dir / inputs["topology"]).resolve())
    conf_path = str((system_dir / inputs["configuration"]).resolve())

    system = System.from_files(topo_path, conf_path)
    assert system.n_atoms == 9

    params = InputParameters.nvt(0.002, 10, 300.0)
    sim = Simulation(system, params)

    e0 = sim.total_energy
    sim.step(10)
    e10 = sim.total_energy

    assert np.isfinite(e0), f"step-0 energy not finite: {e0}"
    assert np.isfinite(e10), f"step-10 energy not finite: {e10}"


def test_run_matches_step_loop():
    """sim.run(steps, ene_freq) must match an equivalent sim.step(1) loop.

    Columns mirror the .tre energy-block layout (newen_io::energy::EnergyFrame):
    [time, kinetic, potential, total, volume, pressure,
     bond, angle, improper, dihedral, lj, coulomb].
    Energy-only columns (indices 0,1,2,3,6,7,8,9,10,11) are checked exactly against
    `sim.energies`; volume/pressure (4,5) aren't exposed on `Energy` so are only
    checked for internal consistency between the two run paths.
    """
    system_dir = REF_DIR / "water_216_box"
    topo_path = str(system_dir / "water_216_box.topo")
    conf_path = str(system_dir / "water_216_box.conf")
    input_path = str(system_dir / "water_216_box.in")

    energy_cols = [0, 1, 2, 3, 6, 7, 8, 9, 10, 11]

    sim_run = Simulation(topo_path, conf_path, input_path)
    frames = sim_run.run(10, ene_freq=1)
    assert frames.shape == (11, 13)

    sim_step = Simulation(topo_path, conf_path, input_path)

    def row(sim):
        e = sim.energies
        return [
            sim.time - sim.dt,  # the energies are those of the step just taken
            e.kinetic,
            e.potential,
            e.total,
            e.bond,
            e.angle,
            e.improper,
            e.dihedral,
            e.lj,
            e.coulomb,
        ]

    manual = [row(sim_step)]
    for _ in range(10):
        sim_step.step(1)
        manual.append(row(sim_step))
    manual = np.array(manual)

    assert np.allclose(frames[:, energy_cols], manual)


def test_start_time_stamps_every_row():
    """STEP T = 10 ps: `Simulation.time` and the energy rows carry T + step*DT, the times of
    gromosXX's `.tre` (`md.cc:241`)."""
    system_dir = REF_DIR / "aladip_vacuum_fep_slowgrowth_t10"
    lines = (system_dir / "expected" / "energies.tre").read_text().splitlines()
    expected = [float(lines[i + 1].split()[1]) for i, x in enumerate(lines) if x == "TIMESTEP"]
    sim = _create_simulation(system_dir)
    assert sim.time == expected[1]
    frames = sim.run(4, ene_freq=1)
    assert np.allclose(frames[:, 0], expected[:5], rtol=0, atol=1e-12)


def test_energy_timeseries():
    """EnergyTimeseries wraps a run() array: column access, block averaging, to_dataframe()."""
    system_dir = REF_DIR / "water_216_box"
    topo_path = str(system_dir / "water_216_box.topo")
    conf_path = str(system_dir / "water_216_box.conf")
    input_path = str(system_dir / "water_216_box.in")

    sim = Simulation(topo_path, conf_path, input_path)
    frames = sim.run(10, ene_freq=2)
    ts = EnergyTimeseries(frames)

    assert len(ts) == frames.shape[0]
    assert np.array_equal(ts.total, frames[:, 3])
    assert np.array_equal(ts.time, frames[:, 0])
    assert ts.time[0] == 0.0

    mean, error = ts.block_average("total", block_size=2)
    assert np.isfinite(mean)
    assert error >= 0.0

    with pytest.raises(ValueError):
        ts.block_average("not_a_component", block_size=2)

    # Default dataframe backend is polars.
    df = ts.to_dataframe()
    assert df.columns == list(newen.timeseries.COLUMNS)
    assert df["total"].to_list() == list(frames[:, 3])

    df_dict = ts.to_dataframe(backend="dict")
    assert set(df_dict.keys()) >= {"time", "total", "kinetic", "potential"}

    with pytest.raises(ValueError):
        ts.to_dataframe(backend="not_a_backend")

    # Default plot backend is plotly; both backends must return a renderable object.
    fig = ts.plot("total", "bond")
    assert fig.data[0].name == "total"

    ax = ts.plot("total", backend="matplotlib")
    assert ax.get_xlabel() == "time (ps)"

    with pytest.raises(ValueError):
        ts.plot("not_a_component")

    with pytest.raises(ValueError):
        ts.plot("total", backend="not_a_backend")

    with pytest.raises(ValueError):
        EnergyTimeseries(np.zeros((3, 5)))


def test_run_frames_carry_the_time_gromosxx_gives_their_energies():
    """Frame k of run() is step k: its energies and its time match row k of gromosXX's .tre.

    The energies of a step are evaluated at the positions the step starts from, after which the
    step counter advances; the time column used to be read after that advance, one dt late."""
    system_dir = REF_DIR / "pair_lj"
    sim = Simulation(
        str(system_dir / "pair_lj.topo"),
        str(system_dir / "pair_lj.conf"),
        str(system_dir / "pair_lj.in"),
    )
    ts = EnergyTimeseries(sim.run(9, ene_freq=1))
    tre = EnergyTrajectory(str(system_dir / "expected" / "energies.tre"), properties=["totene"])

    assert ts.time[0] == 0.0
    np.testing.assert_allclose(ts.time, tre.times, rtol=0, atol=1e-12)
    np.testing.assert_allclose(ts.total, tre["totene"], rtol=ENERGY_REL_TOL)
