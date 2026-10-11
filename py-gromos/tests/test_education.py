"""gromos.education: the theory functions, the toy systems and the nine experiments' checks."""

import inspect
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from gromos import EnergyTimeseries, Simulation
from gromos.education import (
    Check,
    Result,
    argon_pair,
    diatomic,
    experiments,
    lj_gas,
    theory,
    with_bath,
)
from gromos.education.systems import ARGON_MASS, EPSILON, K_B, SIGMA

REF_DIR = (
    Path(__file__).parent.parent.parent / "crates" / "gromos-md" / "tests" / "gromosXX_references"
)

# -- theory -----------------------------------------------------------------------------------


def test_lj_energy_and_minimum_by_hand():
    c6, c12 = 4 * EPSILON * SIGMA**6, 4 * EPSILON * SIGMA**12
    assert theory.lj_energy(SIGMA, c6, c12) == pytest.approx(0.0, abs=1e-12)
    r_min, eps = theory.lj_minimum(c6, c12)
    assert r_min == pytest.approx(2 ** (1 / 6) * SIGMA, rel=1e-12)
    assert eps == pytest.approx(EPSILON, rel=1e-12)
    assert theory.lj_energy(r_min, c6, c12) == pytest.approx(-EPSILON, rel=1e-12)
    np.testing.assert_allclose(theory.lj_energy([1.0, 2.0], 1.0, 1.0), [0.0, 1 / 4096 - 1 / 64])


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_maxwell_boltzmann_pdf_normalised_with_the_closed_form_mean(dim):
    T, m = 120.0, ARGON_MASS
    v = np.linspace(0.0, 3.0, 300_001)
    pdf = theory.maxwell_boltzmann_speed_pdf(v, T, m, dim)
    assert np.trapezoid(pdf, v) == pytest.approx(1.0, rel=1e-6)
    mean = theory.maxwell_boltzmann_mean_speed(T, m, dim)
    assert np.trapezoid(v * pdf, v) == pytest.approx(mean, rel=1e-6)
    closed = {
        1: math.sqrt(2 * K_B * T / (math.pi * m)),
        2: math.sqrt(math.pi * K_B * T / (2 * m)),
        3: math.sqrt(8 * K_B * T / (math.pi * m)),
    }[dim]
    assert mean == pytest.approx(closed, rel=1e-12)


def test_maxwell_boltzmann_pdf_2d_by_hand():
    a = ARGON_MASS / (K_B * 100.0)
    expected = a * 0.2 * math.exp(-a * 0.04 / 2)
    assert theory.maxwell_boltzmann_speed_pdf(0.2, 100.0, ARGON_MASS, 2) == pytest.approx(expected)


def test_canonical_kinetic_energy_is_a_gamma_distribution():
    T, dof = 120.0, 126
    kt = K_B * T
    ke = np.linspace(1e-6, 4 * dof / 2 * kt, 200_001)
    pdf = theory.canonical_kinetic_energy_pdf(ke, T, dof)
    mean = np.trapezoid(ke * pdf, ke)
    assert np.trapezoid(pdf, ke) == pytest.approx(1.0, rel=1e-6)
    assert mean == pytest.approx(dof / 2 * kt, rel=1e-6)
    sd = math.sqrt(np.trapezoid((ke - mean) ** 2 * pdf, ke))
    assert sd / mean == pytest.approx(theory.canonical_kinetic_energy_relative_width(dof), 1e-5)
    # dof = 2 is an exponential: P(K) = exp(-K/kT) / kT
    assert theory.canonical_kinetic_energy_pdf(kt, T, 2) == pytest.approx(math.exp(-1) / kt)


def test_relative_width_step_limit_and_temperature_by_hand():
    assert theory.canonical_kinetic_energy_relative_width(2) == 1.0
    assert theory.canonical_kinetic_energy_relative_width(200) == pytest.approx(0.1)
    assert theory.harmonic_step_limit(9.8) == pytest.approx(3.1194, abs=1e-4)
    assert theory.temperature_from_kinetic(K_B * 300.0 * 1.5, 3) == pytest.approx(300.0)
    np.testing.assert_allclose(theory.temperature_from_kinetic([0.0, K_B], 2), [0.0, 1.0])


# -- systems ----------------------------------------------------------------------------------


def test_lj_gas_2d_stays_planar_and_reports_the_2d_temperature():
    n = 64
    system, recipe = lj_gas(n=n)
    assert system.n_atoms == n
    sim = Simulation(system, recipe)
    z0 = sim.positions[:, 2].copy()
    sim.step(2000)
    np.testing.assert_array_equal(sim.positions[:, 2], z0)
    np.testing.assert_array_equal(sim.velocities[:, 2], 0.0)
    expected = theory.temperature_from_kinetic(sim.kinetic_energy, 2 * n - 2)
    assert sim.temperature == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(("n", "dim", "density"), [(64, 2, 0.45), (63, 2, 0.45), (125, 3, 0.5)])
def test_lj_gas_equal_speeds_have_one_speed_and_no_momentum(n, dim, density):
    system, _ = lj_gas(n=n, dim=dim, density=density, speeds="equal")
    v = system.velocities
    speed = np.linalg.norm(v, axis=1)
    np.testing.assert_allclose(speed, math.sqrt(dim * EPSILON / ARGON_MASS), rtol=1e-7)
    np.testing.assert_allclose(v.sum(axis=0), 0.0, atol=1e-6)
    if dim == 2:
        np.testing.assert_array_equal(v[:, 2], 0.0)


def test_lj_gas_refuses_a_box_shorter_than_twice_the_cutoff():
    with pytest.raises(ValueError, match="cut-off"):
        lj_gas(n=27, dim=3, density=0.8)


def test_argon_pair_reads_its_lj_parameters_from_the_topology():
    system, _ = argon_pair(0.4)
    top = system.topology
    c6, c12, _, _ = top.lj_parameters(int(top.iac[0]), int(top.iac[1]))
    assert c6 == pytest.approx(4 * EPSILON * SIGMA**6, rel=1e-6)
    assert c12 == pytest.approx(4 * EPSILON * SIGMA**12, rel=1e-6)
    assert np.linalg.norm(system.positions[1] - system.positions[0]) == pytest.approx(0.4)


def test_diatomic_vibrates_with_the_requested_period():
    system, recipe = diatomic()
    sim = Simulation(system, recipe.update(control={"dt": 0.0001}))
    stretch = []
    for _ in range(400):
        stretch.append(np.linalg.norm(sim.positions[1] - sim.positions[0]) - 0.1)
        sim.step(1)
    stretch = np.asarray(stretch)
    up = np.flatnonzero((stretch[:-1] < 0) & (stretch[1:] >= 0))
    np.testing.assert_allclose(np.diff(up) * 0.1, 9.8, atol=0.1)  # fs


@pytest.mark.parametrize(
    ("algorithm", "engine"),
    [("berendsen", "berendsen"), ("nose_hoover", {"nose_hoover_chain": 3})],
)
def test_with_bath_couples_one_bath(algorithm, engine):
    _, recipe = lj_gas(n=16, density=0.2)
    before = recipe.ensemble
    thermostat = with_bath(recipe, 150.0, tau=0.3, algorithm=algorithm).ensemble["thermostat"]
    assert thermostat["algorithm"] == engine
    assert thermostat["baths"] == [{"tau": 0.3, "temperature": 150.0}]
    assert recipe.ensemble == before


def test_with_bath_rejects_an_unknown_algorithm():
    _, recipe = lj_gas(n=16, density=0.2)
    with pytest.raises(ValueError, match="algorithm"):
        with_bath(recipe, 120.0, algorithm="andersen")


@pytest.mark.parametrize(("dt_fs", "stable"), [(3.1, True), (3.2, False)])
def test_diatomic_is_bounded_below_the_step_limit_and_explodes_above(dt_fs, stable):
    assert theory.harmonic_step_limit(9.8) == pytest.approx(3.12, abs=0.01)
    system, recipe = diatomic()
    sim = Simulation(system, recipe.update(control={"dt": dt_fs * 1e-3}))
    total = EnergyTimeseries(sim.run(200, ene_freq=1)).total
    growth = np.ptp(total) / abs(total[0])
    if stable:
        assert growth < 100
    else:
        assert not np.isfinite(growth) or growth > 1e6


# -- Check and Result -------------------------------------------------------------------------


def test_check_band_tolerance_and_messages():
    good = Check("slope", 2.01, 2.0, band=(1.9, 2.1))
    assert good and good.passed
    good.assert_passes()
    bad = Check("slope", 2.3, 2.0, band=(1.9, 2.1))
    assert not bad
    with pytest.raises(AssertionError, match=r"^slope: measured 2\.3, .*FAIL$"):
        bad.assert_passes()
    assert "PASS" in repr(good) and "FAIL" in repr(bad)
    assert Check("x", 1.01, 1.0, rel_tol=0.02)
    assert not Check("x", 1.03, 1.0, rel_tol=0.02)
    assert Check("x", 1.03, 1.0, rel_tol=0.02, abs_tol=0.05)
    assert not Check("x", float("nan"), 1.0, band=(-math.inf, math.inf))
    with pytest.raises(ValueError, match="needs a band"):
        Check("x", 1.0, 1.0)
    with pytest.raises(ValueError, match="no expected"):
        Check("x", 1.0, None, rel_tol=0.1)


# -- the experiments --------------------------------------------------------------------------


def _water(**kwargs) -> Result:
    return experiments.water_structure(
        str(REF_DIR / "water_216_box" / "water_216_box.topo"),
        str(REF_DIR / "water_216_box" / "water_216_box.conf"),
        str(REF_DIR / "water_216_nvt" / "water_216_nvt.in"),
        **kwargs,
    )


EXPERIMENTS = {
    "energy_flow": experiments.energy_flow,
    "time_step_scan": experiments.time_step_scan,
    "energy_landscape": experiments.energy_landscape,
    "thermostat_fingerprint": experiments.thermostat_fingerprint,
    "block_error": experiments.block_error,
    "maxwell_boltzmann": experiments.maxwell_boltzmann,
    "liquid_structure": experiments.liquid_structure,
    "water_structure": _water,
    "heat_bath": experiments.heat_bath,
}


@pytest.fixture(scope="module", params=list(EXPERIMENTS))
def result(request) -> Result:
    return EXPERIMENTS[request.param]()


def test_experiment_checks_pass(result):
    assert result.checks
    for check in result.checks:
        check.assert_passes()
    returned = result.check()
    assert returned == (result.checks[0] if len(result.checks) == 1 else result.checks)
    assert all(isinstance(value, np.ndarray) for value in result.data.values())


def test_experiment_plot_returns_a_figure(result):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    fig = result.plot()
    assert isinstance(fig, Figure)
    assert fig.axes and all(ax.has_data() for ax in fig.axes)
    plt.close(fig)


def test_water_structure_on_the_reference_box():
    r = _water(equilibration_steps=2500, frames=20)
    peak = r.data["r"][np.argmax(r.data["g"])]
    assert 0.26 <= peak <= 0.29
    assert 0.9 <= float(r.data["tail"]) <= 1.1


def test_heat_to_bath_integrates_an_exponential_by_hand():
    t0, delta, tau, dof = 120.0, 100.0, 0.5, 126
    time = np.linspace(0.0, 3.0, 3001)
    heat = experiments._heat_to_bath(time, t0 + delta * np.exp(-time / tau), t0, tau, dof)
    expected = dof * K_B / (2 * tau) * delta * tau * (1 - math.exp(-time[-1] / tau))
    assert heat[0] == 0.0
    assert heat[-1] == pytest.approx(expected, rel=1e-3)


def test_heat_bath_short_run():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = experiments.heat_bath(steps=600, equilibration=200, n=32)
    assert isinstance(r, Result)
    assert len(r.checks) == 8
    assert all(math.isfinite(c.measured) for c in r.checks)
    keys = {"time", "taus", "berendsen", "nose_hoover", "total_energy", "heat_to_bath"}
    assert keys | {"kinetic_only", "tau_fit"} == set(r.data)
    assert r.data["berendsen"].shape == (2, len(r.data["time"]))
    plt.close(r.plot())


def test_experiments_take_their_measurements_from_the_engine():
    source = inspect.getsource(experiments)
    assert "np.histogram(" not in source
    for minimum_image in ("np.round(", "np.rint(", "np.floor(", "% box", "%box"):
        assert minimum_image not in source


def test_education_imports_without_anywidget():
    code = (
        "import sys; sys.modules['anywidget'] = None\n"
        "import gromos.education as e\n"
        "assert callable(e.lj_gas) and e.experiments and e.theory\n"
        "try:\n"
        "    e.view\n"
        "except ImportError:\n"
        "    print('viewer unavailable')\n"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "viewer unavailable"
