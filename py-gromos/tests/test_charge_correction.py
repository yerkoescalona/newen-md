"""
`gromos.analysis.charge_correction` — the analytic correction terms of Kastenholz &
Hünenberger, J. Chem. Phys. 124, 224501 (2006), for a charge change under periodic boundaries.

The paper's own numbers are the oracle: Tables III (LS, CM) and IV (BM, BA) list every term for
Na⁺ (R_I = 0.2 nm) in SPC at its box sizes, L = [(N_w + 1)/ρ]^(1/3) with ρ = 33.427 nm⁻³
(Sec. III). The Rust tests (`gromos_analysis::free_energy::charge_correction`) check every row;
these check a row of each scheme through Python and what the binding adds.
"""

import pytest

from gromos import analysis

RHO = 33.427


def edge(n_water: int) -> float:
    return float(((n_water + 1) / RHO) ** (1 / 3))


# (scheme, N_w, R_C, A1, A2, B, C1, C2, tolerance), kJ/mol. The CM A2 converges 0.07-0.25
# below Table III (see the Rust test `a2_cutoff_coulomb_matches_table_iii_to_a_quarter_kj`),
# hence its tolerance.
ROWS = [
    ("LS", 4, None, 0.0, 0.0, -293.17, -61.89, -3.54, 0.01),
    ("LS", 1024, None, 0.0, 0.0, -61.79, -79.67, -0.02, 0.01),
    ("CM", 2048, 1.6, -42.86, 13.54, -1.49, 0.0, 0.0, 0.1),
    ("CM", 512, 0.8, -85.72, 24.48, -0.53, 0.0, 0.0, 0.26),
    ("RF", 1024, 1.4, -48.98, -24.54, -4.59, -77.75, 0.0, 0.1),
    ("BM", 512, 0.8, -85.72, -42.36, -1.47, -76.76, 0.0, 0.1),
    ("BA", 2048, 1.6, -42.86, -21.52, -2.82, -79.60, 0.0, 0.1),
]


@pytest.mark.parametrize("scheme,nw,rc,a1,a2,b,c1,c2,tol", ROWS)
def test_charge_correction_reproduces_the_papers_tables(scheme, nw, rc, a1, a2, b, c1, c2, tol):
    eps_rf = None if scheme in ("LS", "CM") else 66.6
    c = analysis.charge_correction(1.0, edge(nw), rc, eps_rf, "SPC", scheme=scheme)
    for key, paper in zip(("a1", "a2", "b", "c1", "c2"), (a1, a2, b, c1, c2), strict=True):
        assert c[key] == pytest.approx(paper, abs=tol), key
    assert c["total"] == pytest.approx(c["a1"] + c["a2"] + c["b"] + c["c1"] + c["c2"], abs=1e-12)
    assert c["beyond_pb"] == pytest.approx(c["c1"] + c["c2"], abs=1e-12)
    assert c["quadrupole_trace"] == pytest.approx(0.0082)
    assert c["delta_phi_odl"] == pytest.approx(79.76, abs=0.005)


def test_charge_correction_of_a_change_from_a_charged_state():
    full = analysis.charge_correction(1.0, 3.0, 1.4, 0.0)
    second_half = analysis.charge_correction(0.5, 3.0, 1.4, 0.0, charge_a=0.5)
    for key in ("a1", "a2", "b"):
        assert second_half[key] == pytest.approx(0.75 * full[key], rel=1e-9)
    assert second_half["c1"] == pytest.approx(0.5 * full["c1"], rel=1e-12)
    # A removal is the change undone.
    back = analysis.charge_correction(-1.0, 3.0, 1.4, 0.0, charge_a=1.0)
    assert back["total"] == pytest.approx(-full["total"], rel=1e-12)


def test_charge_correction_other_water_models():
    spc = analysis.charge_correction(1.0, 3.0, 1.4, 66.6)
    tip3p = analysis.charge_correction(1.0, 3.0, 1.4, 66.6, "TIP3P", eps_model=66.6)
    assert tip3p["quadrupole_trace"] == pytest.approx(2 * 0.417 * 0.09572**2)
    assert tip3p["c1"] / spc["c1"] == pytest.approx(tip3p["quadrupole_trace"] / 0.0082)
    for key in ("a1", "a2", "b"):
        assert tip3p[key] == pytest.approx(spc[key], rel=1e-12)
    own = analysis.charge_correction(1.0, 3.0, 1.4, 66.6, quadrupole_trace=0.01)
    assert own["c1"] / spc["c1"] == pytest.approx(0.01 / 0.0082)


@pytest.mark.parametrize(
    "args,kwargs,message",
    [
        ((1.0, 3.0, 1.4, 66.6, "TIP5P"), {}, "unknown solvent"),
        ((1.0, 3.0, 1.4, 66.6, "TIP3P"), {}, "eps_model"),
        ((1.0, 3.0, None, None, "TIP3P"), {"scheme": "LS"}, "Eq. 47"),
        ((1.0, 3.0, 1.4, None), {}, "eps_rf"),
        ((1.0, 3.0, None, 66.6), {}, "cutoff"),
        ((1.0, 3.0, 1.4, 66.6), {"scheme": "GM"}, "unknown electrostatics scheme"),
        ((1.0, 3.0, 0.1, 66.6), {}, "ion radius"),
    ],
)
def test_charge_correction_refuses_what_it_has_no_formula_for(args, kwargs, message):
    with pytest.raises(ValueError, match=message):
        analysis.charge_correction(*args, **kwargs)


@pytest.mark.parametrize(("eps_rf", "want"), [(61.0, -73.825), (66.6, -73.875), (78.5, -73.959)])
def test_rf_self_energy_of_a_unit_charge(eps_rf, want):
    """−½ FPEPSI (1 − C_rf/2)/R_C for q 0 → 1 at R_C = 1.4 nm, κ = 0, ε_CS = 1."""
    assert analysis.rf_self_energy(0.0, 1.0, 1.4, eps_rf) == pytest.approx(want, abs=0.01)
    assert analysis.rf_self_energy([1.0], [0.0], 1.4, eps_rf) == pytest.approx(-want, abs=0.01)


def test_rf_self_energy_sums_atoms_and_refuses_mismatched_states():
    one = analysis.rf_self_energy(0.0, 1.0, 1.4, 61.0)
    two = analysis.rf_self_energy([0.5, -1.0], [1.0, 0.0], 1.4, 61.0)
    assert two == pytest.approx(one * (1.0 - 0.25 - 1.0), rel=1e-12)
    assert analysis.rf_self_energy(1.0, 0.0, 0.0, 61.0) == 0.0
    with pytest.raises(ValueError, match="charges in state A"):
        analysis.rf_self_energy([1.0], [0.0, 0.0], 1.4, 61.0)


def test_corrected_charging_dg_takes_the_self_term_off():
    c = analysis.charge_correction(1.0, 3.0, 1.4, 61.0)
    raw, self_term = -400.0, analysis.rf_self_energy(0.0, 1.0, 1.4, 61.0)
    assert analysis.corrected_charging_dg(raw, self_term, c) == pytest.approx(
        raw - self_term + c["total"], abs=1e-12
    )
    assert analysis.corrected_charging_dg(raw, self_term, c, pb_a=-1.0, pb_b=-78.0) == (
        pytest.approx(raw - self_term - 77.0 + c["beyond_pb"], abs=1e-12)
    )
    with pytest.raises(ValueError, match="both"):
        analysis.corrected_charging_dg(raw, self_term, c, pb_a=0.0)
