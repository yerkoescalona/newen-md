"""
`gromos.analysis.ext_ti_ana`, `bar`, `ext_ti_merge`, `mbar` and `ti_trapezoid` (item 1.8).

The programs' agreement with gromos++ (byte for byte) and of `mbar` with pymbar is pinned on the
Rust side (``crates/gromos-analysis/tests/gromospp_references.rs``, ``mbar_reference.rs``). These
tests check that the Python calls run the same code: the files ``ext_ti_ana`` writes are
gromos++'s, the numbers ``bar`` returns are the ones it prints, and ``mbar`` on two states is
``bar``. By default ``bar`` and ``ext_ti_merge`` print the corrected errors (``chain_fixed/``, the
patched gromos++); ``parity="gromospp"`` prints gromos++ 9cd9f6b3's (``chain/``).
"""

from pathlib import Path

import numpy as np
import pytest

from gromos import analysis

REPO_ROOT = Path(__file__).parent.parent.parent
DATA = REPO_ROOT / "crates" / "gromos-analysis" / "tests" / "data" / "ext_ti_ana"
CHAIN = DATA / "chain"
CHAIN_FIXED = DATA / "chain_fixed"
LAMS = ["0.00", "0.25", "0.50", "0.75", "1.00"]

pytestmark = pytest.mark.skipif(not DATA.exists(), reason="reference fixtures not in the checkout")


def window(lam: str, **kw):
    return analysis.ext_ti_ana(
        en_files=[DATA / "windows" / f"ch4_l{lam}.tre"],
        fr_files=[DATA / "windows" / f"ch4_l{lam}.trg"],
        library=DATA / "ene_ana.md++.lib",
        temp=300,
        nrlambdas=21,
        minlam=0,
        maxlam=1,
        slam=float(lam),
        nlams=1,
        **kw,
    )


def test_ext_ti_ana_writes_gromospps_files(tmp_path):
    lams = [float(x) for x in LAMS]
    r = window(
        "0.25",
        countframes=True,
        bar_data=lams,
        dhdl_data=lams,
        outdir=tmp_path,
    )
    assert r["n_frames"] == 20
    names = sorted(Path(f).name for f in r["files"])
    assert names == ["bar_data_l0.25.dat", "dHdl_data_l0.25.dat", "predicted_l0.25.dat"]
    for name in names:
        assert (tmp_path / name).read_text() == (CHAIN / name).read_text(), name
    p = r["predictions"][0]
    assert p["plam"].shape == (21,)
    # the printed column is the returned prediction, to the 6 decimals printed
    rows = [ln.split() for ln in (tmp_path / "predicted_l0.25.dat").read_text().splitlines()]
    printed = np.array([float(x[1]) for x in rows if not x[0].startswith("#")])
    np.testing.assert_allclose(p["dhdl"], printed, atol=5e-7)
    assert p["bar"]["e_p"].shape == (5, 20)
    assert p["dhdl_data"]["de"].shape == (5, 20)


def test_ext_ti_ana_without_writing_and_with_the_builtin_library(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r = analysis.ext_ti_ana(
        en_files=[DATA / "windows" / "ch4_l0.50.tre"],
        fr_files=[DATA / "windows" / "ch4_l0.50.trg"],
        temp=300,
        nrlambdas=21,
        minlam=0,
        maxlam=1,
        slam=0.5,
        nlams=1,
        write=False,
    )
    assert r["files"] == []
    assert list(tmp_path.iterdir()) == []
    assert r["predictions"][0]["bar"] is None


def test_bar_returns_what_it_prints(tmp_path):
    files = [CHAIN / f"bar_data_l{lam}.dat" for lam in reversed(LAMS)]
    r = analysis.bar(files, 300.0, convergence=1e-10, printdist=tmp_path)
    assert r["stdout"] == (CHAIN_FIXED / "bar.stdout").read_text()
    assert r["lam_i"] == [0.0, 0.25, 0.5, 0.75]
    rows = [ln.split() for ln in r["stdout"].splitlines()[1:5]]
    np.testing.assert_allclose(r["dg"], [float(x[2]) for x in rows], atol=5e-6)
    np.testing.assert_allclose(r["error"], [float(x[4]) for x in rows], atol=5e-6)
    assert r["error"] == r["error_corrected"]
    assert r["error_total"] == r["error_total_corrected"]
    g = analysis.bar(files, 300.0, convergence=1e-10, parity="gromospp")
    assert g["stdout"] == (CHAIN / "bar.stdout").read_text()
    assert g["error_corrected"] == r["error_corrected"] and g["error"] != r["error"]
    assert g["dg"] == r["dg"]
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "dist_0_1.dat",
        "dist_1_2.dat",
        "dist_2_3.dat",
        "dist_3_4.dat",
    ]
    assert (tmp_path / "dist_0_1.dat").read_text() == (CHAIN / "dist_0_1.dat").read_text()


def test_bar_bootstrap_with_a_seed_is_reproducible():
    files = [CHAIN / "bar_data_l0.00.dat", CHAIN / "bar_data_l0.25.dat"]
    a = analysis.bar(files, 300.0, bootstrap=10, seed=7)
    b = analysis.bar(files, 300.0, bootstrap=10, seed=7)
    assert a["bootstrap"] == b["bootstrap"]
    assert a["bootstrap"][0] > 0


def test_the_old_two_list_bar_call_is_refused():
    with pytest.raises(TypeError, match="ext_ti_ana @bar_data"):
        analysis.bar("forward.dat", "reverse.dat", temperature=300)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        analysis.bar(np.zeros(3), np.zeros(3))  # type: ignore[arg-type]


def test_ext_ti_merge_reads_the_count_column_as_the_error():
    files = [CHAIN / f"predicted_l{lam}.dat" for lam in ["0.50", "0.00", "1.00", "0.25", "0.75"]]
    r = analysis.ext_ti_merge(files)
    assert r["stdout"] == (CHAIN_FIXED / "merge.stdout").read_text()
    assert len(r["lambda"]) == 21
    assert r["integral_error"] == r["integral_error_corrected"]
    g = analysis.ext_ti_merge(files, parity="gromospp")
    assert g["stdout"] == (CHAIN / "merge.stdout").read_text()
    assert g["integral_error_corrected"] == r["integral_error_corrected"] > g["integral_error"]
    r2 = analysis.ext_ti_merge(files, noerrors=True)
    assert r2["error"] is None and r2["integral_error"] is None
    assert r2["integral"] == pytest.approx(r["integral"], abs=0)


def test_mbar_of_the_bar_data_files_and_its_two_state_case_is_bar():
    preds = [window(lam, bar_data=[float(x) for x in LAMS], write=False) for lam in LAMS]
    bars = [p["predictions"][0]["bar"] for p in preds]
    states = [float(x) for x in LAMS]
    n_k = [len(b["e_s"]) for b in bars]

    def u_kn(idx):
        rows = []
        for k in idx:
            row = []
            for l_ in idx:
                b = bars[l_]
                if k == l_:
                    row.append(b["e_s"])
                else:
                    row.append(b["e_p"][list(b["plam"]).index(states[k])])
            rows.append(np.concatenate(row))
        return np.array(rows)

    r = analysis.mbar(u_kn(range(5)), n_k, 300.0)
    ref = (
        REPO_ROOT / "crates" / "gromos-analysis" / "tests" / "data" / "mbar" / "pymbar_reference.dat"
    )
    for line in ref.read_text().splitlines():
        if line.startswith("#"):
            continue
        i, j, dg, ddg = line.split()
        assert r["dg"][int(i), int(j)] == pytest.approx(float(dg), abs=0.01)
        assert r["ddg"][int(i), int(j)] == pytest.approx(float(ddg), abs=0.01)

    two = analysis.mbar(u_kn([0, 1]), n_k[:2], 300.0, tolerance=1e-14)
    b = analysis.bar([CHAIN / "bar_data_l0.00.dat", CHAIN / "bar_data_l0.25.dat"], 300.0,
                     convergence=1e-14, max_iterations=1000)
    assert two["dg"][0, 1] == pytest.approx(b["dg"][0], abs=1e-6)


def test_ti_trapezoid_integrates_the_windows():
    trgs = [DATA / "windows" / f"ch4_l{lam}.trg" for lam in LAMS]
    r = analysis.ti_trapezoid(trgs, [float(x) for x in LAMS], skip=2)
    assert r["lambda"] == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert len(set(r["n"])) == 1 and r["n"][0] > 2
    expect = np.trapezoid(r["mean"], r["lambda"])
    assert r["dg"] == pytest.approx(expect, rel=1e-12)
    assert r["dg_error"] > 0


# ---- gromos.free_energy: λ windows run from Python ------------------------------------------

REF = REPO_ROOT / "crates" / "gromos-md" / "tests" / "gromosXX_references"
SHARED = REF / "shared"
# the ch4_water_fep_l* references: 10 steps, dH/dλ every step; λ = 0.5 has no .trg reference
CH4_WINDOWS = {
    0.0: "ch4_water_fep_l000",
    0.25: "ch4_water_fep_l025",
    0.5: None,
    0.75: "ch4_water_fep_l075",
    1.0: "ch4_water_fep_l100",
}


def _trg_dhdl(path: Path) -> list[float]:
    lines = path.read_text().splitlines()
    return [float(lines[i + 1]) for i, ln in enumerate(lines) if ln.strip() == "# totals"]


def _ch4_base():
    from gromos import Recipe

    return (
        Recipe.from_imd(str(REF / "ch4_water_fep" / "ch4_water_fep.in"))
        .with_inputs(pttopo=str(SHARED / "ch4_spc_dummy.ptp"))
        .update(control={"steps": 10}, outputs={"ntwg": 1})
    )


@pytest.mark.skipif(not REF.exists(), reason="reference systems not in the checkout")
def test_lambda_windows_reproduce_the_reference_dhdl(tmp_path):
    from gromos import Recipe, System
    from gromos.free_energy import LambdaSchedule, collect, run_windows

    schedule = LambdaSchedule(list(CH4_WINDOWS))
    recipes = schedule.recipes(_ch4_base())
    # each window's recipe is its reference's .imd (but the title), with its @pttopo
    for (lam, name), r in zip(CH4_WINDOWS.items(), recipes, strict=True):
        assert r.perturbation["lambda"] == lam
        if name:
            ref = (
                Recipe.from_imd(str(REF / name / f"{name}.in"))
                .with_inputs(pttopo=str(SHARED / "ch4_spc_dummy.ptp"))
                .to_dict()
            )
            ours = r.to_dict()
            ref.pop("title"), ours.pop("title")
            assert ours == ref, name

    system = System.from_files(str(SHARED / "ch4_spc.top"), str(SHARED / "ch4_spc.cnf"))
    windows = run_windows(system, recipes, workdir=tmp_path)
    assert windows.lambdas == list(CH4_WINDOWS)
    assert sorted(p.name for p in tmp_path.iterdir()) == [f"{n}.npz" for n in schedule.names]
    trgs, lams = [], []
    for w, (lam, name) in zip(windows, CH4_WINDOWS.items(), strict=True):
        assert len(w.dhdl) == 10
        if not name:
            continue
        trg = REF / name / "expected" / "free_energy.trg"
        expected = _trg_dhdl(trg)
        for i, (a, e) in enumerate(zip(w.dhdl, expected, strict=True)):
            assert abs(a - e) <= max(abs(e) * 1e-6, 1e-8), f"{name} frame {i}: {a} vs {e}"
        assert np.mean(w.dhdl) == pytest.approx(np.mean(expected), rel=1e-7)
        trgs.append(trg)
        lams.append(lam)

    dg, err = windows.ti()
    assert np.isfinite(dg) and err > 0
    assert "ΔG (TI, trapezoid)" in windows.report()
    # the estimator is the one of `analysis.ti_trapezoid` on the reference .trg files
    four = run_windows(system, [r for r, n in zip(recipes, CH4_WINDOWS.values(), strict=True) if n])
    r = analysis.ti_trapezoid(trgs, lams, skip=2)
    assert four.ti(skip=2) == pytest.approx((r["dg"], r["dg_error"]), rel=1e-6)

    # kept windows are read back, not run again
    again = run_windows(system, recipes, workdir=tmp_path, dry_run=True)
    kept = collect(tmp_path)
    for a, b, c in zip(windows, again, kept, strict=True):
        assert np.array_equal(a.dhdl, b.dhdl) and np.array_equal(a.dhdl, c.dhdl)
    assert kept.mean_dhdl() == windows.mean_dhdl()


def test_lambda_schedule_sets_lambda_and_the_interaction_lambdas():
    from gromos.free_energy import LambdaSchedule, run_windows

    rows = [{"kind": "vdw", "g1": 1, "g2": 2, "a": 0.0, "b": 0.0, "c": 0.5, "d": 0.5, "e": 0.0}]
    recipes = LambdaSchedule([0.2, 0.8], interaction_lambdas=rows).recipes(_ch4_base())
    assert [r.perturbation["lambda"] for r in recipes] == [0.2, 0.8]
    assert all(r.perturbation["lambdas"] == {"enabled": True, "rows": rows} for r in recipes)
    with pytest.raises(ValueError, match="outside"):
        LambdaSchedule([0.0, 1.5])
    with pytest.raises(ValueError, match="resolve"):
        LambdaSchedule([0.5, 0.5])
    plain = _ch4_base().update(perturbation={"enabled": 0})
    with pytest.raises(Exception, match="no perturbation"):
        LambdaSchedule([0.5]).recipes(plain)
    # a dry run without a workdir runs nothing
    jobs = run_windows(None, recipes, dry_run=True)  # type: ignore[arg-type]
    assert [w.name for w in jobs] == ["l0.200", "l0.800"]
    assert all(w.energies is None and w.recipe is not None for w in jobs)


def test_apply_perturbation_returns_the_net_charges(tmp_path):
    from gromos import Topology

    top = Topology(str(SHARED / "ch4_spc.top"))
    assert top.apply_perturbation(str(SHARED / "ch4_spc_dummy.ptp")) == (0.0, 0.0)
    ptp = tmp_path / "plus.ptp"
    ptp.write_text(
        "TITLE\ncharge\nEND\nPERTATOMPARAM\n1\n    1 1 CM 17 16.043 1.0 17 16.043 0.0 0.0 0.0\nEND\n"
    )
    a, b = Topology(str(SHARED / "ch4_spc.top")).apply_perturbation(str(ptp))
    assert (a, b) == pytest.approx((1.0, 0.0))
