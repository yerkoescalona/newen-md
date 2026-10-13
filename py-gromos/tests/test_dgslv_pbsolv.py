"""
`gromos.analysis.dgslv_pbsolv` — gromos++ ``dGslv_pbsolv`` called in process.

The function returns the numbers the program computes and also hands back the table it prints
and its ``@verbose`` log, so these tests can check one against the other: every potential in the
returned arrays is the one printed in the table, every free energy the one in the log. The
program's agreement with gromos++ itself is pinned on the Rust side
(``crates/gromos-analysis/tests/gromospp_references.rs``); the grids here are coarse to keep the
tests fast.
"""

import re
from pathlib import Path

import pytest

import gromos
from gromos import analysis

REPO_ROOT = Path(__file__).parent.parent.parent
SHARED = REPO_ROOT / "crates" / "gromos-md" / "tests" / "gromosXX_references" / "shared"
ALADIP_TOP = SHARED / "aladip.topo"
ALADIP_CNF = SHARED / "aladip.conf"

pytestmark = pytest.mark.skipif(
    not ALADIP_TOP.exists(), reason="reference fixtures not in the checkout"
)

FD = ["npbc_slv", "npbc_vac", "pbc_slv", "pbc_vac"]
FFT = ["fft_ls_pbc", "fft_rf_pbc"]


def table(stdout: str) -> list[list[str]]:
    """The rows of the program's table (the header and the gather warning start with '#')."""
    return [line.split() for line in stdout.splitlines() if line and not line.startswith("#")]


def log_value(log: str, prefix: str) -> list[float]:
    return [float(line.split()[-1]) for line in log.splitlines() if line.startswith(prefix)]


def test_c1_at_a_coarse_grid_returns_the_programs_numbers(tmp_path):
    log_file = tmp_path / "c1.log"
    r = analysis.dgslv_pbsolv(
        topo=ALADIP_TOP,
        coord=ALADIP_CNF,
        pbc="r",
        atoms="1:a",
        atoms_to_charge="1:a",
        scheme="LS",
        eps_solv=66.6,
        gridspacing=0.3,
        probe_iac=4,
        verbose=log_file,
    )
    rows = table(r["stdout"])
    assert len(rows) == 12
    assert r["atom"] == list(range(1, 13))
    for i, row in enumerate(rows):
        assert int(row[0]) == i + 1
        assert int(row[1]) == r["residue"][i]
        assert row[2] == r["name"][i]
        assert row[3] == f"{r['charge'][i]:.4f}"
        for column, key in zip(row[4:], FD, strict=True):
            assert column == f"{r[key][i]:.7f}"
    assert r["fft_ls_pbc"] is None and r["dg_fft_rf_pbc"] is None

    # The free energies are the log's (printed there with 6 significant digits).
    assert r["log"] == log_file.read_text()
    assert log_value(r["log"], "# DGRESULT NPBC") == [pytest.approx(r["dg_npbc"], rel=1e-5)]
    assert log_value(r["log"], "# DGRESULT PBC") == [pytest.approx(r["dg_pbc"], rel=1e-5)]
    returned = log_value(r["log"], "# DG RETURNED")
    expected = [r["dg_npbc_vac"], r["dg_npbc_slv"], r["dg_pbc_vac"], r["dg_pbc_slv"]]
    assert returned == pytest.approx(expected, rel=1e-5)
    assert r["dg_npbc"] == r["dg_npbc_slv"] - r["dg_npbc_vac"]
    assert r["documented_correction"] == r["dg_npbc"] - r["dg_pbc"]
    assert set(r["iterations"]) == set(FD)


def test_rf_adds_the_two_fft_solves():
    r = analysis.dgslv_pbsolv(
        topo=ALADIP_TOP,
        coord=ALADIP_CNF,
        atoms="1:a",
        atoms_to_charge="1:4,5",
        scheme="RF",
        eps_solv=66.6,
        eps_rf=66.6,
        rcut=1.4,
        gridspacing=0.2,
        probe_iac=4,
    )
    rows = table(r["stdout"])
    assert len(rows) == 2
    for i, row in enumerate(rows):
        for column, key in zip(row[4:], FD + FFT, strict=True):
            assert column == f"{r[key][i]:.7f}"
    # Each FFT solve ends its log section with its free energy alone on a line (3 decimals).
    lines = r["log"].splitlines()
    assert len([line for line in lines if line.startswith("# Convergence after")]) == 2
    values = [line for line in lines if re.fullmatch(r"-?\d+\.\d+", line)]
    assert values == [f"{r['dg_fft_ls_pbc']:.3f}", f"{r['dg_fft_rf_pbc']:.3f}"]
    assert r["documented_correction"] == pytest.approx(
        r["dg_npbc"] - r["dg_pbc"] + r["dg_fft_ls_pbc"] - r["dg_fft_rf_pbc"], abs=1e-12
    )
    assert set(r["iterations"]) == set(FD + FFT)


def test_eps_rf_one_is_refused():
    with pytest.raises(gromos.ProgramError, match="@epsRF 1"):
        analysis.dgslv_pbsolv(
            topo=ALADIP_TOP,
            coord=ALADIP_CNF,
            atoms="1:a",
            atoms_to_charge="1:a",
            scheme="RF",
            eps_solv=66.6,
            eps_rf=1,
            rcut=1.4,
            gridspacing=0.2,
            probe_iac=4,
        )
