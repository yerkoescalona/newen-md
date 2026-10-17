"""
`newen.EnergyTrajectory` — a `.tre`/`.trg` read through the energy library.

What the class adds over `EnergyTimeseries` (the scalars a running `Simulation`
reports) is the *shape* the library declares: per-bath and per-energy-group tables,
with the layout established before any value is believed
(`docs/src/reference/energy-library.md`).
"""

from pathlib import Path

import numpy as np
import pytest

from newen import EnergyTrajectory

REPO_ROOT = Path(__file__).parent.parent.parent
REF_DIR = REPO_ROOT / "crates" / "newen-md" / "tests" / "gromosXX_references"
TRE = REF_DIR / "aladip_solvated" / "expected" / "energies.tre"
TRG = REF_DIR / "ch4_water_fep_l025" / "expected" / "free_energy.trg"

pytestmark = pytest.mark.skipif(
    not TRE.exists(), reason="gromosXX reference trajectories not in the checkout"
)


def test_reads_a_gromosxx_tre_without_a_library():
    tre = EnergyTrajectory(str(TRE), properties=["totene", "totkin", "totpot"])
    assert len(tre) == 10
    # 2023-04-15 is a version the reader knows, so the layout is checked, not assumed
    assert tre.version == "2023-04-15"
    assert tre.warnings == []
    assert tre.profile_version == 1
    assert tre.times[0] == 0.0
    assert tre["totene"].shape == (10,)
    # gromosXX composes the total as potential + kinetic + special
    special = EnergyTrajectory(str(TRE), properties=["totspecial"])["totspecial"]
    np.testing.assert_allclose(tre["totene"], tre["totkin"] + tre["totpot"] + special, rtol=1e-9)


def test_energy_group_pairs_sum_to_the_total():
    """The re-partition workflow: the per-pair table is the total, split up."""
    tre = EnergyTrajectory(str(TRE), properties=["totlj", "totcrf"])
    nonbonded = tre.table("NONBONDED")
    assert nonbonded.shape == (10, len(tre.group_pairs), 6)
    assert tre.group_pairs == [(0, 0), (0, 1), (1, 1)]  # 2 energy groups, upper triangle
    # the file prints 10 significant digits, so the sum meets the total to about that
    np.testing.assert_allclose(nonbonded[:, :, 0].sum(axis=1), tre["totlj"], atol=1e-8)
    np.testing.assert_allclose(nonbonded[:, :, 1].sum(axis=1), tre["totcrf"], atol=1e-8)


def test_long_form_is_ready_for_a_group_by():
    tre = EnergyTrajectory(str(TRE))
    long = tre.to_long("NONBONDED")
    n = 10 * len(tre.group_pairs) * 6
    assert {k: v.shape for k, v in long.items()} == {
        "frame": (n,),
        "time": (n,),
        "row": (n,),
        "col": (n,),
        "group_i": (n,),
        "group_j": (n,),
        "value": (n,),
    }
    # the solute-solute LJ of the last frame, found by selection rather than by index
    keep = (
        (long["frame"] == 9) & (long["col"] == 0) & (long["group_i"] == 0) & (long["group_j"] == 0)
    )
    assert long["value"][keep] == pytest.approx(tre.last("NONBONDED")[0, 0])


def test_only_the_subblocks_of_this_file_type_are_offered():
    tre = EnergyTrajectory(str(TRE))
    assert "ENER" in tre.subblocks
    assert "FREEENER" not in tre.subblocks  # declared by the library, but for a .trg
    with pytest.raises(KeyError):
        tre.table("FREEENER")


@pytest.mark.skipif(not TRG.exists(), reason="no reference .trg in the checkout")
def test_reads_a_trg_as_free_energy():
    trg = EnergyTrajectory(str(TRG), properties=["dvdl", "totfren"], free_energy=True)
    assert len(trg) == 10
    assert trg.table("FREEENER").shape == (10, 52, 1)
    np.testing.assert_allclose(trg["totfren"], trg.table("FREEENER")[:, 0, 0])


def test_a_property_that_was_not_requested_says_so():
    tre = EnergyTrajectory(str(TRE), properties=["totene"])
    with pytest.raises(KeyError, match="properties="):
        tre["totkin"]


def test_a_file_that_does_not_match_its_layout_is_refused(tmp_path):
    """One value short of the layout it claims: an error, not a shifted read."""
    text = TRE.read_text().splitlines(keepends=True)
    out, values = [], 0
    for line in text:
        stripped = line.strip()
        try:
            float(stripped)
            values += 1
            if values == 52:  # drop the last of ENER's 52 totals in every frame
                continue
        except ValueError:
            pass
        out.append(line)
    broken = tmp_path / "short.tre"
    broken.write_text("".join(out))

    with pytest.raises(ValueError, match="subblock ENER"):
        EnergyTrajectory(str(broken))
