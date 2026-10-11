"""
`gromos.rdf` — one frame of gromos++ `rdf` (`gromos_analysis::rdf`).

Checked on bulk SPC water against a numpy minimum-image reference with gromos++'s
normalisation: density of the other oxygens, (N - 1)/V, since every centre is also
counted around, and shell volume 4 pi r_k^2 dr at the bin centre r_k.

The configuration is `water_216_box`'s `expected/final.conf` (100 MD steps on): the
starting `.conf` is a simple-cubic lattice, whose g(r) is a set of spikes.
"""

from pathlib import Path

import numpy as np
import pytest

from gromos import System, rdf

REPO_ROOT = Path(__file__).parent.parent.parent
WATER = REPO_ROOT / "crates" / "gromos-md" / "tests" / "gromosXX_references" / "water_216_box"
TOPO = WATER / "water_216_box.topo"
CONF = WATER / "expected" / "final.conf"

N_BINS = 45
R_MAX = 0.9

pytestmark = pytest.mark.skipif(
    not CONF.exists(), reason="gromosXX reference systems not in the checkout"
)


@pytest.fixture(scope="module")
def water():
    system = System.from_files(str(TOPO), str(CONF))
    oxygens = list(range(0, system.n_atoms, 3))
    return system, oxygens


def numpy_rdf(positions, oxygens, box):
    """Histogram of minimum-image O-O distances, every pair counted from both centres."""
    o = positions[oxygens]
    lengths = np.asarray(box)
    d = o[:, None, :] - o[None, :, :]
    d -= lengths * np.round(d / lengths)
    r = np.linalg.norm(d, axis=-1)[~np.eye(len(o), dtype=bool)]
    dr = R_MAX / N_BINS
    r = r[r < R_MAX]
    counts = np.bincount((r / dr).astype(int), minlength=N_BINS)[:N_BINS]
    centres = (np.arange(N_BINS) + 0.5) * dr
    density = (len(o) - 1) / np.prod(lengths)
    shell = 4.0 * np.pi * centres**2 * dr
    return centres, counts / (len(o) * density * shell)


def test_every_bin_matches_numpy_minimum_image(water):
    system, oxygens = water
    r, g = rdf(system.positions, oxygens, oxygens, n_bins=N_BINS, r_max=R_MAX, box=system.box)
    r_ref, g_ref = numpy_rdf(system.positions, oxygens, system.box)
    assert r.dtype == np.float64 and g.dtype == np.float64
    np.testing.assert_allclose(r, r_ref, rtol=1e-12)
    np.testing.assert_allclose(g, g_ref, rtol=1e-9, atol=0.0)


def test_first_peak_and_long_range_limit(water):
    system, oxygens = water
    r, g = rdf(system.positions, oxygens, oxygens, n_bins=N_BINS, r_max=R_MAX, box=system.box)
    assert 0.26 <= r[np.argmax(g)] <= 0.29
    assert 0.9 <= g[r > 0.6].mean() <= 1.1


def test_box_is_required(water):
    system, oxygens = water
    with pytest.raises(TypeError):
        rdf(system.positions, oxygens, oxygens, N_BINS, R_MAX)  # type: ignore[call-arg]
