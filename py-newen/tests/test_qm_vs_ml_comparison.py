"""QM-vs-ML comparison, from Python — the direct follow-up to `test_ml_potential.py`.

`crates/newen-md/tests/qm_vs_ml_comparison.rs` already does this in Rust (real `XtbInteraction`
vs a trained `SchNetInteraction`, RMSE check). Until now there was no way to get a real QM
reference value from Python at all — `SchNetPotential` could be attached to a `Simulation`, but
nothing exposed `XtbInteraction`. `XtbPotential` (this session's follow-up to P3.7) closes that.

Skips gracefully if `xtb` isn't on `PATH`. Uses the real trained model from
`scripts/train_qmmm_schnet.py` (`/tmp/trained_water_schnet.pt`) if present — the same file
`qm_vs_ml_comparison.rs` uses — for a real, meaningful RMSE check; falls back to a fresh
untrained toy model (`export_toy_schnet.py`) otherwise, with a much weaker claim (finite,
right-shaped output — not a physics check), same honesty tier as the rest of this repo's `ml`
tests.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
TRAINED_MODEL_PATH = Path("/tmp/trained_water_schnet.pt")

try:
    from newen import SchNetPotential, XtbPotential

    HAS_ML = True
except ImportError:
    HAS_ML = False

requires_ml = pytest.mark.skipif(
    not HAS_ML, reason="extension not built with --features ml (SchNetPotential absent)"
)


def _xtb_available() -> bool:
    return shutil.which("xtb") is not None


# Same held-out water geometry `qm_vs_ml_comparison.rs` uses — a deliberately different, larger
# perturbation than any training trajectory the generator produces.
WATER_POSITIONS_NM = np.array(
    [
        [0.012, -0.008, 0.005],
        [0.0758602 + 0.015, 0.010, 0.0504284 - 0.006],
        [0.0758602 - 0.009, -0.004, -0.0504284 + 0.011],
    ]
)
WATER_ELEMENTS = [8, 1, 1]


@requires_ml
def test_qm_and_ml_potentials_agree_on_held_out_water_geometry(tmp_path):
    if not _xtb_available():
        pytest.skip("xtb not found on PATH")

    xtb = XtbPotential(str(tmp_path / "xtb_work"), WATER_ELEMENTS, gfn=2, charge=0, multiplicity=1)
    e_qm, f_qm = xtb.evaluate(WATER_POSITIONS_NM)
    assert np.all(np.isfinite(f_qm))
    assert f_qm.shape == (3, 3)

    if TRAINED_MODEL_PATH.exists():
        # Real trained model (scripts/train_qmmm_schnet.py) — a genuine physics check, same
        # tolerances qm_vs_ml_comparison.rs uses.
        ml = SchNetPotential(str(TRAINED_MODEL_PATH), 1.0, WATER_ELEMENTS)
        e_ml, f_ml = ml.evaluate(WATER_POSITIONS_NM)

        energy_err = abs(e_qm - e_ml)
        force_rmse = float(np.sqrt(np.mean((f_qm - f_ml) ** 2)))
        print(
            f"QM/MM vs ML/MM (Python): E_qm={e_qm:.3f} E_ml={e_ml:.3f} kJ/mol "
            f"(|diff|={energy_err:.3f}), force RMSE={force_rmse:.3f} kJ/mol/nm"
        )
        # Generous, documented tolerances (matches qm_vs_ml_comparison.rs) — proves the trained
        # model tracks xtb's PES from Python too, not chemical accuracy from a toy-scale,
        # CPU-trained model on ~900 frames.
        assert energy_err < 200.0
        assert force_rmse < 5000.0
    else:
        # No trained model available — fall back to a fresh untrained toy model. Only checks the
        # plumbing (finite, right-shaped output), not agreement with QM — an untrained model has
        # no reason to agree with anything.
        out = tmp_path / "toy_schnet.pt"
        subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "export_toy_schnet.py"), str(out)],
            check=True,
            cwd=REPO_ROOT,
        )
        ml = SchNetPotential(str(out), 1.0, WATER_ELEMENTS)
        e_ml, f_ml = ml.evaluate(WATER_POSITIONS_NM)
        assert np.isfinite(e_ml)
        assert np.all(np.isfinite(f_ml))
        assert f_ml.shape == (3, 3)
