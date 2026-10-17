"""ML potential binding (3.7): `resolve_zone_partition` and
`Simulation(..., ml_potential=, ml_region=, ml_buffer=)`.

Skips gracefully (not a failure) if:
  - the extension wasn't built with `--features ml` (`SchNetPotential`/
    `resolve_zone_partition` simply aren't importable), or
  - the real BuRNN `t_06` tutorial data isn't cloned (gitignored `.local/`,
    same guard every Rust test that touches it already uses), or
  - no exported TorchScript model is available (runs
    `scripts/export_toy_schnet.py` on the fly — same untrained-toy tier
    every other `ml` test in this repo uses; no chemical-accuracy claim).
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
REF_DIR = REPO_ROOT / "crates" / "newen-md" / "tests" / "gromosXX_references"
TUTORIAL_ROOT = REPO_ROOT / ".local" / "gromos_tutorial_livecoms" / "tutorial_files" / "t_06"

try:
    from newen import SchNetPotential, Simulation, Topology, resolve_zone_partition

    HAS_ML = True
except ImportError:
    HAS_ML = False

requires_ml = pytest.mark.skipif(
    not HAS_ML,
    reason="extension not built with --features ml (SchNetPotential/resolve_zone_partition absent)",
)


def _toy_model_path(tmp_path: Path) -> Path:
    """Export a fresh untrained toy SchNetPack model via the repo's own script — same
    architecture `nonbonded/schnet.rs` loads, no chemical-accuracy claim."""
    out = tmp_path / "toy_schnet.pt"
    subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "export_toy_schnet.py"), str(out)],
        check=True,
        cwd=REPO_ROOT,
    )
    assert out.exists()
    return out


@requires_ml
def test_resolve_zone_partition_matches_real_t06_solute_zone():
    """Part A — the exit criterion 3.7 specifies: zone *indices* from a selector string
    must match the real fixture `zone_partition_reference.rs` already validates.

    `t_06`'s real inner zone (the BuRNN QM zone, `QMZONE` in `md_burnn/meoh.qmm`) is exactly the
    6-atom methanol solute — indices 0-5, confirmed by `zone_partition_reference.rs`'s own
    `QM_ZONE_SIZE = 6` constant and `(0..QM_ZONE_SIZE).collect()` inner list. That real fixture's
    *buffer* zone is computed geometrically (atoms within a radius each step), not by a static
    selector string, so it isn't reproducible via `from_selections` and isn't compared here —
    this test scopes to the inner zone, which `"m:a"` (all solute atoms) resolves to exactly,
    since t_06's entire solute *is* the methanol QM zone.
    """
    topo_path = TUTORIAL_ROOT / "topo" / "meoh_54a7.top"
    if not topo_path.exists():
        pytest.skip(f"real BuRNN tutorial not cloned at {topo_path}")

    topo = Topology(str(topo_path))
    inner, buffer, outer = resolve_zone_partition(topo, "m:a", None)

    assert inner == list(range(6)), f"expected the 6-atom methanol solute, got {inner}"
    assert buffer == []
    assert len(outer) == topo.n_atoms - 6


@requires_ml
def test_simulation_with_ml_potential_runs_a_step(tmp_path):
    """Part B — `ml_potential=`/`ml_region=` kwargs actually reach the orchestrator-aware step
    path, not just get parsed. Real small system (`water_single`, 3 atoms, real bond+angle+
    Coulomb — the Rust reference suite's own Level-1 system), toy (untrained) SchNet model
    attached to the whole molecule via a real residue-name selector.

    Narrower than an energy-conservation claim (matching this feature's own "still untrained
    toy model" honesty tier, same as `schnet_nve_loop.rs`/every other `ml` test here) — this
    checks the simulation constructs and steps without error, i.e. the ML term is really wired
    into the step loop.
    """
    from newen import InputParameters, System

    sys_dir = REF_DIR / "water_single"
    if not sys_dir.exists():
        pytest.skip(f"reference system not found at {sys_dir}")

    model_path = _toy_model_path(tmp_path)
    potential = SchNetPotential(str(model_path), 1.0, [8, 1, 1])

    system = System.from_files(
        str(sys_dir / "water_single.topo"), str(sys_dir / "water_single.conf")
    )
    params = InputParameters(str(sys_dir / "water_single.in"))

    sim = Simulation(
        system,
        params,
        ml_potential=potential,
        ml_region="1:res(SOL:a)",
    )
    sim.step(1)
    assert sim.n_atoms == 3


def test_ml_kwargs_require_both_potential_and_region():
    """Even without --features ml, `Simulation` should reject a lopsided ml_potential/ml_region
    pair with a clear error, not a confusing downstream failure — this is a real check
    regardless of build (`resolve_ml_spec` validates the pairing before any feature-gated code
    runs) — but only meaningful once `Simulation`/`System`/`InputParameters` are importable."""
    try:
        from newen import InputParameters, Simulation, System
    except ImportError:
        pytest.skip("newen extension not built")

    sys_dir = REF_DIR / "water_single"
    if not sys_dir.exists():
        pytest.skip(f"reference system not found at {sys_dir}")

    system = System.from_files(
        str(sys_dir / "water_single.topo"), str(sys_dir / "water_single.conf")
    )
    params = InputParameters(str(sys_dir / "water_single.in"))

    with pytest.raises(ValueError, match="ml_potential and ml_region must be given together"):
        Simulation(system, params, ml_region="1:res(SOL:a)")
