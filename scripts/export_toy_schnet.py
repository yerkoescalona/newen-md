#!/usr/bin/env python3
"""Export a small, randomly-initialized SchNet model — built from the real
SchNetPack library — to TorchScript.

Earlier version of this script hand-rolled a single cfconv block by hand to
prove the `PotentialProvider` seam (`crates/newen-forces/src/nonbonded/
schnet.rs`) could load *some* TorchScript model. That was necessary at the
time: GROMOS's own real trained tutorial model (SchNetPack v1, see
`.local/gromos_tutorial_livecoms/tutorial_files/t_06`) fails to
`torch.jit.script` — legacy Python-style control flow SchNetPack v1's own
code, not something fixable from the calling side.

**SchNetPack 2.x fixes this.** Verified directly: a `schnetpack.model.
NeuralNetworkPotential` wrapping `schnetpack.representation.SchNet` +
`schnetpack.atomistic.Forces` (the module that registers the energy→forces
autograd derivative) scripts cleanly with `torch.jit.script`, once
`Forces(calc_forces=True)` is included as an output module (a bare
`Atomwise`-only model fails TorchScript's type inference on an empty
`required_derivatives` list — `Forces` populates it). So: still untrained
(random init, no scientific-accuracy claim), but now the *actual* SchNetPack
architecture and code, not a reimplementation, and forces come from the
model's own internal autograd (the Rust side no longer calls `.backward()`
itself — see `schnet.rs`).

**Input/output contract differs from the old hand-rolled model — Dict[str,
Tensor], not positional args**, matching `schnetpack.properties`:
  - in:  `_positions` [N,3], `_atomic_numbers` [N] long, `_idx_i`/`_idx_j`
    [E] long (host-supplied edges — SchNetPack's `PairwiseDistances` does
    NOT compute its own neighbor list, it expects one as input, matching
    this project's `SpatialIndex` design exactly), `_offsets` [E,3] (the
    periodic-image shift per edge — the same quantity `SpatialIndex::
    neighbor_pairs` already returns), `_n_atoms` [1] long, `_idx_m` [N] long
    (which system each atom belongs to; all zeros for single-system
    inference — SchNetPack batches multiple systems together in general).
  - out: `energy` [1], `forces` [N,3] — both already in the output dict, no
    separate backward pass needed by the caller.

Requires `pip install schnetpack` (pulls in `pytorch-lightning`, heavier
than the previous plain-`torch` script) in the same `torch==2.11.0` venv
used to build newen-forces with `--features ml`. Run from the repo root:

    python3 scripts/export_toy_schnet.py [output_path]
"""

import sys
from pathlib import Path

import torch
import schnetpack as spk


def build_model(cutoff: float = 1.0, n_atom_basis: int = 16, n_interactions: int = 1) -> torch.nn.Module:
    pairwise_distance = spk.atomistic.PairwiseDistances()
    radial_basis = spk.nn.GaussianRBF(n_rbf=12, cutoff=cutoff)
    schnet = spk.representation.SchNet(
        n_atom_basis=n_atom_basis,
        n_interactions=n_interactions,
        radial_basis=radial_basis,
        cutoff_fn=spk.nn.CosineCutoff(cutoff),
    )
    pred_energy = spk.atomistic.Atomwise(n_in=n_atom_basis, output_key="energy")
    pred_forces = spk.atomistic.Forces(calc_forces=True, energy_key="energy", force_key="forces")

    model = spk.model.NeuralNetworkPotential(
        representation=schnet,
        input_modules=[pairwise_distance],
        output_modules=[pred_energy, pred_forces],
    )
    model.eval()
    return model


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("toy_schnet.pt")

    torch.manual_seed(0)
    model = build_model()

    scripted = torch.jit.script(model)
    scripted.save(str(out_path))
    print(f"wrote {out_path} (torch {torch.__version__}, schnetpack {spk.__version__})")


if __name__ == "__main__":
    main()
