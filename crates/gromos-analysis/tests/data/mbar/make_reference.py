"""pymbar's free energies for the five CH4 windows in ../ext_ti_ana/chain (gromos-rs `mbar`'s
reference, item 1.8). Needs pymbar (not a dependency of gromos-rs; `uv pip install pymbar`).

    python make_reference.py > pymbar_reference.dat

The reduced potentials are built as `mbar::reduced_potentials` builds them: every
`bar_data_l<lam>.dat` contributes its frames, at its own lambda the E_s column and at the other
simulated lambdas their E_p columns, all divided by kT with gromos++'s Boltzmann constant.
"""

from pathlib import Path

import numpy as np
from pymbar import MBAR

CHAIN = Path(__file__).resolve().parent.parent / "ext_ti_ana" / "chain"
KBT = 1.380658e-23 * 6.0221367e23 / 1000.0 * 300.0


def read(path):
    rows = [ln.split() for ln in path.read_text().splitlines() if not ln.startswith("#")]
    nr = int(rows[0][0])
    lam = [float(x) for x in rows[1][: nr + 1]]
    data = np.array([[float(x) for x in r[: nr + 1]] for r in rows[2:]])
    return lam[0], lam[1:], data


files = [read(p) for p in sorted(CHAIN.glob("bar_data_l*.dat"))]
files.sort(key=lambda f: f[0])
states = [f[0] for f in files]
u_kn = np.zeros((len(states), sum(len(f[2]) for f in files)))
n_k = np.array([len(f[2]) for f in files])
col = 0
for slam, plam, data in files:
    n = len(data)
    for k, lk in enumerate(states):
        u_kn[k, col : col + n] = (data[:, 0] if lk == slam else data[:, 1 + plam.index(lk)]) / KBT
    col += n
r = MBAR(u_kn, n_k, solver_protocol="robust").compute_free_energy_differences()
print("# pymbar", __import__("pymbar").__version__, "on", ", ".join(f"{s:.2f}" for s in states))
print("#  i  j   DG (kJ/mol)            dDG (kJ/mol)")
for i in range(len(states)):
    for j in range(i + 1, len(states)):
        print(f"{i:3d}{j:3d} {r['Delta_f'][i, j] * KBT:.12e} {r['dDelta_f'][i, j] * KBT:.12e}")
