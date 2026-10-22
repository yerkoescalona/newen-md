# Reference systems for real-world GROMOS input styles

The systems in the levels above were all generated the same way, and they share three properties
that no real GROMOS input has: **one temperature bath**, **one value per line** in the parameter
file, and **molecules that are never wrapped** across the periodic boundary. Running the GROMOS
LiveCoMS tutorial suite against newen on 2026-08-30 turned that blind spot into six defects
(CHANGELOG 0.0.42 and 0.0.43). These systems close it.

Each was produced with the same gromosXX binary as the rest of the suite
(`python3 scripts/regen_gromosXX_references.py <name>`), and each is a *style* the tutorials use.

| system | what it covers | how it would fail |
|---|---|---|
| `aladip_multibath` | `MULTIBATH NBATHS=2` with a `DOFSET` (solute bath + solvent bath), all bonds constrained | before 0.0.42 the run was refused outright ("newen supports one bath"); per-bath degrees of freedom, per-bath kinetic energy and per-bath scaling are all exercised here |
| `aladip_multibath_collapsed` | the *same run*, with `NONBONDED` and `CONSTRAINT` wrapped the way real files write them — `NLRELE` on the `APPAK` line, `TOLA2 = 1e-10`, the five `CONSTRAINT` values on one line | a line-indexed reader puts `TOLA2` in `NSLFEXCL`'s integer slot and loses `NTCP`; gromosXX reads every block as a value stream, and this file proves the layout is legal GROMOS |
| `aladip_wrapped` | the solute **split across the periodic boundary** (3.80 nm between two of its atoms in a 3.767 nm box), charge groups kept whole — exactly how GROMOS writes an equilibrated configuration | bonded terms or SHAKE without the minimum image see 3.8 nm "bonds": SHAKE cannot converge, and the angle energy comes out ~70× too large |

`aladip_multibath` and `aladip_wrapped` are the same physical system in two different images, and
gromosXX gives them **the same energies** — the wrapped configuration is a pure bookkeeping change.
That equality is the property the third row tests: any interaction that skips the minimum image
breaks it.

Two more followed as the tutorials were worked through:

| system | what it covers |
|---|---|
| `aladip_harmonic_covalent` | `COVALENTFORM 1 1 0` — harmonic bonds and harmonic angles instead of the quartic and cosine-harmonic forms. Both were implemented all along; nothing selected between them until 0.0.45, and the block was refused, which is what blocked the GaMD tutorial |
| `aladip_gamd` | a `GAMD` block (FORM 3, the dihedral boost; production mode, fixed parameters) with its `@gamd` file: one interaction group over all atoms |
| `aladip_gamd_dual` | `aladip_gamd` with FORM 1 and ET 0: the potential and the dihedral boost both act on the forces; ΔV is the dihedral one only, as gromosXX assigns it twice |
| `aladip_gamd_cmd_search` | `aladip_gamd` with SEARCH 1, NTIGAMDS 1: a cMD search from scratch, no boost; the dihedral statistics end in the final configuration's `GAMDSTAT` block |

`aladip_wrapped` also writes a velocity trajectory (`WRITETRAJ` NTWV), so `@trv` — accepted and
ignored until 0.0.44 — is exercised by a reference.

The comparison in `test_gromosXX_references.rs` also covers the **trajectory** (frame by frame,
positions, periodic images folded out) and the **final configuration** for every
system (`expected/final.conf`, positions and velocities). That is where a wrong step count shows —
the per-step frames look right either way — and it is how the NSTLIM+1 loop of 0.0.41 was found.
Positions are compared modulo a periodic image, because which image a file records is not physics;
NTB = −1 systems skip the coordinate comparison entirely, since newen writes the rotated
triclinic frame and gromosXX the truncated-octahedron one.
