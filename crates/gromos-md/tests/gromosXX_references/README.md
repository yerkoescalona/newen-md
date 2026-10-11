# gromosXX reference systems

One directory per system: the inputs, and `expected/` holding what the working upstream
(`.local/gromosXX`, md++ 2023-04-15) produced for them. `test_gromosXX_references.rs` runs each and
compares energies, trajectory and final configuration; the Python suite drives the same systems
through `Simulation` so both front-ends are held to one oracle.

**A system counts as passing when its numbers agree with md++ — step-0 energies to the printed
precision and the trajectory within the stated tolerances — not when the binary exits 0.**

Adding one: see `crates/gromos-md/.claude/CONTEXT.md` §"How to add a reference test".
`README-livecoms.md` covers the three systems that came from real LiveCoMS tutorial input.

## Status

Run them: `cargo test -p gromos-md --test test_gromosXX_references`

| Lvl | System           | Atoms | Isolates                              | Status   |
|-----|------------------|-------|---------------------------------------|----------|
| 0   | pair_lj          | 2     | Pure LJ, no PBC                      | **PASS** |
| 0   | pair_lj_mixed    | 2     | LJ combination rules                 | **PASS** |
| 0   | nacl_pair        | 2     | Coulomb + LJ ions                    | **PASS** |
| 1   | water_single     | 3     | bond + angle + intramolecular CRF    | **PASS** |
| 1   | water_single_genvel | 3  | NTIVEL=1 Maxwell-Boltzmann velocity generation | **PASS** |
| 1   | benzene_vacuum   | 12    | aromatic ring + improper + torsion   | **PASS** |
| 1   | nacl_pair_box    | 2     | Coulomb + LJ in PBC with RF (no solvent) | **PASS** |
| 1   | butane_vacuum    | 4     | dihedral + 1-4 LJ interaction        | **PASS** |
| 1   | aladip_vacuum    | 12    | all bonded + exclusions + 1-4        | **PASS** |
| 2   | water_3_box      | 9     | PBC + min image + pairlist + CRF     | **PASS** |
| 2   | nacl_1water_box  | 5     | minimal solute-solvent + SHAKE       | **PASS** |
| 2   | nacl_1water_settle | 5   | SETTLE (analytical rigid water)      | **PASS** |
| 2   | nacl_1water_lincs | 5    | LINCS (solvent)                      | **PASS** |
| 2   | nacl_3water_box  | 11    | multiple solvent + solute-solvent pairlist | **PASS** |
| 2   | water_3_box_twinrange | 9 | twin-range pairlist (RCUTP<RCUTL, NSNB=5) | **PASS** |
| 2   | water_10_box     | 32    | 2 ions + 10 SPC, positions away from cutoff | **PASS** |
| 2   | nacl_3water_cutoff | 11  | nacl_3water near cutoff boundary     | **PASS** |
| 2   | nacl_water_box   | 62    | ion-water RF in PBC                  | **PASS** |
| 2   | nacl_water_box_shifted | 62 | nacl_water_box with perturbed positions | **PASS** |
| 3   | water_216_box    | 648   | bulk NVE, pairlist, virial           | **PASS** |
| 3   | water_216_nve_nobath | 648 | absent MULTIBATH block ⇒ no bath (IMD parser defaults, 3.9 A18) | **PASS** |
| 3   | water_216_box_com| 648   | bulk NVE + COM removal (NTICOM=1, NSCM=10) | **PASS** |
| 3   | water_216_box_com_rot | 648 | COM translation+rotation removal (NTICOM=2, NSCM=-10) | **PASS** |
| 3   | water_216_nvt    | 648   | Berendsen thermostat                 | **PASS** |
| 3   | water_216_nvt_nosehoover | 648 | Nosé-Hoover thermostat (single NHC) | **PASS** |
| 3   | water_216_nvt_nhc_chain | 648 | Nosé-Hoover-Chain (3 chains)        | **PASS** |
| 3   | water_216_npt    | 648   | Berendsen barostat                   | **PASS** |
| 4   | aladip_vacuum_lincs | 12 | LINCS (solute, NTC=2)               | **PASS** |
| 4   | aladip_solvated  | 72    | SHAKE + solute-solvent               | **PASS** |
| 4   | aladip_vacuum_em | 12    | steepest descent EM, vacuum          | **PASS** |
| 4   | aladip_vacuum_em_shake | 12 | SD EM + SHAKE, vacuum             | **PASS** |
| 4   | aladip_solvated_em_noshake | 72 | SD EM, solvated, no SHAKE      | **PASS** |
| 4   | aladip_solvated_em_shake | 72 | SD EM + SHAKE, solvated          | **PASS** |
| 4   | aladip_solvated_em_posres | 72 | SD EM + position restraints     | **PASS** |
| 4   | aladip_solvated_em | 72  | SD EM + SHAKE + posres, solvated    | **PASS** |
| 2   | nacl_1water_distres | 5  | distance restraint on Na-Cl pair (NTDIR=2, CDIR*w0) | **PASS** |
| 4   | ch4_water_fep | 2998 | CH4→dummy in 999 SPC water, λ=0.5, twin-range NB FEP | **PASS** (its `.trg` is empty — NTWG=0; dH/dλ is covered by the rows below) |
| 4   | ch4_water_fep_l000 / l025 / l075 / l100 | 2998 | the same system at λ = 0, 0.25, 0.75, 1; ten steps; dH/dλ total, LJ, CRF, bonded per step | **PASS** (2026-08-29) |
| 4   | meoh_water_fep | 2862 | charged CH3OH→dummy in 953 SPC water, λ=0.5, ten steps, dH/dλ per term (soft-core RF for charged atoms) | **PASS** (2026-08-29, after the kernel fix) |
| 4   | aladip_vacuum_fep_groups | 12 | aladip vacuum FEP (RLAM=0.125) with two energy groups (1-6, 7-12), `.trg` every step: per-group bonded/nonbonded energies, per-term per-group bonded dH/dλ, kinetic dH/dλ of the perturbed masses (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_trg | 12 | aladip vacuum FEP with one energy group and a `.trg` every step: the single-group derivative cells carry the dH/dλ totals (1.8) | **PASS** |
| 4   | ch4_water_fep_lambdas | 2998 | `ch4_water_fep_l025` with `LAMBDAS` on: the vdw/vdw_soft polynomials on the solute-solvent pair (1,2) are what it tests; its crf/crf_soft rows are inert, CH4 being uncharged in both states (charged perturbed atoms under crf rows: `aladip_vacuum_fep_lambdas`, solute–solute only) (1.8) | **PASS** |
| 4   | ch4_water_fep_lambdas_off | 2998 | the same `LAMBDAS` block switched off: rows ignored, identical to `ch4_water_fep_l025` outside the title (1.8) | **PASS** |
| 4   | ch4_water_fep_nlam2 | 2998 | NLAM = 2: `(1−λ)ⁿ`, `λⁿ` on the soft-core nonbonded kernel only (1.8) | **PASS** |
| 4   | ch4_water_fep_nrdgl | 2998 | NRDGL = 1, RLAM = 0: λ = 0.25 comes from the PERTDATA block of `ch4_water_fep_l025`'s final configuration, which is the start configuration here; the final configuration's PERTDATA is compared (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_lambdas | 12 | two energy groups, `LAMBDAS` on: bond/angle/dihedral/improper/mass polynomials on (1,1) and (2,2), vdw/crf (+soft) on every pair; the PERTATOMPAIR's LJ and CRF λ differ, as gromosXX allows with a critical message (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_nlam2 | 12 | aladip vacuum FEP, two energy groups, NLAM = 2: bonded terms and masses linear in λ (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_lambdas_nlam2 | 12 | NLAM = 2 with the `LAMBDAS` rows of `aladip_vacuum_fep_lambdas`: the per-kind `(1−L)ⁿ⁻¹·L′` factors with L_lj ≠ L_crf (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_shake | 12 | NTC = 3: the perturbed bonds become perturbed SHAKE constraints at the λ-mixed length; their dH/dλ in the constraints column (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_shake_lambdas | 12 | the same with bond individual λ per group; group 2's perturbed constraint keeps its length (same b0 in A and B), so only group 1's row acts (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_shake_hblock | 12 | `../shared/aladip_hblock.pttopo`: PERTBONDSTRETCHH 4-5 and PERTBONDSTRETCH 2-4 share atom 4 (gromosXX shakes the H block first), and the group-2 bond 8-10 changes length 0.133 → 0.153 nm under its own bond row (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_distres | 12 | perturbed distance restraints (NTDIR = 2, `../shared/check/aladip.distrest`): energy and dH/dλ per group (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_distres_lambdas | 12 | the same with `distanceres` individual λ on (1,1) and (2,2) (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_distres_tavg | 12 | NTDIR = −2, TAUDIR 0.1: the r⁻³ running average in energy, force and dH/dλ of perturbed and unperturbed restraints, and the DISRESEXPAVE block of the final configuration (1.8) | **PASS** |
| 4   | aladip_solvated_em_posres_groups | 72 | SD EM + position restraints with three energy groups (1-4, 5-12, water): restraint energy booked per group of the restrained atom (1.8) | **PASS** |
| 4   | aladip_harmonic_covalent_groups | 72 | `COVALENTFORM` harmonic bonds/angles with three energy groups (1-4, 5-12, water): bonded energies booked by first atom across solute groups (1.8) | **PASS** |
| 4   | aladip_multibath_swapped | 72 | `DOFSET 12 2 2  72 1 1`, baths at 300 K and 250 K: atoms map to baths by DOFSET range, not by bath number (gromosXX `in_bath`), dynamics included (1.8) | **PASS** |

**45 of 45 tests pass** (2026-08-30; nothing ignored). Plus gromosXX's own regression oracle,
ported: `tests/gromos_check_suite.rs` evaluates the `check/{aladip,aladip_unperturbed,
aladip_special}.t.cc` per-term energies (bonded, perturbed, soft-core, both reaction-field
variants, atomic cutoff, distance restraints, λ-derivatives) — 7 pass against the *current*
gromosXX binary to 1e-7 and against the suite's hard-coded table with its own δ; the dihedral
restraint value is `#[ignore]`d until 1.6 lands. It caught NSLFEXCL being parsed and never
consumed (fixed 0.0.35).

`gromos-tools` (2026-08-30): nine reference tests against gromos++ output — `make_top` (2), `com_top`,
`red_top`, `make_pt_top`, `sim_box`, `copy_box`, `inbox`, `ion` (`crates/gromos-tools/tests/`, the exact
gromos++ command in each fixture README). `gromos-md` (2026-08-30): 48 reference systems — the 45 physics levels plus three for real-world
input styles (`README-livecoms.md`), each also compared on its final configuration.
`gromos-analysis` (2026-08-30): fifteen reference tests against gromos++ output — `tser`, `mdf`, `dg_ener`, `dfmult`, `matrix_overlap`, `bilayer_dist`, `bilayer_oparam`, `jval`, `edyn`, `gca`, `eds_update_1`, `eds_update_2`, `jepot`, `pocket`, `dfgrid` (`crates/gromos-analysis/tests/data/README.md`); gromos-tools adds `explode`, `duplicate`, `pt_top`, `make_sasa_top`.

---
