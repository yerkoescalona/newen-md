# gromosXX reference systems

One directory per system: the inputs, and `expected/` holding what the working upstream
(`.local/gromosXX`, md++ 2023-04-15) produced for them. `test_gromosXX_references.rs` runs each and
compares energies, trajectory and final configuration; the Python suite drives the same systems
through `Simulation` so both front-ends are held to one oracle.

**A system counts as passing when its numbers agree with md++ — step-0 energies to the printed
precision and the trajectory within the stated tolerances — not when the binary exits 0.**

Adding one: see `crates/gromos-md/.claude/CONTEXT.md` §"How to add a reference test".
An optional `compare` key in `input.toml` (`frame0_eds_states`, `frame0_eds_states_vr`) restricts
the comparison to frame 0's EDS end-state energies (and V_R), for references whose forces come
from a defective gromosXX branch.
`README-livecoms.md` covers the three systems that came from real LiveCoMS tutorial input.
An optional `[recipe_overrides]` table in `input.toml` sets recipe fields no `.imd` block carries
(`"perturbation.atom_pair_kernel" = "gromosxx-shifted"`, the A-EDS `"aeds.energy_columns"`,
`"aeds.prevalence"`, `"aeds.target_state"`); both runners apply it (`md @recipe_set key=value`,
`Recipe.update`).

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
| 4   | ch4_water_fep_precalc | 2998 | `PRECALCLAM 11 0 1` (extended TI): the A/B LJ energies and derivatives at eleven λ points in the `# precalclam` rows of the `.tre` and `.trg`, AB_kinetic of the perturbed atom (1.8) | **PASS** |
| 4   | meoh_water_fep_precalc | 2862 | the same for charged perturbed atoms: the CRF pair, RF self and RF exclusion ext terms (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_precalc | 12 | AB_bond/AB_angle/AB_improper (regular and soft-core), AB_kinetic, the ABdih row and the PERTATOMPAIR's constant A/B energies at every λ point (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_pairpert_precalc | 12 | `../shared/aladip_pairpert.pttopo`: PERTATOMPAIR rows 4-10 (normal/excluded), 1-11 (excluded/1-4), 4-7 (1-4/normal) with perturbed atoms, next to the unperturbed 7-10: main energies and the ext rows of the filled state; run with `atom_pair_kernel = "gromosxx-shifted"`, it pins gromosXX's argument-shift defect (CHANGELOG 0.0.55) (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_pairpert_precalc_fixed | 12 | the same run, expected output from native 247c95c4 + D1 patch (`tests/gromosXX_references/patches/D1_pair_kernel_args.patch`: the pair kernel's outputs in their own slots, Coulomb scaling 1, the 1-4 scaling for 1-4 states): the intended kernel, gromos-rs's default (1.8, 0.0.56) | **PASS** |
| 4   | aladip_vacuum_fep_lambdas_precalc | 12 | `LAMBDAS` on and `PRECALCLAM 6 0.2 0.7`: the grid ignores the individual-λ polynomials, except the current mass λ in AB_kinetic (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_shake_precalc | 12 | perturbed SHAKE (NTC = 3) with PRECALCLAM: its contribution to the AB_bond derivatives only (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_distres_precalc | 12 | perturbed distance restraints with PRECALCLAM: AB_disres energies and derivatives (1.8) | **PASS** |
| 4   | ch4_water_fep_precalc_twinrange | 2998 | twin range, pairlist every 5 steps, NTWE = NTWG = 2: the long-range perturbed pairs (main energies and ext rows) are computed on pairlist steps and reused, the ext rows only when that step is also an output step (gromosXX quirk) (1.8) | **PASS** |
| 4   | nacl_water_fep / _precalc / _lambdas | 62 | a charge change: Na⁺ → neutral dummy (`../shared/nacl_na_dummy.ptp`, charge and IAC, α 0.5/0.5, ΔQ = −1) in `nacl_water_box` (20 SPC, NLRELE 1, RCRF 0.9, EPSRF 61), RLAM 0.5, 20 steps, `.trg` every step: plain TI, `PRECALCLAM 11 0 1`, and `LAMBDAS` with crf ≠ vdw rows on (1,1) and (1,2). They pin the perturbed RF self term and the RF excluded Na–Cl pair, which move only when the charge does; both are booked in the (1,1) crf cells of the `# nonbonded` rows (dH/dλ −89.58 kJ/mol at frame 0) and in the `A_e_crf`/`B_e_crf` `# precalclam` columns, not in columns of their own (1.8) | **PASS** |
| 4   | nacl_water_fep_14 | 62 | `nacl_water_fep_precalc` with the Na–Cl pair a 1-4 pair instead of an exclusion (`../shared/nacl_water_box_14.topo`): the perturbed 1-4 term (CS6/CS12) of a charge change, its dH/dλ and ext rows (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_slowgrowth | 12 | slow growth, RLAM 0.125, DLAMT 2/ps, 100 steps (λ 0.125 → 0.525), `.trg` every step: λ += DLAMT·Δt after every step, the λ-mixed masses (kinetic dH/dλ, leap-frog) and every perturbed bonded and nonbonded term following it; the `.trg` prints the λ the step ran at, the `.tre` mass the masses at the next λ, the final PERTDATA 0.525 (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_slowgrowth_shake | 12 | the same with NTC = 3 and two energy groups: the perturbed constraint lengths re-evaluated at the moving λ (1.8) | **PASS** |
| 4   | aladip_vacuum_fep_slowgrowth_cont | 12 | NRDGL = 1 from `aladip_vacuum_fep_slowgrowth`'s native final configuration (λ 0.525 from PERTDATA, RLAM 0 ignored), DLAMT 6/ps, 50 steps: λ runs past 1 to 1.125 unbounded, as in gromosXX (1.8) | **PASS** |
| 4   | ch4_water_fep_slowgrowth | 2998 | `ch4_water_fep_l025` with DLAMT 5/ps, 50 steps (λ 0.25 → 0.75), energies and `.trg` every step, positions and forces every 25: soft-core pairs at a moving λ in a periodic box with solvent SHAKE (1.8) | **PASS** |
| 4   | aladip_vacuum_em_slowgrowth | 12 | `aladip_vacuum_em` with `aladip.pttopo`, RLAM 0.125, DLAMT 2/ps: λ moves on every minimisation step but the one the minimiser stops at (final PERTDATA 0.173); the repeated converged `.trg` frame carries the previous step's λ (1.8) | **PASS** |
| 4   | aladip_solvated_epsrf0 | 72 | EPSRF = 0: a conducting boundary, C_RF = −1 (1.8) | **PASS** |
| 4   | aladip_solvated_rcrf | 72 | RCRF = 1.2 with RCUTL = 1.4: the reaction-field constants come from RCRF (1.8) | **PASS** |
| 4   | aladip_solvated_em_posres_groups | 72 | SD EM + position restraints with three energy groups (1-4, 5-12, water): restraint energy booked per group of the restrained atom (1.8) | **PASS** |
| 4   | aladip_harmonic_covalent_groups | 72 | `COVALENTFORM` harmonic bonds/angles with three energy groups (1-4, 5-12, water): bonded energies booked by first atom across solute groups (1.8) | **PASS** |
| 4   | aladip_multibath_swapped | 72 | `DOFSET 12 2 2  72 1 1`, baths at 300 K and 250 K: atoms map to baths by DOFSET range, not by bath number (gromosXX `in_bath`), dynamics included (1.8) | **PASS** |
| 4   | aladip_solvated_npt_shake | 72 | Berendsen barostat (VIRIAL 2) with NTC 3 and solvent SHAKE: the constraint virial enters only the atomic virial; 20 steps, every `VOLUMEPRESSURE03` compared (1.8) | **PASS** |
| 4   | aladip_multibath_rottrans | 72 | `aladip_multibath` + `ROTTRANS 1 12`: 3 COM and 3 IR degrees of freedom leave the solute bath, the system is centred on the solute (1.8) | **PASS** |
| 4   | aladip_vacuum_eds | 12 | EDS FORM 1 (s 0.5) over two end states of CB, N, H (`../shared/aladip_eds2.mpt`): per-state LJ+RF pairs, 1-4 and RF excluded/self terms, V_R (slot 35, special energy), the mixed forces and the `# eds` rows (1.8) | **PASS** |
| 4   | aladip_vacuum_eds3 | 12 | EDS FORM 1 over three end states (`../shared/aladip_eds3.mpt`): gromosXX's running pairwise log-sum-exp in state order (1.8) | **PASS** |
| 4   | aladip_eds_solv | 72 | EDS FORM 1 in 20 SPC, two baths, NTC 3, twin range with NSNB 5: EDS-solvent pairs short and long range (long-range reused between pairlist steps), an EDS-EDS excluded pair (1,4) and 1-4 pair (9,11), the RF self term taken out with the topology charge (H11's state charges differ from it) (1.8) | **PASS** |
| 4   | two_waters_eds | 3000 | LiveCoMS t_05 inputs (`2h2o.top`, `2h2o.cnf`, `two_waters.ptp`, softness columns ignored): EDS FORM 1 over four end states with dummy types across two solute waters, three energy groups (1.8) | **PASS** |
| 4   | aladip_eds_solv_ti | 72 | TI-EDS: `aladip_eds_solv` with EDS atoms N4, H11 (two states) and perturbed atoms CB1, C8, N10, CB12 (`PERTATOMPARAM`, own α), RLAM 0.3, one energy group, `.trg` every step: the soft-core A/B pair of each state for the EDS–perturbed pairlist pairs (EDS atom first 4-12 and second 1-11), the 1-4 pair 4-10 and the excluded pairs 4-1, 4-8, 11-8, 11-10, 11-12 (soft RF), ∂V_i/∂λ in the `.trg`'s `# eds` rows and Σ p_i ∂V_i/∂λ in its slot 35. gromos-rs's EDS atom owns the excluded and 1-4 EDS–perturbed pairs (gromosXX's perturbed atom): the same terms in another summation order, so the stated tolerance is energy_rel 1e-10 / force_abs 1e-8; measured: every printed `.tre`/`.trg` number equal, forces within 1e-9 (1.8) | **PASS** |
| 4   | aladip_eds_solv_ti_groups | 72 | the same with two energy groups (1-5, 6-72: EDS and perturbed atoms on both sides) and `LAMBDAS` on for vdw/crf and their softness: each EDS–perturbed pair takes the λ of its group cell (1.8) | **PASS** |
| 4   | aladip_vacuum_eds_multis_fixed | 12 | EDS FORM 2 (s per state pair) over three end states, expected output from native 247c95c4 + D2 patch (`patches/D2_eds_weights.patch`): V_R over all N(N−1)/2 pairs and the weights ∂V_R/∂V_i, gromos-rs's formulas (1.8) | **PASS** |
| 4   | aladip_vacuum_eds_pairs_fixed | 12 | EDS FORM 3 (pairs (1,2) 0.3, (2,3) 0.1), native + D2 patch: the documented weights (1.8) | **PASS** |
| 4   | aladip_eds_solv_ti_multis_fixed | 72 | `aladip_eds_solv_ti` under FORM 2, native + D2 patch: Σ p_i ∂V_i/∂λ in slot 35 (the released gromosXX books 0 for FORMS 2/3, `eds.cc:603-604, 704-705`); with two states FORM 2 equals FORM 1 (1.8) | **PASS** |
| 4   | aladip_vacuum_eds_multis | 12 | the FORM 2 run on the released gromosXX: only frame 0's end-state energies are compared — its force weights are defective (zero-seeded log-sum-exp, missing `− max`, pairs (0, j≥2) skipped, which also makes its V_R wrong for N ≥ 3), so its forces and trajectory are not those of its V_R (`compare` in `input.toml`) (1.8) | **PASS** |
| 4   | aladip_vacuum_eds_pairs | 12 | the FORM 3 run on the released gromosXX: frame 0's end-state energies and V_R (right for any N); forces and trajectory not compared for the same reason (1.8) | **PASS** |
| 4   | aladip_vacuum_aeds_below / _inside / _above | 12 | A-EDS FORM 1 over the two end states of `../shared/aladip_eds2.mpt`, EIR 0 10, (EMAX, EMIN) = (50, −10) / (40, −60) / (−40, −80): V_mix below, inside and above the window — V_mix (slot 34), V_R* (35), E_max/E_min (36/37), the offset column, the EDS forces and virial scaled by the boost factor inside (1.8) | **PASS** |
| 4   | aladip_vacuum_aeds_f2 / _f4 / _f5 / _f6 | 12 | A-EDS searches over three gentle end states (`../shared/aladip_aeds3.mpt`, N4/H5 charges), EIR 0 0 −3, from Maxwell-Boltzmann velocities, 300 steps (FORM 5: 250): FORM 2 offsets (τ 20 → 50 steps), FORM 4 offsets + E_max/E_min (all four E_min branches, 16 round trips), FORM 5 and 6 prevalence-driven offsets; slots 34-39 and the offsets every frame, and the final AEDSSEARCH block; run with the overrides reproducing gromosXX's D16 (and D25/D32 where the FORM reaches them) (1.8) | **PASS** |
| 4   | aladip_vacuum_aeds_sigma0 | 12 | A-EDS FORM 3, BMAXTYPE 2, EIR −20 0: the never-visited default target state 1 has no variance (b_max = 0); RESTREMIN 1 holds E_min at its fake mean 0 (with RESTREMIN 0 gromosXX divides by zero; refused); run with the D16/D32 reproducing overrides (1.8) | **PASS** |
| 4   | aladip_vacuum_aeds_cont1 / _cont2 | 12 | an A-EDS FORM 4 search in two jobs of 100 steps: job 2 starts from job 1's native final configuration, NTIAEDSS 0, its AEDSSEARCH block overriding the input EMAX 50, EMIN −50, EIR 9 9 9; run with the D16/D32 reproducing overrides (1.8) | **PASS** |
| 4   | aladip_vacuum_aeds_f6_ovf | 12 | A-EDS FORM 6 with EIR 0 6000 0, 100 steps: e^{−β(V_2 − E_2^R)} overflows every step, the released gromosXX's prevalences are NaN and its offsets stay at the step-0 centring (−2000, 4000, −2000); run with the reproducing overrides (1.8 D25) | **PASS** |
| 4   | aladip_vacuum_aeds_f2_fixed / _f4_fixed / _f5_fixed / _f6_fixed / _f6_ovf_fixed / _sigma0_fixed / _cont1_fixed / _cont2_fixed | 12 | the same runs, expected output from native 247c95c4 + `patches/D16_aeds_energy_columns.patch`, `D25_aeds_prevalence_lse.patch`, `D32_aeds_target_state.patch`: E_max, E_min and the offsets of a frame are the ones its V_R* used, log-sum-exp prevalences (the offsets of `_f6_ovf_fixed` move), a sampled state as the E_min target (`_sigma0_fixed`: G and σ_G of state 2, not the never-visited state 1's 0); gromos-rs's defaults (1.8, 0.0.63) | **PASS** |

**119 of 120 registered tests pass** (2026-10-16; `aladip_gamd` ignored). Plus gromosXX's own regression oracle,
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
