# Porting guidelines — faithful to GROMOS, without its bugs

How to translate gromosXX/gromos++ into gromos-rs so the result is *both* faithful (matches GROMOS
where it should) and *correct* (doesn't inherit 30 years of accreted C++ defects). These are
guidelines, not tasks — the live work is on the project board, and the A/B/C classification below is
the discipline itself, not a summary of it.

## The bit-for-bit invariant

The gromosXX reference suite (`crates/gromos-md/tests/gromosXX_references/`) is the oracle. For any
**wired, tested** path:

- Reproducing GROMOS output to full precision is the definition of correct. Tolerances are tight
  (force ~1e-6, energy ~1e-8) by design.
- A change that alters that output is a **regression until proven otherwise** — never a silent
  improvement. If gromosXX is genuinely wrong and you want the corrected behavior, that is a
  *conscious, documented decision* (and a candidate for a `--gromos-compat` vs `--corrected` split),
  not an undiscussed divergence.

## The oracle only covers wired paths

This is the core risk. The reference suite catches divergence **only** for features that are wired
and tested. Every *unwired* path (e.g. triclinic boundaries, grid pairlists, Nosé-Hoover, GaMD,
lattice-sum/PME) is un-audited surface where two failure modes hide:

- **Inherited C++ bugs** — you faithfully port a defect.
- **Silent divergence** — your "cleaner" reimplementation disagrees with GROMOS and nothing notices
  until the feature is wired.

(A live example already in the tree: the triclinic nearest-image — see 1.4.)

FEP (fixed λ, `LAMBDAS`, `PRECALCLAM`, slow growth), EDS and A-EDS left this list in 1.8: they run
through the recipe and are pinned by native and patched-oracle references.

So the standing rule: **write the gromosXX reference test *before* wiring an unwired path, not
after.** The test is what converts un-audited surface into oracle-covered surface.

## Before porting a subsystem

1. **Grep the C++ for the authors' own flags:**
   `grep -rniE 'bug|fixme|wrong|hack' interaction/ algorithm/ math/`
   The gromosXX developers annotate their own defects and uncertainties. Two minutes here catches the
   issues they already knew about.
2. **Classify every suspicious thing you find** into one of three buckets, and treat each differently:
   - **(A) Genuine bug** → reproduce it first, pinned by a native reference and flagged (CHANGELOG
     "Reproduced gromosXX defects", an entry below, a code comment with file:line), so the divergence
     is demonstrated bit-for-bit; correct it in a later, separate commit that keeps the reference as
     documentation of the bug and gates the correction some other way. Where the defect is
     mathematical, state it as a theorem in the `math` repository (`OAI.Physics.Gromos.Divergences`).
     Never port one silently, and never skip straight to the correction.
   - **(B) Intentional GROMOS quirk** → reproduce it *exactly* and **document it as a deliberate
     decision** (in the style of CHANGELOG's "BONDANGLETYPE intentionally unsupported" entries), so it
     is never later mistaken for a porting mistake.
   - **(C) Latent divergence** (our code already differs, nothing tests it yet) → write the reference
     test, then reconcile.

## Second-source the subtle physics

For anything the authors flagged as uncertain, or anything intrinsically delicate (reaction-field
self-terms, virial, Ewald/lattice-sum, perturbed/soft-core terms): **derive it independently from
the GROMOS book / primary paper, implement to that, then diff against the C++.** Where they disagree,
investigate — do not assume the C++ is right. Transcribing uncertain C++ line-by-line is how you
inherit a bug you can't see.

Free leverage: gromosXX's own `check/*.t.cc` files hard-code per-term reference energies (including
perturbed/soft-core variants). They are a second source of truth independent of running the md
binary — use them to validate the very terms you were told not to trust.

## Reproduce-or-correct, made explicit

When GROMOS behavior and "physically correct" behavior diverge, you may eventually want both: a
faithful compatibility mode and a corrected mode. Design so the two can coexist behind a flag rather
than forcing an irreversible choice. Until then, **faithful is the default** and any deviation is
named and documented.

### The list

Each entry: what gromosXX does, which bucket, what gromos-rs does.

- **PRECALCLAM NRLAM = 1** (A). gromosXX accepts it and divides by `NRLAM − 1`: every λ_k is
  NaN. gromos-rs refuses it at build time (`check_precalclam`, `gromos-run/src/plan.rs`).
- **PRECALCLAM long-range ext energies** (B, D34). gromosXX computes the long-range perturbed
  pairs only on pairlist steps and keeps them in `m_longrange_storage` until the next one; their
  ext (λ_k) energies are filled only when that pairlist step is also an NTWG step, and are zero
  on the output steps after one that is not (`perturbed_nonbonded_set.cc:87-96, 280-300`).
  Reproduced (`Forcefield::longrange_perturbed`), pinned by `ch4_water_fep_precalc_twinrange`.
  BAR/MBAR from those rows is biased unless NTWG is a multiple of NSNB (2.058 vs 2.550 kJ/mol
  ⟨∂H/∂λ⟩ at λ = 0, NSNB 5, NTWG 2); `PrepareNote::PrecalclamLongRangeGaps` warns.
- **PRECALCLAM grid in μ** (convention of the published method: de Ruiter & Oostenbrink, JCTC
  2016, 12, 4476, p. 7 and Eqs. 17/19; not a quirk). The grid is in the individual μ_k (no NLAM
  exponent, no `LAMBDAS` polynomial, no L′; softness = μ_k), so `ext_ti_ana` can apply any
  parameterisation afterwards; the kinetic term keeps the simulated mass μ (`slam`) because
  momenta are conserved (`temperature_calculation.cc:200-224` ↔
  `gromos-integrators/src/algorithms/temperature.rs:220-250`). Reproduced.
- **PRECALCLAM current-λ leaks** (B). Inside the μ grid: the improper's ζ wrapped around the
  current-λ q0 and not re-wrapped, perturbed SHAKE adding to the `AB_bond` derivatives only, the
  1-4 ext call taking no Coulomb scaling, and the RF self ext term carrying no self-term
  correction. All reproduced.
- **PRECALCLAM derivative signs** (convention of the published method: JCTC 2016 Eqs. 9/12; not a
  quirk). The derivative record holds `A_de = +∂A/∂λ` but `B_de = −∂B/∂λ` for the nonbonded A/B
  energies, the B-state softness derivative being taken in `1 − μ`
  (`perturbed_nonbonded_term.cc:614-716` ↔ `lj_crf_soft_interaction_ext`,
  `gromos-forces/src/nonbonded/perturbed.rs:348-357`). Reproduced.
- **PRECALCLAM bond mixing** (B). The grid mixes `(1 − λ_k)·b0_A + λ_k·b0_B`
  (`perturbed_quartic_bond_interaction.cc:173-196` and the other bond forms), the main term
  `r0_A + L·Δr0`. Both forms kept, so the main term stays bit-identical to earlier versions and
  the grid to gromosXX.
- **Reaction field at EPSRF = 0** (C, reconciled). gromosXX sets `crf = −1` (conducting
  boundary; `nonbonded_term.cc:58-71`, `perturbed_nonbonded_term.cc:51-61`); gromos-rs used the
  ε formula. Now `crf = −1` in every kernel (`CRFParameters::new`), pinned by
  `aladip_solvated_epsrf0`.
- **Reaction field from RCRF** (C, reconciled). gromosXX builds every RF constant from RCRF
  (`nonbonded.rf_cutoff`, `in_parameter.cc:4606`) and accepts RCRF ≠ RCUTL silently
  (`check_parameter.cc` has no rule); gromos-rs built them from RCUTL. Now from RCRF, silently,
  pinned by `aladip_solvated_rcrf`.
- **PERTATOMPAIR argument shift** (A; corrected in 0.0.56, the reference pins the defect under
  `gromosxx-shifted`; patched oracle `tests/gromosXX_references/patches/D1_pair_kernel_args.patch`).
  `perturbed_nonbonded_pair.cc:469, 677, 992, 1200` call `eds_pert_lj_crf_interaction` with the
  outputs shifted by one slot against `perturbed_nonbonded_term.h:65-72`: for a pair with a
  perturbed atom and a normal or 1-4 state type the energies receive the kernel's dH/dλ, the
  force is `(f + e_lj + e_crf)·r`, dH/dλ is 0, and `eps = 0`, `coulomb_scaling = 0` drop the 1/r
  Coulomb term. gromos-rs computes the intended kernel by default (outputs in their slots,
  Coulomb scaling 1, the 1-4 scaling for 1-4 states), pinned by
  `aladip_vacuum_fep_pairpert_precalc_fixed` (native 247c95c4 + D1 patch); the recipe's
  `perturbation.atom_pair_kernel = "gromosxx-shifted"` reproduces the defect
  (`atom_pair_state`, `gromos-forces/src/nonbonded/perturbed.rs`), pinned by
  `aladip_vacuum_fep_pairpert_precalc`. Lean (`math`, `OAI.Physics.Gromos.Divergences`):
  `shiftedPairForce_eq_pairForce_iff`, `shifted_ti_shortfall`, `crfSoftEnergy_shift`.
- **EDS multi-s / pair-s weights** (A; corrected without a reproducing mode, decided in 1.8:
  the released branch's forces are not the gradient of any energy it reports, and the patch
  reproduces it). gromosXX 247c95c4 seeds the per-state log-prefactors with 0 instead of −∞
  (`eds.cc:518, 620`), drops the `− max` of the first-pair update (`:550-551, 647-648`: a + ln 2
  for a ≥ 0), multi-s never visits the pairs (0, j ≥ 2) (`:557`, so its V_R is wrong for N ≥ 3)
  and neither form books Σ p_i ∂V_i/∂λ (`:603-604, 704-705`). gromos-rs implements the documented
  formulas (`parameter.h:624-658`; `reference_hamiltonian`, `gromos-integrators/src/eds.rs`),
  what native + `patches/D2_eds_weights.patch` computes, pinned by
  `aladip_vacuum_eds_multis_fixed`/`_pairs_fixed`; the released binary pins only the end-state
  energies (and FORM 3's V_R) at frame 0 (`aladip_vacuum_eds_multis`/`_pairs`), and
  `gromos-run/tests/eds_forces.rs` checks F = −∇V_R for all three forms. FORM 1 is gromosXX's,
  bit for bit, including its running pairwise log-sum-exp in state order (`eds.cc:409-424`).
- **EDS end-state bookkeeping** (B), reproduced (`gromos-forces/src/nonbonded/eds.rs`,
  `Topology::apply_eds_perturbation`): `MPERTATOM` moves an EDS atom's excluded and 1-4 pairs off
  the standard lists, the first-listed atom owning an EDS–EDS pair, but they stay off the
  pairlist — `all_exclusion` was built before (`in_perturbation.cc:1757-1777`,
  `in_topology.cc:2111`); the standard RF self term books the topology charge and the EDS loop
  takes it back out of `crf_energy[g][g]` (`perturbed_nonbonded_innerloop.cc:3264-3293`); the
  end-state virials of the 1-4 and RF-excluded terms are kept for every virial type, the pairs'
  only with a virial (`perturbed_nonbonded_set.cc:306-318, 367-376`); the long-range EDS pairs
  of a twin-range run are computed on pairlist steps only. Pinned by the four FORM 1 EDS
  references.
- **EDS inputs gromosXX does not check** (refused): FORM 3 pair states outside 1..N
  (`eds.cc:635, 660` index with them), an `MPERTATOM` IAC of 0 (unsigned, wraps;
  `in_perturbation.cc:1721-1727`). Not ported, refused by name: the shifted reaction field
  (NLRELE 4), REPLICA_EDS; SD and the coarse-grained and virtual-atom variants need blocks
  the recipe does not model.
- **A-EDS inputs gromosXX crashes on or computes nonsense with** (A, refused, not reproduced;
  `check_aeds`, `gromos-run/src/plan.rs`, and `emax_emin_update`): EMIN > EMAX (warning, then
  EIR unread and a segfault, `in_parameter.cc:4084-4088`; D14); BSTEPS < ASTEPS in FORMS 2/4/5/6
  (unsigned `bsteps − asteps`, τ ≈ 4·10⁹ steps, offsets frozen, `eds.cc:186`; D13 — the pure
  `memory_time` reproduces the wrap); BMAXTYPE 2 with a target state without variance, unless
  the RESTREMIN clamp catches it (division by b_max = 0 at `eds.cc:314`, E_min = −∞, V_R* NaN;
  D15). Pinned by unit tests and `gromos-md/tests/aeds_runs.rs`.
- **A-EDS defects corrected by default** (A; corrected in 0.0.63, each reproduced by a recipe
  switch; patched oracle `tests/gromosXX_references/patches/D16_*`, `D25_*`, `D32_*`;
  `gromos-integrators/src/eds.rs`):
  - the energy columns (`eds_emax`, `eds_emin`, the offsets) written after the search update, so
    a frame's `eds_vmix`/`eds_vr` belong to the previous frame's printed parameters
    (`eds.cc:368-372`; D16). Default: the parameters the frame's V_R* used;
    `aeds.energy_columns = "gromosxx-post-update"` reproduces.
  - the FORMS 5/6 prevalences as plain exponentials (`:339-343`, NaN and frozen offsets once
    −β(V_i − E_i^R) > 709.78; D25). Default: log-sum-exp; `aeds.prevalence =
    "gromosxx-plain-exp"` reproduces (`aladip_vacuum_aeds_f6_ovf`).
  - a never-visited state 1 as the default E_min target (`:252-261`; D32, its mean 0 a fake G,
    σ 0 the D15 division). Default: the lowest mean among the states visited more than once,
    state 1 only when none is; `aeds.target_state = "gromosxx-state1-default"` reproduces.
  The native search references set the reproducing switches; the `aladip_vacuum_aeds_*_fixed`
  references pin the defaults against the patched binary; `gromos-run/tests/aeds_replay.rs`
  replays both. Lean statements: `~/.cache/gromos-rs-item-1.8/lean_aeds_statements.md` (not yet
  in `math`).
- **A-EDS bookkeeping** (B, reproduced behaviour documented in GROMOS Vol. 2, not defects;
  `gromos-integrators/src/eds.rs`): the FORMS 5/6 offset step is ASTEPS·Δt (`:355`; D33, not the
  dossier's D17, which is the extended-TI 1-4 scaling) — a
  time if ASTEPS is a step count, an energy if it is the rate Vol. 4 and the tutorial describe;
  τ ramps over the current run only (`sim.steps()` from 0), as Vol. 2 documents ("t_tot is the
  total simulation time of the run") while Perthold 2020 ramps over the whole search, for which
  the input carries no total length; the mean of V_i − E_i^R averages over the offsets of every
  step (`:246`), Vol. 2's \overline{V_i − E_i^R} with E_i^R moving; the boost scales only the
  EDS part of forces, virial and dV/dλ (`:137-176`). Pinned by the `aladip_vacuum_aeds_*`
  references (whole block, 9-decimal AEDSSEARCH) and `gromos-run/tests/aeds_replay.rs`
  (gromosXX's V_i(t) through the pure functions, ≤ 1e-9 relative).
- **TI-EDS inputs gromosXX computes wrongly** (A, refused, not reproduced; `check_ti_eds`,
  `gromos-run/src/prepare.rs`): an atom in both `PERTATOMPARAM` and `MPERTATOM` gets pair code
  6, which no case handles, so its pairs are dropped silently (`perturbed_nonbonded_innerloop.cc:54-65,
  641-643`; dossier D18); a `PERTATOMPAIR` pair with an EDS atom is computed with the topology
  type and charge in every end state (`perturbed_nonbonded_pair.cc` has no EDS branch; D19).
  PRECALCLAM and slow growth (DLAMT ≠ 0) stay refused with EDS (gromosXX locks both against EDS,
  `check_parameter.cc:4271, 2528-2530`), NSCALE ≠ 0 with or without it.
- **TI-EDS pair ownership** (deliberate difference, same terms): gromosXX's perturbed atom owns
  an excluded or 1-4 EDS–perturbed pair (`PERTATOMPARAM` takes the pair before `MPERTATOM` looks,
  `in_perturbation.cc:1363-1391, 1748-1777`; case 4 at `perturbed_nonbonded_innerloop.cc:2277-2341,
  3696-3743`), gromos-rs's EDS atom owns every EDS–perturbed pair (`Topology::apply_eds_perturbation`,
  `EdsAtoms`, `gromos-forces/src/nonbonded/eds.rs`): the same kernel per state, summed in another
  order, so agreement is to rounding, not bit for bit (pinned at energy_rel 1e-10 /
  force_abs 1e-8 by `aladip_eds_solv_ti`, `_groups`, `_multis_fixed`; every printed number
  equal). FORMS 2 and 3 book Σ p_i ∂V_i/∂λ as FORM 1 does (the released gromosXX books 0, D2).
- **Slow growth** (DLAMT ≠ 0; reproduced, `gromos_run::advance_lambda`, called after every step
  by `md` and `Simulation`). gromosXX adds DLAMT·Δt to λ at the end of each step and re-evaluates
  the individual λ and the perturbed masses (`slow_growth.cc:43-53`, pushed after the barostat,
  `create_md_sequence.cc:279-286`); not on the step the minimiser stops at. Its bookkeeping (B):
  the `.trg` prints the λ the step ran at (`topo.old_lambda()`, `out_configuration.cc:1726`), the
  `.tre` mass the masses at the next λ (`:1786`, written after `Slow_Growth`), PERTDATA the λ
  after the last step; the repeated converged `.trg` frame of a minimisation prints the previous
  step's λ next to derivatives taken at the current one. λ has no upper bound (A): RLAM is
  checked to [0, 1] (`in_parameter.cc:878`), the running λ and a PERTDATA λ are not, so a run
  extrapolates past the end state with no warning. Reproduced (`aladip_vacuum_fep_slowgrowth_cont`
  runs to λ = 1.125), with a warning gromosXX does not give: `md` and `Simulation` report the
  step that takes λ out of [0, 1] (`advance_lambda`'s return) and a run that starts outside it
  (`PrepareNote::LambdaOutOfRange`); never a refusal. gromosXX's own locks are refused (`check_slow_growth`, `plan.rs`): LAMBDAS
  with NTIL on, PRECALCLAM, EDS, A-EDS, lattice sums, NLRELE 4, conjugate gradients, GAMD,
  NTPOR 3. Not ported (A): the `.omd` "CUMULATIVE DG" sums DLAMT·∂H/∂λ without the Δt
  (`average.cc:425, 663`), a ΔG 1/Δt too large (−615.5 instead of −1.231 kJ/mol over 20 steps of
  `aladip_vacuum_fep_slowgrowth`'s run); gromos-rs writes no `.omd` averages.
- **`dGslv_pbsolv` constants** (checked in 1.8 against gromos-core; no difference that changes
  a number, nothing to reconcile). 1/(4πε₀) is the topology's `PHYSICALCONSTANTS` value, as in
  `md` (138.9354 in every reference topology, `gromos_core::units::four_pi_eps_i`); only PQR
  input falls back to gromos++'s `gmath/Physics.cc:53` 138.935484611
  (`pb::GROMOSPP_FOUR_PI_EPS_I`, relative 6·10⁻⁷ from gromos-core's, kept for byte identity).
  k_B and N_A are not used: the NPBC boundary has κ = 0 (`PB_Parameters.cc:69`, no ionic
  strength, so no kT) and the energies come out in kJ/mol through 1/(4πε₀); gromos++'s own
  pair (`Physics.cc:51-52`, k_B = 1.380658·10⁻²³ × 6.0221367·10²³ / 1000 = 0.0083145110 vs
  gromos-core's gromosXX 0.00831441) never reaches the program. The water defaults are
  gromos++'s, not gromos-core's (which has none): `@epsNPBC` 78.4 (and the FFT solver's
  default ε 78.4), `@probeRAD` 0.14 nm, the hydrogen radius 0.05 nm, fixed because gromos++
  reads `@radH` but does not list it (`dGslv_pbsolv.cc:448, 454, 498`; `PB_Parameters.cc`
  `hradius`); `@epsSOLV` has no default and should be the
  solvent model's permittivity (66.6 for SPC in the gromos++ examples), `@epsRF`/`@rcut` the
  run's EPSRF/RCRF. Radii come from the topology's LJ parameters with the probe's IAC, the
  engine's own table. Pinned by `dgslv_pbsolv_unit_charges_in_spc_match_gromospp`.
- **Constraint virial under the molecular virial** (C, ours, reconciled). gromosXX adds the
  SHAKE/SETTLE virial only for the atomic virial (`shake.h:380`, `settle.cc:397`); gromos-rs added
  it under every virial type, so the barostat saw −1500 instead of +28 kJ mol⁻¹ nm⁻³. Now gated
  on `atomic_virial` (`shake_buffered`, `settle` in `gromos-integrators/src/constraints.rs`),
  pinned by `aladip_solvated_npt_shake`. Lean: `intragroup_pair_virial_zero`,
  `constraint_virial_pressure_shift`.
- **ROTTRANS degrees of freedom** (C, ours, reconciled). gromosXX removes 3 COM and 3 internal
  dof from atom 1's baths (`topology.cc:1370-1400`); gromos-rs did not. Now in `bath_dofs`
  (`gromos-run/src/dof.rs`), pinned by `aladip_multibath_rottrans`.
- **Per-bath temperature rows** (C, ours, reconciled). gromosXX's `VOLUMEPRESSURE03` prints
  `T_total`/`T_com`/`T_ir`/`scale` per bath (`out_configuration.cc:3270-3346`); gromos-rs printed
  the system temperature on every row. Now per bath (`volume_pressure` in
  `gromos-md/src/bin/md.rs`, written by `gromos-io/src/energy.rs`), compared whole-block by the
  reference test.
- **Virial booking quirks** (B), reproduced in the force-field `apply`
  (`gromos-integrators/src/algorithms/forcefield.rs`), `calculate_bonded_forces`
  (`gromos-forces/src/bonded/mod.rs`) and `settle` (`gromos-integrators/src/constraints.rs`),
  pinned by `VOLUMEPRESSURE03`:
  - bonded and perturbed bonded virial added under every virial type, the gate commented out
    (`quartic_bond_interaction.cc:93`, `perturbed_quartic_bond_interaction.cc:125`);
  - perturbed 1-4, RF-excluded and atom-pair virials unconditional
    (`perturbed_nonbonded_innerloop.cc:2247-2250, 3941-3946`, `perturbed_nonbonded_pair.cc:1559-1561`),
    plain pairs gated (`nonbonded_set.cc:428-430`);
  - SETTLE gated on the atomic virial with no vacuum split (`settle.cc:397`); LINCS adds no
    constraint virial.
- **Thermostat output quirks** (B), reproduced: `NDFMIN` not subtracted from the COM/IR dof
  (`multibath.cc:139-149`; `bath_dofs`, `gromos-run/src/dof.rs`); a bath's `scale` printed 0
  without a thermostat (`multibath.h:48`, `berendsen_thermostat.cc:134`); the kinetic-tensor rows
  printed unrotated while the trace uses the rotated tensor (`out_configuration.cc:3331-3344`);
  the last two in `volume_pressure` (`gromos-md/src/bin/md.rs`).
- **Not ported, no reference**: position-restraint and distance/dihedral-restraint virials
  (`position_restraint_interaction.cc:89`, `distance_restraint_interaction.cc:194`); baths with
  dof ≤ 0 losing their coupling (`topology.cc:1405-1411`).

## In one line

Tests before wiring; derivation before transcription; grep before porting; document every deliberate
quirk. Faithfulness and correctness are both required — the discipline above is how you get both
instead of trading one for the other.
