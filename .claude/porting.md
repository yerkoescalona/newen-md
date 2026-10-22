# Porting guidelines — faithful to GROMOS, without its bugs

How to translate gromosXX/gromos++ into newen so the result is *both* faithful (matches GROMOS
where it should) and *correct* (doesn't inherit 30 years of accreted C++ defects). These are
guidelines, not tasks — the live work is on the project board, and the A/B/C classification below is
the discipline itself, not a summary of it.

## The bit-for-bit invariant

The gromosXX reference suite (`crates/newen-md/tests/gromosXX_references/`) is the oracle. For any
**wired, tested** path:

- Matching the oracle to full precision is the definition of correct. Tolerances are tight
  (force ~1e-6, energy ~1e-8) by design.
- A change that alters that output is a **regression until proven otherwise** — never a silent
  improvement.
- The oracle is gromosXX as released, except where a known defect has been corrected: there it is
  gromosXX with that defect's patch applied (`tests/gromosXX_references/patches/`). newen is gated
  on the patched binary and carries no mode that reproduces the defect (see "One behaviour per
  feature" below).

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
   - **(A) Genuine bug** → newen implements the intended behaviour, and only that. Write the
     smallest patch to gromosXX (or gromos++) that corrects the defect, check that the patched and
     released binaries differ only in what the defect touches, regenerate the reference with the
     patched binary (`<system>_fixed`, the patch named in its `input.toml`) and gate newen on it.
     Record the defect in the dossier in `~/Development/gromos-divergences` (the reproducer against
     released gromosXX, the patch, a Lean theorem where the defect is mathematical) and add an
     entry to the list below. Never port one silently, never correct one without a patched oracle,
     and never add an option that reproduces it.
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

## One behaviour per feature

newen has one behaviour for every feature: the intended one. Where gromosXX is right, that is
gromosXX's behaviour, its conventions and quirks included (bucket B). Where gromosXX is wrong, it is
the corrected behaviour, gated by a patched oracle (bucket A). There is no compatibility mode, no
`--gromos-compat` flag and no recipe value that selects a known defect. Each such switch doubles the
behaviours the code can have, most combinations go untested, and the next reader — person or agent —
can no longer tell which behaviour is intended, and copies the pattern.

The proof that a defect exists (for the gromosXX maintainers) lives in `gromos-divergences`, which
runs released gromosXX against the patched one. newen keeps only the patch and the `_fixed`
reference. A released-binary reference stays in newen only for what the defect does not touch.

Code comments explain the physics. A corrected defect is named by its dossier id (`D41`); the
gromosXX file:line map belongs in the entry below and in the dossier, not in the code.

Bucket-A entries below still marked "Reproduced" predate this rule: each is open work, to be
corrected behind a patched oracle.

### The list

Each entry: what gromosXX does, which bucket, what newen does.

- **PRECALCLAM NRLAM = 1** (A). gromosXX accepts it and divides by `NRLAM − 1`: every λ_k is
  NaN. newen refuses it at build time (`check_precalclam`, `newen-run/src/plan.rs`).
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
  `newen-integrators/src/algorithms/temperature.rs:220-250`). Reproduced.
- **PRECALCLAM current-λ leaks** (B). Inside the μ grid: the improper's ζ wrapped around the
  current-λ q0 and not re-wrapped, perturbed SHAKE adding to the `AB_bond` derivatives only, the
  1-4 ext call taking no Coulomb scaling, and the RF self ext term carrying no self-term
  correction. All reproduced.
- **PRECALCLAM derivative signs** (convention of the published method: JCTC 2016 Eqs. 9/12; not a
  quirk). The derivative record holds `A_de = +∂A/∂λ` but `B_de = −∂B/∂λ` for the nonbonded A/B
  energies, the B-state softness derivative being taken in `1 − μ`
  (`perturbed_nonbonded_term.cc:614-716` ↔ `lj_crf_soft_interaction_ext`,
  `newen-forces/src/nonbonded/perturbed.rs:348-357`). Reproduced.
- **PRECALCLAM bond mixing** (B). The grid mixes `(1 − λ_k)·b0_A + λ_k·b0_B`
  (`perturbed_quartic_bond_interaction.cc:173-196` and the other bond forms), the main term
  `r0_A + L·Δr0`. Both forms kept, so the main term stays bit-identical to earlier versions and
  the grid to gromosXX.
- **Reaction field at EPSRF = 0** (C, reconciled). gromosXX sets `crf = −1` (conducting
  boundary; `nonbonded_term.cc:58-71`, `perturbed_nonbonded_term.cc:51-61`); newen used the
  ε formula. Now `crf = −1` in every kernel (`CRFParameters::new`), pinned by
  `aladip_solvated_epsrf0`.
- **Reaction field from RCRF** (C, reconciled). gromosXX builds every RF constant from RCRF
  (`nonbonded.rf_cutoff`, `in_parameter.cc:4606`) and accepts RCRF ≠ RCUTL silently
  (`check_parameter.cc` has no rule); newen built them from RCUTL. Now from RCRF, silently,
  pinned by `aladip_solvated_rcrf`.
- **PERTATOMPAIR argument shift** (A, D1; corrected in 0.0.56; patched oracle
  `tests/gromosXX_references/patches/D1_pair_kernel_args.patch`).
  `perturbed_nonbonded_pair.cc:469, 677, 992, 1200` call `eds_pert_lj_crf_interaction` with the
  outputs shifted by one slot against `perturbed_nonbonded_term.h:65-72`: for a pair with a
  perturbed atom and a normal or 1-4 state type the energies receive the kernel's dH/dλ, the
  force is `(f + e_lj + e_crf)·r`, dH/dλ is 0, and `eps = 0`, `coulomb_scaling = 0` drop the 1/r
  Coulomb term. newen computes the intended kernel (outputs in their slots, Coulomb scaling 1,
  the 1-4 scaling for 1-4 states; `atom_pair_state`, `newen-forces/src/nonbonded/perturbed.rs`),
  pinned by `aladip_vacuum_fep_pairpert_precalc_fixed` (native 247c95c4 + D1 patch). The
  released-binary reproducer lives in `gromos-divergences`. Lean (`math`,
  `OAI.Physics.Gromos.Divergences`): `shiftedPairForce_eq_pairForce_iff`, `shifted_ti_shortfall`,
  `crfSoftEnergy_shift`.
- **EDS multi-s / pair-s weights** (A; corrected without a reproducing mode, decided in 1.8:
  the released branch's forces are not the gradient of any energy it reports, and the patch
  reproduces it). gromosXX 247c95c4 seeds the per-state log-prefactors with 0 instead of −∞
  (`eds.cc:518, 620`), drops the `− max` of the first-pair update (`:550-551, 647-648`: a + ln 2
  for a ≥ 0), multi-s never visits the pairs (0, j ≥ 2) (`:557`, so its V_R is wrong for N ≥ 3)
  and neither form books Σ p_i ∂V_i/∂λ (`:603-604, 704-705`). newen implements the documented
  formulas (`parameter.h:624-658`; `reference_hamiltonian`, `newen-integrators/src/eds.rs`),
  what native + `patches/D2_eds_weights.patch` computes, pinned by
  `aladip_vacuum_eds_multis_fixed`/`_pairs_fixed`; the released binary pins only the end-state
  energies (and FORM 3's V_R) at frame 0 (`aladip_vacuum_eds_multis`/`_pairs`), and
  `newen-run/tests/eds_forces.rs` checks F = −∇V_R for all three forms. FORM 1 is gromosXX's,
  bit for bit, including its running pairwise log-sum-exp in state order (`eds.cc:409-424`).
- **EDS end-state bookkeeping** (B), reproduced (`newen-forces/src/nonbonded/eds.rs`,
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
  `check_aeds`, `newen-run/src/plan.rs`, and `emax_emin_update`): EMIN > EMAX (warning, then
  EIR unread and a segfault, `in_parameter.cc:4084-4088`; D14); BSTEPS < ASTEPS in FORMS 2/4/5/6
  (unsigned `bsteps − asteps`, τ ≈ 4·10⁹ steps, offsets frozen, `eds.cc:186`; D13 — the pure
  `memory_time` reproduces the wrap); BMAXTYPE 2 with a target state without variance, unless
  the RESTREMIN clamp catches it (division by b_max = 0 at `eds.cc:314`, E_min = −∞, V_R* NaN;
  D15). Pinned by unit tests and `newen-md/tests/aeds_runs.rs` (D15 by the unit test only: the
  `md` input that reached it needs D32's never-visited target, see below).
- **A-EDS defects corrected** (A; corrected in 0.0.63; patched oracle
  `tests/gromosXX_references/patches/D16_*`, `D25_*`, `D32_*`; `newen-integrators/src/eds.rs`):
  - the energy columns (`eds_emax`, `eds_emin`, the offsets) written after the search update, so
    a frame's `eds_vmix`/`eds_vr` belong to the previous frame's printed parameters
    (`eds.cc:368-372`; D16). newen writes the parameters the frame's V_R* used.
  - the FORMS 5/6 prevalences as plain exponentials (`:339-343`, NaN and frozen offsets once
    −β(V_i − E_i^R) > 709.78; D25). newen uses a log-sum-exp (`aladip_vacuum_aeds_f6_ovf_fixed`).
  - a never-visited state 1 as the default E_min target (`:252-261`; D32, its mean 0 a fake G,
    σ 0 the D15 division). newen takes the lowest mean among the states visited more than once,
    state 1 only when none is.
  The `aladip_vacuum_aeds_*_fixed` references pin them against the patched binary;
  `newen-run/tests/aeds_replay.rs` replays them. The released-binary reproducers live in
  `gromos-divergences`. Lean statements: `~/.cache/gromos-rs-item-1.8/lean_aeds_statements.md`
  (not yet in `math`).
- **A-EDS bookkeeping** (B, reproduced behaviour documented in GROMOS Vol. 2, not defects;
  `newen-integrators/src/eds.rs`): the FORMS 5/6 offset step is ASTEPS·Δt (`:355`; D33, not the
  dossier's D17, which is the extended-TI 1-4 scaling) — a
  time if ASTEPS is a step count, an energy if it is the rate Vol. 4 and the tutorial describe;
  τ ramps over the current run only (`sim.steps()` from 0), as Vol. 2 documents ("t_tot is the
  total simulation time of the run") while Perthold 2020 ramps over the whole search, for which
  the input carries no total length; the mean of V_i − E_i^R averages over the offsets of every
  step (`:246`), Vol. 2's \overline{V_i − E_i^R} with E_i^R moving; the boost scales only the
  EDS part of forces, virial and dV/dλ (`:137-176`). Pinned by the `aladip_vacuum_aeds_*`
  references (whole block, 9-decimal AEDSSEARCH) and `newen-run/tests/aeds_replay.rs`
  (gromosXX's V_i(t) through the pure functions, ≤ 1e-9 relative).
- **TI-EDS inputs gromosXX computes wrongly** (A, refused, not reproduced; `check_ti_eds`,
  `newen-run/src/prepare.rs`): an atom in both `PERTATOMPARAM` and `MPERTATOM` gets pair code
  6, which no case handles, so its pairs are dropped silently (`perturbed_nonbonded_innerloop.cc:54-65,
  641-643`; dossier D18); a `PERTATOMPAIR` pair with an EDS atom is computed with the topology
  type and charge in every end state (`perturbed_nonbonded_pair.cc` has no EDS branch; D19).
  PRECALCLAM and slow growth (DLAMT ≠ 0) stay refused with EDS (gromosXX locks both against EDS,
  `check_parameter.cc:4271, 2528-2530`), NSCALE ≠ 0 with or without it.
- **TI-EDS pair ownership** (deliberate difference, same terms): gromosXX's perturbed atom owns
  an excluded or 1-4 EDS–perturbed pair (`PERTATOMPARAM` takes the pair before `MPERTATOM` looks,
  `in_perturbation.cc:1363-1391, 1748-1777`; case 4 at `perturbed_nonbonded_innerloop.cc:2277-2341,
  3696-3743`), newen's EDS atom owns every EDS–perturbed pair (`Topology::apply_eds_perturbation`,
  `EdsAtoms`, `newen-forces/src/nonbonded/eds.rs`): the same kernel per state, summed in another
  order, so agreement is to rounding, not bit for bit (pinned at energy_rel 1e-10 /
  force_abs 1e-8 by `aladip_eds_solv_ti`, `_groups`, `_multis_fixed`; every printed number
  equal). FORMS 2 and 3 book Σ p_i ∂V_i/∂λ as FORM 1 does (the released gromosXX books 0, D2).
- **Slow growth** (DLAMT ≠ 0; reproduced, `newen_run::advance_lambda`, called after every step
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
  `aladip_vacuum_fep_slowgrowth`'s run); newen writes no `.omd` averages.
- **`dGslv_pbsolv` constants** (checked in 1.8 against newen-core; no difference that changes
  a number, nothing to reconcile). 1/(4πε₀) is the topology's `PHYSICALCONSTANTS` value, as in
  `md` (138.9354 in every reference topology, `newen_core::units::four_pi_eps_i`); only PQR
  input falls back to gromos++'s `gmath/Physics.cc:53` 138.935484611
  (`pb::GROMOSPP_FOUR_PI_EPS_I`, relative 6·10⁻⁷ from newen-core's, kept for byte identity).
  k_B and N_A are not used: the NPBC boundary has κ = 0 (`PB_Parameters.cc:69`, no ionic
  strength, so no kT) and the energies come out in kJ/mol through 1/(4πε₀); gromos++'s own
  pair (`Physics.cc:51-52`, k_B = 1.380658·10⁻²³ × 6.0221367·10²³ / 1000 = 0.0083145110 vs
  newen-core's gromosXX 0.00831441) never reaches the program. The water defaults are
  gromos++'s, not newen-core's (which has none): `@epsNPBC` 78.4 (and the FFT solver's
  default ε 78.4), `@probeRAD` 0.14 nm, the hydrogen radius 0.05 nm, fixed because gromos++
  reads `@radH` but does not list it (`dGslv_pbsolv.cc:448, 454, 498`; `PB_Parameters.cc`
  `hradius`); `@epsSOLV` has no default and should be the
  solvent model's permittivity (66.6 for SPC in the gromos++ examples), `@epsRF`/`@rcut` the
  run's EPSRF/RCRF. Radii come from the topology's LJ parameters with the probe's IAC, the
  engine's own table. Pinned by `dgslv_pbsolv_unit_charges_in_spc_match_gromospp`.
- **Constraint virial under the molecular virial** (C, ours, reconciled). gromosXX adds the
  SHAKE/SETTLE virial only for the atomic virial (`shake.h:380`, `settle.cc:397`); newen added
  it under every virial type, so the barostat saw −1500 instead of +28 kJ mol⁻¹ nm⁻³. Now gated
  on `atomic_virial` (`shake_buffered`, `settle` in `newen-integrators/src/constraints.rs`),
  pinned by `aladip_solvated_npt_shake`. Lean: `intragroup_pair_virial_zero`,
  `constraint_virial_pressure_shift`.
- **ROTTRANS degrees of freedom** (C, ours, reconciled). gromosXX removes 3 COM and 3 internal
  dof from atom 1's baths (`topology.cc:1370-1400`); newen did not. Now in `bath_dofs`
  (`newen-run/src/dof.rs`), pinned by `aladip_multibath_rottrans`.
- **Per-bath temperature rows** (C, ours, reconciled). gromosXX's `VOLUMEPRESSURE03` prints
  `T_total`/`T_com`/`T_ir`/`scale` per bath (`out_configuration.cc:3270-3346`); newen printed
  the system temperature on every row. Now per bath (`volume_pressure` in
  `newen-md/src/bin/md.rs`, written by `newen-io/src/energy.rs`), compared whole-block by the
  reference test.
- **Virial booking quirks** (B), reproduced in the force-field `apply`
  (`newen-integrators/src/algorithms/forcefield.rs`), `calculate_bonded_forces`
  (`newen-forces/src/bonded/mod.rs`) and `settle` (`newen-integrators/src/constraints.rs`),
  pinned by `VOLUMEPRESSURE03`:
  - bonded and perturbed bonded virial added under every virial type, the gate commented out
    (`quartic_bond_interaction.cc:93`, `perturbed_quartic_bond_interaction.cc:125`);
  - perturbed 1-4, RF-excluded and atom-pair virials unconditional
    (`perturbed_nonbonded_innerloop.cc:2247-2250, 3941-3946`, `perturbed_nonbonded_pair.cc:1559-1561`),
    plain pairs gated (`nonbonded_set.cc:428-430`);
  - SETTLE gated on the atomic virial with no vacuum split (`settle.cc:397`); LINCS adds no
    constraint virial.
- **Thermostat coupling edge inputs** (C, ours, documented). gromosXX couples a bath whenever
  `tau != -1` (`berendsen_thermostat.cc:99`, `nosehoover_thermostat.cc:123,179`), so any other
  negative τ is coupled with a negative time constant; newen treats every `tau < 0` as uncoupled.
  gromosXX has no guard on `dof` (`2·ekin/(dof·k_B)` is inf or NaN for a 0-dof bath, which it
  decouples at setup, `topology.cc:1405-1411`; ported as `decouple_empty_baths`,
  `newen-run/src/plan.rs`); newen's kernels take `T_free = 0` there and fall back to `T₀`, and
  `kinetic_temperature` prints 0 for every `dof ≤ 0` (gromosXX `out_configuration.cc:3293`
  tests `dof != 0`). Both are in the Lean models (`spec/NewenSpec/Thermostat.lean`); no
  reference reaches either input.
- **Multi-bath Nosé–Hoover** (C, ours, reconciled). newen built the thermostat from the first
  bath only (`NoseHooverThermostat::new_single_bath`/`new_chain_bath` with bath 0's T₀, τ and
  dof) and refused several Nosé–Hoover baths at plan validation; with the refusal lifted it
  indexed a scale it never had for the second bath. Now one ζ chain per bath as gromosXX
  (`nosehoover_thermostat.cc:120-230`); pinned by `aladip_multibath_nosehoover`,
  `aladip_multibath_nhc`.
- **Thermostat output quirks** (B), reproduced: `NDFMIN` not subtracted from the COM/IR dof
  (`multibath.cc:139-149`; `bath_dofs`, `newen-run/src/dof.rs`); a bath's `scale` printed 0
  without a thermostat (`multibath.h:48`, `berendsen_thermostat.cc:134`); the kinetic-tensor rows
  printed unrotated while the trace uses the rotated tensor (`out_configuration.cc:3331-3344`);
  the last two in `volume_pressure` (`newen-md/src/bin/md.rs`).
- **Nosé–Hoover variables of a continued run** (A, D41; corrected in 0.0.69; patched oracle
  `tests/gromosXX_references/patches/D41_nh_restart.patch`).
  gromosXX writes and reads `NHCVARIABLES` for chains only (`out_configuration.cc:606`,
  `in_configuration.cc:1156`) although plain Nosé–Hoover integrates `zeta[0]`
  (`nosehoover_thermostat.cc:142-143`), and refuses NTINHT = 1 without a chain
  (`check_parameter.cc:92`), so a continued plain-NH run restarts from ζ = 0; for a chain it reads
  `zeta.size()` = 1 value per bath before `NoseHoover_Thermostat::init` sizes the chain
  (`in_configuration.cc:3444`, `md.cc:150, 167`, `nosehoover_thermostat.cc:69`), so links 2..M
  restart from 0. newen wrote and read nothing until 0.0.69 (NTINHT parsed, unused). newen
  writes every variable of every bath in gromosXX's block format, plain NH included, and
  restores them when NTINHT = 0 (`restored_zeta`,
  `newen-integrators/src/algorithms/nosehoover_thermostat.rs`; `prepare_system` parses the block
  into `Configuration::nose_hoover`); NTINHT = 1 starts from zero, plain NH included; NTINHT = 0
  without a block starts from zero too, with a warning (gromosXX refuses that for a chain).
  Pinned by `water_216_nvt_nosehoover_cont_fixed` / `water_216_nvt_nhc_cont_fixed` (patched);
  `newen-md/tests/nh_restart.rs` checks that 100 + 100 steps follow the 200-step run as closely
  as weak coupling does. The released-binary reproducer lives in `gromos-divergences`. Lean
  (`gromos-divergences`): `nhScale_restart_eq_iff`, `nhScale_restart_sub`,
  `nhcReadTranscribed_eq_iff`.
- **NDFMIN spread over the baths through an `int` total** (A, D40; corrected; patched oracle
  `tests/gromosXX_references/patches/D40_ndfmin_total.patch`). gromosXX removes NDFMIN from each
  bath in proportion `dof/dof_total`, but sums `dof_total` in an `int` (`multibath.cc:145-152`);
  a DOFSET boundary that cuts an odd number of constraints leaves half-integer baths (each
  constraint takes ½ from the IR bath of each end, `topology.cc:1105-1106`), the total is
  truncated and more than NDFMIN is removed (14.5 and 10.5 with NDFMIN 3: 12.69 and 9.19 instead
  of 12.76 and 9.24). newen sums the total exactly and spreads nothing when it is not positive
  (`bath_dofs`, `newen-run/src/dof.rs`). Pinned by `aladip_vacuum_multibath_ndfmin3_half_fixed`
  (patched) next to `aladip_vacuum_multibath_ndfmin3` (integer baths, released binary). The
  released-binary reproducer lives in `gromos-divergences`. Lean (`gromos-divergences`):
  `dofTotalTranscribed_le_iff`, `ndfminSpread_removes_more`.
- **Perturbed RF self term, "Chris: CHECK!"** (false alarm; not a quirk). The comment at
  `perturbed_nonbonded_term.cc:749-751` (`gromosXX-patched`; `:596-598` in the non-ext
  variant) worries that the `selfterm_correction` adds `+A_q·FPEPSI·crf_cut` with no B-state
  term and no matching `de_rf`. It needs neither: the term is λ-independent, it only cancels the
  state-A self energy the unperturbed excluded loop already booked, and the B state enters
  through the λ-weighted `A_e_rf`/`B_e_rf` above it, whose derivative is in `de_rf`. The net
  self energy is −½·FPEPSI·crf_cut·[(1−λ)ⁿq_A² + λⁿq_B²], and its TI integral is exactly
  −½·FPEPSI·crf_cut·(q_B² − q_A²) (`perturbed_self_energy_correction` ↔
  `newen_analysis::free_energy::rf_self_energy`, test
  `rf_self_energy_is_the_ti_integral_of_the_perturbed_self_term`). Reproduced as is.
- **Soft improper ζ wrap** (C, ours, dossier candidate). newen's soft improper caller wraps ζ into
  `[ζ0 − π, ζ0 + π]` before the soft kernel, as the plain and perturbed impropers do;
  gromosXX's `perturbed_soft_improper_interaction.cc:134-138` does not wrap, so a ζ a turn away
  from ζ0 gets the energy of the long way round there. newen keeps the wrap.
- **LJEXCEPTIONS** (refused, not implemented). gromosXX moves each pair into `all_exclusion` and
  gives it its own C6/C12 (`in_topology.cc:2111`); newen's topology reader refuses a block with
  NEX > 0 ("LJEXCEPTIONS block is not implemented") instead of skipping it as an unknown block.
- **ROTTRANS across baths** (ported). gromosXX refuses roto-translational constraints whose
  atoms `1..=RTCLAST` span more than one (COM, IR) bath pair (`topology.cc:1370-1394`);
  `thermostat_baths` (`newen-run/src/plan.rs`) does the same, with gromosXX's message.
- **Not ported, no reference**: position-restraint and distance/dihedral-restraint virials
  (`position_restraint_interaction.cc:89`, `distance_restraint_interaction.cc:194`).

## In one line

Tests before wiring; derivation before transcription; grep before porting; document every deliberate
quirk; one behaviour per feature — a gromosXX defect is corrected in newen and patched in the oracle,
never kept behind a switch. Faithfulness and correctness are both required — the discipline above is
how you get both instead of trading one for the other.
