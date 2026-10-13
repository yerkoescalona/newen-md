# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [0.0.55] (2026-10-12)

**Two defects found by step 2's gates, fixed as gromosXX does them, and one gromosXX defect
reproduced and flagged** (1.8). The reference test now compares `VOLUMEPRESSURE03` and
`CONSFORCERED` whole-block; two new native references, and the ignored
`aladip_vacuum_fep_pairpert_precalc` now passes; 79 reference systems pass.

### Fixed

- **`PERTATOMPAIR` pairs with a perturbed atom.** The interaction types are read raw (0 excluded,
  1 normal, 2 1-4) as `in_perturbation.cc` does; gromos-rs subtracted one and underflowed on 0.
  Every pair now leaves the exclusion and 1-4 lists at prepare time and is computed once, per
  state, in the atom-pair term with the A/B types, charges and soft core
  (`perturbed_nonbonded_pair.cc:179-298`); its r⊗f goes to the virial and its λ-grid energies to
  `PRECALCLAM`. `pt_top` mapped the B type off by one (type 1 added an exclusion, type 2 nothing);
  it now matches gromos++ `PtTopology.cc:151-170`. Pinned by `aladip_vacuum_fep_pairpert_precalc`
  (no longer ignored: pairs 4-10 1-4/excluded, 1-11 normal/1-4 with both atoms perturbed, 4-7
  excluded/normal); 77 references whole-block at this point. The eleven aladip references that
  carry the old 7-10 pair change only in the `.tre` virial (now r⊗f, as gromosXX) and one 1e-34
  bath value.
- **The constraint virial is added only for the atomic virial** (`shake.h:380`, `settle.cc:397`);
  SHAKE and SETTLE added it under every virial type. Under the molecular virial the barostat saw
  −1500 instead of +28 kJ mol⁻¹ nm⁻³ and shrank the LiveCoMS `t_02` box by 0.09 % per step.
  Lean: `intragroup_pair_virial_zero`, `constraint_virial_pressure_shift`
  (`math`, `OAI.Physics.Gromos.Divergences`). With it, as gromosXX:
  - roto-translational constraints remove 3 COM and 3 internal degrees of freedom from atom 1's
    baths (`topology.cc:1370-1400`), and the `ROTTRANS` reference is centred on the
    chain-gathered COM;
  - the `.tre` writes per-bath `T_total`/`T_com`/`T_ir`/`scale` and gromosXX's
    `VOLUMEPRESSURE03` layout;
  - one pressure, trace/3, for the barostat and the output;
  - `CONSFORCERED` is scaled by 1/dt² once per atom;
  - the perturbed bonded, 1-4, RF-excluded and atom-pair terms contribute their virial.

  Nothing pinned this: `water_216_npt` runs no SHAKE and the test compared `ENERGY03` only. It now
  compares the whole `VOLUMEPRESSURE03` block and `CONSFORCERED`, under which the previous binary
  fails 22 systems; two new references, `aladip_solvated_npt_shake` and
  `aladip_multibath_rottrans`.

### Changed

- `Simulation.pressure` (Python) returns trace(P)/3, the number the barostat couples to and
  gromosXX writes, and is 0 under NVT/NVE (it returned the kinetic term alone).
- `PERTATOMPAIR` interaction codes other than 0/1/2 are a parse error (gromosXX silently zeroes
  the state).

### Reproduced gromosXX defects (flagged, to be corrected in a later version)

- **`eds_pert_lj_crf_interaction` called with its outputs shifted by one slot**
  (`perturbed_nonbonded_pair.cc:469, 677, 992, 1200` against `perturbed_nonbonded_term.h:65-72`).
  For a `PERTATOMPAIR` pair with a perturbed atom and a normal or 1-4 state type, the state's
  energies receive the kernel's dH/dλ, its force is `(f + e_lj + e_crf)·r`, its dH/dλ is 0, and
  `eps = 0` and `coulomb_scaling = 0` are fed into the kernel, so its 1/r Coulomb term vanishes.
  Reproduced in `atom_pair_state` (`gromos-forces/src/nonbonded/perturbed.rs`), pinned by
  `aladip_vacuum_fep_pairpert_precalc`. Lean: `shiftedPairForce_eq_pairForce_iff`,
  `shifted_ti_shortfall`, `crfSoftEnergy_shift` (`math`, `OAI.Physics.Gromos.Divergences`).
- Reproduced quirks, now pinned by `VOLUMEPRESSURE03`:
  - the bonded and perturbed bonded virial is added under every virial type, the gate commented
    out (`quartic_bond_interaction.cc:93`, `perturbed_quartic_bond_interaction.cc:125`);
  - the perturbed 1-4, RF-excluded and atom-pair virials are unconditional
    (`perturbed_nonbonded_innerloop.cc:2247-2250, 3941-3946`, `perturbed_nonbonded_pair.cc:1559-1561`),
    while plain pairs are gated (`nonbonded_set.cc:428-430`);
  - SETTLE is gated on the atomic virial with no vacuum split (`settle.cc:397`); LINCS adds no
    constraint virial;
  - `NDFMIN` is not subtracted from the COM/IR dof (`multibath.cc:139-149`);
  - a bath's `scale` is printed 0 without a thermostat (`multibath.h:48`,
    `berendsen_thermostat.cc:134`);
  - the kinetic-tensor rows are printed unrotated while the trace uses the rotated tensor
    (`out_configuration.cc:3331-3344`).
- Not ported (no reference): the position-restraint and distance/dihedral-restraint virials
  (`position_restraint_interaction.cc:89`, `distance_restraint_interaction.cc:194`); baths with
  dof ≤ 0 losing their coupling (`topology.cc:1405-1411`).

### Verified

- References: 79 passed, 1 ignored (`aladip_gamd`), whole-block, `VOLUMEPRESSURE03` and
  `CONSFORCERED` included; check suite, `fep_null_perturbation` and 563 unit tests green.
- `ENERGY03`, forces, trajectories and final configurations are byte-identical to the previous
  binary on every old reference.
- LiveCoMS `t_02` over 500 steps matches native gromosXX to |ΔV|/V 3e-10 and 2e-5 K per bath.

## [0.0.54] (2026-10-12)

**Extended TI (`PRECALCLAM`): the energies of every perturbed interaction and their λ derivatives
at `NRLAM` extra λ points, written as gromosXX writes them; and the reaction field built as
gromosXX builds it at `EPSRF = 0` and from `RCRF`** (1.8 step 2). Nine new references generated
by native gromosXX pin it: 76 reference systems pass, whole-block.

### Added

- **`PRECALCLAM NRLAM MINLAM MAXLAM`** (`in_parameter.cc:4478-4513`, `check_parameter.cc:261-266`):
  gromosXX's range checks, perturbation on and `NTWE = NTWG > 0` required; refused with lattice
  sum, the conjugate-gradient minimiser, QM/ML terms, GAMD and EDS, as gromosXX does. Allowed
  with `LAMBDAS`, perturbed SHAKE and perturbed distance restraints. `Recipe.perturbation.precalclam`
  (`{nr_lambdas, min_lambda, max_lambda}`, a dict in Python, `None` when absent); `imd_write`
  emits the block.
- On every NTWG step the grid fills, per λ_k, `A/B_lj` and `A/B_crf` per energy-group pair (pairs,
  RF self and exclusions, perturbed 1-4, perturbed atom pairs), `AB_kinetic`, `AB_bond`,
  `AB_angle`, `AB_improper`, `AB_disres`, and the `A_dihedral`/`B_dihedral` scalars; folded and
  summed in `calculate_totals`' order (`energy.cc:397-423`), not part of the potential total.
  `AB_angres`, `AB_dihres`, `AB_disfld` stay 0 (those restraints are not modelled).
- The `.tre` and `.trg` carry `# precalclam` (`# nr_lambdas`, N, the 12 columns) and `# ABdih`
  (`out_configuration.cc:3226-3266`); the `.trg` prints the derivative record. With
  `NRLAM = 0` nothing changes.
- **Nine references**: `ch4_water_fep_precalc`, `meoh_water_fep_precalc` (charge changes: CRF
  pair, RF self and exclusion terms), `aladip_vacuum_fep_precalc` (bonded `AB_*`, `AB_kinetic`,
  `ABdih`), `aladip_vacuum_fep_lambdas_precalc` (`LAMBDAS` with `NRLAM 6 0.2 0.7`),
  `aladip_vacuum_fep_shake_precalc`, `aladip_vacuum_fep_distres_precalc`,
  `ch4_water_fep_precalc_twinrange` (pairlist every 5 steps, output every 2), and
  `aladip_solvated_epsrf0`, `aladip_solvated_rcrf`. `aladip_vacuum_fep_pairpert_precalc` is
  generated and ignored (see Known, deferred).

### Fixed

- **`EPSRF = 0`** (conducting boundary): gromosXX sets `crf = −1` (`nonbonded_term.cc:58-71`,
  `perturbed_nonbonded_term.cc:51-61`); gromos-rs used the ε formula. Step-0 CRF of
  `aladip_solvated_epsrf0` was −56.26, gromosXX −60.87; now equal, in every kernel.
- **The reaction-field constants come from `RCRF`**, not `RCUTL` (`nonbonded.rf_cutoff`,
  `in_parameter.cc:4606`). Step-0 CRF of `aladip_solvated_rcrf` (RCRF 1.2, RCUTL 1.4) was −60.83,
  gromosXX −61.72; now equal. gromosXX accepts `RCRF ≠ RCUTL` silently (`check_parameter.cc` has
  no rule); so does gromos-rs. Every other reference has `RCRF = RCUTL` and `EPSRF ≠ 0` and is
  byte-identical.

### Verified

- References: 76 passed, 2 ignored (`aladip_gamd`, `aladip_vacuum_fep_pairpert_precalc`),
  whole-block; check suite 10 passed, 1 ignored; `fep_null_perturbation` 2 passed; unit crates
  green; the MPI build succeeds.
- `NRLAM = 0`: every `.tre`/`.trg`/`.trf`/`.trc`/final configuration is byte-identical to 0.0.53
  on all 68 previous references.
- gromos++ `ene_ana` prints identical averages on our `ch4_water_fep_precalc` files and
  gromosXX's.
- LiveCoMS tutorial `t_02/DDM/md_TI/md_TI_ASA_Na.imd` runs to completion (50 000 steps, 292 s).

### Measured

Wall time, one thread, 500 steps, interleaved, n = 5, mean ± std (BENCHMARKING.md §9):

| system | 0.0.53 | this version | gromosXX | gromosXX / this |
|---|---|---|---|---|
| `ch4_water_fep_l025`, outputs off | 11.33 ± 0.11 s | 11.30 ± 0.08 s | — | — |
| `ch4_water_fep_l025`, NTWE = NTWG = 1 | — | 11.31 ± 0.06 s | 11.93 ± 0.02 s | 1.055 ± 0.006 |
| `ch4_water_fep_precalc`, `PRECALCLAM 81 0 1`, NTWE = NTWG = 1 | — | 12.15 ± 0.05 s | 13.14 ± 0.02 s | 1.081 ± 0.005 |

Without `PRECALCLAM` nothing changes beyond noise. 81 λ points on every step cost 7 % here and
10 % in gromosXX.

### Deliberate differences from gromosXX

- `NRLAM = 1` is refused: gromosXX divides by `NRLAM − 1` and every λ_k is NaN.
- Conventions of the method (de Ruiter & Oostenbrink, JCTC 2016, 12, 4476; not divergences):
  - the derivative record's signs, `A_de = +∂A/∂λ`, `B_de = −∂B/∂λ`: the B-state softness
    derivative is taken in `1 − μ` (Eqs. 9, 12; `perturbed_nonbonded_term.cc:614-716`);
  - the grid is in the individual μ, not λ, so `LAMBDAS` and `NLAM` are applied afterwards by
    `ext_ti_ana` (p. 7, Eqs. 10–13); `AB_kinetic` keeps the simulated mass μ, `slam`, as momenta
    are conserved (Eqs. 17, 19; `temperature_calculation.cc:200-224`).
- Reproduced quirks (`.claude/porting.md`):
  - the long-range ext energies are refreshed only on pairlist steps that are also output steps,
    and are zero on the output steps after one that is not (`perturbed_nonbonded_set.cc:87-96,
    280-300`);
  - the improper angle is wrapped once, at the current λ, and not re-wrapped per λ_k;
  - perturbed SHAKE books into the `AB_bond` derivative only (`perturbed_shake.cc:229-259`;
    gromosXX's warning, `in_perturbation.cc:327-334`, is a prepare note);
  - two bond-mixing forms are kept for bit-identity: `(1 − λ_k)·b0_A + λ_k·b0_B` on the grid,
    `r0_A + L·Δr0` in the main term.
- The two reaction-field fixes above (`EPSRF = 0` → `crf = −1`; constants from `RCRF`, not
  `RCUTL`, accepted silently as gromosXX does).

### Known, deferred

- `PERTATOMPAIR` pairs with a perturbed atom are wrong in the main energy (`ptp.rs` type
  underflow, list surgery, soft core); pinned by the ignored
  `aladip_vacuum_fep_pairpert_precalc`, fixed next.

## [0.0.53] (2026-10-11)

**Individual λ per interaction and energy-group pair (`LAMBDAS`), the λ exponent done as gromosXX
does it, and the dH/dλ of perturbed SHAKE and perturbed distance restraints** (1.8 step 1). Thirteen
new references generated by native gromosXX pin it, all compared value by value in both energy
blocks: 67 reference systems pass.

### Added

- **`LAMBDAS`** (`in_parameter.cc:4326-4472`): `NTIL on|off|1|0`, rows `NTLI NILG1 NILG2 ALI…ELI`
  with `NTLI` a name or 1–12 (`angleres` is documented by gromosXX and rejected by its parser; so
  here), groups 1-based and bounded by `NEGR`, `NILG2 ≥ NILG1`, a warning for an off-diagonal row
  of a bonded kind, symmetric storage, the last duplicate wins, left-over tokens an error, and
  without `NTG=1` the block is ignored with gromosXX's warning. One table in gromos-core
  (`gromos_core::lambdas::IndividualLambdas`, 14 kinds in gromosXX's order) evaluated by the
  topology's `update_for_lambda` with gromosXX's exact expression (`topology.cc:1442-1455`);
  every perturbed term reads its own kind at `[eg(i)][eg(j)]` (nonbonded, RF self at `[g][g]`,
  1-4, excluded, atom pairs) or `[eg(first atom)][same]` (bonded, restraints, SHAKE, masses,
  kinetic), and multiplies its dH/dλ by dλ_i/dλ. `Recipe.perturbation.lambdas` carries the rows
  (name or number for the kind); `imd_write` emits the block.
- **`NRDGL=1`** reads λ from the configuration's `PERTDATA` block, as gromosXX does
  (`in_configuration.cc:876-897, 2351-2371`), and errors without one; the final configuration
  now carries `PERTDATA` whenever perturbation is on, and `DISRESEXPAVE` with time-averaged
  restraints.
- **Perturbed SHAKE** (`perturbed_shake.cc`): with `NTC=3` the perturbed bonds become constraints
  of length `r0_A + L_bond·(r0_B − r0_A)`, shaken before the unperturbed ones, and their dH/dλ is
  booked per group in the constraints column of the `.trg`. `NTC=2` constrains exactly the
  `BONDH` bonds, which the topology now remembers (`Bond.hydrogen`); a perturbed BONDH bond under
  `NTC=2` is refused as gromosXX does.
- **Perturbed distance restraints** are a full port of
  `perturbed_distance_restraint_interaction.cc` (`r0`/`w0` mixed in L, the `2^(n+m) L^n (1−L)^m`
  prefactor, energy, force and dH/dλ per group), and time averaging (`NTDIR < 0`, the `r⁻³`
  running average) now works for perturbed and unperturbed restraints; before this the average
  was never computed.
- **Thirteen references**: `ch4_water_fep_lambdas` (vdw/vdw_soft rows on the solute–solvent
  pair), `ch4_water_fep_lambdas_off` (`NTIL off`, identical to `ch4_water_fep_l025` outside the
  title), `ch4_water_fep_nlam2`, `ch4_water_fep_nrdgl`, `aladip_vacuum_fep_lambdas`
  (bond/angle/dihedral/improper/mass on (1,1), (2,2), nonbonded on every pair),
  `aladip_vacuum_fep_lambdas_nlam2`, `aladip_vacuum_fep_nlam2`, `aladip_vacuum_fep_shake`,
  `aladip_vacuum_fep_shake_lambdas`, `aladip_vacuum_fep_shake_hblock` (H and non-H perturbed bond
  blocks sharing an atom, a group-2 bond that changes length), `aladip_vacuum_fep_distres`,
  `aladip_vacuum_fep_distres_lambdas`, `aladip_vacuum_fep_distres_tavg`.
- gromosXX's `check/lambdas.t.cc` ported into the check suite: per-group energies of the `on`
  run equal the `off` run at `λ = L`, and `dE_on = dE_off·L′`, for bond, angle, dihedral and
  improper; its `NonBonded` values need `SCALEDINTERACTIONS` and are left out.
- `fep_null_perturbation`: identity `LAMBDAS` rows on a null perturbation stay bit-identical.

### Fixed

- **`NLAM ≠ 1` was applied twice and to the wrong terms.** gromos-rs fed `RLAM^NLAM` to every
  perturbed term and then raised it to `NLAM` again in the soft-core kernel. gromosXX evaluates
  the polynomial on `RLAM`, uses `(1−L)^n` and `L^n` in the nonbonded term only, and keeps bonded
  terms, restraints and masses linear in L (`perturbed_nonbonded_term.cc:322-347`); `powi` became
  `powf` to match `pow`. `NLAM=1` outputs are byte-identical to before.
- The constraint dH/dλ was computed and discarded; the perturbed distance-restraint dH/dλ was
  never computed. Both reach the `.trg` now, and `Simulation.dhdl` returns the full derivative
  (kinetic + potential + special), what the `.trg` writes and `ext_ti_ana` integrates;
  `dhdl_terms` gains `kinetic`, `constraints`, `distanceres`.
- `DLAMT ≠ 0` (slow growth) and `NSCALE ≠ 0` (`SCALEDINTERACTIONS`) were parsed and silently
  ignored; both are refused with a message. The refusals run on every build path, including a
  Python `Simulation(system, recipe, plan=…)`.
- `PERTDISRESSPEC` restraints were applied with `NTG=0`; gromosXX ignores them then.
- `PERTBONDSTRETCHH`, `PERTBONDANGLEH`, `PERTIMPROPERDIHH` and `PERTPROPERDIHH` were read after
  the plain blocks; gromosXX reads them first (`in_perturbation.cc:172-175, 590, 827, 1058`),
  which fixes the SHAKE order when both blocks share an atom.

### Changed

- The plan no longer carries `lambda`/`nlam`: λ lives on the prepared topology (`Topology.lambda`,
  `lambda_exp`, `individual_lambdas`), one source for the force field, the constraints, the
  restraints and the temperature calculation. The two dead `InteractionLambdas`/`LambdaController`
  pairs are gone.

### Deliberate differences from gromosXX

- A perturbed atom pair whose group pair has `L_lj ≠ L_crf`: gromosXX only logs a critical
  message, scales each state's combined LJ+RF force by the LJ factor and finishes the run
  (`perturbed_nonbonded_pair.cc:309-312, 405`). Reproduced (a prepare note says so). With
  `L_lj = L_crf` the previous per-kind sum is kept so NLAM=1 runs stay byte-identical; the last
  bits differ from gromosXX for a charged perturbed 1-4 pair.
- `LINCS` with `NTC=3` and perturbed bonds: gromosXX moves the bonds to constraints LINCS never
  shakes and only warns, leaving them with neither constraint nor force
  (`in_perturbation.cc:265-303`, `create_constraints.cc:139-143`). Refused instead.
- `FORCESCALE 1/2` and `NTDIRA=1` with time-averaged restraints are refused (not implemented);
  gromosXX runs them.
- The `…03` perturbation blocks are still accepted (gromosXX errors on them); "Perturbation of a
  non-existing bond" is still not checked. Both predate this version.

### Still open

- `SCALEDINTERACTIONS` / `NSCALE`, slow growth (`DLAMT`), replica-exchange λ, `FORCESCALE 1/2`.
- `shake_buffered` scales the constraint forces by `1/dt²` once per constraint rather than once per
  atom (gromosXX `perturbed_shake.cc:400-406`); nothing compares `CONSFORCERED` yet.

## [0.0.52] (2026-10-11)

**Energies are booked per energy group the way gromosXX books them, and the reference test now
compares every number of the energy blocks** (1.8 step 0). Before this, the `# nonbonded` rows of
every `.tre` were zero, `# bonded` filled group 1 only, `# special` and the centre-of-mass and
internal kinetic energies were zero, and the `.trg` put all bonded dH/dλ in the bond slot. Nothing
noticed, because the test compared three totals per frame. It now compares all 52 totals and every
row of `# baths`, `# bonded`, `# nonbonded`, `# special`, `# eds`, `# gamd`, `# precalclam` and
`# ABdih`, in the `.tre` and in the `.trg`, against native gromosXX: 54 references pass that way.

### Added

- **Energy groups as gromosXX assigns them** (`in_topology.cc:433-459`): `prepare_system` maps
  every atom, solvent included, to its FORCE group. An NRE list that stops before the last atom
  gets one more group, reported as `PrepareNote::EnergyGroupAdded`; an NRE beyond the last atom
  or not strictly increasing is an error.
- **Every energy and dH/dλ booked into its group cell**: nonbonded pairs, 1-4 pairs, RF excluded
  pairs and the RF self term into `[eg(i)][eg(j)]`, the perturbed corrections into the same cells,
  bonded terms by their first atom (cross-dihedrals by `a`), position and distance restraints by
  the restrained atom. The kernels take an `EnergyGroups` argument and choose a const-generic
  instantiation once per call, so a one-group run compiles to the previous code path.
- **gromosXX's energy layout.** `Energy` carries the per-group arrays and matrices of
  `configuration::Energy` and the centre-of-mass/internal split of the kinetic energy per bath
  (`temperature_calculation.cc:60-113`); dH/dλ is a second `Energy` on
  `State::perturbed_energy_derivatives`, as `conf.current().perturbed_energy_derivatives` is.
- **The kinetic dH/dλ of perturbed masses** (`temperature_calculation.cc:116-237`), and dH/dλ
  per bonded term (bond, angle, improper, dihedral) instead of one bonded sum.
- **One writer for both blocks**: `gromos_io::energy::write_energyred_body` prints
  `_print_energyred_helper` for the `.tre` (`EnergyWriter::write_energy`) and the `.trg`
  (`FreeEnergyWriter::write_derivatives`), including slot 40 of the `.trg` (dH/dλ of the
  potential times the potential, `energy_calculation.cc:87-89`), slots 47/48 (total plus the
  shifted-RF extras) and the one empty GaMD group gromosXX always prints in the `.tre`.
- **Five references** generated by native gromosXX: `aladip_vacuum_fep_groups` (two groups, per-term
  dH/dλ per group, kinetic dH/dλ), `aladip_vacuum_fep_trg` (one group with a `.trg`),
  `aladip_multibath_swapped` (DOFSET `12 2 2  72 1 1`), `aladip_solvated_em_posres_groups` and
  `aladip_harmonic_covalent_groups` (position restraints and harmonic bonded terms split over two
  solute groups).
- `scripts/bench_engines.py --rust-md PATH` compares several Rust builds on the same patched
  input, interleaved.

### Fixed

- **Atoms went to the wrong thermostat bath when DOFSET bath numbers do not rise with atom
  order.** The plan kept one last atom per bath and read the baths as contiguous ranges in bath
  order; with solute in bath 2 and water in bath 1, the solute was scaled by bath 1. Every
  consumer (kinetic energy per bath, the Berendsen and Nosé–Hoover scaling, the kinetic dH/dλ) now
  uses gromosXX's `in_bath` over the DOFSET ranges (`multibath.h:159-171`);
  `aladip_multibath_swapped` pins it, trajectory included.
- A temperature group spread over two baths is refused with gromosXX's message
  (`multibath.cc:161-181`); it used to run, with the bath energies and the scaling disagreeing.
- A TEMPERATUREGROUPS block whose last group ends before the last solute atom is refused when an
  MD run is prepared (`in_topology.cc:278-290`); the shared topology parser still reads it, as
  gromos++ does.
- FORCE now refuses NRE values that are not positive or not strictly increasing
  (`in_parameter.cc:986-993`).
- `Topology::num_solute_atoms()` recursed forever on a topology without instances.
- Under MPI, ranks other than 0 wrote `md.tre`, `md.trc`, `md.trg` and `md.recipe.toml` into
  the working directory. Only rank 0 writes files now.
- The `.trg` comparison now requires the same number of frames as gromosXX; extra frames used to
  pass unseen.

### Changed

- `Energy` no longer has `dhdl_total`, `dhdl_bonded`, `dhdl_lj` and `dhdl_crf`; read
  `State::perturbed_energy_derivatives` instead. `EnergyFrame` and `FreeEnergyFrame` are
  scalar summaries; the engine writes both files from the energy records. The Python API is
  unchanged (`Simulation.dhdl` is the potential part of dH/dλ, as before).

### Deliberate differences from gromosXX

- Per-cell summation order differs from gromosXX (long-range terms are added before the local
  terms, 1-4 pairs are summed separately, perturbed terms are applied as a correction to state A),
  so group cells agree within the reference tolerances, not bit for bit. The totals are summed as
  before and are unchanged on every reference.
- gromosXX's "Multibath not sorted" check never fires (`multibath.cc:189-196` never updates the
  running last atom), so it accepts DOFSET lines whose last atom goes backwards and then counts a
  negative number of degrees of freedom. gromos-rs refuses them.
- `thermostat.cc:66-67` skips a range whose bath scale is exactly 1.0 with `continue`, which also
  skips the update of the range's last atom, so the next range scales the skipped atoms too. This
  is a defect and is not reproduced.

### Measured

Wall time of `md`, one thread, 500 steps, outputs off, interleaved, n = 5, mean ± std:

| system (energy groups) | 0.0.50 | this version |
|---|---|---|
| `ch4_water_fep` (2) | 11.70 ± 0.41 s | 12.13 ± 0.41 s |
| `water_216_box` (1) | 0.314 ± 0.009 s | 0.316 ± 0.009 s |

One group costs nothing measurable. Two groups cost about 4 % ± 5 %, mostly in the long-range
pairs. A per-pair read-modify-write of the group cell cost 14.5 %; a register accumulator that
writes the cell only when it changes brought it down to these figures.

### Still open

- `VOLUMEPRESSURE03` still writes one temperature per bath without the centre-of-mass and internal
  split and without the scaling column.
- The dH/dλ of perturbed distance restraints and of perturbed SHAKE are not computed yet (1.8 step 1).

## [0.0.50] (2026-10-06)

**A first MD run for newcomers: two argon atoms in 30 minutes.** The mdBook tutorial
`docs/src/tutorials/first-md.md` is written, and `py-gromos/notebooks/03_two_argon_atoms.ipynb`
is its runnable version. Each experiment asks for a prediction first, and the notebook asserts
its own numbers, so CI fails if the engine stops reproducing them.

### Added

- **`03_two_argon_atoms.ipynb`** on the `pair_lj` reference system, in four sections: the system
  and the recipe (the step-0 energy against `C12/r¹² − C6/r⁶`); the plan (three edits that
  `Plan.validate()` refuses, with their messages asserted); the time step (the total-energy range
  over 8 ps for dt from 0.0005 to 0.128 ps, log–log slope asserted in [1.9, 2.1], measured 2.01);
  the well depth (the lowest potential energy of the run within 0.1% of `−C6²/(4·C12)`).
  It states its time budget, what it does not teach, and the three most likely errors with
  their fixes. Executed headlessly, it takes about 4 s on a desktop machine.

### Fixed

- **`Simulation.run()` stamped every energy frame one time step late.** A frame's energies are
  evaluated at the positions its step started from, but the time came from the step counter
  after it had advanced, so frame 0 read `t = dt`. It now reads the time gromosXX's `.tre`
  gives the same step; `test_run_frames_carry_the_time_gromosxx_gives_their_energies` holds the
  time column and the total energy to `pair_lj`'s expected `.tre`.
- **`01_load_and_inspect.ipynb` no longer ran.** It cast positions to float32 for `rdf()`, which
  has taken float64 since 0.0.31, and section 7 expected `System.from_files()` to refuse an
  unsolvated topology with a solvated configuration, which it has accepted since 0.0.30 (the
  solvent count comes from the configuration, as in GROMOS). Both cells now show what happens.
- `gromos.pyi` typed `rmsd`, `rdf` and `Frame.positions_array` as float32; the extension takes
  and returns float64. The API reference's `rdf` example also used keyword names that do not
  exist (`group_a`/`group_b` for `group1`/`group2`).
- Learner-facing files no longer cite the retired plan file or its item labels (notebook 01,
  the design mockup, `learning-guide.md`, `pairlist.md`, the FEP tutorial stub).
  `docs/book.toml` pointed at a misspelt GitHub account. `CONTRIBUTING.md` was titled for
  gromosXX, asked for Python 3.9+ (the package needs 3.12+), and linked a style guide that does
  not exist; it now links the architecture and porting guides.

### Changed

- The superseded API design mockup moved from `py-gromos/notebooks/` to
  `py-gromos/docs/development/api-design-mockup.ipynb`. It cannot run (it is a sketch against
  files that do not exist), and the notebook job executes everything in `notebooks/`.

## [0.0.49] (2026-09-14)

**The 33 gromos-tools programs are library functions, and Python can call them.** `gromos.tools`
runs any of them in process — the same Rust code the command-line binary runs, so there is no
binary to find, nothing to ship beside the wheel, and no way for the two to disagree (3.10 steps 0–5).

### Added

- **`gromos.tools`** — the gromos++ system-preparation programs from Python.
  `tools.run("sim_box", topo=…, pbc="r", minwall=[1.0])` returns exactly what the binary prints;
  keywords are the program's own `@key` names (a list is a multi-valued key, `True` a flag,
  `None`/`False` "not given"). `tools.programs()` lists all 33 with their accepted keys, usage and
  **status** — `verified` / `ported` / `stand_in` / `stub` / `no_counterpart`, so a caller can see
  whether they are calling a port of gromos++'s program or a utility that shares its name.
  A failure raises `gromos.ProgramError`, whose `stdout` carries whatever the program had written:
  `check_top` prints its whole report and *then* fails, and that report is the answer.
- **Typed wrappers for the preparation chain**: `make_top`, `com_top`, `red_top`, `pdb2g96`,
  `sim_box`, `build_box`, `ion`, `check_top` take and return `Topology` / `Configuration` /
  `CheckReport` objects, so a workflow needs no file between steps. Any input that names a file
  accepts a path *or* the object. `Topology.from_string` / `Configuration.from_string` /
  `str(obj)` / `.to_g96()` / `.write()` make the object and its file text interchangeable.
- **`System.solvate()` / `System.add_ions()` / `System.check()`** — the same verbs as methods,
  returning a *new* System (an input system is never mutated). `add_ions` does what gromos++'s
  `ion` leaves to the caller: it builds each ion's topology from the `(mtb, ifp)` pair and
  combines it with the solute's, so both halves of the result agree.
- **A program registry** (`gromos_io::args::Program`, `gromos_tools::PROGRAMS`): name, accepted
  keys, usage, status and the `run` function. `gromos_io::run_main` is every binary's whole `main`.
- `gromos_io::normalize(argv)` — `gromos_args()` without the process: `@key` → `--key` and `@f`
  expansion returning a `Result`, so the command line and an in-process caller expand argument
  files identically. `gromos_io::topology::parse_topology` and `coordinate::parse_g96_labeled`
  read from any reader, not only a path.

### Changed

- **25 gromos-analysis programs moved into the crate library too**
  (`src/programs/<category>/<name>.rs`, the same four-part shape and the same category mirroring as
  gromos-tools): the 18 that already parsed with `Arguments`, and the 7 that read `gromos_args()`
  and walked the vector by hand — `atominfo`, `bar`, `ener`, `ext_ti_ana`, `ext_ti_merge`,
  `frameout`, `nhoparam` — whose parsers became `Arguments` in the move. Those 7 now reject an
  unknown `@key`, report a bad value instead of silently substituting a default, and answer
  `--help` with their registered usage; none of them calls `process::exit` any more. The 18
  gromos++ reference tests pass untouched, and `tests/library_parity.rs` holds the library run to
  the binary's stdout for all 25. **The other 57 analysis binaries still have their compute in
  `main()`.**
- **Every gromos-tools program's body moved into the crate library**
  (`src/programs/<category>/<name>.rs`: `pub fn run(argv, out)` plus a `pub const PROGRAM`); each
  `src/bin/<category>/<name>.rs` file is now a one-line `run_main(&PROGRAM)`. The two trees carry
  the same four categories — `topology`, `box`, `conversion`, `utilities` — and
  `program_and_bin_agree_on_category` fails if a body and its shim are filed under different ones.
  `programs/mod.rs` re-exports every program flat, so a caller still writes `programs::make_top`.
  A program no longer reads the process arguments, writes to the process stdout, or exits —
  `just lint` greps for all three (G11).
- **Hand-rolled `env::args()` parsers replaced by `gromos_io::args::Arguments`** in the 22 tools
  that still had them, so they accept `@f` argument files like the rest (T1 for gromos-tools). **Unknown `@key`s are now rejected** where those parsers silently ignored them.
- `make_top` checked for a missing solvent building block *half-way through printing a topology*
  and exited; the check now runs before any output, so a failure never yields a truncated file.
- `ion` took the solvent molecule size from a hard-coded 3 (SPC water); it now reads it from the
  topology's SOLVENTATOM block.
- `Topology` keeps its parsed blocks beside the built form, and `Configuration` keeps residue and
  atom names — the preparation tools work at that level.

### Fixed

- **`atominfo` was built by two crates** (a stub in gromos-tools, the real gromos++ port in
  gromos-analysis), which made `cargo test --workspace` alternate between passing and failing to
  build. The stub is gone; `atominfo` is one program (2.4).

### Tests

- `crates/gromos-tools/tests/library_parity.rs` — for all 33 programs, `run(argv, &mut Vec<u8>)`
  from the library equals the compiled binary's stdout **byte for byte**, and fails exactly when
  the binary exits non-zero, on real input, on `--help` and on the failure paths.
- `crates/gromos-tools/tests/typed_verbs.rs` — each typed verb's result, serialised with the
  program's own writer, equals what `run` writes; `make_top`'s topology reads back from its own
  text; a failed `make_top` writes nothing; `check_top`'s report survives its failure.
- `py-gromos/tests/test_tools.py` — 26 tests carrying those guarantees across the Python boundary,
  including the preparation chain with no files between steps.
- The 13 gromos++ reference tests are untouched and still pass.

## [0.0.48] (2026-09-04)

**LiveCoMS tutorial 2 (double decoupling: ASA in water and bound to PLA2) started.** Its
equilibration input runs on the tutorial's own files with step-0 energies identical to md++
2023-04-15, and a 20-step run tracks md++'s trajectory to 4×10⁻¹⁶ nm. Two things it took: the
charge-group gathering below, and modelling `ROTTRANS`. The TI production runs still need
`PRECALCLAM`.

### Fixed

- **Charge groups are gathered on reading a configuration** (`gromos_core::gather::
  gather_chargegroups`), as gromosXX does in `Configuration::init`. Every GROMOS tool writes a
  configuration with the molecules put *back into the box*, so a charge group straddles the
  boundary; its centre of geometry — which the charge-group pairlist uses as the cutoff unit —
  then lands in empty space and the wrong neighbours are classified as far away. On the LiveCoMS
  tutorial-2 system that put 189 486 kJ/mol of spurious Lennard-Jones energy into the long-range
  solute set and SHAKE failed at step 0. With the gather, the step-0 energies equal md++'s exactly
  and a 20-step run tracks its trajectory to 4×10⁻¹⁶ nm. No reference system had a split charge
  group, which is why nothing caught it; a twin-range chargegroup cutoff (`RCUTP` < `RCUTL`) is
  what makes it fatal rather than merely wrong.

### Added

- **`ROTTRANS` — roto-translational constraints** (`RottransConstraints`, `rottrans.cc`): the six
  rigid-body degrees of freedom of atoms `1..=RTCLAST` constrained after SHAKE, so a solute
  neither drifts nor tumbles once its position restraints are released. `RTC`/`RTCLAST` are read
  and written; the algorithm is in the plan as `rottrans`. It is *not* bit-exact against md++
  yet: the step-0 potential energy matches to the printed precision and the kinetic energy to
  1.6×10⁻⁶ relative, but the reference orientation `init` stores differs slightly. `NTIRTC = 0`
  (continue from the configuration's `ROTTRANSREFPOS`) is unsupported — the block is neither read
  nor written.
- **SHAKE names the constraint it could not satisfy**: `SHAKE failed to converge after 1000
  iterations (atom 1 - atom 2: 2.5806 nm, wants 0.1250 nm)` instead of the iteration count alone
  (`ConstraintResult::worst`). This is what made the charge-group bug findable in one run.
- **`gromos.EnergyTrajectory`** — a `.tre`/`.trg` read from Python through the energy library,
  keeping the shape the library declares rather than only the scalar totals `Simulation.run()`
  reports: `table("NONBONDED")` is `(n_frames, n_pairs, 6)`, `group_pairs` names its rows as the
  energy-group pairs in gromosXX's order, `to_long()` is the tidy form a `group_by` wants, and
  `warnings` says which tier established the layout. The per-pair Lennard-Jones column sums to
  `totlj`, which is the re-partition workflow the library exists for (PLAN.md 2.10).
- **`.trg` files are read through the library too**: `read_free_energy_trajectory` binds the
  built-in one and returns a `FreeEnergyTrajectory` (frames plus layout warnings), so `ext_ti_ana`
  and the FEP reference harness get the tier checks — a `.trg` one value short of the layout it
  claims is an error now, not dH/dλ from the wrong slot. `dhdl_special` is read (it was hard-wired
  to zero). The pre-0.0.33 one-line `FREEENERGY03` block keeps its own path, without checks.

## [0.0.47] (2026-09-03)

**Energy library profile, version 1, implemented** — `docs/src/reference/energy-library.md`, all
of it. Every `.tre`/`.trg` gromos-rs writes now carries a fingerprint of the layout it was written
with, and `ene_ana` refuses a library that disagrees with the file it is applied to instead of
reading every value into the wrong slot. gromos++ and gromosXX read the files unchanged.

### Fixed

- **`.tre`/`.trg` totals written two slots early.** `EnergyWriter` put the special-interaction
  total at `ENER[19]`, SASA at `[20]`, the constraint total at `[22]` and the distance-restraint
  total at `[23]`; the 2023-04-15 layout has them at `[21]`, `[22]`, `[24]` and `[25]` (`[19]` is
  the lattice-sum surface term, `[20]` the self term; `out_configuration.cc:3074-3125`).
  `FreeEnergyWriter` had the same offset for `dhdl_special`. Any `ene_ana` run on a gromos-rs
  file read `totspecial`, `totdisres` and `totconstraint` as zero before. Found while pinning
  the first-frame identity check to gromosXX's slot numbers. `EnergyFrame::special` now
  means what gromosXX's `special_total` means — every special term, restraints included — and is
  the value `ENER[1] = ENER[2] + ENER[3] + ENER[21]` needs.

### Added

- **Layout in code, fingerprinted** (`gromos_io::energy_traj::Schema`). The ENERTRJ and FRENERTRJ
  layouts of `ENEVERSION 2023-04-15` are declared once, in the same grammar as a library; their
  canonical form drops names and numbers the sizes by first occurrence, so a rename keeps the
  fingerprint while a reshape, a shift or a dropped size changes it. The fingerprint of the layout
  in code equals the one computed from md++'s own `data/ene_ana.md++.lib` (tested) and from
  `sha256sum` of the canonical text by hand.
- **Self-describing trajectories.** `md` writes `# energy-schema 1 ENERTRJ sha256:…` and one
  `# energy-layout` line per declaration after `ENEVERSION`, then `# energy-provenance` lines with
  the writer version, the topology's digest and path, and the energy-group ranges. Both upstream
  readers skip `#` lines at column 0, so gromos++'s `ene_ana` output is byte-identical with and
  without them. `md` prints the layout fingerprint in its parameter summary.
- **Reader tiers** (`EnergyTraj::bind`). A self-described file is the authority for its layout
  and a library that does not match it is refused with a per-declaration diff naming both sides
  (`ENERGY03 ENER:   library 50 x 1,   file 52 x 1`); a file with a known `ENEVERSION` but no
  self-description is checked against the built-in layout; anything else is read as the library
  declares, with a warning that says so. The LiveCoMS tutorial's stale `ene_ana.md++.lib` is
  refused at open under both tiers. Structural checks in every tier: a `size` that is not an
  integer, changes between frames or exceeds what is left of the block; a subblock that overruns
  its block; a gromosXX marker comment (`# baths`, …) falling inside a table; and, as a warning on
  the first frame, the total-energy identity, a negative kinetic energy and a non-positive volume.
- **The official library ships in `gromos-io`** (`crates/gromos-io/data/ene_ana.md++.lib`,
  `energy_traj::OFFICIAL_LIBRARY`): generated from the layout in code with md++'s `VARIABLES`
  verbatim, carrying `# energy-schema` lines that let a reader detect an edited schema section —
  editing a `subblock` line makes the file refuse until the line is regenerated. A test keeps the
  checked-in copy byte-equal to the generator (`UPDATE_OFFICIAL_LIBRARY=1` rewrites it). gromos++
  produces identical output with it and with md++'s own. A library with only `VARIABLES` is valid.
- **`ene_ana` without `@library`** uses the built-in one; `@print_library` prints it, and
  `@print_library <trajectory>` generates a library from the file's self-description (or the
  built-in layout for its `ENEVERSION`) for a reader that has none. Warnings — version mismatch
  in gromos++'s words, unverified layout, first-frame identity — go to stderr. The gromos++
  reference cases run against the built-in library now, not the md++ checkout's copy.

## [0.0.46] (2026-08-31)

**LiveCoMS tutorial 1 (GB3, S² order-parameter restraining), end to end on its own files.** The
whole chain up to the production run works from the tutorial's own `.arg` and `.imd` files:
`make_top` → `com_top` → `pdb2g96` → vacuum `ENERGYMIN` → `sim_box` → solvent `ENERGYMIN` with
position restraints → `ion` → `mk_script` → thermalisation → production MD, plus `rmsd`, `rmsf` and
`ene_ana` on the result. `ORDERPARAMRES` — the restraint the tutorial exists to demonstrate — is
still unported and refused by name.

### Fixed

- **`com_top` wrote force-field type blocks in internal units.** `BONDANGLEBENDTYPE`,
  `IMPDIHEDRALTYPE` and `TORSDIHEDRALTYPE` came out in radians and kJ·mol⁻¹·rad⁻² instead of the
  file convention (degrees, per degree²), and lossily — the dihedral phase π was written `3.1`.
  Every topology `com_top` produced therefore had wrong angle, improper and dihedral parameters:
  the GB3 + 2 Na⁺ system started at an improper energy of 1.3×10⁷ kJ/mol and SHAKE failed at step 0.
  `com_top` had hand-rolled its own 340-line writer; it now builds a `ParsedTopology` and uses
  `write_parsed_topology`, which converts. Its `CombinedTopology` mirror struct is gone.
- **The topology reference tests compared angle/dihedral/improper parameter *counts*, never their
  values** — which is why the above survived. They compare the values now, and fail on the old
  `com_top` (`theta0` 0.0274 rad where gromos++ has 1.5708).
- **`mk_script` silently dropped every `MULTIBATH` override.** A data line that is a bare integer —
  `MULTIBATH`'s `ALGORITHM`, `NBATHS` and `DOFSET` each sit alone on one — was taken for a block
  name and reset the data-line counter, so `TEMP0[i]` never reached the temperature line. The GB3
  thermalisation ran all five of its jobs at 60 K instead of 60 → 120 → 180 → 240 → 300 K. A block
  name is now only recognised outside a block. `TAU[i]` is supported too, as in gromos++.
- **The `.cnf` writer numbered residues wrongly.** It wrote `residue_nr + 1` — gromosXX stores the
  number 0-based and adds one on output, ours is already 1-based — so GB3's first residue came out
  as `2 MET`; and every solvent atom carried its index *within* its molecule (2, 3, 4 repeating)
  where gromosXX writes the solvent molecule's number. The column layout and the
  `# first 24 chars ignored` marker now match `out_configuration.cc:982` exactly.
- **`@f` argument files kept inline comments.** gromos++ erases everything from the first `#` on a
  line (`Arguments.cc:170`); we only skipped whole-line comments, so `@bin /path/md # for MPI …`
  passed `#`, `for` and `MPI` through as values.

### Added

- **Transparent gzip on every text input** (`gromos_io::gz`). gromosXX reads its inputs through
  `igzstream` and gromos++ through `Ginstream`, both on zlib's `gzopen`, which reads plain files
  too — so the decision is made on the gzip magic bytes, not the file name. Topologies,
  coordinates, trajectories, energy trajectories, IMD files, position/distance restraints, PDB
  files, job lists and script templates all go through it. The trajectory reader's optional-block
  lookahead used `seek`, which a gzip stream cannot do; it is a one-line pushback now.
- **`mk_script` compresses what it produced**, as gromos++ does (`mk_script.cc:3853`), and picks
  its output flags from the `WRITETRAJ` switches instead of always passing `@trc @tre`. Together
  with the reader change this is why every analysis input in the GROMOS tutorials is a `.gz`.
- **`rmsd` and `rmsf` on the gromos++ command line**, byte-identical to it. Two independent atom
  sets (`@atomsrmsd`/`@atomsrmsf` for the measurement, `@atomsfit` for the superposition, both
  shifted by the *fit* set's centre of geometry), `@ref` as a coordinate file falling back to the
  trajectory's first frame, several `@traj`, `@pbc <type> [gather]`, `@time`; `@prop` for `rmsd`
  and `@outpdb` for `rmsf`. `@reftopo`/`@refpbc` are refused by name rather than ignored. Five
  reference cases against gromos++.
- **`ene_ana` reads the library that describes the trajectory** (`gromos_io::energy_traj`, a port
  of `utils::EnergyTraj` and `gmath::Expression`): `ENERTRJ`/`FRENERTRJ` declarations with
  `block`/`subblock`/`size` and the `matrix_` expansion, `VARIABLES` arithmetic over the resulting
  tables, and `@en_files`/`@fr_files`/`@prop`/`@library [print]`/`@topo`/`@time` on top. stdout and
  every `<prop>.dat` match gromos++ exactly on five reference cases. gromos++'s operator
  associativity is reproduced deliberately — it splits at the rightmost `+`/`-` but the *leftmost*
  `*`/`/`, so `a / b / c` is `a / (b / c)` — and documented at the module head.
- **`Stat::rmsd_strict`/`msd_strict`**: gromos++ does not clamp `⟨x²⟩ − ⟨x⟩²` at zero, so a constant
  series prints `-nan` where our clamped `rmsd` gives `0`. Programs that reproduce gromos++ output
  use the strict pair, as they already did for `ee_strict`.
- **`sim_box` accepts `@f`**: it parsed `env::args()` by hand and rejected the argument-file form
  every tutorial uses. It now uses the shared `Arguments`, with byte-identical output.
- **Energy library profile, version 1, specified** (`docs/src/reference/energy-library.md`): the
  `ene_ana` library grammar as gromos++ defines it, plus the stricter rules gromos-rs will hold
  itself to — a layout fingerprint derived from the schema (names excluded), trajectories and
  libraries that describe themselves in `#` lines both upstream readers ignore, three reader tiers
  that refuse a library disagreeing with the file it is applied to, and the official library
  generated from the layout in code. Motivated by the tutorial's stale library, which carries the
  same `ENEVERSION` as md++'s current one and reads every value into the wrong slot without an
  error. Specification only in this release; implementation is PLAN.md §2.10.

### Changed

- The `gromos` facade crate builds an `rlib` only. Its `cdylib`/`staticlib` wrote
  `target/libgromos.{so,a}` — the same paths py-gromos's extension module writes — and a
  `cargo test --workspace` raced on them, after which the compiler reported "multiple different
  versions of crate `gromos_core`". The suite alternated between passing and failing to build with
  no code change; it is deterministic now. Nothing consumed a C library from the facade.

## [0.0.45] (2026-08-30)

**The first LiveCoMS tutorial runs on gromos-rs**: tutorial 4 (GaMD), all four of its input files,
unmodified. Its classical force field matches gromosXX exactly at frame 0; the GaMD boost itself is
the remaining difference and now has a reference of its own.

### Added

- **`COVALENTFORM` is modelled** (NTBBH / NTBAH / NTBDN): quartic or harmonic bonds, cosine-harmonic
  or harmonic angles, arbitrary or 0°/180° dihedral phase shifts. All six forms were already
  implemented — nothing selected between them, and the block was refused, which is what blocked the
  GaMD tutorial. New reference `aladip_harmonic_covalent` (NTBBH=1, NTBAH=1) generated with
  gromosXX. Note the dihedral default is now the arbitrary-phase form, as GROMOS specifies; for
  force fields whose phase shifts are 0° or 180° (all the reference systems) the two agree exactly.
- **`INNERLOOP` passes through**: it selects which inner-loop kernel gromosXX uses (special solvent
  loops, CUDA device) and carries no physics, so it no longer stops a run.
- **`aladip_gamd` reference** — aladip + 20 SPC with a GAMD block (dual boost, production mode,
  fixed parameters) and its `@gamd` specification file; the harness and the regeneration script
  learned the flag. The test is `#[ignore]`d with the measured reason: gromosXX accumulates the
  dihedral and total forces of each acceleration group during the force calculation, scales them by
  k·(V−E) and reports ΔV = k(V−E)²/2 in the special energy (≈1.23 kJ/mol on this system); our run
  neither applies that scaling the same way nor reports the boost. The reference is the
  specification for that work (PLAN.md 1.9).

## [0.0.44] (2026-08-30)

Continuing down the LiveCoMS list, with a reference for each gap.

### Fixed

- **A trajectory frame now belongs with the energies of the same frame.** gromosXX writes both at
  the same point of a step, so its `.trc` frame 0 is the input configuration; ours wrote the state
  *after* the step, leaving structure and energy one step out of step with each other. Every
  trajectory frame of every reference system is now bit-identical to gromosXX's.
- **Truncated-octahedron output is written in the input's frame.** The engine works in the rotated
  triclinic frame `truncoct_triclinic_box` produces; gromosXX rotates positions and velocities back
  before writing (`out_configuration.cc`, `truncoct_triclinic_rotmat(false)`) and we only did it for
  forces. `aladip_trunc_oct` now compares coordinates too.
- **Every IMD block is read as a value stream**, not just `NONBONDED` and `CONSTRAINT`: `SYSTEM`,
  `STEP`, `BOUNDCOND`, `MULTIBATH` (the algorithm name and the chain length are values too),
  `INITIALISE`, `WRITETRAJ`, `PRINTOUT`, `COMTRANSROT`, `PAIRLIST`, `FORCE`, `POSITIONRES`,
  `DISTANCERES`, `PERTURBATION`, `ENERGYMIN`. `aladip_multibath_collapsed` now has *every* block on
  one line — gromosXX reads it and gives the same trajectory as the expanded file.
  A malformed `INITIALISE` fixture in `io_integration_tests.rs` (eleven values where GROMOS has ten)
  was corrected: the old line-based reader skipped the extra one and the test passed on a default.
- **`@trv` writes a velocity trajectory** instead of being accepted and ignored (`WRITETRAJ` NTWV).

### Added

- The reference harness compares the **trajectory** frame by frame (positions, images folded) and
  the **velocity trajectory**'s frame count, on top of the energies, forces, free energies and final
  configuration it already compared. Truncated-octahedron images are folded with the cube-frame
  rule, since configurations are written in the frame of the input file.
- `aladip_wrapped` writes a velocity trajectory (NTWV), so `@trv` is exercised by a reference.

### Known

- Which half-step a velocity frame carries is not yet gromosXX's: its frame 0 is the input
  velocities, ours is a half-step further on, and from frame 1 the residual is ~1e-4 nm/ps — the
  size of one thermostat scaling. The harness asserts the frames are written, not their values.

## [0.0.43] (2026-08-30)

The reference suite could not have caught the defects of 0.0.42: all 44 systems use one temperature
bath, one value per line, and molecules that are never wrapped. This release closes that blind spot
and fixes the two further defects the new references exposed.

### Added

- **Three reference systems for the input styles real GROMOS files come in**
  (`crates/gromos-md/tests/gromosXX_references/README-livecoms.md`), each generated with the same
  gromosXX binary as the rest of the suite: `aladip_multibath` (`NBATHS=2` + `DOFSET`, all bonds
  constrained), `aladip_multibath_collapsed` (the same run with `NONBONDED`/`CONSTRAINT` wrapped as
  real files write them — `TOLA2 = 1e-10` on the value stream), and `aladip_wrapped` (the solute
  split across the periodic boundary, 3.80 nm between two of its atoms in a 3.767 nm box, charge
  groups whole). gromosXX gives `aladip_multibath` and `aladip_wrapped` **the same energies** — the
  wrap is pure bookkeeping — which is exactly the property an interaction that skips the minimum
  image breaks.
- **The reference harness now compares the final configuration** (`expected/final.conf`, positions
  and velocities) and requires the frame counts to match. That is where a wrong step count shows —
  the per-step frames look right either way — and it is what let the NSTLIM+1 loop hide. Positions
  are compared modulo a periodic image; NTB = −1 systems skip the coordinate comparison, because
  gromos-rs writes the rotated triclinic frame and gromosXX the truncated-octahedron one.

### Fixed

- **Charge groups are put back into the box each step** (`LatticeShift`, gromosXX's
  `Lattice_Shift_Tracker`: solute groups by their centre of geometry, solvent groups by their first
  atom, immediately before the force field). Without it our coordinates drifted out of the box for
  ever, so the written configuration and the charge-group cutoff saw a different image than
  gromosXX's.
- **A constrained bond is no longer also a force term.** With `NTC ≥ 2` gromosXX moves those bonds
  out of the bond list when it creates the constraints, so they contribute no energy and no force;
  we computed them as well, which on `aladip_multibath` (NTC=3) put 19.36 kJ/mol of bond energy into
  a term gromosXX reports as exactly zero. `Topology::constrained_bonds` now carries the set and the
  quartic and harmonic bond loops skip it.

## [0.0.42] (2026-08-30)

Everything in this release came out of running the **GROMOS LiveCoMS tutorial suite**
(`.local/gromos_tutorial_livecoms`, six tutorials, 21 `.imd` inputs) against gromos-rs for the
first time — real GROMOS input files rather than the ones our reference suite generates. Every
defect below was invisible to the 44 reference systems because they share one style: one
temperature bath, one value per line, and molecules that are never wrapped across the box.

### Fixed

- **The MD loop ran NSTLIM+1 steps.** `start()` is `init` + a full `run_step`, so `for step in
  0..=n_steps` integrated one step too many. The written frames were unaffected (which is why the
  reference tests never caught it), but `@fin` — and therefore any continuation from it — was one
  step ahead of gromosXX. Our final configuration is now bit-identical to gromosXX's for the same
  NSTLIM (checked on `water_216_box`).
- **SHAKE ignored the minimum-image convention** (gromosXX: `shake.h::shake_iteration` takes both
  the current and the reference vector as `periodicity.nearest_image`). A solute whose atoms sit on
  opposite sides of the box — i.e. every equilibrated GROMOS configuration — could not converge:
  the tutorials' aladip has two "bonds" 3.3 nm long as stored.
- **The bonded terms ignored it too** (bonds, angles, dihedrals, impropers and their perturbed and
  soft-core variants): 40 vector sites now go through `nearest_image`, as every gromosXX bonded
  interaction does. On the tutorial system this moved the angle energy from 1183 kJ/mol to
  gromosXX's 16.86 — the whole potential energy now matches gromosXX **exactly** at frame 0.
- **`NONBONDED` and `CONSTRAINT` were parsed line by line**, but gromosXX reads every block as a
  value stream (`Block::get_next_parameter`), so the file may wrap the values anywhere. Real inputs
  put `NLRELE` on the same line as `APPAK RCRF EPSRF NSLFEXCL`, which landed `TOLA2 = 1e-10` in
  `NSLFEXCL`'s integer slot and rejected all 21 tutorial files; `CONSTRAINT` on one line silently
  lost NTCP and both tolerances. Both blocks now read through a `Tokens` cursor, error messages name
  the parameter, and `CONSTRAINT` accepts gromosXX's names and numbers (settle is 4, not 3) and
  refuses the algorithms we do not implement instead of mapping them onto another one.
- **Two crates produced a binary named `mdf`** (gromos-analysis's gromos++ port and gromos-md's
  mean-force utility), so they overwrote each other in `target/<profile>/` and the analysis
  reference test intermittently executed the wrong program. gromos-md's binaries are now declared
  explicitly and its utility is `mean_force`. (`atominfo` is still built by both gromos-analysis and
  gromos-tools — the gromos-tools one is a lesser duplicate of a ported program; left for a
  decision.)

### Added

- **Multi-bath temperature coupling** (`MULTIBATH` with NBATHS > 1, weak coupling): per-bath
  degrees of freedom exactly as `Multibath::calculate_degrees_of_freedom` counts them
  (temperature-group COM dof to the COM bath, the rest to the IR bath, constraints subtracted from
  the bath that owns them, NDFMIN spread in proportion), per-bath kinetic energy, and scaling per
  bath range. On the tutorials this reproduces gromosXX's own DOF (24.99 / 7653.01) and per-bath
  temperatures (184.41 K / 298.52 K). Every reference system uses one bath, and that path is
  unchanged. A DOFSET line with COMBATH ≠ IRBATH, and multi-bath Nosé-Hoover, are refused by name.

### Known, not yet fixed (found by the same run)

- gromosXX puts charge groups back into the box; we never wrap. On a system whose molecules are
  wrapped this still moves the trajectory apart (~3e-2 nm after 10 steps on the tutorial system)
  even though frame 0 matches exactly — the charge-group cutoff is the likely path.
- Our `.trc` writes the state *after* the first step as frame 0; gromosXX writes the initial
  configuration. `@trv` is accepted and then ignored.
- The tutorials' remaining blockers are unported features, refused by name: `ROTTRANS`,
  `PRECALCLAM`, `AEDS`, `REPLICA`, `DISTANCEFIELD`, `ORDERPARAMRES`, and `@qmmm`. `GAMD` is applied
  but unvalidated (gromosXX reports a 45.3 kJ/mol boost where we record none).

## [0.0.41] (2026-08-30)

### Added

- **gromos++ ports, batch 3** (reference tests against gromos++ output): `eds_update_1`,
  `eds_update_2` (EDS parameter updates; shared machinery in `gromos_analysis::eds`), `jepot`
  (J-value local-elevation potential, all three input modes), `pocket` (binding-site vectors,
  enclosed volume and area), `dfgrid` (distance-field grid; `cnf` output and real atoms only) and,
  in gromos-tools, `make_sasa_top`. Every gromos++ program of the analysis/tools suite is now
  ported except `cos_dipole`/`cos_epsilon` (need charge-on-spring special trajectories, which
  gromos-md does not produce), `dGslv_pbsolv` (a Poisson–Boltzmann solver — new physics, not an
  analysis port), `prep_bb` (interactive building-block editor) and `prep_xray_le` (X-ray
  local-elevation setup; no X-ray restraints in gromos-md).
- `gromos_io::jvalue`: JVALRESSPEC, JVALUERESEPS and restraint-trajectory (`.trs`) readers,
  shared by `jval` and `jepot`. `gromos_analysis::time::TimeSpec` (`@timespec`/`@timepts` frame
  selection) shared by `jval` and `jepot`; `VectorSpec` accepts `atom(<one atom>)` as a position.
- Knowingly reproduced gromos++ behaviours: `pocket` never gathers (gromos++ parses `@pbc` and
  does not apply it); `dfgrid` keeps gromos++'s open-set scan so that equal-length paths resolve
  identically. One deviation: `pocket @center atom(…)` is evaluated on the first frame — gromos++
  evaluates it before any frame is read, i.e. at the origin.

### Fixed

- `gromos_io::topology::write_parsed_topology` printed the harmonic angle constant CHT and the
  improper constant CQ in the internal unit (per rad²) instead of per degree² as the file format
  requires, and did not follow gromos++'s `OutTopology` layout (`# n` markers, number formats,
  block comments). `pt_top` output was affected. The writer now reproduces gromos++'s layout and
  is compared token by token with gromos++ in the `make_sasa_top` test.
- `TrajectoryReader` accepts configuration files as trajectories, as gromos++'s `InG96` does:
  frames without a `TIMESTEP` block, labelled `POSITION`/`VELOCITY`/`FORCE` blocks and unknown
  trailing blocks (`STOCHINT`, `PERTDATA`, …).

## [0.0.40] (2026-08-30)

### Added

- **gromos++ ports, batch 2** (reference tests against gromos++ output in
  `crates/gromos-analysis/tests/gromospp_references.rs`): `bilayer_dist`, `bilayer_oparam`,
  `jval` (Karplus ³J from dihedral time series; time series, rmsd and averages modes), `edyn`
  (essential dynamics), `gca` (generate conformations by setting bonds, angles, dihedrals).
  Two gromos++ behaviours are reproduced knowingly: `gca @mobile first` moves the first atom in
  the wrong direction (the formula is gromos++'s; `last`, the default, is right). One is not:
  gromos++'s `edyn` projects onto the *covariance-matrix columns* because `edyn.cc` discards the
  eigenvector matrix returned by `diagonaliseSymmetric` (its `EIVEC.out` equals `COVAR.out`);
  ours projects onto the eigenvectors and the test checks that the projection variances equal the
  eigenvalues. `EIVAL`/`EIFLUC`/`COVAR`/`COVATOM` match gromos++.
- `gromos_io::topology::solvate_to_atoms`: the analysis programs read every atom of a frame, as
  gromos++ `select("ALL")` does, so `s:OW`-style selections work on gromosXX trajectories.

## [0.0.39] (2026-08-30)

### Added

- **gromos++ programs missing from gromos-analysis / gromos-tools, batch 1** — each a port of the
  gromos++ source with a reference test against the gromos++ binary's output
  (`crates/gromos-analysis/tests/gromospp_references.rs`, `crates/gromos-tools/tests/`):
  `mdf` (minimum distances), `dg_ener` (exponential-averaging ΔG), `dfmult` (multi-state ΔF with
  gromos++'s log-exponential estimators), `matrix_overlap` (covariance-matrix overlap),
  `explode`, `duplicate`, `pt_top` (state-B topology; `PERTTOPO` output not implemented).
- **`tser` rewritten to gromos++'s semantics**: property specifiers (`d`, `a`, `t`, `tp`, `o`, `op`,
  `pr`, `pa`), `@time`, `@pbc` gathering, `@dist`/`@norm` distributions, `@skip`/`@stride`; output
  identical to gromos++ (titles, averages, `-nan` where gromos++ prints it). The old `tser` was a
  volume/density toy.
- Shared machinery the ports use instead of re-implementing it (the crate rule): `gromos_io::args`
  (gromos++ `Arguments`), `gromos_io::pbc` (`@pbc` and the `cog` gather), `gromos_io::table`
  (numeric column files), `gromos_analysis::{property, distribution, lnexp, time}`,
  `gromos_io::coordinate::format_g96`, `Stat::ee_strict`/`values`.

### Fixed

- **The trajectory reader mis-read gromosXX's GENBOX block** (NTB line, lengths, angles, Euler
  angles, origin): it took the NTB line as the lengths and desynchronised on the next frame —
  every analysis program failed on a real `.trc`. The writer now emits the same five-line block.

## [0.0.38] (2026-08-30)

### Added

- **MPI in `md`** (`cargo build --release --features use-mpi --bin md`, run under `mpirun`): the
  nonbonded pair terms are decomposed across ranks by first atom index (`Forcefield::pair_partition`,
  `RunOptions::pair_partition`); each rank builds the full pairlist itself and keeps its share, and
  one reduce + broadcast per step (`Forcefield::set_nonbonded_reducer`) gives every rank the same
  forces, energies, virial and dH/dλ, so every rank runs the whole integrator and rank 0 writes the
  files. `np=1` is bit-identical to the serial binary; `np=2/4` reproduce the serial energies to all
  printed digits (forces and positions to 1e-9, the summation order). `scripts/bench_mpi.py` times
  it; BENCHMARKING.md §6 has the table. `crates/gromos-run/tests/pair_partition.rs` checks the
  decomposition without MPI. `AlgorithmSequence::find_mut::<T>()` (via `Algorithm::as_any_mut`)
  reaches the built force field; `ForceStorage` carries `dhdl_lj`/`dhdl_crf`.

### Removed

- The dead MPI code: `md_mpi`, `md_mpi_cuda`, `repex_mpi`, `mpi_scaling` binaries and
  `gromos_integrators::{mpi, remd_mpi}` (an API that no longer existed, an undeclared feature, and
  `&[Vec3]` reinterpreted as `&[f32]`). REMD over MPI, when it comes, builds on the new seam.

### Changed

- The force field adds the non-pairlist terms (excluded-pair/self reaction field, 1-4, the
  perturbed corrections that are not pairlist-based) after the pair terms instead of in between;
  the reference suite (45/45), the check-suite port and the null-perturbation test are unchanged
  within their tolerances (the sums differ in the last bits only).

## [0.0.37] (2026-08-30)

### Fixed

- **`make_top` is now a port of gromos++'s `make_top`** (`utils/make_top.h`): begin/end groups
  replace atoms the way gromos++ does (the old code dropped atoms of the residue after an end
  group and skipped atoms between residues), exclusions come from the building blocks (not from
  connectivity), 1-4 pairs from the bond graph (`get14s`), bonded terms ordered as gromos++'s sets,
  several `@build` files, `OutTopology` block names and layout (`IMPDIHEDRALTYPE`,
  `TORSDIHEDRALTYPE` — the old `…TYPECODE` names made gromos++ reject the file), `CROSSDIHEDRAL*`
  and `LJEXCEPTIONS` blocks. For `NH3+ ALA GLY COO-` and a two-MTB methanol the result reads back
  identical to gromos++'s and gromos++'s reader accepts it.
- **`com_top`** wrote IACs 0-based (solute and solvent) and a `SOLUTEATOM` layout our own reader
  read as "no exclusions"; block names as above. Reads back identical to gromos++'s output.
- **The topology reader read `SOLUTEATOM` line-wise**: exclusions beyond the first line (gromos++
  wraps six per line — every protein) were dropped, and a layout with the counts on their own lines
  lost all exclusions silently. It now reads the block as the token stream gromosXX reads.
  `SOLUTEMOLECULES` with several entries per line (gromos++'s layout) was rejected.
- **The trajectory reader rejected multi-line `TITLE` blocks** (gromosXX's own `final configuration`
  files), which is why `copy_box` and `inbox` crashed on them.
- **`red_top`** wrote a comment list instead of a topology; **`copy_box`** dropped the solvent and
  wrote the box inside `POSITION`; **`inbox`** wrote nothing. All three rewritten to gromos++'s
  semantics (`@atoms 1:1-6`, `@dir x:2,y:2`, nearest image to the box centre) with `GENBOX` output.
- **`ion`** chose a different water than gromos++: the potential now includes every atom of every
  other molecule within the cutoff, Coulomb shifted by −1/R_c, as `utils::Energy` computes it.
- **`make_pt_top`** wrote a `PERTURBEDATOM` block no reader knows; now `PERTATOMPARAM` (+ the empty
  bonded blocks) as gromos++ writes.
- **`atominfo`** printed 0-based IACs; now gromos++'s columns.

### Added

- `gromos_io::topology::write_parsed_topology`: a complete writer in gromos++ `OutTopology` layout —
  what it writes reads back identical (used by `red_top`).
- Reference tests against gromos++ output for `make_top` (2), `com_top`, `red_top`, `make_pt_top`,
  `copy_box`, `inbox`, `ion` (`crates/gromos-tools/tests/`), fixtures with the exact gromos++
  commands in their READMEs. The tools audit table lives in `crates/gromos-tools/.claude/CONTEXT.md`.

## [0.0.36] (2026-08-30)

### Fixed

- **`sim_box` solvated differently from gromos++** (323 vs 352 waters, 2.148 vs 2.208 nm box for the
  same methanol and `@minwall 1.0`). Four deviations: the box was the per-axis extent instead of the
  longest solute atom–atom distance; the template was replicated `ceil(box/template)+1` times instead
  of `int(...)+1`, which shifts the lattice; the template was not centred on its centre of geometry;
  and solute hydrogens counted in the clash test (gromos++ ignores atoms of mass 1.008). Now the same
  algorithm, with a GENBOX block on output. `gromos-tools/tests/sim_box_reference.rs` compares box,
  count and every position with gromos++'s output (`tests/data/sim_box/`).
- `crates/gromos/tests/io_integration_tests.rs`: the MULTIBATH fixture had an invented DOFSET layout
  that gromosXX would reject; the energy-writer test still expected the pre-0.0.34 `ENERTRJ` block.

### Changed

- **Code quality (PLAN.md 2.4):** `cargo clippy --workspace --all-targets` is clean (294 warnings
  before). The mechanical fixes were checked bit-for-bit on seven reference systems (all output files
  byte-identical). `too_many_arguments` is allowed for `gromos-forces` and `needless_range_loop` for
  the numeric crates, each with its reason next to the allow. No bare `unwrap()` is left in the
  non-test code of the library crates: user-input parses now return errors (`local_elevation`
  potential files, MOPAC input writing) and the infallible cases say why in an `expect`.
- 52 new unit tests: the SHAKE / SETTLE / LINCS algorithm wrappers, improper dihedrals (planar zero,
  finite-difference forces), temperature, Berendsen thermostat and barostat, COM-motion removal,
  pressure, energy totals, steepest descent, excluded-pair reaction field, CRF constants, IMD
  write→read, posres / distanceres / `.trg` round trips, and `gromos-run`'s degrees of freedom,
  constraint selection, bundle and prepare. Every `src` file of `gromos-integrators`, `gromos-forces`,
  `gromos-io` and `gromos-run` now carries at least one test.
- `EdsBlock` / `GamdBlock` / `ReplicaBlock`: the inherent `to_string` became `render` plus a `Display`
  impl (`.to_string()` keeps working).

## [0.0.35] (2026-08-30)

### Added

- **gromosXX's own regression oracle, ported** — `crates/gromos-md/tests/gromos_check_suite.rs`
  evaluates the per-term energies hard-coded in `md++/src/check/{aladip,aladip_unperturbed,
  aladip_special}.t.cc` (quartic bonds, angles, impropers, dihedrals, their perturbed and soft-core
  variants, the nonbonded term with both reaction-field variants and the atomic cutoff, distance
  restraints, and finite-difference λ-derivatives) on the suite's own inputs (`shared/check/`,
  byte-identical to gromosXX's). Every value is checked against the *current* gromosXX binary
  (gromos-rs matches to ten digits) and against the suite's table with the suite's own δ. The
  dihedral-restraint value is `#[ignore]`d until PLAN.md 1.6 lands.

### Fixed

- **NSLFEXCL was parsed and never consumed.** The reaction-field term of the excluded pairs and the
  self term were always added; gromosXX skips them when NSLFEXCL = 0 (`nonbonded.rf_excluded`,
  also the excluded-state term of perturbed atom pairs). `ForcefieldPlan::rf_excluded` now gates all
  three places. The default is 1, as in gromosXX (`parameter.h`); every reference input already
  says 1, so nothing else moves.
- The MULTIBATH reader assumed one line of `TEMP0 TAU` pairs; gromosXX reads the block as a token
  stream and its own check inputs list one bath per line. Now token-based — both layouts read.

## [0.0.34] (2026-08-30)

### Added

- **`.tre` in gromosXX's native layout** (`TIMESTEP` + `ENERGY03` + `VOLUMEPRESSURE03`; per-bath
  kinetic and per-group LJ/CRF from the engine's own accumulators, the untracked splits written as
  zeros and documented), and frames on gromosXX's schedule (steps 0 … NSTLIM−1). gromos++ `ene_ana`
  with the current `ene_ana.md++.lib` (`.local/gromosXX/md++/data/`) now reads gromos-rs energies
  and free energies and reports the same `totene`/`totkin`/`totpot`/`totlj`/`totcrf`/`dvdl` as for
  gromosXX's files of the same window — the LiveCoMS analysis scripts run unchanged on gromos-rs
  output. `gromos_io::read_energy_trajectory_native` reads the layout; the reference and pairlist
  tests use it.
- **dH/dλ from Python**: `Simulation.dhdl` (kJ/mol, the TI integrand), `Simulation.dhdl_terms`
  (`{"lj", "crf", "bonded"}`), and `dhdl` as the 13th column of `run()` / `EnergyTimeseries`.
- `md` writes the converged energies once more at EM convergence, as gromosXX does, so
  `aladip_vacuum_em` — ignored since the reference was made ("frame count off-by-one") — passes.
  **The Rust reference suite is 45 of 45.**
- `tests/fep_null_perturbation.rs`: a `.ptp` with identical end states and α = 0 must leave energies
  and forces bit-identical to the unperturbed run, and α ≠ 0 must not — the invariance that pinned
  the FEP defects, now a permanent test.

### Fixed — every FEP reference now matches gromosXX

Found by bisecting `aladip_vacuum_fep` against the native gromosXX one perturbation block, atom
and property at a time, ending with *null perturbations* (identical A and B states, α = 0) — a
non-zero Δ there is a bookkeeping bug by construction. `aladip_vacuum_fep`, ignored since 1.7,
passes; so do `meoh_water_fep` and CH4 at five λ. Rust suite 44 passed / 1 ignored (`aladip_vacuum_em`,
an EM frame-count question, not physics); Python suite 344 passed, no xfail left.

- **Perturbed nonbonded electrostatics, two defects found by bisecting `aladip_vacuum_fep` against
  gromosXX with single-block `.ptp` files.** (1) The soft-core pair kernel's reaction-field term
  used `crf/(2·R_c³)` where gromosXX uses `crf/2` before dividing by the softened cutoff³ — an
  extra 1/R_c³ on the r² term of every perturbed pair, invisible with zero charges (CH4) and
  0.16 kJ/mol on methanol. (2) The pairlist correction only soft-cored pairs whose *lower-index*
  atom was perturbed; gromosXX routes a pair to the perturbed list when either atom is. A
  perturbed atom that is not at the start of the topology lost most of its soft-core pairs.
  `meoh_water_fep` now passes (energies, positions, forces, dH/dλ per term, ten steps);
  every single-atom null/charge/type variant of the dipeptide is exact.
- **λ-mixed masses were never applied** in the `md` path (`mass_at_lambda` had no caller):
  `prepare_system` now sets m(λ) = (1−λ)·m_A + λ·m_B for perturbed atoms, as gromosXX's
  `update_for_lambda` does.
- **`PERTATOMPAIR` entries were ignored unless some atom was also perturbed** (the whole perturbed
  nonbonded block was guarded by "any perturbed atom"), their CRF part was never applied, and an
  "excluded" end state dropped the reaction-field term an excluded pair keeps in gromosXX.
- With these, **`aladip_vacuum_fep` passes** — the reference that had been `ignore`d since 1.7 —
  and every FEP reference (CH4 at five λ, methanol, the dipeptide) matches gromosXX. Suite: 44
  passed, 1 ignored (`aladip_vacuum_em`, an EM frame-count question, not physics).

- Perturbed bonded energies (bonds, angles, impropers, dihedrals, and their soft variants) are
  booked in their own `.tre` columns instead of all in `bond_total` — totals were right, columns
  were not (`ForceEnergyLambda` carries the split). In the `aladip_vacuum_fep` bisection the
  angle-, improper- and dihedral-only cases are now exact against gromosXX.
- A perturbation topology whose only content is soft-core bonded terms (`PERTBONDSOFT`,
  `PERTANGLESOFT`, `PERTIMPROPERDIHSOFT`) or perturbed atom pairs counted as *empty*
  (`PerturbedSolute::is_empty` ignored those lists): the regular term was removed from the topology
  and nothing replaced it. Every soft bonded case of the `aladip_vacuum_fep` bisection is now exact
  against gromosXX (the soft-core formulas themselves were already right).

## [0.0.33] (2026-08-29)

### Fixed

- **dH/dλ was never compared against gromosXX.** The reference suite (and `gromos-io`'s `.trg`
  reader, hence `ext_ti_ana`) looked for a one-line `FREEENERGY03` block, while the native binary
  writes `TIMESTEP` + `FREEENERDERIVS03` (`# lambda`, `# totals` in ENERGY03 order, per-bath and
  per-group sections) — so every "dH/dλ tracked" claim since 1.7 rested on zero frames. The reader
  now parses both layouts, the reference test compares the total and the LJ / CRF / bonded parts
  per frame (aligned by time), and the writer emits the native layout with the run's bath and
  energy-group counts (per-group derivatives as zeros).
- The perturbed nonbonded derivative is split into its LJ and CRF parts
  (`PertNBCorrection::dhdl_lj/dhdl_crf`, `Energy::dhdl_lj/dhdl_crf/dhdl_bonded`), so the `.trg`
  carries the same per-term derivatives gromosXX writes.

### Added

- **CH4 λ-range references** `ch4_water_fep_l000/l025/l075/l100` (ten steps each, NTWG=1): the
  proven CH4 → dummy system at λ = 0, 0.25, 0.75, 1 — energies, positions, forces and now dH/dλ
  (total, LJ, CRF, bonded) match gromosXX at every λ. Reference suite: 42 pass, 3 ignored.
- `scripts/ti_ch4.py`: thermodynamic integration of CH4 → dummy through gromosXX and gromos-rs
  from the same inputs (11 windows, NVT 300 K), integrated with the project's `ext_ti_ana`;
  per-window ⟨∂H/∂λ⟩ ± ee (n) for both engines, ΔG and timings, into `bench/work/ti_ch4/`.
- `scripts/regen_gromosXX_references.py <system …>` regenerates only the named systems.

## [0.0.32] (2026-08-29)

### Added (PLAN.md 3.9 step 5 — the first new term through the new door)

- **`meoh_water_fep` reference system** (gromosXX native build, ten steps at λ = 0.5): 54a7 united-atom
  methanol (charged) → dummies in 953 SPC waters, built with gromos++ `make_top` and the project's
  `sim_box`, minimised and equilibrated with gromosXX. It exercises the perturbed reaction-field
  self/excluded-pair terms that the zero-charge `ch4_water_fep` never could — and fails today (CRF off
  by 0.16 kJ/mol at frame 0, LJ exact): registered as `ignore` in the Rust suite and a strict xfail in
  the Python suite until the term is second-sourced (PLAN.md 1.7). The `aladip_vacuum_fep` mismatch was
  bisected block by block against gromosXX; findings in PLAN.md 1.7.

- **`Term("xtb", region=..., elements=[...], gfn=2, charge=0, multiplicity=1, work_dir=None,
  timeout_s=600, coupling="delta")`** — a real GFN-xTB subprocess (`XtbInteraction`) over a region,
  additive on top of the classical force field, no cargo feature. Wiring: the `TermSpec::Xtb`
  variant (`recipe.rs`), its registry lines (`plan.rs`) and one `instantiate` arm (`build.rs`);
  `gromos.terms()` lists it, `_TermKind` in the stubs names it. Each xtb term gets its own
  `work_dir` (`<tmp>/gromos-rs-xtb-term-<index>` unless given) and every xtb call is bounded by
  `timeout_s` (`XtbInteraction::with_timeout`, `qm_subprocess::run_subprocess_with_timeout`: the
  child is killed and the step fails with a clear error instead of hanging).
- **Per-term energies (G10, in memory).** `Energy.term_energies: Vec<(String, f64)>` is filled by
  the orchestrator every step (`ProviderOrchestrator::register_labelled` / `evaluate_with_terms`),
  keyed by the term's registry name (`xtb`, or `xtb:0`/`xtb:1` for a repeated kind), in plan
  order; `Simulation.term_energies -> dict[str, float]`. Their sum is what the terms add to
  `total_energy`; `potential_energy` stays the classical force field alone. Not yet in the `.tre`
  file or in `run()`'s columns.
- Documentation pass for the recipe model: root `README.md` (components, "How a run is described"),
  `py-gromos/README.md` (rewritten), a new concept page `py-gromos/docs/user-guide/recipe.md`
  (in the mkdocs nav), the API reference without the deleted building blocks, `gromos-run`'s crate
  doc, the `md`/`pyo3-gromos`/`parameters` module docs, `.claude/overview.md`/`architecture.md`,
  `FUTURE.md`'s 2026-07 addendum marked resolved, `BENCHMARKING.md`'s parity note, and the design
  mockup notebook marked superseded.
- `py-gromos/tests/test_xtb_term.py` — the physics oracle from the Python side: adding the term
  adds exactly `XtbPotential`'s direct energy and forces of the region (and nothing outside it),
  two terms report separately and add up, NVE through the real step loop holds the Rust oracle's
  thresholds, `coupling="replace"` / a barostat / an `elements` table that does not cover the
  region / a `timeout_s` that expires are each a named error. Skips without `xtb` on PATH.

## [0.0.31] (2026-08-29)

### Removed (PLAN.md 3.9 step 4 — delete the copies)

- The Python **descriptor path**: `resolve_algorithm_sequence` (the second IMD→sequence builder),
  `AlgorithmDescriptor`, and the twelve building-block classes (`Forcefield`, `LeapFrogIntegrator`,
  `LeapFrogVelocity`, `LeapFrogPosition`, `BerendsenThermostat`, `BerendsenBarostat`,
  `ShakeConstraints`, `SteepestDescent`, `TemperatureCalculation`, `PressureCalculation`,
  `EnergyCalculation`, `RemoveCOMMotion`). The MD step is `gromos.Plan` from `recipe.plan(system)`;
  `crates/pyo3-gromos/src/algorithm_sequence.rs` went from 1473 to ~100 lines and builds nothing.
- `gromos.md_runners` (subprocess wrappers around the `md` binary) — `Simulation` is the way to run
  MD from Python. `gromos.analysis` is no longer in the default namespace (`import gromos.analysis`
  still works). The placeholder suites `test_basic.py` / `test_advanced_features.py` (eleven skipped
  classes for a past API) — `test_recipe.py` covers the surface they described.
- `py-gromos`'s `glam` dependency.

### Changed

- `AlgorithmSequence.nve/nvt/npt/minimize/from_parameters(topo, params)` (deprecated) now return the
  `Plan` of `params` — the ensemble comes from the parameters, not from the preset name — and
  `Simulation.from_sequence(topo, conf, params, plan)` takes that `Plan`; `AlgorithmSequence` itself
  cannot be instantiated. `Simulation.recipe` / `.plan` / `.recipe_toml` / `.plan_json` are always
  present (no `None` case left).
- `Vec3`, `Frame`, `rmsd`, `rdf` are f64 (`gromos_core::math::Vec3`), not f32.
- Atom-count and factory-argument errors are `gromos.exceptions.RecipeError` (still a `ValueError`).
- One IMD reader: `gromos_run::read_imd` — the `md` binary, `Recipe.from_imd`, the deprecated
  `InputParameters` and the bundle loader all open parameter files through it.
- `just lint` runs the G6 drift gates after clippy: `AlgorithmSequence::new()` outside `gromos-run`,
  `.push(Box::new(` in the binding, `read_imd_file(`/`parse_imd_str(` outside `gromos-run`/
  `gromos-io`, `process::exit` in `gromos-run` — each grep must print nothing.
- `test_front_end_parity.py`: path C and its `xfail` table are gone; the deprecated
  `AlgorithmSequence`/`from_sequence` names are checked as a translation (`array_equal` against the
  `Recipe` path) instead.

## [0.0.30] (2026-08-29)

### Added

- **`gromos.Recipe`, `Term`, `Algorithm`, `Plan` — the run as data, in Python (PLAN.md 3.9
  step 3).** They *are* the `gromos-run` types (serde ↔ Python through `pythonize`; no Python-side
  copy to drift). `Recipe.from_imd/from_toml/from_json/from_dict/from_bundle`, `Recipe.nve/nvt/npt/
  minimize`, one getter per group (`.control`, `.ensemble`, `.constraints`, …), immutable
  `update(**groups)` (deep merge), `with_term/without_terms/with_inputs/with_execution`,
  `to_imd/save_imd/to_bundle`, `==`, pickling. `recipe.plan(system)` is the validated MD step as an
  editable `Plan` (addressed by index or kind; `insert/insert_after/insert_before/remove/replace/
  validate`; JSON) and `Simulation(system, recipe, plan=plan)` runs it — stage 1 of the same builder,
  not a second one. `Simulation.recipe`/`.plan` objects, `Simulation.from_bundle`. Registries:
  `gromos.terms()`, `gromos.algorithms()` (with the ordering rules), `gromos.build_info()`.
- `gromos.exceptions`: `RecipeError`, `PlanError` (⊂ `ValueError`), `MissingFeatureError`,
  `RunError` (⊂ `RuntimeError`). A typo'd field or an unknown kind is a `RecipeError`, never a
  silent default. `Term("schnet", …)` constructs on a non-`ml` build with `available=False`;
  `Simulation` then raises `MissingFeatureError` (A14).
- Drift guards (G2/G5/G6): `test_front_end_parity.py` compares the `Recipe` front-end and the `Plan`
  front-end against the deprecated `InputParameters` path with `np.array_equal` on every reference
  system (no divergences); `test_every_kind_has_a_parity_case` / `test_pyi_lists_every_kind` fail
  when an algorithm kind is added without a reference case or a stub; `mypy.stubtest` is clean
  against `gromos.pyi` (`MYPYPATH=python python -m mypy.stubtest gromos.gromos --allowlist
  stubtest_allowlist_no_ml.txt`; without the allowlist on an `ml` build).

### Changed

- `test_gromosXX_references.py` builds every reference run as `Simulation(system,
  Recipe.from_imd(...).with_inputs(...))`: the recipe front-end is the gromosXX-guarded one.
- `System(topology, configuration)` / `System.from_files` accept a solute topology with a solvated
  coordinate file (solute plus whole solvent molecules); `prepare_system` solvates from the
  coordinate count, as the `md` binary does.
- `Recipe.__repr__` says NVE when no bath couples (TAU < 0 is gromosXX's "no coupling"; the bath is
  kept so the `.imd` round-trips).
- `gromos.__version__` comes from the installed distribution metadata; notebooks `01`/`02` and the
  Python docs use `Recipe`.

### Deprecated (one release; each emits a `DeprecationWarning` naming its replacement)

- `InputParameters(path)` / `.from_file` / `.nve` / `.nvt` / `.npt` / `.steepest_descent` →
  `Recipe.from_imd` / `Recipe.nve|nvt|npt|minimize`; `Simulation(system, InputParameters)` →
  `Simulation(system, recipe)`.
- `Simulation(..., distrest=, posresspec=, refpos=)` → `recipe.with_inputs(...)`;
  `Simulation(..., ml_potential=, ml_region=, ml_buffer=)` → `recipe.with_term(Term("schnet", ...))`.
- `AlgorithmSequence.nve|nvt|npt|minimize|from_parameters` and `Simulation.from_sequence` →
  `recipe.plan(system)` + `Simulation(system, recipe, plan=plan)`.

## [0.0.29] (2026-08-29)

### Added

- **`RunRecipe` — the run as plain data (PLAN.md 3.9 step 2).** `gromos-run` now builds every
  run from a versioned, serde-serialisable recipe: `.imd` → `RunRecipe::from_imd` →
  `build_plan` → `Vec<AlgorithmSpec>` (fully resolved) → `validate_plan` → `instantiate` (reads
  only the plan) → `start`. The `md` binary and `Simulation` both go through it; a frozen copy of
  the step-1 builder is the test oracle and the recipe path is bit-identical to it on all 41
  reference systems. `md @dump` prints recipe + plan as JSON; `md` writes `<tre>.recipe.toml`
  with the parser diagnostics; `Simulation.recipe_toml` / `.plan_json` / `.diagnostics`.
- `.imd` **writer** (`gromos_io::imd_write::write_imd`) and a lossless, strict parser: every
  field of every modelled block round-trips exactly (multi-line TITLE, NTIRTC/NTISTI, NTWSE/NTWG/
  NTWB, NTPP, NCYC, NTCS0, DOF sets, PRESSURESCALE COUPLE/SEMI, NONBONDED lattice-sum lines;
  NTWV/NTWF are frequencies, not booleans), and a malformed number is an error naming the block
  and token instead of a silent 0. `parse_imd_str` for in-memory text.
- `scripts/roundtrip_imd_gromosXX.py`: feeds every rewritten `.imd` to the real gromosXX —
  **40/40** accepted and reproducing a fresh run of the original exactly. It caught two writer
  defects a Rust round trip could not (fused columns; an invented "off" DISTANCERES block that
  gromosXX validates).
- `validate_plan`: GROMOS ordering invariants as per-kind rules (one `forcefield`, `orchestrator`
  right after it, `energy_calculation` last, steepest descent excludes the dynamics kinds,
  barostat requires the pressure calculation, a term without a virial cannot meet a barostat,
  `coupling=replace` rejected until the zone-aware Forcefield exists).

### Fixed

- **An `.imd` without a MULTIBATH block no longer runs a Berendsen bath**: the parameter defaults
  now mean "block absent" (no bath; `nscm = 0` for a missing COMTRANSROT, as gromosXX's
  `parameter.h`), so `water_216_nve_nobath` passes in both suites and is no longer ignored.
- More than one temperature bath is a named error (`build_plan`) instead of a silent truncation
  to the first bath.

## [0.0.28] (2026-08-29)

### Added

- **`gromos-run` crate — the one run builder (PLAN.md 3.9 step 1).** The IMD → algorithm-
  sequence assembly that the `md` binary and the Python binding each carried (and had let
  drift) now lives once, as a library: `prepare_system` (perturbation topology, truncated-
  octahedron transform, NSM from the coordinate file, validation, initial velocities),
  `build_sequence_from_imd` (the GROMOS step order), `start` (init + step 0), `total_dof`, and a
  `RunError` in place of the binary's `println!`/`process::exit` sites. `md.rs` shrank from 2060
  to 1359 lines; `pyo3-gromos`'s `build_simulation` contains no algorithm construction at all.
- Divergences resolved toward the binary: Python now honours the topology's PHYSICALCONSTANTS
  (`four_pi_eps_i`), takes NSM from the coordinate file, uses the same parallel-kernel policy
  (`ParallelPolicy::Auto` = parallel above 100 atoms), and one DOF formula that also counts
  solute constraints (the binary's old `TODO`; no reference combines NTC>1 with a live
  thermostat, so no reference result moved).
- FEP from Python: `Topology.apply_perturbation(path)` (the `.ptp` merge is now a
  `Topology` method in gromos-core, shared with `md @pttopo`); `ch4_water_fep` joins the Python
  reference suite and passes; `aladip_vacuum_fep` is a strict xfail like the Rust `ignore`.
- `AlgorithmSequence::insert` (gromos-core) so the ML term can be placed after `Forcefield`
  without a second builder.
- `crates/gromos-run/tests/prepare_and_build.rs`: NSM-from-coordinates, solute-constraint DOF,
  perturbation/NTG consistency, named errors.

### Changed

- The descriptor path (`AlgorithmSequence.from_parameters` → `resolve_algorithm_sequence`) uses
  the same parallel policy and physical constants as the shared builder, so the front-end parity
  test stays exact on >100-atom systems; its other gaps (SETTLE, LINCS, Nosé-Hoover, NTICOM=2,
  truncated octahedron, restraints, FEP) remain strict xfails until the recipe plan replaces it
  (PLAN.md 3.9 steps 2–4).
- Suites after the lift: Rust reference 37 passed / 3 ignored; Python 208 passed, 16 skipped,
  18 xfailed.

## [0.0.27] (2026-08-29)

### Added

- **PLAN.md 3.9 step 0 — front-end parity measured before any refactor.**
  `py-gromos/tests/test_front_end_parity.py`: every reference system through `Simulation(...)` and
  through `AlgorithmSequence.from_parameters(...)`, compared with `np.array_equal` (never a
  tolerance) on per-step energies, final positions and forces, plus a `same_path_twice` determinism
  baseline. Result: 27/37 systems bit-identical between the two builders; the 10 divergences
  (SETTLE, LINCS ×2, Nosé-Hoover ×2, NTICOM=2 rotation removal, truncated-octahedron box,
  restraints ×3) are `xfail(strict=True)` entries naming the missing feature, mirrored in PLAN.md
  3.9's divergence table.
- New gromosXX reference `water_216_nve_nobath` (generated with the CODATA-patched
  `.local/gromosXX/md++/BUILD/program/md`, the only build that reproduces the committed references
  bit-for-bit): `water_216_box` with the MULTIBATH block absent. gromosXX runs no temperature
  coupling; gromos-rs's IMD parser keeps `TempBathParameters::default()` and silently thermostats.
  Registered as an ignored test (Rust) / strict xfail (Python) until PLAN.md 3.9 step 2 makes the
  parser presence-aware.
- `scripts/kernel_determinism.py` — run-to-run and thread-count determinism of the `md` binary
  (n = 3): bit-identical run-to-run at a fixed thread count; 1 vs 8 threads differ at ~1e-13
  relative. Recorded in BENCHMARKING.md ("Kernel determinism").

### Changed

- `test_steepest_descent_via_algorithm_sequence` compares energies and positions between the two
  builders, not only algorithm names.
- PLAN.md: finished items condensed to their load-bearing decisions (full notes moved verbatim to
  the new `PLAN_ARCHIVE.md`); 3.8 audit extended with verified engine divergences; 3.9 rewritten as
  a checkable plan (assumptions A1–A19, drift guards G1–G10, steps 0–5 with gates) after a
  five-perspective review — OpenMM, GROMACS, HOOMD-blue, ASE/i-PI/PLUMED, Polars/pyo3 — archived in
  `PLAN_ARCHIVE.md` §3.9-review.

### Known engine defects recorded (not yet fixed)

- An `.imd` without a MULTIBATH block runs a Berendsen bath at 300 K, τ = 0.1, in both the `md`
  binary and Python; every builder also reads only `temp_bath[0]`, so multi-bath inputs are
  truncated. Fix scheduled in PLAN.md 3.9 step 2 (A18).

## [0.0.26] (2026-08-29)

### Performance

Single-core benchmarking against gromosXX 1.6.0 (both `-march=native`) began with gromos-rs at
0.48× (`water_216_box`), 0.77× (`water_1000_spc_gridcell`) and 0.82× (`ch4_water_fep`) of
gromosXX's speed. After this work it is at 0.97×, 1.11× and 0.90× (mean ratios, n = 3–5), and
**1.14× on a 24 000-atom production-regime box** (NSNB=5, twin-range, grid_cell) built with
`sim_box` via `scripts/make_solvated_box.py`. With 8 threads (Phase 3) it is 1.38× faster than
gromosXX's OpenMP build on that box and 1.06× on the 3 000-atom system; before the
parallel-pairlist work it was 0.53× and 0.43× there. See `BENCHMARKING.md`. Every change preserves the reference-suite
output; the nonbonded kernels are bit-identical per pair to the previous scalar code.

- **gromos-core:** `Rectangular::nearest_image` rounds ties-to-even (`rint`, as gromosXX does)
  instead of ties-away-from-zero — the latter has no native x86 encoding and was ~47 % of the
  nonbonded kernel on `water_216_box`.
- **gromos-core:** cell-list pairlist takes its grid spacing from the IMD `PAIRLIST` `SIZE` field
  (`auto` = half the short cutoff, as gromosXX) instead of a hardcoded cutoff-wide cell, and prunes
  cell pairs by true inter-cell distance instead of a cubic shell. On `water_1000` the old sizing
  produced a 3×3×3 grid where every cell neighboured every other, i.e. no pruning at all.
  Pairlist phase 3.15 s → 2.15 s. `SIZE` had been parsed and never used.
- **gromos-core:** `PairlistContainer` stores `(u32, u32)`; the per-step copy into `u32` in the
  forcefield is gone (2.5 s on `ch4_water_fep`) and pairlist write traffic is halved.
- **gromos-core:** `is_excluded_or_14` / `is_excluded` query one side of the (symmetric, sorted)
  lists with gromosXX's `last < j` early-out; 1-4 pairs are now stored symmetrically by the reader.
  Invariants pinned by `gromos-io/tests/exclusion_invariants.rs`.
- **gromos-forces:** LJ/CRF innerloop evaluates four pairs per `f64x4` register (divide and sqrt
  were the arithmetic cost); the flat, charge-group-grouped and solvent kernels now share one
  `process_pair_slice` — the file's "single source of truth" promise now holds.
- **gromos-forces:** parallel kernels size chunks per thread instead of per 1024 pairs, which was
  allocating and merging a full force buffer 61 times per step on `water_216_box`; also fixes a
  latent `par_chunks(0)` panic when there are fewer charge-group pairs than threads.
- **gromos-core / gromos-integrators:** long-range solvent–solvent pairs are stored as one sentinel
  pair per molecule pair and evaluated with a shared periodic shift when that is provably exact
  for the box (`sentinel_long_range_is_exact`; `PairlistContainer::solvent_long_is_sentinel`),
  as the short-range solvent list always was; otherwise the expanded per-atom-image path is kept.
  Verified against the expanded path by `gromos-forces/tests/solvent_longrange_sentinel.rs`.
- **gromos-integrators:** long-range kernels use the `_novirial` variants when no pressure
  coupling reads the virial, matching the short-range path.
- **gromos-core:** pairlist construction runs in parallel — the cell list over cell ranges, the
  standard O(N²) algorithm over charge-group ranges — into private lists appended in order, so
  pair order and all outputs are unchanged. Per-pair debug logging in the standard algorithm's
  hot loop was removed.
- **gromos-integrators:** charge-group pair-group metadata is built in parallel (per-chunk runs
  stitched at chunk boundaries; identical group boundaries).
- **gromos-integrators:** the long-range block dispatches to the parallel kernels when
  `parallel_nonbonded` is set, as the short-range block already did; both phases previously ran
  serial regardless of thread count.
- **bench:** `scripts/bench_engines.py` reports n, mean, standard deviation, median and min;
  `scripts/make_solvated_box.py` builds large solvated boxes through the project's `sim_box`.
- **gromos-integrators:** twin-range long-range solute interactions use the charge-group-grouped
  kernel (one `nearest_image` per CG pair) under the same safety condition as short-range, and
  reuse one force buffer; 1.59 s → 0.60 s on `water_1000`.
- **gromos-integrators / gromos-core:** per-algorithm and per-force-phase timing is opt-in
  (`AlgorithmSequence::enable_timing`); `md @verb 1` prints a `TIMING` block comparable line by
  line with gromosXX's. Previously two `Instant::now()` calls per algorithm per step ran
  unconditionally.

### Fixed

- **gromos-md:** `md` panicked with "remainder with a divisor of zero" whenever `NTWX`, `NTWE` or
  `NTPR` was 0 (GROMOS: "never write").
- **gromos-core:** `PAIRLIST` `TYPE atomic` was parsed and ignored; it now selects the atom-based
  cutoff. `ALGORITHM grid` fell through to the O(N²) algorithm; it now selects the cell list
  (documented deviation: same pair set, different traversal).
- **gromos-core:** `Topology::solvate` left `one_four_pairs` at the solute length; any per-atom
  query on a solvent atom (reachable via `TYPE atomic`) indexed past the end.
- **gromos-md:** `benches/md_bench.rs`, `scripts/benchmark.sh` and
  `scripts/regen_gromosXX_references.py` pointed at non-existent paths.

### Added

- `BENCHMARKING.md`: the engine-vs-engine protocol, assumptions, and measured results;
  `scripts/bench_engines.py`: the harness that produces them.

## [0.0.25] (2026-08-13)

### Features

- **gromos-core:** `SpatialIndex` trait (`spatial_index.rs`) — query-based neighbor service ("pairs within r of this selection"), independent of the MD pairlist's charge-group/twin-range/solvent-sentinel shape so it can serve arbitrary-cutoff consumers (ML radial graphs, QM-zone gathers) later. `ConfigurationSpatialIndex` is the first (brute-force, atom-level, exclusion-free) implementation; carries periodic-image shift vectors per pair so a later periodic ML/QM provider doesn't force a breaking signature change (FUTURE.md Dim 12.5 point 3 / P3)
- **gromos-forces:** `PotentialProvider` trait (`provider.rs`) — the additive-provider shape for force/energy computation (`architecture.md` "provider pattern"; FUTURE.md Dim 12.5), so classical terms, and eventually QM/MM/ML potentials, are all the same kind of thing. `contribute()` takes read-only `Topology`/`Configuration` plus an arbitrary `AtomSelection` region and a `SpatialIndex`, returns a scattered `Contribution` (energy + sparse per-atom forces + virial + extensible `ProviderExtra`) rather than raw `&mut Configuration` access — reviewed and corrected before landing (see Testing) so a provider can't touch atoms outside what it's scoped to
- **gromos-forces:** `LjCrfInteraction` (`nonbonded/interaction.rs`) — the first `PotentialProvider` impl, wrapping the existing `lj_crf_innerloop_novirial` with zero math changes to prove the shape is transparent. Named `*Interaction` (not `*Provider`) to match gromosXX's own `Nonbonded_Interaction`/`QMMM_Interaction` base-class convention. Scoped to solute-solute atom-level pairs only this pass — solvent's charge-group "sentinel pair per molecule pair" pairlist shape doesn't fit the generic `SpatialIndex` query and needs its own adapter, deferred to the wiring pass
- **gromos-forces:** `SchNetInteraction` (`nonbonded/schnet.rs`, feature-gated `ml`, default builds untouched) — the first ML `PotentialProvider`: loads a TorchScript model via `tch`/libtorch and runs it **in-process**, no Python subprocess or embedded interpreter, closing the loop on the whole reason this groundwork exists. Now loads a real `schnetpack.model.NeuralNetworkPotential` (SchNetPack **2.x**, not a hand-rolled reimplementation — see Notes), calling it via its actual `Dict[str, Tensor]` I/O contract; forces are computed by the model's own internal `Forces` module, no manual `.backward()` needed. Consumes the same `SpatialIndex` classical terms use for its radial graph — `SchNetPack`'s own `PairwiseDistances` expects a host-supplied neighbor list too, confirming FUTURE.md P3 (one spatial service, multiple consumers) against a real, current ML library, not just GROMACS's NNPot interface
- **scripts:** `export_toy_schnet.py` — now builds a real (still randomly initialized/untrained) SchNetPack 2 `SchNet` + `Forces` model and exports it via TorchScript, replacing the earlier hand-rolled single cfconv block

### Testing

- **gromos-core:** `spatial_index::tests` — vacuum and periodic (minimum-image + shift-vector) neighbor queries
- **gromos-forces:** `nonbonded::interaction::tests::wrapper_transparency_matches_direct_call` — `LjCrfInteraction` vs. a direct `lj_crf_innerloop_novirial` call, bit-identical (transparency check, not independent oracle coverage on its own)
- **gromos-forces:** `nonbonded::interaction::tests::cross_checks_against_single_point_energy_oracle` — cross-checked against `single_point_energy`'s already reference-suite-validated LJ+CRF total on a 4-atom fixture (no exclusions/1-4, so the only reconciliation needed is gromosXX's unconditional per-atom RF self-energy term, which this pairlist-only provider deliberately excludes) — a genuine second computation path, since `single_point_energy` builds its pairlist via the existing chargegroup-aware `StandardPairlistAlgorithm`, independent of the new `ConfigurationSpatialIndex`
- **gromos-forces:** `nonbonded::schnet::tests::model_forces_match_finite_differences` (feature `ml`) — there is no GROMOS oracle for a neural-network potential (P8), so this is the tolerance-based second test tier the plan called for: central finite differences on the energy surface vs. the model's own forces, agreeing to `<5e-3` absolute (widened from `1e-3`; the real SchNetPack interaction block has more float32-accumulating layers than the old hand-rolled one — verified per-component that the ~1.4e-3 gap is noise, not a bug, before loosening)
- **gromos-md:** `tests/test_provider_reference.rs::lj_crf_interaction_matches_real_gromosxx_energy` — a real gromosXX-reference test, not another internal cross-check: loads the `pair_lj` reference system directly via `gromos-io` readers, runs `LjCrfInteraction` in-process, and matches the actual gromosXX-produced `expected/energies.tre` to the standard `1e-8` relative tolerance. `pair_lj` was picked because it fits this provider's current scope exactly (solute-only, zero charge, no exclusions)
- **gromos-forces:** `tests/schnet_burnn_reference.rs::schnet_provider_on_real_burnn_methanol_qm_zone` (feature `ml`) — runs `SchNetInteraction` on the real BuRNN tutorial's actual equilibrated, solvated methanol configuration (`eq/eq_meoh_5.cnf`), using the real `QMZONE` atom selection/atomic numbers and `RCUTQM` cutoff read directly off the tutorial's own `.qmm`/`.imd` files, not a synthetic fixture. **Caught a real bug**: `SchNetInteraction` assumed every neighbor pair had both endpoints inside the region, which panicked the instant a real region-inside-a-solvated-system geometry was used (`SpatialIndex::neighbor_pairs` deliberately also returns pairs reaching outside the selection, for FUTURE.md P5's future embedding case). Fixed by skipping out-of-region pairs, documented as a current limitation (no electrostatic embedding support yet) rather than silently worked around
- `cargo test -p gromos-md --test test_gromosXX_references` reconfirmed 37/37 unchanged after both the classical and the `ml`-feature work — none of it touches any existing call site

### Features (P2.7 — electrostatic embedding onto MM atoms)

- **gromos-forces:** `Embedding {None, Mechanical, Electrostatic}` + a defaulted `PotentialProvider::embedding()` — how a provider treats atoms *outside* its region is now a declared property rather than an accident. `SchNetInteraction` declares `None` (its `Dict[str,Tensor]` contract has no environment channel) and **rejects** an unsupported scheme with a clear error instead of silently ignoring it — the P2.6 cross-boundary bug was exactly an unstated assumption here
- **gromos-forces:** `ElectrostaticEmbedding` (`nonbonded/embedding.rs`) — the QM↔MM Coulomb coupling term, and the **first provider to put forces on atoms outside its own region** (P5). Implements Poliak et al. 2025's path (c) (QM charges + pairwise Coulomb computed in the MD code; their paths (a)/(b) are adapters over a QM engine's output, not new physics). Validated against a **closed-form two-point-charge oracle** (exact analytic energy *and* force), plus finite differences on an MM atom in both a toy and the real solvated `t_06` system
- **gromos-forces:** `zones.rs` — `Zone {Inner, Buffer, Outer}` / `PairOwner {Provider, Classical, Embedding}` / `ZonePartition::owner()`: the six-pair-class ownership table that prevents double-counting, derived term-by-term from the BuRNN training target rather than assumed. Deliberately an *orchestration* type, not logic hidden inside a provider; degrades to plain QM/MM without a buffer and to fully-classical without zones

### Testing (P2.7)

- **gromos-forces:** `tests/embedding_gather_reference.rs` — confirms assumption A1 on the real `t_06` system (**1363 MM embedding partners** within 1.4 nm of the QM zone; charges, positions and force-output all reachable through today's trait, so **no signature change was needed**), and finite-difference-validates embedding forces on a real MM atom at scale
- **gromos-forces:** `tests/zone_partition_reference.rs` — the ownership table **partitions** the real system's interactions with no gaps (vanished physics) and no overlaps (double-counting): 3513 atoms → inner 6 / buffer 99 / outer 3408, 6659 QM-cutoff pairs cleanly owned
- **gromos-forces:** `zones::tests` (6) and `nonbonded::embedding::tests` (3) — ownership symmetry, exact exclusion set, degradation to plain QM/MM; analytic Coulomb, MM-atom finite differences, and intra-region pairs proven not double-counted

### Features + Testing (P2.8-2a — NN/MM wired into a real step loop)

- **gromos-md:** `ml` feature (forwards to `gromos-forces/ml`) + `tests/schnet_nve_loop.rs` — `SchNetInteraction` driving actual leapfrog dynamics via `gromos-integrators::LeapFrog`, not another single-point evaluation. 500-step NVE trajectory on a small vacuum molecule: kinetic energy goes from a standing start (0) to 0.081 kJ/mol, total energy stays within 0.0047% of its mean. Chose NVE conservation because it's architecture-blind to whether the (still untrained) potential is chemically meaningful while genuinely testing the loop — force eval → integrate → force eval again — the same way the classical suite validates bulk NVE systems like `water_216_box`, just against a self-consistency bound instead of a gromosXX oracle (none exists for an untrained NN, P8)
- `dt=1e-3` was picked by scanning `1e-4..2e-3` and taking the largest step still producing real dynamics with tight conservation — not guessed; tolerances set at ~100x the measured margin

### Notes

- **Two claims made earlier in this work were wrong and are corrected here**, both caught by reading the tutorial's own source rather than trusting a summary: (1) the BuRNN training target is **QM−QM**, not QM−MM — `get_burnn_energy()` subtracts two MOPAC energies, with the vacuum terms as per-species normalisation; (2) the shipped BuRNN model **predicts no charges** — a single `Atomwise` output, `property=energy`, `derivative=forces` — so the charge-output channel initially thought necessary is not needed for this tutorial at all
- **A planned test was dropped rather than faked.** P2.7 Step 5 was scoped as comparing gromos-rs's classical buffer energy against "the MM baseline the training pipeline used". No such baseline exists: no classical term enters BuRNN training anywhere, and `get_reference_energies()` is not merely "not finished" but **never called**. Building it would have manufactured a false-confidence oracle; it was replaced with the partition-consistency check above, which answers the same underlying question ("how would we detect if the two sides' logics differ?") at the level the contract actually lives
- **Nothing is wired yet.** `ZonePartition` is exported but consumed by nothing; `LjCrfInteraction` does not filter by zone; no binary iterates providers. A QM/MM force term can be *validated* today but a QM/MM simulation cannot be *run* — see PLAN.md P2.8, whose first item is the only blocker that changes this
- This is Dim 12.5 groundwork (the classical-side invariants) plus a first proof that the ML seam actually works end to end — not a production ML engine. Design reviewed from three independent lenses (GROMOS-fidelity/porting-discipline, Rust trait/extensibility, ML-integration prior-art) before implementation; corrections from that review are reflected above.
- **Environment finding worth keeping:** getting `tch`/libtorch to actually run in this sandbox took two real fixes, not zero. (1) The cached conda-forge MKL-linked libtorch failed even to `import torch` in Python with "cannot enable executable stack" — a kernel/sandbox restriction on that specific MKL build, unrelated to Rust; the official PyPI CPU wheel doesn't have this problem. (2) `tch` pins an exact libtorch version per release (`torch-sys`'s build script hard-checks it) — `tch = "0.24.0"` requires precisely **libtorch 2.11.0**; nearby versions (2.13 tried first) fail the check outright rather than risk an ABI mismatch. Both are now recorded in `nonbonded/schnet.rs`'s module docs with the exact working recipe, since this is exactly the kind of non-obvious setup friction FUTURE.md's P2 (model-export boundary) warned would exist.
- **Runtime choice for the eventual ML provider, corrected via research:** GROMOS's real SchNet integration ("BuRNN", Poliak et al. 2025, *J. Comput. Chem.*) embeds SchNetPack via `pybind11` — the Python/GIL tax this design avoids. Checked current in-process prior art: OpenMM-ML and GROMACS's NNPot interface both route through libtorch/TorchScript, not ONNX — no production tool examined has a working ONNX export path for SchNet/MACE-class models. `tch` (TorchScript) is therefore the prior-art-backed reference point for a future ML provider, not `ort`/ONNX. No runtime is committed in this pass — `PotentialProvider` is a trait precisely so each model family can pick what fits it when that provider is actually built.
- `crates/gromos-forces/src/qmmm.rs` (file-in/file-out QM engine scaffolding, forks `xtb`, writes XYZ files) is legacy and was never wired into any binary or test — it predates this design and is not the basis for the new provider seam. Left untouched; its fate (repurpose vs. delete) is not decided here.
- **Checked against GROMOS's real NN(QM)/MM tutorial, not just designed from docs.** Cloned `github.com/biomos/gromos_tutorial_livecoms` into `.local/gromos_tutorial_livecoms` (gitignored, same pattern as `.local/gromosXX`) — Tutorial 6 ("BuRNN") is a working NN(QM)/MM example with a real trained SchNetPack model. Confirmed its `.qmm` file's `QMUNIT` block (`QMULEN/QMUENE/QMUFOR/QMUCHR`) is exactly the per-provider unit-conversion approach this design chose, and that `QMZONE`'s `QMZ` column (atomic number) matches the real model's embedding convention (`Embedding(100, 256, padding_idx=0)`) — the model now indexes elements by atomic number to match. **Ran the real trained model end to end** (schnetpack v1 in a torch venv, with `np.int`/`np.float` shimmed back for modern numpy) — real forward+backward pass, output converted through the tutorial's own unit factors into GROMOS units. Then tried `torch.jit.script` on that real model: **it fails** (`schnetpack/nn/neighbors.py` uses a Python-style conditionally-typed return TorchScript's type system rejects) — a concrete, not hypothetical, confirmation of FUTURE.md P2 ("not every architecture exports cleanly")
- **SchNetPack v1 → 2.x, same session.** That scripting failure is a v1 problem, not an inherent SchNet one: verified `schnetpack.model.NeuralNetworkPotential` (real `SchNet` + an `atomistic.Forces` output module, SchNetPack **2.2.0**) scripts cleanly. `export_toy_schnet.py` and `SchNetInteraction` were both rewritten around the real library instead of the earlier hand-rolled stand-in — see Features/Testing above. `Forces` requires including it explicitly as an output module (a bare `Atomwise` alone fails TorchScript's type inference on an empty `required_derivatives` list)
- Explicitly deferred: wiring any MD binary's force-accumulation loop to iterate `Vec<Box<dyn PotentialProvider>>`; a *trained* ML provider (this pass's model is architecturally real but untrained); making `contribute` async/cancellable for external-process providers; a shared typed-units boundary (each provider currently owns its own unit conversion).

## [0.0.24] (2026-07-02)

### Features

- **pyo3-gromos:** `sim.run(steps, ene_freq=100)` — batch MD loop in Rust returning an `(n_frames, 12)` numpy energy array (`[time, kinetic, potential, total, volume, pressure, bond, angle, improper, dihedral, lj, coulomb]`); no `.tre` file round-trip (P3.3)
- **py-gromos:** `EnergyTimeseries` (`python/gromos/timeseries.py`) wraps the `run()` array: named-column access, `block_average()`, `to_dataframe()` (polars/pandas/dict), `plot()` (plotly/matplotlib); defaults configurable via `gromos.timeseries.config` (P3.3)
- **pyo3-gromos:** `InputParameters.nve/nvt/npt` take `constraints="none"|"hbonds"|"allbonds"` (default `"none"`) mapping to GROMOS `NTC`, plus a readable `.constraints` getter — factory-built params can now SHAKE (P3.5 M1)
- **pyo3-gromos:** `InputParameters.steepest_descent()` now actually minimizes when run through `Simulation` — previously silently fell through to leap-frog at `dt=0`. New `SteepestDescent` algorithm-sequence building block and `AlgorithmSequence.minimize()` preset; `from_parameters()` dispatches to it when `NTEM > 0` (P3.5 SD)
- **pyo3-gromos:** `sim.volume` / `sim.pressure` getters, mirroring `sim.temperature`. `pressure` is only physically meaningful under NPT — the virial term is only populated by `PressureCalculation`, which only NPT's sequence includes; under NVE/NVT it returns the kinetic-only term, not zero, and is documented as such

### Fixes

- **gromos-forces:** `calculate_bonded_forces_ntf` summed bond/angle/dihedral/improper into a single scalar with no per-term breakdown; every call site (`algorithms/forcefield.rs` — the main path every `Simulation.step()` uses, plus `replica.rs`, `gamd.rs`) dumped the combined total entirely into `bond_total`, leaving `angle_total`/`dihedral_total`/`improper_total` at zero for every simulation, not just Python. Added `bond_energy`/`angle_energy`/`dihedral_energy`/`improper_energy` fields to `ForceEnergy`, populated where each term is still separate, before combining. Also fixes the GaMD dihedral-boost term in `md.rs`, which read the same always-zero field
- **pyo3-gromos:** `temperature` getter used bare `n_atoms*3` instead of the constraint-aware degrees of freedom the thermostat actually couples to (`3*n_atoms - solvent_constraint_dof - ndfmin`) — the two could silently disagree. Extracted `compute_total_dof()` as the single source of truth for both builders and the getter (P3.5 M2)

### Testing

- **py-gromos:** `test_run_matches_step_loop`, `test_energy_timeseries` — `run()` matches an equivalent `step(1)` loop; `EnergyTimeseries` column access, `block_average`, `to_dataframe`/`plot` backends
- **py-gromos:** `test_factory_constraints_knob`, `test_constrained_system_stable_with_factory_params` — factory `constraints=` sets `NTC` correctly; `aladip_solvated` (flexible solute H-bonds) stays stable under `constraints="hbonds"` and diverges under `"none"` (documented contrast)
- **py-gromos:** `test_steepest_descent_via_simulation`, `test_steepest_descent_via_algorithm_sequence` — EM actually decreases potential energy and converges; the direct and composable EM paths agree
- **py-gromos:** `test_volume_and_pressure_getters` — `sim.volume`/`sim.pressure` match `run()`'s array columns; NVE/NVT volume is exactly fixed, NPT's responds to the barostat
- **gromos-io:** `test_factory_constraints_knob` — `ntc_from_constraints()` mapping and error case

### Documentation

- **py-gromos/notebooks:** replaced `01_understanding_pyo3_bindings.ipynb`, `02_molecular_systems_and_energy.ipynb`, `03_performance_deep_dive.ipynb` (referenced the nonexistent `gromos.State`) with `01_load_and_inspect.ipynb` and `02_short_md.ipynb` against real reference systems — topology/config inspection, RDF, manual `Topology.solvate()` + `System(topo, conf)` composition, native-parameter MD, energy-component breakdown, block-averaging, NVE/NVT/NPT ensemble comparison (temperature/volume/pressure), and energy minimization, all executed with real baked-in plotly/matplotlib output
- **py-gromos:** `pyproject.toml` `notebooks` dependency group (`jupyter`, `ipykernel`, `matplotlib`, `plotly`, `polars`); CI and Makefile use `uv sync --all-groups` instead of ephemeral `--with` packages
- **py-gromos/docs:** `api/reference.md`, `index.md`, `quick-start.md`, `learning-guide.md` updated for `constraints=`, `sim.run()`/`EnergyTimeseries`, `SteepestDescent`/EM, `sim.volume`/`sim.pressure`; removed the stale "planned, not yet implemented" P3.3 section
- **py-gromos:** `md_runners.py` marked deprecated (docstring notice, points at `Simulation`)
- **PLAN.md / FUTURE.md:** P3.1–3.5 marked done with verification notes; audited and parked a composition-pattern redesign (constraints-on-`System`, unified builders) as tech debt with no triggering need yet — recorded as a deliberate divergence bet in `FUTURE.md`, not scheduled

## [0.0.23] (2026-06-29)

### Features

- **pyo3-gromos:** `PySystem` — paired `Topology` + `Configuration` with atom-count validation at construction; `System.from_files(topo, cnf)` convenience loader; `write(path)` for `.cnf` output (P3.1)
- **pyo3-gromos:** `InputParameters` factory methods `nve(dt, steps)`, `nvt(dt, steps, temperature)`, `npt(dt, steps, temperature, pressure)`, `steepest_descent(steps)` — no `.imd` file required (P3.2)
- **pyo3-gromos:** `InputParameters.from_file(path)` staticmethod alias of the constructor
- **pyo3-gromos:** `Simulation(system, params)` two-argument constructor accepting `System` + `InputParameters`; existing three-argument forms unchanged (P3.2)
- **gromos-io:** `ImdParameters::nve / nvt / npt / steepest_descent` factory methods on the Rust struct

### Testing

- **py-gromos:** `test_system_constructor_matches_file_constructor` — `Simulation(System, params)` is bit-identical to `Simulation(Topology, Configuration, params)` on `water_216_nvt` (energy rel=1e-8)
- **py-gromos:** `test_system_factory_workflow` — full P3.2 workflow: `System.from_files` + `InputParameters.nvt` + `Simulation` + `step(10)` on `water_3_box`
- **py-gromos:** `test_factory_nvt_properties` — factory stores `dt`, `nstlim`, `temperature` correctly
- **gromos-io:** unit tests for all four factory methods (`test_factory_nve/nvt/npt/steepest_descent`)
- **py-gromos:** fixed 6 stale `test_advanced_features.py` test classes missing `@pytest.mark.skip`

### Build

- **Makefile (root):** added `.venv`, `requirements`, `build-python`, `build-release`, `test-python`, `docs`, `docs-serve` targets
- **py-gromos/Makefile:** added `docs`, `docs-serve` targets
- **`.gitignore`:** added `py-gromos/site/`

### Documentation

- **py-gromos/docs:** full rewrite of all seven documentation pages to reflect the working API; removed all phantom class references (`State`, `Box`, `LeapFrog`, GaMD/EDS/REMD runners); added accurate status tables, P3.3/FUTURE roadmap sections
- **py-gromos/python/gromos/__init__.py:** stripped to only working names; removed `system_builder` stubs, `md_runners` functions, and `analysis` wrappers from top-level namespace
- **py-gromos/python/gromos/gromos.pyi:** added `InputParameters` factory stubs, `System` class, `Simulation` two-arg `@overload`

## [0.0.22] (2026-06-27)

### Features

- **gromos-core:** `PairlistAlgorithm` enum (`Standard`, `CellList`) with `from_imd()` dispatch
  and `update<BC>()` delegation via match — zero heap allocation on the hot path (9a-0)
- **gromos-core:** `from_imd()` slots match gromosXX exactly: 0=standard (forced), 1=grid
  (ExtendedGrid fallback→Standard), 2=grid_cell (CellList); IMD parser now recognises
  `"grid_cell"` keyword; no auto-heuristic — every value is an explicit instruction
- **gromos-core:** 9 unit tests covering full dispatch matrix (algorithm × box_type × n_atoms × chargegroups)
- **gromos-integrators, gromos-md, pyo3-gromos:** `StandardPairlistAlgorithm` field swapped to
  `PairlistAlgorithm` everywhere; `md.rs` calls `from_imd(imd.algorithm, …)` (9a-0)

### Testing

- **gromos-md:** `test_pairlist_margin.rs` — empirical CellList vs Standard margin test on
  `water_216_box` (100 steps); observed margin = 0.000e0 (bit-identical) (9a-1)
- **gromos-md:** `water_1000_spc_gridcell` reference — 1000 equilibrated SPC molecules
  (3.1057 nm box, vol7 tutorial), gromosXX Grid_Cell_Pairlist (7×7×7 grid); directly
  validates gromos-rs CellList against gromosXX; completes the proof chain

### Reference tests

- 38/38 pass (37 original + `water_1000_spc_gridcell`)

## [0.0.21] (2026-06-27)

### Features

- **gromos-core:** `gather.rs` — PBC molecule gathering: `gather_chain`, `gather_bond`,
  `gather_molecules`, `centre_of_geometry`, `centre_of_mass`; shared primitive for MD engine
  and analysis tools; 5 unit tests
- **gromos-core:** `fit.rs` moved to `gromos-analysis` (analysis-only); `gather.rs` stays in
  `gromos-core` (shared with MD engine)
- **gromos-analysis:** `fit.rs` — Kabsch rotational fit (Horn 1987 quaternion method);
  `superimpose()`, `rmsd()`, `kabsch_rotation()`, `weighted_cog()`; 7 unit tests
- **gromos-analysis:** `rmsd` binary — real Kabsch fit; @atomspec, @ref, @nofit, @pbc
- **gromos-analysis:** `nhoparam` — N-H bond order parameters S²; rotational fit, window
  averaging, `ee()`; port of GROMOS nhoparam algorithm
- **gromos-forces:** `energy.rs` — `single_point_energy()` + `EnergyParams`; computes
  bonded + nonbonded potential without running MD; used by `ener` binary
- **gromos-analysis:** `ener` — real energy recalculation from trajectory; per-frame
  E_bond, E_lj, E_crf, E_pot; optional PBC gathering
- **gromos-analysis:** `bar` — BAR iteration (numerically stable log-sum-exp, GROMOS
  reference, Shirts 2003); bootstrap error estimation
- **gromos-analysis:** `ext_ti_merge` — linear interpolation between λ windows;
  trapezoidal ΔG from merged curve
- **gromos-analysis:** `frameout` — full GROMOS feature parity: PBC gathering (@pbc),
  atom filtering (@include SOLUTE/SOLVENT/ALL), rotational fit (@ref @atomsfit),
  @spec ALL/EVERY/SPEC, @time range, cnf/pdb/trc output, @single; 8 integration tests
- **gromos-io:** `TrajectoryWriter::write_trc_frame()` — write POSITIONRED frames directly
  from positions (no Configuration needed); standard 3-column GROMOS format
- **gromos-io:** `write_pdb_positions()` in `pdb.rs` — write PDB from raw positions +
  optional topology for atom names; re-exported from lib.rs
- **gromos-io:** `TrajectoryReader` now handles both 3-column (standard GROMOS) and legacy
  7-column POSITIONRED formats; GENBOX made optional (vacuum trajectories)
- **gromos-core:** `selection.rs` — `m:` prefix: selects solute molecules by role
  (mirrors `s:` for solvent); improved error messages with full syntax cheatsheet and
  specific hints for common mistakes (`r:`, `mol:`, bare atom names as prefix)
- **gromos-analysis:** tests — 8 `frameout` integration tests + 7 trajectory round-trip
  unit tests in `gromos-io`

### Refactor

- **gromos-analysis:** `fit.rs` relocated from `gromos-core` → `gromos-analysis`;
  `gromos-forces` added as dependency of `gromos-analysis`
- All analysis binaries use only `gromos-io` for I/O — zero hand-rolled BufWriter/File
  in binary entry points
- PLAN.md restructured: Dim 10 section dissolved into Priority 2; priority order updated
  (CellList wiring → code quality → dihedral restraints → Dim 10 Phase 3)

### Reference tests

- 37/37 pass (unchanged); `frameout` integration tests added to `gromos-analysis`

## [0.0.20] (2026-06-26)

### Refactor

- **gromos-core (Dim 10 Phase 2a–3):** Dissolve the solute/solvent split — full instancing model.
  Direct field access everywhere: `topo.moltypes[0].bonds`, no wrapper methods.
  `MoleculeType` carries atoms + all bonded terms; `Atom = MolTypeAtom` (type alias);
  `Solute` reduced to ZST shell; `Solvent` struct and `Vec<Solvent>` field removed.
  `rebuild_flat_arrays()`: iac/mass/charge derived from instance registry (Phase 2a).
  `all_bonds_global()` / `all_angles_global()` / etc.: instance-iterating force-loop
  iterators with local→global index translation — one loop handles flexible solute,
  flexible solvent (future), and any repeated molecule type (Phase 3).
  `Topology::new()` pre-initializes `moltypes[0]` (SOLUTE) + `instances[0]` so no
  guard methods are needed at write sites.
  Follows GROMACS `gmx_moltype_t` / `gmx_molblock_t` convention exactly.
  37/37 reference tests pass byte-identical.

## [0.0.19] (2026-06-24)

### Features

- **gromos-forces:** perturbed bonded forces (P1.7 Step 2) — quartic bond, cos-harmonic
  angle, improper dihedral, proper dihedral; gromosXX-faithful λ-interpolation + dE/dλ;
  `calculate_perturbed_bonded_forces` wired into `Forcefield::apply` step 2b
- **gromos-forces:** perturbed nonbonded corrections (P1.7 Step 3) — soft-core LJ+CRF
  dual-topology via `perturbed_pairlist_correction`; RF self-energy, excluded-pair, 1-4,
  and PERTATOMPAIR corrections; `ch4_water_fep` (CH4→dummy in 999 SPC water, λ=0.5)
  passes to <1e-6 kJ/mol vs gromosXX
- **gromos-io:** fix PERTURBATION block multi-line parser (ALPHLJ/ALPHC silently dropped);
  fix effective alpha = per-atom × global (matches `in_perturbation.cc:1308`)
- **gromos-io:** `FreeEnergyWriter` + `read_free_energy_trajectory` (FREEENERGY03 blocks);
  md binary writes `@trg` when NTG≠0 at NTWE frequency (P1.7 Step 4)
- **gromos-core:** `Stat` — `ave()`, `rmsd()`, `ee()` block averaging, port of gromos++
  `gmath/Stat` (Allen–Tildesley); 4 unit tests
- **gromos-analysis:** `ext_ti_ana` — real TI analysis; reads N `.trg` files, ⟨dH/dλ⟩ ±
  ee() per window, trapezoidal ΔG integration (P1.7 Step 5 + P2.2 stat + P2.3 program)
- **gromos-core:** full gromos++ AtomSpecifier grammar — `a:NAME`, `1:name,name`,
  `1:res(nr:atom)`, `1:res(name:atom)`, `not(spec)`, `minus(spec)`, `;`-union, `all`/`no`;
  routes via `topology.molecules` (no `num_solute_atoms()` threshold, Dim 10 ready);
  30 tests, every index confirmed against gromos++ `atominfo` on aladip topology
- **gromos-core:** per-atom metadata accessors — `atom_name(i)`, `residue_nr(i)`,
  `residue_name(i)`, `molecule_nr(i)` covering solute + solvent uniformly
- **gromos-analysis:** `atominfo` binary — gromos++-compatible TITLE+ATOMS output
- **gromos-core:** Dimension 10 Phase 1 — instancing model alongside legacy structs;
  `Role`, `MolTypeAtom`, `MoleculeType`, `MoleculeInstance`; `Topology::moltypes` +
  `instances` populated by `init_solute_moltype()` and `solvate()`; per-atom accessors
  prefer moltype path; `s:` in AtomSelection uses `role == Role::Solvent`; `promote()`;
  37/37 reference tests byte-identical

### Test infrastructure

- Shared topology/coord/ptp files moved to `tests/gromosXX_references/shared/`;
  all `input.toml` paths updated from `../../` to `../shared/`
- `.trg` free-energy trajectory tracked in `ch4_water_fep` (dH/dλ at 1e-6 rel tol)
- `run_references.py` removed (both copies); Rust runner is the only harness
- 37 tests pass (up from 36); 1 ignored

## [0.0.18] (2026-06-19)

### Features

- **gromos-core:** add `PerturbedAtom` and `PerturbedAtomPair` structs; extend
  `PerturbedSolute` with `atoms` + `atom_pairs` fields; add `is_perturbed: Vec<bool>`
  to `Topology` (P1.7 Step 1 scaffolding)
- **gromos-io:** replace ptp.rs stub with a full gromosXX-faithful `.pttopo` reader
  (`read_pttopo`); parses PERTATOMPARAM, PERTATOMPAIR, PERTBONDSTRETCH(H),
  PERTBONDANGLE(H), PERTIMPROPERDIH(H), PERTPROPERDIH(H); validated against
  `aladip.pttopo` (3 atoms, 1 pair, 2 bonds, 2 angles, 2 impropers, 2 dihedrals)
- **gromos-md:** wire `@pttopo` into md binary — when NTG≠0, read perturbation topology
  and populate `topo.perturbed_solute` + `topo.is_perturbed`
- **gromos-core:** fix `Energies::total()` to include `distanceres_total` (matches
  gromosXX: `total = potential_total + kinetic_total + special_total` where
  `special_total` includes distance restraints)
- **gromos-md:** add `nacl_1water_distres` reference test — instantaneous distance
  restraint (NTDIR=2) on Na-Cl pair; validates total energy accounting and step-0
  forces; 36 of 36 tests pass

### Reference test matrix

- 36 tests pass (up from 34); new: `nacl_1water_distres` (distance restraint)

## [0.0.17] (2026-06-17)

### Features

- **gromos-forces:** port gromosXX distance restraints (PLAN.md §1.6, first priority)
  - `DistanceRestraint` / `DistanceRestraints` — faithful translation of
    `distance_restraint_interaction.cc`: RAH dimensionality encoding
    (`dim_base ∈ {0,10,20,30,40,50,60}` × form `∈ {-1,0,+1}`), harmonic/linear
    switching (gromosXX `r_linear = DIR0 = 0.3 nm`), mode 1 / mode 2 (NTDIR=2 ×w0)
  - `PerturbedDistanceRestraint` / `PerturbedDistanceRestraints` — translation of
    `perturbed_distance_restraint_interaction.cc`: λ-interpolated r0/w0,
    hidden-restraint prefactor `2^(n+m)·λⁿ·(1-λ)ᵐ`
  - Reference test `test_distance_restraint_gromosxx_reference` validates both at
    gromosXX's own hard-coded values: DistanceRestraint = 257.189539 kJ/mol,
    PerturbedDistanceRestraint = 195.899012 kJ/mol (aladip.distrest, tol 1e-3)
- **gromos-core:** add `DistanceRestraintSpec` / `PerturbedDistanceRestraintSpec` to
  `topology.rs`; extend `Topology` with `distance_restraints` /
  `perturbed_distance_restraints` Vec fields
- **gromos-core:** add `distanceres_total: f64` to `Energy` struct
- **gromos-io:** new `distanceres.rs` parser for DISTANCERESSPEC / PERTDISRESSPEC
  blocks (virtual atom type 0 only; unsupported types logged + skipped)
- **gromos-io:** parse DISTANCERES block (`NTDIR NTDIRA CDIR DIR0 TAUDIR FORCESCALE
  VDIR NTWDIR`) and PERTURBATION block (`NTG NRDGL RLAM DLAMT ALPHLJ ALPHC NLAM
  NSCALE`) in `imd.rs`; add `ImdParameters::lambda_and_derivative()` helper
- **gromos-integrators:** wire `DistanceRestraints` / `PerturbedDistanceRestraints`
  into `Forcefield` (step 4c after position restraints); result accumulated into
  `conf.current_mut().energies.distanceres_total`

### Reference test matrix

- 13 restraint tests pass (13 restraint-specific, up from 0 gromosXX-faithful)
- All other workspace tests unchanged

## [0.0.16] (2026-06-14)

### Features

- **gromos-core:** add `CellListPairlistAlgorithm` (the cell-list pairlist) — a charge-group-aware
  linked-cell pairlist, drop-in via the existing `update<BC>()` interface
  - Bins each chargegroup's reference position (solute CG → COG, solvent CG → first-atom
    position, matching `StandardPairlistAlgorithm`'s own distance conventions) into a grid that
    exactly tiles the box, with cells ≥ `long_range_cutoff + skin` and a periodic-wrapped,
    deduplicated 27-cell neighbor search (correct even when `grid_dim` is 1 or 2 along an axis)
  - O(N) for `BoxType::Rectangular`; vacuum/triclinic/truncated-octahedron boxes fall back to
    `StandardPairlistAlgorithm`'s O(N²) path (correct, not yet accelerated)
  - Pairs are classified by `min/max(cg_a, cg_b)` against `num_solute_chargegroups` — both
    chargegroups' roles — avoiding the Martina solute/solvent misclassification bug
    (`extended_grid_pairlist_algorithm.cc:1309`)
  - Replaces the old, unexported, non-chargegroup-aware `GridCellPairlistAlgorithm` stub
- **gromos-core:** add set-equality tests (`test_cell_list_matches_standard_*`) validating that
  `CellListPairlistAlgorithm` produces identical `solute_short/long`/`solvent_short/long` pair
  sets to `StandardPairlistAlgorithm`, covering periodic wrap, the single-cell fallback, and
  mixed solute+solvent with exclusions

### Reference test matrix

- 35 of 35 tests pass (unchanged — `StandardPairlistAlgorithm` remains the only algorithm wired
  into the `md`/`pyo3-gromos` binaries)

## [0.0.15] (2026-06-14)

### Features

- **gromos-core:** port `truncoct_triclinic_rotmat` / `truncoct_triclinic_box` / `truncoct_triclinic`
  (gromosXX `math/transformation.cc`) for NTB=-1 truncated-octahedron boxes
  - `forward`/`!forward` APIs convert the legacy cube-edge BOX block to/from the lower-triangular
    triclinic box vectors and rotate positions/velocities between the cube and triclinic frames
  - Add `Box::truncated_octahedral()`
- **gromos-core:** replace `Triclinic::nearest_image`'s textbook `frac - frac.round()` reduction
  with gromosXX's iterative while-loop z→y→x reduction (`boundary_implementation.cc:285-318`) —
  not equivalent for strongly triclinic cells (FUTURE.md Dim 11, finding #1)
- **gromos-integrators:** `forcefield.rs` periodicity refresh now also matches
  `BoxType::TruncatedOctahedral`
- **gromos-md:** wire NTB=-1 in the `md` binary — convert the BOX block, rotate positions/velocities
  on read, build `Periodicity::Triclinic` from the resulting box, and rotate
  `FREEFORCERED`/`CONSFORCERED` back to the cube frame on output
  (`out_configuration.cc::_print_forcered`)
- **gromos-md:** add `aladip_trunc_oct` gromosXX reference test (truncated octahedron, NTB=-1)

### Bug fixes

- **gromos-core:** fix `truncoct_triclinic_rotmat`'s forward/backward rotation — gromosXX's
  `product(rot, v)` (`gmath.h`) computes `rot^T * v`, not `rot * v`, because `GenericMatrix`'s
  3-`Vec` constructor stores its arguments as rows while `product()` contracts over the matrix's
  first index. The initial port had the two branches exactly swapped. Added a regression test
  pinned to a gromosXX debug-build reference value for `aladip.conf` atom 21.

### Reference test matrix

- 35 of 35 tests pass (up from 34); new: `aladip_trunc_oct` (truncated octahedron / NTB=-1)

## [0.0.14] (2026-06-12)

### Features

- **gromos-integrators:** implement Nosé-Hoover thermostat — single NHC and NHC chain
  - `NoseHooverThermostat` struct with `new_single_bath` (nhc=1) and `new_chain_bath` (nhc≥2)
  - Single NHC: `ζ[0] += dt/τ² · (T_free/T₀ - 1)`, `scale = 1 - ζ[0]·dt`
  - Chain NHC: sequential tail-to-head ζ update (gromosXX `calc_chain_scaling`)
  - Wired in `md.rs` via `imd.temp_bath[0].algorithm` dispatch (0=Berendsen, 1=NHC, N=chain)
- **gromos-io:** fix MULTIBATH algorithm code mapping to match gromosXX exactly
  (0=weak-coupling, 1=nose-hoover, N≥2=nose-hoover-chains); add `nhc_chain` field
- **gromos-md:** add `water_216_nvt_nosehoover` and `water_216_nvt_nhc_chain` reference tests

### Reference test matrix

- 34 of 34 tests pass (up from 32); new: NHC single, NHC chain

## [0.0.13] (2026-06-12)

### Features

- **gromos-integrators:** wire SETTLE, LINCS, and COM rotation removal into the MD algorithm
  sequence; all three are now reference-tested and passing
  - `SettleAlgorithm` / `LincsAlgorithm` wired for solvent (NTCS=settle/lincs) and solute
    (NTCP=lincs); reference tests: `nacl_1water_settle`, `nacl_1water_lincs`, `aladip_vacuum_lincs`
  - COM rotation removal (`RemoveCOMMotion`) fully wired; reference test: `water_216_box_com_rot`

### Bug fixes

- **gromos-integrators:** fix COM rotation removal for PBC systems
  - Use minimum-image (put_into_box) wrapped positions when computing angular momentum and
    inertia tensor, mirroring gromosXX's `gather_chargegroups` init convention
    (`periodicity.cc:175`); without this the position-COM and angular momentum diverged
  - Suppress *periodic* COM rotation removal when the box is periodic, matching gromosXX
    `configuration.cc:555-560` (`param.centreofmass.remove_rot = false` for non-vacuum PBC)
  - Remove leftover debug `eprintln!` that printed all 648 atom positions at step 0
- **gromos-io:** fix `imd.rs` parsing for NTCP/NTCS lincs/settle keywords and LINCS order params

### Reference test matrix

- 32 of 32 tests pass (up from 28); new: SETTLE, LINCS (solvent+solute), COM rotation removal

## [0.0.12] (2026-06-07)

### Features

- **gromos-core:** add `random` module — GSL-compatible `mt19937` RNG (`GslMt19937`) and
  `gsl_ran_gaussian` (Marsaglia polar method), plus `generate_velocities` for NTIVEL=1
  Maxwell-Boltzmann initial velocity generation, matching gromosXX `util::generate_velocities` /
  `math::RandomGeneratorGSL` bit-for-bit
- **gromos-md:** wire NTIVEL=1 into `md` binary's velocity setup (reads NTIVEL/IG/TEMPI from
  INITIALISE block); add `water_single_genvel` reference test verifying generated velocities
  reproduce gromosXX trajectories exactly

## [0.0.11] (2026-06-07)

### Features

- **gromos-integrators:** add steepest descent energy minimization algorithm
- **gromos-md:** add energy minimization reference tests (vacuum + solvated, with/without SHAKE)
- **gromos-io:** add IFP, MTB, posres, mk_script, jobs, and script_template parsers (~2,600 lines)
- **gromos-tools:** add `prep_posres` utility; rewrite `ion`, `sim_box`, `pdb2g96`, `com_top`, `make_top`, and `mk_script` against the new parsers
- **gromos-tools:** implement `build_box` (replaces hardcoded placeholder) — generates a condensed-phase system on a regular grid sized to a target density, matching gromosXX `build_box.cc`
- **core/forces:** update configuration and position restraints handling

### Performance

- **forces/shake:** optimize nonbonded hot path and SHAKE with cache-friendly data structures
  - Replace `HashSet` exclusions with sorted `Vec` + binary search in topology
  - Flatten LJ parameter matrix into a contiguous `LJParamMatrix` (removes double indirection)
  - Add `ShakeBuffers` for precomputed constraint lists and reusable skip arrays
  - Add CG-grouped nonbonded kernel (`CGPairGroup`): compute PBC nearest-image once per charge-group pair
  - Add parallel (`_parallel`) and no-virial (`_novirial`) innerloop variants and `ForceStorage::merge` for thread-local reduction
  - Add Criterion benchmark suite (`gromos-md/benches/md_bench.rs`)
  - Fix duplicate energy/force accumulation in `calculate_dihedral_new_forces`

### Refactor

- **workspace:** remove the `gromos-cli` crate (superseded by the per-crate binary layout from the 0.0.10 restructure)
- **workspace:** condense and reorganize `PLAN.md` into a priority-ordered roadmap (MD engine physics → analysis foundations → py-gromos/education → code quality)

### Bug Fixes

- **ci:** fix Rust workflow steps using an incorrect `working-directory: gromos-rs` (the workspace lives at the repo root); fix coverage lcov output path and the crates.io publish job directory (`crates/gromos`)

### Chores

- **workspace:** clean up workspace `Cargo.toml`, `Makefile`, and update `py-gromos` reference test paths

## [0.0.10] (2026-05-17)

### Refactor

- **workspace:** restructure crates into focused responsibilities
  - `gromos-md`: 8 simulation engine binaries (md, md_mpi, md_mpi_cuda, mdf, remd, repex_mpi, eds, gamd) + integration tests
  - `gromos-tools`: 30 system construction binaries organized in subdirectories (topology/, box/, conversion/, utilities/)
  - `gromos-analysis`: 66 analysis binaries organized in subdirectories (structural/, energy/, distribution/, dynamics/, free_energy/, trajectory/, noe/, clustering/, xray/, special/) + existing library code
  - `gromos-cli`: slimmed to thin unified `gromos` multicall binary (clap only)
  - `gromos` facade: removed `gromos-analysis` re-export to avoid circular dependency; analysis bins use `gromos_core`/`gromos_io` directly
  - All 21 reference tests pass, workspace compiles clean

## [0.0.9](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.8...v0.0.9) (2026-03-29)

### Bug Fixes

- **forces:** fix CRF energy mismatch (~3.82 kJ/mol) for systems with solvent near cutoff boundary
  - Expand solvent long-range pairlist to all atom pairs (was: first-atom only with shared PBC shift)
  - Remove HEAVISIDE truncation from `lj_crf_innerloop` to match gromosXX `#undef XXHEAVISIDE` default
  - Long-range solvent now processed with `lj_crf_innerloop` (per-atom nearest_image) instead of `solvent_innerloop`

### Features

- **tests:** add water_10_box and nacl_water_box_shifted reference systems
  - water_10_box: 2 ions + 10 SPC waters (32 atoms), positions away from cutoff boundaries
  - nacl_water_box_shifted: nacl_water_box with perturbed positions near cutoff
- **tests:** promote nacl_water_box, nacl_water_box_shifted, nacl_3water_cutoff from ignored to active
  - 14 of 19 reference tests now pass (was 11)
- **io:** add verbose logging and md_output.log capture to run_references.py
- **io:** improve force trajectory writer robustness

## [0.0.8](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.7...v0.0.8) (2026-03-08)

### Features

- **core:** add solute/solvent boundary fields to Topology
  - Add `num_solute_chargegroups` field and `num_solute_atoms()` method
  - Set during `build_topology()` for solute/solvent dispatch

- **forces:** gromosXX-compatible nonbonded force architecture
  - Add HEAVISIDE truncation to `lj_crf_innerloop` (skip pairs beyond cutoff²)
  - Add `solvent_innerloop` with shared PBC shift from O-O nearest_image
  - Split `rf_excluded_interactions` into solute/solvent paths:
    - Solute: self-term + excluded pair forces + full RF energy
    - Solvent: no self-term, no forces, only distance-dependent energy
  - Add `cutoff_sq` to CRFParameters

### Refactor

- **pairlist:** separate solute/solvent CG pairs with correct distance metrics
  - Solute CGs: center-of-geometry distance, exclusion checks
  - Solvent CGs: first-atom position distance (gromosXX convention)
  - Intra-CG non-excluded pairs for solute CGs
  - Solvent-solvent stores first-atom pairs only (expanded in innerloop)
  - Add debug logging for CG positions, distances, pair classification

## [0.0.7](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.6...v0.0.7) (2026-03-02)

### Refactor

- **io:** gromosXX-compatible topology init order with solvate() ([60207d9](https://github.com/yerkoescalona/gromos-rs/commit/60207d9584e67c0ae622ed6267f19fd38965467e))
  - Follow gromosXX initialization: read_parameter → read_topology → solvate → read_configuration
  - Add Topology::solvate(nsm) method matching gromosXX topo.solvate(0, nsm)
  - Add SolventAtomTemplate/SolventConstraintTemplate types on Topology
  - build_topology() now stores solvent template without expanding
  - Remove build_topology_with_solvent() — NSM comes from IMD SYSTEM block
  - Reorder md.rs: IMD → topology+solvate → coordinates (was: IMD → coords → topo)
  - Store chargegroup_codes on Topology for later use
  - LJ matrix now includes solvent IAC types at build time
  - All 6 reference tests pass (pair_lj, pair_lj_mixed, nacl_pair, water_single, benzene_vacuum, water_3_box)

## [0.0.6](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.5...v0.0.6) (2026-02-27)

### Features

- **forces:** gromosXX-compatible bonded forces, topology parsing, and SHAKE ([5d409e7](https://github.com/yerkoescalona/gromos-rs/commit/5d409e7da1ffb9528db09d122e68515bda4df6bc))
  - Fix bonded force vector conventions to match gromosXX (v = pos(i) - pos(j))
  - Simplify quartic bond force: avoid unnecessary r division, use r^2 directly
  - Add NTF flag support to calculate_bonded_forces_ntf() for selective force terms
  - Add debug logging for individual bonded energy components
  - RF excluded corrections: add forces (not just energy), rename to rf_excluded_interactions
  - Parse dihedral/improper types and terms from topology (TORSDIHEDRALTYPE, DIHEDRAL, IMPDIHEDRAL)
  - Fix SOLUTEATOM parser: exclusions on same line (col 8+), skip INE14 line
  - Rewrite SHAKE to match gromosXX convergence (tolerance on |constraint_length - r|)
  - Add shake_algorithm.rs module
  - Add ForceStorage to nonbonded for combined force+energy tracking

## [0.0.5](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.4...v0.0.5) (2026-02-27)

### Refactor

- **io:** gromosXX-compatible CLI and IO cleanup ([c0e9ea2](https://github.com/yerkoescalona/gromos-rs/commit/c0e9ea2603942bed73944dab073cd58045966f5d))
  - Move coordinate reader from md.rs to gromos-io (read_coordinates + CoordinateData)
  - Add GENBOX, POSITIONRED, VELOCITYRED block support in coordinate reader
  - Rewrite IMD parser for gromosXX positional format (comment headers + data lines)
  - Add NTF force flags, COMTRANSROT, shake_tol to ImdParameters
  - Refactor md.rs CLI to match gromosXX @-parameters: @topo @conf @input (required), @fin @trc @tre @trf @trv (output), @pttopo @posresspec @refpos @distrest @angrest @dihrest @gamd (input), @verb @print @version @develop (control)
  - All simulation parameters now come from @input .imd/.in file
  - Fix pre-existing f64/usize type errors in energy_binary and trajectory_binary tests
  - 9 new/updated tests for coordinate and IMD parsing
  - pair_lj reference still passes at ~1e-11 precision

## [0.0.4](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.3...v0.0.4) (2026-02-24)

### Refactor

- **core:** implement gromosXX Algorithm/AlgorithmSequence pattern ([c61cd57](https://github.com/yerkoescalona/gromos-rs/commit/c61cd57c5b7b1eb29c0f067fd02038560b923ab9))
  - Add Algorithm trait and AlgorithmSequence in gromos-core with 7 unit tests
  - Implement algorithm wrappers: Forcefield, LeapFrogVelocity, LeapFrogPosition, TemperatureCalculation, EnergyCalculation
  - TemperatureCalculation uses gromosXX averaged kinetic energy formula
  - Refactor md.rs to use AlgorithmSequence instead of manual force+integrate loop
  - Energy convention follows gromosXX: write to current(), exchange_state moves to old(), read from old() for output
  - Energy output uses full f64 scientific notation for exact comparison
  - All 10 steps of pair_lj match gromosXX reference to ~1e-11 precision
  - Fix minor clippy warnings in gamd, remd, integrator, pairlist

## [0.0.3](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.2...v0.0.3) (2026-02-22)

### Refactor

- **forces:** enhance nonbonded interactions with RF corrections and update math structures ([51b8681](https://github.com/yerkoescalona/gromos-rs/commit/51b86819e74ffcbe918b30b00f180b4c5be2c9ae))
- **core:** support f64 instead of f32 ([ade5423](https://github.com/yerkoescalona/gromos-rs/commit/ade5423ca49fc8974dfdc624f71770a4b7c10a92))

### Chores

- remove dead logging module and replace black with ruff ([de003df](https://github.com/yerkoescalona/gromos-rs/commit/de003dfef3066db053255dbc58087739becf87f7))

## [0.0.2](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.1...v0.0.2) (2025-12-21)

### Refactor

- **workspace:** reorganize into crate-based workspace structure ([7d97a54](https://github.com/yerkoescalona/gromos-rs/commit/7d97a54096d902529e878374395b837717e9343f))

## [0.0.1](https://github.com/yerkoescalona/gromos-rs/commits/v0.0.1) (2025-11-22)

### Refactor

- update force writer integration tests for new API ([ec11525](https://github.com/yerkoescalona/gromos-rs/commit/ec11525f4e08aabc8d621a187cd8875afd3b90f6))

### Chores

- clean and update README files ([27f80bc](https://github.com/yerkoescalona/gromos-rs/commit/27f80bc3abcad0df3d19bd20d5a5291ce7f5fb8e))
- add GROMOS and claude entries to .gitignore ([926ff5c](https://github.com/yerkoescalona/gromos-rs/commit/926ff5c41d8658b964a09aa9d93682f7627ba264))

### Initial

- initial commit ([c0a551c](https://github.com/yerkoescalona/gromos-rs/commit/c0a551cfcac12b0a60daaad908d151f689395732))

## Item index — the stable IDs the code cites

Code comments, tests and older entries above cite work by a stable id (`1.5`, `2.8`, `3.9 A16`,
`P2.1` — the `P` prefix is the same id in an older notation, so `P2.1` is `2.1`).

- **Open work** carries its id as the board card's title (`3.10 step 6` is a card name).
- **Finished work** is the table below, and the version entry above it in this file.

| id | what it was | state |
|---|---|---|
| `1.1` | Reproducibility & correctness | ✓ complete |
| `1.2` | Constraints | ✓ complete |
| `1.3` | Thermostat | ✓ complete |
| `1.4` | Boundary | ✓ complete |
| `1.5b` | Real-input compatibility | ✓ 0.0.42–0.0.46 |
| `1.5c` | Tutorial 2 opened | ✓ 0.0.48 |
| `1.7` | FEP / TI | ✓ complete for fixed λ |
| `1.8` | Free-energy completeness | step 0 ✓ 0.0.52 (energy groups, whole-block energy tests); step 1 ✓ 0.0.53 (LAMBDAS, NLAM, perturbed SHAKE/restraint dH/dλ); step 2 ✓ 0.0.54 (PRECALCLAM, RF at EPSRF = 0 / from RCRF); 0.0.55 fixes (PERTATOMPAIR, constraint virial); steps 3–7 open |
| `2.0` | Instancing model, phases 1–2e | ✓ |
| `2.1` | Atom selection | ✓ complete |
| `2.2` | Shared analysis infra | ✓ complete |
| `2.3` | Real program implementations | ✓ complete |
| `2.3b` | Free-energy estimators | core done |
| `2.4` | Code quality | done 2026-08-30 |
| `2.4b` | gromos-tools audit | done 2026-08-30 |
| `2.5` | Stub cleanup | the real implementations landed; the remaining stubs are open |
| `2.6` | The provider seam (QM/MM/ML) | ✓ landed 2026-08-13 |
| `2.7` | Electrostatic embedding onto MM atoms (QM/MM + BuRNN) | ✓ Steps 1–5 |
| `2.8` | Wiring: from "the pieces exist" to "a QM/MM run happens" | ✓ done except 4b and the zone-aware Forcefield |
| `2.9` | The QM/MM → ML/MM pipeline: generate, train, compare | ✓ ran end to end |
| `2.10` | Energy library profile: same grammar, stricter rules, versioned | ✓ implemented 0.0.47 |
| `3.0` | Approach | ✓ |
| `3.1` | System.from_files() | ✓ done — System = Topology + Configuration, atom-count |
| `3.2` | Run a System with native parameters | ✓ done — ImdParameters::default |
| `3.3` | Reporters (energy out without files) | ✓ done — angle/dihedral/improper wired |
| `3.4` | Notebooks | ✓ done — replaced 3 stale notebooks referencing phantom APIs with |
| `3.5` | M1/M2/SD convenience-path fixes | ✓ done |
| `3.6` | Rust↔Python reference-test parity gap | ✓ done, FEP deferred |
| `3.7` | ML potential binding for py-gromos: name zones, don't count atoms | ✓ done |
| `3.8` | py-gromos composition audit | ✓ 2026-08-29 |
| `3.9` | The composition model: one recipe, one builder, three front-ends | ✓ steps 0–5 done |
| `4` | Benchmarking | ✓ Phases 0–4 done |
