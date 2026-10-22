# spec/ — Lean specifications of newen logic

Full explanation (why, how the pieces fit, guarantees, adding a function): [docs/src/verification/index.md](../docs/src/verification/index.md).

Verified so far: `schedule::output_misses_longrange` (the pilot, explained below), the run
schedule (`schedule::output::{due, frames_in, minimisation_stop, final_step}`: which steps write
frames, how many, what a converged minimisation writes on the step it stops at, the step of the
final configuration; `schedule::com::com_removal`: which steps remove the centre-of-mass motion),
`thermo::dof::total_dof`, the six functions of `thermo::bath` (degrees of freedom per
temperature bath, the atom-to-bath mapping, and which baths are empty and decoupled), and the
five functions of `topology` (`pairs::contains_sorted` and `pairs::pair_kind`: excluded, 1-4,
owned or normal pair; `groups::group_of`, `groups::check_last_atoms` and `groups::energy_cell`:
the energy group of an atom, the FORCE `NRE` rule, the energy-matrix cell); and forty-one
floating-point kernels: the Berendsen and Nosé–Hoover thermostat steps of
`crates/newen-integrators`, the Lennard-Jones and reaction-field pair term of `crates/newen-forces`
with its reaction-field constants, excluded-pair and self terms, the seven bonded terms of
`crates/newen-forces` (quartic and harmonic bond, cosine-harmonic and harmonic angle, the two
proper-dihedral forms, the improper dihedral), the dynamics kernels of `crates/newen-integrators`
(the leap-frog velocity and position steps, the kinetic energy of the atoms and of a temperature
group, the temperature, the pressure tensor and scalar pressure, and the Berendsen barostat
factor), the periodic boundaries of `crates/newen-core` (nearest image and put-into-box for
vacuum, rectangular and triclinic boxes, the truncated octahedron being a triclinic box, with the
two box constructors), and the free-energy kernels of `crates/newen-forces`: the λ table, the
soft-core LJ + CRF pair, the soft excluded-pair and self terms, the PRECALCLAM end-state kernels,
and the perturbed bonded terms (see "Floating-point kernels" below and `spec/SPEC.toml`). The pilot:
`output_misses_longrange(nsnb, ntwg)` in `crates/newen-logic/src/schedule/mod.rs`. It answers:
with a twin-range pairlist updated every `nsnb` steps and free-energy output every `ntwg` steps,
does some output step reuse long-range terms from a pairlist step that was *not* an output step?
(Under PRECALCLAM gromosXX then writes λ_k rows without the long-range part; `newen-run` warns
about it.)

The meaning lives in `NewenTheory/` (physics, organised like the physics; it never imports the
code's spec), the model and the bridge theorem in `NewenSpec/` (organised like the crate). Three
statements are checked, each by a different tool, and together they close a loop:

```
             meaning  (Misses: "∃ output step whose pairlist step is not an output step")
              ▲    ▲
   Layer 1    │    │  Layer 3
   Lean proof │    │  exhaustive test, 1 ≤ nsnb, ntwg ≤ 200
              │    │
   model ◄────┼── Rust function
 (missesModel)  Layer 2: Aeneas refinement proof
```

| Layer | Claim | Where | Run | Guarantees | Does not guarantee |
|---|---|---|---|---|---|
| 1. Lean theorem | model = meaning | `NewenTheory/Dynamics/Schedule.lean` (`Misses`), `NewenSpec/Schedule.lean` (`missesModel_iff`) | `lake build` | for **all** naturals, the closed form "nsnb > 1, ntwg > 0, neither divides the other" is exactly `Misses` | that `Misses` is what gromosXX does: that is a reading of `perturbed_nonbonded_set.cc`, checked by humans and by the reference tests |
| 2. Aeneas refinement | Rust = model | `NewenSpec/Extracted/` (generated), `NewenSpec/Refinement.lean` | `scripts/extract.sh`, then `lake build` | for **all** `usize` inputs the Rust function never panics and returns `missesModel` | that Charon/Aeneas translate Rust faithfully (trusted), or anything about callers such as `prepare.rs` |
| 3. Differential test | Rust = meaning, Rust = model, on a grid | `crates/newen-logic/tests/schedule_spec.rs`, `vectors/schedule.json` | `cargo test -p newen-logic` | the compiled Rust agrees with a step-by-step simulation and with the Lean model on 40 000 inputs | anything outside the grid |

Layers 1 + 2 are the proof. Layer 3 needs no Lean and runs with every `cargo test`; it also tests the
*compiled* binary, which the proofs do not see, and its simulation is written independently of
the formula, so a wrong `pairlistStep` definition in Lean would show up there.

## Running

```sh
cd spec && lake build                    # layers 1 and 2: proofs (first run downloads Mathlib's cache)
cd spec && lake exe spec-vectors         # rewrite vectors/*.json from the Lean models
spec/scripts/extract.sh                  # regenerate NewenSpec/Extracted/ after editing the Rust
cargo test --release -p newen-logic     # layer 3
```

`lake build` must end with no `sorry` warning; `#print axioms` on the theorems lists only
`propext`, `Classical.choice`, `Quot.sound` (Lean's standard axioms).

**Toolchain.** Layer 1 uses core Lean only (`omega`, `Nat` lemmas);
`lean NewenTheory/Dynamics/Schedule.lean` checks the theory with no packages at all.
`NewenTheory` and `NewenSpec` are two libraries of one Lake project, on one toolchain. The
Aeneas library needs Mathlib and fixes the toolchain: `lean-toolchain` (v4.35.0) and the commit
in `lakefile.lean` must match the Aeneas binary.
**Install** (no root): download `aeneas-linux-x86_64.tar.gz` from the Aeneas GitHub release whose
tag ends in the commit pinned in `lakefile.lean` and unpack it into `~/.local/aeneas/aeneas` (it
bundles the matching `charon`); then `rustup toolchain install <channel in its rust-toolchain>
--profile minimal -c rustc-dev -c llvm-tools -c rust-src`. `extract.sh` refuses to run if the
binary and the pinned library differ.

## Floating-point kernels

`berendsen_scale` and `nose_hoover_scale` (`crates/newen-integrators/src/algorithms/`),
`CRFParameters::new`, `lj_crf_pair`, `rf_excluded_pair` and `rf_self_term`
(`crates/newen-forces/src/nonbonded/`), `quartic_bond`, `harmonic_bond`, `cosine_angle`,
`harmonic_angle`, `dihedral_chebyshev`, `dihedral_new` and `improper_dihedral`
(`crates/newen-forces/src/bonded/`), `leap_frog_velocity`, `leap_frog_position`,
`kinetic_energy`, `group_kinetic_energy`, `kinetic_temperature`, `pressure_tensor`,
`scalar_pressure` and `berendsen_barostat_scale` (`crates/newen-integrators/src/algorithms/`),
`vacuum_nearest_image`, `Rectangular::new`, `rectangular_nearest_image`,
`rectangular_put_into_box`, `Triclinic::new`, `triclinic_nearest_image` and
`triclinic_put_into_box` (`crates/newen-core/src/math.rs`), and the perturbed kernels of
`crates/newen-forces/src/nonbonded/perturbed.rs` and `src/bonded/perturbed.rs` compute with
`f64`, so Aeneas does not apply (`kind = "float"` in `SPEC.toml`). Three pieces instead:

| Piece | Claim | Where |
|---|---|---|
| ℝ theorems | what the formula means: Berendsen relaxes `T` towards `T₀` by `dt/τ` per step and never overshoots for `dt ≤ τ`; the Nosé–Hoover friction grows exactly when the bath is too hot; the LJ and reaction-field scalar forces are `−E'(r)/r` (`lj_force_is_minus_grad`, `crf_force_is_minus_grad`), the LJ energy vanishes at σ and has its minimum `−C6²/(4 C12)` at `2^(1/6) σ`, the reaction-field energy vanishes at the cutoff, and the excluded-pair and self terms are the pair term without its Coulomb part; the bond forces and the cosine- and harmonic-angle forces (on all three atoms; the harmonic angle for `0 < θ < π`) are `−∇V` (`quarticBond_force_is_minus_grad`, `harmonicBond_force_is_minus_grad`, `cosineAngle_force_is_minus_grad`, `harmonicAngle_force_is_minus_grad`), and so are the Chebyshev dihedral's forces on its outer atoms (`chebyshev_force_is_minus_grad`), a quartic bond matches a harmonic one of force constant `2 K b0²` (`quartic_harmonic_match`), the code's dihedral polynomials are `cos mφ` and their derivatives (`cosMultiple_eq_cos`, `dcosMultiple_is_deriv`), the two dihedral forms agree for `δ ∈ {0, π}` (`dihedral_new_eq_dihedral`), and the angle and dihedral forces sum to zero and exert no torque (`cosineAngle_torque_zero`, `chebyshev_conservation`, `newDihedral_conservation`, `improper_conservation`, …); leap-frog is time-reversible, exact for a constant acceleration and area preserving in one dimension, its half-step kinetic energy is the velocity-Verlet one plus `m dt² \|a\|²/8`, a group's kinetic energy splits into centre of mass and a non-negative rest (König), `T = 2K/(N_df k_B)` inverts equipartition; the pressure is the virial theorem's `2/V · (K − Ξ)` (the ideal gas without forces), and the Berendsen barostat leaves `P₀` fixed, grows the box exactly when `P > P₀` and relaxes the pressure by `dt/τ` under linear response; a nearest image and a put-into-box differ from their input by a lattice vector, lie in their cell (rectangular `[−L/2, L/2]` and `[0, L)`; triclinic fractional coordinates in `(−1/2, 1/2]` and `[0, 1)`), are the only point of the cell that does, and are idempotent; the rectangular image is the shortest, the triclinic one not always (`image_not_shortest`); for the free-energy kernels, the kernel's dV/dλ is the derivative of its energy in λ (`softcore_dvdl_is_deriv` for the soft-core LJ and CRF terms with individual λ, `pert_quartic_dvdl_is_deriv`, `pert_harmonic_dvdl_is_deriv`, `soft_dvdl_is_deriv` for the bonded ones), λ = 0 and λ = 1 are the plain terms of states A and B (`softcore_lambda_zero`, `softcore_lambda_one`, `pert_*_ends`), the soft-core force, excluded pairs included, is `−∂V/∂r / r`, and the Chebyshev polynomials are `cos(mφ)` | `NewenTheory/Dynamics/Thermostats.lean`, `NewenTheory/Interactions/LennardJones.lean`, `NewenTheory/Interactions/ReactionField.lean`, `NewenTheory/Interactions/Bonded.lean`, `NewenTheory/Dynamics/Integration.lean`, `NewenTheory/Dynamics/Pressure.lean`, `NewenTheory/Space/Periodicity.lean`, `NewenTheory/Interactions/SoftCore.lean`, `NewenTheory/Interactions/PerturbedBonded.lean` (import Mathlib) |
| Float model | one definition per kernel, polymorphic in the number type and written operation by operation like the Rust; at ℝ it equals the theory (`berendsenModel_real_eq`, `noseHooverModel_real_eq`, `crfParametersModel_real_eq`, `ljCrfPairModel_real_eq`, `rfExcludedPairModel_real_eq`, `rfSelfTermModel_real_eq`, one `*Model_real_eq` per bonded kernel, and for the dynamics kernels `leapFrogModel_real_eq`, `kineticEnergyModel_real_eq`, `groupKineticModel_real_eq`, `kineticTemperatureModel_real_eq`, `pressureTensorModel_real_eq`, `scalarPressureModel_real_eq`, `berendsenBarostatModel_real_eq`, and for the periodic boundaries `rectangularNearestImageModel_real_eq`, `rectangularPutIntoBoxModel_real_eq`, `triclinicNearestImageModel_real_eq`, `triclinicPutIntoBoxModel_real_eq`, and `perturbedLjCrfPairModel_real_eq` and the other `*_real_eq` of the perturbed kernels; `perturbedLjCrfPairModel_dvdl_is_deriv`: the pair kernel's dV/dλ is the derivative of its energy, with the λ table built by the table kernel), at `Float` it is executed; the triclinic `while` loops are `whileFuel`, bounded by a fuel that the ℝ theorem requires to exceed every fractional coordinate by one | `NewenSpec/Thermostat.lean`, `NewenSpec/Nonbonded.lean`, `NewenSpec/Bonded.lean`, `NewenSpec/Dynamics.lean`, `NewenSpec/Periodicity.lean`, `NewenSpec/SoftCore.lean`, `NewenSpec/PerturbedBonded.lean` (the shared `Sqrt`, `Acos`, `Cos`, `Sin`, `Powf`, `Rint` and `Floor` classes in `NewenSpec/Float.lean`) |
| Differential test | the compiled Rust returns the same 64 bits as the `Float` model on every generated input (6 840 Berendsen rows, 3 036 Nosé–Hoover rows, chains iterated over 1–8 steps; 600 reaction-field constant rows, 6 000 LJ + CRF pair rows, 1 500 excluded-pair rows, 300 self-term rows; 1 482 bond rows, 1 000 angle rows, 710 dihedral rows, each evaluated by every kernel of its kind; 1 200 leap-frog rows, 500 groups of 1–8 atoms for the kinetic energies, 800 temperature rows, 600 pressure rows, 1 500 barostat rows; 4 420 periodic-boundary rows over vacuum, rectangular, truncated-octahedron and triclinic boxes; 2 500 rows of the seven soft-core kernels, 1 900 rows of the six perturbed bonded kernels, whose dV/dλ and forces are also checked against central differences of their own energies; the generator makes as many more as a nightly run wants); for the bonded kernels also a finite-difference check of the forces against the energy on 200 rows per kernel. The finite-difference checks compare the Rust with itself (its forces or dV/dλ against central differences of its own energy): they show the kernel is consistent, not that it computes the theory | `crates/newen-integrators/tests/thermostat_spec.rs`, `vectors/thermostat.json`; `crates/newen-forces/tests/nonbonded_spec.rs`, `vectors/nonbonded.json`; `crates/newen-forces/tests/bonded_spec.rs`, `vectors/bonded.json`; `crates/newen-integrators/tests/dynamics_spec.rs`, `vectors/dynamics.json`; `crates/newen-core/tests/periodicity_spec.rs`, `vectors/periodicity.json`; `crates/newen-forces/tests/softcore_spec.rs`, `vectors/softcore.json`; `crates/newen-forces/tests/perturbed_bonded_spec.rs`, `vectors/perturbed_bonded.json` |

The four-lane copy of the pair kernel, `lj_crf_interaction_x4`, is checked against the scalar
`lj_crf_pair` on the same rows, lane by lane and bit for bit, by a unit test in
`crates/newen-forces/src/nonbonded/innerloops.rs`.

The comparison is exact, not within a tolerance: the kernels use `+ − × ÷` and `sqrt`, which
IEEE 754 rounds correctly, and `floor` and round-to-nearest-even, which are exact, in Rust as in
Lean's `Float` (C `double`), so the same operations in the same order give the same bits. Rust's
`f64::round_ties_even` is C `rint`, which Lean's `Float` lacks (its `Float.round` is C `round`,
halves away from zero, Rust's `f64::round`), so `NewenSpec/Float.lean` binds `rint` itself
(`rintFloat`) and models it over ℝ as `roundEven`. The bonded kernels also call `acos`, `cos`
and `sin`, the Berendsen barostat `powf(x, 1/3)` and the λ table `powf`, which IEEE 754 does not
round correctly; there the bits agree because Lean's `Float.acos`, `Float.cos`, `Float.sin`,
`Float.pow` and Rust's `f64::acos`, `cos`, `sin`, `powf` call the same C library (glibc's `libm`
on Linux), so the check is exact on one machine and says nothing about another `libm`. Replacing
`powf(…, 1.0 / 3.0)` by `cbrt` makes 52 of the 1 500 barostat rows differ. The constants (`k_B`,
`ZERO_TEMPERATURE`, `ZERO_DISTANCE2`, the bonded guards, `π`, the `0.5` and `1/3` literals) and
every input are written as bits. The one exception is NaN: Lean has a single NaN, so a NaN
matches any NaN.

```sh
cd spec && lake exe spec-vectors                              # rewrite every file in vectors/
cargo test -p newen-integrators --test thermostat_spec        # Rust = Float model, bit for bit
cargo test -p newen-integrators --test dynamics_spec          # Rust = Float model, bit for bit
cargo test -p newen-forces --test nonbonded_spec              # Rust = Float model, bit for bit
cargo test -p newen-forces --test bonded_spec                 # Rust = Float model, bit for bit; forces = −∇E by finite differences
cargo test -p newen-core --test periodicity_spec              # Rust = Float model, bit for bit
cargo test -p newen-forces --lib x4_lanes                     # SIMD lanes = scalar kernel, bit for bit
cargo test -p newen-forces --test softcore_spec               # soft-core kernels = Float model; dV/dλ, force by FD
cargo test -p newen-forces --test perturbed_bonded_spec       # perturbed bonded kernels = Float model; dV/dλ by FD
```

Not guaranteed: the distance between the ℝ model and the `Float` model (the rounding error is not
bounded here), and anything outside the generated inputs.

## Adding a function

1. Write it in `crates/newen-logic` in the subset Aeneas translates (safe Rust, integers and
   booleans, no floats, simple loops); end its doc comment with `/// spec: <module>.<name>`.
2. The meaning (a `Prop` over all inputs) in `NewenTheory/<Area>/<Topic>.lean`, placed by its
   physics; in `NewenSpec/<Module>.lean` a model written like the Rust, and a theorem model ↔
   meaning.
3. `scripts/extract.sh`, then a theorem in `Refinement.lean`: extracted function `= ok (model ..)`
   (usually `step*` plus `scalar_tac`/`simp_all`).
4. A test in `crates/newen-logic/tests/` comparing against an independent brute-force meaning and,
   via `SpecVectors.lean`, against the model.
5. A `[[function]]` entry in `SPEC.toml`.
