# spec/ — Lean specifications of newen logic

Full explanation (why, how the pieces fit, guarantees, adding a function): [docs/src/verification/index.md](../docs/src/verification/index.md).

Verified so far: `schedule::output_misses_longrange` (the pilot, explained below),
`thermo::dof::total_dof`, the five functions of `thermo::bath` (degrees of freedom per
temperature bath and the atom-to-bath mapping), and the five functions of `topology`
(`pairs::contains_sorted` and `pairs::pair_kind`: excluded, 1-4, owned or normal pair;
`groups::group_of`, `groups::check_last_atoms` and `groups::energy_cell`: the energy group of an
atom, the FORCE `NRE` rule, the energy-matrix cell); and thirteen floating-point kernels: the
Berendsen and Nosé–Hoover thermostat steps of `crates/newen-integrators`, the Lennard-Jones and
reaction-field pair term of `crates/newen-forces` with its reaction-field constants, excluded-pair
and self terms, and the seven bonded terms of `crates/newen-forces` (quartic and harmonic bond,
cosine-harmonic and harmonic angle, the two proper-dihedral forms, the improper dihedral; see
"Floating-point kernels" below and `spec/SPEC.toml`). The pilot: `output_misses_longrange(nsnb, ntwg)` in
`crates/newen-logic/src/schedule.rs`. It answers: with a twin-range pairlist updated every
`nsnb` steps and free-energy output every `ntwg` steps, does some output step reuse long-range
terms from a pairlist step that was *not* an output step? (Under PRECALCLAM gromosXX then writes
λ_k rows without the long-range part; `newen-run` warns about it.)

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
cd spec && lake exe spec-vectors         # rewrite vectors/schedule.json from the Lean model
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
(`crates/newen-forces/src/nonbonded/`), and `quartic_bond`, `harmonic_bond`, `cosine_angle`,
`harmonic_angle`, `dihedral_chebyshev`, `dihedral_new` and `improper_dihedral`
(`crates/newen-forces/src/bonded/`) compute with `f64`, so Aeneas does not apply
(`kind = "float"` in `SPEC.toml`). Three pieces instead:

| Piece | Claim | Where |
|---|---|---|
| ℝ theorems | what the formula means: Berendsen relaxes `T` towards `T₀` by `dt/τ` per step and never overshoots for `dt ≤ τ`; the Nosé–Hoover friction grows exactly when the bath is too hot; the LJ and reaction-field scalar forces are `−E'(r)/r` (`lj_force_is_minus_grad`, `crf_force_is_minus_grad`), the LJ energy vanishes at σ and has its minimum `−C6²/(4 C12)` at `2^(1/6) σ`, the reaction-field energy vanishes at the cutoff, and the excluded-pair and self terms are the pair term without its Coulomb part; the bond forces and the cosine-angle forces (on all three atoms) are `−∇V` (`quarticBond_force_is_minus_grad`, `harmonicBond_force_is_minus_grad`, `cosineAngle_force_is_minus_grad`), a quartic bond matches a harmonic one of force constant `2 K b0²` (`quartic_harmonic_match`), the code's dihedral polynomials are `cos mφ` and their derivatives (`cosMultiple_eq_cos`, `dcosMultiple_is_deriv`), the two dihedral forms agree for `δ ∈ {0, π}` (`dihedral_new_eq_dihedral`), and the angle and dihedral forces sum to zero and exert no torque (`cosineAngle_torque_zero`, `chebyshev_conservation`, `newDihedral_conservation`, `improper_conservation`, …) | `NewenTheory/Dynamics/Thermostats.lean`, `NewenTheory/Interactions/LennardJones.lean`, `NewenTheory/Interactions/ReactionField.lean`, `NewenTheory/Interactions/Bonded.lean` (import Mathlib) |
| Float model | one definition per kernel, polymorphic in the number type and written operation by operation like the Rust; at ℝ it equals the theory (`berendsenModel_real_eq`, `noseHooverModel_real_eq`, `crfParametersModel_real_eq`, `ljCrfPairModel_real_eq`, `rfExcludedPairModel_real_eq`, `rfSelfTermModel_real_eq`, and one `*Model_real_eq` per bonded kernel), at `Float` it is executed | `NewenSpec/Thermostat.lean`, `NewenSpec/Nonbonded.lean`, `NewenSpec/Bonded.lean` (the shared `Sqrt`, `Acos`, `Cos`, `Sin` classes in `NewenSpec/Float.lean`) |
| Differential test | the compiled Rust returns the same 64 bits as the `Float` model on every generated input (6 840 Berendsen rows, 3 036 Nosé–Hoover rows, chains iterated over 1–8 steps; 600 reaction-field constant rows, 6 000 LJ + CRF pair rows, 1 500 excluded-pair rows, 300 self-term rows; 1 482 bond rows, 1 000 angle rows, 710 dihedral rows, each evaluated by every kernel of its kind; the generator makes as many more as a nightly run wants); for the bonded kernels also a finite-difference check of the forces against the energy on 200 rows per kernel | `crates/newen-integrators/tests/thermostat_spec.rs`, `vectors/thermostat.json`; `crates/newen-forces/tests/nonbonded_spec.rs`, `vectors/nonbonded.json`; `crates/newen-forces/tests/bonded_spec.rs`, `vectors/bonded.json` |

The four-lane copy of the pair kernel, `lj_crf_interaction_x4`, is checked against the scalar
`lj_crf_pair` on the same rows, lane by lane and bit for bit, by a unit test in
`crates/newen-forces/src/nonbonded/innerloops.rs`.

The comparison is exact, not within a tolerance: the kernels use `+ − × ÷` and `sqrt`, which
IEEE 754 rounds correctly, in Rust as in Lean's `Float` (C `double`), so the same operations in
the same order give the same bits. The bonded kernels also call `acos`, `cos` and `sin`, which
IEEE 754 does not round correctly; there the bits agree because Lean's `Float.acos`, `Float.cos`,
`Float.sin` and Rust's `f64::acos`, `cos`, `sin` call the same C library (glibc's `libm` on
Linux), so the check is exact on one machine and says nothing about another `libm`. The constants
(`k_B`, `ZERO_TEMPERATURE`, `ZERO_DISTANCE2`, the bonded guards, `π`) and every input are written
as bits. The one exception is NaN: Lean has a single NaN, so a NaN matches any NaN.

```sh
cd spec && lake exe spec-vectors                              # rewrite vectors/thermostat.json, nonbonded.json, bonded.json
cargo test -p newen-integrators --test thermostat_spec        # Rust = Float model, bit for bit
cargo test -p newen-forces --test nonbonded_spec              # Rust = Float model, bit for bit
cargo test -p newen-forces --test bonded_spec                 # Rust = Float model, bit for bit; forces = −∇E by finite differences
cargo test -p newen-forces --lib x4_lanes                     # SIMD lanes = scalar kernel, bit for bit
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
