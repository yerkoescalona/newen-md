import NewenSpec.Float
import NewenTheory.Interactions.LennardJones
import NewenTheory.Interactions.ReactionField

/-!
# Model of the Lennard-Jones and reaction-field pair kernels of `crates/newen-forces`

Rust (spec keys in `spec/SPEC.toml`):

* `nonbonded.crf_parameters`: `CRFParameters::new` in `src/nonbonded/params.rs`, the reaction-field
  constants of a run;
* `nonbonded.lj_crf_pair`: `lj_crf_pair` in `src/nonbonded/innerloops.rs`, one LJ + CRF pair;
* `nonbonded.rf_excluded_pair`, `nonbonded.rf_self_term`: `rf_excluded_pair` and `rf_self_term` in
  `src/nonbonded/excluded.rs`, the reaction-field terms of excluded solute pairs and of each
  solute atom with itself.

Theory: `NewenTheory/Interactions/LennardJones.lean` and `ReactionField.lean`.

Each kernel is written once, over any type `α` with the operations the Rust uses
(`+ − × ÷`, negation, `sqrt` from `NewenSpec/Float.lean`, `<`, `≤`, `==`, the literals
`0 1 2 6 12`), operation by operation in the Rust order, and used twice, as in
`NewenSpec/Thermostat.lean`:

1. at `ℝ`, where the `*_real_eq` theorems prove it equal to the theory outside the guards;
2. at `Float`, where `lake exe spec-vectors` evaluates it and
   `crates/newen-forces/tests/nonbonded_spec.rs` compares the Rust kernels bit for bit.

The two literals that are not small integers, the zero-distance guard `1e-10` and the `0.5` of
the self term, are parameters, given at `Float` by their bits.

What is not proved: the rounding error between the two instantiations.
-/

namespace NewenSpec.Nonbonded

open NewenTheory.Interactions.LennardJones NewenTheory.Interactions.ReactionField

/-- `CRFParameters`: the reaction-field constants a pair kernel reads. -/
structure CrfParams (α : Type) where
  /-- `(1 − C_RF/2)/R`, the energy shift. -/
  crfCut : α
  /-- `C_RF/(2R³)`, the coefficient of `r²` in the energy. -/
  crf2Cut3i : α
  /-- `C_RF/R³`, the reaction-field part of the scalar force. -/
  crfCut3i : α
  /-- `R²`. -/
  cutoffSq : α

/-! ## 1. Model (mirrors the Rust functions) -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α] [LE α]
  [DecidableRel (α := α) (· < ·)] [DecidableRel (α := α) (· ≤ ·)] [BEq α] [Sqrt α]
  [OfNat α 0] [OfNat α 1] [OfNat α 2] [OfNat α 6] [OfNat α 12]

/-- The constants `CRFParameters::new` derives from `C_RF` and the radius `R` (`cutoff`). -/
def crfDerivedModel (crf cutoff : α) : CrfParams α :=
  let cut3i := 1 / (cutoff * cutoff * cutoff)
  let crfCut3i := crf * cut3i
  let crf2Cut3i := crfCut3i / 2
  let crfCut := (1 - crf / 2) / cutoff
  ⟨crfCut, crf2Cut3i, crfCut3i, cutoff * cutoff⟩

/-- `C_RF` as `CRFParameters::new` computes it: `−1` for the sentinel `ε_rf = 0` (a conducting
boundary), else the Tironi formula with `κR` and `(κR)²` formed first. -/
def crfConstantModel (cutoff epsilon rfEpsilon rfKappa : α) : α :=
  let kappaCut := rfKappa * cutoff
  let kappaCut2 := kappaCut * kappaCut
  if rfEpsilon == 0 then -1
  else
    (2 * (epsilon - rfEpsilon) * (1 + kappaCut) - rfEpsilon * kappaCut2) /
      ((epsilon + 2 * rfEpsilon) * (1 + kappaCut) + rfEpsilon * kappaCut2)

/-- `CRFParameters::new(cutoff, epsilon, rf_epsilon, rf_kappa)`: all zero for `cutoff ≤ 0` (no
reaction field), else the constants derived from `C_RF`. -/
def crfParametersModel (cutoff epsilon rfEpsilon rfKappa : α) : CrfParams α :=
  if cutoff ≤ 0 then ⟨0, 0, 0, 0⟩
  else crfDerivedModel α (crfConstantModel α cutoff epsilon rfEpsilon rfKappa) cutoff

/-- `lj_crf_pair(r2, c6, c12, q_prod, crf)`: `(f_mag, e_lj, e_crf)`, with `ZERO_DISTANCE2` as
`guard`. -/
def ljCrfPairModel (guard r2 c6 c12 q : α) (p : CrfParams α) : α × α × α :=
  if r2 < guard then (0, 0, 0)
  else
    let invR2 := 1 / r2
    let invR6 := invR2 * invR2 * invR2
    let eLj := (c12 * invR6 - c6) * invR6
    let fLj := (12 * c12 * invR6 - 6 * c6) * invR6 * invR2
    let invR := Sqrt.sqrt invR2
    let eCrf := q * (invR - p.crf2Cut3i * r2 - p.crfCut)
    let fCrf := q * (invR * invR2 + p.crfCut3i)
    (fLj + fCrf, eLj, eCrf)

/-- `rf_excluded_pair(r2, q_prod, crf)`: `(f_mag, e)`. -/
def rfExcludedPairModel (r2 q : α) (p : CrfParams α) : α × α :=
  (q * p.crfCut3i, q * ((-p.crf2Cut3i) * r2 - p.crfCut))

/-- `rf_self_term(q, four_pi_eps_i, crf)`, with the literal `0.5` as `half`. -/
def rfSelfTermModel (half q fpe : α) (p : CrfParams α) : α :=
  (-half) * q * q * fpe * p.crfCut

end Model

/-! ## 2. The `Float` instantiation (what `lake exe spec-vectors` runs)

The two decimal literals are given by their bits, so that they are exactly Rust's; the Rust test
asserts `ZERO_DISTANCE2.to_bits()`. -/

/-- `ZERO_DISTANCE2 = 1e-10` (nm²), as `f64` bits. -/
def zeroDistance2Float : Float := Float.ofBits 0x3ddb7cdfd9d7bdbb

/-- The `0.5` of `rf_self_term`, as `f64` bits. -/
def halfFloat : Float := Float.ofBits 0x3fe0000000000000

/-- `CRFParameters::new` in binary64. -/
def crfParametersFloat (cutoff epsilon rfEpsilon rfKappa : Float) : CrfParams Float :=
  crfParametersModel Float cutoff epsilon rfEpsilon rfKappa

/-- `lj_crf_pair` in binary64. -/
def ljCrfPairFloat (r2 c6 c12 q : Float) (p : CrfParams Float) : Float × Float × Float :=
  ljCrfPairModel Float zeroDistance2Float r2 c6 c12 q p

/-- `rf_excluded_pair` in binary64. -/
def rfExcludedPairFloat (r2 q : Float) (p : CrfParams Float) : Float × Float :=
  rfExcludedPairModel Float r2 q p

/-- `rf_self_term` in binary64. -/
def rfSelfTermFloat (q fpe : Float) (p : CrfParams Float) : Float :=
  rfSelfTermModel Float halfFloat q fpe p

/-! ## 3. The ℝ instantiation equals the theory -/

/-- The derived constants over ℝ are the theory's: `(1 − C_RF/2)/R`, `C_RF/(2R³)`, `C_RF/R³`,
`R²`. -/
theorem crfDerivedModel_real (crf R : ℝ) :
    crfDerivedModel ℝ crf R = ⟨(1 - crf / 2) / R, crf / (2 * R ^ 3), crf / R ^ 3, R ^ 2⟩ := by
  simp only [crfDerivedModel, CrfParams.mk.injEq, true_and]
  and_intros <;> ring

/-- Outside the two sentinels (`R ≤ 0`: no reaction field; `ε_rf = 0`: conducting boundary) the
model's `C_RF` is the theory's `crfConstant`, and its constants are derived from it. -/
theorem crfParametersModel_real_eq {R ε εrf κ : ℝ} (hR : 0 < R) (hrf : εrf ≠ 0) :
    crfParametersModel ℝ R ε εrf κ = crfDerivedModel ℝ (crfConstant ε εrf κ R) R := by
  simp [crfParametersModel, crfConstantModel, not_le.mpr hR, hrf, crfConstant, pow_two]

/-- The conducting sentinel `ε_rf = 0` gives `C_RF = −1`, the `ε_rf → ∞` limit of
`crfConstant` (`crf_conducting_limit`). -/
theorem crfParametersModel_real_conducting {R ε κ : ℝ} (hR : 0 < R) :
    crfParametersModel ℝ R ε 0 κ = crfDerivedModel ℝ (-1) R := by
  simp [crfParametersModel, crfConstantModel, not_le.mpr hR]

/-- One pair at distance `r > 0` outside the guard, over ℝ: the model's three outputs are the
theory's scalar force `ljForce + crfForce`, `ljEnergy` and `crfEnergy`, for any `C_RF` and radius
`R ≠ 0` (with `Real.sqrt (1/r²) = 1/r`). -/
theorem ljCrfPairModel_real_eq {guard c6 c12 q crf R r : ℝ} (hr : 0 < r) (hR : R ≠ 0)
    (hg : guard ≤ r ^ 2) :
    ljCrfPairModel ℝ guard (r ^ 2) c6 c12 q (crfDerivedModel ℝ crf R) =
      (ljForce c6 c12 r + crfForce q crf R r, ljEnergy c6 c12 r, crfEnergy q crf R r) := by
  have hsqrt : Real.sqrt (1 / r ^ 2) = 1 / r := by
    rw [one_div, Real.sqrt_inv, Real.sqrt_sq hr.le, one_div]
  have hr0 : r ≠ 0 := hr.ne'
  simp only [ljCrfPairModel, not_lt.mpr hg, ite_false, crfDerivedModel, Sqrt.sqrt, hsqrt,
    ljForce, crfForce, ljEnergy, crfEnergy, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_⟩ <;> field_simp

/-- An excluded pair over ℝ: the model is the theory's `(rfExcludedForce, rfExcludedEnergy)`. -/
theorem rfExcludedPairModel_real_eq (q crf R r : ℝ) :
    rfExcludedPairModel ℝ (r ^ 2) q (crfDerivedModel ℝ crf R) =
      (rfExcludedForce q crf R, rfExcludedEnergy q crf R r) := by
  simp only [rfExcludedPairModel, crfDerivedModel_real, rfExcludedForce, rfExcludedEnergy,
    Prod.mk.injEq, true_and]
  ring

/-- The self term over ℝ, with `half = 1/2`: the theory's `rfSelfEnergy`, half the excluded-pair
energy of the atom with itself (`rf_self_is_half_excluded`). -/
theorem rfSelfTermModel_real_eq (qi fpe crf R : ℝ) :
    rfSelfTermModel ℝ (1 / 2) qi fpe (crfDerivedModel ℝ crf R) = rfSelfEnergy qi fpe crf R := by
  simp only [rfSelfTermModel, crfDerivedModel_real, rfSelfEnergy]
  ring

end NewenSpec.Nonbonded
