module

public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Inv
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.SpecificLimits.Basic
public import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# Coulomb interaction with a reaction field

GROMOS truncates electrostatics at the reaction-field radius `R` and replaces everything beyond it
by a continuum of relative permittivity `ε_rf` and inverse Debye length `κ`, polarised by the
charges inside (Barker & Watts 1973; Tironi, Sperb, Smith & van Gunsteren 1995, J. Chem. Phys.
102, 5451; GROMOS Vol. 2). Inside a medium of permittivity `ε`, a pair with charge product `q`
(which already carries `1/(4πε₀)`) at distance `r ≤ R` then has the energy

  `E(r) = q · (1/r − C_RF r²/(2R³) − (1 − C_RF/2)/R)`

with the reaction-field constant

  `C_RF = (2(ε − ε_rf)(1 + κR) − ε_rf (κR)²) / ((ε + 2ε_rf)(1 + κR) + ε_rf (κR)²)`.

The `r²` term is the reaction field of the continuum, the constant term shifts the energy to zero
at the cutoff (`crf_energy_zero_at_cutoff`). GROMOS writes `ε_rf = 0` for a conducting boundary,
the limit `ε_rf → ∞` where `C_RF = −1` (`crf_conducting_limit`, for every `κ ≥ 0`).

Excluded pairs (bonded neighbours) and each atom with itself get no Coulomb `1/r` term, but they
do polarise the continuum, so they keep the rest of the pair energy
(`rf_excluded_is_crf_minus_coulomb`).

As for Lennard-Jones (`NewenTheory/Interactions/LennardJones.lean`), the "force" of a pair is the
scalar `−E'(r)/r`, multiplied by `r_vec = x_i − x_j` to give the force on atom `i`.
-/

@[expose] public section

namespace NewenTheory.Interactions.ReactionField

open Filter Topology

/-- The reaction-field constant `C_RF` for permittivity `ε` inside, `ε_rf` outside, inverse Debye
length `κ` and radius `R`. -/
noncomputable def crfConstant (ε εrf κ R : ℝ) : ℝ :=
  (2 * (ε - εrf) * (1 + κ * R) - εrf * (κ * R) ^ 2) /
    ((ε + 2 * εrf) * (1 + κ * R) + εrf * (κ * R) ^ 2)

/-- The reaction-field pair energy for charge product `q`, constant `crf`, radius `R`, distance
`r`. -/
noncomputable def crfEnergy (q crf R r : ℝ) : ℝ :=
  q * (1 / r - crf * r ^ 2 / (2 * R ^ 3) - (1 - crf / 2) / R)

/-- The scalar reaction-field force `−E'(r)/r`. -/
noncomputable def crfForce (q crf R r : ℝ) : ℝ := q * (1 / r ^ 3 + crf / R ^ 3)

/-- The energy of an excluded pair: the reaction-field part of the pair energy, without `1/r`. -/
noncomputable def rfExcludedEnergy (q crf R r : ℝ) : ℝ :=
  q * (-(crf * r ^ 2 / (2 * R ^ 3)) - (1 - crf / 2) / R)

/-- The scalar force of an excluded pair. -/
noncomputable def rfExcludedForce (q crf R : ℝ) : ℝ := q * (crf / R ^ 3)

/-- The self term of an atom with charge `qi`, with `fpe = 1/(4πε₀)`. -/
noncomputable def rfSelfEnergy (qi fpe crf R : ℝ) : ℝ :=
  -(1 / 2) * (qi ^ 2 * fpe) * ((1 - crf / 2) / R)

/-- **Force = −dE/dr.** The derivative of the energy at `r ≠ 0` is `−crfForce · r`. -/
theorem crf_force_is_minus_grad (q crf R : ℝ) {r : ℝ} (hr : r ≠ 0) :
    HasDerivAt (crfEnergy q crf R) (-(crfForce q crf R r) * r) r := by
  have h : HasDerivAt (fun x => q * (x⁻¹ - crf / (2 * R ^ 3) * x ^ 2 - (1 - crf / 2) / R))
      (q * (-(r ^ 2)⁻¹ - crf / (2 * R ^ 3) * (↑2 * r ^ (2 - 1)))) r :=
    HasDerivAt.const_mul q (HasDerivAt.sub_const ((1 - crf / 2) / R)
      (HasDerivAt.sub (hasDerivAt_inv hr)
        (HasDerivAt.const_mul (crf / (2 * R ^ 3)) (hasDerivAt_pow 2 r))))
  have hE : crfEnergy q crf R =
      fun x => q * (x⁻¹ - crf / (2 * R ^ 3) * x ^ 2 - (1 - crf / 2) / R) := by
    funext x
    simp only [crfEnergy, one_div]
    ring
  rw [hE]
  refine h.congr_deriv ?_
  simp only [crfForce]
  field_simp
  ring

/-- The pair energy is zero at the cutoff. -/
theorem crf_energy_zero_at_cutoff (q crf : ℝ) {R : ℝ} (hR : R ≠ 0) :
    crfEnergy q crf R R = 0 := by
  simp only [crfEnergy]
  field_simp
  ring

/-- Without ionic strength (`κ = 0`) and with the same permittivity inside and out, there is no
reaction field: `C_RF = 0`, and the energy is the Coulomb energy shifted to zero at `R`. (With
`κ ≠ 0`, `ε_rf = ε` still leaves `C_RF = −ε_rf (κR)² / (3ε_rf + ε_rf (κR)²)`.) -/
theorem crf_eps_rf_eq_eps (ε R : ℝ) :
    crfConstant ε ε 0 R = 0 ∧
      ∀ q r, crfEnergy q (crfConstant ε ε 0 R) R r = q * (1 / r - 1 / R) := by
  have h0 : crfConstant ε ε 0 R = 0 := by simp [crfConstant]
  refine ⟨h0, fun q r => ?_⟩
  rw [h0]
  simp [crfEnergy]

/-- The excluded-pair energy is the pair energy without its Coulomb term `q/r`, and the
excluded-pair force the pair force without `q/r³`. -/
theorem rf_excluded_is_crf_minus_coulomb (q crf R r : ℝ) :
    rfExcludedEnergy q crf R r = crfEnergy q crf R r - q / r ∧
      rfExcludedForce q crf R = crfForce q crf R r - q / r ^ 3 := by
  constructor
  · simp only [rfExcludedEnergy, crfEnergy]
    ring
  · simp only [rfExcludedForce, crfForce]
    ring

/-- The self term is half the excluded-pair energy of the atom with itself (`r = 0`, charge
product `qi²·fpe`): the half because the pair `(i, i)` is one atom, not two. -/
theorem rf_self_is_half_excluded (qi fpe crf R : ℝ) :
    rfSelfEnergy qi fpe crf R = rfExcludedEnergy (qi ^ 2 * fpe) crf R 0 / 2 := by
  simp only [rfSelfEnergy, rfExcludedEnergy]
  ring

/-- The excluded-pair force is the derivative of the excluded-pair energy, as for the full pair:
`HasDerivAt E (−rfExcludedForce · r) r`, at every `r`. -/
theorem rf_excluded_force_is_minus_grad (q crf R r : ℝ) :
    HasDerivAt (rfExcludedEnergy q crf R) (-(rfExcludedForce q crf R) * r) r := by
  have h : HasDerivAt (fun x => q * (-(crf / (2 * R ^ 3) * x ^ 2) - (1 - crf / 2) / R))
      (q * -(crf / (2 * R ^ 3) * (↑2 * r ^ (2 - 1)))) r :=
    HasDerivAt.const_mul q (HasDerivAt.sub_const ((1 - crf / 2) / R)
      (HasDerivAt.neg (HasDerivAt.const_mul (crf / (2 * R ^ 3)) (hasDerivAt_pow 2 r))))
  have hE : rfExcludedEnergy q crf R =
      fun x => q * (-(crf / (2 * R ^ 3) * x ^ 2) - (1 - crf / 2) / R) := by
    funext x
    simp only [rfExcludedEnergy]
    ring
  rw [hE]
  refine h.congr_deriv ?_
  simp only [rfExcludedForce]
  by_cases hR : R = 0
  · simp [hR]
  · field_simp
    ring

/-- A conducting boundary: for every ionic strength `κ ≥ 0` and radius `R > 0`, `C_RF → −1` as
`ε_rf → ∞`. This is the value GROMOS assigns to the sentinel `ε_rf = 0`. -/
theorem crf_conducting_limit (ε : ℝ) {κ R : ℝ} (hκ : 0 ≤ κ) (hR : 0 < R) :
    Tendsto (fun εrf => crfConstant ε εrf κ R) atTop (𝓝 (-1)) := by
  -- numerator `2 ε a − ε_rf c`, denominator `ε a + ε_rf c`, with `a = 1 + κR`, `c = 2a + (κR)²`
  set a := 1 + κ * R
  set c := 2 * a + (κ * R) ^ 2
  have hc : 0 < c := by have : 0 ≤ κ * R := mul_nonneg hκ hR.le; positivity
  have hinv : Tendsto (fun x : ℝ => x⁻¹) atTop (𝓝 0) := tendsto_inv_atTop_zero
  have hn : Tendsto (fun x : ℝ => 2 * ε * a * x⁻¹ - c) atTop (𝓝 (2 * ε * a * 0 - c)) :=
    (hinv.const_mul (2 * ε * a)).sub_const c
  have hd : Tendsto (fun x : ℝ => ε * a * x⁻¹ + c) atTop (𝓝 (ε * a * 0 + c)) :=
    (hinv.const_mul (ε * a)).add_const c
  have hlim := hn.div hd (by simp only [mul_zero, zero_add]; exact hc.ne')
  rw [show (2 * ε * a * 0 - c) / (ε * a * 0 + c) = (-1 : ℝ) by
    rw [mul_zero, mul_zero, zero_sub, zero_add, neg_div, div_self hc.ne']] at hlim
  refine hlim.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with x hx
  have hx' : x ≠ 0 := hx.ne'
  have hnum : 2 * ε * a * x⁻¹ - c = (2 * (ε - x) * a - x * (κ * R) ^ 2) / x := by
    simp only [c]; field_simp; ring
  have hden : ε * a * x⁻¹ + c = ((ε + 2 * x) * a + x * (κ * R) ^ 2) / x := by
    simp only [c]; field_simp; ring
  rw [Pi.div_apply, hnum, hden, div_div_div_cancel_right₀ hx']
  rfl

end NewenTheory.Interactions.ReactionField

end
