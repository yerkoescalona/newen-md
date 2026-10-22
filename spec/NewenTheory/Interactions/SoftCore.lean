module

public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Comp
public import Mathlib.Analysis.Calculus.Deriv.Inv
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.Calculus.Deriv.Pow
public import Mathlib.Analysis.SpecialFunctions.Sqrt
public import NewenTheory.Interactions.LennardJones
public import NewenTheory.Interactions.ReactionField
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# Soft-core Lennard-Jones and reaction-field interactions of a perturbed pair

A free-energy calculation couples two end states A and B of a pair through a parameter `λ`
(GROMOS Vol. 2; Beutler, Mark, van Schaik, Gerber & van Gunsteren 1994, Chem. Phys. Lett. 222,
529). Each state keeps its own Lennard-Jones parameters `C6`, `C12` and charge product `q`, and the
pair energy is `V(λ) = (1 − λ)ⁿ V_A(λ) + λⁿ V_B(λ)` with `n` = NLAM. Where one state is a dummy
(no Lennard-Jones, no charge) the other state's energy would diverge at `r → 0` before its
weight vanishes, so each state is *softened* by the other state's `λ`:

* Lennard-Jones: `r⁶` becomes `r⁶ + α_LJ · s · C12/C6`, with `s = λ²` for A and `s = (1 − λ)²`
  for B;
* reaction field: `r²` becomes `r² + α_CRF · s`, in the Coulomb term and in the reaction-field
  term, where the cutoff `R²` is softened the same way (`R² + α_CRF · s`); the constant shift
  `(1 − C_RF/2)/R` is not.

At `λ = 0` state A is unsoftened and has weight 1, so the pair is the plain `ljEnergy` and
`crfEnergy` of `LennardJones.lean` and `ReactionField.lean` (`softcore_lambda_zero`); at `λ = 1`
the same holds for B (`softcore_lambda_one`).

Thermodynamic integration needs `dV/dλ`. The kernels compute it in closed form, and
`softcore_dvdl_is_deriv_lj` / `softcore_dvdl_is_deriv_crf` prove that the closed form is the
derivative of the energy. They are stated for *individual* λ (GROMOS LAMBDAS): the weights take
`Λ(λ)`, the softness `Λₛ(λ)`, two differentiable functions of the global λ, and `dV/dλ` is the chain
rule `∂V/∂Λ · Λ′ + ∂V/∂Λₛ · Λₛ′`; with `Λ = Λₛ = λ` it is the plain derivative
(`softcore_dvdl_is_deriv_lj_plain`, `softcore_dvdl_is_deriv_crf_plain`).

As in the other interaction files, the force of a pair is the scalar `−∂V/∂r / r`
(`softcore_force_is_minus_grad`). Everything here is over ℝ; how newen computes it in binary64 is
in `NewenSpec/SoftCore.lean`.
-/

@[expose] public section

namespace NewenTheory.Interactions.SoftCore

open NewenTheory.Interactions.LennardJones NewenTheory.Interactions.ReactionField

/-! ## One end state, softened by `s` -/

/-- The Lennard-Jones parameters of one end state: `C6`, `C12`, and the soft-core length
`c126 = C12/C6` (or `0` without dispersion). -/
structure LjState where
  c6 : ℝ
  c12 : ℝ
  c126 : ℝ

/-- The soft-core Lennard-Jones energy of one state at softness `s`. -/
noncomputable def scLjEnergy (αlj : ℝ) (p : LjState) (s r : ℝ) : ℝ :=
  p.c12 / (r ^ 6 + αlj * s * p.c126) ^ 2 - p.c6 / (r ^ 6 + αlj * s * p.c126)

/-- `∂/∂s` of `scLjEnergy`. -/
noncomputable def scLjDs (αlj : ℝ) (p : LjState) (s r : ℝ) : ℝ :=
  -(αlj * p.c126) *
    (2 * p.c12 / (r ^ 6 + αlj * s * p.c126) ^ 3 - p.c6 / (r ^ 6 + αlj * s * p.c126) ^ 2)

/-- The scalar soft-core Lennard-Jones force `−∂E/∂r / r` of one state. -/
noncomputable def scLjForce (αlj : ℝ) (p : LjState) (s r : ℝ) : ℝ :=
  (12 * p.c12 / (r ^ 6 + αlj * s * p.c126) ^ 3 - 6 * p.c6 / (r ^ 6 + αlj * s * p.c126) ^ 2) *
    r ^ 4

/-- The soft-core reaction-field energy of one state with charge product `q` at softness `s`;
`cs` scales the Coulomb term (`1` except for 1-4 pairs). -/
noncomputable def scCrfEnergy (αc cs q crf R s r : ℝ) : ℝ :=
  q * (cs / Real.sqrt (r ^ 2 + αc * s) - crf / 2 / Real.sqrt ((R ^ 2 + αc * s) ^ 3) * r ^ 2 -
    (1 - crf / 2) / R)

/-- `∂/∂s` of `scCrfEnergy`. -/
noncomputable def scCrfDs (αc cs q crf R s r : ℝ) : ℝ :=
  -(αc / 2) * q *
    (cs / Real.sqrt (r ^ 2 + αc * s) / (r ^ 2 + αc * s) -
      3 * (crf / 2 / Real.sqrt ((R ^ 2 + αc * s) ^ 3)) / (R ^ 2 + αc * s) * r ^ 2)

/-- The scalar soft-core reaction-field force `−∂E/∂r / r` of one state. -/
noncomputable def scCrfForce (αc cs q crf R s r : ℝ) : ℝ :=
  q * (cs / Real.sqrt (r ^ 2 + αc * s) / (r ^ 2 + αc * s) +
    crf / Real.sqrt ((R ^ 2 + αc * s) ^ 3))

theorem scLj_hasDerivAt_s (αlj : ℝ) (p : LjState) (r : ℝ) {s : ℝ}
    (hX : r ^ 6 + αlj * s * p.c126 ≠ 0) :
    HasDerivAt (fun s => scLjEnergy αlj p s r) (scLjDs αlj p s r) s := by
  have hX' : HasDerivAt (fun s => r ^ 6 + αlj * s * p.c126) (αlj * p.c126) s := by
    simpa using (((hasDerivAt_id s).const_mul αlj).mul_const p.c126).const_add (r ^ 6)
  have h := ((hX'.pow 2).inv (pow_ne_zero 2 hX)).const_mul p.c12 |>.sub
    ((hX'.inv hX).const_mul p.c6)
  have hE : (fun s => scLjEnergy αlj p s r) =
      fun s => p.c12 * ((r ^ 6 + αlj * s * p.c126) ^ 2)⁻¹ - p.c6 * (r ^ 6 + αlj * s * p.c126)⁻¹ := by
    funext s
    simp [scLjEnergy, div_eq_mul_inv]
  rw [hE]
  refine h.congr_deriv ?_
  simp only [scLjDs, Pi.pow_apply]
  field_simp
  ring

theorem scLj_force_is_minus_grad (αlj : ℝ) (p : LjState) (s : ℝ) {r : ℝ}
    (hX : r ^ 6 + αlj * s * p.c126 ≠ 0) :
    HasDerivAt (fun r => scLjEnergy αlj p s r) (-(scLjForce αlj p s r) * r) r := by
  have hX' : HasDerivAt (fun r => r ^ 6 + αlj * s * p.c126) (6 * r ^ 5) r := by
    simpa using (hasDerivAt_pow 6 r).add_const (αlj * s * p.c126)
  have h := ((hX'.pow 2).inv (pow_ne_zero 2 hX)).const_mul p.c12 |>.sub
    ((hX'.inv hX).const_mul p.c6)
  have hE : (fun r => scLjEnergy αlj p s r) =
      fun r => p.c12 * ((r ^ 6 + αlj * s * p.c126) ^ 2)⁻¹ - p.c6 * (r ^ 6 + αlj * s * p.c126)⁻¹ := by
    funext r
    simp [scLjEnergy, div_eq_mul_inv]
  rw [hE]
  refine h.congr_deriv ?_
  simp only [scLjForce, Pi.pow_apply]
  field_simp
  ring

/-- `1/√(f y)` has derivative `−f′/(2 f √f)` where `f > 0`. -/
theorem hasDerivAt_inv_sqrt {f : ℝ → ℝ} {f' x : ℝ} (hf : HasDerivAt f f' x) (h0 : 0 < f x) :
    HasDerivAt (fun y => 1 / Real.sqrt (f y)) (-(f' / (2 * (Real.sqrt (f x) * f x)))) x := by
  have hs : Real.sqrt (f x) ≠ 0 := (Real.sqrt_pos.mpr h0).ne'
  have h := (hf.sqrt h0.ne').inv hs
  have hE : (fun y => 1 / Real.sqrt (f y)) = fun y => (Real.sqrt (f y))⁻¹ := by
    funext y
    simp
  rw [hE]
  refine h.congr_deriv ?_
  rw [Real.sq_sqrt h0.le]
  field_simp

theorem scCrf_hasDerivAt_s (αc cs q crf R r : ℝ) {s : ℝ} (hX : 0 < r ^ 2 + αc * s)
    (hY : 0 < R ^ 2 + αc * s) :
    HasDerivAt (fun s => scCrfEnergy αc cs q crf R s r) (scCrfDs αc cs q crf R s r) s := by
  have hX' : HasDerivAt (fun s => r ^ 2 + αc * s) αc s := by
    simpa using ((hasDerivAt_id s).const_mul αc).const_add (r ^ 2)
  have hY' : HasDerivAt (fun s => R ^ 2 + αc * s) αc s := by
    simpa using ((hasDerivAt_id s).const_mul αc).const_add (R ^ 2)
  have hY3 : 0 < (R ^ 2 + αc * s) ^ 3 := pow_pos hY 3
  have h1 := (hasDerivAt_inv_sqrt hX' hX).const_mul cs
  have h2 := (hasDerivAt_inv_sqrt (hY'.pow 3) hY3).const_mul (crf / 2 * r ^ 2)
  have h := ((h1.sub h2).sub_const ((1 - crf / 2) / R)).const_mul q
  have hE : (fun s => scCrfEnergy αc cs q crf R s r) = fun s =>
      q * (cs * (1 / Real.sqrt (r ^ 2 + αc * s)) -
        crf / 2 * r ^ 2 * (1 / Real.sqrt ((fun s => (R ^ 2 + αc * s)) s ^ 3)) -
        (1 - crf / 2) / R) := by
    funext s
    simp only [scCrfEnergy]
    ring
  rw [hE]
  refine h.congr_deriv ?_
  have hsX : Real.sqrt (r ^ 2 + αc * s) ≠ 0 := (Real.sqrt_pos.mpr hX).ne'
  have hsY : Real.sqrt ((R ^ 2 + αc * s) ^ 3) ≠ 0 := (Real.sqrt_pos.mpr hY3).ne'
  simp only [scCrfDs, Pi.pow_apply]
  field_simp
  ring

theorem scCrf_force_is_minus_grad (αc cs q crf R s : ℝ) {r : ℝ} (hX : 0 < r ^ 2 + αc * s) :
    HasDerivAt (fun r => scCrfEnergy αc cs q crf R s r) (-(scCrfForce αc cs q crf R s r) * r) r := by
  have hX' : HasDerivAt (fun r => r ^ 2 + αc * s) (2 * r) r := by
    simpa using (hasDerivAt_pow 2 r).add_const (αc * s)
  have h1 := (hasDerivAt_inv_sqrt hX' hX).const_mul cs
  have h2 := (hasDerivAt_pow 2 r).const_mul (crf / 2 / Real.sqrt ((R ^ 2 + αc * s) ^ 3))
  have h := ((h1.sub h2).sub_const ((1 - crf / 2) / R)).const_mul q
  have hE : (fun r => scCrfEnergy αc cs q crf R s r) = fun r =>
      q * (cs * (1 / Real.sqrt (r ^ 2 + αc * s)) -
        crf / 2 / Real.sqrt ((R ^ 2 + αc * s) ^ 3) * r ^ 2 - (1 - crf / 2) / R) := by
    funext r
    simp only [scCrfEnergy]
    ring
  rw [hE]
  refine h.congr_deriv ?_
  have hsX : Real.sqrt (r ^ 2 + αc * s) ≠ 0 := (Real.sqrt_pos.mpr hX).ne'
  simp only [scCrfForce]
  by_cases hY : Real.sqrt ((R ^ 2 + αc * s) ^ 3) = 0
  · simp only [hY, div_zero]
    field_simp
    ring
  · field_simp
    ring

/-- The excluded-pair reaction-field energy of one state at softness `s`: the reaction-field
part without the Coulomb term (`rfExcludedEnergy`, softened). -/
noncomputable def scRfExclEnergy (αc q crf R s r : ℝ) : ℝ :=
  q * (-(crf / 2 / Real.sqrt ((R ^ 2 + αc * s) ^ 3)) * r ^ 2 - (1 - crf / 2) / R)

/-- `∂/∂s` of `scRfExclEnergy`. -/
noncomputable def scRfExclDs (αc q crf R s r : ℝ) : ℝ :=
  (αc / 2) * q * (3 * (crf / 2 / Real.sqrt ((R ^ 2 + αc * s) ^ 3)) / (R ^ 2 + αc * s) * r ^ 2)

theorem scRfExcl_hasDerivAt_s (αc q crf R r : ℝ) {s : ℝ} (hY : 0 < R ^ 2 + αc * s) :
    HasDerivAt (fun s => scRfExclEnergy αc q crf R s r) (scRfExclDs αc q crf R s r) s := by
  have hY' : HasDerivAt (fun s => R ^ 2 + αc * s) αc s := by
    simpa using ((hasDerivAt_id s).const_mul αc).const_add (R ^ 2)
  have hY3 : 0 < (R ^ 2 + αc * s) ^ 3 := pow_pos hY 3
  have h2 := (hasDerivAt_inv_sqrt (hY'.pow 3) hY3).const_mul (crf / 2 * r ^ 2)
  have h := ((h2.neg).sub_const ((1 - crf / 2) / R)).const_mul q
  have hE : (fun s => scRfExclEnergy αc q crf R s r) = fun s =>
      q * (-(crf / 2 * r ^ 2 * (1 / Real.sqrt ((fun s => (R ^ 2 + αc * s)) s ^ 3))) -
        (1 - crf / 2) / R) := by
    funext s
    simp only [scRfExclEnergy]
    ring
  rw [hE]
  refine h.congr_deriv ?_
  have hsY : Real.sqrt ((R ^ 2 + αc * s) ^ 3) ≠ 0 := (Real.sqrt_pos.mpr hY3).ne'
  simp only [scRfExclDs, Pi.pow_apply]
  field_simp
  ring

/-! ## The perturbed pair: weights `(1 − Λ)ⁿ`, `Λⁿ`, softness from `Λₛ` -/

/-- The Lennard-Jones energy of a perturbed pair at weight λ `L` and softness λ `Ls`. -/
noncomputable def pertLjEnergy (n : ℕ) (αlj : ℝ) (A B : LjState) (r L Ls : ℝ) : ℝ :=
  (1 - L) ^ n * scLjEnergy αlj A (Ls ^ 2) r + L ^ n * scLjEnergy αlj B ((1 - Ls) ^ 2) r

/-- Its `dV/dλ` when `L` and `Ls` move with `dL`, `dLs` per unit λ. -/
noncomputable def pertLjDvdl (n : ℕ) (αlj : ℝ) (A B : LjState) (r L Ls dL dLs : ℝ) : ℝ :=
  (1 - L) ^ n * (scLjDs αlj A (Ls ^ 2) r * (2 * Ls * dLs)) +
    L ^ n * (scLjDs αlj B ((1 - Ls) ^ 2) r * (-(2 * (1 - Ls) * dLs))) +
    n * (L ^ (n - 1) * scLjEnergy αlj B ((1 - Ls) ^ 2) r -
      (1 - L) ^ (n - 1) * scLjEnergy αlj A (Ls ^ 2) r) * dL

/-- The scalar Lennard-Jones force of a perturbed pair. -/
noncomputable def pertLjForce (n : ℕ) (αlj : ℝ) (A B : LjState) (r L Ls : ℝ) : ℝ :=
  (1 - L) ^ n * scLjForce αlj A (Ls ^ 2) r + L ^ n * scLjForce αlj B ((1 - Ls) ^ 2) r

/-- The reaction-field energy of a perturbed pair with charge products `qA`, `qB` (without
`1/(4πε₀)`, which is `fpe`). -/
noncomputable def pertCrfEnergy (n : ℕ) (αc cs fpe qA qB crf R r L Ls : ℝ) : ℝ :=
  ((1 - L) ^ n * scCrfEnergy αc cs qA crf R (Ls ^ 2) r +
    L ^ n * scCrfEnergy αc cs qB crf R ((1 - Ls) ^ 2) r) * fpe

/-- Its `dV/dλ`. -/
noncomputable def pertCrfDvdl (n : ℕ) (αc cs fpe qA qB crf R r L Ls dL dLs : ℝ) : ℝ :=
  ((1 - L) ^ n * (scCrfDs αc cs qA crf R (Ls ^ 2) r * (2 * Ls * dLs)) +
    L ^ n * (scCrfDs αc cs qB crf R ((1 - Ls) ^ 2) r * (-(2 * (1 - Ls) * dLs))) +
    n * (L ^ (n - 1) * scCrfEnergy αc cs qB crf R ((1 - Ls) ^ 2) r -
      (1 - L) ^ (n - 1) * scCrfEnergy αc cs qA crf R (Ls ^ 2) r) * dL) * fpe

/-- The scalar reaction-field force of a perturbed pair. -/
noncomputable def pertCrfForce (n : ℕ) (αc cs fpe qA qB crf R r L Ls : ℝ) : ℝ :=
  ((1 - L) ^ n * scCrfForce αc cs qA crf R (Ls ^ 2) r +
    L ^ n * scCrfForce αc cs qB crf R ((1 - Ls) ^ 2) r) * fpe

/-- The reaction-field energy of an excluded perturbed pair. -/
noncomputable def pertRfExclEnergy (n : ℕ) (αc fpe qA qB crf R r L Ls : ℝ) : ℝ :=
  ((1 - L) ^ n * scRfExclEnergy αc qA crf R (Ls ^ 2) r +
    L ^ n * scRfExclEnergy αc qB crf R ((1 - Ls) ^ 2) r) * fpe

/-- Its `dV/dλ`. -/
noncomputable def pertRfExclDvdl (n : ℕ) (αc fpe qA qB crf R r L Ls dL dLs : ℝ) : ℝ :=
  ((1 - L) ^ n * (scRfExclDs αc qA crf R (Ls ^ 2) r * (2 * Ls * dLs)) +
    L ^ n * (scRfExclDs αc qB crf R ((1 - Ls) ^ 2) r * (-(2 * (1 - Ls) * dLs))) +
    n * (L ^ (n - 1) * scRfExclEnergy αc qB crf R ((1 - Ls) ^ 2) r -
      (1 - L) ^ (n - 1) * scRfExclEnergy αc qA crf R (Ls ^ 2) r) * dL) * fpe

/-- The scalar reaction-field force of an excluded perturbed pair (it does not depend on `r`). -/
noncomputable def pertRfExclForce (n : ℕ) (αc fpe qA qB crf R L Ls : ℝ) : ℝ :=
  ((1 - L) ^ n * (qA * (crf / Real.sqrt ((R ^ 2 + αc * Ls ^ 2) ^ 3))) +
    L ^ n * (qB * (crf / Real.sqrt ((R ^ 2 + αc * (1 - Ls) ^ 2) ^ 3)))) * fpe

/-- The reaction-field self term of a perturbed atom with squared charges `qa2`, `qb2`, as the
change from the state-A self term, times `2` (the kernel books half of it): the excluded-pair
energy at `r = 0`, where the softness drops out, weighted, minus state A's. `c = (1 − C_RF/2)/R ·
fpe`. -/
noncomputable def rfSelfChange (n : ℕ) (qa2 qb2 c L : ℝ) : ℝ :=
  (qa2 - ((1 - L) ^ n * qa2 + L ^ n * qb2)) * c

/-! ## `dV/dλ` is the derivative of `V` -/

/-- The weight `(1 − Λ)ⁿ` and its derivative. -/
theorem hasDerivAt_weight_a (n : ℕ) {Λ : ℝ → ℝ} {dΛ lam : ℝ} (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => (1 - Λ l) ^ n) (-(n * (1 - Λ lam) ^ (n - 1) * dΛ)) lam := by
  have h := (hΛ.const_sub 1).pow n
  refine h.congr_deriv ?_
  ring

/-- The weight `Λⁿ` and its derivative. -/
theorem hasDerivAt_weight_b (n : ℕ) {Λ : ℝ → ℝ} {dΛ lam : ℝ} (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => Λ l ^ n) (n * Λ lam ^ (n - 1) * dΛ) lam :=
  hΛ.pow n

/-- The softness `Λₛ²` of state A. -/
theorem hasDerivAt_soft_a {Λs : ℝ → ℝ} {dΛs lam : ℝ} (hΛs : HasDerivAt Λs dΛs lam) :
    HasDerivAt (fun l => Λs l ^ 2) (2 * Λs lam * dΛs) lam := by
  have h := hΛs.pow 2
  refine h.congr_deriv ?_
  simp

/-- The softness `(1 − Λₛ)²` of state B. -/
theorem hasDerivAt_soft_b {Λs : ℝ → ℝ} {dΛs lam : ℝ} (hΛs : HasDerivAt Λs dΛs lam) :
    HasDerivAt (fun l => (1 - Λs l) ^ 2) (-(2 * (1 - Λs lam) * dΛs)) lam := by
  have h := (hΛs.const_sub 1).pow 2
  refine h.congr_deriv ?_
  simp

/-- **dV/dλ of the soft-core Lennard-Jones term is the derivative of its energy**, for
individual λ `Λ` (weights) and `Λₛ` (softness) that are differentiable in λ, wherever both
softened `r⁶` are nonzero. -/
theorem softcore_dvdl_is_deriv_lj (n : ℕ) (αlj : ℝ) (A B : LjState) (r : ℝ)
    {Λ Λs : ℝ → ℝ} {dΛ dΛs lam : ℝ} (hΛ : HasDerivAt Λ dΛ lam) (hΛs : HasDerivAt Λs dΛs lam)
    (hA : r ^ 6 + αlj * Λs lam ^ 2 * A.c126 ≠ 0)
    (hB : r ^ 6 + αlj * (1 - Λs lam) ^ 2 * B.c126 ≠ 0) :
    HasDerivAt (fun l => pertLjEnergy n αlj A B r (Λ l) (Λs l))
      (pertLjDvdl n αlj A B r (Λ lam) (Λs lam) dΛ dΛs) lam := by
  have hEA := (scLj_hasDerivAt_s αlj A r hA).comp lam (hasDerivAt_soft_a hΛs)
  have hEB := (scLj_hasDerivAt_s αlj B r hB).comp lam (hasDerivAt_soft_b hΛs)
  have h := ((hasDerivAt_weight_a n hΛ).mul hEA).add ((hasDerivAt_weight_b n hΛ).mul hEB)
  refine h.congr_deriv ?_
  simp only [pertLjDvdl, Function.comp_apply]
  ring

/-- **dV/dλ of the soft-core reaction-field term is the derivative of its energy**, for
individual λ, wherever the softened `r²` and `R²` of both states are positive. -/
theorem softcore_dvdl_is_deriv_crf (n : ℕ) (αc cs fpe qA qB crf R r : ℝ)
    {Λ Λs : ℝ → ℝ} {dΛ dΛs lam : ℝ} (hΛ : HasDerivAt Λ dΛ lam) (hΛs : HasDerivAt Λs dΛs lam)
    (hXA : 0 < r ^ 2 + αc * Λs lam ^ 2) (hYA : 0 < R ^ 2 + αc * Λs lam ^ 2)
    (hXB : 0 < r ^ 2 + αc * (1 - Λs lam) ^ 2) (hYB : 0 < R ^ 2 + αc * (1 - Λs lam) ^ 2) :
    HasDerivAt (fun l => pertCrfEnergy n αc cs fpe qA qB crf R r (Λ l) (Λs l))
      (pertCrfDvdl n αc cs fpe qA qB crf R r (Λ lam) (Λs lam) dΛ dΛs) lam := by
  have hEA := (scCrf_hasDerivAt_s αc cs qA crf R r hXA hYA).comp lam (hasDerivAt_soft_a hΛs)
  have hEB := (scCrf_hasDerivAt_s αc cs qB crf R r hXB hYB).comp lam (hasDerivAt_soft_b hΛs)
  have h := (((hasDerivAt_weight_a n hΛ).mul hEA).add
    ((hasDerivAt_weight_b n hΛ).mul hEB)).mul_const fpe
  refine h.congr_deriv ?_
  simp only [pertCrfDvdl, Function.comp_apply]
  ring

/-- The same for the reaction-field term of an excluded pair. -/
theorem softcore_dvdl_is_deriv_rf_excluded (n : ℕ) (αc fpe qA qB crf R r : ℝ)
    {Λ Λs : ℝ → ℝ} {dΛ dΛs lam : ℝ} (hΛ : HasDerivAt Λ dΛ lam) (hΛs : HasDerivAt Λs dΛs lam)
    (hYA : 0 < R ^ 2 + αc * Λs lam ^ 2) (hYB : 0 < R ^ 2 + αc * (1 - Λs lam) ^ 2) :
    HasDerivAt (fun l => pertRfExclEnergy n αc fpe qA qB crf R r (Λ l) (Λs l))
      (pertRfExclDvdl n αc fpe qA qB crf R r (Λ lam) (Λs lam) dΛ dΛs) lam := by
  have hEA := (scRfExcl_hasDerivAt_s αc qA crf R r hYA).comp lam (hasDerivAt_soft_a hΛs)
  have hEB := (scRfExcl_hasDerivAt_s αc qB crf R r hYB).comp lam (hasDerivAt_soft_b hΛs)
  have h := (((hasDerivAt_weight_a n hΛ).mul hEA).add
    ((hasDerivAt_weight_b n hΛ).mul hEB)).mul_const fpe
  refine h.congr_deriv ?_
  simp only [pertRfExclDvdl, Function.comp_apply]
  ring

/-- And for the self term: `d/dλ` of `rfSelfChange` is `n((1 − Λ)ⁿ⁻¹ qa2 − Λⁿ⁻¹ qb2) c Λ′`. -/
theorem rf_self_dvdl_is_deriv (n : ℕ) (qa2 qb2 c : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => rfSelfChange n qa2 qb2 c (Λ l))
      (n * ((1 - Λ lam) ^ (n - 1) * qa2 - Λ lam ^ (n - 1) * qb2) * c * dΛ) lam := by
  have h := ((((hasDerivAt_weight_a n hΛ).mul_const qa2).add
    ((hasDerivAt_weight_b n hΛ).mul_const qb2)).const_sub qa2).mul_const c
  refine h.congr_deriv ?_
  ring

/-- With one λ for weights and softness (`Λ = Λₛ = λ`), the Lennard-Jones `dV/dλ` is the plain
derivative. -/
theorem softcore_dvdl_is_deriv_lj_plain (n : ℕ) (αlj : ℝ) (A B : LjState) (r lam : ℝ)
    (hA : r ^ 6 + αlj * lam ^ 2 * A.c126 ≠ 0) (hB : r ^ 6 + αlj * (1 - lam) ^ 2 * B.c126 ≠ 0) :
    HasDerivAt (fun l => pertLjEnergy n αlj A B r l l) (pertLjDvdl n αlj A B r lam lam 1 1) lam :=
  softcore_dvdl_is_deriv_lj n αlj A B r (hasDerivAt_id lam) (hasDerivAt_id lam) hA hB

/-- With one λ, the reaction-field `dV/dλ` is the plain derivative. -/
theorem softcore_dvdl_is_deriv_crf_plain (n : ℕ) (αc cs fpe qA qB crf R r lam : ℝ)
    (hXA : 0 < r ^ 2 + αc * lam ^ 2) (hYA : 0 < R ^ 2 + αc * lam ^ 2)
    (hXB : 0 < r ^ 2 + αc * (1 - lam) ^ 2) (hYB : 0 < R ^ 2 + αc * (1 - lam) ^ 2) :
    HasDerivAt (fun l => pertCrfEnergy n αc cs fpe qA qB crf R r l l)
      (pertCrfDvdl n αc cs fpe qA qB crf R r lam lam 1 1) lam :=
  softcore_dvdl_is_deriv_crf n αc cs fpe qA qB crf R r (hasDerivAt_id lam) (hasDerivAt_id lam)
    hXA hYA hXB hYB

/-- **The central theorem**: for one perturbed pair, the kernel's `dV/dλ` (Lennard-Jones plus
reaction field) is the derivative of its energy, for individual λ `Λ_LJ`, `Λ_LJs`, `Λ_CRF`,
`Λ_CRFs`. -/
theorem softcore_dvdl_is_deriv (n : ℕ) (αlj αc cs fpe : ℝ) (A B : LjState) (qA qB crf R r : ℝ)
    {Λ Λs Λc Λcs : ℝ → ℝ} {dΛ dΛs dΛc dΛcs lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) (hΛs : HasDerivAt Λs dΛs lam)
    (hΛc : HasDerivAt Λc dΛc lam) (hΛcs : HasDerivAt Λcs dΛcs lam)
    (hA : r ^ 6 + αlj * Λs lam ^ 2 * A.c126 ≠ 0)
    (hB : r ^ 6 + αlj * (1 - Λs lam) ^ 2 * B.c126 ≠ 0)
    (hXA : 0 < r ^ 2 + αc * Λcs lam ^ 2) (hYA : 0 < R ^ 2 + αc * Λcs lam ^ 2)
    (hXB : 0 < r ^ 2 + αc * (1 - Λcs lam) ^ 2) (hYB : 0 < R ^ 2 + αc * (1 - Λcs lam) ^ 2) :
    HasDerivAt
      (fun l => pertLjEnergy n αlj A B r (Λ l) (Λs l) +
        pertCrfEnergy n αc cs fpe qA qB crf R r (Λc l) (Λcs l))
      (pertLjDvdl n αlj A B r (Λ lam) (Λs lam) dΛ dΛs +
        pertCrfDvdl n αc cs fpe qA qB crf R r (Λc lam) (Λcs lam) dΛc dΛcs) lam :=
  (softcore_dvdl_is_deriv_lj n αlj A B r hΛ hΛs hA hB).add
    (softcore_dvdl_is_deriv_crf n αc cs fpe qA qB crf R r hΛc hΛcs hXA hYA hXB hYB)

/-! ## Force = −dV/dr at fixed λ -/

/-- **Force = −dV/dr** for the soft-core pair: the scalar force is `−∂V/∂r / r`. -/
theorem softcore_force_is_minus_grad (n : ℕ) (αlj αc cs fpe : ℝ) (A B : LjState)
    (qA qB crf R L Ls Lc Lcs : ℝ) {r : ℝ}
    (hA : r ^ 6 + αlj * Ls ^ 2 * A.c126 ≠ 0) (hB : r ^ 6 + αlj * (1 - Ls) ^ 2 * B.c126 ≠ 0)
    (hXA : 0 < r ^ 2 + αc * Lcs ^ 2) (hXB : 0 < r ^ 2 + αc * (1 - Lcs) ^ 2) :
    HasDerivAt
      (fun r => pertLjEnergy n αlj A B r L Ls + pertCrfEnergy n αc cs fpe qA qB crf R r Lc Lcs)
      (-(pertLjForce n αlj A B r L Ls + pertCrfForce n αc cs fpe qA qB crf R r Lc Lcs) * r) r := by
  have hLA := (scLj_force_is_minus_grad αlj A (Ls ^ 2) hA).const_mul ((1 - L) ^ n)
  have hLB := (scLj_force_is_minus_grad αlj B ((1 - Ls) ^ 2) hB).const_mul (L ^ n)
  have hCA := (scCrf_force_is_minus_grad αc cs qA crf R (Lcs ^ 2) hXA).const_mul ((1 - Lc) ^ n)
  have hCB := (scCrf_force_is_minus_grad αc cs qB crf R ((1 - Lcs) ^ 2) hXB).const_mul (Lc ^ n)
  have h := (hLA.add hLB).add ((hCA.add hCB).mul_const fpe)
  refine h.congr_deriv ?_
  simp only [pertLjForce, pertCrfForce]
  ring

/-! ## End states -/

theorem scLjEnergy_zero (αlj : ℝ) (p : LjState) (r : ℝ) :
    scLjEnergy αlj p 0 r = ljEnergy p.c6 p.c12 r := by
  simp only [scLjEnergy, ljEnergy, mul_zero, zero_mul, add_zero]
  ring

theorem scCrfEnergy_zero (αc q crf : ℝ) {R r : ℝ} (hR : 0 < R) (hr : 0 < r) :
    scCrfEnergy αc 1 q crf R 0 r = crfEnergy q crf R r := by
  have h1 : Real.sqrt (r ^ 2) = r := Real.sqrt_sq hr.le
  have h3 : Real.sqrt ((R ^ 2) ^ 3) = R ^ 3 := by
    rw [show (R ^ 2) ^ 3 = (R ^ 3) ^ 2 by ring]
    exact Real.sqrt_sq (pow_nonneg hR.le 3)
  simp only [scCrfEnergy, crfEnergy, mul_zero, add_zero, h1, h3]
  have := hR.ne'
  field_simp

theorem scRfExclEnergy_zero (αc q crf : ℝ) {R : ℝ} (hR : 0 < R) (r : ℝ) :
    scRfExclEnergy αc q crf R 0 r = rfExcludedEnergy q crf R r := by
  have h3 : Real.sqrt ((R ^ 2) ^ 3) = R ^ 3 := by
    rw [show (R ^ 2) ^ 3 = (R ^ 3) ^ 2 by ring]
    exact Real.sqrt_sq (pow_nonneg hR.le 3)
  simp only [scRfExclEnergy, rfExcludedEnergy, mul_zero, add_zero, h3]
  have := hR.ne'
  field_simp

/-- **λ = 0 is state A**: with `n ≥ 1`, the soft-core pair at `Λ = Λₛ = 0` is the plain
Lennard-Jones and reaction-field pair of state A (`ljEnergy`, `crfEnergy` with charge product
`qA · fpe`), and the excluded pair is `rfExcludedEnergy`. -/
theorem softcore_lambda_zero {n : ℕ} (hn : n ≠ 0) (αlj αc fpe : ℝ) (A B : LjState)
    (qA qB crf : ℝ) {R r : ℝ} (hR : 0 < R) (hr : 0 < r) :
    pertLjEnergy n αlj A B r 0 0 = ljEnergy A.c6 A.c12 r ∧
      pertCrfEnergy n αc 1 fpe qA qB crf R r 0 0 = crfEnergy (qA * fpe) crf R r ∧
      pertRfExclEnergy n αc fpe qA qB crf R r 0 0 = rfExcludedEnergy (qA * fpe) crf R r := by
  have h0 : (0 : ℝ) ^ n = 0 := zero_pow hn
  have hs : (0 : ℝ) ^ 2 = 0 := by norm_num
  refine ⟨?_, ?_, ?_⟩
  · simp only [pertLjEnergy, sub_zero, one_pow, one_mul, h0, zero_mul, add_zero, hs,
      scLjEnergy_zero]
  · simp only [pertCrfEnergy, sub_zero, one_pow, one_mul, h0, zero_mul, add_zero, hs,
      scCrfEnergy_zero αc qA crf hR hr, crfEnergy]
    ring
  · simp only [pertRfExclEnergy, sub_zero, one_pow, one_mul, h0, zero_mul, add_zero, hs,
      scRfExclEnergy_zero αc qA crf hR r, rfExcludedEnergy]
    ring

/-- **λ = 1 is state B**, likewise. -/
theorem softcore_lambda_one {n : ℕ} (hn : n ≠ 0) (αlj αc fpe : ℝ) (A B : LjState)
    (qA qB crf : ℝ) {R r : ℝ} (hR : 0 < R) (hr : 0 < r) :
    pertLjEnergy n αlj A B r 1 1 = ljEnergy B.c6 B.c12 r ∧
      pertCrfEnergy n αc 1 fpe qA qB crf R r 1 1 = crfEnergy (qB * fpe) crf R r ∧
      pertRfExclEnergy n αc fpe qA qB crf R r 1 1 = rfExcludedEnergy (qB * fpe) crf R r := by
  have h0 : (0 : ℝ) ^ n = 0 := zero_pow hn
  have hs : (0 : ℝ) ^ 2 = 0 := by norm_num
  refine ⟨?_, ?_, ?_⟩
  · simp only [pertLjEnergy, sub_self, one_pow, one_mul, h0, zero_mul, zero_add, hs,
      scLjEnergy_zero]
  · simp only [pertCrfEnergy, sub_self, one_pow, one_mul, h0, zero_mul, zero_add, hs,
      scCrfEnergy_zero αc qB crf hR hr, crfEnergy]
    ring
  · simp only [pertRfExclEnergy, sub_self, one_pow, one_mul, h0, zero_mul, zero_add, hs,
      scRfExclEnergy_zero αc qB crf hR r, rfExcludedEnergy]
    ring

/-- The self-term change vanishes at λ = 0 and is `(qa2 − qb2) c` at λ = 1. -/
theorem rfSelfChange_ends {n : ℕ} (hn : n ≠ 0) (qa2 qb2 c : ℝ) :
    rfSelfChange n qa2 qb2 c 0 = 0 ∧ rfSelfChange n qa2 qb2 c 1 = (qa2 - qb2) * c := by
  have h0 : (0 : ℝ) ^ n = 0 := zero_pow hn
  constructor <;> simp [rfSelfChange, h0]

end NewenTheory.Interactions.SoftCore

end
