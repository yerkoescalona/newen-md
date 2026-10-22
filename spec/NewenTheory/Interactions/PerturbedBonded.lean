module

public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Inv
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.Calculus.Deriv.Pow
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.IntervalCases
import Mathlib.Tactic.Ring

/-!
# Perturbed bonded interactions

In a free-energy calculation a bonded term has parameters for two end states A and B, and GROMOS
mixes them linearly in λ (GROMOS Vol. 2): the force constant `K(λ) = (1 − λ) K_A + λ K_B` and the
reference value `b0(λ) = (1 − λ) b0_A + λ b0_B` (`mix`). The term is then evaluated with the
mixed parameters: the quartic bond `¼ K (b² − b0²)²`, the cosine-harmonic angle `½ K (cos θ −
cos θ0)²` and the harmonic improper dihedral `½ K (ζ − ζ0)²` (`harmonicEnergy`, once in `cos θ`
and once in `ζ`). A proper dihedral mixes the energies instead,
`(1 − λ) K_A (1 + cos δ_A cos(m_A φ)) + λ K_B (…)`, since the two states may have different
multiplicities, and `cos(m φ)` is a polynomial in `cos φ` (`chebT`, `chebT_cos`).

Bonds, angles and impropers that exist in only one state use a soft-core form (PERTBONDSOFT and
friends; `softEnergy`): `½ [(1 − λ) K_A / (1 + α λ Δ²) + λ K_B / (1 + α (1 − λ) Δ²)] Δ²`.

The kernels compute `dV/dλ` in closed form, with GROMOS's individual λ: the term takes `L(λ)`
and multiplies by `L′(λ)`. The `*_dvdl_is_deriv` theorems prove that the closed form is the
derivative of the energy (`HasDerivAt`), and the `*_ends` theorems that λ = 0 and λ = 1 give the
plain terms of the two states.
-/

@[expose] public section

namespace NewenTheory.Interactions.PerturbedBonded

/-- Linear mixing of an A-state and a B-state parameter. -/
def mix (a b L : ℝ) : ℝ := (1 - L) * a + L * b

theorem mix_zero (a b : ℝ) : mix a b 0 = a := by simp [mix]

theorem mix_one (a b : ℝ) : mix a b 1 = b := by simp [mix]

/-- The mixed parameter is affine in `L`, so its derivative is `(b − a) · L′`. -/
theorem hasDerivAt_mix (a b : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ} (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => mix a b (Λ l)) ((b - a) * dΛ) lam := by
  have h := ((hΛ.const_sub 1).mul_const a).add (hΛ.mul_const b)
  refine h.congr_deriv ?_
  ring

/-! ## Quartic bond -/

/-- The quartic bond energy `¼ K (b² − b0²)²` at squared length `d2 = b²`. -/
noncomputable def quarticEnergy (K b0 d2 : ℝ) : ℝ := K / 4 * (d2 - b0 ^ 2) ^ 2

/-- The perturbed quartic bond. -/
noncomputable def pertQuarticEnergy (KA KB bA bB d2 L : ℝ) : ℝ :=
  quarticEnergy (mix KA KB L) (mix bA bB L) d2

/-- Its `dV/dλ`: `¼ (K_d Δ² − 4 K b_d b0 Δ) L′` with `Δ = b² − b0²`. -/
noncomputable def pertQuarticDvdl (KA KB bA bB d2 L dL : ℝ) : ℝ :=
  1 / 4 * ((KB - KA) * (d2 - mix bA bB L ^ 2) ^ 2 -
    4 * mix KA KB L * (bB - bA) * mix bA bB L * (d2 - mix bA bB L ^ 2)) * dL

theorem pert_quartic_dvdl_is_deriv (KA KB bA bB d2 : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => pertQuarticEnergy KA KB bA bB d2 (Λ l))
      (pertQuarticDvdl KA KB bA bB d2 (Λ lam) dΛ) lam := by
  have hK := hasDerivAt_mix KA KB hΛ
  have hb := hasDerivAt_mix bA bB hΛ
  have h := (hK.div_const 4).mul (((hb.pow 2).const_sub d2).pow 2)
  refine h.congr_deriv ?_
  simp only [pertQuarticDvdl, Pi.pow_apply]
  ring

/-- The kernel's force factor `f = −K Δ` (the force on atom `i` is `f · r_vec`) is
`−dE/db / b`. -/
theorem quartic_force_is_minus_grad (K b0 b : ℝ) :
    HasDerivAt (fun b => quarticEnergy K b0 (b ^ 2)) (-(-K * (b ^ 2 - b0 ^ 2)) * b) b := by
  have h := ((((hasDerivAt_pow 2 b).sub_const (b0 ^ 2)).pow 2).const_mul (K / 4))
  refine h.congr_deriv ?_
  ring

theorem pert_quartic_ends (KA KB bA bB d2 : ℝ) :
    pertQuarticEnergy KA KB bA bB d2 0 = quarticEnergy KA bA d2 ∧
      pertQuarticEnergy KA KB bA bB d2 1 = quarticEnergy KB bB d2 := by
  simp [pertQuarticEnergy, mix_zero, mix_one]

/-! ## Harmonic in a coordinate: the cosine-harmonic angle and the improper dihedral -/

/-- `½ K (x − x0)²`: the angle term with `x = cos θ`, the improper term with `x = ζ`. -/
noncomputable def harmonicEnergy (K x0 x : ℝ) : ℝ := K / 2 * (x - x0) ^ 2

noncomputable def pertHarmonicEnergy (KA KB xA xB x L : ℝ) : ℝ :=
  harmonicEnergy (mix KA KB L) (mix xA xB L) x

/-- Its `dV/dλ`: `½ (K_d Δ² − 2 K x_d Δ) L′` with `Δ = x − x0(λ)`. -/
noncomputable def pertHarmonicDvdl (KA KB xA xB x L dL : ℝ) : ℝ :=
  1 / 2 * ((KB - KA) * (x - mix xA xB L) ^ 2 -
    2 * mix KA KB L * (xB - xA) * (x - mix xA xB L)) * dL

theorem pert_harmonic_dvdl_is_deriv (KA KB xA xB x : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => pertHarmonicEnergy KA KB xA xB x (Λ l))
      (pertHarmonicDvdl KA KB xA xB x (Λ lam) dΛ) lam := by
  have hK := hasDerivAt_mix KA KB hΛ
  have hx := hasDerivAt_mix xA xB hΛ
  have h := (hK.div_const 2).mul ((hx.const_sub x).pow 2)
  refine h.congr_deriv ?_
  simp only [pertHarmonicDvdl, Pi.pow_apply]
  ring

/-- The kernel's `df = −K (x − x0)` is `−dE/dx`. -/
theorem harmonic_df_is_minus_grad (K x0 x : ℝ) :
    HasDerivAt (fun x => harmonicEnergy K x0 x) (-(-K * (x - x0))) x := by
  have h := (((hasDerivAt_id x).sub_const x0).pow 2).const_mul (K / 2)
  refine h.congr_deriv ?_
  simp
  ring

theorem pert_harmonic_ends (KA KB xA xB x : ℝ) :
    pertHarmonicEnergy KA KB xA xB x 0 = harmonicEnergy KA xA x ∧
      pertHarmonicEnergy KA KB xA xB x 1 = harmonicEnergy KB xB x := by
  simp [pertHarmonicEnergy, mix_zero, mix_one]

/-! ## Proper dihedral -/

/-- `cos(m φ)` as the polynomial `T_m(cos φ)` the kernels evaluate, for `m = 1 … 6`; GROMOS
returns `0` for `m = 0` and any other `m` (where `cos(0 · φ)` is `1`). -/
def chebT : ℕ → ℝ → ℝ
  | 1, c => c
  | 2, c => 2 * c ^ 2 - 1
  | 3, c => 4 * c ^ 3 - 3 * c
  | 4, c => 8 * c ^ 4 - 8 * c ^ 2 + 1
  | 5, c => 16 * c ^ 5 - 20 * c ^ 3 + 5 * c
  | 6, c => 32 * c ^ 6 - 48 * c ^ 4 + 18 * c ^ 2 - 1
  | _, _ => 0

/-- `d T_m / d(cos φ)`, the second value the kernels return. -/
def chebDT : ℕ → ℝ → ℝ
  | 1, _ => 1
  | 2, c => 4 * c
  | 3, c => 12 * c ^ 2 - 3
  | 4, c => 32 * c ^ 3 - 16 * c
  | 5, c => 80 * c ^ 4 - 60 * c ^ 2 + 5
  | 6, c => 192 * c ^ 5 - 192 * c ^ 3 + 36 * c
  | _, _ => 0

/-- `cos(x + θ) = 2 cos θ cos x − cos(x − θ)`, the recurrence behind `chebT`. -/
theorem cos_add_rec (x θ : ℝ) : Real.cos (x + θ) = 2 * Real.cos θ * Real.cos x - Real.cos (x - θ) := by
  rw [Real.cos_add, Real.cos_sub]
  ring

/-- **`chebT m (cos φ) = cos(m φ)`** for `m = 1 … 6`. -/
theorem chebT_cos (θ : ℝ) {m : ℕ} (h1 : 1 ≤ m) (h6 : m ≤ 6) :
    chebT m (Real.cos θ) = Real.cos (m * θ) := by
  set c := Real.cos θ with hc
  have c2 : Real.cos (2 * θ) = 2 * c ^ 2 - 1 := by
    rw [show 2 * θ = θ + θ by ring, cos_add_rec, sub_self, Real.cos_zero]; ring
  have c3 : Real.cos (3 * θ) = 4 * c ^ 3 - 3 * c := by
    rw [show 3 * θ = 2 * θ + θ by ring, cos_add_rec, show 2 * θ - θ = θ by ring, c2]; ring
  have c4 : Real.cos (4 * θ) = 8 * c ^ 4 - 8 * c ^ 2 + 1 := by
    rw [show 4 * θ = 3 * θ + θ by ring, cos_add_rec, show 3 * θ - θ = 2 * θ by ring, c3, c2]
    ring
  have c5 : Real.cos (5 * θ) = 16 * c ^ 5 - 20 * c ^ 3 + 5 * c := by
    rw [show 5 * θ = 4 * θ + θ by ring, cos_add_rec, show 4 * θ - θ = 3 * θ by ring, c4, c3]
    ring
  have c6 : Real.cos (6 * θ) = 32 * c ^ 6 - 48 * c ^ 4 + 18 * c ^ 2 - 1 := by
    rw [show 6 * θ = 5 * θ + θ by ring, cos_add_rec, show 5 * θ - θ = 4 * θ by ring, c5, c4]
    ring
  interval_cases m <;> simp [chebT, c2, c3, c4, c5, c6, hc]

theorem chebT_hasDerivAt (m : ℕ) (c : ℝ) : HasDerivAt (chebT m) (chebDT m c) c := by
  match m with
  | 0 => exact hasDerivAt_const c 0
  | 1 => exact hasDerivAt_id c
  | 2 =>
    have h := ((hasDerivAt_pow 2 c).const_mul 2).sub_const 1
    exact h.congr_deriv (by simp [chebDT]; ring)
  | 3 =>
    have h := ((hasDerivAt_pow 3 c).const_mul 4).sub ((hasDerivAt_id c).const_mul 3)
    exact h.congr_deriv (by simp [chebDT]; ring)
  | 4 =>
    have h := (((hasDerivAt_pow 4 c).const_mul 8).sub ((hasDerivAt_pow 2 c).const_mul 8)).add_const 1
    exact h.congr_deriv (by simp [chebDT]; ring)
  | 5 =>
    have h := (((hasDerivAt_pow 5 c).const_mul 16).sub ((hasDerivAt_pow 3 c).const_mul 20)).add
      ((hasDerivAt_id c).const_mul 5)
    exact h.congr_deriv (by simp [chebDT]; ring)
  | 6 =>
    have h := ((((hasDerivAt_pow 6 c).const_mul 32).sub ((hasDerivAt_pow 4 c).const_mul 48)).add
      ((hasDerivAt_pow 2 c).const_mul 18)).sub_const 1
    exact h.congr_deriv (by simp [chebDT]; ring)
  | k + 7 => exact hasDerivAt_const c 0

/-- One dihedral state: `K (1 + cos δ · T_m(cos φ))`. -/
noncomputable def dihedralEnergy (K cospd : ℝ) (m : ℕ) (c : ℝ) : ℝ := K * (1 + cospd * chebT m c)

/-- The perturbed dihedral mixes the two states' energies; its `dV/dλ` is `(E_B − E_A) L′`. -/
theorem pert_dihedral_dvdl_is_deriv (EA EB : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => mix EA EB (Λ l)) ((EB - EA) * dΛ) lam :=
  hasDerivAt_mix EA EB hΛ

/-! ## Soft-core harmonic terms -/

/-- The soft-core term in the coordinate `x` with softness `α`. -/
noncomputable def softEnergy (KA KB xA xB α x L : ℝ) : ℝ :=
  1 / 2 * ((1 - L) * KA / (1 + α * L * (x - mix xA xB L) ^ 2) +
    L * KB / (1 + α * (1 - L) * (x - mix xA xB L) ^ 2)) * (x - mix xA xB L) ^ 2

/-- Its `dV/dλ`, gromosXX's `e_lambda` formula. -/
noncomputable def softDvdl (KA KB xA xB α x L dL : ℝ) : ℝ :=
  let d := x - mix xA xB L
  let sA := 1 + α * L * d ^ 2
  let sB := 1 + α * (1 - L) * d ^ 2
  let st1 := 1 + α * d ^ 2
  let st2 := -2 * α * L * (1 - L) * d * (xB - xA)
  (1 / 2 * d ^ 2 * (KA / sA ^ 2 * (-st1 - st2) + KB / sB ^ 2 * (st1 - st2)) -
    ((1 - L) * KA / sA + L * KB / sB) * d * (xB - xA)) * dL

/-- `∂E/∂x = K_eff · Δ` with `K_eff = (1 − λ) K_A / S_A² + λ K_B / S_B²`. -/
noncomputable def softGrad (KA KB xA xB α x L : ℝ) : ℝ :=
  let d := x - mix xA xB L
  ((1 - L) * KA / (1 + α * L * d ^ 2) ^ 2 + L * KB / (1 + α * (1 - L) * d ^ 2) ^ 2) * d

theorem soft_dvdl_is_deriv (KA KB xA xB α x : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam)
    (hA : 1 + α * Λ lam * (x - mix xA xB (Λ lam)) ^ 2 ≠ 0)
    (hB : 1 + α * (1 - Λ lam) * (x - mix xA xB (Λ lam)) ^ 2 ≠ 0) :
    HasDerivAt (fun l => softEnergy KA KB xA xB α x (Λ l))
      (softDvdl KA KB xA xB α x (Λ lam) dΛ) lam := by
  have hd := (hasDerivAt_mix xA xB hΛ).const_sub x
  have hd2 := hd.pow 2
  have hsA := ((hΛ.const_mul α).mul hd2).const_add 1
  have hsB := (((hΛ.const_sub 1).const_mul α).mul hd2).const_add 1
  have hwA := (hΛ.const_sub 1).mul_const KA
  have hwB := hΛ.mul_const KB
  have h := (((hwA.div hsA hA).add (hwB.div hsB hB)).const_mul (1 / 2)).mul hd2
  refine h.congr_deriv ?_
  dsimp only [softDvdl]
  simp only [Pi.pow_apply, Pi.add_apply, Pi.div_apply, Pi.mul_apply]
  field_simp
  ring

/-- One state of the soft term: `d/dΔ [½ w K Δ² / (1 + k Δ²)] = w K Δ / (1 + k Δ²)²`. -/
theorem soft_state_grad (w K k d : ℝ) (hs : 1 + k * d ^ 2 ≠ 0) :
    1 / 2 * ((0 * (1 + k * d ^ 2) - w * K * (k * (2 * d))) / (1 + k * d ^ 2) ^ 2) * d ^ 2 +
      1 / 2 * (w * K / (1 + k * d ^ 2)) * (2 * d) = w * K / (1 + k * d ^ 2) ^ 2 * d := by
  field_simp
  ring

theorem soft_grad_is_deriv (KA KB xA xB α L : ℝ) {x : ℝ}
    (hA : 1 + α * L * (x - mix xA xB L) ^ 2 ≠ 0)
    (hB : 1 + α * (1 - L) * (x - mix xA xB L) ^ 2 ≠ 0) :
    HasDerivAt (fun x => softEnergy KA KB xA xB α x L) (softGrad KA KB xA xB α x L) x := by
  have hd := (hasDerivAt_id' x).sub_const (mix xA xB L)
  have hd2 := hd.pow 2
  have hsA := (hd2.const_mul (α * L)).const_add 1
  have hsB := (hd2.const_mul (α * (1 - L))).const_add 1
  have h := ((((hasDerivAt_const x ((1 - L) * KA)).div hsA hA).add
    ((hasDerivAt_const x (L * KB)).div hsB hB)).const_mul (1 / 2)).mul hd2
  refine h.congr_deriv ?_
  simp only [softGrad, Pi.pow_apply, Pi.add_apply, Pi.div_apply]
  generalize x - mix xA xB L = d at hA hB ⊢
  have eA := soft_state_grad (1 - L) KA (α * L) d hA
  have eB := soft_state_grad L KB (α * (1 - L)) d hB
  norm_num
  linear_combination eA + eB

theorem soft_ends (KA KB xA xB α x : ℝ) :
    softEnergy KA KB xA xB α x 0 = harmonicEnergy KA xA x ∧
      softEnergy KA KB xA xB α x 1 = harmonicEnergy KB xB x := by
  constructor <;> simp only [softEnergy, harmonicEnergy, mix_zero, mix_one] <;> ring

end NewenTheory.Interactions.PerturbedBonded

end
