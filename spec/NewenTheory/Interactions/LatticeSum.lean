module

public import Mathlib.Analysis.Calculus.Deriv.Polynomial
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
public import Mathlib.Analysis.SpecialFunctions.ExpDeriv
public import Mathlib.Analysis.SpecialFunctions.Sqrt
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
public import Mathlib.Data.Matrix.Mul
public import Mathlib.LinearAlgebra.CrossProduct
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring

/-!
# Lattice-sum electrostatics (GROMOS Vol. 2)

In a periodic system the Coulomb sum `Σ q_i q_j / r` over all images converges only
conditionally. GROMOS (Vol. 2, its reference HU99.4) splits it with a charge-shaping function
`γ` of width `a`: every point charge is screened by a shaped charge of the opposite sign, and
the shaped charges are added back in reciprocal space. The electrostatic energy is then
(Vol. 2, `chg_fe_all`)

  `E = E_γ + E_η + E_A + E_slf + E_srf`:

* `E_η`, the **real-space** sum: pairs at distance `r` interact through `η(r/a)/r`, where the
  switch function `η` (`eta_function_def`) falls from `η(0) = 1` to zero
  (`realSpacePair`; excluded pairs keep only `(η − 1)/r`, `realSpaceExcluded`);
* `E_γ`, the **reciprocal-space** sum over the wave vectors `k ≠ 0` with weight `γ̂(ak)/k²`
  (`e_gamma_ewald`, `kSpace`), where `γ̂` is the Fourier transform of the shape
  (`latsum_gamma_ft`);
* `E_A` and `E_slf`, constants built from the total charge `S = Σ q_i`, the sum of squares
  `S̃² = Σ q_i²` and the three constants `A1, A2, A3` (`e_a_def`, `e_slf`; `aTerm`, `selfTerm`);
* `E_srf`, the **surface** term of the box dipole `M` in a dielectric of permittivity `ε_LS`
  (`e_srf`, `surface`).

That the five terms add up to the lattice sum of `1/r` is the definition of the method (Vol. 2,
`chg_fe_all`, and the derivation before it); it is not proved here. What is proved: the switch
functions and their derivatives of the twelve shapes of Vol. 2 (`eta_hasDerivAt`), their
compact support (`eta_compact_support`), `η(0) = 1` (`eta_zero_eq_one`), that the
truncated-polynomial shapes carry unit charge (`chargeShape_normalized`, whence `γ̂(0) = 1`),
that the tabulated `γ̂'` is the derivative of the tabulated `γ̂` (`gammaHat_hasDerivAt`), that
`A3` is `−η'(0)/a` (`a3Tilde_eq_neg_etaDeriv_zero`), and that the forces of the real-space,
reciprocal-space and surface terms are minus the gradients of their energies. For the Ewald sum:
the reciprocal basis is dual to the box (`latticeVector_reciprocal_dot`), so moving an atom by a
lattice vector changes no `C(k)² + S(k)²` (`structureFactorSq_lattice_shift`); the set of wave
vectors `kSet` is symmetric (`kSet_symmetric`) and `C² + S²` even in `k`
(`structureFactorSq_neg`), which is the ½ of the reciprocal-space energy (`kSpace_half_space`);
the k-space force is minus the gradient over `kSet` (`ewald_total_force_is_minus_grad`); and the
shells of the integer lattice partition the cube, each the union of gromosXX's half shell and its
mirror image, so `A2` summed by doubled half shells is the lattice sum (`sum_cube_erase_zero`,
`shell_eq_face_union_neg`, `a2Shells_eq_cube`).

Every energy here is in units where the charges carry no `1/(4πε₀)`; `ε0` is the permittivity
of vacuum (gromosXX `four_pi_eps_i = 1/(4πε₀)`, `eps0_i = 1/ε₀`), and the relative permittivity
inside the box is 1 (Vol. 2's `ε_ls`).

Mathlib has no error function, so `erf` is defined here by its integral; its derivative is
proved (`erfc_hasDerivAt`), its limit `erf x → 1` is not needed and not proved.
-/

@[expose] public section

namespace NewenTheory.Interactions.LatticeSum

open Polynomial Real

/-! ## The error function -/

/-- `erf x = 2/√π · ∫₀ˣ exp(−t²) dt`. -/
noncomputable def erf (x : ℝ) : ℝ := 2 / √π * ∫ t in (0 : ℝ)..x, exp (-t ^ 2)

/-- `erfc x = 1 − erf x`. -/
noncomputable def erfc (x : ℝ) : ℝ := 1 - erf x

theorem erfc_zero : erfc 0 = 1 := by simp [erfc, erf]

/-- `erfc' = −2/√π · exp(−x²)`, at every `x`. -/
theorem erfc_hasDerivAt (x : ℝ) : HasDerivAt erfc (-(2 / √π) * exp (-x ^ 2)) x := by
  have hi : HasDerivAt (fun u => ∫ t in (0 : ℝ)..u, exp (-t ^ 2)) (exp (-x ^ 2)) x :=
    ((show Continuous fun t : ℝ => exp (-t ^ 2) by fun_prop).integral_hasStrictDerivAt 0 x).hasDerivAt
  have h := (hi.const_mul (2 / √π)).const_sub 1
  exact h.congr_deriv (by ring)

/-! ## The charge shapes -/

/-- The order `N` of a truncated-polynomial shape, `NSHAPE = 0 … 10`. -/
inductive TPOrder
  | o0 | o1 | o2 | o3 | o4 | o5 | o6 | o7 | o8 | o9 | o10
  deriving DecidableEq

/-- The charge-shaping functions of GROMOS (Vol. 2, table `ls_charge_shaping`): the Gaussian
(`NSHAPE = −1`) and the truncated polynomials of order `N = 0 … 10`. -/
inductive Shape
  | gaussian
  | poly (N : TPOrder)
  deriving DecidableEq

/-- `N` as a natural number. -/
def TPOrder.toNat : TPOrder → ℕ
  | .o0 => 0 | .o1 => 1 | .o2 => 2 | .o3 => 3 | .o4 => 4 | .o5 => 5 | .o6 => 6 | .o7 => 7
  | .o8 => 8 | .o9 => 9 | .o10 => 10

/-- The LONGRANGE switch `NSHAPE` of a shape. -/
def Shape.nshape : Shape → ℤ
  | .gaussian => -1
  | .poly N => N.toNat

/-- The shape of `NSHAPE`, if there is one. -/
def Shape.ofNshape (n : ℤ) : Option Shape :=
  [Shape.gaussian, .poly .o0, .poly .o1, .poly .o2, .poly .o3, .poly .o4, .poly .o5, .poly .o6,
    .poly .o7, .poly .o8, .poly .o9, .poly .o10].find? (·.nshape = n)

/-- The Heaviside step, `H(x) = 1` for `x ≥ 0` (Vol. 2's convention). -/
noncomputable def heaviside (x : ℝ) : ℝ := if 0 ≤ x then 1 else 0

/-- The switch function `η(ξ)` of a truncated polynomial without its factor `H(1 − ξ)`, as a
polynomial (Vol. 2, table `ls_switch_function`). -/
noncomputable def etaPoly : TPOrder → ℝ[X]
  | .o0 => C (1 / 2) * (1 - X) ^ 2 * (X + 2)
  | .o1 => (1 - X) ^ 3 * (X + 1)
  | .o2 => C (1 / 2) * (1 - X) ^ 4 * (3 * X + 2)
  | .o3 => C (1 / 4) * (1 - X) ^ 4 * (4 * X ^ 2 + 7 * X + 4)
  | .o4 => C (1 / 8) * (1 - X) ^ 5 * (15 * X ^ 2 + 19 * X + 8)
  | .o5 => (1 - X) ^ 6 * (3 * X ^ 2 + 3 * X + 1)
  | .o6 => C (1 / 16) * (1 - X) ^ 6 * (35 * X ^ 3 + 66 * X ^ 2 + 51 * X + 16)
  | .o7 => C (1 / 8) * (1 - X) ^ 7 * (32 * X ^ 3 + 49 * X ^ 2 + 31 * X + 8)
  | .o8 => C (1 / 16) * (1 - X) ^ 8 * (105 * X ^ 3 + 136 * X ^ 2 + 73 * X + 16)
  | .o9 => C (1 / 32) * (1 - X) ^ 8 * (160 * X ^ 4 + 335 * X ^ 3 + 312 * X ^ 2 + 151 * X + 32)
  | .o10 =>
    C (1 / 128) * (1 - X) ^ 9 * (1155 * X ^ 4 + 2075 * X ^ 3 + 1665 * X ^ 2 + 697 * X + 128)

/-- `η'(ξ)` of a truncated polynomial without its factor `H(1 − ξ)`: minus the column `−η'(ξ)`
of Vol. 2, table `ls_switch_prime`. -/
noncomputable def etaPolyDeriv : TPOrder → ℝ → ℝ
  | .o0, ξ => -(3 / 2 * (1 - ξ) * (ξ + 1))
  | .o1, ξ => -(2 * (1 - ξ) ^ 2 * (2 * ξ + 1))
  | .o2, ξ => -(5 / 2 * (1 - ξ) ^ 3 * (3 * ξ + 1))
  | .o3, ξ => -(3 / 4 * (1 - ξ) ^ 3 * (8 * ξ ^ 2 + 9 * ξ + 3))
  | .o4, ξ => -(21 / 8 * (1 - ξ) ^ 4 * (5 * ξ ^ 2 + 4 * ξ + 1))
  | .o5, ξ => -(3 * (1 - ξ) ^ 5 * (8 * ξ ^ 2 + 5 * ξ + 1))
  | .o6, ξ => -(9 / 16 * (1 - ξ) ^ 5 * (35 * ξ ^ 3 + 47 * ξ ^ 2 + 25 * ξ + 5))
  | .o7, ξ => -(5 / 8 * (1 - ξ) ^ 6 * (64 * ξ ^ 3 + 69 * ξ ^ 2 + 30 * ξ + 5))
  | .o8, ξ => -(55 / 16 * (1 - ξ) ^ 7 * (21 * ξ ^ 3 + 19 * ξ ^ 2 + 7 * ξ + 1))
  | .o9, ξ => -(15 / 32 * (1 - ξ) ^ 7 * (128 * ξ ^ 4 + 203 * ξ ^ 3 + 141 * ξ ^ 2 + 49 * ξ + 7))
  | .o10, ξ => -(65 / 128 * (1 - ξ) ^ 8 * (231 * ξ ^ 4 + 312 * ξ ^ 3 + 186 * ξ ^ 2 + 56 * ξ + 7))

/-- `π γ(ξ)` of a truncated polynomial without its factor `H(1 − ξ)` (Vol. 2, table
`ls_charge_shaping`; the shaped charge density of width `a` is `a⁻³ γ(r/a)`). -/
noncomputable def piGammaPoly : TPOrder → ℝ → ℝ
  | .o0, _ => 3 / 4
  | .o1, ξ => 3 * (1 - ξ)
  | .o2, ξ => 15 / 2 * (1 - ξ) ^ 2
  | .o3, ξ => 15 / 4 * (1 - ξ) ^ 2 * (1 + 2 * ξ)
  | .o4, ξ => 105 / 16 * (1 - ξ) ^ 3 * (3 * ξ + 1)
  | .o5, ξ => 21 / 2 * (1 - ξ) ^ 4 * (4 * ξ + 1)
  | .o6, ξ => 63 / 8 * (1 - ξ) ^ 4 * (5 * ξ ^ 2 + 4 * ξ + 1)
  | .o7, ξ => 45 / 4 * (1 - ξ) ^ 5 * (8 * ξ ^ 2 + 5 * ξ + 1)
  | .o8, ξ => 165 / 32 * (1 - ξ) ^ 6 * (35 * ξ ^ 2 + 18 * ξ + 3)
  | .o9, ξ => 165 / 64 * (1 - ξ) ^ 6 * (64 * ξ ^ 3 + 69 * ξ ^ 2 + 30 * ξ + 5)
  | .o10, ξ => 2145 / 128 * (1 - ξ) ^ 7 * (21 * ξ ^ 3 + 19 * ξ ^ 2 + 7 * ξ + 1)

/-- The switch function `η(ξ)`, `ξ = r/a` (Vol. 2, `eta_function_def`, table
`ls_switch_function`): `erfc ξ` for the Gaussian, the polynomial times `H(1 − ξ)` otherwise. -/
noncomputable def eta : Shape → ℝ → ℝ
  | .gaussian, ξ => erfc ξ
  | .poly N, ξ => (etaPoly N).eval ξ * heaviside (1 - ξ)

/-- `η'(ξ)` (Vol. 2, table `ls_switch_prime`). -/
noncomputable def etaDeriv : Shape → ℝ → ℝ
  | .gaussian, ξ => -(2 / √π) * exp (-ξ ^ 2)
  | .poly N, ξ => etaPolyDeriv N ξ * heaviside (1 - ξ)

/-- `π γ(ξ)`, the charge-shaping function times `π` (Vol. 2, table `ls_charge_shaping`). -/
noncomputable def piGamma : Shape → ℝ → ℝ
  | .gaussian, ξ => exp (-ξ ^ 2) / √π
  | .poly N, ξ => piGammaPoly N ξ * heaviside (1 - ξ)

/-- `simp only` set that differentiates and evaluates the polynomials above. -/
macro "poly_simp" : tactic => `(tactic| simp only [derivative_mul, derivative_pow,
  derivative_sub, derivative_add, derivative_neg, derivative_X, derivative_C, derivative_one, derivative_zero,
  derivative_ofNat, eval_mul, eval_pow, eval_sub, eval_add, eval_neg, eval_X, eval_C, eval_one,
  eval_zero, eval_natCast, eval_ofNat, map_natCast, zero_mul, mul_zero, add_zero, zero_add,
  sub_zero, Nat.cast_ofNat])

/-- The tabulated `η'` is the derivative of the tabulated `η`. -/
theorem etaPoly_derivative (N : TPOrder) (ξ : ℝ) :
    (etaPoly N).derivative.eval ξ = etaPolyDeriv N ξ := by
  cases N <;> simp only [etaPoly, etaPolyDeriv] <;> poly_simp <;> ring

/-- `η'' = 4 ξ · π γ(ξ)`: the tabulated `η` and `γ` agree. Differentiating `eta_function_def`,
`η(ξ) = 4π ∫_ξ^∞ ρ (ρ − ξ) γ(ρ) dρ`, twice gives `η''(ξ) = 4π ξ γ(ξ)`. -/
theorem etaPoly_derivative_two (N : TPOrder) (ξ : ℝ) :
    (etaPoly N).derivative.derivative.eval ξ = 4 * ξ * piGammaPoly N ξ := by
  cases N <;> simp only [etaPoly, piGammaPoly] <;> poly_simp <;> ring

theorem etaPoly_eval_zero (N : TPOrder) : (etaPoly N).eval 0 = 1 := by
  cases N <;> simp only [etaPoly] <;> poly_simp <;> norm_num

theorem etaPoly_eval_one (N : TPOrder) : (etaPoly N).eval 1 = 0 := by
  cases N <;> simp only [etaPoly] <;> poly_simp <;> norm_num

theorem etaPolyDeriv_one (N : TPOrder) : etaPolyDeriv N 1 = 0 := by
  cases N <;> simp [etaPolyDeriv]

/-- A polynomial that vanishes with its derivative at 1, cut off by `H(1 − ξ)`, is
differentiable everywhere, with the cut-off derivative. -/
theorem hasDerivAt_mul_heaviside (p : ℝ[X]) (h0 : p.eval 1 = 0)
    (h1 : p.derivative.eval 1 = 0) (ξ : ℝ) :
    HasDerivAt (fun x => p.eval x * heaviside (1 - x))
      (p.derivative.eval ξ * heaviside (1 - ξ)) ξ := by
  rcases lt_trichotomy ξ 1 with hξ | rfl | hξ
  · have hev : (fun x => p.eval x * heaviside (1 - x)) =ᶠ[nhds ξ] fun x => p.eval x := by
      filter_upwards [Iio_mem_nhds hξ] with x hx
      simp [heaviside, sub_nonneg.mpr (le_of_lt hx)]
    rw [show heaviside (1 - ξ) = 1 by simp [heaviside, sub_nonneg.mpr hξ.le], mul_one]
    exact (Polynomial.hasDerivAt p ξ).congr_of_eventuallyEq hev
  · have hle : HasDerivWithinAt (fun x => p.eval x * heaviside (1 - x))
        (p.derivative.eval 1 * heaviside (1 - 1)) (Set.Iic 1) 1 := by
      rw [show heaviside (1 - 1) = 1 by simp [heaviside], mul_one]
      refine (Polynomial.hasDerivAt p 1).hasDerivWithinAt.congr (fun x hx => ?_) ?_
      · simp [heaviside, sub_nonneg.mpr (Set.mem_Iic.mp hx)]
      · simp [heaviside]
    have hge : HasDerivWithinAt (fun x => p.eval x * heaviside (1 - x))
        (p.derivative.eval 1 * heaviside (1 - 1)) (Set.Ici 1) 1 := by
      rw [h1, zero_mul]
      refine (hasDerivAt_const (1 : ℝ) (0 : ℝ)).hasDerivWithinAt.congr (fun x hx => ?_) ?_
      · rcases (Set.mem_Ici.mp hx).lt_or_eq with hx | hx
        · simp [heaviside, not_le.mpr (sub_neg.mpr hx)]
        · simp [← hx, h0]
      · simp [h0]
    have h := hle.union hge
    rwa [Set.Iic_union_Ici, hasDerivWithinAt_univ] at h
  · have hev : (fun x => p.eval x * heaviside (1 - x)) =ᶠ[nhds ξ] fun _ => 0 := by
      filter_upwards [Ioi_mem_nhds hξ] with x hx
      simp [heaviside, not_le.mpr (sub_neg.mpr (Set.mem_Ioi.mp hx))]
    rw [show heaviside (1 - ξ) = 0 by simp [heaviside, not_le.mpr (sub_neg.mpr hξ)], mul_zero]
    exact (hasDerivAt_const ξ (0 : ℝ)).congr_of_eventuallyEq hev

/-- **`η'` is the derivative of `η`**, for every shape and every `ξ`: the Gaussian through
`erfc' = −2/√π · exp(−ξ²)`, the truncated polynomials across the edge of their support too
(they vanish there with their first derivative). -/
theorem eta_hasDerivAt (s : Shape) (ξ : ℝ) : HasDerivAt (eta s) (etaDeriv s ξ) ξ := by
  cases s with
  | gaussian => exact erfc_hasDerivAt ξ
  | poly N =>
    have h := hasDerivAt_mul_heaviside (etaPoly N) (etaPoly_eval_one N)
      (by rw [etaPoly_derivative, etaPolyDeriv_one]) ξ
    rw [etaPoly_derivative] at h
    exact h

/-- The truncated polynomials and their derivatives vanish for `ξ ≥ 1` (`r ≥ a`): the
real-space sum is finite. -/
theorem eta_compact_support (N : TPOrder) {ξ : ℝ} (hξ : 1 ≤ ξ) :
    eta (.poly N) ξ = 0 ∧ etaDeriv (.poly N) ξ = 0 := by
  rcases hξ.lt_or_eq with hξ | rfl
  · simp [eta, etaDeriv, heaviside, not_le.mpr (sub_neg.mpr hξ)]
  · simp [eta, etaDeriv, etaPoly_eval_one, etaPolyDeriv_one]

/-- `η(0) = 1` for every shape: at contact the real-space term is the full Coulomb term. -/
theorem eta_zero_eq_one (s : Shape) : eta s 0 = 1 := by
  cases s with
  | gaussian => exact erfc_zero
  | poly N => simp [eta, heaviside, etaPoly_eval_zero]

/-- **Unit charge.** The truncated-polynomial shapes are normalised (Vol. 2,
`ewa_shapefct_norm`): `4π ∫₀^∞ ξ² γ(ξ) dξ = 4 ∫₀¹ ξ² · πγ(ξ) dξ = 1`. Proof: `ξ η'(ξ) − η(ξ)` is
an antiderivative of `ξ η'' = 4 ξ² · πγ`, and `η(0) = 1`, `η(1) = η'(1) = 0`. -/
theorem chargeShape_normalized (N : TPOrder) :
    ∫ ξ in (0 : ℝ)..1, 4 * ξ ^ 2 * piGamma (.poly N) ξ = 1 := by
  set p := etaPoly N
  set q : ℝ[X] := X * p.derivative - p
  have hq : ∀ x, q.derivative.eval x = 4 * x ^ 2 * piGammaPoly N x := by
    intro x
    simp only [q, derivative_mul, derivative_sub, derivative_X, one_mul, eval_sub, eval_add,
      eval_mul, eval_X]
    rw [etaPoly_derivative_two]
    ring
  have hcongr : Set.EqOn (fun ξ => 4 * ξ ^ 2 * piGamma (.poly N) ξ)
      (fun x => q.derivative.eval x) (Set.uIcc 0 1) := by
    intro x hx
    rw [Set.uIcc_of_le zero_le_one] at hx
    simp [piGamma, heaviside, sub_nonneg.mpr hx.2, hq]
  rw [intervalIntegral.integral_congr hcongr,
    intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => Polynomial.hasDerivAt q x)
      ((Polynomial.differentiable _).continuous.intervalIntegrable _ _)]
  simp only [q, eval_sub, eval_mul, eval_X, etaPoly_derivative, etaPolyDeriv_one,
    etaPoly_eval_one, p, etaPoly_eval_zero]
  ring

/-! ## Fourier coefficients -/

/-- `(A(κ) + B(κ) cos κ + D(κ) sin κ) / κⁿ`, the form of every tabulated `γ̂` and `γ̂'` of a
truncated polynomial. -/
noncomputable def fourierForm (A B D : ℝ[X]) (n : ℕ) (κ : ℝ) : ℝ :=
  (A.eval κ + B.eval κ * cos κ + D.eval κ * sin κ) / κ ^ n

/-- The three polynomials of `κ^(N+3) γ̂(κ) = A + B cos κ + D sin κ` (Vol. 2, table `ls_fcoef`). -/
noncomputable def fcoef : TPOrder → ℝ[X] × ℝ[X] × ℝ[X]
  | .o0 => (0, 3 * -X, 3)
  | .o1 => (12 * 2, 12 * -2, 12 * -X)
  | .o2 => (60 * (2 * X), 60 * X, 60 * -3)
  | .o3 => (90 * 8, 90 * (X ^ 2 - 8), 90 * (-5 * X))
  | .o4 => (630 * (8 * X), 630 * (7 * X), 630 * (X ^ 2 - 15))
  | .o5 => (5040 * (4 * (X ^ 2 - 6)), 5040 * -(X ^ 2 - 24), 5040 * (9 * X))
  | .o6 => (7560 * (48 * X), 7560 * (-(X ^ 2 - 57) * X), 7560 * (3 * (4 * X ^ 2 - 35)))
  | .o7 => (75600 * (24 * (X ^ 2 - 8)), 75600 * (-3 * (5 * X ^ 2 - 64)),
      75600 * (-(X ^ 2 - 87) * X))
  | .o8 => (831600 * (8 * X * (X ^ 2 - 24)), 831600 * ((X ^ 2 - 123) * X),
      831600 * -(18 * X ^ 2 - 315))
  | .o9 => (1247400 * (192 * (X ^ 2 - 10)), 1247400 * (X ^ 4 - 207 * X ^ 2 + 1920),
      1247400 * (-(22 * X ^ 2 - 975) * X))
  | .o10 => (16216200 * (64 * X * (X ^ 2 - 30)), 16216200 * ((26 * X ^ 2 - 1545) * X),
      16216200 * (X ^ 4 - 285 * X ^ 2 + 3465))

/-- The three polynomials of `κ^(N+4) γ̂'(κ)` (Vol. 2, table `ls_fcoef_prime`). -/
noncomputable def fcoefPrime : TPOrder → ℝ[X] × ℝ[X] × ℝ[X]
  | .o0 => (0, 3 * (3 * X), 3 * (X ^ 2 - 3))
  | .o1 => (12 * -8, 12 * -(X ^ 2 - 8), 12 * (5 * X))
  | .o2 => (60 * (-8 * X), 60 * (-7 * X), 60 * -(X ^ 2 - 15))
  | .o3 => (90 * -48, 90 * -(9 * X ^ 2 - 48), 90 * (-(X ^ 2 - 33) * X))
  | .o4 => (630 * (-48 * X), 630 * ((X ^ 2 - 57) * X), 630 * (-3 * (4 * X ^ 2 - 35)))
  | .o5 => (5040 * (-24 * (X ^ 2 - 8)), 5040 * (3 * (5 * X ^ 2 - 64)),
      5040 * ((X ^ 2 - 87) * X))
  | .o6 => (7560 * (-384 * X), 7560 * (3 * (6 * X ^ 2 - 187) * X),
      7560 * (X ^ 4 - 141 * X ^ 2 + 945))
  | .o7 => (75600 * (-192 * (X ^ 2 - 10)), 75600 * -(X ^ 4 - 207 * X ^ 2 + 1920),
      75600 * ((22 * X ^ 2 - 975) * X))
  | .o8 => (831600 * (-64 * X * (X ^ 2 - 30)), 831600 * (-(26 * X ^ 2 - 1545) * X),
      831600 * -(X ^ 4 - 285 * X ^ 2 + 3465))
  | .o9 => (1247400 * (-1920 * (X ^ 2 - 12)), 1247400 * (-15 * (2 * X ^ 4 - 203 * X ^ 2 + 1536)),
      1247400 * (-(X ^ 4 - 405 * X ^ 2 + 12645) * X))
  | .o10 => (16216200 * (-640 * X * (X ^ 2 - 36)),
      16216200 * ((X ^ 4 - 545 * X ^ 2 + 22005) * X),
      16216200 * (-5 * (7 * X ^ 4 - 936 * X ^ 2 + 9009)))

/-- The Fourier coefficient `γ̂(κ)` of the shape at `κ = ak` (Vol. 2, `latsum_gamma_ft`, table
`ls_fcoef`): `exp(−κ²/4)` for the Gaussian; for a truncated polynomial the tabulated closed form,
and `1` at `κ = 0`, as Vol. 2 defines it there (the unit charge, `chargeShape_normalized`). -/
noncomputable def gammaHat : Shape → ℝ → ℝ
  | .gaussian, κ => exp (-κ ^ 2 / 4)
  | .poly N, κ =>
    if κ = 0 then 1 else fourierForm (fcoef N).1 (fcoef N).2.1 (fcoef N).2.2 (N.toNat + 3) κ

/-- `γ̂'(κ)` (Vol. 2, table `ls_fcoef_prime`), `0` at `κ = 0` for a truncated polynomial. -/
noncomputable def gammaHatDeriv : Shape → ℝ → ℝ
  | .gaussian, κ => -(1 / 2) * κ * exp (-κ ^ 2 / 4)
  | .poly N, κ =>
    if κ = 0 then 0
    else fourierForm (fcoefPrime N).1 (fcoefPrime N).2.1 (fcoefPrime N).2.2 (N.toNat + 4) κ

/-- `γ̂(0) = 1` for every shape: the Gaussian's `exp 0`; for the truncated polynomials Vol. 2's
value at `k = 0`, the total charge of the normalised shape (`chargeShape_normalized`). -/
theorem gammaHat_zero_eq_one (s : Shape) : gammaHat s 0 = 1 := by
  cases s <;> simp [gammaHat]

/-- The derivative of a `fourierForm` is a `fourierForm` with one power of `κ` more. -/
theorem fourierForm_hasDerivAt (A B D : ℝ[X]) (n : ℕ) {κ : ℝ} (hκ : κ ≠ 0) :
    HasDerivAt (fourierForm A B D (n + 1))
      (fourierForm (X * derivative A - C ((n : ℝ) + 1) * A)
        (X * derivative B + X * D - C ((n : ℝ) + 1) * B)
        (X * derivative D - X * B - C ((n : ℝ) + 1) * D) (n + 2) κ) κ := by
  have hnum : HasDerivAt (fun x => A.eval x + B.eval x * cos x + D.eval x * sin x)
      ((derivative A).eval κ + ((derivative B).eval κ * cos κ + B.eval κ * -sin κ) +
        ((derivative D).eval κ * sin κ + D.eval κ * cos κ)) κ :=
    ((Polynomial.hasDerivAt A κ).add ((Polynomial.hasDerivAt B κ).mul (hasDerivAt_cos κ))).add
      ((Polynomial.hasDerivAt D κ).mul (hasDerivAt_sin κ))
  have h := hnum.div (hasDerivAt_pow (n + 1) κ) (pow_ne_zero _ hκ)
  refine h.congr_deriv ?_
  simp only [fourierForm, eval_sub, eval_add, eval_mul, eval_X, eval_C]
  field_simp
  push_cast
  ring

/-- **The tabulated `γ̂'` is the derivative of the tabulated `γ̂`**, for every shape, at every
`κ ≠ 0`. This checks table `ls_fcoef_prime` against `ls_fcoef`, the two places where gromosXX
`latticesum.h` marks earlier versions "wrong!!!" included. -/
theorem gammaHat_hasDerivAt (s : Shape) {κ : ℝ} (hκ : κ ≠ 0) :
    HasDerivAt (gammaHat s) (gammaHatDeriv s κ) κ := by
  cases s with
  | gaussian =>
    have h1 : HasDerivAt (fun x : ℝ => -x ^ 2 / 4) (-(2 * κ) / 4) κ := by
      simpa using ((hasDerivAt_pow 2 κ).neg).div_const 4
    refine h1.exp.congr_deriv ?_
    simp only [gammaHatDeriv]
    ring
  | poly N =>
    have hev : gammaHat (.poly N) =ᶠ[nhds κ]
        fourierForm (fcoef N).1 (fcoef N).2.1 (fcoef N).2.2 (N.toNat + 2 + 1) := by
      filter_upwards [eventually_ne_nhds hκ] with x hx
      simp [gammaHat, hx]
    refine ((fourierForm_hasDerivAt _ _ _ _ hκ).congr_of_eventuallyEq hev).congr_deriv ?_
    simp only [gammaHatDeriv, hκ, ite_false]
    cases N <;> simp only [fcoef, fcoefPrime, fourierForm, TPOrder.toNat] <;> poly_simp <;>
      field_simp <;> ring

/-! ## The lattice-sum energy -/

/-- A position or wave vector. -/
abbrev Vec3 := Fin 3 → ℝ

/-- The real-space interaction of a pair with charge product `q` at distance `r` (Vol. 2,
`e_eta_approx`, first sum): `q η(r/a) / r`. -/
noncomputable def realSpacePair (q a : ℝ) (s : Shape) (r : ℝ) : ℝ := q * eta s (r / a) / r

/-- Its scalar force `−E'(r)/r` in the kernel's form `(q η/r − q η'/a) / r²` (Vol. 2,
`e_eta_approx_force`). -/
noncomputable def realSpacePairForce (q a : ℝ) (s : Shape) (r : ℝ) : ℝ :=
  (q * eta s (r / a) / r - q * etaDeriv s (r / a) / a) / r ^ 2

/-- The real-space term of an excluded pair (Vol. 2, `e_eta_approx`, second sum):
`q (η(r/a) − 1) / r`, the shaped-charge interaction the reciprocal sum includes, without the
Coulomb term the exclusion removes. -/
noncomputable def realSpaceExcluded (q a : ℝ) (s : Shape) (r : ℝ) : ℝ :=
  q * (eta s (r / a) - 1) / r

/-- Its scalar force, `(q (η − 1)/r − q η'/a) / r²`. -/
noncomputable def realSpaceExcludedForce (q a : ℝ) (s : Shape) (r : ℝ) : ℝ :=
  (q * (eta s (r / a) - 1) / r - q * etaDeriv s (r / a) / a) / r ^ 2

/-- **Force = −dE/dr** for the real-space pair, at `r ≠ 0`. -/
theorem realSpacePair_force_is_minus_grad (q : ℝ) {a : ℝ} (s : Shape) {r : ℝ} (hr : r ≠ 0)
    (ha : a ≠ 0) :
    HasDerivAt (realSpacePair q a s) (-(realSpacePairForce q a s r) * r) r := by
  have he : HasDerivAt (fun x => eta s (x / a)) (etaDeriv s (r / a) * (1 / a)) r := by
    have := (eta_hasDerivAt s (r / a)).comp r ((hasDerivAt_id r).div_const a)
    simpa [Function.comp_def, one_div] using this
  have h := (he.const_mul q).div (hasDerivAt_id r) hr
  refine h.congr_deriv ?_
  simp only [realSpacePairForce, id]
  field_simp
  ring

/-- **Force = −dE/dr** for the excluded-pair term, at `r ≠ 0`. -/
theorem realSpaceExcluded_force_is_minus_grad (q : ℝ) {a : ℝ} (s : Shape) {r : ℝ} (hr : r ≠ 0)
    (ha : a ≠ 0) :
    HasDerivAt (realSpaceExcluded q a s) (-(realSpaceExcludedForce q a s r) * r) r := by
  have he : HasDerivAt (fun x => eta s (x / a) - 1) (etaDeriv s (r / a) * (1 / a)) r := by
    have := ((eta_hasDerivAt s (r / a)).comp r ((hasDerivAt_id r).div_const a)).sub_const 1
    simpa [Function.comp_def, one_div] using this
  have h := (he.const_mul q).div (hasDerivAt_id r) hr
  refine h.congr_deriv ?_
  simp only [realSpaceExcludedForce, id]
  field_simp
  ring

/-! ## The wave vectors of a box -/

/-- The volume of the box with edges `a, b, c`, `(a × b) · c` (positive for GROMOS's
right-handed boxes). -/
noncomputable def boxVolume (a b c : Vec3) : ℝ := crossProduct a b ⬝ᵥ c

/-- The reciprocal basis of the box (Vol. 2 after `e_gamma_ewald`; gromosXX
`KSpace_Utils::l_to_k_matrix`): `2π/V · (b × c, c × a, a × b)`. -/
noncomputable def reciprocal (a b c : Vec3) : Fin 3 → Vec3 :=
  ![(2 * π / boxVolume a b c) • crossProduct b c, (2 * π / boxVolume a b c) • crossProduct c a,
    (2 * π / boxVolume a b c) • crossProduct a b]

/-- The vector of the integer triple `l` in the basis `B`: `l₀ B₀ + l₁ B₁ + l₂ B₂`. With
`B = reciprocal a b c` it is the wave vector of `l`; with `B = ![a, b, c]`, a lattice vector of
the box. -/
noncomputable def latticeVector (B : Fin 3 → Vec3) (l : Fin 3 → ℤ) : Vec3 := ∑ j, (l j : ℝ) • B j

theorem latticeVector_neg (B : Fin 3 → Vec3) (l : Fin 3 → ℤ) :
    latticeVector B (-l) = -latticeVector B l := by
  simp [latticeVector, Finset.sum_neg_distrib]

/-- **The reciprocal basis is dual to the box**: the wave vector `k` of `l` has `k · a = 2π l₀`,
`k · b = 2π l₁`, `k · c = 2π l₂`. So `exp(i k · r)` has the period of the box, and an atom moved
by a lattice vector leaves the reciprocal-space sum unchanged
(`structureFactorSq_lattice_shift`). -/
theorem latticeVector_reciprocal_dot {a b c : Vec3} (hV : boxVolume a b c ≠ 0)
    (l : Fin 3 → ℤ) :
    latticeVector (reciprocal a b c) l ⬝ᵥ a = 2 * π * l 0 ∧
      latticeVector (reciprocal a b c) l ⬝ᵥ b = 2 * π * l 1 ∧
      latticeVector (reciprocal a b c) l ⬝ᵥ c = 2 * π * l 2 := by
  set V := boxVolume a b c with hVdef
  have hk : ∀ w : Vec3, latticeVector (reciprocal a b c) l ⬝ᵥ w =
      (l 0 : ℝ) * (2 * π / V) * (crossProduct b c ⬝ᵥ w) +
        (l 1 : ℝ) * (2 * π / V) * (crossProduct c a ⬝ᵥ w) +
        (l 2 : ℝ) * (2 * π / V) * (crossProduct a b ⬝ᵥ w) := by
    intro w
    simp only [latticeVector, reciprocal, Fin.sum_univ_three, add_dotProduct, smul_dotProduct,
      smul_eq_mul, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.head_cons, Matrix.tail_cons, ← hVdef]
    ring
  have triple : ∀ u v w : Vec3, crossProduct u v ⬝ᵥ w =
      (u 1 * v 2 - u 2 * v 1) * w 0 + (u 2 * v 0 - u 0 * v 2) * w 1 +
        (u 0 * v 1 - u 1 * v 0) * w 2 := fun u v w => by
    simp [cross_apply, dotProduct, Fin.sum_univ_three]
  have hVe : V = crossProduct a b ⬝ᵥ c := rfl
  refine ⟨?_, ?_, ?_⟩ <;> rw [hk] <;> simp only [triple] at hVe ⊢ <;> field_simp <;>
    rw [hVe] <;> ring

/-- The integer triples of the Ewald sum, `|l_j| ≤ NK_j` (NONBONDED NKX, NKY, NKZ). -/
def lBox (nk : Fin 3 → ℕ) : Finset (Fin 3 → ℤ) :=
  Fintype.piFinset fun j => Finset.Icc (-(nk j : ℤ)) (nk j)

theorem mem_lBox {nk : Fin 3 → ℕ} {l : Fin 3 → ℤ} :
    l ∈ lBox nk ↔ ∀ j, -(nk j : ℤ) ≤ l j ∧ l j ≤ nk j := by
  simp [lBox, Fintype.mem_piFinset]

theorem neg_mem_lBox {nk : Fin 3 → ℕ} {l : Fin 3 → ℤ} : -l ∈ lBox nk ↔ l ∈ lBox nk := by
  simp only [mem_lBox, Pi.neg_apply]
  exact forall_congr' fun j => by omega

open Classical in
/-- The wave vectors of the Ewald sum (gromosXX `configuration::calculate_k_space`,
`kspace.cc:48-126`): `k = latticeVector B l` for `l` in `lBox nk` with `k ≠ 0` and `|k|² ≤ KCUT`.
gromosXX compares `|k|²`, not `|k|`, with KCUT; Vol. 2 cuts at `|k| ≤ k_c`. -/
noncomputable def kSet (nk : Fin 3 → ℕ) (kcut : ℝ) (B : Fin 3 → Vec3) : Finset Vec3 :=
  ((lBox nk).filter fun l => latticeVector B l ⬝ᵥ latticeVector B l ≤ kcut ∧
    latticeVector B l ⬝ᵥ latticeVector B l ≠ 0).image (latticeVector B)

/-- **The set of wave vectors is symmetric**: `k` is in it exactly when `−k` is. -/
theorem kSet_symmetric (nk : Fin 3 → ℕ) (kcut : ℝ) (B : Fin 3 → Vec3) (k : Vec3) :
    k ∈ kSet nk kcut B ↔ -k ∈ kSet nk kcut B := by
  have one : ∀ k, k ∈ kSet nk kcut B → -k ∈ kSet nk kcut B := by
    intro k hk
    simp only [kSet, Finset.mem_image, Finset.mem_filter] at hk ⊢
    obtain ⟨l, ⟨hl, hc, h0⟩, rfl⟩ := hk
    exact ⟨-l, ⟨neg_mem_lBox.mpr hl, by simpa [latticeVector_neg] using hc,
      by simpa [latticeVector_neg] using h0⟩, latticeVector_neg B l⟩
  exact ⟨one k, fun h => by simpa using one (-k) h⟩

section Reciprocal

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- `C(k) = Σ_i q_i cos(k · r_i)` (Vol. 2, `c_and_s_fcts_def`). -/
noncomputable def cosSum (q : ι → ℝ) (x : ι → Vec3) (k : Vec3) : ℝ :=
  ∑ i, q i * cos (k ⬝ᵥ x i)

/-- `S(k) = Σ_i q_i sin(k · r_i)` (Vol. 2, `c_and_s_fcts_def`). -/
noncomputable def sinSum (q : ι → ℝ) (x : ι → Vec3) (k : Vec3) : ℝ :=
  ∑ i, q i * sin (k ⬝ᵥ x i)

/-- `|Σ_i q_i exp(i k · r_i)|² = C(k)² + S(k)²`. -/
noncomputable def structureFactorSq (q : ι → ℝ) (x : ι → Vec3) (k : Vec3) : ℝ :=
  cosSum q x k ^ 2 + sinSum q x k ^ 2

/-- The reciprocal-space energy over a finite set `K` of wave vectors `k ≠ 0` (Vol. 2,
`e_gamma_ewald`): `(2 ε₀ V)⁻¹ Σ_k γ̂(a|k|)/k² · (C(k)² + S(k)²)`. -/
noncomputable def kSpace (ε0 V a : ℝ) (s : Shape) (K : Finset Vec3) (q : ι → ℝ)
    (x : ι → Vec3) : ℝ :=
  (2 * ε0 * V)⁻¹ * ∑ k ∈ K, gammaHat s (a * √(k ⬝ᵥ k)) / (k ⬝ᵥ k) * structureFactorSq q x k

/-- The reciprocal-space force on atom `i` (Vol. 2, `ewald_reci_force`):
`(ε₀ V)⁻¹ q_i Σ_k γ̂(a|k|)/k² · k · (C(k) sin(k · r_i) − S(k) cos(k · r_i))`. -/
noncomputable def kSpaceForce (ε0 V a : ℝ) (s : Shape) (K : Finset Vec3) (q : ι → ℝ)
    (x : ι → Vec3) (i : ι) : Vec3 :=
  fun α => (ε0 * V)⁻¹ * q i * ∑ k ∈ K, gammaHat s (a * √(k ⬝ᵥ k)) / (k ⬝ᵥ k) * k α *
    (cosSum q x k * sin (k ⬝ᵥ x i) - sinSum q x k * cos (k ⬝ᵥ x i))

/-- Atom `i` moved by `t e`, every other atom in place. -/
def displace (x : ι → Vec3) (i : ι) (e : Vec3) (t : ℝ) : ι → Vec3 :=
  Function.update x i (x i + t • e)

/-- A sum over the atoms in which only the term of atom `i` sees the displacement. -/
theorem sum_displace {β : Type*} [AddCommGroup β] (f : ι → Vec3 → β) (x : ι → Vec3) (i : ι)
    (e : Vec3) (t : ℝ) :
    ∑ j, f j (displace x i e t j) = ∑ j, f j (x j) - f i (x i) + f i (x i + t • e) := by
  rw [← Finset.add_sum_erase _ _ (Finset.mem_univ i), ← Finset.add_sum_erase _ _ (Finset.mem_univ i)]
  have : ∑ j ∈ Finset.univ.erase i, f j (displace x i e t j) =
      ∑ j ∈ Finset.univ.erase i, f j (x j) :=
    Finset.sum_congr rfl fun j hj => by
      simp [displace, Function.update_of_ne (Finset.ne_of_mem_erase hj)]
  rw [this]
  simp only [displace, Function.update_self]
  abel

/-- **The gradient of `C(k)² + S(k)²`**, per wave vector: moving atom `i` along `e` changes it at
the rate `2 q_i (k · e) (S(k) cos(k · r_i) − C(k) sin(k · r_i))`; with `e` a unit vector along a
coordinate axis, that is one component of the gradient. -/
theorem structureFactorSq_hasDerivAt (q : ι → ℝ) (x : ι → Vec3) (k : Vec3) (i : ι) (e : Vec3) :
    HasDerivAt (fun t => structureFactorSq q (displace x i e t) k)
      (-(2 * q i * (k ⬝ᵥ e) *
        (cosSum q x k * sin (k ⬝ᵥ x i) - sinSum q x k * cos (k ⬝ᵥ x i)))) 0 := by
  have hc : ∀ t, cosSum q (displace x i e t) k =
      cosSum q x k - q i * cos (k ⬝ᵥ x i) + q i * cos (k ⬝ᵥ x i + t * (k ⬝ᵥ e)) := by
    intro t
    simp only [cosSum]
    rw [sum_displace (fun j v => q j * cos (k ⬝ᵥ v)), dotProduct_add, dotProduct_smul,
      smul_eq_mul]
  have hs : ∀ t, sinSum q (displace x i e t) k =
      sinSum q x k - q i * sin (k ⬝ᵥ x i) + q i * sin (k ⬝ᵥ x i + t * (k ⬝ᵥ e)) := by
    intro t
    simp only [sinSum]
    rw [sum_displace (fun j v => q j * sin (k ⬝ᵥ v)), dotProduct_add, dotProduct_smul,
      smul_eq_mul]
  have harg : HasDerivAt (fun t : ℝ => k ⬝ᵥ x i + t * (k ⬝ᵥ e)) (k ⬝ᵥ e) 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).mul_const (k ⬝ᵥ e)).const_add (k ⬝ᵥ x i)
  have hC : HasDerivAt (fun t => cosSum q (displace x i e t) k)
      (q i * (-sin (k ⬝ᵥ x i + 0 * (k ⬝ᵥ e)) * (k ⬝ᵥ e))) 0 := by
    simp only [hc]
    exact (harg.cos.const_mul (q i)).const_add _
  have hS : HasDerivAt (fun t => sinSum q (displace x i e t) k)
      (q i * (cos (k ⬝ᵥ x i + 0 * (k ⬝ᵥ e)) * (k ⬝ᵥ e))) 0 := by
    simp only [hs]
    exact (harg.sin.const_mul (q i)).const_add _
  have h := (hC.pow 2).add (hS.pow 2)
  refine (h.congr_deriv ?_).congr_of_eventuallyEq (Filter.Eventually.of_forall fun t => rfl)
  simp only [hc, hs, zero_mul, add_zero]
  ring

/-- **Force = −∇E** for the reciprocal-space energy: moving atom `i` along `e` changes `kSpace`
at the rate `−kSpaceForce i · e`. -/
theorem kSpace_force_is_minus_grad (ε0 V a : ℝ) (s : Shape) (K : Finset Vec3) (q : ι → ℝ)
    (x : ι → Vec3) (i : ι) (e : Vec3) :
    HasDerivAt (fun t => kSpace ε0 V a s K q (displace x i e t))
      (-(kSpaceForce ε0 V a s K q x i ⬝ᵥ e)) 0 := by
  have hF : kSpaceForce ε0 V a s K q x i ⬝ᵥ e = (ε0 * V)⁻¹ * q i * ∑ k ∈ K,
      gammaHat s (a * √(k ⬝ᵥ k)) / (k ⬝ᵥ k) * (k ⬝ᵥ e) *
        (cosSum q x k * sin (k ⬝ᵥ x i) - sinSum q x k * cos (k ⬝ᵥ x i)) := by
    have hdot : ∀ k : Vec3, k ⬝ᵥ e = ∑ α, k α * e α := fun _ => rfl
    simp only [hdot]
    simp only [kSpaceForce, Finset.mul_sum, Finset.sum_mul]
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun α _ => by ring
  have h := (HasDerivAt.fun_sum (u := K) (fun k _ =>
    (structureFactorSq_hasDerivAt q x k i e).const_mul
      (gammaHat s (a * √(k ⬝ᵥ k)) / (k ⬝ᵥ k)))).const_mul (2 * ε0 * V)⁻¹
  refine h.congr_deriv ?_
  have hsum : ∑ k ∈ K, gammaHat s (a * √(k ⬝ᵥ k)) / (k ⬝ᵥ k) *
      -(2 * q i * (k ⬝ᵥ e) * (cosSum q x k * sin (k ⬝ᵥ x i) - sinSum q x k * cos (k ⬝ᵥ x i))) =
      -(2 * q i) * ∑ k ∈ K, gammaHat s (a * √(k ⬝ᵥ k)) / (k ⬝ᵥ k) * (k ⬝ᵥ e) *
        (cosSum q x k * sin (k ⬝ᵥ x i) - sinSum q x k * cos (k ⬝ᵥ x i)) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun k _ => by ring
  rw [hF, hsum]
  ring

omit [DecidableEq ι] in
/-- **Translation invariance**: shifting every atom by the same vector `d` leaves
`C(k)² + S(k)²`, hence the reciprocal-space energy, unchanged. -/
theorem kSpace_translation_invariant (q : ι → ℝ) (x : ι → Vec3) (k d : Vec3) :
    structureFactorSq q (fun j => x j + d) k = structureFactorSq q x k := by
  have hc : cosSum q (fun j => x j + d) k =
      cosSum q x k * cos (k ⬝ᵥ d) - sinSum q x k * sin (k ⬝ᵥ d) := by
    simp only [cosSum, sinSum, dotProduct_add, cos_add, Finset.sum_mul, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  have hs : sinSum q (fun j => x j + d) k =
      sinSum q x k * cos (k ⬝ᵥ d) + cosSum q x k * sin (k ⬝ᵥ d) := by
    simp only [cosSum, sinSum, dotProduct_add, sin_add, Finset.sum_mul, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  simp only [structureFactorSq, hc, hs]
  linear_combination (cosSum q x k ^ 2 + sinSum q x k ^ 2) * sin_sq_add_cos_sq (k ⬝ᵥ d)

omit [DecidableEq ι] in
theorem cosSum_neg (q : ι → ℝ) (x : ι → Vec3) (k : Vec3) : cosSum q x (-k) = cosSum q x k := by
  simp [cosSum, neg_dotProduct, Real.cos_neg]

omit [DecidableEq ι] in
theorem sinSum_neg (q : ι → ℝ) (x : ι → Vec3) (k : Vec3) : sinSum q x (-k) = -sinSum q x k := by
  simp [sinSum, neg_dotProduct, Real.sin_neg, Finset.sum_neg_distrib]

omit [DecidableEq ι] in
/-- **`C(k)² + S(k)²` is even in `k`**: `C(−k) = C(k)`, `S(−k) = −S(k)`. -/
theorem structureFactorSq_neg (q : ι → ℝ) (x : ι → Vec3) (k : Vec3) :
    structureFactorSq q x (-k) = structureFactorSq q x k := by
  simp [structureFactorSq, cosSum_neg, sinSum_neg]

open Pointwise in
/-- A sum over a set and its mirror image of an even function is twice the sum over the set. -/
theorem sum_union_neg_eq_two_mul {β : Type*} [AddCommGroup β] [DecidableEq β] {H : Finset β}
    (hd : Disjoint H (-H)) (f : β → ℝ) (hf : ∀ x, f (-x) = f x) :
    ∑ x ∈ H ∪ -H, f x = 2 * ∑ x ∈ H, f x := by
  rw [Finset.sum_union hd, ← Finset.image_neg_eq_neg,
    Finset.sum_image (fun x _ y _ h => neg_injective h)]
  simp only [hf]
  ring

omit [DecidableEq ι] in
open Pointwise in
/-- **The ½ of `kSpace`**: over a half `H` of the wave vectors and its mirror image `−H`, the
reciprocal-space energy is `(ε₀V)⁻¹ Σ_{k ∈ H}`, the half-space sum without the ½. -/
theorem kSpace_half_space (ε0 V a : ℝ) (s : Shape) {H : Finset Vec3} (hd : Disjoint H (-H))
    (q : ι → ℝ) (x : ι → Vec3) :
    kSpace ε0 V a s (H ∪ -H) q x =
      (ε0 * V)⁻¹ * ∑ k ∈ H, gammaHat s (a * √(k ⬝ᵥ k)) / (k ⬝ᵥ k) * structureFactorSq q x k := by
  rw [kSpace, sum_union_neg_eq_two_mul hd]
  · ring
  · intro k
    simp [structureFactorSq_neg, neg_dotProduct, dotProduct_neg]

/-- **Force = −∇E for the Ewald sum**: over the wave vectors of `kSet`, moving atom `i` along `e`
changes the reciprocal-space energy at the rate `−F_i · e`, `F_i` the sum of the per-wave-vector
forces (`kSpace_force_is_minus_grad`, `structureFactorSq_hasDerivAt`). -/
theorem ewald_total_force_is_minus_grad (ε0 V a : ℝ) (s : Shape) (nk : Fin 3 → ℕ) (kcut : ℝ)
    (B : Fin 3 → Vec3) (q : ι → ℝ) (x : ι → Vec3) (i : ι) (e : Vec3) :
    HasDerivAt (fun t => kSpace ε0 V a s (kSet nk kcut B) q (displace x i e t))
      (-(kSpaceForce ε0 V a s (kSet nk kcut B) q x i ⬝ᵥ e)) 0 :=
  kSpace_force_is_minus_grad ε0 V a s _ q x i e

/-- **Periodicity**: moving one atom by a lattice vector `n₀ a + n₁ b + n₂ c` of the box leaves
`C(k)² + S(k)²` unchanged at every wave vector `k` of the reciprocal lattice
(`latticeVector_reciprocal_dot`). gromosXX's wave vectors in a triclinic box are not of this
form (it multiplies `l` by the transpose of the reciprocal basis, `kspace.cc:88`). -/
theorem structureFactorSq_lattice_shift {a b c : Vec3} (hV : boxVolume a b c ≠ 0)
    (l n : Fin 3 → ℤ) (q : ι → ℝ) (x : ι → Vec3) (i : ι) :
    structureFactorSq q (Function.update x i (x i + latticeVector ![a, b, c] n))
        (latticeVector (reciprocal a b c) l) =
      structureFactorSq q x (latticeVector (reciprocal a b c) l) := by
  set k := latticeVector (reciprocal a b c) l
  obtain ⟨ha, hb, hc⟩ := latticeVector_reciprocal_dot hV l
  have hshift : k ⬝ᵥ latticeVector ![a, b, c] n =
      ((l 0 * n 0 + l 1 * n 1 + l 2 * n 2 : ℤ) : ℝ) * (2 * π) := by
    simp only [latticeVector, Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, dotProduct_add, dotProduct_smul,
      smul_eq_mul]
    rw [ha, hb, hc]
    push_cast
    ring
  have hterm : ∀ j, k ⬝ᵥ Function.update x i (x i + latticeVector ![a, b, c] n) j =
      k ⬝ᵥ x j ∨ k ⬝ᵥ Function.update x i (x i + latticeVector ![a, b, c] n) j =
        k ⬝ᵥ x j + ((l 0 * n 0 + l 1 * n 1 + l 2 * n 2 : ℤ) : ℝ) * (2 * π) := by
    intro j
    by_cases hj : j = i
    · subst hj
      right
      simp [dotProduct_add, hshift]
    · left
      simp [Function.update_of_ne hj]
  have hc' : cosSum q (Function.update x i (x i + latticeVector ![a, b, c] n)) k = cosSum q x k := by
    refine Finset.sum_congr rfl fun j _ => ?_
    rcases hterm j with h | h <;> rw [h]
    rw [Real.cos_add_int_mul_two_pi]
  have hs' : sinSum q (Function.update x i (x i + latticeVector ![a, b, c] n)) k = sinSum q x k := by
    refine Finset.sum_congr rfl fun j _ => ?_
    rcases hterm j with h | h <;> rw [h]
    rw [Real.sin_add_int_mul_two_pi]
  simp only [structureFactorSq, hc', hs']

/-! ## Constant terms -/

/-- The box dipole `M = Σ_i q_i (r_i − r_c)`, `r_c` the centre of the box (Vol. 2,
`latsum_dip`). -/
noncomputable def dipole (q : ι → ℝ) (x : ι → Vec3) (c : Vec3) : Vec3 := ∑ i, q i • (x i - c)

/-- The surface term (Vol. 2, `e_srf`): `|M|² / (2 (2ε_LS + 1) ε₀ V)`, for a surrounding
dielectric of relative permittivity `ε_LS` (GROMOS input `EPSLS = 0` stands for `ε_LS = ∞`,
conducting boundaries, where the term vanishes). -/
noncomputable def surface (ε0 V εLS : ℝ) (q : ι → ℝ) (x : ι → Vec3) (c : Vec3) : ℝ :=
  dipole q x c ⬝ᵥ dipole q x c / (2 * (2 * εLS + 1) * ε0 * V)

/-- The surface force on atom `i` (Vol. 2, `e_srf_forces`): `−q_i M / ((2ε_LS + 1) ε₀ V)`. -/
noncomputable def surfaceForce (ε0 V εLS : ℝ) (q : ι → ℝ) (x : ι → Vec3) (c : Vec3) (i : ι) :
    Vec3 :=
  -(q i / ((2 * εLS + 1) * ε0 * V)) • dipole q x c

/-- **Force = −∇E** for the surface term: moving atom `i` along `e` changes it at the rate
`−surfaceForce i · e`. -/
theorem surface_force_is_minus_grad (ε0 V εLS : ℝ) (q : ι → ℝ) (x : ι → Vec3) (c : Vec3)
    (i : ι) (e : Vec3) :
    HasDerivAt (fun t => surface ε0 V εLS q (displace x i e t) c)
      (-(surfaceForce ε0 V εLS q x c i ⬝ᵥ e)) 0 := by
  set M := dipole q x c
  have hM : ∀ t, dipole q (displace x i e t) c = M + (t * q i) • e := by
    intro t
    simp only [dipole, M]
    rw [sum_displace (fun j v => q j • (v - c))]
    ext α
    simp only [Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
    ring
  have hsq : ∀ t, dipole q (displace x i e t) c ⬝ᵥ dipole q (displace x i e t) c =
      M ⬝ᵥ M + t * (2 * q i * (M ⬝ᵥ e)) + t ^ 2 * (q i ^ 2 * (e ⬝ᵥ e)) := by
    intro t
    rw [hM, add_dotProduct, dotProduct_add, dotProduct_add, dotProduct_smul, smul_dotProduct,
      smul_dotProduct, dotProduct_smul, dotProduct_comm e M]
    simp only [smul_eq_mul]
    ring
  have hpoly : HasDerivAt (fun t : ℝ => M ⬝ᵥ M + t * (2 * q i * (M ⬝ᵥ e)) +
      t ^ 2 * (q i ^ 2 * (e ⬝ᵥ e))) (2 * q i * (M ⬝ᵥ e)) 0 := by
    have := (((hasDerivAt_id (0 : ℝ)).mul_const (2 * q i * (M ⬝ᵥ e))).const_add (M ⬝ᵥ M)).fun_add
      ((hasDerivAt_pow 2 (0 : ℝ)).mul_const (q i ^ 2 * (e ⬝ᵥ e)))
    simpa using this
  have h := hpoly.div_const (2 * (2 * εLS + 1) * ε0 * V)
  simp only [surface]
  simp only [hsq]
  refine h.congr_deriv ?_
  simp only [surfaceForce, smul_dotProduct, smul_eq_mul]
  by_cases h0 : (2 * εLS + 1) * ε0 * V = 0
  · have : 2 * (2 * εLS + 1) * ε0 * V = 0 := by linear_combination 2 * h0
    simp [this, h0]
  · field_simp
    ring

/-- The `A`-term (Vol. 2, `e_a_def` with the computed `Ã2` in place of `A2`):
`(A1 S² − (A1 + Ã2) S̃²) / (8π ε₀)`, `S = Σ q_i`, `S̃² = Σ q_i²`. -/
noncomputable def aTerm (ε0 A1 A2t : ℝ) (q : ι → ℝ) : ℝ :=
  (A1 * (∑ i, q i) ^ 2 - (A1 + A2t) * ∑ i, q i ^ 2) / (8 * π * ε0)

/-- The self term (Vol. 2, `e_slf`): `(A1 + A2 + A3) S̃² / (8π ε₀)`. -/
noncomputable def selfTerm (ε0 A1 A2 A3 : ℝ) (q : ι → ℝ) : ℝ :=
  (A1 + A2 + A3) * (∑ i, q i ^ 2) / (8 * π * ε0)

end Reciprocal

/-- `−V π⁻¹ a⁻² A1` per shape (Vol. 2, table `ls_acst`, first column; gromosXX
`nonbonded_outerloop.cc:2291-2348`). -/
noncomputable def a1Tilde : Shape → ℝ
  | .gaussian => 1
  | .poly .o0 => 2 / 5 | .poly .o1 => 4 / 15 | .poly .o2 => 4 / 21 | .poly .o3 => 3 / 14
  | .poly .o4 => 1 / 6 | .poly .o5 => 2 / 15 | .poly .o6 => 8 / 55 | .poly .o7 => 4 / 33
  | .poly .o8 => 4 / 39 | .poly .o9 => 10 / 91 | .poly .o10 => 2 / 21

/-- `−a A3` per shape (Vol. 2, table `ls_acst`, second column). -/
noncomputable def a3Tilde : Shape → ℝ
  | .gaussian => 2 / √π
  | .poly .o0 => 3 / 2 | .poly .o1 => 2 | .poly .o2 => 5 / 2 | .poly .o3 => 9 / 4
  | .poly .o4 => 21 / 8 | .poly .o5 => 3 | .poly .o6 => 45 / 16 | .poly .o7 => 25 / 8
  | .poly .o8 => 55 / 16 | .poly .o9 => 105 / 32 | .poly .o10 => 455 / 128

/-- `A1 = −π a² Ã1 / V` (Vol. 2, `latsum_a1_def`). -/
noncomputable def a1 (s : Shape) (a V : ℝ) : ℝ := -(π * a ^ 2 / V) * a1Tilde s

/-- `A3 = −Ã3 / a` (Vol. 2, `latsum_a3_def`). -/
noncomputable def a3 (s : Shape) (a : ℝ) : ℝ := -(a3Tilde s / a)

/-- The table's `−a A3` is `−η'(0)`: in `latsum_a3_def` the image `n = 0` contributes
`lim (η(r/a) − 1)/r = η'(0)/a` (the other images vanish for a truncated polynomial of width at
most half the box). -/
theorem a3Tilde_eq_neg_etaDeriv_zero (s : Shape) : a3Tilde s = -etaDeriv s 0 := by
  cases s with
  | gaussian => simp [a3Tilde, etaDeriv]
  | poly N => cases N <;> simp [a3Tilde, etaDeriv, etaPolyDeriv, heaviside] <;> norm_num

/-- The Wigner constant of the simple cubic lattice, `ξ_EW = −2.83729748`: `A1 + A2 + A3 = ξ_EW / L`
in a cube of edge `L` (gromosXX `nonbonded_outerloop.cc:2379`). Vol. 2 prints `−2.83792748`
(after `e_slf`), a transposition of two digits. -/
noncomputable def wignerCube : ℝ := -2.83729748

/-! ## `A2` by shells of the integer lattice -/

/-- The integer triples with `|l_j| ≤ n`. -/
def cube (n : ℕ) : Finset (Fin 3 → ℤ) := lBox fun _ => n

/-- Shell `n ≥ 1` of the integer lattice: `max_j |l_j| = n`. -/
def shell (n : ℕ) : Finset (Fin 3 → ℤ) := cube n \ cube (n - 1)

/-- gromosXX's half of shell `n` (`nonbonded_outerloop.cc:2405-2459`): the face `l₀ = n`, the
face `l₁ = n` without the edges on `l₀ = ±n`, the face `l₂ = n` without those on `l₀ = ±n` or
`l₁ = ±n`. -/
def face (n : ℕ) : Finset (Fin 3 → ℤ) :=
  (cube n).filter fun l => l 0 = n ∨ (l 1 = n ∧ -(n : ℤ) < l 0 ∧ l 0 < n) ∨
    (l 2 = n ∧ -(n : ℤ) < l 0 ∧ l 0 < n ∧ -(n : ℤ) < l 1 ∧ l 1 < n)

theorem mem_cube {n : ℕ} {l : Fin 3 → ℤ} : l ∈ cube n ↔ ∀ j, -(n : ℤ) ≤ l j ∧ l j ≤ n :=
  mem_lBox

theorem cube_mono {m n : ℕ} (h : m ≤ n) : cube m ⊆ cube n := fun l hl => by
  rw [mem_cube] at hl ⊢
  exact fun j => ⟨by have := (hl j).1; omega, by have := (hl j).2; omega⟩

/-- **The shells partition the lattice**: a sum over the cube `|l_j| ≤ N` without `l = 0` is the
sum over shells `1 … N`. -/
theorem sum_cube_erase_zero (f : (Fin 3 → ℤ) → ℝ) (N : ℕ) :
    ∑ l ∈ (cube N).erase 0, f l = ∑ n ∈ Finset.Icc 1 N, ∑ l ∈ shell n, f l := by
  induction N with
  | zero =>
    have : (cube 0).erase 0 = ∅ := by
      ext l
      simp only [Finset.mem_erase, mem_cube, Finset.notMem_empty, iff_false, not_and]
      intro hne h
      apply hne
      funext j
      have := h j
      simp only [Pi.zero_apply]
      omega
    simp [this]
  | succ N ih =>
    have hsub : (cube N).erase 0 ⊆ (cube (N + 1)).erase 0 :=
      Finset.erase_subset_erase 0 (cube_mono (Nat.le_succ N))
    have h0 : (0 : Fin 3 → ℤ) ∈ cube N := by
      rw [mem_cube]
      intro j
      simp
    have hdiff : (cube (N + 1)).erase 0 \ (cube N).erase 0 = shell (N + 1) := by
      ext l
      simp only [Finset.mem_sdiff, Finset.mem_erase, shell, Nat.add_sub_cancel]
      constructor
      · rintro ⟨⟨hne, h1⟩, h2⟩
        exact ⟨h1, fun h => h2 ⟨hne, h⟩⟩
      · rintro ⟨h1, h2⟩
        refine ⟨⟨?_, h1⟩, fun h => h2 h.2⟩
        rintro rfl
        exact h2 h0
    rw [← Finset.sum_sdiff hsub, hdiff, ih, Finset.sum_Icc_succ_top (by omega), add_comm]

open Pointwise in
/-- **A shell is gromosXX's half shell and its mirror image**, which do not overlap. -/
theorem shell_eq_face_union_neg {n : ℕ} (hn : 1 ≤ n) :
    shell n = face n ∪ -face n ∧ Disjoint (face n) (-face n) := by
  have hmemneg : ∀ l, l ∈ -face n ↔ -l ∈ face n := fun l => by
    rw [Finset.mem_neg]
    constructor
    · rintro ⟨m, hm, rfl⟩
      simpa using hm
    · intro h
      exact ⟨-l, h, neg_neg l⟩
  have three : ∀ (P : Fin 3 → Prop), (∀ j, P j) ↔ P 0 ∧ P 1 ∧ P 2 := fun P =>
    ⟨fun h => ⟨h 0, h 1, h 2⟩, fun ⟨h0, h1, h2⟩ j => by fin_cases j <;> assumption⟩
  refine ⟨?_, ?_⟩
  · ext l
    rw [Finset.mem_union, hmemneg]
    simp only [shell, face, Finset.mem_sdiff, Finset.mem_filter, mem_cube, three, Pi.neg_apply,
      Nat.cast_sub hn, Nat.cast_one]
    omega
  · rw [Finset.disjoint_left]
    intro l h1 h2
    rw [hmemneg] at h2
    simp only [face, Finset.mem_filter, mem_cube, three, Pi.neg_apply] at h1 h2
    omega

/-- The term of `l` in `A2`: `γ̂(a|k|)/k²`, `k = latticeVector B l`. -/
noncomputable def a2Term (s : Shape) (a : ℝ) (B : Fin 3 → Vec3) (l : Fin 3 → ℤ) : ℝ :=
  gammaHat s (a * √(latticeVector B l ⬝ᵥ latticeVector B l)) /
    (latticeVector B l ⬝ᵥ latticeVector B l)

/-- `A2` over the shells `1 … N` as gromosXX sums it (`nonbonded_outerloop.cc:2387-2522`): twice
the half shell, `4π/V` at the end (Vol. 2, `latsum_a2_calc`). -/
noncomputable def a2Shells (s : Shape) (a : ℝ) (B : Fin 3 → Vec3) (V : ℝ) (N : ℕ) : ℝ :=
  4 * π / V * ∑ n ∈ Finset.Icc 1 N, 2 * ∑ l ∈ face n, a2Term s a B l

/-- **gromosXX's shells are the lattice sum**: `A2` over `N` shells is `4π/V Σ γ̂/k²` over every
`l ≠ 0` with `|l_j| ≤ N`. -/
theorem a2Shells_eq_cube (s : Shape) (a : ℝ) (B : Fin 3 → Vec3) (V : ℝ) (N : ℕ) :
    a2Shells s a B V N = 4 * π / V * ∑ l ∈ (cube N).erase 0, a2Term s a B l := by
  rw [a2Shells, sum_cube_erase_zero]
  congr 1
  refine Finset.sum_congr rfl fun n hn => ?_
  obtain ⟨hshell, hd⟩ := shell_eq_face_union_neg (Finset.mem_Icc.mp hn).1
  rw [hshell, sum_union_neg_eq_two_mul hd]
  intro l
  simp [a2Term, latticeVector_neg]

/-- `A2` in a cube of edge `L` from the Wigner constant (gromosXX
`nonbonded_outerloop.cc:2376-2385`): `A1 + A2 + A3 = ξ_EW / L`. -/
noncomputable def a2Cube (W L A1 A3 : ℝ) : ℝ := W / L - A1 - A3

theorem a1_add_a2Cube_add_a3 (W L A1 A3 : ℝ) : A1 + a2Cube W L A1 A3 + A3 = W / L := by
  simp only [a2Cube]
  ring

end NewenTheory.Interactions.LatticeSum

end
