module

public import Mathlib.Analysis.Fourier.ZMod
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Sinc
public import NewenTheory.Interactions.LatticeSum
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring
import Mathlib.Tactic.IntervalCases

/-!
# P3M (GROMOS Vol. 2, particle–particle particle–mesh)

P3M replaces the reciprocal-space sum of the Ewald method by a sum over a grid (Vol. 2,
`assign_p3m` … `force_from_field_p3m`): each charge is spread onto the grid points near it with an
assignment function `w_p` of order `p` (a B-spline), the grid charge density is Fourier
transformed by a discrete Fourier transform, multiplied by the optimal influence function `Ĝ`
(`optim_influ`), and the field comes back on the grid by `ik` or by a finite-difference operator
(`finite_diff`, `exact_diff`) and is interpolated to the charges with `w_p` again.

What is defined and proved here:

* the assignment functions `w p` for `p = 1 … 5` as gromosXX writes them
  (`charge_assignment_1d`), non-negative (`w_nonneg`), zero for `|x| ≥ p/2` (`w_support`), and a
  partition of unity: the weights a charge gives the grid add up to one
  (`w_partition_of_unity`), except, for `p = 1`, at a charge exactly half-way between two grid
  points, where the strict `|x| < 1/2` of the code gives both points weight zero
  (`w1_half_integer`);
* the discrete Fourier transform of FFTW's convention on `ZMod N` (Mathlib's `ZMod.dft`, forward
  with `e^{−2πi jk/N}`, both directions unnormalised): the backward transform of the forward one
  is `N` times the original (`dft_inv`), a shift multiplies the transform by a phase
  (`dft_shift`), and a real grid has a Hermitian transform (`dft_neg_eq_conj`);
* `sincPow x p = (sin(x/2)/(x/2))^p`, the transform of `w_p` (Vol. 2 `wei_product`), `1` at `0`
  (`sincPow_zero`), Mathlib's `sinc (x/2) ^ p` (`sincPow_eq_sinc`), so continuous there
  (`sincPow_tendsto_zero`);
* the coefficients of the finite-difference operators (Vol. 2 `tab:fd_coef`) and what they are
  for: the operator of order `q` differentiates every polynomial of degree `≤ 2q` exactly
  (`fd_exact_on_polynomials`, from `fd_exact_1` … `fd_exact_5`);
* the optimal influence function (Vol. 2 `optim_influ`) over a finite set of alias vectors, and
  for a single alias with `ik` differentiation its value `γ̂(ak)/(k² P̂(k)²)`
  (`influence_single_alias`);
* the k-space energy `½ V_G ε₀⁻¹ Σ_l Ĝ_l |ρ̂_l|²/N` and that the complex form
  `Σ Ĝ_l ρ̂_l ρ̂_{−l}` it comes from is real for a real charge density (`p3m_energy_real`).

That P3M approximates the Ewald sum is the method's motivation (Hockney and Eastwood; Vol. 2) and
is not proved here; the tests compare the two numerically.
-/

@[expose] public section

namespace NewenTheory.Interactions.P3M

open Real
open NewenTheory.Interactions.LatticeSum (Vec3 Shape gammaHat)

/-! ## The assignment functions -/

/-- The one-dimensional charge-assignment function `w_p(x)` of order `p`, as gromosXX's
`Lattice_Sum::charge_assignment_1d` (`latticesum.h:508-565`) and Vol. 2's table of weights; zero
for every other `p`. -/
noncomputable def w (p : ℕ) (x : ℝ) : ℝ :=
  match p with
  | 1 => if |x| < 1 / 2 then 1 else 0
  | 2 => if |x| < 1 then 1 - |x| else 0
  | 3 =>
    if |x| < 1 / 2 then 3 / 4 - x ^ 2
    else if |x| < 3 / 2 then (2 * |x| - 3) ^ 2 / 8 else 0
  | 4 =>
    if |x| < 1 then (3 * |x| ^ 3 - 6 * x ^ 2 + 4) / 6
    else if |x| < 2 then -(|x| - 2) ^ 3 / 6 else 0
  | 5 =>
    if |x| < 1 / 2 then (48 * x ^ 4 - 120 * x ^ 2 + 115) / 192
    else if |x| < 3 / 2 then -(16 * x ^ 4 - 80 * |x| ^ 3 + 120 * x ^ 2 - 20 * |x| - 55) / 96
    else if |x| < 5 / 2 then (2 * |x| - 5) ^ 4 / 384 else 0
  | _ => 0

theorem w_nonneg {p : ℕ} (x : ℝ) : 0 ≤ w p x := by
  have ha := abs_nonneg x
  have h2 : x ^ 2 = |x| ^ 2 := (sq_abs x).symm
  have h4 : x ^ 4 = |x| ^ 4 := (Even.pow_abs ⟨2, rfl⟩ x).symm
  rcases p with _ | _ | _ | _ | _ | _ | p <;> simp only [w, h2, h4] <;> (try split_ifs) <;>
    (try push Not at *) <;>
    first
    | positivity
    | nlinarith
    -- order 4, `1 ≤ |x| < 2`: `−(|x| − 2)³ ≥ 0`
    | nlinarith [pow_nonneg (by linarith : (0 : ℝ) ≤ 2 - |x|) 3]
    -- order 5, `1/2 ≤ |x| < 3/2`
    | nlinarith [mul_nonneg (by linarith : (0 : ℝ) ≤ |x| - 1 / 2) (by linarith : (0 : ℝ) ≤ 3 / 2 - |x|),
        mul_nonneg (mul_nonneg (by linarith : (0 : ℝ) ≤ |x| - 1 / 2)
          (by linarith : (0 : ℝ) ≤ 3 / 2 - |x|)) (sq_nonneg (|x| - 1)),
        sq_nonneg (|x| - 1)]

theorem w_support {p : ℕ} {x : ℝ} (hx : (p : ℝ) / 2 ≤ |x|) : w p x = 0 := by
  rcases p with _ | _ | _ | _ | _ | _ | p <;> norm_num at hx <;> simp only [w] <;> (try split_ifs) <;>
    first | rfl | linarith

/-- `if a < b` with `a < b` known. -/
theorem if_lt_of_lt {a b x y : ℝ} (h : a < b) : (if a < b then x else y) = x :=
  ite_eq_left_iff.mpr fun h' => absurd h h'

/-- `if a < b` with `b ≤ a` known. -/
theorem if_lt_of_le {a b x y : ℝ} (h : b ≤ a) : (if a < b then x else y) = y :=
  ite_eq_right_iff.mpr fun h' => absurd h' (not_lt.mpr h)

/-- Decide every `if a < b` of the goal by `linarith`, outermost first. -/
macro "decide_ifs" : tactic =>
  `(tactic| repeat (first | rw [if_lt_of_lt (by linarith)] | rw [if_lt_of_le (by linarith)]))

/-- `|t + j|` for `t ∈ [0, 1)` and `j = 3 … −3`. -/
theorem abs_window {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) :
    |t + 3| = t + 3 ∧ |t + 2| = t + 2 ∧ |t + 1| = t + 1 ∧ |t| = t ∧ |t - 1| = 1 - t ∧
      |t - 2| = 2 - t ∧ |t - 3| = 3 - t := by
  refine ⟨abs_of_nonneg (by linarith), abs_of_nonneg (by linarith), abs_of_nonneg (by linarith),
    abs_of_nonneg h0, ?_, ?_, ?_⟩ <;> rw [abs_of_neg (by linarith)] <;> ring

/-- `t ∈ [0, 1)` is `0`, in `(0, 1/2)`, `1/2` or in `(1/2, 1)`: inside the open pieces every
branch of `w p (t + j)` is decided by `linarith`, at the two points by arithmetic. -/
theorem window_cases {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) :
    t = 0 ∨ (0 < t ∧ t < 1 / 2) ∨ t = 1 / 2 ∨ (1 / 2 < t ∧ t < 1) := by
  rcases h0.eq_or_lt with h | h
  · exact Or.inl h.symm
  rcases lt_trichotomy t (1 / 2) with h' | h' | h'
  · exact Or.inr (Or.inl ⟨h, h'⟩)
  · exact Or.inr (Or.inr (Or.inl h'))
  · exact Or.inr (Or.inr (Or.inr ⟨h', h1⟩))

/-- The seven weights of a charge at `t ∈ [0, 1)` (grid units) from the grid points `−3 … 3`
add up to one, `p = 2 … 5`. -/
theorem w_window {p : ℕ} (hp2 : 2 ≤ p) (hp5 : p ≤ 5) {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) :
    w p (t + 3) + w p (t + 2) + w p (t + 1) + w p t + w p (t - 1) + w p (t - 2) + w p (t - 3) =
      1 := by
  obtain ⟨a0, a1, a2, a3, a4, a5, a6⟩ := abs_window h0 h1
  rcases window_cases h0 h1 with rfl | ⟨hl, hr⟩ | rfl | ⟨hl, hr⟩ <;> interval_cases p <;>
    simp only [w, a0, a1, a2, a3, a4, a5, a6] <;> decide_ifs <;> ring

/-- The same for `p = 1`, away from `t = 1/2`. -/
theorem w1_window {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) (hh : t ≠ 1 / 2) :
    w 1 (t + 3) + w 1 (t + 2) + w 1 (t + 1) + w 1 t + w 1 (t - 1) + w 1 (t - 2) + w 1 (t - 3) =
      1 := by
  obtain ⟨a0, a1, a2, a3, a4, a5, a6⟩ := abs_window h0 h1
  rcases window_cases h0 h1 with rfl | ⟨hl, hr⟩ | rfl | ⟨hl, hr⟩
  all_goals first
    | exact absurd rfl hh
    | (simp only [w, a0, a1, a2, a3, a4, a5, a6]; decide_ifs; ring)

/-- The window sum around `⌊x⌋` is the window sum of the fractional part `t = x − ⌊x⌋`. -/
theorem window_sum_eq (p : ℕ) (x : ℝ) :
    ∑ k ∈ Finset.range 7, w p (x - ((⌊x⌋ - 3 + k : ℤ) : ℝ)) =
      let t := x - ⌊x⌋
      w p (t + 3) + w p (t + 2) + w p (t + 1) + w p t + w p (t - 1) + w p (t - 2) +
        w p (t - 3) := by
  simp only [Finset.sum_range_succ, Finset.sum_range_zero]
  push_cast
  ring_nf

/-- **Partition of unity**: the seven grid points `⌊x⌋ − 3 … ⌊x⌋ + 3` around a charge at `x`
(grid units) carry weights adding up to one, for `p = 2 … 5` (no other point carries any,
`w_support`). -/
theorem w_partition_of_unity {p : ℕ} (hp2 : 2 ≤ p) (hp5 : p ≤ 5) (x : ℝ) :
    ∑ k ∈ Finset.range 7, w p (x - ((⌊x⌋ - 3 + k : ℤ) : ℝ)) = 1 := by
  rw [window_sum_eq]
  exact w_window hp2 hp5 (Int.fract_nonneg x) (Int.fract_lt_one x)

/-- The same for `p = 1` (nearest grid point), away from the half-integers. -/
theorem w1_partition_of_unity {x : ℝ} (hx : ∀ n : ℤ, x ≠ n + 1 / 2) :
    ∑ k ∈ Finset.range 7, w 1 (x - ((⌊x⌋ - 3 + k : ℤ) : ℝ)) = 1 := by
  rw [window_sum_eq]
  refine w1_window (Int.fract_nonneg x) (Int.fract_lt_one x) fun h => hx ⌊x⌋ ?_
  have : x - ⌊x⌋ = 1 / 2 := h
  linarith

/-- At a half-integer, `w 1` gives no grid point any weight: gromosXX's strict `|x| < 1/2`
loses a charge that sits exactly half-way between two grid points. -/
theorem w1_half_integer (n : ℤ) :
    ∑ k ∈ Finset.range 7, w 1 ((n : ℝ) + 1 / 2 - ((n - 3 + k : ℤ) : ℝ)) = 0 := by
  simp only [Finset.sum_range_succ, Finset.sum_range_zero]
  push_cast
  ring_nf
  norm_num [w, abs_of_nonneg, abs_of_nonpos]

/-! ## The discrete Fourier transform -/

open ZMod in
/-- FFTW's backward transform (unnormalised, `e^{+2πi jk/N}`) through Mathlib's forward one:
`Σ_k Ψ_k e^{2πi jk/N} = 𝓕 Ψ (−j)`. -/
noncomputable def dftBackward {N : ℕ} [NeZero N] (Ψ : ZMod N → ℂ) (j : ZMod N) : ℂ :=
  𝓕 Ψ (-j)

open ZMod in
/-- **Inversion**: the backward transform of the forward transform is `N` times the grid. -/
theorem dft_inv {N : ℕ} [NeZero N] (Φ : ZMod N → ℂ) :
    dftBackward (𝓕 Φ) = fun j => (N : ℂ) * Φ j := by
  funext j
  simp [dftBackward, dft_dft]

open ZMod in
/-- `stdAddChar` lands on the unit circle: its conjugate is its inverse. -/
theorem stdAddChar_conj {N : ℕ} [NeZero N] (x : ZMod N) :
    (starRingEnd ℂ) (stdAddChar x) = stdAddChar (-x) := by
  rw [AddChar.map_neg_eq_inv, stdAddChar_apply, ← Circle.coe_inv, Circle.coe_inv_eq_conj]

open ZMod in
/-- **Shift**: moving the grid by `a` multiplies its transform by `e^{2πi ak/N}`. -/
theorem dft_shift {N : ℕ} [NeZero N] (Φ : ZMod N → ℂ) (a k : ZMod N) :
    𝓕 (fun j => Φ (j + a)) k = (stdAddChar (a * k) : ℂ) * 𝓕 Φ k := by
  rw [dft_apply, dft_apply, Finset.mul_sum]
  refine Fintype.sum_equiv (Equiv.addRight a) _ _ fun j => ?_
  simp only [Equiv.coe_addRight, smul_eq_mul]
  rw [← mul_assoc, ← AddChar.map_add_eq_mul]
  congr 2
  ring

open ZMod in
/-- A real grid has a Hermitian transform: `ρ̂(−k) = conj ρ̂(k)`. -/
theorem dft_neg_eq_conj {N : ℕ} [NeZero N] (ρ : ZMod N → ℝ) (k : ZMod N) :
    𝓕 (fun j => (ρ j : ℂ)) (-k) = (starRingEnd ℂ) (𝓕 (fun j => (ρ j : ℂ)) k) := by
  rw [dft_apply, dft_apply, map_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp only [smul_eq_mul, map_mul, Complex.conj_ofReal, stdAddChar_conj, mul_neg, neg_neg]

/-! ## The transform of the assignment function -/

/-- `P̂_p(x) = (sin(x/2) / (x/2))^p`, `1` at `x = 0` (gromosXX
`fourier_charge_assignment_1d`, `latticesum.h:583-589`; Vol. 2 `wei_product`). -/
noncomputable def sincPow (x : ℝ) (p : ℕ) : ℝ := if x = 0 then 1 else (2 / x * sin (x / 2)) ^ p

theorem sincPow_zero (p : ℕ) : sincPow 0 p = 1 := by simp [sincPow]

/-- `sincPow x p = sinc(x/2)^p` with Mathlib's `Real.sinc` (`sin x / x`, `1` at `0`). -/
theorem sincPow_eq_sinc (x : ℝ) (p : ℕ) : sincPow x p = sinc (x / 2) ^ p := by
  by_cases hx : x = 0
  · simp [sincPow, hx]
  · have : x / 2 ≠ 0 := by positivity
    simp only [sincPow, hx, ite_false, sinc_of_ne_zero this]
    congr 1
    field_simp

/-- The value at `0` is the limit: `sincPow` is continuous there. -/
theorem sincPow_tendsto_zero (p : ℕ) :
    Filter.Tendsto (fun x => sincPow x p) (nhds 0) (nhds 1) := by
  have h : Continuous fun x : ℝ => sinc (x / 2) ^ p := by fun_prop
  simpa [sincPow_eq_sinc] using h.tendsto 0

/-! ## The finite-difference operators -/

/-- The coefficients `c_q(j)`, `j = 1 … q`, of the finite-difference operator of order `q`
(Vol. 2 `tab:fd_coef`; gromosXX `finite_difference_constants`, `latticesum.h:598-633`); zero
elsewhere. -/
def fdCoeff : ℕ → ℕ → ℚ
  | 1, 1 => 1
  | 2, 1 => 4 / 3
  | 2, 2 => -1 / 3
  | 3, 1 => 3 / 2
  | 3, 2 => -3 / 5
  | 3, 3 => 1 / 10
  | 4, 1 => 8 / 5
  | 4, 2 => -4 / 5
  | 4, 3 => 8 / 35
  | 4, 4 => -1 / 35
  | 5, 1 => 5 / 3
  | 5, 2 => -20 / 21
  | 5, 3 => 5 / 14
  | 5, 4 => -5 / 63
  | 5, 5 => 1 / 126
  | _, _ => 0

/-- The finite-difference derivative of order `q` of `f` at `x` (grid spacing 1):
`Σ_{j=1}^{q} c_q(j) (f(x + j) − f(x − j)) / (2j)`. -/
noncomputable def fdDeriv (q : ℕ) (f : ℝ → ℝ) (x : ℝ) : ℝ :=
  ∑ j ∈ Finset.Icc 1 q, (fdCoeff q j : ℝ) * (f (x + j) - f (x - j)) / (2 * j)

/-- The index sets `{1, …, q}` of the five operators as explicit sets. -/
theorem Icc_one_eq :
    Finset.Icc 1 1 = {1} ∧ Finset.Icc 1 2 = {1, 2} ∧ Finset.Icc 1 3 = {1, 2, 3} ∧
      Finset.Icc 1 4 = {1, 2, 3, 4} ∧ Finset.Icc 1 5 = {1, 2, 3, 4, 5} := by
  decide

/-- Order 1 (`c = 1`): exact up to degree 2. -/
theorem fd_exact_1 {n : ℕ} (hn : n ≤ 2) (x : ℝ) :
    fdDeriv 1 (fun y => y ^ n) x = n * x ^ (n - 1) := by
  interval_cases n
  · simp [fdDeriv]
  · simp [fdDeriv, fdCoeff]
  · simp [fdDeriv, fdCoeff]
    ring

/-- Order 2: exact up to degree 4. -/
theorem fd_exact_2 {n : ℕ} (hn : n ≤ 4) (x : ℝ) :
    fdDeriv 2 (fun y => y ^ n) x = n * x ^ (n - 1) := by
  interval_cases n <;> simp [fdDeriv, Icc_one_eq.2.1, Finset.sum_insert, fdCoeff] <;> ring

/-- Order 3: exact up to degree 6. -/
theorem fd_exact_3 {n : ℕ} (hn : n ≤ 6) (x : ℝ) :
    fdDeriv 3 (fun y => y ^ n) x = n * x ^ (n - 1) := by
  interval_cases n <;> simp [fdDeriv, Icc_one_eq.2.2.1, Finset.sum_insert, fdCoeff] <;> ring

/-- Order 4: exact up to degree 8. -/
theorem fd_exact_4 {n : ℕ} (hn : n ≤ 8) (x : ℝ) :
    fdDeriv 4 (fun y => y ^ n) x = n * x ^ (n - 1) := by
  interval_cases n <;> simp [fdDeriv, Icc_one_eq.2.2.2.1, Finset.sum_insert, fdCoeff] <;> ring

/-- Order 5: exact up to degree 10. -/
theorem fd_exact_5 {n : ℕ} (hn : n ≤ 10) (x : ℝ) :
    fdDeriv 5 (fun y => y ^ n) x = n * x ^ (n - 1) := by
  interval_cases n <;> simp [fdDeriv, Icc_one_eq.2.2.2.2, Finset.sum_insert, fdCoeff] <;> ring

/-- **What the table encodes**: the operator of order `q = 1 … 5` differentiates `x ↦ xⁿ`
exactly for every `n ≤ 2q`, so every polynomial of degree `≤ 2q`. -/
theorem fd_exact_on_polynomials {q : ℕ} (hq1 : 1 ≤ q) (hq5 : q ≤ 5) {n : ℕ} (hn : n ≤ 2 * q)
    (x : ℝ) : fdDeriv q (fun y => y ^ n) x = n * x ^ (n - 1) := by
  interval_cases q
  · exact fd_exact_1 hn x
  · exact fd_exact_2 hn x
  · exact fd_exact_3 hn x
  · exact fd_exact_4 hn x
  · exact fd_exact_5 hn x

/-! ## The optimal influence function -/

/-- `|k|²`. -/
noncomputable def norm2 (k : Vec3) : ℝ := k ⬝ᵥ k

/-- The optimal influence function of Vol. 2 (`optim_influ`) at a wave vector with
differentiation operator `D` (`k` itself for `ik`), over the finite set `K` of its alias vectors
`k_m = k + 2π H⁻ᵀ m`, with `P̂(k_m)` the transform of the three-dimensional assignment function:

  `Ĝ = D · Σ_m k_m γ̂(a|k_m|)/|k_m|² P̂(k_m)² / (|D|² (Σ_m P̂(k_m)²)²)`. -/
noncomputable def influence (s : Shape) (a : ℝ) (pHat : Vec3 → ℝ) (D : Vec3) (K : Finset Vec3) :
    ℝ :=
  (D ⬝ᵥ ∑ k ∈ K, (gammaHat s (a * √(norm2 k)) / norm2 k * pHat k ^ 2) • k) /
    (norm2 D * (∑ k ∈ K, pHat k ^ 2) ^ 2)

/-- **A single alias with `ik` differentiation** (NALIAS = 0, NFDORD = 0): `Ĝ(k) =
γ̂(a|k|) / (|k|² P̂(k)²)`, the Ewald weight `γ̂/k²` divided by what the assignment smooths away. -/
theorem influence_single_alias (s : Shape) (a : ℝ) (pHat : Vec3 → ℝ) {k : Vec3}
    (hk : norm2 k ≠ 0) (hP : pHat k ≠ 0) :
    influence s a pHat k {k} = gammaHat s (a * √(norm2 k)) / (norm2 k * pHat k ^ 2) := by
  simp only [influence, Finset.sum_singleton, dotProduct_smul, smul_eq_mul]
  rw [show k ⬝ᵥ k = norm2 k from rfl]
  field_simp

/-! ## The k-space energy -/

open ZMod in
/-- The k-space energy of P3M on a grid of `N` points (Vol. 2 `e_gamma_p3m`; gromosXX
`calculate_potential_and_energy`): `½ V_G ε₀⁻¹ Σ_l Ĝ_l |ρ̂_l|² / N`, with the unnormalised forward
transform `ρ̂` of the grid charge density. -/
noncomputable def p3mEnergy {N : ℕ} [NeZero N] (vCell eps0Inv : ℝ) (ghat : ZMod N → ℝ)
    (ρ : ZMod N → ℝ) : ℝ :=
  1 / 2 * vCell * eps0Inv * ∑ l, ghat l * ‖𝓕 (fun j => (ρ j : ℂ)) l‖ ^ 2 / N

open ZMod in
/-- **The energy is real**: for a real charge density the complex sum `Σ_l Ĝ_l ρ̂_l ρ̂_{−l}` it
is written from is the real `Σ_l Ĝ_l |ρ̂_l|²`. -/
theorem p3m_energy_real {N : ℕ} [NeZero N] (ghat : ZMod N → ℝ) (ρ : ZMod N → ℝ) :
    (∑ l, (ghat l : ℂ) * 𝓕 (fun j => (ρ j : ℂ)) l * 𝓕 (fun j => (ρ j : ℂ)) (-l)) =
      ((∑ l, ghat l * ‖𝓕 (fun j => (ρ j : ℂ)) l‖ ^ 2 : ℝ) : ℂ) := by
  push_cast
  refine Finset.sum_congr rfl fun l _ => ?_
  rw [dft_neg_eq_conj, mul_assoc, Complex.mul_conj, Complex.normSq_eq_norm_sq]
  push_cast
  ring

end NewenTheory.Interactions.P3M

end
