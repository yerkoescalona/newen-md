module

public import Mathlib.Analysis.Calculus.FDeriv.Extend
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Comp
public import Mathlib.Analysis.Complex.Trigonometric
public import Mathlib.Topology.Algebra.InfiniteSum.Basic
public import Mathlib.Data.Matrix.Mul
public import NewenTheory.Interactions.P3M
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring
import Mathlib.Tactic.IntervalCases

/-!
# Smooth particle-mesh Ewald (Essmann et al., J. Chem. Phys. 103, 8577 (1995))

SPME is newen's third lattice-sum method (NONBONDED NLRELE 5, newen only: gromosXX has no
SPME). It is the P3M pipeline of `P3M.lean` with three substitutions: the charges are spread with
the cardinal B-spline `M_n` of order `n` (Essmann eq. 4.1), the influence function is
`γ̂(a|k|)/k² · B(m)` with the Euler exponential-spline factors `B` (eqs. 4.4, 4.7, 4.8), and the
forces are the derivatives of the splines against the convolution of the grid (eq. 4.9).

What is defined and proved here:

* the cardinal B-spline `bspline n` by its recursion (`M_1 = χ_[0,1)`,
  `M_n(u) = u/(n−1) M_{n−1}(u) + (n−u)/(n−1) M_{n−1}(u−1)`): zero outside `[0, n)`
  (`bspline_support`), non-negative (`bspline_nonneg`), a partition of unity for every `n ≥ 1`,
  over a window of `n` grid points (`bspline_partition_of_unity`) and over all of `ℤ`
  (`bspline_tsum_eq_one`), continuous for `n ≥ 2` (`bspline_continuous`), and differentiable
  with `M_n' = M_{n−1}(u) − M_{n−1}(u − 1)`: for `n ≥ 3` at every `u` (`bspline_deriv`), for
  `n = 2` away from the integers, its knots (`bspline_deriv_two`);
* the B-splines are P3M's assignment functions shifted to their centre: `w p x = M_p(x + p/2)`
  for `p = 2 … 5` and every `x`, and for `p = 1` everywhere but at `x = −1/2`, where gromosXX's
  strict `|x| < 1/2` gives weight zero and `M_1(0) = 1` (`bspline_eq_p3m_w`);
* the Euler factor `|b(m)|² = 1/|Σ_{j=0}^{n−2} M_n(j+1) e^{2πimj/K}|²` (`eulerFactor`), positive
  where its denominator does not vanish (`eulerFactor_pos`), `1` at `m = 0`
  (`eulerFactor_zero`);
* the mesh energy `½ Σ_{g,g'} Q(g) θ(g,g') Q(g')` of the spline-assigned charges `Q` with a
  symmetric kernel `θ`, and **`spme_force_is_minus_grad`**: for `n ≥ 3` the force of Essmann
  eq. 4.9, the derivatives of the splines against `θ ⋆ Q`, is minus the derivative of the mesh
  energy with respect to every coordinate of every charge, at every position. This is what SPME
  has over P3M: its force is the exact gradient of its energy;
* the k-space energy of `P3M.lean`, `½ V_G ε₀⁻¹ Σ_l Ĝ_l |ρ̂_l|²/N` on a grid `ZMod N`, is such a
  mesh energy, with the kernel `θ(g,g') = V_G ε₀⁻¹/N Σ_l Ĝ_l Re e^{2πi l(g'−g)/N}`
  (`p3mEnergy_eq_meshEnergy`), which is symmetric (`p3mKernel_symm`).

That SPME approximates the Ewald sum is the method's motivation and is not proved here; the
tests compare the two numerically.
-/

@[expose] public section

namespace NewenTheory.Interactions.SPME

open Real Filter Topology Set

/-! ## The cardinal B-splines -/

/-- The cardinal B-spline `M_n` of order `n` (Essmann eq. 4.1 and the recursion after it);
`M_0 = 0`. -/
noncomputable def bspline : ℕ → ℝ → ℝ
  | 0, _ => 0
  | 1, u => if 0 ≤ u ∧ u < 1 then 1 else 0
  | n + 2, u => u / (n + 1) * bspline (n + 1) u + (n + 2 - u) / (n + 1) * bspline (n + 1) (u - 1)

theorem bspline_one (u : ℝ) : bspline 1 u = if 0 ≤ u ∧ u < 1 then 1 else 0 := rfl

theorem bspline_succ_succ (n : ℕ) (u : ℝ) :
    bspline (n + 2) u =
      u / (n + 1) * bspline (n + 1) u + (n + 2 - u) / (n + 1) * bspline (n + 1) (u - 1) := rfl

/-- `M_1` as two comparisons, the form the case analyses below decide. -/
theorem bspline_one_eq_ite (u : ℝ) :
    bspline 1 u = if u < 0 then 0 else if u < 1 then 1 else 0 := by
  rw [bspline_one]
  by_cases h0 : u < 0
  · simp [h0, not_le.mpr h0]
  · simp [h0, not_lt.mp h0]

/-- **Support**: `M_n` vanishes outside `[0, n)`. -/
theorem bspline_support : ∀ (n : ℕ) {u : ℝ}, (u < 0 ∨ (n : ℝ) ≤ u) → bspline n u = 0
  | 0, _, _ => rfl
  | 1, u, h => by
    rw [bspline_one, ite_eq_right]
    rintro ⟨h1, h2⟩
    rcases h with h | h
    · linarith
    · norm_num at h; linarith
  | n + 2, u, h => by
    have h1 : bspline (n + 1) u = 0 := bspline_support (n + 1) (by
      rcases h with h | h
      · exact Or.inl h
      · right; push_cast at h ⊢; linarith)
    have h2 : bspline (n + 1) (u - 1) = 0 := bspline_support (n + 1) (by
      rcases h with h | h
      · left; linarith
      · right; push_cast at h ⊢; linarith)
    rw [bspline_succ_succ, h1, h2]
    ring

/-- `M_n ≥ 0`. -/
theorem bspline_nonneg : ∀ (n : ℕ) (u : ℝ), 0 ≤ bspline n u
  | 0, _ => le_rfl
  | 1, u => by rw [bspline_one]; split_ifs <;> norm_num
  | n + 2, u => by
    by_cases h0 : u < 0
    · rw [bspline_support (n + 2) (Or.inl h0)]
    by_cases h1 : ((n + 2 : ℕ) : ℝ) ≤ u
    · rw [bspline_support (n + 2) (Or.inr h1)]
    push Not at h0 h1
    push_cast at h1
    rw [bspline_succ_succ]
    have a := bspline_nonneg (n + 1) u
    have b := bspline_nonneg (n + 1) (u - 1)
    have c1 : 0 ≤ u / (n + 1) := by positivity
    have c2 : 0 ≤ (n + 2 - u) / (n + 1) := div_nonneg (by linarith) (by positivity)
    positivity

/-- The `n` values of `M_n` a charge at fractional offset `t ∈ [0, 1)` gives its grid points add
up to one. -/
theorem bspline_window : ∀ (n : ℕ), 1 ≤ n → ∀ {t : ℝ}, 0 ≤ t → t < 1 →
    ∑ j ∈ Finset.range n, bspline n (t + j) = 1
  | 0, h, _, _, _ => absurd h (by norm_num)
  | 1, _, t, h0, h1 => by
    rw [Finset.sum_range_one, Nat.cast_zero, add_zero, bspline_one, ite_eq_left ⟨h0, h1⟩]
  | n + 2, _, t, h0, h1 => by
    have ih := bspline_window (n + 1) (by omega) h0 h1
    simp only [bspline_succ_succ, Finset.sum_add_distrib]
    -- the last term of the first sum and the first of the second vanish (support)
    have hA : ∑ j ∈ Finset.range (n + 2), (t + j) / (n + 1) * bspline (n + 1) (t + j) =
        ∑ j ∈ Finset.range (n + 1), (t + j) / (n + 1) * bspline (n + 1) (t + j) := by
      rw [Finset.sum_range_succ, bspline_support (n + 1) (Or.inr (by push_cast; linarith)),
        mul_zero, add_zero]
    have hB : ∑ j ∈ Finset.range (n + 2), (n + 2 - (t + j)) / (n + 1) * bspline (n + 1) (t + j - 1) =
        ∑ j ∈ Finset.range (n + 1), (n + 1 - (t + j)) / (n + 1) * bspline (n + 1) (t + j) := by
      rw [Finset.sum_range_succ', Nat.cast_zero, add_zero,
        bspline_support (n + 1) (u := t - 1) (Or.inl (by linarith)), mul_zero, add_zero]
      refine Finset.sum_congr rfl fun j _ => ?_
      have : t + ((j + 1 : ℕ) : ℝ) - 1 = t + j := by push_cast; ring
      rw [this]
      push_cast
      ring
    rw [hA, hB, ← Finset.sum_add_distrib, ← ih]
    refine Finset.sum_congr rfl fun j _ => ?_
    field_simp
    ring

/-- **Partition of unity**: the grid points `⌊u⌋ − j`, `j = 0 … n − 1`, carry weights
`M_n(u − (⌊u⌋ − j))` adding up to one, for every order `n ≥ 1` and every `u` (no other point
carries any, `bspline_support`). -/
theorem bspline_partition_of_unity {n : ℕ} (hn : 1 ≤ n) (u : ℝ) :
    ∑ j ∈ Finset.range n, bspline n (u - ((⌊u⌋ - j : ℤ) : ℝ)) = 1 := by
  have h : ∀ j : ℕ, u - ((⌊u⌋ - j : ℤ) : ℝ) = Int.fract u + j := fun j => by
    rw [Int.fract]; push_cast; ring
  simp only [h]
  exact bspline_window n hn (Int.fract_nonneg u) (Int.fract_lt_one u)

/-- The same over all of `ℤ`: `Σ_k M_n(u − k) = 1`. -/
theorem bspline_tsum_eq_one {n : ℕ} (hn : 1 ≤ n) (u : ℝ) :
    ∑' k : ℤ, bspline n (u - k) = 1 := by
  set s := (Finset.range n).image fun j : ℕ => ⌊u⌋ - (j : ℤ)
  rw [tsum_eq_sum (s := s)]
  · rw [Finset.sum_image (fun a _ b _ h => by simpa using h)]
    exact bspline_partition_of_unity hn u
  · intro k hk
    simp only [s, Finset.mem_image, Finset.mem_range, not_exists, not_and] at hk
    apply bspline_support
    by_cases hlt : ⌊u⌋ < k
    · left
      have : (⌊u⌋ : ℝ) + 1 ≤ k := by exact_mod_cast hlt
      linarith [Int.lt_floor_add_one u]
    · right
      push Not at hlt
      have hk' : k ≤ ⌊u⌋ - n := by
        by_contra hc
        push Not at hc
        exact hk (⌊u⌋ - k).toNat (by omega) (by omega)
      have : (k : ℝ) ≤ ⌊u⌋ - n := by exact_mod_cast hk'
      linarith [Int.floor_le u]


/-! ## Continuity and the derivative -/

/-- `M_2` is the hat function `max 0 (1 − |u − 1|)`. -/
theorem bspline_two (u : ℝ) : bspline 2 u = max 0 (1 - |u - 1|) := by
  rw [bspline_succ_succ, bspline_one_eq_ite, bspline_one_eq_ite]
  split_ifs <;> push_cast <;>
  first
    | (exfalso; linarith)
    | (rw [abs_of_neg (by linarith), max_eq_left (by linarith)]; ring)
    | (rw [abs_of_neg (by linarith), max_eq_right (by linarith)]; ring)
    | (rw [abs_of_nonneg (by linarith), max_eq_left (by linarith)]; ring)
    | (rw [abs_of_nonneg (by linarith), max_eq_right (by linarith)]; ring)

/-- `M_n` is continuous for `n ≥ 2`. -/
theorem bspline_continuous {n : ℕ} (hn : 2 ≤ n) : Continuous (bspline n) := by
  induction n, hn using Nat.le_induction with
  | base =>
    have : bspline 2 = fun u => max 0 (1 - |u - 1|) := funext bspline_two
    rw [this]
    fun_prop
  | succ n hn ih =>
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    have : bspline (m + 1 + 1) = fun u =>
        u / (m + 1) * bspline (m + 1) u + (m + 2 - u) / (m + 1) * bspline (m + 1) (u - 1) :=
      funext (bspline_succ_succ m)
    rw [this]
    have hc : Continuous fun u : ℝ => bspline (m + 1) (u - 1) := ih.comp (by fun_prop)
    fun_prop

/-- Near a point that is not an integer `M_1` is constant. -/
theorem bspline_one_hasDerivAt {u : ℝ} (hu : ∀ z : ℤ, u ≠ z) : HasDerivAt (bspline 1) 0 u := by
  have hlt : (⌊u⌋ : ℝ) < u := lt_of_le_of_ne (Int.floor_le u) (fun h => hu ⌊u⌋ h.symm)
  have hgt : u < (⌊u⌋ : ℝ) + 1 := Int.lt_floor_add_one u
  have hev : bspline 1 =ᶠ[𝓝 u] fun _ => bspline 1 u := by
    filter_upwards [Ioo_mem_nhds hlt hgt] with y hy
    have hfy : ⌊y⌋ = ⌊u⌋ := by
      rw [Int.floor_eq_iff]
      exact ⟨hy.1.le, hy.2⟩
    have key : ∀ x : ℝ, bspline 1 x = if ⌊x⌋ = 0 then 1 else 0 := fun x => by
      rw [bspline_one]
      congr 1
      exact propext (Int.floor_eq_zero_iff.symm.trans (by simp [Set.mem_Ico]))
    rw [key y, key u, hfy]
  exact (hasDerivAt_const u (bspline 1 u)).congr_of_eventuallyEq hev

/-- An integer translate of a point that is not an integer is not one either. -/
theorem not_int_sub_one {u : ℝ} (hu : ∀ z : ℤ, u ≠ z) : ∀ z : ℤ, u - 1 ≠ z := fun z h =>
  hu (z + 1) (by push_cast; linarith)

/-- One step of the recursion differentiated: if `M = M_{n+1}` has the derivatives `D₁` at `u`
and `D₂` at `u − 1`, `M_{n+2}` has `(M(u) + u D₁ − M(u − 1) + (n + 2 − u) D₂)/(n + 1)`. -/
theorem bspline_succ_succ_hasDerivAt {n : ℕ} {u D₁ D₂ : ℝ}
    (h₁ : HasDerivAt (bspline (n + 1)) D₁ u) (h₂ : HasDerivAt (bspline (n + 1)) D₂ (u - 1)) :
    HasDerivAt (bspline (n + 2))
      (1 / (n + 1) * bspline (n + 1) u + u / (n + 1) * D₁ +
        (-1 / (n + 1) * bspline (n + 1) (u - 1) + (n + 2 - u) / (n + 1) * D₂)) u := by
  have e : bspline (n + 2) = fun y =>
      y / (n + 1) * bspline (n + 1) y + (n + 2 - y) / (n + 1) * bspline (n + 1) (y - 1) :=
    funext (bspline_succ_succ n)
  rw [e]
  have hy : HasDerivAt (fun y : ℝ => y / (n + 1)) (1 / (n + 1)) u := by
    simpa using (hasDerivAt_id u).div_const ((n : ℝ) + 1)
  have hy' : HasDerivAt (fun y : ℝ => (n + 2 - y) / (n + 1)) (-1 / (n + 1)) u := by
    simpa using ((hasDerivAt_id u).const_sub ((n : ℝ) + 2)).div_const ((n : ℝ) + 1)
  have h₂' : HasDerivAt (fun y => bspline (n + 1) (y - 1)) D₂ u := by
    simpa using h₂.comp_sub_const u 1
  exact (hy.mul h₁).add (hy'.mul h₂')

/-- **The derivative away from the integers**: for `n ≥ 2`,
`M_n'(u) = M_{n−1}(u) − M_{n−1}(u − 1)` at every `u` that is not an integer. -/
theorem bspline_hasDerivAt_of_not_int {n : ℕ} (hn : 2 ≤ n) {u : ℝ} (hu : ∀ z : ℤ, u ≠ z) :
    HasDerivAt (bspline n) (bspline (n - 1) u - bspline (n - 1) (u - 1)) u := by
  induction n, hn using Nat.le_induction generalizing u with
  | base =>
    have h := bspline_succ_succ_hasDerivAt (n := 0) (bspline_one_hasDerivAt hu)
      (bspline_one_hasDerivAt (not_int_sub_one hu))
    convert h using 1
    norm_num
    ring
  | succ n hn ih =>
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
    have h := bspline_succ_succ_hasDerivAt (n := m + 1) (ih hu) (ih (not_int_sub_one hu))
    convert h using 1
    simp only [Nat.add_sub_cancel, show m + 2 - 1 = m + 1 by omega, bspline_succ_succ]
    push_cast
    field_simp
    ring

/-- `f` has the derivative `g x` at `x` if it has `g y` at every `y ≠ x` closer than `1` and `f`,
`g` are continuous at `x` (Mathlib's `hasDerivAt_of_hasDerivAt_of_ne`, on a neighbourhood). -/
theorem hasDerivAt_of_punctured {f g : ℝ → ℝ} {x : ℝ}
    (hd : ∀ y, y ≠ x → |y - x| < 1 → HasDerivAt f (g y) y) (hf : ContinuousAt f x)
    (hg : ContinuousAt g x) : HasDerivAt f (g x) x := by
  have A : HasDerivWithinAt f (g x) (Ici x) x := by
    have hs : Ioo x (x + 1) ∈ 𝓝[>] x := Ioo_mem_nhdsGT (by linarith)
    have hd' : ∀ y ∈ Ioo x (x + 1), HasDerivAt f (g y) y := fun y hy =>
      hd y (ne_of_gt hy.1) (by rw [abs_of_pos (by linarith [hy.1])]; linarith [hy.2])
    apply hasDerivWithinAt_Ici_of_tendsto_deriv
      (fun y hy => (hd' y hy).differentiableAt.differentiableWithinAt)
      hf.continuousWithinAt hs
    have : Tendsto g (𝓝[>] x) (𝓝 (g x)) := tendsto_inf_left hg
    apply this.congr'
    filter_upwards [hs] with y hy
    exact (hd' y hy).deriv.symm
  have B : HasDerivWithinAt f (g x) (Iic x) x := by
    have hs : Ioo (x - 1) x ∈ 𝓝[<] x := Ioo_mem_nhdsLT (by linarith)
    have hd' : ∀ y ∈ Ioo (x - 1) x, HasDerivAt f (g y) y := fun y hy =>
      hd y (ne_of_lt hy.2) (by rw [abs_of_neg (by linarith [hy.2])]; linarith [hy.1])
    apply hasDerivWithinAt_Iic_of_tendsto_deriv
      (fun y hy => (hd' y hy).differentiableAt.differentiableWithinAt)
      hf.continuousWithinAt hs
    have : Tendsto g (𝓝[<] x) (𝓝 (g x)) := tendsto_inf_left hg
    apply this.congr'
    filter_upwards [hs] with y hy
    exact (hd' y hy).deriv.symm
  simpa using B.union A

/-- **The derivative**: for `n ≥ 3`, `M_n'(u) = M_{n−1}(u) − M_{n−1}(u − 1)` at every `u`, the
knots included (`M_n` is `C^{n−2}`). -/
theorem bspline_deriv {n : ℕ} (hn : 3 ≤ n) (u : ℝ) :
    HasDerivAt (bspline n) (bspline (n - 1) u - bspline (n - 1) (u - 1)) u := by
  by_cases hu : ∀ z : ℤ, u ≠ z
  · exact bspline_hasDerivAt_of_not_int (by omega) hu
  push Not at hu
  obtain ⟨z, rfl⟩ := hu
  have hc : Continuous (bspline (n - 1)) := bspline_continuous (by omega)
  have hc' : Continuous fun y : ℝ => bspline (n - 1) y - bspline (n - 1) (y - 1) :=
    hc.sub (hc.comp (continuous_id.sub continuous_const))
  refine hasDerivAt_of_punctured (fun y hy hyd => ?_)
    (bspline_continuous (n := n) (by omega)).continuousAt hc'.continuousAt
  refine bspline_hasDerivAt_of_not_int (by omega) fun k hk => hy ?_
  subst hk
  have h1 : |(k : ℝ) - z| < 1 := hyd
  have h2 : |k - z| < 1 := by exact_mod_cast h1
  have : k = z := by rw [abs_lt] at h2; omega
  rw [this]

/-- For `n = 2`, away from the knots. -/
theorem bspline_deriv_two {u : ℝ} (hu : ∀ z : ℤ, u ≠ z) :
    HasDerivAt (bspline 2) (bspline 1 u - bspline 1 (u - 1)) u :=
  bspline_hasDerivAt_of_not_int le_rfl hu


/-! ## P3M's assignment functions are the B-splines -/

open NewenTheory.Interactions.P3M (w)

/-- Every `u` is below `0`, above `p`, or at or just above a multiple `k/2` of one half,
`k ≤ 2p`. -/
theorem half_cases (p : ℕ) (u : ℝ) :
    u < 0 ∨ (p : ℝ) < u ∨
      ∃ k : ℕ, k ≤ 2 * p ∧ (u = k / 2 ∨ ((k : ℝ) / 2 < u ∧ u < k / 2 + 1 / 2)) := by
  rcases lt_or_ge u 0 with h0 | h0
  · exact Or.inl h0
  rcases lt_or_ge (p : ℝ) u with hp | hp
  · exact Or.inr (Or.inl hp)
  refine Or.inr (Or.inr ⟨⌊2 * u⌋₊, ?_, ?_⟩)
  · have := Nat.floor_le_floor (a := 2 * u) (b := 2 * p) (by linarith)
    rw [show (2 * (p : ℝ)) = ((2 * p : ℕ) : ℝ) by push_cast; ring, Nat.floor_natCast] at this
    exact this
  · have h1 := Nat.floor_le (a := 2 * u) (by linarith)
    have h2 := Nat.lt_floor_add_one (2 * u)
    rcases h1.eq_or_lt with h | h
    · left; linarith
    · right; constructor <;> linarith

/-- Rewrite `|x|` by the sign of `x`, decide every comparison, and close by `ring`. -/
macro "w_region" : tactic => `(tactic| first
  | (exfalso; linarith)
  | (rw [abs_of_nonneg (by linarith)]; decide_ifs; ring)
  | (rw [abs_of_neg (by linarith)]; decide_ifs; ring))

/-- `w p (u − p/2) = M_p(u)` region by region (`half_cases`). -/
macro "w_cases" p:term "at" u:term : tactic => `(tactic| (
  simp only [w, bspline_succ_succ, Nat.reduceAdd, zero_add, bspline_one_eq_ite]
  push_cast
  rcases half_cases $p $u with h | h | ⟨k, hk, h | ⟨h1, h2⟩⟩
  · w_region
  · norm_num at h; w_region
  · interval_cases k <;> norm_num at h <;> subst h <;> w_region
  · interval_cases k <;> norm_num at h1 h2 <;> w_region))

theorem w_two (u : ℝ) : w 2 (u - 2 / 2) = bspline 2 u := by w_cases 2 at u
theorem w_three (u : ℝ) : w 3 (u - 3 / 2) = bspline 3 u := by w_cases 3 at u
theorem w_four (u : ℝ) : w 4 (u - 4 / 2) = bspline 4 u := by w_cases 4 at u
theorem w_five (u : ℝ) : w 5 (u - 5 / 2) = bspline 5 u := by w_cases 5 at u

/-- **P3M's `w_p` is the B-spline of order `p` centred on its grid point**:
`w p x = M_p(x + p/2)` for `p = 2 … 5` and every `x`. -/
theorem bspline_eq_p3m_w {p : ℕ} (hp2 : 2 ≤ p) (hp5 : p ≤ 5) (x : ℝ) :
    w p x = bspline p (x + p / 2) := by
  obtain ⟨u, rfl⟩ : ∃ u, x = u - p / 2 := ⟨x + p / 2, by ring⟩
  rw [sub_add_cancel]
  interval_cases p
  · exact_mod_cast w_two u
  · exact_mod_cast w_three u
  · exact_mod_cast w_four u
  · exact_mod_cast w_five u

/-- For `p = 1` everywhere but at `x = −1/2`. -/
theorem bspline_eq_p3m_w_one {x : ℝ} (hx : x ≠ -1 / 2) : w 1 x = bspline 1 (x + 1 / 2) := by
  simp only [w, bspline_one]
  rcases lt_trichotomy x (-1 / 2) with h | h | h
  · rw [ite_eq_right (fun h' => absurd h' (by rw [abs_of_neg (by linarith)]; linarith)),
      ite_eq_right (fun h' => absurd h'.1 (by linarith))]
  · exact absurd h hx
  · by_cases h2 : x < 1 / 2
    · rw [ite_eq_left (show |x| < 1 / 2 by rw [abs_lt]; constructor <;> linarith),
        ite_eq_left ⟨by linarith, by linarith⟩]
    · rw [ite_eq_right (fun h' => absurd h' (by rw [abs_of_nonneg (by linarith)]; linarith)),
        ite_eq_right (fun h' => absurd h'.2 (by linarith))]

/-- At `x = −1/2` the two differ: gromosXX's strict `|x| < 1/2` gives the point weight zero. -/
theorem p3m_w_one_half : w 1 (-1 / 2) = 0 ∧ bspline 1 (-1 / 2 + 1 / 2) = 1 := by
  norm_num [w, bspline_one]

/-! ## The Euler exponential-spline factors -/

/-- `M_n(0) = 0` for `n ≥ 2`. -/
theorem bspline_at_zero {n : ℕ} (hn : 2 ≤ n) : bspline n 0 = 0 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  rw [bspline_succ_succ, bspline_support (m + 1) (u := 0 - 1) (Or.inl (by norm_num))]
  ring

/-- The denominator `Σ_{j=0}^{n−2} M_n(j+1) e^{2πimj/K}` of Essmann's `b(m)` (eq. 4.4); the
numerator `e^{2πi(n−1)m/K}` has modulus one. -/
noncomputable def eulerDenom (n K : ℕ) (m : ℤ) : ℂ :=
  ∑ j ∈ Finset.range (n - 1),
    (bspline n (j + 1) : ℂ) * Complex.exp (((2 * π * m * j / K : ℝ) : ℂ) * Complex.I)

/-- The Euler exponential-spline factor `|b(m)|²` of an axis of `K` grid points (Essmann
eqs. 4.4, 4.7). -/
noncomputable def eulerFactor (n K : ℕ) (m : ℤ) : ℝ := 1 / Complex.normSq (eulerDenom n K m)

/-- Positive wherever the denominator does not vanish (for every `m` when `n` is even; for odd
`n` it vanishes at `2m = K`, which this file does not prove). -/
theorem eulerFactor_pos {n K : ℕ} {m : ℤ} (h : eulerDenom n K m ≠ 0) : 0 < eulerFactor n K m :=
  one_div_pos.mpr (Complex.normSq_pos.mpr h)

/-- At `m = 0` the denominator is the partition of unity at the grid points: `1`. -/
theorem eulerDenom_zero {n : ℕ} (hn : 2 ≤ n) (K : ℕ) : eulerDenom n K 0 = 1 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  have hw := bspline_window (m + 1) (by omega) (t := 0) le_rfl one_pos
  rw [Finset.sum_range_succ', Nat.cast_zero, zero_add, bspline_at_zero hn, add_zero] at hw
  have hw' : ∑ j ∈ Finset.range m, bspline (m + 1) ((j : ℝ) + 1) = 1 := by
    refine Eq.trans ?_ hw
    refine Finset.sum_congr rfl fun j _ => ?_
    push_cast
    ring_nf
  have h1 : ∀ j : ℕ, (bspline (m + 1) ((j : ℝ) + 1) : ℂ) *
      Complex.exp (((2 * π * ((0 : ℤ) : ℝ) * j / K : ℝ) : ℂ) * Complex.I) =
      ((bspline (m + 1) ((j : ℝ) + 1) : ℝ) : ℂ) := fun j => by simp
  rw [eulerDenom, Nat.add_sub_cancel]
  simp only [h1]
  rw [← Complex.ofReal_sum, hw', Complex.ofReal_one]

theorem eulerFactor_zero {n : ℕ} (hn : 2 ≤ n) (K : ℕ) : eulerFactor n K 0 = 1 := by
  simp [eulerFactor, eulerDenom_zero hn]

/-- `|b(m)|²` by the real and imaginary parts of its denominator, as the code sums them. -/
theorem eulerFactor_eq_re_im (n K : ℕ) (m : ℤ) :
    eulerFactor n K m =
      1 / ((∑ j ∈ Finset.range (n - 1), bspline n (j + 1) * cos (2 * π * m * j / K)) ^ 2 +
        (∑ j ∈ Finset.range (n - 1), bspline n (j + 1) * sin (2 * π * m * j / K)) ^ 2) := by
  rw [eulerFactor, Complex.normSq_apply, eulerDenom, Complex.re_sum, Complex.im_sum]
  simp only [Complex.re_ofReal_mul, Complex.im_ofReal_mul, Complex.exp_ofReal_mul_I_re,
    Complex.exp_ofReal_mul_I_im]
  ring

/-! ## The mesh energy and its gradient -/

section Mesh

variable {ι : Type*} [Fintype ι] {N : ℕ}

/-- The weight a charge at grid coordinates `u` gives the (unwrapped) grid point `c`:
`W(u, c) = Π_β M_n(u_β − c_β)` (Essmann eq. 4.6). -/
noncomputable def weight (n : ℕ) (u : Fin 3 → ℝ) (c : Fin 3 → ℤ) : ℝ := ∏ β, bspline n (u β - c β)

/-- The charge on the grid point `g` of the charges `q` at grid coordinates `u`, each grid point
the images `S g` of itself (the periodic wrap): `Q(g) = Σ_i q_i Σ_{c ∈ S g} W(u_i, c)`. -/
noncomputable def meshCharge (n : ℕ) (S : ι → Finset (Fin 3 → ℤ)) (q : Fin N → ℝ)
    (u : Fin N → Fin 3 → ℝ) (g : ι) : ℝ :=
  ∑ i, q i * ∑ c ∈ S g, weight n (u i) c

/-- The SPME reciprocal energy on a finite grid, `½ Σ_{g,g'} Q(g) θ(g,g') Q(g')`, with the
convolution kernel `θ` (`p3mEnergy_eq_meshEnergy`: the k-space energy
`½ V_cell ε₀⁻¹ Σ_l Ĝ_l |ρ̂_l|²/N` is this with `θ = p3mKernel` and `ρ = Q/V_cell`). -/
noncomputable def meshEnergy (θ : ι → ι → ℝ) (Q : ι → ℝ) : ℝ :=
  1 / 2 * ∑ g, ∑ g', Q g * θ g g' * Q g'

/-- `M_n'(x) = M_{n−1}(x) − M_{n−1}(x − 1)`. -/
noncomputable def bsplineDeriv (n : ℕ) (x : ℝ) : ℝ := bspline (n - 1) x - bspline (n - 1) (x - 1)

/-- The derivative of `W(u + t a, c)` at `t = 0`:
`Σ_β (Π_{γ ≠ β} M_n(u_γ − c_γ)) · M_n'(u_β − c_β) a_β`. -/
noncomputable def weightDeriv (n : ℕ) (u a : Fin 3 → ℝ) (c : Fin 3 → ℤ) : ℝ :=
  ∑ β, (∏ γ ∈ Finset.univ.erase β, bspline n (u γ - c γ)) * (bsplineDeriv n (u β - c β) * a β)

/-- The force of Essmann eq. 4.9 on charge `i` along the direction whose grid coordinates move
by `a` (for the Cartesian component `α`, `a_β = ∂u_β/∂r_α`): minus the charge times the
derivatives of its splines against the convolution `θ ⋆ Q`. -/
noncomputable def spmeForce (n : ℕ) (S : ι → Finset (Fin 3 → ℤ)) (θ : ι → ι → ℝ) (q : Fin N → ℝ)
    (u : Fin N → Fin 3 → ℝ) (a : Fin 3 → ℝ) (i : Fin N) : ℝ :=
  -(q i * ∑ g, (∑ c ∈ S g, weightDeriv n (u i) a c) * ∑ g', θ g g' * meshCharge n S q u g')

/-- Charge `i` moved by `t a` in grid coordinates. -/
def moved (u : Fin N → Fin 3 → ℝ) (i : Fin N) (a : Fin 3 → ℝ) (t : ℝ) : Fin N → Fin 3 → ℝ :=
  fun i' β => u i' β + if i' = i then t * a β else 0

theorem moved_zero (u : Fin N → Fin 3 → ℝ) (i : Fin N) (a : Fin 3 → ℝ) : moved u i a 0 = u := by
  funext i' β
  simp [moved]

theorem moved_self (u : Fin N → Fin 3 → ℝ) (i : Fin N) (a : Fin 3 → ℝ) (t : ℝ) :
    moved u i a t i = fun β => u i β + t * a β := by
  funext β
  show u i β + (if i = i then t * a β else 0) = _
  simp

theorem moved_ne (u : Fin N → Fin 3 → ℝ) {i i' : Fin N} (a : Fin 3 → ℝ) (t : ℝ) (h : i' ≠ i) :
    moved u i a t i' = u i' := by
  funext β
  show u i' β + (if i' = i then t * a β else 0) = _
  simp [h]

/-- The weight of a moving charge: the product rule over the three axes and the chain rule
through `bspline_deriv`. -/
theorem weight_hasDerivAt {n : ℕ} (hn : 3 ≤ n) (u a : Fin 3 → ℝ) (c : Fin 3 → ℤ) :
    HasDerivAt (fun t => weight n (fun β => u β + t * a β) c) (weightDeriv n u a c) 0 := by
  have hf : ∀ β ∈ (Finset.univ : Finset (Fin 3)),
      HasDerivAt (fun t => bspline n (u β + t * a β - c β)) (bsplineDeriv n (u β - c β) * a β) 0 := by
    intro β _
    have hin : HasDerivAt (fun t : ℝ => u β + t * a β - c β) (a β) 0 := by
      simpa using (((hasDerivAt_id (0 : ℝ)).mul_const (a β)).const_add (u β)).sub_const
        ((c β : ℤ) : ℝ)
    have hout := bspline_deriv hn (u β - c β)
    exact HasDerivAt.comp_of_eq (hh₂ := hout) (hh := hin) (hy := by simp)
  have h := HasDerivAt.fun_finsetProd hf
  simp only [zero_mul, add_zero, smul_eq_mul] at h
  exact h

omit [Fintype ι] in
/-- The grid charge as charge `i` moves. -/
theorem meshCharge_hasDerivAt {n : ℕ} (hn : 3 ≤ n) (S : ι → Finset (Fin 3 → ℤ)) (q : Fin N → ℝ)
    (u : Fin N → Fin 3 → ℝ) (i : Fin N) (a : Fin 3 → ℝ) (g : ι) :
    HasDerivAt (fun t => meshCharge n S q (moved u i a t) g)
      (q i * ∑ c ∈ S g, weightDeriv n (u i) a c) 0 := by
  unfold meshCharge
  have h : ∀ i' ∈ (Finset.univ : Finset (Fin N)),
      HasDerivAt (fun t => q i' * ∑ c ∈ S g, weight n (moved u i a t i') c)
        (if i' = i then q i * ∑ c ∈ S g, weightDeriv n (u i) a c else 0) 0 := by
    intro i' _
    by_cases hi : i' = i
    · subst hi
      simp only [↓reduceIte]
      refine HasDerivAt.const_mul _ (HasDerivAt.fun_sum fun c _ => ?_)
      simp only [moved_self]
      exact weight_hasDerivAt hn (u i') a c
    · simp only [hi, ↓reduceIte, moved_ne u a _ hi]
      exact hasDerivAt_const _ _
  have hs := HasDerivAt.fun_sum h
  simpa using hs

/-- **The SPME force is minus the gradient of the mesh energy** (Essmann eq. 4.9): for B-splines
of order `n ≥ 3` and a symmetric kernel `θ`, the derivative of the mesh energy as charge `i`
moves along `a` (in grid coordinates) is minus the force `spmeForce`, at every configuration
(the knots included), for every set of periodic images `S`. -/
theorem spme_force_is_minus_grad_grid {n : ℕ} (hn : 3 ≤ n) (S : ι → Finset (Fin 3 → ℤ))
    {θ : ι → ι → ℝ} (hθ : ∀ g g', θ g g' = θ g' g) (q : Fin N → ℝ) (u : Fin N → Fin 3 → ℝ)
    (i : Fin N) (a : Fin 3 → ℝ) :
    HasDerivAt (fun t => meshEnergy θ (meshCharge n S q (moved u i a t)))
      (-spmeForce n S θ q u a i) 0 := by
  set Q := fun t g => meshCharge n S q (moved u i a t) g
  set Q' := fun g => q i * ∑ c ∈ S g, weightDeriv n (u i) a c
  have hQ : ∀ g, HasDerivAt (fun t => Q t g) (Q' g) 0 := fun g =>
    meshCharge_hasDerivAt hn S q u i a g
  have hE : HasDerivAt (fun t => 1 / 2 * ∑ g, ∑ g', Q t g * θ g g' * Q t g')
      (1 / 2 * ∑ g, ∑ g', (Q' g * θ g g' * Q 0 g' + Q 0 g * θ g g' * Q' g')) 0 := by
    refine HasDerivAt.const_mul _ (HasDerivAt.fun_sum fun g _ => HasDerivAt.fun_sum fun g' _ => ?_)
    have := ((hQ g).mul_const (θ g g')).mul (hQ g')
    convert this using 1
  have h0 : ∀ g, Q 0 g = meshCharge n S q u g := fun g => by simp [Q, moved_zero]
  -- the second double sum is the first with g, g' exchanged (θ symmetric)
  have hswap : ∑ g, ∑ g', meshCharge n S q u g * θ g g' * Q' g' =
      ∑ g, ∑ g', Q' g * θ g g' * meshCharge n S q u g' := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun g _ => Finset.sum_congr rfl fun g' _ => ?_
    rw [hθ g' g]
    ring
  have key : ∑ g, ∑ g', Q' g * θ g g' * meshCharge n S q u g' = -spmeForce n S θ q u a i := by
    rw [spmeForce, neg_neg, Finset.mul_sum]
    refine Finset.sum_congr rfl fun g _ => ?_
    rw [Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun g' _ => ?_
    simp only [Q']
    ring
  unfold meshEnergy
  convert hE using 1
  simp only [h0, Finset.sum_add_distrib, hswap, key]
  ring

/-- Grid coordinates `u_i = A r_i` (`A = K H⁻ᵀ`: `A β α = ∂u_β/∂r_α`). -/
def gridCoords (A : Matrix (Fin 3) (Fin 3) ℝ) (r : Fin N → Fin 3 → ℝ) : Fin N → Fin 3 → ℝ :=
  fun i => A.mulVec (r i)

/-- **`spme_force_is_minus_grad`**: in Cartesian coordinates. Moving charge `i` by `t` along the
axis `α` moves its grid coordinates along the column `A_{·α}`; the derivative of the mesh energy
at `t = 0` is minus Essmann's force `spmeForce` for that column, for every `n ≥ 3`, every
configuration and every symmetric kernel. -/
theorem spme_force_is_minus_grad {n : ℕ} (hn : 3 ≤ n) (S : ι → Finset (Fin 3 → ℤ))
    {θ : ι → ι → ℝ} (hθ : ∀ g g', θ g g' = θ g' g) (q : Fin N → ℝ)
    (A : Matrix (Fin 3) (Fin 3) ℝ) (r : Fin N → Fin 3 → ℝ) (i : Fin N) (α : Fin 3) :
    HasDerivAt
      (fun t : ℝ => meshEnergy θ (meshCharge n S q (gridCoords A
        (fun i' => if i' = i then r i' + t • Pi.single α 1 else r i'))))
      (-spmeForce n S θ q (gridCoords A r) (fun β => A β α) i) 0 := by
  have e : ∀ t : ℝ, gridCoords A (fun i' => if i' = i then r i' + t • Pi.single α 1 else r i') =
      moved (gridCoords A r) i (fun β => A β α) t := fun t => by
    funext i' β
    by_cases hi : i' = i
    · subst hi
      show (Matrix.mulVec A (if i' = i' then r i' + t • Pi.single α 1 else r i')) β =
        (Matrix.mulVec A (r i')) β + (if i' = i' then t * A β α else 0)
      simp only [↓reduceIte]
      rw [Matrix.mulVec_add, Matrix.mulVec_smul, Matrix.mulVec_single_one]
      simp [Matrix.col_apply]
    · show (Matrix.mulVec A (if i' = i then r i' + t • Pi.single α 1 else r i')) β =
        (Matrix.mulVec A (r i')) β + (if i' = i then t * A β α else 0)
      simp only [hi, ↓reduceIte, add_zero]
  simp only [e]
  exact spme_force_is_minus_grad_grid hn S hθ q (gridCoords A r) i (fun β => A β α)

end Mesh

/-! ## P3M's k-space energy is a mesh energy with a symmetric kernel -/

section Transform

open ZMod

/-- `|ρ̂_l|²` of a real grid as a double sum over the grid: `Σ_{g,g'} ρ_g ρ_{g'} Re e^{2πi l(g'−g)/N}`. -/
theorem norm_sq_dft {N : ℕ} [NeZero N] (ρ : ZMod N → ℝ) (l : ZMod N) :
    ‖𝓕 (fun j => (ρ j : ℂ)) l‖ ^ 2 =
      ∑ g, ∑ g', ρ g * ρ g' * (stdAddChar (l * (g' - g)) : ℂ).re := by
  have h := Complex.mul_conj (𝓕 (fun j => (ρ j : ℂ)) l)
  rw [Complex.normSq_eq_norm_sq] at h
  have hre := congrArg Complex.re h
  rw [Complex.ofReal_re] at hre
  rw [← hre, dft_apply, map_sum, Finset.sum_mul_sum, Complex.re_sum]
  refine Finset.sum_congr rfl fun g _ => ?_
  rw [Complex.re_sum]
  refine Finset.sum_congr rfl fun g' _ => ?_
  simp only [smul_eq_mul, map_mul, Complex.conj_ofReal, NewenTheory.Interactions.P3M.stdAddChar_conj]
  have : (stdAddChar (-(g * l)) : ℂ) * (stdAddChar (-(-(g' * l))) : ℂ) =
      stdAddChar (l * (g' - g)) := by
    rw [← AddChar.map_add_eq_mul]
    congr 1
    ring
  rw [show (stdAddChar (-(g * l)) : ℂ) * (ρ g : ℂ) * ((stdAddChar (-(-(g' * l))) : ℂ) * (ρ g' : ℂ)) =
      ((ρ g * ρ g' : ℝ) : ℂ) * ((stdAddChar (-(g * l)) : ℂ) * (stdAddChar (-(-(g' * l))) : ℂ)) by
    push_cast; ring, this, Complex.re_ofReal_mul]

/-- The kernel of P3M's (and SPME's) k-space energy on a grid of `N` points:
`θ(g,g') = V_G ε₀⁻¹/N Σ_l Ĝ_l Re e^{2πi l(g'−g)/N}`. -/
noncomputable def p3mKernel {N : ℕ} [NeZero N] (vCell eps0Inv : ℝ) (ghat : ZMod N → ℝ)
    (g g' : ZMod N) : ℝ :=
  vCell * eps0Inv / N * ∑ l, ghat l * (stdAddChar (l * (g' - g)) : ℂ).re

/-- The kernel is symmetric. -/
theorem p3mKernel_symm {N : ℕ} [NeZero N] (vCell eps0Inv : ℝ) (ghat : ZMod N → ℝ) (g g' : ZMod N) :
    p3mKernel vCell eps0Inv ghat g g' = p3mKernel vCell eps0Inv ghat g' g := by
  unfold p3mKernel
  congr 1
  refine Finset.sum_congr rfl fun l _ => ?_
  congr 1
  rw [show l * (g - g') = -(l * (g' - g)) by ring, ← NewenTheory.Interactions.P3M.stdAddChar_conj,
    Complex.conj_re]

/-- **P3M's k-space energy is the mesh energy of its kernel**: `½ V_G ε₀⁻¹ Σ_l Ĝ_l |ρ̂_l|²/N =
½ Σ_{g,g'} ρ_g θ(g,g') ρ_{g'}`, so `spme_force_is_minus_grad` applies to it with
`θ = p3mKernel` (symmetric by `p3mKernel_symm`) and `ρ = Q/V_cell`. -/
theorem p3mEnergy_eq_meshEnergy {N : ℕ} [NeZero N] (vCell eps0Inv : ℝ) (ghat ρ : ZMod N → ℝ) :
    NewenTheory.Interactions.P3M.p3mEnergy vCell eps0Inv ghat ρ =
      meshEnergy (p3mKernel vCell eps0Inv ghat) ρ := by
  unfold NewenTheory.Interactions.P3M.p3mEnergy meshEnergy p3mKernel
  simp only [norm_sq_dft, Finset.mul_sum, Finset.sum_div]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun g _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun g' _ => ?_
  rw [Finset.sum_mul, Finset.mul_sum]
  refine Finset.sum_congr rfl fun l _ => ?_
  ring

end Transform

end NewenTheory.Interactions.SPME

end
