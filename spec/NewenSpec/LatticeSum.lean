import NewenSpec.Float
import NewenTheory.Interactions.LatticeSum
import NewenTheory.Interactions.LennardJones

/-!
# Model of the lattice-sum kernels of `crates/newen-forces`

Rust (spec keys in `spec/SPEC.toml`), all in `src/nonbonded/lattice_sum.rs`:

* `lattice_sum.shape_switch`: `shape_switch`, the switch function `η(ξ)` and `η'(ξ)` of a charge
  shape (gromosXX `charge_shape_switch`);
* `lattice_sum.shape_fourier`: `shape_fourier`, its Fourier coefficient `γ̂(κ)` and `γ̂'(κ)`
  (gromosXX `charge_shape_fourier`);
* `lattice_sum.lj_ls_pair`: `lj_ls_pair`, one LJ + real-space lattice-sum pair (gromosXX
  `lj_ls_interaction`);
* `lattice_sum.ls_excluded_pair`: `ls_excluded_pair`, the real-space term of an excluded pair
  (gromosXX `ls_excluded_interaction`);
* `lattice_sum.self_constants`: `self_constants`, the constants `A1` and `A3` (gromosXX
  `ls_self_outerloop`);
* the Ewald sum (NLRELE 2), section 4 below: `lattice_sum.reciprocal_box`, `ewald.k_vectors`,
  `ewald.kspace`, `lattice_sum.a2_cubic`, `lattice_sum.a2_shells`, `lattice_sum.self_and_a_terms`,
  `lattice_sum.surface_term` (`lattice_sum.rs`, `ewald.rs`).

Theory: `NewenTheory/Interactions/LatticeSum.lean` (and `LennardJones.lean`).

Each kernel is written once, over any type `α` with the operations the Rust uses (`+ − × ÷`,
negation, `<`, `sqrt`, `exp`, `erfc`, `cos`, `sin` from `NewenSpec/Float.lean`), operation by
operation in the Rust order, which is gromosXX's. The integer literals (up to `16216200`) come
from `[NatCast α]` through Mathlib's numerals; at `Float` that is `Float.ofNat`, exact below
`2⁵³`, which `lake exe spec-vectors` checks. The fractions of the code (`1.0 / 128.0`,
`2.0 / 5.0`, …) are written as divisions, as in the C. `π` is a parameter, given at `Float` by its
bits. Used twice, as in `NewenSpec/Nonbonded.lean`:

1. at `ℝ`, where the `*_real_eq` theorems prove it equal to the theory: for every shape and
   every `ξ`, every `κ ≠ 0`, every `r > 0` and width `a ≠ 0`. The ℝ instance of `Erfc` is the
   theory's `erfc` (Mathlib has none), defined by its integral;
2. at `Float`, where `lake exe spec-vectors` evaluates it and
   `crates/newen-forces/tests/lattice_sum_spec.rs` (and `ewald_spec.rs` for section 4) compares
   the Rust kernels bit for bit.

What is not proved: the rounding error between the two instantiations, and that C `erfc` is
within some ulp of the theory's `erfc`.
-/

namespace NewenSpec.LatticeSum

open NewenTheory.Interactions.LatticeSum NewenTheory.Interactions.LennardJones

/-- Over ℝ the complementary error function is the theory's. -/
noncomputable instance : Erfc ℝ := ⟨NewenTheory.Interactions.LatticeSum.erfc⟩

/-! ## 1. Model (mirrors the Rust functions) -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α]
  [DecidableRel (α := α) (· < ·)] [Sqrt α] [Exp α] [Erfc α] [Cos α] [Sin α]
  [OfNat α 0] [OfNat α 1] [NatCast α]

/-- The truncated polynomial of order `N` at `ξ` with `m = 1 − ξ ≥ 0`: `(η, η')`. -/
def tpSwitchModel (N : TPOrder) (xi m : α) : α × α :=
  match N with
  | .o0 =>
    let m2 := m * m
    (1 / 2 * m2 * (xi + 2), 3 / 2 * (xi * xi - 1))
  | .o1 =>
    let m2 := m * m
    let m3 := m2 * m
    (m3 * (xi + 1), -2 * m2 * (1 + 2 * xi))
  | .o2 =>
    let m2 := m * m
    let m3 := m2 * m
    let m4 := m2 * m2
    (1 / 2 * m4 * (3 * xi + 2), -5 / 2 * m3 * (1 + 3 * xi))
  | .o3 =>
    let m2 := m * m
    let m3 := m2 * m
    let m4 := m2 * m2
    (1 / 4 * m4 * ((4 * xi + 7) * xi + 4), -3 / 4 * m3 * ((8 * xi + 9) * xi + 3))
  | .o4 =>
    let m2 := m * m
    let m4 := m2 * m2
    let m5 := m4 * m
    (1 / 8 * m5 * ((15 * xi + 19) * xi + 8), -21 / 8 * m4 * (1 + xi * (4 + 5 * xi)))
  | .o5 =>
    let m2 := m * m
    let m4 := m2 * m2
    let m5 := m4 * m
    let m6 := m4 * m2
    (m6 * ((3 * xi + 3) * xi + 1), -3 * m5 * (1 + xi * (5 + 8 * xi)))
  | .o6 =>
    let m2 := m * m
    let m5 := m2 * m2 * m
    let m6 := m5 * m
    (1 / 16 * m6 * (((35 * xi + 66) * xi + 51) * xi + 16),
      -9 / 16 * m5 * (5 + xi * (25 + xi * (47 + 35 * xi))))
  | .o7 =>
    let m2 := m * m
    let m6 := m2 * m2 * m2
    let m7 := m6 * m
    (1 / 8 * m7 * (((32 * xi + 49) * xi + 31) * xi + 8),
      -5 / 8 * m6 * (5 + xi * (30 + xi * (69 + 64 * xi))))
  | .o8 =>
    let m3 := m * m * m
    let m7 := m3 * m3 * m
    let m8 := m7 * m
    (1 / 16 * m8 * (((105 * xi + 136) * xi + 73) * xi + 16),
      -55 / 16 * m7 * (1 + 3 * xi) * (1 + xi * (4 + 7 * xi)))
  | .o9 =>
    let m3 := m * m * m
    let m7 := m3 * m3 * m
    let m8 := m7 * m
    (1 / 32 * m8 * ((((160 * xi + 335) * xi + 312) * xi + 151) * xi + 32),
      -15 / 32 * m7 * (7 + xi * (49 + xi * (141 + xi * (203 + 128 * xi)))))
  | .o10 =>
    let m2 := m * m
    let m4 := m2 * m2
    let m8 := m4 * m4
    let m9 := m8 * m
    (1 / 128 * m9 * ((((1155 * xi + 2075) * xi + 1665) * xi + 697) * xi + 128),
      -65 / 128 * m8 * (7 + xi * (56 + 3 * xi * (62 + xi * (104 + 77 * xi)))))

/-- `shape_switch(shape, xi)`: `(η, η')`, with `π` as `pi`. -/
def shapeSwitchModel (pi : α) (s : Shape) (xi : α) : α × α :=
  match s with
  | .gaussian => (Erfc.erfc xi, Exp.exp (-xi * xi) * -2 / Sqrt.sqrt pi)
  | .poly N =>
    let m := 1 - xi
    if m < 0 then (0, 0) else tpSwitchModel α N xi m

/-- The truncated polynomial of order `N` at `κ`: `(γ̂, γ̂')`. -/
def tpFourierModel (N : TPOrder) (kappa : α) : α × α :=
  match N with
  | .o0 =>
    let k2 := kappa * kappa
    let k3 := k2 * kappa
    let s := Sin.sin kappa
    let c := Cos.cos kappa
    let kc := kappa * c
    let k4 := k2 * k2
    (3 * (s - kc) / k3, 3 * (3 * kc + (k2 - 3) * s) / k4)
  | .o1 =>
    let k2 := kappa * kappa
    let k4 := k2 * k2
    let c := Cos.cos kappa
    let ks := kappa * Sin.sin kappa
    let k5 := k4 * kappa
    (12 * (2 - 2 * c - ks) / k4, 12 * (-8 - (k2 - 8) * c + 5 * ks) / k5)
  | .o2 =>
    let k2 := kappa * kappa
    let k5 := k2 * k2 * kappa
    let c := Cos.cos kappa
    let s := Sin.sin kappa
    let k6 := k5 * kappa
    (60 * (2 * kappa + kappa * c - 3 * s) / k5,
      60 * (-8 * kappa - 7 * kappa * c - (k2 - 15) * s) / k6)
  | .o3 =>
    let k2 := kappa * kappa
    let k6 := k2 * k2 * k2
    let c := Cos.cos kappa
    let ks := kappa * Sin.sin kappa
    let k7 := k6 * kappa
    (90 * (8 + (k2 - 8) * c - 5 * ks) / k6,
      90 * (-48 - (9 * k2 - 48) * c - (k2 - 33) * ks) / k7)
  | .o4 =>
    let k2 := kappa * kappa
    let k7 := k2 * k2 * k2 * kappa
    let kc := kappa * Cos.cos kappa
    let s := Sin.sin kappa
    let k8 := k7 * kappa
    (630 * (8 * kappa + 7 * kc + (k2 - 15) * s) / k7,
      630 * (-48 * kappa + (k2 - 57) * kc - 3 * (4 * k2 - 35) * s) / k8)
  | .o5 =>
    let k2 := kappa * kappa
    let k4 := k2 * k2
    let k8 := k4 * k4
    let c := Cos.cos kappa
    let ks := kappa * Sin.sin kappa
    let k9 := k8 * kappa
    (5040 * (4 * (k2 - 6) - (k2 - 24) * c + 9 * ks) / k8,
      5040 * (-24 * (k2 - 8) + 3 * (5 * k2 - 64) * c + (k2 - 87) * ks) / k9)
  | .o6 =>
    let k2 := kappa * kappa
    let k4 := k2 * k2
    let k9 := k4 * k4 * kappa
    let kc := kappa * Cos.cos kappa
    let s := Sin.sin kappa
    let k10 := k9 * kappa
    (7560 * (48 * kappa - (k2 - 57) * kc + 3 * (4 * k2 - 35) * s) / k9,
      7560 * (-384 * kappa + 3 * (6 * k2 - 187) * kc + (k4 - 141 * k2 + 945) * s) / k10)
  | .o7 =>
    let k2 := kappa * kappa
    let k4 := k2 * k2
    let k10 := k4 * k4 * k2
    let c := Cos.cos kappa
    let ks := kappa * Sin.sin kappa
    let k11 := k10 * kappa
    (75600 * (24 * (k2 - 8) - 3 * (5 * k2 - 64) * c - (k2 - 87) * ks) / k10,
      75600 * (-192 * (k2 - 10) - (k4 - 207 * k2 + 1920) * c + (22 * k2 - 975) * ks) / k11)
  | .o8 =>
    let k2 := kappa * kappa
    let k4 := k2 * k2
    let k11 := k4 * k4 * k2 * kappa
    let kc := kappa * Cos.cos kappa
    let s := Sin.sin kappa
    let k12 := k11 * kappa
    (831600 * (8 * kappa * (k2 - 24) + (k2 - 123) * kc - (18 * k2 - 315) * s) / k11,
      831600 * (-64 * kappa * (k2 - 30) - (26 * k2 - 1545) * kc - (k4 - 285 * k2 + 3465) * s) /
        k12)
  | .o9 =>
    let k2 := kappa * kappa
    let k4 := k2 * k2
    let k12 := k4 * k4 * k4
    let c := Cos.cos kappa
    let ks := kappa * Sin.sin kappa
    let k13 := k12 * kappa
    (1247400 * (192 * (k2 - 10) + (k4 - 207 * k2 + 1920) * c - (22 * k2 - 975) * ks) / k12,
      1247400 * (-1920 * (k2 - 12) - 15 * (2 * k4 - 203 * k2 + 1536) * c -
        (k4 - 405 * k2 + 12645) * ks) / k13)
  | .o10 =>
    let k2 := kappa * kappa
    let k4 := k2 * k2
    let k13 := k4 * k4 * k4 * kappa
    let kc := kappa * Cos.cos kappa
    let s := Sin.sin kappa
    let k14 := k13 * kappa
    (16216200 * (64 * kappa * (k2 - 30) + (26 * k2 - 1545) * kc + (k4 - 285 * k2 + 3465) * s) /
        k13,
      16216200 * (-640 * kappa * (k2 - 36) + (k4 - 545 * k2 + 22005) * kc -
        5 * (7 * k4 - 936 * k2 + 9009) * s) / k14)

/-- `shape_fourier(shape, kappa)`: `(γ̂, γ̂')`. -/
def shapeFourierModel (s : Shape) (kappa : α) : α × α :=
  match s with
  | .gaussian =>
    let k2Quarter := kappa * kappa * -(1 / 4)
    (Exp.exp k2Quarter, -(1 / 2) * kappa * Exp.exp k2Quarter)
  | .poly N => tpFourierModel α N kappa

/-- `lj_ls_pair(r2, c6, c12, q_eps, shape, inv_width)`: `(f_mag, e_lj, e_ls)`. -/
def ljLsPairModel (pi r2 c6 c12 qEps : α) (s : Shape) (invWidth : α) : α × α × α :=
  let dist := Sqrt.sqrt r2
  let dist2i := 1 / r2
  let dist6i := dist2i * dist2i * dist2i
  let disti := 1 / dist
  let c12Dist6i := c12 * dist6i
  let aiDist := invWidth * dist
  let eLj := (c12Dist6i - c6) * dist6i
  let sw := shapeSwitchModel α pi s aiDist
  let eLs := qEps * disti * sw.1
  let fLs := dist2i * (eLs - qEps * sw.2 * invWidth)
  ((c12Dist6i + c12Dist6i - c6) * 6 * dist6i * dist2i + fLs, eLj, eLs)

/-- `ls_excluded_pair(r2, q_eps, shape, inv_width)`: `(f_mag, e_ls)`. -/
def lsExcludedPairModel (pi r2 qEps : α) (s : Shape) (invWidth : α) : α × α :=
  let dist := Sqrt.sqrt r2
  let dist2i := 1 / r2
  let disti := 1 / dist
  let aiDist := invWidth * dist
  let sw := shapeSwitchModel α pi s aiDist
  let eLs := qEps * disti * (sw.1 - 1)
  (dist2i * (eLs - qEps * sw.2 * invWidth), eLs)

/-- The table of `ls_self_outerloop`: `(Ã1, Ã3)`. -/
def selfTableModel (pi : α) (s : Shape) : α × α :=
  match s with
  | .gaussian => (1, 2 * 1 / Sqrt.sqrt pi)
  | .poly .o0 => (2 / 5, 3 / 2)
  | .poly .o1 => (4 / 15, 2)
  | .poly .o2 => (4 / 21, 5 / 2)
  | .poly .o3 => (3 / 14, 9 / 4)
  | .poly .o4 => (1 / 6, 21 / 8)
  | .poly .o5 => (2 / 15, 3)
  | .poly .o6 => (8 / 55, 45 / 16)
  | .poly .o7 => (4 / 33, 25 / 8)
  | .poly .o8 => (4 / 39, 55 / 16)
  | .poly .o9 => (10 / 91, 105 / 32)
  | .poly .o10 => (2 / 21, 455 / 128)

/-- `self_constants(shape, width, volume)`: `(A1, A3)`. -/
def selfConstantsModel (pi : α) (s : Shape) (width volume : α) : α × α :=
  let t := selfTableModel α pi s
  (t.1 * (-pi * width * width / volume), t.2 * (-1 / width))

end Model

/-! ## 2. The `Float` instantiation (what `lake exe spec-vectors` runs)

`π` is given by its bits, so that it is exactly Rust's `std::f64::consts::PI`; the Rust test
asserts it. -/

/-- `π` as `f64` bits (the double nearest to π, gromosXX `math::Pi`). -/
def piFloat : Float := Float.ofBits 0x400921fb54442d18

def shapeSwitchFloat (s : Shape) (xi : Float) : Float × Float :=
  shapeSwitchModel Float piFloat s xi

def shapeFourierFloat (s : Shape) (kappa : Float) : Float × Float :=
  shapeFourierModel Float s kappa

def ljLsPairFloat (r2 c6 c12 qEps : Float) (s : Shape) (invWidth : Float) :
    Float × Float × Float :=
  ljLsPairModel Float piFloat r2 c6 c12 qEps s invWidth

def lsExcludedPairFloat (r2 qEps : Float) (s : Shape) (invWidth : Float) : Float × Float :=
  lsExcludedPairModel Float piFloat r2 qEps s invWidth

def selfConstantsFloat (s : Shape) (width volume : Float) : Float × Float :=
  selfConstantsModel Float piFloat s width volume

/-! ## 3. The ℝ instantiation equals the theory -/

/-- The switch function over ℝ, with `pi = π`: the theory's `(η(ξ), η'(ξ))` for every shape and
every `ξ`, the support edge and beyond included. -/
theorem shapeSwitchModel_real_eq (s : Shape) (ξ : ℝ) :
    shapeSwitchModel ℝ Real.pi s ξ = (eta s ξ, etaDeriv s ξ) := by
  cases s with
  | gaussian =>
    simp only [shapeSwitchModel, Erfc.erfc, Exp.exp, Sqrt.sqrt, eta, etaDeriv,
      show -ξ * ξ = -ξ ^ 2 by ring, Prod.mk.injEq, true_and]
    ring
  | poly N =>
    by_cases h : 1 - ξ < 0
    · simp [shapeSwitchModel, eta, etaDeriv, heaviside, h, not_le.mpr h]
    · have hH : heaviside (1 - ξ) = 1 := by simp [heaviside, not_lt.mp h]
      simp only [shapeSwitchModel, h, ite_false, eta, etaDeriv, hH, mul_one]
      cases N <;> simp only [tpSwitchModel, etaPoly, etaPolyDeriv, Prod.mk.injEq] <;> poly_simp <;>
        constructor <;> ring

/-- The Fourier coefficient over ℝ: the theory's `(γ̂(κ), γ̂'(κ))` at every `κ ≠ 0` (at `κ = 0`
the polynomial shapes divide by zero; Vol. 2 defines `γ̂(0) = 1`). -/
theorem shapeFourierModel_real_eq (s : Shape) {κ : ℝ} (hκ : κ ≠ 0) :
    shapeFourierModel ℝ s κ = (gammaHat s κ, gammaHatDeriv s κ) := by
  cases s with
  | gaussian =>
    simp only [shapeFourierModel, Exp.exp, gammaHat, gammaHatDeriv,
      show κ * κ * -(1 / 4) = -κ ^ 2 / 4 by ring]
  | poly N =>
    simp only [shapeFourierModel, gammaHat, gammaHatDeriv, hκ, ite_false]
    cases N <;>
      simp only [tpFourierModel, fcoef, fcoefPrime, fourierForm, TPOrder.toNat, Sin.sin, Cos.cos,
        Prod.mk.injEq] <;>
      poly_simp <;> constructor <;> field_simp <;> ring

/-- One pair at distance `r > 0` with shape width `a ≠ 0`, over ℝ, `inv_width = 1/a`: the
theory's scalar force `ljForce + realSpacePairForce`, `ljEnergy` and `realSpacePair`. -/
theorem ljLsPairModel_real_eq (c6 c12 q : ℝ) (s : Shape) {a r : ℝ} (hr : 0 < r) (ha : a ≠ 0) :
    ljLsPairModel ℝ Real.pi (r ^ 2) c6 c12 q s (1 / a) =
      (ljForce c6 c12 r + realSpacePairForce q a s r, ljEnergy c6 c12 r,
        realSpacePair q a s r) := by
  have hr0 : r ≠ 0 := hr.ne'
  simp only [ljLsPairModel, Sqrt.sqrt, Real.sqrt_sq hr.le, show 1 / a * r = r / a by ring,
    shapeSwitchModel_real_eq, ljForce, ljEnergy, realSpacePairForce, realSpacePair,
    Prod.mk.injEq]
  refine ⟨?_, ?_, ?_⟩
  · field_simp
    ring
  · field_simp
  · field_simp

/-- An excluded pair over ℝ: the theory's `(realSpaceExcludedForce, realSpaceExcluded)`. -/
theorem lsExcludedPairModel_real_eq (q : ℝ) (s : Shape) {a r : ℝ} (hr : 0 < r) (ha : a ≠ 0) :
    lsExcludedPairModel ℝ Real.pi (r ^ 2) q s (1 / a) =
      (realSpaceExcludedForce q a s r, realSpaceExcluded q a s r) := by
  have hr0 : r ≠ 0 := hr.ne'
  simp only [lsExcludedPairModel, Sqrt.sqrt, Real.sqrt_sq hr.le, show 1 / a * r = r / a by ring,
    shapeSwitchModel_real_eq, realSpaceExcludedForce, realSpaceExcluded, Prod.mk.injEq]
  refine ⟨?_, ?_⟩ <;> field_simp

/-- The self-term constants over ℝ: the theory's `(A1, A3)`. -/
theorem selfConstantsModel_real_eq (s : Shape) (a V : ℝ) :
    selfConstantsModel ℝ Real.pi s a V = (a1 s a V, a3 s a) := by
  cases s with
  | gaussian =>
    simp only [selfConstantsModel, selfTableModel, Sqrt.sqrt, a1, a3, a1Tilde, a3Tilde,
      Prod.mk.injEq]
    constructor <;> ring
  | poly N =>
    cases N <;> simp only [selfConstantsModel, selfTableModel, a1, a3, a1Tilde, a3Tilde,
      Prod.mk.injEq] <;> constructor <;> ring

/-! ## 4. The wave vectors, the Ewald sum, `A2`, the self, A and surface terms

Rust (spec keys in `spec/SPEC.toml`): `lattice_sum.reciprocal_box` (`reciprocal_box`),
`lattice_sum.a2_cubic`, `lattice_sum.a2_shells`, `lattice_sum.self_and_a_terms`
(`self_and_a_terms`, `charge_sums`, `self_and_a_virials`), `lattice_sum.surface_term` in
`lattice_sum.rs`; `ewald.k_vectors` and `ewald.kspace` in `ewald.rs`.

A `[f64; 3]` or `glam::DVec3` is a `V3`, a `Sym3` a `Sym`, a `KVector` a `KVec`. The loops of the
Rust are folds from the left in the same order; `-nk..=nk` is `intRange nk`. Two tests are
written with `<` where the Rust writes `==`: `k2 != 0.0` is `0 < k2` (a sum of squares is never
negative, and a NaN fails `k2 ≤ KCUT` first), and `eps_ls == 0.0` is `¬(eps_ls < 0) ∧
¬(0 < eps_ls)` (the same for every double but NaN). The k-space forces: Rust adds the term of
every wave vector to every atom's force in one pass; the model folds each atom's force over the
wave vectors, recomputing `C(k)`, `S(k)` and `sin`, `cos` of that atom: the same operations on
the same values in the same order per atom, hence the same bits. -/

/-- A 3-vector (`[f64; 3]`, `glam::DVec3`). -/
structure V3 (α : Type) where
  x : α
  y : α
  z : α

/-- A symmetric tensor `[xx, xy, xz, yy, yz, zz]` (`lattice_sum::Sym3`). -/
structure Sym (α : Type) where
  xx : α
  xy : α
  xz : α
  yy : α
  yz : α
  zz : α

/-- One wave vector with what is cached for it (`ewald::KVector`). -/
structure KVec (α : Type) where
  k : V3 α
  k2 : α
  k2i : α
  ak : α
  gammaHat : α
  k2iGammaHat : α
  akGammaHatPrime : α

/-- `-n..=n`. -/
def intRange (n : ℕ) : List ℤ := (List.range (2 * n + 1)).map fun i : ℕ => (i : ℤ) - n

/-- gromosXX's half of shell `n` (`a2_shells`): `(n, a, b)` with `|a|, |b| ≤ n`, then `(a, n, b)`
with `|a| ≤ n − 1`, `|b| ≤ n`, then `(a, b, n)` with `|a|, |b| ≤ n − 1`. -/
def faceTriples (n : ℕ) : List (ℤ × ℤ × ℤ) :=
  ((intRange n).flatMap fun a => (intRange n).map fun b => ((n : ℤ), a, b)) ++
    ((intRange (n - 1)).flatMap fun a => (intRange n).map fun b => (a, (n : ℤ), b)) ++
    ((intRange (n - 1)).flatMap fun a => (intRange (n - 1)).map fun b => (a, b, (n : ℤ)))

section Model2

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α] [LE α]
  [DecidableRel (α := α) (· < ·)] [DecidableRel (α := α) (· ≤ ·)] [Sqrt α] [Exp α] [Erfc α]
  [Cos α] [Sin α] [Fabs α] [OfNat α 0] [OfNat α 1] [NatCast α] [IntCast α]

variable {α}

def V3.dot (u v : V3 α) : α := u.x * v.x + u.y * v.y + u.z * v.z

def V3.cross (u v : V3 α) : V3 α :=
  ⟨u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x⟩

def V3.scale (v : V3 α) (s : α) : V3 α := ⟨v.x * s, v.y * s, v.z * s⟩

def V3.zero : V3 α := ⟨0, 0, 0⟩

def Sym.zero : Sym α := ⟨0, 0, 0, 0, 0, 0⟩

/-- `sym_tensor(v1, v2)`. -/
def symTensor (v1 v2 : V3 α) : Sym α :=
  ⟨v1.x * v2.x, v1.x * v2.y, v1.x * v2.z, v1.y * v2.y, v1.y * v2.z, v1.z * v2.z⟩

/-- `m[i] += t[i]`. -/
def Sym.add (m t : Sym α) : Sym α :=
  ⟨m.xx + t.xx, m.xy + t.xy, m.xz + t.xz, m.yy + t.yy, m.yz + t.yz, m.zz + t.zz⟩

/-- `add_to_diagonal(m, d)`. -/
def Sym.addDiag (m : Sym α) (d : α) : Sym α := ⟨m.xx + d, m.xy, m.xz, m.yy + d, m.yz, m.zz + d⟩

/-- `m.map(|x| x * s)`. -/
def Sym.mulRight (m : Sym α) (s : α) : Sym α :=
  ⟨m.xx * s, m.xy * s, m.xz * s, m.yy * s, m.yz * s, m.zz * s⟩

/-- `m.map(|x| s * x)`. -/
def Sym.mulLeft (s : α) (m : Sym α) : Sym α :=
  ⟨s * m.xx, s * m.xy, s * m.xz, s * m.yy, s * m.yz, s * m.zz⟩

variable (α)

/-- `reciprocal_box(edges, rectangular)`: `(V, b₁, b₂, b₃)`. -/
def reciprocalBoxModel (pi : α) (a b c : V3 α) (rectangular : Bool) :
    α × V3 α × V3 α × V3 α :=
  let norm := fun v : V3 α => Sqrt.sqrt (v.x * v.x + v.y * v.y + v.z * v.z)
  let volume := if rectangular then norm a * norm b * norm c else (V3.cross a b).dot c
  let s := 2 * pi / volume
  (volume, (V3.cross b c).scale s, (V3.cross c a).scale s, (V3.cross a b).scale s)

/-- `lattice_vector(recip, l)`. -/
def latticeVectorModel (u v w : V3 α) (l : ℤ × ℤ × ℤ) : V3 α :=
  let x : α := l.1
  let y : α := l.2.1
  let z : α := l.2.2
  ⟨u.x * x + v.x * y + w.x * z, u.y * x + v.y * y + w.y * z, u.z * x + v.z * y + w.z * z⟩

/-- The wave vector of `l` with its cached quantities, if it is kept. -/
def kVectorOf (nkcut : α) (u v w : V3 α) (s : Shape) (width : α) (l : ℤ × ℤ × ℤ) :
    Option (KVec α) :=
  let k := latticeVectorModel α u v w l
  let k2 := k.x * k.x + k.y * k.y + k.z * k.z
  if k2 ≤ nkcut ∧ 0 < k2 then
    let ak := Sqrt.sqrt k2 * width
    let k2i := 1 / k2
    let g := shapeFourierModel α s ak
    some ⟨k, k2, k2i, ak, g.1, k2i * g.1, ak * g.2⟩
  else none

/-- `k_vectors(nk, kcut, recip, shape, width)`. -/
def kVectorsModel (nk : ℕ × ℕ × ℕ) (kcut : α) (u v w : V3 α) (s : Shape) (width : α) :
    List (KVec α) :=
  (intRange nk.1).flatMap fun lx => (intRange nk.2.1).flatMap fun ly =>
    (intRange nk.2.2).filterMap fun lz => kVectorOf α kcut u v w s width (lx, ly, lz)

/-- `C(k)` and `S(k)`: two loops over the atoms in order. -/
def cosSinModel (atoms : List (V3 α × α)) (k : V3 α) : α × α :=
  (atoms.foldl (fun acc a => acc + a.2 * Cos.cos (k.dot a.1)) 0,
    atoms.foldl (fun acc a => acc + a.2 * Sin.sin (k.dot a.1)) 0)

/-- The k-space force on the atom at `r` before its factor `q · eps0_inv / V`. -/
def kForceSumModel (atoms : List (V3 α × α)) (ks : List (KVec α)) (r : V3 α) : V3 α :=
  ks.foldl (fun f kv =>
    let cs := cosSinModel α atoms kv.k
    let rk := kv.k.dot r
    let t := kv.k2iGammaHat * (cs.1 * Sin.sin rk - cs.2 * Cos.cos rk)
    ⟨f.x + kv.k.x * t, f.y + kv.k.y * t, f.z + kv.k.z * t⟩) V3.zero

/-- `ewald_kspace(positions, charges, kvecs, volume, eps0_inv)`:
`(energy, forces, Ã2, Ξ, dÃ2)`. -/
def ewaldKspaceModel (pi : α) (atoms : List (V3 α × α)) (ks : List (KVec α))
    (volume eps0Inv : α) : α × List (V3 α) × α × Sym α × Sym α :=
  let epsVi := eps0Inv / volume
  let energy := ks.foldl (fun e kv =>
    let cs := cosSinModel α atoms kv.k
    e + kv.k2iGammaHat * (cs.1 * cs.1 + cs.2 * cs.2)) 0
  let a2t := ks.foldl (fun acc kv => acc + kv.k2iGammaHat) 0
  let vir := ks.foldl (fun (vs : Sym α × Sym α) kv =>
    let cs := cosSinModel α atoms kv.k
    let ef := cs.1 * cs.1 + cs.2 * cs.2
    let iso := (kv.akGammaHatPrime - 2 * kv.gammaHat) * (kv.k2i * kv.k2i)
    let term := symTensor ⟨iso * kv.k.x, iso * kv.k.y, iso * kv.k.z⟩ kv.k
    (vs.1.add term, (vs.2.add (Sym.mulLeft ef term)).addDiag (ef * kv.k2iGammaHat)))
    (Sym.zero, Sym.zero)
  let forces := atoms.map fun a =>
    let f := kForceSumModel α atoms ks a.1
    let s := a.2 * epsVi
    (⟨f.x * s, f.y * s, f.z * s⟩ : V3 α)
  let fourPiV := 4 * pi / volume
  (energy * (epsVi * (1 / 2)), forces, a2t * fourPiV,
    vir.2.mulRight (-(1 / 4) * epsVi), Sym.mulLeft fourPiV (vir.1.addDiag a2t))

/-- `a2_cubic(wigner, edge, a1, a3)`. -/
def a2CubicModel (wigner edge a1 a3 : α) : α × α :=
  (wigner / edge - a1 - a3, wigner / (3 * edge) - a1)

/-- One half shell of `a2_shells`: the sum `Σ γ̂/k²` and the derivative sum, from `sg`. -/
def a2ShellModel (s : Shape) (width : α) (u v w : V3 α) (n : ℕ) (sg : Sym α) : α × Sym α :=
  (faceTriples n).foldl (fun (acc : α × Sym α) l =>
    let k := latticeVectorModel α u v w l
    let k2 := k.x * k.x + k.y * k.y + k.z * k.z
    let ak := Sqrt.sqrt k2 * width
    let g := shapeFourierModel α s ak
    let iso := 2 * (ak * g.2 - 2 * g.1) / (k2 * k2)
    (acc.1 + g.1 / k2, acc.2.add (symTensor ⟨iso * k.x, iso * k.y, iso * k.z⟩ k))) (0, sg)

/-- The shell loop of `a2_shells` from shell `lmax + 1` on, with `fuel` shells left:
`(A2 before 4π/V, shells, converged, derivative sum)`. -/
def a2ShellsLoop (s : Shape) (width : α) (u v w : V3 α) (tola2 : α) :
    ℕ → ℕ → α → Sym α → α × ℕ × Bool × Sym α
  | 0, lmax, a2, sg => (a2, lmax, false, sg)
  | fuel + 1, lmax, a2, sg =>
    let sh := a2ShellModel α s width u v w (lmax + 1) sg
    let term := sh.1 * 2
    let a2' := a2 + term
    let sg' := sh.2.addDiag term
    if tola2 < Fabs.fabs (term / a2') then
      match fuel with
      | 0 => (a2', lmax + 1, false, sg')
      | f + 1 => a2ShellsLoop s width u v w tola2 (f + 1) (lmax + 1) a2' sg'
    else (a2', lmax + 1, true, sg')

/-- `a2_shells(shape, width, recip, volume, tola2, max_shells)` for `max_shells ≥ 1`:
`(A2, shells, converged, derivative)`. -/
def a2ShellsModel (pi : α) (s : Shape) (width : α) (u v w : V3 α) (volume tola2 : α)
    (maxShells : ℕ) : α × ℕ × Bool × Sym α :=
  let r := a2ShellsLoop α s width u v w tola2 maxShells 0 0 Sym.zero
  let sc := 4 * pi / volume
  (r.1 * sc, r.2.1, r.2.2.1, r.2.2.2.mulRight sc)

/-- `charge_sums(charges)`: `(Σ q, Σ q²)`. -/
def chargeSumsModel (qs : List α) : α × α :=
  (qs.foldl (fun acc q => acc + q) 0, qs.foldl (fun acc q => acc + q * q) 0)

/-- `self_and_a_terms(a1, a2, a3, a2_tilde, sum_q, sum_q2, eps0_inv)`. -/
def selfATermsModel (pi a1 a2 a3 a2t sumQ sumQ2 eps0Inv : α) : α × α :=
  let s2 := sumQ * sumQ
  let eightPi := 8 * pi
  ((a1 + a2 + a3) * sumQ2 * eps0Inv / eightPi, (a1 * s2 - (a1 + a2t) * sumQ2) * eps0Inv / eightPi)

/-- `self_and_a_virials(a1, sum_q, sum_q2, d_a2, d_a2_tilde, four_pi_eps_i)`. -/
def selfAVirialsModel (a1 sumQ sumQ2 : α) (dA2 dA2t : Sym α) (fpe : α) : Sym α × Sym α :=
  let s2 := sumQ * sumQ
  let sv := (Sym.zero : Sym α).addDiag (-(1 / 4) * fpe * a1 * sumQ2)
  let cv := (Sym.zero : Sym α).addDiag (-(1 / 4) * fpe * a1 * (s2 - sumQ2))
  let f := -(1 / 4) * sumQ2 * fpe
  (⟨sv.xx + f * dA2.xx, sv.xy + f * dA2.xy, sv.xz + f * dA2.xz, sv.yy + f * dA2.yy,
      sv.yz + f * dA2.yz, sv.zz + f * dA2.zz⟩,
    ⟨cv.xx - f * dA2t.xx, cv.xy - f * dA2t.xy, cv.xz - f * dA2t.xz, cv.yy - f * dA2t.yy,
      cv.yz - f * dA2t.yz, cv.zz - f * dA2t.zz⟩)

/-- The box dipole `Σ q (r − c)`. -/
def dipoleModel (atoms : List (V3 α × α)) (c : V3 α) : V3 α :=
  atoms.foldl (fun m a => ⟨m.x + a.2 * (a.1.x - c.x), m.y + a.2 * (a.1.y - c.y),
    m.z + a.2 * (a.1.z - c.z)⟩) V3.zero

/-- `surface_term(positions, charges, centre, eps_ls, volume, eps0_inv)`:
`(energy, forces, Ξ)`. -/
def surfaceModel (atoms : List (V3 α × α)) (c : V3 α) (epsLs volume eps0Inv : α) :
    α × List (V3 α) × Sym α :=
  if ¬(epsLs < 0) ∧ ¬(0 < epsLs) then (0, atoms.map fun _ => V3.zero, Sym.zero)
  else
    let m := dipoleModel α atoms c
    let p := eps0Inv / ((epsLs * 2 + 1) * volume)
    let abs2 := m.x * m.x + m.y * m.y + m.z * m.z
    let iso := -(1 / 4) * p * abs2
    let hp := 1 / 2 * p
    let t := symTensor ⟨hp * m.x, hp * m.y, hp * m.z⟩ m
    (1 / 2 * abs2 * p,
      atoms.map fun a => (⟨-p * a.2 * m.x, -p * a.2 * m.y, -p * a.2 * m.z⟩ : V3 α),
      ⟨iso + t.xx, t.xy, t.xz, iso + t.yy, t.yz, iso + t.zz⟩)

end Model2

/-! ## 5. The Ewald kernels at `Float` -/

def reciprocalBoxFloat (a b c : V3 Float) (rectangular : Bool) :
    Float × V3 Float × V3 Float × V3 Float :=
  reciprocalBoxModel Float piFloat a b c rectangular

def kVectorsFloat (nk : ℕ × ℕ × ℕ) (kcut : Float) (u v w : V3 Float) (s : Shape)
    (width : Float) : List (KVec Float) :=
  kVectorsModel Float nk kcut u v w s width

def ewaldKspaceFloat (atoms : List (V3 Float × Float)) (ks : List (KVec Float))
    (volume eps0Inv : Float) : Float × List (V3 Float) × Float × Sym Float × Sym Float :=
  ewaldKspaceModel Float piFloat atoms ks volume eps0Inv

def a2CubicFloat (wigner edge a1 a3 : Float) : Float × Float :=
  a2CubicModel Float wigner edge a1 a3

def a2ShellsFloat (s : Shape) (width : Float) (u v w : V3 Float) (volume tola2 : Float)
    (maxShells : ℕ) : Float × ℕ × Bool × Sym Float :=
  a2ShellsModel Float piFloat s width u v w volume tola2 maxShells

def chargeSumsFloat (qs : List Float) : Float × Float := chargeSumsModel Float qs

def selfATermsFloat (a1 a2 a3 a2t sumQ sumQ2 eps0Inv : Float) : Float × Float :=
  selfATermsModel Float piFloat a1 a2 a3 a2t sumQ sumQ2 eps0Inv

def selfAVirialsFloat (a1 sumQ sumQ2 : Float) (dA2 dA2t : Sym Float) (fpe : Float) :
    Sym Float × Sym Float :=
  selfAVirialsModel Float a1 sumQ sumQ2 dA2 dA2t fpe

def surfaceFloat (atoms : List (V3 Float × Float)) (c : V3 Float) (epsLs volume eps0Inv : Float) :
    Float × List (V3 Float) × Sym Float :=
  surfaceModel Float atoms c epsLs volume eps0Inv


/-! ## 6. The Ewald kernels at ℝ: the model is the theory -/

/-- A model vector over ℝ as the theory's. -/
def V3.toVec (v : V3 ℝ) : Vec3 := ![v.x, v.y, v.z]

theorem V3.dot_eq (u v : V3 ℝ) : u.dot v = u.toVec ⬝ᵥ v.toVec := by
  simp [V3.dot, V3.toVec, dotProduct, Fin.sum_univ_three]

theorem V3.cross_toVec (u v : V3 ℝ) : (u.cross v).toVec = crossProduct u.toVec v.toVec := by
  ext i
  fin_cases i <;> simp [V3.cross, V3.toVec, cross_apply]

theorem V3.scale_toVec (v : V3 ℝ) (s : ℝ) : (v.scale s).toVec = s • v.toVec := by
  ext i
  fin_cases i <;> simp [V3.scale, V3.toVec, mul_comm]

/-- The box volume and the reciprocal basis over ℝ (triclinic flag): the theory's `boxVolume` and
`reciprocal`. -/
theorem reciprocalBoxModel_real_eq (a b c : V3 ℝ) :
    (reciprocalBoxModel ℝ Real.pi a b c false).1 = boxVolume a.toVec b.toVec c.toVec ∧
      (reciprocalBoxModel ℝ Real.pi a b c false).2.1.toVec =
        reciprocal a.toVec b.toVec c.toVec 0 ∧
      (reciprocalBoxModel ℝ Real.pi a b c false).2.2.1.toVec =
        reciprocal a.toVec b.toVec c.toVec 1 ∧
      (reciprocalBoxModel ℝ Real.pi a b c false).2.2.2.toVec =
        reciprocal a.toVec b.toVec c.toVec 2 := by
  simp only [reciprocalBoxModel, Bool.false_eq_true, ite_false, V3.scale_toVec, V3.cross_toVec,
    V3.dot_eq, reciprocal, boxVolume, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, and_self]

/-- In a rectangular box `diag(Lx, Ly, Lz)` with positive edges, the rectangular volume
`|a| |b| |c|` is `(a × b) · c`. -/
theorem reciprocalBoxModel_rectangular_volume {lx ly lz : ℝ} (hx : 0 < lx) (hy : 0 < ly)
    (hz : 0 < lz) :
    (reciprocalBoxModel ℝ Real.pi ⟨lx, 0, 0⟩ ⟨0, ly, 0⟩ ⟨0, 0, lz⟩ true).1 =
      boxVolume (V3.toVec ⟨lx, 0, 0⟩) (V3.toVec ⟨0, ly, 0⟩) (V3.toVec ⟨0, 0, lz⟩) := by
  simp only [reciprocalBoxModel, ite_true, boxVolume, V3.toVec, cross_apply, dotProduct,
    Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons]
  simp only [Sqrt.sqrt, mul_zero, zero_mul, add_zero, zero_add, sub_zero,
    Real.sqrt_mul_self hx.le, Real.sqrt_mul_self hy.le, Real.sqrt_mul_self hz.le]

/-- An integer triple as the theory's `Fin 3 → ℤ`. -/
def toFin (l : ℤ × ℤ × ℤ) : Fin 3 → ℤ := ![l.1, l.2.1, l.2.2]

theorem toFin_injective : Function.Injective toFin := by
  rintro ⟨a, b, c⟩ ⟨a', b', c'⟩ h
  have h0 := congrFun h 0
  have h1 := congrFun h 1
  have h2 := congrFun h 2
  simp only [toFin, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons] at h0 h1 h2
  rw [h0, h1, h2]

/-- The wave vector of `l` over ℝ is the theory's `latticeVector`. -/
theorem latticeVectorModel_real_eq (u v w : V3 ℝ) (l : ℤ × ℤ × ℤ) :
    (latticeVectorModel ℝ u v w l).toVec = latticeVector ![u.toVec, v.toVec, w.toVec] (toFin l) := by
  ext i
  fin_cases i <;>
    simp [latticeVectorModel, latticeVector, toFin, V3.toVec, Fin.sum_univ_three] <;> ring

theorem mem_intRange {n : ℕ} {z : ℤ} : z ∈ intRange n ↔ -(n : ℤ) ≤ z ∧ z ≤ n := by
  simp only [intRange, List.mem_map, List.mem_range]
  constructor
  · rintro ⟨i, hi, rfl⟩
    omega
  · intro h
    exact ⟨(z + n).toNat, by omega, by omega⟩

theorem intRange_nodup (n : ℕ) : (intRange n).Nodup := by
  refine List.Nodup.map (fun i j h => ?_) List.nodup_range
  omega

/-- A sum loop over ℝ from `0` is the list sum. -/
theorem foldl_add_eq_sum {β : Type} (f : β → ℝ) (l : List β) (init : ℝ) :
    l.foldl (fun acc b => acc + f b) init = init + (l.map f).sum := by
  induction l generalizing init with
  | nil => simp
  | cons b l ih => simp [ih]; ring

theorem foldl_v3_add {β : Type} (l : List β) (gx gy gz : β → ℝ) (init : V3 ℝ) :
    l.foldl (fun f b => (⟨f.x + gx b, f.y + gy b, f.z + gz b⟩ : V3 ℝ)) init =
      ⟨init.x + (l.map gx).sum, init.y + (l.map gy).sum, init.z + (l.map gz).sum⟩ := by
  induction l generalizing init with
  | nil => simp
  | cons b l ih =>
    simp only [List.foldl_cons, ih, List.map_cons, List.sum_cons, V3.mk.injEq]
    refine ⟨?_, ?_, ?_⟩ <;> ring

/-- What `k_vectors` caches for a wave vector over ℝ: `k2 = |k|² ≠ 0`, `k2i = 1/k2`,
`ak = √k2 · a`, `γ̂(ak)`, `k2i · γ̂`, `ak · γ̂'(ak)`. -/
structure KVec.WF (s : Shape) (width : ℝ) (kv : KVec ℝ) : Prop where
  k2 : kv.k2 = kv.k.toVec ⬝ᵥ kv.k.toVec
  k2_ne : kv.k2 ≠ 0
  k2i : kv.k2i = 1 / kv.k2
  ak : kv.ak = √kv.k2 * width
  gammaHat : kv.gammaHat = NewenTheory.Interactions.LatticeSum.gammaHat s kv.ak
  k2iGammaHat : kv.k2iGammaHat = kv.k2i * kv.gammaHat
  akGammaHatPrime : kv.akGammaHatPrime = kv.ak * gammaHatDeriv s kv.ak

theorem k2_model_eq (k : V3 ℝ) : k.x * k.x + k.y * k.y + k.z * k.z = k.toVec ⬝ᵥ k.toVec := by
  simp [V3.toVec, dotProduct, Fin.sum_univ_three]

theorem toVec_dot_self_nonneg (k : V3 ℝ) : 0 ≤ k.toVec ⬝ᵥ k.toVec := by
  rw [← k2_model_eq]
  exact add_nonneg (add_nonneg (mul_self_nonneg _) (mul_self_nonneg _)) (mul_self_nonneg _)

/-- `kVectorOf` over ℝ keeps the wave vector of `l` exactly when the theory's `kSet` filter does. -/
theorem kVectorOf_real_map (kcut : ℝ) (u v w : V3 ℝ) (s : Shape) (width : ℝ) (l : ℤ × ℤ × ℤ) :
    (kVectorOf ℝ kcut u v w s width l).map (fun kv => kv.k.toVec) =
      if latticeVector ![u.toVec, v.toVec, w.toVec] (toFin l) ⬝ᵥ
            latticeVector ![u.toVec, v.toVec, w.toVec] (toFin l) ≤ kcut ∧
          latticeVector ![u.toVec, v.toVec, w.toVec] (toFin l) ⬝ᵥ
            latticeVector ![u.toVec, v.toVec, w.toVec] (toFin l) ≠ 0
        then some (latticeVector ![u.toVec, v.toVec, w.toVec] (toFin l)) else none := by
  have hpos := toVec_dot_self_nonneg (latticeVectorModel ℝ u v w l)
  simp only [kVectorOf, k2_model_eq, latticeVectorModel_real_eq] at hpos ⊢
  split_ifs with h1 h2 h2
  · simp [latticeVectorModel_real_eq]
  · exact absurd ⟨h1.1, h1.2.ne'⟩ h2
  · exact absurd ⟨h2.1, lt_of_le_of_ne hpos (Ne.symm h2.2)⟩ h1
  · rfl

/-- **The model's wave vectors are the theory's `kSet`**: over ℝ, the set of the `k` that
`kVectorsModel` returns is `kSet` of the box of triples `|l_j| ≤ NK_j`, KCUT and the
reciprocal basis. -/
theorem kVectorsModel_real_toFinset (nk : ℕ × ℕ × ℕ) (kcut : ℝ) (u v w : V3 ℝ) (s : Shape)
    (width : ℝ) :
    ((kVectorsModel ℝ nk kcut u v w s width).map fun kv => kv.k.toVec).toFinset =
      kSet ![nk.1, nk.2.1, nk.2.2] kcut ![u.toVec, v.toVec, w.toVec] := by
  classical
  ext k
  have hmem : ∀ l : ℤ × ℤ × ℤ, toFin l ∈ lBox ![nk.1, nk.2.1, nk.2.2] ↔
      l.1 ∈ intRange nk.1 ∧ l.2.1 ∈ intRange nk.2.1 ∧ l.2.2 ∈ intRange nk.2.2 := by
    rintro ⟨a, b, c⟩
    simp only [mem_lBox, mem_intRange]
    constructor
    · intro h
      exact ⟨by simpa [toFin] using h 0, by simpa [toFin] using h 1, by simpa [toFin] using h 2⟩
    · rintro ⟨ha, hb, hc⟩ j
      fin_cases j <;> simpa [toFin]
  simp only [List.mem_toFinset, List.mem_map, kVectorsModel, List.mem_flatMap,
    List.mem_filterMap, kSet, Finset.mem_image, Finset.mem_filter]
  constructor
  · rintro ⟨kv, ⟨lx, hx, ly, hy, lz, hz, hkv⟩, rfl⟩
    have := kVectorOf_real_map kcut u v w s width (lx, ly, lz)
    rw [hkv] at this
    split_ifs at this with h
    · refine ⟨toFin (lx, ly, lz), ⟨(hmem _).mpr ⟨hx, hy, hz⟩, h⟩, ?_⟩
      simpa using this.symm
    · simp at this
  · rintro ⟨l, ⟨hl, h⟩, rfl⟩
    have hl' : toFin (l 0, l 1, l 2) = l := by
      funext j
      fin_cases j <;> rfl
    rw [← hl'] at hl h ⊢
    obtain ⟨hx, hy, hz⟩ := (hmem _).mp hl
    have := kVectorOf_real_map kcut u v w s width (l 0, l 1, l 2)
    simp only [ne_eq, h.1, h.2, not_false_eq_true, and_self, ite_true] at this
    obtain ⟨kv, hkv, hk⟩ := Option.map_eq_some_iff.mp this
    exact ⟨kv, ⟨l 0, hx, l 1, hy, l 2, hz, hkv⟩, hk⟩

/-- **Every cached wave vector is consistent** (`KVec.WF`), for a width `a ≠ 0`. -/
theorem kVectorsModel_real_wf (nk : ℕ × ℕ × ℕ) (kcut : ℝ) (u v w : V3 ℝ) (s : Shape)
    {width : ℝ} (hw : width ≠ 0) :
    ∀ kv ∈ kVectorsModel ℝ nk kcut u v w s width, kv.WF s width := by
  intro kv hkv
  simp only [kVectorsModel, List.mem_flatMap, List.mem_filterMap] at hkv
  obtain ⟨lx, -, ly, -, lz, -, hkv⟩ := hkv
  simp only [kVectorOf] at hkv
  split_ifs at hkv with h
  obtain rfl := (Option.some.inj hkv).symm
  have hk2 : (latticeVectorModel ℝ u v w (lx, ly, lz)).x * (latticeVectorModel ℝ u v w (lx, ly, lz)).x +
      (latticeVectorModel ℝ u v w (lx, ly, lz)).y * (latticeVectorModel ℝ u v w (lx, ly, lz)).y +
      (latticeVectorModel ℝ u v w (lx, ly, lz)).z * (latticeVectorModel ℝ u v w (lx, ly, lz)).z ≠ 0 :=
    h.2.ne'
  have hak : √((latticeVectorModel ℝ u v w (lx, ly, lz)).x * (latticeVectorModel ℝ u v w (lx, ly, lz)).x +
      (latticeVectorModel ℝ u v w (lx, ly, lz)).y * (latticeVectorModel ℝ u v w (lx, ly, lz)).y +
      (latticeVectorModel ℝ u v w (lx, ly, lz)).z * (latticeVectorModel ℝ u v w (lx, ly, lz)).z) *
        width ≠ 0 :=
    mul_ne_zero (Real.sqrt_pos.mpr h.2).ne' hw
  exact ⟨k2_model_eq _, hk2, rfl, rfl, by simp [Sqrt.sqrt, shapeFourierModel_real_eq s hak], rfl,
    by simp [Sqrt.sqrt, shapeFourierModel_real_eq s hak]⟩

section KSpaceReal

variable {n : ℕ} (pos : Fin n → V3 ℝ) (q : Fin n → ℝ)

/-- The atoms of a kernel call over ℝ. -/
def atomsOf : List (V3 ℝ × ℝ) := List.ofFn fun i => (pos i, q i)

theorem cosSinModel_real_eq (k : V3 ℝ) :
    cosSinModel ℝ (atomsOf pos q) k =
      (cosSum q (fun i => (pos i).toVec) k.toVec, sinSum q (fun i => (pos i).toVec) k.toVec) := by
  simp only [cosSinModel, atomsOf, foldl_add_eq_sum, zero_add, List.map_ofFn, List.sum_ofFn,
    cosSum, sinSum, Function.comp_def, V3.dot_eq, Cos.cos, Sin.sin]

/-- **The k-space energy over ℝ** is the theory's sum over the cached wave vectors, with
`ε₀ = 1/eps0_inv`; over a list without repeated `k`, it is `kSpace` over their set. -/
theorem ewaldKspaceModel_real_energy (s : Shape) (width : ℝ) (ks : List (KVec ℝ))
    (hwf : ∀ kv ∈ ks, kv.WF s width) (V e : ℝ) :
    (ewaldKspaceModel ℝ Real.pi (atomsOf pos q) ks V e).1 =
      (2 * e⁻¹ * V)⁻¹ * (ks.map fun kv => gammaHat s (width * √(kv.k.toVec ⬝ᵥ kv.k.toVec)) /
        (kv.k.toVec ⬝ᵥ kv.k.toVec) * structureFactorSq q (fun i => (pos i).toVec) kv.k.toVec).sum := by
  simp only [ewaldKspaceModel]
  rw [foldl_add_eq_sum, zero_add]
  have : ∀ kv ∈ ks, kv.k2iGammaHat * ((cosSinModel ℝ (atomsOf pos q) kv.k).1 *
      (cosSinModel ℝ (atomsOf pos q) kv.k).1 + (cosSinModel ℝ (atomsOf pos q) kv.k).2 *
        (cosSinModel ℝ (atomsOf pos q) kv.k).2) =
      gammaHat s (width * √(kv.k.toVec ⬝ᵥ kv.k.toVec)) / (kv.k.toVec ⬝ᵥ kv.k.toVec) *
        structureFactorSq q (fun i => (pos i).toVec) kv.k.toVec := by
    intro kv hkv
    obtain ⟨h2, -, hi, ha, hg, hig, -⟩ := hwf kv hkv
    rw [cosSinModel_real_eq, hig, hg, ha, hi, h2, mul_comm (√_) width]
    simp only [structureFactorSq]
    ring
  rw [List.map_congr_left this]
  simp only [mul_inv_rev, inv_inv]
  ring

theorem ewaldKspaceModel_real_eq_kSpace (s : Shape) (width : ℝ) (ks : List (KVec ℝ))
    (hwf : ∀ kv ∈ ks, kv.WF s width) (hnd : (ks.map fun kv => kv.k.toVec).Nodup) (V e : ℝ) :
    (ewaldKspaceModel ℝ Real.pi (atomsOf pos q) ks V e).1 =
      kSpace e⁻¹ V width s (ks.map fun kv => kv.k.toVec).toFinset q fun i => (pos i).toVec := by
  classical
  rw [ewaldKspaceModel_real_energy pos q s width ks hwf, kSpace, List.sum_toFinset _ hnd,
    List.map_map]
  rfl

/-- **The k-space forces over ℝ** are the theory's `kSpaceForce` over the cached wave vectors. -/
theorem ewaldKspaceModel_real_forces (s : Shape) (width : ℝ) (ks : List (KVec ℝ))
    (hwf : ∀ kv ∈ ks, kv.WF s width) (hnd : (ks.map fun kv => kv.k.toVec).Nodup) (V e : ℝ) :
    (ewaldKspaceModel ℝ Real.pi (atomsOf pos q) ks V e).2.1.map V3.toVec =
      List.ofFn fun i => kSpaceForce e⁻¹ V width s (ks.map fun kv => kv.k.toVec).toFinset q
        (fun i => (pos i).toVec) i := by
  classical
  simp only [ewaldKspaceModel, atomsOf, List.map_ofFn, Function.comp_def]
  congr 1
  funext i
  simp only [kForceSumModel, V3.zero]
  rw [foldl_v3_add]
  ext α
  have hsum : ∀ f : Vec3 → ℝ, ∑ k ∈ (ks.map fun kv => kv.k.toVec).toFinset, f k =
      (ks.map fun kv => f kv.k.toVec).sum := fun f => by
    rw [List.sum_toFinset _ hnd, List.map_map]
    rfl
  have hterm : ∀ (j : Fin 3), ∀ kv ∈ ks, kv.k.toVec j * (kv.k2iGammaHat *
      ((cosSinModel ℝ (List.ofFn fun i => (pos i, q i)) kv.k).1 * Sin.sin (kv.k.dot (pos i)) -
        (cosSinModel ℝ (List.ofFn fun i => (pos i, q i)) kv.k).2 * Cos.cos (kv.k.dot (pos i)))) =
      gammaHat s (width * √(kv.k.toVec ⬝ᵥ kv.k.toVec)) / (kv.k.toVec ⬝ᵥ kv.k.toVec) *
        kv.k.toVec j * (cosSum q (fun i => (pos i).toVec) kv.k.toVec *
          Real.sin (kv.k.toVec ⬝ᵥ (pos i).toVec) - sinSum q (fun i => (pos i).toVec) kv.k.toVec *
            Real.cos (kv.k.toVec ⬝ᵥ (pos i).toVec)) := by
    intro j kv hkv
    obtain ⟨h2, -, hi, ha, hg, hig, -⟩ := hwf kv hkv
    have hcs := cosSinModel_real_eq pos q kv.k
    simp only [atomsOf] at hcs
    rw [hcs, hig, hg, ha, hi, h2, mul_comm (√_) width, V3.dot_eq]
    simp only [Sin.sin, Cos.cos]
    ring
  simp only [kSpaceForce]
  rw [hsum]
  fin_cases α
  · have h := List.map_congr_left (hterm 0)
    simp only [V3.toVec, Fin.zero_eta, Fin.isValue, Matrix.cons_val_zero, zero_add] at h ⊢
    rw [h]
    simp only [mul_inv_rev, inv_inv]
    ring
  · have h := List.map_congr_left (hterm 1)
    simp only [V3.toVec, Fin.mk_one, Fin.isValue, Matrix.cons_val_one, Matrix.cons_val_zero,
      zero_add] at h ⊢
    rw [h]
    simp only [mul_inv_rev, inv_inv]
    ring
  · have h := List.map_congr_left (hterm 2)
    simp only [V3.toVec, Fin.reduceFinMk, Fin.isValue, Matrix.cons_val_two, Matrix.tail_cons,
      Matrix.head_cons, zero_add] at h ⊢
    rw [h]
    simp only [mul_inv_rev, inv_inv]
    ring

/-- `Ã2` over ℝ: `4π/V Σ γ̂/k²` over the cached wave vectors. -/
theorem ewaldKspaceModel_real_a2Tilde (s : Shape) (width : ℝ) (ks : List (KVec ℝ))
    (hwf : ∀ kv ∈ ks, kv.WF s width) (V e : ℝ) :
    (ewaldKspaceModel ℝ Real.pi (atomsOf pos q) ks V e).2.2.1 =
      4 * Real.pi / V * (ks.map fun kv => gammaHat s (width * √(kv.k.toVec ⬝ᵥ kv.k.toVec)) /
        (kv.k.toVec ⬝ᵥ kv.k.toVec)).sum := by
  simp only [ewaldKspaceModel]
  rw [foldl_add_eq_sum, zero_add]
  have : ∀ kv ∈ ks, kv.k2iGammaHat =
      gammaHat s (width * √(kv.k.toVec ⬝ᵥ kv.k.toVec)) / (kv.k.toVec ⬝ᵥ kv.k.toVec) := by
    intro kv hkv
    obtain ⟨h2, -, hi, ha, hg, hig, -⟩ := hwf kv hkv
    rw [hig, hg, ha, hi, h2, mul_comm (√_) width]
    ring
  rw [List.map_congr_left this]
  ring

/-- `(Σ q, Σ q²)` over ℝ. -/
theorem chargeSumsModel_real_eq :
    chargeSumsModel ℝ (List.ofFn q) = (∑ i, q i, ∑ i, q i ^ 2) := by
  have h1 := foldl_add_eq_sum (fun x : ℝ => x) (List.ofFn q) 0
  simp only [chargeSumsModel, h1, foldl_add_eq_sum, zero_add, List.map_ofFn, List.sum_ofFn,
    Function.comp_def, sq]

/-- **The self and A terms over ℝ** are the theory's `selfTerm` and `aTerm`, with
`ε₀ = 1/eps0_inv`, `S = Σ q`, `S̃² = Σ q²`. -/
theorem selfATermsModel_real_eq (A1 A2 A3 A2t e : ℝ) :
    selfATermsModel ℝ Real.pi A1 A2 A3 A2t (∑ i, q i) (∑ i, q i ^ 2) e =
      (selfTerm e⁻¹ A1 A2 A3 q, aTerm e⁻¹ A1 A2t q) := by
  simp only [selfATermsModel, selfTerm, aTerm, Prod.mk.injEq, div_eq_mul_inv, mul_inv_rev,
    inv_inv]
  constructor <;> ring

/-- The dipole over ℝ is the theory's `dipole`. -/
theorem dipoleModel_real_eq (c : V3 ℝ) :
    (dipoleModel ℝ (atomsOf pos q) c).toVec = dipole q (fun i => (pos i).toVec) c.toVec := by
  simp only [dipoleModel, atomsOf, V3.zero]
  rw [foldl_v3_add]
  ext j
  fin_cases j <;>
    simp [V3.toVec, dipole, List.map_ofFn, List.sum_ofFn, Function.comp_def, Finset.sum_apply]

/-- **The surface term over ℝ**: zero for EPSLS = 0 (conducting boundary), otherwise the
theory's `surface` and `surfaceForce` with `ε₀ = 1/eps0_inv`. -/
theorem surfaceModel_real_eq (c : V3 ℝ) {ε : ℝ} (hε : ε ≠ 0) (V e : ℝ) :
    (surfaceModel ℝ (atomsOf pos q) c ε V e).1 =
        surface e⁻¹ V ε q (fun i => (pos i).toVec) c.toVec ∧
      (surfaceModel ℝ (atomsOf pos q) c ε V e).2.1.map V3.toVec =
        List.ofFn fun i => surfaceForce e⁻¹ V ε q (fun i => (pos i).toVec) c.toVec i := by
  have hc : ¬(¬ε < 0 ∧ ¬0 < ε) := fun h => hε (le_antisymm (not_lt.mp h.2) (not_lt.mp h.1))
  have hM := dipoleModel_real_eq pos q c
  simp only [surfaceModel, hc, ↓reduceIte]
  refine ⟨?_, ?_⟩
  · rw [surface, ← hM, ← k2_model_eq]
    simp only [div_eq_mul_inv, mul_inv_rev, inv_inv]
    ring
  · simp only [atomsOf, List.map_ofFn, Function.comp_def]
    congr 1
    funext i
    rw [surfaceForce, ← hM]
    have hv : ∀ (t : ℝ) (m : V3 ℝ), (⟨t * m.x, t * m.y, t * m.z⟩ : V3 ℝ).toVec = t • m.toVec := by
      intro t m
      ext j
      fin_cases j <;> rfl
    rw [hv]
    congr 1
    simp only [div_eq_mul_inv, mul_inv_rev, inv_inv]
    ring

theorem surfaceModel_real_zero (c : V3 ℝ) (V e : ℝ) :
    surfaceModel ℝ (atomsOf pos q) c 0 V e = (0, (atomsOf pos q).map fun _ => V3.zero, Sym.zero) := by
  simp [surfaceModel]

end KSpaceReal

/-- `A2` in a cube over ℝ is the theory's `a2Cube`. -/
theorem a2CubicModel_real_eq (W L A1 A3 : ℝ) : (a2CubicModel ℝ W L A1 A3).1 = a2Cube W L A1 A3 :=
  rfl


/-! ### `A2` by shells over ℝ -/

theorem mem_faceTriples {n : ℕ} (hn : 1 ≤ n) {t : ℤ × ℤ × ℤ} :
    t ∈ faceTriples n ↔
      (t.1 = n ∧ -(n : ℤ) ≤ t.2.1 ∧ t.2.1 ≤ n ∧ -(n : ℤ) ≤ t.2.2 ∧ t.2.2 ≤ n) ∨
        (t.2.1 = n ∧ -(n : ℤ) < t.1 ∧ t.1 < n ∧ -(n : ℤ) ≤ t.2.2 ∧ t.2.2 ≤ n) ∨
        (t.2.2 = n ∧ -(n : ℤ) < t.1 ∧ t.1 < n ∧ -(n : ℤ) < t.2.1 ∧ t.2.1 < n) := by
  obtain ⟨a, b, c⟩ := t
  have hc : ((n - 1 : ℕ) : ℤ) = n - 1 := by omega
  simp only [faceTriples, List.mem_append, List.mem_flatMap, List.mem_map, mem_intRange, hc,
    Prod.mk.injEq]
  constructor
  · rintro ((⟨x, hx, y, hy, h1, h2, h3⟩ | ⟨x, hx, y, hy, h1, h2, h3⟩) | ⟨x, hx, y, hy, h1, h2, h3⟩) <;>
      subst h1 h2 h3 <;> omega
  · rintro (⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩)
    · exact Or.inl (Or.inl ⟨b, ⟨h2.1, h2.2.1⟩, c, ⟨h2.2.2.1, h2.2.2.2⟩, h1.symm, rfl, rfl⟩)
    · exact Or.inl (Or.inr ⟨a, by omega, c, ⟨h2.2.2.1, h2.2.2.2⟩, rfl, h1.symm, rfl⟩)
    · exact Or.inr ⟨a, by omega, b, by omega, rfl, rfl, h1.symm⟩

theorem faceTriples_nodup (n : ℕ) (hn : 1 ≤ n) : (faceTriples n).Nodup := by
  have part : ∀ (m k : ℕ) (f : ℤ → ℤ → ℤ × ℤ × ℤ), (∀ a b a' b', f a b = f a' b' → a = a' ∧ b = b') →
      ((intRange m).flatMap fun a => (intRange k).map fun b => f a b).Nodup := by
    intro m k f hf
    rw [List.nodup_flatMap]
    refine ⟨fun a _ => (intRange_nodup k).map fun b b' h => (hf a b a b' h).2, ?_⟩
    refine (intRange_nodup m).pairwise_of_forall_ne fun a _ a' _ hne => ?_
    rw [Function.onFun, List.disjoint_left]
    intro x hx hx'
    obtain ⟨b, -, rfl⟩ := List.mem_map.mp hx
    obtain ⟨b', -, h⟩ := List.mem_map.mp hx'
    exact hne (hf a' b' a b h).1.symm
  have hc : ((n - 1 : ℕ) : ℤ) = n - 1 := by omega
  simp only [faceTriples, List.nodup_append, List.mem_append, List.mem_flatMap, List.mem_map,
    mem_intRange, hc]
  refine ⟨⟨part _ _ _ ?_, part _ _ _ ?_, ?_⟩, part _ _ _ ?_, ?_⟩
  · intro a b a' b' h
    simp only [Prod.mk.injEq] at h
    exact ⟨h.2.1, h.2.2⟩
  · intro a b a' b' h
    simp only [Prod.mk.injEq] at h
    exact ⟨h.1, h.2.2⟩
  · rintro _ ⟨x, hx, y, hy, rfl⟩ _ ⟨x', hx', y', hy', rfl⟩ h
    simp only [Prod.mk.injEq] at h
    omega
  · intro a b a' b' h
    simp only [Prod.mk.injEq] at h
    exact ⟨h.1, h.2.1⟩
  · rintro _ (⟨x, hx, y, hy, rfl⟩ | ⟨x, hx, y, hy, rfl⟩) _ ⟨x', hx', y', hy', rfl⟩ h <;>
      simp only [Prod.mk.injEq] at h <;> omega

/-- **gromosXX's enumeration of a half shell is the theory's `face`**, each triple once. -/
theorem faceTriples_toFinset {n : ℕ} (hn : 1 ≤ n) :
    ((faceTriples n).map toFin).toFinset = face n ∧ ((faceTriples n).map toFin).Nodup := by
  classical
  refine ⟨?_, (faceTriples_nodup n hn).map toFin_injective⟩
  ext l
  have hl : toFin (l 0, l 1, l 2) = l := by
    funext j
    fin_cases j <;> rfl
  have three : ∀ (P : Fin 3 → Prop), (∀ j, P j) ↔ P 0 ∧ P 1 ∧ P 2 := fun P =>
    ⟨fun h => ⟨h 0, h 1, h 2⟩, fun ⟨h0, h1, h2⟩ j => by fin_cases j <;> assumption⟩
  simp only [List.mem_toFinset, List.mem_map, face, Finset.mem_filter, mem_cube, three]
  constructor
  · rintro ⟨t, ht, rfl⟩
    have := (mem_faceTriples hn).mp ht
    simp only [toFin, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.head_cons, Matrix.tail_cons]
    omega
  · intro h
    refine ⟨(l 0, l 1, l 2), (mem_faceTriples hn).mpr ?_, hl⟩
    simp only
    omega

/-- The first component of a fold that adds `f` to it and does anything to the second. -/
theorem foldl_fst_add {β σ : Type} (l : List β) (f : β → ℝ) (g : ℝ × σ → β → σ) (a : ℝ) (s : σ) :
    (l.foldl (fun (acc : ℝ × σ) x => (acc.1 + f x, g acc x)) (a, s)).1 = a + (l.map f).sum := by
  induction l generalizing a s with
  | nil => simp
  | cons x l ih => simp only [List.foldl_cons, ih, List.map_cons, List.sum_cons]; ring

/-- **One half shell over ℝ** is the sum of `a2Term` over the theory's `face`, for a width `a ≠ 0`
and a basis that sends no triple `l ≠ 0` to `k = 0`. -/
theorem a2ShellModel_real_eq (s : Shape) {width : ℝ} (hw : width ≠ 0) (u v w : V3 ℝ)
    (hB : ∀ l : Fin 3 → ℤ, l ≠ 0 → latticeVector ![u.toVec, v.toVec, w.toVec] l ⬝ᵥ
      latticeVector ![u.toVec, v.toVec, w.toVec] l ≠ 0)
    {n : ℕ} (hn : 1 ≤ n) (sg : Sym ℝ) :
    (a2ShellModel ℝ s width u v w n sg).1 =
      ∑ l ∈ face n, a2Term s width ![u.toVec, v.toVec, w.toVec] l := by
  classical
  obtain ⟨hset, hnd⟩ := faceTriples_toFinset hn
  rw [← hset, List.sum_toFinset _ hnd, List.map_map]
  simp only [a2ShellModel]
  rw [foldl_fst_add, zero_add]
  refine congrArg List.sum (List.map_congr_left fun t ht => ?_)
  have hne : toFin t ≠ 0 := by
    intro h0
    have := (mem_faceTriples hn).mp ht
    have h1 := congrFun h0 0
    have h2 := congrFun h0 1
    have h3 := congrFun h0 2
    simp only [toFin, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.head_cons, Matrix.tail_cons, Pi.zero_apply] at h1 h2 h3
    omega
  have hk := hB _ hne
  have hk2 := k2_model_eq (latticeVectorModel ℝ u v w t)
  rw [latticeVectorModel_real_eq] at hk2
  have hpos : 0 < latticeVector ![u.toVec, v.toVec, w.toVec] (toFin t) ⬝ᵥ
      latticeVector ![u.toVec, v.toVec, w.toVec] (toFin t) :=
    lt_of_le_of_ne (by rw [← hk2]; exact add_nonneg (add_nonneg (mul_self_nonneg _)
      (mul_self_nonneg _)) (mul_self_nonneg _)) (Ne.symm hk)
  have hak : √(latticeVector ![u.toVec, v.toVec, w.toVec] (toFin t) ⬝ᵥ
      latticeVector ![u.toVec, v.toVec, w.toVec] (toFin t)) * width ≠ 0 :=
    mul_ne_zero (Real.sqrt_pos.mpr hpos).ne' hw
  simp only [Function.comp_apply, hk2, Sqrt.sqrt, a2Term]
  rw [shapeFourierModel_real_eq s hak, mul_comm width]

/-- The shell loop over ℝ: what it adds to `a2` is the shell sums `× 2` of the shells it ran. -/
theorem a2ShellsLoop_real (s : Shape) {width : ℝ} (hw : width ≠ 0) (u v w : V3 ℝ)
    (hB : ∀ l : Fin 3 → ℤ, l ≠ 0 → latticeVector ![u.toVec, v.toVec, w.toVec] l ⬝ᵥ
      latticeVector ![u.toVec, v.toVec, w.toVec] l ≠ 0) (tola2 : ℝ) :
    ∀ (fuel lmax : ℕ) (a2 : ℝ) (sg : Sym ℝ),
      (a2ShellsLoop ℝ s width u v w tola2 fuel lmax a2 sg).1 =
          a2 + ∑ m ∈ Finset.Ico (lmax + 1) ((a2ShellsLoop ℝ s width u v w tola2 fuel lmax a2 sg).2.1 + 1),
            (∑ l ∈ face m, a2Term s width ![u.toVec, v.toVec, w.toVec] l) * 2 ∧
        lmax ≤ (a2ShellsLoop ℝ s width u v w tola2 fuel lmax a2 sg).2.1 := by
  intro fuel
  induction fuel with
  | zero => intro lmax a2 sg; simp [a2ShellsLoop]
  | succ fuel ih =>
    intro lmax a2 sg
    have hsh := a2ShellModel_real_eq s hw u v w hB (Nat.le_add_left 1 lmax) sg
    have one : a2 + (a2ShellModel ℝ s width u v w (lmax + 1) sg).1 * 2 =
          a2 + ∑ m ∈ Finset.Ico (lmax + 1) (lmax + 1 + 1),
            (∑ l ∈ face m, a2Term s width ![u.toVec, v.toVec, w.toVec] l) * 2 := by
      rw [Nat.Ico_succ_singleton, Finset.sum_singleton, hsh]
    rw [a2ShellsLoop.eq_def]
    dsimp only
    split_ifs with h
    · cases fuel with
      | zero =>
        dsimp only
        exact ⟨one, by omega⟩
      | succ f =>
        dsimp only
        obtain ⟨h1, h2⟩ := ih (lmax + 1) (a2 + (a2ShellModel ℝ s width u v w (lmax + 1) sg).1 * 2)
          ((a2ShellModel ℝ s width u v w (lmax + 1) sg).2.addDiag
            ((a2ShellModel ℝ s width u v w (lmax + 1) sg).1 * 2))
        refine ⟨?_, le_trans (Nat.le_succ _) h2⟩
        rw [h1, Finset.sum_eq_sum_Ico_succ_bot (Nat.lt_succ_of_le h2), hsh]
        ring
    · dsimp only
      exact ⟨one, by omega⟩

/-- **`A2` by shells over ℝ is the theory's `a2Shells`** over the shells the model summed
(`a2Shells_eq_cube`: the lattice sum over the cube of those shells), for a width `a ≠ 0` and a
basis that sends no triple `l ≠ 0` to `k = 0`. -/
theorem a2ShellsModel_real_eq (s : Shape) {width : ℝ} (hw : width ≠ 0) (u v w : V3 ℝ)
    (hB : ∀ l : Fin 3 → ℤ, l ≠ 0 → latticeVector ![u.toVec, v.toVec, w.toVec] l ⬝ᵥ
      latticeVector ![u.toVec, v.toVec, w.toVec] l ≠ 0) (V tola2 : ℝ) (maxShells : ℕ) :
    (a2ShellsModel ℝ Real.pi s width u v w V tola2 maxShells).1 =
      a2Shells s width ![u.toVec, v.toVec, w.toVec] V
        (a2ShellsModel ℝ Real.pi s width u v w V tola2 maxShells).2.1 := by
  obtain ⟨h1, -⟩ := a2ShellsLoop_real s hw u v w hB tola2 maxShells 0 0 Sym.zero
  simp only [a2ShellsModel, a2Shells]
  rw [h1, zero_add, Finset.Ico_add_one_right_eq_Icc, Finset.mul_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun m _ => ?_
  ring

end NewenSpec.LatticeSum
