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
  `ls_self_outerloop`).

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
   `crates/newen-forces/tests/lattice_sum_spec.rs` compares the Rust kernels bit for bit.

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

end NewenSpec.LatticeSum
