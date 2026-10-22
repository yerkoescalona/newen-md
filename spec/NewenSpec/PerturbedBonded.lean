import NewenSpec.Float
import NewenSpec.Bonded
import NewenTheory.Interactions.PerturbedBonded

/-!
# Model of the perturbed bonded kernels of `crates/newen-forces`

Rust, all in `src/bonded/perturbed.rs` (spec keys in `spec/SPEC.toml`):

* `bonded.perturbed_quartic_bond_term`: `perturbed_quartic_bond_term`;
* `bonded.perturbed_angle_term`: `perturbed_angle_term` (cosine-harmonic);
* `bonded.perturbed_improper_term`: `perturbed_improper_term`, with the wrap of `ζ` into
  `[ζ0(λ) − π, ζ0(λ) + π]` (the theory's `wrapZeta`, `wrapModel_real_eq`);
* `bonded.chebyshev`, `bonded.perturbed_dihedral_term`: `cos(m φ)` as a polynomial in `cos φ`,
  and the proper dihedral that mixes the two states' energies;
* `bonded.soft_perturbed_term`: `soft_perturbed_term`, the soft-core form of the bond, the angle
  and the improper (PERTBONDSOFT, PERTANGLESOFT, PERTIMPROPERDIHSOFT) and of their PRECALCLAM rows.

Theory: `NewenTheory/Interactions/PerturbedBonded.lean`. As in `NewenSpec/SoftCore.lean`, each
kernel is written once over any number type, operation by operation like the Rust, proven equal
to the theory at `ℝ` and executed at `Float` by `lake exe spec-vectors`;
`crates/newen-forces/tests/perturbed_bonded_spec.rs` compares the Rust bit for bit. The decimal
literals `0.25`, `0.5` and `π` are parameters, given at `Float` by their bits. `x.powi(k)` is
written as the multiplications compiler-rt's `__powidf2` and LLVM's expansion perform
(`x · (x · x)` for `k = 3`, `x² · x²` for `4`, `x · x⁴` for `5`, `x² · x⁴` for `6`).

What is not proved: the rounding error between the two instantiations.
-/

namespace NewenSpec.PerturbedBonded

open NewenTheory.Interactions.PerturbedBonded
open NewenTheory.Interactions.Bonded (wrapZeta wrapZeta_sub_eq_int_mul wrapZeta_eq_of_abs_lt)
open Filter Topology

/-! ## 1. Model -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α]
  [DecidableRel (α := α) (· < ·)] [OfNat α 1] [OfNat α 2] [OfNat α 4]

/-- `perturbed_quartic_bond_term(dist2, k_a, k_b, r0_a, r0_b, lambda, lambda_derivative)`:
`(energy, f, dE/dλ)`, with the literal `0.25` as `quarter`. -/
def quarticModel (quarter dist2 kA kB r0A r0B lam dlam : α) : α × α × α :=
  let k := (1 - lam) * kA + lam * kB
  let r0 := (1 - lam) * r0A + lam * r0B
  let kD := kB - kA
  let bD := r0B - r0A
  let r02 := r0 * r0
  let delta := dist2 - r02
  (quarter * k * delta * delta, (-k) * delta,
    quarter * dlam * (kD * delta * delta - 4 * k * bD * r0 * delta))

/-- `perturbed_angle_term(cost, k_a, k_b, cos0a, cos0b, lambda, lambda_derivative)`:
`(energy, df, dE/dλ)`, with the literal `0.5` as `half`. -/
def angleModel (half cost kA kB cos0A cos0B lam dlam : α) : α × α × α :=
  let k := (1 - lam) * kA + lam * kB
  let cos0 := (1 - lam) * cos0A + lam * cos0B
  let kD := kB - kA
  let cdD := cos0B - cos0A
  let dCos := cost - cos0
  (half * k * dCos * dCos, (-k) * dCos,
    half * dlam * (kD * dCos * dCos - 2 * k * cdD * dCos))

/-- `while ζ < ζ0 − π { ζ += 2π }`, at most `fuel` times. -/
def wrapUp (pi q0 : α) : Nat → α → α
  | 0, z => z
  | n + 1, z => if z < q0 - pi then wrapUp pi q0 n (z + 2 * pi) else z

/-- `while ζ > ζ0 + π { ζ -= 2π }`, at most `fuel` times. -/
def wrapDown (pi q0 : α) : Nat → α → α
  | 0, z => z
  | n + 1, z => if q0 + pi < z then wrapDown pi q0 n (z - 2 * pi) else z

/-- `wrap_around(zeta, q0)`; the 64 turns of fuel are far more than any angle needs
(`wrapModel_real_eq` covers every `|ζ − ζ0| ≤ 129 π`). -/
def wrapModel (pi zeta q0 : α) : α := wrapDown α pi q0 64 (wrapUp α pi q0 64 zeta)

/-- `perturbed_improper_term(zeta, k_a, k_b, q0_a, q0_b, lambda, lambda_derivative)`:
`(energy, −K·Δ, dE/dλ, ζ)`. -/
def improperModel (half pi zeta kA kB q0A q0B lam dlam : α) : α × α × α × α :=
  let k := (1 - lam) * kA + lam * kB
  let q0 := (1 - lam) * q0A + lam * q0B
  let kD := kB - kA
  let qD := q0B - q0A
  let zetaAdj := wrapModel α pi zeta q0
  let dZ := zetaAdj - q0
  (half * k * dZ * dZ, (-k) * dZ,
    half * dlam * ((-2) * k * qD * dZ + kD * dZ * dZ), zetaAdj)

/-- `soft_perturbed_term(x, k_a, k_b, x0_a, x0_b, alpha, lambda, lambda_derivative)`:
`(energy, K_eff·Δ, dE/dλ)`. -/
def softModel (half x kA kB x0A x0B alpha lam dlam : α) : α × α × α :=
  let x0 := (1 - lam) * x0A + lam * x0B
  let x0D := x0B - x0A
  let diff := x - x0
  let diff2 := diff * diff
  let sA := 1 + alpha * lam * diff2
  let sB := 1 + alpha * (1 - lam) * diff2
  let sA2 := sA * sA
  let sB2 := sB * sB
  let kEff := (1 - lam) * kA / sA2 + lam * kB / sB2
  let kSoft := (1 - lam) * kA / sA + lam * kB / sB
  let st1 := 1 + alpha * diff2
  let st2 := (-2) * alpha * lam * (1 - lam) * diff * x0D
  (half * kSoft * diff2, kEff * diff,
    dlam * (half * diff2 * (kA / sA2 * ((-st1) - st2) + kB / sB2 * (st1 - st2)) -
      kSoft * diff * x0D))

end Model

section Chebyshev

variable (α : Type) [Add α] [Sub α] [Mul α] [OfNat α 0] [OfNat α 1] [OfNat α 2] [OfNat α 3]
  [OfNat α 4] [OfNat α 5] [OfNat α 8] [OfNat α 12] [OfNat α 16] [OfNat α 18] [OfNat α 20]
  [OfNat α 32] [OfNat α 36] [OfNat α 48] [OfNat α 60] [OfNat α 80] [OfNat α 192]

/-- `x.powi(3)`. -/
def powi3 (x : α) : α := x * x * x
/-- `x.powi(4)`. -/
def powi4 (x : α) : α := x * x * (x * x)
/-- `x.powi(5)`. -/
def powi5 (x : α) : α := x * (x * x * (x * x))
/-- `x.powi(6)`. -/
def powi6 (x : α) : α := x * x * (x * x * (x * x))

/-- `chebyshev(m, cos_phi)`: `(cos(m φ), d cos(m φ)/d cos φ)`. -/
def chebyshevModel (m : Int) (c : α) : α × α :=
  if m == 1 then (c, 1)
  else if m == 2 then (2 * c * c - 1, 4 * c)
  else if m == 3 then (4 * powi3 α c - 3 * c, 12 * c * c - 3)
  else if m == 4 then (8 * powi4 α c - 8 * c * c + 1, 32 * powi3 α c - 16 * c)
  else if m == 5 then
    (16 * powi5 α c - 20 * powi3 α c + 5 * c, 80 * powi4 α c - 60 * c * c + 5)
  else if m == 6 then
    (32 * powi6 α c - 48 * powi4 α c + 18 * c * c - 1,
      192 * powi5 α c - 192 * powi3 α c + 36 * c)
  else (0, 0)

/-- `perturbed_dihedral_term(cos_phi, (k_a, cospd_a, m_a), (k_b, cospd_b, m_b), lambda,
lambda_derivative)`: `(e_a, e_b, g_a, g_b, energy, dE/dλ)`. -/
def dihedralModel [Neg α] (c kA cospdA : α) (mA : Int) (kB cospdB : α) (mB : Int) (lam dlam : α) :
    α × α × α × α × α × α :=
  let (cmA, dcmA) := chebyshevModel α mA c
  let eA := kA * (1 + cospdA * cmA)
  let gA := (-kA) * cospdA * dcmA
  let (cmB, dcmB) := chebyshevModel α mB c
  let eB := kB * (1 + cospdB * cmB)
  let gB := (-kB) * cospdB * dcmB
  (eA, eB, gA, gB, (1 - lam) * eA + lam * eB, dlam * (eB - eA))

end Chebyshev

/-! ## 2. The `Float` instantiation -/

/-- The literal `0.25`, as bits. -/
def quarterFloat : Float := Float.ofBits 0x3fd0000000000000
/-- The literal `0.5`, as bits. -/
def halfFloat : Float := Float.ofBits 0x3fe0000000000000
/-- `std::f64::consts::PI`, as bits. -/
def piFloat : Float := Float.ofBits 0x400921fb54442d18

def quarticFloat (dist2 kA kB r0A r0B lam dlam : Float) : Float × Float × Float :=
  quarticModel Float quarterFloat dist2 kA kB r0A r0B lam dlam

def angleFloat (cost kA kB cos0A cos0B lam dlam : Float) : Float × Float × Float :=
  angleModel Float halfFloat cost kA kB cos0A cos0B lam dlam

def improperFloat (zeta kA kB q0A q0B lam dlam : Float) : Float × Float × Float × Float :=
  improperModel Float halfFloat piFloat zeta kA kB q0A q0B lam dlam

def softFloat (x kA kB x0A x0B alpha lam dlam : Float) : Float × Float × Float :=
  softModel Float halfFloat x kA kB x0A x0B alpha lam dlam

def chebyshevFloat (m : Int) (c : Float) : Float × Float := chebyshevModel Float m c

def dihedralFloat (c kA cospdA : Float) (mA : Int) (kB cospdB : Float) (mB : Int) (lam dlam : Float) :
    Float × Float × Float × Float × Float × Float :=
  dihedralModel Float c kA cospdA mA kB cospdB mB lam dlam

/-! ## 3. The ℝ instantiation equals the theory -/

/-- The quartic bond kernel over ℝ: the mixed energy, the force factor `−K(λ)(b² − b0(λ)²)`,
and `pertQuarticDvdl`. -/
theorem quarticModel_real_eq (d2 kA kB r0A r0B L dL : ℝ) :
    quarticModel ℝ (1 / 4) d2 kA kB r0A r0B L dL =
      (pertQuarticEnergy kA kB r0A r0B d2 L, -mix kA kB L * (d2 - mix r0A r0B L ^ 2),
        pertQuarticDvdl kA kB r0A r0B d2 L dL) := by
  simp only [quarticModel, pertQuarticEnergy, quarticEnergy, pertQuarticDvdl, mix, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_⟩ <;> ring

/-- The cosine-harmonic angle kernel over ℝ. -/
theorem angleModel_real_eq (cost kA kB cA cB L dL : ℝ) :
    angleModel ℝ (1 / 2) cost kA kB cA cB L dL =
      (pertHarmonicEnergy kA kB cA cB cost L, -mix kA kB L * (cost - mix cA cB L),
        pertHarmonicDvdl kA kB cA cB cost L dL) := by
  simp only [angleModel, pertHarmonicEnergy, harmonicEnergy, pertHarmonicDvdl, mix,
    Prod.mk.injEq]
  refine ⟨?_, ?_, ?_⟩ <;> first | trivial | ring

/-- The loops of `wrap_around` are those of the unperturbed improper kernel
(`NewenSpec.Bonded.wrapUp`, `wrapDown`) with bounds `ζ0 ∓ π` and step `2π`. -/
theorem wrapUp_eq_bonded (pi q0 : ℝ) :
    ∀ (n : Nat) (z : ℝ), wrapUp ℝ pi q0 n z = NewenSpec.Bonded.wrapUp ℝ n (q0 - pi) (2 * pi) z
  | 0, _ => rfl
  | n + 1, z => by
    simp only [wrapUp, NewenSpec.Bonded.wrapUp, wrapUp_eq_bonded pi q0 n]

theorem wrapDown_eq_bonded (pi q0 : ℝ) :
    ∀ (n : Nat) (z : ℝ), wrapDown ℝ pi q0 n z = NewenSpec.Bonded.wrapDown ℝ n (q0 + pi) (2 * pi) z
  | 0, _ => rfl
  | n + 1, z => by
    simp only [wrapDown, NewenSpec.Bonded.wrapDown, wrapDown_eq_bonded pi q0 n]

/-- **`wrap_around` over ℝ is `wrapZeta`**, for every `ζ` within `129 π` of `ζ0` (the 64 turns
of fuel, `NewenSpec.Bonded.wrap_real_eq_wrapZeta`). -/
theorem wrapModel_real_eq {zeta q0 : ℝ} (h : |zeta - q0| ≤ 129 * Real.pi) :
    wrapModel ℝ Real.pi zeta q0 = wrapZeta q0 zeta := by
  rw [wrapModel, wrapUp_eq_bonded, wrapDown_eq_bonded]
  exact NewenSpec.Bonded.wrap_real_eq_wrapZeta 64 (by norm_num; linarith)

/-- The improper kernel over ℝ, for every `ζ` within `129 π` of `ζ0(λ)`: the mixed energy,
`−K(λ)·Δ` and `pertHarmonicDvdl` at the wrapped angle `wrapZeta ζ0(λ) ζ`, which it also returns. -/
theorem improperModel_real_eq (zeta kA kB qA qB L dL : ℝ)
    (hfuel : |zeta - mix qA qB L| ≤ 129 * Real.pi) :
    improperModel ℝ (1 / 2) Real.pi zeta kA kB qA qB L dL =
      (pertHarmonicEnergy kA kB qA qB (wrapZeta (mix qA qB L) zeta) L,
        -mix kA kB L * (wrapZeta (mix qA qB L) zeta - mix qA qB L),
        pertHarmonicDvdl kA kB qA qB (wrapZeta (mix qA qB L) zeta) L dL,
        wrapZeta (mix qA qB L) zeta) := by
  have hw : wrapModel ℝ Real.pi zeta ((1 - L) * qA + L * qB) = wrapZeta (mix qA qB L) zeta :=
    wrapModel_real_eq hfuel
  simp only [improperModel, hw, pertHarmonicEnergy, harmonicEnergy, pertHarmonicDvdl, mix,
    Prod.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> first | trivial | ring

/-- The soft-core kernel over ℝ: `softEnergy`, `softGrad`, `softDvdl`. -/
theorem softModel_real_eq (x kA kB xA xB a L dL : ℝ) :
    softModel ℝ (1 / 2) x kA kB xA xB a L dL =
      (softEnergy kA kB xA xB a x L, softGrad kA kB xA xB a x L, softDvdl kA kB xA xB a x L dL) := by
  simp only [softModel, softEnergy, softGrad, softDvdl, mix, Prod.mk.injEq]
  simp only [div_eq_mul_inv, mul_inv, ← inv_pow]
  ring_nf
  simp only [and_self]

/-- `chebyshev` over ℝ is `(chebT m, chebDT m)` for `m = 1 … 6`, hence `cos(m φ)` at
`c = cos φ` (`chebT_cos`), and `(0, 0)` for `m = 0`. -/
theorem chebyshevModel_real_eq (c : ℝ) {m : ℕ} (h1 : 1 ≤ m) (h6 : m ≤ 6) :
    chebyshevModel ℝ (m : ℤ) c = (chebT m c, chebDT m c) := by
  obtain rfl | rfl | rfl | rfl | rfl | rfl : m = 1 ∨ m = 2 ∨ m = 3 ∨ m = 4 ∨ m = 5 ∨ m = 6 := by
    omega
  all_goals
    simp [chebyshevModel, chebT, chebDT, powi3, powi4, powi5, powi6]
    try (first | (constructor <;> ring) | ring)

theorem chebyshevModel_zero (c : ℝ) : chebyshevModel ℝ 0 c = (0, 0) := by
  simp [chebyshevModel]

/-- The dihedral kernel over ℝ mixes the two states' `dihedralEnergy` and returns
`(E_B − E_A) L′`, the derivative of the mixture (`pert_dihedral_dvdl_is_deriv`). -/
theorem dihedralModel_real_eq (c kA cA kB cB L dL : ℝ) {mA mB : ℕ} (hA1 : 1 ≤ mA) (hA6 : mA ≤ 6)
    (hB1 : 1 ≤ mB) (hB6 : mB ≤ 6) :
    let o := dihedralModel ℝ c kA cA mA kB cB mB L dL
    o.2.2.2.2.1 = mix (dihedralEnergy kA cA mA c) (dihedralEnergy kB cB mB c) L ∧
      o.2.2.2.2.2 = (dihedralEnergy kB cB mB c - dihedralEnergy kA cA mA c) * dL := by
  simp only [dihedralModel, chebyshevModel_real_eq c hA1 hA6, chebyshevModel_real_eq c hB1 hB6,
    dihedralEnergy, mix]
  refine ⟨?_, ?_⟩ <;> first | trivial | ring

/-! ## 4. The kernels' `dE/dλ` are the derivatives of their energies -/

/-- The quartic bond kernel, with `L(λ)` differentiable. -/
theorem quarticModel_dvdl_is_deriv (d2 kA kB r0A r0B : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => (quarticModel ℝ (1 / 4) d2 kA kB r0A r0B (Λ l) dΛ).1)
      (quarticModel ℝ (1 / 4) d2 kA kB r0A r0B (Λ lam) dΛ).2.2 lam := by
  simp only [quarticModel_real_eq]
  exact pert_quartic_dvdl_is_deriv kA kB r0A r0B d2 hΛ

/-- The angle kernel. -/
theorem angleModel_dvdl_is_deriv (cost kA kB cA cB : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) :
    HasDerivAt (fun l => (angleModel ℝ (1 / 2) cost kA kB cA cB (Λ l) dΛ).1)
      (angleModel ℝ (1 / 2) cost kA kB cA cB (Λ lam) dΛ).2.2 lam := by
  simp only [angleModel_real_eq]
  exact pert_harmonic_dvdl_is_deriv kA kB cA cB cost hΛ

/-- **The improper kernel's `dE/dλ` is the derivative of its energy** on every λ-interval where
the wrapped angle does not jump: `ζ0(λ)` moves with λ, and the wrap jumps by `2π` where `ζ` meets
`ζ0(λ) ± π`. The hypothesis `hjump` keeps the wrapped angle strictly inside
`(ζ0(λ) − π, ζ0(λ) + π)`; there it stays the same turn of `ζ` for all λ near `lam`. -/
theorem improperModel_dvdl_is_deriv (zeta kA kB qA qB : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam)
    (hfuel : |zeta - mix qA qB (Λ lam)| < 129 * Real.pi)
    (hjump : |wrapZeta (mix qA qB (Λ lam)) zeta - mix qA qB (Λ lam)| < Real.pi) :
    HasDerivAt (fun l => (improperModel ℝ (1 / 2) Real.pi zeta kA kB qA qB (Λ l) dΛ).1)
      (improperModel ℝ (1 / 2) Real.pi zeta kA kB qA qB (Λ lam) dΛ).2.2.1 lam := by
  have hq : ContinuousAt (fun l => mix qA qB (Λ l)) lam := (hasDerivAt_mix qA qB hΛ).continuousAt
  set w := wrapZeta (mix qA qB (Λ lam)) zeta
  obtain ⟨k, hk⟩ := wrapZeta_sub_eq_int_mul (mix qA qB (Λ lam)) zeta
  have ev1 : ∀ᶠ l in 𝓝 lam, |zeta - mix qA qB (Λ l)| < 129 * Real.pi :=
    (continuousAt_const.sub hq).abs.eventually_lt continuousAt_const hfuel
  have ev2 : ∀ᶠ l in 𝓝 lam, |w - mix qA qB (Λ l)| < Real.pi :=
    (continuousAt_const.sub hq).abs.eventually_lt continuousAt_const hjump
  have heq : (fun l => (improperModel ℝ (1 / 2) Real.pi zeta kA kB qA qB (Λ l) dΛ).1)
      =ᶠ[𝓝 lam] fun l => pertHarmonicEnergy kA kB qA qB w (Λ l) := by
    filter_upwards [ev1, ev2] with l h1 h2
    rw [improperModel_real_eq _ _ _ _ _ _ _ h1.le]
    have hzk : zeta + 2 * Real.pi * k = w := by linarith
    rw [wrapZeta_eq_of_abs_lt k (by rw [hzk]; exact h2), hzk]
  rw [improperModel_real_eq _ _ _ _ _ _ _ hfuel.le]
  exact (pert_harmonic_dvdl_is_deriv kA kB qA qB w hΛ).congr_of_eventuallyEq heq

/-- The soft-core kernel, wherever both softening denominators are nonzero. -/
theorem softModel_dvdl_is_deriv (x kA kB xA xB a : ℝ) {Λ : ℝ → ℝ} {dΛ lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam)
    (hA : 1 + a * Λ lam * (x - mix xA xB (Λ lam)) ^ 2 ≠ 0)
    (hB : 1 + a * (1 - Λ lam) * (x - mix xA xB (Λ lam)) ^ 2 ≠ 0) :
    HasDerivAt (fun l => (softModel ℝ (1 / 2) x kA kB xA xB a (Λ l) dΛ).1)
      (softModel ℝ (1 / 2) x kA kB xA xB a (Λ lam) dΛ).2.2 lam := by
  simp only [softModel_real_eq]
  exact soft_dvdl_is_deriv kA kB xA xB a x hΛ hA hB

end NewenSpec.PerturbedBonded
