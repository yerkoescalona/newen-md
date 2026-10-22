module

public import Mathlib.Analysis.SpecialFunctions.Sqrt
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse
public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Inv
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.IntervalCases

/-!
# The bonded terms

GROMOS (Vol. 2) holds a molecule together with four kinds of covalent term:

* a **bond** between atoms `i`, `j` at distance `b = |x_i − x_j|`: the quartic form
  `V = K/4 · (b² − b0²)²` (the GROMOS default; no square root) or the harmonic `V = K/2 · (b − b0)²`;
* a **bond angle** `θ` between `r_ij = x_i − x_j` and `r_kj = x_k − x_j`: the cosine-harmonic form
  `V = K/2 · (cos θ − cos θ0)²` (the default) or the harmonic `V = K/2 · (θ − θ0)²`;
* a **proper dihedral** `φ` of four atoms `i j k l`: `V = K · (1 + cos δ · cos mφ)` with
  `δ ∈ {0, π}` and `cos mφ` a polynomial in `cos φ` (COVALENTFORM NTBDN = 1), or
  `V = K · (1 + cos(mφ − δ))` for any `δ` (NTBDN = 0, the default);
* an **improper dihedral** `ζ`, which keeps a group planar or chiral: `V = K/2 · (ζ − ζ0)²`.

The two bond forms agree near the minimum when `K_harm = 2 K_quartic b0²`
(`quartic_harmonic_match`), and the two dihedral forms agree on every GROMOS topology because
`δ` is 0 or π there (`dihedral_new_eq_dihedral`).

A force is minus the gradient of the energy, and the forces of an isolated term sum to zero and
exert no torque. Here these are theorems about the force expressions newen evaluates (the same
expressions as `crates/newen-forces/src/bonded/`, over ℝ): `*_force_is_minus_grad`,
`*_forces_sum_zero`, `*_torque_zero`. How the kernels compute them in floating point is in
`NewenSpec/Bonded.lean`.

Vectors are `V3`, three real components, with `dot`, `cross` and `len` written out, so that
every identity below is a statement about nine (or six) real numbers.
-/

@[expose] public section

namespace NewenTheory.Interactions.Bonded

/-! ## Vectors of ℝ³ -/

/-- A vector of ℝ³: a bond vector (nm) or a force (kJ mol⁻¹ nm⁻¹). -/
@[ext] structure V3 where
  x : ℝ
  y : ℝ
  z : ℝ

namespace V3

instance : Add V3 := ⟨fun u v => ⟨u.x + v.x, u.y + v.y, u.z + v.z⟩⟩
instance : Sub V3 := ⟨fun u v => ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩⟩
instance : Neg V3 := ⟨fun u => ⟨-u.x, -u.y, -u.z⟩⟩
instance : Zero V3 := ⟨⟨0, 0, 0⟩⟩
instance : SMul ℝ V3 := ⟨fun a u => ⟨a * u.x, a * u.y, a * u.z⟩⟩

@[simp] theorem add_x (u v : V3) : (u + v).x = u.x + v.x := rfl
@[simp] theorem add_y (u v : V3) : (u + v).y = u.y + v.y := rfl
@[simp] theorem add_z (u v : V3) : (u + v).z = u.z + v.z := rfl
@[simp] theorem sub_x (u v : V3) : (u - v).x = u.x - v.x := rfl
@[simp] theorem sub_y (u v : V3) : (u - v).y = u.y - v.y := rfl
@[simp] theorem sub_z (u v : V3) : (u - v).z = u.z - v.z := rfl
@[simp] theorem neg_x (u : V3) : (-u).x = -u.x := rfl
@[simp] theorem neg_y (u : V3) : (-u).y = -u.y := rfl
@[simp] theorem neg_z (u : V3) : (-u).z = -u.z := rfl
@[simp] theorem zero_x : (0 : V3).x = 0 := rfl
@[simp] theorem zero_y : (0 : V3).y = 0 := rfl
@[simp] theorem zero_z : (0 : V3).z = 0 := rfl
@[simp] theorem smul_x (a : ℝ) (u : V3) : (a • u).x = a * u.x := rfl
@[simp] theorem smul_y (a : ℝ) (u : V3) : (a • u).y = a * u.y := rfl
@[simp] theorem smul_z (a : ℝ) (u : V3) : (a • u).z = a * u.z := rfl

/-- `u · v`. -/
def dot (u v : V3) : ℝ := u.x * v.x + u.y * v.y + u.z * v.z

/-- `u × v`. -/
def cross (u v : V3) : V3 :=
  ⟨u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y⟩

@[simp] theorem cross_x (u v : V3) : (cross u v).x = u.y * v.z - v.y * u.z := rfl
@[simp] theorem cross_y (u v : V3) : (cross u v).y = u.z * v.x - v.z * u.x := rfl
@[simp] theorem cross_z (u v : V3) : (cross u v).z = u.x * v.y - v.x * u.y := rfl

/-- `|u|`. -/
noncomputable def len (u : V3) : ℝ := Real.sqrt (dot u u)

theorem dot_self_nonneg (u : V3) : 0 ≤ dot u u := by
  unfold dot; nlinarith [mul_self_nonneg u.x, mul_self_nonneg u.y, mul_self_nonneg u.z]

theorem dot_self_pos {u : V3} (hu : u ≠ 0) : 0 < dot u u := by
  rcases (dot_self_nonneg u).eq_or_lt with h | h
  · exfalso
    apply hu
    simp only [dot] at h
    have hx : u.x * u.x = 0 := by
      nlinarith [mul_self_nonneg u.x, mul_self_nonneg u.y, mul_self_nonneg u.z]
    have hy : u.y * u.y = 0 := by
      nlinarith [mul_self_nonneg u.x, mul_self_nonneg u.y, mul_self_nonneg u.z]
    have hz : u.z * u.z = 0 := by
      nlinarith [mul_self_nonneg u.x, mul_self_nonneg u.y, mul_self_nonneg u.z]
    ext <;> simp [mul_self_eq_zero.mp hx, mul_self_eq_zero.mp hy, mul_self_eq_zero.mp hz]
  · exact h

theorem len_sq (u : V3) : len u ^ 2 = dot u u := Real.sq_sqrt (dot_self_nonneg u)

theorem len_nonneg (u : V3) : 0 ≤ len u := Real.sqrt_nonneg _

theorem len_pos {u : V3} (hu : u ≠ 0) : 0 < len u := Real.sqrt_pos.mpr (dot_self_pos hu)

/-- Lagrange's identity: `|u|²|v|² − (u · v)² = |u × v|²`. -/
theorem lagrange (u v : V3) : dot u u * dot v v - dot u v ^ 2 = dot (cross u v) (cross u v) := by
  simp only [dot, cross_x, cross_y, cross_z]; ring

/-- Cauchy–Schwarz: the cosine of the angle between two nonzero vectors lies in `[−1, 1]`, so
the clamps of the kernels never act over ℝ. -/
theorem cos_between_mem (u v : V3) :
    -1 ≤ dot u v / (len u * len v) ∧ dot u v / (len u * len v) ≤ 1 := by
  have hcs : dot u v ^ 2 ≤ dot u u * dot v v := by
    have := lagrange u v
    have := dot_self_nonneg (cross u v)
    linarith
  have habs : |dot u v| ≤ len u * len v := by
    rw [len, len, ← Real.sqrt_mul (dot_self_nonneg u)]
    exact Real.abs_le_sqrt hcs
  rcases (mul_nonneg (len_nonneg u) (len_nonneg v)).eq_or_lt with h0 | hpos
  · rw [← h0, div_zero]; norm_num
  · rw [abs_le] at habs
    constructor
    · rw [le_div_iff₀ hpos]; linarith
    · rw [div_le_iff₀ hpos]; linarith

end V3

open V3

/-- Write every vector expression as its components, keeping `u.x * u.x + u.y * u.y + u.z * u.z`
in the order `dot` gives it (so that `field_simp` matches it against a hypothesis). -/
macro "v3_expand" : tactic =>
  `(tactic| simp only [dot, V3.add_x, V3.add_y, V3.add_z, V3.sub_x, V3.sub_y, V3.sub_z,
    V3.neg_x, V3.neg_y, V3.neg_z, V3.smul_x, V3.smul_y, V3.smul_z, V3.cross_x, V3.cross_y,
    V3.cross_z, V3.zero_x, V3.zero_y, V3.zero_z])

/-! ## Bonds -/

/-- The quartic bond energy at length `b` (GROMOS Vol. 2): `K/4 · (b² − b0²)²`. -/
noncomputable def quarticBond (k b0 b : ℝ) : ℝ := k / 4 * (b ^ 2 - b0 ^ 2) ^ 2

/-- The harmonic bond energy at length `b`: `K/2 · (b − b0)²`. -/
noncomputable def harmonicBond (k b0 b : ℝ) : ℝ := k / 2 * (b - b0) ^ 2

/-- The force on atom `i` of a term that depends on `b = |r|` only, `r = x_i − x_j`, with
`dE/db = dEdb`: `−dEdb · r/b`. Atom `j` gets the opposite. -/
noncomputable def radialForce (dEdb : ℝ) (r : V3) : V3 := (-dEdb / len r) • r

/-- `dV/db` of the quartic bond. -/
theorem quarticBond_hasDerivAt (k b0 b : ℝ) :
    HasDerivAt (quarticBond k b0) (k * b * (b ^ 2 - b0 ^ 2)) b := by
  have h := (((hasDerivAt_pow 2 b).sub_const (b0 ^ 2)).pow 2).const_mul (k / 4)
  refine (h.congr_deriv ?_)
  push_cast
  ring

/-- `dV/db` of the harmonic bond. -/
theorem harmonicBond_hasDerivAt (k b0 b : ℝ) :
    HasDerivAt (harmonicBond k b0) (k * (b - b0)) b := by
  have h := (((hasDerivAt_id b).sub_const b0).pow 2).const_mul (k / 2)
  refine (h.congr_deriv ?_)
  simp only [id]
  push_cast
  ring

/-- A term `E(|r|)` with `E'(|r|) = dEdb`: along any direction `e`, the derivative of the energy
at `r ≠ 0` is `−radialForce dEdb r · e`. So `radialForce` is `−∇E`. -/
theorem radialForce_is_minus_grad {E : ℝ → ℝ} {dEdb : ℝ} {r : V3}
    (hE : HasDerivAt E dEdb (len r)) (hr : r ≠ 0) (e : V3) :
    HasDerivAt (fun t : ℝ => E (len (r + t • e))) (-dot (radialForce dEdb r) e) 0 := by
  have hd := dot_self_pos hr
  -- `|r + t e|² = r·r + 2t (r·e) + t² (e·e)`
  have hq : HasDerivAt (fun t : ℝ => dot (r + t • e) (r + t • e)) (2 * dot r e) 0 := by
    have h1 := ((hasDerivAt_id (0 : ℝ)).const_mul (2 * dot r e)).add
      (((hasDerivAt_id (0 : ℝ)).pow 2).const_mul (dot e e))
    have hfun : (fun t : ℝ => dot (r + t • e) (r + t • e)) =
        fun t => dot r r + ((2 * dot r e) * id t + dot e e * id t ^ 2) := by
      funext t; simp only [dot, add_x, add_y, add_z, smul_x, smul_y, smul_z, id]; ring
    rw [hfun]
    exact (h1.const_add (dot r r)).congr_deriv (by simp)
  have hq0 : dot (r + (0 : ℝ) • e) (r + (0 : ℝ) • e) = dot r r := by
    simp only [dot, add_x, add_y, add_z, smul_x, smul_y, smul_z]; ring
  have hn : HasDerivAt (fun t : ℝ => len (r + t • e))
      (2 * dot r e / (2 * Real.sqrt (dot r r))) 0 := by
    have h := hq.sqrt (by rw [hq0]; exact hd.ne')
    rw [hq0] at h
    exact h
  have hE' : HasDerivAt E dEdb (len (r + (0 : ℝ) • e)) := by
    have h0 : len (r + (0 : ℝ) • e) = len r := by unfold len; rw [hq0]
    rw [h0]; exact hE
  refine (hE'.comp (0 : ℝ) hn).congr_deriv ?_
  have hs : Real.sqrt (dot r r) ≠ 0 := (Real.sqrt_pos.mpr hd).ne'
  rw [radialForce, show len r = Real.sqrt (dot r r) from rfl]
  simp only [dot, smul_x, smul_y, smul_z] at hs ⊢
  field_simp
  ring

/-- **The quartic bond force is `−∇V`.** The kernel's force on `i`, `−K (b² − b0²) r`, is
`radialForce (dV/db) r`, and along every direction it is minus the derivative of the energy. -/
theorem quarticBond_force_is_minus_grad (k b0 : ℝ) {r : V3} (hr : r ≠ 0) :
    radialForce (k * len r * (len r ^ 2 - b0 ^ 2)) r = (-(k * (dot r r - b0 ^ 2))) • r ∧
      ∀ e : V3, HasDerivAt (fun t : ℝ => quarticBond k b0 (len (r + t • e)))
        (-dot (radialForce (k * len r * (len r ^ 2 - b0 ^ 2)) r) e) 0 := by
  refine ⟨?_, radialForce_is_minus_grad (quarticBond_hasDerivAt k b0 (len r)) hr⟩
  have hn := (len_pos hr).ne'
  rw [radialForce, ← len_sq r]
  congr 1
  field_simp

/-- **The harmonic bond force is `−∇V`**: `radialForce (K (b − b0)) r`, the kernel's
`−K (b − b0)/b · r`. -/
theorem harmonicBond_force_is_minus_grad (k b0 : ℝ) {r : V3} (hr : r ≠ 0) (e : V3) :
    HasDerivAt (fun t : ℝ => harmonicBond k b0 (len (r + t • e)))
      (-dot (radialForce (k * (len r - b0)) r) e) 0 :=
  radialForce_is_minus_grad (harmonicBond_hasDerivAt k b0 (len r)) hr e

/-- **Quartic and harmonic bonds match at the minimum** (GROMOS Vol. 2): the curvature of the
quartic bond at `b0` is `2 K b0²`, so a quartic bond of force constant `K` behaves near `b0` like
a harmonic bond of force constant `K_harm = 2 K b0²`. -/
theorem quartic_harmonic_match (k b0 : ℝ) :
    deriv (deriv (quarticBond k b0)) b0 = 2 * k * b0 ^ 2 ∧
      deriv (deriv (harmonicBond (2 * k * b0 ^ 2) b0)) b0 = 2 * k * b0 ^ 2 := by
  have hq : deriv (quarticBond k b0) = fun b => k * b * (b ^ 2 - b0 ^ 2) :=
    funext fun b => (quarticBond_hasDerivAt k b0 b).deriv
  have hh : deriv (harmonicBond (2 * k * b0 ^ 2) b0) = fun b => 2 * k * b0 ^ 2 * (b - b0) :=
    funext fun b => (harmonicBond_hasDerivAt _ b0 b).deriv
  constructor
  · rw [hq]
    have h := ((hasDerivAt_id b0).const_mul k).mul ((hasDerivAt_pow 2 b0).sub_const (b0 ^ 2))
    have h' : HasDerivAt (fun b => k * b * (b ^ 2 - b0 ^ 2)) (2 * k * b0 ^ 2) b0 := by
      refine (h.congr_deriv ?_)
      simp only [id]; push_cast; ring
    exact h'.deriv
  · rw [hh]
    have h := ((hasDerivAt_id b0).sub_const b0).const_mul (2 * k * b0 ^ 2)
    exact (h.congr_deriv (by simp)).deriv

/-! ## Bond angles -/

/-- `cos θ` of the angle between `r_ij` and `r_kj`. -/
noncomputable def cosAngle (rij rkj : V3) : ℝ := dot rij rkj / (len rij * len rkj)

/-- The cosine-harmonic angle energy at `c = cos θ`: `K/2 · (c − cos θ0)²`. -/
noncomputable def cosineAngle (k cos0 c : ℝ) : ℝ := k / 2 * (c - cos0) ^ 2

/-- The harmonic angle energy at `θ`: `K/2 · (θ − θ0)²`. -/
noncomputable def harmonicAngle (k θ0 θ : ℝ) : ℝ := k / 2 * (θ - θ0) ^ 2

/-- The forces `(f_i, f_j, f_k)` of a three-atom term from `f_i` and `f_k`: `f_j = −(f_i + f_k)`. -/
def threeBody (fi fk : V3) : V3 × V3 × V3 := (fi, -(fi + fk), fk)

/-- The cosine-angle forces, as the kernel writes them, with `df = −K (cos θ − cos θ0)`:
`f_i = df/d_ij · (r_kj/d_kj − cos θ · r_ij/d_ij)`, `f_k` the same with `i` and `k` swapped. -/
noncomputable def cosineAngleForces (k cos0 : ℝ) (rij rkj : V3) : V3 × V3 × V3 :=
  let dij := len rij
  let dkj := len rkj
  let c := cosAngle rij rkj
  let df := -k * (c - cos0)
  threeBody ((df / dij) • ((1 / dkj) • rkj - (c / dij) • rij))
    ((df / dkj) • ((1 / dij) • rij - (c / dkj) • rkj))

/-- The harmonic-angle forces, as the kernel writes them away from `sin θ = 0`, with
`ff = K (θ − θ0)/sin θ`: `f_i = ff/d_ij · (r_kj/d_kj − cos θ · r_ij/d_ij)`. -/
noncomputable def harmonicAngleForces (k θ0 : ℝ) (rij rkj : V3) : V3 × V3 × V3 :=
  let dij := len rij
  let dkj := len rkj
  let c := cosAngle rij rkj
  let θ := Real.arccos c
  let ff := k * (θ - θ0) / Real.sin θ
  threeBody ((ff / dij) • ((1 / dkj) • rkj - (c / dij) • rij))
    ((ff / dkj) • ((1 / dij) • rij - (c / dkj) • rkj))

/-- The torque of three-atom forces about any origin, with `x_i = x_j + r_ij`,
`x_k = x_j + r_kj`, is `r_ij × f_i + r_kj × f_k`. -/
theorem threeBody_torque (xj rij rkj fi fk : V3) :
    let f := threeBody fi fk
    cross (xj + rij) f.1 + cross xj f.2.1 + cross (xj + rkj) f.2.2 =
      cross rij fi + cross rkj fk := by
  ext <;> simp [threeBody] <;> ring

/-- Three-atom forces sum to zero. -/
theorem threeBody_sum_zero (fi fk : V3) :
    let f := threeBody fi fk
    f.1 + f.2.1 + f.2.2 = 0 := by
  ext <;> simp [threeBody]

/-- **The cosine-angle forces sum to zero.** -/
theorem cosineAngle_forces_sum_zero (k cos0 : ℝ) (rij rkj : V3) :
    let f := cosineAngleForces k cos0 rij rkj
    f.1 + f.2.1 + f.2.2 = 0 :=
  threeBody_sum_zero _ _

/-- **The cosine-angle forces exert no torque**, about any origin `x_j` with `x_i = x_j + r_ij`,
`x_k = x_j + r_kj`. -/
theorem cosineAngle_torque_zero (k cos0 : ℝ) (xj rij rkj : V3) :
    let f := cosineAngleForces k cos0 rij rkj
    cross (xj + rij) f.1 + cross xj f.2.1 + cross (xj + rkj) f.2.2 = 0 := by
  simp only [cosineAngleForces]
  rw [threeBody_torque]
  ext <;> simp <;> ring

/-- The harmonic-angle forces sum to zero and exert no torque. -/
theorem harmonicAngle_conservation (k θ0 : ℝ) (xj rij rkj : V3) :
    let f := harmonicAngleForces k θ0 rij rkj
    f.1 + f.2.1 + f.2.2 = 0 ∧
      cross (xj + rij) f.1 + cross xj f.2.1 + cross (xj + rkj) f.2.2 = 0 := by
  refine ⟨threeBody_sum_zero _ _, ?_⟩
  simp only [harmonicAngleForces]
  rw [threeBody_torque]
  ext <;> simp <;> ring

/-! ## Proper and improper dihedrals -/

/-- `r_im`, the part of `r_ij` perpendicular to `r_kj`. -/
noncomputable def perpI (rij rkj : V3) : V3 := rij - (dot rij rkj / dot rkj rkj) • rkj

/-- `r_ln`, the part of `−r_kl = x_l − x_k` perpendicular to `r_kj`. -/
noncomputable def perpL (rkj rkl : V3) : V3 := (dot rkl rkj / dot rkj rkj) • rkj - rkl

/-- `cos φ` of the dihedral `i j k l`: the cosine of the angle between `r_im` and `r_ln`. -/
noncomputable def cosDihedral (rij rkj rkl : V3) : ℝ := cosAngle (perpI rij rkj) (perpL rkj rkl)

/-- The dihedral angle `φ ∈ [−π, π]`, negative when `r_ij · (r_kj × r_kl) < 0`. -/
noncomputable def signedDihedral (rij rkj rkl : V3) : ℝ :=
  if dot rij (cross rkj rkl) < 0 then -Real.arccos (cosDihedral rij rkj rkl)
  else Real.arccos (cosDihedral rij rkj rkl)

/-- The improper dihedral angle `ζ ∈ [−π, π]`: the angle between the plane normals
`r_mj = r_ij × r_kj` and `r_nk = r_kj × r_kl`, negative when `r_ij · r_nk < 0`. -/
noncomputable def improperAngle (rij rkj rkl : V3) : ℝ :=
  if dot rij (cross rkj rkl) < 0 then -Real.arccos (cosAngle (cross rij rkj) (cross rkj rkl))
  else Real.arccos (cosAngle (cross rij rkj) (cross rkj rkl))

/-- `cos mφ` as the polynomial in `c = cos φ` the code evaluates: the Chebyshev polynomial
`T_m(c)` for `1 ≤ m ≤ 6` (`cosMultiple_eq_cos`). For `m = 0` the code (newen and gromosXX
`dihedral_interaction.cc`) takes `0`, not `T_0 = 1`; larger `m` are not supported. -/
noncomputable def cosMultiple : ℕ → ℝ → ℝ
  | 1, c => c
  | 2, c => 2 * c ^ 2 - 1
  | 3, c => 4 * c ^ 3 - 3 * c
  | 4, c => 8 * c ^ 4 - 8 * c ^ 2 + 1
  | 5, c => 16 * c ^ 5 - 20 * c ^ 3 + 5 * c
  | 6, c => 32 * c ^ 6 - 48 * c ^ 4 + 18 * c ^ 2 - 1
  | _, _ => 0

/-- The code's `d cos mφ / d cos φ`, the derivative of `cosMultiple m` (`dcosMultiple_is_deriv`). -/
noncomputable def dcosMultiple : ℕ → ℝ → ℝ
  | 1, _ => 1
  | 2, c => 4 * c
  | 3, c => 12 * c ^ 2 - 3
  | 4, c => 32 * c ^ 3 - 16 * c
  | 5, c => 80 * c ^ 4 - 60 * c ^ 2 + 5
  | 6, c => 192 * c ^ 5 - 192 * c ^ 3 + 36 * c
  | _, _ => 0

/-- The dihedral energy of COVALENTFORM NTBDN = 1 at `c = cos φ`: `K (1 + cos δ · cos mφ)`. -/
noncomputable def dihedral (k cosδ : ℝ) (m : ℕ) (c : ℝ) : ℝ := k * (1 + cosδ * cosMultiple m c)

/-- The dihedral energy of NTBDN = 0 at `φ`: `K (1 + cos(mφ − δ))`. -/
noncomputable def dihedralNew (k δ m φ : ℝ) : ℝ := k * (1 + Real.cos (m * φ - δ))

/-- The improper dihedral energy at `ζ`: `K/2 · (ζ − ζ0)²`. -/
noncomputable def improper (k ζ0 ζ : ℝ) : ℝ := k / 2 * (ζ - ζ0) ^ 2

/-- `cos((n + 2)φ) = 2 cos φ cos((n + 1)φ) − cos(nφ)`, the Chebyshev recursion. -/
theorem cos_add_two_mul (n φ : ℝ) :
    Real.cos ((n + 2) * φ) = 2 * Real.cos φ * Real.cos ((n + 1) * φ) - Real.cos (n * φ) := by
  have h1 : (n + 2) * φ = (n + 1) * φ + φ := by ring
  have h2 : n * φ = (n + 1) * φ - φ := by ring
  rw [h1, h2, Real.cos_add, Real.cos_sub]
  ring

/-- **The code's polynomials are `cos mφ`**: for `1 ≤ m ≤ 6`,
`cosMultiple m (cos φ) = cos(mφ)`. -/
theorem cosMultiple_eq_cos {m : ℕ} (h1 : 1 ≤ m) (h6 : m ≤ 6) (φ : ℝ) :
    cosMultiple m (Real.cos φ) = Real.cos (m * φ) := by
  have c2 := Real.cos_two_mul φ
  have c3 := Real.cos_three_mul φ
  have c4 : Real.cos (4 * φ) = 2 * Real.cos φ * Real.cos (3 * φ) - Real.cos (2 * φ) := by
    have h := cos_add_two_mul 2 φ; norm_num at h; exact h
  have c5 : Real.cos (5 * φ) = 2 * Real.cos φ * Real.cos (4 * φ) - Real.cos (3 * φ) := by
    have h := cos_add_two_mul 3 φ; norm_num at h; exact h
  have c6 : Real.cos (6 * φ) = 2 * Real.cos φ * Real.cos (5 * φ) - Real.cos (4 * φ) := by
    have h := cos_add_two_mul 4 φ; norm_num at h; exact h
  interval_cases m
  · simp [cosMultiple]
  · simp only [cosMultiple]; push_cast; rw [c2]
  · simp only [cosMultiple]; push_cast; rw [c3]
  · simp only [cosMultiple]; push_cast; rw [c4, c3, c2]; ring
  · simp only [cosMultiple]; push_cast; rw [c5, c4, c3, c2]; ring
  · simp only [cosMultiple]; push_cast; rw [c6, c5, c4, c3, c2]; ring

/-- `a · cⁿ` has derivative `a · n · cⁿ⁻¹`. -/
private theorem hasDerivAt_monomial (a : ℝ) (n : ℕ) (c : ℝ) :
    HasDerivAt (fun x : ℝ => a * x ^ n) (a * (n * c ^ (n - 1))) c :=
  (hasDerivAt_pow n c).const_mul a

/-- **The code's derivative polynomial is the derivative of `cosMultiple m`**, for every `m`. -/
theorem dcosMultiple_is_deriv (m : ℕ) (c : ℝ) :
    HasDerivAt (cosMultiple m) (dcosMultiple m c) c := by
  match m with
  | 0 => exact hasDerivAt_const c (0 : ℝ)
  | 1 => exact hasDerivAt_id c
  | 2 =>
    have h := (hasDerivAt_monomial 2 2 c).sub_const 1
    exact h.congr_deriv (by simp [dcosMultiple]; ring)
  | 3 =>
    have h := (hasDerivAt_monomial 4 3 c).sub (hasDerivAt_monomial 3 1 c)
    refine (h.congr_deriv (by simp [dcosMultiple]; ring)).congr_of_eventuallyEq ?_
    exact Filter.Eventually.of_forall fun x => by simp [cosMultiple]
  | 4 =>
    have h := ((hasDerivAt_monomial 8 4 c).sub (hasDerivAt_monomial 8 2 c)).add_const 1
    exact h.congr_deriv (by simp [dcosMultiple]; ring)
  | 5 =>
    have h := ((hasDerivAt_monomial 16 5 c).sub (hasDerivAt_monomial 20 3 c)).add
      (hasDerivAt_monomial 5 1 c)
    refine (h.congr_deriv (by simp [dcosMultiple]; ring)).congr_of_eventuallyEq ?_
    exact Filter.Eventually.of_forall fun x => by simp [cosMultiple]
  | 6 =>
    have h := (((hasDerivAt_monomial 32 6 c).sub (hasDerivAt_monomial 48 4 c)).add
      (hasDerivAt_monomial 18 2 c)).sub_const 1
    exact h.congr_deriv (by simp [dcosMultiple]; ring)
  | n + 7 => exact hasDerivAt_const c (0 : ℝ)

/-- With `δ = 0` or `δ = π`, `cos(mφ − δ) = cos δ · cos(mφ)`. -/
theorem dihedral_new_eq_dihedral_at_phase_0_or_pi {δ : ℝ} (hδ : δ = 0 ∨ δ = Real.pi)
    (k m φ : ℝ) : k * (1 + Real.cos (m * φ - δ)) = k * (1 + Real.cos δ * Real.cos (m * φ)) := by
  rcases hδ with h | h <;> simp [h, Real.cos_sub_pi]

/-- **Why the two dihedral forms agree on every GROMOS topology**: with `δ ∈ {0, π}` and
`1 ≤ m ≤ 6`, `K (1 + cos(mφ − δ)) = K (1 + cos δ · cosMultiple m (cos φ))`. -/
theorem dihedral_new_eq_dihedral {δ : ℝ} (hδ : δ = 0 ∨ δ = Real.pi) {m : ℕ} (h1 : 1 ≤ m)
    (h6 : m ≤ 6) (k φ : ℝ) :
    dihedralNew k δ m φ = dihedral k (Real.cos δ) m (Real.cos φ) := by
  rw [dihedralNew, dihedral, cosMultiple_eq_cos h1 h6,
    dihedral_new_eq_dihedral_at_phase_0_or_pi hδ]

/-- The forces `(f_i, f_j, f_k, f_l)` of a four-atom term from `f_i` and `f_l`, as both dihedral
kernels and the improper kernel complete them: `f_j = (r_ij·r_kj/|r_kj|² − 1) f_i −
(r_kl·r_kj/|r_kj|²) f_l`, `f_k = −(f_i + f_j + f_l)`. -/
noncomputable def fourBody (rij rkj rkl fi fl : V3) : V3 × V3 × V3 × V3 :=
  let fj := (dot rij rkj / dot rkj rkj - 1) • fi - (dot rkl rkj / dot rkj rkj) • fl
  (fi, fj, -(fi + fj + fl), fl)

/-- The forces of the Chebyshev dihedral (NTBDN = 1) as the kernel writes them:
`k_i = −K cos δ · T'_m(cos φ)/d_im`, `f_i = k_i (r_ln/d_ln − cos φ · r_im/d_im)`, `f_l` the same
with `i` and `l` swapped. -/
noncomputable def chebyshevForces (k cosδ : ℝ) (m : ℕ) (rij rkj rkl : V3) :
    V3 × V3 × V3 × V3 :=
  let rim := perpI rij rkj
  let rln := perpL rkj rkl
  let dim := len rim
  let dln := len rln
  let c := cosAngle rim rln
  let ki := -k * cosδ * dcosMultiple m c / dim
  let kl := -k * cosδ * dcosMultiple m c / dln
  fourBody rij rkj rkl ((ki / dln) • rln - (ki * c / dim) • rim)
    ((kl / dim) • rim - (kl * c / dln) • rln)

/-- The forces of the dihedral of NTBDN = 0 as the kernel writes them (away from a linear bond
angle): `k_i = K m sin(mφ − δ)`, `f_i = k_i |r_kj|/|r_mj|² · r_mj`,
`f_l = −k_i |r_kj|/|r_nk|² · r_nk`. -/
noncomputable def newDihedralForces (k δ m : ℝ) (rij rkj rkl : V3) : V3 × V3 × V3 × V3 :=
  let rmj := cross rij rkj
  let rnk := cross rkj rkl
  let ki := k * m * Real.sin (m * signedDihedral rij rkj rkl - δ)
  fourBody rij rkj rkl ((ki * len rkj / dot rmj rmj) • rmj) ((-ki * len rkj / dot rnk rnk) • rnk)

/-- The forces of the improper dihedral at angle `ζ` (the kernel's `ζ`, moved into
`[ζ0 − π, ζ0 + π]`) as the kernel writes them: `a = −K (ζ − ζ0) |r_kj|`, `f_i = a/|r_mj|² · r_mj`,
`f_l = −a/|r_nk|² · r_nk`. -/
noncomputable def improperForces (k ζ0 ζ : ℝ) (rij rkj rkl : V3) : V3 × V3 × V3 × V3 :=
  let rmj := cross rij rkj
  let rnk := cross rkj rkl
  let a := -k * (ζ - ζ0) * len rkj
  fourBody rij rkj rkl ((a / dot rmj rmj) • rmj) ((-a / dot rnk rnk) • rnk)

/-- Four-atom forces sum to zero. -/
theorem fourBody_sum_zero (rij rkj rkl fi fl : V3) :
    let f := fourBody rij rkj rkl fi fl
    f.1 + f.2.1 + f.2.2.1 + f.2.2.2 = 0 := by
  ext <;> simp [fourBody] <;> ring

/-- The torque of four-atom forces about any origin, with `x_i = x_j + r_ij`, `x_k = x_j + r_kj`,
`x_l = x_k − r_kl`, is `r_im × f_i + r_ln × f_l`. -/
theorem fourBody_torque (xj rij rkj rkl fi fl : V3) :
    let f := fourBody rij rkj rkl fi fl
    cross (xj + rij) f.1 + cross xj f.2.1 + cross (xj + rkj) f.2.2.1 +
        cross (xj + rkj - rkl) f.2.2.2 =
      cross (perpI rij rkj) fi + cross (perpL rkj rkl) fl := by
  ext <;> simp [fourBody, perpI, perpL] <;> ring

/-- Forces along the plane normals, `f_i = b/|r_mj|² · r_mj` and `f_l = −b/|r_nk|² · r_nk`, exert
no torque: `r_im × r_mj = −(|r_mj|²/|r_kj|²) r_kj` and `r_ln × r_nk = −(|r_nk|²/|r_kj|²) r_kj`
(Lagrange's identity). -/
theorem normal_forces_torque_zero {rij rkj rkl : V3} (hkj : dot rkj rkj ≠ 0)
    (hmj : dot (cross rij rkj) (cross rij rkj) ≠ 0) (hnk : dot (cross rkj rkl) (cross rkj rkl) ≠ 0)
    (b : ℝ) :
    cross (perpI rij rkj) ((b / dot (cross rij rkj) (cross rij rkj)) • cross rij rkj) +
        cross (perpL rkj rkl) ((-b / dot (cross rkj rkl) (cross rkj rkl)) • cross rkj rkl) = 0 := by
  -- `|r_kj|²` as an atom `D` for `field_simp`, replaced by its components afterwards
  have hm : cross (perpI rij rkj) (cross rij rkj) =
      (-(dot (cross rij rkj) (cross rij rkj) / dot rkj rkj)) • rkj := by
    have key : ∀ D : ℝ, D = dot rkj rkj → D ≠ 0 →
        cross (rij - (dot rij rkj / D) • rkj) (cross rij rkj) =
          (-(dot (cross rij rkj) (cross rij rkj) / D)) • rkj := by
      intro D hD hD0
      ext <;> v3_expand <;> field_simp <;> (simp only [dot] at hD; subst hD; ring)
    exact key _ rfl hkj
  have hn : cross (perpL rkj rkl) (cross rkj rkl) =
      (-(dot (cross rkj rkl) (cross rkj rkl) / dot rkj rkj)) • rkj := by
    have key : ∀ D : ℝ, D = dot rkj rkj → D ≠ 0 →
        cross ((dot rkl rkj / D) • rkj - rkl) (cross rkj rkl) =
          (-(dot (cross rkj rkl) (cross rkj rkl) / D)) • rkj := by
      intro D hD hD0
      ext <;> v3_expand <;> field_simp <;> (simp only [dot] at hD; subst hD; ring)
    exact key _ rfl hkj
  have hs (u v : V3) (a : ℝ) : cross u (a • v) = a • cross u v := by
    ext <;> simp <;> ring
  rw [hs, hs, hm, hn]
  ext <;> simp <;> field_simp <;> ring

/-- **The Chebyshev-dihedral forces sum to zero and exert no torque.** -/
theorem chebyshev_conservation (k cosδ : ℝ) (m : ℕ) (xj rij rkj rkl : V3) :
    let f := chebyshevForces k cosδ m rij rkj rkl
    f.1 + f.2.1 + f.2.2.1 + f.2.2.2 = 0 ∧
      cross (xj + rij) f.1 + cross xj f.2.1 + cross (xj + rkj) f.2.2.1 +
        cross (xj + rkj - rkl) f.2.2.2 = 0 := by
  refine ⟨fourBody_sum_zero _ _ _ _ _, ?_⟩
  simp only [chebyshevForces]
  rw [fourBody_torque]
  ext <;> simp <;> ring

/-- **The NTBDN = 0 dihedral forces sum to zero and exert no torque** (away from a linear bond
angle, where the kernel drops `f_i` or `f_l`). -/
theorem newDihedral_conservation (k δ m : ℝ) (xj : V3) {rij rkj rkl : V3}
    (hkj : dot rkj rkj ≠ 0) (hmj : dot (cross rij rkj) (cross rij rkj) ≠ 0)
    (hnk : dot (cross rkj rkl) (cross rkj rkl) ≠ 0) :
    let f := newDihedralForces k δ m rij rkj rkl
    f.1 + f.2.1 + f.2.2.1 + f.2.2.2 = 0 ∧
      cross (xj + rij) f.1 + cross xj f.2.1 + cross (xj + rkj) f.2.2.1 +
        cross (xj + rkj - rkl) f.2.2.2 = 0 := by
  refine ⟨fourBody_sum_zero _ _ _ _ _, ?_⟩
  simp only [newDihedralForces]
  rw [fourBody_torque, neg_mul]
  exact normal_forces_torque_zero hkj hmj hnk _

/-- **The improper-dihedral forces sum to zero and exert no torque** (away from a linear bond
angle). -/
theorem improper_conservation (k ζ0 ζ : ℝ) (xj : V3) {rij rkj rkl : V3}
    (hkj : dot rkj rkj ≠ 0) (hmj : dot (cross rij rkj) (cross rij rkj) ≠ 0)
    (hnk : dot (cross rkj rkl) (cross rkj rkl) ≠ 0) :
    let f := improperForces k ζ0 ζ rij rkj rkl
    f.1 + f.2.1 + f.2.2.1 + f.2.2.2 = 0 ∧
      cross (xj + rij) f.1 + cross xj f.2.1 + cross (xj + rkj) f.2.2.1 +
        cross (xj + rkj - rkl) f.2.2.2 = 0 := by
  refine ⟨fourBody_sum_zero _ _ _ _ _, ?_⟩
  simp only [improperForces]
  rw [fourBody_torque]
  exact normal_forces_torque_zero hkj hmj hnk _

/-! ## The cosine-angle force is minus the gradient of its energy -/

theorem line_zero (u a : V3) : u + (0 : ℝ) • a = u := by ext <;> simp

/-- `|u + t a|` has the derivative `u · a/|u|` at `t = 0`, for `u ≠ 0`. -/
theorem hasDerivAt_len_line {u : V3} (hu : u ≠ 0) (a : V3) :
    HasDerivAt (fun t : ℝ => len (u + t • a)) (dot u a / len u) 0 := by
  have hq : HasDerivAt (fun t : ℝ => dot (u + t • a) (u + t • a)) (2 * dot u a) 0 := by
    have h1 := ((hasDerivAt_id (0 : ℝ)).const_mul (2 * dot u a)).add
      (((hasDerivAt_id (0 : ℝ)).pow 2).const_mul (dot a a))
    have hfun : (fun t : ℝ => dot (u + t • a) (u + t • a)) =
        fun t => dot u u + ((2 * dot u a) * id t + dot a a * id t ^ 2) := by
      funext t; v3_expand; simp only [id]; ring
    rw [hfun]
    exact (h1.const_add (dot u u)).congr_deriv (by simp)
  have h := hq.sqrt (by rw [line_zero]; exact (dot_self_pos hu).ne')
  rw [line_zero] at h
  refine h.congr_deriv ?_
  have hl : len u ≠ 0 := (len_pos hu).ne'
  rw [show Real.sqrt (dot u u) = len u from rfl]
  field_simp

/-- `(u + t a) · (v + t b)` has the derivative `a · v + u · b` at `t = 0`. -/
theorem hasDerivAt_dot_lines (u v a b : V3) :
    HasDerivAt (fun t : ℝ => dot (u + t • a) (v + t • b)) (dot a v + dot u b) 0 := by
  have h1 := ((hasDerivAt_id (0 : ℝ)).const_mul (dot a v + dot u b)).add
    (((hasDerivAt_id (0 : ℝ)).pow 2).const_mul (dot a b))
  have hfun : (fun t : ℝ => dot (u + t • a) (v + t • b)) =
      fun t => dot u v + ((dot a v + dot u b) * id t + dot a b * id t ^ 2) := by
    funext t; v3_expand; simp only [id]; ring
  rw [hfun]
  exact (h1.const_add (dot u v)).congr_deriv (by simp)

/-- The derivative of the cosine-angle energy when `r_ij` moves along `a` and `r_kj` along `b`. -/
theorem cosineAngle_hasDerivAt_lines (k cos0 : ℝ) {u v : V3} (hu : u ≠ 0) (hv : v ≠ 0)
    (a b : V3) :
    HasDerivAt (fun t : ℝ => cosineAngle k cos0 (cosAngle (u + t • a) (v + t • b)))
      (k / 2 * (2 * (cosAngle u v - cos0)) *
        (((dot a v + dot u b) * (len u * len v) -
            dot u v * (dot u a / len u * len v + len u * (dot v b / len v))) /
          (len u * len v) ^ 2)) 0 := by
  have hd : len u * len v ≠ 0 := mul_ne_zero (len_pos hu).ne' (len_pos hv).ne'
  have hden := (hasDerivAt_len_line hu a).mul (hasDerivAt_len_line hv b)
  have hg := HasDerivAt.div (hasDerivAt_dot_lines u v a b) hden
    (by simp only [Pi.mul_apply, line_zero]; exact hd)
  simp only [Pi.mul_apply, line_zero] at hg
  have hE := (((hg.sub_const cos0).pow 2).const_mul (k / 2))
  refine hE.congr_deriv ?_
  simp only [Pi.div_apply, Pi.mul_apply, line_zero, cosAngle]
  push_cast
  ring

/-- **The cosine-angle forces are `−∇V`.** Moving atom `i` along `e` (`r_ij + t e`), atom `k`
(`r_kj + t e`) or atom `j` (`r_ij − t e`, `r_kj − t e`), the energy changes at the rate
`−f_i · e`, `−f_k · e`, `−f_j · e`. -/
theorem cosineAngle_force_is_minus_grad (k cos0 : ℝ) {rij rkj : V3} (hij : rij ≠ 0)
    (hkj : rkj ≠ 0) (e : V3) :
    let f := cosineAngleForces k cos0 rij rkj
    HasDerivAt (fun t : ℝ => cosineAngle k cos0 (cosAngle (rij + t • e) rkj)) (-dot f.1 e) 0 ∧
      HasDerivAt (fun t : ℝ => cosineAngle k cos0 (cosAngle rij (rkj + t • e)))
        (-dot f.2.2 e) 0 ∧
      HasDerivAt (fun t : ℝ => cosineAngle k cos0 (cosAngle (rij - t • e) (rkj - t • e)))
        (-dot f.2.1 e) 0 := by
  have hi := (len_pos hij).ne'
  have hk := (len_pos hkj).ne'
  refine ⟨?_, ?_, ?_⟩
  · have h := cosineAngle_hasDerivAt_lines k cos0 hij hkj e 0
    refine (h.congr_of_eventuallyEq (Filter.Eventually.of_forall fun t => ?_)).congr_deriv ?_
    · simp only
      congr 2; ext <;> simp
    · simp only [cosineAngleForces, threeBody, cosAngle]
      v3_expand
      field_simp
      ring
  · have h := cosineAngle_hasDerivAt_lines k cos0 hij hkj 0 e
    refine (h.congr_of_eventuallyEq (Filter.Eventually.of_forall fun t => ?_)).congr_deriv ?_
    · simp only
      congr 2; ext <;> simp
    · simp only [cosineAngleForces, threeBody, cosAngle]
      v3_expand
      field_simp
      ring
  · have h := cosineAngle_hasDerivAt_lines k cos0 hij hkj (-e) (-e)
    refine (h.congr_of_eventuallyEq (Filter.Eventually.of_forall fun t => ?_)).congr_deriv ?_
    · simp only
      congr 2 <;> ext <;> simp <;> ring
    · simp only [cosineAngleForces, threeBody, cosAngle]
      v3_expand
      field_simp
      ring

end NewenTheory.Interactions.Bonded

end
