import NewenSpec.Float
import NewenTheory.Interactions.Bonded

/-!
# Model of the bonded kernels of `crates/newen-forces`

Rust (spec keys in `spec/SPEC.toml`), all in `src/bonded/`:

* `bonded.quartic_bond`, `bonded.harmonic_bond`: `quartic_bond`, `harmonic_bond` in `bonds.rs`;
* `bonded.cosine_angle`, `bonded.harmonic_angle`: `cosine_angle`, `harmonic_angle` in `angles.rs`;
* `bonded.dihedral_chebyshev`, `bonded.dihedral_new`: `dihedral_chebyshev`, `dihedral_new` in
  `dihedrals.rs` (COVALENTFORM NTBDN = 1 and 0);
* `bonded.improper_dihedral`: `improper_dihedral` in `improper.rs`.

Theory: `NewenTheory/Interactions/Bonded.lean`.

Each kernel is written once, over any type `α` with the operations the Rust uses (`+ − × ÷`,
negation, `sqrt`, `acos`, `cos`, `sin` from `NewenSpec/Float.lean`, `<`, the small integer
literals), operation by operation in the Rust order, on a `Vec3 α` whose `dot`, `cross`, `length`
and scalar `*`, `/` are glam's (`x·x' + y·y' + z·z'` left to right, no fused multiply-add). It is
used twice, as in `NewenSpec/Nonbonded.lean`:

1. at `ℝ`, where the `*_real_eq` theorems prove it equal to the theory's energy and force
   expressions outside the guards;
2. at `Float`, where `lake exe spec-vectors` evaluates it and
   `crates/newen-forces/tests/bonded_spec.rs` compares the Rust kernels bit for bit.

The literals that are not small integers are parameters, given at `Float` by their bits: the
guards `1e-10` (`ZERO_DISTANCE`, `ZERO_DISTANCE2`, `ZERO_SINE`, `NEAR_LINEAR`, all the same
double), `0.25`, `0.5` and `π`.

What is not proved: the rounding error between the two instantiations.
-/

namespace NewenSpec.Bonded

open NewenTheory.Interactions.Bonded

/-- A vector as glam's `DVec3`. -/
structure Vec3 (α : Type) where
  x : α
  y : α
  z : α

namespace Vec3

variable {α : Type}

instance [Add α] : Add (Vec3 α) := ⟨fun u v => ⟨u.x + v.x, u.y + v.y, u.z + v.z⟩⟩
instance [Sub α] : Sub (Vec3 α) := ⟨fun u v => ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩⟩
instance [Neg α] : Neg (Vec3 α) := ⟨fun u => ⟨-u.x, -u.y, -u.z⟩⟩
/-- `u * s`: every component times the scalar. -/
instance [Mul α] : HMul (Vec3 α) α (Vec3 α) := ⟨fun u s => ⟨u.x * s, u.y * s, u.z * s⟩⟩
/-- `u / s`: every component divided by the scalar. -/
instance [Div α] : HDiv (Vec3 α) α (Vec3 α) := ⟨fun u s => ⟨u.x / s, u.y / s, u.z / s⟩⟩

/-- `Vec3::ZERO`. -/
def zero [OfNat α 0] : Vec3 α := ⟨0, 0, 0⟩

/-- `u.dot(v)`: `(x·x' + y·y') + z·z'`. -/
def dot [Add α] [Mul α] (u v : Vec3 α) : α := u.x * v.x + u.y * v.y + u.z * v.z

/-- `u.cross(v)`, glam's component formulas. -/
def cross [Sub α] [Mul α] (u v : Vec3 α) : Vec3 α :=
  ⟨u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y⟩

/-- `u.length()`: `sqrt(u · u)`. -/
def length [Add α] [Mul α] [Sqrt α] (u : Vec3 α) : α := Sqrt.sqrt (dot u u)

end Vec3

/-! ## 1. Model (mirrors the Rust functions) -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α]
  [DecidableRel (α := α) (· < ·)] [Sqrt α] [Acos α] [Cos α] [Sin α]
  [OfNat α 0] [OfNat α 1] [OfNat α 2] [OfNat α 3] [OfNat α 4] [OfNat α 5] [OfNat α 8]
  [OfNat α 12] [OfNat α 16] [OfNat α 18] [OfNat α 20] [OfNat α 32] [OfNat α 36] [OfNat α 48]
  [OfNat α 60] [OfNat α 80] [OfNat α 192]

open Vec3

/-- `x.clamp(-1.0, 1.0)`: `if x < −1 { −1 }`, then `if x > 1 { 1 }` (a NaN passes through). -/
def clampUnit (x : α) : α :=
  let x := if x < -1 then -1 else x
  if 1 < x then 1 else x

/-- The comparison `x.abs() < t`, with `|x|` as `if x < 0 { −x } else { x }`: the same truth value
as C `fabs` for every input, `−0` and NaN included. -/
def absLt (x t : α) : Prop := (if x < 0 then -x else x) < t

instance (x t : α) : Decidable (absLt α x t) := by unfold absLt; infer_instance

/-- `quartic_bond(r, k, b0)`: `(energy, f_i)`, with the literal `0.25` as `quarter`. -/
def quarticBondModel (quarter : α) (r : Vec3 α) (k b0 : α) : α × Vec3 α :=
  let r2 := r.dot r
  let dr2 := r2 - b0 * b0
  (quarter * k * dr2 * dr2, r * ((-k) * dr2))

/-- `harmonic_bond(r, k, b0)`: `(energy, f_i)`, with `ZERO_DISTANCE` as `guard` and `0.5` as
`half`. -/
def harmonicBondModel (guard half : α) (r : Vec3 α) (k b0 : α) : α × Vec3 α :=
  let b := r.length
  if b < guard then (0, zero)
  else
    let dr := b - b0
    (half * k * dr * dr, r * ((-k) * dr / b))

/-- `cosine_angle(r_ij, r_kj, k, cos0)`: `(energy, f_i, f_j, f_k)`. -/
def cosineAngleModel (guard half : α) (rij rkj : Vec3 α) (k cos0 : α) :
    α × Vec3 α × Vec3 α × Vec3 α :=
  let dij := rij.length
  let dkj := rkj.length
  if dij < guard ∨ dkj < guard then (0, zero, zero, zero)
  else
    let cost := clampUnit α (rij.dot rkj / (dij * dkj))
    let dCos := cost - cos0
    let energy := half * k * dCos * dCos
    let df := (-k) * dCos
    let fi := (rkj * (1 / dkj) - rij * (cost / dij)) * (df / dij)
    let fk := (rij * (1 / dij) - rkj * (cost / dkj)) * (df / dkj)
    (energy, fi, -(fi + fk), fk)

/-- `harmonic_angle(r_ij, r_kj, k, theta0)`: `(energy, f_i, f_j, f_k)`, with `ZERO_DISTANCE` as
`guard`, `ZERO_SINE` as `sineGuard` and `π` as `pi`. -/
def harmonicAngleModel (guard sineGuard half pi : α) (rij rkj : Vec3 α) (k theta0 : α) :
    α × Vec3 α × Vec3 α × Vec3 α :=
  let nothing : α × Vec3 α × Vec3 α × Vec3 α := (0, zero, zero, zero)
  let dij := rij.length
  let dkj := rkj.length
  if dij < guard ∨ dkj < guard then nothing
  else
    let cosTheta := clampUnit α (rij.dot rkj / (dij * dkj))
    let theta := Acos.acos cosTheta
    let sinTheta := Sin.sin theta
    let forces (ki kk energy : α) : α × Vec3 α × Vec3 α × Vec3 α :=
      let unitIj := rij / dij
      let unitKj := rkj / dkj
      let fi := (unitKj - unitIj * cosTheta) * ki
      let fk := (unitIj - unitKj * cosTheta) * kk
      (energy, fi, -(fi + fk), fk)
    if absLt α sinTheta sineGuard then
      if pi + sineGuard < theta0 ∨ theta0 < pi - sineGuard then nothing
      else forces ((-k) / dij) ((-k) / dkj) (k * (1 + cosTheta))
    else
      let ff := k * (theta - theta0) / sinTheta
      let dTheta := theta - theta0
      forces (ff / dij) (ff / dkj) (half * k * dTheta * dTheta)

/-- `cos_multiple(m, c)`: `(cos mφ, d cos mφ / d cos φ)` at `c = cos φ`, `none` for `m > 6`. -/
def cosMultipleModel (m : Nat) (c : α) : Option (α × α) :=
  let c2 := c * c
  let c3 := c2 * c
  let c4 := c3 * c
  match m with
  | 0 => some (0, 0)
  | 1 => some (c, 1)
  | 2 => some (2 * c2 - 1, 4 * c)
  | 3 => some (4 * c3 - 3 * c, 12 * c2 - 3)
  | 4 => some (8 * c4 - 8 * c2 + 1, 32 * c3 - 16 * c)
  | 5 => some (16 * c4 * c - 20 * c3 + 5 * c, 80 * c4 - 60 * c2 + 5)
  | 6 => some (32 * c4 * c2 - 48 * c4 + 18 * c2 - 1, 192 * c4 * c - 192 * c3 + 36 * c)
  | _ => none

/-- `dihedral_chebyshev(r_ij, r_kj, r_kl, k, cos_delta, m)`: `(energy, f_i, f_j, f_k, f_l)`,
with `ZERO_DISTANCE` as `guard` and `ZERO_DISTANCE2` as `guard2`. -/
def dihedralChebyshevModel (guard guard2 : α) (rij rkj rkl : Vec3 α) (k cosDelta : α)
    (m : Nat) : α × Vec3 α × Vec3 α × Vec3 α × Vec3 α :=
  let nothing : α × Vec3 α × Vec3 α × Vec3 α × Vec3 α := (0, zero, zero, zero, zero)
  let dkj2 := rkj.dot rkj
  if dkj2 < guard2 then nothing
  else
    let frim := rij.dot rkj / dkj2
    let frln := rkl.dot rkj / dkj2
    let rim := rij - rkj * frim
    let rln := rkj * frln - rkl
    let dim := rim.length
    let dln := rln.length
    if dim < guard ∨ dln < guard then nothing
    else
      let cosPhi := clampUnit α (rim.dot rln / (dim * dln))
      match cosMultipleModel α m cosPhi with
      | none => nothing
      | some (cosMPhi, dCosMPhi) =>
        let energy := k * (1 + cosDelta * cosMPhi)
        let ki := (-k) * cosDelta * dCosMPhi / dim
        let kl := (-k) * cosDelta * dCosMPhi / dln
        let fi := rln * (ki / dln) - rim * (ki * cosPhi / dim)
        let fl := rim * (kl / dim) - rln * (kl * cosPhi / dln)
        let fj := fi * (frim - 1) - fl * frln
        (energy, fi, fj, -(fi + fj + fl), fl)

/-- `dihedral_new(r_ij, r_kj, r_kl, k, delta, m)`: `(energy, f_i, f_j, f_k, f_l)`, with
`ZERO_DISTANCE` as `guard` and `NEAR_LINEAR` as `nearLinear`. -/
def dihedralNewModel (guard nearLinear : α) (rij rkj rkl : Vec3 α) (k delta m : α) :
    α × Vec3 α × Vec3 α × Vec3 α × Vec3 α :=
  let nothing : α × Vec3 α × Vec3 α × Vec3 α × Vec3 α := (0, zero, zero, zero, zero)
  let rmj := rij.cross rkj
  let rnk := rkj.cross rkl
  let dkj2 := rkj.dot rkj
  let dmj2 := rmj.dot rmj
  let dnk2 := rnk.dot rnk
  let dkj := Sqrt.sqrt dkj2
  if dkj < guard then nothing
  else
    let frim := rij.dot rkj / dkj2
    let frln := rkl.dot rkj / dkj2
    let rim := rij - rkj * frim
    let rln := rkj * frln - rkl
    let dim := rim.length
    let dln := rln.length
    if dim < guard ∨ dln < guard then nothing
    else
      let cosPhi := clampUnit α (rim.dot rln / (dim * dln))
      let phi := Acos.acos cosPhi
      let phi := if rij.dot rnk < 0 then -phi else phi
      let x := m * phi - delta
      let energy := k * (1 + Cos.cos x)
      let ki := k * m * Sin.sin x
      let kl := -ki
      let fi := if nearLinear * dkj2 < dmj2 then rmj * (ki * dkj / dmj2) else zero
      let fl := if nearLinear * dkj2 < dnk2 then rnk * (kl * dkj / dnk2) else zero
      let fj := fi * (frim - 1) - fl * frln
      (energy, fi, fj, -(fi + fj + fl), fl)

/-- `while z < lo { z += step }`, at most `fuel` times. -/
def wrapUp : Nat → α → α → α → α
  | 0, _, _, z => z
  | fuel + 1, lo, step, z => if z < lo then wrapUp fuel lo step (z + step) else z

/-- `while z > hi { z -= step }`, at most `fuel` times. -/
def wrapDown : Nat → α → α → α → α
  | 0, _, _, z => z
  | fuel + 1, hi, step, z => if hi < z then wrapDown fuel hi step (z - step) else z

/-- `improper_dihedral(r_ij, r_kj, r_kl, k, zeta0)`: `(energy, f_i, f_j, f_k, f_l)`, with
`ZERO_DISTANCE2` as `guard2`, `NEAR_LINEAR` as `nearLinear`, `0.5` as `half`, `π` as `pi`, and the
two `while` loops run at most `fuel` times each. -/
def improperModel (guard2 nearLinear half pi : α) (fuel : Nat) (rij rkj rkl : Vec3 α)
    (k zeta0 : α) : α × Vec3 α × Vec3 α × Vec3 α × Vec3 α :=
  let rmj := rij.cross rkj
  let rnk := rkj.cross rkl
  let dkj2 := rkj.dot rkj
  let dmj2 := rmj.dot rmj
  let dnk2 := rnk.dot rnk
  let dkj := Sqrt.sqrt dkj2
  let dmj := Sqrt.sqrt dmj2
  let dnk := Sqrt.sqrt dnk2
  if dmj < guard2 ∨ dnk < guard2 then (0, zero, zero, zero, zero)
  else
    let acs := clampUnit α (rmj.dot rnk / (dmj * dnk))
    let zeta := Acos.acos acs
    let zeta := if rij.dot rnk < 0 then -zeta else zeta
    let zeta := wrapUp α fuel (zeta0 - pi) (2 * pi) zeta
    let zeta := wrapDown α fuel (zeta0 + pi) (2 * pi) zeta
    let dZeta := zeta - zeta0
    let energy := half * k * dZeta * dZeta
    let ki := (-k) * dZeta * dkj
    let kl := -ki
    let ki := if dmj2 < nearLinear * dkj2 then 0 else ki / dmj2
    let kl := if dnk2 < nearLinear * dkj2 then 0 else kl / dnk2
    let kj1 := rij.dot rkj / dkj2 - 1
    let kj2 := rkl.dot rkj / dkj2
    let fi := rmj * ki
    let fl := rnk * kl
    let fj := fi * kj1 - fl * kj2
    (energy, fi, fj, -(fi + fj + fl), fl)

end Model

/-! ## 2. The `Float` instantiation (what `lake exe spec-vectors` runs)

The decimal literals are given by their bits, so that they are exactly Rust's; the Rust test
asserts the bits of `ZERO_DISTANCE`, `ZERO_DISTANCE2`, `ZERO_SINE`, `NEAR_LINEAR` and `PI`. -/

/-- `1e-10`, the value of all four guards, as `f64` bits. -/
def tinyFloat : Float := Float.ofBits 0x3ddb7cdfd9d7bdbb

/-- `0.25`. -/
def quarterFloat : Float := Float.ofBits 0x3fd0000000000000

/-- `0.5`. -/
def halfFloat : Float := Float.ofBits 0x3fe0000000000000

/-- `std::f64::consts::PI`. -/
def piFloat : Float := Float.ofBits 0x400921fb54442d18

/-- How often `improperFloat` may turn ζ by 2π in each direction; one turn suffices for
`|ζ0| ≤ π`, since `ζ ∈ [−π, π]`. -/
def wrapFuel : Nat := 64

def quarticBondFloat (r : Vec3 Float) (k b0 : Float) : Float × Vec3 Float :=
  quarticBondModel Float quarterFloat r k b0

def harmonicBondFloat (r : Vec3 Float) (k b0 : Float) : Float × Vec3 Float :=
  harmonicBondModel Float tinyFloat halfFloat r k b0

def cosineAngleFloat (rij rkj : Vec3 Float) (k cos0 : Float) :
    Float × Vec3 Float × Vec3 Float × Vec3 Float :=
  cosineAngleModel Float tinyFloat halfFloat rij rkj k cos0

def harmonicAngleFloat (rij rkj : Vec3 Float) (k theta0 : Float) :
    Float × Vec3 Float × Vec3 Float × Vec3 Float :=
  harmonicAngleModel Float tinyFloat tinyFloat halfFloat piFloat rij rkj k theta0

def dihedralChebyshevFloat (rij rkj rkl : Vec3 Float) (k cosDelta : Float) (m : Nat) :
    Float × Vec3 Float × Vec3 Float × Vec3 Float × Vec3 Float :=
  dihedralChebyshevModel Float tinyFloat tinyFloat rij rkj rkl k cosDelta m

def dihedralNewFloat (rij rkj rkl : Vec3 Float) (k delta m : Float) :
    Float × Vec3 Float × Vec3 Float × Vec3 Float × Vec3 Float :=
  dihedralNewModel Float tinyFloat tinyFloat rij rkj rkl k delta m

def improperFloat (rij rkj rkl : Vec3 Float) (k zeta0 : Float) :
    Float × Vec3 Float × Vec3 Float × Vec3 Float × Vec3 Float :=
  improperModel Float tinyFloat tinyFloat halfFloat piFloat wrapFuel rij rkj rkl k zeta0

/-! ## 3. The ℝ instantiation equals the theory

Outside the guards, every model at ℝ returns the theory's energy (`quarticBond`, `cosineAngle`,
`dihedral`, …) at the theory's geometry (`len`, `cosAngle`, `cosDihedral`, `signedDihedral`,
`improperAngle`) and the theory's force expressions (`radialForce`, `cosineAngleForces`,
`chebyshevForces`, …), about which `NewenTheory/Interactions/Bonded.lean` proves
force = −∇energy and the conservation laws. Over ℝ the clamps never act (`cos_between_mem`). -/

open NewenTheory.Interactions.Bonded.V3 (len len_sq cos_between_mem)

namespace Vec3

/-- A model vector over ℝ as a theory vector. -/
def toV3 (u : Vec3 ℝ) : V3 := ⟨u.x, u.y, u.z⟩

@[simp] theorem toV3_x (u : Vec3 ℝ) : u.toV3.x = u.x := rfl
@[simp] theorem toV3_y (u : Vec3 ℝ) : u.toV3.y = u.y := rfl
@[simp] theorem toV3_z (u : Vec3 ℝ) : u.toV3.z = u.z := rfl
@[simp] theorem add_x' (u v : Vec3 ℝ) : (u + v).x = u.x + v.x := rfl
@[simp] theorem add_y' (u v : Vec3 ℝ) : (u + v).y = u.y + v.y := rfl
@[simp] theorem add_z' (u v : Vec3 ℝ) : (u + v).z = u.z + v.z := rfl
@[simp] theorem sub_x' (u v : Vec3 ℝ) : (u - v).x = u.x - v.x := rfl
@[simp] theorem sub_y' (u v : Vec3 ℝ) : (u - v).y = u.y - v.y := rfl
@[simp] theorem sub_z' (u v : Vec3 ℝ) : (u - v).z = u.z - v.z := rfl
@[simp] theorem neg_x' (u : Vec3 ℝ) : (-u).x = -u.x := rfl
@[simp] theorem neg_y' (u : Vec3 ℝ) : (-u).y = -u.y := rfl
@[simp] theorem neg_z' (u : Vec3 ℝ) : (-u).z = -u.z := rfl
@[simp] theorem mul_x' (u : Vec3 ℝ) (s : ℝ) : (u * s).x = u.x * s := rfl
@[simp] theorem mul_y' (u : Vec3 ℝ) (s : ℝ) : (u * s).y = u.y * s := rfl
@[simp] theorem mul_z' (u : Vec3 ℝ) (s : ℝ) : (u * s).z = u.z * s := rfl
@[simp] theorem div_x' (u : Vec3 ℝ) (s : ℝ) : (u / s).x = u.x / s := rfl
@[simp] theorem div_y' (u : Vec3 ℝ) (s : ℝ) : (u / s).y = u.y / s := rfl
@[simp] theorem div_z' (u : Vec3 ℝ) (s : ℝ) : (u / s).z = u.z / s := rfl
@[simp] theorem zero_x' : (zero : Vec3 ℝ).x = 0 := rfl
@[simp] theorem zero_y' : (zero : Vec3 ℝ).y = 0 := rfl
@[simp] theorem zero_z' : (zero : Vec3 ℝ).z = 0 := rfl
@[simp] theorem cross_x' (u v : Vec3 ℝ) : (cross u v).x = u.y * v.z - v.y * u.z := rfl
@[simp] theorem cross_y' (u v : Vec3 ℝ) : (cross u v).y = u.z * v.x - v.z * u.x := rfl
@[simp] theorem cross_z' (u v : Vec3 ℝ) : (cross u v).z = u.x * v.y - v.x * u.y := rfl

theorem dot_eq (u v : Vec3 ℝ) : dot u v = V3.dot u.toV3 v.toV3 := rfl
theorem length_eq (u : Vec3 ℝ) : length u = len u.toV3 := rfl
theorem sqrt_dot_eq (u : Vec3 ℝ) : Sqrt.sqrt (dot u u) = len u.toV3 := rfl
theorem toV3_cross (u v : Vec3 ℝ) : (cross u v).toV3 = V3.cross u.toV3 v.toV3 := rfl

end Vec3

open Vec3

/-- The outputs of a two-atom model with its vector as a theory vector. -/
def out2 (p : ℝ × Vec3 ℝ) : ℝ × V3 := (p.1, p.2.toV3)

/-- The outputs of a three-atom model with its vectors as theory vectors. -/
def out3 (p : ℝ × Vec3 ℝ × Vec3 ℝ × Vec3 ℝ) : ℝ × V3 × V3 × V3 :=
  (p.1, p.2.1.toV3, p.2.2.1.toV3, p.2.2.2.toV3)

/-- The outputs of a four-atom model with its vectors as theory vectors. -/
def out4 (p : ℝ × Vec3 ℝ × Vec3 ℝ × Vec3 ℝ × Vec3 ℝ) : ℝ × V3 × V3 × V3 × V3 :=
  (p.1, p.2.1.toV3, p.2.2.1.toV3, p.2.2.2.1.toV3, p.2.2.2.2.toV3)

/-- Over ℝ the clamp is the identity on `[−1, 1]`. -/
theorem clampUnit_real {c : ℝ} (h : -1 ≤ c ∧ c ≤ 1) : clampUnit ℝ c = c := by
  simp [clampUnit, not_lt.mpr h.1, not_lt.mpr h.2]

/-- Over ℝ, `absLt x t` is `|x| < t`. -/
theorem absLt_real (x t : ℝ) : absLt ℝ x t ↔ |x| < t := by
  unfold absLt
  split_ifs with h
  · rw [abs_of_neg h]
  · rw [abs_of_nonneg (not_lt.mp h)]

/-- For `m ≤ 6`, the model's polynomials are the theory's `cosMultiple` and `dcosMultiple`. -/
theorem cosMultipleModel_real {m : Nat} (hm : m ≤ 6) (c : ℝ) :
    cosMultipleModel ℝ m c = some (cosMultiple m c, dcosMultiple m c) := by
  match m, hm with
  | 0, _ | 1, _ => rfl
  | 2, _ | 3, _ | 4, _ | 5, _ | 6, _ =>
    simp only [cosMultipleModel, cosMultiple, dcosMultiple, Option.some.injEq, Prod.mk.injEq]
    refine ⟨?_, ?_⟩ <;> first | trivial | ring
  | n + 7, h => exact absurd h (by omega)

/-- The two `while` loops do nothing when `ζ` is already in `[lo, hi]`. -/
theorem wrap_real_of_mem (fuel : Nat) {lo hi step z : ℝ} (hlo : lo ≤ z) (hhi : z ≤ hi) :
    wrapDown ℝ fuel hi step (wrapUp ℝ fuel lo step z) = z := by
  have hu : wrapUp ℝ fuel lo step z = z := by
    cases fuel <;> simp [wrapUp, not_lt.mpr hlo]
  rw [hu]
  cases fuel <;> simp [wrapDown, not_lt.mpr hhi]

/-- A quartic bond over ℝ, `0.25` as `1/4`: the theory's energy `quarticBond` at `b = |r|` and the
force `radialForce (dV/db) r`, which is `−∇V` (`quarticBond_force_is_minus_grad`). -/
theorem quarticBondModel_real_eq (r : Vec3 ℝ) (k b0 : ℝ) (hr : r.toV3 ≠ 0) :
    out2 (quarticBondModel ℝ (1 / 4) r k b0) =
      (quarticBond k b0 (len r.toV3),
        radialForce (k * len r.toV3 * (len r.toV3 ^ 2 - b0 ^ 2)) r.toV3) := by
  rw [(quarticBond_force_is_minus_grad k b0 hr).1]
  simp only [quarticBondModel, out2, Prod.mk.injEq, quarticBond, len_sq, dot_eq]
  constructor
  · ring
  · ext <;> simp <;> ring

/-- A harmonic bond over ℝ at `|r|` not below the guard: the theory's `harmonicBond` and
`radialForce (K (b − b0)) r` (`harmonicBond_force_is_minus_grad`). -/
theorem harmonicBondModel_real_eq {guard : ℝ} (r : Vec3 ℝ) (k b0 : ℝ) (hg : guard ≤ len r.toV3) :
    out2 (harmonicBondModel ℝ guard (1 / 2) r k b0) =
      (harmonicBond k b0 (len r.toV3), radialForce (k * (len r.toV3 - b0)) r.toV3) := by
  simp only [harmonicBondModel, length_eq, not_lt.mpr hg, ite_false, out2, Prod.mk.injEq,
    harmonicBond, radialForce]
  constructor
  · ring
  · ext <;> simp <;> ring

/-- A cosine-harmonic angle over ℝ with both bonds not below the guard: the theory's
`cosineAngle` at `cos θ` and `cosineAngleForces`, which are `−∇V`
(`cosineAngle_force_is_minus_grad`) and conserve momentum and angular momentum. -/
theorem cosineAngleModel_real_eq {guard : ℝ} (rij rkj : Vec3 ℝ) (k cos0 : ℝ)
    (hij : guard ≤ len rij.toV3) (hkj : guard ≤ len rkj.toV3) :
    out3 (cosineAngleModel ℝ guard (1 / 2) rij rkj k cos0) =
      (cosineAngle k cos0 (cosAngle rij.toV3 rkj.toV3),
        cosineAngleForces k cos0 rij.toV3 rkj.toV3) := by
  have hc := clampUnit_real (cos_between_mem rij.toV3 rkj.toV3)
  simp only [cosineAngleModel, length_eq, dot_eq, not_or.mpr ⟨not_lt.mpr hij,
    not_lt.mpr hkj⟩, ite_false, cosAngle, hc, out3, cosineAngleForces, threeBody, Prod.mk.injEq,
    cosineAngle]
  refine ⟨by ring, ?_, ?_, ?_⟩ <;> ext <;> simp <;> ring

/-- A harmonic angle over ℝ with both bonds not below the guard and `|sin θ|` not below the sine
guard: the theory's `harmonicAngle` at `θ = arccos(cos θ)` and `harmonicAngleForces`. -/
theorem harmonicAngleModel_real_eq {guard sineGuard pi : ℝ} (rij rkj : Vec3 ℝ) (k theta0 : ℝ)
    (hij : guard ≤ len rij.toV3) (hkj : guard ≤ len rkj.toV3)
    (hs : sineGuard ≤ |Real.sin (Real.arccos (cosAngle rij.toV3 rkj.toV3))|) :
    out3 (harmonicAngleModel ℝ guard sineGuard (1 / 2) pi rij rkj k theta0) =
      (harmonicAngle k theta0 (Real.arccos (cosAngle rij.toV3 rkj.toV3)),
        harmonicAngleForces k theta0 rij.toV3 rkj.toV3) := by
  have hc : clampUnit ℝ (V3.dot rij.toV3 rkj.toV3 / (len rij.toV3 * len rkj.toV3)) =
      cosAngle rij.toV3 rkj.toV3 := clampUnit_real (cos_between_mem rij.toV3 rkj.toV3)
  have hs' : ¬ absLt ℝ (Real.sin (Real.arccos (cosAngle rij.toV3 rkj.toV3))) sineGuard := by
    rw [absLt_real]; exact not_lt.mpr hs
  simp only [harmonicAngleModel, length_eq, dot_eq, not_or.mpr ⟨not_lt.mpr hij,
    not_lt.mpr hkj⟩, ite_false, hc, Acos.acos, Sin.sin, hs', out3,
    harmonicAngleForces, threeBody, Prod.mk.injEq, harmonicAngle]
  refine ⟨by ring, ?_, ?_, ?_⟩ <;> ext <;> simp <;> ring

/-- The model's projections are the theory's `perpI` and `perpL`. -/
theorem toV3_perpI (rij rkj : Vec3 ℝ) :
    (rij - rkj * (dot rij rkj / dot rkj rkj)).toV3 = perpI rij.toV3 rkj.toV3 := by
  ext <;> simp [perpI, dot_eq] <;> ring

theorem toV3_perpL (rkj rkl : Vec3 ℝ) :
    (rkj * (dot rkl rkj / dot rkj rkj) - rkl).toV3 = perpL rkj.toV3 rkl.toV3 := by
  ext <;> simp [perpL, dot_eq] <;> ring

/-- The Chebyshev dihedral over ℝ, `m ≤ 6`, outside the guards (`|r_kj|²`, `|r_im|`, `|r_ln|`):
the theory's `dihedral` at `cos φ` and `chebyshevForces`, which conserve momentum and angular
momentum (`chebyshev_conservation`). With `δ ∈ {0, π}` and `1 ≤ m`, the energy is
`K (1 + cos(mφ − δ))` (`dihedral_new_eq_dihedral`). -/
theorem dihedralChebyshevModel_real_eq {guard guard2 : ℝ} (rij rkj rkl : Vec3 ℝ) (k cosDelta : ℝ)
    {m : Nat} (hm : m ≤ 6) (hkj : guard2 ≤ V3.dot rkj.toV3 rkj.toV3)
    (him : guard ≤ len (perpI rij.toV3 rkj.toV3)) (hln : guard ≤ len (perpL rkj.toV3 rkl.toV3)) :
    out4 (dihedralChebyshevModel ℝ guard guard2 rij rkj rkl k cosDelta m) =
      (dihedral k cosDelta m (cosDihedral rij.toV3 rkj.toV3 rkl.toV3),
        chebyshevForces k cosDelta m rij.toV3 rkj.toV3 rkl.toV3) := by
  have hc : clampUnit ℝ (V3.dot (perpI rij.toV3 rkj.toV3) (perpL rkj.toV3 rkl.toV3) /
      (len (perpI rij.toV3 rkj.toV3) * len (perpL rkj.toV3 rkl.toV3))) =
      cosDihedral rij.toV3 rkj.toV3 rkl.toV3 := clampUnit_real (cos_between_mem _ _)
  have hdot : dot (rij - rkj * (dot rij rkj / dot rkj rkj)) (rkj * (dot rkl rkj / dot rkj rkj) - rkl)
      = V3.dot (perpI rij.toV3 rkj.toV3) (perpL rkj.toV3 rkl.toV3) := by
    rw [dot_eq, toV3_perpI, toV3_perpL]
  have hkj' : ¬ dot rkj rkj < guard2 := not_lt.mpr hkj
  have hg : ¬ (len (perpI rij.toV3 rkj.toV3) < guard ∨ len (perpL rkj.toV3 rkl.toV3) < guard) :=
    not_or.mpr ⟨not_lt.mpr him, not_lt.mpr hln⟩
  simp only [dihedralChebyshevModel, hkj', ite_false, length_eq, hdot]
  simp only [toV3_perpI, toV3_perpL]
  simp only [hg, ite_false, hc, cosMultipleModel_real hm]
  simp only [out4, chebyshevForces, fourBody, dihedral, Prod.mk.injEq, dot_eq]
  refine ⟨trivial, ?_, ?_, ?_, ?_⟩ <;> ext <;> simp [perpI, perpL, cosDihedral] <;> ring

/-- The NTBDN = 0 dihedral over ℝ outside the guards (`|r_kj|`, `|r_im|`, `|r_ln|`, neither bond
angle linear): the theory's `dihedralNew` at the signed `φ` and `newDihedralForces`, which conserve
momentum and angular momentum (`newDihedral_conservation`). -/
theorem dihedralNewModel_real_eq {guard nearLinear : ℝ} (rij rkj rkl : Vec3 ℝ) (k delta m : ℝ)
    (hkj : guard ≤ len rkj.toV3)
    (him : guard ≤ len (perpI rij.toV3 rkj.toV3)) (hln : guard ≤ len (perpL rkj.toV3 rkl.toV3))
    (hmj : nearLinear * V3.dot rkj.toV3 rkj.toV3 <
      V3.dot (V3.cross rij.toV3 rkj.toV3) (V3.cross rij.toV3 rkj.toV3))
    (hnk : nearLinear * V3.dot rkj.toV3 rkj.toV3 <
      V3.dot (V3.cross rkj.toV3 rkl.toV3) (V3.cross rkj.toV3 rkl.toV3)) :
    out4 (dihedralNewModel ℝ guard nearLinear rij rkj rkl k delta m) =
      (dihedralNew k delta m (signedDihedral rij.toV3 rkj.toV3 rkl.toV3),
        newDihedralForces k delta m rij.toV3 rkj.toV3 rkl.toV3) := by
  have hc : clampUnit ℝ (V3.dot (perpI rij.toV3 rkj.toV3) (perpL rkj.toV3 rkl.toV3) /
      (len (perpI rij.toV3 rkj.toV3) * len (perpL rkj.toV3 rkl.toV3))) =
      cosDihedral rij.toV3 rkj.toV3 rkl.toV3 := clampUnit_real (cos_between_mem _ _)
  have hdot : dot (rij - rkj * (dot rij rkj / dot rkj rkj)) (rkj * (dot rkl rkj / dot rkj rkj) - rkl)
      = V3.dot (perpI rij.toV3 rkj.toV3) (perpL rkj.toV3 rkl.toV3) := by
    rw [dot_eq, toV3_perpI, toV3_perpL]
  have hkj' : ¬ len rkj.toV3 < guard := not_lt.mpr hkj
  have hg : ¬ (len (perpI rij.toV3 rkj.toV3) < guard ∨ len (perpL rkj.toV3 rkl.toV3) < guard) :=
    not_or.mpr ⟨not_lt.mpr him, not_lt.mpr hln⟩
  have hmj' : nearLinear * dot rkj rkj < dot (cross rij rkj) (cross rij rkj) := hmj
  have hnk' : nearLinear * dot rkj rkj < dot (cross rkj rkl) (cross rkj rkl) := hnk
  simp only [dihedralNewModel, sqrt_dot_eq, hkj', ite_false, length_eq, hdot]
  simp only [toV3_perpI, toV3_perpL]
  simp only [hg, ite_false, hc, hmj', hnk', ite_true, Acos.acos, Cos.cos, Sin.sin]
  simp only [dot_eq, toV3_cross]
  simp only [out4, newDihedralForces, fourBody, dihedralNew, signedDihedral, Prod.mk.injEq]
  -- the model's `if` and the theory's carry different `Decidable` instances: split both
  split_ifs with h <;> simp only [h, ite_true, ite_false] <;>
    refine ⟨trivial, ?_, ?_, ?_, ?_⟩ <;> ext <;> simp <;> ring

/-- The improper dihedral over ℝ outside the guards (both normals not below the guard, neither
bond angle linear), with `ζ` already in `[ζ0 − π, ζ0 + π]` (always so for `ζ0 = 0`): the theory's
`improper` at `ζ` and `improperForces`, which conserve momentum and angular momentum
(`improper_conservation`). -/
theorem improperModel_real_eq {guard2 nearLinear : ℝ} (fuel : Nat) (rij rkj rkl : Vec3 ℝ)
    (k zeta0 : ℝ)
    (hmj : guard2 ≤ len (V3.cross rij.toV3 rkj.toV3))
    (hnk : guard2 ≤ len (V3.cross rkj.toV3 rkl.toV3))
    (hlo : zeta0 - Real.pi ≤ improperAngle rij.toV3 rkj.toV3 rkl.toV3)
    (hhi : improperAngle rij.toV3 rkj.toV3 rkl.toV3 ≤ zeta0 + Real.pi)
    (hmj2 : nearLinear * V3.dot rkj.toV3 rkj.toV3 ≤
      V3.dot (V3.cross rij.toV3 rkj.toV3) (V3.cross rij.toV3 rkj.toV3))
    (hnk2 : nearLinear * V3.dot rkj.toV3 rkj.toV3 ≤
      V3.dot (V3.cross rkj.toV3 rkl.toV3) (V3.cross rkj.toV3 rkl.toV3)) :
    out4 (improperModel ℝ guard2 nearLinear (1 / 2) Real.pi fuel rij rkj rkl k zeta0) =
      (improper k zeta0 (improperAngle rij.toV3 rkj.toV3 rkl.toV3),
        improperForces k zeta0 (improperAngle rij.toV3 rkj.toV3 rkl.toV3)
          rij.toV3 rkj.toV3 rkl.toV3) := by
  have hc : clampUnit ℝ (V3.dot (V3.cross rij.toV3 rkj.toV3) (V3.cross rkj.toV3 rkl.toV3) /
      (len (V3.cross rij.toV3 rkj.toV3) * len (V3.cross rkj.toV3 rkl.toV3))) =
      cosAngle (V3.cross rij.toV3 rkj.toV3) (V3.cross rkj.toV3 rkl.toV3) :=
    clampUnit_real (cos_between_mem _ _)
  have hg : ¬ (len (V3.cross rij.toV3 rkj.toV3) < guard2 ∨
      len (V3.cross rkj.toV3 rkl.toV3) < guard2) :=
    not_or.mpr ⟨not_lt.mpr hmj, not_lt.mpr hnk⟩
  have hmj2' : ¬ dot (cross rij rkj) (cross rij rkj) < nearLinear * dot rkj rkj :=
    not_lt.mpr hmj2
  have hnk2' : ¬ dot (cross rkj rkl) (cross rkj rkl) < nearLinear * dot rkj rkj :=
    not_lt.mpr hnk2
  have hζ : (if V3.dot rij.toV3 (V3.cross rkj.toV3 rkl.toV3) < 0 then
      -Real.arccos (cosAngle (V3.cross rij.toV3 rkj.toV3) (V3.cross rkj.toV3 rkl.toV3))
      else Real.arccos (cosAngle (V3.cross rij.toV3 rkj.toV3) (V3.cross rkj.toV3 rkl.toV3))) =
      improperAngle rij.toV3 rkj.toV3 rkl.toV3 := by
    unfold improperAngle; split_ifs <;> rfl
  simp only [improperModel, sqrt_dot_eq, toV3_cross]
  simp only [hg, ite_false, hmj2', hnk2']
  simp only [dot_eq, toV3_cross, hc, Acos.acos]
  simp only [hζ, wrap_real_of_mem fuel hlo hhi]
  simp only [out4, improperForces, fourBody, improper, Prod.mk.injEq]
  refine ⟨by ring, ?_, ?_, ?_, ?_⟩ <;> ext <;> simp <;> ring

end NewenSpec.Bonded
