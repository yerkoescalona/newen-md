import NewenSpec.Float
import NewenTheory.Space.Periodicity

/-!
# Model of the periodic-boundary kernels of `crates/newen-core`

Rust (spec keys in `spec/SPEC.toml`), all in `crates/newen-core/src/math.rs`:

* `periodicity.vacuum_nearest_image`: `vacuum_nearest_image`, the plain difference;
* `periodicity.rectangular_new`, `periodicity.rectangular_nearest_image`,
  `periodicity.rectangular_put_into_box`: `Rectangular::new` (half edges and reciprocals),
  `rectangular_nearest_image` and `rectangular_put_into_box`;
* `periodicity.triclinic_new`, `periodicity.triclinic_nearest_image`,
  `periodicity.triclinic_put_into_box`: `Triclinic::new` (glam's `DMat3::inverse`),
  `triclinic_nearest_image` (a z → y → x reduction of the fractional coordinates, as
  gromosXX's comments intend; gromosXX itself reduces into a brick, the same image below half the
  smallest height: `NewenTheory…Box.image_shortest_of_lt_half_height`) and
  `triclinic_put_into_box`.
  The truncated octahedron is a triclinic box (`truncoct_triclinic_box`).

Theory: `NewenTheory/Space/Periodicity.lean`.

Each kernel is written once over any type `α` with the operations the Rust uses (`+ − × ÷`,
negation, `≤`, `rint` and `floor` from `NewenSpec/Float.lean`, the literal `1`, and glam's vector
and matrix operations in glam's order), and used twice, as in `NewenSpec/Nonbonded.lean`: at ℝ,
where the `*_real_eq` theorems prove it equal to the theory, and at `Float`, where
`lake exe spec-vectors` evaluates it and `crates/newen-core/tests/periodicity_spec.rs` compares the
Rust bit for bit. The literal `0.5` is a parameter, given at `Float` by its bits.

The triclinic reduction is six `while` loops. Lean needs every recursion to end, so each loop is
`whileFuel`, which stops after `fuel` turns; the theorems hold once `fuel` exceeds every
fractional coordinate by one, and the generator runs with a fuel far above what its inputs need
(it checks that halving the fuel changes no bit).

What is not proved: the rounding error between the two instantiations. In particular the `Float`
image can leave its cell by a rounding error, and the rectangular image computes `x · (1/L)`, not
`x / L`, so at `Float` the half-box ties are those of the rounded reciprocal.
-/

namespace NewenSpec.Periodicity

open NewenTheory.Space.Periodicity NewenTheory.Space.Periodicity.Vec3

/-- glam's `DMat3`: three columns. -/
structure Mat3 (α : Type) where
  xAxis : Vec3 α
  yAxis : Vec3 α
  zAxis : Vec3 α

/-- `while p x do x := f x`, at most `fuel` times. -/
def whileFuel {β : Type} (p : β → Prop) [DecidablePred p] (f : β → β) : Nat → β → β
  | 0, x => x
  | n + 1, x => if p x then whileFuel p f n (f x) else x

/-! ## 1. Model (mirrors the Rust functions) -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LE α]
  [DecidableRel (α := α) (· ≤ ·)] [Rint α] [Floor α] [OfNat α 1]

/-- glam `DVec3 * DVec3`, component by component. -/
def mulC (u v : Vec3 α) : Vec3 α := ⟨u.x * v.x, u.y * v.y, u.z * v.z⟩

/-- glam `DVec3 * f64`. -/
def scale (u : Vec3 α) (s : α) : Vec3 α := ⟨u.x * s, u.y * s, u.z * s⟩

/-- glam `DVec3::dot`: `(x·x' + y·y') + z·z'`. -/
def dot (u v : Vec3 α) : α := u.x * v.x + u.y * v.y + u.z * v.z

/-- glam `DVec3::cross`. -/
def cross (u v : Vec3 α) : Vec3 α :=
  ⟨u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y⟩

/-- glam `DMat3 * DVec3`: `(x_axis·v.x + y_axis·v.y) + z_axis·v.z`. -/
def mulVec (m : Mat3 α) (v : Vec3 α) : Vec3 α :=
  scale α m.xAxis v.x + scale α m.yAxis v.y + scale α m.zAxis v.z

/-- glam `DMat3::inverse`: the cross products of the columns over the determinant, transposed. -/
def inverseModel (m : Mat3 α) : Mat3 α :=
  let tmp0 := cross α m.yAxis m.zAxis
  let tmp1 := cross α m.zAxis m.xAxis
  let tmp2 := cross α m.xAxis m.yAxis
  let det := dot α m.zAxis tmp2
  let invDet := 1 / det
  let c0 := scale α tmp0 invDet
  let c1 := scale α tmp1 invDet
  let c2 := scale α tmp2 invDet
  ⟨⟨c0.x, c1.x, c2.x⟩, ⟨c0.y, c1.y, c2.y⟩, ⟨c0.z, c1.z, c2.z⟩⟩

/-- `vacuum_nearest_image(ri, rj)`. -/
def vacuumNearestImageModel (ri rj : Vec3 α) : Vec3 α := ri - rj

/-- `Rectangular::new(box_size)`: the half edges (`box_size * 0.5`) and the reciprocals
(`box_size.recip()`). -/
def rectangularNewModel (half : α) (box : Vec3 α) : Vec3 α × Vec3 α :=
  (scale α box half, ⟨1 / box.x, 1 / box.y, 1 / box.z⟩)

/-- `rectangular_nearest_image(box_size, inv_box, ri, rj)`. -/
def rectangularNearestImageModel (box inv ri rj : Vec3 α) : Vec3 α :=
  let r := ri - rj
  let k := mulC α r inv
  let n : Vec3 α := ⟨Rint.rint k.x, Rint.rint k.y, Rint.rint k.z⟩
  r - mulC α box n

/-- `rectangular_put_into_box(box_size, inv_box, pos)`. -/
def rectangularPutIntoBoxModel (box inv pos : Vec3 α) : Vec3 α :=
  let k := mulC α pos inv
  pos - mulC α box ⟨Floor.floor k.x, Floor.floor k.y, Floor.floor k.z⟩

/-- The `nim_y` of `triclinic_nearest_image`: `nim.y − (c.y / c.z)·nim.z`. -/
def nimY (c v : Vec3 α) : α := v.y - c.y / c.z * v.z

/-- The `nim_x` of `triclinic_nearest_image`:
`nim.x − b.x·(nim.y − c.y·nim.z / c.z) / b.y − c.x·nim.z / c.z`. -/
def nimX (b c v : Vec3 α) : α := v.x - b.x * (v.y - c.y * v.z / c.z) / b.y - c.x * v.z / c.z

/-- `triclinic_nearest_image(box_matrix, ri, rj)`, each `while` loop at most `fuel` times. -/
def triclinicNearestImageModel (half : α) (fuel : Nat) (m : Mat3 α) (ri rj : Vec3 α) :
    Vec3 α :=
  let a := m.xAxis
  let b := m.yAxis
  let c := m.zAxis
  let nim := ri - rj
  let nim := whileFuel (fun v => half * c.z ≤ v.z) (· - c) fuel nim
  let nim := whileFuel (fun v => v.z ≤ (-half) * c.z) (· + c) fuel nim
  let nim := whileFuel (fun v => half * b.y ≤ nimY α c v) (· - b) fuel nim
  let nim := whileFuel (fun v => nimY α c v ≤ (-half) * b.y) (· + b) fuel nim
  let nim := whileFuel (fun v => half * a.x ≤ nimX α b c v) (· - a) fuel nim
  whileFuel (fun v => nimX α b c v ≤ (-half) * a.x) (· + a) fuel nim

/-- `triclinic_put_into_box(box_matrix, inv_box_matrix, pos)`. -/
def triclinicPutIntoBoxModel (m inv : Mat3 α) (pos : Vec3 α) : Vec3 α :=
  let frac := mulVec α inv pos
  let wrapped := frac - ⟨Floor.floor frac.x, Floor.floor frac.y, Floor.floor frac.z⟩
  mulVec α m wrapped

end Model

/-! ## 2. The `Float` instantiation (what `lake exe spec-vectors` runs) -/

/-- The `0.5` of the kernels, as `f64` bits. -/
def pointFiveFloat : Float := Float.ofBits 0x3fe0000000000000

/-- The loop bound of the `Float` triclinic reduction: far above the turns the generated inputs
need (fewer than a hundred per loop). -/
def fuelFloat : Nat := 100000

def vacuumNearestImageFloat (ri rj : Vec3 Float) : Vec3 Float :=
  vacuumNearestImageModel Float ri rj

def rectangularNewFloat (box : Vec3 Float) : Vec3 Float × Vec3 Float :=
  rectangularNewModel Float pointFiveFloat box

/-- `Rectangular::new(box)` then `nearest_image(ri, rj)`. -/
def rectangularNearestImageFloat (box ri rj : Vec3 Float) : Vec3 Float :=
  rectangularNearestImageModel Float box (rectangularNewFloat box).2 ri rj

/-- `Rectangular::new(box)` then `put_into_box(pos)`. -/
def rectangularPutIntoBoxFloat (box pos : Vec3 Float) : Vec3 Float :=
  rectangularPutIntoBoxModel Float box (rectangularNewFloat box).2 pos

def triclinicNewFloat (m : Mat3 Float) : Mat3 Float := inverseModel Float m

def triclinicNearestImageFloat (fuel : Nat) (m : Mat3 Float) (ri rj : Vec3 Float) : Vec3 Float :=
  triclinicNearestImageModel Float pointFiveFloat fuel m ri rj

/-- `Triclinic::new(m)` then `put_into_box(pos)`. -/
def triclinicPutIntoBoxFloat (m : Mat3 Float) (pos : Vec3 Float) : Vec3 Float :=
  triclinicPutIntoBoxModel Float m (triclinicNewFloat m) pos

/-! ## 3. The loops -/

section Loops

variable {β : Type} (p : β → Prop) [DecidablePred p] (f : β → β)

/-- What every turn of a loop keeps holds after the loop. -/
theorem whileFuel_invariant (Q : β → Prop) (hQ : ∀ x, Q x → p x → Q (f x)) :
    ∀ (n : Nat) (x : β), Q x → Q (whileFuel p f n x) := by
  intro n
  induction n with
  | zero => intro x hx; exact hx
  | succ n ih =>
    intro x hx
    simp only [whileFuel]
    split_ifs with hp
    · exact ih _ (hQ x hx hp)
    · exact hx

/-- A loop "while `t ≤ μ x`, lower `μ` by one" ends below `t` if it has enough turns. -/
theorem whileFuel_down (μ : β → ℝ) (t : ℝ) (hp : ∀ x, p x ↔ t ≤ μ x)
    (hf : ∀ x, p x → μ (f x) = μ x - 1) :
    ∀ (n : Nat) (x : β), μ x < t + n → μ (whileFuel p f n x) < t := by
  intro n
  induction n with
  | zero => intro x hx; simpa [whileFuel] using hx
  | succ n ih =>
    intro x hx
    simp only [whileFuel]
    split_ifs with h
    · apply ih
      rw [hf x h]
      push_cast at hx
      linarith
    · exact lt_of_not_ge fun h' => h ((hp x).mpr h')

/-- A loop "while `μ x ≤ t`, raise `μ` by one" ends above `t` if it has enough turns. -/
theorem whileFuel_up (μ : β → ℝ) (t : ℝ) (hp : ∀ x, p x ↔ μ x ≤ t)
    (hf : ∀ x, p x → μ (f x) = μ x + 1) :
    ∀ (n : Nat) (x : β), t - n < μ x → t < μ (whileFuel p f n x) := by
  intro n
  induction n with
  | zero => intro x hx; simpa [whileFuel] using hx
  | succ n ih =>
    intro x hx
    simp only [whileFuel]
    split_ifs with h
    · apply ih
      rw [hf x h]
      push_cast at hx
      linarith
    · exact lt_of_not_ge fun h' => h ((hp x).mpr h')

/-- One axis of the reduction: a loop that lowers the coordinate `μ` while it is at least `1/2`,
then one that raises it while it is at most `−1/2`, leave it in `(−1/2, 1/2]`. -/
theorem whileFuel_pair {pd pu : β → Prop} [DecidablePred pd] [DecidablePred pu]
    {fd fu : β → β} (μ : β → ℝ)
    (hpd : ∀ x, pd x ↔ 1 / 2 ≤ μ x) (hfd : ∀ x, pd x → μ (fd x) = μ x - 1)
    (hpu : ∀ x, pu x ↔ μ x ≤ -(1 / 2)) (hfu : ∀ x, pu x → μ (fu x) = μ x + 1)
    (n : Nat) (x : β) (hn : |μ x| + 1 ≤ n) :
    μ (whileFuel pu fu n (whileFuel pd fd n x)) ∈ Set.Ioc (-(1 / 2)) (1 / 2) := by
  have habs := abs_le.mp (le_refl |μ x|)
  have hn0 : (0 : ℝ) ≤ n := by linarith [abs_nonneg (μ x)]
  set y := whileFuel pd fd n x
  have hy_hi : μ y < 1 / 2 := whileFuel_down pd fd μ (1 / 2) hpd hfd n x (by linarith)
  have hy_lo : -(1 / 2) - n < μ y :=
    whileFuel_invariant pd fd (fun v => -(1 / 2) - n < μ v)
      (fun v _ hp => by rw [hfd v hp]; linarith [(hpd v).mp hp]) n x (by linarith)
  constructor
  · exact whileFuel_up pu fu μ (-(1 / 2)) hpu hfu n y (by linarith)
  · exact whileFuel_invariant pu fu (fun v => μ v ≤ 1 / 2)
      (fun v _ hp => by rw [hfu v hp]; linarith [(hpu v).mp hp]) n y hy_hi.le

end Loops

/-! ## 4. The ℝ instantiation equals the theory -/

/-- The box matrix of a theory box: its edges as columns. -/
def boxMat (B : Box) : Mat3 ℝ := ⟨B.a, B.b, B.c⟩

/-- Vacuum: the plain difference. -/
theorem vacuumNearestImageModel_real_eq (ri rj : Vec3 ℝ) :
    vacuumNearestImageModel ℝ ri rj = ri - rj := rfl

/-- `Rectangular::new` over ℝ: half edges and reciprocals. -/
theorem rectangularNewModel_real_eq (L : Vec3 ℝ) :
    rectangularNewModel ℝ (1 / 2) L = (⟨L.x / 2, L.y / 2, L.z / 2⟩, ⟨1 / L.x, 1 / L.y, 1 / L.z⟩) := by
  simp only [rectangularNewModel, scale, Prod.mk.injEq, and_true]
  ext <;> simp <;> ring

/-- The rectangular nearest image over ℝ is `rectImage` in every component: `x · (1/L)` is
`x / L`, and `rint` is `roundEven`. -/
theorem rectangularNearestImageModel_real_eq (L ri rj : Vec3 ℝ) :
    rectangularNearestImageModel ℝ L (rectangularNewModel ℝ (1 / 2) L).2 ri rj =
      rectImageVec L (ri - rj) := by
  simp only [rectangularNearestImageModel, rectangularNewModel, mulC, Rint.rint, rectImageVec,
    rectImage]
  ext <;> simp [div_eq_mul_inv]

/-- The rectangular put-into-box over ℝ is `rectWrap` in every component. -/
theorem rectangularPutIntoBoxModel_real_eq (L pos : Vec3 ℝ) :
    rectangularPutIntoBoxModel ℝ L (rectangularNewModel ℝ (1 / 2) L).2 pos =
      ⟨rectWrap L.x pos.x, rectWrap L.y pos.y, rectWrap L.z pos.z⟩ := by
  simp only [rectangularPutIntoBoxModel, rectangularNewModel, mulC, Floor.floor, rectWrap]
  ext <;> simp [div_eq_mul_inv]

/-- glam's inverse over ℝ, applied to `r`, gives the fractional coordinates of a lower-triangular
box. -/
theorem mulVec_inverseModel_real {B : Box} (h : B.LowerTriangular) (r : Vec3 ℝ) :
    mulVec ℝ (inverseModel ℝ (boxMat B)) r = B.frac r := by
  have ha := h.ax_pos.ne'
  have hb := h.by_pos.ne'
  have hc := h.cz_pos.ne'
  simp only [mulVec, inverseModel, boxMat, cross, dot, scale, Box.frac, h.ay, h.az, h.bz]
  ext <;> simp <;> field_simp <;> ring

/-- glam's `DMat3 * DVec3` with the box edges as columns is `Box.point`. -/
theorem mulVec_toMat3 (B : Box) (s : Vec3 ℝ) : mulVec ℝ (boxMat B) s = B.point s := by
  simp only [mulVec, boxMat, scale, Box.point]
  ext <;> simp <;> ring

/-- The triclinic put-into-box over ℝ is `Box.wrap`: into the positive cell. -/
theorem triclinicPutIntoBoxModel_real_eq {B : Box} (h : B.LowerTriangular) (pos : Vec3 ℝ) :
    triclinicPutIntoBoxModel ℝ (boxMat B) (inverseModel ℝ (boxMat B)) pos = B.wrap pos := by
  simp only [triclinicPutIntoBoxModel, mulVec_inverseModel_real h, mulVec_toMat3, Box.wrap,
    Floor.floor]
  rw [Box.point_sub, Box.point_frac h]

/-- `nim_y` over ℝ is `b_y` times the fractional coordinate along `b`. -/
theorem nimY_real {B : Box} (h : B.LowerTriangular) (v : Vec3 ℝ) :
    nimY ℝ B.c v = (B.frac v).y * B.b.y := by
  have hb := h.by_pos.ne'
  have hc := h.cz_pos.ne'
  simp only [nimY, Box.frac]
  field_simp

/-- `nim_x` over ℝ is `a_x` times the fractional coordinate along `a`. -/
theorem nimX_real {B : Box} (h : B.LowerTriangular) (v : Vec3 ℝ) :
    nimX ℝ B.b B.c v = (B.frac v).x * B.a.x := by
  have ha := h.ax_pos.ne'
  have hb := h.by_pos.ne'
  have hc := h.cz_pos.ne'
  simp only [nimX, Box.frac]
  field_simp

/-- The triclinic nearest image over ℝ, with `half = 1/2` and enough turns for every loop, is
`Box.image`: the image of `ri − rj` in the centred cell. -/
theorem triclinicNearestImageModel_real_eq {B : Box} (h : B.LowerTriangular) {ri rj : Vec3 ℝ}
    {fuel : Nat} (hx : |(B.frac (ri - rj)).x| + 1 ≤ fuel)
    (hy : |(B.frac (ri - rj)).y| + 1 ≤ fuel) (hz : |(B.frac (ri - rj)).z| + 1 ≤ fuel) :
    triclinicNearestImageModel ℝ (1 / 2) fuel (boxMat B) ri rj = B.image (ri - rj) := by
  have hax := h.ax_pos
  have hby := h.by_pos
  have hcz := h.cz_pos
  have fa := Box.frac_a h
  have fb := Box.frac_b h
  have fc := Box.frac_c h
  -- the fractional coordinates move by one per turn, along the loop's own axis only
  have dsub : ∀ v w : Vec3 ℝ, B.frac (v - w) = B.frac v - B.frac w := Box.frac_sub B
  have dadd : ∀ v w : Vec3 ℝ, B.frac (v + w) = B.frac v + B.frac w := Box.frac_add B
  set r := ri - rj
  set fz : Vec3 ℝ → ℝ := fun v => (B.frac v).z
  set fy : Vec3 ℝ → ℝ := fun v => (B.frac v).y
  set fx : Vec3 ℝ → ℝ := fun v => (B.frac v).x
  -- the loop conditions in fractional coordinates
  have cz_le : ∀ v : Vec3 ℝ, 1 / 2 * B.c.z ≤ v.z ↔ 1 / 2 ≤ fz v := by
    intro v
    simp only [fz, Box.frac]
    rw [le_div_iff₀ hcz]
  have le_cz : ∀ v : Vec3 ℝ, v.z ≤ -(1 / 2) * B.c.z ↔ fz v ≤ -(1 / 2) := by
    intro v
    simp only [fz, Box.frac]
    rw [div_le_iff₀ hcz]
  have by_le : ∀ v : Vec3 ℝ, 1 / 2 * B.b.y ≤ nimY ℝ B.c v ↔ 1 / 2 ≤ fy v := by
    intro v
    rw [nimY_real h]
    exact mul_le_mul_iff_of_pos_right hby
  have le_by : ∀ v : Vec3 ℝ, nimY ℝ B.c v ≤ -(1 / 2) * B.b.y ↔ fy v ≤ -(1 / 2) := by
    intro v
    rw [nimY_real h]
    exact mul_le_mul_iff_of_pos_right hby
  have ax_le : ∀ v : Vec3 ℝ, 1 / 2 * B.a.x ≤ nimX ℝ B.b B.c v ↔ 1 / 2 ≤ fx v := by
    intro v
    rw [nimX_real h]
    exact mul_le_mul_iff_of_pos_right hax
  have le_ax : ∀ v : Vec3 ℝ, nimX ℝ B.b B.c v ≤ -(1 / 2) * B.a.x ↔ fx v ≤ -(1 / 2) := by
    intro v
    rw [nimX_real h]
    exact mul_le_mul_iff_of_pos_right hax
  -- the six loops
  simp only [triclinicNearestImageModel, boxMat]
  set v1 := whileFuel (fun v => 1 / 2 * B.c.z ≤ v.z) (· - B.c) fuel r
  set v2 := whileFuel (fun v => v.z ≤ (-(1 / 2)) * B.c.z) (· + B.c) fuel v1
  set v3 := whileFuel (fun v => 1 / 2 * B.b.y ≤ nimY ℝ B.c v) (· - B.b) fuel v2
  set v4 := whileFuel (fun v => nimY ℝ B.c v ≤ (-(1 / 2)) * B.b.y) (· + B.b) fuel v3
  set v5 := whileFuel (fun v => 1 / 2 * B.a.x ≤ nimX ℝ B.b B.c v) (· - B.a) fuel v4
  set v6 := whileFuel (fun v => nimX ℝ B.b B.c v ≤ (-(1 / 2)) * B.a.x) (· + B.a) fuel v5
  -- a coordinate that no loop of a stretch touches keeps its value through it
  have keep : ∀ (μ : Vec3 ℝ → ℝ) (p : Vec3 ℝ → Prop) [DecidablePred p] (f : Vec3 ℝ → Vec3 ℝ),
      (∀ v, μ (f v) = μ v) → ∀ v, μ (whileFuel p f fuel v) = μ v := by
    intro μ p _ f hf v
    exact whileFuel_invariant p f (fun w => μ w = μ v) (fun w hw _ => by rw [hf, hw]) fuel v rfl
  have lat : ∀ (p : Vec3 ℝ → Prop) [DecidablePred p] (f : Vec3 ℝ → Vec3 ℝ),
      (∀ v, B.Lattice (r - v) → B.Lattice (r - f v)) →
        ∀ v, B.Lattice (r - v) → B.Lattice (r - whileFuel p f fuel v) := by
    intro p _ f hf v hv
    exact whileFuel_invariant p f (fun w => B.Lattice (r - w)) (fun w hw _ => hf w hw) fuel v hv
  have lat_sub : ∀ (e : Vec3 ℝ), B.Lattice e → ∀ v, B.Lattice (r - v) → B.Lattice (r - (v - e)) :=
    fun e he v hv => by
      rw [show r - (v - e) = (r - v) + e by ext <;> simp <;> ring]
      exact hv.add he
  have lat_add : ∀ (e : Vec3 ℝ), B.Lattice e → ∀ v, B.Lattice (r - v) → B.Lattice (r - (v + e)) :=
    fun e he v hv => by
      rw [show r - (v + e) = (r - v) + -e by ext <;> simp <;> ring]
      exact hv.add he.neg
  -- each turn's effect on the three coordinates
  have zc : ∀ v, fz (v - B.c) = fz v - 1 := fun v => by simp [fz, dsub, fc]
  have zc' : ∀ v, fz (v + B.c) = fz v + 1 := fun v => by simp [fz, dadd, fc]
  have yc : ∀ v, fy (v - B.c) = fy v := fun v => by simp [fy, dsub, fc]
  have yc' : ∀ v, fy (v + B.c) = fy v := fun v => by simp [fy, dadd, fc]
  have xc : ∀ v, fx (v - B.c) = fx v := fun v => by simp [fx, dsub, fc]
  have xc' : ∀ v, fx (v + B.c) = fx v := fun v => by simp [fx, dadd, fc]
  have zb : ∀ v, fz (v - B.b) = fz v := fun v => by simp [fz, dsub, fb]
  have zb' : ∀ v, fz (v + B.b) = fz v := fun v => by simp [fz, dadd, fb]
  have yb : ∀ v, fy (v - B.b) = fy v - 1 := fun v => by simp [fy, dsub, fb]
  have yb' : ∀ v, fy (v + B.b) = fy v + 1 := fun v => by simp [fy, dadd, fb]
  have xb : ∀ v, fx (v - B.b) = fx v := fun v => by simp [fx, dsub, fb]
  have xb' : ∀ v, fx (v + B.b) = fx v := fun v => by simp [fx, dadd, fb]
  have za : ∀ v, fz (v - B.a) = fz v := fun v => by simp [fz, dsub, fa]
  have za' : ∀ v, fz (v + B.a) = fz v := fun v => by simp [fz, dadd, fa]
  have ya : ∀ v, fy (v - B.a) = fy v := fun v => by simp [fy, dsub, fa]
  have ya' : ∀ v, fy (v + B.a) = fy v := fun v => by simp [fy, dadd, fa]
  have xa : ∀ v, fx (v - B.a) = fx v - 1 := fun v => by simp [fx, dsub, fa]
  have xa' : ∀ v, fx (v + B.a) = fx v + 1 := fun v => by simp [fx, dadd, fa]
  -- the coordinates at the start of each axis' pair of loops are those of `r`
  have y2 : fy v2 = fy r := by rw [keep fy _ _ yc', keep fy _ _ yc]
  have x4 : fx v4 = fx r := by
    rw [keep fx _ _ xb', keep fx _ _ xb, keep fx _ _ xc', keep fx _ _ xc]
  -- each pair of loops puts its coordinate into (−1/2, 1/2]
  have hz2 : fz v2 ∈ Set.Ioc (-(1 / 2)) (1 / 2) :=
    whileFuel_pair fz (fun v => cz_le v) (fun v _ => zc v)
      (fun v => le_cz v) (fun v _ => zc' v) fuel r hz
  have hy4 : fy v4 ∈ Set.Ioc (-(1 / 2)) (1 / 2) :=
    whileFuel_pair fy (fun v => by_le v) (fun v _ => yb v)
      (fun v => le_by v) (fun v _ => yb' v) fuel v2 (by rw [y2]; exact hy)
  have hx6 : fx v6 ∈ Set.Ioc (-(1 / 2)) (1 / 2) :=
    whileFuel_pair fx (fun v => ax_le v) (fun v _ => xa v)
      (fun v => le_ax v) (fun v _ => xa' v) fuel v4 (by rw [x4]; exact hx)
  have hz6 : fz v6 = fz v2 := by
    rw [keep fz _ _ za', keep fz _ _ za, keep fz _ _ zb', keep fz _ _ zb]
  have hy6 : fy v6 = fy v4 := by rw [keep fy _ _ ya', keep fy _ _ ya]
  have hcell : B.InCell v6 :=
    ⟨hx6, by rw [show (B.frac v6).y = fy v6 from rfl, hy6]; exact hy4,
      by rw [show (B.frac v6).z = fz v6 from rfl, hz6]; exact hz2⟩
  have hlat : B.Lattice (r - v6) := by
    have h0 : B.Lattice (r - r) := by
      rw [show r - r = 0 by ext <;> simp]; exact Box.lattice_zero
    exact lat _ _ (lat_add _ Box.lattice_a) _ (lat _ _ (lat_sub _ Box.lattice_a) _
      (lat _ _ (lat_add _ Box.lattice_b) _ (lat _ _ (lat_sub _ Box.lattice_b) _
      (lat _ _ (lat_add _ Box.lattice_c) _ (lat _ _ (lat_sub _ Box.lattice_c) _ h0)))))
  exact (Box.image_eq_of h hcell hlat).symm

/-- The triclinic nearest image over ℝ is idempotent once it has two turns per loop: the image
lies in the centred cell, whose fractional coordinates are at most `1/2` in size. -/
theorem triclinicNearestImageModel_real_idem {B : Box} (h : B.LowerTriangular) {r : Vec3 ℝ}
    {fuel : Nat} (hfuel : 2 ≤ fuel) (hc : B.InCell r) :
    triclinicNearestImageModel ℝ (1 / 2) fuel (boxMat B) r 0 = r := by
  have hf : (2 : ℝ) ≤ fuel := by exact_mod_cast hfuel
  have bound : ∀ s ∈ Set.Ioc (-(1 / 2) : ℝ) (1 / 2), |s| + 1 ≤ fuel := fun s hs => by
    have := abs_le.mpr (⟨by linarith [hs.1], hs.2⟩ : -(1 / 2 : ℝ) ≤ s ∧ s ≤ 1 / 2)
    linarith
  have e : r - 0 = r := by ext <;> simp
  rw [triclinicNearestImageModel_real_eq h (by rw [e]; exact bound _ hc.1)
    (by rw [e]; exact bound _ hc.2.1) (by rw [e]; exact bound _ hc.2.2), e,
    Box.image_of_inCell h hc]

end NewenSpec.Periodicity
