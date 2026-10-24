module

public import Mathlib.Algebra.Order.Round
public import Mathlib.Algebra.Order.Archimedean.Real.Basic
public import Mathlib.Analysis.Real.Sqrt
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# Periodic boundaries: the lattice, its cells and the nearest image

A periodic system is one box repeated through space. The box is spanned by three edge vectors
`a`, `b`, `c`; the copies of a point `r` are the points `r + i·a + j·b + k·c` for integers
`i, j, k`, so two points are images of each other exactly when they differ by a **lattice vector**
(`Box.Lattice`). Every image of a point carries the same physics; a program picks one of them
with a rule, and the rule is a choice of **cell**: a region that holds exactly one image of every
point.

GROMOS writes the box with `a` along `x` and `b` in the `xy` plane (`a = (aₓ, 0, 0)`,
`b = (bₓ, b_y, 0)`, `c = (cₓ, c_y, c_z)`, `Box.LowerTriangular`). Then the fractional coordinates
`s` of a point (`r = sₐ·a + s_b·b + s_c·c`) follow by back substitution (`Box.frac`), and two cells
are used:

* the **centred cell**, fractional coordinates in `(−1/2, 1/2]` (`Box.InCell`), where the
  nearest image of a pair lives (`Box.image`, newen's triclinic `nearest_image`);
* the **positive cell**, fractional coordinates in `[0, 1)` (`Box.InPositiveCell`), where
  newen's `put_into_box` puts a point (`Box.wrap`).

For both: the result differs from the input by a lattice vector, lies in its cell, is the only
point that does (`image_eq_of`, `wrap_eq_of`), and so is unchanged when applied again.

The centred cell is a parallelepiped. For a rectangular box it is the set of shortest images,
but not for a skewed box: there a lattice vector can make the image shorter
(`image_not_shortest`). gromosXX's triclinic `nearest_image` is meant to be this reduction, but
it reads `Box(i, j)` as component `i` of edge `j`, so its skew terms vanish and it reduces into the
brick `|x| < aₓ/2, |y| < b_y/2, |z| < c_z/2` instead. Below half the smallest height of the box
both return the vector itself, the unique shortest image (`image_shortest_of_lt_half_height`);
beyond it the two can choose different images, neither of them always the shortest.

The rectangular box is treated component by component: `rectImage L x = x − L·roundEven(x/L)`,
with `roundEven` the round-to-nearest-even of IEEE 754 (C `rint`; Mathlib's `round` sends halves
up). It lies in the closed interval `[−L/2, L/2]` (both ends are reached: `L/2` stays, `3L/2`
becomes `−L/2`), is unchanged when applied again, and is the shortest image
(`rectImage_shortest`, `rectImageVec_shortest`). gromosXX's rectangular loop
(`while x > L/2: x −= L`, then `while x < −L/2: x += L`) also ends in `[−L/2, L/2]` but keeps
`3L/2` at `L/2`: at `±3L/2, ±7L/2, …` the two choose images a box length apart, of the same
length.
-/

@[expose] public section

namespace NewenTheory.Space.Periodicity

/-! ## Vectors -/

/-- A point or displacement in space, as its three Cartesian components. -/
@[ext] structure Vec3 (α : Type) where
  x : α
  y : α
  z : α

namespace Vec3

variable {α : Type}

instance [Add α] : Add (Vec3 α) := ⟨fun u v => ⟨u.x + v.x, u.y + v.y, u.z + v.z⟩⟩
instance [Sub α] : Sub (Vec3 α) := ⟨fun u v => ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩⟩
instance [Neg α] : Neg (Vec3 α) := ⟨fun u => ⟨-u.x, -u.y, -u.z⟩⟩
instance [Zero α] : Zero (Vec3 α) := ⟨⟨0, 0, 0⟩⟩
/-- A number times a vector. -/
instance [Mul α] : SMul α (Vec3 α) := ⟨fun s u => ⟨s * u.x, s * u.y, s * u.z⟩⟩

@[simp] theorem add_x [Add α] (u v : Vec3 α) : (u + v).x = u.x + v.x := rfl
@[simp] theorem add_y [Add α] (u v : Vec3 α) : (u + v).y = u.y + v.y := rfl
@[simp] theorem add_z [Add α] (u v : Vec3 α) : (u + v).z = u.z + v.z := rfl
@[simp] theorem sub_x [Sub α] (u v : Vec3 α) : (u - v).x = u.x - v.x := rfl
@[simp] theorem sub_y [Sub α] (u v : Vec3 α) : (u - v).y = u.y - v.y := rfl
@[simp] theorem sub_z [Sub α] (u v : Vec3 α) : (u - v).z = u.z - v.z := rfl
@[simp] theorem neg_x [Neg α] (u : Vec3 α) : (-u).x = -u.x := rfl
@[simp] theorem neg_y [Neg α] (u : Vec3 α) : (-u).y = -u.y := rfl
@[simp] theorem neg_z [Neg α] (u : Vec3 α) : (-u).z = -u.z := rfl
@[simp] theorem zero_x [Zero α] : (0 : Vec3 α).x = 0 := rfl
@[simp] theorem zero_y [Zero α] : (0 : Vec3 α).y = 0 := rfl
@[simp] theorem zero_z [Zero α] : (0 : Vec3 α).z = 0 := rfl
@[simp] theorem smul_x [Mul α] (s : α) (u : Vec3 α) : (s • u).x = s * u.x := rfl
@[simp] theorem smul_y [Mul α] (s : α) (u : Vec3 α) : (s • u).y = s * u.y := rfl
@[simp] theorem smul_z [Mul α] (s : α) (u : Vec3 α) : (s • u).z = s * u.z := rfl

/-- The squared length `|u|²`. -/
def normSq (u : Vec3 ℝ) : ℝ := u.x ^ 2 + u.y ^ 2 + u.z ^ 2

end Vec3

open Vec3

/-! ## The rectangular box, component by component -/

/-- Round to the nearest integer, halves to the even one: IEEE 754 roundToIntegralTiesToEven,
C `rint`. -/
noncomputable def roundEven (x : ℝ) : ℤ :=
  if x - ⌊x⌋ = 1 / 2 then (if Even ⌊x⌋ then ⌊x⌋ else ⌊x⌋ + 1) else round x

/-- `roundEven` is a nearest integer. -/
theorem abs_sub_roundEven_le (x : ℝ) : |x - roundEven x| ≤ 1 / 2 := by
  unfold roundEven
  split_ifs with h he
  · rw [h, abs_of_pos (by norm_num)]
  · push_cast
    rw [show x - (⌊x⌋ + 1) = (x - ⌊x⌋) - 1 by ring, h]
    norm_num
  · exact abs_sub_round x

/-- `roundEven` keeps `[−1/2, 1/2]` at `0`, both halves included (`0` is even). -/
theorem roundEven_eq_zero {x : ℝ} (h : |x| ≤ 1 / 2) : roundEven x = 0 := by
  have hl := (abs_le.mp h).1
  have hu := (abs_le.mp h).2
  unfold roundEven
  split_ifs with ht he
  · have h0 : ⌊x⌋ = 0 ∨ ⌊x⌋ = -1 := by
      have h1 : (⌊x⌋ : ℝ) ≤ x := Int.floor_le x
      have h2 : x < ⌊x⌋ + 1 := Int.lt_floor_add_one x
      have h3 : -1 < ⌊x⌋ + 1 := by exact_mod_cast (show (-1 : ℝ) < ⌊x⌋ + 1 by linarith)
      have h4 : ⌊x⌋ < 1 := by exact_mod_cast (show (⌊x⌋ : ℝ) < 1 by linarith)
      omega
    rcases h0 with h0 | h0
    · exact h0
    · exact absurd (h0 ▸ he) (by decide)
  · have h0 : ⌊x⌋ = 0 ∨ ⌊x⌋ = -1 := by
      have h1 : (⌊x⌋ : ℝ) ≤ x := Int.floor_le x
      have h2 : x < ⌊x⌋ + 1 := Int.lt_floor_add_one x
      have h3 : -1 < ⌊x⌋ + 1 := by exact_mod_cast (show (-1 : ℝ) < ⌊x⌋ + 1 by linarith)
      have h4 : ⌊x⌋ < 1 := by exact_mod_cast (show (⌊x⌋ : ℝ) < 1 by linarith)
      omega
    rcases h0 with h0 | h0
    · exact absurd (h0 ▸ he) (by decide)
    · omega
  · rw [round_eq, Int.floor_eq_iff]
    have hne : x ≠ 1 / 2 := by
      rintro rfl
      exact ht (by norm_num [Int.floor_eq_iff])
    have hne' : x ≠ -(1 / 2) := by
      rintro rfl
      exact ht (by norm_num [Int.floor_eq_iff])
    constructor
    · push_cast
      linarith [lt_of_le_of_ne hl (Ne.symm hne')]
    · push_cast
      linarith [lt_of_le_of_ne hu hne]

/-- The nearest image of a separation `x` along a rectangular edge of length `L`, as newen
computes it: `x` minus the nearest multiple of `L`, halves to the even multiple. -/
noncomputable def rectImage (L x : ℝ) : ℝ := x - L * roundEven (x / L)

/-- A rectangular coordinate moved into `[0, L)`. -/
noncomputable def rectWrap (L x : ℝ) : ℝ := x - L * ⌊x / L⌋

/-- **Lattice**: the image differs from `x` by a whole number of edges. -/
theorem rectImage_lattice (L x : ℝ) : ∃ n : ℤ, x - rectImage L x = n * L :=
  ⟨roundEven (x / L), by unfold rectImage; ring⟩

/-- **In the cell**: the image lies in `[−L/2, L/2]`. -/
theorem rectImage_mem {L : ℝ} (hL : 0 < L) (x : ℝ) : |rectImage L x| ≤ L / 2 := by
  have h := abs_sub_roundEven_le (x / L)
  have e : rectImage L x = L * (x / L - roundEven (x / L)) := by
    unfold rectImage; field_simp
  rw [e, abs_mul, abs_of_pos hL]
  nlinarith

/-- **Idempotent**: the image of an image is itself. -/
theorem rectImage_idem {L : ℝ} (hL : 0 < L) (x : ℝ) :
    rectImage L (rectImage L x) = rectImage L x := by
  have h : |rectImage L x / L| ≤ 1 / 2 := by
    rw [abs_div, abs_of_pos hL, div_le_iff₀ hL]
    linarith [rectImage_mem hL x]
  conv_lhs => rw [rectImage, roundEven_eq_zero h]
  simp

/-- **Shortest**: no other image of `x` is closer to the origin. -/
theorem rectImage_shortest {L : ℝ} (hL : 0 < L) (x : ℝ) (n : ℤ) :
    |rectImage L x| ≤ |x - n * L| := by
  set t := x / L
  have ht : x = t * L := by simp [t, hL.ne']
  have e1 : rectImage L x = (t - roundEven t) * L := by
    unfold rectImage; rw [ht]; simp only [t]; field_simp
  have e2 : x - n * L = (t - n) * L := by rw [ht]; ring
  rw [e1, e2, abs_mul, abs_mul]
  refine mul_le_mul_of_nonneg_right ?_ (abs_nonneg L)
  have h1 := abs_sub_roundEven_le t
  by_cases hn : n = roundEven t
  · rw [hn]
  · have h2 : (1 : ℝ) ≤ |(n : ℝ) - roundEven t| := by
      have : (1 : ℤ) ≤ |n - roundEven t| := Int.one_le_abs (sub_ne_zero.mpr hn)
      exact_mod_cast this
    have h3 : |(n : ℝ) - roundEven t| ≤ |t - n| + |t - roundEven t| := by
      calc |(n : ℝ) - roundEven t| = |(t - roundEven t) - (t - n)| := by ring_nf
        _ ≤ |t - roundEven t| + |t - n| := abs_sub _ _
        _ = |t - n| + |t - roundEven t| := add_comm _ _
    linarith

/-- The wrapped coordinate differs from `x` by a whole number of edges. -/
theorem rectWrap_lattice (L x : ℝ) : ∃ n : ℤ, x - rectWrap L x = n * L :=
  ⟨⌊x / L⌋, by unfold rectWrap; ring⟩

/-- The wrapped coordinate lies in `[0, L)`. -/
theorem rectWrap_mem {L : ℝ} (hL : 0 < L) (x : ℝ) : 0 ≤ rectWrap L x ∧ rectWrap L x < L := by
  have e : rectWrap L x = L * (x / L - ⌊x / L⌋) := by unfold rectWrap; field_simp
  have h1 := Int.floor_le (x / L)
  have h2 := Int.lt_floor_add_one (x / L)
  rw [e]
  constructor
  · exact mul_nonneg hL.le (by linarith)
  · nlinarith

/-- Wrapping is idempotent. -/
theorem rectWrap_idem {L : ℝ} (hL : 0 < L) (x : ℝ) :
    rectWrap L (rectWrap L x) = rectWrap L x := by
  obtain ⟨h0, h1⟩ := rectWrap_mem hL x
  have hf : ⌊rectWrap L x / L⌋ = 0 := by
    rw [Int.floor_eq_iff]
    constructor
    · push_cast; exact div_nonneg h0 hL.le
    · push_cast; rw [div_lt_iff₀ hL]; linarith
  conv_lhs => rw [rectWrap, hf]
  simp

/-- The rectangular image of a separation `r`, component by component, in the box of edges
`L.x, L.y, L.z`. -/
noncomputable def rectImageVec (L r : Vec3 ℝ) : Vec3 ℝ :=
  ⟨rectImage L.x r.x, rectImage L.y r.y, rectImage L.z r.z⟩

/-- **Shortest** in space: the image is no longer than any other image `r − (i Lₓ, j L_y, k L_z)`
of `r`. -/
theorem rectImageVec_shortest {L : Vec3 ℝ} (hx : 0 < L.x) (hy : 0 < L.y) (hz : 0 < L.z)
    (r : Vec3 ℝ) (i j k : ℤ) :
    normSq (rectImageVec L r) ≤ normSq (r - ⟨i * L.x, j * L.y, k * L.z⟩) := by
  have sq_le : ∀ {u v : ℝ}, |u| ≤ |v| → u ^ 2 ≤ v ^ 2 := fun h => sq_le_sq.mpr h
  unfold normSq rectImageVec
  simp only [sub_x, sub_y, sub_z]
  have := sq_le (rectImage_shortest hx r.x i)
  have := sq_le (rectImage_shortest hy r.y j)
  have := sq_le (rectImage_shortest hz r.z k)
  linarith

/-! ## The triclinic box -/

/-- A periodic box: its edge vectors, the columns of GROMOS's box matrix. -/
structure Box where
  a : Vec3 ℝ
  b : Vec3 ℝ
  c : Vec3 ℝ

namespace Box

/-- GROMOS's convention: `a` along `x`, `b` in the `xy` plane, the diagonal positive. -/
structure LowerTriangular (B : Box) : Prop where
  ay : B.a.y = 0
  az : B.a.z = 0
  bz : B.b.z = 0
  ax_pos : 0 < B.a.x
  by_pos : 0 < B.b.y
  cz_pos : 0 < B.c.z

/-- The point with fractional coordinates `s`: `sₐ·a + s_b·b + s_c·c`. -/
def point (B : Box) (s : Vec3 ℝ) : Vec3 ℝ := s.x • B.a + s.y • B.b + s.z • B.c

/-- The fractional coordinates of `r` in a lower-triangular box, by back substitution: `c` alone
reaches `z`, then `b` the rest of `y`, then `a` the rest of `x`. -/
noncomputable def frac (B : Box) (r : Vec3 ℝ) : Vec3 ℝ :=
  let sc := r.z / B.c.z
  let sb := (r.y - sc * B.c.y) / B.b.y
  ⟨(r.x - sb * B.b.x - sc * B.c.x) / B.a.x, sb, sc⟩

/-- The integer combinations of the edges: the vectors between images of one point. -/
def Lattice (B : Box) (w : Vec3 ℝ) : Prop := ∃ i j k : ℤ, w = B.point ⟨i, j, k⟩

/-- In the centred cell: every fractional coordinate in `(−1/2, 1/2]`. -/
def InCell (B : Box) (r : Vec3 ℝ) : Prop :=
  (B.frac r).x ∈ Set.Ioc (-(1 / 2)) (1 / 2) ∧ (B.frac r).y ∈ Set.Ioc (-(1 / 2)) (1 / 2) ∧
    (B.frac r).z ∈ Set.Ioc (-(1 / 2)) (1 / 2)

/-- In the positive cell: every fractional coordinate in `[0, 1)`. -/
def InPositiveCell (B : Box) (r : Vec3 ℝ) : Prop :=
  (B.frac r).x ∈ Set.Ico 0 1 ∧ (B.frac r).y ∈ Set.Ico 0 1 ∧ (B.frac r).z ∈ Set.Ico 0 1

/-- The nearest image in the centred cell: from each fractional coordinate `s` subtract the
integer `⌈s − 1/2⌉`, the one that brings it into `(−1/2, 1/2]`. -/
noncomputable def image (B : Box) (r : Vec3 ℝ) : Vec3 ℝ :=
  let s := B.frac r
  r - B.point ⟨(⌈s.x - 1 / 2⌉ : ℝ), (⌈s.y - 1 / 2⌉ : ℝ), (⌈s.z - 1 / 2⌉ : ℝ)⟩

/-- The image in the positive cell: subtract `⌊s⌋` from each fractional coordinate. -/
noncomputable def wrap (B : Box) (r : Vec3 ℝ) : Vec3 ℝ :=
  let s := B.frac r
  r - B.point ⟨(⌊s.x⌋ : ℝ), (⌊s.y⌋ : ℝ), (⌊s.z⌋ : ℝ)⟩

variable {B : Box}

theorem point_sub (B : Box) (s t : Vec3 ℝ) : B.point (s - t) = B.point s - B.point t := by
  ext <;> simp [point] <;> ring

theorem point_add (B : Box) (s t : Vec3 ℝ) : B.point (s + t) = B.point s + B.point t := by
  ext <;> simp [point] <;> ring

/-- `frac` undoes `point`. -/
theorem frac_point (h : B.LowerTriangular) (s : Vec3 ℝ) : B.frac (B.point s) = s := by
  have ha := h.ax_pos.ne'
  have hb := h.by_pos.ne'
  have hc := h.cz_pos.ne'
  ext <;> simp [frac, point, h.ay, h.az, h.bz] <;> field_simp <;> ring

/-- `point` undoes `frac`. -/
theorem point_frac (h : B.LowerTriangular) (r : Vec3 ℝ) : B.point (B.frac r) = r := by
  have ha := h.ax_pos.ne'
  have hb := h.by_pos.ne'
  have hc := h.cz_pos.ne'
  ext <;> simp [frac, point, h.ay, h.az, h.bz] <;> field_simp <;> ring

/-- `frac` is linear. -/
theorem frac_sub (B : Box) (u v : Vec3 ℝ) : B.frac (u - v) = B.frac u - B.frac v := by
  ext <;> simp [frac] <;> ring

theorem frac_add (B : Box) (u v : Vec3 ℝ) : B.frac (u + v) = B.frac u + B.frac v := by
  ext <;> simp [frac] <;> ring

theorem frac_a (h : B.LowerTriangular) : B.frac B.a = ⟨1, 0, 0⟩ := by
  have e : B.point ⟨1, 0, 0⟩ = B.a := by ext <;> simp [point]
  rw [← e, frac_point h]

theorem frac_b (h : B.LowerTriangular) : B.frac B.b = ⟨0, 1, 0⟩ := by
  have e : B.point ⟨0, 1, 0⟩ = B.b := by ext <;> simp [point]
  rw [← e, frac_point h]

theorem frac_c (h : B.LowerTriangular) : B.frac B.c = ⟨0, 0, 1⟩ := by
  have e : B.point ⟨0, 0, 1⟩ = B.c := by ext <;> simp [point]
  rw [← e, frac_point h]

theorem lattice_zero : B.Lattice 0 := ⟨0, 0, 0, by ext <;> simp [point]⟩

theorem Lattice.add {u v : Vec3 ℝ} (hu : B.Lattice u) (hv : B.Lattice v) : B.Lattice (u + v) := by
  obtain ⟨i, j, k, rfl⟩ := hu
  obtain ⟨i', j', k', rfl⟩ := hv
  exact ⟨i + i', j + j', k + k', by rw [← point_add]; congr 1; ext <;> simp⟩

theorem Lattice.neg {u : Vec3 ℝ} (hu : B.Lattice u) : B.Lattice (-u) := by
  obtain ⟨i, j, k, rfl⟩ := hu
  exact ⟨-i, -j, -k, by ext <;> simp [point] <;> ring⟩

theorem lattice_a : B.Lattice B.a := ⟨1, 0, 0, by ext <;> simp [point]⟩
theorem lattice_b : B.Lattice B.b := ⟨0, 1, 0, by ext <;> simp [point]⟩
theorem lattice_c : B.Lattice B.c := ⟨0, 0, 1, by ext <;> simp [point]⟩

/-- A lattice vector has integer fractional coordinates. -/
theorem frac_of_lattice (h : B.LowerTriangular) {w : Vec3 ℝ} (hw : B.Lattice w) :
    ∃ i j k : ℤ, B.frac w = ⟨i, j, k⟩ := by
  obtain ⟨i, j, k, rfl⟩ := hw
  exact ⟨i, j, k, frac_point h _⟩

/-- Two numbers of `(−1/2, 1/2]` that differ by an integer are equal. -/
theorem eq_of_Ioc {u v : ℝ} {n : ℤ} (hu : u ∈ Set.Ioc (-(1 / 2)) (1 / 2))
    (hv : v ∈ Set.Ioc (-(1 / 2)) (1 / 2)) (h : u - v = n) : u = v := by
  obtain ⟨hu1, hu2⟩ := hu
  obtain ⟨hv1, hv2⟩ := hv
  have h1 : (n : ℝ) < 1 := by linarith
  have h2 : (-1 : ℝ) < n := by linarith
  have : n = 0 := by
    have h1' : n < 1 := by exact_mod_cast h1
    have h2' : -1 < n := by exact_mod_cast h2
    omega
  subst this
  simp only [Int.cast_zero] at h
  linarith

/-- Two numbers of `[0, 1)` that differ by an integer are equal. -/
theorem eq_of_Ico {u v : ℝ} {n : ℤ} (hu : u ∈ Set.Ico 0 1) (hv : v ∈ Set.Ico 0 1)
    (h : u - v = n) : u = v := by
  obtain ⟨hu1, hu2⟩ := hu
  obtain ⟨hv1, hv2⟩ := hv
  have : n = 0 := by
    have h1' : n < 1 := by exact_mod_cast (show (n : ℝ) < 1 by linarith)
    have h2' : -1 < n := by exact_mod_cast (show (-1 : ℝ) < n by linarith)
    omega
  subst this
  simp only [Int.cast_zero] at h
  linarith

/-- **One image per cell**: two points of the centred cell that differ by a lattice vector are
the same point. -/
theorem inCell_unique (h : B.LowerTriangular) {u v : Vec3 ℝ} (hu : B.InCell u)
    (hv : B.InCell v) (huv : B.Lattice (u - v)) : u = v := by
  obtain ⟨i, j, k, hijk⟩ := frac_of_lattice h huv
  rw [frac_sub] at hijk
  have hx : (B.frac u).x = (B.frac v).x :=
    eq_of_Ioc hu.1 hv.1 (by simpa using congrArg Vec3.x hijk)
  have hy : (B.frac u).y = (B.frac v).y :=
    eq_of_Ioc hu.2.1 hv.2.1 (by simpa using congrArg Vec3.y hijk)
  have hz : (B.frac u).z = (B.frac v).z :=
    eq_of_Ioc hu.2.2 hv.2.2 (by simpa using congrArg Vec3.z hijk)
  rw [← point_frac h u, ← point_frac h v]
  congr 1
  ext <;> assumption

/-- One image per positive cell. -/
theorem inPositiveCell_unique (h : B.LowerTriangular) {u v : Vec3 ℝ} (hu : B.InPositiveCell u)
    (hv : B.InPositiveCell v) (huv : B.Lattice (u - v)) : u = v := by
  obtain ⟨i, j, k, hijk⟩ := frac_of_lattice h huv
  rw [frac_sub] at hijk
  have hx : (B.frac u).x = (B.frac v).x :=
    eq_of_Ico hu.1 hv.1 (by simpa using congrArg Vec3.x hijk)
  have hy : (B.frac u).y = (B.frac v).y :=
    eq_of_Ico hu.2.1 hv.2.1 (by simpa using congrArg Vec3.y hijk)
  have hz : (B.frac u).z = (B.frac v).z :=
    eq_of_Ico hu.2.2 hv.2.2 (by simpa using congrArg Vec3.z hijk)
  rw [← point_frac h u, ← point_frac h v]
  congr 1
  ext <;> assumption

/-- **Lattice**: the image differs from `r` by a lattice vector. -/
theorem image_lattice (B : Box) (r : Vec3 ℝ) : B.Lattice (r - B.image r) := by
  refine ⟨⌈(B.frac r).x - 1 / 2⌉, ⌈(B.frac r).y - 1 / 2⌉, ⌈(B.frac r).z - 1 / 2⌉, ?_⟩
  ext <;> simp [image]

theorem centred_mem (s : ℝ) : s - ⌈s - 1 / 2⌉ ∈ Set.Ioc (-(1 / 2)) (1 / 2) := by
  have h1 := Int.le_ceil (s - 1 / 2)
  have h2 := Int.ceil_lt_add_one (s - 1 / 2)
  constructor <;> linarith

/-- **In the cell**: the image lies in the centred cell. -/
theorem image_inCell (h : B.LowerTriangular) (r : Vec3 ℝ) : B.InCell (B.image r) := by
  have e : B.frac (B.image r) = ⟨(B.frac r).x - ⌈(B.frac r).x - 1 / 2⌉,
      (B.frac r).y - ⌈(B.frac r).y - 1 / 2⌉, (B.frac r).z - ⌈(B.frac r).z - 1 / 2⌉⟩ := by
    simp only [image, frac_sub, frac_point h]
    rfl
  unfold InCell
  rw [e]
  exact ⟨centred_mem _, centred_mem _, centred_mem _⟩

/-- The image is the only point of the centred cell that is an image of `r`. -/
theorem image_eq_of (h : B.LowerTriangular) {r v : Vec3 ℝ} (hv : B.InCell v)
    (hrv : B.Lattice (r - v)) : B.image r = v := by
  refine inCell_unique h (image_inCell h r) hv ?_
  have : B.image r - v = (r - v) + -(r - B.image r) := by ext <;> simp
  rw [this]
  exact hrv.add (image_lattice B r).neg

/-- **Idempotent**: a point of the centred cell is its own image, so `image (image r) = image r`. -/
theorem image_of_inCell (h : B.LowerTriangular) {r : Vec3 ℝ} (hr : B.InCell r) :
    B.image r = r :=
  image_eq_of h hr (by rw [show r - r = 0 by ext <;> simp]; exact lattice_zero)

theorem image_idem (h : B.LowerTriangular) (r : Vec3 ℝ) : B.image (B.image r) = B.image r :=
  image_of_inCell h (image_inCell h r)

theorem wrap_lattice (B : Box) (r : Vec3 ℝ) : B.Lattice (r - B.wrap r) := by
  refine ⟨⌊(B.frac r).x⌋, ⌊(B.frac r).y⌋, ⌊(B.frac r).z⌋, ?_⟩
  ext <;> simp [wrap]

theorem positive_mem (s : ℝ) : s - ⌊s⌋ ∈ Set.Ico 0 1 := by
  have h1 := Int.floor_le s
  have h2 := Int.lt_floor_add_one s
  constructor <;> linarith

/-- The wrapped point lies in the positive cell. -/
theorem wrap_inPositiveCell (h : B.LowerTriangular) (r : Vec3 ℝ) :
    B.InPositiveCell (B.wrap r) := by
  have e : B.frac (B.wrap r) = ⟨(B.frac r).x - ⌊(B.frac r).x⌋, (B.frac r).y - ⌊(B.frac r).y⌋,
      (B.frac r).z - ⌊(B.frac r).z⌋⟩ := by
    simp only [wrap, frac_sub, frac_point h]
    rfl
  unfold InPositiveCell
  rw [e]
  exact ⟨positive_mem _, positive_mem _, positive_mem _⟩

theorem wrap_eq_of (h : B.LowerTriangular) {r v : Vec3 ℝ} (hv : B.InPositiveCell v)
    (hrv : B.Lattice (r - v)) : B.wrap r = v := by
  refine inPositiveCell_unique h (wrap_inPositiveCell h r) hv ?_
  have : B.wrap r - v = (r - v) + -(r - B.wrap r) := by ext <;> simp
  rw [this]
  exact hrv.add (wrap_lattice B r).neg

theorem wrap_idem (h : B.LowerTriangular) (r : Vec3 ℝ) : B.wrap (B.wrap r) = B.wrap r :=
  wrap_eq_of h (wrap_inPositiveCell h r)
    (by rw [show B.wrap r - B.wrap r = 0 by ext <;> simp]; exact lattice_zero)

/-- **Not the shortest image**: in a skewed box the centred cell can hold an image for which a
lattice vector gives a shorter one. Box `a = (2, 0, 0)`, `b = (1, 1, 0)`, `c = (0, 0, 1)`: the
point `(6/5, 2/5, 0)` has fractional coordinates `(2/5, 2/5, 0)`, so it is its own image, of
squared length `8/5`, while its image `r − a = (−4/5, 2/5, 0)` has squared length `4/5`. -/
theorem image_not_shortest :
    ∃ (B : Box) (r w : Vec3 ℝ), B.LowerTriangular ∧ B.Lattice w ∧
      normSq (B.image r - w) < normSq (B.image r) := by
  let B : Box := ⟨⟨2, 0, 0⟩, ⟨1, 1, 0⟩, ⟨0, 0, 1⟩⟩
  have hB : B.LowerTriangular := ⟨rfl, rfl, rfl, by norm_num, by norm_num, by norm_num⟩
  have hr : B.InCell ⟨6 / 5, 2 / 5, 0⟩ := by
    refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩⟩ <;> simp [B, frac] <;> norm_num
  refine ⟨B, ⟨6 / 5, 2 / 5, 0⟩, B.a, hB, lattice_a, ?_⟩
  rw [image_of_inCell hB hr]
  simp [B, normSq]
  norm_num

/-! ## Below half the smallest height both reductions agree

The **height** of the box along `a` is the distance between the two faces spanned by `b` and `c`,
`h_a = V / |b × c|` (likewise `h_b = V / |c × a|`, `h_c = V / |a × b|`). A fractional coordinate
is the dot product with a reciprocal vector, `s_a = u_a · r` with `u_a = (b × c) / V`, so
`|s_a| ≤ |r| / h_a`. A vector shorter than half the smallest height therefore has every
fractional coordinate inside `(−1/2, 1/2)`: it is its own centred-cell image
(`image_shortest_of_lt_half_height`), and every other image is longer. In a lower-triangular box
it also lies in the brick `|x| < aₓ/2, |y| < b_y/2, |z| < c_z/2` (`h_a ≤ aₓ`, `h_b ≤ b_y`,
`h_c = c_z`), the cell gromosXX's triclinic `nearest_image` reduces into (it reads
`Box(i, j)` as component `i` of edge `j`, so its skew terms `b_z`, `a_y`, `a_z` vanish); so the
two reductions return the same, shortest, image for every such pair.
-/

/-- The dot product. -/
def dot (u v : Vec3 ℝ) : ℝ := u.x * v.x + u.y * v.y + u.z * v.z

/-- Cauchy–Schwarz for three components. -/
theorem dot_sq_le (u v : Vec3 ℝ) : dot u v ^ 2 ≤ normSq u * normSq v := by
  simp only [dot, normSq]
  nlinarith [sq_nonneg (u.y * v.z - u.z * v.y), sq_nonneg (u.z * v.x - u.x * v.z),
    sq_nonneg (u.x * v.y - u.y * v.x)]

/-- The cross product. -/
def cross (u v : Vec3 ℝ) : Vec3 ℝ :=
  ⟨u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x⟩

/-- The reciprocal vector of `a` in a lower-triangular box: `(Box.frac r).x = recA · r`. -/
noncomputable def recA (B : Box) : Vec3 ℝ :=
  ⟨1 / B.a.x, -B.b.x / (B.a.x * B.b.y),
    (B.b.x * B.c.y - B.b.y * B.c.x) / (B.a.x * B.b.y * B.c.z)⟩

/-- The reciprocal vector of `b`: `(Box.frac r).y = recB · r`. -/
noncomputable def recB (B : Box) : Vec3 ℝ := ⟨0, 1 / B.b.y, -B.c.y / (B.b.y * B.c.z)⟩

/-- The reciprocal vector of `c`: `(Box.frac r).z = recC · r`. -/
noncomputable def recC (B : Box) : Vec3 ℝ := ⟨0, 0, 1 / B.c.z⟩

/-- The heights of the box: `1 / |u|` for each reciprocal vector `u` (`= V / |b × c|`, …). -/
noncomputable def heightA (B : Box) : ℝ := (Real.sqrt (normSq B.recA))⁻¹
noncomputable def heightB (B : Box) : ℝ := (Real.sqrt (normSq B.recB))⁻¹
noncomputable def heightC (B : Box) : ℝ := (Real.sqrt (normSq B.recC))⁻¹

/-- The smallest height. -/
noncomputable def minHeight (B : Box) : ℝ := min (min B.heightA B.heightB) B.heightC

/-- `u_a = (b × c) / V`, `u_b = (c × a) / V`, `u_c = (a × b) / V`, with `V = aₓ b_y c_z`. -/
theorem rec_eq_cross (h : B.LowerTriangular) :
    B.recA = (1 / (B.a.x * B.b.y * B.c.z)) • cross B.b B.c ∧
      B.recB = (1 / (B.a.x * B.b.y * B.c.z)) • cross B.c B.a ∧
      B.recC = (1 / (B.a.x * B.b.y * B.c.z)) • cross B.a B.b := by
  have ha := h.ax_pos.ne'
  have hb := h.by_pos.ne'
  have hc := h.cz_pos.ne'
  refine ⟨?_, ?_, ?_⟩ <;> ext <;> simp [recA, recB, recC, cross, h.ay, h.az, h.bz] <;>
    field_simp

theorem frac_eq_dot (h : B.LowerTriangular) (r : Vec3 ℝ) :
    B.frac r = ⟨dot B.recA r, dot B.recB r, dot B.recC r⟩ := by
  have ha := h.ax_pos.ne'
  have hb := h.by_pos.ne'
  have hc := h.cz_pos.ne'
  ext <;> simp [frac, dot, recA, recB, recC] <;> field_simp <;> ring

theorem normSq_nonneg (u : Vec3 ℝ) : 0 ≤ normSq u := by
  simp only [normSq]; positivity

/-- `|r| < h / 2` with `h = 1 / |u|`, squared: `4 |r|² |u|² < 1`. -/
theorem sq_of_lt_half_height {u r : Vec3 ℝ} (hu : 0 < normSq u)
    (hr : Real.sqrt (normSq r) < (Real.sqrt (normSq u))⁻¹ / 2) :
    4 * (normSq r * normSq u) < 1 := by
  have su := Real.sqrt_pos.mpr hu
  have hr' : Real.sqrt (normSq r) * Real.sqrt (normSq u) < 1 / 2 := by
    have m := mul_lt_mul_of_pos_right hr su
    have e : (Real.sqrt (normSq u))⁻¹ / 2 * Real.sqrt (normSq u) = 1 / 2 := by
      field_simp
    linarith
  have e : normSq r * normSq u = (Real.sqrt (normSq r) * Real.sqrt (normSq u)) ^ 2 := by
    rw [mul_pow, Real.sq_sqrt (normSq_nonneg r), Real.sq_sqrt (normSq_nonneg u)]
  rw [e]
  have h0 : 0 ≤ Real.sqrt (normSq r) * Real.sqrt (normSq u) := by positivity
  nlinarith [hr']

theorem recA_pos (h : B.LowerTriangular) : 0 < normSq B.recA := by
  have : 0 < (1 / B.a.x) ^ 2 := by have := h.ax_pos; positivity
  simp only [normSq, recA]; nlinarith [sq_nonneg (-B.b.x / (B.a.x * B.b.y)),
    sq_nonneg ((B.b.x * B.c.y - B.b.y * B.c.x) / (B.a.x * B.b.y * B.c.z))]

theorem recB_pos (h : B.LowerTriangular) : 0 < normSq B.recB := by
  have : 0 < (1 / B.b.y) ^ 2 := by have := h.by_pos; positivity
  simp only [normSq, recB]; nlinarith [sq_nonneg (-B.c.y / (B.b.y * B.c.z))]

theorem recC_pos (h : B.LowerTriangular) : 0 < normSq B.recC := by
  have : 0 < (1 / B.c.z) ^ 2 := by have := h.cz_pos; positivity
  simp only [normSq, recC]; nlinarith

/-- A fractional coordinate `u · r` with `4 |r|² |u|² < 1` lies in `(−1/2, 1/2)`. -/
theorem dot_mem_of_sq {u r : Vec3 ℝ} (hs : 4 * (normSq r * normSq u) < 1) :
    dot u r ∈ Set.Ioc (-(1 / 2)) (1 / 2) := by
  have cs := dot_sq_le u r
  have : dot u r ^ 2 < 1 / 4 := by nlinarith
  constructor <;> nlinarith [sq_nonneg (dot u r - 1 / 2), sq_nonneg (dot u r + 1 / 2)]

/-- A non-zero lattice vector `w` is at least one height long along some reciprocal vector:
`(u · w)² ≥ 1` for `u` one of `u_a, u_b, u_c`. -/
theorem lattice_dot_ge (h : B.LowerTriangular) {w : Vec3 ℝ} (hw : B.Lattice w) (hw0 : w ≠ 0) :
    1 ≤ dot B.recA w ^ 2 ∨ 1 ≤ dot B.recB w ^ 2 ∨ 1 ≤ dot B.recC w ^ 2 := by
  obtain ⟨i, j, k, hijk⟩ := frac_of_lattice h hw
  rw [frac_eq_dot h] at hijk
  have hx : dot B.recA w = i := congrArg Vec3.x hijk
  have hy : dot B.recB w = j := congrArg Vec3.y hijk
  have hz : dot B.recC w = k := congrArg Vec3.z hijk
  have one_le : ∀ n : ℤ, n ≠ 0 → (1 : ℝ) ≤ (n : ℝ) ^ 2 := by
    intro n hn
    have : (1 : ℤ) ≤ n ^ 2 := by nlinarith [sq_nonneg n, Int.one_le_abs hn, sq_abs n]
    exact_mod_cast this
  by_cases hi : i = 0
  · by_cases hj : j = 0
    · by_cases hk : k = 0
      · exfalso; apply hw0
        rw [← point_frac h w, frac_eq_dot h, hx, hy, hz, hi, hj, hk]
        ext <;> simp [point]
      · exact Or.inr (Or.inr (hz ▸ one_le k hk))
    · exact Or.inr (Or.inl (hy ▸ one_le j hj))
  · exact Or.inl (hx ▸ one_le i hi)

/-- If `4 |d|² |u|² < 1` and `(u · w)² ≥ 1`, then `d − w` is longer than `d`. -/
theorem longer_of {u d w : Vec3 ℝ} (hd : 4 * (normSq d * normSq u) < 1) (hw : 1 ≤ dot u w ^ 2) :
    normSq d < normSq (d - w) := by
  have csw := dot_sq_le u w
  have csd := dot_sq_le d w
  have hu : 0 < normSq u := by
    by_contra hu; push Not at hu
    have := normSq_nonneg u
    have : normSq u = 0 := le_antisymm hu this
    nlinarith [normSq_nonneg w]
  -- |w|² ≥ 1/|u|² > 4|d|²
  have hw2 : 4 * normSq d < normSq w := by
    have h1 : 1 ≤ normSq u * normSq w := le_trans hw csw
    nlinarith [normSq_nonneg d, normSq_nonneg w]
  -- (d·w)² ≤ |d|²|w|² < |w|⁴/4, so 2 d·w < |w|²
  have h2 : 2 * dot d w < normSq w := by
    have hwpos : 0 < normSq w := by nlinarith [normSq_nonneg d]
    have : (2 * dot d w) ^ 2 < normSq w ^ 2 := by nlinarith [normSq_nonneg d]
    nlinarith [sq_nonneg (2 * dot d w + normSq w), abs_lt_of_sq_lt_sq' this hwpos.le]
  have e : normSq (d - w) = normSq d - 2 * dot d w + normSq w := by
    simp only [normSq, dot, sub_x, sub_y, sub_z]; ring
  rw [e]; linarith

/-- **Below half the smallest height the image is the vector itself, and the unique shortest
one**: in a lower-triangular box, `|d| < h_min / 2` gives `image d = d`, and every other image
`d − w` (`w` a non-zero lattice vector) is strictly longer. gromosXX's brick reduction returns
the same vector (see the section note), so the two agree on every such pair. -/
theorem image_shortest_of_lt_half_height (h : B.LowerTriangular) {d : Vec3 ℝ}
    (hd : Real.sqrt (normSq d) < B.minHeight / 2) :
    B.image d = d ∧ ∀ w : Vec3 ℝ, B.Lattice w → w ≠ 0 → normSq d < normSq (d - w) := by
  have hA : Real.sqrt (normSq d) < B.heightA / 2 := by
    have : B.minHeight ≤ B.heightA := le_trans (min_le_left _ _) (min_le_left _ _)
    linarith
  have hB : Real.sqrt (normSq d) < B.heightB / 2 := by
    have : B.minHeight ≤ B.heightB := le_trans (min_le_left _ _) (min_le_right _ _)
    linarith
  have hC : Real.sqrt (normSq d) < B.heightC / 2 := by
    have : B.minHeight ≤ B.heightC := min_le_right _ _
    linarith
  have sA := sq_of_lt_half_height (recA_pos h) hA
  have sB := sq_of_lt_half_height (recB_pos h) hB
  have sC := sq_of_lt_half_height (recC_pos h) hC
  refine ⟨image_of_inCell h ?_, fun w hw hw0 => ?_⟩
  · unfold InCell
    rw [frac_eq_dot h]
    exact ⟨dot_mem_of_sq sA, dot_mem_of_sq sB, dot_mem_of_sq sC⟩
  · rcases lattice_dot_ge h hw hw0 with h1 | h1 | h1
    · exact longer_of sA h1
    · exact longer_of sB h1
    · exact longer_of sC h1

end Box

end NewenTheory.Space.Periodicity
