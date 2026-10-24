import NewenSpec.P3M
import NewenTheory.Interactions.SPME

/-!
# Model of the SPME kernels of `crates/newen-forces` (NLRELE 5, newen only, experimental)

Rust (spec keys in `spec/SPEC.toml`), all in `src/nonbonded/spme.rs`:

* `spme.bspline`, `spme.bspline_and_deriv`: the cardinal B-spline `M_n` by its recursion and
  `M_n' = M_{n−1}(u) − M_{n−1}(u − 1)`;
* `spme.spline_weights`: the `n` weights and derivatives a charge gives its grid points, the
  recursion on the whole table at once;
* `spme.euler_spline_factor`, `spme.influence_at`: the Euler factor `|b(m)|²` and the SPME
  influence function `Ĝ = γ̂/k² · B` with its derivative tensor `Γ̂`;
* `spme.assign`, `spme.assign_squared_charges`, `spme.forces`: the grid charge density, the
  squared-charge grid of `Ã2`, and the forces from the splines' derivatives against the
  convolution grid.

Theory: `NewenTheory/Interactions/SPME.lean`.

Each kernel is written once over any number type `α`, operation by operation in the Rust order,
as `NewenSpec/P3M.lean` does: at `Float` it is what `lake exe spec-vectors` evaluates and
`crates/newen-forces/tests/spme_spec.rs` compares bit for bit; at `ℝ` the `*_real_eq` theorems
tie the per-point kernels to the theory. The transforms (PhastFT) are not modelled: the k-space
energy and the convolution are compared within a tolerance, against the Lean Ewald model and P3M.
-/

namespace NewenSpec.SPME

open NewenSpec.LatticeSum NewenSpec.P3M
open NewenTheory.Interactions.LatticeSum (Shape)

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α] [LE α]
  [DecidableRel (α := α) (· < ·)] [DecidableRel (α := α) (· ≤ ·)] [Sqrt α] [Exp α] [Erfc α]
  [Cos α] [Sin α] [Floor α] [ToInt α] [OfNat α 0] [OfNat α 1] [NatCast α] [IntCast α]

/-- `bspline(n, u)`. -/
def bsplineModel : ℕ → α → α
  | 0, _ => 0
  | 1, u => if 0 ≤ u ∧ u < 1 then 1 else 0
  | n + 2, u =>
    let k : α := ((n + 2 : ℕ) : α)
    let a := bsplineModel (n + 1) u
    let b := bsplineModel (n + 1) (u - 1)
    (u * a + (k - u) * b) / (k - 1)

/-- `bspline_and_deriv(n, u)`. -/
def bsplineDerivModel (n : ℕ) (u : α) : α × α :=
  (bsplineModel α n u, bsplineModel α (n - 1) u - bsplineModel α (n - 1) (u - 1))

/-- One order `k` of `spline_weights`' table on a list of `n` entries: per `j < k`,
`x = w + j`, `t_j = (x · t_j + (k − x) · t_{j−1}) / (k − 1)` (`t_{−1} = 0`), the old values on the
right (the Rust loop runs `j` downwards in place). -/
def tableStep (w : α) (n k : ℕ) (t : List α) : List α :=
  let kf : α := (k : α)
  let div := kf - 1
  (List.range n).map fun j =>
    if j < k then
      let x := w + (j : α)
      let below := if 0 < j then t.getD (j - 1) 0 else 0
      (x * t.getD j 0 + (kf - x) * below) / div
    else t.getD j 0

/-- The table after the orders `2 … r + 1`: `(1, 0, …, 0)` for `r = 0`. -/
def tableModel (w : α) (n : ℕ) : ℕ → List α
  | 0 => (List.range n).map fun j => if j = 0 then 1 else 0
  | r + 1 => tableStep α w n (r + 2) (tableModel w n r)

/-- `dθ_0 = t_0`, `dθ_j = t_j − t_{j−1}`. -/
def diffModel (n : ℕ) (t : List α) : List α :=
  (List.range n).map fun j => if j = 0 then t.getD 0 0 else t.getD j 0 - t.getD (j - 1) 0

/-- `spline_weights(n, u)`: `(θ, dθ)` with `w = u − floor u`. -/
def splineWeightsModel (n : ℕ) (u : α) : List α × List α :=
  let w := u - Floor.floor u
  (tableModel α w n (n - 1), diffModel α n (tableModel α w n (n - 2)))

/-- `euler_spline_factor(n, k_alpha, m)`, with `π` as `pi`. -/
def eulerModel (pi : α) (n K : ℕ) (m : ℤ) : α :=
  let twoPi := 2 * pi
  let term := fun (j : ℕ) => (bsplineModel α n ((j + 1 : ℕ) : α),
    twoPi * (m : α) * (j : α) / (K : α))
  let re := (List.range (n - 1)).foldl (fun acc j =>
    let t := term j
    acc + t.1 * Cos.cos t.2) 0
  let im := (List.range (n - 1)).foldl (fun acc j =>
    let t := term j
    acc + t.1 * Sin.sin t.2) 0
  1 / (re * re + im * im)

/-- `influence_from(k, b, shape, width)`: `(Ĝ, Γ̂)`. -/
def influenceFromModel (k : V3 α) (b : α) (s : Shape) (width : α) : α × Sym α :=
  let k2 := k.x * k.x + k.y * k.y + k.z * k.z
  let k2i := 1 / k2
  let ak := width * Sqrt.sqrt k2
  let f := shapeFourierModel α s ak
  let g := k2i * f.1 * b
  let sc := (ak * f.2 - 2 * f.1) * (k2i * k2i) * b
  (g, symTensor ⟨sc * k.x, sc * k.y, sc * k.z⟩ k)

/-- `spme_influence_at(l, grid, n, shape, width, reciprocal)`. -/
def influenceAtModel (pi : α) (l : ℤ × ℤ × ℤ) (g : ℕ × ℕ × ℕ) (n : ℕ) (s : Shape) (width : α)
    (recip : M3 α) : α × Sym α :=
  let k := latticeVectorModel α recip.r0 recip.r1 recip.r2 l
  let b := eulerModel α pi n g.1 l.1 * eulerModel α pi n g.2.1 l.2.1 *
    eulerModel α pi n g.2.2 l.2.2
  influenceFromModel α k b s width

/-- The base grid point `floor u` and the weights and derivatives of one atom along each axis. -/
def atomSplines (m : Mesh α) (n : ℕ) (r : V3 α) :
    (ℤ × ℤ × ℤ) × (List α × List α) × (List α × List α) × (List α × List α) :=
  let u := m.hInv.tmulVec r
  ((ToInt.toInt (Floor.floor u.x), ToInt.toInt (Floor.floor u.y), ToInt.toInt (Floor.floor u.z)),
    splineWeightsModel α n u.x, splineWeightsModel α n u.y, splineWeightsModel α n u.z)

/-- `spme_assign(positions, charges, mesh, n)`. -/
def assignModel (m : Mesh α) (n : ℕ) (atoms : List (V3 α × α)) : Array α :=
  atoms.foldl (fun (rho : Array α) a =>
    let sp := atomSplines α m n a.1
    let (bx, by_, bz) := sp.1
    let wx := sp.2.1.1
    let wy := sp.2.2.1.1
    let wz := sp.2.2.2.1
    (List.range n).foldl (fun rho (jx : ℕ) =>
      let px := bx - (jx : ℤ)
      (List.range n).foldl (fun rho (jy : ℕ) =>
        let py := by_ - (jy : ℤ)
        (List.range n).foldl (fun rho (jz : ℕ) =>
          let pz := bz - (jz : ℤ)
          let w := wx.getD jx 0 * wy.getD jy 0 * wz.getD jz 0 / m.cellVolume
          rho.modify (gridIndex m.grid px py pz) (· + w * a.2)) rho) rho) rho)
    (Array.replicate m.len 0)

/-- `overlaps(theta)`: `f(d) = Σ_j θ_j θ_{j−d}`, `d = −(n − 1) … n − 1`. -/
def overlapsModel (theta : List α) : List α :=
  let n := theta.length
  (fullStencil n).map fun d =>
    (List.range n).foldl (fun s (j : ℕ) =>
      let i : ℤ := (j : ℤ) - d
      if 0 ≤ i ∧ i < n then s + theta.getD j 0 * theta.getD i.toNat 0 else s) 0

/-- `spme_assign_squared_charges(positions, charges, mesh, n)`. -/
def assignSquaredModel (m : Mesh α) (n : ℕ) (atoms : List (V3 α × α)) : Array α :=
  let cvi := 1 / m.cellVolume
  atoms.foldl (fun (grid : Array α) a =>
    let q2 := a.2 * a.2
    let sp := atomSplines α m n a.1
    let fx := overlapsModel α sp.2.1.1
    let fy := overlapsModel α sp.2.2.1.1
    let fz := overlapsModel α sp.2.2.2.1
    let b : ℤ := ((n - 1 : ℕ) : ℤ)
    (fullStencil n).foldl (fun grid (dx : ℤ) =>
      let wx := fx.getD (dx + b).toNat 0
      (fullStencil n).foldl (fun grid (dy : ℤ) =>
        let wy := fy.getD (dy + b).toNat 0
        (fullStencil n).foldl (fun grid (dz : ℤ) =>
          let wz := fz.getD (dz + b).toNat 0
          let w := wx * wy * wz * cvi
          grid.modify (gridIndex m.grid dx dy dz) (· + w * q2)) grid) grid) grid)
    (Array.replicate m.len 0)

/-- `spme_forces` for one atom `(r, q)` and the convolution grid `conv`. -/
def forceModel (m : Mesh α) (n : ℕ) (conv : Array α) (eps0Inv : α) (a : V3 α × α) : V3 α :=
  let sp := atomSplines α m n a.1
  let (bx, by_, bz) := sp.1
  let (wx, dx) := sp.2.1
  let (wy, dy) := sp.2.2.1
  let (wz, dz) := sp.2.2.2
  let g := (List.range n).foldl (fun g (jx : ℕ) =>
    let px := bx - (jx : ℤ)
    (List.range n).foldl (fun g (jy : ℕ) =>
      let py := by_ - (jy : ℤ)
      (List.range n).foldl (fun (g : V3 α) (jz : ℕ) =>
        let pz := bz - (jz : ℤ)
        let c := conv.getD (gridIndex m.grid px py pz) 0
        ⟨g.x + dx.getD jx 0 * wy.getD jy 0 * wz.getD jz 0 * c,
          g.y + wx.getD jx 0 * dy.getD jy 0 * wz.getD jz 0 * c,
          g.z + wx.getD jx 0 * wy.getD jy 0 * dz.getD jz 0 * c⟩) g) g) V3.zero
  let grad := m.hInv.mulVec g
  let s := -(a.2 * eps0Inv) / (m.len : α)
  ⟨s * grad.x, s * grad.y, s * grad.z⟩

end Model

/-! ## The `Float` instantiation (what `lake exe spec-vectors` runs) -/

def bsplineFloat (n : ℕ) (u : Float) : Float := bsplineModel Float n u
def bsplineDerivFloat (n : ℕ) (u : Float) : Float × Float := bsplineDerivModel Float n u
def splineWeightsFloat (n : ℕ) (u : Float) : List Float × List Float :=
  splineWeightsModel Float n u
def eulerFloat (n K : ℕ) (m : ℤ) : Float := eulerModel Float piFloat n K m
def influenceAtFloat (l : ℤ × ℤ × ℤ) (g : ℕ × ℕ × ℕ) (n : ℕ) (s : Shape) (width : Float)
    (recip : M3 Float) : Float × Sym Float := influenceAtModel Float piFloat l g n s width recip
def assignFloat (m : Mesh Float) (n : ℕ) (atoms : List (V3 Float × Float)) : Array Float :=
  assignModel Float m n atoms
def assignSquaredFloat (m : Mesh Float) (n : ℕ) (atoms : List (V3 Float × Float)) :
    Array Float := assignSquaredModel Float m n atoms
def forceFloat (m : Mesh Float) (n : ℕ) (conv : Array Float) (eps0Inv : Float)
    (a : V3 Float × Float) : V3 Float := forceModel Float m n conv eps0Inv a

/-- `spme_influence_function` in its loop order, with the Euler factors of every storage
position: `(Ĝ, Γ̂)`. -/
def influenceFunctionFloat (g : ℕ × ℕ × ℕ) (n : ℕ) (s : Shape) (width : Float)
    (recip : M3 Float) : Array Float × Array (Sym Float) :=
  let ex := (List.range g.1).map fun i => eulerFloat n g.1 (waveIndex i g.1)
  let ey := (List.range g.2.1).map fun i => eulerFloat n g.2.1 (waveIndex i g.2.1)
  let ez := (List.range g.2.2).map fun i => eulerFloat n g.2.2 (waveIndex i g.2.2)
  let pts := (List.range g.1).flatMap fun x => (List.range g.2.1).flatMap fun y =>
    (List.range g.2.2).map fun z => (x, y, z)
  let vals := pts.map fun (x, y, z) =>
    let l := (waveIndex x g.1, waveIndex y g.2.1, waveIndex z g.2.2)
    if l.1 = 0 ∧ l.2.1 = 0 ∧ l.2.2 = 0 then ((0 : Float), (Sym.zero : Sym Float))
    else
      let k := latticeVectorModel Float recip.r0 recip.r1 recip.r2 l
      let b := ex.getD x 0 * ey.getD y 0 * ez.getD z 0
      influenceFromModel Float k b s width
  (vals.toArray.map (·.1), vals.toArray.map (·.2))

/-! ## The `ℝ` instantiation is the theory -/

section Real

open NewenTheory.Interactions.SPME
open NewenTheory.Interactions.LatticeSum (gammaHat)

/-- `bspline` over ℝ is the theory's `M_n`, for every order and every `u`. -/
theorem bsplineModel_real_eq : ∀ (n : ℕ) (u : ℝ), bsplineModel ℝ n u = bspline n u
  | 0, _ => rfl
  | 1, u => by simp only [bsplineModel, bspline_one]
  | n + 2, u => by
    rw [bsplineModel, bspline_succ_succ, bsplineModel_real_eq (n + 1) u,
      bsplineModel_real_eq (n + 1) (u - 1)]
    have h : ((n + 2 : ℕ) : ℝ) - 1 = (n : ℝ) + 1 := by push_cast; ring
    have hn : (n : ℝ) + 1 ≠ 0 := by positivity
    simp only [h]
    push_cast
    field_simp

/-- `bspline_and_deriv` over ℝ: `(M_n(u), M_n'(u))` with the theory's `bsplineDeriv`. -/
theorem bsplineDerivModel_real_eq (n : ℕ) (u : ℝ) :
    bsplineDerivModel ℝ n u = (bspline n u, bsplineDeriv n u) := by
  simp [bsplineDerivModel, bsplineModel_real_eq, bsplineDeriv]

theorem getD_map_range {β : Type} (n j : ℕ) (f : ℕ → β) (d : β) (h : j < n) :
    ((List.range n).map f).getD j d = f j := by
  simp [List.getD_eq_getElem?_getD, h]

theorem length_tableModel (w : ℝ) (n : ℕ) : ∀ r, (tableModel ℝ w n r).length = n
  | 0 => by simp [tableModel]
  | r + 1 => by simp [tableModel, tableStep]

/-- **The table is the recursion**: after the orders up to `r + 1 ≤ n`, entry `j` is
`M_{r+1}(w + j)` for `j ≤ r` and `0` above, for `w ∈ [0, 1)`. -/
theorem tableModel_real_eq {w : ℝ} (h0 : 0 ≤ w) (h1 : w < 1) (n : ℕ) :
    ∀ r, r + 1 ≤ n → ∀ j < n,
      (tableModel ℝ w n r).getD j 0 = if j ≤ r then bspline (r + 1) (w + j) else 0
  | 0, _, j, hj => by
    rw [tableModel, getD_map_range n j _ 0 hj]
    rcases Nat.eq_zero_or_pos j with rfl | hpos
    · simp [bspline_one, h0, h1]
    · simp [show j ≠ 0 by omega, show ¬ j ≤ 0 by omega]
  | r + 1, hr, j, hj => by
    have ih := tableModel_real_eq h0 h1 n r (by omega)
    rw [tableModel, tableStep, getD_map_range n j _ 0 hj]
    dsimp only
    by_cases hjk : j < r + 2
    · simp only [hjk, show j ≤ r + 1 by omega, ↓reduceIte]
      -- the old entry `j` and the one below it, as values of `M_{r+1}`
      have hold : (tableModel ℝ w n r).getD j 0 = bspline (r + 1) (w + j) := by
        rw [ih j hj]
        split_ifs with hle
        · rfl
        · have hrj : (r : ℝ) + 1 ≤ j := by exact_mod_cast (by omega : r + 1 ≤ j)
          exact (bspline_support (r + 1) (Or.inr (by push_cast; linarith))).symm
      have hbelow : (if 0 < j then (tableModel ℝ w n r).getD (j - 1) 0 else 0) =
          bspline (r + 1) (w + j - 1) := by
        split_ifs with hj0
        · rw [ih (j - 1) (by omega)]
          simp only [show j - 1 ≤ r by omega, ↓reduceIte]
          congr 1
          rw [Nat.cast_sub (by omega)]
          push_cast
          ring
        · rw [show j = 0 by omega]
          exact (bspline_support (r + 1) (Or.inl (by push_cast; linarith))).symm
      rw [hold, hbelow, show r + 1 + 1 = r + 2 by omega, bspline_succ_succ]
      have hr0 : (r : ℝ) + 1 ≠ 0 := by positivity
      have h2 : ((r + 2 : ℕ) : ℝ) - 1 = (r : ℝ) + 1 := by push_cast; ring
      rw [h2]
      push_cast
      field_simp
    · simp only [hjk, show ¬ j ≤ r + 1 by omega, ↓reduceIte]
      rw [ih j hj]
      simp only [show ¬ j ≤ r by omega, ↓reduceIte]

/-- **`spline_weights` over ℝ** (`n ≥ 2`): the weight `j` is `M_n(w + j)` and the derivative
`j` is `M_n'(w + j) = M_{n−1}(w + j) − M_{n−1}(w + j − 1)`, with `w = u − ⌊u⌋`. -/
theorem splineWeightsModel_real_eq {n : ℕ} (hn : 2 ≤ n) (u : ℝ) {j : ℕ} (hj : j < n) :
    (splineWeightsModel ℝ n u).1.getD j 0 = bspline n (Int.fract u + j) ∧
      (splineWeightsModel ℝ n u).2.getD j 0 = bsplineDeriv n (Int.fract u + j) := by
  have hw0 := Int.fract_nonneg u
  have hw1 := Int.fract_lt_one u
  have hf : u - Floor.floor u = Int.fract u := rfl
  simp only [splineWeightsModel, hf]
  constructor
  · rw [tableModel_real_eq hw0 hw1 n (n - 1) (by omega) j hj, show n - 1 + 1 = n by omega]
    simp only [show j ≤ n - 1 by omega, ↓reduceIte]
  · rw [diffModel, getD_map_range n j _ 0 hj, bsplineDeriv]
    have ht := tableModel_real_eq hw0 hw1 n (n - 2) (by omega)
    have hn1 : n - 2 + 1 = n - 1 := by omega
    -- entry `j` of the table of order `n − 1`, zero at `j = n − 1` (support)
    have hent : ∀ i < n, (tableModel ℝ (Int.fract u) n (n - 2)).getD i 0 =
        bspline (n - 1) (Int.fract u + i) := fun i hi => by
      rw [ht i hi, hn1]
      split_ifs with hle
      · rfl
      · refine (bspline_support (n - 1) (Or.inr ?_)).symm
        have : ((n - 1 : ℕ) : ℝ) ≤ i := by exact_mod_cast (by omega : n - 1 ≤ i)
        linarith
    split_ifs with hj0
    · subst hj0
      rw [hent 0 (by omega), bspline_support (n - 1) (u := Int.fract u + (0 : ℕ) - 1)
        (Or.inl (by push_cast; linarith))]
      ring
    · rw [hent j hj, hent (j - 1) (by omega), Nat.cast_sub (by omega)]
      push_cast
      ring_nf

/-- A fold that adds `f j` is the sum. -/
theorem foldl_add_eq_sum {β : Type} (l : List β) (f : β → ℝ) (a : ℝ) :
    l.foldl (fun acc j => acc + f j) a = a + (l.map f).sum := by
  induction l generalizing a with
  | nil => simp
  | cons b l ih => simp [ih]; ring

theorem sum_range_list (n : ℕ) (f : ℕ → ℝ) :
    ((List.range n).map f).sum = ∑ j ∈ Finset.range n, f j := by
  induction n with
  | zero => simp
  | succ n ih => simp [List.range_succ, Finset.sum_range_succ, ih]

/-- `euler_spline_factor` over ℝ is the theory's `|b(m)|²`. -/
theorem eulerModel_real_eq (n K : ℕ) (m : ℤ) :
    eulerModel ℝ Real.pi n K m = eulerFactor n K m := by
  rw [eulerFactor_eq_re_im, sq, sq]
  simp only [eulerModel, foldl_add_eq_sum, zero_add, sum_range_list, bsplineModel_real_eq,
    Cos.cos, Sin.sin]
  push_cast
  rfl

/-- `influence_from` over ℝ: `Ĝ = γ̂(a|k|)/|k|² · b` for `k ≠ 0` and `a ≠ 0`. -/
theorem influenceFromModel_real_eq (k : V3 ℝ) (b : ℝ) (s : Shape) {a : ℝ} (ha : a ≠ 0)
    (hk : 0 < k.x * k.x + k.y * k.y + k.z * k.z) :
    (influenceFromModel ℝ k b s a).1 =
      gammaHat s (a * Real.sqrt (k.x * k.x + k.y * k.y + k.z * k.z)) /
        (k.x * k.x + k.y * k.y + k.z * k.z) * b := by
  have hak : a * Real.sqrt (k.x * k.x + k.y * k.y + k.z * k.z) ≠ 0 :=
    mul_ne_zero ha (Real.sqrt_ne_zero'.mpr hk)
  simp only [influenceFromModel, Sqrt.sqrt, shapeFourierModel_real_eq s hak]
  field_simp

/-- `spme_influence_at` over ℝ: the Gaussian-type weight `γ̂(a|k_l|)/|k_l|²` times the theory's
Euler factors of the three axes. -/
theorem influenceAtModel_real_eq (l : ℤ × ℤ × ℤ) (g : ℕ × ℕ × ℕ) (n : ℕ) (s : Shape) {a : ℝ}
    (ha : a ≠ 0) (recip : M3 ℝ) (k : V3 ℝ)
    (hkl : k = latticeVectorModel ℝ recip.r0 recip.r1 recip.r2 l)
    (hk : 0 < k.x * k.x + k.y * k.y + k.z * k.z) :
    (influenceAtModel ℝ Real.pi l g n s a recip).1 =
      gammaHat s (a * Real.sqrt (k.x * k.x + k.y * k.y + k.z * k.z)) /
        (k.x * k.x + k.y * k.y + k.z * k.z) *
        (eulerFactor n g.1 l.1 * eulerFactor n g.2.1 l.2.1 * eulerFactor n g.2.2 l.2.2) := by
  simp only [influenceAtModel, ← hkl, eulerModel_real_eq]
  exact influenceFromModel_real_eq k _ s ha hk

end Real

end NewenSpec.SPME
