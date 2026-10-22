import NewenSpec.LatticeSum
import NewenTheory.Interactions.P3M

/-!
# Model of the P3M kernels of `crates/newen-forces` (NLRELE 3)

Rust (spec keys in `spec/SPEC.toml`), all in `src/nonbonded/p3m.rs`:

* `p3m.assignment_1d`, `p3m.nearest_grid_point`, `p3m.assignment_weights`: the assignment
  function `w_p`, the centre of its stencil and the stencil's weights (gromosXX
  `charge_assignment_1d`);
* `p3m.mesh`: the grid matrix `H`, `H⁻¹` and the cell volume (gromosXX `math::inverse`, `det`);
* `p3m.assign_charges`, `p3m.selfterm_fp`, `p3m.assign_squared_charges`: the grid charge density
  and the squared-charge grid of `Ã2`;
* `p3m.fd_coefficients`, `p3m.sinc_power`, `p3m.influence_at`, `p3m.quality_at`,
  `p3m.influence_box_correction`, `p3m.force_error`: the optimal influence function per wave
  vector, its quality, its first-order box correction and the rms force error;
* `p3m.kspace_energy`, `p3m.potential_hat`, `p3m.field_ik`, `p3m.field_fd`,
  `p3m.interpolate_forces`, `p3m.a2_tilde`: what is computed on the transformed grids.

Theory: `NewenTheory/Interactions/P3M.lean`.

Each kernel is written once over any number type `α`, operation by operation in the Rust order,
which is gromosXX's, as `NewenSpec/LatticeSum.lean` does: at `Float` it is what
`lake exe spec-vectors` evaluates and `crates/newen-forces/tests/p3m_spec.rs` compares bit for
bit; at `ℝ` the `*_real_eq` theorems tie it to the theory. Grids are `Array`s in the storage
order `z + N_z (y + N_y x)`; a grid point's integer coordinates come from `ToInt` (Rust
`as i64`, truncation).

The FFT is not modelled bit for bit: PhastFT's butterflies are not gromosXX's, and neither is
the Lean code's. `dftFloat` evaluates the discrete Fourier transform by its definition, `O(N²)`,
for small grids, and the Rust FFT is compared with it within a tolerance, as is the whole
k-space energy (`p3mKspaceFloat`).
-/

namespace NewenSpec.P3M

open NewenSpec.LatticeSum
open NewenTheory.Interactions.LatticeSum (Shape)

/-- Truncation toward zero to an integer (Rust `x as i64`): `Float.toInt64` at `Float`. -/
class ToInt (α : Type) where
  toInt : α → ℤ

instance : ToInt Float := ⟨fun x => x.toInt64.toInt⟩

/-- A 3 × 3 matrix by rows (`[[f64; 3]; 3]`). -/
structure M3 (α : Type) where
  r0 : V3 α
  r1 : V3 α
  r2 : V3 α

/-- The stencil offsets of order `p`: `−(p − 1)/2 ..= p/2`. -/
def stencil (p : ℕ) : List ℤ :=
  (List.range p).map fun i => (i : ℤ) - (((p - 1) / 2 : ℕ) : ℤ)

/-- `rem_euclid` of a grid coordinate and the storage index `z + N_z (y + N_y x)`. -/
def gridIndex (g : ℕ × ℕ × ℕ) (x y z : ℤ) : ℕ :=
  ((z % g.2.2) + g.2.2 * ((y % g.2.1) + g.2.1 * (x % g.1))).toNat

/-- `wave_index(i, n)`. -/
def waveIndex (i n : ℕ) : ℤ := if n / 2 < i then (i : ℤ) - n else i

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α] [LE α]
  [DecidableRel (α := α) (· < ·)] [DecidableRel (α := α) (· ≤ ·)] [Sqrt α] [Exp α] [Erfc α]
  [Cos α] [Sin α] [Fabs α] [Powf α] [Rint α] [Floor α] [ToInt α] [OfNat α 0] [OfNat α 1]
  [NatCast α] [IntCast α]

/-- `assignment_1d(order, xi)`: `w_p(ξ)`. -/
def assignmentModel (p : ℕ) (xi : α) : α :=
  let absf := Fabs.fabs xi
  match p with
  | 1 => if absf < 1 / 2 then 1 else 0
  | 2 => if absf < 1 then 1 - absf else 0
  | 3 =>
    if absf < 1 / 2 then 3 / 4 - xi * xi
    else if absf < 3 / 2 then
      let term := 2 * absf - 3
      1 / 8 * term * term
    else 0
  | 4 =>
    if absf < 1 then
      let absf2 := absf * absf
      1 / 6 * (3 * absf2 * absf - 6 * absf2 + 4)
    else if absf < 2 then
      let term := absf - 2; (-1 / 6) * term * term * term
    else 0
  | 5 =>
    if absf < 1 / 2 then
      let xi2 := xi * xi
      let xi4 := xi2 * xi2
      1 / 192 * (48 * xi4 - 120 * xi2 + 115)
    else if absf < 3 / 2 then
      let xi2 := xi * xi
      let xi4 := xi2 * xi2
      let absf3 := absf * absf * absf; (-1 / 96) * (16 * xi4 - 80 * absf3 + 120 * xi2 - 20 * absf - 55)
    else if absf < 5 / 2 then
      let tmp := 2 * absf - 5
      let tmp2 := tmp * tmp
      1 / 384 * tmp2 * tmp2
    else 0
  | _ => 0

/-- `nearest_grid_point(order, u)` as a number: `rint u` for odd `p`, `u` truncated toward zero
for even `p`. -/
def nearestModel (p : ℕ) (u : α) : α :=
  if p % 2 = 1 then Rint.rint u
  else if u < 0 then -Floor.floor (-u) else Floor.floor u

/-- `nearest_grid_point` as an integer. -/
def nearestIntModel (p : ℕ) (u : α) : ℤ := ToInt.toInt (nearestModel α p u)

/-- `assignment_weights(order, u)`: `w_p(u − (n₀ + d))` per stencil offset `d`. -/
def weightsModel (p : ℕ) (u : α) : List α :=
  let n0 := nearestIntModel α p u
  (stencil p).map fun d => assignmentModel α p (u - ((n0 + d : ℤ) : α))

/-- `selfterm_fp(order, n, s)` for `|n|`. -/
def selftermFpModel (p n : ℕ) (s : α) : α :=
  let s2 := s * s
  match p, n with
  | 1, _ => 1
  | 2, 0 => 1 / 2 * (1 + 4 * s2)
  | 2, 1 => 1 / 4 * (1 - 4 * s2)
  | 3, 0 => 1 / 32 * (19 + s2 * (-24 + s2 * 48))
  | 3, 1 => 1 / 16 * (3 + s2 * (8 - s2 * 16))
  | 3, 2 =>
    let term := 1 - 4 * s2
    1 / 64 * term * term
  | 4, 0 => 1 / 576 * (265 + s2 * (204 + s2 * (-528 + s2 * 320)))
  | 4, 1 => 1 / 2304 * (575 + s2 * (-564 + s2 * (1488 + s2 * (-960))))
  | 4, 2 => 1 / 1152 * (23 + s2 * (84 + s2 * (-240 + s2 * 192)))
  | 4, 3 =>
    let term := 1 - 4 * s2
    1 / 2304 * term * term * term
  | 5, 0 => 1 / 73728 * (32227 + s2 * (-9520 + s2 * (28960 + s2 * (-29440 + s2 * 8960))))
  | 5, 1 => 1 / 18432 * (4389 + s2 * (1792 + s2 * (-5472 + s2 * (5632 + s2 * (-1792)))))
  | 5, 2 => 1 / 36864 * (1559 + s2 * (-1456 + s2 * (4512 + s2 * (-4864 + s2 * 1792))))
  | 5, 3 => 1 / 18432 * (19 + s2 * (128 + s2 * (-416 + s2 * (512 + s2 * (-256)))))
  | 5, 4 =>
    let term := 1 - 4 * s2
    let term2 := term * term
    1 / 147456 * term2 * term2
  | _, _ => 0

/-- `fd_coefficients(nfdord)`: five entries, unused ones zero. -/
def fdCoefficientsModel (q : ℕ) : List α :=
  match q with
  | 1 => [1, 0, 0, 0, 0]
  | 2 => [4 / 3, -1 / 3, 0, 0, 0]
  | 3 => [3 / 2, -3 / 5, 1 / 10, 0, 0]
  | 4 => [8 / 5, -4 / 5, 8 / 35, -1 / 35, 0]
  | 5 => [5 / 3, -20 / 21, 5 / 14, -5 / 63, 1 / 126]
  | _ => [0, 0, 0, 0, 0]

/-- `sinc_power(x, p)`: `1` at `x = 0` (written `¬(x < 0) ∧ ¬(0 < x)`, the same for every double
but NaN), else `pow(2/x · sin(x/2), p)`. -/
def sincPowModel (x : α) (p : ℕ) : α :=
  if ¬(x < 0) ∧ ¬(0 < x) then 1 else Powf.powf (2 / x * Sin.sin (x / 2)) (p : α)

variable {α}

/-- `det3(m)`. -/
def M3.det (m : M3 α) : α :=
  m.r0.y * m.r1.z * m.r2.x - m.r0.z * m.r1.y * m.r2.x + m.r0.z * m.r1.x * m.r2.y -
    m.r0.x * m.r1.z * m.r2.y - m.r0.y * m.r1.x * m.r2.z + m.r0.x * m.r1.y * m.r2.z

/-- `inverse3(m)`. -/
def M3.inv (m : M3 α) : M3 α :=
  let d := 1 / m.det
  ⟨⟨(m.r1.y * m.r2.z - m.r1.z * m.r2.y) * d, (m.r0.z * m.r2.y - m.r0.y * m.r2.z) * d,
      (m.r0.y * m.r1.z - m.r0.z * m.r1.y) * d⟩,
    ⟨(m.r1.z * m.r2.x - m.r1.x * m.r2.z) * d, (m.r0.x * m.r2.z - m.r0.z * m.r2.x) * d,
      (m.r0.z * m.r1.x - m.r0.x * m.r1.z) * d⟩,
    ⟨(m.r1.x * m.r2.y - m.r1.y * m.r2.x) * d, (m.r0.y * m.r2.x - m.r0.x * m.r2.y) * d,
      (m.r0.x * m.r1.y - m.r0.y * m.r1.x) * d⟩⟩

/-- `transpose3(m)`. -/
def M3.transpose (m : M3 α) : M3 α :=
  ⟨⟨m.r0.x, m.r1.x, m.r2.x⟩, ⟨m.r0.y, m.r1.y, m.r2.y⟩, ⟨m.r0.z, m.r1.z, m.r2.z⟩⟩

/-- `product(m, v)`: `m v`, row by row. -/
def M3.mulVec (m : M3 α) (v : V3 α) : V3 α := ⟨m.r0.dot v, m.r1.dot v, m.r2.dot v⟩

/-- `product_t(m, v)`: `mᵀ v`, gromosXX's `math::product(Matrix, Vec)`. -/
def M3.tmulVec (m : M3 α) (v : V3 α) : V3 α :=
  ⟨m.r0.x * v.x + m.r1.x * v.y + m.r2.x * v.z, m.r0.y * v.x + m.r1.y * v.y + m.r2.y * v.z,
    m.r0.z * v.x + m.r1.z * v.y + m.r2.z * v.z⟩

variable (α)

/-- A grid in a box: `(grid, H, H⁻¹, V_cell)`, `Mesh::new(edges, grid)`. -/
structure Mesh (α : Type) where
  grid : ℕ × ℕ × ℕ
  h : M3 α
  hInv : M3 α
  cellVolume : α

/-- `Mesh::new(edges, grid)`: rows `edges_i / N_i`. -/
def meshModel (edges : M3 α) (g : ℕ × ℕ × ℕ) : Mesh α :=
  let nx : α := (g.1 : α)
  let ny : α := (g.2.1 : α)
  let nz : α := (g.2.2 : α)
  let h : M3 α := ⟨⟨edges.r0.x / nx, edges.r0.y / nx, edges.r0.z / nx⟩,
    ⟨edges.r1.x / ny, edges.r1.y / ny, edges.r1.z / ny⟩,
    ⟨edges.r2.x / nz, edges.r2.y / nz, edges.r2.z / nz⟩⟩
  ⟨g, h, h.inv, h.transpose.det⟩

/-- The number of grid points. -/
def Mesh.len {α : Type} (m : Mesh α) : ℕ := m.grid.1 * m.grid.2.1 * m.grid.2.2

/-- `assign_charges(mesh, order, positions, charges)`. -/
def assignChargesModel (m : Mesh α) (p : ℕ) (atoms : List (V3 α × α)) : Array α :=
  atoms.foldl (fun (rho : Array α) a =>
    let u := m.hInv.tmulVec a.1
    let nx := nearestIntModel α p u.x
    let ny := nearestIntModel α p u.y
    let nz := nearestIntModel α p u.z
    (stencil p).foldl (fun rho dx =>
      let px := nx + dx
      let wx := assignmentModel α p (u.x - (px : α))
      (stencil p).foldl (fun rho dy =>
        let py := ny + dy
        let wy := assignmentModel α p (u.y - (py : α))
        (stencil p).foldl (fun rho dz =>
          let pz := nz + dz
          let wz := assignmentModel α p (u.z - (pz : α))
          let w := wx * wy * wz / m.cellVolume
          rho.modify (gridIndex m.grid px py pz) (· + w * a.2)) rho) rho) rho)
    (Array.replicate m.len 0)

/-- `−(p − 1) ..= p − 1`. -/
def fullStencil (p : ℕ) : List ℤ := (List.range (2 * p - 1)).map fun i => (i : ℤ) - (p - 1 : ℕ)

/-- `assign_squared_charges(mesh, order, positions, charges)`. -/
def assignSquaredModel (m : Mesh α) (p : ℕ) (atoms : List (V3 α × α)) : Array α :=
  let cvi := 1 / m.cellVolume
  atoms.foldl (fun (grid : Array α) a =>
    let q2 := a.2 * a.2
    let u := m.hInv.tmulVec a.1
    let off := fun (c : α) =>
      let n0 := nearestIntModel α p c
      if p % 2 = 1 then c - (n0 : α) else c - (n0 : α) - 1 / 2
    let sx := off u.x
    let sy := off u.y
    let sz := off u.z
    (fullStencil p).foldl (fun grid dx =>
      let wx := selftermFpModel α p dx.natAbs sx
      (fullStencil p).foldl (fun grid dy =>
        let wy := selftermFpModel α p dy.natAbs sy
        (fullStencil p).foldl (fun grid dz =>
          let wz := selftermFpModel α p dz.natAbs sz
          let w := wx * wy * wz * cvi
          grid.modify (gridIndex m.grid dx dy dz) (· + w * q2)) grid) grid) grid)
    (Array.replicate m.len 0)

/-- What the influence function of a grid needs (`InfluenceInputs`). -/
structure Inputs (α : Type) where
  order : ℕ
  nfdord : ℕ
  nalias : ℕ
  shape : Shape
  width : α
  mesh : Mesh α
  recip : M3 α

/-- `d_hat(inp, k_l)`: `(D̂, |D̂|²)`. -/
def dHatModel (inp : Inputs α) (k : V3 α) : V3 α × α :=
  let d :=
    if inp.nfdord = 0 then k
    else
      let c := fdCoefficientsModel α inp.nfdord
      let th := inp.mesh.h.mulVec k
      let s := (List.range inp.nfdord).foldl (fun (d : V3 α) i =>
        let jf : α := ((i + 1 : ℕ) : α)
        let f := c.getD i 0 / jf
        ⟨d.x + Sin.sin (jf * th.x) * f, d.y + Sin.sin (jf * th.y) * f,
          d.z + Sin.sin (jf * th.z) * f⟩) V3.zero
      inp.mesh.hInv.mulVec s
  (d, d.x * d.x + d.y * d.y + d.z * d.z)

/-- The alias sums: `(numerator, denominator, Σγ̂²/k², Γ̂ numerator)`. -/
structure AliasSums (α : Type) where
  num : V3 α
  den : α
  k2if2 : α
  gnum : Sym α

/-- One alias vector's terms added to the sums (the body of `alias_sums`' loop). -/
def aliasStep (pi : α) (inp : Inputs α) (kl d : V3 α) (abs2d : α) (withGamma : Bool)
    (acc : AliasSums α) (m : ℤ × ℤ × ℤ) : AliasSums α :=
  let twoPi := 2 * pi
  let hm := inp.mesh.hInv.mulVec ⟨(m.1 : α), (m.2.1 : α), (m.2.2 : α)⟩
  let k : V3 α := ⟨kl.x + twoPi * hm.x, kl.y + twoPi * hm.y, kl.z + twoPi * hm.z⟩
  let k2 := k.x * k.x + k.y * k.y + k.z * k.z
  let k2i := 1 / k2
  let ak := inp.width * Sqrt.sqrt k2
  let f := shapeFourierModel α inp.shape ak
  let th := inp.mesh.h.mulVec k
  let ph := sincPowModel α th.x inp.order * sincPowModel α th.y inp.order *
    sincPowModel α th.z inp.order
  let ph2 := ph * ph
  let k2iF := k2i * f.1
  let sv := k2iF * ph2
  let num : V3 α := ⟨acc.num.x + k.x * sv, acc.num.y + k.y * sv, acc.num.z + k.z * sv⟩
  let gnum :=
    if withGamma then
      let kd := k.x * d.x + k.y * d.y + k.z * d.z
      let kxd := symTensor k d
      let dxk := symTensor d k
      let kxk := symTensor k k
      let dxd := symTensor d d
      let d2i := 1 / abs2d
      let twoKd := 2 * kd
      let pre := ph2 * k2i
      let deriv := ak * f.2 * k2i * kd
      let term := fun (a b c e : α) => pre * ((a + b - twoKd * (k2i * c + d2i * e)) * f.1 + deriv * c)
      (⟨acc.gnum.xx + term kxd.xx dxk.xx kxk.xx dxd.xx, acc.gnum.xy + term kxd.xy dxk.xy kxk.xy dxd.xy,
        acc.gnum.xz + term kxd.xz dxk.xz kxk.xz dxd.xz, acc.gnum.yy + term kxd.yy dxk.yy kxk.yy dxd.yy,
        acc.gnum.yz + term kxd.yz dxk.yz kxk.yz dxd.yz,
        acc.gnum.zz + term kxd.zz dxk.zz kxk.zz dxd.zz⟩ : Sym α)
    else acc.gnum
  ⟨num, acc.den + ph2, acc.k2if2 + f.1 * k2iF, gnum⟩

/-- `alias_sums(inp, k_l, D̂, |D̂|², with_gamma)`. -/
def aliasSumsModel (pi : α) (inp : Inputs α) (kl d : V3 α) (abs2d : α) (withGamma : Bool) :
    AliasSums α :=
  let ms := (intRange inp.nalias).flatMap fun m1 => (intRange inp.nalias).flatMap fun m2 =>
    (intRange inp.nalias).map fun m3 => (m1, m2, m3)
  let s := ms.foldl (aliasStep α pi inp kl d abs2d withGamma) ⟨V3.zero, 0, 0, Sym.zero⟩
  if inp.order = 1 then { s with den := 1 } else s

/-- `influence_at(inp, l, with_gamma)`: `(Ĝ, Γ̂, q)`, with `eps` the code's `1e-12`. -/
def influenceAtModel (pi eps : α) (inp : Inputs α) (l : ℤ × ℤ × ℤ) (withGamma : Bool) :
    α × Sym α × α :=
  let kl := latticeVectorModel α inp.recip.r0 inp.recip.r1 inp.recip.r2 l
  let dd := dHatModel α inp kl
  let s := aliasSumsModel α pi inp kl dd.1 dd.2 withGamma
  if dd.2 < eps then (0, Sym.zero, s.k2if2)
  else
    let numerator := dd.1.x * s.num.x + dd.1.y * s.num.y + dd.1.z * s.num.z
    let denominator := dd.2 * s.den * s.den
    let g := numerator / denominator
    (g, if withGamma then
        ⟨s.gnum.xx / denominator, s.gnum.xy / denominator, s.gnum.xz / denominator,
          s.gnum.yy / denominator, s.gnum.yz / denominator, s.gnum.zz / denominator⟩
      else Sym.zero,
      s.k2if2 - numerator * g)

/-- `quality_at(inp, l, ghat)`: the two terms added to `Q`. -/
def qualityAtModel (pi eps : α) (inp : Inputs α) (l : ℤ × ℤ × ℤ) (g : α) : α × α :=
  let kl := latticeVectorModel α inp.recip.r0 inp.recip.r1 inp.recip.r2 l
  let dd := dHatModel α inp kl
  let s := aliasSumsModel α pi inp kl dd.1 dd.2 false
  if inp.nfdord ≠ 0 ∧ dd.2 < eps then (s.k2if2, 0)
  else
    (g * g * dd.2 * s.den * s.den,
      s.k2if2 - 2 * g * (dd.1.x * s.num.x + dd.1.y * s.num.y + dd.1.z * s.num.z))

/-- The full matrix of a symmetric tensor. -/
def symMatrix (m : Sym α) : M3 α := ⟨⟨m.xx, m.xy, m.xz⟩, ⟨m.xy, m.yy, m.yz⟩, ⟨m.xz, m.yz, m.zz⟩⟩

/-- `influence_box_correction(ghat, gammahat, tl0_inv, delta_tl)`. -/
def boxCorrectionModel (g : α) (gam : Sym α) (t0 dl : M3 α) : α :=
  let gm := symMatrix α gam
  let row := fun (r : V3 α) => (⟨0 + r.x * t0.r0.x + r.y * t0.r1.x + r.z * t0.r2.x,
    0 + r.x * t0.r0.y + r.y * t0.r1.y + r.z * t0.r2.y,
    0 + r.x * t0.r0.z + r.y * t0.r1.z + r.z * t0.r2.z⟩ : V3 α)
  let p0 := row gm.r0
  let p1 := row gm.r1
  let p2 := row gm.r2
  let t0' := 0 + p0.x * dl.r0.x + p0.y * dl.r1.x + p0.z * dl.r2.x
  let t1' := 0 + p1.x * dl.r0.y + p1.y * dl.r1.y + p1.z * dl.r2.y
  let t2' := 0 + p2.x * dl.r0.z + p2.y * dl.r1.z + p2.z * dl.r2.z
  g - (t0' + t1' + t2')

/-- `force_error_estimate(q, sum_q2, volume, n_atoms, fpe)`. -/
def forceErrorModel (q sumQ2 volume : α) (nAtoms : ℕ) (fpe : α) : α :=
  fpe * sumQ2 * Sqrt.sqrt (q / (volume * (nAtoms : α)))

/-- `Γ̂ + Ĝ 1`. -/
def withDiag (t : Sym α) (g : α) : Sym α := ⟨t.xx + g, t.xy, t.xz, t.yy + g, t.yz, t.zz + g⟩

/-- `kspace_energy(ghat, gammahat, rho_re, rho_im, cell_volume, eps0_inv)` over the grid points
`(Ĝ, Γ̂, Re ρ̂, Im ρ̂)` in storage order: `(E, Ξ)`. -/
def kspaceEnergyModel (pts : List (α × Sym α × α × α)) (withGamma : Bool) (cellVolume eps0Inv : α) :
    α × Sym α :=
  let n : α := (pts.length : α)
  let r := pts.foldl (fun (acc : α × Sym α) pt =>
    let abs2 := (pt.2.2.1 * pt.2.2.1 + pt.2.2.2 * pt.2.2.2) / n
    let e := acc.1 + pt.1 * abs2
    if withGamma then
      let t := withDiag α pt.2.1 pt.1
      (e, ⟨acc.2.xx + t.xx * abs2, acc.2.xy + t.xy * abs2, acc.2.xz + t.xz * abs2,
        acc.2.yy + t.yy * abs2, acc.2.yz + t.yz * abs2, acc.2.zz + t.zz * abs2⟩)
    else (e, acc.2)) (0, Sym.zero)
  let s := -(1 / 4) * eps0Inv * cellVolume
  (1 / 2 * r.1 * cellVolume * eps0Inv, r.2.mulRight s)

/-- `potential_hat` at one grid point: `((eps0_inv Ĝ) Re ρ̂ / √N, (eps0_inv Ĝ) Im ρ̂ / √N)`. -/
def potentialHatModel (n : ℕ) (g re im eps0Inv : α) : α × α :=
  let sqrtN := Sqrt.sqrt (n : α)
  let f := eps0Inv * g
  (f * re / sqrtN, f * im / sqrtN)

/-- `field_ik` at one grid point `(x, y, z)` with `φ̂ = (re, im)`: `[(Re Ê_α, Im Ê_α)]`. -/
def fieldIkModel (m : Mesh α) (recip : M3 α) (x y z : ℕ) (re im : α) : V3 α × V3 α :=
  let sqrtN := Sqrt.sqrt (m.len : α)
  let l := (waveIndex x m.grid.1, waveIndex y m.grid.2.1, waveIndex z m.grid.2.2)
  let k := latticeVectorModel α recip.r0 recip.r1 recip.r2 l
  let mre := im / sqrtN
  let mim := -re / sqrtN
  (⟨k.x * mre, k.y * mre, k.z * mre⟩, ⟨k.x * mim, k.y * mim, k.z * mim⟩)

/-- `field_fd` at one grid point `(x, y, z)` of the potential grid `phi`. -/
def fieldFdModel (m : Mesh α) (q : ℕ) (phi : Array α) (x y z : ℤ) : V3 α :=
  let sqrtN := Sqrt.sqrt (m.len : α)
  let c := fdCoefficientsModel α q
  let at_ := fun (a b e : ℤ) => phi.getD (gridIndex m.grid a b e) 0
  let g := (List.range q).foldl (fun (g : V3 α) (i : ℕ) =>
    let j : ℤ := (i : ℤ) + 1
    let twoJ := 2 * (j : α)
    let dx := (at_ (x - j) y z - at_ (x + j) y z) / twoJ
    let dy := (at_ x (y - j) z - at_ x (y + j) z) / twoJ
    let dz := (at_ x y (z - j) - at_ x y (z + j)) / twoJ
    let f := c.getD i 0 / sqrtN
    ⟨g.x + f * dx, g.y + f * dy, g.z + f * dz⟩) V3.zero
  m.hInv.mulVec g

/-- `interpolate_forces` for one atom `(r, q)` and the field grids `ex, ey, ez`. -/
def interpolateModel (m : Mesh α) (p : ℕ) (ex ey ez : Array α) (a : V3 α × α) : V3 α :=
  let u := m.hInv.tmulVec a.1
  let nx := nearestIntModel α p u.x
  let ny := nearestIntModel α p u.y
  let nz := nearestIntModel α p u.z
  let e := (stencil p).foldl (fun e dx =>
    let px := nx + dx
    let wx := assignmentModel α p ((px : α) - u.x)
    (stencil p).foldl (fun e dy =>
      let py := ny + dy
      let wy := assignmentModel α p ((py : α) - u.y)
      (stencil p).foldl (fun (e : V3 α) dz =>
        let pz := nz + dz
        let wz := assignmentModel α p ((pz : α) - u.z)
        let w := wx * wy * wz
        let i := gridIndex m.grid px py pz
        ⟨e.x + w * ex.getD i 0, e.y + w * ey.getD i 0, e.z + w * ez.getD i 0⟩) e) e) V3.zero
  ⟨a.2 * e.x, a.2 * e.y, a.2 * e.z⟩

/-- `a2_tilde_p3m(ghat, gammahat, rq_re, sum_q2)` over the grid points `(Ĝ, Γ̂, Re R̂)`. -/
def a2TildeModel (pi : α) (pts : List (α × Sym α × α)) (withGamma : Bool) (sumQ2 : α) :
    α × Sym α :=
  let n : α := (pts.length : α)
  let r := pts.foldl (fun (acc : α × Sym α) pt =>
    let s := acc.1 + pt.1 * pt.2.2
    if withGamma then
      let t := withDiag α pt.2.1 pt.1
      (s, ⟨acc.2.xx + t.xx * pt.2.2, acc.2.xy + t.xy * pt.2.2, acc.2.xz + t.xz * pt.2.2,
        acc.2.yy + t.yy * pt.2.2, acc.2.yz + t.yz * pt.2.2, acc.2.zz + t.zz * pt.2.2⟩)
    else (s, acc.2)) (0, Sym.zero)
  let fourPi := 4 * pi
  let s := fourPi / n / sumQ2
  (fourPi * r.1 / n / sumQ2, Sym.mulLeft s r.2)

end Model

/-! ## The `Float` instantiation (what `lake exe spec-vectors` runs) -/

/-- gromosXX's `math::epsilon = 0.000000000001` by its bits. -/
def epsFloat : Float := Float.ofBits 0x3d719799812dea11

def assignmentFloat (p : ℕ) (x : Float) : Float := assignmentModel Float p x
def nearestIntFloat (p : ℕ) (u : Float) : ℤ := nearestIntModel Float p u
def weightsFloat (p : ℕ) (u : Float) : List Float := weightsModel Float p u
def selftermFpFloat (p n : ℕ) (s : Float) : Float := selftermFpModel Float p n s
def fdCoefficientsFloat (q : ℕ) : List Float := fdCoefficientsModel Float q
def sincPowFloat (x : Float) (p : ℕ) : Float := sincPowModel Float x p
def meshFloat (edges : M3 Float) (g : ℕ × ℕ × ℕ) : Mesh Float := meshModel Float edges g
def assignChargesFloat (m : Mesh Float) (p : ℕ) (atoms : List (V3 Float × Float)) : Array Float :=
  assignChargesModel Float m p atoms
def assignSquaredFloat (m : Mesh Float) (p : ℕ) (atoms : List (V3 Float × Float)) : Array Float :=
  assignSquaredModel Float m p atoms
def influenceAtFloat (inp : Inputs Float) (l : ℤ × ℤ × ℤ) (withGamma : Bool) :
    Float × Sym Float × Float :=
  influenceAtModel Float piFloat epsFloat inp l withGamma
def qualityAtFloat (inp : Inputs Float) (l : ℤ × ℤ × ℤ) (g : Float) : Float × Float :=
  qualityAtModel Float piFloat epsFloat inp l g
def boxCorrectionFloat (g : Float) (gam : Sym Float) (t0 dl : M3 Float) : Float :=
  boxCorrectionModel Float g gam t0 dl
def forceErrorFloat (q sumQ2 volume : Float) (nAtoms : ℕ) (fpe : Float) : Float :=
  forceErrorModel Float q sumQ2 volume nAtoms fpe
def kspaceEnergyFloat (pts : List (Float × Sym Float × Float × Float)) (withGamma : Bool)
    (cellVolume eps0Inv : Float) : Float × Sym Float :=
  kspaceEnergyModel Float pts withGamma cellVolume eps0Inv
def potentialHatFloat (n : ℕ) (g re im eps0Inv : Float) : Float × Float :=
  potentialHatModel Float n g re im eps0Inv
def fieldIkFloat (m : Mesh Float) (recip : M3 Float) (x y z : ℕ) (re im : Float) :
    V3 Float × V3 Float := fieldIkModel Float m recip x y z re im
def fieldFdFloat (m : Mesh Float) (q : ℕ) (phi : Array Float) (x y z : ℤ) : V3 Float :=
  fieldFdModel Float m q phi x y z
def interpolateFloat (m : Mesh Float) (p : ℕ) (ex ey ez : Array Float) (a : V3 Float × Float) :
    V3 Float := interpolateModel Float m p ex ey ez a
def a2TildeFloat (pts : List (Float × Sym Float × Float)) (withGamma : Bool) (sumQ2 : Float) :
    Float × Sym Float := a2TildeModel Float piFloat pts withGamma sumQ2

/-- The influence function on the whole grid in gromosXX's order (`influence_function`):
`(Ĝ, Γ̂, Q)`. -/
def influenceFunctionFloat (inp : Inputs Float) (volume : Float) :
    Array Float × Array (Sym Float) × Float :=
  let g := inp.mesh.grid
  let ny2 : ℤ := (g.2.1 / 2 : ℕ)
  let nz2 : ℤ := (g.2.2 / 2 : ℕ)
  let ls := fun (lo hi : ℤ) => (List.range (hi - lo + 1).toNat).map fun i => lo + i
  let r := (List.range g.1).foldl (fun (acc : Array Float × Array (Sym Float) × Float) x =>
    let l0 := waveIndex x g.1
    let plane := (ls (1 - ny2) ny2).foldl (fun (acc : Array Float × Array (Sym Float) × Float) l1 =>
      (ls (1 - nz2) nz2).foldl (fun (acc : Array Float × Array (Sym Float) × Float) l2 =>
        if l0 = 0 ∧ l1 = 0 ∧ l2 = 0 then acc
        else
          let pt := influenceAtFloat inp (l0, l1, l2) true
          let i := gridIndex g l0 l1 l2
          (acc.1.set! i pt.1, acc.2.1.set! i pt.2.1, acc.2.2 + pt.2.2)) acc) (acc.1, acc.2.1, 0)
    (plane.1, plane.2.1, acc.2.2 + plane.2.2))
    (Array.replicate inp.mesh.len 0, Array.replicate inp.mesh.len Sym.zero, 0)
  (r.1, r.2.1, r.2.2 * (16 * piFloat * piFloat / volume))

/-- The discrete Fourier transform of a grid by its definition, `O(N²)`, with the sign `sgn` of
the exponent (−1 forward, +1 backward, both unnormalised): the twiddle of each axis from a table
of `cos`, `sin` of `2π m / N_axis`, `m = j k mod N_axis`. -/
def dftFloat (g : ℕ × ℕ × ℕ) (re im : Array Float) (sgn : Float) : Array Float × Array Float :=
  let tw := fun (n : ℕ) => (Array.range n).map fun m =>
    let a := 2 * piFloat * m.toFloat / n.toFloat
    (Float.cos a, sgn * Float.sin a)
  let (tx, ty, tz) := (tw g.1, tw g.2.1, tw g.2.2)
  let cmul := fun (a b : Float × Float) => (a.1 * b.1 - a.2 * b.2, a.1 * b.2 + a.2 * b.1)
  let out := (List.range g.1).flatMap fun kx => (List.range g.2.1).flatMap fun ky =>
    (List.range g.2.2).map fun kz =>
      (List.range g.1).foldl (fun acc jx =>
        let wx := tx.getD (kx * jx % g.1) (1, 0)
        (List.range g.2.1).foldl (fun acc jy =>
          let wxy := cmul wx (ty.getD (ky * jy % g.2.1) (1, 0))
          (List.range g.2.2).foldl (fun acc jz =>
            let w := cmul wxy (tz.getD (kz * jz % g.2.2) (1, 0))
            let i := jz + g.2.2 * (jy + g.2.1 * jx)
            let v := cmul w (re.getD i 0, im.getD i 0)
            (acc.1 + v.1, acc.2 + v.2)) acc) acc) (0, 0)
  (out.toArray.map (·.1), out.toArray.map (·.2))

/-! ## The `ℝ` instantiation is the theory -/

section Real

open NewenTheory.Interactions.P3M

/-- Over ℝ an integer-valued number truncates to itself (only used through `nearestIntModel`). -/
noncomputable instance : ToInt ℝ := ⟨fun x => if x < 0 then -⌊-x⌋ else ⌊x⌋⟩

/-- `assignment_1d` over ℝ is the theory's `w p`, for every order and every `x`. -/
theorem assignmentModel_real_eq (p : ℕ) (x : ℝ) : assignmentModel ℝ p x = w p x := by
  match p with
  | 0 => rfl
  | 1 => simp [assignmentModel, w, Fabs.fabs]
  | 2 => simp [assignmentModel, w, Fabs.fabs]
  | 3 => simp only [assignmentModel, w, Fabs.fabs]; split_ifs <;> ring
  | 4 => simp only [assignmentModel, w, Fabs.fabs]; split_ifs <;> (try rw [← sq_abs x]) <;> ring
  | 5 => simp only [assignmentModel, w, Fabs.fabs]; split_ifs <;> ring
  | p + 6 => rfl

/-- `sinc_power` over ℝ is the theory's `sincPow`. -/
theorem sincPowModel_real_eq (x : ℝ) (p : ℕ) : sincPowModel ℝ x p = sincPow x p := by
  unfold sincPowModel sincPow
  by_cases h : x = 0
  · simp [h]
  · have hx : ¬(¬x < 0 ∧ ¬0 < x) := fun ⟨a, b⟩ => h (le_antisymm (not_lt.mp b) (not_lt.mp a))
    simp only [hx, h, ite_false]
    simp [Powf.powf, Sin.sin, Real.rpow_natCast]

/-- `fd_coefficients` over ℝ is the theory's table `fdCoeff`. -/
theorem fdCoefficientsModel_real_eq {q j : ℕ} (hq : q ≤ 5) (hj1 : 1 ≤ j) (hj5 : j ≤ 5) :
    (fdCoefficientsModel ℝ q).getD (j - 1) 0 = (fdCoeff q j : ℝ) := by
  interval_cases q <;> interval_cases j <;> norm_num [fdCoefficientsModel, fdCoeff]

/-- The first component of a fold that adds `f b` to it and leaves the second alone. -/
theorem foldl_fst_add_const {β σ : Type} (f : β → ℝ) (l : List β) (a : ℝ) (s : σ) :
    l.foldl (fun (acc : ℝ × σ) b => (acc.1 + f b, acc.2)) (a, s) = (a + (l.map f).sum, s) := by
  induction l generalizing a with
  | nil => simp
  | cons b l ih => simp [ih]; ring

/-- `kspace_energy` over ℝ (without the virial): `½ (Σ_l Ĝ_l |ρ̂_l|²/N) V_G ε₀⁻¹`. -/
theorem kspaceEnergyModel_real_eq (pts : List (ℝ × Sym ℝ × ℝ × ℝ)) (vc e : ℝ) :
    (kspaceEnergyModel ℝ pts false vc e).1 =
      1 / 2 * (pts.map fun pt => pt.1 * ((pt.2.2.1 * pt.2.2.1 + pt.2.2.2 * pt.2.2.2) /
        (pts.length : ℝ))).sum * vc * e := by
  simp only [kspaceEnergyModel, Bool.false_eq_true, ite_false]
  rw [foldl_fst_add_const]
  ring

open ZMod in
/-- **The model's k-space energy is the theory's** `p3mEnergy` when it is given the transform of
a real grid (one axis of `n + 1` points; `ZMod (n + 1)` is `Fin (n + 1)`): the fold over the
points in order is `½ V_G ε₀⁻¹ Σ_l Ĝ_l |ρ̂_l|²/N`. -/
theorem kspaceEnergyModel_real_p3mEnergy {n : ℕ} (ghat ρ : ZMod (n + 1) → ℝ) (vc e : ℝ) :
    (kspaceEnergyModel ℝ (List.ofFn fun l : Fin (n + 1) =>
        (ghat l, (Sym.zero : Sym ℝ), (𝓕 (fun j => (ρ j : ℂ)) l).re,
          (𝓕 (fun j => (ρ j : ℂ)) l).im)) false vc e).1 = p3mEnergy vc e ghat ρ := by
  rw [kspaceEnergyModel_real_eq, List.length_ofFn, List.map_ofFn, List.sum_ofFn, p3mEnergy]
  have h : ∀ z : ℂ, z.re * z.re + z.im * z.im = ‖z‖ ^ 2 := fun z => by
    rw [Complex.sq_norm, Complex.normSq_apply]
  simp only [Function.comp_def, h]
  rw [show ((n + 1 : ℕ) : ℝ) = ((n + 1 : ℕ) : ℝ) from rfl]
  have hs : (∑ l : Fin (n + 1), ghat l * (‖𝓕 (fun j => (ρ j : ℂ)) l‖ ^ 2 / ((n + 1 : ℕ) : ℝ))) =
      ∑ l : ZMod (n + 1), ghat l * ‖𝓕 (fun j => (ρ j : ℂ)) l‖ ^ 2 / ((n + 1 : ℕ) : ℝ) := by
    refine Finset.sum_congr rfl fun l _ => ?_
    ring
  rw [hs]
  push_cast
  ring

/-- The transform of the three-dimensional assignment function at `v` for the grid matrix `h`
(rows `h_α`): `Π_α P̂_p(h_α · v)`. -/
noncomputable def pHatOf (h : M3 ℝ) (p : ℕ) (v : NewenTheory.Interactions.LatticeSum.Vec3) : ℝ :=
  sincPow (h.r0.toVec ⬝ᵥ v) p * sincPow (h.r1.toVec ⬝ᵥ v) p * sincPow (h.r2.toVec ⬝ᵥ v) p

theorem intRange_zero : intRange 0 = [0] := by simp [intRange]

/-- **One alias, `ik` differentiation** (NALIAS = 0, NFDORD = 0, NASORD ≥ 2): `influence_at`
over ℝ is the theory's optimal influence function (Vol. 2 `optim_influ`) over the single alias
`k_l` itself, for a wave vector with `|k_l|² ≥ eps > 0` and a width `a ≠ 0`; so
`Ĝ = γ̂(a|k|)/(|k|² P̂(k)²)` where `P̂(k) ≠ 0` (`influence_single_alias`). -/
theorem influenceAtModel_real_single {eps : ℝ} (heps : 0 < eps) (inp : Inputs ℝ)
    (hna : inp.nalias = 0) (hfd : inp.nfdord = 0) (hord : inp.order ≠ 1) (hw : inp.width ≠ 0)
    (l : ℤ × ℤ × ℤ) (k : V3 ℝ) (hk : k = latticeVectorModel ℝ inp.recip.r0 inp.recip.r1 inp.recip.r2 l)
    (hkk : eps ≤ k.x * k.x + k.y * k.y + k.z * k.z) :
    (influenceAtModel ℝ Real.pi eps inp l false).1 =
      influence inp.shape inp.width (pHatOf inp.mesh.h inp.order) k.toVec {k.toVec} := by
  have hkk0 : 0 < k.x * k.x + k.y * k.y + k.z * k.z := lt_of_lt_of_le heps hkk
  have hn2 : norm2 k.toVec = k.x * k.x + k.y * k.y + k.z * k.z := by
    simp [norm2, V3.toVec, dotProduct, Fin.sum_univ_three]
  have hak : inp.width * Real.sqrt (k.x * k.x + k.y * k.y + k.z * k.z) ≠ 0 :=
    mul_ne_zero hw (Real.sqrt_ne_zero'.mpr hkk0)
  have hf := shapeFourierModel_real_eq inp.shape hak
  simp only [influenceAtModel, ← hk, dHatModel, hfd, ite_true, aliasSumsModel, hna,
    intRange_zero, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil,
    List.append_nil, List.foldl_cons, List.foldl_nil, hord, ite_false, not_lt.mpr hkk]
  simp only [aliasStep, M3.mulVec, V3.dot, Int.cast_zero, mul_zero, add_zero, Sqrt.sqrt, hf,
    sincPowModel_real_eq, Bool.false_eq_true, ite_false, V3.zero]
  simp only [influence, Finset.sum_singleton, pHatOf, hn2, ← V3.dot_eq, V3.dot]
  simp only [V3.toVec, dotProduct, Fin.sum_univ_three, Pi.smul_apply, smul_eq_mul,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons]
  congr 1 <;> ring

/-- `force_error_estimate` over ℝ: `fpe · S̃² · √(Q/(V N))`. -/
theorem forceErrorModel_real_eq (q sumQ2 V : ℝ) (n : ℕ) (fpe : ℝ) :
    forceErrorModel ℝ q sumQ2 V n fpe = fpe * sumQ2 * Real.sqrt (q / (V * n)) := rfl

end Real

end NewenSpec.P3M
