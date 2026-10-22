import Aeneas
import NewenSpec.Schedule
import NewenSpec.Thermo.Dof
import NewenSpec.Thermo.Bath
import NewenSpec.Topology.Pairs
import NewenSpec.Topology.Groups
import NewenSpec.Extracted.NewenLogic

/-!
# The extracted Rust computes the Lean model

`NewenSpec/Extracted/NewenLogic.lean` is what Aeneas made of `crates/newen-logic`. For
`crates/newen-logic/src/schedule.rs` it is a pure Lean function over `Usize` (the machine integer,
`0 ≤ n < 2^64`) returning `Result Bool`, where `fail` would stand for a Rust panic (here: `%` by
zero) and `ok b` for a normal return of `b`.

The first theorems below say: for every pair of `usize` inputs the Rust function returns
normally (it never panics) and returns exactly `missesModel`, hence (with
`Schedule.missesModel_iff`) exactly `Misses` of `NewenTheory/Dynamics/Schedule.lean`.
`thermo/dof.rs` (`total_dof`, over `U64` and `I64`), `thermo/bath.rs` (five functions over
slices) and `topology/` (pair kinds and group lookups, five functions) have their own sections
at the end.
-/

open Aeneas Aeneas.Std Aeneas.Std.WP Result

namespace NewenSpec.Refinement

open newen_logic.schedule (output_misses_longrange)

/-- A computation whose every outcome is `v` returns `v`. -/
theorem eq_ok_of_spec {α : Type} {m : Result α} {v : α} (h : m ⦃ b => b = v ⦄) : m = ok v := by
  cases m with
  | ret r => rw [spec_ok] at h; rw [h]
  | vis i k => rw [spec_vis] at h; exact h.elim
  | div => rw [spec_div] at h; exact h.elim

/-- Aeneas' style of statement (a Hoare triple): the Rust function returns `missesModel`. -/
theorem output_misses_longrange_spec (nsnb ntwg : Usize) :
    output_misses_longrange nsnb ntwg ⦃ b => b = Schedule.missesModel nsnb.val ntwg.val ⦄ := by
  unfold output_misses_longrange Schedule.missesModel
  -- `step*` walks the function: it splits the `if`s and, before each `%`, proves the divisor is
  -- nonzero (so `%` cannot panic) and records the remainder's value.
  step*
  all_goals
    have e1 : ¬ nsnb.val ≤ 1 := by scalar_tac
    have e2 : ntwg.val ≠ 0 := by scalar_tac
    simp only [e1, e2, ite_false]
  · rw [Bool.eq_iff_iff]
    simp_all [bne_iff_ne, UScalar.eq_equiv]
  · simp_all [UScalar.eq_equiv]

/-- **Rust = model**: for all `usize` inputs, `output_misses_longrange` returns `missesModel`. -/
theorem output_misses_longrange_eq (nsnb ntwg : Usize) :
    output_misses_longrange nsnb ntwg = ok (Schedule.missesModel nsnb.val ntwg.val) :=
  eq_ok_of_spec (output_misses_longrange_spec nsnb ntwg)

/-- **Rust = meaning**: the Rust function returns `true` exactly when some output step reuses a
pairlist step that was not an output step. -/
theorem output_misses_longrange_iff (nsnb ntwg : Usize) :
    output_misses_longrange nsnb ntwg = ok true ↔
      NewenTheory.Dynamics.Schedule.Misses nsnb.val ntwg.val := by
  rw [output_misses_longrange_eq, Result.ok.injEq, Schedule.missesModel_iff]

/-! ## `thermo::dof::total_dof`

The extracted function converts the three `u64` counts to `i64` (`UScalar.hcast`, a wrapping
cast exactly like Rust's `as`) and then multiplies and subtracts in `i64`, where Aeneas turns an
overflow into `fail` (a Rust panic with overflow checks on). The theorems below assume
`3 · n_atoms + solvent_constraints + solute_constraints + |ndfmin| ≤ i64::MAX`. That is enough
for every cast to keep its value and for every intermediate result to fit, so the function
returns normally with the model's value. Outside that bound the proof says nothing: depending on
the input, an operation overflows (`fail`), or a count `≥ 2^63` wraps in the cast. No caller can
reach that region: a count of atoms or constraints is at most the length of a Rust collection,
far below `2^63 / 3`. -/

open newen_logic.thermo.dof (total_dof)

/-- No `i64` overflow anywhere in `total_dof`: the bound the refinement assumes. -/
def TotalDofInRange (nAtoms solventConstraints soluteConstraints : Nat) (ndfmin : Int) : Prop :=
  3 * (nAtoms : Int) + solventConstraints + soluteConstraints + ndfmin.natAbs ≤ I64.max

/-- Rust's `x as i64` keeps the value of a `u64` that fits in `i64`. -/
theorem hcast_I64_val (x : U64) (h : (x.val : Int) ≤ I64.max) :
    (UScalar.hcast .I64 x).val = (x.val : Int) := by
  have hs := UScalar.hcast_inBounds_spec .I64 x (by simpa using h)
  simpa [lift, spec_ok] using hs

/-- Aeneas' style of statement (a Hoare triple): within the bound, the Rust function returns
`totalDofModel`. -/
theorem total_dof_spec (nAtoms solventConstraints soluteConstraints : U64) (ndfmin : I64)
    (h : TotalDofInRange nAtoms.val solventConstraints.val soluteConstraints.val ndfmin.val) :
    total_dof nAtoms solventConstraints soluteConstraints ndfmin ⦃ r =>
      r.val = Thermo.Dof.totalDofModel nAtoms.val solventConstraints.val soluteConstraints.val
        ndfmin.val ⦄ := by
  unfold TotalDofInRange at h
  unfold total_dof Thermo.Dof.totalDofModel
  -- Each cast keeps its value; `step*` then proves each `*` and `-` stays within `i64`.
  have hn := hcast_I64_val nAtoms (by omega)
  have hs := hcast_I64_val solventConstraints (by omega)
  have hc := hcast_I64_val soluteConstraints (by omega)
  have hf := ndfmin.hBounds
  step*
  scalar_tac

/-- A computation whose every outcome satisfies `P` returns a value satisfying `P`. -/
theorem exists_ok_of_spec {α : Type} {m : Result α} {P : α → Prop} (h : m ⦃ b => P b ⦄) :
    ∃ b, m = ok b ∧ P b := by
  cases m with
  | ret r => rw [spec_ok] at h; exact ⟨r, rfl, h⟩
  | vis i k => rw [spec_vis] at h; exact h.elim
  | div => rw [spec_div] at h; exact h.elim

/-- **Rust = model**: within the bound, `total_dof` returns normally, with `totalDofModel`. -/
theorem total_dof_eq (nAtoms solventConstraints soluteConstraints : U64) (ndfmin : I64)
    (h : TotalDofInRange nAtoms.val solventConstraints.val soluteConstraints.val ndfmin.val) :
    ∃ r, total_dof nAtoms solventConstraints soluteConstraints ndfmin = ok r ∧
      r.val = Thermo.Dof.totalDofModel nAtoms.val solventConstraints.val soluteConstraints.val
        ndfmin.val :=
  exists_ok_of_spec (total_dof_spec nAtoms solventConstraints soluteConstraints ndfmin h)

open NewenTheory.Dynamics.DegreesOfFreedom in
/-- **Rust = meaning**: called as `newen-run` calls it, for solute molecules `solute` and `k`
copies of the solvent molecule `solvent`, `total_dof` returns the degrees of freedom of that
system. -/
theorem total_dof_eq_dof (solute : List Molecule) (solvent : Molecule) (k : Nat)
    (nAtoms solventConstraints soluteConstraints : U64) (ndfmin : I64)
    (hn : nAtoms.val = atoms solute + k * solvent.atoms)
    (hs : solventConstraints.val = k * solvent.constraints)
    (hc : soluteConstraints.val = constraints solute)
    (h : TotalDofInRange nAtoms.val solventConstraints.val soluteConstraints.val ndfmin.val) :
    ∃ r, total_dof nAtoms solventConstraints soluteConstraints ndfmin = ok r ∧
      r.val = dof (solute ++ List.replicate k solvent) ndfmin.val := by
  obtain ⟨r, hr, hv⟩ := total_dof_eq nAtoms solventConstraints soluteConstraints ndfmin h
  refine ⟨r, hr, ?_⟩
  rw [hv, hn, hs, hc, Thermo.Dof.totalDofModel_eq_dof]

/-- **Overflow is not silent**: an atom count that fits in `i64` but whose `3 · n_atoms` does not
makes `total_dof` fail (panic) instead of returning a wrapped number. -/
theorem total_dof_overflow (nAtoms solventConstraints soluteConstraints : U64) (ndfmin : I64)
    (h1 : (nAtoms.val : Int) ≤ I64.max) (h2 : I64.max < 3 * (nAtoms.val : Int)) (r : I64) :
    total_dof nAtoms solventConstraints soluteConstraints ndfmin ≠ ok r := by
  unfold total_dof
  have hn := hcast_I64_val nAtoms h1
  have hm := IScalar.mul_equiv 3#i64 (UScalar.hcast .I64 nAtoms)
  simp only [lift, bind_ok]
  generalize hx : UScalar.hcast IScalarTy.I64 nAtoms = x at *
  generalize hmul : (3#i64 * x) = m at *
  simp only [HMul.hMul] at hmul
  cases m with
  | ret z =>
    rw [hmul, Result.match.ok] at hm
    scalar_tac
  | vis i k => simp only [bind_vis]; exact vis_not_ok
  | div =>
    rw [hmul, Result.match.div] at hm
    exact hm.elim

/-! ## `thermo::bath`: degrees of freedom per bath

The five functions loop over slices; Aeneas turns each `while` loop into a `loop` over the loop
state, and each proof below runs the loop with `loop.spec_decr_nat`: an invariant saying that
what is accumulated so far plus what the model gives for the rest of the slice (`drop i`) is the
model's value for the whole slice, and a measure (`len − i`) that shrinks.

The integer arithmetic is `i64` with overflow checks, so each proof assumes that no operation
overflows (`BathHalfDofInRange`, `SpreadInRange`). For the ranges `newen-run` passes (they cover
the `N` atoms once, `Covers`) the range loop adds `Σ_r n_r = N` atoms and `Σ_r g_r ≤ G` groups,
so `BathHalfDofInRange` holds whenever `6·(N + 2·G) + 2·C + 12 ≤ i64::MAX` (`N`, `G`, `C` atoms,
groups and constraints, each at most a Rust collection's length) and no `last_atom` is
`usize::MAX`, i.e. for any system below `10^17` atoms. Every atom, group and constraint end is
counted in one bath, so the half counts that reach `ndfmin_spread` have
`Σ_b |h_b| ≤ 6·N + 12·G + 2·C + 12`; `SpreadInRange` needs
`(Σ_b |h_b| + 2·|ndfmin| + 2)² ≤ i64::MAX`, which holds for systems below about `10^8` atoms
(`√i64::MAX ≈ 3·10^9`). -/

section Bath

open newen_logic.thermo.bath
open NewenSpec.Thermo.Bath

/-- The extracted `DofRange` as the model's `Range`. -/
def toRange (d : DofRange) : Range := ⟨d.last_atom.val, d.com_bath.val, d.ir_bath.val⟩

/-- A slice of DOFSET ranges as the model sees it. -/
def rangesOf (s : Slice DofRange) : List Range := s.val.map toRange

/-- A slice of atom indices as naturals. -/
def natsOf (s : Slice Usize) : List Nat := s.val.map (·.val)

/-- A slice of constraints as pairs of atoms. -/
def consOf (s : Slice Constraint) : List (Nat × Nat) := s.val.map fun c => (c.i.val, c.j.val)

/-- A slice of `i64` as integers. -/
def intsOf (s : Slice I64) : List Int := s.val.map (·.val)

/-- `baths_of`'s result as the model's. -/
def pairVal (o : Option (Usize × Usize)) : Option (Nat × Nat) := o.map fun p => (p.1.val, p.2.val)

/-- `HalfDof` as the model's pair `(com, ir)`. -/
def halfVal (h : HalfDof) : Int × Int := (h.com.val, h.ir.val)

theorem bathsOfModel_drop (rs : List Range) (i : Nat) (h : i < rs.length) (a : Nat) :
    bathsOfModel (rs.drop i) a =
      if a ≤ rs[i].lastAtom then some (rs[i].comBath, rs[i].irBath)
      else bathsOfModel (rs.drop (i + 1)) a := by
  rw [List.drop_eq_getElem_cons h]; rfl

@[step]
theorem baths_of_spec (ranges : Slice DofRange) (atom : Usize) :
    baths_of ranges atom ⦃ o => pairVal o = bathsOfModel (rangesOf ranges) atom.val ⦄ := by
  unfold baths_of baths_of_loop
  apply loop.spec_decr_nat (fun r => ranges.length - r.val)
    (fun r => r.val ≤ ranges.length ∧
      bathsOfModel ((rangesOf ranges).drop r.val) atom.val = bathsOfModel (rangesOf ranges) atom.val)
  · intro r ⟨hr, hinv⟩
    unfold baths_of_loop.body
    step*
    · rename_i hlt hle
      have hl : r.val < (rangesOf ranges).length := by simp [rangesOf]; scalar_tac
      rw [← hinv, bathsOfModel_drop _ _ hl]
      have hle' : atom.val ≤ ((rangesOf ranges)[r.val]'hl).lastAtom := by
        simp [rangesOf, toRange, ← dr_post]; scalar_tac
      simp only [hle', ite_true]
      simp [pairVal, rangesOf, toRange, dr_post]
    · rename_i hlt _ hle
      have hl : r.val < (rangesOf ranges).length := by simp [rangesOf]; scalar_tac
      refine ⟨by scalar_tac, ?_, by scalar_tac⟩
      rw [← hinv, bathsOfModel_drop _ _ hl]
      have hle' : ¬ atom.val ≤ ((rangesOf ranges)[r.val]'hl).lastAtom := by
        simp [rangesOf, toRange, ← dr_post]; scalar_tac
      simp only [hle', ite_false, r1_post]
    · rename_i hge
      have : r.val = ranges.length := by scalar_tac
      rw [← hinv, this]
      simp [pairVal, rangesOf, bathsOfModel]
  · simp

theorem groupsStartingInModel_drop (gs : List Nat) (i : Nat) (h : i < gs.length) (f l : Nat) :
    groupsStartingInModel (gs.drop i) f l =
      (if f ≤ gs[i] ∧ gs[i] ≤ l then 1 else 0) + groupsStartingInModel (gs.drop (i + 1)) f l := by
  rw [List.drop_eq_getElem_cons h]; rfl

theorem groupsStartingInModel_le (gs : List Nat) (f l : Nat) :
    groupsStartingInModel gs f l ≤ gs.length := by
  induction gs with
  | nil => simp [groupsStartingInModel]
  | cons g gs ih => simp only [groupsStartingInModel, List.length_cons]; split <;> omega

@[step]
theorem groups_starting_in_spec (gf : Slice Usize) (first last : Usize) :
    groups_starting_in gf first last ⦃ n =>
      n.val = groupsStartingInModel (natsOf gf) first.val last.val ⦄ := by
  unfold groups_starting_in groups_starting_in_loop
  apply loop.spec_decr_nat (fun p => gf.length - p.2.val)
    (fun p => p.2.val ≤ gf.length ∧ p.1.val ≤ p.2.val ∧
      p.1.val + groupsStartingInModel ((natsOf gf).drop p.2.val) first.val last.val
        = groupsStartingInModel (natsOf gf) first.val last.val)
  · intro ⟨count, g⟩ ⟨hg, hc, hinv⟩
    unfold groups_starting_in_loop.body
    simp only at hg hc hinv ⊢
    have hlen := gf.length_ineq
    have hmax : Usize.max ≤ U64.max := by
      rcases Usize.bounds_eq with h | h <;> rw [h]; scalar_tac
    step*
    · rename_i hlt
      have hl : g.val < (natsOf gf).length := by simp [natsOf]; scalar_tac
      have hd := groupsStartingInModel_drop (natsOf gf) g.val hl first.val last.val
      have hgi : (natsOf gf)[g.val]'hl = i1.val := by simp [natsOf, i1_post]
      rw [hgi] at hd
      by_cases h1 : first ≤ i1 <;> by_cases h2 : i1 ≤ last <;> simp only [h1, h2, ite_true, ite_false]
      all_goals
        step*
        have hu1 : first.val ≤ i1.val ↔ first ≤ i1 := Iff.rfl
        have hu2 : i1.val ≤ last.val ↔ i1 ≤ last := Iff.rfl
        simp only [← hu1, ← hu2] at h1 h2
        simp only [h1, h2, and_true, and_false, ite_true, ite_false] at hd
        rw [g1_post]
        refine ⟨by scalar_tac, by scalar_tac, by scalar_tac, by scalar_tac⟩
    · rename_i hge
      have : g.val = gf.length := by scalar_tac
      rw [← hinv, this]
      simp [natsOf, groupsStartingInModel]
  · simp

theorem hcast_I64_val' {ty : UScalarTy} (x : UScalar ty) (h : (x.val : Int) ≤ I64.max) :
    (UScalar.hcast .I64 x).val = (x.val : Int) := by
  have hs := UScalar.hcast_inBounds_spec .I64 x (by simpa using h)
  simpa [lift, spec_ok] using hs

theorem i64_min_eq : I64.min = -I64.max - 1 := by scalar_tac

theorem eq_ok_of_add (x y : I64) (h1 : I64.min ≤ x.val + y.val) (h2 : x.val + y.val ≤ I64.max) :
    ∃ z : I64, x + y = ok z ∧ z.val = x.val + y.val := by
  have := I64.add_spec h1 h2
  cases e : x + y with
  | ret z => rw [e, spec_ok] at this; exact ⟨z, rfl, this⟩
  | vis i k => rw [e, spec_vis] at this; exact this.elim
  | div => rw [e, spec_div] at this; exact this.elim

/-- `Σ_r (n_r + 2·g_r)`: what the range loop of `bath_half_dof` can add up, in units of 6. -/
def rangesWeight (gf : List Nat) : Nat → List Range → Nat
  | _, [] => 0
  | start, r :: rs => (r.lastAtom + 1 - start) + 2 * groupsStartingInModel gf start r.lastAtom
      + rangesWeight gf (r.lastAtom + 1) rs

def BathHalfDofInRange (ranges : List Range) (gf : List Nat) (cs : List (Nat × Nat)) : Prop :=
  6 * (rangesWeight gf 0 ranges : Int) + 2 * cs.length + 12 ≤ I64.max
    ∧ ∀ r, r ∈ ranges → r.lastAtom < Usize.max

theorem rangesHalfModel_drop (bath : Nat) (gf : List Nat) (rs : List Range) (i : Nat)
    (h : i < rs.length) (start : Nat) :
    rangesHalfModel bath gf start (rs.drop i) =
      let r := rs[i]
      let atoms : Int := ((r.lastAtom + 1 - start : Nat) : Int)
      let groups : Int := (groupsStartingInModel gf start r.lastAtom : Int)
      let rest := rangesHalfModel bath gf (r.lastAtom + 1) (rs.drop (i + 1))
      ((if r.comBath = bath then 6 * groups else 0) + rest.1,
       (if r.irBath = bath then 6 * atoms - 6 * groups else 0) + rest.2) := by
  rw [List.drop_eq_getElem_cons h]; rfl

theorem rangesWeight_drop (gf : List Nat) (rs : List Range) (i : Nat) (h : i < rs.length)
    (start : Nat) :
    rangesWeight gf start (rs.drop i) = (rs[i].lastAtom + 1 - start)
      + 2 * groupsStartingInModel gf start rs[i].lastAtom
      + rangesWeight gf (rs[i].lastAtom + 1) (rs.drop (i + 1)) := by
  rw [List.drop_eq_getElem_cons h]; rfl

set_option maxHeartbeats 4000000 in
theorem bath_half_dof_loop0_spec (bath : Usize) (ranges : Slice DofRange) (gf : Slice Usize)
    (cs : List (Nat × Nat)) (hb : BathHalfDofInRange (rangesOf ranges) (natsOf gf) cs) :
    bath_half_dof_loop0 bath ranges gf 0#i64 0#i64 0#usize 0#usize ⦃ p =>
      p.1.val = (rangesHalfModel bath.val (natsOf gf) 0 (rangesOf ranges)).1
      ∧ p.2.val = (rangesHalfModel bath.val (natsOf gf) 0 (rangesOf ranges)).2
      ∧ p.1.val.natAbs + p.2.val.natAbs ≤ 6 * rangesWeight (natsOf gf) 0 (rangesOf ranges) ⦄ := by
  obtain ⟨hw, hlast⟩ := hb
  unfold bath_half_dof_loop0
  apply loop.spec_decr_nat (fun p => ranges.length - p.2.2.2.val)
    (fun p => p.2.2.2.val ≤ ranges.length ∧
      p.1.val + (rangesHalfModel bath.val (natsOf gf) p.2.2.1.val ((rangesOf ranges).drop p.2.2.2.val)).1
        = (rangesHalfModel bath.val (natsOf gf) 0 (rangesOf ranges)).1 ∧
      p.2.1.val + (rangesHalfModel bath.val (natsOf gf) p.2.2.1.val ((rangesOf ranges).drop p.2.2.2.val)).2
        = (rangesHalfModel bath.val (natsOf gf) 0 (rangesOf ranges)).2 ∧
      p.1.val.natAbs + p.2.1.val.natAbs
        + 6 * rangesWeight (natsOf gf) p.2.2.1.val ((rangesOf ranges).drop p.2.2.2.val)
        ≤ 6 * rangesWeight (natsOf gf) 0 (rangesOf ranges))
  · intro ⟨com, ir, start, r⟩ ⟨hr, hcom, hir, hwt⟩
    unfold bath_half_dof_loop0.body
    simp only at hr hcom hir hwt ⊢
    step*
    · rename_i hlt
      have hl : r.val < (rangesOf ranges).length := by simp [rangesOf]; scalar_tac
      have hR : (rangesOf ranges)[r.val]'hl = toRange dr := by simp [rangesOf, dr_post]
      have hd := rangesHalfModel_drop bath.val (natsOf gf) (rangesOf ranges) r.val hl start.val
      have hwd := rangesWeight_drop (natsOf gf) (rangesOf ranges) r.val hl start.val
      have hlast' := hlast _ (List.getElem_mem hl)
      rw [hR] at hd hwd hlast'
      simp only [toRange] at hd hwd hlast'
      rw [hd] at hcom hir
      rw [hwd] at hwt
      have hgl := groupsStartingInModel_le (natsOf gf) start.val dr.last_atom.val
      have hlt' : r.val < ranges.length := by scalar_tac
      have hmin := i64_min_eq
      -- the range's atom count is `last − start + 1`, or 0 when the range is empty or backwards
      split
      on_goal 1 =>
        simp only at hcom hir
        step*
        have hxv := hcast_I64_val' x (by scalar_tac)
        generalize UScalar.hcast IScalarTy.I64 x = A at *
        have hA : A.val = ((dr.last_atom.val + 1 - start.val : Nat) : Int) := by
          rw [hxv]; scalar_tac
      on_goal 2 =>
        simp only at hcom hir
        step*
        generalize hA0 : (0#i64) = A at *
        have hA : A.val = ((dr.last_atom.val + 1 - start.val : Nat) : Int) := by
          rw [← hA0]; simp; scalar_tac
      all_goals
        have hgv : groups.val
            = (groupsStartingInModel (natsOf gf) start.val dr.last_atom.val : Int) := by
          rw [groups_post, hcast_I64_val' i1 (by scalar_tac), i1_post]
        have hcb' : (dr.com_bath = bath) ↔ dr.com_bath.val = bath.val := UScalar.eq_equiv _ _
        have hib' : (dr.ir_bath = bath) ↔ dr.ir_bath.val = bath.val := UScalar.eq_equiv _ _
        by_cases hcb : dr.com_bath.val = bath.val <;> by_cases hib : dr.ir_bath.val = bath.val <;>
          simp only [hcb, hib, hcb', hib', ite_true, ite_false] at hcom hir ⊢ <;> step*
        all_goals clear hd hwd hR dr_post hlast
        all_goals first
          | omega
          | (rw [start1_post, r1_post]; refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> omega)
    · rename_i hge
      have : r.val = ranges.length := by scalar_tac
      have hnil : (rangesOf ranges).drop r.val = [] := by simp [rangesOf, this]
      rw [hnil] at hcom hir hwt
      simp only [rangesHalfModel, rangesWeight] at hcom hir hwt
      refine ⟨by omega, by omega, by omega⟩
  · simp

theorem endsModel_drop (bath : Nat) (rs : List Range) (cs : List (Nat × Nat)) (i : Nat)
    (h : i < cs.length) :
    endsModel bath rs (cs.drop i) = endInModel bath rs cs[i].1 + endInModel bath rs cs[i].2
      + endsModel bath rs (cs.drop (i + 1)) := by
  rw [List.drop_eq_getElem_cons h]; rfl

/-- What one `if let Some((_, ir)) = baths_of(..)` adds to the constraint-end count. -/
def endOf (bath : Nat) : Option (Usize × Usize) → Int
  | none => 0
  | some (_, y) => if y.val = bath then 1 else 0

theorem endInModel_of (bath : Nat) (rs : List Range) (a : Nat) (o : Option (Usize × Usize))
    (h : pairVal o = bathsOfModel rs a) : endInModel bath rs a = endOf bath o := by
  unfold endInModel
  rw [← h]
  cases o with
  | none => rfl
  | some p => rfl

theorem endInModel_le (bath : Nat) (rs : List Range) (a : Nat) :
    0 ≤ endInModel bath rs a ∧ endInModel bath rs a ≤ 1 := by
  unfold endInModel
  split
  · split <;> omega
  · omega

@[step]
theorem end_in_bath_spec (ranges : Slice DofRange) (atom bath : Usize) :
    end_in_bath ranges atom bath ⦃ e =>
      e.val = endInModel bath.val (rangesOf ranges) atom.val ⦄ := by
  unfold end_in_bath
  step as ⟨o, o_post⟩
  rw [endInModel_of bath.val _ _ o o_post]
  rcases o with _ | ⟨x, y⟩
  · simp [endOf]
  · by_cases hy : y = bath
    · have hyv : y.val = bath.val := by rw [hy]
      simp [endOf, hy]
    · have hyv : ¬ y.val = bath.val := fun h => hy (UScalar.val_eq_imp _ _ h)
      simp [endOf, hy, hyv]

theorem bath_half_dof_loop1_spec (bath : Usize) (ranges : Slice DofRange) (cs : Slice Constraint)
    (hc : 2 * (cs.length : Int) ≤ I64.max) :
    bath_half_dof_loop1 bath ranges cs 0#i64 0#usize ⦃ e =>
      e.val = endsModel bath.val (rangesOf ranges) (consOf cs)
      ∧ 0 ≤ e.val ∧ e.val ≤ 2 * cs.length ⦄ := by
  unfold bath_half_dof_loop1
  apply loop.spec_decr_nat (fun p => cs.length - p.2.val)
    (fun p => p.2.val ≤ cs.length ∧ 0 ≤ p.1.val ∧ p.1.val ≤ 2 * (p.2.val : Int) ∧
      p.1.val + endsModel bath.val (rangesOf ranges) ((consOf cs).drop p.2.val)
        = endsModel bath.val (rangesOf ranges) (consOf cs))
  · intro ⟨ends, c⟩ ⟨hcl, h0, h2, hinv⟩
    unfold bath_half_dof_loop1.body
    simp only at hcl h0 h2 hinv ⊢
    split
    · rename_i hlt
      step as ⟨c1, c1_post⟩
      have hl : c.val < (consOf cs).length := by simp [consOf]; scalar_tac
      have hd := endsModel_drop bath.val (rangesOf ranges) (consOf cs) c.val hl
      have hci : ((consOf cs)[c.val]'hl) = (c1.i.val, c1.j.val) := by simp [consOf, c1_post]
      rw [hci] at hd
      simp only at hd
      have hlt' : c.val < cs.length := by scalar_tac
      have hmin := i64_min_eq
      have hd1 := endInModel_le bath.val (rangesOf ranges) c1.i.val
      have hd2 := endInModel_le bath.val (rangesOf ranges) c1.j.val
      step*
      all_goals scalar_tac
    · rename_i hge
      have : c.val = cs.length := by scalar_tac
      have hnil : (consOf cs).drop c.val = [] := by simp [consOf, this]
      rw [hnil] at hinv
      simp only [endsModel] at hinv
      simp only [spec_ok]
      refine ⟨by omega, by omega, by omega⟩
  · simp

/-- The roto-translational part of `bath_half_dof`: 6 halves off each bath of atom 0. -/
theorem rottrans_part_spec (com ir1 : I64) (c0 i0 bath : Usize)
    (hc : I64.min + 6 ≤ com.val) (hi : I64.min + 6 ≤ ir1.val) :
    (do
      let com1 ← if c0 = bath then com - 6#i64 else ok com
      if i0 = bath then do
          let ir2 ← ir1 - 6#i64
          ok ({ com := com1, ir := ir2 } : HalfDof)
        else ok ({ com := com1, ir := ir1 } : HalfDof)) ⦃ h =>
      h.com.val = (if c0.val = bath.val then com.val - 6 else com.val)
      ∧ h.ir.val = (if i0.val = bath.val then ir1.val - 6 else ir1.val) ⦄ := by
  have hc0 := UScalar.eq_equiv c0 bath
  have hi0 := UScalar.eq_equiv i0 bath
  by_cases e1 : c0.val = bath.val <;> by_cases e2 : i0.val = bath.val <;>
    simp only [e1, e2, ↓reduceIte] <;> split <;> rename_i h1 <;>
    (try exact absurd (hc0.mp h1) e1) <;> (try exact absurd (hc0.mpr e1) h1) <;> step*

theorem bath_half_dof_spec (bath : Usize) (ranges : Slice DofRange) (gf : Slice Usize)
    (cs : Slice Constraint) (rottrans : Bool)
    (hb : BathHalfDofInRange (rangesOf ranges) (natsOf gf) (consOf cs)) :
    bath_half_dof bath ranges gf cs rottrans ⦃ h =>
      (h.com.val, h.ir.val)
        = bathHalfDofModel bath.val (rangesOf ranges) (natsOf gf) (consOf cs) rottrans ⦄ := by
  have hlen : (consOf cs).length = cs.length := by simp [consOf]
  have hw := hb.1
  have hc : 2 * (cs.length : Int) ≤ I64.max := by rw [hlen] at hw; omega
  have hmin := i64_min_eq
  unfold bath_half_dof
  step with bath_half_dof_loop0_spec bath ranges gf (consOf cs) hb as ⟨com, ir, hcom, hir, habs⟩
  step with bath_half_dof_loop1_spec bath ranges cs hc as ⟨ends, he, he0, he2⟩
  step as ⟨ir1, hir1⟩
  unfold bathHalfDofModel
  simp only
  cases rottrans
  · simp only [Bool.false_eq_true, ite_false, spec_ok]
    simp only [Prod.mk.injEq]
    omega
  · simp only [ite_true]
    step*
    · rename_i ho
      rw [← o_post, ho]
      simp only [pairVal, Option.map_none, Prod.mk.injEq]
      omega
    · rename_i ho
      rw [← o_post, ho]
      obtain ⟨c0, i0⟩ := p
      simp only [pairVal, Option.map_some]
      apply spec_mono (rottrans_part_spec com ir1 c0 i0 bath (by scalar_tac) (by scalar_tac))
      intro h ⟨h1, h2⟩
      rw [h1, h2]
      by_cases e1 : c0.val = bath.val <;> by_cases e2 : i0.val = bath.val <;>
        simp only [e1, e2, ↓reduceIte, Prod.mk.injEq] <;> omega

theorem bath_half_dofs_spec (nBaths : Usize) (ranges : Slice DofRange) (gf : Slice Usize)
    (cs : Slice Constraint) (rottrans : Bool)
    (hb : BathHalfDofInRange (rangesOf ranges) (natsOf gf) (consOf cs)) :
    bath_half_dofs nBaths ranges gf cs rottrans ⦃ v =>
      v.val.map halfVal
        = bathHalfDofsModel nBaths.val (rangesOf ranges) (natsOf gf) (consOf cs) rottrans ⦄ := by
  unfold bath_half_dofs bath_half_dofs_loop
  apply loop.spec_decr_nat (fun p => nBaths.val - p.2.val)
    (fun p => p.2.val ≤ nBaths.val ∧ p.1.val.length = p.2.val ∧
      p.1.val.map halfVal = (List.range p.2.val).map (fun b =>
        bathHalfDofModel b (rangesOf ranges) (natsOf gf) (consOf cs) rottrans))
  · intro ⟨halves, b⟩ ⟨hbn, hlen, hinv⟩
    unfold bath_half_dofs_loop.body
    simp only at hbn hlen hinv ⊢
    split
    · rename_i hlt
      step with bath_half_dof_spec b ranges gf cs rottrans hb as ⟨h, hh⟩
      step as ⟨halves1, hpush⟩
      step as ⟨b1, hb1⟩
      refine ⟨by scalar_tac, by simp [hpush, hlen, hb1], ?_, by scalar_tac⟩
      rw [hpush, hb1, List.map_append, hinv, List.range_succ, List.map_append]
      simp [halfVal, hh]
    · rename_i hge
      simp only [spec_ok]
      have : b.val = nBaths.val := by scalar_tac
      rw [hinv, this]
      rfl
  · simp

/-- `Σ_b |h_b|`. -/
def absSum (h : List Int) : Nat := (h.map Int.natAbs).sum

/-- No `i64` overflow anywhere in `ndfmin_spread`. -/
def SpreadInRange (h : List Int) (ndfmin : Int) : Prop :=
  ((absSum h + 2 * ndfmin.natAbs + 2) * (absSum h + 2 * ndfmin.natAbs + 2) : Nat) ≤ I64.max

theorem absSum_drop (h : List Int) (i : Nat) (hi : i < h.length) :
    absSum (h.drop i) = h[i].natAbs + absSum (h.drop (i + 1)) := by
  rw [List.drop_eq_getElem_cons hi]; simp only [absSum, List.map_cons, List.sum_cons]

theorem sum_drop (h : List Int) (i : Nat) (hi : i < h.length) :
    (h.drop i).sum = h[i] + (h.drop (i + 1)).sum := by
  rw [List.drop_eq_getElem_cons hi, List.sum_cons]

theorem le_absSum (h : List Int) (i : Nat) (hi : i < h.length) : h[i].natAbs ≤ absSum h := by
  induction h generalizing i with
  | nil => simp at hi
  | cons x xs ih =>
    cases i with
    | zero => simp [absSum]
    | succ i =>
      have := ih i (by simp at hi; omega)
      simp only [absSum, List.map_cons, List.sum_cons, List.getElem_cons_succ] at this ⊢
      omega

theorem ndfmin_spread_loop0_spec (hs : Slice I64) (hb : 2 * (absSum (intsOf hs) : Int) ≤ I64.max) :
    ndfmin_spread_loop0 hs 0#i64 0#usize ⦃ t =>
      t.val = (intsOf hs).sum ∧ t.val.natAbs ≤ absSum (intsOf hs) ⦄ := by
  unfold ndfmin_spread_loop0
  apply loop.spec_decr_nat (fun p => hs.length - p.2.val)
    (fun p => p.2.val ≤ hs.length ∧ p.1.val + ((intsOf hs).drop p.2.val).sum = (intsOf hs).sum ∧
      p.1.val.natAbs + absSum ((intsOf hs).drop p.2.val) ≤ absSum (intsOf hs))
  · intro ⟨total, b⟩ ⟨hbl, hsum, habs⟩
    unfold ndfmin_spread_loop0.body
    simp only at hbl hsum habs ⊢
    have hmin := i64_min_eq
    split
    · rename_i hlt
      have hl : b.val < (intsOf hs).length := by simp [intsOf]; scalar_tac
      have hd := sum_drop (intsOf hs) b.val hl
      have ha := absSum_drop (intsOf hs) b.val hl
      step with Slice.index_usize_spec as ⟨x, hx⟩
      have hxv : (intsOf hs)[b.val]'hl = x.val := by simp [intsOf, hx]
      rw [hxv] at hd ha
      step*
      all_goals scalar_tac
    · rename_i hge
      simp only [spec_ok]
      have : b.val = hs.length := by scalar_tac
      have hnil : (intsOf hs).drop b.val = [] := by simp [intsOf, this]
      rw [hnil] at hsum habs
      have h0 : absSum ([] : List Int) = 0 := rfl
      rw [h0] at habs
      simp only [List.sum_nil, Int.add_zero] at hsum
      omega
  · simp

theorem le_absSum_of_mem (h : List Int) (x : Int) (hx : x ∈ h) : x.natAbs ≤ absSum h := by
  induction h with
  | nil => simp at hx
  | cons y ys ih =>
    simp only [absSum, List.map_cons, List.sum_cons] at ih ⊢
    rcases List.mem_cons.mp hx with e | e
    · subst e; omega
    · have := ih e; omega

theorem mul_in_range (x f : Int) (M : Nat) (hx : x.natAbs ≤ M) (hf : f.natAbs ≤ M)
    (hM : ((M * M : Nat) : Int) ≤ I64.max) : I64.min ≤ x * f ∧ x * f ≤ I64.max := by
  have h : (x * f).natAbs ≤ M * M := by rw [Int.natAbs_mul]; exact Nat.mul_le_mul hx hf
  have := i64_min_eq
  generalize x * f = p at *
  omega

theorem ndfmin_spread_loop1_spec (hs : Slice I64) (factor : I64)
    (hp : ∀ x, x ∈ intsOf hs → I64.min ≤ x * factor.val ∧ x * factor.val ≤ I64.max) :
    ndfmin_spread_loop1 hs factor (alloc.vec.Vec.new I64) 0#usize ⦃ v =>
      v.val.map (·.val) = (intsOf hs).map (· * factor.val) ⦄ := by
  unfold ndfmin_spread_loop1
  apply loop.spec_decr_nat (fun p => hs.length - p.2.val)
    (fun p => p.2.val ≤ hs.length ∧ p.1.val.length = p.2.val ∧
      p.1.val.map (·.val) ++ ((intsOf hs).drop p.2.val).map (· * factor.val)
        = (intsOf hs).map (· * factor.val))
  · intro ⟨nums, k⟩ ⟨hk, hlen, hinv⟩
    unfold ndfmin_spread_loop1.body
    simp only at hk hlen hinv ⊢
    split
    · rename_i hlt
      have hl : k.val < (intsOf hs).length := by simp [intsOf]; scalar_tac
      step with Slice.index_usize_spec as ⟨x, hx⟩
      have hxv : (intsOf hs)[k.val]'hl = x.val := by simp [intsOf, hx]
      have hpx := hp x.val (by rw [← hxv]; exact List.getElem_mem hl)
      step as ⟨y, hy⟩
      step as ⟨nums1, hpush⟩
      step as ⟨k1, hk1⟩
      refine ⟨by scalar_tac, by simp [hpush, hlen, hk1], ?_, by scalar_tac⟩
      rw [hpush, hk1, ← hinv, List.drop_eq_getElem_cons hl, hxv]
      simp [hy]
    · rename_i hge
      simp only [spec_ok]
      have : k.val = hs.length := by scalar_tac
      have hnil : (intsOf hs).drop k.val = [] := by simp [intsOf, this]
      rw [hnil] at hinv
      simpa using hinv
  · simp

theorem ndfmin_spread_spec (hs : Slice I64) (ndfmin : I64)
    (hb : SpreadInRange (intsOf hs) ndfmin.val) :
    ndfmin_spread hs ndfmin ⦃ s =>
      (s.numerators.val.map (·.val), s.denominator.val) = ndfminSpreadModel (intsOf hs) ndfmin.val ⦄ := by
  unfold SpreadInRange at hb
  generalize hA : absSum (intsOf hs) = A at hb
  have hM2 : 2 * (A + 2 * ndfmin.val.natAbs + 2) ≤ (A + 2 * ndfmin.val.natAbs + 2) * (A + 2 * ndfmin.val.natAbs + 2) :=
    Nat.mul_le_mul_right _ (by omega)
  have hmin := i64_min_eq
  unfold ndfmin_spread
  step with ndfmin_spread_loop0_spec hs (by omega) as ⟨total, htot, htabs⟩
  rw [hA] at htabs
  generalize hMdef : A + 2 * ndfmin.val.natAbs + 2 = M at hb hM2
  have hxs : ∀ x, x ∈ intsOf hs → x.natAbs ≤ M := fun x hx => by
    have := le_absSum_of_mem _ x hx; omega
  have hbm := hb
  generalize hMM : M * M = MM at hb hM2
  simp only [ndfminSpreadModel, ← htot]
  split
  · rename_i hpos
    have hpos' : 0 < total.val := by scalar_tac
    step as ⟨i, i_post⟩
    step as ⟨f1, f1_post⟩
    · scalar_tac
    step as ⟨d1, d1_post⟩
    have hf : f1.val.natAbs ≤ M := by scalar_tac
    step with ndfmin_spread_loop1_spec hs f1
      (fun x hx => mul_in_range x f1.val M (hxs x hx) hf hbm) as ⟨v, hv⟩
    simp only [hpos', ↓reduceIte, hv, f1_post, d1_post, i_post]
    exact ⟨trivial, trivial⟩
  · rename_i hpos
    have hpos' : ¬ 0 < total.val := by scalar_tac
    step with ndfmin_spread_loop1_spec hs 1#i64
      (fun x hx => mul_in_range x (1#i64).val M (hxs x hx) (by simp; omega) hbm) as ⟨v, hv⟩
    simp only [hpos', ↓reduceIte, hv]
    simp

/-! ### Rust = model, and Rust = meaning -/

/-- **Rust = model**: `baths_of` returns normally, with `bathsOfModel`. -/
theorem baths_of_eq (ranges : Slice DofRange) (atom : Usize) :
    ∃ o, baths_of ranges atom = ok o ∧ pairVal o = bathsOfModel (rangesOf ranges) atom.val :=
  exists_ok_of_spec (baths_of_spec ranges atom)

/-- **Rust = model**: `groups_starting_in` returns normally, with `groupsStartingInModel`. -/
theorem groups_starting_in_eq (gf : Slice Usize) (first last : Usize) :
    ∃ n, groups_starting_in gf first last = ok n
      ∧ n.val = groupsStartingInModel (natsOf gf) first.val last.val :=
  exists_ok_of_spec (groups_starting_in_spec gf first last)

/-- **Rust = model**: within the bound, `bath_half_dof` returns normally, with
`bathHalfDofModel`. -/
theorem bath_half_dof_eq (bath : Usize) (ranges : Slice DofRange) (gf : Slice Usize)
    (cs : Slice Constraint) (rottrans : Bool)
    (hb : BathHalfDofInRange (rangesOf ranges) (natsOf gf) (consOf cs)) :
    ∃ h, bath_half_dof bath ranges gf cs rottrans = ok h
      ∧ halfVal h = bathHalfDofModel bath.val (rangesOf ranges) (natsOf gf) (consOf cs) rottrans :=
  exists_ok_of_spec (bath_half_dof_spec bath ranges gf cs rottrans hb)

/-- **Rust = model**: within the bound, `bath_half_dofs` returns normally, with
`bathHalfDofsModel`. -/
theorem bath_half_dofs_eq (nBaths : Usize) (ranges : Slice DofRange) (gf : Slice Usize)
    (cs : Slice Constraint) (rottrans : Bool)
    (hb : BathHalfDofInRange (rangesOf ranges) (natsOf gf) (consOf cs)) :
    ∃ v, bath_half_dofs nBaths ranges gf cs rottrans = ok v
      ∧ v.val.map halfVal
        = bathHalfDofsModel nBaths.val (rangesOf ranges) (natsOf gf) (consOf cs) rottrans :=
  exists_ok_of_spec (bath_half_dofs_spec nBaths ranges gf cs rottrans hb)

/-- **Rust = model**: within the bound, `ndfmin_spread` returns normally, with
`ndfminSpreadModel`. -/
theorem ndfmin_spread_eq (hs : Slice I64) (ndfmin : I64)
    (hb : SpreadInRange (intsOf hs) ndfmin.val) :
    ∃ s, ndfmin_spread hs ndfmin = ok s
      ∧ (s.numerators.val.map (·.val), s.denominator.val)
        = ndfminSpreadModel (intsOf hs) ndfmin.val :=
  exists_ok_of_spec (ndfmin_spread_spec hs ndfmin hb)

open NewenTheory.Dynamics.DegreesOfFreedom in
/-- **Rust = meaning**: for ranges that cover the `n` atoms in order, `bath_half_dof` returns
the theory's half degrees of freedom of the bath, atom by atom. -/
theorem bath_half_dof_eq_theory {n : Nat} (bath : Usize) (ranges : Slice DofRange)
    (gf : Slice Usize) (cs : Slice Constraint) (rottrans : Bool)
    (hcov : Covers 0 (rangesOf ranges) n)
    (hb : BathHalfDofInRange (rangesOf ranges) (natsOf gf) (consOf cs)) :
    ∃ h, bath_half_dof bath ranges gf cs rottrans = ok h
      ∧ halfVal h = (comHalf (systemOf n (rangesOf ranges) (natsOf gf) (consOf cs) rottrans) bath.val,
          irHalf (systemOf n (rangesOf ranges) (natsOf gf) (consOf cs) rottrans) bath.val) := by
  obtain ⟨h, hr, hv⟩ := bath_half_dof_eq bath ranges gf cs rottrans hb
  exact ⟨h, hr, hv.trans (bathHalfDofModel_eq hcov _ _ _ _)⟩

open NewenTheory.Dynamics.DegreesOfFreedom in
/-- **Rust = meaning**: `ndfmin_spread` returns the theory's spread, which removes exactly
`ndfmin` in all (`spread_sum`). -/
theorem ndfmin_spread_eq_theory (hs : Slice I64) (ndfmin : I64)
    (hb : SpreadInRange (intsOf hs) ndfmin.val) :
    ∃ s, ndfmin_spread hs ndfmin = ok s
      ∧ (s.numerators.val.map (·.val), s.denominator.val)
        = (spreadNumerators (intsOf hs) ndfmin.val, spreadDenominator (intsOf hs)) := by
  obtain ⟨s, hr, hv⟩ := ndfmin_spread_eq hs ndfmin hb
  exact ⟨s, hr, hv.trans (ndfminSpreadModel_eq _ _)⟩

end Bath

/-! ## `topology`: pair kinds and group lookups

The loops of `contains_sorted`, `group_of` and `check_last_atoms` run by `loop.spec_decr_nat`
as in the `Bath` section: the model of the rest of the slice (`drop i`, or the descent search
from `i`) equals the model of the whole. Only `energy_cell` does arithmetic that can overflow;
its refinement assumes the cell of two in-range groups fits in `usize` (`EnergyCellInRange`),
which holds for fewer than `2^32` energy groups. -/

section Topology

open newen_logic.topology.pairs
open newen_logic.topology.groups
open NewenSpec.Topology.Pairs
open NewenSpec.Topology.Groups

theorem scanModel_drop (l : List Nat) (k : Nat) (h : k < l.length) (j : Nat) :
    scanModel (l.drop k) j = if l[k] = j then true else scanModel (l.drop (k + 1)) j := by
  rw [List.drop_eq_getElem_cons h]; rfl

@[step]
theorem contains_sorted_loop_spec (list : Slice Usize) (j n k : Usize) (hn : n.val = list.length)
    (hk : k.val ≤ n.val) :
    contains_sorted_loop list j n k ⦃ b => b = scanModel ((natsOf list).drop k.val) j.val ⦄ := by
  unfold contains_sorted_loop
  apply loop.spec_decr_nat (fun k' => n.val - k'.val)
    (fun k' => k'.val ≤ n.val ∧
      scanModel ((natsOf list).drop k'.val) j.val = scanModel ((natsOf list).drop k.val) j.val)
  · intro k' ⟨hk', hinv⟩
    unfold contains_sorted_loop.body
    have hlen : (natsOf list).length = list.length := by simp [natsOf]
    step*
    · rename_i hlt heq
      have hl : k'.val < (natsOf list).length := by scalar_tac
      rw [← hinv, scanModel_drop _ _ hl]
      have e : (natsOf list)[k'.val]'hl = j.val := by simp [natsOf, ← i_post, heq]
      simp [e]
    · rename_i hlt _ hne
      have hl : k'.val < (natsOf list).length := by scalar_tac
      refine ⟨by scalar_tac, ?_, by scalar_tac⟩
      rw [← hinv, scanModel_drop _ _ hl, k1_post]
      have e : ¬ (natsOf list)[k'.val]'hl = j.val := by
        simp only [natsOf, List.getElem_map, ← i_post]
        intro h; exact hne (UScalar.eq_of_val_eq h)
      simp [e]
    · rename_i hge
      have e : k'.val = (natsOf list).length := by scalar_tac
      rw [← hinv, e, List.drop_length]
      rfl
  · exact ⟨hk, rfl⟩

@[step]
theorem contains_sorted_spec (list : Slice Usize) (j : Usize) :
    contains_sorted list j ⦃ b => b = containsSortedModel (natsOf list) j.val ⦄ := by
  unfold contains_sorted
  have hlen : (natsOf list).length = list.length := by simp [natsOf]
  step*
  · have h0 : (natsOf list).length = 0 := by scalar_tac
    simp [containsSortedModel, h0]
  · rename_i hn _ hlt
    have h0 : ¬ (natsOf list).length = 0 := by scalar_tac
    have hl : (natsOf list)[(natsOf list).length - 1]'(by omega) = i1.val := by
      simp [natsOf, i1_post, i_post]
    unfold containsSortedModel
    rw [dite_eq_right h0, ite_eq_left (by rw [hl]; scalar_tac)]
  · rename_i hn _ hge
    have h0 : ¬ (natsOf list).length = 0 := by scalar_tac
    have hl : (natsOf list)[(natsOf list).length - 1]'(by omega) = i1.val := by
      simp [natsOf, i1_post, i_post]
    unfold containsSortedModel
    rw [dite_eq_right h0, ite_eq_right (by rw [hl]; scalar_tac), b_post, List.drop_zero]

/-- The Rust `PairKind` as the theory's. -/
def kindVal : newen_logic.topology.pairs.PairKind → NewenTheory.Topology.Pairs.PairKind
  | .Normal => .normal
  | .Excluded => .excluded
  | .OneFour => .oneFour
  | .AtomPair => .atomPair

@[step]
theorem pair_kind_spec (e f a : Slice Usize) (j : Usize) :
    pair_kind e f a j ⦃ k => kindVal k = pairKindModel (natsOf e) (natsOf f) (natsOf a) j.val ⦄ := by
  unfold pair_kind
  step*
  all_goals simp_all [pairKindModel, kindVal]

theorem groupOfLoop_drop (L : List Nat) (g : Nat) (h : g < L.length) (a : Nat) :
    groupOfLoop (L.drop g) a g =
      if a ≤ L[g] then some g else groupOfLoop (L.drop (g + 1)) a (g + 1) := by
  rw [List.drop_eq_getElem_cons h]; rfl

@[step]
theorem group_of_spec (L : Slice Usize) (atom : Usize) :
    group_of L atom ⦃ o => o.map (·.val) = groupOfModel (natsOf L) atom.val ⦄ := by
  unfold group_of group_of_loop
  have hlen : (natsOf L).length = L.length := by simp [natsOf]
  apply loop.spec_decr_nat (fun g => L.length - g.val)
    (fun g => g.val ≤ L.length ∧
      groupOfLoop ((natsOf L).drop g.val) atom.val g.val = groupOfModel (natsOf L) atom.val)
  · intro g ⟨hg, hinv⟩
    unfold group_of_loop.body
    step*
    · rename_i hlt hle
      have hl : g.val < (natsOf L).length := by scalar_tac
      rw [← hinv, groupOfLoop_drop _ _ hl]
      have hle' : atom.val ≤ (natsOf L)[g.val]'hl := by
        simp only [natsOf, List.getElem_map, ← i1_post]; scalar_tac
      simp [hle']
    · rename_i hlt _ hle
      have hl : g.val < (natsOf L).length := by scalar_tac
      refine ⟨by scalar_tac, ?_, by scalar_tac⟩
      rw [← hinv, groupOfLoop_drop _ _ hl]
      have hle' : ¬ atom.val ≤ (natsOf L)[g.val]'hl := by
        simp only [natsOf, List.getElem_map, ← i1_post]; scalar_tac
      simp only [hle', ite_false, g1_post]
    · rename_i hge
      have e : g.val = (natsOf L).length := by scalar_tac
      rw [← hinv, e, List.drop_length]
      rfl
  · simp [groupOfModel]

@[step]
theorem check_last_atoms_loop_spec (L : Slice Usize) (n g : Usize) (hn : n.val = L.length)
    (hg1 : 1 ≤ g.val) (hg : g.val ≤ n.val) :
    check_last_atoms_loop L n g ⦃ r => r.val = firstDescentFrom (natsOf L) g.val ⦄ := by
  unfold check_last_atoms_loop
  have hlen : (natsOf L).length = L.length := by simp [natsOf]
  apply loop.spec_decr_nat (fun g' => n.val - g'.val)
    (fun g' => 1 ≤ g'.val ∧ g'.val ≤ n.val ∧
      firstDescentFrom (natsOf L) g'.val = firstDescentFrom (natsOf L) g.val)
  · intro g' ⟨hg1', hg', hinv⟩
    unfold check_last_atoms_loop.body
    step*
    · rename_i hlt _ _ _ hinc
      have hl : g'.val < (natsOf L).length := by scalar_tac
      have hinc' : (natsOf L)[g'.val - 1]'(by omega) < (natsOf L)[g'.val]'hl := by
        simp only [natsOf, List.getElem_map]
        have e1 : (L.val[g'.val - 1]'(by scalar_tac)) = i1 := by simp [i1_post, i_post]
        rw [e1, ← i2_post]; scalar_tac
      refine ⟨by scalar_tac, by scalar_tac, ?_, by scalar_tac⟩
      rw [← hinv]
      conv_rhs => rw [firstDescentFrom]
      rw [dite_eq_left hl, ite_eq_left hinc', g1_post]
    · rename_i hlt _ _ hinc
      have hl : g'.val < (natsOf L).length := by scalar_tac
      have hinc' : ¬ (natsOf L)[g'.val - 1]'(by omega) < (natsOf L)[g'.val]'hl := by
        simp only [natsOf, List.getElem_map]
        have e1 : (L.val[g'.val - 1]'(by scalar_tac)) = i1 := by simp [i1_post, i_post]
        rw [e1, ← i2_post]; scalar_tac
      rw [← hinv, firstDescentFrom, dite_eq_left hl, ite_eq_right hinc']
    · rename_i hge
      have hl : ¬ g'.val < (natsOf L).length := by scalar_tac
      rw [← hinv, firstDescentFrom, dite_eq_right hl]
  · exact ⟨hg1, hg, rfl⟩

/-- The Rust `LastAtomsCheck` as the model's. -/
def checkVal : newen_logic.topology.groups.LastAtomsCheck → NewenSpec.Topology.Groups.LastAtomsCheck
  | .Ok => .ok
  | .Empty => .empty
  | .NotIncreasing g => .notIncreasing g.val
  | .LastNotFinal => .lastNotFinal

@[step]
theorem check_last_atoms_spec (L : Slice Usize) (n : Usize) :
    check_last_atoms L n ⦃ r => checkVal r = checkLastAtomsModel (natsOf L) n.val ⦄ := by
  unfold check_last_atoms
  have hlen : (natsOf L).length = L.length := by simp [natsOf]
  step*
  · have h0 : (natsOf L).length = 0 := by scalar_tac
    simp [checkLastAtomsModel, h0, checkVal]
  · rename_i h0 hlt
    have h0' : ¬ (natsOf L).length = 0 := by scalar_tac
    have hd : firstDescentFrom (natsOf L) 1 < (natsOf L).length := by scalar_tac
    simp [checkLastAtomsModel, h0', hd, checkVal, g_post]
  · rename_i h0 hge hn0
    have h0' : ¬ (natsOf L).length = 0 := by scalar_tac
    have hd : ¬ firstDescentFrom (natsOf L) 1 < (natsOf L).length := by scalar_tac
    simp [checkLastAtomsModel, h0', hd, checkVal]
  · rename_i h0 _ hge hn0 _ _ hne
    have h0' : ¬ (natsOf L).length = 0 := by scalar_tac
    have hd : ¬ firstDescentFrom (natsOf L) 1 < (natsOf L).length := by scalar_tac
    have hn0' : ¬ n.val = 0 := by scalar_tac
    have hl : (natsOf L)[(natsOf L).length - 1]'(by omega) = i1.val := by
      simp [natsOf, i1_post, i_post]
    have hne' : (natsOf L)[(natsOf L).length - 1]'(by omega) ≠ n.val - 1 := by
      rw [hl]; simp only [bne_iff_ne, ne_eq] at hne
      intro h; exact hne (UScalar.eq_of_val_eq (by scalar_tac))
    simp [checkLastAtomsModel, h0', hd, hn0', hne', checkVal]
  · rename_i h0 _ hge hn0 _ _ heq
    have h0' : ¬ (natsOf L).length = 0 := by scalar_tac
    have hd : ¬ firstDescentFrom (natsOf L) 1 < (natsOf L).length := by scalar_tac
    have hn0' : ¬ n.val = 0 := by scalar_tac
    have hl : (natsOf L)[(natsOf L).length - 1]'(by omega) = i1.val := by
      simp [natsOf, i1_post, i_post]
    have heq' : (natsOf L)[(natsOf L).length - 1]'(by omega) = n.val - 1 := by
      rw [hl]; simp only [bne_iff_ne, ne_eq, Decidable.not_not] at heq
      rw [heq]; scalar_tac
    simp [checkLastAtomsModel, h0', hd, hn0', heq', checkVal]

/-- No `usize` overflow in `energy_cell`: the cell of two groups in range fits. -/
def EnergyCellInRange (i j n : Nat) : Prop := i < n → j < n → i * n + j ≤ Usize.max

@[step]
theorem energy_cell_spec (i j n : Usize) (hb : EnergyCellInRange i.val j.val n.val) :
    energy_cell i j n ⦃ o => o.map (·.val) = energyCellModel i.val j.val n.val ⦄ := by
  unfold energy_cell
  step*
  all_goals first
    | (have := hb (by scalar_tac) (by scalar_tac)
       scalar_tac)
    | simp_all [energyCellModel]

/-! ### Rust = model, and Rust = meaning -/

/-- **Rust = model**: `contains_sorted` returns normally, with `containsSortedModel`. -/
theorem contains_sorted_eq (list : Slice Usize) (j : Usize) :
    contains_sorted list j = ok (containsSortedModel (natsOf list) j.val) :=
  eq_ok_of_spec (contains_sorted_spec list j)

/-- **Rust = meaning**: for an ascending slice, `contains_sorted` is membership. -/
theorem contains_sorted_iff (list : Slice Usize) (j : Usize)
    (h : (natsOf list).Pairwise (· ≤ ·)) :
    contains_sorted list j = ok true ↔ j.val ∈ natsOf list := by
  rw [contains_sorted_eq, Result.ok.injEq, containsSortedModel_iff h]

/-- **Rust = model**: `pair_kind` returns normally, with `pairKindModel`. -/
theorem pair_kind_eq (e f a : Slice Usize) (j : Usize) :
    ∃ k, pair_kind e f a j = ok k
      ∧ kindVal k = pairKindModel (natsOf e) (natsOf f) (natsOf a) j.val :=
  exists_ok_of_spec (pair_kind_spec e f a j)

/-- **Rust = meaning**: for ascending slices, `pair_kind` is the theory's kind by membership. -/
theorem pair_kind_eq_theory (e f a : Slice Usize) (j : Usize)
    (he : (natsOf e).Pairwise (· ≤ ·)) (hf : (natsOf f).Pairwise (· ≤ ·))
    (ha : (natsOf a).Pairwise (· ≤ ·)) :
    ∃ k, pair_kind e f a j = ok k
      ∧ kindVal k = NewenTheory.Topology.Pairs.pairKind (natsOf e) (natsOf f) (natsOf a) j.val := by
  obtain ⟨k, hr, hv⟩ := pair_kind_eq e f a j
  exact ⟨k, hr, hv.trans (pairKindModel_eq he hf ha j.val)⟩

/-- **Rust = model**: `group_of` returns normally, with `groupOfModel`. -/
theorem group_of_eq (L : Slice Usize) (atom : Usize) :
    ∃ o, group_of L atom = ok o ∧ o.map (·.val) = groupOfModel (natsOf L) atom.val :=
  exists_ok_of_spec (group_of_spec L atom)

/-- **Rust = meaning**: `group_of` returns the least group whose last atom reaches the atom. -/
theorem group_of_eq_theory (L : Slice Usize) (atom : Usize) :
    ∃ o, group_of L atom = ok o
      ∧ o.map (·.val) = NewenTheory.Topology.Groups.groupOf (natsOf L) atom.val := by
  obtain ⟨o, hr, hv⟩ := group_of_eq L atom
  exact ⟨o, hr, hv.trans (groupOfModel_eq _ _)⟩

/-- **Rust = model**: `check_last_atoms` returns normally, with `checkLastAtomsModel`. -/
theorem check_last_atoms_eq (L : Slice Usize) (n : Usize) :
    ∃ r, check_last_atoms L n = ok r ∧ checkVal r = checkLastAtomsModel (natsOf L) n.val :=
  exists_ok_of_spec (check_last_atoms_spec L n)

/-- **Rust = meaning**: `check_last_atoms` returns `Ok` exactly for lists of last atoms that
partition the `n` atoms into non-empty consecutive groups. -/
theorem check_last_atoms_ok_iff (L : Slice Usize) (n : Usize) :
    check_last_atoms L n = ok .Ok ↔ NewenTheory.Topology.Groups.LastAtomsOk (natsOf L) n.val := by
  obtain ⟨r, hr, hv⟩ := check_last_atoms_eq L n
  rw [hr, Result.ok.injEq, ← checkLastAtomsModel_iff, ← hv]
  cases r <;> simp [checkVal]

/-- **Rust = model**: within the bound, `energy_cell` returns normally, with
`energyCellModel`. -/
theorem energy_cell_eq (i j n : Usize) (hb : EnergyCellInRange i.val j.val n.val) :
    ∃ o, energy_cell i j n = ok o ∧ o.map (·.val) = energyCellModel i.val j.val n.val :=
  exists_ok_of_spec (energy_cell_spec i j n hb)

/-- **Rust = meaning**: within the bound, `energy_cell` returns the row-major cell of two
in-range groups, and none otherwise. -/
theorem energy_cell_eq_theory (i j n : Usize) (hb : EnergyCellInRange i.val j.val n.val) :
    ∃ o, energy_cell i j n = ok o
      ∧ o.map (·.val) = if i.val < n.val ∧ j.val < n.val
          then some (NewenTheory.Topology.Groups.energyCell n.val i.val j.val) else none := by
  obtain ⟨o, hr, hv⟩ := energy_cell_eq i j n hb
  exact ⟨o, hr, hv.trans (energyCellModel_eq _ _ _)⟩

end Topology


end NewenSpec.Refinement
