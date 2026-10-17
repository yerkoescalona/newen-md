import Aeneas
import NewenSpec.Schedule
import NewenSpec.Thermo.Dof
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
`thermo/dof.rs` (`total_dof`, over `U64` and `I64`) has its own section at the end.
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

end NewenSpec.Refinement
