import Aeneas
import NewenSpec.Schedule
import NewenSpec.Extracted.NewenLogic

/-!
# The extracted Rust computes the Lean model

`NewenSpec/Extracted/NewenLogic.lean` is what Aeneas made of
`crates/newen-logic/src/schedule.rs`: a pure Lean function over `Usize` (the machine integer,
`0 ≤ n < 2^64`) returning `Result Bool`, where `fail` would stand for a Rust panic (here: `%` by
zero) and `ok b` for a normal return of `b`.

The theorems below say: for every pair of `usize` inputs the Rust function returns normally
(it never panics) and returns exactly `missesModel`, hence (with `Schedule.missesModel_iff`)
exactly `Misses` of `NewenTheory/Dynamics/Schedule.lean`.
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

end NewenSpec.Refinement
