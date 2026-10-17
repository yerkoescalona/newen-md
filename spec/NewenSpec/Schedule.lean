import NewenTheory.Dynamics.Schedule

/-!
# Model of `crates/newen-logic/src/schedule.rs`

Rust: `crates/newen-logic/src/schedule.rs`, `output_misses_longrange`
(spec key `schedule.output_misses_longrange` in `spec/SPEC.toml`).
Theory: `NewenTheory/Dynamics/Schedule.lean` (`pairlistStep`, `RowsComplete`, `Misses`).

This file has two parts:

1. the *model*: `missesModel`, a closed form written line by line like the Rust function;
2. the *bridge*: `missesModel_iff`, proving the model means `Misses` of the theory.

Only core Lean is used (`omega` and `Nat` lemmas), no Mathlib.
-/

namespace NewenSpec.Schedule

open NewenTheory.Dynamics.Schedule

/-! ## 1. Model (mirrors the Rust function) -/

/-- Closed form of `output_misses_longrange`, written as the Rust code is:
`nsnb ≤ 1` or `ntwg = 0` first (no reuse, or no output), then "neither divides the other". -/
def missesModel (nsnb ntwg : Nat) : Bool :=
  if nsnb ≤ 1 then false
  else if ntwg = 0 then false
  else ntwg % nsnb != 0 && nsnb % ntwg != 0

/-! ## 2. Bridge: the model means `Misses` -/

/-- The model in words: `nsnb > 1`, `ntwg > 0`, and neither divides the other. -/
theorem missesModel_eq_true_iff (nsnb ntwg : Nat) :
    missesModel nsnb ntwg = true ↔
      1 < nsnb ∧ 0 < ntwg ∧ ¬ nsnb ∣ ntwg ∧ ¬ ntwg ∣ nsnb := by
  unfold missesModel
  by_cases h1 : nsnb ≤ 1
  · simp only [h1, ite_true, Bool.false_eq_true, false_iff]; omega
  · by_cases h0 : ntwg = 0
    · simp only [h1, h0, ite_true, ite_false, Bool.false_eq_true, false_iff]; omega
    · simp only [h1, h0, ite_false, Nat.dvd_iff_mod_eq_zero, bne_iff_ne, Bool.and_eq_true]
      omega

/-- **The model means what the Rust doc comment says**, for every `nsnb` and `ntwg`.

The edge cases agree without hypotheses: `nsnb = 1` (every step is a pairlist step,
`pairlistStep_one`), `nsnb = 0` (the pairlist step is always `0`, an output step) and `ntwg = 0`
(the only "output step" is `0`, itself a pairlist step) all give no miss, and the model says
`false` for each. -/
theorem missesModel_iff (nsnb ntwg : Nat) :
    missesModel nsnb ntwg = true ↔ Misses nsnb ntwg := by
  rw [missesModel_eq_true_iff, misses_iff_not_complete]
  by_cases hs : 1 < nsnb
  · by_cases hw : 0 < ntwg
    · rw [rows_complete_iff (by omega) hw]
      constructor
      · intro ⟨_, _, h1, h2⟩ h
        cases h with
        | inl h => exact h1 h
        | inr h => exact h2 h
      · intro h
        exact ⟨hs, hw, fun hd => h (Or.inl hd), fun hd => h (Or.inr hd)⟩
    · -- ntwg = 0: only step 0 is an "output step", and it is its own pairlist step.
      have h0 : ntwg = 0 := by omega
      subst h0
      constructor
      · intro ⟨_, h, _⟩; omega
      · intro h
        apply False.elim (h _)
        intro s hs'
        rw [Nat.eq_zero_of_zero_dvd hs']
        simp [pairlistStep]
  · -- nsnb ≤ 1: the pairlist step is `s` (nsnb = 1) or `0` (nsnb = 0); both are output steps.
    constructor
    · intro ⟨h, _⟩; omega
    · intro h
      apply False.elim (h _)
      intro s hs'
      unfold pairlistStep
      cases Nat.lt_or_ge nsnb 1 with
      | inl h0 => rw [show nsnb = 0 by omega]; simp
      | inr h1 => rw [show nsnb = 1 by omega, Nat.one_mul, Nat.div_one]; exact hs'

end NewenSpec.Schedule
