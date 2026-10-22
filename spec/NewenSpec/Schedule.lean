import NewenTheory.Dynamics.Schedule

/-!
# Models of `crates/newen-logic/src/schedule/`

Rust: `crates/newen-logic/src/schedule/mod.rs` (`output_misses_longrange`),
`schedule/output.rs` (`due`, `frames_in`, `minimisation_stop`, `final_step`) and `schedule/com.rs`
(`com_removal`); spec keys `schedule.*` in `spec/SPEC.toml`.
Theory: `NewenTheory/Dynamics/Schedule.lean` (`pairlistStep`, `RowsComplete`, `Misses`; `Due`,
`framesIn`, the top-of-loop and end-of-step frames, `finalStep`, `comRemoval`).

For each function this file has

1. the *model*, written line by line like the Rust function;
2. the *bridge*, proving the model means what the theory says.

The Aeneas refinements (Rust = model) are in `NewenSpec/Refinement.lean`.

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


/-! ## `schedule::output` -/

/-- `due`: `every = 0` first (output off), then "`step` is a multiple of `every`". -/
def dueModel (step every : Nat) : Bool :=
  if every = 0 then false else step % every == 0

/-- **`due` means `Due`.** -/
theorem dueModel_iff (step every : Nat) : dueModel step every = true ↔ Due every step := by
  unfold dueModel Due
  by_cases h : every = 0
  · simp [h]
  · simp [h]; omega

/-- `frames_in`: `every = 0` first, then the quotient, plus one for a partial last stretch. -/
def framesInModel (n every : Nat) : Nat :=
  if every = 0 then 0 else if n % every = 0 then n / every else n / every + 1

/-- **`frames_in` counts the due steps of `0 .. n − 1`** (`framesIn`, `framesIn_eq_card`). -/
theorem framesInModel_eq (n every : Nat) : framesInModel n every = framesIn n every := by
  unfold framesInModel
  by_cases h0 : every = 0
  · simp [h0, framesIn_zero]
  · have he : 0 < every := by omega
    rw [framesIn_eq_ceil n every he]
    have hdm := Nat.div_add_mod n every
    have hr := Nat.mod_lt n he
    simp only [h0, ite_false]
    generalize n / every = q at hdm
    by_cases hm : n % every = 0
    · simp only [hm, ite_true]
      symm
      apply Nat.div_eq_of_lt_le
      · rw [Nat.mul_comm]; omega
      · rw [Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega
    · simp only [hm, ite_false]
      symm
      apply Nat.div_eq_of_lt_le
      · rw [Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega
      · rw [Nat.add_mul, Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega

/-- `StopFrames`: the frames written at the stop step. -/
structure Stop where
  positions : Bool
  forces : Bool
  energies : Bool
  freeEnergy : Bool
deriving DecidableEq, Repr

/-- `minimisation_stop`: no stop unless converged after `nmin`; then the positions of `step` and
the lagged outputs of `step − 1`, each when due. -/
def minimisationStopModel (step nmin : Nat) (converged : Bool) (ntwx ntwf ntwe ntwg : Nat) :
    Option Stop :=
  if converged then
    if step ≤ nmin then none
    else some ⟨dueModel step ntwx, dueModel (step - 1) ntwf, dueModel (step - 1) ntwe,
      dueModel (step - 1) ntwg⟩
  else none

/-- A minimisation stops exactly when it has converged after `nmin` steps. -/
theorem minimisationStopModel_isSome (step nmin : Nat) (converged : Bool)
    (ntwx ntwf ntwe ntwg : Nat) :
    (minimisationStopModel step nmin converged ntwx ntwf ntwe ntwg).isSome = true
      ↔ converged = true ∧ nmin < step := by
  unfold minimisationStopModel
  cases converged <;> by_cases h : step ≤ nmin <;> simp [h] <;> omega

/-- **The stop writes gromosXX's frames.** An end-of-step writer that, on top of the steps
`0 .. step − 1`, writes the frames `minimisation_stop` names at the stop writes, for positions,
forces, energies and free energy each, exactly the frames gromosXX writes at the top of steps
`0 .. step` and in its final write. -/
theorem minimisationStopModel_frames {step nmin : Nat} {converged : Bool}
    {ntwx ntwf ntwe ntwg : Nat} {st : Stop}
    (h : minimisationStopModel step nmin converged ntwx ntwf ntwe ntwg = some st) :
    endOfStepPositions ntwx step st.positions = topOfLoopPositions ntwx step
      ∧ endOfStepLagged ntwf step st.forces = topOfLoopLagged ntwf step
      ∧ endOfStepLagged ntwe step st.energies = topOfLoopLagged ntwe step
      ∧ endOfStepLagged ntwg step st.freeEnergy = topOfLoopLagged ntwg step := by
  unfold minimisationStopModel at h
  cases converged
  · simp at h
  · by_cases hs : step ≤ nmin
    · simp [hs] at h
    · simp only [hs, ite_false, ite_true, Option.some.injEq] at h
      subst h
      have e : ∀ e s, dueModel s e = decide (Due e s) := fun e s => by
        rw [Bool.eq_iff_iff, dueModel_iff, decide_eq_true_iff]
      simp only [e]
      exact ⟨stop_positions_eq _ _, stop_lagged_eq _ _, stop_lagged_eq _ _, stop_lagged_eq _ _⟩

/-- `final_step`. -/
def finalStepModel (n : Nat) (stoppedAt : Option Nat) : Nat :=
  match stoppedAt with
  | some step => step
  | none => n

/-- **`final_step` is the theory's `finalStep`.** -/
theorem finalStepModel_eq (n : Nat) (stoppedAt : Option Nat) :
    finalStepModel n stoppedAt = finalStep n stoppedAt := by
  cases stoppedAt <;> rfl

/-! ## `schedule::com` -/

/-- `com_removal`: step 0 first, then `|nscm|`, then the multiple test. -/
def comRemovalModel (step : Nat) (nticom nscm : Int) (periodic : Bool) : Bool × Bool :=
  if step = 0 then (decide (nticom ≥ 1), decide (nticom ≥ 2))
  else
    let every := if nscm < 0 then -nscm else nscm
    if every = 0 then (false, false)
    else if step % every.toNat ≠ 0 then (false, false)
    else (true, decide (nscm < 0) && !periodic)

/-- **`com_removal` is the theory's `comRemoval`.** -/
theorem comRemovalModel_eq (step : Nat) (nticom nscm : Int) (periodic : Bool) :
    comRemovalModel step nticom nscm periodic = comRemoval step nticom nscm periodic := by
  unfold comRemovalModel comRemoval Due
  by_cases h0 : step = 0
  · simp [h0]
  · have ha : (if nscm < 0 then -nscm else nscm).toNat = nscm.natAbs := by split <;> omega
    simp only [h0, ite_false]
    by_cases hz : (if nscm < 0 then -nscm else nscm) = 0
    · have : nscm.natAbs = 0 := by split at hz <;> omega
      simp [hz, this]
    · have hpos : 0 < nscm.natAbs := by split at hz <;> omega
      simp only [hz, ite_false, ha]
      by_cases hm : step % nscm.natAbs = 0
      · cases periodic <;> simp [hm, hpos]
      · simp [hm]

end NewenSpec.Schedule
