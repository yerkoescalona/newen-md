/-!
# Pairlist and free-energy output schedules

Rust: `crates/gromos-logic/src/schedule.rs`, `output_misses_longrange`
(spec key `schedule.output_misses_longrange` in `spec/SPEC.toml`).

An MD run updates its twin-range pairlist every `nsnb` steps and writes free-energy output every
`ntwg` steps. Between pairlist updates the long-range interactions are reused from the last
pairlist step. Under PRECALCLAM gromosXX adds the long-range λ_k energies to an output row only if
that pairlist step was itself an output step
(`perturbed_nonbonded_set.cc:86-96`, `perturbed_nonbonded_innerloop.cc:315`).

This file has three parts:

1. the *meaning*: `pairlistStep`, `RowsComplete`, `Misses`, stated over all steps;
2. the *model*: `missesModel`, a closed form written line by line like the Rust function;
3. the *bridge*: `rows_complete_iff` and `missesModel_iff`, proving the model means (1).

Only core Lean is used (`omega` and `Nat` lemmas), no Mathlib.
-/

namespace GromosSpec.Schedule

/-! ## 1. Meaning -/

/-- The pairlist step whose long-range interactions step `s` reuses: the last multiple of `nsnb`
at or before `s`, `nsnb * ⌊s / nsnb⌋`. -/
def pairlistStep (nsnb s : Nat) : Nat := nsnb * (s / nsnb)

/-- Every output step (a multiple of `ntwg`) reuses a pairlist step that was also an output step,
so every written row carries its long-range part. -/
def RowsComplete (nsnb ntwg : Nat) : Prop :=
  ∀ s, ntwg ∣ s → ntwg ∣ pairlistStep nsnb s

/-- Some output step reuses a pairlist step that was not an output step: that row misses its
long-range part. This is what the Rust doc comment says the function decides. -/
def Misses (nsnb ntwg : Nat) : Prop :=
  ∃ s, ntwg ∣ s ∧ ¬ ntwg ∣ pairlistStep nsnb s

/-- With `nsnb = 1` every step is a pairlist step, so nothing is reused and no row can miss. -/
theorem pairlistStep_one (s : Nat) : pairlistStep 1 s = s := by
  simp [pairlistStep]

/-! ## 2. Model (mirrors the Rust function) -/

/-- Closed form of `output_misses_longrange`, written as the Rust code is:
`nsnb ≤ 1` or `ntwg = 0` first (no reuse, or no output), then "neither divides the other". -/
def missesModel (nsnb ntwg : Nat) : Bool :=
  if nsnb ≤ 1 then false
  else if ntwg = 0 then false
  else ntwg % nsnb != 0 && nsnb % ntwg != 0

/-! ## 3. Bridge: the model means `Misses` -/

/-- Every output step reuses an output pairlist step iff `nsnb ∣ ntwg` or `ntwg ∣ nsnb`. -/
theorem rows_complete_iff {nsnb ntwg : Nat} (hs : 0 < nsnb) (hw : 0 < ntwg) :
    RowsComplete nsnb ntwg ↔ nsnb ∣ ntwg ∨ ntwg ∣ nsnb := by
  unfold RowsComplete pairlistStep
  constructor
  · intro h
    apply Classical.byContradiction
    intro hn
    have hsw : ¬ nsnb ∣ ntwg := fun hd => hn (Or.inl hd)
    have hws : ¬ ntwg ∣ nsnb := fun hd => hn (Or.inr hd)
    have hr : ntwg % nsnb ≠ 0 := fun h0 => hsw (Nat.dvd_of_mod_eq_zero h0)
    have hr' : nsnb % ntwg ≠ 0 := fun h0 => hws (Nat.dvd_of_mod_eq_zero h0)
    have hne : nsnb ≠ ntwg := fun he => hsw (he ▸ Nat.dvd_refl nsnb)
    rcases Nat.lt_or_gt_of_ne hne with hlt | hgt
    · -- nsnb < ntwg: the output step `ntwg` reuses `p = nsnb ⌊ntwg/nsnb⌋`, with `0 < p < ntwg`.
      have hp := h ntwg (Nat.dvd_refl ntwg)
      have hdm := Nat.div_add_mod ntwg nsnb
      have hrl := Nat.mod_lt ntwg hs
      have hpos : 0 < nsnb * (ntwg / nsnb) :=
        Nat.mul_pos hs (Nat.div_pos (Nat.le_of_lt hlt) hs)
      have hle := Nat.le_of_dvd hpos hp
      omega
    · -- ntwg < nsnb: the first output step after `nsnb` lies in `(nsnb, nsnb + ntwg]` and
      -- reuses `p = nsnb`, which is not an output step.
      have hdm := Nat.div_add_mod nsnb ntwg
      have hrl := Nat.mod_lt nsnb hw
      have hs1 : nsnb < ntwg * (nsnb / ntwg + 1) := by
        rw [Nat.mul_add, Nat.mul_one]; omega
      have hs2 : ntwg * (nsnb / ntwg + 1) < 2 * nsnb := by
        rw [Nat.mul_add, Nat.mul_one]; omega
      have hq : ntwg * (nsnb / ntwg + 1) / nsnb = 1 :=
        Nat.div_eq_of_lt_le (by omega) (by omega)
      have hp := h (ntwg * (nsnb / ntwg + 1)) (Nat.dvd_mul_right ntwg _)
      rw [hq, Nat.mul_one] at hp
      exact hws hp
  · intro hd s hs'
    cases hd with
    | inl hd => rw [Nat.mul_div_cancel' (Nat.dvd_trans hd hs')]; exact hs'
    | inr hd => exact Nat.dvd_trans hd (Nat.dvd_mul_right nsnb _)

/-- `Misses` is exactly the failure of `RowsComplete`. -/
theorem misses_iff_not_complete (nsnb ntwg : Nat) :
    Misses nsnb ntwg ↔ ¬ RowsComplete nsnb ntwg := by
  unfold Misses RowsComplete
  constructor
  · intro ⟨s, hs, hn⟩ h
    exact hn (h s hs)
  · intro h
    apply Classical.byContradiction
    intro hne
    exact h fun s hs => Classical.byContradiction fun hn => hne ⟨s, hs, hn⟩

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

end GromosSpec.Schedule
