/-!
# Pairlist steps and output steps

A GROMOS run with a twin-range cutoff evaluates the pairs within the short cutoff `RCUTP` on every
step. The pairs between `RCUTP` and the long cutoff `RCUTL` are evaluated only on a *pairlist
step*, when the pairlist is rebuilt, every `NSNB` steps (PAIRLIST block): steps `0, NSNB, 2 NSNB,
...`. Every step in between reuses the long-range interactions of the last pairlist step at or
before it.

An *output step* is a step on which something is written. Free-energy output (the λ derivatives
and, under PRECALCLAM, the energies at the extra λ values λ_k) is written every `NTWG` steps
(WRITETRAJ block): steps `0, NTWG, 2 NTWG, ...`.

Under PRECALCLAM, gromosXX adds the long-range part of the λ_k energies on a pairlist step only if
that step is also an output step. An output step then writes complete λ_k rows only if the
pairlist step it reuses was itself an output step. `RowsComplete` says every output step does;
`Misses` says some output step does not.

Only core Lean is used (`omega` and `Nat` lemmas), no Mathlib.
-/

namespace GromosTheory.Dynamics.Schedule

/-- The pairlist step whose long-range interactions step `s` reuses: the last multiple of `nsnb`
at or before `s`, `nsnb * ⌊s / nsnb⌋`. -/
def pairlistStep (nsnb s : Nat) : Nat := nsnb * (s / nsnb)

/-- Every output step (a multiple of `ntwg`) reuses a pairlist step that was also an output step,
so every written row carries its long-range part. -/
def RowsComplete (nsnb ntwg : Nat) : Prop :=
  ∀ s, ntwg ∣ s → ntwg ∣ pairlistStep nsnb s

/-- Some output step reuses a pairlist step that was not an output step: that row misses its
long-range part. -/
def Misses (nsnb ntwg : Nat) : Prop :=
  ∃ s, ntwg ∣ s ∧ ¬ ntwg ∣ pairlistStep nsnb s

/-- With `nsnb = 1` every step is a pairlist step, so nothing is reused and no row can miss. -/
theorem pairlistStep_one (s : Nat) : pairlistStep 1 s = s := by
  simp [pairlistStep]

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

end GromosTheory.Dynamics.Schedule
