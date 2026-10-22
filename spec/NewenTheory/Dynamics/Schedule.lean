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

The later sections: which steps write output frames and how many a run writes (`Due`,
`framesIn`), what a minimisation writes on the step it stops at (gromosXX writes at the top of each
step, an end-of-step writer at its end; both write the same frames), the step the final
configuration carries (`finalStep`), and on which steps the centre-of-mass motion is removed
(`comRemoval`).

Only core Lean is used (`omega` and `Nat` lemmas), no Mathlib.
-/

namespace NewenTheory.Dynamics.Schedule

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


/-! ## Output steps

Steps are 0-based: step `s` is the `s`-th step of the run, and its frame describes the
configuration at its start. An output written every `every` steps (WRITETRAJ's `NTWX`, `NTWF`,
`NTWE`, `NTWG`) writes a frame on the multiples of `every` (gromosXX `steps % every == 0`);
`every = 0` turns it off. -/

/-- A frame of an output written every `every` steps is due at step `s`. -/
def Due (every s : Nat) : Prop := 0 < every ∧ s % every = 0

instance (every s : Nat) : Decidable (Due every s) := inferInstanceAs (Decidable (_ ∧ _))

theorem due_iff_dvd (every s : Nat) : Due every s ↔ 0 < every ∧ every ∣ s := by
  simp [Due, Nat.dvd_iff_mod_eq_zero]

/-- The due steps among `steps`, in order. -/
def dueSteps (every : Nat) (steps : List Nat) : List Nat :=
  steps.filter (fun s => decide (Due every s))

/-- How many frames the steps `0 .. n − 1` write, counted step by step. -/
def framesIn (n every : Nat) : Nat :=
  match n with
  | 0 => 0
  | n + 1 => framesIn n every + if Due every n then 1 else 0

/-- **`framesIn` counts the due steps.** -/
theorem framesIn_eq_card (n every : Nat) :
    framesIn n every = (dueSteps every (List.range n)).length := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [framesIn, ih, dueSteps, dueSteps, List.range_succ, List.filter_append,
      List.length_append]
    by_cases h : Due every n <;> simp [h]

/-- One more step adds a frame exactly when it is due: `⌈(n + 1)/e⌉ = ⌈n/e⌉ + [e ∣ n]`. -/
theorem ceil_succ (n e : Nat) (he : 0 < e) :
    (n + 1 + e - 1) / e = (n + e - 1) / e + if n % e = 0 then 1 else 0 := by
  have hdm := Nat.div_add_mod n e
  have hr := Nat.mod_lt n he
  generalize n / e = q at hdm
  generalize n % e = r at hdm hr ⊢
  have hhi : (n + 1 + e - 1) / e = q + 1 := by
    apply Nat.div_eq_of_lt_le
    · rw [Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega
    · rw [Nat.add_mul, Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega
  rw [hhi]
  by_cases h0 : r = 0
  · have hlo : (n + e - 1) / e = q := by
      apply Nat.div_eq_of_lt_le
      · rw [Nat.mul_comm]; omega
      · rw [Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega
    simp [hlo, h0]
  · have hlo : (n + e - 1) / e = q + 1 := by
      apply Nat.div_eq_of_lt_le
      · rw [Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega
      · rw [Nat.add_mul, Nat.add_mul, Nat.one_mul, Nat.mul_comm]; omega
    simp [hlo, h0]

/-- **A run of `n` steps writes `⌈n / every⌉` frames** (none when `every = 0`). -/
theorem framesIn_eq_ceil (n every : Nat) (he : 0 < every) :
    framesIn n every = (n + every - 1) / every := by
  induction n with
  | zero =>
    simp only [framesIn, Nat.zero_add]
    exact (Nat.div_eq_of_lt (by omega)).symm
  | succ n ih =>
    rw [framesIn, ih, ceil_succ n every he]
    simp [Due, he]

theorem framesIn_zero (n : Nat) : framesIn n 0 = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => simp [framesIn, ih, Due]

/-! ## The step a minimisation stops at

gromosXX writes at the *top* of each step `s`: the positions of `s`, and for `s ≥ 1` the forces,
energies and free-energy frame of `s − 1` (the step just completed; *lagged* outputs). A
minimisation that stops at step `k` has gone through the top of steps `0 .. k`; it then writes its
final frames, the lagged outputs of `k − 1` once more when due (`out_configuration.cc:624-641`).

An *end-of-step* writer writes the positions and the lagged outputs of a step at its end. It
completes steps `0 .. k − 1`; at the stop it writes the positions of `k` and the lagged outputs of
`k − 1` if the flags `atStop` say so. With the flags `Due ntwx k` and `Due every (k − 1)` both
writers produce the same frames, in the same order, for every output. -/

/-- gromosXX: the position frames of a minimisation stopping at `k`. -/
def topOfLoopPositions (ntwx k : Nat) : List Nat := dueSteps ntwx (List.range (k + 1))

/-- gromosXX: the lagged frames (forces, energies, free energy) of a minimisation stopping at `k`:
those of `s − 1` at the top of each step `1 ≤ s ≤ k`, then the final one of `k − 1`. -/
def topOfLoopLagged (every k : Nat) : List Nat :=
  ((List.range (k + 1)).filter (fun s => decide (1 ≤ s ∧ Due every (s - 1)))).map (· - 1)
    ++ if Due every (k - 1) then [k - 1] else []

/-- An end-of-step writer: the position frames of steps `0 .. k − 1`, then that of `k` if
`atStop`. -/
def endOfStepPositions (ntwx k : Nat) (atStop : Bool) : List Nat :=
  dueSteps ntwx (List.range k) ++ if atStop then [k] else []

/-- An end-of-step writer: the lagged frames of steps `0 .. k − 1`, then that of `k − 1` again
if `atStop`. -/
def endOfStepLagged (every k : Nat) (atStop : Bool) : List Nat :=
  dueSteps every (List.range k) ++ if atStop then [k - 1] else []

/-- The top-of-loop lagged frames up to step `k` are the due steps below `k`. -/
theorem lagged_shift (every k : Nat) :
    ((List.range (k + 1)).filter (fun s => decide (1 ≤ s ∧ Due every (s - 1)))).map (· - 1)
      = dueSteps every (List.range k) := by
  induction k with
  | zero => simp [List.range_succ, dueSteps]
  | succ k ih =>
    rw [List.range_succ, List.filter_append, List.map_append, ih, dueSteps, dueSteps,
      List.range_succ (n := k), List.filter_append]
    by_cases h : Due every k <;> simp [h]

/-- **The stop writes gromosXX's position frames.** -/
theorem stop_positions_eq (ntwx k : Nat) :
    endOfStepPositions ntwx k (decide (Due ntwx k)) = topOfLoopPositions ntwx k := by
  unfold endOfStepPositions topOfLoopPositions dueSteps
  rw [List.range_succ, List.filter_append]
  by_cases h : Due ntwx k <;> simp [h]

/-- **The stop writes gromosXX's lagged frames** (forces, energies, free energy). -/
theorem stop_lagged_eq (every k : Nat) :
    endOfStepLagged every k (decide (Due every (k - 1))) = topOfLoopLagged every k := by
  unfold endOfStepLagged topOfLoopLagged
  rw [lagged_shift]
  by_cases h : Due every (k - 1) <;> simp [h]

/-! ## The final step

The final configuration carries the step counter after the loop (gromosXX `sim.steps()`,
`out_configuration.cc:519`): `n` after a run of `n` steps, `k` for a minimisation that stopped at
step `k` (it did not complete step `k`). -/

/-- The step the final configuration carries. -/
def finalStep (n : Nat) : Option Nat → Nat
  | some k => k
  | none => n

/-- The step counter after a loop over the steps `0 .. n − 1` that, at the top of step `s`,
stops if `stop s` and otherwise integrates step `s` and moves on: from step `s`, with `m` steps
left. -/
def loopCounter (stop : Nat → Bool) : (s m : Nat) → Nat
  | s, 0 => s
  | s, m + 1 => if stop s then s else loopCounter stop (s + 1) m

theorem loopCounter_eq (stop : Nat → Bool) (s m : Nat) :
    loopCounter stop s m = finalStep (s + m) ((List.range' s m).find? stop) := by
  induction m generalizing s with
  | zero => simp [loopCounter, finalStep]
  | succ m ih =>
    rw [loopCounter, List.range'_succ, List.find?_cons]
    cases h : stop s
    · simp only [Bool.false_eq_true, ite_false, ih]
      congr 1
      omega
    · simp [finalStep]

/-- **The final step is the number of steps integrated.** A loop over `n` steps that stops at
the top of the first step `k` with `stop k` has integrated `k` steps, and `finalStep` of that
first stop is `k`; a loop that never stops has integrated `n` steps, and `finalStep n none` is
`n`. -/
theorem finalStep_eq_loopCounter (stop : Nat → Bool) (n : Nat) :
    finalStep n ((List.range n).find? stop) = loopCounter stop 0 n := by
  rw [loopCounter_eq, Nat.zero_add, List.range_eq_range']

/-! ## Centre-of-mass motion removal

Step 0 follows INITIALISE's `NTICOM` (≥ 1: translation, ≥ 2: also rotation). Every later step
follows COMTRANSROT's `NSCM`: on the multiples of `|NSCM|`, the translation when `NSCM ≠ 0`, the
rotation too when `NSCM < 0`, except in a periodic box, where the rotation is never removed after
step 0 (gromosXX `remove_com_motion.cc:240-262`, `in_parameter.cc:1306-1317`,
`configuration.cc:555-560`). -/

/-- `(translation, rotation)` removed at step `s`. -/
def comRemoval (s : Nat) (nticom nscm : Int) (periodic : Bool) : Bool × Bool :=
  if s = 0 then (decide (1 ≤ nticom), decide (2 ≤ nticom))
  else (decide (Due nscm.natAbs s), decide (Due nscm.natAbs s ∧ nscm < 0 ∧ periodic = false))

/-- **Step 0 follows `NTICOM`**, whatever `NSCM` and the box. -/
theorem com_step0 (nticom nscm : Int) (periodic : Bool) :
    comRemoval 0 nticom nscm periodic = (decide (1 ≤ nticom), decide (2 ≤ nticom)) := by
  simp [comRemoval]

/-- **Later steps follow `NSCM`**: the translation on the multiples of `|NSCM| > 0`, the rotation
on the same steps when `NSCM < 0` and the box is not periodic. -/
theorem com_period (s : Nat) (hs : 0 < s) (nticom nscm : Int) (periodic : Bool) :
    ((comRemoval s nticom nscm periodic).1 = true ↔ nscm ≠ 0 ∧ nscm.natAbs ∣ s)
    ∧ ((comRemoval s nticom nscm periodic).2 = true
        ↔ nscm < 0 ∧ periodic = false ∧ nscm.natAbs ∣ s) := by
  have hs' : s ≠ 0 := by omega
  simp only [comRemoval, hs', ite_false, decide_eq_true_eq, due_iff_dvd]
  constructor
  · constructor
    · intro ⟨h0, hd⟩; exact ⟨by omega, hd⟩
    · intro ⟨h0, hd⟩; exact ⟨by omega, hd⟩
  · constructor
    · intro ⟨⟨_, hd⟩, hn, hp⟩; exact ⟨hn, hp, hd⟩
    · intro ⟨hn, hp, hd⟩; exact ⟨⟨by omega, hd⟩, hn, hp⟩

end NewenTheory.Dynamics.Schedule
