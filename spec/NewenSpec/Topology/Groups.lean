import NewenTheory.Topology.Groups
import NewenSpec.Thermo.Bath

/-!
# Model of `crates/newen-logic/src/topology/groups.rs`

Rust: `crates/newen-logic/src/topology/groups.rs`: `group_of`, `check_last_atoms`,
`energy_cell` (spec keys `topology.groups.*` in `spec/SPEC.toml`).
Theory: `NewenTheory/Topology/Groups.lean` (`groupOf`, `LastAtomsOk`, `energyCell`).

1. The *models*, written as the Rust is: `group_of` walks the groups from the first;
   `check_last_atoms` looks for the first entry that does not exceed its predecessor, then checks
   the last entry.
2. The *bridges*: `group_of` is the theory's least index (`groupOfModel_eq`),
   `check_last_atoms` accepts exactly the lists that partition the atoms
   (`checkLastAtomsModel_iff`) and reports each refusal in the Rust's check order
   (`checkLastAtomsModel_empty_iff`, `_notIncreasing_iff`, `_lastNotFinal_iff`), `energy_cell` is the row-major cell (`energyCellModel_eq`); and
   `thermo.bath.baths_of` is the same lookup on the DOFSET ranges (`bathsOf_eq_groupOf`).

The Aeneas refinements (Rust = model) are in `NewenSpec/Refinement.lean`. Core Lean only.
-/

namespace NewenSpec.Topology.Groups

open NewenTheory.Topology.Groups

/-! ## 1. Models (mirror the Rust functions) -/

/-- The loop of `group_of` from group `g` on, over the remaining last atoms. -/
def groupOfLoop : List Nat → Nat → Nat → Option Nat
  | [], _, _ => none
  | l :: ls, atom, g => if atom ≤ l then some g else groupOfLoop ls atom (g + 1)

/-- `group_of`: the first group whose last atom reaches `atom`. -/
def groupOfModel (L : List Nat) (atom : Nat) : Option Nat := groupOfLoop L atom 0

/-- The result of `check_last_atoms` (`LastAtomsCheck`). -/
inductive LastAtomsCheck where
  | ok
  | empty
  | notIncreasing (g : Nat)
  | lastNotFinal
  deriving DecidableEq, Repr

/-- The loop of `check_last_atoms`: from index `g` on, the first index whose entry does not
exceed its predecessor, or the length. -/
def firstDescentFrom (L : List Nat) (g : Nat) : Nat :=
  if h : g < L.length then
    if L[g - 1] < L[g] then firstDescentFrom L (g + 1) else g
  else g
termination_by L.length - g

/-- `check_last_atoms`: empty, else the first descent, else the last entry against `n − 1`. -/
def checkLastAtomsModel (L : List Nat) (n : Nat) : LastAtomsCheck :=
  if h : L.length = 0 then .empty
  else if firstDescentFrom L 1 < L.length then .notIncreasing (firstDescentFrom L 1)
  else if n = 0 then .lastNotFinal
  else if L[L.length - 1] ≠ n - 1 then .lastNotFinal
  else .ok

/-- `energy_cell`: the row-major cell when both groups are in range. -/
def energyCellModel (i j n : Nat) : Option Nat :=
  if i < n then (if j < n then some (i * n + j) else none) else none

/-! ## 2. Bridges -/

theorem groupOfLoop_eq (L : List Nat) (a g : Nat) :
    groupOfLoop L a g = (groupOf L a).map (· + g) := by
  induction L generalizing g with
  | nil => simp [groupOfLoop, groupOf]
  | cons l ls ih =>
    unfold groupOf at ih ⊢
    rw [List.findIdx?_cons]
    by_cases h : a ≤ l
    · simp [groupOfLoop, h]
    · simp only [groupOfLoop, h, ite_false, decide_false, Bool.false_eq_true, ih, Option.map_map]
      congr 1
      funext i
      simp only [Function.comp_apply]
      omega

/-- **`group_of` is the theory's least index.** -/
theorem groupOfModel_eq (L : List Nat) (a : Nat) : groupOfModel L a = groupOf L a := by
  simp [groupOfModel, groupOfLoop_eq]

theorem firstDescentFrom_ge (L : List Nat) (g : Nat) : g ≤ firstDescentFrom L g := by
  rw [firstDescentFrom]
  by_cases h : g < L.length
  · simp only [h, dite_true]
    split
    · have := firstDescentFrom_ge L (g + 1); omega
    · exact Nat.le_refl g
  · simp only [h, dite_false]; exact Nat.le_refl g
termination_by L.length - g

/-- The loop runs to the end exactly when every entry from `g` on exceeds its predecessor. -/
theorem firstDescentFrom_spec (L : List Nat) (g : Nat) (hg : 1 ≤ g) :
    L.length ≤ firstDescentFrom L g ↔ ∀ k (hk : k < L.length), g ≤ k → L[k - 1] < L[k] := by
  rw [firstDescentFrom]
  by_cases h : g < L.length
  · simp only [h, dite_true]
    by_cases hd : L[g - 1] < L[g]
    · simp only [hd, ite_true]
      rw [firstDescentFrom_spec L (g + 1) (by omega)]
      constructor
      · intro H k hk hgk
        by_cases e : k = g
        · subst e; exact hd
        · exact H k hk (by omega)
      · intro H k hk hgk
        exact H k hk (by omega)
    · simp only [hd, ite_false]
      constructor
      · intro H; omega
      · intro H; exact absurd (H g h (Nat.le_refl g)) hd
  · simp only [h, dite_false]
    constructor
    · intro _ k hk hgk; omega
    · intro _; omega
termination_by L.length - g

/-- **`check_last_atoms` accepts exactly the lists that partition the atoms.** -/
theorem checkLastAtomsModel_iff (L : List Nat) (n : Nat) :
    checkLastAtomsModel L n = .ok ↔ LastAtomsOk L n := by
  have spec := firstDescentFrom_spec L 1 (Nat.le_refl 1)
  unfold checkLastAtomsModel
  constructor
  · intro hok
    by_cases h0 : L.length = 0
    · simp [h0] at hok
    · simp only [h0, dite_false] at hok
      by_cases hd : firstDescentFrom L 1 < L.length
      · simp [hd] at hok
      · simp only [hd, ite_false] at hok
        by_cases hn : n = 0
        · simp [hn] at hok
        · simp only [hn, ite_false] at hok
          by_cases hl : L[L.length - 1] = n - 1
          · have hinc := spec.mp (by omega)
            refine ⟨by omega, fun g hg => ?_, ?_, by omega⟩
            · have := hinc (g + 1) hg (by omega)
              simpa using this
            · rw [List.getLast?_eq_getElem?, List.getElem?_eq_getElem (by omega), hl]
          · simp [hl] at hok
  · intro h
    have hl := h.last_eq
    have hlen := h.nonempty
    have hn := h.atoms
    have hd : ¬ firstDescentFrom L 1 < L.length := by
      have := spec.mpr (fun k hk hk1 => by
        have := h.increasing (k - 1) (by omega)
        simpa [show k - 1 + 1 = k by omega] using this)
      omega
    simp [show L.length ≠ 0 by omega, hd, show n ≠ 0 by omega, hl]

/-- Every index the loop passes is an ascent. -/
theorem firstDescentFrom_asc (L : List Nat) (g : Nat) :
    ∀ k (hk : k < L.length), g ≤ k → k < firstDescentFrom L g → L[k - 1] < L[k] := by
  intro k hk hgk hkd
  rw [firstDescentFrom] at hkd
  by_cases h : g < L.length
  · simp only [h, dite_true] at hkd
    by_cases hd : L[g - 1] < L[g]
    · simp only [hd, ite_true] at hkd
      by_cases e : k = g
      · subst e; exact hd
      · exact firstDescentFrom_asc L (g + 1) k hk (by omega) hkd
    · simp only [hd, ite_false] at hkd; omega
  · simp only [h, dite_false] at hkd; omega
termination_by L.length - g

/-- Where the loop stops inside the list, the entry does not exceed its predecessor. -/
theorem firstDescentFrom_desc (L : List Nat) (g d : Nat) (e : firstDescentFrom L g = d)
    (hd : d < L.length) : ¬ L[d - 1] < L[d] := by
  rw [firstDescentFrom] at e
  by_cases hg : g < L.length
  · simp only [hg, dite_true] at e
    by_cases ha : L[g - 1] < L[g]
    · simp only [ha, ite_true] at e
      exact firstDescentFrom_desc L (g + 1) d e hd
    · simp only [ha, ite_false] at e
      subst e
      exact ha
  · simp only [hg, dite_false] at e
    omega
termination_by L.length - g

/-- **`check_last_atoms` reports `Empty` exactly for the empty list.** -/
theorem checkLastAtomsModel_empty_iff (L : List Nat) (n : Nat) :
    checkLastAtomsModel L n = .empty ↔ L = [] := by
  unfold checkLastAtomsModel
  by_cases h0 : L.length = 0
  · simp [List.length_eq_zero_iff.mp h0]
  · have : L ≠ [] := fun e => h0 (by simp [e])
    simp only [h0, dite_false, this, iff_false]
    repeat' split
    all_goals simp

/-- **`check_last_atoms` reports `NotIncreasing g` exactly when `g` is the first index whose entry
does not exceed its predecessor** (`1 ≤ g`, `L[g − 1] ≥ L[g]`, ascents before it). -/
theorem checkLastAtomsModel_notIncreasing_iff (L : List Nat) (n g : Nat) :
    checkLastAtomsModel L n = .notIncreasing g ↔
      ∃ hg : g < L.length, 1 ≤ g ∧ L[g] ≤ L[g - 1] ∧
        ∀ k (hk : k < g), 1 ≤ k → L[k - 1] < L[k] := by
  have hasc := firstDescentFrom_asc L 1
  have hge := firstDescentFrom_ge L 1
  unfold checkLastAtomsModel
  by_cases h0 : L.length = 0
  · simp only [h0, dite_true, reduceCtorEq, false_iff]
    intro ⟨hg, _⟩; omega
  · simp only [h0, dite_false]
    by_cases hd : firstDescentFrom L 1 < L.length
    · have hdesc := firstDescentFrom_desc L 1 _ rfl hd
      simp only [hd, ite_true, LastAtomsCheck.notIncreasing.injEq]
      constructor
      · rintro rfl
        exact ⟨hd, hge, Nat.le_of_not_lt hdesc,
          fun k hk hk1 => hasc k (by omega) hk1 hk⟩
      · rintro ⟨hg, hg1, hle, hbefore⟩
        apply Nat.le_antisymm
        · apply Nat.le_of_not_lt
          intro hlt
          exact absurd (hasc g hg hg1 hlt) (Nat.not_lt_of_le hle)
        · apply Nat.le_of_not_lt
          intro hlt
          exact hdesc (hbefore _ hlt hge)
    · simp only [hd, ite_false]
      constructor
      · intro h
        revert h
        repeat' split
        all_goals simp
      · rintro ⟨hg, hg1, hle, _⟩
        exact absurd (hasc g hg hg1 (by omega)) (Nat.not_lt_of_le hle)

/-- **`check_last_atoms` reports `LastNotFinal` exactly for a non-empty, strictly increasing
list whose last entry is not `n − 1`** (or `n = 0`): the third check, after the first two
passed. -/
theorem checkLastAtomsModel_lastNotFinal_iff (L : List Nat) (n : Nat) :
    checkLastAtomsModel L n = .lastNotFinal ↔
      0 < L.length ∧ (∀ k (hk : k < L.length), 1 ≤ k → L[k - 1] < L[k]) ∧
        (n = 0 ∨ L.getLast? ≠ some (n - 1)) := by
  have spec := firstDescentFrom_spec L 1 (Nat.le_refl 1)
  unfold checkLastAtomsModel
  by_cases h0 : L.length = 0
  · simp [h0]
  · simp only [h0, dite_false]
    have hlast : L.getLast? = some L[L.length - 1] := by
      rw [List.getLast?_eq_getElem?, List.getElem?_eq_getElem]
    by_cases hd : firstDescentFrom L 1 < L.length
    · simp only [hd, ite_true, reduceCtorEq, false_iff]
      intro ⟨_, hinc, _⟩
      have := spec.mpr hinc
      omega
    · have hinc := spec.mp (by omega)
      simp only [hd, ite_false]
      constructor
      · intro h
        refine ⟨by omega, hinc, ?_⟩
        by_cases hn : n = 0
        · exact Or.inl hn
        · right
          rw [hlast]
          intro e
          simp only [Option.some.injEq] at e
          simp [hn, e] at h
      · rintro ⟨_, _, hn | hl⟩
        · simp [hn]
        · have : L[L.length - 1] ≠ n - 1 := fun e => hl (by rw [hlast, e])
          by_cases hn : n = 0 <;> simp [hn, this]

/-- **`energy_cell` is the row-major cell of two groups in range.** -/
theorem energyCellModel_eq (i j n : Nat) :
    energyCellModel i j n = if i < n ∧ j < n then some (energyCell n i j) else none := by
  unfold energyCellModel energyCell
  by_cases hi : i < n <;> by_cases hj : j < n <;> simp [hi, hj]

open NewenSpec.Thermo.Bath in
/-- **`thermo.bath.baths_of` is the same lookup**: the baths of an atom are those of the range
`groupOf` picks from the ranges' last atoms. -/
theorem bathsOf_eq_groupOf (rs : List Range) (a : Nat) :
    bathsOfModel rs a
      = (groupOf (rs.map (·.lastAtom)) a).bind fun g => rs[g]?.map fun r => (r.comBath, r.irBath) := by
  induction rs with
  | nil => simp [bathsOfModel, groupOf]
  | cons r rs ih =>
    unfold groupOf at ih ⊢
    rw [List.map_cons, List.findIdx?_cons]
    by_cases h : a ≤ r.lastAtom
    · simp [bathsOfModel, h]
    · simp [bathsOfModel, h, ih, Option.bind_map]
      rfl

end NewenSpec.Topology.Groups
