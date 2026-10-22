import NewenTheory.Topology.Pairs

/-!
# Model of `crates/newen-logic/src/topology/pairs.rs`

Rust: `crates/newen-logic/src/topology/pairs.rs`: `contains_sorted`, `pair_kind` (spec keys
`topology.pairs.*` in `spec/SPEC.toml`).
Theory: `NewenTheory/Topology/Pairs.lean` (`pairKind`, by membership).

1. The *models*, written as the Rust is: an empty list or a last entry below `j` answers
   `false` at once (gromosXX's early-out), otherwise the list is scanned.
2. The *bridges*: for ascending lists the early-out loses nothing, so the model is membership
   (`containsSortedModel_iff`) and the pair kind is the theory's (`pairKindModel_eq`).

The Aeneas refinements (Rust = model) are in `NewenSpec/Refinement.lean`. Core Lean only.
-/

namespace NewenSpec.Topology.Pairs

open NewenTheory.Topology.Pairs

/-! ## 1. Models (mirror the Rust functions) -/

/-- The loop of `contains_sorted`: is some entry equal to `j`? -/
def scanModel : List Nat → Nat → Bool
  | [], _ => false
  | x :: xs, j => if x = j then true else scanModel xs j

/-- `contains_sorted`: `false` for an empty list or a last entry below `j`, else the scan. -/
def containsSortedModel (list : List Nat) (j : Nat) : Bool :=
  if h : list.length = 0 then false
  else if list[list.length - 1] < j then false
  else scanModel list j

/-- `pair_kind`: excluded, else 1-4, else owned pair, else normal. -/
def pairKindModel (E F A : List Nat) (j : Nat) : PairKind :=
  if containsSortedModel E j then .excluded
  else if containsSortedModel F j then .oneFour
  else if containsSortedModel A j then .atomPair
  else .normal

/-! ## 2. Bridges: the models are membership -/

theorem scanModel_iff (list : List Nat) (j : Nat) : scanModel list j = true ↔ j ∈ list := by
  induction list with
  | nil => simp [scanModel]
  | cons x xs ih =>
    by_cases e : x = j
    · subst e; simp [scanModel]
    · simp only [scanModel, e, ite_false, ih, List.mem_cons]
      constructor
      · exact Or.inr
      · rintro (h | h)
        · exact absurd h.symm e
        · exact h

/-- In an ascending list every entry is at most the last. -/
theorem le_last {list : List Nat} (h : list.Pairwise (· ≤ ·)) {x : Nat} (hx : x ∈ list) :
    x ≤ list[list.length - 1]'(by
      have := List.length_pos_of_mem hx; omega) := by
  obtain ⟨i, hi, rfl⟩ := List.mem_iff_getElem.mp hx
  by_cases e : i = list.length - 1
  · subst e; exact Nat.le_refl _
  · exact List.pairwise_iff_getElem.mp h i (list.length - 1) hi (by omega) (by omega)

/-- **For an ascending list, `contains_sorted` is membership**: the early-out never answers
`false` for an entry of the list. -/
theorem containsSortedModel_iff {list : List Nat} (h : list.Pairwise (· ≤ ·)) (j : Nat) :
    containsSortedModel list j = true ↔ j ∈ list := by
  unfold containsSortedModel
  by_cases h0 : list.length = 0
  · simp only [h0, dite_true, Bool.false_eq_true, false_iff]
    intro hj
    have := List.length_pos_of_mem hj
    omega
  · simp only [h0, dite_false]
    by_cases hl : list[list.length - 1] < j
    · simp only [hl, ite_true, Bool.false_eq_true, false_iff]
      intro hj
      have := le_last h hj
      omega
    · simp only [hl, ite_false]
      exact scanModel_iff list j

/-- **For ascending lists, `pair_kind` is the theory's pair kind.** -/
theorem pairKindModel_eq {E F A : List Nat} (hE : E.Pairwise (· ≤ ·)) (hF : F.Pairwise (· ≤ ·))
    (hA : A.Pairwise (· ≤ ·)) (j : Nat) : pairKindModel E F A j = pairKind E F A j := by
  have e1 := containsSortedModel_iff hE j
  have e2 := containsSortedModel_iff hF j
  have e3 := containsSortedModel_iff hA j
  unfold pairKindModel pairKind
  by_cases m1 : j ∈ E <;> by_cases m2 : j ∈ F <;> by_cases m3 : j ∈ A <;>
    simp_all

end NewenSpec.Topology.Pairs
