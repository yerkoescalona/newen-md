/-!
# Groups of consecutive atoms

GROMOS partitions the atoms into consecutive groups by listing the last atom of each group: the
energy groups of the FORCE block (`NRE`), the temperature-bath ranges of `MULTIBATH` (DOFSET), the
solute molecules, the temperature and pressure groups. Group `g` holds the atoms after group
`g − 1`'s last atom up to its own, so an atom belongs to the *first* group whose last atom
reaches it (`groupOf`). The list describes a partition of `n` atoms into non-empty groups when it
is non-empty, strictly increasing, and ends at atom `n − 1` (`LastAtomsOk`, the rule the FORCE
block's `NRE` must obey).

Energies are booked per pair of energy groups, in an `n × n` matrix stored row by row
(`energyCell`).

`thermo.bath.baths_of` (`NewenSpec/Thermo/Bath.lean`, `bathsOfModel`) is the same rule on the
DOFSET ranges: `NewenSpec.Topology.Groups.bathsOf_eq_groupOf` states it.

Indices are 0-based. Core Lean only.
-/

namespace NewenTheory.Topology.Groups

/-- The group of `atom`: the least `g` with `atom ≤ L[g]`, none if no group reaches it. -/
def groupOf (L : List Nat) (atom : Nat) : Option Nat :=
  L.findIdx? fun l => decide (atom ≤ l)

/-- `groupOf` is the least index whose last atom reaches the atom. -/
theorem groupOf_eq_some_iff (L : List Nat) (a g : Nat) :
    groupOf L a = some g ↔ ∃ h : g < L.length, a ≤ L[g] ∧ ∀ k (_ : k < g), L[k] < a := by
  unfold groupOf
  rw [List.findIdx?_eq_some_iff_getElem]
  simp [Nat.not_le]

/-- The list of last atoms partitions `n` atoms into non-empty consecutive groups: at least one
group, strictly increasing, the last group ending at atom `n − 1`. -/
structure LastAtomsOk (L : List Nat) (n : Nat) : Prop where
  nonempty : 0 < L.length
  increasing : ∀ g (_ : g + 1 < L.length), L[g] < L[g + 1]
  final : L.getLast? = some (n - 1)
  atoms : 0 < n

namespace LastAtomsOk

variable {L : List Nat} {n : Nat}

theorem lt (h : LastAtomsOk L n) {a b : Nat} (hab : a < b) (hb : b < L.length) : L[a] < L[b] := by
  induction b with
  | zero => omega
  | succ b ih =>
    have hs := h.increasing b hb
    by_cases e : a = b
    · subst e; exact hs
    · exact Nat.lt_trans (ih (by omega) (by omega)) hs

theorem le (h : LastAtomsOk L n) {a b : Nat} (hab : a ≤ b) (hb : b < L.length) : L[a] ≤ L[b] := by
  by_cases e : a = b
  · subst e; exact Nat.le_refl _
  · exact Nat.le_of_lt (h.lt (by omega) hb)

theorem last_eq (h : LastAtomsOk L n) : L[L.length - 1]'(by have := h.nonempty; omega) = n - 1 := by
  have e := h.final
  rw [List.getLast?_eq_getElem?, List.getElem?_eq_getElem (by have := h.nonempty; omega)] at e
  exact Option.some.inj e

/-- Every group ends at an atom below `n`. -/
theorem lt_atoms (h : LastAtomsOk L n) (g : Nat) (hg : g < L.length) : L[g] < n := by
  have := h.le (show g ≤ L.length - 1 by omega) (by omega)
  have := h.last_eq
  have := h.atoms
  omega

end LastAtomsOk

/-- **Every atom has a group, and nothing beyond the atoms has one.** -/
theorem groupOf_some_iff {L : List Nat} {n : Nat} (h : LastAtomsOk L n) (a : Nat) :
    (groupOf L a).isSome ↔ a < n := by
  constructor
  · intro hs
    cases e : groupOf L a with
    | none => simp [e] at hs
    | some g =>
      obtain ⟨hg, hag, _⟩ := (groupOf_eq_some_iff L a g).mp e
      exact Nat.lt_of_le_of_lt hag (h.lt_atoms g hg)
  · intro ha
    have hl := h.nonempty
    have hlast := h.last_eq
    unfold groupOf
    rw [List.findIdx?_isSome, List.any_eq_true]
    exact ⟨L[L.length - 1], List.getElem_mem _, by simp; omega⟩

/-- **Groups are contiguous**: the atoms of group `g` are exactly those after group `g − 1`'s
last atom up to group `g`'s. -/
theorem groupOf_contiguous {L : List Nat} {n : Nat} (h : LastAtomsOk L n) {g : Nat}
    (hg : g < L.length) (a : Nat) :
    groupOf L a = some g ↔ (∀ hpos : 0 < g, L[g - 1] < a) ∧ a ≤ L[g] := by
  rw [groupOf_eq_some_iff]
  constructor
  · rintro ⟨_, hag, hk⟩
    exact ⟨fun hpos => hk (g - 1) (by omega), hag⟩
  · rintro ⟨hlo, hag⟩
    refine ⟨hg, hag, fun k hk => ?_⟩
    have hle := h.le (show k ≤ g - 1 by omega) (by omega)
    have := hlo (by omega)
    omega

/-- **Later atoms are in the same or a later group** (for any list). -/
theorem groupOf_mono {L : List Nat} {a b ga gb : Nat} (hab : a ≤ b)
    (ha : groupOf L a = some ga) (hb : groupOf L b = some gb) : ga ≤ gb := by
  obtain ⟨_, _, hka⟩ := (groupOf_eq_some_iff L a ga).mp ha
  obtain ⟨hgb, hbl, _⟩ := (groupOf_eq_some_iff L b gb).mp hb
  refine Nat.le_of_not_lt fun hlt => ?_
  have := hka gb hlt
  omega

/-- The row-major cell `[i][j]` of an `n × n` energy matrix. -/
def energyCell (n i j : Nat) : Nat := i * n + j

theorem energyCell_lt {n i j : Nat} (hi : i < n) (hj : j < n) : energyCell n i j < n * n := by
  unfold energyCell
  have : (i + 1) * n ≤ n * n := Nat.mul_le_mul_right n hi
  rw [Nat.succ_mul] at this
  omega

/-- **Distinct group pairs book into distinct cells.** -/
theorem energyCell_injective {n i j i' j' : Nat} (hj : j < n) (hj' : j' < n)
    (e : energyCell n i j = energyCell n i' j') : i = i' ∧ j = j' := by
  unfold energyCell at e
  have key : ∀ {x y u v : Nat}, u < n → x < y → x * n + u < y * n + v := by
    intro x y u v hu hxy
    have : (x + 1) * n ≤ y * n := Nat.mul_le_mul_right n hxy
    rw [Nat.succ_mul] at this
    omega
  rcases Nat.lt_trichotomy i i' with hlt | heq | hgt
  · exact absurd e (Nat.ne_of_lt (key hj hlt))
  · subst heq; exact ⟨rfl, by omega⟩
  · exact absurd e.symm (Nat.ne_of_lt (key hj' hgt))

end NewenTheory.Topology.Groups
