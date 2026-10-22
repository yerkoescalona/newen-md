/-!
# Kinds of nonbonded pairs

A force field does not treat every pair of atoms alike. Atoms one or two bonds apart (1-2, 1-3)
are *excluded*: the bond and angle terms already describe them, so they get no Lennard-Jones or
Coulomb term (in GROMOS only the reaction-field correction of excluded pairs). Atoms three bonds
apart (1-4) interact with special, softer 1-4 parameters. A pair listed in a perturbation
topology's `PERTATOMPAIR` block, or a pair of an EDS atom, is computed by that atom's own loop.
All three stay off the pairlist; every other pair is *normal*.

Each atom `i` carries three lists of partners: `E i` (excluded), `F i` (1-4) and `A i` (owned
pairs). The kind of a pair is decided by membership, with exclusion first: a pair that is both
excluded and 1-4 is excluded.

Core Lean only.
-/

namespace NewenTheory.Topology.Pairs

/-- The kind of a pair of atoms. -/
inductive PairKind where
  | normal
  | excluded
  | oneFour
  | atomPair
  deriving DecidableEq, Repr

/-- The kind of the pair `(i, j)` from atom `i`'s three partner lists, by priority. -/
def pairKind (E F A : List Nat) (j : Nat) : PairKind :=
  if j ∈ E then .excluded
  else if j ∈ F then .oneFour
  else if j ∈ A then .atomPair
  else .normal

/-- A partner relation stored per atom is symmetric: `j` is a partner of `i` exactly when `i` is
a partner of `j`. -/
def Symmetric (L : Nat → List Nat) : Prop := ∀ i j, j ∈ L i ↔ i ∈ L j

theorem pairKind_excluded_iff (E F A : List Nat) (j : Nat) :
    pairKind E F A j = .excluded ↔ j ∈ E := by
  unfold pairKind
  by_cases hE : j ∈ E <;> by_cases hF : j ∈ F <;> by_cases hA : j ∈ A <;> simp [hE, hF, hA]

theorem pairKind_oneFour_iff (E F A : List Nat) (j : Nat) :
    pairKind E F A j = .oneFour ↔ j ∉ E ∧ j ∈ F := by
  unfold pairKind
  by_cases hE : j ∈ E <;> by_cases hF : j ∈ F <;> by_cases hA : j ∈ A <;> simp [hE, hF, hA]

theorem pairKind_atomPair_iff (E F A : List Nat) (j : Nat) :
    pairKind E F A j = .atomPair ↔ j ∉ E ∧ j ∉ F ∧ j ∈ A := by
  unfold pairKind
  by_cases hE : j ∈ E <;> by_cases hF : j ∈ F <;> by_cases hA : j ∈ A <;> simp [hE, hF, hA]

/-- A pair is normal, i.e. on the pairlist, exactly when it is in none of the three lists. -/
theorem pairKind_normal_iff (E F A : List Nat) (j : Nat) :
    pairKind E F A j = .normal ↔ j ∉ E ∧ j ∉ F ∧ j ∉ A := by
  unfold pairKind
  by_cases hE : j ∈ E <;> by_cases hF : j ∈ F <;> by_cases hA : j ∈ A <;> simp [hE, hF, hA]

/-- **The kind of a pair does not depend on which atom is asked**, when each list is
symmetric: `(i, j)` and `(j, i)` are the same pair. -/
theorem pairKind_symm {E F A : Nat → List Nat} (hE : Symmetric E) (hF : Symmetric F)
    (hA : Symmetric A) (i j : Nat) :
    pairKind (E i) (F i) (A i) j = pairKind (E j) (F j) (A j) i := by
  simp only [pairKind, hE i j, hF i j, hA i j]

end NewenTheory.Topology.Pairs
