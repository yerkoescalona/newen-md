/-!
# Kinetic degrees of freedom

A system of `N` point particles in three dimensions has `3 N` coordinates. Every independent
holonomic constraint (a fixed bond length, a rigid water geometry) removes one of them. Overall
motion that the simulation does not let the particles exchange with a heat bath removes more: a
fixed or removed centre-of-mass translation removes 3, and GROMOS takes the number to remove as
`NDFMIN` (BOUNDCOND block). What is left is the number of kinetic degrees of freedom `N_df`.

This number matters because of equipartition: in equilibrium every quadratic kinetic degree of
freedom carries `½ k_B T` on average, so `⟨K⟩ = ½ N_df k_B T`. The instantaneous temperature
`T = 2 K / (N_df k_B)` and every thermostat that drives `K` towards a target therefore depend on
counting `N_df` correctly. Equipartition itself is not proved here; only the counting is.

A system is described molecule by molecule: each `Molecule` has its atoms and its constraints.
The count is the sum over molecules of `3 · atoms − constraints`, minus `NDFMIN` for the whole
system. It is an integer (`Int`) on purpose: more constraints than coordinates is not a physical
system, and the count then becomes negative rather than meaningless.

Only core Lean is used (`omega` and `Int` lemmas), no Mathlib.
-/

namespace NewenTheory.Dynamics.DegreesOfFreedom

/-- The number of coordinates of `n` point particles in three dimensions. -/
def coordinates (n : Nat) : Nat := 3 * n

/-- A molecule as far as counting goes: its atoms and the independent holonomic constraints
imposed on them. -/
structure Molecule where
  atoms : Nat
  constraints : Nat

/-- The degrees of freedom one molecule keeps: its coordinates less its constraints. -/
def Molecule.dof (m : Molecule) : Int := (coordinates m.atoms : Int) - m.constraints

/-- All atoms of a collection of molecules. -/
def atoms : List Molecule → Nat
  | [] => 0
  | m :: ms => m.atoms + atoms ms

/-- All constraints of a collection of molecules. -/
def constraints : List Molecule → Nat
  | [] => 0
  | m :: ms => m.constraints + constraints ms

/-- The sum over molecules of the degrees of freedom each keeps. -/
def internalDof : List Molecule → Int
  | [] => 0
  | m :: ms => m.dof + internalDof ms

/-- The kinetic degrees of freedom of a system of molecules from which `ndfmin` degrees of
freedom of overall motion are removed. Physically `ndfmin ≥ 0`; it is an `Int` so that the count
is defined for every input. -/
def dof (ms : List Molecule) (ndfmin : Int) : Int := internalDof ms - ndfmin

/-- The molecule-by-molecule sum is the global count: `3 N` minus all constraints. -/
theorem internalDof_eq (ms : List Molecule) :
    internalDof ms = (coordinates (atoms ms) : Int) - constraints ms := by
  induction ms with
  | nil => simp [internalDof, atoms, constraints, coordinates]
  | cons m ms ih =>
    simp only [internalDof, atoms, constraints, Molecule.dof, coordinates, ih]
    omega

/-- **Additivity**: the degrees of freedom of two disjoint subsystems add up. -/
theorem internalDof_append (ms₁ ms₂ : List Molecule) :
    internalDof (ms₁ ++ ms₂) = internalDof ms₁ + internalDof ms₂ := by
  induction ms₁ with
  | nil => simp [internalDof]
  | cons m ms ih =>
    simp only [List.cons_append, internalDof, ih]
    omega

/-- Additivity with the overall motion removed once, from the whole system. -/
theorem dof_append (ms₁ ms₂ : List Molecule) (ndfmin : Int) :
    dof (ms₁ ++ ms₂) ndfmin = internalDof ms₁ + internalDof ms₂ - ndfmin := by
  rw [dof, internalDof_append]

/-- **Identical molecules**: `k` copies of one molecule keep `k` times its degrees of freedom.
This is how a solvent stored as one template molecule repeated `k` times is counted. -/
theorem internalDof_replicate (k : Nat) (m : Molecule) :
    internalDof (List.replicate k m) = (k : Int) * m.dof := by
  induction k with
  | zero => simp [internalDof]
  | succ k ih =>
    rw [List.replicate_succ, internalDof, ih, Int.natCast_succ, Int.add_mul, Int.one_mul,
      Int.add_comm]

/-- The atoms of `k` identical molecules. -/
theorem atoms_replicate (k : Nat) (m : Molecule) :
    atoms (List.replicate k m) = k * m.atoms := by
  induction k with
  | zero => simp [atoms]
  | succ k ih => rw [List.replicate_succ, atoms, ih, Nat.succ_mul, Nat.add_comm]

/-- The constraints of `k` identical molecules. -/
theorem constraints_replicate (k : Nat) (m : Molecule) :
    constraints (List.replicate k m) = k * m.constraints := by
  induction k with
  | zero => simp [constraints]
  | succ k ih => rw [List.replicate_succ, constraints, ih, Nat.succ_mul, Nat.add_comm]

/-- The atoms of two subsystems add up. -/
theorem atoms_append (ms₁ ms₂ : List Molecule) : atoms (ms₁ ++ ms₂) = atoms ms₁ + atoms ms₂ := by
  induction ms₁ with
  | nil => simp [atoms]
  | cons m ms ih => simp only [List.cons_append, atoms, ih]; omega

/-- The constraints of two subsystems add up. -/
theorem constraints_append (ms₁ ms₂ : List Molecule) :
    constraints (ms₁ ++ ms₂) = constraints ms₁ + constraints ms₂ := by
  induction ms₁ with
  | nil => simp [constraints]
  | cons m ms ih => simp only [List.cons_append, constraints, ih]; omega

/-- **Physicality**: the count is non-negative exactly when the constraints and the removed
overall motion together do not exceed the `3 N` coordinates. -/
theorem dof_nonneg_iff (ms : List Molecule) (ndfmin : Int) :
    0 ≤ dof ms ndfmin ↔ (constraints ms : Int) + ndfmin ≤ coordinates (atoms ms) := by
  rw [dof, internalDof_eq]
  omega

end NewenTheory.Dynamics.DegreesOfFreedom
