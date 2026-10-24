import NewenTheory.Dynamics.DegreesOfFreedom

/-!
# Model of `crates/newen-logic/src/thermo/dof.rs`

Rust: `crates/newen-logic/src/thermo/dof.rs`, `total_dof`
(spec key `thermo.dof.total_dof` in `spec/SPEC.toml`).
Theory: `NewenTheory/Dynamics/DegreesOfFreedom.lean` (`Molecule`, `internalDof`, `dof`).

This file has two parts:

1. the *model*: `totalDofModel`, the Rust expression over `Nat` and `Int`;
2. the *bridge*: `totalDofModel_eq_dof`, proving that what the caller passes in (the solute's
   atoms and constraint pairs, `k` copies of the solvent template) makes the model the theory's
   count for that system.

The Aeneas refinement (Rust = model) is `Refinement.total_dof_eq` in `NewenSpec/Refinement.lean`.

Only core Lean is used (`omega` and `Int` lemmas), no Mathlib.
-/

namespace NewenSpec.Thermo.Dof

open NewenTheory.Dynamics.DegreesOfFreedom

/-! ## 1. Model (mirrors the Rust function) -/

/-- `total_dof`, written as the Rust code is: `3 · n_atoms − solvent_constraints −
solute_constraints − ndfmin`, in signed arithmetic. -/
def totalDofModel (nAtoms solventConstraints soluteConstraints : Nat) (ndfmin : Int) : Int :=
  3 * (nAtoms : Int) - solventConstraints - soluteConstraints - ndfmin

/-! ## 2. Bridge: the model is the theory's count -/

/-- For any system, the model given its atom count and its constraints split into two parts is
the theory's count. -/
theorem totalDofModel_eq_dof_of_split (ms : List Molecule) (c₁ c₂ : Nat)
    (h : c₁ + c₂ = constraints ms) (ndfmin : Int) :
    totalDofModel (atoms ms) c₁ c₂ ndfmin = dof ms ndfmin := by
  rw [dof, internalDof_eq, ← h]
  simp only [totalDofModel, coordinates]
  omega

/-- **The model is what the caller means.** `newen-run`'s `total_dof` passes the number of atoms,
`k · (constraints of the solvent template)` for `k` solvent molecules, and the number of solute
constraint pairs. For a system of solute molecules followed by `k` identical solvent molecules
that is exactly the theory's degrees of freedom of the whole system. -/
theorem totalDofModel_eq_dof (solute : List Molecule) (solvent : Molecule) (k : Nat)
    (ndfmin : Int) :
    totalDofModel (atoms solute + k * solvent.atoms) (k * solvent.constraints)
        (constraints solute) ndfmin =
      dof (solute ++ List.replicate k solvent) ndfmin := by
  have ha : atoms solute + k * solvent.atoms = atoms (solute ++ List.replicate k solvent) := by
    rw [atoms_append, atoms_replicate]
  rw [ha]
  apply totalDofModel_eq_dof_of_split
  rw [constraints_append, constraints_replicate, Nat.add_comm]

end NewenSpec.Thermo.Dof
