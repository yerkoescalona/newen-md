module

public import Mathlib.Analysis.Real.Sqrt

/-!
# What the floating-point models share

The models of the floating-point kernels (`NewenSpec/Thermostat.lean`, `NewenSpec/Nonbonded.lean`)
are written once over any number type `α` and used at `Float` (executed) and at `ℝ` (proven
about). `+ − × ÷ <` come from the standard classes; the square root has none that covers both, so
it is this class. A Lean `module`, as the theory files that import Mathlib.
-/

@[expose] public section

namespace NewenSpec

/-- The square root a kernel calls: `Float.sqrt` (C `sqrt`, correctly rounded) or `Real.sqrt`. -/
class Sqrt (α : Type) where
  sqrt : α → α

instance : Sqrt Float := ⟨Float.sqrt⟩

noncomputable instance : Sqrt ℝ := ⟨Real.sqrt⟩

end NewenSpec

end
