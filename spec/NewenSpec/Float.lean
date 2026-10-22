module

public import Mathlib.Analysis.Real.Sqrt
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse

/-!
# What the floating-point models share

The models of the floating-point kernels (`NewenSpec/Thermostat.lean`, `NewenSpec/Nonbonded.lean`,
`NewenSpec/Bonded.lean`) are written once over any number type `α` and used at `Float` (executed)
and at `ℝ` (proven about). `+ − × ÷ <` come from the standard classes; the square root and the
three transcendental functions of the bonded kernels have none that covers both, so they are these
classes. A Lean `module`, as the theory files that import Mathlib.

`sqrt` is correctly rounded (IEEE 754); `acos`, `cos` and `sin` are not, so `Float` and Rust's
`f64` agree on them only because both call the same C library (`Float.acos` is C `acos`, Rust's
`f64::acos` the platform `acos`, glibc's on Linux).
-/

@[expose] public section

namespace NewenSpec

/-- The square root a kernel calls: `Float.sqrt` (C `sqrt`, correctly rounded) or `Real.sqrt`. -/
class Sqrt (α : Type) where
  sqrt : α → α

instance : Sqrt Float := ⟨Float.sqrt⟩

noncomputable instance : Sqrt ℝ := ⟨Real.sqrt⟩

/-- The arc cosine a kernel calls: `Float.acos` (C `acos`) or `Real.arccos`. -/
class Acos (α : Type) where
  acos : α → α

/-- The cosine a kernel calls: `Float.cos` (C `cos`) or `Real.cos`. -/
class Cos (α : Type) where
  cos : α → α

/-- The sine a kernel calls: `Float.sin` (C `sin`) or `Real.sin`. -/
class Sin (α : Type) where
  sin : α → α

instance : Acos Float := ⟨Float.acos⟩
instance : Cos Float := ⟨Float.cos⟩
instance : Sin Float := ⟨Float.sin⟩

noncomputable instance : Acos ℝ := ⟨Real.arccos⟩
noncomputable instance : Cos ℝ := ⟨Real.cos⟩
noncomputable instance : Sin ℝ := ⟨Real.sin⟩

end NewenSpec

end
