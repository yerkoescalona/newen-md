module

public import Mathlib.Analysis.Real.Sqrt
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse
public import Mathlib.Analysis.SpecialFunctions.Pow.Real
public import NewenTheory.Space.Periodicity

/-!
# What the floating-point models share

The models of the floating-point kernels (`NewenSpec/Thermostat.lean`, `NewenSpec/Nonbonded.lean`,
`NewenSpec/Bonded.lean`, `NewenSpec/Dynamics.lean`, `NewenSpec/Periodicity.lean`,
`NewenSpec/SoftCore.lean`, `NewenSpec/PerturbedBonded.lean`) are written once over any number type
`α` and used at `Float` (executed) and at `ℝ` (proven about). `+ − × ÷ <` come from the standard
classes; the square root, the three transcendental functions of the bonded kernels, the power
function, the absolute value and the two roundings to an integer have none that covers both, so
they are the classes below. A Lean `module`, as the theory files that import Mathlib.

`sqrt` is correctly rounded (IEEE 754); `acos`, `cos`, `sin` and `pow` are not, so `Float` and
Rust's `f64` agree on them only because both call the same C library (`Float.acos` is C `acos`,
Rust's `f64::acos` the platform `acos`, glibc's on Linux).
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

/-- The power function a kernel calls (Rust `f64::powf`): `Float.pow`, the C library's `pow`,
which IEEE 754 does not require to be correctly rounded, or `Real.rpow`. -/
class Powf (α : Type) where
  powf : α → α → α

instance : Powf Float := ⟨Float.pow⟩

noncomputable instance : Powf ℝ := ⟨fun x y => x ^ y⟩

/-- C `rint`: the nearest integer, halves to the even one (IEEE 754 roundToIntegralTiesToEven in
the default rounding mode). Exact, so every correct implementation returns the same bits, the
sign of zero included; Rust's `f64::round_ties_even` is one. Lean's own `Float.round` is C
`round` (halves away from zero, Rust's `f64::round`), a different function. Like `Float.floor`
it has no logical model: it runs only in compiled code (`lake exe spec-vectors`). -/
@[extern "rint"] opaque rintFloat : Float → Float

/-- Round to the nearest integer, halves to the even one: `rintFloat` (C `rint`, Rust
`f64::round_ties_even`) or `roundEven` over ℝ. -/
class Rint (α : Type) where
  rint : α → α

instance : Rint Float := ⟨rintFloat⟩

noncomputable instance : Rint ℝ := ⟨fun x => (NewenTheory.Space.Periodicity.roundEven x : ℝ)⟩

/-- The largest integer not above: `Float.floor` (C `floor`, Rust `f64::floor`, exact) or `⌊x⌋`
over ℝ. -/
class Floor (α : Type) where
  floor : α → α

instance : Floor Float := ⟨Float.floor⟩

noncomputable instance : Floor ℝ := ⟨fun x => ((⌊x⌋ : ℤ) : ℝ)⟩

/-- `x.abs()`: `Float.abs` (exact: clears the sign bit) or `|x|`. -/
class Fabs (α : Type) where
  fabs : α → α

instance : Fabs Float := ⟨Float.abs⟩

noncomputable instance : Fabs ℝ := ⟨fun x => |x|⟩

end NewenSpec

end
