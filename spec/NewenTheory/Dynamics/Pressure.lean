module

public import Mathlib.Analysis.SpecialFunctions.Pow.Real
public import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Pressure and the Berendsen barostat

**The virial theorem** (Clausius 1870; GROMOS Vol. 2) gives the pressure of a system of particles
in volume `V` from its kinetic energy and its virial. GROMOS writes it as a tensor,

```text
P = 2/V · (K − Ξ),        K = ½ Σ m v ⊗ v,        Ξ = −½ Σ r ⊗ f,
```

and the scalar pressure is a third of its trace, `P = 2/(3V) · (tr K − tr Ξ)`
(`scalarPressure_virial`). Without forces it is the ideal gas law (`ideal_gas`).

**The Berendsen barostat** (Berendsen et al. 1984, J. Chem. Phys. 81, 3684) lets the pressure
relax towards `P₀` with time constant `τ`, `dP/dt = (P₀ − P)/τ`, by scaling every coordinate and
the box by `μ`, so the volume by `μ³`. For a system of isothermal compressibility
`β = −(1/V) dV/dP` that needs `μ³ = 1 − β dt/τ (P₀ − P)`:

* the reference pressure is a fixed point, `μ = 1` (`barostat_fixed`);
* the box grows exactly when the pressure is above the reference (`barostat_grows_iff`), and
  shrinks exactly when it is below (`barostat_shrinks_iff`);
* the volume is scaled by the radicand (`barostat_volume`), and if the pressure responds to the
  volume as the compressibility says, the new pressure has moved the fraction `dt/τ` of the way
  to `P₀`, as in the Berendsen thermostat (`barostat_relax`).

Everything here is over ℝ; newen's floating-point version is `NewenSpec/Dynamics.lean`.
-/

@[expose] public section

namespace NewenTheory.Dynamics.Pressure

open Matrix

/-! ## The virial theorem -/

/-- The pressure tensor `2/V · (K − Ξ)` of the kinetic-energy tensor `K` and the virial `Ξ`. -/
noncomputable def pressureTensor (K Ξ : Matrix (Fin 3) (Fin 3) ℝ) (V : ℝ) :
    Matrix (Fin 3) (Fin 3) ℝ :=
  (2 / V) • (K - Ξ)

/-- The scalar pressure, a third of the trace of the tensor. -/
noncomputable def scalarPressure (P : Matrix (Fin 3) (Fin 3) ℝ) : ℝ := P.trace / 3

/-- The scalar virial theorem: `P = 2/(3V) · (tr K − tr Ξ)`. -/
theorem scalarPressure_virial (K Ξ : Matrix (Fin 3) (Fin 3) ℝ) (V : ℝ) :
    scalarPressure (pressureTensor K Ξ V) = 2 / (3 * V) * (K.trace - Ξ.trace) := by
  simp only [scalarPressure, pressureTensor, trace_smul, trace_sub, smul_eq_mul]
  field_simp

/-- **The ideal gas.** No virial and the kinetic energy of `N` free particles at temperature `T`
(equipartition, `tr K = 3/2 · N k_B T`): `P V = N k_B T`. -/
theorem ideal_gas {K : Matrix (Fin 3) (Fin 3) ℝ} {V N kB T : ℝ} (hV : V ≠ 0)
    (hK : K.trace = 3 / 2 * N * kB * T) :
    scalarPressure (pressureTensor K 0 V) * V = N * kB * T := by
  rw [scalarPressure_virial, trace_zero, sub_zero, hK]
  field_simp

/-! ## The Berendsen barostat -/

/-- The Berendsen scaling factor of coordinates and box for pressure `P`, reference `P₀`,
compressibility `β`, coupling time `τ` and time step `dt`. -/
noncomputable def barostatScale (P P₀ β τ dt : ℝ) : ℝ :=
  (1 - β * dt / τ * (P₀ - P)) ^ (1 / 3 : ℝ)

/-- The volume is scaled by `μ³ = 1 − β dt/τ (P₀ − P)`, when that is not negative. -/
theorem barostat_volume {P P₀ β τ dt : ℝ} (hr : 0 ≤ 1 - β * dt / τ * (P₀ - P)) :
    barostatScale P P₀ β τ dt ^ 3 = 1 - β * dt / τ * (P₀ - P) := by
  unfold barostatScale
  rw [← Real.rpow_natCast, ← Real.rpow_mul hr]
  norm_num

/-- At the reference pressure nothing moves. -/
theorem barostat_fixed (P β τ dt : ℝ) : barostatScale P P β τ dt = 1 := by
  simp [barostatScale]

/-- With a positive coupling `β dt/τ` the box grows exactly when the pressure is above the
reference. -/
theorem barostat_grows_iff {P P₀ β τ dt : ℝ} (hc : 0 < β * dt / τ)
    (hr : 0 ≤ 1 - β * dt / τ * (P₀ - P)) :
    1 < barostatScale P P₀ β τ dt ↔ P₀ < P := by
  unfold barostatScale
  rw [show (1 < (1 - β * dt / τ * (P₀ - P)) ^ (1 / 3 : ℝ)) ↔
      (1 : ℝ) ^ (1 / 3 : ℝ) < (1 - β * dt / τ * (P₀ - P)) ^ (1 / 3 : ℝ) by rw [Real.one_rpow],
    Real.rpow_lt_rpow_iff zero_le_one hr (by norm_num)]
  constructor <;> intro h <;> nlinarith

/-- With a positive coupling the box shrinks exactly when the pressure is below the reference. -/
theorem barostat_shrinks_iff {P P₀ β τ dt : ℝ} (hc : 0 < β * dt / τ)
    (hr : 0 ≤ 1 - β * dt / τ * (P₀ - P)) :
    barostatScale P P₀ β τ dt < 1 ↔ P < P₀ := by
  unfold barostatScale
  rw [show ((1 - β * dt / τ * (P₀ - P)) ^ (1 / 3 : ℝ) < 1) ↔
      (1 - β * dt / τ * (P₀ - P)) ^ (1 / 3 : ℝ) < (1 : ℝ) ^ (1 / 3 : ℝ) by rw [Real.one_rpow],
    Real.rpow_lt_rpow_iff hr zero_le_one (by norm_num)]
  constructor <;> intro h <;> nlinarith

/-- The pressure at volume `V'` of a system at pressure `P` and volume `V` with compressibility
`β`, to first order: `P − (V' − V)/(β V)`. -/
noncomputable def linearResponse (P β V V' : ℝ) : ℝ := P - (V' - V) / (β * V)

/-- **Relaxation.** If the pressure responds to the volume as the compressibility says, one
barostat step moves it the fraction `dt/τ` of the way to `P₀`: `P' = P + dt/τ · (P₀ − P)`. -/
theorem barostat_relax {P P₀ β τ dt V : ℝ} (hβ : β ≠ 0) (hV : V ≠ 0)
    (hr : 0 ≤ 1 - β * dt / τ * (P₀ - P)) :
    linearResponse P β V (barostatScale P P₀ β τ dt ^ 3 * V) = P + dt / τ * (P₀ - P) := by
  rw [barostat_volume hr, linearResponse]
  field_simp
  ring

end NewenTheory.Dynamics.Pressure

end
