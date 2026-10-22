module

public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Inv
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# The Lennard-Jones pair potential

Two atoms at distance `r` interact through `E(r) = C12/r¹² − C6/r⁶` (GROMOS writes the
Lennard-Jones potential with the two coefficients `C12 = 4εσ¹²` and `C6 = 4εσ⁶`, GROMOS Vol. 2).
The `r⁻¹²` term is the repulsion of overlapping electron clouds, the `r⁻⁶` term the dispersion
attraction.

The force on atom `i` is `−∇ᵢ E = −E'(r) · r_vec / r`, with `r_vec = x_i − x_j`. newen and gromosXX
compute the scalar `f = −E'(r)/r` once per pair and multiply it by `r_vec`, so the "force" of a
pair below is that scalar, `ljForce`, and `lj_force_is_minus_grad` says that it is the
derivative it should be.

Everything here is over ℝ; how newen computes it in floating point is in
`NewenSpec/Nonbonded.lean`. A Lean `module` (Mathlib's import system), as
`NewenTheory/Dynamics/Thermostats.lean`.
-/

@[expose] public section

namespace NewenTheory.Interactions.LennardJones

/-- The Lennard-Jones energy of a pair at distance `r`. -/
noncomputable def ljEnergy (c6 c12 r : ℝ) : ℝ := c12 / r ^ 12 - c6 / r ^ 6

/-- The scalar Lennard-Jones force of a pair at distance `r`, `−E'(r)/r`: the force on atom `i` is
this times `r_vec = x_i − x_j`. -/
noncomputable def ljForce (c6 c12 r : ℝ) : ℝ := 12 * c12 / r ^ 14 - 6 * c6 / r ^ 8

/-- **Force = −dE/dr.** The derivative of the energy at `r ≠ 0` is `−ljForce · r`, that is,
`ljForce = −E'(r)/r`. -/
theorem lj_force_is_minus_grad (c6 c12 : ℝ) {r : ℝ} (hr : r ≠ 0) :
    HasDerivAt (ljEnergy c6 c12) (-(ljForce c6 c12 r) * r) r := by
  have h : HasDerivAt (fun x => c12 * (x ^ 12)⁻¹ - c6 * (x ^ 6)⁻¹)
      (c12 * (-(↑12 * r ^ (12 - 1)) / (r ^ 12) ^ 2) - c6 * (-(↑6 * r ^ (6 - 1)) / (r ^ 6) ^ 2)) r :=
    HasDerivAt.sub
      (HasDerivAt.const_mul c12 (HasDerivAt.inv (hasDerivAt_pow 12 r) (pow_ne_zero 12 hr)))
      (HasDerivAt.const_mul c6 (HasDerivAt.inv (hasDerivAt_pow 6 r) (pow_ne_zero 6 hr)))
  have hE : ljEnergy c6 c12 = fun x => c12 * (x ^ 12)⁻¹ - c6 * (x ^ 6)⁻¹ := by
    funext x
    simp [ljEnergy, div_eq_mul_inv]
  rw [hE]
  refine h.congr_deriv ?_
  simp only [ljForce]
  field_simp
  ring

/-- The energy is zero at `r = σ`, where `σ⁶ = C12/C6`. -/
theorem lj_zero_at_sigma {c6 c12 r : ℝ} (hc12 : 0 < c12)
    (h : r ^ 6 = c12 / c6) : ljEnergy c6 c12 r = 0 := by
  have h12 : r ^ 12 = (c12 / c6) ^ 2 := by rw [← h]; ring
  rw [ljEnergy, h12, h]
  field_simp
  ring

/-- The minimum: at `r⁶ = 2·C12/C6` (`r = 2^(1/6) σ`) the force vanishes, so does the
derivative of the energy, and the energy is the well depth `−C6²/(4·C12)` (`−ε`). -/
theorem lj_minimum {c6 c12 r : ℝ} (hc6 : 0 < c6) (hc12 : 0 < c12) (hr : 0 < r)
    (h : r ^ 6 = 2 * c12 / c6) :
    ljForce c6 c12 r = 0 ∧ HasDerivAt (ljEnergy c6 c12) 0 r ∧
      ljEnergy c6 c12 r = -(c6 ^ 2 / (4 * c12)) := by
  have hc6' := hc6.ne'
  have hc12' := hc12.ne'
  have h8 : r ^ 8 = r ^ 2 * (2 * c12 / c6) := by rw [← h]; ring
  have h14 : r ^ 14 = r ^ 2 * (2 * c12 / c6) ^ 2 := by rw [← h]; ring
  have h12 : r ^ 12 = (2 * c12 / c6) ^ 2 := by rw [← h]; ring
  have hr2 : r ^ 2 ≠ 0 := pow_ne_zero 2 hr.ne'
  have hF : ljForce c6 c12 r = 0 := by
    rw [ljForce, h8, h14]
    field_simp
    ring
  refine ⟨hF, ?_, ?_⟩
  · simpa [hF] using lj_force_is_minus_grad c6 c12 hr.ne'
  · rw [ljEnergy, h12, h]
    field_simp
    ring

end NewenTheory.Interactions.LennardJones

end
