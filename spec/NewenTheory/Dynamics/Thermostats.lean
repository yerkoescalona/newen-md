module

public import Mathlib.Analysis.Real.Sqrt
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Thermostats: Berendsen and Nosé–Hoover

A thermostat drives the kinetic energy `K` of a bath towards the value equipartition assigns to
the reference temperature `T₀`, `K₀ = ½ N_df k_B T₀`. Both thermostats here act once per step by
multiplying the velocities of the bath by a factor `λ`, so the temperature
`T = 2 K / (N_df k_B)` becomes `λ² T`.

**Berendsen weak coupling** (Berendsen, Postma, van Gunsteren, DiNola & Haak 1984,
J. Chem. Phys. 81, 3684) chooses `λ` so that the temperature relaxes towards `T₀` with time
constant `τ`, `dT/dt = (T₀ − T)/τ`, integrated over one step `dt`:
`λ² = 1 + dt/τ · (T₀/T − 1)`. It does not sample the canonical ensemble; it is a relaxation.

**Nosé–Hoover** (Nosé 1984, Mol. Phys. 52, 255; Hoover 1985, Phys. Rev. A 31, 1695) adds a
friction `ζ` with its own equation of motion, `dζ/dt = (T/T₀ − 1)/τ²`: a hot bath increases the
friction, a cold one decreases it, and `λ = 1 − ζ dt` is the velocity damping of one step.
**Nosé–Hoover chains** (Martyna, Klein & Tuckerman 1992, J. Chem. Phys. 97, 2635) thermostat
`ζ` itself with a second friction, that one with a third, and so on.

Everything here is over ℝ; how newen computes the same numbers in floating point is in
`NewenSpec/Thermostat.lean`. Lean's division by zero is zero, so a formula with `x / 0` is
defined but means nothing; the theorems carry the hypotheses that keep the divisions honest.

This file imports Mathlib (`Real.sqrt`, `field_simp`, `nlinarith`), and is therefore a Lean
`module` (Mathlib's import system); every definition is exposed to the files importing it.
-/

@[expose] public section

namespace NewenTheory.Dynamics.Thermostats

/-! ## Berendsen -/

/-- The Berendsen velocity scaling factor for temperature `T`, reference temperature `T₀`, time
step `dt` and coupling time `τ`. -/
noncomputable def berendsenScale (T T₀ dt τ : ℝ) : ℝ :=
  Real.sqrt (1 + dt / τ * (T₀ / T - 1))

/-- After one step the temperature has moved towards `T₀` by the fraction `dt/τ` of the gap:
`λ² T = T + dt/τ · (T₀ − T)`. Needs a non-negative radicand, which `0 < dt ≤ τ` and `T₀ ≥ 0`
guarantee (`berendsen_between`); `τ ≠ 0` is not needed for the identity. -/
theorem berendsen_relax {T T₀ dt τ : ℝ} (hT : 0 < T) (hr : 0 ≤ 1 + dt / τ * (T₀ / T - 1)) :
    berendsenScale T T₀ dt τ ^ 2 * T = T + dt / τ * (T₀ - T) := by
  unfold berendsenScale
  rw [Real.sq_sqrt hr]
  field_simp

/-- A bath already at the reference temperature is left alone. -/
theorem berendsen_fixed {T T₀ dt τ : ℝ} (hT : T ≠ 0) (h : T = T₀) :
    berendsenScale T T₀ dt τ = 1 := by
  subst h
  simp [berendsenScale, div_self hT]

/-- The deviation from `T₀` shrinks by the factor `1 − dt/τ` per step. -/
theorem berendsen_contraction {T T₀ dt τ : ℝ} (hT : 0 < T)
    (hr : 0 ≤ 1 + dt / τ * (T₀ / T - 1)) :
    berendsenScale T T₀ dt τ ^ 2 * T - T₀ = (1 - dt / τ) * (T - T₀) := by
  rw [berendsen_relax hT hr]
  ring

/-- With `0 < dt ≤ τ` the new temperature lies between the old one and `T₀`: Berendsen
coupling never overshoots. (With `dt = τ` it lands on `T₀`.) -/
theorem berendsen_between {T T₀ dt τ : ℝ} (hT : 0 < T) (hT₀ : 0 ≤ T₀) (hdt : 0 < dt)
    (hdτ : dt ≤ τ) :
    min T T₀ ≤ berendsenScale T T₀ dt τ ^ 2 * T ∧
      berendsenScale T T₀ dt τ ^ 2 * T ≤ max T T₀ := by
  have hτ : 0 < τ := lt_of_lt_of_le hdt hdτ
  have hf0 : 0 ≤ dt / τ := (div_pos hdt hτ).le
  have hf1 : dt / τ ≤ 1 := (div_le_one hτ).mpr hdτ
  have hq : 0 ≤ T₀ / T := div_nonneg hT₀ hT.le
  have hr : 0 ≤ 1 + dt / τ * (T₀ / T - 1) := by nlinarith
  rw [berendsen_relax hT hr]
  rcases le_total T T₀ with h | h
  · rw [min_eq_left h, max_eq_right h]
    constructor <;> nlinarith
  · rw [min_eq_right h, max_eq_left h]
    constructor <;> nlinarith

/-! ## Nosé–Hoover, one friction -/

/-- One Nosé–Hoover step: the friction `ζ` is advanced by `dt/τ² · (T/T₀ − 1)`, and the
velocity scaling factor of the step is `1 − ζ' dt` with the new friction `ζ'`. Returns
`(ζ', λ)`. -/
noncomputable def nhStep (ζ T T₀ dt τ : ℝ) : ℝ × ℝ :=
  let ζ' := ζ + dt / τ ^ 2 * (T / T₀ - 1)
  (ζ', 1 - ζ' * dt)

/-- At the reference temperature with no friction nothing happens. -/
theorem nh_fixed {ζ T T₀ dt τ : ℝ} (hT₀ : T₀ ≠ 0) (h : T = T₀ ∧ ζ = 0) :
    (nhStep ζ T T₀ dt τ).1 = 0 ∧ (nhStep ζ T T₀ dt τ).2 = 1 := by
  obtain ⟨rfl, rfl⟩ := h
  simp [nhStep, div_self hT₀]

/-- Negative feedback: the friction grows exactly when the bath is hotter than `T₀`. -/
theorem nh_feedback {ζ T T₀ dt τ : ℝ} (hdt : 0 < dt) (hτ : τ ≠ 0) (hT₀ : 0 < T₀) :
    T > T₀ ↔ (nhStep ζ T T₀ dt τ).1 > ζ := by
  have hc : 0 < dt / τ ^ 2 := div_pos hdt (by positivity)
  simp only [nhStep, gt_iff_lt]
  constructor
  · intro h
    have : 0 < T / T₀ - 1 := by rw [sub_pos, one_lt_div hT₀]; exact h
    nlinarith
  · intro h
    have h' : 0 < dt / τ ^ 2 * (T / T₀ - 1) := by linarith
    have : 0 < T / T₀ - 1 := (pos_iff_pos_of_mul_pos h').mp hc
    rwa [sub_pos, one_lt_div hT₀] at this

/-- The step cools (`λ < 1`) exactly when the new friction is positive. -/
theorem nh_scale_lt_one_iff {ζ T T₀ dt τ : ℝ} (hdt : 0 < dt) :
    (nhStep ζ T T₀ dt τ).2 < 1 ↔ (nhStep ζ T T₀ dt τ).1 > 0 := by
  simp only [nhStep, gt_iff_lt]
  constructor
  · intro h
    have : 0 < (ζ + dt / τ ^ 2 * (T / T₀ - 1)) * dt := by linarith
    exact (pos_iff_pos_of_mul_pos this).mpr hdt
  · intro h
    have := mul_pos h hdt
    linarith

/-! ## Nosé–Hoover chains

Link `0` is the friction on the particles; link `i ≥ 1` is a friction on link `i − 1`, with
"mass" `Q_0 = τ²` and `Q_i = τ²/N_df`. Over one step the links are advanced from the last one
down to link 0, each using the value its upper neighbour has *already* been given this step and
the *old* value of its lower neighbour (the order gromosXX uses):

```text
ζ_{n−1} += (Q_{n−2} ζ_{n−2}² − 1/N_df) / Q_{n−1} · dt
ζ_i     += ((Q_{i−1} ζ_{i−1}² − 1/N_df) / Q_i − ζ_i ζ_{i+1}) · dt      i = n−2 … 1
ζ_0     += ((T/T₀ − 1) / Q_0 − ζ_1 ζ_0) · dt
```

Only this transcription is given for chains, with one theorem (`nhc_single_link`). The fixed
point of a chain is not `ζ = 0`: link `i ≥ 1` is driven towards `Q_{i−1} ζ_{i−1}² = 1/N_df`,
the equipartition of the thermostat variable below it (its "kinetic energy" `½ Q ζ²` carries
`½ k_B T₀` per degree of freedom, in reduced units), so at equilibrium the lower links fluctuate
around a non-zero mean square. -/

/-- Links `1 …` of a chain, given the old value `zb` and mass `Qb` of the link below the first
one, and the mass `Q = τ²/N_df` of every upper link. -/
noncomputable def nhcUpper (N Q dt : ℝ) : ℝ → ℝ → List ℝ → List ℝ
  | _, _, [] => []
  | zb, Qb, [z] => [z + (Qb * zb ^ 2 - 1 / N) / Q * dt]
  | zb, Qb, z :: z' :: rest =>
    let up := nhcUpper N Q dt z Q (z' :: rest)
    (z + ((Qb * zb ^ 2 - 1 / N) / Q - z * up.headD 0) * dt) :: up

/-- One step of a Nosé–Hoover chain `ζ = [ζ_0, …]` with `N` degrees of freedom: the new chain
and the velocity scaling factor `1 − ζ_0' dt`. -/
noncomputable def nhcStep (ζ : List ℝ) (T T₀ dt τ N : ℝ) : List ℝ × ℝ :=
  match ζ with
  | [] => ([], 1)
  | z0 :: rest =>
    let up := nhcUpper N (τ ^ 2 / N) dt z0 (τ ^ 2) rest
    let z0' := z0 + ((T / T₀ - 1) / τ ^ 2 - up.headD 0 * z0) * dt
    (z0' :: up, 1 - z0' * dt)

/-- A chain of one link is plain Nosé–Hoover. (gromosXX never runs one: it asserts `nhc > 1`
and uses `nhStep` for one link.) -/
theorem nhc_single_link (ζ T T₀ dt τ N : ℝ) :
    nhcStep [ζ] T T₀ dt τ N = ([(nhStep ζ T T₀ dt τ).1], (nhStep ζ T T₀ dt τ).2) := by
  simp only [nhcStep, nhcUpper, nhStep, List.headD_nil]
  congr 1 <;> [congr 1; skip] <;> ring

end NewenTheory.Dynamics.Thermostats

end
