module

public import Mathlib.Analysis.Calculus.Deriv.Add
public import Mathlib.Analysis.Calculus.Deriv.Mul
public import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Module
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# The leap-frog integrator, kinetic energy and temperature

**Leap-frog** (Hockney 1970; GROMOS Vol. 2) advances positions at whole steps and velocities at
half steps. Its state at time `t` is the position `x(t)` and the velocity `v(t − dt/2)` of the
half step before; with the acceleration `a(x) = f(x)/m`, one step is a kick and a drift:

```text
v(t + dt/2) = v(t − dt/2) + a(x(t)) dt
x(t + dt)   = x(t) + v(t + dt/2) dt
```

Three properties make it the standard integrator of molecular dynamics, and they are proven
below: it is **time-reversible** (`leapFrog_reversible`), it is **exact for a constant
acceleration** (`leapFrog_exact_constant`, `leapFrog_iterate_exact`), and in one dimension its
**Jacobian determinant is 1** (`leapFrog_jacobian`): it conserves phase-space volume, the property
behind its good long-time energy behaviour.

**Kinetic energy.** With velocities known only at half steps, GROMOS takes the kinetic energy at
time `t` as the average of the two half-step energies, `½ m (|v(t−dt/2)|² + |v(t+dt/2)|²)/2`
(`halfStepKinetic`). It is not the energy of the averaged velocity, the velocity-Verlet
`v(t)`: it exceeds it by exactly `m dt² |a|²/8` (`halfStep_vs_verlet`), a second-order term.

**Centre of mass.** For a group of atoms the kinetic energy splits into that of the centre of mass
and the rest (König's theorem, `koenig`). GROMOS books the two parts to different baths; with the
half-step average the rest is never negative (`group_internal_nonneg`).

**Temperature.** Equipartition gives every degree of freedom `½ k_B T` on average, so `N_df`
degrees of freedom hold `K = ½ N_df k_B T`; the temperature of a kinetic energy is the inverse,
`T = 2K/(N_df k_B)` (`temperatureOf_iff`).

Everything here is over ℝ (positions and velocities in any real vector space, norms in any real
inner product space); newen's floating-point version is `NewenSpec/Dynamics.lean`. A Lean `module`,
as `NewenTheory/Dynamics/Thermostats.lean`.
-/

@[expose] public section

namespace NewenTheory.Dynamics.Integration

/-- Three-dimensional Euclidean space, where an atom's position and velocity live. -/
abbrev R3 := EuclideanSpace ℝ (Fin 3)

/-! ## Leap-frog -/

section LeapFrog

variable {E : Type*} [AddCommGroup E] [Module ℝ E]

/-- One leap-frog step in the acceleration field `a`: `(x(t), v(t − dt/2)) ↦ (x(t + dt),
v(t + dt/2))`, the kick `v' = v + a(x) dt`, then the drift `x' = x + v' dt`. -/
def leapFrog (a : E → E) (dt : ℝ) (s : E × E) : E × E :=
  let v' := s.2 + dt • a s.1
  (s.1 + dt • v', v')

/-- Reversing time at `(x(t), v(t − dt/2))`: the half step before `t` on the reversed trajectory
is the half step after `t` on this one, with the opposite sign, `−v(t + dt/2) = −(v + a(x) dt)`. -/
def reverse (a : E → E) (dt : ℝ) (s : E × E) : E × E :=
  (s.1, -(s.2 + dt • a s.1))

/-- Reversing twice is the identity. -/
theorem reverse_reverse (a : E → E) (dt : ℝ) (s : E × E) : reverse a dt (reverse a dt s) = s := by
  obtain ⟨x, v⟩ := s
  simp only [reverse, Prod.mk.injEq, true_and]
  module

/-- **Time reversibility.** A step, a reversal and a step again land on the starting position
with the velocity reversed: the leap-frog trajectory run backwards is the leap-frog trajectory. -/
theorem leapFrog_reversible (a : E → E) (dt : ℝ) (s : E × E) :
    leapFrog a dt (reverse a dt (leapFrog a dt s)) = reverse a dt s := by
  obtain ⟨x, v⟩ := s
  have hx : x + dt • (v + dt • a x) + dt • (-(v + dt • a x + dt • a (x + dt • (v + dt • a x))) +
      dt • a (x + dt • (v + dt • a x))) = x := by module
  simp only [leapFrog, reverse, Prod.mk.injEq]
  exact ⟨hx, by module⟩

/-- The same statement as an involution: step, then reverse, done twice, is the identity. -/
theorem reverse_comp_leapFrog_involutive (a : E → E) (dt : ℝ) :
    Function.Involutive (reverse a dt ∘ leapFrog a dt) := by
  intro s
  simp only [Function.comp_apply, leapFrog_reversible, reverse_reverse]

/-- **Exact for constant acceleration.** On the trajectory `x(t) = x₀ + t u + t²/2 a₀` with the
velocity `x'(t) = u + t a₀` taken at the half step before, one step lands on the trajectory again,
a whole step later. -/
theorem leapFrog_exact_constant (x₀ u a₀ : E) (t dt : ℝ) :
    leapFrog (fun _ => a₀) dt (x₀ + t • u + (t ^ 2 / 2) • a₀, u + (t - dt / 2) • a₀) =
      (x₀ + (t + dt) • u + ((t + dt) ^ 2 / 2) • a₀, u + (t + dt / 2) • a₀) := by
  simp only [leapFrog, Prod.mk.injEq]
  constructor <;> module

/-- `n` steps under constant acceleration are exact: they land on `x(t + n dt)`. -/
theorem leapFrog_iterate_exact (x₀ u a₀ : E) (dt : ℝ) (n : ℕ) :
    (leapFrog (fun _ => a₀) dt)^[n] (x₀, u - (dt / 2) • a₀) =
      (x₀ + (n * dt) • u + ((n * dt) ^ 2 / 2) • a₀, u + (n * dt - dt / 2) • a₀) := by
  induction n with
  | zero =>
    simp only [Function.iterate_zero, id_eq, Nat.cast_zero, zero_mul, Prod.mk.injEq]
    constructor <;> module
  | succ n ih =>
    rw [Function.iterate_succ_apply', ih, leapFrog_exact_constant]
    simp only [Nat.cast_succ, Prod.mk.injEq]
    constructor <;> module

/-- **Phase-space volume (one dimension).** The partial derivatives of the step map at `(x, v)`,
for an acceleration with derivative `a'` at `x`, and the determinant of that Jacobian, which is
`1`: leap-frog is area preserving. -/
theorem leapFrog_jacobian (a : ℝ → ℝ) {a' x v dt : ℝ} (ha : HasDerivAt a a' x) :
    HasDerivAt (fun y => (leapFrog a dt (y, v)).1) (1 + dt ^ 2 * a') x ∧
      HasDerivAt (fun w => (leapFrog a dt (x, w)).1) dt v ∧
      HasDerivAt (fun y => (leapFrog a dt (y, v)).2) (dt * a') x ∧
      HasDerivAt (fun w => (leapFrog a dt (x, w)).2) 1 v ∧
      (1 + dt ^ 2 * a') * 1 - dt * (dt * a') = 1 := by
  simp only [leapFrog, smul_eq_mul]
  refine ⟨?_, ?_, ?_, ?_, by ring⟩
  · have h := (hasDerivAt_id' x).add (((ha.const_mul dt).const_add v).const_mul dt)
    convert h using 1
    ring
  · have h := (hasDerivAt_id' v).add_const (dt * a x) |>.const_mul dt |>.const_add x
    convert h using 1
    ring
  · exact (ha.const_mul dt).const_add v
  · exact (hasDerivAt_id' v).add_const (dt * a x)

end LeapFrog

/-! ## Kinetic energy -/

section Kinetic

variable {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F]

/-- The kinetic energy at time `t` of an atom of mass `m` from its two half-step velocities, as
GROMOS takes it: the average of the two half-step energies. -/
noncomputable def halfStepKinetic (m : ℝ) (vOld vNew : F) : ℝ :=
  1 / 2 * m * ((‖vNew‖ ^ 2 + ‖vOld‖ ^ 2) / 2)

/-- The average of two kinetic energies is the energy of the average velocity plus the energy
of half the difference (the parallelogram law). -/
theorem halfStep_average (m : ℝ) (vOld vNew : F) :
    halfStepKinetic m vOld vNew =
      1 / 2 * m * ‖(1 / 2 : ℝ) • (vNew + vOld)‖ ^ 2 +
        1 / 2 * m * ‖(1 / 2 : ℝ) • (vNew - vOld)‖ ^ 2 := by
  simp only [halfStepKinetic, norm_smul, mul_pow, norm_add_sq_real, norm_sub_sq_real]
  norm_num
  ring

/-- **Against velocity Verlet.** After the kick `v(t + dt/2) = v(t − dt/2) + a dt`, the average of
the half-step energies is the energy of the full-step velocity `v(t) = v(t − dt/2) + a dt/2` (the
one velocity Verlet uses) plus `m dt² |a|²/8`: equal to first order in `dt`, never smaller. -/
theorem halfStep_vs_verlet (m dt : ℝ) (vOld a : F) :
    halfStepKinetic m vOld (vOld + dt • a) =
      1 / 2 * m * ‖vOld + (dt / 2) • a‖ ^ 2 + m * dt ^ 2 / 8 * ‖a‖ ^ 2 := by
  rw [halfStep_average]
  have h1 : (1 / 2 : ℝ) • (vOld + dt • a + vOld) = vOld + (dt / 2) • a := by module
  have h2 : (1 / 2 : ℝ) • (vOld + dt • a - vOld) = (dt / 2) • a := by module
  rw [h1, h2, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs]
  ring

/-- Total mass of a list of `(mass, velocity)` atoms. -/
def totalMass (l : List (ℝ × F)) : ℝ := (l.map Prod.fst).sum

/-- Total momentum `Σ m v`. -/
def momentum (l : List (ℝ × F)) : F := (l.map fun p => p.1 • p.2).sum

/-- The centre-of-mass velocity `Σ m v / M`. -/
noncomputable def comVelocity (l : List (ℝ × F)) : F := (totalMass l)⁻¹ • momentum l

/-- The kinetic energy `Σ ½ m |v|²`. -/
noncomputable def kinetic (l : List (ℝ × F)) : ℝ := (l.map fun p => 1 / 2 * p.1 * ‖p.2‖ ^ 2).sum

/-- The atoms seen from a frame moving with velocity `w`. -/
def shift (w : F) (l : List (ℝ × F)) : List (ℝ × F) := l.map fun p => (p.1, p.2 - w)

/-- The kinetic energy seen from a frame moving with any velocity `w`. -/
theorem kinetic_shift (w : F) (l : List (ℝ × F)) :
    kinetic (shift w l) =
      kinetic l - inner ℝ (momentum l) w + 1 / 2 * totalMass l * ‖w‖ ^ 2 := by
  induction l with
  | nil => simp [kinetic, shift, momentum, totalMass]
  | cons p l ih =>
    simp only [shift, kinetic, momentum, totalMass, List.map_cons, List.sum_cons] at ih ⊢
    rw [ih, norm_sub_sq_real, inner_add_left, real_inner_smul_left]
    ring

/-- **König's theorem.** The kinetic energy of a group is that of its centre of mass,
`½ M |V|²`, plus the kinetic energy relative to the centre of mass. -/
theorem koenig (l : List (ℝ × F)) (hM : totalMass l ≠ 0) :
    kinetic l = 1 / 2 * totalMass l * ‖comVelocity l‖ ^ 2 + kinetic (shift (comVelocity l) l) := by
  have hp : momentum l = totalMass l • comVelocity l := by
    simp only [comVelocity, smul_smul, mul_inv_cancel₀ hM, one_smul]
  rw [kinetic_shift, hp, real_inner_smul_left, real_inner_self_eq_norm_sq]
  ring

omit [InnerProductSpace ℝ F] in
/-- The kinetic energy of non-negative masses is not negative. -/
theorem kinetic_nonneg (l : List (ℝ × F)) (hm : ∀ p ∈ l, 0 ≤ p.1) : 0 ≤ kinetic l := by
  induction l with
  | nil => simp [kinetic]
  | cons p l ih =>
    simp only [kinetic, List.map_cons, List.sum_cons, List.mem_cons] at hm ih ⊢
    have h0 : 0 ≤ p.1 := hm p (Or.inl rfl)
    have := ih fun q hq => hm q (Or.inr hq)
    positivity

/-- A group's atoms as `(m, v_old, v_new)`: its kinetic energy at time `t`, the half-step average
of every atom, summed. -/
noncomputable def groupKinetic (g : List (ℝ × F × F)) : ℝ :=
  (g.map fun a => halfStepKinetic a.1 a.2.1 a.2.2).sum

/-- The atoms with their averaged velocities `(v_old + v_new)/2`, whose centre of mass GROMOS
books to the COM bath. -/
noncomputable def averaged (g : List (ℝ × F × F)) : List (ℝ × F) :=
  g.map fun a => (a.1, (1 / 2 : ℝ) • (a.2.2 + a.2.1))

/-- The half kicks `(v_new − v_old)/2`. -/
noncomputable def halfKicks (g : List (ℝ × F × F)) : List (ℝ × F) :=
  g.map fun a => (a.1, (1 / 2 : ℝ) • (a.2.2 - a.2.1))

/-- **The split GROMOS books.** The kinetic energy of a group at time `t` is the centre-of-mass
energy of the averaged velocities, plus the averaged velocities relative to their centre of mass,
plus the half kicks. The first term goes to the COM bath, the other two to the internal one. -/
theorem group_split (g : List (ℝ × F × F)) (hM : totalMass (averaged g) ≠ 0) :
    groupKinetic g =
      1 / 2 * totalMass (averaged g) * ‖comVelocity (averaged g)‖ ^ 2 +
        (kinetic (shift (comVelocity (averaged g)) (averaged g)) + kinetic (halfKicks g)) := by
  have hsplit : groupKinetic g = kinetic (averaged g) + kinetic (halfKicks g) := by
    clear hM
    induction g with
    | nil => simp [groupKinetic, kinetic, averaged, halfKicks]
    | cons a g ih =>
      simp only [groupKinetic, kinetic, averaged, halfKicks, List.map_cons, List.sum_cons] at ih ⊢
      rw [ih, halfStep_average]
      ring
  rw [hsplit, koenig _ hM]
  ring

/-- The internal part GROMOS books (`total − com`) is never negative for non-negative masses. -/
theorem group_internal_nonneg (g : List (ℝ × F × F)) (hM : totalMass (averaged g) ≠ 0)
    (hm : ∀ a ∈ g, 0 ≤ a.1) :
    0 ≤ groupKinetic g - 1 / 2 * totalMass (averaged g) * ‖comVelocity (averaged g)‖ ^ 2 := by
  rw [group_split g hM, add_sub_cancel_left]
  have hav : ∀ p ∈ averaged g, 0 ≤ p.1 := by
    intro p hp
    obtain ⟨a, ha, rfl⟩ := List.mem_map.mp hp
    exact hm a ha
  have h1 : ∀ p ∈ shift (comVelocity (averaged g)) (averaged g), 0 ≤ p.1 := by
    intro p hp
    obtain ⟨q, hq, rfl⟩ := List.mem_map.mp hp
    exact hav q hq
  have h2 : ∀ p ∈ halfKicks g, 0 ≤ p.1 := by
    intro p hp
    obtain ⟨a, ha, rfl⟩ := List.mem_map.mp hp
    exact hm a ha
  exact add_nonneg (kinetic_nonneg _ h1) (kinetic_nonneg _ h2)

end Kinetic

/-! ## Temperature -/

/-- The temperature of the kinetic energy `K` held by `N` degrees of freedom. -/
noncomputable def temperatureOf (K N kB : ℝ) : ℝ := 2 * K / (N * kB)

/-- **Equipartition, inverted.** `T = 2K/(N k_B)` exactly when `K = ½ N k_B T`. -/
theorem temperatureOf_iff {K N kB T : ℝ} (h : N * kB ≠ 0) :
    temperatureOf K N kB = T ↔ K = 1 / 2 * N * kB * T := by
  unfold temperatureOf
  rw [div_eq_iff h]
  constructor <;> intro h' <;> linarith

end NewenTheory.Dynamics.Integration

end
