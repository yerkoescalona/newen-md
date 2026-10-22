import NewenSpec.Float
import NewenTheory.Dynamics.Integration
import NewenTheory.Dynamics.Pressure

/-!
# Model of the dynamics kernels of `crates/newen-integrators`

Rust (spec keys in `spec/SPEC.toml`):

* `dynamics.leap_frog_velocity`, `dynamics.leap_frog_position`: `leap_frog_velocity` and
  `leap_frog_position` in `src/algorithms/leap_frog.rs`, one atom's leap-frog step;
* `dynamics.kinetic_energy`, `dynamics.group_kinetic_energy`, `dynamics.kinetic_temperature`:
  `kinetic_energy`, `group_kinetic_energy` and `kinetic_temperature` in
  `src/algorithms/temperature.rs`, the kinetic energy of the atoms, of one temperature group with
  its centre-of-mass part, and the temperature of a kinetic energy;
* `dynamics.pressure_tensor`, `dynamics.scalar_pressure`: `pressure_tensor` and
  `scalar_pressure` in `src/algorithms/pressure_calculation.rs`;
* `dynamics.berendsen_barostat_scale`: `berendsen_barostat_scale` in
  `src/algorithms/berendsen_barostat.rs`.

Theory: `NewenTheory/Dynamics/Integration.lean` and `NewenTheory/Dynamics/Pressure.lean`.

Each kernel is written once, over any type `α` with the operations the Rust uses (`+ − × ÷`, `<`,
`powf` from `NewenSpec/Float.lean`, the literals `0 1 2 3`), operation by operation in the Rust
order, and used twice, as in `NewenSpec/Thermostat.lean`:

1. at `ℝ`, where the `*_real_eq` theorems prove it equal to the theory;
2. at `Float`, where `lake exe spec-vectors` evaluates it and
   `crates/newen-integrators/tests/dynamics_spec.rs` compares the Rust kernels bit for bit.

A `glam::DVec3` is a `V3`, its operations component by component as glam writes them
(`length_squared` is `x·x + y·y + z·z`, left to right). A sum over atoms is a Rust loop
`acc += term` from `0.0`, so the model folds the list from the left (`sumModel`). The two
literals that are not small integers, the `0.5` of the kinetic energies and the pressure, and the
`1.0 / 3.0` of the barostat, are parameters, given at `Float` by their bits.

The barostat calls `powf`, which IEEE 754 does not require to be correctly rounded: Rust's
`f64::powf` and Lean's `Float.pow` both call the C library's `pow`, so the bits agree when both
run against the same `libm`.

What is not proved: the rounding error between the two instantiations.
-/

namespace NewenSpec.Dynamics

open NewenTheory.Dynamics.Integration NewenTheory.Dynamics.Pressure

/-- A `glam::DVec3`. -/
structure V3 (α : Type) where
  x : α
  y : α
  z : α

/-- A newen `Mat3` (`glam::DMat3`) by its columns `x_axis`, `y_axis`, `z_axis`; newen keeps
GROMOS's row `i` as column `i`. -/
structure M3 (α : Type) where
  c0 : V3 α
  c1 : V3 α
  c2 : V3 α

/-- One atom as the kinetic-energy kernels read it: mass, new and old half-step velocity. -/
structure Atom (α : Type) where
  m : α
  vNew : V3 α
  vOld : V3 α

/-- `DVec3::length_squared`: `x·x + y·y + z·z`. -/
def V3.lengthSq {α : Type} [Add α] [Mul α] (a : V3 α) : α := a.x * a.x + a.y * a.y + a.z * a.z

/-! ## 1. Model (mirrors the Rust functions) -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [LT α] [DecidableRel (α := α) (· < ·)]
  [Powf α]
  [OfNat α 0] [OfNat α 1] [OfNat α 2] [OfNat α 3]

/-- `leap_frog_velocity(v, f, inverse_mass, dt)`: `v + f·m⁻¹·dt`, component by component. -/
def leapFrogVelocityModel (v f : V3 α) (im dt : α) : V3 α :=
  ⟨v.x + f.x * im * dt, v.y + f.y * im * dt, v.z + f.z * im * dt⟩

/-- `leap_frog_position(x, v, dt)`: `x + v·dt`, component by component. -/
def leapFrogPositionModel (x v : V3 α) (dt : α) : V3 α :=
  ⟨x.x + v.x * dt, x.y + v.y * dt, x.z + v.z * dt⟩

/-- A Rust sum loop, `acc += f(b)` from `0` over `l` in order. -/
def sumModel {β : Type} (f : β → α) (l : List β) : α := l.foldl (fun acc b => acc + f b) 0

/-- `kinetic_energy(mass, vel_new, vel_old)`: `(Σ ½·m·(|v_new|² + |v_old|²)/2, Σ ½·m·|v_new|²)`,
with the literal `0.5` as `half`. -/
def kineticEnergyModel (half : α) (atoms : List (Atom α)) : α × α :=
  (sumModel α (fun a => half * a.m * (a.vNew.lengthSq + a.vOld.lengthSq) / 2) atoms,
    sumModel α (fun a => half * a.m * a.vNew.lengthSq) atoms)

/-- `0.5 * (new_v + old_v)` of one atom. -/
def averagedModel (half : α) (a : Atom α) : V3 α :=
  ⟨half * (a.vNew.x + a.vOld.x), half * (a.vNew.y + a.vOld.y), half * (a.vNew.z + a.vOld.z)⟩

/-- `com_v += m * v` over the atoms, then `com_v /= tot_mass`. -/
def comModel (atoms : List (Atom α)) (vel : Atom α → V3 α) (M : α) : V3 α :=
  ⟨sumModel α (fun a => a.m * (vel a).x) atoms / M, sumModel α (fun a => a.m * (vel a).y) atoms / M,
    sumModel α (fun a => a.m * (vel a).z) atoms / M⟩

/-- `group_kinetic_energy(mass, vel_new, vel_old)`: `(com, total, new_com, new_total)`, with the
literal `0.5` as `half`. -/
def groupKineticModel (half : α) (atoms : List (Atom α)) : α × α × α × α :=
  let M := sumModel α (·.m) atoms
  let com := comModel α atoms (averagedModel α half) M
  let newCom := comModel α atoms (·.vNew) M
  (half * M * com.lengthSq,
    sumModel α (fun a => a.m * (a.vNew.lengthSq + a.vOld.lengthSq) / 2) atoms * half,
    half * M * newCom.lengthSq,
    sumModel α (fun a => a.m * a.vNew.lengthSq) atoms * half)

/-- `kinetic_temperature(ekin, dof, k_boltzmann)`: `2·E_kin/k_B/N_df`, `0` for `N_df ≤ 0` (and
for a NaN `N_df`), as the thermostats. -/
def kineticTemperatureModel (ekin dof kB : α) : α :=
  if 0 < dof then 2 * ekin / kB / dof else 0

/-- `pressure_tensor(kinetic, virial, volume)`: `(K + W·½)·(2/V)` element by element, with the
literal `0.5` as `half`. -/
def pressureTensorModel (half : α) (K W : M3 α) (V : α) : M3 α :=
  let factor := 2 / V
  let col (k w : V3 α) : V3 α :=
    ⟨(k.x + w.x * half) * factor, (k.y + w.y * half) * factor, (k.z + w.z * half) * factor⟩
  ⟨col K.c0 W.c0, col K.c1 W.c1, col K.c2 W.c2⟩

/-- `scalar_pressure(p)`: `(p_xx + p_yy + p_zz)/3`. -/
def scalarPressureModel (p : M3 α) : α := (p.c0.x + p.c1.y + p.c2.z) / 3

/-- `berendsen_barostat_scale(pressure, pressure0, compressibility, tau, dt)`:
`powf(1 − β·dt/τ·(P₀ − P), 1/3)`, with the literal `1.0 / 3.0` as `third`. -/
def berendsenBarostatModel (third pressure pressure0 compressibility tau dt : α) : α :=
  Powf.powf (1 - compressibility * dt / tau * (pressure0 - pressure)) third

end Model

/-! ## 2. The `Float` instantiation (what `lake exe spec-vectors` runs)

The two decimal literals are given by their bits, so that they are exactly Rust's. -/

/-- The `0.5` of the kinetic energies and of the pressure tensor, as `f64` bits. -/
def halfFloat : Float := Float.ofBits 0x3fe0000000000000

/-- The `1.0 / 3.0` of the barostat, as `f64` bits (Rust folds the division at compile time). -/
def thirdFloat : Float := Float.ofBits 0x3fd5555555555555

def leapFrogVelocityFloat (v f : V3 Float) (im dt : Float) : V3 Float :=
  leapFrogVelocityModel Float v f im dt

def leapFrogPositionFloat (x v : V3 Float) (dt : Float) : V3 Float :=
  leapFrogPositionModel Float x v dt

def kineticEnergyFloat (atoms : List (Atom Float)) : Float × Float :=
  kineticEnergyModel Float halfFloat atoms

def groupKineticFloat (atoms : List (Atom Float)) : Float × Float × Float × Float :=
  groupKineticModel Float halfFloat atoms

def kineticTemperatureFloat (ekin dof kB : Float) : Float :=
  kineticTemperatureModel Float ekin dof kB

def pressureTensorFloat (K W : M3 Float) (V : Float) : M3 Float :=
  pressureTensorModel Float halfFloat K W V

def scalarPressureFloat (p : M3 Float) : Float := scalarPressureModel Float p

def berendsenBarostatFloat (pressure pressure0 compressibility tau dt : Float) : Float :=
  berendsenBarostatModel Float thirdFloat pressure pressure0 compressibility tau dt

/-! ## 3. The ℝ instantiation equals the theory -/

/-- A `V3` over ℝ as a point of `ℝ × ℝ × ℝ` (the leap-frog step needs only the vector space). -/
def V3.toProd (a : V3 ℝ) : ℝ × ℝ × ℝ := (a.x, a.y, a.z)

/-- A `V3` over ℝ as a Euclidean vector (the kinetic energies need its norm). -/
noncomputable def V3.toE (a : V3 ℝ) : R3 := !₂[a.x, a.y, a.z]

theorem V3.norm_toE_sq (a : V3 ℝ) : ‖a.toE‖ ^ 2 = a.lengthSq := by
  rw [V3.toE, EuclideanSpace.norm_sq_eq]
  simp [Fin.sum_univ_three, V3.lengthSq, sq]

/-- One leap-frog step over ℝ, the velocity kernel then the position kernel, is the theory's
`leapFrog` with the acceleration `m⁻¹·f`. -/
theorem leapFrogModel_real_eq (x v f : V3 ℝ) (im dt : ℝ) :
    ((leapFrogPositionModel ℝ x (leapFrogVelocityModel ℝ v f im dt) dt).toProd,
        (leapFrogVelocityModel ℝ v f im dt).toProd) =
      leapFrog (fun _ => im • f.toProd) dt (x.toProd, v.toProd) := by
  simp only [leapFrogPositionModel, leapFrogVelocityModel, V3.toProd, leapFrog, Prod.smul_mk,
    Prod.mk_add_mk, smul_eq_mul, Prod.mk.injEq]
  refine ⟨⟨?_, ?_, ?_⟩, ?_, ?_, ?_⟩ <;> ring

/-- A sum loop over ℝ is the list sum. -/
theorem sumModel_real {β : Type} (f : β → ℝ) (l : List β) : sumModel ℝ f l = (l.map f).sum := by
  have h : ∀ acc : ℝ, l.foldl (fun acc b => acc + f b) acc = acc + (l.map f).sum := by
    induction l with
    | nil => simp
    | cons b l ih => intro acc; simp [ih]; ring
  simp [sumModel, h]

/-- The atoms of a kernel call as the theory's `(m, v_old, v_new)`. -/
noncomputable def theoryAtoms (atoms : List (Atom ℝ)) :
    List (ℝ × R3 × R3) :=
  atoms.map fun a => (a.m, a.vOld.toE, a.vNew.toE)

/-- The atoms with their new velocities, `(m, v_new)`. -/
noncomputable def newVelocities (atoms : List (Atom ℝ)) : List (ℝ × R3) :=
  atoms.map fun a => (a.m, a.vNew.toE)

/-- The kinetic energy over ℝ, with `half = 1/2`: the theory's half-step average summed over the
atoms, and the kinetic energy of the new velocities. -/
theorem kineticEnergyModel_real_eq (atoms : List (Atom ℝ)) :
    kineticEnergyModel ℝ (1 / 2) atoms =
      (groupKinetic (theoryAtoms atoms), kinetic (newVelocities atoms)) := by
  simp only [kineticEnergyModel, sumModel_real, groupKinetic, kinetic, theoryAtoms, newVelocities,
    List.map_map, Prod.mk.injEq]
  refine ⟨?_, ?_⟩ <;> congr 1 <;> apply List.map_congr_left <;> intro a _
  · simp only [Function.comp_apply, halfStepKinetic, V3.norm_toE_sq]
    ring
  · simp only [Function.comp_apply, V3.norm_toE_sq]

/-- `Σ m·u` component by component, as a Euclidean vector, is the theory's momentum. -/
theorem toE_momentum (atoms : List (Atom ℝ)) (vel : Atom ℝ → V3 ℝ) :
    (⟨sumModel ℝ (fun a => a.m * (vel a).x) atoms, sumModel ℝ (fun a => a.m * (vel a).y) atoms,
        sumModel ℝ (fun a => a.m * (vel a).z) atoms⟩ : V3 ℝ).toE =
      momentum (atoms.map fun a => (a.m, (vel a).toE)) := by
  simp only [sumModel_real, momentum, List.map_map]
  induction atoms with
  | nil => ext i; fin_cases i <;> simp [V3.toE]
  | cons a atoms ih =>
    simp only [List.map_cons, List.sum_cons, Function.comp_apply] at ih ⊢
    rw [← ih]
    ext i
    fin_cases i <;> simp [V3.toE]

/-- The centre-of-mass velocity of the model over ℝ is the theory's. -/
theorem comModel_real_eq (atoms : List (Atom ℝ)) (vel : Atom ℝ → V3 ℝ) :
    (comModel ℝ atoms vel (sumModel ℝ (·.m) atoms)).toE =
      comVelocity (atoms.map fun a => (a.m, (vel a).toE)) := by
  have hM : sumModel ℝ (·.m) atoms = totalMass (atoms.map fun a => (a.m, (vel a).toE)) := by
    simp [sumModel_real, totalMass, List.map_map, Function.comp_def]
  rw [comVelocity, ← toE_momentum, hM]
  ext i
  fin_cases i <;> simp [comModel, V3.toE, div_eq_inv_mul]

/-- One temperature group over ℝ, with `half = 1/2`: the centre-of-mass energy of the averaged
velocities and the half-step kinetic energy of the theory's `group_split`, and König's
centre-of-mass energy and kinetic energy of the new velocities. -/
theorem groupKineticModel_real_eq (atoms : List (Atom ℝ)) :
    groupKineticModel ℝ (1 / 2) atoms =
      (1 / 2 * totalMass (averaged (theoryAtoms atoms)) *
          ‖comVelocity (averaged (theoryAtoms atoms))‖ ^ 2,
        groupKinetic (theoryAtoms atoms),
        1 / 2 * totalMass (newVelocities atoms) * ‖comVelocity (newVelocities atoms)‖ ^ 2,
        kinetic (newVelocities atoms)) := by
  have hav : averaged (theoryAtoms atoms) =
      atoms.map fun a => (a.m, (averagedModel ℝ (1 / 2) a).toE) := by
    simp only [averaged, theoryAtoms, List.map_map]
    apply List.map_congr_left
    intro a _
    simp only [Function.comp_apply, Prod.mk.injEq, true_and]
    ext i
    fin_cases i <;> simp [V3.toE, averagedModel, mul_add]
  have hM : ∀ vel : Atom ℝ → V3 ℝ,
      sumModel ℝ (·.m) atoms = totalMass (atoms.map fun a => (a.m, (vel a).toE)) := by
    intro vel
    simp [sumModel_real, totalMass, List.map_map, Function.comp_def]
  rw [hav]
  simp only [groupKineticModel, Prod.mk.injEq]
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [← comModel_real_eq, V3.norm_toE_sq, ← hM]
  · simp only [sumModel_real, groupKinetic, theoryAtoms, List.map_map, ← List.sum_map_mul_right]
    congr 1
    apply List.map_congr_left
    intro a _
    simp [halfStepKinetic, V3.norm_toE_sq]
    ring
  · rw [newVelocities, ← comModel_real_eq, V3.norm_toE_sq, ← hM]
  · simp only [sumModel_real, kinetic, newVelocities, List.map_map, ← List.sum_map_mul_right]
    congr 1
    apply List.map_congr_left
    intro a _
    simp [V3.norm_toE_sq]
    ring

/-- The temperature over ℝ, with degrees of freedom, is the theory's `2K/(N k_B)`. -/
theorem kineticTemperatureModel_real_eq {ekin dof kB : ℝ} (hdof : 0 < dof) :
    kineticTemperatureModel ℝ ekin dof kB = temperatureOf ekin dof kB := by
  simp [kineticTemperatureModel, hdof, temperatureOf, div_div, mul_comm]

/-- A `M3` over ℝ as the GROMOS matrix: row `i` is newen's column `i`. -/
def M3.toMatrix (p : M3 ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  !![p.c0.x, p.c0.y, p.c0.z; p.c1.x, p.c1.y, p.c1.z; p.c2.x, p.c2.y, p.c2.z]

/-- The pressure tensor over ℝ, with `half = 1/2`, is the theory's `2/V · (K − Ξ)` with the
GROMOS virial `Ξ = −½ W` of newen's stored `W = Σ r⊗f`. -/
theorem pressureTensorModel_real_eq (K W : M3 ℝ) (V : ℝ) :
    (pressureTensorModel ℝ (1 / 2) K W V).toMatrix =
      pressureTensor K.toMatrix (-(1 / 2 : ℝ) • W.toMatrix) V := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [pressureTensorModel, M3.toMatrix, pressureTensor] <;> ring

/-- The scalar pressure over ℝ is the theory's, a third of the trace. -/
theorem scalarPressureModel_real_eq (p : M3 ℝ) :
    scalarPressureModel ℝ p = scalarPressure p.toMatrix := by
  simp [scalarPressureModel, scalarPressure, Matrix.trace_fin_three, M3.toMatrix]

/-- The barostat factor over ℝ, with `third = 1/3`, is the theory's `barostatScale`, and its
cube is the volume factor `1 − β·dt/τ·(P₀ − P)` (`barostat_volume`). The radicand must not be
negative (`hr`): below it Rust's `powf` and `Float.pow` return NaN, while ℝ's `rpow` returns a real
number that is not a cube root. -/
theorem berendsenBarostatModel_real_eq {P P₀ β τ dt : ℝ} (hr : 0 ≤ 1 - β * dt / τ * (P₀ - P)) :
    berendsenBarostatModel ℝ (1 / 3) P P₀ β τ dt = barostatScale P P₀ β τ dt ∧
      berendsenBarostatModel ℝ (1 / 3) P P₀ β τ dt ^ 3 = 1 - β * dt / τ * (P₀ - P) :=
  ⟨rfl, barostat_volume hr⟩

end NewenSpec.Dynamics
