import NewenSpec.Float
import NewenTheory.Dynamics.Thermostats

/-!
# Model of the thermostat kernels of `crates/newen-integrators`

Rust: `crates/newen-integrators/src/algorithms/berendsen_thermostat.rs`, `berendsen_scale`, and
`nosehoover_thermostat.rs`, `nose_hoover_scale` (spec keys `thermostat.berendsen_scale` and
`thermostat.nose_hoover_scale` in `spec/SPEC.toml`).
Theory: `NewenTheory/Dynamics/Thermostats.lean` (`berendsenScale`, `nhStep`, `nhcStep`).

One definition per kernel, written over any type `α` with the operations the Rust uses
(`+ − × ÷`, `sqrt` from `NewenSpec/Float.lean`, `<`, `==`, the literals `0 1 2`), decision by decision and operation by
operation in the Rust order. It is used twice (the pattern of LeanLJ, Ugwuanyi et al. 2025):

1. at `ℝ`, where `berendsenModel_real_eq` and `noseHooverModel_real_eq` prove it equal to the
   theory once the fallbacks are ruled out;
2. at `Float` (IEEE 754 binary64; `+ − × ÷ sqrt` correctly rounded, like Rust's `f64`), where
   `lake exe spec-vectors` evaluates it and `crates/newen-integrators/tests/thermostat_spec.rs`
   compares the Rust kernels bit for bit.

What is not proved: the distance between the two instantiations (rounding). The ℝ theorems say
what the formula means; the vectors say the Rust computes the same floating-point formula.

Mathlib is used through the theory (`Real.sqrt`).
-/

namespace NewenSpec.Thermostat

open NewenTheory.Dynamics.Thermostats

/-! ## 1. Model (mirrors the Rust functions) -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [LT α] [DecidableRel (α := α) (· < ·)]
  [BEq α] [Sqrt α] [OfNat α 0] [OfNat α 1] [OfNat α 2]

/-- The bath temperature both kernels start from: `2·E_kin / (N_df·k_B)`, `0` without degrees of
freedom, and the reference temperature when that is below `eps` (gromosXX's divide-by-zero
measure). -/
def freeTempModel (kB eps ekin dof temperature : α) : α :=
  let freeTemp := if 0 < dof then 2 * ekin / (dof * kB) else 0
  if freeTemp < eps then temperature else freeTemp

/-- `berendsen_scale(ekin, dof, temperature, tau, dt)`, with `k_B` and `ZERO_TEMPERATURE` as
`kB` and `eps`. -/
def berendsenModel (kB eps ekin dof temperature tau dt : α) : α :=
  if tau < 0 then 1
  else
    let freeTemp := freeTempModel α kB eps ekin dof temperature
    if freeTemp < eps then 1
    else if tau == 0 then Sqrt.sqrt (temperature / freeTemp)
    else Sqrt.sqrt (1 + dt / tau * (temperature / freeTemp - 1))

/-- The chain links `1 …` of `nose_hoover_scale`, from the tail down: `zb` is the old value of
the link below the first one, `tvb` its `tv`, `tvUp = τ·τ/N_df` the `tv` of every upper link. -/
def chainUpperModel (dof tvUp dt : α) : α → α → List α → List α
  | _, _, [] => []
  | zb, tvb, [z] => [z + (tvb * zb * zb - 1 / dof) / tvUp * dt]
  | zb, tvb, z :: z' :: rest =>
    let up := chainUpperModel dof tvUp dt z tvUp (z' :: rest)
    (z + ((tvb * zb * zb - 1 / dof) / tvUp - z * up.headD 0) * dt) :: up

/-- `nose_hoover_scale(zeta, ekin, dof, temperature, tau, dt)`: the new chain and the scaling
factor. An empty chain is outside the Rust function's domain (`apply` refuses it). -/
def noseHooverModel (kB eps : α) (zeta : List α) (ekin dof temperature tau dt : α) :
    List α × α :=
  if tau < 0 then (zeta, 1)
  else
    let freeTemp := freeTempModel α kB eps ekin dof temperature
    match zeta with
    | [] => ([], 1)
    | [z0] =>
      let z0' := z0 + dt / (tau * tau) * (freeTemp / temperature - 1)
      ([z0'], 1 - z0' * dt)
    | z0 :: z1 :: rest =>
      let tau2 := tau * tau
      let up := chainUpperModel α dof (tau2 / dof) dt z0 tau2 (z1 :: rest)
      let z0' := z0 + ((freeTemp / temperature - 1) / tau2 - up.headD 0 * z0) * dt
      (z0' :: up, 1 - z0' * dt)

/-- `τ < 0` is no coupling, for every `α`. -/
theorem berendsenModel_uncoupled {kB eps ekin dof temperature tau dt : α} (h : tau < 0) :
    berendsenModel α kB eps ekin dof temperature tau dt = 1 := by
  simp [berendsenModel, h]

omit [BEq α] [Sqrt α] in
/-- `τ < 0` is no coupling and leaves the chain alone, for every `α`. -/
theorem noseHooverModel_uncoupled {kB eps ekin dof temperature tau dt : α} {zeta : List α}
    (h : tau < 0) :
    noseHooverModel α kB eps zeta ekin dof temperature tau dt = (zeta, 1) := by
  simp [noseHooverModel, h]

end Model

/-! ## 2. The `Float` instantiation (what `lake exe spec-vectors` runs)

The constants are given by their bits, not as decimal literals, so that they are exactly Rust's:
Lean's `OfScientific Float` is not promised to round like Rust's float parser. The Rust test
asserts the same bit patterns. -/

/-- `newen_core::units::kB = 0.00831441` (kJ mol⁻¹ K⁻¹, the gromosXX value), as `f64` bits. -/
def kBFloat : Float := Float.ofBits 0x3f8107253848c6c3

/-- `ZERO_TEMPERATURE = 1e-12` (K, gromosXX `math::epsilon`), as `f64` bits. -/
def zeroTemperatureFloat : Float := Float.ofBits 0x3d719799812dea11

/-- `berendsen_scale` in binary64. -/
def berendsenFloat (ekin dof temperature tau dt : Float) : Float :=
  berendsenModel Float kBFloat zeroTemperatureFloat ekin dof temperature tau dt

/-- `nose_hoover_scale` in binary64. -/
def noseHooverFloat (zeta : List Float) (ekin dof temperature tau dt : Float) :
    List Float × Float :=
  noseHooverModel Float kBFloat zeroTemperatureFloat zeta ekin dof temperature tau dt

/-! ## 3. The ℝ instantiation equals the theory -/

/-- The bath temperature over ℝ, when there are degrees of freedom and it is not below `eps`. -/
theorem freeTempModel_real {kB eps ekin dof temperature : ℝ} (hdof : 0 < dof)
    (hft : eps ≤ 2 * ekin / (dof * kB)) :
    freeTempModel ℝ kB eps ekin dof temperature = 2 * ekin / (dof * kB) := by
  simp [freeTempModel, hdof, not_lt.mpr hft]

/-- Coupled Berendsen (`τ > 0`) over ℝ is the theory's `berendsenScale` of the bath
temperature `2·E_kin / (N_df·k_B)`, once the fallbacks are ruled out. -/
theorem berendsenModel_real_eq {kB eps ekin dof temperature tau dt : ℝ} (hτ : 0 < tau)
    (hdof : 0 < dof) (hft : eps ≤ 2 * ekin / (dof * kB)) :
    berendsenModel ℝ kB eps ekin dof temperature tau dt =
      berendsenScale (2 * ekin / (dof * kB)) temperature dt tau := by
  simp [berendsenModel, freeTempModel_real hdof hft, not_lt.mpr hτ.le, not_lt.mpr hft,
    hτ.ne', berendsenScale, Sqrt.sqrt]

/-- `τ = 0` over ℝ: instantaneous scaling to the reference temperature,
`λ = sqrt(T₀ / T)`. -/
theorem berendsenModel_real_instant {kB eps ekin dof temperature dt : ℝ} (hdof : 0 < dof)
    (hft : eps ≤ 2 * ekin / (dof * kB)) :
    berendsenModel ℝ kB eps ekin dof temperature 0 dt =
      Real.sqrt (temperature / (2 * ekin / (dof * kB))) := by
  simp [berendsenModel, freeTempModel_real hdof hft, not_lt.mpr hft, Sqrt.sqrt]

/-- One coupled Nosé–Hoover link (`τ ≥ 0`) over ℝ is the theory's `nhStep` of the bath
temperature, once the fallbacks are ruled out. -/
theorem noseHooverModel_real_eq {kB eps ekin dof temperature tau dt zeta : ℝ} (hτ : 0 ≤ tau)
    (hdof : 0 < dof) (hft : eps ≤ 2 * ekin / (dof * kB)) :
    noseHooverModel ℝ kB eps [zeta] ekin dof temperature tau dt =
      ([(nhStep zeta (2 * ekin / (dof * kB)) temperature dt tau).1],
        (nhStep zeta (2 * ekin / (dof * kB)) temperature dt tau).2) := by
  simp [noseHooverModel, freeTempModel_real hdof hft, not_lt.mpr hτ, nhStep, pow_two]

/-- The chain links over ℝ are the theory's `nhcUpper`. -/
theorem chainUpperModel_real_eq (dof tvUp dt : ℝ) :
    ∀ (zb tvb : ℝ) (zs : List ℝ),
      chainUpperModel ℝ dof tvUp dt zb tvb zs = nhcUpper dof tvUp dt zb tvb zs
  | _, _, [] => by simp [chainUpperModel, nhcUpper]
  | zb, tvb, [z] => by simp [chainUpperModel, nhcUpper, pow_two, mul_assoc]
  | zb, tvb, z :: z' :: rest => by
    simp only [chainUpperModel, nhcUpper, chainUpperModel_real_eq dof tvUp dt z tvUp (z' :: rest),
      pow_two, mul_assoc]

/-- A coupled Nosé–Hoover chain of two or more links (`τ ≥ 0`) over ℝ is the theory's
`nhcStep` of the bath temperature, once the fallbacks are ruled out. -/
theorem noseHooverModel_real_eq_chain {kB eps ekin dof temperature tau dt z0 z1 : ℝ}
    {rest : List ℝ} (hτ : 0 ≤ tau) (hdof : 0 < dof) (hft : eps ≤ 2 * ekin / (dof * kB)) :
    noseHooverModel ℝ kB eps (z0 :: z1 :: rest) ekin dof temperature tau dt =
      nhcStep (z0 :: z1 :: rest) (2 * ekin / (dof * kB)) temperature dt tau dof := by
  simp only [noseHooverModel, freeTempModel_real hdof hft, not_lt.mpr hτ, ite_false, nhcStep,
    chainUpperModel_real_eq, pow_two]

end NewenSpec.Thermostat
