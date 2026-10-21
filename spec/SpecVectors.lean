import NewenSpec.Schedule
import NewenSpec.Thermo.Dof
import NewenSpec.Thermo.Bath
import NewenSpec.Thermostat
import NewenSpec.Nonbonded

/-! `lake exe spec-vectors`: evaluate the Lean models on the grids the Rust tests sweep (and the
thermostat and nonbonded models on random `Float` inputs) and write them to `vectors/`, so
`cargo test` can compare the Rust functions against Lean without Lean. -/

open NewenSpec.Schedule
open NewenSpec.Thermo.Dof
open NewenSpec.Thermo.Bath
open NewenSpec.Thermostat
open NewenSpec.Nonbonded

/-- Grid bound of the schedule grid, equal to `GRID` in
`crates/newen-logic/tests/schedule_spec.rs`. -/
def grid : Nat := 200

/-- Row `nsnb` as a string of `0`/`1`, one character per `ntwg = 1..grid`. -/
def row (nsnb : Nat) : String :=
  String.ofList ((List.range grid).map fun i => if missesModel nsnb (i + 1) then '1' else '0')

def scheduleJson : String :=
  let rows := (List.range grid).map fun i => s!"    \"{row (i + 1)}\""
  "{\n" ++
  "  \"spec\": \"schedule.output_misses_longrange\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Schedule.missesModel)\",\n" ++
  s!"  \"layout\": \"rows[nsnb - 1][ntwg - 1] for 1 <= nsnb, ntwg <= {grid}; 1 = true\",\n" ++
  "  \"rows\": [\n" ++ ",\n".intercalate rows ++ "\n  ]\n}\n"

/-- Count bound of the `total_dof` grid, equal to `GRID` in
`crates/newen-logic/tests/dof_spec.rs`. -/
def dofGrid : Nat := 10

/-- The `ndfmin` values of the `total_dof` grid: `-3 ..= 6`, negative ones included. -/
def ndfmins : List Int := (List.range 10).map fun (i : Nat) => (i : Int) - 3

/-- Row `(nAtoms, solvent, solute)`: `totalDofModel` for every `ndfmin` in `ndfmins`. -/
def dofRow (nAtoms solvent solute : Nat) : String :=
  "[" ++ ", ".intercalate (ndfmins.map fun f => toString (totalDofModel nAtoms solvent solute f))
    ++ "]"

def dofJson : String :=
  let counts := List.range (dofGrid + 1)
  let rows := counts.flatMap fun n => counts.flatMap fun s => counts.map fun c =>
    s!"    {dofRow n s c}"
  "{\n" ++
  "  \"spec\": \"thermo.dof.total_dof\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Thermo.Dof.totalDofModel)\",\n" ++
  s!"  \"layout\": \"rows[(n_atoms * {dofGrid + 1} + solvent_constraints) * {dofGrid + 1} + " ++
  s!"solute_constraints][ndfmin + 3] for 0 <= n_atoms, solvent_constraints, " ++
  s!"solute_constraints <= {dofGrid} and -3 <= ndfmin <= 6\",\n" ++
  "  \"rows\": [\n" ++ ",\n".intercalate rows ++ "\n  ]\n}\n"

/-- The subsets of `xs` with at most `k` elements, each in the order of `xs`. -/
def subsetsUpTo {α : Type} : List α → Nat → List (List α)
  | [], _ => [[]]
  | _ :: _, 0 => [[]]
  | x :: xs, k + 1 => subsetsUpTo xs (k + 1) ++ (subsetsUpTo xs k).map (x :: ·)

/-- DOFSET splits of `n` atoms into 1–3 ranges: the last atoms, ascending, ending at `n − 1`. -/
def splits (n : Nat) : List (List Nat) :=
  (subsetsUpTo (List.range (n - 1)) 2).map (· ++ [n - 1])

/-- Temperature groups of `n` atoms by their first atoms: ascending, from 0, at most 3. -/
def groupSets (n : Nat) : List (List Nat) :=
  (subsetsUpTo ((List.range n).drop 1) 2).map (0 :: ·)

/-- Constraint sets of `n` atoms: at most 2 of the pairs `i < j`. -/
def constraintSets (n : Nat) : List (List (Nat × Nat)) :=
  subsetsUpTo ((List.range n).flatMap fun i =>
    ((List.range n).filter (i < ·)).map fun j => (i, j)) 2

/-- The `bath.json` systems `(n, ranges, nBaths, groupFirsts, constraints)`: every `n ≤ 6`, every
split, 1–3 baths with up to three bath assignments `(com, ir) = ((r + a) % nBaths,
(2r + a + 1) % nBaths)` of range `r`, and, for the `k`-th of them, group set `k` and constraint
set `7k` (cyclically). -/
def bathSystems : List (Nat × List Range × Nat × List Nat × List (Nat × Nat)) :=
  let base := (List.range 6).flatMap fun n0 => (splits (n0 + 1)).flatMap fun lasts =>
    (List.range 3).flatMap fun b0 =>
      (List.range (if b0 = 0 then 1 else 3)).map fun a => (n0 + 1, lasts, b0 + 1, a)
  base.zipIdx.map fun ((n, lasts, nb, a), k) =>
    let ranges := lasts.zipIdx.map fun (l, r) => (⟨l, (r + a) % nb, (2 * r + a + 1) % nb⟩ : Range)
    let gs := groupSets n
    let cs := constraintSets n
    (n, ranges, nb, gs.getD (k % gs.length) [], cs.getD ((7 * k) % cs.length) [])

/-- One row: the inputs, then `bathHalfDofsModel` and `ndfminSpreadModel` of them. -/
def bathRow (n : Nat) (ranges : List Range) (nb : Nat) (gs : List Nat) (cs : List (Nat × Nat))
    (rottrans : Bool) (ndfmin : Int) : List Int :=
  let halves := bathHalfDofsModel nb ranges gs cs rottrans
  let spread := ndfminSpreadModel (halves.map fun h => h.1 + h.2) ndfmin
  [(n : Int), nb, if rottrans then 1 else 0, ndfmin, ranges.length]
    ++ ranges.flatMap (fun r => [(r.lastAtom : Int), r.comBath, r.irBath])
    ++ [(gs.length : Int)] ++ gs.map Int.ofNat
    ++ [(cs.length : Int)] ++ cs.flatMap (fun c => [(c.1 : Int), c.2])
    ++ halves.flatMap (fun h => [h.1, h.2]) ++ spread.1 ++ [spread.2]

def bathJson : String :=
  let rows := bathSystems.flatMap fun (n, ranges, nb, gs, cs) =>
    [false, true].flatMap fun rt => ([0, 3, 6] : List Int).map fun f =>
      "    [" ++ ", ".intercalate ((bathRow n ranges nb gs cs rt f).map toString) ++ "]"
  "{\n" ++
  "  \"spec\": \"thermo.bath.bath_half_dofs, thermo.bath.ndfmin_spread\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Thermo.Bath.bathHalfDofsModel, " ++
  "ndfminSpreadModel)\",\n" ++
  "  \"layout\": \"one system per row: n_atoms, n_baths, rottrans (0/1), ndfmin, n_ranges, " ++
  "(last_atom, com_bath, ir_bath) per range, n_groups, the group first atoms, n_constraints, " ++
  "(i, j) per constraint; then (com, ir) half dof per bath, the n_baths numerators, the " ++
  "denominator\",\n" ++
  "  \"rows\": [\n" ++ ",\n".intercalate rows ++ "\n  ]\n}\n"

/-! ## Thermostats: random inputs as `Float`, every number written as its 64 bits -/

/-- SplitMix64 (Steele, Lea & Flood 2014): the generator state and one draw. -/
def splitMix (s : UInt64) : UInt64 × UInt64 :=
  let s := s + 0x9E3779B97F4A7C15
  let z := (s ^^^ (s >>> 30)) * 0xBF58476D1CE4E5B9
  let z := (z ^^^ (z >>> 27)) * 0x94D049BB133111EB
  (z ^^^ (z >>> 31), s)

abbrev Gen := StateM UInt64

def draw : Gen UInt64 := modifyGet splitMix

/-- Uniform in `[0, 1)`: 52 random mantissa bits under the exponent of 1, minus 1. -/
def unit : Gen Float := do
  return Float.ofBits ((0x3FF0000000000000 : UInt64) ||| ((← draw) >>> 12)) - 1

def uniform (a b : Float) : Gen Float := do return a + (b - a) * (← unit)

/-- Log-uniform in `[a, b]` (`exp` and `log` only make inputs; they are written as bits). -/
def logUniform (a b : Float) : Gen Float := do
  return Float.exp (Float.log a + (Float.log b - Float.log a) * (← unit))

def below (n : Nat) : Gen Nat := do return ((← draw) % n.toUInt64).toNat

def pick (xs : List Float) : Gen Float := do return xs.getD (← below xs.length) 0

/-- 16 hex digits of the bits of `x`, quoted. -/
def hexBits (x : Float) : String :=
  let b := x.toBits
  let digit (i : Nat) : Char :=
    "0123456789abcdef".toList.getD ((b >>> (4 * (15 - i)).toUInt64) &&& 15).toNat '0'
  "\"0x" ++ String.ofList ((List.range 16).map digit) ++ "\""

def hexRow (xs : List Float) : String := "    [" ++ ", ".intercalate (xs.map hexBits) ++ "]"

def thermoDofs : List Float := [1, 2, 3, 5.5, 11.5, 12.76, 100, 6000, 1e5]
def thermoTemperatures : List Float := [0, 1, 300, 1000]
def thermoDts : List Float := [1e-4, 5e-4, 2e-3, 1e-2]

/-- A `dof`: one of `thermoDofs`, random in `[1, 1e5]`, or `0`. -/
def drawDof : Gen Float := do
  match ← below (thermoDofs.length + 2) with
  | 0 => return 0
  | 1 => uniform 1 1e5
  | i => return thermoDofs.getD (i - 2) 0

def drawTemperature : Gen Float := do
  match ← below (thermoTemperatures.length + 1) with
  | 0 => uniform 0 2000
  | i => return thermoTemperatures.getD (i - 1) 0

/-- A kinetic energy whose bath temperature is log-uniform in `[1e-14, 1e6]` K, or `0`. -/
def drawEkin (dof : Float) : Gen Float := do
  if (← below 16) == 0 then return 0
  let t ← logUniform 1e-14 1e6
  return t * (if dof > 0 then dof else 1) * kBFloat / 2

def berendsenTaus (dt : Float) : List Float := [-1, 0, dt, 0.1, 1, 10]

def berendsenRow (ekin dof temperature tau dt : Float) : List Float :=
  [ekin, dof, temperature, tau, dt, berendsenFloat ekin dof temperature tau dt]

/-- Every fixed `(dt, dof, temperature, tau)`, `dof = 0` included, with four kinetic energies
(the first `0`); then random rows with every parameter drawn from its family or at random. -/
def berendsenRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for dt in thermoDts do
    for dof in thermoDofs ++ [0] do
      for temperature in thermoTemperatures do
        for tau in berendsenTaus dt do
          for k in List.range 4 do
            let ekin ← if k == 0 then pure 0 else drawEkin dof
            rows := rows.push (berendsenRow ekin dof temperature tau dt)
  for _ in List.range 3000 do
    let dt ← pick thermoDts
    let dof ← drawDof
    let temperature ← drawTemperature
    let taus := berendsenTaus dt
    let tau ← match ← below (taus.length + 1) with
      | 0 => uniform dt 10
      | i => pure (taus.getD (i - 1) 0)
    let ekin ← drawEkin dof
    rows := rows.push (berendsenRow ekin dof temperature tau dt)
  return rows.toList

def noseHooverTaus : List Float := [-1, 0.05, 0.1, 1]

/-- One chain run: random parameters and starting `ζ` in `[-50, 50]`, then 1–8 steps, each with
a new kinetic energy and the `ζ` the previous step left; one row per step. -/
def noseHooverRun : Gen (List String) := do
  let nhc := [1, 2, 3, 5].getD (← below 4) 1
  let dt ← pick thermoDts
  let dof ← drawDof
  let temperature ← drawTemperature
  let tau ← match ← below (noseHooverTaus.length + 1) with
    | 0 => uniform 0.01 2
    | i => pure (noseHooverTaus.getD (i - 1) 0)
  let mut zeta := #[]
  for _ in List.range nhc do
    zeta := zeta.push (← uniform (-50) 50)
  let steps := 1 + (← below 8)
  let mut rows := #[]
  for _ in List.range steps do
    let ekin ← drawEkin dof
    let (zeta', scale) := noseHooverFloat zeta.toList ekin dof temperature tau dt
    rows := rows.push ("    [" ++ ", ".intercalate ([ekin, dof, temperature, tau, dt].map hexBits
      ++ [toString nhc] ++ (zeta.toList ++ zeta' ++ [scale]).map hexBits) ++ "]")
    zeta := zeta'.toArray
  return rows.toList

def noseHooverRows : Gen (List String) := do
  let mut rows := #[]
  for _ in List.range 700 do
    rows := rows ++ (← noseHooverRun).toArray
  return rows.toList

/-- The rows of `vectors/thermostat.json`, from a fixed seed. -/
def thermostatRows : List (List Float) × List String :=
  (do return (← berendsenRows, ← noseHooverRows) : Gen _).run' 0x6e6577656e2d6d64

def thermostatJson (berendsen : List (List Float)) (noseHoover : List String) : String :=
  "{\n" ++
  "  \"spec\": \"thermostat.berendsen_scale, thermostat.nose_hoover_scale\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Thermostat.berendsenFloat, " ++
  "noseHooverFloat; SplitMix64, fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits); nhc is a plain integer\",\n" ++
  "  \"layout_berendsen\": \"ekin, dof, temperature, tau, dt, scale\",\n" ++
  "  \"layout_nose_hoover\": \"ekin, dof, temperature, tau, dt, nhc, zeta_in[0..nhc), " ++
  "zeta_out[0..nhc), scale; consecutive rows of one run carry zeta_out into zeta_in\",\n" ++
  "  \"berendsen\": [\n" ++ ",\n".intercalate (berendsen.map hexRow) ++ "\n  ],\n" ++
  "  \"nose_hoover\": [\n" ++ ",\n".intercalate noseHoover ++ "\n  ]\n}\n"

/-! ## Nonbonded: the LJ + CRF pair kernels on random `Float` inputs, written as bits -/

def crfCutoffs : List Float := [0, 0.8, 0.9, 1.0, 1.4]
def crfRfEpsilons : List Float := [0, 1, 54, 61, 62, 78.5]
def crfKappas : List Float := [0, 0.5, 1]

/-- One of `xs`, or uniform in `[lo, hi)`, each with the same probability. -/
def drawFrom (xs : List Float) (lo hi : Float) : Gen Float := do
  match ← below (xs.length + 1) with
  | 0 => uniform lo hi
  | i => return xs.getD (i - 1) 0

/-- The four inputs of `CRFParameters::new`: `cutoff`, `epsilon`, `rf_epsilon`, `rf_kappa`. -/
def drawCrfInputs : Gen (List Float) := do
  let cutoff ← drawFrom crfCutoffs 0.5 2
  let epsilon ← drawFrom [1] 1 4
  let rfEpsilon ← drawFrom crfRfEpsilons 1 100
  let rfKappa ← drawFrom crfKappas 0 3
  return [cutoff, epsilon, rfEpsilon, rfKappa]

/-- `CRFParameters::new` of the four inputs. -/
def crfOf (inputs : List Float) : CrfParams Float :=
  crfParametersFloat (inputs.getD 0 0) (inputs.getD 1 0) (inputs.getD 2 0) (inputs.getD 3 0)

def crfParametersRow (inputs : List Float) : List Float :=
  let p := crfOf inputs
  inputs ++ [p.crfCut, p.crf2Cut3i, p.crfCut3i, p.cutoffSq]

/-- Every fixed `(cutoff, rf_epsilon, rf_kappa)` with `epsilon = 1` and with a random `epsilon`,
then random rows. -/
def crfParametersRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for cutoff in crfCutoffs do
    for rfEpsilon in crfRfEpsilons do
      for rfKappa in crfKappas do
        rows := rows.push (crfParametersRow [cutoff, 1, rfEpsilon, rfKappa])
        rows := rows.push (crfParametersRow [cutoff, ← uniform 1 4, rfEpsilon, rfKappa])
  for _ in List.range 420 do
    rows := rows.push (crfParametersRow (← drawCrfInputs))
  return rows.toList

/-- `0` one time in eight, else `draw`. -/
def orZero (draw : Gen Float) : Gen Float := do
  if (← below 8) == 0 then return 0 else draw

/-- A squared distance: log-uniform in `[1e-4, 4]` nm², or below the guard (`0` or log-uniform in
`[1e-16, 9e-11]`), exactly the guard, or exactly `cutoff²`. -/
def drawR2 (cutoff : Float) : Gen Float := do
  match ← below 80 with
  | 0 => return 0
  | 1 | 2 => logUniform 1e-16 9e-11
  | 3 => return zeroDistance2Float
  | 4 | 5 | 6 => return cutoff * cutoff
  | _ => logUniform 1e-4 4

def ljCrfPairRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 6000 do
    let inputs ← drawCrfInputs
    let r2 ← drawR2 (inputs.getD 0 0)
    let c6 ← orZero (logUniform 1e-5 1e-2)
    let c12 ← orZero (logUniform 1e-8 1e-5)
    let q ← orZero (uniform (-150) 150)
    let (f, eLj, eCrf) := ljCrfPairFloat r2 c6 c12 q (crfOf inputs)
    rows := rows.push ([r2, c6, c12, q] ++ inputs ++ [f, eLj, eCrf])
  return rows.toList

def rfExcludedPairRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 1500 do
    let inputs ← drawCrfInputs
    let r2 ← drawR2 (inputs.getD 0 0)
    let q ← orZero (uniform (-150) 150)
    let (f, e) := rfExcludedPairFloat r2 q (crfOf inputs)
    rows := rows.push ([r2, q] ++ inputs ++ [f, e])
  return rows.toList

/-- A charge in `[-2, 2]` e and `four_pi_eps_i` (138.9354, 1, or uniform in `[1, 200]`). -/
def rfSelfTermRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let inputs ← drawCrfInputs
    let q ← orZero (uniform (-2) 2)
    let fpe ← drawFrom [138.9354, 1] 1 200
    rows := rows.push ([q, fpe] ++ inputs ++ [rfSelfTermFloat q fpe (crfOf inputs)])
  return rows.toList

/-- The four row sets of `vectors/nonbonded.json`, from a fixed seed. -/
def nonbondedRows : List (List Float) × List (List Float) × List (List Float) ×
    List (List Float) :=
  (do return (← crfParametersRows, ← ljCrfPairRows, ← rfExcludedPairRows, ← rfSelfTermRows) :
    Gen _).run' 0x6e6577656e2d6e62

def nonbondedJson (crf lj excl selfRows : List (List Float)) : String :=
  let part (name : String) (rows : List (List Float)) : String :=
    s!"  \"{name}\": [\n" ++ ",\n".intercalate (rows.map hexRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"nonbonded.crf_parameters, nonbonded.lj_crf_pair, nonbonded.rf_excluded_pair, " ++
  "nonbonded.rf_self_term\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Nonbonded.crfParametersFloat, " ++
  "ljCrfPairFloat, rfExcludedPairFloat, rfSelfTermFloat; SplitMix64, fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits)\",\n" ++
  "  \"layout_crf_parameters\": \"cutoff, epsilon, rf_epsilon, rf_kappa, crf_cut, crf_2cut3i, " ++
  "crf_cut3i, cutoff_sq\",\n" ++
  "  \"layout_lj_crf_pair\": \"r2, c6, c12, q_prod, cutoff, epsilon, rf_epsilon, rf_kappa, " ++
  "f_mag, e_lj, e_crf; the kernel reads CRFParameters::new(cutoff, epsilon, rf_epsilon, " ++
  "rf_kappa)\",\n" ++
  "  \"layout_rf_excluded_pair\": \"r2, q_prod, cutoff, epsilon, rf_epsilon, rf_kappa, f_mag, " ++
  "e\",\n" ++
  "  \"layout_rf_self_term\": \"q, four_pi_eps_i, cutoff, epsilon, rf_epsilon, rf_kappa, e\",\n" ++
  part "crf_parameters" crf ++ ",\n" ++ part "lj_crf_pair" lj ++ ",\n" ++
  part "rf_excluded_pair" excl ++ ",\n" ++ part "rf_self_term" selfRows ++ "\n}\n"

def main : IO Unit := do
  IO.FS.createDirAll "vectors"
  IO.FS.writeFile "vectors/schedule.json" scheduleJson
  IO.println s!"wrote vectors/schedule.json ({grid} x {grid})"
  IO.FS.writeFile "vectors/dof.json" dofJson
  IO.println s!"wrote vectors/dof.json ({(dofGrid + 1) ^ 3} rows x {ndfmins.length})"
  IO.FS.writeFile "vectors/bath.json" bathJson
  IO.println s!"wrote vectors/bath.json ({bathSystems.length * 6} rows)"
  -- the literals the models use must be the exact doubles 0, 1, 2, 6, 12 (and 0.5)
  unless (0 : Float).toBits == 0 && (1 : Float).toBits == 0x3FF0000000000000 &&
      (2 : Float).toBits == 0x4000000000000000 && (6 : Float).toBits == 0x4018000000000000 &&
      (12 : Float).toBits == 0x4028000000000000 && halfFloat == 1 / 2 do
    throw (IO.userError "Float literals 0, 1, 2, 6, 12 are not exact")
  let (berendsen, noseHoover) := thermostatRows
  IO.FS.writeFile "vectors/thermostat.json" (thermostatJson berendsen noseHoover)
  IO.println (s!"wrote vectors/thermostat.json ({berendsen.length} Berendsen rows, " ++
    s!"{noseHoover.length} Nosé–Hoover rows)")
  let (crf, lj, excl, selfRows) := nonbondedRows
  IO.FS.writeFile "vectors/nonbonded.json" (nonbondedJson crf lj excl selfRows)
  IO.println (s!"wrote vectors/nonbonded.json ({crf.length} crf_parameters, {lj.length} " ++
    s!"lj_crf_pair, {excl.length} rf_excluded_pair, {selfRows.length} rf_self_term rows)")
