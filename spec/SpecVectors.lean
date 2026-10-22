import NewenSpec.Schedule
import NewenSpec.Thermo.Dof
import NewenSpec.Thermo.Bath
import NewenSpec.Thermostat
import NewenSpec.Nonbonded
import NewenSpec.Bonded

/-! `lake exe spec-vectors`: evaluate the Lean models on the grids the Rust tests sweep (and the
thermostat, nonbonded and bonded models on random `Float` inputs) and write them to `vectors/`, so
`cargo test` can compare the Rust functions against Lean without Lean. -/

open NewenSpec.Schedule
open NewenSpec.Thermo.Dof
open NewenSpec.Thermo.Bath
open NewenSpec.Thermostat
open NewenSpec.Nonbonded
open NewenSpec.Bonded (Vec3 tinyFloat quarterFloat piFloat quarticBondFloat harmonicBondFloat
  cosineAngleFloat harmonicAngleFloat dihedralChebyshevFloat dihedralNewFloat improperFloat)

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

/-! ## Bonded: the bonded kernels on random `Float` inputs, written as bits

Three row sets share their inputs between the kernels of one geometry: `bonds` (quartic and
harmonic bond), `angles` (cosine-harmonic and harmonic angle), `dihedrals` (Chebyshev and NTBDN = 0
dihedral, improper dihedral). Each set has random rows first, then rows at the edges: zero and
guard lengths, collinear and planar geometries, `cos φ = ±1`, `|sin θ| < 1e-10`, a linear bond angle
inside a dihedral, ζ next to ζ0 ± π. -/

def vec (x y z : Float) : Vec3 Float := ⟨x, y, z⟩

def vecBits (v : Vec3 Float) : List Float := [v.x, v.y, v.z]

/-- A vector with components uniform in `[−a, a]`. -/
def drawVec (a : Float) : Gen (Vec3 Float) := do
  return vec (← uniform (-a) a) (← uniform (-a) a) (← uniform (-a) a)

/-- A vector of random direction and length `len`. -/
def drawVecOfLength (len : Float) : Gen (Vec3 Float) := do
  let d ← drawVec 1
  let n := Float.sqrt (d.x * d.x + d.y * d.y + d.z * d.z)
  return d * (len / n)

/-- A bond-like vector: one time in four uniform in `[−0.5, 0.5]³`, else of length within 20 % of
`b0` (in `[0.09, 0.2]` nm when `b0 = 0`). -/
def drawBondVec (b0 : Float) : Gen (Vec3 Float) := do
  if (← below 4) == 0 then drawVec 0.5
  else
    let b0 ← if b0 == 0 then uniform 0.09 0.2 else pure b0
    drawVecOfLength (b0 * (← uniform 0.8 1.2))

/-- `degrees · π/180`. -/
def rad (deg : Float) : Float := deg * piFloat / 180

/-- One `bonds` row: `r, k_quartic, b0, k_harmonic`, then the outputs of both kernels. -/
def bondRow (r : Vec3 Float) (kq b0 kh : Float) : List Float :=
  let (eq, fq) := quarticBondFloat r kq b0
  let (eh, fh) := harmonicBondFloat r kh b0
  vecBits r ++ [kq, b0, kh] ++ [eq] ++ vecBits fq ++ [eh] ++ vecBits fh

def bondRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 1450 do
    let kq ← logUniform 1e5 2e7
    let b0 ← uniform 0.09 0.2
    let kh ← logUniform 1e5 5e5
    rows := rows.push (bondRow (← drawBondVec b0) kq b0 kh)
  let edges : List (Vec3 Float) :=
    [vec 0 0 0, vec 1e-11 0 0, vec 0 (-5e-11) 5e-11, vec 0 0 1e-10, vec tinyFloat 0 0,
      vec 0.1 0 0, vec 0 (-0.153) 0, vec 0.5 (-0.5) 0.5]
  for r in edges do
    for b0 in [0.1, 0.153] do
      rows := rows.push (bondRow r 1.57e7 b0 (2 * 1.57e7 * b0 * b0))
      rows := rows.push (bondRow r (← logUniform 1e5 2e7) b0 (← logUniform 1e5 5e5))
  return rows.toList

/-- One `angles` row: `r_ij, r_kj, k_cosine, cos0, k_harmonic, theta0`, then the outputs of both
kernels; `cos0 = cos theta0` for the random rows. -/
def angleRow (rij rkj : Vec3 Float) (kc cos0 kh theta0 : Float) : List Float :=
  let (ec, fic, fjc, fkc) := cosineAngleFloat rij rkj kc cos0
  let (eh, fih, fjh, fkh) := harmonicAngleFloat rij rkj kh theta0
  vecBits rij ++ vecBits rkj ++ [kc, cos0, kh, theta0] ++
    [ec] ++ vecBits fic ++ vecBits fjc ++ vecBits fkc ++
    [eh] ++ vecBits fih ++ vecBits fjh ++ vecBits fkh

/-- Force constants and a reference angle in `[90°, 125°]` (or `θ0`): `k_cosine, cos0, k_harmonic,
theta0`. -/
def drawAngleParams (theta0 : Option Float := none) : Gen (Float × Float × Float × Float) := do
  let theta0 ← match theta0 with
    | some t => pure t
    | none => do let d ← uniform 90 125; pure (rad d)
  return (← uniform 200 800, Float.cos theta0, ← uniform 100 800, theta0)

def angleRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 950 do
    let (kc, cos0, kh, theta0) ← drawAngleParams
    rows := rows.push (angleRow (← drawBondVec 0) (← drawBondVec 0) kc cos0 kh theta0)
  for i in List.range 50 do
    let rij ← drawBondVec 0
    let s ← uniform 0.5 2
    let a ← uniform 0.09 0.2
    -- collinear (θ = 0 and θ = π; along an axis cos θ = ±1 exactly, so |sin θ| < 1e-10),
    -- perpendicular, a zero bond, a bond at the guard
    let (rij, rkj) := match i % 8 with
      | 0 => (rij, rij * s)
      | 1 => (rij, rij * (-s))
      | 2 | 3 => (vec a 0 0, vec (-a * s) 0 0)
      | 4 => (vec 0 a 0, vec 0 (a * s) 0)
      | 5 => (rij, vec (-rij.y) rij.x 0)
      | 6 => (rij, vec 0 0 0)
      | _ => (rij, vec 0 0 tinyFloat)
    let (kc, cos0, kh, theta0) ← drawAngleParams (if i % 8 == 2 then some piFloat else none)
    rows := rows.push (angleRow rij rkj kc cos0 kh theta0)
  return rows.toList

/-- One `dihedrals` row: `r_ij, r_kj, r_kl, k, cos_delta, delta, m, k_improper, zeta0`, then the
outputs of the Chebyshev, the NTBDN = 0 and the improper kernel (`m` as a float). -/
def dihedralRow (rij rkj rkl : Vec3 Float) (k delta : Float) (m : Nat) (ki zeta0 : Float) :
    List Float :=
  let cosDelta := Float.cos delta
  let (ec, fic, fjc, fkc, flc) := dihedralChebyshevFloat rij rkj rkl k cosDelta m
  let (en, fin, fjn, fkn, fln) := dihedralNewFloat rij rkj rkl k delta m.toFloat
  let (ei, fii, fji, fki, fli) := improperFloat rij rkj rkl ki zeta0
  let out (e : Float) (fi fj fk fl : Vec3 Float) :=
    [e] ++ vecBits fi ++ vecBits fj ++ vecBits fk ++ vecBits fl
  vecBits rij ++ vecBits rkj ++ vecBits rkl ++ [k, cosDelta, delta, m.toFloat, ki, zeta0] ++
    out ec fic fjc fkc flc ++ out en fin fjn fkn fln ++ out ei fii fji fki fli

/-- `k, delta, m, k_improper, zeta0`: `δ ∈ {0, π}` two times in three, else uniform in
`[−π, π]`; `m` in `0..6` (and the unsupported 7 one time in sixteen); `ζ0 ∈ {0, 35.26°}` three
times in four, else uniform in `[−π, π]`. -/
def drawDihedralParams : Gen (Float × Float × Nat × Float × Float) := do
  let k ← uniform 0.1 50
  let delta ← match ← below 3 with
    | 0 => pure 0
    | 1 => pure piFloat
    | _ => uniform (-piFloat) piFloat
  let m ← if (← below 16) == 0 then pure 7 else below 7
  let zeta0 ← match ← below 4 with
    | 0 => pure 0
    | 1 | 2 => pure (rad 35.26)
    | _ => uniform (-piFloat) piFloat
  return (k, delta, m, ← uniform 100 400, zeta0)

def dihedralRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 650 do
    let (k, delta, m, ki, zeta0) ← drawDihedralParams
    rows := rows.push
      (dihedralRow (← drawBondVec 0) (← drawBondVec 0) (← drawBondVec 0) k delta m ki zeta0)
  for i in List.range 60 do
    let rij ← drawBondVec 0
    let rkj ← drawBondVec 0
    let rkl ← drawBondVec 0
    let a ← uniform (-1.5) 1.5
    let b ← uniform (-1.5) 1.5
    let s ← uniform 0.5 2
    let perp := (rij.cross rkj) * (1e-6 / Float.sqrt ((rij.cross rkj).dot (rij.cross rkj)))
    let (rij, rkj, rkl) := match i % 10 with
      | 0 => (rij, rkj, rij * a + rkj * b)          -- planar: cos φ = ±1
      | 1 => (rij, rkj, rij * (-s) + rkj * b)       -- planar, the other side
      | 2 => (rkj * s, rkj, rkl)                    -- i, j, k collinear
      | 3 => (rij, rkj, rkj * s)                    -- j, k, l collinear
      | 4 => (rkj * s + perp, rkj, rkl)             -- bond angle within 1e-5 rad of 0 or π
      | 5 => (rij, vec 0 0 0, rkl)                  -- r_kj = 0
      | 6 => (rij, rkj * 1e-6, rkl)                 -- |r_kj|² below the guard
      | 7 => (rij, rkj, rij * (-1 : Float) + rkj * b)  -- trans: ζ, φ at ±π
      | 8 => (rij, rkj, rij * (-1 : Float) + rkj * b + perp * (10 : Float))  -- just off trans
      | _ => (rij, rkj, rkl)
    let (k, delta, m, ki, zeta0) ← drawDihedralParams
    rows := rows.push (dihedralRow rij rkj rkl k delta m ki zeta0)
  return rows.toList

/-- The three row sets of `vectors/bonded.json`, from a fixed seed. -/
def bondedRows : List (List Float) × List (List Float) × List (List Float) :=
  (do return (← bondRows, ← angleRows, ← dihedralRows) : Gen _).run' 0x6e6577656e2d6264

def bondedJson (bonds angles dihedrals : List (List Float)) : String :=
  let part (name : String) (rows : List (List Float)) : String :=
    s!"  \"{name}\": [\n" ++ ",\n".intercalate (rows.map hexRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"bonded.quartic_bond, bonded.harmonic_bond, bonded.cosine_angle, " ++
  "bonded.harmonic_angle, bonded.dihedral_chebyshev, bonded.dihedral_new, " ++
  "bonded.improper_dihedral\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Bonded.quarticBondFloat, " ++
  "harmonicBondFloat, cosineAngleFloat, harmonicAngleFloat, dihedralChebyshevFloat, " ++
  "dihedralNewFloat, improperFloat; SplitMix64, fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits); m is a float\",\n" ++
  "  \"layout_bonds\": \"r (3), k_quartic, b0, k_harmonic; quartic_bond energy, f_i (3); " ++
  "harmonic_bond energy, f_i (3)\",\n" ++
  "  \"layout_angles\": \"r_ij (3), r_kj (3), k_cosine, cos0, k_harmonic, theta0; cosine_angle " ++
  "energy, f_i, f_j, f_k (3 each); harmonic_angle energy, f_i, f_j, f_k\",\n" ++
  "  \"layout_dihedrals\": \"r_ij (3), r_kj (3), r_kl (3), k, cos_delta, delta, m, k_improper, " ++
  "zeta0; dihedral_chebyshev(k, cos_delta, m) energy, f_i, f_j, f_k, f_l (3 each); " ++
  "dihedral_new(k, delta, m) energy, f_i, f_j, f_k, f_l; improper_dihedral(k_improper, zeta0) " ++
  "energy, f_i, f_j, f_k, f_l\",\n" ++
  part "bonds" bonds ++ ",\n" ++ part "angles" angles ++ ",\n" ++
  part "dihedrals" dihedrals ++ "\n}\n"

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
      (12 : Float).toBits == 0x4028000000000000 && NewenSpec.Nonbonded.halfFloat == 1 / 2 &&
      NewenSpec.Bonded.halfFloat == 1 / 2 && quarterFloat == 1 / 4 do
    throw (IO.userError "Float literals 0, 1, 2, 6, 12 are not exact")
  -- the other literals of the bonded models, and 1e-10 and π as the nearest doubles
  let lits : List Nat := [3, 4, 5, 8, 16, 18, 20, 32, 36, 48, 60, 80, 192, 180]
  unless lits.all (fun n => (Float.ofNat n).toUInt64.toNat == n) &&
      piFloat == 3.141592653589793 &&
      tinyFloat == 1e-10 do
    throw (IO.userError "bonded literals are not the expected doubles")
  let (berendsen, noseHoover) := thermostatRows
  IO.FS.writeFile "vectors/thermostat.json" (thermostatJson berendsen noseHoover)
  IO.println (s!"wrote vectors/thermostat.json ({berendsen.length} Berendsen rows, " ++
    s!"{noseHoover.length} Nosé–Hoover rows)")
  let (crf, lj, excl, selfRows) := nonbondedRows
  IO.FS.writeFile "vectors/nonbonded.json" (nonbondedJson crf lj excl selfRows)
  IO.println (s!"wrote vectors/nonbonded.json ({crf.length} crf_parameters, {lj.length} " ++
    s!"lj_crf_pair, {excl.length} rf_excluded_pair, {selfRows.length} rf_self_term rows)")
  let (bonds, angles, dihedrals) := bondedRows
  IO.FS.writeFile "vectors/bonded.json" (bondedJson bonds angles dihedrals)
  IO.println (s!"wrote vectors/bonded.json ({bonds.length} bonds, {angles.length} angles, " ++
    s!"{dihedrals.length} dihedrals rows)")
