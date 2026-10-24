import NewenSpec.Schedule
import NewenSpec.Thermo.Dof
import NewenSpec.Thermo.Bath
import NewenSpec.Topology.Pairs
import NewenSpec.Topology.Groups
import NewenSpec.Thermostat
import NewenSpec.Nonbonded
import NewenSpec.Bonded
import NewenSpec.Dynamics
import NewenSpec.Periodicity
import NewenSpec.SoftCore
import NewenSpec.PerturbedBonded
import NewenSpec.LatticeSum
import NewenSpec.P3M
import NewenSpec.SPME

/-! `lake exe spec-vectors`: evaluate the Lean models on the grids the Rust tests sweep (and the
thermostat, nonbonded, bonded, dynamics, periodicity, soft-core, perturbed bonded, lattice-sum and
Ewald models on random `Float` inputs) and write them to `vectors/`, so `cargo test` can compare the Rust functions against Lean without Lean. -/

open NewenSpec.Schedule
open NewenSpec.Thermo.Dof
open NewenSpec.Thermo.Bath
open NewenSpec.Topology.Pairs
open NewenSpec.Topology.Groups
open NewenSpec.Thermostat
open NewenSpec.Nonbonded
open NewenSpec.Bonded (tinyFloat quarterFloat piFloat quarticBondFloat harmonicBondFloat
  cosineAngleFloat harmonicAngleFloat dihedralChebyshevFloat dihedralNewFloat improperFloat)
open NewenSpec.Dynamics (V3 M3 Atom leapFrogVelocityFloat leapFrogPositionFloat kineticEnergyFloat
  groupKineticFloat kineticTemperatureFloat pressureTensorFloat scalarPressureFloat
  berendsenBarostatFloat thirdFloat)
open NewenSpec.Periodicity

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

/-! ## Output steps, the minimisation stop, the final step, COM removal (`schedule2.json`)

Grid bounds equal to the constants in `crates/newen-logic/tests/schedule_spec.rs`. -/

/-- `due`, `frames_in`: `0 ≤ every ≤ outputEvery`, `0 ≤ step, n_steps ≤ outputSteps`. -/
def outputEvery : Nat := 64
def outputSteps : Nat := 256
/-- `minimisation_stop`: `0 ≤ step < 16`, `0 ≤ nmin < 4`, `0 ≤ ntwx, ntwf, ntwe, ntwg < 4`. -/
def stopSteps : Nat := 16
def stopNmin : Nat := 4
def stopEvery : Nat := 4
/-- `final_step`: `0 ≤ n_steps, stopped_at ≤ finalSteps`. -/
def finalSteps : Nat := 20
/-- `com_removal`: `0 ≤ step ≤ comSteps`, `-1 ≤ nticom ≤ 3`, `-comNscm ≤ nscm ≤ comNscm`. -/
def comSteps : Nat := 60
def comNscm : Nat := 6

def bit (b : Bool) : Char := if b then '1' else '0'

def hexDigit (n : Nat) : Char := "0123456789abcdef".toList.getD n '?'

def quoted (s : String) : String := "    \"" ++ s ++ "\""

def intRow (xs : List Int) : String := "    [" ++ ", ".intercalate (xs.map toString) ++ "]"

def dueRows : List String := (List.range (outputEvery + 1)).map fun every =>
  quoted (String.ofList ((List.range (outputSteps + 1)).map fun s => bit (dueModel s every)))

def framesInRows : List String := (List.range (outputEvery + 1)).map fun every =>
  intRow ((List.range (outputSteps + 1)).map fun n => (framesInModel n every : Int))

/-- `-` for no stop, else the hex digit `positions + 2·forces + 4·energies + 8·free_energy`. -/
def stopChar : Option Stop → Char
  | none => '-'
  | some s => hexDigit (s.positions.toNat + 2 * s.forces.toNat + 4 * s.energies.toNat
      + 8 * s.freeEnergy.toNat)

def stopRows : List String :=
  (List.range stopSteps).flatMap fun step => (List.range stopNmin).flatMap fun nmin =>
    [false, true].map fun converged =>
      let ws := List.range stopEvery
      quoted (String.ofList (ws.flatMap fun x => ws.flatMap fun f => ws.flatMap fun e =>
        ws.map fun g => stopChar (minimisationStopModel step nmin converged x f e g)))

def finalRows : List String := (List.range (finalSteps + 1)).map fun n =>
  intRow ((finalStepModel n none : Int) ::
    (List.range (finalSteps + 1)).map fun k => (finalStepModel n (some k) : Int))

def comRows : List String :=
  ((List.range 5).map fun (i : Nat) => (i : Int) - 1).flatMap fun nticom =>
    ((List.range (2 * comNscm + 1)).map fun (i : Nat) => (i : Int) - comNscm).flatMap fun nscm =>
      [false, true].map fun periodic =>
        quoted (String.ofList ((List.range (comSteps + 1)).map fun s =>
          let p := comRemovalModel s nticom nscm periodic
          hexDigit (p.1.toNat + 2 * p.2.toNat)))

def jsonSection (name : String) (rows : List String) : String :=
  s!"  \"{name}\": [\n" ++ ",\n".intercalate rows ++ "\n  ]"

def schedule2Json : String :=
  "{\n" ++
  "  \"spec\": \"schedule.output.due, schedule.output.frames_in, " ++
  "schedule.output.minimisation_stop, schedule.output.final_step, schedule.com.com_removal\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Schedule.dueModel, framesInModel, " ++
  "minimisationStopModel, finalStepModel, comRemovalModel)\",\n" ++
  s!"  \"layout_due\": \"due[every][step] for 0 <= every <= {outputEvery}, " ++
  s!"0 <= step <= {outputSteps}; 1 = true\",\n" ++
  s!"  \"layout_frames_in\": \"frames_in[every][n_steps] for 0 <= every <= {outputEvery}, " ++
  s!"0 <= n_steps <= {outputSteps}\",\n" ++
  s!"  \"layout_minimisation_stop\": \"minimisation_stop[(step * {stopNmin} + nmin) * 2 + " ++
  s!"converged][((ntwx * {stopEvery} + ntwf) * {stopEvery} + ntwe) * {stopEvery} + ntwg] for " ++
  s!"0 <= step < {stopSteps}, 0 <= nmin < {stopNmin}, converged 0/1, 0 <= ntwx, ntwf, ntwe, " ++
  s!"ntwg < {stopEvery}; '-' = None, else the hex digit positions + 2 forces + 4 energies + " ++
  "8 free_energy\",\n" ++
  s!"  \"layout_final_step\": \"final_step[n_steps] = [None, Some(0), ..., " ++
  s!"Some({finalSteps})] for 0 <= n_steps <= {finalSteps}\",\n" ++
  s!"  \"layout_com_removal\": \"com_removal[((nticom + 1) * {2 * comNscm + 1} + nscm + " ++
  s!"{comNscm}) * 2 + periodic][step] for -1 <= nticom <= 3, -{comNscm} <= nscm <= {comNscm}, " ++
  s!"periodic 0/1, 0 <= step <= {comSteps}; digit translation + 2 rotation\",\n" ++
  jsonSection "due" dueRows ++ ",\n" ++ jsonSection "frames_in" framesInRows ++ ",\n" ++
  jsonSection "minimisation_stop" stopRows ++ ",\n" ++ jsonSection "final_step" finalRows ++
  ",\n" ++ jsonSection "com_removal" comRows ++ "\n}\n"

/-! ## Empty baths (`bath_empty.json`) -/

/-- Every list of half counts in `-2 ..= 2` of length `0 ..= 4`. -/
def halfLists : List (List Int) :=
  (List.range 5).flatMap fun len =>
    (List.range len).foldl (fun acc _ => acc.flatMap fun l =>
      ((List.range 5).map fun (i : Nat) => (i : Int) - 2).map (l ++ [·])) [[]]

def bathEmptyJson : String :=
  let rows := halfLists.map fun h => intRow (h ++ (emptyBathsModel h).map fun b => (b.toNat : Int))
  "{\n" ++
  "  \"spec\": \"thermo.bath.empty_baths\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Thermo.Bath.emptyBathsModel)\",\n" ++
  "  \"layout\": \"one list per row: the n half counts (every list over -2..=2 of length " ++
  "0..=4), then the n flags (1 = empty)\",\n" ++
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

/-! ## Topology: pair kinds and group lookups, flat integer rows -/

/-- Every ascending subset of `0 .. n − 1`, by bitmask. -/
def allSubsets (n : Nat) : List (List Nat) :=
  (List.range (2 ^ n)).map fun m => (List.range n).filter fun b => m / 2 ^ b % 2 = 1

/-- Every list of exactly `k` entries from `vals`, in any order, repeats allowed. -/
def tuples (vals : List Nat) : Nat → List (List Nat)
  | 0 => [[]]
  | k + 1 => (tuples vals k).flatMap fun t => vals.map (· :: t)

def ints (xs : List Nat) : List Int := xs.map Int.ofNat

/-- `PairKind` as the Rust enum's declaration order. -/
def kindCode : NewenTheory.Topology.Pairs.PairKind → Int
  | .normal => 0
  | .excluded => 1
  | .oneFour => 2
  | .atomPair => 3

/-- `LastAtomsCheck` as `(code, index)`: Ok 0, Empty 1, NotIncreasing 2 with its index,
LastNotFinal 3. -/
def checkCode : NewenSpec.Topology.Groups.LastAtomsCheck → List Int
  | .ok => [0, 0]
  | .empty => [1, 0]
  | .notIncreasing g => [2, g]
  | .lastNotFinal => [3, 0]

def optCode : Option Nat → Int
  | some g => g
  | none => -1

/-- Tag 0: one ascending list, then `contains_sorted` for `j = 0 ..= jMax`. -/
def containsRow (l : List Nat) (jMax : Nat) : List Int :=
  [(0 : Int), l.length] ++ ints l ++ (List.range (jMax + 1)).map fun j =>
    if containsSortedModel l j then 1 else 0

/-- Tag 1: three ascending lists, then `pair_kind` for `j = 0 ..= jMax`. -/
def pairRow (e f a : List Nat) (jMax : Nat) : List Int :=
  [(1 : Int), e.length] ++ ints e ++ [(f.length : Int)] ++ ints f ++ [(a.length : Int)] ++ ints a
    ++ (List.range (jMax + 1)).map fun j => kindCode (pairKindModel e f a j)

/-- Tag 2: `n_atoms` and a list of last atoms, then `check_last_atoms`, then `group_of` for
every atom `0 ..= n_atoms` (−1 for none). -/
def groupRow (n : Nat) (l : List Nat) : List Int :=
  [(2 : Int), n, l.length] ++ ints l ++ checkCode (checkLastAtomsModel l n)
    ++ (List.range (n + 1)).map fun a => optCode (groupOfModel l a)

/-- Tag 3: `eg_i, eg_j, n_groups`, then `energy_cell` (−1 for none). -/
def cellRow (i j n : Nat) : List Int := [3, i, j, n, optCode (energyCellModel i j n)]

/-- The `topology.json` rows. -/
def topologyRows : List (List Int) :=
  let subs8 := allSubsets 8
  let subs4 := allSubsets 4
  let pick (k : Nat) : List Nat := subs8.getD (k % subs8.length) []
  (subs8.map fun l => containsRow l 8)
    ++ (subs4.flatMap fun e => subs4.flatMap fun f => subs4.map fun a => pairRow e f a 4)
    ++ ((List.range subs8.length).map fun k => pairRow (pick k) (pick (7 * k + 3)) (pick (13 * k + 5)) 8)
    ++ ((List.range 9).flatMap fun n0 => (allSubsets 9).map fun l => groupRow (n0 + 1) l)
    ++ ((List.range 5).flatMap fun n0 => ((List.range 4).flatMap (tuples (List.range 5))).map
        fun l => groupRow (n0 + 1) l)
    ++ ((List.range 5).flatMap fun n => (List.range 6).flatMap fun i => (List.range 6).map
        fun j => cellRow i j n)

def topologyJson : String :=
  let rows := topologyRows.map fun r => "    [" ++ ", ".intercalate (r.map toString) ++ "]"
  "{\n" ++
  "  \"spec\": \"topology.pairs.contains_sorted, topology.pairs.pair_kind, " ++
  "topology.groups.group_of, topology.groups.check_last_atoms, topology.groups.energy_cell\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Topology.Pairs.containsSortedModel, " ++
  "pairKindModel; NewenSpec.Topology.Groups.groupOfModel, checkLastAtomsModel, " ++
  "energyCellModel)\",\n" ++
  "  \"layout\": \"one case per row, by its first entry. 0: len, the ascending list, then " ++
  "contains_sorted (0/1) for j = 0..=8. 1: three ascending lists (exclusions, one_four, " ++
  "atom_pair) each as len then entries, then pair_kind (0 Normal, 1 Excluded, 2 OneFour, " ++
  "3 AtomPair) for j = 0..=jmax, jmax = 4 or 8 (the rest of the row). 2: n_atoms, len, the " ++
  "last atoms, check_last_atoms as (code, index) with 0 Ok, 1 Empty, 2 NotIncreasing(index), " ++
  "3 LastNotFinal, then group_of (-1 = None) for atom = 0..=n_atoms. 3: eg_i, eg_j, n_groups, " ++
  "energy_cell (-1 = None)\",\n" ++
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

/-- A negative τ other than `−1`: `−0.5`, `−2` or uniform in `[−10, −1e-3)`. -/
def drawNegativeTau : Gen Float := do
  match ← below 3 with
  | 0 => return -0.5
  | 1 => return -2
  | _ => uniform (-10) (-1e-3)

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

/-- Rows with a negative τ other than `−1` (decoupled as well), drawn after all other rows. -/
def berendsenNegativeTauRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 40 do
    let dt ← pick thermoDts
    let dof ← drawDof
    let temperature ← drawTemperature
    let tau ← drawNegativeTau
    let ekin ← drawEkin dof
    rows := rows.push (berendsenRow ekin dof temperature tau dt)
  return rows.toList

def noseHooverTaus : List Float := [-1, 0.05, 0.1, 1]

/-- One chain run: random parameters (τ given, or drawn) and starting `ζ` in `[-50, 50]`, then
1–8 steps, each with a new kinetic energy and the `ζ` the previous step left; one row per step. -/
def noseHooverRun (tauGiven : Option (Gen Float) := none) : Gen (List String) := do
  let nhc := [1, 2, 3, 5].getD (← below 4) 1
  let dt ← pick thermoDts
  let dof ← drawDof
  let temperature ← drawTemperature
  let tau ← match tauGiven with
    | some drawTau => drawTau
    | none => match ← below (noseHooverTaus.length + 1) with
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

/-- Runs with a negative τ other than `−1` (ζ and the velocities left alone), drawn after all
other rows. -/
def noseHooverNegativeTauRows : Gen (List String) := do
  let mut rows := #[]
  for _ in List.range 10 do
    rows := rows ++ (← noseHooverRun (some drawNegativeTau)).toArray
  return rows.toList

/-- The rows of `vectors/thermostat.json`, from a fixed seed. -/
def thermostatRows : List (List Float) × List String :=
  (do
    let (b, n) := (← berendsenRows, ← noseHooverRows)
    return (b ++ (← berendsenNegativeTauRows), n ++ (← noseHooverNegativeTauRows)) :
    Gen _).run' 0x6e6577656e2d6d64

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

section Bonded
open NewenSpec.Bonded (Vec3)

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

/-- θ ≈ 0 with θ0 = π: the linear branch at its other end; drawn after all other rows. -/
def angleZeroRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for i in List.range 6 do
    let a ← uniform 0.09 0.2
    let s ← uniform 0.5 2
    let rij ← drawBondVec 0
    let (rij, rkj) := if i % 2 == 0 then (vec 0 0 a, vec 0 0 (a * s)) else (rij, rij * s)
    let (kc, cos0, kh, theta0) ← drawAngleParams (some piFloat)
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

/-- The improper angle ζ of a geometry before the wrap, as the kernel forms it. -/
def rawZeta (rij rkj rkl : Vec3 Float) : Float :=
  let mj := rij.cross rkj
  let nk := rkj.cross rkl
  let c := mj.dot nk / (Float.sqrt (mj.dot mj) * Float.sqrt (nk.dot nk))
  let z := Float.acos (if c < -1 then -1 else if 1 < c then 1 else c)
  if rij.dot nk < 0 then -z else z

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

/-- j, k, l within 1e-7 nm of a line with i off it (NEAR_LINEAR on the l side drops f_l), and
ζ0 beyond ±π; drawn after all other rows. -/
def dihedralEdgeRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 5 do
    let rij ← drawBondVec 0
    let rkj ← drawBondVec 0
    let s ← uniform 0.5 2
    let side := rkj.cross rij
    let side := side * (1e-7 / Float.sqrt (side.dot side))
    let (k, delta, m, ki, zeta0) ← drawDihedralParams
    rows := rows.push (dihedralRow rij rkj (rkj * s + side) k delta m ki zeta0)
  -- |ζ0| > π, |ζ − ζ0| in (π, 3π) (one turn of the wrap) or in (3π, 4π) (two turns)
  for i in List.range 8 do
    let rij ← drawBondVec 0
    let rkj ← drawBondVec 0
    let rkl ← drawBondVec 0
    let (k, delta, m, ki, _) ← drawDihedralParams
    let d ← if i % 2 == 0 then uniform (piFloat + 0.1) (3 * piFloat)
      else uniform (3 * piFloat + 0.1) (4 * piFloat - 0.1)
    let z := rawZeta rij rkj rkl
    let zeta0 := if z < 0 then z - d else z + d
    rows := rows.push (dihedralRow rij rkj rkl k delta m ki zeta0)
  return rows.toList

/-- The three row sets of `vectors/bonded.json`, from a fixed seed. -/
def bondedRows : List (List Float) × List (List Float) × List (List Float) :=
  (do
    let (b, a, d) := (← bondRows, ← angleRows, ← dihedralRows)
    return (b, a ++ (← angleZeroRows), d ++ (← dihedralEdgeRows)) :
    Gen _).run' 0x6e6577656e2d6264

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

end Bonded

/-! ## Dynamics: leap-frog, kinetic energy, temperature, pressure, barostat -/

def v3 (a : V3 Float) : List Float := [a.x, a.y, a.z]

def m3 (p : M3 Float) : List Float := v3 p.c0 ++ v3 p.c1 ++ v3 p.c2

def drawV3 (lo hi : Float) : Gen (V3 Float) := do
  return ⟨← uniform lo hi, ← uniform lo hi, ← uniform lo hi⟩

def dynDts : List Float := [5e-4, 1e-3, 2e-3, 4e-3]
def atomMasses : List Float := [1.008, 12.011, 14.0067, 15.9994, 32.06, 39.948]

/-- A mass in u: one of `atomMasses`, uniform in `[0.5, 200]`, or `0` one time in 16. -/
def drawMass : Gen Float := do
  if (← below 16) == 0 then return 0
  drawFrom atomMasses 0.5 200

/-- A velocity: uniform in `[-3, 3]` nm/ps per component, or zero one time in ten. -/
def drawVelocity : Gen (V3 Float) := do
  if (← below 10) == 0 then return ⟨0, 0, 0⟩
  drawV3 (-3) 3

/-- A force: uniform in `[-2000, 2000]` kJ mol⁻¹ nm⁻¹ per component, up to `1e6`, or zero. -/
def drawForce : Gen (V3 Float) := do
  match ← below 8 with
  | 0 => return ⟨0, 0, 0⟩
  | 1 => drawV3 (-1e6) 1e6
  | _ => drawV3 (-2000) 2000

/-- One atom per row: `v, f, m⁻¹, dt, x`, then the new velocity, then the position step from the
new velocity. The inverse mass is `1/m` as `Topology::compute_inverse_masses` forms it (`0` for
a massless atom). -/
def leapFrogRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 1200 do
    let m ← drawMass
    let im := if m > 0 then 1 / m else 0
    let dt ← drawFrom dynDts 1e-4 5e-3
    let v ← drawVelocity
    let f ← drawForce
    let x ← drawV3 (-10) 10
    let v' := leapFrogVelocityFloat v f im dt
    let x' := leapFrogPositionFloat x v' dt
    rows := rows.push (v3 v ++ v3 f ++ [im, dt] ++ v3 x ++ v3 v' ++ v3 x')
  return rows.toList

/-- One group of 1–8 atoms per row. The old velocity is the new one less a leap-frog kick, or
independent; one group in 40 has only massless atoms (a centre of mass of `0/0`). -/
def kineticRows : Gen (List String) := do
  let mut rows := #[]
  for _ in List.range 500 do
    let n := 1 + (← below 8)
    let massless := (← below 40) == 0
    let mut atoms := #[]
    for _ in List.range n do
      let m ← if massless then pure 0 else drawMass
      let vNew ← drawVelocity
      let vOld ← if (← below 2) == 0 then drawVelocity else do
        let dt ← pick dynDts
        let kick ← drawV3 (-2000) 2000
        let im := if m > 0 then 1 / m else 0
        pure ⟨vNew.x - kick.x * im * dt, vNew.y - kick.y * im * dt, vNew.z - kick.z * im * dt⟩
      atoms := atoms.push (⟨m, vNew, vOld⟩ : Atom Float)
    let (eKin, eKinNew) := kineticEnergyFloat atoms.toList
    let (com, total, newCom, newTotal) := groupKineticFloat atoms.toList
    let fields := atoms.toList.flatMap fun a => [a.m] ++ v3 a.vNew ++ v3 a.vOld
    rows := rows.push ("    [" ++ ", ".intercalate ([toString n] ++
      (fields ++ [eKin, eKinNew, com, total, newCom, newTotal]).map hexBits) ++ "]")
  return rows.toList

/-- `ekin, dof, k_B, T`: `dof` among the GROMOS-sized counts, `0` and a negative one, `k_B` the
GROMOS value, the CODATA value or near them. -/
def temperatureRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 800 do
    let ekin ← orZero (logUniform 1e-6 1e7)
    let dof ← drawFrom [0, 1, 3, 5.5, 12.76, 100, 6000, 1e5, -3] 1 1e5
    let kB ← drawFrom [kBFloat, 0.00831446261815324] 0.008 0.0085
    rows := rows.push [ekin, dof, kB, kineticTemperatureFloat ekin dof kB]
  return rows.toList

/-- `K, W, V` (tensors by columns), the pressure tensor, its scalar pressure. `K` has a
non-negative diagonal; `V` is log-uniform in `[1, 1e4]` nm³, or `0` one time in 40. -/
def pressureRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 600 do
    let k ← (List.range 9).mapM fun i =>
      if i % 4 == 0 then uniform 0 1e4 else uniform (-1e3) 1e3
    let w ← (List.range 9).mapM fun _ => uniform (-1e5) 1e5
    let col (l : List Float) (j : Nat) : V3 Float :=
      ⟨l.getD (3 * j) 0, l.getD (3 * j + 1) 0, l.getD (3 * j + 2) 0⟩
    let K : M3 Float := ⟨col k 0, col k 1, col k 2⟩
    let W : M3 Float := ⟨col w 0, col w 1, col w 2⟩
    let V ← if (← below 40) == 0 then pure 0 else logUniform 1 1e4
    let P := pressureTensorFloat K W V
    rows := rows.push (m3 K ++ m3 W ++ [V] ++ m3 P ++ [scalarPressureFloat P])
  return rows.toList

/-- `P, P₀, β, τ, dt, μ` in GROMOS units (kJ mol⁻¹ nm⁻³; 1 bar = 0.06102): one row in ten at
`P = P₀`, one in twenty with `|P|` up to `1e5` (a negative radicand, `μ` NaN), one in thirty
with `τ = 0`. In GROMOS units the coupling term `β dt/τ (P₀ − P)` is far below the last bit of
`1`, so the result hides its rounding; a quarter of the rows are stressed instead, with
`β ∈ [0.1, 1]`, `τ ∈ [dt/2, 4 dt]` and `|P|` up to `2`, a term of order `1`. -/
def barostatRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 1500 do
    let dt ← drawFrom dynDts 1e-4 5e-3
    let p0 ← drawFrom [0.06102, 0, 1] (-1) 1
    let stressed := (← below 4) == 0
    let p ← if stressed then uniform (-2) 2 else match ← below 20 with
      | 0 | 1 => pure p0
      | 2 => uniform (-1e5) 1e5
      | _ => drawFrom [0, 0.06102, 1, -1] (-50) 50
    let beta ← if stressed then uniform 0.1 1 else drawFrom [4.575e-4, 7.513e-4] 1e-5 1e-3
    let tau ← if stressed then pure (dt * (← uniform 0.5 4))
      else if (← below 30) == 0 then pure 0 else drawFrom [0.1, 0.5, 1, 2] 0.05 5
    rows := rows.push [p, p0, beta, tau, dt, berendsenBarostatFloat p p0 beta tau dt]
  return rows.toList

/-- The row sets of `vectors/dynamics.json`, from a fixed seed. -/
def dynamicsRows : List (List Float) × List String × List (List Float) × List (List Float) ×
    List (List Float) :=
  (do return (← leapFrogRows, ← kineticRows, ← temperatureRows, ← pressureRows,
    ← barostatRows) : Gen _).run' 0x6e6577656e2d6479

def dynamicsJson (leap : List (List Float)) (kin : List String)
    (temp pres baro : List (List Float)) : String :=
  let part (name : String) (rows : List String) : String :=
    s!"  \"{name}\": [\n" ++ ",\n".intercalate rows ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"dynamics.leap_frog_velocity, dynamics.leap_frog_position, " ++
  "dynamics.kinetic_energy, dynamics.group_kinetic_energy, dynamics.kinetic_temperature, " ++
  "dynamics.pressure_tensor, dynamics.scalar_pressure, dynamics.berendsen_barostat_scale\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Dynamics.*Float; SplitMix64, " ++
  "fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits); n is a plain integer\",\n" ++
  "  \"layout_leap_frog\": \"v(3), f(3), inverse_mass, dt, x(3), v_new(3), x_new(3); x_new is " ++
  "leap_frog_position(x, v_new, dt)\",\n" ++
  "  \"layout_kinetic\": \"n, then (mass, v_new(3), v_old(3)) per atom, then kinetic_energy " ++
  "(average, new) and group_kinetic_energy (com, total, new_com, new_total)\",\n" ++
  "  \"layout_temperature\": \"ekin, dof, k_boltzmann, temperature\",\n" ++
  "  \"layout_pressure\": \"kinetic(9), virial(9), volume, pressure(9), scalar_pressure; " ++
  "tensors by columns (x_axis, y_axis, z_axis), each x, y, z\",\n" ++
  "  \"layout_barostat\": \"pressure, pressure0, compressibility, tau, dt, mu\",\n" ++
  part "leap_frog" (leap.map hexRow) ++ ",\n" ++ part "kinetic" kin ++ ",\n" ++
  part "temperature" (temp.map hexRow) ++ ",\n" ++ part "pressure" (pres.map hexRow) ++ ",\n" ++
  part "barostat" (baro.map hexRow) ++ "\n}\n"

/-! ## Periodicity: nearest image and put-into-box on random `Float` inputs, written as bits -/

section Periodicity
open NewenTheory.Space.Periodicity (Vec3)

def vecList (v : Vec3 Float) : List Float := [v.x, v.y, v.z]

def matList (m : Mat3 Float) : List Float := vecList m.xAxis ++ vecList m.yAxis ++ vecList m.zAxis

/-- The neighbours of a positive finite `x` one ulp up and down. -/
def ulpUp (x : Float) : Float := Float.ofBits (x.toBits + 1)
def ulpDown (x : Float) : Float := Float.ofBits (x.toBits - 1)

def nanFloat : Float := 0 / 0
def infFloat : Float := 1 / 0

/-- Box edges: whole, binary-exact and not. -/
def rectEdges : List Float := [1, 2, 3, 4, 2.5, 3.1, 1.8653, 3.767055681, 6.2]

/-- A rectangular box: cubic from `rectEdges`, cubic random, or each edge random in `[1.5, 8]`. -/
def drawRectBox : Gen (Vec3 Float) := do
  match ← below 3 with
  | 0 => let l ← pick rectEdges; return ⟨l, l, l⟩
  | 1 => let l ← uniform 1.5 8; return ⟨l, l, l⟩
  | _ => return ⟨← uniform 1.5 8, ← uniform 1.5 8, ← uniform 1.5 8⟩

/-- One component of a separation along an edge `l`: inside, a few boxes out, fifty boxes out,
an odd multiple of `l/2` (a tie), one ulp either side of `l/2`, a whole multiple of `l`, or zero
of either sign. -/
def drawSep (l : Float) : Gen Float := do
  let half := pointFiveFloat * l
  let sign : Float := if (← below 2) == 0 then 1 else -1
  match ← below 10 with
  | 0 => uniform (-half) half
  | 1 => uniform (-l) l
  | 2 => uniform (-3 * l) (3 * l)
  | 3 => uniform (-50 * l) (50 * l)
  | 4 => return sign * half
  | 5 => return sign * half * [1, 3, 5, 7].getD (← below 4) 1
  | 6 => return sign * (if (← below 2) == 0 then ulpUp half else ulpDown half)
  | 7 => return ((← below 7).toFloat - 3) * l
  | 8 => return if (← below 2) == 0 then 0 else -0.0
  | _ => uniform (-2 * l) (2 * l)

/-- `ri, rj` with the separation `sep`: `(sep, 0)`, `(0, −sep)`, or `rj` in the box and
`ri = rj + sep` (rounded). -/
def splitSep (sep inBox : Vec3 Float) : Gen (Vec3 Float × Vec3 Float) := do
  match ← below 3 with
  | 0 => return (sep, ⟨0, 0, 0⟩)
  | 1 => return (⟨0, 0, 0⟩, ⟨-sep.x, -sep.y, -sep.z⟩)
  | _ => return (inBox + sep, inBox)

/-- One vector in thirty gets a NaN or an infinite component (`inf = false`: NaN only). -/
def spoil (inf : Bool) (v : Vec3 Float) : Gen (Vec3 Float) := do
  if (← below 30) != 0 then return v
  let bad ← if inf then pick [nanFloat, infFloat, -infFloat] else pure nanFloat
  match ← below 3 with
  | 0 => return { v with x := bad }
  | 1 => return { v with y := bad }
  | _ => return { v with z := bad }

def vacuumRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 200 do
    let ri ← spoil true ⟨← drawSep 3, ← drawSep 3, ← drawSep 3⟩
    let rj ← spoil true ⟨← uniform (-5) 5, ← uniform (-5) 5, ← uniform (-5) 5⟩
    rows := rows.push (vecList ri ++ vecList rj ++ vecList (vacuumNearestImageFloat ri rj))
  return rows.toList

def rectangularNewRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for l in rectEdges ++ [0.1, 1e-3, 1e3, 0.7] do
    let (h, inv) := rectangularNewFloat ⟨l, l, l⟩
    rows := rows.push ([l, l, l] ++ vecList h ++ vecList inv)
  for _ in List.range 107 do
    let b ← drawRectBox
    let (h, inv) := rectangularNewFloat b
    rows := rows.push (vecList b ++ vecList h ++ vecList inv)
  return rows.toList

def rectangularNearestImageRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 1200 do
    let b ← drawRectBox
    let sep : Vec3 Float := ⟨← drawSep b.x, ← drawSep b.y, ← drawSep b.z⟩
    let inBox : Vec3 Float := ⟨← uniform 0 b.x, ← uniform 0 b.y, ← uniform 0 b.z⟩
    let (ri, rj) ← splitSep sep inBox
    let ri ← spoil true ri
    rows := rows.push (vecList b ++ vecList ri ++ vecList rj ++
      vecList (rectangularNearestImageFloat b ri rj))
  return rows.toList

/-- A position along an edge `l`: inside, around, far, a whole multiple of `l`, one ulp below
`l`, or zero of either sign. -/
def drawPos (l : Float) : Gen Float := do
  match ← below 8 with
  | 0 => uniform 0 l
  | 1 => uniform (-l) (2 * l)
  | 2 => uniform (-50 * l) (50 * l)
  | 3 => return ((← below 7).toFloat - 3) * l
  | 4 => return ulpDown l
  | 5 => return if (← below 2) == 0 then 0 else -0.0
  | 6 => uniform (-3 * l) 0
  | _ => uniform (-l) l

def rectangularPutIntoBoxRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 600 do
    let b ← drawRectBox
    let pos ← spoil true ⟨← drawPos b.x, ← drawPos b.y, ← drawPos b.z⟩
    rows := rows.push (vecList b ++ vecList pos ++ vecList (rectangularPutIntoBoxFloat b pos))
  return rows.toList

/-- `truncoct_triclinic_box(diag(l), true)`, operation by operation: the truncated octahedron of
cube edge `l` as a lower-triangular triclinic box. -/
def truncoctBoxFloat (l : Float) : Mat3 Float :=
  let d := pointFiveFloat * Float.sqrt 3 * l
  ⟨⟨d, 0, 0⟩, ⟨d / 3, 2 * Float.sqrt 2 / 3 * d, 0⟩,
    ⟨-d / 3, Float.sqrt 2 / 3 * d, Float.sqrt 6 / 3 * d⟩⟩

def truncoctEdges : List Float := [3.767055681, 4, 5.5, 7.4]

/-- A triclinic box: rectangular, a truncated octahedron, a reduced triclinic box
(`|bₓ|, |cₓ| ≤ aₓ/2`, `|c_y| ≤ b_y/2`), a strongly skewed one, or a monoclinic one. -/
def drawTriclinicBox : Gen (Mat3 Float) := do
  match ← below 5 with
  | 0 => return ⟨⟨← uniform 2 8, 0, 0⟩, ⟨0, ← uniform 2 8, 0⟩, ⟨0, 0, ← uniform 2 8⟩⟩
  | 1 => return truncoctBoxFloat (← drawFrom truncoctEdges 3 8)
  | 2 =>
    let ax ← uniform 2 8
    let by_ ← uniform 2 8
    return ⟨⟨ax, 0, 0⟩, ⟨← uniform (-ax / 2) (ax / 2), by_, 0⟩,
      ⟨← uniform (-ax / 2) (ax / 2), ← uniform (-by_ / 2) (by_ / 2), ← uniform 2 8⟩⟩
  | 3 =>
    let ax ← uniform 2 8
    let by_ ← uniform 2 8
    return ⟨⟨ax, 0, 0⟩, ⟨← uniform (-ax) ax, by_, 0⟩,
      ⟨← uniform (-ax) ax, ← uniform (-by_) by_, ← uniform 2 8⟩⟩
  | _ => return ⟨⟨← uniform 2 8, 0, 0⟩, ⟨0, ← uniform 2 8, 0⟩, ⟨← uniform (-2) 2, 0, ← uniform 2 8⟩⟩

/-- A fractional coordinate: inside the cell, a few cells out, forty cells out, a half-integer
(a face of the cell, as far as rounding lets it be), an integer, or zero. -/
def drawFracCoord : Gen Float := do
  match ← below 8 with
  | 0 => uniform (-0.5) 0.5
  | 1 => uniform (-1.5) 1.5
  | 2 => uniform (-4) 4
  | 3 => uniform (-40) 40
  | 4 => return [1, 3, 5, 9].getD (← below 4) 1 * (if (← below 2) == 0 then 0.5 else -0.5)
  | 5 => return (← below 7).toFloat - 3
  | 6 => return 0
  | _ => uniform (-0.6) 0.6

/-- The point with fractional coordinates `s` in the box `m` (glam's matrix product). -/
def pointFloat (m : Mat3 Float) (s : Vec3 Float) : Vec3 Float := mulVec Float m s

/-- Half a box vector or a sum of halves, exactly: a separation on a face of the cell. -/
def drawFace (m : Mat3 Float) : Gen (Vec3 Float) := do
  let h (v : Vec3 Float) : Vec3 Float := scale Float v pointFiveFloat
  let n (v : Vec3 Float) : Vec3 Float := ⟨-v.x, -v.y, -v.z⟩
  let faces := [h m.zAxis, n (h m.zAxis), h m.yAxis, n (h m.yAxis), h m.xAxis, n (h m.xAxis),
    h m.xAxis + h m.yAxis, h m.xAxis + h m.yAxis + h m.zAxis]
  return faces.getD (← below faces.length) ⟨0, 0, 0⟩

def triclinicNewRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 180 do
    let m ← drawTriclinicBox
    rows := rows.push (matList m ++ matList (triclinicNewFloat m))
  for _ in List.range 19 do
    let m : Mat3 Float := ⟨⟨← uniform (-3) 3, ← uniform (-3) 3, ← uniform (-3) 3⟩,
      ⟨← uniform (-3) 3, ← uniform (-3) 3, ← uniform (-3) 3⟩,
      ⟨← uniform (-3) 3, ← uniform (-3) 3, ← uniform (-3) 3⟩⟩
    rows := rows.push (matList m ++ matList (triclinicNewFloat m))
  let singular : Mat3 Float := ⟨⟨1, 2, 3⟩, ⟨2, 4, 6⟩, ⟨0, 1, 1⟩⟩
  rows := rows.push (matList singular ++ matList (triclinicNewFloat singular))
  return rows.toList

/-- The triclinic nearest image at `fuelFloat` turns per loop, refused if half the fuel gives other
bits, that is if some loop came near its bound. -/
def triclinicNearestImageChecked (m : Mat3 Float) (ri rj : Vec3 Float) :
    Except String (Vec3 Float) :=
  let full := triclinicNearestImageFloat fuelFloat m ri rj
  let half := triclinicNearestImageFloat (fuelFloat / 2) m ri rj
  if (vecList full).map Float.toBits == (vecList half).map Float.toBits ||
      (vecList full).any Float.isNaN then .ok full
  else .error s!"triclinic loop near its fuel bound at {vecList ri} {vecList rj}"

def triclinicNearestImageRows : Gen (Except String (List (List Float))) := do
  let mut rows := #[]
  for _ in List.range 1500 do
    let m ← drawTriclinicBox
    let sep ← if (← below 8) == 0 then drawFace m
      else pure (pointFloat m ⟨← drawFracCoord, ← drawFracCoord, ← drawFracCoord⟩)
    let inBox := pointFloat m ⟨← unit, ← unit, ← unit⟩
    let (ri, rj) ← splitSep sep inBox
    let ri ← spoil false ri
    match triclinicNearestImageChecked m ri rj with
    | .ok out => rows := rows.push (matList m ++ vecList ri ++ vecList rj ++ vecList out)
    | .error e => return .error e
  return .ok rows.toList

def triclinicPutIntoBoxRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 600 do
    let m ← drawTriclinicBox
    let pos ← if (← below 6) == 0 then drawFace m
      else pure (pointFloat m ⟨← drawFracCoord, ← drawFracCoord, ← drawFracCoord⟩)
    let pos ← spoil true pos
    rows := rows.push (matList m ++ vecList pos ++ vecList (triclinicPutIntoBoxFloat m pos))
  return rows.toList

def truncoctRows : List (List Float) :=
  (truncoctEdges ++ [3, 6.123, 8]).map fun l => [l] ++ matList (truncoctBoxFloat l)

/-- The row sets of `vectors/periodicity.json`, from a fixed seed. -/
def periodicityRows : Except String (List (String × List (List Float))) :=
  (do
    let vac ← vacuumRows
    let rnew ← rectangularNewRows
    let rnim ← rectangularNearestImageRows
    let rput ← rectangularPutIntoBoxRows
    let tnew ← triclinicNewRows
    let tnim ← triclinicNearestImageRows
    let tput ← triclinicPutIntoBoxRows
    return tnim.map fun tnim =>
      [("vacuum_nearest_image", vac), ("rectangular_new", rnew),
        ("rectangular_nearest_image", rnim), ("rectangular_put_into_box", rput),
        ("truncoct_box", truncoctRows), ("triclinic_new", tnew),
        ("triclinic_nearest_image", tnim), ("triclinic_put_into_box", tput)] :
    Gen _).run' 0x6e6577656e2d7062

def periodicityJson (parts : List (String × List (List Float))) : String :=
  let part (p : String × List (List Float)) : String :=
    s!"  \"{p.1}\": [\n" ++ ",\n".intercalate (p.2.map hexRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"periodicity.vacuum_nearest_image, periodicity.rectangular_new, " ++
  "periodicity.rectangular_nearest_image, periodicity.rectangular_put_into_box, " ++
  "periodicity.triclinic_new, periodicity.triclinic_nearest_image, " ++
  "periodicity.triclinic_put_into_box\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.Periodicity, the *Float definitions; " ++
  "SplitMix64, fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits); a vector is x, y, z; a matrix is its columns a, b, c\",\n" ++
  "  \"layout_vacuum_nearest_image\": \"ri, rj, nearest_image\",\n" ++
  "  \"layout_rectangular_new\": \"box_size, half_box, inv_box\",\n" ++
  "  \"layout_rectangular_nearest_image\": \"box_size, ri, rj, nearest_image\",\n" ++
  "  \"layout_rectangular_put_into_box\": \"box_size, pos, put_into_box\",\n" ++
  "  \"layout_truncoct_box\": \"cube edge L, truncoct_triclinic_box(diag(L), true): the " ++
  "truncated octahedra among the triclinic boxes\",\n" ++
  "  \"layout_triclinic_new\": \"box_matrix, inv_box_matrix\",\n" ++
  "  \"layout_triclinic_nearest_image\": \"box_matrix, ri, rj, nearest_image\",\n" ++
  "  \"layout_triclinic_put_into_box\": \"box_matrix, pos, put_into_box\",\n" ++
  ",\n".intercalate (parts.map part) ++ "\n}\n"

end Periodicity

/-! ## Soft-core: the perturbed nonbonded kernels on random `Float` inputs, written as bits -/

def scLambdas : List Float := [0, 0.25, 0.5, 0.75, 1]
def scAlphas : List Float := [0, 0.5, 1.0, 1.51]

/-- A λ: one of `scLambdas` or uniform in `[0, 1)`. -/
def drawLam : Gen Float := drawFrom scLambdas 0 1

/-- An individual λ: one of `scLambdas` or uniform in `[-0.1, 1.1)` (a λ schedule may leave
`[0, 1]`). -/
def drawIndividualLam : Gen Float := drawFrom scLambdas (-0.1) 1.1

/-- An α: one of `scAlphas` or uniform in `[0, 2)`. -/
def drawAlpha : Gen Float := drawFrom scAlphas 0 2

/-- The inputs of `PerturbedLambdaParams::from_lambda`: half the time one λ with derivatives 1
(no individual λ), else four λ and four derivatives in `[-2, 2)`; NLAM 1, 2, 3 or 5. -/
def drawTable : Gen (List Float × Nat) := do
  let n := [1, 2, 3, 5].getD (← below 4) 1
  if (← below 2) == 0 then
    let l ← drawLam
    return ([l, l, l, l, 1, 1, 1, 1], n)
  else
    let mut xs := #[]
    for _ in List.range 4 do
      xs := xs.push (← drawIndividualLam)
    for _ in List.range 4 do
      xs := xs.push (← uniform (-2) 2)
    return (xs.toList, n)

def tableOf (t : List Float × Nat) : NewenSpec.SoftCore.LambdaParams Float :=
  let x (i : Nat) := t.1.getD i 0
  NewenSpec.SoftCore.lambdaParamsFloat (x 0) (x 1) (x 2) (x 3) (x 4) (x 5) (x 6) (x 7)
    (Float.ofNat t.2) (Float.ofNat (t.2 - 1))

/-- The table inputs as row fields: eight floats as bits, NLAM as a plain integer. -/
def tableFields (t : List Float × Nat) : List String := t.1.map hexBits ++ [toString t.2]

/-- A `C6`: `0` one time in nine, tiny (`0 < |C6| ≤ 1e-10`, either sign) one time in six,
negative one time in eighteen, else log-uniform in `[1e-5, 1e-2]`. -/
def drawC6 : Gen Float := do
  match ← below 18 with
  | 0 | 1 => return 0
  | 2 => logUniform 1e-12 1e-10
  | 3 => return 1e-10
  | 4 => do return -(← logUniform 1e-12 1e-10)
  | 5 => do return -(← logUniform 1e-5 1e-2)
  | _ => logUniform 1e-5 1e-2

/-- One end state `(C6, C12, q)`: a dummy (all zero), uncharged, or charged. -/
def drawState : Gen (Float × Float × Float) := do
  match ← below 3 with
  | 0 => return (0, 0, 0)
  | k =>
    let c6 ← drawC6
    let c12 ← logUniform 1e-8 1e-5
    let q ← if k == 1 then pure 0 else uniform (-2) 2
    return (c6, c12, q)

def drawFpe : Gen Float := drawFrom [138.9354, 1] 1 200

/-- `perturbed_lj_crf_pair` rows. -/
def perturbedPairRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 800 do
    let inputs ← drawCrfInputs
    let r2 ← drawR2 (inputs.getD 0 0)
    let (aC6, aC12, aQ) ← drawState
    let (bC6, bC12, bQ) ← drawState
    let αlj ← drawAlpha
    let αc ← drawAlpha
    let fpe ← drawFpe
    let cs ← drawFrom [1, 0.5, 0.8333333333333334] 0 1
    let t ← drawTable
    let (f, eLj, eCrf, deLj, deCrf) := NewenSpec.SoftCore.perturbedLjCrfPairFloat r2 aC6 aC12 bC6
      bC12 aQ bQ αlj αc (tableOf t) (crfOf inputs) fpe cs
    rows := rows.push (([r2, aC6, aC12, bC6, bC12, aQ, bQ, αlj, αc] ++ inputs ++ [fpe, cs]).map
      hexBits ++ tableFields t ++ [f, eLj, eCrf, deLj, deCrf].map hexBits)
  -- overlapping atoms (r² below 1e-10) with a soft core in both states
  for _ in List.range 60 do
    let inputs ← drawCrfInputs
    let r2 ← if (← below 4) == 0 then pure 0 else logUniform 1e-16 9e-11
    let aC6 ← logUniform 1e-5 1e-2
    let aC12 ← logUniform 1e-8 1e-5
    let bC6 ← logUniform 1e-5 1e-2
    let bC12 ← logUniform 1e-8 1e-5
    let aQ ← uniform (-2) 2
    let bQ ← uniform (-2) 2
    let αlj ← uniform 0.1 2
    let αc ← uniform 0.1 2
    let fpe ← drawFpe
    let t ← drawTable
    let (f, eLj, eCrf, deLj, deCrf) := NewenSpec.SoftCore.perturbedLjCrfPairFloat r2 aC6 aC12 bC6
      bC12 aQ bQ αlj αc (tableOf t) (crfOf inputs) fpe 1
    rows := rows.push (([r2, aC6, aC12, bC6, bC12, aQ, bQ, αlj, αc] ++ inputs ++ [fpe, 1]).map
      hexBits ++ tableFields t ++ [f, eLj, eCrf, deLj, deCrf].map hexBits)
  return rows.toList

/-- Rows of a kernel `(r2, a_q, b_q, alpha_crf, lp, crf, fpe) → 3 floats`. -/
def rfTableRows (count : Nat)
    (kernel : Float → Float → Float → Float → NewenSpec.SoftCore.LambdaParams Float →
      NewenSpec.Nonbonded.CrfParams Float → Float → Float × Float × Float) :
    Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range count do
    let inputs ← drawCrfInputs
    let r2 ← drawR2 (inputs.getD 0 0)
    let aQ ← orZero (uniform (-2) 2)
    let bQ ← orZero (uniform (-2) 2)
    let αc ← drawAlpha
    let fpe ← drawFpe
    let t ← drawTable
    let (f, e, de) := kernel r2 aQ bQ αc (tableOf t) (crfOf inputs) fpe
    rows := rows.push (([r2, aQ, bQ, αc] ++ inputs ++ [fpe]).map hexBits ++ tableFields t ++
      [f, e, de].map hexBits)
  return rows.toList

/-- `perturbed_rf_self_term` rows: squared charges in `[0, 4)`. -/
def perturbedSelfRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 200 do
    let inputs ← drawCrfInputs
    let qa ← orZero (uniform (-2) 2)
    let qb ← orZero (uniform (-2) 2)
    let fpe ← drawFpe
    let t ← drawTable
    let (e, de) := NewenSpec.SoftCore.perturbedRfSelfFloat (qa * qa) (qb * qb) (tableOf t)
      (crfOf inputs) fpe
    rows := rows.push (([qa * qa, qb * qb] ++ inputs ++ [fpe]).map hexBits ++ tableFields t ++
      [e, de].map hexBits)
  return rows.toList

/-- `PerturbedLambdaParams::from_lambda` rows: the inputs, then the sixteen fields. -/
def lambdaTableRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let t ← drawTable
    let p := tableOf t
    rows := rows.push (tableFields t ++ ([p.aLjN, p.bLjN, p.aCrfN, p.bCrfN, p.aLjN1, p.bLjN1,
      p.aCrfN1, p.bCrfN1, p.aLjs, p.bLjs, p.aLjs2, p.bLjs2, p.aCrfs, p.bCrfs, p.aCrfs2,
      p.bCrfs2]).map hexBits)
  return rows.toList

/-- `lj_crf_soft_interaction_ext` rows. -/
def ljCrfExtRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 400 do
    let inputs ← drawCrfInputs
    let r2 ← drawR2 (inputs.getD 0 0)
    let (aC6, aC12, aQ) ← drawState
    let (bC6, bC12, bQ) ← drawState
    let αlj ← drawAlpha
    let αc ← drawAlpha
    let fpe ← drawFpe
    let lam ← drawLam
    let ((e1, e2, e3, e4), (d1, d2, d3, d4)) := NewenSpec.SoftCore.ljCrfSoftExtFloat r2 aC6 aC12
      bC6 bC12 aQ bQ αlj αc (crfOf inputs) fpe lam
    rows := rows.push (([r2, aC6, aC12, bC6, bC12, aQ, bQ, αlj, αc] ++ inputs ++
      [fpe, lam, e1, e2, e3, e4, d1, d2, d3, d4]).map hexBits)
  return rows.toList

/-- `rf_soft_interaction_ext` rows. -/
def rfExtRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 200 do
    let inputs ← drawCrfInputs
    let r2 ← drawR2 (inputs.getD 0 0)
    let aQ ← orZero (uniform (-2) 2)
    let bQ ← orZero (uniform (-2) 2)
    let αc ← drawAlpha
    let fpe ← drawFpe
    let lam ← drawLam
    let ((e1, e2), (d1, d2)) := NewenSpec.SoftCore.rfSoftExtFloat r2 aQ bQ αc (crfOf inputs) fpe lam
    rows := rows.push (([r2, aQ, bQ, αc] ++ inputs ++ [fpe, lam, e1, e2, d1, d2]).map hexBits)
  return rows.toList

structure SoftCoreRows where
  lambdaParams : List (List String)
  pair : List (List String)
  rfSoft : List (List String)
  rfExcluded : List (List String)
  rfSelf : List (List String)
  ljCrfExt : List (List String)
  rfExt : List (List String)

/-- The row sets of `vectors/softcore.json`, from a fixed seed. -/
def softCoreRows : SoftCoreRows :=
  (do return {
      lambdaParams := ← lambdaTableRows, pair := ← perturbedPairRows,
      rfSoft := ← rfTableRows 300 NewenSpec.SoftCore.rfSoftFloat,
      rfExcluded := ← rfTableRows 300 NewenSpec.SoftCore.perturbedRfExcludedFloat,
      rfSelf := ← perturbedSelfRows, ljCrfExt := ← ljCrfExtRows, rfExt := ← rfExtRows } :
    Gen _).run' 0x6e6577656e2d7363

def fieldRow (fields : List String) : String := "    [" ++ ", ".intercalate fields ++ "]"

def softCoreJson (r : SoftCoreRows) : String :=
  let part (name : String) (rows : List (List String)) : String :=
    s!"  \"{name}\": [\n" ++ ",\n".intercalate (rows.map fieldRow) ++ "\n  ]"
  let table := "lj, ljs, crf, crfs, dlj, dljs, dcrf, dcrfs, nlam"
  "{\n" ++
  "  \"spec\": \"nonbonded.perturbed_lambda_params, nonbonded.perturbed_lj_crf_pair, " ++
  "nonbonded.rf_soft_interaction, nonbonded.perturbed_rf_excluded_pair, " ++
  "nonbonded.perturbed_rf_self_term, nonbonded.lj_crf_soft_interaction_ext, " ++
  "nonbonded.rf_soft_interaction_ext\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.SoftCore.*Float; SplitMix64, fixed " ++
  "seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits); nlam is a plain integer\",\n" ++
  s!"  \"layout_lambda_params\": \"{table}; a_lj_n, b_lj_n, a_crf_n, b_crf_n, a_lj_n_1, " ++
  "b_lj_n_1, a_crf_n_1, b_crf_n_1, a_ljs, b_ljs, a_ljs2, b_ljs2, a_crfs, b_crfs, a_crfs2, " ++
  "b_crfs2\",\n" ++
  "  \"layout_perturbed_lj_crf_pair\": \"r2, a_c6, a_c12, b_c6, b_c12, a_q, b_q, alpha_lj, " ++
  s!"alpha_crf, cutoff, epsilon, rf_epsilon, rf_kappa, four_pi_eps_i, coulomb_scaling, {table}; " ++
  "force, e_lj, e_crf, de_lj, de_crf\",\n" ++
  "  \"layout_rf_soft_interaction\": \"r2, a_q, b_q, alpha_crf, cutoff, epsilon, rf_epsilon, " ++
  s!"rf_kappa, four_pi_eps_i, {table}; force, e_rf, de_rf\",\n" ++
  "  \"layout_perturbed_rf_excluded_pair\": \"as layout_rf_soft_interaction\",\n" ++
  "  \"layout_perturbed_rf_self_term\": \"qa2, qb2, cutoff, epsilon, rf_epsilon, rf_kappa, " ++
  s!"four_pi_eps_i, {table}; e_rf, de_rf\",\n" ++
  "  \"layout_lj_crf_soft_interaction_ext\": \"r2, a_c6, a_c12, b_c6, b_c12, a_q, b_q, " ++
  "alpha_lj, alpha_crf, cutoff, epsilon, rf_epsilon, rf_kappa, four_pi_eps_i, lam; A_e_lj, " ++
  "B_e_lj, A_e_crf, B_e_crf, A_de_lj, B_de_lj, A_de_crf, B_de_crf\",\n" ++
  "  \"layout_rf_soft_interaction_ext\": \"r2, a_q, b_q, alpha_crf, cutoff, epsilon, " ++
  "rf_epsilon, rf_kappa, four_pi_eps_i, lam; A_e_rf, B_e_rf, A_de_rf, B_de_rf\",\n" ++
  part "lambda_params" r.lambdaParams ++ ",\n" ++
  part "perturbed_lj_crf_pair" r.pair ++ ",\n" ++
  part "rf_soft_interaction" r.rfSoft ++ ",\n" ++
  part "perturbed_rf_excluded_pair" r.rfExcluded ++ ",\n" ++
  part "perturbed_rf_self_term" r.rfSelf ++ ",\n" ++
  part "lj_crf_soft_interaction_ext" r.ljCrfExt ++ ",\n" ++
  part "rf_soft_interaction_ext" r.rfExt ++ "\n}\n"

/-! ## Perturbed bonded: the kernels on random `Float` inputs, written as bits -/

/-- A λ derivative: `1` half the time (no individual λ), else uniform in `[-2, 2)`. -/
def drawDlam : Gen Float := do
  if (← below 2) == 0 then return 1 else uniform (-2) 2

/-- A force constant: `0` one time in eight (the term is absent in that state), else uniform. -/
def drawK (lo hi : Float) : Gen Float := orZero (uniform lo hi)

def piF : Float := NewenSpec.PerturbedBonded.piFloat

def pertQuarticRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let d2 ← logUniform 0.002 0.1
    let kA ← drawK 1e6 2e7
    let kB ← drawK 1e6 2e7
    let r0A ← uniform 0.09 0.2
    let r0B ← uniform 0.09 0.2
    let lam ← drawLam
    let dlam ← drawDlam
    let (e, f, de) := NewenSpec.PerturbedBonded.quarticFloat d2 kA kB r0A r0B lam dlam
    rows := rows.push ([d2, kA, kB, r0A, r0B, lam, dlam, e, f, de].map hexBits)
  return rows.toList

def pertAngleRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let cost ← uniform (-1) 1
    let kA ← drawK 100 800
    let kB ← drawK 100 800
    let cA ← uniform (-1) 0
    let cB ← uniform (-1) 0
    let lam ← drawLam
    let dlam ← drawDlam
    let (e, df, de) := NewenSpec.PerturbedBonded.angleFloat cost kA kB cA cB lam dlam
    rows := rows.push ([cost, kA, kB, cA, cB, lam, dlam, e, df, de].map hexBits)
  return rows.toList

/-- Impropers: `ζ` in `[-π, π)`, one time in eight in `[-3π, 3π)` (several turns to wrap). -/
def pertImproperRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let zeta ← if (← below 8) == 0 then uniform (-3 * piF) (3 * piF) else uniform (-piF) piF
    let kA ← drawK 0 400
    let kB ← drawK 0 400
    let qA ← drawFrom [0, 0.6154797086703874, -0.6154797086703874] (-piF) piF
    let qB ← drawFrom [0, 0.6154797086703874, -0.6154797086703874] (-piF) piF
    let lam ← drawLam
    let dlam ← drawDlam
    let (e, mkd, de, z) := NewenSpec.PerturbedBonded.improperFloat zeta kA kB qA qB lam dlam
    rows := rows.push ([zeta, kA, kB, qA, qB, lam, dlam, e, mkd, de, z].map hexBits)
  return rows.toList

/-- `chebyshev(m, c)` for `m = 0 … 7` (7: unsupported, `(0, 0)`), `c` in `[-1, 1)` or `±1`, `0`. -/
def pertChebyshevRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let m := (← below 8)
    let c ← drawFrom [-1, 0, 1] (-1) 1
    let (cm, dcm) := NewenSpec.PerturbedBonded.chebyshevFloat (m : Int) c
    rows := rows.push ([toString m] ++ [c, cm, dcm].map hexBits)
  return rows.toList

def pertDihedralRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let c ← uniform (-1) 1
    let kA ← drawK 0 10
    let kB ← drawK 0 10
    let pA ← drawFrom [-1, 1] (-1) 1
    let pB ← drawFrom [-1, 1] (-1) 1
    let mA := (← below 7)
    let mB := (← below 7)
    let lam ← drawLam
    let dlam ← drawDlam
    let (eA, eB, gA, gB, e, de) := NewenSpec.PerturbedBonded.dihedralFloat c kA pA (mA : Int) kB pB
      (mB : Int) lam dlam
    rows := rows.push ([c, kA, pA].map hexBits ++ [toString mA] ++ [kB, pB].map hexBits ++
      [toString mB] ++ [lam, dlam, eA, eB, gA, gB, e, de].map hexBits)
  return rows.toList

/-- Soft terms in a bond length, `cos θ` or `ζ`; one state's force constant may be `0` (the term
is absent there, and the kernel's caller gives it the other state's reference value). -/
def pertSoftRows : Gen (List (List String)) := do
  let mut rows := #[]
  for _ in List.range 400 do
    let (x, x0A, x0B, k) : Float × Float × Float × Float ← match ← below 3 with
      | 0 => do pure (← uniform 0.05 0.3, ← uniform 0.09 0.2, ← uniform 0.09 0.2, 4e5)
      | 1 => do pure (← uniform (-1) 1, ← uniform (-1) 0, ← uniform (-1) 0, 800)
      | _ => do pure (← uniform (-piF) piF, ← uniform (-1) 1, ← uniform (-1) 1, 400)
    let kA ← drawK 0 k
    let kB ← drawK 0 k
    let alpha ← drawFrom scAlphas 0 1000
    let lam ← drawLam
    let dlam ← drawDlam
    let (e, g, de) := NewenSpec.PerturbedBonded.softFloat x kA kB x0A x0B alpha lam dlam
    rows := rows.push ([x, kA, kB, x0A, x0B, alpha, lam, dlam, e, g, de].map hexBits)
  return rows.toList

structure PerturbedBondedRows where
  chebyshev : List (List String)
  quartic : List (List String)
  angle : List (List String)
  improper : List (List String)
  dihedral : List (List String)
  soft : List (List String)

/-- The row sets of `vectors/perturbed_bonded.json`, from a fixed seed. -/
def perturbedBondedRows : PerturbedBondedRows :=
  (do return {
      chebyshev := ← pertChebyshevRows, quartic := ← pertQuarticRows, angle := ← pertAngleRows,
      improper := ← pertImproperRows, dihedral := ← pertDihedralRows, soft := ← pertSoftRows } :
    Gen _).run' 0x6e6577656e2d7062

def perturbedBondedJson (r : PerturbedBondedRows) : String :=
  let part (name : String) (rows : List (List String)) : String :=
    s!"  \"{name}\": [\n" ++ ",\n".intercalate (rows.map fieldRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"bonded.chebyshev, bonded.perturbed_quartic_bond_term, " ++
  "bonded.perturbed_angle_term, bonded.perturbed_improper_term, " ++
  "bonded.perturbed_dihedral_term, bonded.soft_perturbed_term\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.PerturbedBonded.*Float; SplitMix64, " ++
  "fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits); multiplicities are plain integers\",\n" ++
  "  \"layout_chebyshev\": \"m, cos_phi; cos(m phi), d cos(m phi) / d cos_phi\",\n" ++
  "  \"layout_quartic\": \"dist2, k_a, k_b, r0_a, r0_b, lambda, lambda_derivative; energy, f, " ++
  "de_dl\",\n" ++
  "  \"layout_angle\": \"cost, k_a, k_b, cos0a, cos0b, lambda, lambda_derivative; energy, df, " ++
  "de_dl\",\n" ++
  "  \"layout_improper\": \"zeta, k_a, k_b, q0_a, q0_b, lambda, lambda_derivative; energy, " ++
  "minus_k_dzeta, de_dl, zeta_wrapped\",\n" ++
  "  \"layout_dihedral\": \"cos_phi, k_a, cospd_a, m_a, k_b, cospd_b, m_b, lambda, " ++
  "lambda_derivative; e_a, e_b, g_a, g_b, energy, de_dl\",\n" ++
  "  \"layout_soft\": \"x, k_a, k_b, x0_a, x0_b, alpha, lambda, lambda_derivative; energy, " ++
  "k_eff_diff, de_dl\",\n" ++
  part "chebyshev" r.chebyshev ++ ",\n" ++ part "quartic" r.quartic ++ ",\n" ++
  part "angle" r.angle ++ ",\n" ++ part "improper" r.improper ++ ",\n" ++
  part "dihedral" r.dihedral ++ ",\n" ++ part "soft" r.soft ++ "\n}\n"

/-! ## Lattice sum: the charge-shape and real-space kernels on grids and random `Float` inputs

`vectors/lattice_sum.json`: every shape (`NSHAPE = −1 … 10`) at `ξ = 0`, `1`, the doubles next to
`1` on either side, a log grid in `[10⁻³, 3]`, random points in `[0, 3]` and a few large and
negative ones (`erfc` underflows near 27); `κ = 0`, a log grid in `[10⁻³, 40]`, random points in
`[0, 40]`; pair rows with `r²` log-uniform in `[10⁻⁴, 4]` nm² or at the support edge `r = a`,
widths `a ∈ {0.4, 0.8, 1.4}` and the `C6`, `C12`, charge families of `nonbonded.json`; the
self-term constants of every shape for several widths and volumes. -/

section LatticeSumRows

open NewenTheory.Interactions.LatticeSum (Shape TPOrder)
open NewenSpec.LatticeSum

def lsShapes : List Shape :=
  [.gaussian, .poly .o0, .poly .o1, .poly .o2, .poly .o3, .poly .o4, .poly .o5, .poly .o6,
    .poly .o7, .poly .o8, .poly .o9, .poly .o10]

def nshapeFloat (s : Shape) : Float := Float.ofInt s.nshape

/-- `n` points log-spaced from `a` to `b`. -/
def logGrid (a b : Float) (n : Nat) : List Float :=
  (List.range n).map fun i =>
    Float.exp (Float.log a + (Float.log b - Float.log a) * i.toFloat / (n - 1).toFloat)

def lsWidths : List Float := [0.4, 0.8, 1.4]

/-- The doubles next to 1: `1 − 2⁻⁵³` and `1 + 2⁻⁵²`. -/
def belowOne : Float := Float.ofBits 0x3fefffffffffffff
def aboveOne : Float := Float.ofBits 0x3ff0000000000001

def shapeSwitchRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for s in lsShapes do
    let mut xs := [0, 1, belowOne, aboveOne, -0.25, 5, 10, 26.5, 27.3, 30] ++ logGrid 1e-3 3 120
    for _ in List.range 60 do
      xs := xs ++ [← uniform 0 3]
    for xi in xs do
      let (eta, dEta) := shapeSwitchFloat s xi
      rows := rows.push [nshapeFloat s, xi, eta, dEta]
  return rows.toList

def shapeFourierRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for s in lsShapes do
    let mut ks := [0, 1, NewenSpec.LatticeSum.piFloat, 2 * NewenSpec.LatticeSum.piFloat, 40] ++ logGrid 1e-3 40 120
    for _ in List.range 60 do
      ks := ks ++ [← uniform 0 40]
    for kappa in ks do
      let (g, dg) := shapeFourierFloat s kappa
      rows := rows.push [nshapeFloat s, kappa, g, dg]
  return rows.toList

/-- A squared distance for a shape of width `a`: log-uniform in `[1e-4, 4]` nm², or at the
support edge, `a²` and the doubles `a²(1 ± 2⁻⁴⁰)` around it. -/
def drawLsR2 (width : Float) : Gen Float := do
  match ← below 20 with
  | 0 => return width * width
  | 1 => return width * width * (1 - 9.094947017729282e-13)
  | 2 => return width * width * (1 + 9.094947017729282e-13)
  | _ => logUniform 1e-4 4

def lsShapeAndWidth : Gen (Shape × Float × Float) := do
  let s := lsShapes.getD (← below lsShapes.length) .gaussian
  let width ← drawFrom lsWidths 0.2 1.5
  return (s, width, 1 / width)

def ljLsPairRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 2400 do
    let (s, width, invWidth) ← lsShapeAndWidth
    let r2 ← drawLsR2 width
    let c6 ← orZero (logUniform 1e-5 1e-2)
    let c12 ← orZero (logUniform 1e-8 1e-5)
    let q ← orZero (uniform (-150) 150)
    let (f, eLj, eLs) := ljLsPairFloat r2 c6 c12 q s invWidth
    rows := rows.push [r2, c6, c12, q, nshapeFloat s, width, invWidth, f, eLj, eLs]
  return rows.toList

def lsExcludedPairRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 1200 do
    let (s, width, invWidth) ← lsShapeAndWidth
    let r2 ← drawLsR2 width
    let q ← orZero (uniform (-150) 150)
    let (f, e) := lsExcludedPairFloat r2 q s invWidth
    rows := rows.push [r2, q, nshapeFloat s, width, invWidth, f, e]
  return rows.toList

def selfConstantsRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for s in lsShapes do
    for width in lsWidths ++ [← uniform 0.2 1.5] do
      for volume in [1, 27, 125, ← uniform 1 1000] do
        let (a1, a3) := selfConstantsFloat s width volume
        rows := rows.push [nshapeFloat s, width, volume, a1, a3]
  return rows.toList

structure LatticeSumRows where
  switch : List (List Float)
  fourier : List (List Float)
  pair : List (List Float)
  excluded : List (List Float)
  selfConstants : List (List Float)

def latticeSumRows : LatticeSumRows :=
  (do return ⟨← shapeSwitchRows, ← shapeFourierRows, ← ljLsPairRows, ← lsExcludedPairRows,
      ← selfConstantsRows⟩ : Gen _).run' 0x6e6577656e2d6c73

def latticeSumJson (r : LatticeSumRows) : String :=
  let part (name : String) (rows : List (List Float)) : String :=
    s!"  \"{name}\": [\n" ++ ",\n".intercalate (rows.map hexRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"lattice_sum.shape_switch, lattice_sum.shape_fourier, lattice_sum.lj_ls_pair, " ++
  "lattice_sum.ls_excluded_pair, lattice_sum.self_constants\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.LatticeSum.shapeSwitchFloat, " ++
  "shapeFourierFloat, ljLsPairFloat, lsExcludedPairFloat, selfConstantsFloat; SplitMix64, " ++
  "fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits (f64::to_bits, " ++
  "Float.toBits); nshape is NSHAPE as a float (-1 Gaussian, 0..10 polynomial)\",\n" ++
  "  \"layout_shape_switch\": \"nshape, xi, eta, d_eta\",\n" ++
  "  \"layout_shape_fourier\": \"nshape, kappa, gamma_hat, gamma_hat_prime\",\n" ++
  "  \"layout_lj_ls_pair\": \"r2, c6, c12, q_eps, nshape, width, inv_width = 1/width, f_mag, " ++
  "e_lj, e_ls\",\n" ++
  "  \"layout_ls_excluded_pair\": \"r2, q_eps, nshape, width, inv_width = 1/width, f_mag, " ++
  "e_ls\",\n" ++
  "  \"layout_self_constants\": \"nshape, width, volume, a1, a3\",\n" ++
  part "shape_switch" r.switch ++ ",\n" ++ part "shape_fourier" r.fourier ++ ",\n" ++
  part "lj_ls_pair" r.pair ++ ",\n" ++ part "ls_excluded_pair" r.excluded ++ ",\n" ++
  part "self_constants" r.selfConstants ++ "\n}\n"

/-- The integer literals of the lattice-sum models, each an exact double through `NatCast`. -/
def lsLiterals : List Nat :=
  [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 19, 21, 22, 24, 25, 26, 30, 31, 32, 33, 35,
    36, 39, 45, 47, 48, 49, 51, 55, 56, 57, 60, 62, 64, 65, 66, 69, 73, 77, 87, 90, 91, 104, 105,
    123, 128, 136, 141, 151, 160, 187, 192, 203, 207, 285, 312, 315, 335, 384, 405, 455, 545,
    630, 640, 697, 936, 945, 975, 1155, 1536, 1545, 1665, 1920, 2075, 3465, 5040, 7560, 9009,
    12645, 22005, 75600, 831600, 1247400, 16216200]

end LatticeSumRows

/-! ## Ewald: the wave vectors, the k-space sum, `A2`, the self, A and surface terms

`vectors/ewald.json`: random small systems (2 to 12 atoms, neutral or charged) in cubic,
rectangular and triclinic boxes (GROMOS's lower-triangular form, angles 60°–120°), `NK` up to 3,
KCUT that cuts the cube of integer triples and KCUT that does not, the twelve shapes, widths
0.4–1.4, `EPSLS` 0 and 61, and `1/ε₀` both as gromosXX's literal and as `4π · 138.9354`. Every row
holds its inputs and the model's outputs; lists of variable length carry their length. -/

section EwaldRows

open NewenTheory.Interactions.LatticeSum (Shape TPOrder)
open NewenSpec.LatticeSum

def v3List (v : NewenSpec.LatticeSum.V3 Float) : List Float := [v.x, v.y, v.z]

def symList (m : Sym Float) : List Float := [m.xx, m.xy, m.xz, m.yy, m.yz, m.zz]

/-- `1/ε₀`: gromosXX's literal `eps0_i = 1745.9137` or, with a topology, `4π · 138.9354` as
gromosXX computes it (`138.9354 * 4.0 * π`). -/
def eps0Invs : List Float := [1745.9137, 138.9354 * 4 * NewenSpec.LatticeSum.piFloat]

/-- A box: `(rectangular, a, b, c)`; cubic, rectangular or triclinic. -/
def drawBox : Gen (Bool × NewenSpec.LatticeSum.V3 Float × NewenSpec.LatticeSum.V3 Float × NewenSpec.LatticeSum.V3 Float) := do
  let la ← uniform 1.2 3
  match ← below 3 with
  | 0 => return (true, ⟨la, 0, 0⟩, ⟨0, la, 0⟩, ⟨0, 0, la⟩)
  | 1 => return (true, ⟨la, 0, 0⟩, ⟨0, ← uniform 1.2 3, 0⟩, ⟨0, 0, ← uniform 1.2 3⟩)
  | _ =>
    let lb ← uniform 1.2 3
    let lc ← uniform 1.2 3
    let deg := NewenSpec.LatticeSum.piFloat / 180
    let al := (← uniform 60 120) * deg
    let be := (← uniform 60 120) * deg
    let ga := (← uniform 60 120) * deg
    let cx := lc * Float.cos be
    let cy := lc * (Float.cos al - Float.cos be * Float.cos ga) / Float.sin ga
    let cz := Float.sqrt (lc * lc - cx * cx - cy * cy)
    return (false, ⟨la, 0, 0⟩, ⟨lb * Float.cos ga, lb * Float.sin ga, 0⟩, ⟨cx, cy, cz⟩)

def drawShape : Gen Shape := do return lsShapes.getD (← below lsShapes.length) .gaussian

def drawNk : Gen (ℕ × ℕ × ℕ) := do
  return ((← below 3) + 1, (← below 3) + 1, (← below 3) + 1)

/-- KCUT that cuts (`|k|² ≤ KCUT` inside the cube of triples) or that does not. -/
def drawKcut : Gen Float := do
  match ← below 3 with
  | 0 => return 1e6
  | _ => uniform 5 200

/-- `n` atoms in the box (and a little outside it), with charges summing to zero or not. -/
def drawAtoms (a b c : NewenSpec.LatticeSum.V3 Float) : Gen (List (NewenSpec.LatticeSum.V3 Float × Float)) := do
  let n := (← below 11) + 2
  let neutral := (← below 2) == 0
  let mut atoms : Array (NewenSpec.LatticeSum.V3 Float × Float) := #[]
  for _ in List.range n do
    let fa ← uniform (-0.1) 1.1
    let fb ← uniform (-0.1) 1.1
    let fc ← uniform (-0.1) 1.1
    let r : NewenSpec.LatticeSum.V3 Float := ⟨fa * a.x + fb * b.x + fc * c.x, fa * a.y + fb * b.y + fc * c.y,
      fa * a.z + fb * b.z + fc * c.z⟩
    atoms := atoms.push (r, ← uniform (-1) 1)
  if neutral then
    let total := (atoms.toList.take (n - 1)).foldl (fun s a => s + a.2) 0
    atoms := atoms.modify (n - 1) fun a => (a.1, -total)
  return atoms.toList

def kvecList (kv : KVec Float) : List Float :=
  v3List kv.k ++ [kv.k2, kv.k2i, kv.ak, kv.gammaHat, kv.k2iGammaHat, kv.akGammaHatPrime]

def natF (n : ℕ) : Float := Float.ofNat n

def reciprocalBoxRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 120 do
    let (rect, a, b, c) ← drawBox
    let (v, u, w, x) := reciprocalBoxFloat a b c rect
    rows := rows.push ([if rect then 1 else 0] ++ v3List a ++ v3List b ++ v3List c ++ [v] ++
      v3List u ++ v3List w ++ v3List x)
  return rows.toList

/-- One case of the sum over wave vectors: the box, its reciprocal basis and the k-space
parameters. -/
structure KCase where
  rect : Bool
  a : NewenSpec.LatticeSum.V3 Float
  b : NewenSpec.LatticeSum.V3 Float
  c : NewenSpec.LatticeSum.V3 Float
  volume : Float
  u : NewenSpec.LatticeSum.V3 Float
  v : NewenSpec.LatticeSum.V3 Float
  w : NewenSpec.LatticeSum.V3 Float
  nk : ℕ × ℕ × ℕ
  kcut : Float
  shape : Shape
  width : Float

def drawKCase : Gen KCase := do
  let (rect, a, b, c) ← drawBox
  let (volume, u, v, w) := reciprocalBoxFloat a b c rect
  return ⟨rect, a, b, c, volume, u, v, w, ← drawNk, ← drawKcut, ← drawShape,
    ← uniform 0.4 1.4⟩

def kCaseList (k : KCase) : List Float :=
  [natF k.nk.1, natF k.nk.2.1, natF k.nk.2.2, k.kcut, nshapeFloat k.shape, k.width] ++
    v3List k.u ++ v3List k.v ++ v3List k.w

def kVectorRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 24 do
    let k ← drawKCase
    let ks := kVectorsFloat k.nk k.kcut k.u k.v k.w k.shape k.width
    rows := rows.push (kCaseList k ++ [natF ks.length] ++ ks.flatMap kvecList)
  return rows.toList

def kspaceRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 160 do
    let k ← drawKCase
    let eps0Inv ← pick eps0Invs
    let atoms ← drawAtoms k.a k.b k.c
    let ks := kVectorsFloat k.nk k.kcut k.u k.v k.w k.shape k.width
    let (e, forces, a2t, vir, da2t) := ewaldKspaceFloat atoms ks k.volume eps0Inv
    rows := rows.push (kCaseList k ++ [k.volume, eps0Inv, natF atoms.length] ++
      atoms.flatMap (fun a => v3List a.1 ++ [a.2]) ++ [e, a2t] ++ symList vir ++ symList da2t ++
      forces.flatMap v3List)
  return rows.toList

def a2CubicRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 200 do
    let s ← drawShape
    let edge ← uniform 1.2 6
    let width ← uniform 0.4 1.4
    let (a1, a3) := selfConstantsFloat s width (edge * edge * edge)
    let (a2, d) := a2CubicFloat (-2.83729748) edge a1 a3
    rows := rows.push [-2.83729748, edge, a1, a3, a2, d]
  return rows.toList

def a2ShellsRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 48 do
    let k ← drawKCase
    let tola2 ← pick [1e-4, 1e-6, 1e-8, 1e-10]
    let maxShells := (← below 12) + 1
    let (a2, shells, conv, d) := a2ShellsFloat k.shape k.width k.u k.v k.w k.volume tola2 maxShells
    rows := rows.push ([nshapeFloat k.shape, k.width] ++ v3List k.u ++ v3List k.v ++ v3List k.w ++
      [k.volume, tola2, natF maxShells, a2, natF shells, if conv then 1 else 0] ++ symList d)
  return rows.toList

def chargeSumRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 100 do
    let n := (← below 12) + 1
    let mut qs := #[]
    for _ in List.range n do
      qs := qs.push (← uniform (-2) 2)
    let (sq, sq2) := chargeSumsFloat qs.toList
    rows := rows.push ([natF n] ++ qs.toList ++ [sq, sq2])
  return rows.toList

def selfATermsRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 300 do
    let a1 ← uniform (-1) 0
    let a2 ← uniform (-3) 3
    let a3 ← uniform (-8) 0
    let a2t ← uniform (-3) 3
    let sumQ ← orZero (uniform (-5) 5)
    let sumQ2 ← uniform 0 50
    let eps0Inv ← pick eps0Invs
    let (e1, e2) := selfATermsFloat a1 a2 a3 a2t sumQ sumQ2 eps0Inv
    rows := rows.push [a1, a2, a3, a2t, sumQ, sumQ2, eps0Inv, e1, e2]
  return rows.toList

def drawSym : Gen (Sym Float) := do
  return ⟨← uniform (-5) 5, ← uniform (-5) 5, ← uniform (-5) 5, ← uniform (-5) 5,
    ← uniform (-5) 5, ← uniform (-5) 5⟩

def selfAVirialRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 200 do
    let a1 ← uniform (-1) 0
    let sumQ ← orZero (uniform (-5) 5)
    let sumQ2 ← uniform 0 50
    let d ← drawSym
    let dt ← drawSym
    let fpe ← pick [138.9354, 138.935458]
    let (sv, cv) := selfAVirialsFloat a1 sumQ sumQ2 d dt fpe
    rows := rows.push ([a1, sumQ, sumQ2] ++ symList d ++ symList dt ++ [fpe] ++ symList sv ++
      symList cv)
  return rows.toList

def surfaceRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 160 do
    let (rect, a, b, c) ← drawBox
    let (volume, _, _, _) := reciprocalBoxFloat a b c rect
    let atoms ← drawAtoms a b c
    let centre : NewenSpec.LatticeSum.V3 Float := ⟨a.x / 2 + b.x / 2 + c.x / 2, a.y / 2 + b.y / 2 + c.y / 2,
      a.z / 2 + b.z / 2 + c.z / 2⟩
    let epsLs ← pick [0, 61, 1, 78.5]
    let eps0Inv ← pick eps0Invs
    let (e, forces, vir) := surfaceFloat atoms centre epsLs volume eps0Inv
    rows := rows.push (v3List centre ++ [epsLs, volume, eps0Inv, natF atoms.length] ++
      atoms.flatMap (fun at_ => v3List at_.1 ++ [at_.2]) ++ [e] ++ symList vir ++
      forces.flatMap v3List)
  return rows.toList

structure EwaldRows where
  reciprocal : List (List Float)
  kvectors : List (List Float)
  kspace : List (List Float)
  a2Cubic : List (List Float)
  a2Shells : List (List Float)
  chargeSums : List (List Float)
  selfA : List (List Float)
  selfAVirial : List (List Float)
  surface : List (List Float)

def ewaldRows : EwaldRows :=
  (do return ⟨← reciprocalBoxRows, ← kVectorRows, ← kspaceRows, ← a2CubicRows, ← a2ShellsRows,
      ← chargeSumRows, ← selfATermsRows, ← selfAVirialRows, ← surfaceRows⟩ : Gen _).run'
    0x6e6577656e2d6577

def ewaldJson (r : EwaldRows) : String :=
  let part (name : String) (rows : List (List Float)) : String :=
    s!"  \"{name}\": [\n" ++ ",\n".intercalate (rows.map hexRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"lattice_sum.reciprocal_box, ewald.k_vectors, ewald.kspace, lattice_sum.a2_cubic, " ++
  "lattice_sum.a2_shells, lattice_sum.self_and_a_terms, lattice_sum.surface_term\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.LatticeSum.reciprocalBoxFloat, " ++
  "kVectorsFloat, ewaldKspaceFloat, a2CubicFloat, a2ShellsFloat, chargeSumsFloat, " ++
  "selfATermsFloat, selfAVirialsFloat, surfaceFloat; SplitMix64, fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits; counts, flags " ++
  "and nshape are floats; a vector is x, y, z; a symmetric tensor xx, xy, xz, yy, yz, zz\",\n" ++
  "  \"layout_reciprocal_box\": \"rectangular, a, b, c, volume, b1, b2, b3\",\n" ++
  "  \"layout_k_vectors\": \"nkx, nky, nkz, kcut, nshape, width, b1, b2, b3, count, then per " ++
  "wave vector k, k2, k2i, ak, gamma_hat, k2i_gamma_hat, ak_gamma_hat_prime\",\n" ++
  "  \"layout_kspace\": \"nkx, nky, nkz, kcut, nshape, width, b1, b2, b3, volume, eps0_inv, " ++
  "n_atoms, per atom r, q; energy, a2_tilde, virial, a2_tilde_derivative, per atom force " ++
  "(the wave vectors are k_vectors of the first fields)\",\n" ++
  "  \"layout_a2_cubic\": \"wigner, edge, a1, a3, a2, d_a2\",\n" ++
  "  \"layout_a2_shells\": \"nshape, width, b1, b2, b3, volume, tola2, max_shells, a2, shells, " ++
  "converged, derivative\",\n" ++
  "  \"layout_charge_sums\": \"n, the charges, sum_q, sum_q2\",\n" ++
  "  \"layout_self_and_a_terms\": \"a1, a2, a3, a2_tilde, sum_q, sum_q2, eps0_inv, ls_self, " ++
  "ls_a_term\",\n" ++
  "  \"layout_self_and_a_virials\": \"a1, sum_q, sum_q2, d_a2, d_a2_tilde, four_pi_eps_i, " ++
  "self_virial, constant_virial\",\n" ++
  "  \"layout_surface_term\": \"centre, eps_ls, volume, eps0_inv, n_atoms, per atom r, q; " ++
  "energy, virial, per atom force\",\n" ++
  part "reciprocal_box" r.reciprocal ++ ",\n" ++ part "k_vectors" r.kvectors ++ ",\n" ++
  part "kspace" r.kspace ++ ",\n" ++ part "a2_cubic" r.a2Cubic ++ ",\n" ++
  part "a2_shells" r.a2Shells ++ ",\n" ++ part "charge_sums" r.chargeSums ++ ",\n" ++
  part "self_and_a_terms" r.selfA ++ ",\n" ++ part "self_and_a_virials" r.selfAVirial ++ ",\n" ++
  part "surface_term" r.surface ++ "\n}\n"

end EwaldRows

/-! ## P3M: the per-point kernels, and the transforms within a tolerance

`vectors/p3m.json`: the assignment functions, the stencil weights and centres, `f_p`, the
finite-difference coefficients and `P̂_p` on grids of their arguments; grid matrices of random
boxes; the influence function at single wave vectors (grids 4³ and 8³, every order, NFDORD 0–3,
NALIAS 0–2, four shapes, widths 0.4 and 1.4 nm, cubic, rectangular and triclinic boxes) and on
whole 4³ grids, its quality, its box correction, the force error; charge and squared-charge
assignment, the energy, potential, fields, `Ã2` and the interpolated forces on small grids. All
of these are compared bit for bit. Two sections are compared within a tolerance: `dft` (the
direct transform `dftFloat` of grids whose entries are `sin(0.37 i + seed)` and
`cos(0.91 i + 2 seed)`, so only the outputs are stored) and `p3m_energy` (the whole k-space
energy on an 8³ grid through `dftFloat`). -/

section P3MRows

open NewenTheory.Interactions.LatticeSum (Shape TPOrder)
open NewenSpec.LatticeSum (Sym reciprocalBoxFloat piFloat)
open NewenSpec.P3M

abbrev LV3 := NewenSpec.LatticeSum.V3
abbrev PM3 := NewenSpec.P3M.M3

def m3List (m : PM3 Float) : List Float := v3List m.r0 ++ v3List m.r1 ++ v3List m.r2

def intF (z : ℤ) : Float := Float.ofInt z

def p3mShapes : List Shape := [.gaussian, .poly .o0, .poly .o3, .poly .o7]

/-- The doubles on either side of a positive `x`. -/
def p3mNeighbours (x : Float) : List Float :=
  [x, Float.ofBits (x.toBits - 1), Float.ofBits (x.toBits + 1)]

def p3mSpecials : List Float :=
  let pos := [0.5, 1, 1.5, 2, 2.5].flatMap p3mNeighbours
  [0] ++ pos ++ pos.map (fun x => -x)

def assignmentRows : List (List Float) :=
  (List.range 5).flatMap fun i =>
    let p := i + 1
    let xs := (List.range 97).map (fun j => j.toFloat / 16 - 3) ++ p3mSpecials
    xs.map fun x => [natF p, x, assignmentFloat p x]

def padFive (xs : List Float) : List Float := (xs ++ [0, 0, 0, 0, 0]).take 5

def weightsRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for i in List.range 5 do
    let p := i + 1
    let mut us := (List.range 65).map (fun j => j.toFloat / 32 - 1) ++
      [0.5, 1.5, 2.5, 31.5, -0.5, 3, 31]
    for _ in List.range 20 do
      us := us ++ [← uniform 0 32]
    for u in us do
      rows := rows.push ([natF p, u, intF (nearestIntFloat p u)] ++ padFive (weightsFloat p u))
  return rows.toList

def selftermRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for i in List.range 5 do
    let p := i + 1
    for n in List.range p do
      let mut ss := (List.range 17).map (fun j => j.toFloat / 16 - 0.5)
      for _ in List.range 3 do
        ss := ss ++ [← uniform (-0.5) 0.5]
      for sv in ss do
        rows := rows.push [natF p, natF n, sv, selftermFpFloat p n sv]
  return rows.toList

def fdRows : List (List Float) :=
  (List.range 6).map fun q => [natF q] ++ fdCoefficientsFloat q

def sincRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for i in List.range 5 do
    let p := i + 1
    let mut xs := (List.range 97).map (fun j => j.toFloat / 4 - 12)
    for _ in List.range 5 do
      xs := xs ++ [← uniform (-40) 40]
    for x in xs do
      rows := rows.push [natF p, x, sincPowFloat x p]
  return rows.toList

def drawGridSize : Gen ℕ := do return [2, 4, 8, 16, 32].getD (← below 5) 4

def boxM3 (a b c : LV3 Float) : PM3 Float := ⟨a, b, c⟩

def boxList (rect : Bool) (a b c : LV3 Float) : List Float :=
  [if rect then 1 else 0] ++ v3List a ++ v3List b ++ v3List c

def meshRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 60 do
    let (rect, a, b, c) ← drawBox
    let g := (← drawGridSize, ← drawGridSize, ← drawGridSize)
    let m := meshFloat (boxM3 a b c) g
    rows := rows.push (boxList rect a b c ++ [natF g.1, natF g.2.1, natF g.2.2] ++ m3List m.h ++
      m3List m.hInv ++ [m.cellVolume])
  return rows.toList

/-- The inputs of the influence function for a box and a grid, and the box volume. -/
def mkInputs (order nfdord nalias : ℕ) (sh : Shape) (width : Float) (rect : Bool)
    (a b c : LV3 Float) (g : ℕ × ℕ × ℕ) : Inputs Float × Float :=
  let (vol, u, v, w) := reciprocalBoxFloat a b c rect
  (⟨order, nfdord, nalias, sh, width, meshFloat (boxM3 a b c) g, ⟨u, v, w⟩⟩, vol)

/-- One influence case: the grid, the orders, the shape, the box. -/
structure ICase where
  g : ℕ × ℕ × ℕ
  order : ℕ
  nfdord : ℕ
  nalias : ℕ
  shape : Shape
  width : Float
  rect : Bool
  a : LV3 Float
  b : LV3 Float
  c : LV3 Float

def ICase.inputs (k : ICase) : Inputs Float × Float :=
  mkInputs k.order k.nfdord k.nalias k.shape k.width k.rect k.a k.b k.c k.g

def ICase.list (k : ICase) : List Float :=
  [natF k.g.1, natF k.g.2.1, natF k.g.2.2, natF k.order, natF k.nfdord, natF k.nalias,
    nshapeFloat k.shape, k.width] ++ boxList k.rect k.a k.b k.c

def drawICase (grids : List ℕ) (maxAlias : ℕ) : Gen ICase := do
  let n := grids.getD (← below grids.length) 4
  let (rect, a, b, c) ← drawBox
  return ⟨(n, n, n), (← below 5) + 1, ← below 4, ← below (maxAlias + 1),
    p3mShapes.getD (← below p3mShapes.length) .gaussian, ← pick [0.4, 1.4], rect, a, b, c⟩

/-- A wave index in `(−n/2, n/2]`, the triple not zero. -/
def drawL (g : ℕ × ℕ × ℕ) : Gen (ℤ × ℤ × ℤ) := do
  let one := fun (n : ℕ) => do return ((← below n : ℕ) : ℤ) - ((n / 2 : ℕ) : ℤ) + 1
  let l := (← one g.1, ← one g.2.1, ← one g.2.2)
  return if l = (0, 0, 0) then (1, 0, 0) else l

def lList (l : ℤ × ℤ × ℤ) : List Float := [intF l.1, intF l.2.1, intF l.2.2]

def influenceRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 250 do
    let k ← drawICase [4, 8] 2
    let l ← drawL k.g
    let wg := (← below 2) == 0
    let (g, gam, q) := influenceAtFloat k.inputs.1 l wg
    rows := rows.push (k.list ++ lList l ++ [if wg then 1 else 0, g] ++ symList gam ++ [q])
  return rows.toList

def qualityRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 80 do
    let k ← drawICase [4, 8] 2
    let l ← drawL k.g
    let g := (influenceAtFloat k.inputs.1 l false).1 * (← uniform 0.9 1.1)
    let (f, s) := qualityAtFloat k.inputs.1 l g
    rows := rows.push (k.list ++ lList l ++ [g, f, s])
  return rows.toList

def drawM3 (lo hi : Float) : Gen (PM3 Float) := do
  return ⟨⟨← uniform lo hi, ← uniform lo hi, ← uniform lo hi⟩,
    ⟨← uniform lo hi, ← uniform lo hi, ← uniform lo hi⟩,
    ⟨← uniform lo hi, ← uniform lo hi, ← uniform lo hi⟩⟩

def boxCorrectionRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 100 do
    let g ← uniform (-2) 2
    let gam ← drawSym
    let (_, a, b, c) ← drawBox
    let t0 := (boxM3 a b c).inv
    let dl ← drawM3 (-0.05) 0.05
    rows := rows.push ([g] ++ symList gam ++ m3List t0 ++ m3List dl ++
      [boxCorrectionFloat g gam t0 dl])
  return rows.toList

def forceErrorRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 50 do
    let q ← logUniform 1e-12 1e-2
    let s2 ← uniform 0.5 500
    let v ← uniform 1 30
    let n := (← below 3000) + 1
    let fpe ← pick [138.9354, 138.935458]
    rows := rows.push [q, s2, v, natF n, fpe, forceErrorFloat q s2 v n fpe]
  return rows.toList

def influenceFunctionRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 6 do
    let k ← drawICase [4] 1
    let (inp, vol) := k.inputs
    let (gh, gam, q) := influenceFunctionFloat inp vol
    rows := rows.push (k.list ++ [natF gh.size] ++ gh.toList ++ gam.toList.flatMap symList ++ [q])
  return rows.toList

/-- Atoms with fractional coordinates in `[0, 1)` (in the box). -/
def drawAtomsIn (a b c : LV3 Float) (maxN : ℕ) : Gen (List (LV3 Float × Float)) := do
  let n := (← below (maxN - 1)) + 2
  let mut atoms := #[]
  for _ in List.range n do
    let fa ← unit
    let fb ← unit
    let fc ← unit
    let r : LV3 Float := ⟨fa * a.x + fb * b.x + fc * c.x, fa * a.y + fb * b.y + fc * c.y,
      fa * a.z + fb * b.z + fc * c.z⟩
    atoms := atoms.push (r, ← uniform (-1) 1)
  return atoms.toList

def atomsList (atoms : List (LV3 Float × Float)) : List Float :=
  [natF atoms.length] ++ atoms.flatMap fun at_ => v3List at_.1 ++ [at_.2]

def drawSmallGrid : Gen (ℕ × ℕ × ℕ) := do
  return [(4, 4, 4), (8, 4, 2), (2, 4, 8)].getD (← below 3) (4, 4, 4)

def assignRows (squared : Bool) : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 30 do
    let g ← drawSmallGrid
    let p := (← below 5) + 1
    let (rect, a, b, c) ← drawBox
    let atoms ← drawAtomsIn a b c 8
    let m := meshFloat (boxM3 a b c) g
    let out := if squared then assignSquaredFloat m p atoms else assignChargesFloat m p atoms
    rows := rows.push ([natF g.1, natF g.2.1, natF g.2.2, natF p] ++ boxList rect a b c ++
      atomsList atoms ++ out.toList)
  return rows.toList

def p3mKspaceRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 8 do
    let n := 64
    let wg := (← below 2) == 0
    let vc ← uniform 1e-3 1e-2
    let e0 ← pick eps0Invs
    let mut pts := #[]
    for _ in List.range n do
      pts := pts.push (← uniform 0 1, ← drawSym, ← uniform (-1) 1, ← uniform (-1) 1)
    let (e, vir) := kspaceEnergyFloat pts.toList wg vc e0
    let phi := pts.toList.flatMap fun pt =>
      let (r, i) := potentialHatFloat n pt.1 pt.2.2.1 pt.2.2.2 e0
      [r, i]
    rows := rows.push ([natF n, if wg then 1 else 0, vc, e0] ++
      pts.toList.flatMap (fun pt => [pt.1] ++ symList pt.2.1 ++ [pt.2.2.1, pt.2.2.2]) ++
      [e] ++ symList vir ++ phi)
  return rows.toList

def a2Rows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 6 do
    let n := 64
    let wg := (← below 2) == 0
    let s2 ← uniform 0.5 50
    let mut pts := #[]
    for _ in List.range n do
      pts := pts.push (← uniform 0 1, ← drawSym, ← uniform (-1) 1)
    let (a2t, d) := a2TildeFloat pts.toList wg s2
    rows := rows.push ([natF n, if wg then 1 else 0, s2] ++
      pts.toList.flatMap (fun pt => [pt.1] ++ symList pt.2.1 ++ [pt.2.2]) ++ [a2t] ++ symList d)
  return rows.toList

/-- Every grid point `(x, y, z)` in storage order. -/
def gridPoints (g : ℕ × ℕ × ℕ) : List (ℕ × ℕ × ℕ) :=
  (List.range g.1).flatMap fun x => (List.range g.2.1).flatMap fun y =>
    (List.range g.2.2).map fun z => (x, y, z)

def fieldIkRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 4 do
    let g ← drawSmallGrid
    let (rect, a, b, c) ← drawBox
    let (_, u, v, w) := reciprocalBoxFloat a b c rect
    let m := meshFloat (boxM3 a b c) g
    let mut vals := #[]
    for _ in gridPoints g do
      vals := vals.push (← uniform (-1) 1, ← uniform (-1) 1)
    let out := (gridPoints g).zip vals.toList |>.flatMap fun (pt, ri) =>
      let (er, ei) := fieldIkFloat m ⟨u, v, w⟩ pt.1 pt.2.1 pt.2.2 ri.1 ri.2
      v3List er ++ v3List ei
    rows := rows.push ([natF g.1, natF g.2.1, natF g.2.2] ++ boxList rect a b c ++
      vals.toList.flatMap (fun ri => [ri.1, ri.2]) ++ out)
  return rows.toList

def fieldFdRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 10 do
    let g ← drawSmallGrid
    let q := (← below 5) + 1
    let (rect, a, b, c) ← drawBox
    let m := meshFloat (boxM3 a b c) g
    let mut phi := #[]
    for _ in gridPoints g do
      phi := phi.push (← uniform (-1) 1)
    let out := (gridPoints g).flatMap fun pt =>
      v3List (fieldFdFloat m q phi pt.1 pt.2.1 pt.2.2)
    rows := rows.push ([natF g.1, natF g.2.1, natF g.2.2, natF q] ++ boxList rect a b c ++
      phi.toList ++ out)
  return rows.toList

def interpolateRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 15 do
    let g ← drawSmallGrid
    let p := (← below 5) + 1
    let (rect, a, b, c) ← drawBox
    let m := meshFloat (boxM3 a b c) g
    let atoms ← drawAtomsIn a b c 6
    let mut ex := #[]
    let mut ey := #[]
    let mut ez := #[]
    for _ in gridPoints g do
      ex := ex.push (← uniform (-5) 5)
      ey := ey.push (← uniform (-5) 5)
      ez := ez.push (← uniform (-5) 5)
    let forces := atoms.flatMap fun at_ => v3List (interpolateFloat m p ex ey ez at_)
    rows := rows.push ([natF g.1, natF g.2.1, natF g.2.2, natF p] ++ boxList rect a b c ++
      atomsList atoms ++ ex.toList ++ ey.toList ++ ez.toList ++ forces)
  return rows.toList

/-- The input grid of a `dft` row: `sin(0.37 i + seed)`, `cos(0.91 i + 2 seed)`. -/
def dftInput (n : ℕ) (seed : Float) : Array Float × Array Float :=
  ((Array.range n).map fun i => Float.sin (0.37 * i.toFloat + seed),
    (Array.range n).map fun i => Float.cos (0.91 * i.toFloat + 2 * seed))

def dftRows : List (List Float) :=
  [((4, 4, 4), 1.0, -1.0), ((4, 4, 4), 2.0, 1.0), ((8, 8, 8), 1.0, -1.0), ((8, 8, 8), 3.0, 1.0),
    ((8, 4, 2), 4.0, -1.0), ((16, 16, 16), 1.0, -1.0)].map fun (g, seed, sgn) =>
    let n := g.1 * g.2.1 * g.2.2
    let (re, im) := dftInput n seed
    let (ore, oim) := dftFloat g re im sgn
    [natF g.1, natF g.2.1, natF g.2.2, seed, sgn] ++ ore.toList ++ oim.toList

def p3mEnergyRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 8 do
    let k ← drawICase [8] 1
    let (inp, vol) := k.inputs
    let atoms ← drawAtomsIn k.a k.b k.c 12
    let e0 ← pick eps0Invs
    let (gh, _, _) := influenceFunctionFloat inp vol
    let rho := assignChargesFloat inp.mesh k.order atoms
    let (re, im) := dftFloat k.g rho (Array.replicate rho.size 0) (-1)
    let pts := (List.range rho.size).map fun i =>
      (gh.getD i 0, (NewenSpec.LatticeSum.Sym.zero : Sym Float), re.getD i 0, im.getD i 0)
    let (e, _) := kspaceEnergyFloat pts false inp.mesh.cellVolume e0
    rows := rows.push (k.list ++ [e0] ++ atomsList atoms ++ [e])
  return rows.toList

structure P3MRows where
  weights : List (List Float)
  selfterm : List (List Float)
  sinc : List (List Float)
  mesh : List (List Float)
  influence : List (List Float)
  quality : List (List Float)
  boxCorrection : List (List Float)
  forceError : List (List Float)
  influenceFunction : List (List Float)
  assign : List (List Float)
  assignSquared : List (List Float)
  kspace : List (List Float)
  a2 : List (List Float)
  fieldIk : List (List Float)
  fieldFd : List (List Float)
  interpolate : List (List Float)
  p3mEnergy : List (List Float)

def p3mRows : P3MRows :=
  (do return ⟨← weightsRows, ← selftermRows, ← sincRows, ← meshRows, ← influenceRows,
      ← qualityRows, ← boxCorrectionRows, ← forceErrorRows, ← influenceFunctionRows,
      ← assignRows false, ← assignRows true, ← p3mKspaceRows, ← a2Rows, ← fieldIkRows,
      ← fieldFdRows, ← interpolateRows, ← p3mEnergyRows⟩ : Gen _).run' 0x6e6577656e2d7033

def p3mParts (r : P3MRows) : List (String × List (List Float)) :=
  [("assignment_1d", assignmentRows), ("weights", r.weights), ("selfterm_fp", r.selfterm),
    ("fd_coefficients", fdRows), ("sinc_power", r.sinc), ("mesh", r.mesh),
    ("influence_at", r.influence), ("quality_at", r.quality), ("box_correction", r.boxCorrection),
    ("force_error", r.forceError), ("influence_function", r.influenceFunction),
    ("assign_charges", r.assign), ("assign_squared_charges", r.assignSquared),
    ("kspace_energy", r.kspace), ("a2_tilde", r.a2), ("field_ik", r.fieldIk),
    ("field_fd", r.fieldFd), ("interpolate_forces", r.interpolate), ("dft", dftRows),
    ("p3m_energy", r.p3mEnergy)]

def p3mJson (parts : List (String × List (List Float))) : String :=
  let part (p : String × List (List Float)) : String :=
    s!"  \"{p.1}\": [\n" ++ ",\n".intercalate (p.2.map hexRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"p3m.assignment_1d, p3m.nearest_grid_point, p3m.assignment_weights, p3m.selfterm_fp, " ++
  "p3m.fd_coefficients, p3m.sinc_power, p3m.mesh, p3m.influence_at, p3m.quality_at, " ++
  "p3m.influence_box_correction, p3m.force_error, p3m.assign_charges, " ++
  "p3m.assign_squared_charges, p3m.kspace_energy, p3m.potential_hat, p3m.a2_tilde, p3m.field_ik, " ++
  "p3m.field_fd, p3m.interpolate_forces (bit for bit); p3m.fft3, p3m.kspace (within a tolerance)\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.P3M, the *Float functions; SplitMix64, fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits; counts, orders, flags, " ++
  "nshape and integer coordinates are floats; a vector is x, y, z; a matrix row by row; a symmetric " ++
  "tensor xx, xy, xz, yy, yz, zz; a box is rectangular, a, b, c; grids in storage order z + Nz (y + Ny x)\",\n" ++
  "  \"layout_assignment_1d\": \"p, x, w_p(x)\",\n" ++
  "  \"layout_weights\": \"p, u, nearest grid point, five weights (unused zero)\",\n" ++
  "  \"layout_selfterm_fp\": \"p, |n|, s, f_p\",\n" ++
  "  \"layout_fd_coefficients\": \"q, five coefficients\",\n" ++
  "  \"layout_sinc_power\": \"p, x, P_p(x)\",\n" ++
  "  \"layout_mesh\": \"box, nx, ny, nz, H, H_inv, cell volume\",\n" ++
  "  \"layout_influence_at\": \"case (nx, ny, nz, order, nfdord, nalias, nshape, width, box), l, " ++
  "with_gamma, ghat, gammahat, q (reciprocal basis and mesh from reciprocal_box and Mesh::new of the box)\",\n" ++
  "  \"layout_quality_at\": \"case, l, ghat, the two terms\",\n" ++
  "  \"layout_box_correction\": \"ghat, gammahat, tl0_inv, delta_tl, corrected ghat\",\n" ++
  "  \"layout_force_error\": \"q, sum_q2, volume, n_atoms, fpe, error\",\n" ++
  "  \"layout_influence_function\": \"case, N, ghat per point, gammahat per point, Q\",\n" ++
  "  \"layout_assign_charges\": \"nx, ny, nz, order, box, n_atoms, per atom r, q; grid\",\n" ++
  "  \"layout_assign_squared_charges\": \"as assign_charges\",\n" ++
  "  \"layout_kspace_energy\": \"N, with_gamma, cell volume, eps0_inv, per point ghat, gammahat, re, im; " ++
  "energy, virial, per point the potential_hat re, im\",\n" ++
  "  \"layout_a2_tilde\": \"N, with_gamma, sum_q2, per point ghat, gammahat, re R; a2_tilde, derivative\",\n" ++
  "  \"layout_field_ik\": \"nx, ny, nz, box, per point phi re, im; per point Re E, Im E\",\n" ++
  "  \"layout_field_fd\": \"nx, ny, nz, q, box, per point phi; per point E\",\n" ++
  "  \"layout_interpolate_forces\": \"nx, ny, nz, order, box, atoms, Ex, Ey, Ez grids; per atom force\",\n" ++
  "  \"layout_dft\": \"nx, ny, nz, seed, sign of the exponent; re, im of the transform (input sin(0.37 i + seed), " ++
  "cos(0.91 i + 2 seed))\",\n" ++
  "  \"layout_p3m_energy\": \"case (grid 8^3), eps0_inv, atoms; k-space energy through dftFloat\",\n" ++
  ",\n".intercalate (parts.map part) ++ "\n}\n"

end P3MRows

/-! ## SPME: the per-point kernels and the grid kernels, bit for bit

`vectors/spme.json`: `M_n` and `M_n'` by the recursion for `n = 2 … 12` on a grid of `u` (the
knots and their neighbouring doubles included), the spline table of `spline_weights` for the same
orders at random and special `u`, the Euler factors of even orders on axes of 4 … 32 points, the
influence function `(Ĝ, Γ̂)` at single wave vectors (grids 4³ and 8³, orders 4 … 12, four shapes,
widths 0.4 and 1.4 nm, cubic, rectangular and triclinic boxes) and on whole 4³ and 8³ grids, the
charge and squared-charge assignment and the forces on random convolution grids, on 4³ grids. All
of these are compared bit for bit; `spme_energy` (the whole k-space energy on 8³ grids through
`dftFloat`) within a tolerance. -/

section SPMERows

open NewenTheory.Interactions.LatticeSum (Shape TPOrder)
open NewenSpec.LatticeSum (Sym reciprocalBoxFloat piFloat)

abbrev SV3 := NewenSpec.LatticeSum.V3
abbrev SM3 := NewenSpec.P3M.M3

def spmeOrders : List ℕ := [4, 6, 8, 10, 12]

def spmeShapes : List Shape := [.gaussian, .poly .o0, .poly .o3, .poly .o7]

def spmeM3List (m : SM3 Float) : List Float :=
  v3List m.r0 ++ v3List m.r1 ++ v3List m.r2

def spmeBoxList (rect : Bool) (a b c : SV3 Float) : List Float :=
  [if rect then 1 else 0] ++ v3List a ++ v3List b ++ v3List c

def spmeIntF (z : ℤ) : Float := Float.ofInt z

/-- The doubles at and on either side of `x`. -/
def spmeNeighbours (x : Float) : List Float :=
  if x == 0 then [0, Float.ofBits 1, -Float.ofBits 1]
  else [x, Float.ofBits (x.toBits - 1), Float.ofBits (x.toBits + 1)]

def bsplineRows : List (List Float) :=
  (List.range 11).flatMap fun i =>
    let n := i + 2
    let grid := (List.range (8 * (n + 2) + 1)).map fun j => j.toFloat / 8 - 1
    let knots := (List.range (n + 1)).flatMap fun k => spmeNeighbours k.toFloat
    (grid ++ knots).map fun u =>
      let d := NewenSpec.SPME.bsplineDerivFloat n u
      [natF n, u, NewenSpec.SPME.bsplineFloat n u, d.1, d.2]

def splineWeightRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for i in List.range 11 do
    let n := i + 2
    let mut us : List Float := [0, 0.5, 1, 3, 31, 31.999999, -0.25, -1, -3.5] ++
      spmeNeighbours 7 ++ spmeNeighbours 16
    for _ in List.range 16 do
      us := us ++ [← uniform 0 32]
    for _ in List.range 4 do
      us := us ++ [← uniform (-4) 0]
    for u in us do
      let (t, d) := NewenSpec.SPME.splineWeightsFloat n u
      rows := rows.push ([natF n, u] ++ t ++ d)
  return rows.toList

def eulerRows : List (List Float) :=
  spmeOrders.flatMap fun n => [4, 8, 16, 32].flatMap fun k =>
    (List.range k).map fun i =>
      let m := NewenSpec.P3M.waveIndex i k
      [natF n, natF k, spmeIntF m, NewenSpec.SPME.eulerFloat n k m]

/-- One SPME influence case: the grid, the order, the shape, the box. -/
structure SCase where
  g : ℕ × ℕ × ℕ
  n : ℕ
  shape : Shape
  width : Float
  rect : Bool
  a : SV3 Float
  b : SV3 Float
  c : SV3 Float

def SCase.recip (k : SCase) : SM3 Float :=
  let (_, u, v, w) := reciprocalBoxFloat k.a k.b k.c k.rect
  ⟨u, v, w⟩

def SCase.list (k : SCase) : List Float :=
  [natF k.g.1, natF k.g.2.1, natF k.g.2.2, natF k.n, nshapeFloat k.shape, k.width] ++
    spmeBoxList k.rect k.a k.b k.c

def drawSCase (grids : List ℕ) : Gen SCase := do
  let g := grids.getD (← below grids.length) 4
  let (rect, a, b, c) ← drawBox
  return ⟨(g, g, g), spmeOrders.getD (← below spmeOrders.length) 4,
    spmeShapes.getD (← below spmeShapes.length) .gaussian, ← pick [0.4, 1.4], rect, a, b, c⟩

def spmeInfluenceRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 150 do
    let k ← drawSCase [4, 8]
    let l ← drawL k.g
    let (g, gam) := NewenSpec.SPME.influenceAtFloat l k.g k.n k.shape k.width k.recip
    rows := rows.push (k.list ++ lList l ++ [g] ++ symList gam)
  return rows.toList

def spmeInfluenceFunctionRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for i in List.range 6 do
    let k ← drawSCase [if i < 4 then 4 else 8]
    let (g, gam) := NewenSpec.SPME.influenceFunctionFloat k.g k.n k.shape k.width k.recip
    rows := rows.push (k.list ++ g.toList ++ gam.toList.flatMap symList)
  return rows.toList

def spmeAssignRows (squared : Bool) : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 20 do
    let g ← drawSmallGrid
    let n := spmeOrders.getD (← below spmeOrders.length) 4
    let (rect, a, b, c) ← drawBox
    let atoms ← drawAtomsIn a b c 6
    let m := NewenSpec.P3M.meshFloat ⟨a, b, c⟩ g
    let out := if squared then NewenSpec.SPME.assignSquaredFloat m n atoms
      else NewenSpec.SPME.assignFloat m n atoms
    rows := rows.push ([natF g.1, natF g.2.1, natF g.2.2, natF n] ++ spmeBoxList rect a b c ++
      atomsList atoms ++ out.toList)
  return rows.toList

def spmeForceRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 20 do
    let g ← drawSmallGrid
    let n := spmeOrders.getD (← below spmeOrders.length) 4
    let (rect, a, b, c) ← drawBox
    let atoms ← drawAtomsIn a b c 6
    let e0 ← pick eps0Invs
    let m := NewenSpec.P3M.meshFloat ⟨a, b, c⟩ g
    let mut conv := #[]
    for _ in gridPoints g do
      conv := conv.push (← uniform (-1) 1)
    let forces := atoms.flatMap fun at_ => v3List (NewenSpec.SPME.forceFloat m n conv e0 at_)
    rows := rows.push ([natF g.1, natF g.2.1, natF g.2.2, natF n] ++ spmeBoxList rect a b c ++
      atomsList atoms ++ [e0] ++ conv.toList ++ forces)
  return rows.toList

/-- The whole k-space energy of SPME on an 8³ grid through the direct transform `dftFloat`
(compared within a tolerance: the Rust uses PhastFT). -/
def spmeEnergyRows : Gen (List (List Float)) := do
  let mut rows := #[]
  for _ in List.range 6 do
    let k ← drawSCase [8]
    let atoms ← drawAtomsIn k.a k.b k.c 12
    let e0 ← pick eps0Invs
    let (gh, _) := NewenSpec.SPME.influenceFunctionFloat k.g k.n k.shape k.width k.recip
    let m := NewenSpec.P3M.meshFloat ⟨k.a, k.b, k.c⟩ k.g
    let rho := NewenSpec.SPME.assignFloat m k.n atoms
    let (re, im) := NewenSpec.P3M.dftFloat k.g rho (Array.replicate rho.size 0) (-1)
    let pts := (List.range rho.size).map fun i =>
      (gh.getD i 0, (NewenSpec.LatticeSum.Sym.zero : Sym Float), re.getD i 0, im.getD i 0)
    let (e, _) := NewenSpec.P3M.kspaceEnergyFloat pts false m.cellVolume e0
    rows := rows.push (k.list ++ [e0] ++ atomsList atoms ++ [e])
  return rows.toList

structure SPMERows where
  energy : List (List Float)
  weights : List (List Float)
  influence : List (List Float)
  influenceFunction : List (List Float)
  assign : List (List Float)
  assignSquared : List (List Float)
  forces : List (List Float)

def spmeRows : SPMERows :=
  (do return ⟨← spmeEnergyRows, ← splineWeightRows, ← spmeInfluenceRows, ← spmeInfluenceFunctionRows,
      ← spmeAssignRows false, ← spmeAssignRows true, ← spmeForceRows⟩ : Gen _).run'
    0x6e6577656e2d7370

def spmeParts (r : SPMERows) : List (String × List (List Float)) :=
  [("bspline", bsplineRows), ("spline_weights", r.weights), ("euler_spline_factor", eulerRows),
    ("influence_at", r.influence), ("influence_function", r.influenceFunction),
    ("assign", r.assign), ("assign_squared_charges", r.assignSquared), ("forces", r.forces),
    ("spme_energy", r.energy)]

def spmeJson (parts : List (String × List (List Float))) : String :=
  let part (p : String × List (List Float)) : String :=
    s!"  \"{p.1}\": [\n" ++ ",\n".intercalate (p.2.map hexRow) ++ "\n  ]"
  "{\n" ++
  "  \"spec\": \"spme.bspline, spme.bspline_and_deriv, spme.spline_weights, " ++
  "spme.euler_spline_factor, spme.influence_at, spme.assign, spme.assign_squared_charges, " ++
  "spme.forces (bit for bit); spme.kspace (within a tolerance)\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (NewenSpec.SPME, the *Float functions; SplitMix64, fixed seed)\",\n" ++
  "  \"encoding\": \"every float is the hex string of its IEEE 754 binary64 bits; counts, orders, " ++
  "nshape and integer coordinates are floats; a vector is x, y, z; a symmetric tensor xx, xy, xz, " ++
  "yy, yz, zz; a box is rectangular, a, b, c; grids in storage order z + Nz (y + Ny x)\",\n" ++
  "  \"layout_bspline\": \"n, u, bspline(n, u), bspline_and_deriv(n, u)\",\n" ++
  "  \"layout_spline_weights\": \"n, u, n weights, n derivatives\",\n" ++
  "  \"layout_euler_spline_factor\": \"n, K, m, |b(m)|^2\",\n" ++
  "  \"layout_influence_at\": \"case (nx, ny, nz, n, nshape, width, box), l, ghat, gammahat " ++
  "(reciprocal basis from reciprocal_box of the box)\",\n" ++
  "  \"layout_influence_function\": \"case; ghat per point, gammahat per point\",\n" ++
  "  \"layout_assign\": \"nx, ny, nz, n, box, n_atoms, per atom r, q; grid\",\n" ++
  "  \"layout_assign_squared_charges\": \"as assign\",\n" ++
  "  \"layout_forces\": \"nx, ny, nz, n, box, atoms, eps0_inv, convolution grid; per atom force\",\n" ++
  "  \"layout_spme_energy\": \"case (grid 8^3), eps0_inv, atoms; k-space energy through dftFloat\",\n" ++
  ",\n".intercalate (parts.map part) ++ "\n}\n"

end SPMERows

def main : IO Unit := do
  IO.FS.createDirAll "vectors"
  IO.FS.writeFile "vectors/schedule.json" scheduleJson
  IO.println s!"wrote vectors/schedule.json ({grid} x {grid})"
  IO.FS.writeFile "vectors/dof.json" dofJson
  IO.println s!"wrote vectors/dof.json ({(dofGrid + 1) ^ 3} rows x {ndfmins.length})"
  IO.FS.writeFile "vectors/bath.json" bathJson
  IO.println s!"wrote vectors/bath.json ({bathSystems.length * 6} rows)"
  IO.FS.writeFile "vectors/schedule2.json" schedule2Json
  IO.println (s!"wrote vectors/schedule2.json ({dueRows.length} due, {framesInRows.length} " ++
    s!"frames_in, {stopRows.length} minimisation_stop, {finalRows.length} final_step, " ++
    s!"{comRows.length} com_removal rows)")
  IO.FS.writeFile "vectors/bath_empty.json" bathEmptyJson
  IO.println s!"wrote vectors/bath_empty.json ({halfLists.length} rows)"
  IO.FS.writeFile "vectors/topology.json" topologyJson
  IO.println s!"wrote vectors/topology.json ({topologyRows.length} rows)"
  -- the literals the models use must be the exact doubles 0, 1, 2, 6, 12 (and 0.5)
  unless (0 : Float).toBits == 0 && (1 : Float).toBits == 0x3FF0000000000000 &&
      (2 : Float).toBits == 0x4000000000000000 && (6 : Float).toBits == 0x4018000000000000 &&
      (12 : Float).toBits == 0x4028000000000000 && NewenSpec.Nonbonded.halfFloat == 1 / 2 &&
      NewenSpec.Bonded.halfFloat == 1 / 2 && quarterFloat == 1 / 4 &&
      (3 : Float).toBits == 0x4008000000000000 && NewenSpec.Dynamics.halfFloat == 1 / 2 &&
      thirdFloat == 1 / 3 do
    throw (IO.userError "Float literals 0, 1, 2, 3, 6, 12 are not exact")
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
  let (leap, kin, temp, pres, baro) := dynamicsRows
  IO.FS.writeFile "vectors/dynamics.json" (dynamicsJson leap kin temp pres baro)
  IO.println (s!"wrote vectors/dynamics.json ({leap.length} leap_frog, {kin.length} kinetic, " ++
    s!"{temp.length} temperature, {pres.length} pressure, {baro.length} barostat rows)")
  unless pointFiveFloat == 1 / 2 && (3 : Float).toBits == 0x4008000000000000 do
    throw (IO.userError "Float literals 0.5, 3 are not exact")
  match periodicityRows with
  | .error e => throw (IO.userError e)
  | .ok parts =>
    IO.FS.writeFile "vectors/periodicity.json" (periodicityJson parts)
    IO.println (s!"wrote vectors/periodicity.json (" ++
      ", ".intercalate (parts.map fun p => s!"{p.2.length} {p.1}") ++ " rows)")
  let sc := softCoreRows
  IO.FS.writeFile "vectors/softcore.json" (softCoreJson sc)
  IO.println (s!"wrote vectors/softcore.json ({sc.lambdaParams.length} lambda_params, " ++
    s!"{sc.pair.length} perturbed_lj_crf_pair, {sc.rfSoft.length} rf_soft_interaction, " ++
    s!"{sc.rfExcluded.length} perturbed_rf_excluded_pair, {sc.rfSelf.length} " ++
    s!"perturbed_rf_self_term, {sc.ljCrfExt.length} lj_crf_soft_interaction_ext, " ++
    s!"{sc.rfExt.length} rf_soft_interaction_ext rows)")
  let pb := perturbedBondedRows
  IO.FS.writeFile "vectors/perturbed_bonded.json" (perturbedBondedJson pb)
  IO.println (s!"wrote vectors/perturbed_bonded.json ({pb.chebyshev.length} chebyshev, " ++
    s!"{pb.quartic.length} quartic, {pb.angle.length} angle, {pb.improper.length} improper, " ++
    s!"{pb.dihedral.length} dihedral, {pb.soft.length} soft rows)")
  unless lsLiterals.all (fun n => (Float.ofNat n).toUInt64.toNat == n &&
        ((n : Float) : Float).toBits == (Float.ofNat n).toBits) &&
      (1 / 2 : Float) == 0.5 && (1 / 128 : Float) == 0.0078125 &&
      NewenSpec.LatticeSum.piFloat == 3.141592653589793 do
    throw (IO.userError "lattice-sum literals are not the expected doubles")
  let ls := latticeSumRows
  IO.FS.writeFile "vectors/lattice_sum.json" (latticeSumJson ls)
  IO.println (s!"wrote vectors/lattice_sum.json ({ls.switch.length} shape_switch, " ++
    s!"{ls.fourier.length} shape_fourier, {ls.pair.length} lj_ls_pair, " ++
    s!"{ls.excluded.length} ls_excluded_pair, {ls.selfConstants.length} self_constants rows)")
  unless (List.range 41).all (fun i =>
      let z : Int := (i : Int) - 20
      (Float.ofInt z).toBits == (if z < 0 then -Float.ofNat z.natAbs else Float.ofNat z.natAbs).toBits &&
        (Float.ofNat z.natAbs).toUInt64.toNat == z.natAbs) do
    throw (IO.userError "Float.ofInt is not exact on the integer triples")
  let ew := ewaldRows
  IO.FS.writeFile "vectors/ewald.json" (ewaldJson ew)
  IO.println (s!"wrote vectors/ewald.json ({ew.reciprocal.length} reciprocal_box, " ++
    s!"{ew.kvectors.length} k_vectors, {ew.kspace.length} kspace, {ew.a2Cubic.length} a2_cubic, " ++
    s!"{ew.a2Shells.length} a2_shells, {ew.chargeSums.length} charge_sums, " ++
    s!"{ew.selfA.length} self_and_a_terms, {ew.selfAVirial.length} self_and_a_virials, " ++
    s!"{ew.surface.length} surface_term rows)")
  unless NewenSpec.P3M.epsFloat == 0.000000000001 && (1 / 10 : Float) == 0.1 &&
      [3, 6, 8, 16, 32, 35, 63, 64, 96, 126, 192, 384, 576, 1152, 2304, 18432, 36864, 73728,
        147456, 32227, 29440].all (fun n => (Float.ofNat n).toUInt64.toNat == n) do
    throw (IO.userError "P3M literals are not the expected doubles")
  let parts := p3mParts p3mRows
  IO.FS.writeFile "vectors/p3m.json" (p3mJson parts)
  IO.println (s!"wrote vectors/p3m.json (" ++
    ", ".intercalate (parts.map fun p => s!"{p.2.length} {p.1}") ++ " rows)")
  let sparts := spmeParts spmeRows
  IO.FS.writeFile "vectors/spme.json" (spmeJson sparts)
  IO.println (s!"wrote vectors/spme.json (" ++
    ", ".intercalate (sparts.map fun p => s!"{p.2.length} {p.1}") ++ " rows)")
