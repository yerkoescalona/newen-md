import NewenSpec.Schedule
import NewenSpec.Thermo.Dof
import NewenSpec.Thermo.Bath

/-! `lake exe spec-vectors`: evaluate the Lean models on the grid the Rust tests sweep and write
them to `vectors/`, so `cargo test` can compare the Rust functions against Lean without Lean. -/

open NewenSpec.Schedule
open NewenSpec.Thermo.Dof
open NewenSpec.Thermo.Bath

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

def main : IO Unit := do
  IO.FS.createDirAll "vectors"
  IO.FS.writeFile "vectors/schedule.json" scheduleJson
  IO.println s!"wrote vectors/schedule.json ({grid} x {grid})"
  IO.FS.writeFile "vectors/dof.json" dofJson
  IO.println s!"wrote vectors/dof.json ({(dofGrid + 1) ^ 3} rows x {ndfmins.length})"
  IO.FS.writeFile "vectors/bath.json" bathJson
  IO.println s!"wrote vectors/bath.json ({bathSystems.length * 6} rows)"
