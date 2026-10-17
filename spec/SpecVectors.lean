import GromosSpec.Schedule

/-! `lake exe spec-vectors`: evaluate the Lean models on the grid the Rust tests sweep and write
them to `vectors/`, so `cargo test` can compare the Rust functions against Lean without Lean. -/

open GromosSpec.Schedule

/-- Grid bound, equal to `GRID` in `crates/gromos-logic/tests/schedule_spec.rs`. -/
def grid : Nat := 200

/-- Row `nsnb` as a string of `0`/`1`, one character per `ntwg = 1..grid`. -/
def row (nsnb : Nat) : String :=
  String.ofList ((List.range grid).map fun i => if missesModel nsnb (i + 1) then '1' else '0')

def scheduleJson : String :=
  let rows := (List.range grid).map fun i => s!"    \"{row (i + 1)}\""
  "{\n" ++
  "  \"spec\": \"schedule.output_misses_longrange\",\n" ++
  "  \"generator\": \"lake exe spec-vectors (GromosSpec.Schedule.missesModel)\",\n" ++
  s!"  \"layout\": \"rows[nsnb - 1][ntwg - 1] for 1 <= nsnb, ntwg <= {grid}; 1 = true\",\n" ++
  "  \"rows\": [\n" ++ ",\n".intercalate rows ++ "\n  ]\n}\n"

def main : IO Unit := do
  IO.FS.createDirAll "vectors"
  IO.FS.writeFile "vectors/schedule.json" scheduleJson
  IO.println s!"wrote vectors/schedule.json ({grid} x {grid})"
