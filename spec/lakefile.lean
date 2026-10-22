import Lake
open Lake DSL

-- The Aeneas Lean library: the `Result` monad and the machine integers (`Std.Usize`, ...) the
-- extracted Rust is written in. It brings Mathlib along and fixes the toolchain
-- (`lean-toolchain` must equal the one in Aeneas' `backends/lean`). Keep the commit equal to
-- the Aeneas binary `scripts/extract.sh` runs (`aeneas -version`).
require aeneas from git
  "https://github.com/AeneasVerif/aeneas" @ "e62dd83d133d9043e7d3628096838b7430f36565" / "backends/lean"

package «newen-spec»

/-- The physics and bookkeeping of GROMOS, organised like the physics. Imports neither
`NewenSpec` nor Aeneas. -/
@[default_target] lean_lib NewenTheory

/-- What newen computes: models of `crates/newen-logic`, bridges to `NewenTheory`, and the
Aeneas refinement proofs. -/
@[default_target] lean_lib NewenSpec

/-- Writes `vectors/schedule.json`, `vectors/dof.json`, `vectors/bath.json` and
`vectors/topology.json` (the Lean models on the grids the Rust tests sweep),
`vectors/thermostat.json`, `vectors/nonbonded.json`, `vectors/bonded.json` and
`vectors/dynamics.json` (the thermostat, LJ + reaction-field, bonded, and leap-frog,
kinetic-energy, temperature, pressure and barostat models at `Float` on random inputs, as bits). -/
lean_exe «spec-vectors» where
  root := `SpecVectors
