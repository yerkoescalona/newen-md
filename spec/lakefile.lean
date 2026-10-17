import Lake
open Lake DSL

-- The Aeneas Lean library: the `Result` monad and the machine integers (`Std.Usize`, ...) the
-- extracted Rust is written in. It brings Mathlib along and fixes the toolchain
-- (`lean-toolchain` must equal the one in Aeneas' `backends/lean`). Keep the commit equal to
-- the Aeneas binary `scripts/extract.sh` runs (`aeneas -version`).
require aeneas from git
  "https://github.com/AeneasVerif/aeneas" @ "e62dd83d133d9043e7d3628096838b7430f36565" / "backends/lean"

package «gromos-spec»

/-- The physics and bookkeeping of GROMOS, organised like the physics. Imports neither
`GromosSpec` nor Aeneas. -/
@[default_target] lean_lib GromosTheory

/-- What gromos-rs computes: models of `crates/gromos-logic`, bridges to `GromosTheory`, and the
Aeneas refinement proofs. -/
@[default_target] lean_lib GromosSpec

/-- Writes `vectors/schedule.json`: the Lean model on the grid the Rust test sweeps. -/
lean_exe «spec-vectors» where
  root := `SpecVectors
