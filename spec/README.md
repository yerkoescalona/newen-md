# spec/ — Lean specifications of newen logic

Full explanation (why, how the pieces fit, guarantees, adding a function): [docs/src/verification/index.md](../docs/src/verification/index.md).

Verified so far: `schedule::output_misses_longrange` (the pilot, explained below) and
`thermo::dof::total_dof` (see `spec/SPEC.toml`). The pilot: `output_misses_longrange(nsnb, ntwg)` in
`crates/newen-logic/src/schedule.rs`. It answers: with a twin-range pairlist updated every
`nsnb` steps and free-energy output every `ntwg` steps, does some output step reuse long-range
terms from a pairlist step that was *not* an output step? (Under PRECALCLAM gromosXX then writes
λ_k rows without the long-range part; `newen-run` warns about it.)

The meaning lives in `NewenTheory/` (physics, organised like the physics; it never imports the
code's spec), the model and the bridge theorem in `NewenSpec/` (organised like the crate). Three
statements are checked, each by a different tool, and together they close a loop:

```
             meaning  (Misses: "∃ output step whose pairlist step is not an output step")
              ▲    ▲
   Layer 1    │    │  Layer 3
   Lean proof │    │  exhaustive test, 1 ≤ nsnb, ntwg ≤ 200
              │    │
   model ◄────┼── Rust function
 (missesModel)  Layer 2: Aeneas refinement proof
```

| Layer | Claim | Where | Run | Guarantees | Does not guarantee |
|---|---|---|---|---|---|
| 1. Lean theorem | model = meaning | `NewenTheory/Dynamics/Schedule.lean` (`Misses`), `NewenSpec/Schedule.lean` (`missesModel_iff`) | `lake build` | for **all** naturals, the closed form "nsnb > 1, ntwg > 0, neither divides the other" is exactly `Misses` | that `Misses` is what gromosXX does: that is a reading of `perturbed_nonbonded_set.cc`, checked by humans and by the reference tests |
| 2. Aeneas refinement | Rust = model | `NewenSpec/Extracted/` (generated), `NewenSpec/Refinement.lean` | `scripts/extract.sh`, then `lake build` | for **all** `usize` inputs the Rust function never panics and returns `missesModel` | that Charon/Aeneas translate Rust faithfully (trusted), or anything about callers such as `prepare.rs` |
| 3. Differential test | Rust = meaning, Rust = model, on a grid | `crates/newen-logic/tests/schedule_spec.rs`, `vectors/schedule.json` | `cargo test -p newen-logic` | the compiled Rust agrees with a step-by-step simulation and with the Lean model on 40 000 inputs | anything outside the grid |

Layers 1 + 2 are the proof. Layer 3 needs no Lean and runs with every `cargo test`; it also tests the
*compiled* binary, which the proofs do not see, and its simulation is written independently of
the formula, so a wrong `pairlistStep` definition in Lean would show up there.

## Running

```sh
cd spec && lake build                    # layers 1 and 2: proofs (first run downloads Mathlib's cache)
cd spec && lake exe spec-vectors         # rewrite vectors/schedule.json from the Lean model
spec/scripts/extract.sh                  # regenerate NewenSpec/Extracted/ after editing the Rust
cargo test --release -p newen-logic     # layer 3
```

`lake build` must end with no `sorry` warning; `#print axioms` on the theorems lists only
`propext`, `Classical.choice`, `Quot.sound` (Lean's standard axioms).

**Toolchain.** Layer 1 uses core Lean only (`omega`, `Nat` lemmas);
`lean NewenTheory/Dynamics/Schedule.lean` checks the theory with no packages at all.
`NewenTheory` and `NewenSpec` are two libraries of one Lake project, on one toolchain. The
Aeneas library needs Mathlib and fixes the toolchain: `lean-toolchain` (v4.35.0) and the commit
in `lakefile.lean` must match the Aeneas binary.
**Install** (no root): download `aeneas-linux-x86_64.tar.gz` from the Aeneas GitHub release whose
tag ends in the commit pinned in `lakefile.lean` and unpack it into `~/.local/aeneas/aeneas` (it
bundles the matching `charon`); then `rustup toolchain install <channel in its rust-toolchain>
--profile minimal -c rustc-dev -c llvm-tools -c rust-src`. `extract.sh` refuses to run if the
binary and the pinned library differ.

## Adding a function

1. Write it in `crates/newen-logic` in the subset Aeneas translates (safe Rust, integers and
   booleans, no floats, simple loops); end its doc comment with `/// spec: <module>.<name>`.
2. The meaning (a `Prop` over all inputs) in `NewenTheory/<Area>/<Topic>.lean`, placed by its
   physics; in `NewenSpec/<Module>.lean` a model written like the Rust, and a theorem model ↔
   meaning.
3. `scripts/extract.sh`, then a theorem in `Refinement.lean`: extracted function `= ok (model ..)`
   (usually `step*` plus `scalar_tac`/`simp_all`).
4. A test in `crates/newen-logic/tests/` comparing against an independent brute-force meaning and,
   via `SpecVectors.lean`, against the model.
5. A `[[function]]` entry in `SPEC.toml`.
