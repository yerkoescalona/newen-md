# Verifying newen against its theory

This chapter explains how parts of newen are checked against a mathematical specification
written in the Lean theorem prover, and how the Rust code is tied to that specification. You do
not need to know Lean or Aeneas to read it. By the end you should know why we do this, how the
pieces fit together, how to run them, what is and is not guaranteed, and how to add a function of
your own.

Everything described here lives in two places in the repository:

- `crates/newen-logic`: the Rust functions that carry a specification;
- `spec/`: the Lean project that states and proves what those functions mean (`spec/README.md` is
  the short operational summary of this chapter).

## Why

newen is a port of gromosXX. The obvious way to check a port is to run both programs on the
same input and compare the output, and newen does that on every reference system in
`crates/newen-md/tests/gromosXX_references/`. That test answers one question very well: *does
newen do what gromosXX does?*

It cannot answer a second question: *does gromosXX do what the equations say?* While porting we
found 39 places where gromosXX disagrees with the equations it was meant to implement. Few of them
are deep. Most are bookkeeping:

- arguments shifted by one position;
- a value written at the wrong step;
- items matched in the wrong order;
- a missing time step Δt in a running sum;
- a long-range term missing on some output steps.

Each of these produces plausible numbers. A comparison against gromosXX passes, because newen
faithfully reproduces the defect. A comparison against gromosXX only proves that we copy it,
defects included.

What catches this kind of error is a statement of what the code should do that is independent of
any code: a specification. The findings are collected, with reproducers and proposed fixes, in
the divergence dossier of the gromos-divergences repository. Each has a number (D1, D2, …); this
chapter uses one of them, D34, as its worked example.

Where gromosXX is a defect against the equations, newen computes the equations, and the reference
it is checked against comes from gromosXX with that defect patched. The defect itself is shown in
the dossier, not kept as an option in newen. Where the difference is a GROMOS convention rather
than a defect, newen keeps gromosXX's behaviour and warns when it matters. The specification is
what tells us where those places are.

## The idea in one picture

Every verified function has three descriptions, and every pair of neighbouring descriptions is
connected by a check done with a different tool:

```text
   MEANING            what the physics requires, as a Lean proposition (Prop)
   "there is an       stated over all inputs, with no algorithm in it
    output step
    that ..."
       ^
       |  Lean theorem: model <=> meaning, for all natural numbers
       |  (checked by the Lean kernel, `lake build`)
       v
   MODEL              a Lean function, written the way the Rust is written
   missesModel        (same branches, same arithmetic)
       ^
       |  Aeneas: the Rust function is translated into Lean (Charon + Aeneas),
       |  and a Lean refinement proof shows the translation returns the model
       |  for every machine input, and never panics (`lake build`)
       v
   RUST               the function newen-run actually calls
   output_misses_longrange

   and, closing the loop from the other side:

   RUST  <----  exhaustive / differential test (`cargo test`)  ---->  MEANING, MODEL
                compiled binary vs. a step-by-step simulation
                and vs. values the Lean model computed
```

The proofs cover all inputs but look at Lean code. The test looks at the compiled binary but only
at a finite grid of inputs. Together they leave very little room for a mistake to hide.

The Lean side is two libraries of one project. The meaning lives in `NewenTheory`
(`spec/NewenTheory/`): the physics and bookkeeping of a GROMOS run, organised like the physics
(this pilot's meaning is in `Dynamics/Schedule.lean`), and it never imports the code's
specification. The model and the bridge theorem linking it to the meaning live in `NewenSpec`
(`spec/NewenSpec/`), one file per module of `crates/newen-logic`, next to the refinement proofs.
The same definitions can then judge newen here and gromosXX in the divergence dossier.

## Two kinds of code, two techniques

Not all code can be verified the same way. newen splits it in two.

**Integer decision logic.** Which step writes output, how many degrees of freedom a system has,
which atom belongs to which temperature bath, which index in one list matches which in another.
These are small functions over integers and booleans. They are where most of the 39 defects live,
and they are exactly what current tools can prove. They go into `crates/newen-logic`: safe Rust,
no floating point, no dependencies, the subset that Charon and Aeneas translate into Lean. A Lean
proof then shows that the Rust function equals the model for every input.

**Floating-point kernels.** Lennard-Jones and reaction-field pair terms, thermostats. These compute
with `f64`. Proving that a floating-point kernel computes its real-number equation, within a
rounding bound, is a research problem: Boldo et al. (2013) needed a long paper and machine-checked
proofs of both the method error and the rounding error for a single C program solving the one-dimensional
wave equation. For these kernels we use the technique AWS uses for its Cedar authorization
language: write an executable model in Lean, prove properties of the model, and test the Rust
kernel against the model on many random inputs (differential random testing).

| | Integer decision logic | Floating-point kernels |
|---|---|---|
| Examples | output schedules, degrees of freedom, bath mapping, index matching | LJ, reaction field, Berendsen, Nosé–Hoover |
| Rust lives in | `crates/newen-logic` | the existing kernel crates |
| Lean has | meaning, model, theorem model ⇔ meaning | meaning, executable model, theorems about the model |
| Rust ⇔ model | proof, for every input (Charon + Aeneas + refinement theorem) | differential random testing |
| Status | working (one function, the pilot below) | working for the Berendsen and Nosé–Hoover thermostats (below); LJ and reaction field planned |

## Worked example: `output_misses_longrange`

The pilot function decides one question about a run's schedule. We go through it from the physics
to the last proof.

### The physics

GROMOS uses a twin-range cutoff (see [Pairlist & Neighbor Search](../algorithms/pairlist.md)).
Pairs within the short cutoff `RCUTP` are evaluated every step. Pairs between `RCUTP` and the long
cutoff `RCUTL` are evaluated only when the pairlist is rebuilt, every `NSNB` steps, and their
contribution is reused until the next rebuild.

In a free-energy run, the PRECALCLAM block asks gromosXX to compute, on every free-energy output
step (every `NTWG` steps), the energies at a list of extra λ values λ_k, for BAR and MBAR. gromosXX
adds the long-range part of those λ_k energies only if the pairlist step that produced it was
itself an output step (`perturbed_nonbonded_set.cc:86-96`, `perturbed_nonbonded_innerloop.cc:315`).
An output step that reuses a pairlist step which was *not* an output step writes λ_k rows without
their long-range part. This is D34 in the dossier.

Take `NSNB = 5`, `NTWG = 2`. Output step `s` reuses pairlist step `5 * floor(s / 5)`:

```text
step          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
pairlist      P              P              P              P              P
output        W     W     W     W     W     W     W     W     W     W     W
reuses        0     0     0     5     5    10    10    10    15    15    20
complete    yes   yes   yes    NO    NO   yes   yes   yes    NO    NO   yes
```

Steps 6, 8, 16 and 18 reuse pairlist steps 5 and 15, which are odd and so not output steps; their
rows miss the long-range term. In the dossier's methane-in-water reproducer this is about
1.25 kJ/mol on two of every five frames, and every BAR or MBAR estimate built from those rows is
biased. With `NTWG = 10` instead, every output step reuses an output pairlist step and nothing is
missed.

newen reproduces this behaviour (so its output matches gromosXX bit for bit) and warns before
the run. The warning asks one integer question: *for this `NSNB` and `NTWG`, does some output step
miss its long-range part?* That question is the pilot function.

### The Rust code

`crates/newen-logic/src/schedule.rs`:

```rust
/// Whether some free-energy output step (every `ntwg` steps) follows a pairlist step (every
/// `nsnb` steps) that was not itself an output step.
///
/// With a twin-range pairlist the long-range interactions are evaluated only on pairlist steps
/// and reused until the next one. Under PRECALCLAM gromosXX adds the long-range λ_k energies to
/// the written rows only when that pairlist step was also an output step, so such an output step
/// writes λ_k rows without their long-range part. That happens on some output step exactly when
/// `nsnb > 1`, `ntwg > 0`, and neither of `nsnb` and `ntwg` divides the other. `ntwg == 0` means
/// no output; `nsnb <= 1` means every step is a pairlist step, so nothing is reused.
///
/// spec: schedule.output_misses_longrange
// `%` rather than `is_multiple_of`: plain arithmetic keeps the Lean extraction self-explanatory.
#[allow(clippy::manual_is_multiple_of)]
pub fn output_misses_longrange(nsnb: usize, ntwg: usize) -> bool {
    if nsnb <= 1 || ntwg == 0 {
        return false;
    }
    ntwg % nsnb != 0 && nsnb % ntwg != 0
}
```

Three lines of logic and a doc comment making a claim: the answer is "yes" exactly when neither
number divides the other. The line `/// spec: schedule.output_misses_longrange` is the key that
connects this function to its entry in `spec/SPEC.toml`.

The caller in `crates/newen-run/src/prepare.rs` adds the conditions that are not integer logic
(PRECALCLAM on, twin range present) and leaves the schedule question to `newen-logic`:

```rust
pub fn precalclam_longrange_gaps(imd: &ImdParameters) -> Option<(usize, usize)> {
    let precalc = imd.ntg != 0 && imd.precalclam.is_some_and(|b| b.nr_lambdas > 0);
    let twin_range = imd.rcutp < imd.rcutl - 1e-10;
    let (nsnb, ntwg) = (imd.nsnb, imd.ntwg);
    (precalc && twin_range && newen_logic::schedule::output_misses_longrange(nsnb, ntwg))
        .then_some((nsnb, ntwg))
}
```

### The meaning, in Lean

`spec/NewenTheory/Dynamics/Schedule.lean` states what the question *means*, without any formula
and without any reference to newen:

```lean
/-- The pairlist step whose long-range interactions step `s` reuses: the last multiple of `nsnb`
at or before `s`, `nsnb * ⌊s / nsnb⌋`. -/
def pairlistStep (nsnb s : Nat) : Nat := nsnb * (s / nsnb)

/-- Every output step (a multiple of `ntwg`) reuses a pairlist step that was also an output step,
so every written row carries its long-range part. -/
def RowsComplete (nsnb ntwg : Nat) : Prop :=
  ∀ s, ntwg ∣ s → ntwg ∣ pairlistStep nsnb s

/-- Some output step reuses a pairlist step that was not an output step: that row misses its
long-range part. -/
def Misses (nsnb ntwg : Nat) : Prop :=
  ∃ s, ntwg ∣ s ∧ ¬ ntwg ∣ pairlistStep nsnb s
```

Read aloud: "`pairlistStep nsnb s` is `nsnb` times `s` divided by `nsnb`, rounded down", that is,
the last pairlist step at or before `s`. "`Misses nsnb ntwg` holds when there exists a step `s`
such that `ntwg` divides `s` (it is an output step) and `ntwg` does not divide the pairlist step
that `s` reuses (that pairlist step was not an output step)."

`RowsComplete` is the healthy case, read the same way: "for every step `s`, if `s` is an output
step then its pairlist step is an output step too." `Misses` is exactly its negation (the file
proves that as `misses_iff_not_complete`). `Misses` is what the Rust doc comment says the function
decides.

This is a direct transcription of the physics in the previous section. It says nothing about
divisibility of `nsnb` and `ntwg`, and it ranges over infinitely many steps, so you cannot run it.
It is the thing we want to be true.

### The model and the theorem `missesModel_iff`

Next comes the model, in `spec/NewenSpec/Schedule.lean`: a Lean function written branch for
branch like the Rust.

```lean
/-- Closed form of `output_misses_longrange`, written as the Rust code is:
`nsnb ≤ 1` or `ntwg = 0` first (no reuse, or no output), then "neither divides the other". -/
def missesModel (nsnb ntwg : Nat) : Bool :=
  if nsnb ≤ 1 then false
  else if ntwg = 0 then false
  else ntwg % nsnb != 0 && nsnb % ntwg != 0
```

And the theorem that connects model and meaning:

```lean
theorem missesModel_iff (nsnb ntwg : Nat) :
    missesModel nsnb ntwg = true ↔ Misses nsnb ntwg := by
```

Read aloud: "for every pair of natural numbers `nsnb` and `ntwg`, the model returns `true` if and
only if some output step misses its long-range part." There are no hypotheses: it holds for zero,
for one, for every size.

The proof (the lines after `by`) rests on a lemma of the theory, `rows_complete_iff`: for
positive `nsnb` and `ntwg`, every output step reuses an output pairlist step exactly when `nsnb`
divides `ntwg` or `ntwg` divides `nsnb`. When neither divides the other, the proof exhibits the
failing step: if `nsnb < ntwg` it is step `ntwg` itself; if `ntwg < nsnb` it is the first output
step after `nsnb`.

The small values needed care, and this is a typical place for specification work to pay off. The
divisibility argument only makes sense for positive periods, so the edge cases are proven
separately, against the meaning rather than the formula:

- `nsnb = 1`: every step is a pairlist step, so `pairlistStep 1 s = s` (the lemma
  `pairlistStep_one`) and every output step reuses itself. The formula alone would also say "no"
  here (1 divides everything), but only the meaning shows *why*: nothing is reused at all.
- `nsnb = 0`: in Lean `s / 0 = 0`, so the pairlist step is always step 0, which is an output step.
  In Rust, `ntwg % 0` would panic. The `nsnb <= 1` guard is what keeps the Rust from dividing by
  zero, and the theorem confirms that answering `false` there agrees with the meaning.
- `ntwg = 0`: the only multiple of 0 is 0, so the only "output step" is step 0, which is its own
  pairlist step. Again `false`, and again the guard prevents `nsnb % 0` in Rust.

So the guards in the Rust are not just defensive programming: each one is a case of the theorem.

### The extracted Lean

Charon reads the compiled Rust (its MIR, the compiler's mid-level representation) and Aeneas
turns it into a pure Lean function. `spec/scripts/extract.sh` writes the result to
`spec/NewenSpec/Extracted/NewenLogic.lean`; you never edit that file by hand.

```lean
def schedule.output_misses_longrange
  (nsnb : Std.Usize) (ntwg : Std.Usize) : Result Bool := do
  if nsnb <= 1#usize
  then ok false
  else
    if ntwg = 0#usize
    then ok false
    else
      let i ← ntwg % nsnb
      if i != 0#usize
      then let i1 ← nsnb % ntwg
           ok (i1 != 0#usize)
      else ok false
```

Line by line against the Rust:

| Rust | Extracted Lean | Why it looks like that |
|---|---|---|
| `nsnb: usize, ntwg: usize` | `(nsnb : Std.Usize) (ntwg : Std.Usize)` | `Usize` is a machine integer, `0 ≤ n < 2^64`, not an unbounded `Nat` |
| `-> bool` | `: Result Bool` | a Rust function can panic; `Result` is either `ok b` (normal return) or a failure (panic) |
| `if nsnb <= 1 \|\| ntwg == 0` | two nested `if`s | `\|\|` short-circuits, so in the compiled code it is control flow |
| `return false;` | `ok false` | a normal return of `false` |
| `ntwg % nsnb` | `let i ← ntwg % nsnb` | `%` can panic (division by zero), so it is a `Result` too, and `←` takes its value if it succeeded |
| `... != 0 && nsnb % ntwg != 0` | `if i != 0#usize then ... else ok false` | `&&` short-circuits: the second `%` runs only if the first test passed |

Nothing in the extracted function knows about pairlists. It is what the compiler made of the Rust,
in a language Lean can reason about.

### The refinement theorem `output_misses_longrange_iff`

`spec/NewenSpec/Refinement.lean` proves that the extracted function computes the model:

```lean
/-- **Rust = model**: for all `usize` inputs, `output_misses_longrange` returns `missesModel`. -/
theorem output_misses_longrange_eq (nsnb ntwg : Usize) :
    output_misses_longrange nsnb ntwg = ok (Schedule.missesModel nsnb.val ntwg.val) :=
  eq_ok_of_spec (output_misses_longrange_spec nsnb ntwg)
```

Read aloud: "for every pair of `usize` inputs, the Rust function returns normally (`ok`), and what
it returns is the model applied to the numeric values of the inputs." Returning `ok` for every
input means it never panics: neither `%` ever divides by zero.

Chaining this with `missesModel_iff` gives the statement we wanted from the start:

```lean
/-- **Rust = meaning**: the Rust function returns `true` exactly when some output step reuses a
pairlist step that was not an output step. -/
theorem output_misses_longrange_iff (nsnb ntwg : Usize) :
    output_misses_longrange nsnb ntwg = ok true ↔
      NewenTheory.Dynamics.Schedule.Misses nsnb.val ntwg.val := by
  rw [output_misses_longrange_eq, Result.ok.injEq, Schedule.missesModel_iff]
```

Read aloud: "the Rust function returns `true` exactly when the physical problem happens; it never
panics." The proof is one line of rewriting: replace the Rust function by the model
(`output_misses_longrange_eq`), strip the `ok` on both sides, replace the model by the meaning
(`missesModel_iff`).

The work happens in `output_misses_longrange_spec`, whose proof is mostly the Aeneas tactic
`step*`. It walks through the extracted function, splits each `if`, and before each `%` proves the
divisor is nonzero, using the facts the earlier branches established.

### The exhaustive test

`crates/newen-logic/tests/schedule_spec.rs` runs the *compiled* Rust on every
`1 ≤ nsnb, ntwg ≤ 200` (40 000 inputs) and compares it with two independent references:

- `misses_by_simulation`, which walks the MD steps one by one, the way a run does, remembering the
  last pairlist step and whether it was an output step. It does not use the closed form, so it is
  a second, executable transcription of the meaning.
- `spec/vectors/schedule.json`, the Lean model evaluated on the same grid by
  `lake exe spec-vectors` and checked in, so the test needs no Lean.

```rust
/// Walks steps `0..lcm(nsnb, ntwg)` (the schedule repeats after that): remember the last pairlist
/// step and whether it was an output step; report an output step that reuses one that was not.
fn misses_by_simulation(nsnb: usize, ntwg: usize) -> bool {
    let period = nsnb / gcd(nsnb, ntwg) * ntwg;
    let mut pairlist_was_output = true;
    for step in 0..period {
        let is_output = step % ntwg == 0;
        if step % nsnb == 0 {
            pairlist_was_output = is_output;
        }
        if is_output && !pairlist_was_output {
            return true;
        }
    }
    false
}
```

This test does what the proofs cannot: it sees the binary the compiler produced, so a translation
error in Charon or Aeneas would show up here. And it runs with every `cargo test`, for contributors
who never install Lean.

It does catch bugs. Replace the last line of the Rust function with the simpler-looking
`ntwg % nsnb != 0` (warn whenever `NTWG` is not a multiple of `NSNB`) and both tests fail,
listing `(nsnb, ntwg)` pairs such as `(5, 1)`: with output on every step, each pairlist step is an
output step and nothing is missed, but the simpler check would warn.

## Floating-point kernels: the thermostats

The Berendsen and Nosé–Hoover thermostats compute one scaling factor per temperature bath and
step. They are `f64` code, so they stay in `crates/newen-integrators` and Aeneas plays no part.
The three pieces follow LeanLJ (one definition, proven about over ℝ and executed over `Float`)
and Cedar (the Rust tested against the executable model). We go through Berendsen; Nosé–Hoover
has the same three pieces.

### The Rust kernel

`crates/newen-integrators/src/algorithms/berendsen_thermostat.rs`, the whole decision that
`BerendsenThermostat::apply` makes for one bath, in gromosXX's order of operations:

```rust
pub fn berendsen_scale(ekin: f64, dof: f64, temperature: f64, tau: f64, dt: f64) -> f64 {
    if tau < 0.0 {
        return 1.0;
    }
    let mut free_temp = if dof > 0.0 {
        2.0 * ekin / (dof * kB)
    } else {
        0.0
    };
    if free_temp < ZERO_TEMPERATURE {
        free_temp = temperature;
    }
    if free_temp < ZERO_TEMPERATURE {
        1.0
    } else if tau == 0.0 {
        (temperature / free_temp).sqrt()
    } else {
        (1.0 + dt / tau * (temperature / free_temp - 1.0)).sqrt()
    }
}
```

`ZERO_TEMPERATURE = 1e-12` is gromosXX's `math::epsilon`, the one threshold both thermostats use.

### The meaning, over ℝ

`spec/NewenTheory/Dynamics/Thermostats.lean` states what the factor is for. Berendsen et al.
(1984) let the temperature relax towards `T₀` with time constant `τ`; over one step that is
`λ² = 1 + dt/τ · (T₀/T − 1)`:

```lean
noncomputable def berendsenScale (T T₀ dt τ : ℝ) : ℝ :=
  Real.sqrt (1 + dt / τ * (T₀ / T - 1))

theorem berendsen_relax (hT : 0 < T) (hr : 0 ≤ 1 + dt / τ * (T₀ / T - 1)) :
    berendsenScale T T₀ dt τ ^ 2 * T = T + dt / τ * (T₀ - T)

theorem berendsen_between (hT : 0 < T) (hT₀ : 0 ≤ T₀) (hdt : 0 < dt) (hdτ : dt ≤ τ) :
    min T T₀ ≤ berendsenScale T T₀ dt τ ^ 2 * T ∧
      berendsenScale T T₀ dt τ ^ 2 * T ≤ max T T₀
```

Read aloud: the new temperature `λ² T` has moved the fraction `dt/τ` of the way to `T₀`, and with
`dt ≤ τ` it never overshoots. `berendsen_contraction` (the gap shrinks by `1 − dt/τ`) and
`berendsen_fixed` (a bath at `T₀` is left alone) complete it. For Nosé–Hoover, `nh_feedback`
says the friction grows exactly when the bath is hotter than `T₀`.

### One model, two number types

`spec/NewenSpec/Thermostat.lean` transcribes the Rust once, for any type with the operations the
kernel uses, and takes `k_B` and the threshold as parameters:

```lean
def berendsenModel (kB eps ekin dof temperature tau dt : α) : α :=
  if tau < 0 then 1
  else
    let freeTemp := freeTempModel α kB eps ekin dof temperature
    if freeTemp < eps then 1
    else if tau == 0 then Sqrt.sqrt (temperature / freeTemp)
    else Sqrt.sqrt (1 + dt / tau * (temperature / freeTemp - 1))
```

At `α = ℝ` it is proven equal to the meaning once the fallbacks are ruled out:

```lean
theorem berendsenModel_real_eq (hτ : 0 < tau) (hdof : 0 < dof)
    (hft : eps ≤ 2 * ekin / (dof * kB)) :
    berendsenModel ℝ kB eps ekin dof temperature tau dt =
      berendsenScale (2 * ekin / (dof * kB)) temperature dt tau
```

At `α = Float` the same definition runs. `lake exe spec-vectors` draws random inputs from a
fixed-seed generator (typical and random degrees of freedom, temperatures, `dt` and `τ`, with
`τ < 0`, `τ = 0`, `dof = 0` and bath temperatures from `1e-14` to `1e6` K) and writes each input and
each result as the 16 hex digits of its bits to `spec/vectors/thermostat.json`. Nosé–Hoover chains
of 1, 2, 3 and 5 links are iterated for 1–8 steps, so the state carried from one step to the next
is covered.

### The bit-exact test

`crates/newen-integrators/tests/thermostat_spec.rs` runs the compiled kernels on every row and
compares `to_bits()`: no tolerance. That is possible because IEEE 754 requires `+ − × ÷` and
`sqrt` to be correctly rounded, the exact result rounded to the nearest double. Rust's `f64` and
Lean's `Float` (a C `double`) both obey it, so the same operations on the same inputs in the same
order give the same bits. The model must therefore follow the Rust's order exactly: the test
fails if the Rust writes `dt * (T₀/T − 1) / τ` instead of `dt / τ * (T₀/T − 1)` (about one
Berendsen row in thirty then differs, by 1 to several thousand ulp). The constants are written as
bits too, and the test checks `kB.to_bits()` and `ZERO_TEMPERATURE.to_bits()` against the
literals in the Lean file, because Lean's decimal parser is not promised to round like Rust's.
One value has no bits in Lean: its `Float` has a single NaN, so a NaN matches any NaN.

Bit-exact comparison stops working when a kernel calls a function that IEEE 754 does not require
to be correctly rounded, such as `exp`, `log`, `pow` or `erfc` (Ewald sums need `erfc`):
Rust's and C's libraries may differ in the last bit. It also stops working when the two sides
sum in a different order, as a parallel reduction over pairs does. There a tolerance, justified by a rounding analysis, has to replace equality.

What is not guaranteed: how far the `Float` model is from the ℝ model (no rounding bound is
proved), and inputs outside the generated set.

## Where each mistake is caught

Lean does not run your code and look for errors. You write a statement ("the Rust function equals
the formula") and a proof. If the code or the formula is wrong, the proof can no longer be
completed and `lake build` stops at the step that fails. *Detection means a proof breaks* — or,
for the tests, a comparison fails. Each layer catches a different kind of mistake, at a different
moment:

| Level | Compares | Runs when | Catches | Example from this pilot |
|---|---|---|---|---|
| 0. Dossier | a formula as gromosXX transcribed it (written into Lean by hand) vs the formula of the paper or manual | when someone writes the transcription | wrong physics in someone else's code | the 39 rows D1–D39 in the gromos-divergences repository |
| 1. Meaning ⇔ model | `Misses` vs `missesModel` (`missesModel_iff`) | `lake build` | a formula that does not mean what you think | "warn when `NTWG % NSNB != 0`" looks right but is wrong for NTWG = 1 (every step is written, nothing is missed): the theorem cannot be proved, so the mistake is caught before any Rust exists |
| 2. Rust ⇔ model | the Aeneas translation of the Rust vs `missesModel`, for **every** pair of `usize` | `spec/scripts/extract.sh`, then `lake build` | transcription mistakes, possible panics | swapping the two `%` in one place, `nsnb <= 1` written as `nsnb < 1`, or deleting the `ntwg == 0` guard (Aeneas translates `%` by zero as `fail`, so "always returns `ok`" becomes unprovable) |
| 3. Binary ⇔ meaning | the compiled function vs a step-by-step simulation and the Lean vectors, on all pairs up to 200 | `cargo test`, no Lean needed | the same transcription mistakes, and mistranslation by Charon or Aeneas themselves | replacing the last line with `ntwg % nsnb != 0` fails both tests; the mismatches include `(5, 1)` |

A correct refactor also changes the extracted Lean. The proof may then need adjusting, but it can
be re-proved; a wrong change cannot. That asymmetry is the whole mechanism.

What no level catches: if `Misses` itself misdescribes what the physics or gromosXX does, all
levels stay green and agree with a wrong idea. Comparisons with gromosXX references and published
values remain the guard for that. The method turns "is the code right?" into "is the spec
right?", and the spec is ten lines you can read and argue about.

## Reading Lean in five minutes

You can read the files in `spec/` with this legend. It covers exactly the syntax that appears
there.

| Syntax | Meaning |
|---|---|
| `def f (a b : Nat) : Nat := ...` | a function definition: `f` takes two natural numbers and returns one |
| `theorem name (a : Nat) : statement := proof` | a claim, and its proof; Lean refuses to build if the proof is wrong |
| `Nat` | natural numbers 0, 1, 2, ..., unbounded |
| `Usize`, `1#usize` | the Rust machine type `usize`, and the literal 1 of that type; `n.val` is its value as a `Nat` |
| `Bool` | `true` or `false`, a value you can compute |
| `Prop` | a proposition: a statement that may be true or false, possibly about infinitely many things; you prove it, you do not compute it |
| `∀ s, ...` | for every `s` |
| `∃ s, ...` | there exists an `s` |
| `a ∣ b` | `a` divides `b` (typed `\mid`, not the bar `\|`) |
| `a % b`, `a / b` | remainder, and division rounded down; in Lean `a / 0 = 0` and `a % 0 = a` |
| `↔` | if and only if |
| `¬`, `∧`, `∨` | not, and, or (for propositions) |
| `!=`, `&&` | not-equal and and, as booleans (computable) |
| `Result α`, `ok v` | the outcome of a Rust function: `ok v` is a normal return of `v`, a failure is a panic |
| `do`, `let x ← e` | run `e`; if it fails, the whole function fails; otherwise call its value `x` |
| `m ⦃ b => P b ⦄` | Aeneas' notation: "running `m` succeeds and its result `b` satisfies `P b`" |
| `by ...` | the proof script that follows. Tactics (`simp`, `omega`, `rw`, `step*`, `scalar_tac`) tell Lean how to build the proof. You do not need to read them: if `lake build` succeeds, Lean has checked the proof |
| `sorry` | a hole in a proof. Lean accepts it with a warning; a build with a `sorry` warning proves nothing |
| `/-- ... -/`, `--` | doc comment, and line comment |
| `namespace`, `open` | name prefixes, as Rust modules and `use` |
| `#print axioms t` | lists the basic assumptions the proof of `t` depends on |

The last line is a useful sanity check. Put this in a scratch file and run it with
`cd spec && lake env lean <file>`:

```lean
import NewenSpec
#print axioms NewenSpec.Refinement.output_misses_longrange_iff
```

It should list only Lean's three standard axioms:

```text
'NewenSpec.Refinement.output_misses_longrange_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Anything else, in particular `sorryAx`, means part of the proof is missing.

## What is guaranteed and what is not

Each check proves something, relies on something, and leaves something out. Here is the honest
account.

| Check | Guarantees | Trusts | Does not cover |
|---|---|---|---|
| Lean theorem `missesModel_iff` (model ⇔ meaning) | for **all** natural numbers, the model returns `true` exactly when `Misses` holds | the Lean kernel; that `Misses` is a correct reading of the GROMOS behaviour (a human reading of `perturbed_nonbonded_set.cc`, backed by the reference test `ch4_water_fep_precalc_twinrange`) | whether the physics behind `Misses` is the physics you want |
| Aeneas refinement `output_misses_longrange_iff` (Rust ⇔ model) | for **all** `usize` inputs, the Rust function never panics and returns the model's answer | the Lean kernel; that Charon and Aeneas translate the Rust faithfully (both are nightly research tools); the Rust compiler | the callers, such as `precalclam_longrange_gaps` in `newen-run`, which add their own conditions |
| Exhaustive test (Rust ⇔ simulation, Rust ⇔ model) | the compiled binary agrees with an independent simulation and with the Lean model on 40 000 inputs | the test's simulation; the checked-in vector file being current | inputs outside `1 ≤ nsnb, ntwg ≤ 200` (the proofs cover those) |

And at the level of the whole project:

- **Not covered: everything outside `crates/newen-logic`.** The crate holds
  `schedule::output_misses_longrange`, `thermo::dof::total_dof`, and the per-bath degrees of
  freedom with the atom-to-bath mapping in `thermo::bath` (`baths_of`, `groups_starting_in`,
  `bath_half_dof`, `bath_half_dofs`, `ndfmin_spread`); those are proven. The rest of newen is
  checked the way it always was, by reference tests and unit tests.
- **Not covered: floating-point kernels**, until the differential-testing layer for them exists.
- **Not covered: physics that is self-consistent but wrong.** A specification can be internally
  consistent and still describe the wrong thing. The comparisons against gromosXX references and
  against published values stay necessary; verification adds to them, it does not replace them.

## Running it

**Rust only.** The exhaustive test needs nothing but cargo, and runs in the normal test suite:

```sh
cargo test -p newen-logic
```

**Lean.** Install elan, Lean's toolchain manager; `spec/lean-toolchain` pins the Lean version
(`leanprover/lean4:v4.35.0`, the version the Aeneas Lean library requires), and elan fetches it on
first use. Then:

```sh
cd spec && lake build
```

The first build downloads Mathlib's cache (the Aeneas library depends on Mathlib): about 7.8 GB in
`spec/.lake`, which git ignores. Later builds take seconds. The build must end without any
`sorry` warning. To rewrite the vector file the Rust test reads, from the Lean model:

```sh
cd spec && lake exe spec-vectors
```

`lake build` alone needs only elan. You need Aeneas only when you change Rust code in
`newen-logic` and must regenerate the extracted Lean.

**Aeneas and Charon.** Use the prebuilt release, not a source build. The release must be the one
whose commit is pinned in `spec/lakefile.lean`; today that is `nightly-2026.10.17-e62dd83`. It
bundles the matching Charon. Unpack it so the binaries sit in `~/.local/aeneas/aeneas`, the
default location `extract.sh` uses (override with `AENEAS_DIR`):

```sh
mkdir -p ~/.local/aeneas/aeneas
curl -L -o aeneas-linux-x86_64.tar.gz \
  https://github.com/AeneasVerif/aeneas/releases/download/nightly-2026.10.17-e62dd83/aeneas-linux-x86_64.tar.gz
tar -xzf aeneas-linux-x86_64.tar.gz -C ~/.local/aeneas/aeneas
```

Charon is a Rust compiler plugin and needs the Rust nightly named in
`~/.local/aeneas/aeneas/rust-toolchain`, which is `nightly-2026-09-17`, with three components:

```sh
rustup toolchain install nightly-2026-09-17 --profile minimal -c rustc-dev -c llvm-tools -c rust-src
```

Then, from the repository root:

```sh
spec/scripts/extract.sh
```

It runs Charon on `crates/newen-logic`, then Aeneas, and rewrites
`spec/NewenSpec/Extracted/`. It refuses to run if the Aeneas binary and the Aeneas Lean library
pinned in `spec/lakefile.lean` are different commits, because a mismatched pair produces Lean the
library cannot check. After extracting, run `cd spec && lake build` again.

When you upgrade Aeneas, change three things together: the commit in `spec/lakefile.lean`, the
Lean version in `spec/lean-toolchain` (it must equal the one in Aeneas' `backends/lean`), and the
installed binary.

## Adding a new function: the checklist

1. **Decide which kind it is.** Integers and booleans only, and small: it goes into
   `crates/newen-logic` and gets an Aeneas proof. Floating point: it stays in its kernel crate and
   gets a Lean model plus a differential test.

2. **Write the Lean first.** Before the Rust, write:
   - in `spec/NewenTheory/<Area>/<Topic>.lean`, placed by its physics, the meaning: a `Prop` over
     all inputs, transcribed from the physics or the GROMOS manual and containing no algorithm;
   - in `spec/NewenSpec/<Module>.lean`, the model, a Lean function shaped like the Rust you
     intend to write;
   - a theorem `model ↔ meaning` (or, for a float model, the properties the model must have).

   This is verification-guided development, model first: writing the meaning before the code is
   where the specification catches a misreading, before it is ported.

3. **Write the Rust.** In `crates/newen-logic`, in the subset Aeneas translates: safe Rust,
   integers and booleans, no floats, simple loops, no dependencies. End the doc comment with
   `/// spec: <module>.<name>`.

4. **Register it** with a `[[function]]` entry in `spec/SPEC.toml`. The pilot's entry:

   ```toml
   [[function]]
   key = "schedule.output_misses_longrange"
   rust = "crates/newen-logic/src/schedule.rs::output_misses_longrange"
   used_by = ["crates/newen-run/src/prepare.rs::precalclam_longrange_gaps"]
   origin = "D34: gromosXX perturbed_nonbonded_set.cc:86-96 (PRECALCLAM rows without long-range terms)"
   # Layer 1, model = meaning (core Lean, no Mathlib); the meaning lives in NewenTheory
   theory = "NewenTheory.Dynamics.Schedule.Misses"
   model = "NewenSpec.Schedule.missesModel"
   theorem = "NewenSpec.Schedule.missesModel_iff"
   # Layer 2, Rust = model (Aeneas)
   extracted = "NewenSpec.Extracted.NewenLogic: newen_logic.schedule.output_misses_longrange"
   refinement = "NewenSpec.Refinement.output_misses_longrange_eq"
   # Layer 3, Rust = meaning and Rust = model on 1 <= nsnb, ntwg <= 200 (cargo test)
   test = "crates/newen-logic/tests/schedule_spec.rs"
   vectors = "spec/vectors/schedule.json"
   ```

5. **Extract and prove refinement.** Run `spec/scripts/extract.sh`, then add a theorem to
   `spec/NewenSpec/Refinement.lean`: the extracted function `= ok (model ...)`. The pilot's proof
   is the pattern: `step*`, then `scalar_tac` or `simp_all` for what remains. For a float kernel
   this step is the differential test instead: generate random inputs, evaluate the Lean model and
   the Rust kernel, compare within a stated tolerance.

6. **Add the test.** In `crates/newen-logic/tests/`, compare the Rust against an independent
   brute-force version of the meaning and, through vectors written by `spec/SpecVectors.lean`,
   against the model. Check that the test can fail: plant a bug, watch it fail, remove it.

7. **Gate.** `cd spec && lake build` with no `sorry` warning, and `cargo test -p newen-logic`.

**When a refactor breaks a proof.** You changed Rust in `newen-logic`, re-extracted, and
`lake build` now fails in `Refinement.lean`. One of two things happened:

- *The behaviour changed.* Then the exhaustive test usually fails too. Decide which is right, the
  new code or the meaning. If the code is wrong, fix it. If the meaning was wrong, fix the meaning
  and the model, re-prove `model ↔ meaning`, and say so in the commit message: that is a change
  in what newen claims.
- *Only the shape changed* (a branch reordered, a helper extracted). Then the theorem is still
  true and only the proof script needs updating; usually re-running `step*` on the new shape is
  enough.

Never commit a `sorry` to get the build green, and never edit the files in
`spec/NewenSpec/Extracted/` by hand.

## Verified core and experimental zone

Verification is a gate for **promotion**, not for entry.

newen is also a research code. New methods change daily while they are being tried, and
requiring a Lean model and proof before every experiment would stop the experiments. Cedar can
afford model first for every change because its specification, an authorization language, is
stable. Ours is not, while a method is still being developed.

The agreed design, not yet implemented, has two zones:

- **Experimental zone.** New code sits behind an `experimental` switch, and any output it
  produces is marked as experimental in the run's provenance. It needs only the cheap
  self-consistency checks: forces equal −dV/dr by finite difference, dH/dλ equals the
  finite-difference derivative in λ, a null perturbation (identical end states) gives zero, and an
  NVE run conserves energy without drift.
- **Verified core.** Code is promoted when it becomes a default, or goes into a paper or a
  tutorial. Promotion requires a Lean model, the theorem, and the test described in this chapter.

This keeps the cost of verification where its value is: on code that other people will rely on.

## Prior art

None of this is new in kind; newen combines existing pieces.

- **LeanLJ.** Ugwuanyi, Jones, Velkey & Josephson (2025), [arXiv:2505.09095](https://arxiv.org/abs/2505.09095).
  Lennard-Jones energies in Lean, with functions polymorphic over ℝ and `Float`, so the same
  definition is both proven about and executed; the results match NIST benchmark values. A
  candidate to reuse for our LJ model.
- **Bobbin et al. (2024)**, Digital Discovery 3, 264,
  [doi:10.1039/D3DD00077J](https://doi.org/10.1039/D3DD00077J). Formalizing chemical physics in
  Lean, by the same group.
- **Cedar.** Cutler et al. (2024), [arXiv:2407.01688](https://arxiv.org/abs/2407.01688). AWS's
  authorization language: Lean models and proofs, with the Rust implementation tested against
  the models by differential random testing. The proofs found 4 bugs, differential testing 21. Our
  floating-point plan follows this approach.
- **Boldo et al. (2013)**, Journal of Automated Reasoning,
  [arXiv:1112.1795](https://arxiv.org/abs/1112.1795). A complete mechanized proof of a small C
  program solving the wave equation, rounding errors included: the measure of what a full proof
  of floating-point code costs.
- **Aeneas.** Ho & Protzenko, [github.com/AeneasVerif/aeneas](https://github.com/AeneasVerif/aeneas).
  Translates safe Rust, through Charon, into pure functions in Lean, so that Rust code can be
  proven against a specification.

## Roadmap

Planned, in the order chosen:

- **Integer logic, with Aeneas:** degrees of freedom, and the mapping of atoms to temperature
  baths. Done: `thermo::dof::total_dof` for the whole system, and `thermo::bath` for the
  degrees of freedom per bath (in half degrees of freedom, so a constraint shared between two
  baths stays an integer; the NDFMIN spread as an exact rational, where gromosXX truncates the
  total, D40) and the atom-to-bath mapping that the thermostats and the temperature calculation
  read.
- **Floating-point kernels, with a Lean model and a differential test:** the Berendsen and
  Nosé–Hoover thermostats (see [Thermostats & Barostats](../algorithms/thermostats.md)) and the
  Lennard-Jones and reaction-field pair term (see
  [Nonbonded Interactions](../algorithms/nonbonded.md)). Done: the thermostats
  (`berendsen_scale`, `nose_hoover_scale`; ℝ theorems, a `Float` model, a bit-exact differential
  test). Planned: the Lennard-Jones and reaction-field pair term.
- **Continuous integration (planned):** on every pull request, the cargo tests and lint, which
  include the exhaustive test; a Lean job (`lake build`) only when `spec/` or
  `crates/newen-logic` changes, since the Mathlib cache is large; and nightly long runs of the
  differential tests with many more random inputs than a pull request can afford.
