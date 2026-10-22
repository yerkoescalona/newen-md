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

**Floating-point kernels.** Lennard-Jones and reaction-field pair terms, thermostats, the
integrator and the pressure. These compute with `f64`. Proving that a floating-point kernel
computes its real-number equation, within a rounding bound, is a research problem: Boldo et al. (2013) needed a long paper and machine-checked
proofs of both the method error and the rounding error for a single C program solving the one-dimensional
wave equation. For these kernels we use the technique AWS uses for its Cedar authorization
language: write an executable model in Lean, prove properties of the model, and test the Rust
kernel against the model on many random inputs (differential random testing).

| | Integer decision logic | Floating-point kernels |
|---|---|---|
| Examples | output schedules, degrees of freedom, bath mapping, index matching | LJ, reaction field, bonded terms, Berendsen, Nosé–Hoover, leap-frog, kinetic energy, pressure, periodic boundaries |
| Rust lives in | `crates/newen-logic` | the existing kernel crates |
| Lean has | meaning, model, theorem model ⇔ meaning | meaning, executable model, theorems about the model |
| Rust ⇔ model | proof, for every input (Charon + Aeneas + refinement theorem) | differential random testing |
| Status | working (one function, the pilot below) | working for the Berendsen and Nosé–Hoover thermostats, the Lennard-Jones and reaction-field pair term, the bonded terms, the leap-frog, kinetic-energy, temperature, pressure and Berendsen barostat kernels, and the periodic boundaries (below) |

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

`crates/newen-logic/src/schedule/mod.rs`:

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
Rust's and C's libraries may differ in the last bit. (The Berendsen barostat calls `pow`; both
sides reach the same C `pow` there, see below.) It also stops working when the two sides
sum in a different order, as a parallel reduction over pairs does. There a tolerance, justified by a rounding analysis, has to replace equality.

What is not guaranteed: how far the `Float` model is from the ℝ model (no rounding bound is
proved), and inputs outside the generated set.

### The Lennard-Jones and reaction-field pair term

The same three pieces cover the pair term of the nonbonded loops in `crates/newen-forces`:
`CRFParameters::new` (the reaction-field constants of a run), `lj_crf_pair` (one LJ + CRF pair at
squared distance `r2`; `lj_crf_interaction` is `lj_crf_pair(r.length_squared(), …)`), and
`rf_excluded_pair` and `rf_self_term` (the reaction-field terms of excluded solute pairs and of each
solute atom with itself).

**The meaning.** `spec/NewenTheory/Interactions/LennardJones.lean` defines
`ljEnergy c6 c12 r = c12/r¹² − c6/r⁶` and the scalar force the kernels compute, `ljForce`, which
multiplies `r_vec` to give the force on atom `i`. The theorem that ties them together reads:

```lean
theorem lj_force_is_minus_grad (c6 c12 : ℝ) {r : ℝ} (hr : r ≠ 0) :
    HasDerivAt (ljEnergy c6 c12) (-(ljForce c6 c12 r) * r) r
```

Read aloud: at every `r ≠ 0` the derivative of the energy is `−ljForce · r`, so `ljForce` is
`−E'(r)/r`: force = −dE/dr. `lj_zero_at_sigma` (the energy vanishes at `r⁶ = c12/c6`) and
`lj_minimum` (at `r⁶ = 2 c12/c6` the derivative vanishes and the energy is `−c6²/(4 c12)`) complete
it. `ReactionField.lean` does the same for the reaction field of Tironi et al. (1995):
`crfConstant ε ε_rf κ R`, `crfEnergy q crf R r = q (1/r − crf r²/(2R³) − (1 − crf/2)/R)` with
`crf_force_is_minus_grad`, `crf_energy_zero_at_cutoff` (the energy is zero at `r = R`),
`crf_eps_rf_eq_eps` (`κ = 0` and `ε_rf = ε`: no reaction field, a Coulomb energy shifted to zero at
`R`), `crf_conducting_limit` (`C_RF → −1` as `ε_rf → ∞`, which is why GROMOS writes `ε_rf = 0` for a
conducting boundary and newen sets `C_RF = −1` there), and `rf_excluded_is_crf_minus_coulomb` (an
excluded pair keeps the pair term without its `q/r` and `q/r³`), with `rf_self_is_half_excluded`
(the self term is half the excluded term of the atom with itself at `r = 0`).

**One model.** `spec/NewenSpec/Nonbonded.lean` writes the four kernels once over any number type,
operation by operation like the Rust, with the guard `ZERO_DISTANCE2 = 1e-10` and the `0.5` of the
self term as parameters given at `Float` by their bits. At ℝ, `ljCrfPairModel_real_eq` proves the
pair model returns `(ljForce + crfForce, ljEnergy, crfEnergy)` for `r > 0` outside the guard (the
model takes `sqrt(1/r²)`, which is `1/r` there), `crfParametersModel_real_eq` that its `C_RF` is
`crfConstant` outside the two sentinels (`R ≤ 0`, `ε_rf = 0`), and `rfExcludedPairModel_real_eq`,
`rfSelfTermModel_real_eq` the same for the excluded-pair and self terms.

**The bit-exact test.** `lake exe spec-vectors` writes `spec/vectors/nonbonded.json`: 600 rows of
reaction-field constants, 6 000 pair rows, 1 500 excluded-pair rows and 300 self-term rows, over
cutoffs `0, 0.8, 0.9, 1.0, 1.4` nm and random ones, `ε_rf ∈ {0, 1, 54, 61, 62, 78.5}` and random,
`κ ∈ {0, 0.5, 1}` and random, `r²` log-uniform in `[10⁻⁴, 4]` nm² with rows below the guard, at it
and at exactly `R²`, and GROMOS-sized `C6`, `C12` and charge products, zeros included. Every pair row
carries the four inputs of `CRFParameters::new`, so the Rust rebuilds the parameters and the pair
rows test both. `crates/newen-forces/tests/nonbonded_spec.rs` requires the same 64 bits on every
row, and `ZERO_DISTANCE2.to_bits()` to equal the Lean constant.

The production loops evaluate four pairs at a time in `lj_crf_interaction_x4`, a copy of the
kernel on `wide::f64x4`. A unit test in `innerloops.rs` runs it on the same pair rows, four at a
time, and requires every lane to return the scalar `lj_crf_pair`'s bits. It passes with
`target-cpu=native` and with `target-cpu=generic`, in debug and release builds; moving one
parenthesis in the SIMD force line makes 4 940 lanes differ.

The model follows newen's Rust, not gromosXX, and the order differs in two places, both within
the reference tolerances:

- gromosXX computes the LJ force as `(c12/r⁶ + c12/r⁶ − c6) · 6 · r⁻⁶ · r⁻²`
  (`nonbonded_term.cc:241`), newen as `(12 · c12 · r⁻⁶ − 6 · c6) · r⁻⁶ · r⁻²`, the textbook
  derivative. Written gromosXX's way, 1 267 of the 6 000 pair rows differ, by 1 or 2 ulp.
- gromosXX forms `(κR)²` as `κ · R · κ · R`, newen as `(κ · R) · (κ · R)`.

The self term agrees bit for bit outside underflow: gromosXX books `0.5 ·` the excluded-pair
energy at `r = 0` with charge product `q_i²`, newen `−0.5 · q_i · q_i · (4πε₀)⁻¹ · crf_cut`, and a
factor 0.5 is exact.

```sh
cargo test -p newen-forces --test nonbonded_spec    # Rust = Float model, bit for bit
cargo test -p newen-forces --lib x4_lanes           # SIMD lanes = scalar kernel, bit for bit
```

### The bonded terms

The covalent terms of `crates/newen-forces/src/bonded/` have one pure kernel per functional form,
called by the loops that walk the topology (which keep the bookkeeping: energy groups, virial,
skipped constrained bonds):

| Kernel | Energy | gromosXX |
|---|---|---|
| `quartic_bond(r, k, b0)` | `K/4 · (b² − b0²)²` | `quartic_bond_interaction.cc` |
| `harmonic_bond(r, k, b0)` | `K/2 · (b − b0)²` | `harmonic_bond_interaction.cc` |
| `cosine_angle(r_ij, r_kj, k, cos0)` | `K/2 · (cos θ − cos θ0)²` | `angle_interaction.cc` |
| `harmonic_angle(r_ij, r_kj, k, theta0)` | `K/2 · (θ − θ0)²` | `harm_angle_interaction.cc` |
| `dihedral_chebyshev(r_ij, r_kj, r_kl, k, cos_delta, m)` | `K (1 + cos δ · cos mφ)`, NTBDN = 1 | `dihedral_interaction.cc` |
| `dihedral_new(r_ij, r_kj, r_kl, k, delta, m)` | `K (1 + cos(mφ − δ))`, NTBDN = 0 | `dihedral_new_interaction.cc` |
| `improper_dihedral(r_ij, r_kj, r_kl, k, zeta0)` | `K/2 · (ζ − ζ0)²` | `improper_dihedral_interaction.cc` |

Each returns the energy and the forces on its atoms; its doc comment gives the formula in the
order the code evaluates it.

**The meaning.** `spec/NewenTheory/Interactions/Bonded.lean` writes vectors of ℝ³ out by components
(`V3`, with `dot`, `cross`, `len`) and defines the energies and the force expressions the kernels
evaluate. The theorems:

- `quarticBond_force_is_minus_grad`, `harmonicBond_force_is_minus_grad`: along every direction
  `e`, the derivative of `V(|r + t e|)` at `t = 0` is `−f_i · e`, so the bond force is `−∇V`;
  `quartic_harmonic_match`: the curvature of the quartic bond at `b0` is `2 K b0²`, the GROMOS
  relation `K_harm = 2 K_quartic b0²` (Vol. 2).
- `cosMultiple_eq_cos`: for `1 ≤ m ≤ 6` the code's polynomial in `cos φ` is `cos(mφ)` (by the
  recursion `cos((n+2)φ) = 2 cos φ cos((n+1)φ) − cos(nφ)`); `dcosMultiple_is_deriv`: its second
  polynomial is the derivative of the first; `dihedral_new_eq_dihedral`: with `δ ∈ {0, π}` the
  NTBDN = 0 and NTBDN = 1 energies are equal, which is why the two forms agree on every GROMOS
  topology.
- `cosineAngle_forces_sum_zero`, `cosineAngle_torque_zero`, `harmonicAngle_conservation`,
  `chebyshev_conservation`, `newDihedral_conservation`, `improper_conservation`: the forces of an
  angle or a dihedral sum to zero and exert no torque about any origin. For the two forms whose
  forces lie along the plane normals (NTBDN = 0 and the improper) this needs Lagrange's identity
  `|r_ij|²|r_kj|² − (r_ij · r_kj)² = |r_ij × r_kj|²` and holds away from a linear bond angle, where
  the kernels drop a force on purpose.
- `cosineAngle_force_is_minus_grad`: moving atom `i`, `j` or `k` along any `e`, the
  cosine-harmonic energy changes at the rate `−f · e` of that atom's force. The other angle and the
  dihedrals have the conservation theorems but no gradient theorem; the finite-difference test below
  covers them.

**One model.** `spec/NewenSpec/Bonded.lean` writes the seven kernels once over any number type, on
a `Vec3 α` with glam's `dot` (`x·x' + y·y' + z·z'`, left to right), `cross` and `length`. The bonded
kernels also call `acos`, `cos` and `sin`, so `NewenSpec/Float.lean` gains the classes `Acos`, `Cos`
and `Sin` next to `Sqrt` (`Float.acos`, … at `Float`; `Real.arccos`, … at ℝ). The guards (all
`1e-10`: `ZERO_DISTANCE`, `ZERO_DISTANCE2`, `ZERO_SINE`, `NEAR_LINEAR`), `0.25`, `0.5` and `π` are
parameters given at `Float` by their bits; the two `while` loops that move ζ into `[ζ0 − π, ζ0 + π]`
run with fuel (one turn suffices for `|ζ0| ≤ π`). At ℝ, one `*Model_real_eq` per kernel proves
that, outside the guards, the model returns the theory's energy at the theory's geometry and the
theory's force expressions; the clamps of `cos θ`, `cos φ`, `cos ζ` to `[−1, 1]` never act there
(`cos_between_mem`, Cauchy–Schwarz). The improper's theorem assumes ζ already in
`[ζ0 − π, ζ0 + π]`, always true for `ζ0 = 0`.

**The tests.** `lake exe spec-vectors` writes `spec/vectors/bonded.json` (under 2 MB): 1 482 bond
rows, 1 000 angle rows and 710 dihedral rows, each row evaluated by every kernel of its kind
(quartic and harmonic bond; cosine and harmonic angle; Chebyshev, NTBDN = 0 and improper
dihedral). Random rows have components in `[−0.5, 0.5]` nm or bond lengths within 20 % of `b0`,
GROMOS-sized parameters (`K_quartic ∈ [10⁵, 2·10⁷]` kJ mol⁻¹ nm⁻⁴, `b0 ∈ [0.09, 0.2]` nm,
`θ0 ∈ [90°, 125°]`, dihedral `K ∈ [0.1, 50]` kJ mol⁻¹, `m ∈ 0..6` and the unsupported 7,
`δ ∈ {0, π}` or random, `ζ0 ∈ {0, 35.26°}` or random); the edge rows hit every guard: zero and
guard lengths, collinear angles with `cos θ = ±1` exactly (`|sin θ| < 10⁻¹⁰`, with `θ0 = π` and
without), planar dihedrals, a bond angle within 10⁻⁵ rad of 180° inside a dihedral, `|r_kj|²`
below the guard, ζ at ±π. `crates/newen-forces/tests/bonded_spec.rs` requires the same 64 bits as
the Lean model on every row, and the bits of the four guards and of `std::f64::consts::PI`. It also
checks the physics independently of the bits: on 200 well-conditioned rows per kernel, the
central difference of the kernel's energy (`h = 10⁻⁶` nm) along every coordinate of every atom is
minus the returned force, to `10⁻⁵` of the row's largest force component plus `10⁻⁶`.

`acos`, `cos` and `sin` are not correctly rounded by IEEE 754. The bits agree because Lean's
`Float` functions and Rust's `f64` functions call the same C library (glibc's `libm`); on a
machine with another `libm` the Lean vectors would have to be regenerated there. The test passes
in debug, release and `target-cpu=generic` builds.

Two changes to the Rust came with the kernels; neither moves a gromosXX reference:

- **`powi` replaced by gromosXX's products.** The Chebyshev dihedral wrote `cos³φ … cos⁶φ` with
  `powi`, whose rounding is the compiler's business; gromosXX forms `c² = c·c`, `c³ = c²·c`,
  `c⁴ = c³·c` and writes `16 · c⁴ · c` for `c⁵`. The kernel now does the same. On the 637
  dihedral rows with `m ≤ 6` outside the guards, the `powi` form gives other bits for 138 (all with
  `m ≥ 3`), and at 71 % of 200 001 evenly spaced `cos φ ∈ [−1, 1]` for at least one `m`. No gromosXX reference used
  NTBDN = 1, so `butane_vacuum_ntbdn1` (`butane_vacuum` with `COVALENTFORM 0 0 1`) was added first
  and passes before and after.
- **`cos θ0` once per angle type.** The cosine-angle loop took `theta0.cos()` for every angle;
  gromosXX stores `cos θ0` when it reads the topology. The loop now computes it once per type
  before the walk and passes it to `cosine_angle`; the bits are the same.

The Rust follows newen where it differs from gromosXX in the order of operations, all within the
reference tolerances (each listed in the `origin` of its `SPEC.toml` entry): the cosine-angle
force is `(r_kj · (1/d_kj) − r_ij · (cos θ/d_ij)) · (df/d_ij)` where gromosXX writes
`df/d_ij · (r_kj/d_kj − r_ij/d_ij · cos θ)`; the Chebyshev force groups `k_i` with the scalars; the
NTBDN = 0 force divides the scalar `k_i · d_kj` by `|r_mj|²` before scaling `r_mj`. One oddity is
shared with gromosXX and stated in the theory: for `m = 0` the Chebyshev form takes `cos mφ = 0`,
not `cos 0 = 1`, so its energy is `K` where the NTBDN = 0 form gives `K (1 + cos δ)`; the force is
zero in both.

```sh
cargo test -p newen-forces --test bonded_spec    # Rust = Float model, bit for bit; forces = −∇E by finite differences
cargo test -p newen-md --test test_gromosXX_references butane_vacuum_ntbdn1
```

### The integrator, temperature and pressure

The same three pieces cover the per-step dynamics of `crates/newen-integrators`. Each production
algorithm calls a pure kernel, extracted without changing one operation:

| Kernel | Called by | Computes |
|---|---|---|
| `leap_frog_velocity(v, f, m⁻¹, dt)` | `LeapFrogVelocity::apply` | `v + f · m⁻¹ · dt` per component |
| `leap_frog_position(x, v, dt)` | `LeapFrogPosition::apply` | `x + v · dt` per component |
| `kinetic_energy(m, v_new, v_old)` | `TemperatureCalculation::apply` | `Σ ½ m (‖v_new‖² + ‖v_old‖²)/2` and `Σ ½ m ‖v_new‖²`, in atom order |
| `group_kinetic_energy(m, v_new, v_old)` | `TemperatureCalculation::apply`, per temperature group | the centre-of-mass energy of the averaged velocities, the half-step total, and the same two for the new velocities |
| `kinetic_temperature(K, N_df, k_B)` | `md` (bath temperatures of the energy trajectory, the energy check) | `2K / k_B / N_df`, 0 for `N_df = 0` |
| `pressure_tensor(K, W, V)` | `PressureCalculation::apply` | `(K + W · ½) · (2/V)` element by element, `W = Σ r⊗f` |
| `scalar_pressure(P)` | `BerendsenBarostat::apply` | `(P_xx + P_yy + P_zz)/3` |
| `berendsen_barostat_scale(P, P₀, β, τ, dt)` | `BerendsenBarostat::apply` | `powf(1 − β·dt/τ·(P₀ − P), 1/3)` |

**The meaning.** `spec/NewenTheory/Dynamics/Integration.lean` states the leap-frog step for an
acceleration field in any real vector space, `leapFrog a dt (x, v) = (x + dt·v', v')` with
`v' = v + dt·a(x)`, and proves the three properties that make it the integrator of molecular
dynamics:

```lean
theorem leapFrog_reversible (a : E → E) (dt : ℝ) (s : E × E) :
    leapFrog a dt (reverse a dt (leapFrog a dt s)) = reverse a dt s

theorem leapFrog_iterate_exact (x₀ u a₀ : E) (dt : ℝ) (n : ℕ) :
    (leapFrog (fun _ => a₀) dt)^[n] (x₀, u - (dt / 2) • a₀) =
      (x₀ + (n * dt) • u + ((n * dt) ^ 2 / 2) • a₀, u + (n * dt - dt / 2) • a₀)

theorem leapFrog_jacobian (a : ℝ → ℝ) {a' x v dt : ℝ} (ha : HasDerivAt a a' x) :
    HasDerivAt (fun y => (leapFrog a dt (y, v)).1) (1 + dt ^ 2 * a') x ∧ … ∧
      (1 + dt ^ 2 * a') * 1 - dt * (dt * a') = 1
```

Read aloud: a step, a velocity reversal and a step bring the system back to where it started,
with the velocity reversed (`reverse` flips the velocity of the *other* half step,
`−(v + a(x) dt)`, because the state carries the half-step velocity before `x`; step-then-reverse
is an involution, `reverse_comp_leapFrog_involutive`); under a constant acceleration `n` steps
land exactly on the parabola; and in one dimension the Jacobian of the step has determinant 1, so
phase-space area is conserved.

GROMOS takes the kinetic energy at time `t` as the average of the two half-step energies.
`halfStep_vs_verlet` says what that is, exactly: after the kick `v₊ = v₋ + a dt`,

```lean
halfStepKinetic m v₋ (v₋ + dt • a) = 1/2 * m * ‖v₋ + (dt/2) • a‖ ^ 2 + m * dt ^ 2 / 8 * ‖a‖ ^ 2
```

the energy of the full-step velocity that velocity Verlet would use, plus a second-order term that
is never negative. For a temperature group, `group_split` decomposes the half-step energy into the
centre-of-mass energy of the averaged velocities (what GROMOS books to the COM bath), the averaged
velocities relative to their centre of mass, and the half kicks (König's theorem, `koenig`, plus
the parallelogram law), so `group_internal_nonneg`: the part booked to the internal bath is never
negative. `temperatureOf_iff` says `T = 2K/(N_df k_B)` is equipartition `K = ½ N_df k_B T`
inverted.

`spec/NewenTheory/Dynamics/Pressure.lean` states the virial theorem in GROMOS's tensor form,
`P = 2/V · (K − Ξ)` with `Ξ = −½ Σ r⊗f`; `scalarPressure_virial` gives the scalar form
`P = 2/(3V) · (tr K − tr Ξ)` and `ideal_gas` the ideal gas law without forces. For the barostat
factor `μ = (1 − β dt/τ (P₀ − P))^{1/3}`: `barostat_fixed` (`P = P₀` gives `μ = 1`),
`barostat_grows_iff` and `barostat_shrinks_iff` (with a positive coupling the box grows exactly
when `P > P₀` and shrinks exactly when `P < P₀`), `barostat_volume` (the volume is scaled by the
radicand), and `barostat_relax`: if the pressure responds to the volume as the compressibility
says, `P' = P − (V' − V)/(βV)`, one step moves the pressure the fraction `dt/τ` of the way to `P₀`,
the barostat's analogue of `berendsen_relax`.

**One model.** `spec/NewenSpec/Dynamics.lean` writes the eight kernels once over any number type,
a `glam::DVec3` as a structure of three numbers with `length_squared` as glam writes it
(`x·x + y·y + z·z`), and every sum over atoms as the Rust loop `acc += term` from zero (a left
fold). The `0.5` literal and the `1.0 / 3.0` of the barostat are parameters, given at `Float` by
their bits. At ℝ, `leapFrogModel_real_eq` proves the two leap-frog kernels are the theory's step
with acceleration `m⁻¹ f`, `kineticEnergyModel_real_eq` and `groupKineticModel_real_eq` that the
kinetic energies are the theory's sums (velocities as vectors of `EuclideanSpace ℝ (Fin 3)`, the
centre of mass `comVelocity`), `kineticTemperatureModel_real_eq` that the temperature is
`temperatureOf`, and `pressureTensorModel_real_eq`, `scalarPressureModel_real_eq`,
`berendsenBarostatModel_real_eq` the same for the pressure and the barostat.

**The bit-exact test.** `lake exe spec-vectors` writes `spec/vectors/dynamics.json`: 1 200
leap-frog rows (velocities in `[−3, 3]` nm/ps, forces to `10⁶` kJ mol⁻¹ nm⁻¹, GROMOS masses,
random ones and massless atoms with `m⁻¹ = 0`), 500 groups of 1–8 atoms for both kinetic-energy
kernels (one in 40 massless, a centre of mass of `0/0`), 800 temperature rows (`N_df = 0` and a
negative count included), 600 pressure rows (`V = 0` included) and 1 500 barostat rows, a quarter
of them with a coupling term of order 1, because at GROMOS values `β dt/τ (P₀ − P)` is far below
the last bit of 1 and the result hides how the term was rounded.
`crates/newen-integrators/tests/dynamics_spec.rs` requires the same 64 bits on every row, the
position kernel fed the model's new velocity so that each kernel is tested on its own. Planting a
reassociation fails it: `v + f · (m⁻¹ · dt)` makes 331 of the 1 200 leap-frog rows differ,
`2K/(N_df k_B)` 216 of the 800 temperature rows, `com_v *= 1/M` instead of `com_v /= M` 137 of
the 500 groups, summing the atoms in reverse 208 of them, `K·(2/V) + W·½·(2/V)` 442 of the 600
pressure rows, `β · (dt/τ)` 17 of the 1 500 barostat rows. Rewrites that only scale by a power of
two (`½ m ((a + b)/2)` for `(½ m (a + b))/2`, `1/V · 2` for `2/V`) are exact and pass, as they
should.

The barostat is the one kernel that calls a function IEEE 754 does not round correctly,
`powf(x, 1/3)`. Rust's `f64::powf` and Lean's `Float.pow` both call the C library's `pow` (the
Lean executable links the system `libm` dynamically, as Rust's standard library does), so the
comparison is kept strict and holds where both run against the same `libm`. It does catch a
change of function: `cbrt` in place of `powf(…, 1.0 / 3.0)` makes 52 of the 1 500 rows differ (by
1 ulp, and for a negative radicand `cbrt` returns a negative number where `pow` returns NaN).

The model follows newen, not gromosXX, and differs from it in these places:

- gromosXX's velocity step is `v + f · dt / m` (`leap_frog.cc`), newen's `v + f · m⁻¹ · dt` with
  `m⁻¹` precomputed: differences in the last bits. A massless atom has `m⁻¹ = 0` in newen and
  does not move; gromosXX divides by zero.
- gromosXX multiplies the pressure tensor by `tot_cg_factor` (1 without coarse graining, an exact
  factor) and divides by any volume; newen skips the update for `V < 10⁻³⁰`.
- The temperatures written out are `2K/k_B/N_df` (gromosXX's `out_configuration.cc` order), the
  thermostats' `2K/(N_df k_B)` (gromosXX's thermostat order); the two can differ in the last bit.

```sh
cargo test -p newen-integrators --test dynamics_spec    # Rust = Float model, bit for bit
```

### Periodic boundaries

Every pair interaction starts from the separation of two atoms in a periodic box, so the
nearest-image function is called more often than any other kernel in newen. The three pieces
cover `crates/newen-core/src/math.rs`: `vacuum_nearest_image`, `Rectangular::new` with
`rectangular_nearest_image` and `rectangular_put_into_box`, and `Triclinic::new` with
`triclinic_nearest_image` and `triclinic_put_into_box`. The `BoundaryCondition` methods of
`Vacuum`, `Rectangular` and `Triclinic`, and so `Periodicity`, call these functions. A truncated
octahedron is read into a triclinic box (`truncoct_triclinic_box`), so it runs the triclinic
kernels.

**The meaning.** `spec/NewenTheory/Space/Periodicity.lean` describes the lattice of a box.
Two points are images of each other when they differ by an integer combination of the box
vectors `a, b, c` (`Box.Lattice`). A **cell** is a region that holds exactly one image of every
point. In GROMOS's lower-triangular box (`a` along `x`, `b` in the `xy` plane) the fractional
coordinates of a point follow by back substitution (`Box.frac`). Two cells are used:

- the centred cell, every fractional coordinate in `(−1/2, 1/2]`. The nearest image
  (`Box.image`) lives there.
- the positive cell, every coordinate in `[0, 1)`. newen's put-into-box (`Box.wrap`) puts a
  point there.

For both, the theory proves that the result differs from the input by a lattice vector
(`image_lattice`), lies in its cell (`image_inCell`), is the only point of the cell that does
(`inCell_unique`, `image_eq_of`), and so is unchanged when applied again (`image_idem`). The
same four hold for `wrap`.

The rectangular kernel is stated component by component:
`rectImage L x = x − L · roundEven(x / L)`, where `roundEven` rounds halves to the even integer as
IEEE 754 does (Mathlib's `round` sends halves up). It is a lattice translate, lies in
`[−L/2, L/2]`, is idempotent, and is the **shortest** image:

```lean
theorem rectImage_shortest {L : ℝ} (hL : 0 < L) (x : ℝ) (n : ℤ) :
    |rectImage L x| ≤ |x - n * L|
```

Read aloud: no other image of `x` is closer to the origin. `rectImageVec_shortest` lifts this to
the squared length of the vector. The triclinic reduction promises less, and the theory says so:

```lean
theorem image_not_shortest :
    ∃ (B : Box) (r w : Vec3 ℝ), B.LowerTriangular ∧ B.Lattice w ∧
      normSq (B.image r - w) < normSq (B.image r)
```

Read aloud: there is a box and a point whose image in the centred cell can be shortened by a
lattice vector. In the box `a = (2, 0, 0)`, `b = (1, 1, 0)`, `c = (0, 0, 1)` the point
`(6/5, 2/5, 0)` is its own image, with squared length `8/5`, while `r − a` has `4/5`. What
gromosXX's triclinic `nearest_image` guarantees is the centred cell, not the shortest vector.

**One model.** `spec/NewenSpec/Periodicity.lean` writes each kernel over any number type, with
glam's vector and matrix operations in glam's order: `DMat3 * DVec3` is
`(x_axis·v.x + y_axis·v.y) + z_axis·v.z`, and `DMat3::inverse` is the cross products of the
columns over the determinant, transposed. Two operations have no common Lean class, so
`NewenSpec/Float.lean` adds them:

- `Floor`: `Float.floor` (C `floor`, which is what Rust's `f64::floor` computes) and `⌊x⌋` over ℝ.
- `Rint`: Rust's `f64::round_ties_even` is C `rint`, halves to even. Lean's `Float.round` is C
  `round`, halves away from zero, which is Rust's `f64::round` and a different function. So
  `Float.lean` binds C `rint` itself (`@[extern "rint"] opaque rintFloat`) and models it over ℝ
  as `roundEven`.

Both roundings are exact, so any correct implementation gives the same bits, the sign of zero
included.

The triclinic kernel is six `while` loops: subtract `c` while the fractional coordinate along
`c` is at least `1/2`, add it while it is at most `−1/2`, then the same for `b` and for `a`. A
Lean function must terminate, so each loop is `whileFuel`, which stops after `fuel` turns.
`triclinicNearestImageModel_real_eq` proves that over ℝ, with `0.5` as `1/2` and a fuel larger
than every fractional coordinate plus one, the six loops compute `Box.image`. The proof tracks
the three fractional coordinates through the loops: each loop moves only its own coordinate, by
one per turn. The other model theorems are `rectangularNearestImageModel_real_eq` (the model is
`rectImage` in each component, because `x · (1/L) = x / L` and `rint` is `roundEven`),
`rectangularPutIntoBoxModel_real_eq`, `mulVec_inverseModel_real` (glam's inverse applied to a
point gives its fractional coordinates) and `triclinicPutIntoBoxModel_real_eq` (the model is
`Box.wrap`).

**The bit-exact test.** `lake exe spec-vectors` writes `spec/vectors/periodicity.json`: 4 420
rows over cubic, anisotropic, truncated-octahedron, reduced triclinic, strongly skewed and
monoclinic boxes. The separations and positions include:

- values inside the box, a few boxes out and fifty boxes out, and negative values;
- exact odd multiples of half an edge (the ties of `rint`), and one ulp on either side of half
  an edge;
- whole multiples of the edge, and zeros of both signs;
- separations exactly on the faces of the triclinic cell (half a box vector, and sums of
  halves);
- NaN and infinite components. The triclinic nearest image gets NaN only, because an infinite
  separation never leaves the Rust `while` loops.

Each row carries the box, not its reciprocals or inverse. So the Rust test builds
`Rectangular::new` or `Triclinic::new` from the box and calls the kernel through `Periodicity`,
the way the force loops do. The generator runs the triclinic model with 100 000 turns per loop,
and refuses a row whose bits change at half that. A separate section checks that the truncated
octahedra among the boxes are the ones `truncoct_triclinic_box` builds. Replacing
`round_ties_even` by `round` in the rectangular kernel, or `>=` by `>` in the first triclinic
loop, makes the test fail.

The model follows newen's Rust. Where it differs from gromosXX:

- **Rectangular nearest image.** gromosXX loops: while `x > L/2` subtract `L`, while `x < −L/2`
  add `L` (`boundary_implementation.cc:224-261`). newen subtracts `L · rint(x · (1/L))` in one
  step. Over ℝ both end in `[−L/2, L/2]`, but at `x = ±3L/2, ±7L/2, …` they pick images a box
  length apart (gromosXX keeps `L/2`, newen goes to `−L/2`; both are equally long). In
  floating point, the reciprocal and the single multiply round differently from gromosXX's
  repeated subtraction.
- **Put-into-box.** gromosXX's `Periodicity::put_into_box` is the nearest image of the origin,
  the centred cell (`periodicity.cc:36-41`). newen's `put_into_box` uses the positive cell,
  `[0, L)` and fractional `[0, 1)`. `gather.rs` reproduces gromosXX's choice with
  `nearest_image(v, 0)` where it matters.
- **Triclinic nearest image.** It is gromosXX's, loop for loop and operation for operation.

```sh
cargo test -p newen-core --test periodicity_spec    # Rust = Float model, bit for bit
```

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
  `schedule::output_misses_longrange`, the run schedule (`schedule::output::due`, `frames_in`,
  `minimisation_stop`, `final_step`; `schedule::com::com_removal`), `thermo::dof::total_dof`,
  the per-bath degrees of freedom with the atom-to-bath mapping in `thermo::bath` (`baths_of`,
  `groups_starting_in`, `bath_half_dof`, `bath_half_dofs`, `ndfmin_spread`, `empty_baths`), and
  the topology bookkeeping in `topology` (`pairs::contains_sorted`, `pairs::pair_kind`,
  `groups::group_of`, `groups::check_last_atoms`, `groups::energy_cell`); those are proven. Their
  callers in newen-core, newen-run and newen-forces (building the partner lists sorted and
  symmetric, which the pair-kind theorems assume) are checked by unit tests. The pairlist rebuild
  counter and the refusal matrix are not there yet. The rest of newen is checked the way it always
  was, by reference tests and unit tests.
- **Floating-point kernels: covered on generated inputs only.** The thermostat steps, the
  Lennard-Jones and reaction-field pair term, the bonded terms, the leap-frog, kinetic-energy,
  temperature, pressure and barostat kernels, and the periodic boundaries equal their Lean models
  bit for bit on the rows `lake exe spec-vectors` writes; no rounding bound is proved, and inputs
  outside those rows are not checked. The bonded kernels' agreement also assumes the same C
  `acos`, `cos` and `sin` on both sides, the barostat's the same C `pow`.
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
   rust = "crates/newen-logic/src/schedule/mod.rs::output_misses_longrange"
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
- **Run schedule: output frames, minimisation stop, COM removal, bath decoupling, done;
  pairlist counter and the refusal matrix planned.** `schedule::output` (which steps write a
  frame, how many a run writes, the frames a converged minimisation writes on the step it stops
  at, proven equal to gromosXX's top-of-loop and final writes; the final configuration's step),
  `schedule::com::com_removal` (NTICOM at step 0, NSCM after it, no periodic rotation removal in
  a periodic box) and `thermo::bath::empty_baths` (the baths gromosXX decouples, decided before
  the NDFMIN spread).
- **Topology bookkeeping: pair kinds and group lookups, done; index conversions and
  perturbed-atom matching planned.** Done: `topology::pairs` (whether a pair is excluded, 1-4,
  owned by a `PERTATOMPAIR`/EDS atom, or on the pairlist) and `topology::groups` (the energy
  group of an atom, the FORCE `NRE` rule, the energy-matrix cell).
- **Floating-point kernels, with a Lean model and a differential test:** the Berendsen and
  Nosé–Hoover thermostats (see [Thermostats & Barostats](../algorithms/thermostats.md)) and the
  Lennard-Jones and reaction-field pair term (see
  [Nonbonded Interactions](../algorithms/nonbonded.md)). Done: the thermostats
  (`berendsen_scale`, `nose_hoover_scale`), the pair term (`CRFParameters::new`,
  `lj_crf_pair`, `rf_excluded_pair`, `rf_self_term`, with the SIMD lanes checked against the scalar
  kernel), the bonded terms (`quartic_bond`, `harmonic_bond`, `cosine_angle`,
  `harmonic_angle`, `dihedral_chebyshev`, `dihedral_new`, `improper_dihedral`, with a
  finite-difference check of every force), and the dynamics kernels (`leap_frog_velocity`,
  `leap_frog_position`, `kinetic_energy`, `group_kinetic_energy`, `kinetic_temperature`,
  `pressure_tensor`, `scalar_pressure`, `berendsen_barostat_scale`); each with ℝ theorems, a
  `Float` model and a bit-exact differential test. Not yet covered in `newen-integrators`: the
  constraint algorithms (SHAKE, SETTLE, LINCS), centre-of-mass removal, the anisotropic barostat
  modes gromosXX has, and the older `integrator.rs` and `barostats.rs` implementations, which the
  MD sequence does not use.
- **Periodic boundaries, the same way** (see [Periodic boundaries](#periodic-boundaries)). Done:
  the nearest image and the put-into-box of every box type (vacuum, rectangular, triclinic; the
  truncated octahedron is a triclinic box) and the two box constructors, with the lattice, cell,
  uniqueness and idempotence theorems, the rectangular image proven shortest, and the triclinic
  one proven not always shortest. Not covered: image code outside `newen-core`, such as the
  `[−L/2, L/2]` wrapping in `remove_com_motion.rs`, which uses its own `f64::round` and `x / L`.
- **Continuous integration (planned):** on every pull request, the cargo tests and lint, which
  include the exhaustive test; a Lean job (`lake build`) only when `spec/` or
  `crates/newen-logic` changes, since the Mathlib cache is large; and nightly long runs of the
  differential tests with many more random inputs than a pull request can afford.
