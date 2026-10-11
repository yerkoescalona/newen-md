# Porting guidelines — faithful to GROMOS, without its bugs

How to translate gromosXX/gromos++ into gromos-rs so the result is *both* faithful (matches GROMOS
where it should) and *correct* (doesn't inherit 30 years of accreted C++ defects). These are
guidelines, not tasks — the live work is on the project board, and the A/B/C classification below is
the discipline itself, not a summary of it.

## The bit-for-bit invariant

The gromosXX reference suite (`crates/gromos-md/tests/gromosXX_references/`) is the oracle. For any
**wired, tested** path:

- Reproducing GROMOS output to full precision is the definition of correct. Tolerances are tight
  (force ~1e-6, energy ~1e-8) by design.
- A change that alters that output is a **regression until proven otherwise** — never a silent
  improvement. If gromosXX is genuinely wrong and you want the corrected behavior, that is a
  *conscious, documented decision* (and a candidate for a `--gromos-compat` vs `--corrected` split),
  not an undiscussed divergence.

## The oracle only covers wired paths

This is the core risk. The reference suite catches divergence **only** for features that are wired
and tested. Every *unwired* path (e.g. triclinic boundaries, FEP, grid pairlists, Nosé-Hoover,
EDS/GaMD, lattice-sum/PME) is un-audited surface where two failure modes hide:

- **Inherited C++ bugs** — you faithfully port a defect.
- **Silent divergence** — your "cleaner" reimplementation disagrees with GROMOS and nothing notices
  until the feature is wired.

(A live example already in the tree: the triclinic nearest-image — see 1.4.)

So the standing rule: **write the gromosXX reference test *before* wiring an unwired path, not
after.** The test is what converts un-audited surface into oracle-covered surface.

## Before porting a subsystem

1. **Grep the C++ for the authors' own flags:**
   `grep -rniE 'bug|fixme|wrong|hack' interaction/ algorithm/ math/`
   The gromosXX developers annotate their own defects and uncertainties. Two minutes here catches the
   issues they already knew about.
2. **Classify every suspicious thing you find** into one of three buckets, and treat each differently:
   - **(A) Genuine bug** → do *not* port it; note it so a future reader knows it was deliberate.
   - **(B) Intentional GROMOS quirk** → reproduce it *exactly* and **document it as a deliberate
     decision** (in the style of CHANGELOG's "BONDANGLETYPE intentionally unsupported" entries), so it
     is never later mistaken for a porting mistake.
   - **(C) Latent divergence** (our code already differs, nothing tests it yet) → write the reference
     test, then reconcile.

## Second-source the subtle physics

For anything the authors flagged as uncertain, or anything intrinsically delicate (reaction-field
self-terms, virial, Ewald/lattice-sum, perturbed/soft-core terms): **derive it independently from
the GROMOS book / primary paper, implement to that, then diff against the C++.** Where they disagree,
investigate — do not assume the C++ is right. Transcribing uncertain C++ line-by-line is how you
inherit a bug you can't see.

Free leverage: gromosXX's own `check/*.t.cc` files hard-code per-term reference energies (including
perturbed/soft-core variants). They are a second source of truth independent of running the md
binary — use them to validate the very terms you were told not to trust.

## Reproduce-or-correct, made explicit

When GROMOS behavior and "physically correct" behavior diverge, you may eventually want both: a
faithful compatibility mode and a corrected mode. Design so the two can coexist behind a flag rather
than forcing an irreversible choice. Until then, **faithful is the default** and any deviation is
named and documented.

### The list

Each entry: what gromosXX does, which bucket, what gromos-rs does.

- **PRECALCLAM NRLAM = 1** (A). gromosXX accepts it and divides by `NRLAM − 1`: every λ_k is
  NaN. gromos-rs refuses it at build time (`check_precalclam`, `gromos-run/src/plan.rs`).
- **PRECALCLAM long-range ext energies** (B). gromosXX computes the long-range perturbed pairs
  only on pairlist steps and keeps them in `m_longrange_storage` until the next one; their ext
  (λ_k) energies are filled only when that pairlist step is also an NTWG step, and are zero on
  the output steps after one that is not (`perturbed_nonbonded_set.cc:87-96, 280-300`).
  Reproduced (`Forcefield::longrange_perturbed`), pinned by `ch4_water_fep_precalc_twinrange`.
- **PRECALCLAM current-λ leaks** (B). The grid uses the raw λ_k everywhere, except: the kinetic
  term's `slam` (the current mass λ), the improper's ζ wrapped around the current-λ q0 and not
  re-wrapped, perturbed SHAKE adding to the `AB_bond` derivatives only, the 1-4 ext call taking
  no Coulomb scaling, and the RF self ext term carrying no self-term correction. All reproduced.

## In one line

Tests before wiring; derivation before transcription; grep before porting; document every deliberate
quirk. Faithfulness and correctness are both required — the discipline above is how you get both
instead of trading one for the other.
