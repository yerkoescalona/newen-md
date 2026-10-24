# Energy Library Profile

**Profile version 1 — specified 2026-09-03.**

The energy library (`ene_ana.md++.lib`) is the file that tells `ene_ana` how to read an energy
trajectory (`.tre`) or a free-energy trajectory (`.trg`), and which quantities to derive from
it. GROMOS deliberately made it a *file* rather than code: a researcher who re-partitions a
system into different energy groups, recomputes the energies with a new topology, and wants to
reconstruct a log with different sections, can describe the new file without recompiling
anything. That flexibility is the point of the design and this profile keeps all of it.

What the profile adds is rigor. The upstream format has one weakness: nothing ties a library to
the trajectory it is applied to. A library that declares the wrong shape for a block reads every
value into the wrong slot and reports plausible-looking nonsense — a total kinetic energy of
−52 kJ/mol — with no error. The only guard, the `ENEVERSION` string, is hand-typed in both files
and two different layouts in the wild carry the same one.

The profile therefore states:

> **Same grammar, stricter rules.** Every file newen writes is a valid GROMOS file that
> gromos++ and gromosXX read unchanged; every file newen reads is checked against a
> fingerprint of the layout it was written with, and a library that disagrees with that
> fingerprint is refused with a diff — never applied.

## Status

| Part | State |
|---|---|
| §1 Grammar, reader, expression evaluation | implemented (`newen_io::energy_traj`, 0.0.46), byte-identical to gromos++ on the reference cases |
| §2 Canonical form, fingerprint, self-description lines | implemented (0.0.47): `Schema::fingerprint`, `md` writes R3 and R8 lines, the fingerprint of the layout in code equals md++'s own library |
| §3 Reader tiers and refusals | implemented (0.0.47): `EnergyTraj::bind`; the LiveCoMS tutorial's stale library is refused at open with the diff |
| §4 Official library from `newen-io` | implemented (0.0.47): `crates/newen-io/data/ene_ana.md++.lib`, generated from the layout in code, `ene_ana` uses it without `@library`, `@print_library` prints it |

Implementation landed in 0.0.47 (item 2.10; see CHANGELOG.md).

## 1. The grammar — unchanged from gromos++

### 1.1 File structure

A library is a sequence of GROMOS blocks, each a name line, content, and `END`. Lines whose
first non-blank character is `#` are comments; text from `#` to end of line is stripped
everywhere (gromos++ `Ginstream.cc:165-172`).

| Block | Content |
|---|---|
| `TITLE` | free text |
| `ENEVERSION` | one token, the layout version the library was written for (`2023-04-15`) |
| `ENERTRJ` | schema of an energy trajectory (`@en_files`) |
| `FRENERTRJ` | schema of a free-energy trajectory (`@fr_files`) |
| `VARIABLES` | derived quantities, one `name = expression` per line |

gromos++ ignores blocks it does not know (`ene_ana.cc:404-423`); so does newen.

### 1.2 Schema declarations

Inside `ENERTRJ` / `FRENERTRJ`, three keywords describe how the values in a trajectory block are
consumed, in order:

```
block    NAME               the next GROMOS block in the trajectory must be NAME
size     NAME               read one integer from the stream, bind it as NAME
subblock NAME ROWS COLS     read ROWS × COLS values into the table NAME
```

`ROWS` and `COLS` are literal integers, a bound `size` name, or `matrix_<size>` which expands to
`N(N+1)/2` — the upper triangle of a symmetric group–group matrix. A `size` name refers to one
slot: declaring `size NUM_BATHS` in `ENERGY03` and again in `VOLUMEPRESSURE03` re-reads the
same slot, so a later subblock may use a size bound in an earlier block.

The values of a block are one whitespace-separated stream; declarations consume it left to
right. Fewer values than declared is an error (`Not enough values in NAME`); values left over
when the declarations are exhausted is an error (`Block definition does not agree with
trajectory data for NAME`). **A declaration that consumes the right number of values in the
wrong shape is not an error** — that is the gap §2 closes.

### 1.3 Derived quantities

`VARIABLES` holds one definition per line, `name = expression`. An expression is a sequence of
whitespace-separated tokens — operators must be surrounded by blanks, as in the shipped library:

| Token | Meaning |
|---|---|
| `TABLE[i]`, `TABLE[i][j]` | a value of a subblock, **1-based**; `[i]` on a one-column table |
| `MASS`, `NUMMOL` | constants from `@topo`: total mass, number of solute molecules |
| `BOLTZ` | Boltzmann's constant in GROMOS units, built in |
| `TIME[2]` | the current time (`TIME` is the 2 × 1 subblock of `TIMESTEP`) |
| number | a literal |
| `+ - * /` | binary operators |
| `sin cos exp log` | functions, applied to the whole expression that follows |
| `name` | another variable, evaluated on demand |

Precedence reproduces `gmath::Expression` exactly: an expression is split at the **rightmost**
`+` or `-`, else at the **leftmost** `*` or `/`. Addition and subtraction therefore associate
left, multiplication and division associate *right* — `a / b / c` is `a / (b / c)`. Use blanks
and reordering rather than parentheses; the grammar has none.

### 1.4 The canonical md++ layout — `ENEVERSION 2023-04-15`

The layout gromosXX has written since 2023 and newen's `md` writes (`newen_io::energy`):

`ENERTRJ` (`.tre`)

| `block` | declaration | rows × cols | meaning |
|---|---|---|---|
| `TIMESTEP` | `subblock TIME 2 1` | 2 × 1 | step, time |
| `ENERGY03` | `subblock ENER 52 1` | 52 × 1 | totals: total, kinetic, potential, covalent, bond, angle, improper, dihedral, crossdihedral, nonbonded, LJ, CRF, … |
| | `size NUM_BATHS` | | |
| | `subblock KINENER NUM_BATHS 3` | N_b × 3 | per bath: total, centre-of-mass, internal kinetic energy |
| | `size NUM_ENERGY_GROUPS` | | |
| | `subblock BONDED NUM_ENERGY_GROUPS 5` | N_g × 5 | per group: bond, angle, improper, dihedral, crossdihedral |
| | `subblock NONBONDED matrix_NUM_ENERGY_GROUPS 6` | N_g(N_g+1)/2 × 6 | per group pair: LJ, CRF, lattice-sum total, lattice-sum k-space, shift extra (orig, phys) |
| | `subblock SPECIAL NUM_ENERGY_GROUPS 13` | N_g × 13 | per group: constraints, position / distance / distance-field / angle / dihedral restraints, SASA, SASA volume, RDC, … |
| | `size NUM_EDS_STATES` / `subblock EDS NUM_EDS_STATES 6` | N_eds × 6 | |
| | `size NUM_GAMD_GROUPS` / `subblock GAMD NUM_GAMD_GROUPS 5` | N_gamd × 5 | |
| | `size NUM_LAMBDAS` / `subblock PRECALCLAM NUM_LAMBDAS 12` | N_λ × 12 | |
| | `subblock ABDIH 1 2` | 1 × 2 | |
| `VOLUMEPRESSURE03` | `subblock MASS 1 1` | 1 × 1 | |
| | `size NUM_BATHS` / `subblock TEMPERATURE NUM_BATHS 4` | N_b × 4 | per bath: total, COM, internal temperature, scaling factor |
| | `subblock VOLUME 10 1` | 10 × 1 | volume, box vectors |
| | `subblock PRESSURE 30 1` | 30 × 1 | pressure, virial and kinetic-energy tensors |

`FRENERTRJ` (`.trg`) is the same shape with `FREEENERDERIVS03` in place of `ENERGY03`, an
extra leading `subblock RLAM 1 1` (λ), the subblocks prefixed `FREE`, and no
`VOLUMEPRESSURE03`.

Any layout newen knows how to write or check is listed in this section; a new upstream
layout gets a new row group under its own `ENEVERSION`, never an edit to an existing one.

## 2. Profile rules

Normative; "must" is a requirement on newen, "may" an option left open. Rationale follows
each rule in one sentence.

### R1 — Compatibility

Every trajectory and library newen writes must be read unchanged by gromos++ `ene_ana` and,
for trajectories, be indistinguishable to gromosXX from its own output. All profile additions
are lines whose first character is `#`.

*Why:* the reference suite — gromosXX and gromos++ run on the same inputs — is the project's
oracle, and the LiveCoMS tutorials are the proof that real GROMOS workflows run; a new block
would break gromos++, which reads blocks strictly in library order (`EnergyTraj.cc:454-457`,
"Block X expected. Got Y").

### R2 — Canonical form and fingerprint

The **canonical form** of a schema (`ENERTRJ` or `FRENERTRJ`) is the text produced by:

1. one line per declaration, in order, keyword and operands separated by one blank;
2. subblock names dropped: `subblock ENER 52 1` → `subblock 52 1`;
3. size names replaced by `$k`, numbering distinct names by first occurrence from 1:
   `size NUM_BATHS` → `size $1`, `subblock KINENER NUM_BATHS 3` → `subblock $1 3`,
   `matrix_NUM_ENERGY_GROUPS` → `matrix($2)`;
4. block names kept;
5. lines joined with `\n`, terminated by `\n`, no other whitespace.

The **fingerprint** is `sha256:` followed by the lowercase hex SHA-256 of that text.

For the 2023-04-15 `ENERTRJ` layout the canonical form begins:

```
block TIMESTEP
subblock 2 1
block ENERGY03
subblock 52 1
size $1
subblock $1 3
size $2
subblock $2 5
subblock matrix($2) 6
subblock $2 13
```

*Why:* two libraries that read the bytes identically get the same fingerprint, two that read
them differently cannot; names are dropped because renaming a table or a size changes nothing
about which byte lands where, and the profile must leave that free. SHA-256 so that anyone can
verify by hand with `sha256sum`.

### R3 — Trajectory self-description

A trajectory written by newen must carry, once, immediately after the `END` of its
`ENEVERSION` block and before the first `TIMESTEP`:

```
# energy-schema 1 ENERTRJ sha256:<64 hex>
# energy-layout block TIMESTEP
# energy-layout subblock TIME 2 1
# energy-layout block ENERGY03
# energy-layout subblock ENER 52 1
# energy-layout size NUM_BATHS
…
```

- `energy-schema` names the profile version, the file type, and the fingerprint of the layout
  that follows.
- one `energy-layout` line per declaration, **with names** — this is the full schema, from
  which a library can be generated (§4).
- the `ENEVERSION` block is written exactly as gromosXX writes it and is not replaced.

The fingerprint must equal the fingerprint of the canonical form of the `energy-layout` lines;
a reader that finds them unequal must refuse the file as corrupt.

*Why:* the file becomes the authority on its own layout; a library can no longer claim
something about a file that the file can contradict. Written once because the layout cannot
change within a file (a `size` may, its shape may not).

### R4 — Library self-description

A library may carry, anywhere outside a block:

```
# energy-schema 1 ENERTRJ sha256:<64 hex>
# energy-schema 1 FRENERTRJ sha256:<64 hex>
```

If present, each must equal the fingerprint of the library's own declaration of that file type.
A library newen generates always carries them. A reader that finds one unequal must refuse
the library: *"schema section edited after generation"*.

*Why:* this catches the tutorial's failure mode exactly — a copied library whose schema drifted
while its `ENEVERSION` still claims the original — as an error at open, not as wrong numbers
at the end.

### R5 — Reader tiers

Before reading a frame, the reader establishes the file's layout and checks the library against
it. Exactly one tier applies, the first that matches:

| Tier | Condition | Authority for the layout | Library schema … |
|---|---|---|---|
| **a** | file carries `energy-schema` (R3) | the file | must match the file's fingerprint → else **error with a diff** |
| **b** | no self-description; `ENEVERSION` is one newen knows (§1.4) | the built-in layout for that version | must match its fingerprint → else **error with a diff** |
| **c** | neither | the library | accepted with the warning *"layout unverified"*; R6 checks are the only guard |

The diff names each differing declaration on both sides:

```
energy library does not describe this trajectory
  ENERGY03 ENER:   library 26 x 2,   file 52 x 1
  ENERGY03 SPECIAL: library NUM_ENERGY_GROUPS x 12,   file NUM_ENERGY_GROUPS x 13
  library:    tutorial/ene_ana.md++.lib  (ENEVERSION 2023-04-15)
  trajectory: md.tre  (ENEVERSION 2023-04-15, written by newen md 0.0.47)
```

There is no option to override tiers **a** or **b**. The layout is a fact about the bytes; the
only thing legitimately different between libraries is `VARIABLES`, which no tier examines.

*Why:* tier **a** is the guarantee for our own files; tier **b** is the best available for
gromosXX's, bounded by upstream's version string not being unique — a `.tre` from an older
md++ that also says `2023-04-15` but wrote 50 totals is caught only by R6. Tier **c** keeps
the reconstruct workflow open for files written by anything else.

### R6 — Structural checks, every tier

Independent of the library's origin, the reader must refuse a frame when:

- a `size` value is not a non-negative integer, or differs from the same size in the first
  frame;
- the number of values a subblock would consume exceeds the values remaining in the block
  (checked before allocating);
- the block carries marker comments (`# totals`, `# baths`, … — gromosXX writes them in every
  `ENERGY03`) and a marker falls *inside* a declared subblock's run of values rather than
  between two declarations.

And on the first frame after binding, it must **warn** when:

- `ENER[1] ≠ ENER[2] + ENER[3] + ENER[21]` beyond 1e-6 relative — md++ composes the total
  energy as potential + kinetic + special (`energy.cc`, `Energy::calculate_totals`), and
  writes the special total at slot 21 of `ENER`;
- `ENER[2] < 0`;
- `VOLUME[1] ≤ 0` when the block is present.

The warning says why it matters for the tier in force: under **c** "the library's layout is
unverified and may mislabel the values"; under **a**/**b** "the layout is verified, so the
writer composed these totals inconsistently". It is a warning, not a refusal, because a
verified layout is a fact about where the values sit while the identity is a fact about how
the writer computed them: a file whose slots are proven right must stay readable even if its
writer summed badly (newen' own `.tre` files before 0.0.47 wrote the special total two
slots early and would otherwise be unreadable).

*Why:* the first two turn an allocation or a silent shift into an error; the marker check
catches every mis-sized subblock in a gromosXX file; the invariants are the only check that
catches a **reshape** (same value count, wrong slots) in a file with no fingerprint, and they
hold for any partitioning of the system.

### R7 — The official library

`newen-io` provides the library for every layout it knows (§1.4): the schema sections are
generated from the layout definition in code, the `VARIABLES` section is md++'s. `ene_ana`
uses it when `@library` is absent. A copy is checked in at `crates/newen-io/data/` for people
who need the file — for gromos++, for editing — and a test asserts the copy equals what the
code generates, so the two cannot drift.

A user library may:

- carry a full schema (gromos++-compatible) — checked by R5;
- carry **only** `VARIABLES` — the schema then comes from tier **a** or **b**; such a library
  is newen-only, since gromos++ needs `ENERTRJ`.

*Why:* the thing anyone can edit is never the thing the reader trusts; the source of truth is
code, and the file is derived from it and verified against it.

### R8 — Provenance

A trajectory written by newen should carry, after the layout lines:

```
# energy-provenance writer newen md 0.0.47
# energy-provenance topology sha256:<hex>   <path as given>
# energy-provenance energy-groups 1-63 64-3000
```

Readers ignore `energy-provenance` keys they do not know. `ene_ana` reports writer and
topology in its header and in every R5 diff.

*Why:* the re-partition workflow's real question, a year later, is "which grouping was this
file written with?" — today nothing anywhere records it.

### R9 — Versioning

- `energy-schema N`: the **profile** version — the rules in this chapter. It changes only when
  the canonical form, the self-description syntax, or a reader obligation changes. A reader
  must refuse a profile version it does not know.
- `ENEVERSION`: the **upstream layout** version, unchanged in meaning and written exactly as
  gromosXX writes it.
- A new upstream layout adds a new entry to §1.4 and to the code, under its own `ENEVERSION`;
  existing entries are never edited.

*Why:* the two versions answer different questions — "how do I verify this file?" and "which
bytes are in it?" — and conflating them is how the upstream check became a no-op.

### R10 — What stays free

The profile does not constrain, and must not be extended to constrain:

- the number of energy groups, baths, EDS states, λ points — these are `size` values in the
  stream, read per file;
- subblock and size **names** — dropped from the canonical form on purpose;
- everything in `VARIABLES`;
- new file types or blocks, given a full schema for them (tier **c**).

The one thing removed is the ability to be silently wrong.

## 3. Refusals

| Condition | Message (first line) | Rule |
|---|---|---|
| library fingerprint ≠ file fingerprint | `energy library does not describe this trajectory` | R5a |
| library fingerprint ≠ built-in for `ENEVERSION` | `energy library does not match the <version> layout` | R5b |
| file's `energy-schema` ≠ hash of its `energy-layout` lines | `trajectory self-description is corrupt` | R3 |
| library's `energy-schema` ≠ hash of its own schema | `library schema section edited after generation` | R4 |
| unknown `energy-schema` profile version | `energy-schema <N> is newer than this reader` | R9 |
| `size` not an integer / negative / changed | `size NAME: expected an integer, got …` | R6 |
| subblock overruns the block | `subblock NAME: needs R x C values, block has …` | R6 |
| marker inside a subblock | `"# baths" falls inside subblock ENER` | R6 |
| `size` in a later frame differs from the first | `size NAME: 3 in this frame, 2 in the first` | R6 |
| `size` larger than what is left of the block | `size NAME: 40 rows, but only 13 values left in the block` | R6 |
| no way to know the layout and no library | `no way to establish the layout` | R5 |

Warnings, never refusals: tier **c** (`layout unverified: <file> (ENEVERSION …) is read as
<library> declares it`); the first-frame invariants (`ENERTRJ first frame: ENER[1] is not
ENER[2] + ENER[3] + ENER[21]: … vs …; …`, R6); `ENEVERSION` present in one of library and file
but not the other, or different (gromos++'s `Version number mismatch!`/`missing!` text, kept).
Warnings go to stderr; the analysis output is unchanged.

## 4. Workflows

**Analyse a newen or gromosXX run** — no library argument:

```sh
ene_ana @en_files md.tre @prop totene totkin totpot @topo sys.top
```

**Add your own derived quantities** — a `VARIABLES`-only library, or generate the full one and
edit its `VARIABLES`:

```sh
ene_ana @print_library > mine.lib          # official library for the built-in layout
ene_ana @print_library md.tre > mine.lib   # the schema md.tre says it was written with
ene_ana @en_files md.tre @library mine.lib @prop solute_intra
```

The generated file carries R4 fingerprints; editing a `VARIABLES` line leaves them valid,
editing a `subblock` line makes the file refuse until the fingerprint line is regenerated —
which is the point.

**Re-partition and reconstruct.** Recompute with a topology whose energy groups are the ones you
now want; the writer records the new `NUM_ENERGY_GROUPS` in the stream and the grouping in
`energy-provenance`; the same library reads the result. Only if the recomputation produces a
differently *shaped* file — a tool of your own — do you need a schema, and then `@print_library
newfile.tre` writes it for you from the file's self-description.

## 5. Decisions

Options considered for "who owns the official library":

| Option | Verdict |
|---|---|
| Ship the file, locate at runtime (gromos++'s `data/` model) | editable in place, nothing verifies it — the status quo this profile replaces |
| Embed the file text in the binary | immutable without a rebuild, but the writer and the embedded text can still disagree |
| **Layout defined in code; library, fingerprints and writer check derived from it** | chosen — writer, reader and library cannot disagree, and the file is a verified artifact |
| Signed distribution | the threat is accident and drift, not an adversary; a fingerprint line already makes any edit visible |

Other decisions: comments rather than a new block (R1, forced by gromos++'s strict block order);
names outside the fingerprint (R2, so that renaming stays free); SHA-256 rather than a short
in-tree hash (verifiable by hand); no override for tiers **a**/**b** (there is no legitimate
case); tier **c** kept (the reconstruct workflow needs it).

The derived-quantity language is deliberately **not** extended beyond gromos++'s. A library
that uses a function gromos++ lacks runs on one implementation only; anyone who needs more
than `+ - * / sin cos exp log` over the tables should take the tables into Python.
