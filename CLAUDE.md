# newen

Project guidance lives under `.claude/`. Imported here so it loads automatically:

- `@.claude/overview.md` — workspace router: crate map, prime directives, global decisions, how to verify.
- `@.claude/architecture.md` — how to shape new and refactored code.
- `@.claude/porting.md` — how to translate gromosXX faithfully without inheriting its defects.

When working inside a specific crate, also read its stage contract:
`crates/<crate>/.claude/CONTEXT.md` — job, inputs/outputs, current status, key files, crate rules.

When touching code with a `/// spec:` line or anything in `crates/newen-logic` or `spec/`, read
`docs/src/verification/index.md` (the Lean specification: what is proven, what is tested, how to
add a function).
