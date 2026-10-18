**GaMD reads gromosXX's GAMD block and runs gromosXX's algorithm.** A one-parameter toggle sweep
found that newen read the `GAMD` block in a layout of its own (`SEARCH FORM THRESH` /
`SIGMA0 ×2` / `K, E ×4`, FORM 1 = dihedral, 3 = dual): on a real gromosXX input no line matched,
the run went on with newen's defaults, and no boost was applied. Its boost, applied after the
step, was not gromosXX's either.

### Fixed

- **The GAMD block** (`GamdBlock::from_lines`, `newen-io/src/input/gamd_block.rs`) is read as
  gromosXX's `read_GAMD` reads it (`in_parameter.cc:4148-4320`): one token stream, GAMD,
  SEARCH FORM THRESH NTIGAMDS, AGROUPS, IGROUPS, DIHSTD TOTSTD, then IGROUPS values each of ED,
  ET, KD, KT, then EQSTEPS and GWINDOW, every parameter with its range and gromosXX's message
  (`GAMD block: FORM parameter is 4, should be:  one of: 1,2,3`, "reached END before",
  "Left-over parameters"). FORM 1 = dual boost, 2 = potential-energy boost, 3 = dihedral boost;
  SEARCH 0 = fixed parameters, 1 = cMD search, 2 = GaMD search; THRESH 1 = lower, 2 = upper bound.
  A malformed block now stops `md` (it was a warning, and the run went on without GaMD). `md`
  reads it from the parsed `.imd` (`raw_blocks`), not by re-reading the file. `render` writes
  gromosXX's layout and parses back unchanged, also through `write_imd`.
- **The `@gamd` file** (`GamdGroups::parse`): GAMDATOMS and GAMDGROUPS with gromosXX's checks and
  messages (`in_gamd.cc:76-238`); a repeated GAMDGROUPS pair keeps its first entry, as
  `std::map::insert` does.
- **The algorithm** (`newen_integrators::gamd::Gamd`, run at the end of `Forcefield::apply`
  through `Forcefield::gamd`, gromosXX's place right after the force field,
  `create_md_sequence.cc:143-146`): gromosXX's `GAMD::init` and `GAMD::apply` (`gamd.cc`). V_D is
  the proper and improper dihedral energy, V_T the bond, angle and nonbonded energy; the
  statistics (Welford over `stepsdone`, EQSTEPS, GWINDOW resets) for SEARCH 1/2, E and k from
  them for SEARCH 2 and for KD = KT = 0 at start, then except in a cMD search the boost
  ΔV = k(V − E)²/2 for V < E with F += k(V − E)·F_term and the same scaling of the term's virial
  (the nonbonded part only with a virial, `nonbonded_set.cc:431-441`). ΔV is booked in the
  special energy, ENER[51] and the `# gamd` rows. The old after-step boost in `md` (and its
  `gamd_stats.dat` / `gamd_boost.dat`) is gone.
- **The `# gamd` section** of the `.tre`: IGROUPS + 1 rows of ED, ET, KD, KT, ΔV per interaction
  group (`out_configuration.cc:3202-3218`) whenever a GAMD block is read, GAMD on or off
  (`configuration.cc:422-444`); it was always one zero row. `write_energyred_body` takes the
  rows (`EnergyWriter::set_gamd_rows`) and prints their ΔV sum as ENER[51].
- **GAMDSTAT**: a search (SEARCH 1/2) appends the statistics to the final configuration
  (`out_configuration.cc:619-621, 4013-4034`, `%15.9f`), and with NTIGAMDS = 0 starts from the
  input configuration's block (`in_configuration.cc:1458-1507`): "no GAMDSTAT block in
  configuration." without one, as gromosXX.

### Reproduced gromosXX defects (flagged, not corrected)

- **FORM 1 reports the dihedral ΔV only.** `GAMD::apply` assigns `gamd_DV` for the potential
  boost and then again for the dihedral boost (`gamd.cc:364, 381`): when both act, both scale the
  forces but ΔV (special energy, ENER[51], `# gamd`) is the dihedral one, so the reported energy
  is not the one whose gradient the forces are. Pinned by `aladip_gamd_dual`.
- **The reaction-field excluded solute pairs leave the accelerated group.**
  `RF_excluded_interaction_innerloop` declares `gamdi = 0` and sets a shadowing `gamdi` inside the
  `if (gamd)` of the self term (`nonbonded_innerloop.cc:1694, 1711`), so every excluded pair is
  looked up with the key (0, group of j) (`:1721-1724`): with GAMDGROUPS not listing (0, g) its
  energy and forces go to the never-boosted group 0, while the self term stays in g. On
  `aladip_gamd_dual` this moves 50.19 kJ/mol out of V_T (−85.78 instead of −35.59 kJ/mol), a
  2.41× larger potential boost. newen books those pairs to the group gromosXX does
  (`Gamd::excluded_group`, `GamdGroups::excluded_pair_group`).

### Changed

- **Refused instead of run wrongly** (`md`, `gamd_setup`): GAMD = 1 with a GAMDATOMS/GAMDGROUPS
  assignment that books the interactions to more than one interaction group (the nonbonded
  pair terms are not split by group), with a perturbation (NTG = 1), with NLRELE ≠ 1, or with
  QM/ML terms; gromosXX runs them. An `@gamd` atom range past the last atom is refused
  (gromosXX's `iatom > num_atoms` check after the 0-based shift lets FINATOM = N + 1 through and
  writes past the array, `in_gamd.cc:138`).

### Added

- References (native gromosXX): `aladip_gamd` is no longer ignored (FORM 3, the dihedral boost);
  `aladip_gamd_dual` (FORM 1, ET 0: both boosts act) and `aladip_gamd_cmd_search` (SEARCH 1,
  NTIGAMDS 1: statistics, no boost, `GAMDSTAT`). The comparator checks the final
  configuration's GAMDSTAT block (to 1e-8). 125 of 125 pass.
