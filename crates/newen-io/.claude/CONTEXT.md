# newen-io — stage contract

## Job
L0 I/O service. All file format parsing and writing. The only place in the workspace where
GROMOS file formats are read or written.

## Inputs (consumes from)
`newen-core` types only.

## Outputs (public API)
Parsers + writers for every GROMOS file format used in the workspace.

## Status
- Topology parser ✓: all blocks (SOLUTEATOM, SOLVENTCONSTR, LJPARAMETERS, PRESSUREGROUPS, ...)
- Coordinate reader/writer ✓: POSITION, POSITIONRED, VELOCITY (labelled or three columns),
  VELOCITYRED, GENBOX, PERTDATA, AEDSSEARCH (`parse_aedssearch`, `G96Writer::write_aedssearch`),
  NHCVARIABLES (`parse_nhcvariables`, `G96Writer::write_nhcvariables`, gromosXX's `%15.9f` rows)
- IMD: the `AEDS` block (`AedsBlock`, `parse_aeds`, written by `write_imd`); ENERGY03 slots 34-39
  and the `# eds` offset column carry the A-EDS fields
- Trajectory ✓: `TrajectoryReader` (3-column and legacy 7-column POSITIONRED; GENBOX optional for vacuum),
  `TrajectoryWriter::write_frame` (via `write_trc_frame`), `write_trc_frame(step, time, positions, box)`
- Energy/force writers ✓: ENERTRJ, FREEFORCERED/CONSFORCERED
- PDB: `write_pdb_positions(path, title, positions, box, topology)` ✓
- G96: `write_g96(path, title, positions, velocities, box, topology)` ✓
- Free energy: `FreeEnergyWriter` (FREEENERGY03), `read_free_energy_trajectory` ✓
- Energy library profile ✓ (`docs/src/reference/energy-library.md`): `energy_traj.rs` holds the
  2023-04-15 ENERTRJ/FRENERTRJ layouts in code, fingerprints them (`Schema::fingerprint`), reads
  ene_ana `.lib` files and `.tre`/`.trg` blocks against them (`EnergyTraj::bind`, tiers a/b/c), and
  generates the official library `data/ene_ana.md++.lib` (`OFFICIAL_LIBRARY`, `include_str!`).
  Writers self-describe (`# energy-schema`/`# energy-layout`) and record provenance.
- MTB, IFP, IMD, PTP parsers ✓; position restraints (.por/.rpr) ✓

## Key files
```
src/topology.rs    — topology file parser (~1200 LOC; P4 split candidate)
src/trajectory.rs  — TrajectoryReader + TrajectoryWriter (write_trc_frame)
src/g96.rs         — write_g96, G96Writer
src/pdb.rs         — write_pdb_positions
src/free_energy.rs — FreeEnergyWriter, read_free_energy_trajectory
src/energy_traj.rs — Schema, EnergyTraj (library-driven .tre/.trg reader), official library
data/ene_ana.md++.lib — generated; `UPDATE_OFFICIAL_LIBRARY=1 cargo test -p newen-io --lib energy_traj` rewrites it
src/imd.rs         — IMD parameter file parser
```

## Decisions logged
- **POSITIONRED format:** standard is 3-column (x y z); reader accepts both 3-col and legacy 7-col.
  `write_trc_frame` always writes 3-column (standard GROMOS format).
- **GENBOX optional:** vacuum trajectories omit GENBOX; reader returns `Vec3::ZERO` when missing.
- **PERTURBATION block is multi-line:** splits NTG/NRDGL/RLAM/DLAMT on line 1 and ALPHLJ/ALPHC/NLAM/NSCALE
  on line 2. Always combine `data_lines` before tokenising. Bug: `data_lines.first()` silently dropped line 2.
- **Energy output:** full f64 scientific notation for exact comparison.
- **NLRELE values (2026-10-23):** −1 reaction field (Vol. 4's "LSERF compatibility mode"; gromosXX
  reads it exactly as 1, `in_parameter.cc:4641-4648`), 0 none, 1 reaction field, 2 Ewald, 3 P3M,
  4 shifted reaction field (`in_parameter.cc:4663-4674`), 5 SPME (newen only, experimental). Why 5
  for SPME: gromosXX's own 4 is the shifted reaction field, and its reader refuses everything
  outside −1..4 (`in_parameter.cc:4604`), so a gromosXX input can never be read as SPME. Vol. 4's
  notes "abs(NLRELE) = 3,4" (`vol4.tex:7785`, `:7792`) are leftovers of an older reading (in
  the dossier's manual and comment issues, no id); gromosXX's enum has `el_spme`
  (`simulation/parameter.h:480`), which nothing sets; NSPORD is read into `spme_bspline`
  (`in_parameter.cc:4624`) and used nowhere. newen reads NSPORD into
  `NonbondedExtra::nspord` as before; only NLRELE 5 uses it.
- **Energy layout lives in code, not in the data file:** `data/ene_ana.md++.lib` is generated from
  `energy_traj::LAYOUTS` and a test keeps the checked-in copy byte-equal. A layout change is a new
  ENEVERSION, never an edit of the 2023-04-15 entry; a library whose `# energy-schema` line no longer
  matches its blocks is refused.

## Crate-specific rules
- **ALL parsing/writing lives here.** No file format code in binaries or other crates.
- **All analysis binaries** use `write_g96`, `write_pdb_positions`, `TrajectoryWriter::write_trc_frame`
  — never hand-roll `BufWriter`/`writeln!` for GROMOS formats in binary entry points.
