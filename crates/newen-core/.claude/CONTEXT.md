# newen-core — stage contract

## Job
L0 data core. Owns `(Topology, Configuration)` shared by everything. No physics, no I/O.

## Inputs (consumes from)
`newen-logic` only (the proven integer decisions: pair kinds, `due` for `PrecalcLambdas::active`).

## Outputs (public API)
`Topology`, `Configuration`, `Vec3`/`Mat3`, `Periodicity` (Vacuum/Rectangular/Triclinic),
`Algorithm` trait + `AlgorithmSequence`, `AtomSelection`, `Stat`,
`PairlistAlgorithm` (enum: Standard/CellList), `StandardPairlistAlgorithm`, `CellListPairlistAlgorithm`,
`gather_chain`, `gather_bond`, `gather_molecules` (PBC gathering),
`SpatialIndex` trait + `ConfigurationSpatialIndex` (query-based neighbor service for the
`PotentialProvider` seam, 2.6).

## Status
- Boundary: Vacuum ✓, Rectangular ✓, Triclinic ✓ (NTB=-1 truncated-octahedron, and since
  2026-10-22 GENBOX's NTB 2 cell with angles: `genbox_box` builds the box vectors,
  `genbox_lengths_angles` writes them back with the `acos` argument clamped (D57),
  `box_heights` gives the cell heights the run's cutoff check uses)
- Boundary kernels checked against a Lean model (2026-10-22): `vacuum_nearest_image`,
  `Rectangular::new`, `rectangular_nearest_image`, `rectangular_put_into_box`, `Triclinic::new`,
  `triclinic_nearest_image`, `triclinic_put_into_box` (pure functions in `math.rs`, called by the
  `BoundaryCondition` impls; `/// spec: periodicity.*`). ℝ theorems in
  `spec/NewenTheory/Space/Periodicity.lean`, model in `spec/NewenSpec/Periodicity.lean`, bit-exact
  test `tests/periodicity_spec.rs` on `spec/vectors/periodicity.json`. Keep the operation order
  of these functions: any change shows up as differing bits there; regenerate the vectors only
  after changing the Lean model to match.
- Topology: the instancing model complete — `moltypes[0]` = solute, `moltypes[1..]` = solvent types,
  `instances[k].role` = Solute/Solvent; flat arrays derived from instances; `Solvent` struct removed;
  `Solute` is a ZST shell (atoms moved to moltypes). Direct field access: `topo.moltypes[0].bonds`.
- AtomSelection: full newen grammar ✓ — `a:name`, `1:res(N:CA)`, `not()`, `minus()`, `;`-union,
  `m:` (solute by role), `s:` (solvent by role), syntax-error hints. `a:` is every *solute*
  molecule, never solvent (gromos++ `AtomSpecifier`), so on a solvated system `a:OW` is empty
- PBC gathering: `gather.rs` — `gather_chain`, `gather_bond`, `gather_molecules` ✓
- Statistics: `stat.rs` — `Stat` with `ave()`, `rmsd()`, `ee()` (block averaging) ✓
- Perturbation net charge (1.8, 0.0.59): `NetCharge` (A/B, `delta()`, per molecule);
  `Topology::apply_perturbation` returns the solute's, `perturbed_molecule_charges` per molecule
- Pair tests: `Topology::is_excluded`, `pair_kind_of`, `is_excluded_or_14` (every pairlist
  builder) call `newen_logic::topology::pairs::{contains_sorted, pair_kind}` and re-export
  `PairKind`; `exclusions` and `one_four_pairs` are sorted and deduplicated.
- Pairlist: `PairlistAlgorithm` enum wired everywhere (9a-0 ✓); `CellListPairlistAlgorithm` bit-identical to Standard (9a-1 ✓, margin=0); `water_1000_spc_gridcell` reference validates against gromosXX Grid_Cell_Pairlist
- `SpatialIndex`: query-based neighbor service (`neighbor_pairs`, periodic-image shift per pair),
  independent of the MD pairlist's charge-group/twin-range/solvent-sentinel shape. Used by every
  `PotentialProvider` in `newen-forces` (LJ+CRF, electrostatic embedding, SchNet) — not by the
  production MD pairlist itself, which still uses `PairlistAlgorithm` above.

## Key files
```
src/topology.rs   — Topology, MoleculeType, MoleculeInstance, Role
src/selection.rs  — AtomSelection, full grammar + error hints
src/gather.rs     — PBC molecule gathering
src/stat.rs       — Stat, block-averaging ee()
src/pairlist.rs   — StandardPairlistAlgorithm, CellListPairlistAlgorithm
src/math.rs       — Vec3, Periodicity, the boundary kernels, truncoct_triclinic_*
tests/periodicity_spec.rs — boundary kernels = Lean model, bit for bit
src/spatial_index.rs — SpatialIndex trait, ConfigurationSpatialIndex
```

## Crate-specific rules
- **No physics.** Energy/force belongs in newen-forces.
- **Instancing convention:** `moltypes[0]` is always the solute (GROMOS `moltype[0]` convention).
  Direct field access: `topo.moltypes[0].bonds.push(...)` — no accessor methods.
- **Role-based dispatch:** `is_solvent_atom(i)`, `role_of_atom(i)` — never check `i >= num_solute_atoms()`.
- **gather.rs:** `nearest_image(ri, rj)` returns `ri − rj` (displacement); gathered position = `rj + nearest_image(ri, rj)`.

## The instancing model

`topology.rs` implements this; the comments there assume it rather than restate it. Phases 1–2e are
built; the bonded loops over instances are the open part (item `2.0`).

**The problem it solves.** GROMOS's solute/solvent split fuses two orthogonal concerns into one
rigid boundary:

- *Representation* — store a molecule's parameters once, not once per copy. 216 waters should not
  be 216 copies of the water topology.
- *Role* — "this molecule is the thing I am studying" vs "this is bulk".

gromosXX fuses them as an index threshold over one array (`is_solvent(i)` is `i >= num_solute_atoms()`),
which is rigid but uniform. gromos++ fuses them harder, in *separate containers of different
classes* (`mol(i)` → `Molecule`, `sol(i)` → `Solvent`), which is why every analysis program has two
code paths and `AtomSpecifier` bifurcates into `m:` and `s:` namespaces.

**The separation.** A `MoleculeType` holds the parameters once; a `MoleculeInstance` is one
placement of it, carrying its own role. Then:

- Analysis iterates **one** atom array and filters by attribute — the dual-code-path problem
  disappears structurally rather than by convention.
- "Treat this one water as solute" is `promote()`: de-instance that molecule into an explicit,
  editable one while the other 215 stay instanced and compact. In gromosXX this means rebuilding the
  topology; here it is local. This is the operation the QM/MM zone work needs too — moving a water
  into the QM region *is* de-instancing.
- The solvent fast path survives: the moltype carries the "regular/rigid, use the vectorised 3-site
  kernel" hint, so gromosXX's `solvent_innerloop`/SPC-table speed is keyed on *moltype* rather than
  on the solute/solvent boundary.

**The cost, stated honestly.** Instancing adds one indirection (`atom → instance → moltype →
parameter`) per parameter lookup. The mitigation is to keep the genuinely hot per-atom parameters
(charge, IAC) materialised in a flat array *as a rebuildable cache* derived from the moltype — dedup
as the source of truth, flat array as the accelerator — while bonds, constraints and exclusions stay
in the moltype.
