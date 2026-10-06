# Pairlist & Neighbor Search

The pairlist (neighbor list) is the data structure that answers the question: *which pairs of
atoms are close enough to interact?* Because nonbonded interactions have a finite cutoff radius
`r_c`, only pairs with `|r_ij| ≤ r_c` contribute to the force. The pairlist pre-computes those
pairs so the inner force loop does not have to test every possible pair on every MD step.

## Why GROMOS uses charge groups, not atoms

Most MD engines (GROMACS, AMBER) use atoms as the cutoff unit. GROMOS uses **charge groups**.

A charge group is a set of covalently bonded atoms whose partial charges sum to zero (or to a
small integer). Water (SPC/E) is one charge group of three atoms. A carbonyl C=O is one charge
group of two atoms.

The reason is the **reaction-field (RF) self-term**. When an atom pair is cut off by the
spherical cutoff, a reaction-field correction is applied to approximate the electrostatic
contribution from everything beyond the cutoff. This correction is derived assuming a uniform
dielectric *outside* the cutoff sphere. If you cut off *atoms* individually, you can split a
charge group across the boundary — some atoms are inside, some outside — leaving a net charge
imbalance at the cutoff surface that the RF formula was not designed to handle. The resulting
energy drift is systematic and grows with system size.

By cutting off *charge groups* as a unit, GROMOS ensures that the total charge inside the cutoff
sphere changes only in neutral increments. The RF self-term is then exact to the order of the
approximation.

> **Implementation note.** The cutoff test is performed on the **center of geometry (COG)** of
> each charge group. If two CG COGs are within the cutoff, *all* atom pairs between those two
> CGs are added to the pairlist — regardless of the individual atom-atom distances. This is the
> GROMOS standard and is reproduced exactly in gromos-rs.

## Twin-range cutoffs

GROMOS supports two cutoff radii:

- **Short-range cutoff** `r_s` (IMD `RCUTP`): pairs within this distance are evaluated every
  step.
- **Long-range cutoff** `r_l` (IMD `RCUTL`): pairs between `r_s` and `r_l` are evaluated only
  on pairlist-update steps (every NSNB steps) and their contribution is reused ("twin-range"
  or "charge-group shifted") between updates.

When `r_s == r_l` (the common case), there is no twin-range and the `solute_long` /
`solvent_long` lists are always empty.

`PairlistContainer` separates pairs into four lists:

| List | Content |
|------|---------|
| `solute_short` | Solute–solute pairs within `r_s` (exclusions removed) |
| `solute_long` | Solute–solute pairs in `(r_s, r_l]` |
| `solvent_short` | Solvent–solvent and solute–solvent pairs within `r_s` |
| `solvent_long` | Solvent–solvent and solute–solvent pairs in `(r_s, r_l]` |

The split matters for the inner loop: solvent molecules are rigid and identical, so the
`solvent_innerloop` kernel can use a specialized fast path that expands the first-atom sentinel
pair to all atom pairs within the molecule type.

## The skin distance

The pairlist is not rebuilt every step (that would dominate runtime). Instead it is rebuilt
every `NSNB` steps (IMD `NSNB`, typically 5). To avoid missing pairs that drift into the cutoff
between rebuilds, the pairlist is built with cutoff `r_c + skin` where `skin` is a small extra
distance (typically 0.1 nm). Any pair that enters the true cutoff before the next rebuild was
already in the list.

The skin introduces a trade-off: larger skin → safer (fewer missed pairs) but larger pairlist →
more force evaluations. The displacement-triggered rebuild criterion resolves this
more precisely: rebuild when `2 · max_displacement > skin`.

The factor of **2** is load-bearing. Two atoms can each move `max_displacement` toward each
other, closing their separation by `2 · max_displacement`. A rebuild triggered on factor-1
would miss pairs under convergent motion.

## Algorithm complexity

### StandardPairlistAlgorithm — O(N²)

Tests every CG pair. Scales quadratically with system size. Correct for all box types
(rectangular, triclinic, vacuum) and all system sizes. Kept as the always-correct reference
oracle in tests.

```
for each CG pair (i, j) with i < j:
    compute r_ij = nearest_image(COG_i, COG_j)
    if |r_ij| ≤ r_l: add atom pairs to appropriate list
```

### CellListPairlistAlgorithm — O(N)

Bins CG COGs into a spatial grid of cells sized ≥ `r_l + skin`. Each CG only tests pairs
against CGs in its own cell and the 26 neighboring cells. For a uniform distribution of N
atoms, the number of pairs per cell is constant → O(N) total.

```
build grid, bin each CG COG into cell(cg)
for each cell_a:
    for each neighbor cell_b of cell_a (including itself):
        for each CG pair (cg1 ∈ cell_a, cg2 ∈ cell_b) with cg1 < cg2:
            compute r and classify exactly as StandardPairlistAlgorithm
```

The CG pair classification logic (`process_cg_pair`) is identical in both algorithms —
CellList only changes *which* pairs are tested, not how they are classified. This is the
invariant that the set-equality tests verify.

Falls back to `StandardPairlistAlgorithm` for triclinic and vacuum boxes (axis-aligned cells
are not periodicity-safe without extra stencil machinery).

## Validation strategy

Every change to the pairlist must pass two independent checks:

1. **Set-equality test** (unit test, no MD): the pair set produced by `CellListPairlistAlgorithm`
   is identical to the set produced by `StandardPairlistAlgorithm` after normalizing order. This
   runs on every reference topology and is the ground-truth correctness check for the algorithm.

2. **Energy regression test** (integration test, full MD): energies from the reference test suite
   match GROMOS to within `ENERGY_REL_TOL = 1e-8` on every step for all 37 systems. A change to
   the pairlist that passes set-equality but changes *iteration order* will change floating-point
   summation order and may cause marginal failures here — the acceptable margin must be measured
   empirically before lowering the auto-select threshold below the largest reference system.

## Cell-list roadmap

The cell-list cards (`9b`–`9e`) on the board carry the implementation sequence. The sub-dimensions in order:

- **9a-0**: Add `PairlistAlgorithm` enum; wire dispatch; all reference systems keep Standard → zero float change.
- **9a-1**: Measure Standard vs CellList energy margin over 100 steps; lower threshold if safe.
- **9d**: `ChargeGroupTable` — COG cache, cg-of-atom lookup, energy-neutral (COGs are bit-identical).
- **9e**: `SpatialIndex` service — general neighbor query API (prerequisite for QM/MM).
- **9c**: Parallel cell build + displacement-triggered rebuild.
- **9b**: Spatial reordering of atom arrays for cache coherence (last, most invasive).
