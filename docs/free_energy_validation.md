# Free-energy validation: the hydration free energy of Na⁺

**Purpose.** This page checks the whole charge-change chain of newen on one system: the λ
windows, the TI, BAR and MBAR estimators, the reaction-field self term, `dGslv_pbsolv` and the
Kastenholz–Hünenberger correction terms. The result is compared with the methodology-independent
hydration free energy of Na⁺ that Kastenholz and Hünenberger derived. This is a **check of the
chain at the level of a few kJ/mol**. It is not a converged number: under a 15-minute wall-time
budget each window ran for 20 ps.

Raw outputs: `~/.cache/gromos-rs-item-1.8/fe11/` (`README.md` lists every file and command). Every
number below is in `box3.0/analysis.json`, `box3.0/pb.json` or `box3.0/correction.json`, and
`box3.0/tables.md` renders them.

## System

- One Na⁺ (atom 1 of `crates/newen-md/tests/gromosXX_references/nacl_water_box/nacl_water_box.topo`,
  IAC 34, +1 e) in 878 SPC. The box is cubic, 3.0 nm at the start (2635 atoms; `sim_box` from the
  LiveCoMS t_01 SPC box) and 2.9984 ± 0.0022 nm on average under NPT (mean ± std over the 33
  window means).
- The perturbation is `shared/nacl_na_dummy.ptp` (the FE-9 reference): Na⁺ (IAC 34, +1) becomes a
  dummy (IAC 19, 0). The mass is kept and α_LJ = α_CRF = 0.5. Charge and Lennard-Jones
  interactions are switched off together along one λ.
- The Na⁺–OW parameters are C12 = 1.886450e-7 kJ mol⁻¹ nm¹² and C6 = 4.190170e-4 kJ mol⁻¹ nm⁶.
  Kastenholz and Hünenberger used 2.303628e-7 and 4.298760e-4 (their Sec. III). **The ion is
  therefore not the same model.**

## Protocol

| stage | settings |
|---|---|
| minimisation | `na_em.imd` |
| equilibration | `na_eq.imd`: state A, 100 ps NPT, 298.15 K / 1 atm (Berendsen 0.1 / 0.5 ps), from an earlier session (`prepare.out`, 286 s) |
| windows | `na_window.imd`: λ = 0, 0.1, …, 1 (11 windows), each started from `eq.cnf`; NPT as above; 2 fs; SHAKE; reaction field RCRF 1.4 nm, EPSRF 61, EPSCS 1, κ 0; pair list every 5 steps (NSNB 5), RCUTP = RCUTL = 1.4; NTWG = NTWE = 50 (a multiple of NSNB); `PRECALCLAM 11 0 1` |
| length | **10 000 steps = 20 ps per window**, n = 3 replicas (seed 210185 + 1000·r + i for replica r, window i) |
| analysis | the first 4 ps (40 of 200 frames, 20 %) of each window dropped, so 160 frames are kept; TI by trapezoid, BAR and MBAR over the `PRECALCLAM` rows (`newen.free_energy.collect`) |
| PB | `newen.analysis.dgslv_pbsolv` on `eq.cnf` for both end states: `@pbc r`, `@atoms 1:a`, `@atomsTOcharge 1:a`, `@epsSOLV 66.6`, `@probeIAC 5`, `@gridspacing 0.1`, `@NPBCsize 5.0`, RF with `@epsRF 61 @rcut 1.4`, and LS for comparison |

The window length was set from a timing probe: 33 windows × 1500 steps, all at once, took 96 s
(`probe33/`), about 64 ms per step per window with 33 single-threaded processes on 8 cores /
16 threads. Production was 33 windows × 10 000 steps in 672 s. FE-11 used 810 s of wall time in
all, counting the probe and a 42-s start that was abandoned.

Committed inputs: `docs/free_energy_validation/na_em.imd`, `na_eq.imd`, `na_window.imd` and
`scripts/na_hydration.py`. Run them as follows:

```bash
cd py-newen
W=~/.cache/gromos-rs-item-1.8/fe11/box3.0
uv run python ../scripts/na_hydration.py build    --workdir $W --box 3.0
uv run python ../scripts/na_hydration.py prepare  --workdir $W
RAYON_NUM_THREADS=1 uv run python ../scripts/na_hydration.py windows --workdir $W \
    --replicas 0,1,2 --steps 10000 --workers 33
uv run python ../scripts/na_hydration.py analyse  --workdir $W --skip-ps 4
uv run python ../scripts/na_hydration.py pb       --workdir $W
uv run python ../scripts/na_hydration.py correct  --workdir $W
```

## Raw ΔG(Na⁺ → dummy)

| estimator | n | mean ± std (kJ/mol) | replicas | mean single-run error |
|---|--:|--:|---|--:|
| TI | 3 | 384.24 ± 1.92 | 383.72, 382.63, 386.36 | 1.27 |
| BAR | 3 | 388.31 ± 2.39 | 387.89, 386.15, 390.88 | 1.95 |
| MBAR | 3 | 388.13 ± 2.27 | 387.88, 386.00, 390.51 | 1.35 |

The std is the sample std (n − 1) over replicas. The single-run error is the error each estimator
reports, averaged over the replicas.

BAR and MBAR agree to 0.2 kJ/mol. TI is 4.1 kJ/mol lower in every replica (3.5 to 4.5). This is
the quadrature error of the trapezoid on 11 points: ⟨∂H/∂λ⟩ rises from 743 to 1059 kJ/mol
between λ = 0 and 0.1 and then falls steeply (table below), so TI needs a finer grid at the
charged end. **BAR and MBAR are the reference values here.**

### ⟨∂H/∂λ⟩ per window

| λ | mean ± std (n = 3), kJ/mol | mean block error | frames kept per replica |
|---|--:|--:|--:|
| 0.0 | 743.25 ± 5.79 | 3.09 | 160 |
| 0.1 | 1058.70 ± 6.06 | 6.39 | 160 |
| 0.2 | 881.56 ± 6.11 | 6.79 | 160 |
| 0.3 | 563.75 ± 10.51 | 5.30 | 160 |
| 0.4 | 343.31 ± 6.02 | 4.13 | 160 |
| 0.5 | 214.69 ± 6.20 | 3.04 | 160 |
| 0.6 | 140.14 ± 2.40 | 2.49 | 160 |
| 0.7 | 98.43 ± 5.14 | 2.56 | 160 |
| 0.8 | 71.49 ± 7.66 | 2.01 | 160 |
| 0.9 | 65.28 ± 3.23 | 1.12 | 160 |
| 1.0 | 66.77 ± 1.08 | 0.89 | 160 |

## dGslv_pbsolv on the end states

| state | NPBC | PBC | FFT LS | FFT RF | correction RF | correction LS |
|---|--:|--:|--:|--:|--:|--:|
| A (Na⁺) | −416.63 | −351.27 | −337.65 | −325.93 | −77.08 | −65.36 |
| B (dummy) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

All values are in kJ/mol. The RF correction is CB/NPBC − LS/PBC + LS/PBC[FFT] − RF/PBC[FFT], which
gives PB(B) − PB(A) = +77.08 kJ/mol for the change 1 → 0.

## The chain, term by term

The charge goes from q_A = +1 to q_B = 0 (kJ/mol; `correction.json`; L = 2.9984 nm, R_C = 1.4 nm,
EPSRF 61, R_I = 0.2 nm, SPC):

| term | value | source |
|---|--:|---|
| ΔG_self (RF self term held by the raw ΔG) | +73.825 | `rf_self_energy(1, 0, 1.4, 61)` = −½·FPEPSI·(1 − C_rf/2)/R_C·(0 − 1) |
| PB(B) − PB(A) (A1 + A2 + B for the real solute) | +77.083 | `dgslv_pbsolv`, RF |
| A1 / A2 / B, analytic sphere | +48.984 / +24.435 / +5.386 | `charge_correction(-1, 2.9984, 1.4, 61, charge_a=1)` |
| C1 | +77.587 | same |
| C2 | 0 | same (cut-off scheme) |
| ΔG_cmp (standard state, 298 K) | +7.95 | K&H Eq. 13 (`literature/k2.txt`, line 383) |

The analytic A1 + A2 + B for the sphere is 78.805, against 77.083 from PB on the real solute.

The corrected change and the hydration free energy:

ΔG_cor(1 → 0) = [ΔG_raw − ΔG_self] + PB(B) − PB(A) + C1 + C2
             = [ΔG_raw − 73.825] + 77.083 + 77.587 + 0

ΔG°_slv(Na⁺) = −ΔG_cor(1 → 0) + ΔG_cmp

The raw ΔG removes the Lennard-Jones interactions as well as the charge, so it already contains
the cavity term that Kastenholz and Hünenberger add separately (their ΔG_cav = 5.67 kJ/mol for
their ion). Only ΔG_cmp is added here.

| estimator | ΔG_cor(1 → 0), PB route | ΔG°_slv(Na⁺), PB route | ΔG°_slv(Na⁺), analytic sphere | mean single-run error |
|---|--:|--:|--:|--:|
| TI | 465.08 ± 1.92 | −457.13 ± 1.92 | −458.85 ± 1.92 | 1.27 |
| BAR | 469.15 ± 2.39 | −461.20 ± 2.39 | −462.92 ± 2.39 | 1.95 |
| MBAR | 468.97 ± 2.27 | −461.02 ± 2.27 | −462.75 ± 2.27 | 1.35 |

All values are mean ± std over n = 3, in kJ/mol. Every term except the raw ΔG is a constant, so
the std is the std of the raw ΔG. The analytic route is [raw − ΔG_self] + A1 + A2 + B + C1 + C2
+ ΔG_cmp, the same formula with the sphere's A1 + A2 + B in place of PB.

## Comparison with Kastenholz and Hünenberger

The target is the methodology-independent ΔG°_slv(Na⁺) = **−466.5 ± 1.4 kJ/mol**. This is the
mean of the corrected values of their simulations with N_w ≥ 512 over every boundary condition
and scheme. The std of those simulations is 11.6 kJ/mol, and ±1.4 is the error of the mean
(J. Chem. Phys. 124, 224501 (2006); `literature/k2.txt`, lines 1224–1234).

- **BAR/MBAR: −461.2 ± 2.4 and −461.0 ± 2.3 kJ/mol.** These lie 5.3 and 5.5 kJ/mol above the
  target. The standard error of the difference is √(2.4²/3 + 1.4²) ≈ 2.0 kJ/mol, so the gap is
  about 2.7 standard errors. Counting only our error bars, this is a disagreement at the few-kJ/mol
  level. It is half the 11.6 kJ/mol spread of the paper's own corrected values across methods.
  It is within the precision this run was meant to check.
- The gap has three known contributions that this run does not resolve:
  - The Na⁺ Lennard-Jones parameters differ from the paper's (C12 18 % and C6 3 % smaller).
  - Every window is 20 ps and starts from the state-A configuration, so each gets only 4 ps to
    relax before sampling.
  - The PB and analytic routes differ by 1.7 kJ/mol (A1 + A2 + B for the sphere against the real
    solute).
- TI (−457.1 ± 1.9) is a further 4 kJ/mol above BAR/MBAR because of its quadrature error, not
  because of the chain.
- **The self term is not optional.** Without the subtraction, every corrected value moves by
  −73.8 kJ/mol (BAR: −535.0), 69 kJ/mol beyond the target. With it, the chain lands within a few
  kJ/mol. This is the check this run was designed for. Subtracting the self term and also adding
  the analytic A1 on top of PB would count the polarisation twice (49 kJ/mol).

**Neutral control.** CH₄ → dummy in 999 SPC (BENCHMARKING.md §7, 11 × 20 ps, n = 3) gives a raw
ΔG of −8.37 ± 2.11 (BAR) and −8.43 ± 2.14 (MBAR) kJ/mol, against about −8 kJ/mol from
experiment. With no charge change, ΔG_self, PB and C1 are all zero, so that number needs no
correction. The engine and estimators are thus right at the 2-kJ/mol level for a neutral solute,
and the Na⁺ result above tests what the charge adds.

## What a converged study would add

- Windows of 100 ps or more, with the first 20 ps dropped.
- A denser λ grid near the charged end, or separate charge and LJ stages.
- The paper's Na⁺ parameters.
- Several box sizes and cut-offs to show that the corrected value no longer depends on them.
