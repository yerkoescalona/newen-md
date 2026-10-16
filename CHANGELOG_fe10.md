### Added (1.8 FE-10 — the analytic companions of `dGslv_pbsolv`)

- **`gromos_analysis::free_energy::charge_correction`**: the correction terms of Kastenholz &
  Hünenberger, J. Chem. Phys. 124, 224501 (2006) (doi:10.1063/1.2201698) for a charge change
  q_A → q_B on a spherical ion in a cubic box, each added to the raw ΔG: A1 (Born beyond the
  cut-off, Eq. 25), A2 (polarisation inside the cut-off: the paper's radial Fredholm equation,
  Eqs. 26–29 with kernels A3/A6, solved by product integration and Richardson-extrapolated),
  B (Eq. 32 for lattice sums, the fit of Eq. 33 for cut-offs), C1 (summation scheme, from the
  solvent's quadrupole-moment trace γ, Eqs. 34–42), C2 (Eq. 46, and the SPC fit of Eq. 47 for
  lattice sums). Schemes LS (tinfoil), CM, BM (= GROMOS RF with charge groups; `RF`), BA;
  EPSRF 0 is conducting, as in the engine. Quadratic terms scale with q_B² − q_A², linear ones
  with q_B − q_A. γ is Eq. 35 over the model's charges: SPC 0.0082, SPC/E 0.008476, TIP3P
  0.0076414, TIP4P 0.0092948 e nm². `ChargeCorrection::beyond_pb` = C1 + C2 is what to add on
  top of `dGslv_pbsolv`, whose correction already holds the polarisation terms (B for LS,
  A1 + A2 + B for RF); `total` is all five, for a sphere without the PB solves.
- Python: `gromos.analysis.charge_correction(delta_q, box_edge, cutoff, eps_rf, solvent="SPC",
  *, scheme="RF", charge_a=0.0, ion_radius=0.2, eps_solvent=78.0, eps_model=None,
  density=33.427, quadrupole_trace=None, four_pi_eps_i=138.9354)` → dict `a1, a2, b, c1, c2,
  total, beyond_pb, delta_phi_odl, quadrupole_trace` (kJ/mol), over the private binding
  `_charge_correction` (stub in `gromos.pyi`). Refused with `ValueError`: an unknown solvent or
  scheme (GM/GA are not covered), C2 for a lattice sum in anything but SPC, A2 without a model
  permittivity (tabulated for SPC only, 66.6), a cut-off not beyond the ion radius.
- Checked against the paper's Tables III and IV (Na⁺ in SPC, R_I = 0.2 nm, L = [(N_w + 1)/ρ]^⅓):
  A1, C1, and B, C1, C2 of all ten lattice-sum boxes to 0.01 kJ/mol; B of all 34 CM/BM/BA rows
  to 0.02; A2 for BM within 0.06. **Not reproduced to 0.1:** A2 for CM, where our converged
  solve lies 0.07–0.25 kJ/mol below Table III (24.23 vs 24.48 at R_C = 0.8 nm, 13.47 vs 13.54
  at 1.6), the gap shrinking with R_C; the same solver unconverged (spacing R_C/80) lands within
  0.12, so it is most likely the table's discretisation. That test runs at 0.26 kJ/mol and says
  so. The corrected ΔG°_slv (ΔG_chg + the five terms + 13.62) of all 44 LS/CM/BM/BA rows comes
  out within 0.06 kJ/mol of the tables, CM within 0.26 (its A2); six rows are pinned in the test.
  `py-gromos/tests/test_charge_correction.py`: a row of each scheme, the charge-change scaling,
  other water models, the refusals.
- User guide: "Correcting a charge change" in `docs/user-guide/free-energy.md` — the five
  terms, which of them `dGslv_pbsolv` provides, and the chain raw ΔG (TI/BAR/MBAR) +
  `dGslv_pbsolv` + `charge_correction(...)["beyond_pb"]`.
- Paper oddities found on the way (not ported): Table III prints L = 3.13 nm for the LS row with
  N_w = 2048 (the box is 3.94 nm; its B and C1 are 3.94's), Table IV prints A1 = −557.15 for
  −57.15 (BA, N_w = 1024, R_C = 1.2), and Table V lists C1 = −79.76 (= −ΔΦ_ODL) for every GM and
  GA row although Eqs. 43–44 depend on R_C and R_I (about −72.7 kJ/mol at R_C = 0.8 nm as we
  read Eq. 43).
