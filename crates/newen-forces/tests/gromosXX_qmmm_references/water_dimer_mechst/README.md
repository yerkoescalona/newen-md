# Water dimer, mechanical embedding with constant QM charges (`NTQMMM=-1`)

Real gromosXX QM/MM input and output data, extracted from the published dataset accompanying:

> Poliak, P. et al. "A Robust and Versatile QM/MM Interface for Molecular Dynamics in GROMOS."
> *J. Comput. Chem.* 2025, 46, e70053. https://doi.org/10.1002/jcc.70053

Data repository: https://doi.org/10.5281/zenodo.14549978 (CC BY 4.0). Files here are a small
excerpt of `data-gromos_qmmm.zip`'s `WATER_DIMER/` directory, kept under this license with
attribution.

## Files

- `water_dimer.top` — `WATER_DIMER/topo/water.top`: one solute water (the would-be QM zone) +
  one `NSM=1` solvent water (the would-be MM zone), standard SPC/54A7 parameters.
- `water_dimer.cnf` — `WATER_DIMER/eq/eq_2waters_am1_mechst_5.cnf`: the real equilibrated starting
  configuration fed into `md_2waters_am1_mechst_1.imd` (see that run's own `@conf` argument in the
  archived `.run` script).
- `water_dimer.imd` — `WATER_DIMER/md_NVE/md_2waters_am1_mechst_1.imd`: the real input parameters
  (`PAIRLIST`/`NONBONDED`/`QMMM` blocks) for that run.
- `expected/am1_mechst_nonb.out` — `WATER_DIMER/md_NVE/ana/ene_ana/am1_mechst_nonb.out`: the real
  `ene_ana`-derived nonbonded-energy time series gromosXX produced for this run.

## Why this specific run is usable without a QM engine

`QMMM` block: `NTQMMM=-1` — mechanical embedding with **constant** QM charges. Per gromosXX's own
`QMMM_Interaction::modify_exclusions` (`interaction/qmmm/qmmm_interaction.cc`), mechanical
embedding leaves every QM-MM pair's exclusions untouched — the classical nonbonded loop runs
exactly as if there were no QM/MM at all, using the topology's ordinary charges. So the `nonb`
column is the plain classical LJ+CRF total for this system, reproducible with no QM engine, no
per-step charge derivation, and no QM-zone atom selection needed.

The `am1`/`om2`/`om3`/`pm3`/`elst`/`mech` (non-`mechst`) variants in the full archive all need a
semi-empirical QM program (MNDO or similar; not installed, not readily installable here) to
reproduce their QM-dependent energy terms — only `mechst` was extracted for that reason.

## Raw `.tre`/`.trc` trajectories were not archived

The Zenodo deposit keeps only inputs (`.cnf`/`.imd`/`.run`) and `ene_ana`-post-processed energy
time series (`.out`), not the raw gromosXX `.tre` output — so the comparison in
`water_dimer_qmmm_mechst_reference.rs` is against `am1_mechst_nonb.out`'s `t=0` row, not a raw
`ENERGY03` block the way the rest of this repo's reference suite compares.
