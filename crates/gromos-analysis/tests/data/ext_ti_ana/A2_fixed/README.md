# Case A2 through the patched gromos++ ext_ti_ana (9cd9f6b3 + D5), item 1.8

Command: case A2 of `../README.md`, as `the regeneration scripts kept outside the repository PROGRAMS A2_fixed` runs it (the
same script with the unpatched program reproduces `../A2/` byte for byte). Binary
<gromos++ 9cd9f6b3 + patches, BUILD/programs/>ext_ti_ana.

Only `predicted_l0.125.dat` changes: with `@lamdih 0 0 1 0 0` (Lambda_dih = lambda^2) the dihedral term is now
multiplied by dLambda/dlambda = 2 lambda (D5). The point at lambda = 0.5 is unchanged (factor 1).
`bar_data_l0.125.dat` and stdout are byte-identical to `../A2/`. The predicted values equal
~/Development/gromos-divergences/data/e3_dih.csv `pred_corrected` (variant a2_quad) to printed precision.
