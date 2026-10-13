# Chain fixture through the patched gromos++ (9cd9f6b3 + D4 + D5 + D8), item 1.8

Binaries: <gromos++ 9cd9f6b3 + patches, BUILD/programs/>{ext_ti_ana,bar,ext_ti_merge}
(worktree of gromos++ on branch `divergences-fix`; patches `crates/gromos-md/tests/gromosXX_references/patches/D4_bar_error.patch`,
`D5_ext_ti_ana_dlam.patch`, `D8_ext_ti_merge_loop.patch`; build recipe `crates/gromos-md/tests/gromosXX_references/patches/BUILD.md`).
Commands: exactly those of `../README.md` (chain), as `the regeneration scripts kept outside the repository PROGRAMS chain_fixed`
run them; the same script with the unpatched programs reproduces `../chain/` byte for byte.

Against `../chain/` only two files change:
- `bar.stdout`: the four window errors and the total (D4: no second division by d^2):
  2.11698, 2.40196, 1.49959, 0.44748, total 3.56370 (was 2.56637, 2.41098, 1.23029, 0.54694, 3.76986).
  DG, the overlap integrals and `dist_*.dat` are unchanged.
- `merge.stdout`: the error of the integral (D8: the loop now includes point N-2): 2.302129 (was 2.248844);
  the integral -51.651947 and the 21 points are unchanged. The "error" column is still the @countframes
  frame count read as an error (gromos++ quirk), so only the change is meaningful.
The ext_ti_ana outputs are unchanged: this chain uses linear lambda dependences, where D5 multiplies by
dlam = 1.

Values agree with the corrected values of ~/Development/gromos-divergences/data/e3_bar.csv (err_undivided)
and e3_merge.csv (err_full_python) to printed precision (`the regeneration scripts kept outside the repository`).
