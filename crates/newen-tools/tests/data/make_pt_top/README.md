# make_pt_top reference (gromos++, .local/gromosPlsPls, 2026-08-30)

`shared/` is `crates/newen-md/tests/gromosXX_references/shared/`.

ch4_B.top is shared/ch4_spc.top with atom 1 changed to IAC 18, charge 0.1; ch4_to_B.gromospp.ptp: gromos++ `make_pt_top @topo shared/ch4_spc.top ch4_B.top @softpar 1.51 0.5`.
