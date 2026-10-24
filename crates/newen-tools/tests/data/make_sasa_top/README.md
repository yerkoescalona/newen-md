# make_sasa_top reference (gromos++)

`aladip_sasa.gromospp.top`: `make_sasa_top @topo shared/aladip.topo @sasaspec aladip.sasaspec`;
`aladip_sasa_noH.gromospp.top`: the same with `@noH`. `aladip.sasaspec` is a synthetic SASASPEC
block covering the six atom types of aladip. The whole output (topology in gromos++'s
`OutTopology` layout plus the SASAPARAMETERS block) is compared token by token.
