# explode reference (gromos++)

`meoh2.top`: gromos++ `com_top @topo shared/meoh_spc.top shared/meoh_spc.top @param 1 @solv 1`; `meoh2.cnf`: the methanol solute twice, 1 nm apart, two waters.
`explode_2_1.5.gromospp.cnf`: `explode @topo meoh2.top @pos meoh2.cnf @nsm 2 @dist 1.5`.
