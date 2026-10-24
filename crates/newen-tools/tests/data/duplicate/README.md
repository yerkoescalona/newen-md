# duplicate reference (gromos++)

`meoh2_dup.cnf`: `explode/meoh2.top` coordinates with molecule 2 an exact copy of molecule 1 and a duplicated water.
`report.gromospp.out`: `duplicate @topo ../explode/meoh2.top @pos meoh2_dup.cnf @pbc r`; `write.gromospp.cnf`: the same with `@write`.
