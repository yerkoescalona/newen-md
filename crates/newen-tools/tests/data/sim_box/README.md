# sim_box reference (gromos++)

`meoh.cnf`: united-atom methanol (3 atoms, the solute of `newen-md/tests/gromosXX_references/shared/meoh_spc.cnf`).
`meoh_spc_minwall1.gromospp.cnf`: what gromos++ `sim_box` (`.local/gromosPlsPls`, 2026-08-30) writes for

    sim_box @topo ../../../newen-md/tests/gromosXX_references/shared/meoh_spc.top @pbc r \
            @pos meoh.cnf \
            @solvent ../../../newen-md/tests/gromosXX_references/water_1000_spc_gridcell/water_1000_spc.cnf \
            @minwall 1.0

(gather warnings stripped): 352 SPC molecules in a 2.207681950 nm cubic box. `tests/sim_box_reference.rs`
requires our `sim_box` to produce the same box, the same count and the same positions.
