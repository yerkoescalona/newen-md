# make_top reference (gromos++)

Force-field files from the LiveCoMS GROMOS tutorial (`t_06/topo`: GROMOS 54A7 `54a7.mtb`, `54a7.ifp`,
and `MeOH_exp_H.mtb`, a methanol building block with explicit hydrogens). The `*.gromospp.top` files
are what gromos++ `make_top` (`.local/gromosPlsPls`, 2026-08-30) writes for

    make_top @build 54a7.mtb @param 54a7.ifp @seq NH3+ ALA GLY COO- @solv H2O   → nh3_ala_gly_coo.gromospp.top
    make_top @build 54a7.mtb MeOH_exp_H.mtb @param 54a7.ifp @seq MeOH @solv H2O → meoh_exp_h.gromospp.top

`tests/make_top_reference.rs` requires our `make_top` to produce topologies that read back identical
(atoms, exclusions, 1-4 pairs, every bonded term, every parameter block, the solvent).
