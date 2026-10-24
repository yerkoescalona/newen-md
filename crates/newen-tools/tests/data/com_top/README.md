# com_top reference (gromos++)

`aladip_ch4.gromospp.top` is what gromos++ `com_top` (`.local/gromosPlsPls`, 2026-08-30) writes for

    com_top @topo ../../../newen-md/tests/gromosXX_references/shared/aladip.topo \
                  ../../../newen-md/tests/gromosXX_references/shared/ch4_spc.top @param 1 @solv 1

`tests/com_top_reference.rs` requires our `com_top` to produce a topology that reads back identical.
