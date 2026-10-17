# pdb2g96 fixture

`aladip.pdb`: the solute of `newen-md/tests/gromosXX_references/shared/aladip.conf` written as PDB
(nm → Å, PDB's fixed columns). Generated from that configuration, not from gromos++ — there is no
gromos++ reference for `pdb2g96` yet, so this fixture only serves `tests/library_parity.rs`: it
gives the program real input on which the library `run` and the binary must produce the same bytes.

A gromos++ reference would be

    pdb2g96 @topo ../../../newen-md/tests/gromosXX_references/shared/aladip.topo @pdb aladip.pdb
