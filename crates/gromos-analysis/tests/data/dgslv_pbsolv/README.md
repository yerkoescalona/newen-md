# dGslv_pbsolv reference outputs (gromos++, 2026-10-11)

The gromos++ build ships only the "no FFTW" stub of `dGslv_pbsolv`: the `#ifdef HAVE_LIBFFTW3`
at `programs/dGslv_pbsolv.cc:108` comes before `#include "../config.h"`, and line 109 is a
deliberate syntax error. The `pb` library itself is built with FFTW. The reference binary was
built from a patched copy of the program linked against the built library:

    PP=.local/gromosPlsPls/gromos++
    sed -e '108s/.*/#if 1/' -e '109d' $PP/programs/dGslv_pbsolv.cc > dGslv_pbsolv.cc
    g++ -std=c++11 -O3 -funroll-loops -fno-gcse -DNDEBUG -DHAVE_CONFIG_H \
      -iquote $PP/programs -iquote $PP/BUILD/programs -I$PP/src \
      dGslv_pbsolv.cc $PP/BUILD/src/.libs/libgromos.a -lfftw3 -lgsl -lgslcblas -lz -lm -o dGslv_pbsolv

`shared/` is `crates/gromos-md/tests/gromosXX_references/shared/`. `<case>.stdout` is the
program's stdout, `<case>.log` the `@verbose` log.

| case | gromos++ command |
|------|------------------|
| c1 | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @pbc r @atoms 1:a @atomsTOcharge 1:a @schemeELEC LS @epsSOLV 66.6 @gridspacing 0.1 @probeIAC 4 @verbose c1.log` |
| sigma | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @atoms 1:a @atomsTOcharge 1:4,5 @schemeELEC LS @epsSOLV 66.6 @gridspacing 0.1 @NPBCsize 5.0 @increasegrid 2 2 2 @rminORsigma 1 @radscal 1.1 @probeIAC 4 @verbose sigma.log` |
| order | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @pbc r @atoms 1:5,3,1 1:7 @atomsTOcharge 1:9,2 @schemeELEC ls @epsSOLV 66.6 @epsNPBC 66.6 @gridspacing 0.13 @probeIAC 5 @probeRAD 0.12 @maxiter 40 @verbose order.log` |
| pqr | `dGslv_pbsolv @pqr aladip.pqr @coordinates 3.767055681 3.767055681 3.767055681 @nogridpoints 24 @gridspacing 0.15 @atoms 1-12 @atomsTOcharge 2,3,5 @schemeELEC LS @epsSOLV 66.6 @radscal 1.1 @verbose pqr.log` |

What each case pins:

- **c1**: the four finite-difference solves with default settings, the gather warning on
  stdout and the COG gather, the hydrogen-radius branch, the PBC neutralisation and average
  removal, the NPBC Debye–Hückel boundary. Grids: PBC 38³, NPBC 78³.
- **sigma**: σ radii and `@radscal`, `@NPBCsize` (NPBC 50³) and `@increasegrid` on the PBC grid
  only (40³), a read-out set different from the charged set, no `@pbc` (no gather, no warning).
- **order**: `@atoms` in the order written (5, 3, 1, 7: the charge-spreading order and the read-out
  rows follow it), a lower-case `@schemeELEC`, another probe type and probe radius, and
  `@maxiter 40`, which three of the four solves reach (the `NO CONVERGENCE` log lines).
- **pqr**: PQR input, a synthetic file made from `aladip.conf` (positions ×10, in Å), the topology
  charges and hand-picked radii in Å (0 for both hydrogens), with 10- and 11-field ATOM lines and
  HETATM lines; specifiers without a molecule, the spacing taken from `@nogridpoints` (gromos++
  ignores the value of `@gridspacing` except for the PBC grid size), the solvent-first solve order.
  gromos++ exits 1 after printing the table and writes a `four_pi_eps_i` warning to stderr (not
  stored); ours exits 0 with the same warning.

`tests/gromospp_references.rs` runs ours with the same arguments and requires stdout and the log
to be byte-identical.
