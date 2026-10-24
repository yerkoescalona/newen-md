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

`shared/` is `crates/newen-md/tests/gromosXX_references/shared/`. `<case>.stdout` is the
program's stdout, `<case>.log` the `@verbose` log.

| case | gromos++ command |
|------|------------------|
| c1 | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @pbc r @atoms 1:a @atomsTOcharge 1:a @schemeELEC LS @epsSOLV 66.6 @gridspacing 0.1 @probeIAC 4 @verbose c1.log` |
| sigma | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @atoms 1:a @atomsTOcharge 1:4,5 @schemeELEC LS @epsSOLV 66.6 @gridspacing 0.1 @NPBCsize 5.0 @increasegrid 2 2 2 @rminORsigma 1 @radscal 1.1 @probeIAC 4 @verbose sigma.log` |
| order | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @pbc r @atoms 1:5,3,1 1:7 @atomsTOcharge 1:9,2 @schemeELEC ls @epsSOLV 66.6 @epsNPBC 66.6 @gridspacing 0.13 @probeIAC 5 @probeRAD 0.12 @maxiter 40 @verbose order.log` |
| pqr | `dGslv_pbsolv @pqr aladip.pqr @coordinates 3.767055681 3.767055681 3.767055681 @nogridpoints 24 @gridspacing 0.15 @atoms 1-12 @atomsTOcharge 2,3,5 @schemeELEC LS @epsSOLV 66.6 @radscal 1.1 @verbose pqr.log` |
| c2 | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @atoms 1:a @atomsTOcharge 1:4,5 @schemeELEC RF @epsSOLV 66.6 @epsRF 66.6 @rcut 1.4 @cubesFFT 4 @gridspacing 0.1 @NPBCsize 5.0 @increasegrid 2 2 2 @rminORsigma 1 @radscal 1.1 @probeIAC 4 @verbose c2.log` |
| c3b | `dGslv_pbsolv @topo shared/aladip.topo @coord shared/aladip.conf @pbc r cog @atoms 1:1-4 @atomsTOcharge 1:2-4 @schemeELEC RF @epsSOLV 78.4 @epsRF 0 @rcut 0.9 @gridspacing 0.12 @probeIAC 4 @verbose c3b.log` (gromos++ wrote it as `c3.log`) |

The 1.8 step-4b gate cases (2026-10-16, the same binary, `~/.cache/gromos-rs-item-1.8/pb_oracle/`),
each with `@schemeELEC LS` and with `@schemeELEC RF @epsRF 78.5 @rcut 1.4` (`ch4_water_fep`'s
EPSRF and RCRF), all with

    @coord shared/ch4_spc.cnf @pbc r @atoms 1:a @atomsTOcharge 1:a @epsSOLV 66.6
    @gridspacing 0.1 @NPBCsize 5.0 @probeIAC 5

| case | `@topo` | system |
|------|---------|--------|
| na_ls, na_rf | `na.top` | one Na⁺ (`t_02/topo/Na.top` of the LiveCoMS tutorial, 54A8, charge +1) at the CH4's place in `ch4_spc.cnf`'s cubic box of 999 SPC (3.106 nm) |
| ch4q_ls, ch4q_rf | `ch4q.top` | `ch4_water_fep`'s CH4 in state B of `ch4_plus.ptp`: charge +1, everything else state A. `ch4_plus.ptp` is `newen.free_energy.Perturbation.mutate(Topology("shared/ch4_spc.top"), "1:1", charge_b=1.0, alpha_lj=0.0, alpha_crf=0.0)`, `ch4q.top` is `pt_top @topo shared/ch4_spc.top @pttopo ch4_plus.ptp @type TOPO` (newen; gromos++ reads it) |

Grids 32³ (PBC) and 50³ (NPBC), probe IAC 5 = OW. The gromos++ free energies (kJ/mol; the
`# DGRESULT` lines and, for RF, the value that ends each FFT section of the log), and the
correction they make (CB/NPBC − LS/PBC, + LS/PBC[FFT] − RF/PBC[FFT] for RF):

| case | NPBC | PBC | FFT LS | FFT RF | correction LS | correction RF |
|------|-----:|----:|-------:|-------:|--------------:|--------------:|
| na | −296.866 | −235.993 | −168.810 | −151.944 | −60.873 | −77.739 |
| ch4q | −296.549 | −235.623 | −206.065 | −189.005 | −60.926 | −77.986 |

Ours: stdout and every `DGRESULT` byte-identical; the RF logs differ only in the sign of a few
printed zeros of the FFT section (`-0.000` / `0.000`), which the test allows.
`py-newen/tests/test_free_energy_charge.py` takes the `ch4q` numbers as the oracle of the
correction it applies through `newen.analysis.dgslv_pbsolv` (0.01 kJ/mol).

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

- **c2**: `@schemeELEC RF` on the `sigma` settings: after the four FD solves, the FFT solve with
  the lattice-sum boundary (`@atoms` centred on the 40³ grid, hat charges in both alias lists:
  six heavy charged atoms with α1, the two hydrogens with α2) and the one with the reaction-field
  boundary (εRF = 66.6 in the vacuum field, εs in the tensors, `@rcut 1.4`, `@cubesFFT 4`), the
  six-column table and the log's switch to `fixed`/3 decimals at the first `# GRIDCHECK` line.
  63 and 55 iterations.
- **c3b**: a net-charged subset (`@atoms 1:1-4`, read out `1:2-4`) whose grid centre still comes
  from all twelve solute atoms (so the centred atoms do not sit at the grid centre), a COG
  gather, a conducting reaction field (`@epsRF 0`, srf = 1), no small-ion list, a 32³ grid at
  0.12 nm. 43 and 25 iterations.

`tests/gromospp_references.rs` runs ours with the same arguments. stdout and the log of the
finite-difference solves are byte-identical. The FFT solves use rustfft where gromos++ uses
FFTW, so the test lets their table columns differ by 1e-7 and the numbers of their log section
(after `# DGRESULT PBC`) by one unit in the last printed digit, with the same iteration counts;
on the machine they were made on, both FFT cases came out byte-identical too.
