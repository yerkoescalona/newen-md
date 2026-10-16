### Changed (1.8 FE-12 — housekeeping)

- **Perturbation topologies are checked as gromosXX checks them.** The pre-2005 block names
  `PERTBOND03`, `PERTBANGLE03`, `PERTIMPDIHEDRAL03`, `PERTATOMPAIR03` and `PERTATOM03` were read as
  their successors; `parse_pttopo` now refuses them with gromosXX's message ("The PERTBOND03 block
  was renamed to PERTBONDSTRETCH.", `in_perturbation.cc:164-170` and the four like it).
  `PERTDIHEDRAL03`, which gromosXX never reads and so refuses as "block PERTDIHEDRAL03 not
  supported!" (`:1809-1813`), is refused with that message. No checked-in `.ptp` uses any of them.
- **A perturbed bond must exist.** `Topology::apply_perturbation` returns
  `Result<NetCharge, String>` and refuses, before changing the topology, a `PERTBONDSTRETCH(H)` or
  `PERTBONDSOFT` bond that is not a bond of the topology: "Perturbation of a non-existing bond in
  PERTBONDSTRETCH block." (`in_perturbation.cc:232-247, 271-286, 420-425`), with the atom numbers
  appended. As gromosXX's `bondMatcher` does, the pair is matched in file order (`2 1` does not
  find bond `1 2`). `md`/`RunRecipe` report it as an inconsistent input, Python's
  `Topology.apply_perturbation` and `Perturbation` raise. Every reference `.ptp` passes the check
  (its perturbed bonds are all topology bonds, in topology order).
- **The Python package takes the workspace version.** `py-gromos/pyproject.toml` declared
  `version = "0.1.0"` against the workspace's 0.0.64; it now declares `dynamic = ["version"]` and
  maturin reads it from `py-gromos/Cargo.toml` (`version.workspace = true`), so
  `importlib.metadata.version("gromos")` and `gromos.__version__` give 0.0.64. `uv.lock` drops
  the static version.
- **CI runs on pushes to `develop`** as well as `main` and `ci` (`ci.yml`, `test-rust.yml`,
  `test-python.yml`; pull requests already ran for every target branch). `release.yml` unchanged.
- The free-energy user guide says that slow growth (DLAMT ≠ 0) runs and warns when λ leaves
  [0, 1], that NSCALE ≠ 0 is refused, and that the two checks above refuse a `.ptp`.

- **Soft-core terms made from nothing.** A `PERTBONDSOFT`, `PERTANGLESOFT` or
  `PERTIMPROPERDIHSOFT` line with state-A type 0 (a bond, angle or improper that exists only in
  state B) underflowed `type − 1` (a panic in debug, a wrapped index in release). Type 0 now
  means "absent in that state" in either column, as in gromosXX (`in_perturbation.cc:435-449,
  765-777, 1003-1015`): `SoftBond`/`SoftAngle`/`SoftImproper::a_type` is `Option<usize>` like
  `b_type`, and the soft kernels give the absent state K = 0 and the other state's equilibrium
  value, the dummy type gromosXX appends (`perturbed_soft_bond_interaction.cc:249-263` and the
  angle and improper equivalents). Type 0 in both columns is refused ("Conversion from bond type
  0 to bond type 0 in PERTBONDSOFT block."). A soft bond being made is exempt from the
  non-existing-bond check, as gromosXX's `t_A > 0` condition makes it, and the reader logs
  gromosXX's "You are making / breaking a bond" warnings. No reference uses type 0, so no
  output changes.
- **Every `.ptp` block the reader does not apply is refused** ("block X not supported!",
  gromosXX `in_perturbation.cc:1803-1820`). Before, such a block was dropped silently. This
  covers `PERTCONSTRAINT03` and `PERTPOLPARAM`, which gromosXX reads and gromos-rs does not
  model, and replaces the special case for `PERTDIHEDRAL03`. `LAMBDADEP` gets gromosXX's own
  message ("… is deprecated. Use LAMBDAS block in the input file"). `SCALEDINTERACTIONS` still
  passes: NSCALE ≠ 0 is refused by the run, and `check/lambdas.pttopo` (the check-suite input,
  run with NSCALE forced to 0) carries the block. No other checked-in `.ptp`/`.pttopo`/`.mpt`
  has a block outside the read set.

### Not done (left for later)

- gromosXX's other existence checks: perturbing a non-existing angle, improper or proper dihedral
  (`in_perturbation.cc:643, 882, 1103`) is still not refused.
- `SCALEDINTERACTIONS` is accepted with NSCALE = 0, where gromosXX errors (its message
  "Scaling turned on but no SCALEDINTERACTIONS block." sits in the scaling-off branch,
  `in_perturbation.cc:1421-1427`).
