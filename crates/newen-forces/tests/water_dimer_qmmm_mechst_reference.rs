//! A real gromosXX QM/MM nonbonded energy, not a synthetic fixture or a self-consistency check.
//!
//! Every other QM/MM-seam test in this suite validates against a closed-form oracle, finite
//! differences, or gromosXX *source code* read directly — none against a real gromosXX QM/MM
//! *run's output*, because that needs either a QM engine (not installed here) or a raw `.tre`
//! nobody has archived. This test closes that specific gap for the one case that needs neither:
//! mechanical embedding with constant QM charges (`NTQMMM=-1`) leaves the classical nonbonded
//! loop completely untouched by QM/MM (see `gromosXX_qmmm_references/water_dimer_mechst/
//! README.md` for the exact source-code citation), so its reported nonbonded energy is plain
//! classical LJ+CRF — reproducible with nothing but `newen-io` + `newen-forces`.
//!
//! Fixture: `WATER_DIMER/` from the dataset accompanying Poliak et al. 2025 (*J. Comput. Chem.*
//! 46:e70053), the paper describing GROMOS's mainline QM/MM interface — a different, larger real
//! system than the BuRNN tutorial's `t_06` used elsewhere in this suite. Small excerpt checked in
//! directly (not gitignored under `.local/`, unlike the multi-GB tutorial/gromosXX clones) under
//! CC BY 4.0 with attribution; see the fixture directory's `README.md` for provenance and license.

use std::fs;
use std::path::PathBuf;

use newen_core::math::Vec3;
use newen_forces::energy::{single_point_energy, EnergyParams};
use newen_io::coordinate::read_coordinates;
use newen_io::imd::read_imd_file;
use newen_io::topology::{build_topology, read_topology_file};

const ENERGY_REL_TOL: f64 = 1e-8;

fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/gromosXX_qmmm_references/water_dimer_mechst")
}

/// The `t=0` row of the real `ene_ana` nonbonded-energy time series — the value at exactly the
/// starting configuration this test also loads.
fn first_frame_nonb(path: &std::path::Path) -> f64 {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    for line in content.lines() {
        let t = line.trim();
        if t.starts_with('#') || t.is_empty() {
            continue;
        }
        let fields: Vec<&str> = t.split_whitespace().collect();
        assert_eq!(
            fields.len(),
            2,
            "expected 'time nonb' rows in {}",
            path.display()
        );
        let time: f64 = fields[0].parse().expect("time column");
        if time == 0.0 {
            return fields[1].parse().expect("nonb column");
        }
    }
    panic!("no t=0 row found in {}", path.display());
}

#[test]
fn mechanical_embedding_constant_charge_nonb_matches_real_gromosxx() {
    let root = fixture_root();

    let mut topo = build_topology(
        read_topology_file(root.join("water_dimer.top")).expect("read water_dimer.top"),
    );
    let coords = read_coordinates(root.join("water_dimer.cnf")).expect("read water_dimer.cnf");
    let imd = read_imd_file(root.join("water_dimer.imd")).expect("read water_dimer.imd");

    // NSM=1 in the real .imd's SYSTEM block: one solvent water molecule.
    let n_solute = topo.num_solute_atoms();
    let atoms_per_solvent = topo.solvent_atom_template.len();
    topo.solvate((coords.positions.len() - n_solute) / atoms_per_solvent);

    // One solute water (would-be QM zone) + NSM=1 solvent water (would-be MM zone) — see
    // README.md. Real values straight from the real .imd, not hand-transcribed.
    let params = EnergyParams {
        cutoff: imd.rcrf,
        epsilon: 1.0,
        rf_epsilon: imd.epsrf,
        rf_kappa: imd.appak,
        pairlist_freq: 1,
        ntf: [true, true, true, true],
        atoms_per_solvent: 3,
        quartic_bonds: true,
    };

    let result = single_point_energy(&topo, &coords.positions, Vec3::ZERO, &params);
    let nonb = result.lj + result.crf;

    let expected = first_frame_nonb(&root.join("expected/am1_mechst_nonb.out"));
    let rel_diff = (nonb - expected).abs() / expected.abs();
    assert!(
        rel_diff < ENERGY_REL_TOL,
        "nonb {nonb} vs real gromosXX QM/MM run's reported nonb {expected} (rel diff {rel_diff})"
    );
    eprintln!("computed nonb = {nonb} kJ/mol vs real gromosXX QM/MM run's {expected} kJ/mol (rel diff {rel_diff:e})");
}
