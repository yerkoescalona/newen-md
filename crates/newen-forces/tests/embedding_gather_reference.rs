//! ** P2.7 Step 1** — falsify or confirm assumption **A1**: does today's
//! [`PotentialProvider`] signature already carry everything classic electrostatic embedding
//! needs, without a trait change?
//!
//! Electrostatic embedding (Poliak et al. 2025, *J. Comput. Chem.* 46:e70053) requires exactly
//! four things, and this test checks each is reachable from what `contribute()` already
//! receives — no new channel, no new parameter:
//!
//! | EE needs | reachable today via |
//! |---|---|
//! | which MM atoms are near the QM zone | `neigh.neighbor_pairs(region, cutoff)` — already returns pairs with only *one* endpoint in `region` |
//! | their partial charges | `topo.charge` |
//! | their positions | `conf.current().pos` |
//! | a place to put forces on them | `Contribution.forces: Vec<(usize, Vec3)>` — arbitrary global indices |
//!
//! The first row is the load-bearing one, and it is the exact set `SchNetInteraction` currently
//! throws away (`schnet.rs`, the cross-boundary skip added in P2.6) — so this test also pins
//! down what that skip is discarding, before P2.7 Step 2 makes the policy explicit.
//!
//! Deliberately **read-only**: it uses the real fixture and today's API only, and changes no
//! production code. If it passes, A1 holds and Steps 2-5 need no signature change. Fixture is
//! the same real BuRNN tutorial system as `schnet_burnn_reference.rs` (gitignored `.local/`
//! clone; skips gracefully when absent).

use newen_core::configuration::Configuration;
use newen_core::math::{Periodicity, Rectangular};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::{ConfigurationSpatialIndex, SpatialIndex};
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};
use std::path::PathBuf;

/// Real `RCUTQM` from the tutorial's `md_burnn/md.imd` QMMM block.
const RCUTQM: f64 = 1.4;
/// Real `QMZONE` from `md_burnn/meoh.qmm`: the 6 methanol atoms, which are the system's first 6.
const QM_ZONE_SIZE: usize = 6;

fn tutorial_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../.local/gromos_tutorial_livecoms/tutorial_files/t_06")
}

#[test]
fn a1_electrostatic_embedding_inputs_are_reachable_from_todays_trait() {
    let root = tutorial_root();
    let topo_path = root.join("topo/meoh_54a7.top");
    let cnf_path = root.join("eq/eq_meoh_5.cnf");

    if !topo_path.exists() || !cnf_path.exists() {
        eprintln!(
            "skipping: real BuRNN tutorial not cloned at {} \
             (git clone github.com/biomos/gromos_tutorial_livecoms into .local/)",
            root.display()
        );
        return;
    }

    // ── Build the real system exactly as the md binary does: read topology, then solvate ──
    let mut topo = build_topology(read_topology_file(&topo_path).expect("read t_06 topology"));
    let coords = read_coordinates(&cnf_path).expect("read t_06 eq_meoh_5.cnf");

    let n_atoms = coords.positions.len();
    let n_solute = topo.num_solute_atoms();
    let atoms_per_solvent = topo.solvent_atom_template.len();
    assert!(
        atoms_per_solvent > 0,
        "t_06 topology should carry an SPC solvent template"
    );
    assert_eq!(
        (n_atoms - n_solute) % atoms_per_solvent,
        0,
        "coordinate count should be solute + whole solvent molecules"
    );
    topo.solvate((n_atoms - n_solute) / atoms_per_solvent);
    assert_eq!(
        topo.num_atoms(),
        n_atoms,
        "solvated topology must match the coordinate file"
    );

    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions.clone();
    let periodicity = Periodicity::Rectangular(Rectangular::new(coords.box_dims));
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
    let region = AtomSelection::from_indices((0..QM_ZONE_SIZE).collect(), n_atoms).unwrap();

    // ── EE input 1: which MM atoms are within the QM cutoff ──────────────────────────────
    // `neighbor_pairs` returns pairs where *at least one* endpoint is in `region` — so the
    // cross-boundary subset is exactly the QM↔MM pair set electrostatic embedding acts on.
    let pairs = index.neighbor_pairs(&region, RCUTQM);
    let in_region = |i: usize| i < QM_ZONE_SIZE;

    let mut mm_neighbors: Vec<usize> = pairs
        .iter()
        .filter_map(|&(i, j, _shift)| match (in_region(i), in_region(j)) {
            (true, false) => Some(j),
            (false, true) => Some(i),
            _ => None, // intra-QM pair: not an embedding partner
        })
        .collect();
    mm_neighbors.sort_unstable();
    mm_neighbors.dedup();

    assert!(
        !mm_neighbors.is_empty(),
        "expected MM atoms within {RCUTQM} nm of the QM zone in a solvated system — \
         if this is empty, the cross-boundary pairs EE depends on are not reachable"
    );
    assert!(
        mm_neighbors.iter().all(|&i| i >= QM_ZONE_SIZE),
        "embedding partners must all be outside the QM zone"
    );

    // ── EE input 2 + 3: their charges and positions ───────────────────────────────────────
    let embedding_charges: Vec<(usize, f64, newen_core::math::Vec3)> = mm_neighbors
        .iter()
        .map(|&i| (i, topo.charge[i], conf.current().pos[i]))
        .collect();

    assert_eq!(embedding_charges.len(), mm_neighbors.len());
    assert!(
        embedding_charges.iter().any(|&(_, q, _)| q.abs() > 1e-12),
        "embedding partners should carry non-zero partial charges (SPC water is charged) — \
         a QM zone embedded in neutral-charge atoms would not be polarized at all"
    );
    // SPC water charge magnitudes are ~0.41 e (H) / ~0.82 e (O); anything wildly outside that
    // means we are reading the wrong atoms, not merely a different force field.
    assert!(
        embedding_charges.iter().all(|&(_, q, _)| q.abs() < 2.0),
        "partial charges outside a physically sane range — likely an index misalignment"
    );

    // ── EE output: forces on those MM atoms already fit `Contribution.forces` ─────────────
    // `Vec<(usize, Vec3)>` over arbitrary *global* indices, so no new output channel is
    // needed to return embedding forces onto atoms outside the region. Construct one to prove
    // the shape accepts it (values are placeholders — Step 3 computes real ones).
    let embedding_forces: Vec<(usize, newen_core::math::Vec3)> = embedding_charges
        .iter()
        .map(|&(i, _, _)| (i, newen_core::math::Vec3::ZERO))
        .collect();
    let contribution = newen_forces::provider::Contribution {
        energy: 0.0,
        forces: embedding_forces,
        virial: newen_core::math::Mat3::ZERO,
        extra: newen_forces::provider::ProviderExtra::default(),
    };
    assert!(
        contribution.forces.iter().all(|&(i, _)| i >= QM_ZONE_SIZE),
        "Contribution.forces must be able to carry forces on atoms outside the region"
    );

    eprintln!(
        "A1 holds: {} MM embedding partners within {RCUTQM} nm of the {QM_ZONE_SIZE}-atom QM zone, \
         charges and positions readable, forces representable — no trait change needed.",
        mm_neighbors.len()
    );
}

/// ** P2.7 Step 3 at scale.** The unit tests in `nonbonded/embedding.rs` validate
/// [`ElectrostaticEmbedding`] against a closed-form two-charge oracle; this runs the same
/// provider on the real solvated BuRNN system — a 6-atom QM zone embedded in ~1.4k MM point
/// charges across a periodic box — and finite-difference-checks a force on a real **MM** atom.
///
/// That combination is the actual Step 3 exit criterion: forces reaching atoms *outside* the
/// region (P5), on real data, under periodic boundary conditions.
#[test]
fn embedding_forces_on_real_mm_atoms_match_finite_differences() {
    use newen_core::math::Vec3;
    use newen_forces::nonbonded::ElectrostaticEmbedding;
    use newen_forces::provider::PotentialProvider;

    let root = tutorial_root();
    let topo_path = root.join("topo/meoh_54a7.top");
    let cnf_path = root.join("eq/eq_meoh_5.cnf");
    if !topo_path.exists() || !cnf_path.exists() {
        eprintln!(
            "skipping: real BuRNN tutorial not cloned at {}",
            root.display()
        );
        return;
    }

    let mut topo = build_topology(read_topology_file(&topo_path).expect("read t_06 topology"));
    let coords = read_coordinates(&cnf_path).expect("read t_06 eq_meoh_5.cnf");
    let n_atoms = coords.positions.len();
    let n_solute = topo.num_solute_atoms();
    topo.solvate((n_atoms - n_solute) / topo.solvent_atom_template.len());

    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions.clone();
    let periodicity = Periodicity::Rectangular(Rectangular::new(coords.box_dims));
    let region = AtomSelection::from_indices((0..QM_ZONE_SIZE).collect(), n_atoms).unwrap();

    // Stand-in QM charges on the methanol zone. A real run refreshes these from the QM engine
    // each step; the force path being validated here is identical either way (Poliak path (c)).
    let mut qm_charges = vec![0.0; n_atoms];
    for (i, q) in [0.27, -0.64, 0.04, 0.04, 0.04, 0.41].iter().enumerate() {
        qm_charges[i] = *q;
    }
    let mut provider = ElectrostaticEmbedding::new(RCUTQM, qm_charges);

    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
    let contribution = provider
        .contribute(&region, &topo, &conf, &index)
        .expect("embedding should run on the real solvated system");

    assert!(contribution.energy.is_finite());

    // ** P2.8-3, real-data half.** `trace(virial) == energy` is an *exact* per-pair
    // identity for a `1/r` potential (`nonbonded/embedding.rs`'s
    // `virial_trace_equals_energy_for_a_pure_coulomb_pair` proves it on one pair; summed over
    // pairs it's still exact by linearity). Checking it here — thousands of pairs, real periodic
    // shift vectors — is a much sharper test of the PBC bookkeeping than the vacuum unit tests
    // can be: a wrong minimum-image convention would perturb `r_vec` in the force term but not
    // (consistently) in the energy term, breaking this identity even though energy and forces
    // might each look individually plausible.
    let virial_trace =
        contribution.virial.x_axis.x + contribution.virial.y_axis.y + contribution.virial.z_axis.z;
    assert!(
        (virial_trace - contribution.energy).abs() < 1e-8 * contribution.energy.abs().max(1.0),
        "trace(virial) = {virial_trace} vs energy = {} on the real periodic t_06 system",
        contribution.energy
    );

    let mm_forced: Vec<usize> = contribution
        .forces
        .iter()
        .map(|&(i, _)| i)
        .filter(|&i| i >= QM_ZONE_SIZE)
        .collect();
    assert!(
        !mm_forced.is_empty(),
        "embedding must put forces on MM atoms outside the QM zone"
    );

    // Pick the MM atom carrying the largest force — the most sensitive FD probe.
    let (probe, probe_force) = contribution
        .forces
        .iter()
        .filter(|&&(i, _)| i >= QM_ZONE_SIZE)
        .max_by(|a, b| a.1.length().partial_cmp(&b.1.length()).unwrap())
        .copied()
        .unwrap();

    let mut energy_at = |positions: &[Vec3]| -> f64 {
        let mut c2 = Configuration::new(n_atoms, 1, 1);
        c2.current_mut().pos = positions.to_vec();
        let idx = ConfigurationSpatialIndex::new(&c2, &periodicity);
        provider
            .contribute(&region, &topo, &c2, &idx)
            .unwrap()
            .energy
    };

    let h = 1e-6;
    for axis in 0..3 {
        let delta = match axis {
            0 => Vec3::new(h, 0.0, 0.0),
            1 => Vec3::new(0.0, h, 0.0),
            _ => Vec3::new(0.0, 0.0, h),
        };
        let mut plus = coords.positions.clone();
        let mut minus = coords.positions.clone();
        plus[probe] += delta;
        minus[probe] -= delta;

        let fd = -(energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
        let analytic = match axis {
            0 => probe_force.x,
            1 => probe_force.y,
            _ => probe_force.z,
        };
        assert!(
            (fd - analytic).abs() < 1e-3,
            "real MM atom {probe} axis {axis}: finite-diff {fd} vs provider {analytic}"
        );
    }

    eprintln!(
        "Step 3 holds on real data: E_embed = {:.6} kJ/mol, forces on {} MM atoms, \
         FD-verified on MM atom {probe}.",
        contribution.energy,
        mm_forced.len()
    );
}
