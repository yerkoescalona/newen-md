//! `ProviderOrchestrator` against real gromosXX references — the P2.8-2 exit criterion.
//!
//! P2.8-2 originally called for "classical-only configuration reproduces today's 37/37
//! bit-for-bit." That premise doesn't survive contact with P2.6's own documented scope limit:
//! `LjCrfInteraction` is solute-solute atom-level only (no exclusions, no solvent charge-group
//! pairlist shape — see its module docs), so 33 of the 37 reference systems (anything with
//! solvent, or bonded terms this provider never computed in the first place) are out of reach
//! for *any* orchestrator built only from today's providers — that's a scope gap in the provider
//! roster, not something an orchestrator could paper over. Checked before writing this file,
//! same discipline as P2.7 Step 5's premise check.
//!
//! What *is* reachable, and is what "the orchestrator is transparent when only one provider
//! exists" actually needs proving against a real external oracle: the four reference systems
//! that are pure two-atom nonbonded, no bonds, no solvent — `pair_lj`, `pair_lj_mixed` (both
//! zero-charge, pure LJ), `nacl_pair` (vacuum, real ion charges, so CRF is genuinely exercised),
//! `nacl_pair_box` (periodic, real charges — exercises both PBC wrapping and CRF through the
//! orchestrator at once). For each, a `ProviderOrchestrator` holding exactly one
//! `LjCrfInteraction` must reproduce gromosXX's real reported potential energy to the same
//! `1e-8` relative tolerance `test_gromosXX_references.rs` uses everywhere else, *and* match a
//! direct (non-orchestrated) call to the same provider bit-for-bit — proving the orchestrator's
//! bookkeeping doesn't perturb the one-provider case, the actual content of "transparent."

use std::fs;
use std::path::{Path, PathBuf};

use newen_core::configuration::{Box as SimBox, Configuration};
use newen_core::math::{Periodicity, Rectangular, Vacuum};
use newen_core::selection::AtomSelection;
use newen_core::spatial_index::ConfigurationSpatialIndex;
use newen_forces::nonbonded::LjCrfInteraction;
use newen_forces::orchestrator::ProviderOrchestrator;
use newen_forces::provider::PotentialProvider;
use newen_io::coordinate::read_coordinates;
use newen_io::topology::{build_topology, read_topology_file};

const ENERGY_REL_TOL: f64 = 1e-8;

fn ref_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

/// Same extraction as `test_provider_reference.rs` — kept file-local on purpose, see that
/// file's comment on why this small helper isn't shared across test binaries.
fn first_frame_potential(path: &Path) -> f64 {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut in_energy = false;
    let mut vals: Vec<f64> = Vec::new();
    for line in content.lines() {
        let t = line.trim();
        match t {
            "ENERGY03" => {
                in_energy = true;
                vals.clear();
            },
            "END" if in_energy => break,
            _ if in_energy && !t.starts_with('#') && !t.is_empty() => {
                if let Ok(v) = t.parse::<f64>() {
                    vals.push(v);
                }
            },
            _ => {},
        }
    }
    assert!(
        vals.len() >= 3,
        "expected ENERGY03 total/kinetic/potential in {}",
        path.display()
    );
    vals[2]
}

/// One system: directory name under `gromosXX_references/`, the `LjCrfInteraction` cutoff to
/// use, and the reaction-field parameters. `nacl_pair`/`nacl_pair_box` use `RCRF` (the reaction-
/// field cutoff from the `.in`/`.imd` `NONBONDED` block), not `RCUTP` (the pairlist search
/// cutoff) — `LjCrfInteraction` has one `cutoff` field that serves both roles (a known scope
/// limitation: no twin-range support), so this picks the value that matches gromosXX's actual
/// CRF formula. Safe to widen the geometric search this way: `RCRF >= RCUTP` on every system
/// here, so the single pair each system has is still found, and using the *larger* cutoff for a
/// 2-atom system can't spuriously add another pair.
struct RefCase {
    name: &'static str,
    cutoff: f64,
    rf_epsilon: f64,
    rf_kappa: f64,
}

const CASES: &[RefCase] = &[
    // RCUTP=0.8 in the .in file; zero charge makes RCRF irrelevant, so RCUTP is fine here
    // (matches the precedent already validated in test_provider_reference.rs).
    RefCase {
        name: "pair_lj",
        cutoff: 0.8,
        rf_epsilon: 1.0,
        rf_kappa: 0.0,
    },
    RefCase {
        name: "pair_lj_mixed",
        cutoff: 0.8,
        rf_epsilon: 1.0,
        rf_kappa: 0.0,
    },
    // RCUTP=0.8 but RCRF=1.4 — real ion charges, so this one actually exercises CRF.
    RefCase {
        name: "nacl_pair",
        cutoff: 1.4,
        rf_epsilon: 1.0,
        rf_kappa: 0.0,
    },
    // RCUTP=RCUTL=RCRF=0.9, periodic (2 nm cubic box) — exercises PBC wrapping too.
    RefCase {
        name: "nacl_pair_box",
        cutoff: 0.9,
        rf_epsilon: 1.0,
        rf_kappa: 0.0,
    },
];

#[test]
fn single_provider_orchestrator_matches_real_gromosxx_energy_on_every_pure_pair_system() {
    for case in CASES {
        let sys_dir = ref_root().join(case.name);
        let topo_data =
            read_topology_file(sys_dir.join(format!("{}.topo", case.name))).expect("read topology");
        let topo = build_topology(topo_data);
        let coords = read_coordinates(sys_dir.join(format!("{}.conf", case.name)))
            .expect("read coordinates");

        let n_atoms = topo.num_atoms();
        let mut conf = Configuration::new(n_atoms, 1, 1);
        conf.current_mut().pos = coords.positions.clone();

        let periodicity = if coords.box_type == 0 {
            Periodicity::Vacuum(Vacuum)
        } else {
            // `LjCrfInteraction::contribute` derives its own periodicity from
            // `conf.current().box_config`, not from whatever `Periodicity` the caller built —
            // `Configuration::new` otherwise defaults to `Box::vacuum()`. Must set this or a
            // periodic case like nacl_pair_box silently uses unwrapped distances (found and
            // documented while building the P2.8-1 exit test, same trap).
            conf.current_mut().box_config =
                SimBox::rectangular(coords.box_dims.x, coords.box_dims.y, coords.box_dims.z);
            Periodicity::Rectangular(Rectangular::new(coords.box_dims))
        };
        let index = ConfigurationSpatialIndex::new(&conf, &periodicity);
        let region = AtomSelection::all(n_atoms);

        // ── Direct call: the provider on its own, no orchestrator involved ────────────────
        let direct = LjCrfInteraction::new(case.cutoff, 1.0, case.rf_epsilon, case.rf_kappa)
            .contribute(&region, &topo, &conf, &index)
            .unwrap_or_else(|e| panic!("{}: direct contribute failed: {e}", case.name));

        // ── Through the orchestrator: exactly one registered provider ──────────────────────
        let mut orchestrator = ProviderOrchestrator::new();
        orchestrator.register(
            AtomSelection::all(n_atoms),
            Box::new(LjCrfInteraction::new(
                case.cutoff,
                1.0,
                case.rf_epsilon,
                case.rf_kappa,
            )),
        );
        let orchestrated = orchestrator
            .evaluate(&topo, &conf, &index)
            .unwrap_or_else(|e| panic!("{}: orchestrator evaluate failed: {e}", case.name));

        // Transparency: orchestrating one provider must not perturb its result at all.
        assert_eq!(
            orchestrated.energy, direct.energy,
            "{}: orchestrator energy diverged from a direct provider call",
            case.name
        );

        // The real oracle: gromosXX's own reported potential energy for this system. No bonds,
        // no exclusions, no restraints on any of these four — but for charged systems (nacl_pair,
        // nacl_pair_box) gromosXX's total also includes the unconditional per-atom RF self-energy
        // (`rf_excluded_interactions`, added regardless of exclusions — every solute atom gets
        // it), which `LjCrfInteraction` deliberately does not compute (it's a separate GROMOS
        // term, not part of the pairlist LJ+CRF loop — same accounting `interaction.rs`'s own
        // `cross_checks_against_single_point_energy_oracle` test already establishes). Zero for
        // the two uncharged systems, so this correction is safe to apply uniformly.
        let crf = newen_forces::nonbonded::CRFParameters::new(
            case.cutoff,
            1.0,
            case.rf_epsilon,
            case.rf_kappa,
        );
        let self_energy: f64 = topo
            .charge
            .iter()
            .map(|q| -0.5 * q * q * newen_core::units::four_pi_eps_i * crf.crf_cut)
            .sum();

        let expected = first_frame_potential(&sys_dir.join("expected/energies.tre"));
        let rel_diff = (orchestrated.energy - (expected - self_energy)).abs() / expected.abs();
        assert!(
            rel_diff < ENERGY_REL_TOL,
            "{}: orchestrated energy {} vs real gromosXX energy {} minus self-energy {} = {} \
             (rel diff {})",
            case.name,
            orchestrated.energy,
            expected,
            self_energy,
            expected - self_energy,
            rel_diff
        );
    }
}

/// The orchestrator itself, independent of any reference system: registering zero providers
/// must give exactly zero, and this must not panic or require a nonempty registry.
#[test]
fn empty_orchestrator_on_a_real_topology_is_zero() {
    let sys_dir = ref_root().join("pair_lj");
    let topo =
        build_topology(read_topology_file(sys_dir.join("pair_lj.topo")).expect("read topology"));
    let coords = read_coordinates(sys_dir.join("pair_lj.conf")).expect("read coordinates");
    let n_atoms = topo.num_atoms();
    let mut conf = Configuration::new(n_atoms, 1, 1);
    conf.current_mut().pos = coords.positions;
    let periodicity = Periodicity::Vacuum(Vacuum);
    let index = ConfigurationSpatialIndex::new(&conf, &periodicity);

    let mut orchestrator = ProviderOrchestrator::new();
    let result = orchestrator.evaluate(&topo, &conf, &index).unwrap();
    assert_eq!(result.energy, 0.0);
    assert!(result.forces.is_empty());
}
