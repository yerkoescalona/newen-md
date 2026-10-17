//! gromosXX's own regression oracle, ported: `md++/src/check/{aladip,aladip_unperturbed,
//! aladip_special}.t.cc` hard-code per-term energies for the alanine dipeptide in 20 SPC waters
//! (`shared/check/`, byte-identical inputs). Each term is evaluated on its own at the starting
//! configuration — independent of the `md` binary's trajectory comparison, and a *second source*
//! for the bonded, perturbed, soft-core and reaction-field terms (, cross-cutting).
//!
//! What the C++ does (`check_forcefield.cc`): for every interaction term, zero the energies,
//! call the term alone, compare `potential_total + special_total` to the hard-coded value
//! (`CHECK_APPROX_EQUAL`, δ = 1e-3 … 1e-8), then a finite-difference check of the λ-derivative
//! (ε = 0.001). Here: the bonded terms through the `newen-forces` functions the force field
//! calls, the nonbonded totals through the assembled force field at step 0 (`lj_total + crf_total`
//! is the nonbonded term; the bonded totals must equal the sum of the split), with `NSLFEXCL`
//! toggled for the two reaction-field variants.
//!
//! `lambdas.t.cc` (individual λ, LAMBDAS on vs off) is ported for its bonded terms and its
//! individual-λ loop; its NonBonded values need SCALEDINTERACTIONS (NSCALE = 1), which newen
//! does not implement. Not portable here: `scaling.t.cc`, `aladip_ls.t.cc` (lattice sum), the CG
//! function, and the dihedral / angle / X-ray / local elevation / order-parameter restraints
//! (1.6 — parked or next).

use std::path::{Path, PathBuf};

use newen_core::configuration::Configuration;
use newen_core::Topology;
use newen_forces::bonded::perturbed::{
    calculate_perturbed_angle_forces, calculate_perturbed_bond_forces,
    calculate_perturbed_dihedral_forces, calculate_perturbed_improper_dihedral_forces,
    calculate_soft_angle_forces, calculate_soft_bond_forces, calculate_soft_improper_forces,
};
use newen_forces::bonded::{calculate_bonded_forces_ntf, calculate_perturbed_bonded_forces};
use newen_io::coordinate::read_coordinates;
use newen_io::imd::ImdParameters;
use newen_io::topology::{build_topology, read_topology_file};
use newen_run::{
    build_sequence_from_imd, prepare_system, read_imd, start, Coordinates, Prepared, RunInputs,
    RunOptions,
};

fn shared() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references/shared")
}

/// `check_forcefield.cc`'s `CHECK_APPROX_EQUAL`: relative deviation, the suite's own δ is
/// 1e-3 for every energy. Its hard-coded table predates the current gromosXX: the binary built
/// from `.local/gromosXX` (2026-08) reproduces the bonded values to 1e-9 but the nonbonded ones
/// only to ~3e-6 (e.g. NonBonded −50.1969469 now vs −50.196817 in the table). So each value is
/// checked twice: against the current binary's number (`NOW`, tight — newen matches it to
/// ten digits) and against the table with the suite's δ (`SUITE`).
fn close(actual: f64, expected: f64, rel: f64, what: &str) {
    let diff = if expected != 0.0 {
        ((actual - expected) / expected).abs()
    } else {
        actual.abs()
    };
    assert!(
        diff <= rel,
        "{what}: got {actual:.9}, expected {expected:.9} (rel diff {diff:.2e}, δ {rel:.0e})"
    );
}

/// Tolerances: current binary (values printed by `check_interaction` with `setprecision(12)`),
/// and the suite's hard-coded table with its own δ.
const NOW: f64 = 1e-7;
const SUITE: f64 = 1e-3;

/// Check one term against both sources.
fn term(actual: f64, now: f64, suite: f64, what: &str) {
    close(actual, now, NOW, &format!("{what} vs current gromosXX"));
    close(
        actual,
        suite,
        SUITE,
        &format!("{what} vs check-suite table"),
    );
}

struct System {
    imd: ImdParameters,
    prepared: Prepared,
    inputs: RunInputs,
}

fn load(input: &str, pttopo: bool, distrest: bool) -> System {
    load_with(input, pttopo.then_some("aladip.pttopo"), distrest)
}

/// [`load`] with the perturbation topology `pttopo` (relative to `shared/`).
fn load_with(input: &str, pttopo: Option<&str>, distrest: bool) -> System {
    let s = shared();
    let parsed = read_topology_file(s.join("aladip.topo")).expect("topology");
    let physical_constants = parsed.physical_constants;
    let topology = build_topology(parsed);
    let coords: Coordinates = read_coordinates(s.join("aladip.conf"))
        .expect("coordinates")
        .into();
    let imd = read_imd(s.join("check").join(input)).expect("imd");
    let inputs = RunInputs {
        pttopo: pttopo.map(|p| s.join(p)),
        distrest: distrest.then(|| s.join("check/aladip.distrest")),
        ..Default::default()
    };
    let prepared = prepare_system(&imd, topology, physical_constants, coords, &inputs)
        .expect("prepare_system");
    System {
        imd,
        prepared,
        inputs,
    }
}

/// Energies of the starting configuration: build the force field for `imd` (optionally with
/// NSLFEXCL and RLAM overridden), run `init` + step 0, read the accumulated totals.
fn step0_energies(sys: &System, nslfexcl: Option<i32>, rlam: Option<f64>) -> newen_core::Energy {
    step0_state(sys, nslfexcl, rlam).energies
}

/// The old state after step 0: its energies and their λ derivatives.
fn step0_state(sys: &System, nslfexcl: Option<i32>, rlam: Option<f64>) -> newen_core::State {
    let mut imd = sys.imd.clone();
    if let Some(v) = nslfexcl {
        imd.nslfexcl = v;
    }
    if let Some(l) = rlam {
        imd.rlam = l;
    }
    // The check inputs couple solute and solvent to separate baths; the recipe supports one
    // bath (MULTIBATH NBATHS=1). The energies of the starting configuration do not depend on
    // the thermostat, so evaluate them with one bath over all atoms.
    if let Some(bath) = imd.temp_bath.first_mut() {
        bath.num_bath_groups = 1;
        bath.temp0.truncate(1);
        bath.tau.truncate(1);
        let last = bath.dof_sets.last().map(|d| d[0]).unwrap_or(0);
        bath.dof_sets = vec![[last, 1, 1]];
    }
    // prepare again: `Prepared` owns the configuration that `start` advances
    let s = shared();
    let parsed = read_topology_file(s.join("aladip.topo")).expect("topology");
    let physical_constants = parsed.physical_constants;
    let topology = build_topology(parsed);
    let coords: Coordinates = read_coordinates(s.join("aladip.conf"))
        .expect("coordinates")
        .into();
    let mut prepared = prepare_system(&imd, topology, physical_constants, coords, &sys.inputs)
        .expect("prepare_system");
    let built = build_sequence_from_imd(&imd, &prepared, &sys.inputs, &RunOptions::default())
        .expect("build");
    let mut seq = built.sequence;
    let state = newen_core::algorithm::SimulationState::new(imd.dt, imd.nstlim);
    start(
        &mut seq,
        &prepared.topology,
        &mut prepared.configuration,
        &state,
    )
    .expect("start");
    prepared.configuration.old().clone()
}

/// The check suite's reference values are single molecules in vacuum: no box, no minimum image.
fn vacuum() -> newen_core::math::Periodicity {
    newen_core::math::Periodicity::Vacuum(newen_core::math::Vacuum)
}

fn bonded_split(topo: &Topology, conf: &Configuration, lambda: f64) -> [f64; 12] {
    let unpert = calculate_bonded_forces_ntf(
        topo,
        conf,
        &vacuum(),
        true,
        newen_forces::bonded::CovalentForm::default(),
        true,
        true,
        true,
        true,
        &mut newen_forces::bonded::BondedGroupEnergies::default(),
    );
    let dl = 1.0; // NLAM = 1 → dλ/dλ = 1
    let pb = calculate_perturbed_bond_forces(topo, conf, &vacuum(), lambda, dl);
    let sb = calculate_soft_bond_forces(topo, conf, &vacuum(), lambda, dl);
    let pa = calculate_perturbed_angle_forces(topo, conf, &vacuum(), lambda, dl);
    let sa = calculate_soft_angle_forces(topo, conf, &vacuum(), lambda, dl);
    let pi = calculate_perturbed_improper_dihedral_forces(topo, conf, &vacuum(), lambda, dl);
    let si = calculate_soft_improper_forces(topo, conf, &vacuum(), lambda, dl);
    let pd = calculate_perturbed_dihedral_forces(topo, conf, &vacuum(), lambda, dl);
    [
        unpert.bond_energy,
        pb.energy,
        sb.energy,
        unpert.angle_energy,
        pa.energy,
        sa.energy,
        unpert.improper_energy,
        pi.energy,
        si.energy,
        unpert.dihedral_energy,
        pd.energy,
        // λ-derivatives, summed over the perturbed bonded terms
        pb.lambda_derivative
            + sb.lambda_derivative
            + pa.lambda_derivative
            + sa.lambda_derivative
            + pi.lambda_derivative
            + si.lambda_derivative
            + pd.lambda_derivative,
    ]
}

// ── aladip.t.cc: perturbed, λ = 0.125, ALPHLJ = ALPHC = 1 ────────────────────────────────

#[test]
fn aladip_bonded_terms_match_the_check_suite() {
    let sys = load("aladip.in", true, false);
    let lambda = sys.imd.rlam;
    let e = bonded_split(&sys.prepared.topology, &sys.prepared.configuration, lambda);
    term(e[0], 18.0538111135, 18.053811, "QuarticBond");
    term(e[1], 1.14956775938, 1.149568, "PerturbedQuarticBond");
    term(e[2], 0.00128067341101, 0.001281, "PerturbedSoftBond");
    term(e[3], 11.8704063662, 11.870406, "Angle");
    term(e[4], 0.714818210072, 0.714818, "PerturbedAngle");
    term(e[5], 0.262309698988, 0.262310, "PerturbedSoftAngle");
    term(e[6], 0.809250083318, 0.809250, "ImproperDihedral");
    term(e[7], 2.6427799241, 2.642780, "PerturbedImproperDihedral");
    term(e[8], 0.136206551853, 0.136207, "PerturbedSoftImproper");
    term(e[9], 2.25520636037, 2.255206, "Dihedral");
    term(e[10], 13.3146020649, 13.314602, "PerturbedDihedral");

    // The force field books the same split into its totals.
    let st = step0_state(&sys, None, None);
    let tot = &st.energies;
    close(
        tot.bond_total,
        e[0] + e[1] + e[2],
        1e-9,
        "bond_total = unperturbed + perturbed + soft",
    );
    close(tot.angle_total, e[3] + e[4] + e[5], 1e-9, "angle_total");
    close(
        tot.improper_total,
        e[6] + e[7] + e[8],
        1e-9,
        "improper_total",
    );
    close(tot.dihedral_total, e[9] + e[10], 1e-9, "dihedral_total");
    // dH/dλ per term, as gromosXX books it (soft terms with their kind)
    let d = &st.perturbed_energy_derivatives;
    let (topo, conf) = (&sys.prepared.topology, &sys.prepared.configuration);
    let v = vacuum();
    let dl = |f: fn(&Topology, &Configuration, &newen_core::math::Periodicity, f64, f64) -> _| {
        let r: newen_forces::bonded::ForceEnergyLambda = f(topo, conf, &v, lambda, 1.0);
        r.lambda_derivative
    };
    let bond = dl(calculate_perturbed_bond_forces) + dl(calculate_soft_bond_forces);
    let angle = dl(calculate_perturbed_angle_forces) + dl(calculate_soft_angle_forces);
    let improper =
        dl(calculate_perturbed_improper_dihedral_forces) + dl(calculate_soft_improper_forces);
    let dihedral = dl(calculate_perturbed_dihedral_forces);
    close(d.bond_total, bond, 1e-9, "dH/dλ of the bonds");
    close(d.angle_total, angle, 1e-9, "dH/dλ of the angles");
    close(d.improper_total, improper, 1e-9, "dH/dλ of the impropers");
    close(d.dihedral_total, dihedral, 1e-9, "dH/dλ of the dihedrals");
    close(
        d.bond_total + d.angle_total + d.improper_total + d.dihedral_total,
        e[11],
        1e-9,
        "dH/dλ of the bonded terms",
    );
}

#[test]
fn aladip_bonded_lambda_derivatives_are_finite_differences() {
    // check_lambda_derivative: (E(λ+ε) − E(λ−ε)) / 2ε ≈ dE/dλ, ε = 0.001, δ = 0.001
    let sys = load("aladip.in", true, false);
    let (topo, conf) = (&sys.prepared.topology, &sys.prepared.configuration);
    let lambda = sys.imd.rlam;
    let eps = 0.001;
    let energy = |l: f64| {
        let e = bonded_split(topo, conf, l);
        e[1] + e[2] + e[4] + e[5] + e[7] + e[8] + e[10]
    };
    let fd = (energy(lambda + eps) - energy(lambda - eps)) / (2.0 * eps);
    let analytic = bonded_split(topo, conf, lambda)[11];
    assert!(
        (fd - analytic).abs() < 1e-3 * analytic.abs().max(1.0),
        "bonded dH/dλ: finite difference {fd:.6} vs analytic {analytic:.6}"
    );
}

#[test]
fn aladip_nonbonded_matches_the_check_suite_for_both_reaction_field_variants() {
    let sys = load("aladip.in", true, false);
    // rf_excluded = 0: "NonBonded"; rf_excluded = 1 (NSLFEXCL): "NonBonded_newRF"
    let e0 = step0_energies(&sys, Some(0), None);
    term(
        e0.lj_total + e0.crf_total,
        -50.1969469,
        -50.196817,
        "NonBonded (NSLFEXCL=0)",
    );
    let e1 = step0_energies(&sys, Some(1), None);
    term(
        e1.lj_total + e1.crf_total,
        -84.092632735,
        -84.092443,
        "NonBonded_newRF (NSLFEXCL=1)",
    );
}

#[test]
fn aladip_nonbonded_lambda_derivative_is_a_finite_difference() {
    let sys = load("aladip.in", true, false);
    let lambda = sys.imd.rlam;
    let eps = 0.001;
    let nb = |l: f64| {
        let e = step0_energies(&sys, Some(1), Some(l));
        e.lj_total + e.crf_total
    };
    let fd = (nb(lambda + eps) - nb(lambda - eps)) / (2.0 * eps);
    let d = step0_state(&sys, Some(1), None).perturbed_energy_derivatives;
    let analytic = d.lj_total + d.crf_total;
    assert!(
        (fd - analytic).abs() < 1e-3 * analytic.abs().max(1.0),
        "nonbonded dH/dλ: finite difference {fd:.6} vs analytic {analytic:.6}"
    );
}

// ── aladip_unperturbed.t.cc: no perturbation topology ────────────────────────────────────

#[test]
fn aladip_unperturbed_terms_match_the_check_suite() {
    let sys = load("aladip_unperturbed.in", false, false);
    let e = bonded_split(&sys.prepared.topology, &sys.prepared.configuration, 0.0);
    term(e[0], 19.3560714867, 19.356071, "QuarticBond");
    term(e[3], 12.6704131861, 12.670413, "Angle");
    term(e[6], 1.48547050344, 1.485471, "ImproperDihedral");
    term(e[9], 6.29262384561, 6.292624, "Dihedral");
    let e0 = step0_energies(&sys, Some(0), None);
    term(
        e0.lj_total + e0.crf_total,
        -50.8922163577,
        -50.892127,
        "NonBonded (NSLFEXCL=0)",
    );
    let e1 = step0_energies(&sys, Some(1), None);
    term(
        e1.lj_total + e1.crf_total,
        -67.6141086463,
        -67.613934,
        "NonBonded_newRF (NSLFEXCL=1)",
    );
}

#[test]
fn aladip_unperturbed_atomic_cutoff_matches_the_check_suite() {
    // check_atomic_cutoff: PAIRLIST TYPE = atomic instead of chargegroup
    let mut sys = load("aladip_unperturbed.in", false, false);
    sys.imd.type_ = 1; // PAIRLIST TYPE 1 = atomic
    let e = step0_energies(&sys, Some(0), None);
    // `check_atomic_cutoff` does not go through `check_interaction`; only the table value.
    close(
        e.lj_total + e.crf_total,
        -50.791941,
        SUITE,
        "NonBonded_atomic vs check-suite table",
    );
}

// ── aladip_special.t.cc: restraints ──────────────────────────────────────────────────────

#[test]
fn aladip_distance_restraints_match_the_check_suite() {
    // aladip_distres.in = aladip_special.in without the restraint blocks newen does not
    // model (angle, dihedral, X-ray, local elevation, order parameter — 1.6).
    // `check_interaction` evaluates the unperturbed (257.189539) and the perturbed
    // (195.899012) distance restraints separately; here they land in one total.
    let sys = load("aladip_distres.in", true, true);
    let e = step0_energies(&sys, None, None);
    close(
        e.distanceres_total,
        257.189539 + 195.899012,
        SUITE,
        "DistanceRestraint + PerturbedDistanceRestraint vs check-suite table",
    );
}

#[test]
#[ignore = "dihedral restraints are not implemented (1.6 — next); check value 2127.910749"]
fn aladip_dihedral_restraints_match_the_check_suite() {
    let sys = load("aladip_special.in", false, false);
    let e = step0_energies(&sys, None, None);
    close(
        e.special_total,
        2127.910749,
        SUITE,
        "DihedralRestraint vs check-suite table",
    );
}

#[test]
fn aladip_perturbed_distance_restraint_lambda_derivative_is_a_finite_difference() {
    // check_lambda_derivative on the special terms: (E(λ+ε) − E(λ−ε)) / 2ε, ε = δ = 0.001
    let sys = load("aladip_distres.in", true, true);
    let lambda = sys.imd.rlam;
    let eps = 0.001;
    let e = |l: f64| step0_energies(&sys, None, Some(l)).distanceres_total;
    let fd = (e(lambda + eps) - e(lambda - eps)) / (2.0 * eps);
    let analytic = step0_state(&sys, None, None)
        .perturbed_energy_derivatives
        .distanceres_total;
    assert!(
        (fd - analytic).abs() < 1e-3 * analytic.abs().max(1.0),
        "distance-restraint dH/dλ: finite difference {fd:.6} vs analytic {analytic:.6}"
    );
}

// ── lambdas.t.cc: individual λ, LAMBDAS on vs off, RLAM = 0.125 ──────────────────────────
//
// `lambdas.{on,off}.in` (byte-identical copies in `shared/check/`) set NSCALE = 1, and
// `lambdas.pttopo` carries a SCALEDINTERACTIONS block. The scaling only touches the nonbonded
// pairs it lists; newen does not implement it and refuses NSCALE ≠ 0, so both inputs are
// evaluated with NSCALE forced to 0. The bonded terms and the individual-λ loop below do not
// depend on it; the suite's NonBonded and NonBonded_newRF values do, and are left out. Its
// kinetic check compares two zeros (`check_forcefield` runs no temperature calculation); the mass
// L′ of the kinetic dH/dλ is judged by the `aladip_vacuum_fep_lambdas` reference instead.

fn load_lambdas(input: &str) -> System {
    let mut sys = load_with(input, Some("check/lambdas.pttopo"), false);
    sys.imd.nscale = 0;
    sys
}

#[test]
fn lambdas_bonded_terms_of_the_on_run_match_the_check_suite() {
    // check_forcefield on the "on" run: every perturbed term at its own individual λ
    let sys = load_lambdas("lambdas.on.in");
    let (topo, conf) = (&sys.prepared.topology, &sys.prepared.configuration);
    let unpert = calculate_bonded_forces_ntf(
        topo,
        conf,
        &vacuum(),
        true,
        newen_forces::bonded::CovalentForm::default(),
        true,
        true,
        true,
        true,
        &mut newen_forces::bonded::BondedGroupEnergies::default(),
    );
    let pert = calculate_perturbed_bonded_forces(
        topo,
        conf,
        &vacuum(),
        &mut newen_forces::bonded::BondedGroupEnergies::default(),
        &mut newen_forces::bonded::BondedGroupEnergies::default(),
        None,
    );
    let table = |actual: f64, expected: f64, what: &str| {
        close(
            actual,
            expected,
            SUITE,
            &format!("{what} vs check-suite table"),
        )
    };
    table(unpert.bond_energy, 18.055276, "QuarticBond");
    table(pert.bond_energy, 1.271962, "PerturbedQuarticBond");
    table(unpert.angle_energy, 12.170290, "Angle");
    table(pert.angle_energy, 0.500162, "PerturbedAngle");
    table(unpert.improper_energy, 0.965060, "ImproperDihedral");
    table(pert.improper_energy, 0.560988, "PerturbedImproperDihedral");
    table(unpert.dihedral_energy, 2.255206, "Dihedral");
    table(pert.dihedral_energy, 4.626301, "PerturbedDihedral");
}

#[test]
fn lambdas_individual_lambdas_equal_the_off_run_at_their_own_lambda() {
    // The individual-λ loop of lambdas.t.cc: for every interaction kind and energy group whose
    // polynomial is not the identity, the "off" run at RLAM = L gives the "on" run's per-group
    // energy, and dE_on = dE_off · L′ (δ = 1e-7). As there, a term whose dE_on is zero (no such
    // interaction in that group: the restraint rows, group 2) is not counted as checked.
    use newen_core::lambdas::{InteractionLambda as K, IDENTITY_COEFFICIENTS};
    let on = load_lambdas("lambdas.on.in");
    let off = load_lambdas("lambdas.off.in");
    let lam = on.imd.rlam;
    assert_eq!(lam, off.imd.rlam, "on and off runs are at different λ");
    let block = on.imd.lambdas.as_ref().expect("LAMBDAS block");
    assert!(block.enabled);
    assert!(!off.imd.lambdas.as_ref().expect("LAMBDAS block").enabled);
    let s_on = step0_state(&on, None, None);
    let n_groups = on.prepared.topology.energy_groups.len();
    assert_eq!(n_groups, 2);

    let mut checked = Vec::new();
    for kind in [K::Bond, K::Angle, K::Improper, K::Dihedral, K::Disres] {
        for n1 in 0..n_groups {
            // gromosXX's coefficient matrices: the identity, overwritten by the rows in order
            let [a, b, c, d, e] = block
                .rows
                .iter()
                .rfind(|r| r.kind == kind && r.g1 == n1 + 1 && r.g2 == n1 + 1)
                .map_or(IDENTITY_COEFFICIENTS, |r| r.coefficients);
            if [a, b, c, d, e] == IDENTITY_COEFFICIENTS {
                continue;
            }
            let l_int =
                a * lam * lam * lam * lam + b * lam * lam * lam + c * lam * lam + d * lam + e;
            let dl_int = 4.0 * a * lam * lam * lam + 3.0 * b * lam * lam + 2.0 * c * lam + d;
            let s_off = step0_state(&off, None, Some(l_int));
            let pick = |s: &newen_core::State| {
                let (e, d) = (&s.energies, &s.perturbed_energy_derivatives);
                match kind {
                    K::Bond => (e.bond_energy[n1], d.bond_energy[n1]),
                    K::Angle => (e.angle_energy[n1], d.angle_energy[n1]),
                    K::Improper => (e.improper_energy[n1], d.improper_energy[n1]),
                    K::Dihedral => (e.dihedral_energy[n1], d.dihedral_energy[n1]),
                    K::Disres => (e.distanceres_energy[n1], d.distanceres_energy[n1]),
                    _ => unreachable!(),
                }
            };
            let ((e_on, de_on), (e_off, de_off)) = (pick(&s_on), pick(&s_off));
            if de_on == 0.0 {
                continue;
            }
            let what = format!("{} (group {})", kind.gromos_name(), n1 + 1);
            close(e_on, e_off, 1e-7, &format!("E_on = E_off(L), {what}"));
            close(
                de_on,
                de_off * dl_int,
                1e-7,
                &format!("dE_on = dE_off(L)·L′, {what}"),
            );
            checked.push(what);
        }
    }
    // gromosXX checks the same four: every perturbed bonded term of lambdas.pttopo is in group 1
    assert_eq!(
        checked,
        [
            "bond (group 1)",
            "angle (group 1)",
            "improper (group 1)",
            "dihedral (group 1)"
        ],
        "individual-λ checks run"
    );
}
