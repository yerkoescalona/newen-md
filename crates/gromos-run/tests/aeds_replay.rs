//! A-EDS search replay (1.8): the pure update functions of `gromos_integrators::eds`, fed with the
//! end-state energies V_i(t) gromosXX printed for the `aladip_vacuum_aeds_f*` references, must
//! give the V_mix and V_R* it printed for every frame (computed with the parameters of the frame
//! before: the columns are written after the update) and the E_max, E_min, G, σ_G and offsets
//! it printed after the update, and end in its final AEDSSEARCH block. This pins the search
//! bookkeeping independently of the MD. The `_fixed` references (gromosXX with the D16, D25
//! and D32 patches) replay under the intended rules, their columns before the update.

use std::path::{Path, PathBuf};

use gromos_integrators::eds::{
    boosted_reference, initialise_search, reference_hamiltonian, search_step, AedsColumns,
    AedsForm, AedsParams, AedsPrevalence, AedsSearch, AedsTarget, EdsForm,
};
use gromos_io::coordinate::{parse_aedssearch, read_coordinates};
use gromos_io::imd::read_imd_file;
use gromos_io::topology::read_topology_file;

fn refs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
}

/// One printed frame: ENER[1..] and the `# eds` rows.
struct Frame {
    ener: Vec<f64>,
    eds: Vec<Vec<f64>>,
}

impl Frame {
    /// ENER[k], 1-based as in the ene_ana library.
    fn e(&self, k: usize) -> f64 {
        self.ener[k - 1]
    }
}

fn frames(path: &Path) -> Vec<Frame> {
    let text = std::fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = text.lines().map(str::trim).collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        if lines[i] == "ENERGY03" {
            let mut f = Frame {
                ener: Vec::new(),
                eds: Vec::new(),
            };
            i += 1;
            while lines[i] != "END" {
                if lines[i] == "# totals" {
                    i += 1;
                    while !lines[i].starts_with('#') {
                        f.ener.push(lines[i].parse().unwrap());
                        i += 1;
                    }
                    continue;
                }
                if lines[i] == "# eds" {
                    let n: usize = lines[i + 2].parse().unwrap();
                    for k in 0..n {
                        f.eds.push(
                            lines[i + 4 + k]
                                .split_whitespace()
                                .map(|x| x.parse().unwrap())
                                .collect(),
                        );
                    }
                    i += 4 + n;
                    continue;
                }
                i += 1;
            }
            out.push(f);
        }
        i += 1;
    }
    out
}

/// Relative difference with an absolute floor at the printed precision.
fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs().max(1.0)
}

/// gromosXX 247c95c4 as released.
const RELEASED: (AedsColumns, AedsPrevalence, AedsTarget) = (
    AedsColumns::GromosxxPostUpdate,
    AedsPrevalence::GromosxxPlainExp,
    AedsTarget::GromosxxState1Default,
);
/// gromosXX with the D16, D25 and D32 patches: gromos-rs's defaults.
const PATCHED: (AedsColumns, AedsPrevalence, AedsTarget) = (
    AedsColumns::Intended,
    AedsPrevalence::Intended,
    AedsTarget::Intended,
);

/// Replay one reference under `rules`; returns the largest relative deviation seen.
fn replay(name: &str, rules: (AedsColumns, AedsPrevalence, AedsTarget)) -> f64 {
    let dir = refs().join(name);
    let imd = read_imd_file(dir.join(format!("{name}.in"))).unwrap();
    let a = imd.aeds.clone().expect("AEDS block");
    let kb = read_topology_file(refs().join("shared/aladip.topo"))
        .unwrap()
        .physical_constants
        .kB;
    let beta = 1.0 / (kb * imd.temp_bath[0].temp0[0]);
    let params = AedsParams {
        form: AedsForm::from_gromos(a.form).unwrap(),
        bmax_type: a.bmaxtype as u32,
        bmax: a.bmax,
        asteps: a.asteps,
        bsteps: a.bsteps,
        columns: rules.0,
        prevalence: rules.1,
        target: rules.2,
    };
    let mut s = AedsSearch::new(a.emax, a.emin, a.eir.clone(), !a.restremin);
    let printed = frames(&dir.join("expected/energies.tre"));
    assert_eq!(printed.len(), imd.nstlim);
    let n = a.num_states;
    let (mut w, mut scratch) = (vec![0.0; n], vec![0.0; n]);
    // the printed values have ten significant digits; the replay feeds them back in, so the
    // memory averages carry that round-off along
    let tol = 1e-9;
    let mut worst = 0.0f64;
    let mut check = |what: &str, k: usize, ours: f64, theirs: f64| {
        worst = worst.max((ours - theirs).abs() / theirs.abs().max(1.0));
        assert!(
            close(ours, theirs, tol),
            "{name} frame {k}: {what} {ours} vs gromosXX {theirs}"
        );
    };
    for (k, f) in printed.iter().enumerate() {
        let vi: Vec<f64> = f.eds.iter().map(|r| r[0]).collect();
        let v_mix = reference_hamiltonian(&EdsForm::SingleS { s: 1.0 }, beta, &s.eir, &vi, &mut w);
        if s.initialise {
            initialise_search(params.form, &mut s, v_mix, beta);
        }
        let (v_r, _) = boosted_reference(v_mix, s.emax, s.emin);
        check("eds_vmix", k, v_mix, f.e(34));
        check("eds_vr", k, v_r, f.e(35));
        let before = (s.emax, s.emin, s.eir.clone());
        let globmin = search_step(
            &params,
            &mut s,
            &vi,
            v_mix,
            v_r,
            beta,
            k,
            imd.nstlim,
            imd.dt,
            &mut scratch,
        )
        .unwrap();
        s.initialise = false;
        let (emax, emin, eir) = match params.columns {
            AedsColumns::Intended => before,
            AedsColumns::GromosxxPostUpdate => (s.emax, s.emin, s.eir.clone()),
        };
        check("eds_emax", k, emax, f.e(36));
        check("eds_emin", k, emin, f.e(37));
        if let Some((g, fluc)) = globmin {
            check("eds_globmin", k, g, f.e(38));
            check("eds_globminfluc", k, fluc, f.e(39));
        }
        for (i, r) in f.eds.iter().enumerate() {
            check(&format!("e{}r", i + 1), k, eir[i], r[3]);
        }
    }
    let fin = read_coordinates(dir.join("expected/final.conf")).unwrap();
    let theirs = parse_aedssearch(fin.aedssearch.as_ref().unwrap(), n).unwrap();
    assert_eq!(
        (s.emaxcounts, s.oldstate, s.fullemin),
        (theirs.emaxcounts, theirs.oldstate, theirs.fullemin),
        "{name}: AEDSSEARCH counters"
    );
    assert_eq!(s.visited, theirs.visited, "{name}");
    assert_eq!(s.visitcounts, theirs.visitcounts, "{name}");
    assert_eq!(s.framecounts, theirs.framecounts, "{name}");
    for (what, ours, theirs) in [
        ("emax", vec![s.emax], vec![theirs.emax]),
        ("emin", vec![s.emin], vec![theirs.emin]),
        ("searchemax", vec![s.searchemax], vec![theirs.searchemax]),
        ("eir", s.eir.clone(), theirs.eir.clone()),
        ("lnexpde", s.lnexpde.clone(), theirs.lnexpde.clone()),
        ("statefren", s.statefren.clone(), theirs.statefren.clone()),
        ("avgenergy", s.avgenergy.clone(), theirs.avgenergy.clone()),
        (
            "eiravgenergy",
            s.eiravgenergy.clone(),
            theirs.eiravgenergy.clone(),
        ),
        ("bigs", s.bigs.clone(), theirs.bigs.clone()),
        (
            "stdevenergy",
            s.stdevenergy.clone(),
            theirs.stdevenergy.clone(),
        ),
    ] {
        for (a, b) in ours.iter().zip(&theirs) {
            assert!(
                close(*a, *b, tol),
                "{name}: AEDSSEARCH {what} {a} vs gromosXX {b}"
            );
        }
    }
    worst
}

#[test]
fn form2_offset_search_replays_gromosxx() {
    eprintln!("f2 worst {:.1e}", replay("aladip_vacuum_aeds_f2", RELEASED));
}

#[test]
fn form4_full_search_replays_gromosxx() {
    eprintln!("f4 worst {:.1e}", replay("aladip_vacuum_aeds_f4", RELEASED));
}

#[test]
fn form5_advanced_search_replays_gromosxx() {
    eprintln!("f5 worst {:.1e}", replay("aladip_vacuum_aeds_f5", RELEASED));
}

#[test]
fn form6_advanced_offsets_replay_gromosxx() {
    eprintln!("f6 worst {:.1e}", replay("aladip_vacuum_aeds_f6", RELEASED));
}

#[test]
fn form3_zero_variance_clamp_replays_gromosxx() {
    eprintln!(
        "sigma0 worst {:.1e}",
        replay("aladip_vacuum_aeds_sigma0", RELEASED)
    );
}

#[test]
fn form6_overflow_freezes_the_released_offsets() {
    eprintln!(
        "f6_ovf worst {:.1e}",
        replay("aladip_vacuum_aeds_f6_ovf", RELEASED)
    );
}

/// The patched oracle (D16, D25, D32) under gromos-rs's default rules: the columns are the ones
/// the frame's V_mix and V_R* used, so the replay's pre-update values must match.
#[test]
fn patched_searches_replay_under_the_intended_rules() {
    for name in [
        "aladip_vacuum_aeds_f2_fixed",
        "aladip_vacuum_aeds_f4_fixed",
        "aladip_vacuum_aeds_f5_fixed",
        "aladip_vacuum_aeds_f6_fixed",
        "aladip_vacuum_aeds_f6_ovf_fixed",
        "aladip_vacuum_aeds_sigma0_fixed",
    ] {
        eprintln!("{name} worst {:.1e}", replay(name, PATCHED));
    }
}
