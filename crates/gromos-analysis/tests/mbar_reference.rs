//! `mbar` has no gromos++ counterpart: it is checked against pymbar on the five CH4 windows of
//! the ext_ti_ana chain (`tests/data/mbar/make_reference.py` made the numbers), and against
//! `bar` for two states, where MBAR is BAR.

use gromos_analysis::programs::{bar, mbar};
use std::path::{Path, PathBuf};

fn chain(file: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/ext_ti_ana/chain")
        .join(file)
        .to_string_lossy()
        .into_owned()
}

fn reference() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mbar/pymbar_reference.dat")
}

const LAMS: [&str; 5] = ["0.00", "0.25", "0.50", "0.75", "1.00"];

#[test]
fn mbar_matches_pymbar_on_the_five_ch4_windows() {
    let files: Vec<_> = LAMS
        .iter()
        .map(|l| bar::read_bar_file(&chain(&format!("bar_data_l{l}.dat"))).unwrap())
        .collect();
    let kbt = gromos_analysis::GROMOSPP_KB * 300.0;
    let (lam, u_kn, n_k) = mbar::reduced_potentials(&files, kbt).unwrap();
    assert_eq!(lam, vec![0.0, 0.25, 0.5, 0.75, 1.0]);
    let r = mbar::mbar(&u_kn, &n_k, 10000, 1e-10).unwrap();
    let text = std::fs::read_to_string(reference()).unwrap();
    let mut pairs = 0;
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let t: Vec<&str> = line.split_whitespace().collect();
        let (i, j): (usize, usize) = (t[0].parse().unwrap(), t[1].parse().unwrap());
        let (dg, ddg): (f64, f64) = (t[2].parse().unwrap(), t[3].parse().unwrap());
        let ours = r.delta_f[i][j] * kbt;
        let ours_err = r.d_delta_f[i][j] * kbt;
        assert!(
            (ours - dg).abs() < 0.01,
            "ΔG {i}->{j}: {ours} vs pymbar {dg}"
        );
        assert!(
            (ours_err - ddg).abs() < 0.01,
            "δΔG {i}->{j}: {ours_err} vs pymbar {ddg}"
        );
        pairs += 1;
    }
    assert_eq!(pairs, 10);
}

#[test]
fn mbar_of_two_states_is_bar() {
    for pair in LAMS.windows(2) {
        let files: Vec<_> = pair
            .iter()
            .map(|l| bar::read_bar_file(&chain(&format!("bar_data_l{l}.dat"))).unwrap())
            .collect();
        let opts = bar::BarOptions {
            temp: 300.0,
            max_iterations: 1000,
            convergence: 1e-14,
            bootstrap: 0,
            seed: None,
            printdist: None,
        };
        let b = bar::bar(&files, &opts, &mut std::io::sink()).unwrap();
        let kbt = gromos_analysis::GROMOSPP_KB * 300.0;
        let (_, u_kn, n_k) = mbar::reduced_potentials(&files, kbt).unwrap();
        let m = mbar::mbar(&u_kn, &n_k, 100000, 1e-14).unwrap();
        let dg_mbar = m.delta_f[0][1] * kbt;
        assert!(
            (dg_mbar - b.intervals[0].dg).abs() < 1e-6,
            "{pair:?}: MBAR {dg_mbar} vs BAR {}",
            b.intervals[0].dg
        );
    }
}

#[test]
fn bar_bootstrap_is_reproducible_with_a_seed() {
    let files: Vec<_> = ["0.00", "0.25"]
        .iter()
        .map(|l| bar::read_bar_file(&chain(&format!("bar_data_l{l}.dat"))).unwrap())
        .collect();
    let opts = bar::BarOptions {
        temp: 300.0,
        max_iterations: 500,
        convergence: 1e-10,
        bootstrap: 20,
        seed: Some(12345),
        printdist: None,
    };
    let mut out_a = Vec::new();
    let mut out_b = Vec::new();
    let a = bar::bar(&files, &opts, &mut out_a).unwrap();
    let b = bar::bar(&files, &opts, &mut out_b).unwrap();
    assert_eq!(out_a, out_b);
    let boot = a.intervals[0].bootstrap.unwrap();
    assert!(boot.is_finite() && boot > 0.0);
    assert_eq!(a, b);
    // the corrected analytical error leaves out gromos++'s extra 1/(ln<f>)² factor
    assert!(a.intervals[0].error_corrected != a.intervals[0].error);
}
