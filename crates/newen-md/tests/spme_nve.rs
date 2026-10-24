//! SPME (NLRELE 5, newen only, experimental) in a run: 200 NVE steps of the 216 waters, whose
//! total-energy drift must stay below twice that of the same run with P3M (NLRELE 3). gromosXX
//! has no SPME, so there is no gromosXX reference; the run is compared with newen's own P3M.

use std::path::{Path, PathBuf};
use std::process::Command;

fn refs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("newen_spme_nve_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Run `md` on the 216 waters with `imd` (text) from `conf`; returns its standard output.
fn md(dir: &Path, imd: &str, conf: &Path) -> String {
    let input = dir.join("run.imd");
    std::fs::write(&input, imd).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_md"))
        .arg("@topo")
        .arg(refs().join("water_216_box/water_216_box.topo"))
        .arg("@conf")
        .arg(conf)
        .arg("@input")
        .arg(&input)
        .arg("@fin")
        .arg(dir.join("final.cnf"))
        .arg("@tre")
        .arg(dir.join("ener.tre"))
        .arg("@trc")
        .arg(dir.join("traj.trc"))
        .output()
        .unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout).to_string();
    assert!(
        out.status.success(),
        "{stdout}{}",
        String::from_utf8_lossy(&out.stderr)
    );
    stdout
}

/// The first number of every ENERGY03 block: the total energy per step.
fn totals(path: &Path) -> Vec<f64> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = Vec::new();
    let mut inside = false;
    for l in text.lines().map(str::trim) {
        if l == "ENERGY03" {
            inside = true;
        } else if inside && !l.starts_with('#') {
            out.push(l.split_whitespace().next().unwrap().parse::<f64>().unwrap());
            inside = false;
        }
    }
    out
}

/// The P3M reference input with `NLRELE`, `steps` steps, the energies every step, velocities
/// from the configuration (NTIVEL 0) and the bath coupled with `tau` (−1: not coupled, NVE).
fn input(nlrele: i32, steps: usize, tau: f64) -> String {
    let base = std::fs::read_to_string(
        refs().join("water_216_box_p3m_gauss_fixed/water_216_box_p3m_gauss_fixed.in"),
    )
    .unwrap();
    let base = base.replace(
        "#   NLRELE\n         3\n",
        &format!("#   NLRELE\n         {nlrele}\n"),
    );
    assert!(base.contains(&format!("         {nlrele}\n#    APPAK")));
    let base = base.replace(
        "       10       0.0     0.002",
        &format!("{steps:>9}       0.0     0.002"),
    );
    assert!(base.contains("    300.0  -1.0\n"));
    base.replace("    300.0  -1.0\n", &format!("    300.0  {tau}\n"))
}

/// The slope of the least-squares line through `e` against the step (kJ mol⁻¹ per step), and
/// the largest deviation from the first value.
fn drift(e: &[f64]) -> (f64, f64) {
    let n = e.len() as f64;
    let mx = (n - 1.0) / 2.0;
    let my = e.iter().sum::<f64>() / n;
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (i, y) in e.iter().enumerate() {
        let dx = i as f64 - mx;
        sxy += dx * (y - my);
        sxx += dx * dx;
    }
    let dev = e.iter().fold(0.0f64, |m, y| m.max((y - e[0]).abs()));
    (sxy / sxx, dev)
}

#[test]
fn spme_conserves_energy_at_least_as_well_as_p3m_within_a_factor_two() {
    // the generated box relaxes for hundreds of steps: equilibrate first (P3M, weak coupling
    // to 300 K, τ = 0.1 ps), then compare the two methods from the same configuration
    let eq = scratch("equilibrate");
    md(
        &eq,
        &input(3, 1000, 0.1),
        &refs().join("water_216_box/water_216_box.conf"),
    );
    let start = eq.join("final.cnf");
    let dir = scratch("spme");
    let out = md(&dir, &input(5, 200, -1.0), &start);
    assert!(
        out.contains(
            "SPME is experimental: no gromosXX oracle; checked against the Lean Ewald model and P3M"
        ),
        "{out}"
    );
    let spme = totals(&dir.join("ener.tre"));
    let dir_p = scratch("p3m");
    let out_p = md(&dir_p, &input(3, 200, -1.0), &start);
    assert!(!out_p.contains("SPME is experimental"));
    let p3m = totals(&dir_p.join("ener.tre"));
    assert_eq!(spme.len(), 200);
    assert_eq!(p3m.len(), 200);
    let (s_spme, dev_spme) = drift(&spme);
    let (s_p3m, dev_p3m) = drift(&p3m);
    let ps = 0.002;
    println!(
        "200 NVE steps of 216 SPC waters (dt 2 fs): total-energy drift SPME (NSPORD 4, 32³) \
         {:.3} kJ/mol/ps (largest deviation {dev_spme:.3} kJ/mol), P3M (NASORD 3, ik, 32³) \
         {:.3} kJ/mol/ps (largest deviation {dev_p3m:.3} kJ/mol); E_tot(0) SPME {:.4}, P3M {:.4}",
        s_spme / ps,
        s_p3m / ps,
        spme[0],
        p3m[0]
    );
    assert!(
        s_spme.abs() < 2.0 * s_p3m.abs(),
        "SPME drift {s_spme:e} vs P3M {s_p3m:e} per step"
    );
}
