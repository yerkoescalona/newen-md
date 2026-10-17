//! The real-space ingredients of gromos++'s FFT Poisson solver: the solute occupancy
//! (`pb::FFTInsideOutside`), the hat charge-shaping function (`FFTChargeShapingFunction_Hat`) and
//! the two vacuum fields (`FFTVacuumField_LS`, Ewald-type in k-space, and `FFTVacuumField_RF`,
//! reaction-field-type in real space).
//!
//! Every expression keeps gromos++'s operation order: the cubelet tests are strict comparisons a
//! one-ulp change can flip. The grid-point loops run in parallel over x slices; each point still
//! sees the atoms in gromos++'s order, so the results do not depend on the thread count.

use super::fft::{k_component, FftGrid};
use super::fft_poisson::StreamLog;
use super::PbSystem;
use gromos_core::Vec3;
use rayon::prelude::*;
use rustfft::num_complex::Complex64;

/// gromos++ `PB_Parameters::tiny_real`.
pub(crate) const TINY: f64 = 0.0000000001;

/// Minimum image of `cellori − pos` on the grid's own period `n·d` (the expressions of
/// `FFTInsideOutside.cc:139-156` and `FFTVacuumField_RF.cc:213-235`).
#[inline]
fn nearest(cellori: [f64; 3], pos: Vec3, n: [usize; 3], d: [f64; 3]) -> [f64; 3] {
    let mut nim = [cellori[0] - pos.x, cellori[1] - pos.y, cellori[2] - pos.z];
    for u in 0..3 {
        let boxh = n[u] as f64 * d[u] / 2.0;
        let bx = n[u] as f64 * d[u];
        if nim[u] > boxh {
            nim[u] -= bx;
        } else if nim[u] < -boxh {
            nim[u] += bx;
        }
    }
    nim
}

/// `FFTInsideOutside::inside_sol` (`FFTInsideOutside.cc:46-267`): the fraction of each grid cell's
/// `(2·nc+1)³` cubelets that lie inside any atom of `atoms` (radius `radius[a]`, zero radii
/// skipped), counted as a union over atoms; values within `TINY` of 1 become exactly 1.
pub fn inside_sol(
    grid: &FftGrid,
    nc: usize,
    positions: &[Vec3],
    atoms: &[usize],
    radius: &[f64],
) -> Vec<f64> {
    let gridn = grid.n;
    let gridd = grid.dr;
    let dr = 0.5 * (gridd[0] * gridd[0] + gridd[1] * gridd[1] + gridd[2] * gridd[2]).sqrt();
    let nc = nc as i64;
    let ce = 2 * nc + 1;
    let cei = 1.0 / ce as f64;
    let ddx = gridd[0] * cei;
    let ddy = gridd[1] * cei;
    let ddz = gridd[2] * cei;
    let ctoi = cei * cei * cei;
    let cecece = (ce * ce * ce) as usize;

    let mut inside = vec![0.0; grid.ngr3];
    let plane = gridn[1] * gridn[2];
    inside
        .par_chunks_mut(plane)
        .enumerate()
        .for_each(|(ii, slice)| {
            let mut cear = vec![false; cecece];
            for jj in 0..gridn[1] {
                for kk in 0..gridn[2] {
                    cear.fill(false);
                    let mut cin = 0usize;
                    let cellori = [
                        ii as f64 * gridd[0],
                        jj as f64 * gridd[1],
                        kk as f64 * gridd[2],
                    ];
                    for &a in atoms {
                        let rr = radius[a];
                        if rr.abs() < TINY {
                            continue;
                        }
                        let rplus = rr + dr;
                        let rminus = if rr > dr { rr - dr } else { 0.0 };
                        let r2plus = rplus * rplus;
                        let r2minus = rminus * rminus;
                        let r2 = rr * rr;
                        let nim = nearest(cellori, positions[a], gridn, gridd);
                        let distance2 = nim[0] * nim[0] + nim[1] * nim[1] + nim[2] * nim[2];
                        if distance2 >= r2plus + TINY {
                            continue;
                        }
                        if (distance2 < r2minus + TINY && distance2 > TINY) || distance2 <= TINY {
                            cin = cecece;
                            cear.fill(true);
                        } else {
                            let mut c = 0usize;
                            for xx in -nc..=nc {
                                let dx = nim[0] + xx as f64 * ddx;
                                let dx2 = dx * dx;
                                for yy in -nc..=nc {
                                    let dy = nim[1] + yy as f64 * ddy;
                                    let dy2 = dy * dy;
                                    for zz in -nc..=nc {
                                        let dz = nim[2] + zz as f64 * ddz;
                                        let dz2 = dz * dz;
                                        let dist2c = dx2 + dy2 + dz2;
                                        if dist2c < r2 + TINY && !cear[c] {
                                            cin += 1;
                                            cear[c] = true;
                                        }
                                        c += 1;
                                    }
                                }
                            }
                        }
                    }
                    let mut v = cin as f64 * ctoi;
                    if v > 1.0 || (v - 1.0).abs() < TINY {
                        v = 1.0;
                    }
                    slice[jj * gridn[2] + kk] = v;
                }
            }
        });
    inside
}

/// `FFTInsideOutside::integr_inside` (`:273-298`): every value below `1 − TINY` becomes 0; logs
/// the count of points left and the integral of the fractional values.
pub fn integr_inside(inside: &mut [f64], log: &mut StreamLog) {
    let mut integral = 0.0;
    let mut in_points = inside.len() as i64;
    for v in inside.iter_mut() {
        integral += *v;
        if *v < 1.0 - TINY {
            *v = 0.0;
            in_points -= 1;
        }
    }
    let n = inside.len();
    log.line(format!(
        "# FFTInsideOutside::integr_inside: gridpoints: {n}"
    ));
    log.line(format!(
        "# FFTInsideOutside::integr_inside: inPoints_noB: {in_points}"
    ));
    let f = log.fmt;
    log.line(format!(
        "# FFTInsideOutside::integr_inside: integral {}",
        f.f(integral)
    ));
    log.line(format!(
        "# FFTInsideOutside::integr_inside: normalized integral {}",
        f.f(integral / (1.0 * n as f64))
    ));
}

/// The k-space potential of a unit hat charge of radius `alpha` plus its neutralising background,
/// `1/(ε₀k²)·12(2 − 2cos αk − αk sin αk)/(αk)⁴`, and 0 for `k² < TINY`
/// (`FFTChargeShapingFunction_Hat.cc:50-56`).
#[inline]
pub fn hat(x: f64, y: f64, z: f64, alpha: f64, eps0: f64) -> f64 {
    let k2 = x * x + y * y + z * z;
    if k2 >= TINY {
        let alphak = alpha * k2.sqrt();
        1.0 / (eps0 * k2) * 12.0 * (2.0 - 2.0 * alphak.cos() - alphak * alphak.sin())
            / (alphak * alphak * alphak * alphak)
    } else {
        0.0
    }
}

/// The charged atoms of `@atoms` split by radius (`FFTVacuumField_LS::makeatomlist`,
/// `FFTVacuumField_LS.cc:129-181`): above `alpha1 + TINY` they are shaped with `alpha1`, the
/// others with `alpha2`. Entries are positions in `sys.atoms`.
pub struct IonLists {
    pub big: Vec<usize>,
    pub small: Vec<usize>,
}

impl IonLists {
    pub fn new(
        sys: &PbSystem,
        alpha1: f64,
        alpha2: f64,
        log: &mut StreamLog,
    ) -> Result<Self, String> {
        let mut big = Vec::new();
        let mut small = Vec::new();
        for (ion, &a) in sys.atoms.iter().enumerate() {
            let r = sys.radii[a];
            let charged = sys.charges[a].abs() > TINY;
            if r > alpha1 + TINY {
                if charged {
                    big.push(ion);
                }
            } else if charged {
                if r < alpha2 {
                    let f = log.fmt;
                    log.line(format!("Ion {ion} : radius {}", f.f(r)));
                    return Err("FFTVacuumField_LS: Ion too small. Exiting.".into());
                }
                small.push(ion);
            }
        }
        let f = log.fmt;
        for (n, list) in [(1, &big), (2, &small)] {
            log.line(format!("# IONLIST {n}: NUMBER OF ENTRIES: {}", list.len()));
            for (e, &ion) in list.iter().enumerate() {
                log.line(format!(
                    "# list {n} : entry {e} has ion_list{n}[ {e} ] = {ion} and radius {}",
                    f.f(sys.radii[sys.atoms[ion]])
                ));
            }
        }
        Ok(IonLists { big, small })
    }
}

/// `FFTVacuumField_LS::positionIndependentVacField` (`:243-304`): per grid k-vector,
/// `F_c(k) = Σ_n −(k_c + n_c·ka_c)·ψ(k + n∘ka)` over the aliases `n ∈ [−N, N]³`.
fn position_independent(
    grid: &FftGrid,
    ka: [f64; 3],
    nalias: i64,
    alpha: f64,
    eps0: f64,
) -> [Vec<f64>; 3] {
    let [nx, ny, nz] = grid.n;
    let plane = ny * nz;
    let slices: Vec<[Vec<f64>; 3]> = (0..nx)
        .into_par_iter()
        .map(|i| {
            let mut out = [vec![0.0; plane], vec![0.0; plane], vec![0.0; plane]];
            let kx = k_component(i, nx, grid.dk[0]);
            for j in 0..ny {
                let ky = k_component(j, ny, grid.dk[1]);
                for k in 0..nz {
                    let kz = k_component(k, nz, grid.dk[2]);
                    let mut fld = [0.0f64; 3];
                    for iax in -nalias..=nalias {
                        for iay in -nalias..=nalias {
                            for iaz in -nalias..=nalias {
                                let ax = kx + iax as f64 * ka[0];
                                let ay = ky + iay as f64 * ka[1];
                                let az = kz + iaz as f64 * ka[2];
                                let psi = hat(ax, ay, az, alpha, eps0);
                                fld[0] += -ax * psi;
                                fld[1] += -ay * psi;
                                fld[2] += -az * psi;
                            }
                        }
                    }
                    let m = k + nz * j;
                    for c in 0..3 {
                        out[c][m] = fld[c];
                    }
                }
            }
            out
        })
        .collect();
    let mut dest = [
        Vec::with_capacity(grid.ngr3),
        Vec::with_capacity(grid.ngr3),
        Vec::with_capacity(grid.ngr3),
    ];
    for s in slices {
        for c in 0..3 {
            dest[c].extend_from_slice(&s[c]);
        }
    }
    dest
}

/// `updateMultipliers` (`:362-394`): `Σ q sin(k·r)` (real slot) and `Σ q cos(k·r)` (imaginary
/// slot) over the list, in list order.
fn multipliers(
    grid: &FftGrid,
    sys: &PbSystem,
    positions: &[Vec3],
    list: &[usize],
) -> Vec<Complex64> {
    let [nx, ny, nz] = grid.n;
    let plane = ny * nz;
    let mut mult = vec![Complex64::default(); grid.ngr3];
    mult.par_chunks_mut(plane).enumerate().for_each(|(i, s)| {
        let kx = k_component(i, nx, grid.dk[0]);
        for j in 0..ny {
            let ky = k_component(j, ny, grid.dk[1]);
            for k in 0..nz {
                let kz = k_component(k, nz, grid.dk[2]);
                let mut sum_im = 0.0;
                let mut sum_re = 0.0;
                for &ion in list {
                    let a = sys.atoms[ion];
                    let p = positions[a];
                    let kr = kx * p.x + ky * p.y + kz * p.z;
                    let charge = sys.charges[a];
                    sum_im += charge * kr.cos();
                    sum_re += charge * kr.sin();
                }
                s[k + nz * j] = Complex64::new(sum_re, sum_im);
            }
        }
    });
    mult
}

/// `FFTVacuumField_LS::calcVacField` (`:186-238`) with `recyclefield`: the Ewald-type vacuum field
/// of the charged `@atoms` in k-space, `−i k ψ(k) Σ q e^{−ik·r}` (aliased), into `v` (zeroed).
#[allow(clippy::too_many_arguments)]
pub fn ls_vacuum_field(
    grid: &FftGrid,
    sys: &PbSystem,
    positions: &[Vec3],
    lists: &IonLists,
    alpha: [f64; 2],
    nalias: [i64; 2],
    eps0: f64,
    pi: f64,
    v: &mut [Vec<Complex64>; 3],
    log: &mut StreamLog,
) {
    log.line("# FFTVacuumField_LS::calcVacField :  Calculating vac field from scratch ...");
    let pi2 = 2.0 * pi;
    let ka = [pi2 / grid.dr[0], pi2 / grid.dr[1], pi2 / grid.dr[2]];
    let zero = || {
        [
            vec![0.0; grid.ngr3],
            vec![0.0; grid.ngr3],
            vec![0.0; grid.ngr3],
        ]
    };
    let slices = |log: &mut StreamLog| {
        for i in 0..grid.n[0] {
            log.line(format!(
                "# FFTVacuumField_LS::positionIndependentVacField : slice {i}"
            ));
        }
    };

    let store1 = if lists.big.is_empty() {
        zero()
    } else {
        let f = log.fmt;
        log.line(format!(
            "# Computing position-independent part of the vacuum field due to big ions, with alpha {}",
            f.f(alpha[0])
        ));
        slices(log);
        position_independent(grid, ka, nalias[0], alpha[0], eps0)
    };
    let f = log.fmt;
    for m in [0, 100] {
        for (c, name) in ["Vx1_store", "Vy1_store", "Vz1_store"].iter().enumerate() {
            log.line(format!(
                "# after positionIndependentVacField: {name}[{m}] {}",
                f.f(store1[c][m])
            ));
        }
    }
    let store2 = if lists.small.is_empty() {
        zero()
    } else {
        log.line(format!(
            "# Computing position-independent part of the vacuum field due to small ions, with alpha {}",
            f.f(alpha[1])
        ));
        slices(log);
        position_independent(grid, ka, nalias[1], alpha[1], eps0)
    };

    log.line("# FFTVacuumField_LS::recyclefield: Recycling vacuum field");
    for c in v.iter_mut() {
        c.fill(Complex64::default());
    }
    // complexFromDouble puts the store into both slots; the products are summed per slot.
    for (store, list) in [(&store1, &lists.big), (&store2, &lists.small)] {
        let mult = multipliers(grid, sys, positions, list);
        for c in 0..3 {
            for (l, d) in v[c].iter_mut().enumerate() {
                let s = store[c][l];
                d.re += s * mult[l].re;
                d.im += s * mult[l].im;
            }
        }
    }
    for m in [0, 100] {
        for (c, name) in ["fldx_k", "fldy_k", "fldz_k"].iter().enumerate() {
            log.line(format!(
                "# after recyclefield: {name}[{m}] {}",
                f.f(super::fft::slot(&v[c], m))
            ));
        }
    }
}

/// `FFTVacuumField_RF::calcVacFieldVincent` (`FFTVacuumField_RF.cc:61-307`): the real-space field
/// of the charged `@atoms` with a reaction-field correction (`srf = 2(ε−1)/(2ε+1)`, 1 for ε = 0)
/// inside the cutoff `rcut`, the cutoff sphere resolved in cubelets at its surface. A grid point
/// on top of an atom has its field reset to zero, dropping earlier atoms' contributions there too.
pub fn rf_vacuum_field(
    grid: &FftGrid,
    sys: &PbSystem,
    positions: &[Vec3],
    eps_rf: f64,
    rcut: f64,
    v: &mut [Vec<Complex64>; 3],
    log: &mut StreamLog,
) -> Result<(), String> {
    let f = log.fmt;
    log.line(format!(
        "# computing vacuum field with epsilon {}",
        f.f(eps_rf)
    ));
    log.line(
        "# Warning: The vacuum field on a grid point hosting a particle will be set to zero ...",
    );
    let rcut_3 = 1.0 / (rcut * rcut * rcut);
    let srf = if eps_rf.abs() < TINY {
        1.0
    } else {
        2.0 * (eps_rf - 1.0) / (2.0 * eps_rf + 1.0)
    };
    let frfs = srf * rcut_3;
    let [drx, dry, drz] = grid.dr;
    let dr = 0.5 * (drx * drx + dry * dry + drz * drz).sqrt();
    if rcut < dr {
        return Err(
            "FFTVacuumField_RF: Cutoff too small for this grid (cutoff<dr). Exiting.".into(),
        );
    }
    if rcut > grid.len[0] || rcut > grid.len[1] || rcut > grid.len[2] {
        return Err(
            "FFTVacuumField_RF: Cutoff too large for this grid (cutoff>gridlength). Exiting."
                .into(),
        );
    }
    let nc = grid.ncubes as i64;
    let ddx = drx / (2.0 * nc as f64 + 1.0);
    let ddy = dry / (2.0 * nc as f64 + 1.0);
    let ddz = drz / (2.0 * nc as f64 + 1.0);
    let ce = 2 * nc + 1;
    let cto = ce * ce * ce;
    let inverse_cubelets = 1.0 / cto as f64;
    let rplus = rcut + dr;
    let rminus = rcut - dr;
    let r2plus = rplus * rplus;
    let r2minus = rminus * rminus;
    let r2 = rcut * rcut;
    let scaled: Vec<f64> = sys
        .atoms
        .iter()
        .map(|&a| sys.fpepsi * sys.charges[a])
        .collect();

    let [nx, ny, nz] = grid.n;
    let plane = ny * nz;
    let gridd = grid.dr;
    let field: Vec<[f64; 3]> = (0..nx)
        .into_par_iter()
        .flat_map_iter(|ii| {
            let mut out = vec![[0.0f64; 3]; plane];
            for jj in 0..ny {
                for kk in 0..nz {
                    let fld = &mut out[kk + nz * jj];
                    let cellori = [
                        ii as f64 * gridd[0],
                        jj as f64 * gridd[1],
                        kk as f64 * gridd[2],
                    ];
                    for (aa, &a) in sys.atoms.iter().enumerate() {
                        if scaled[aa].abs() < TINY {
                            continue;
                        }
                        let nim = nearest(cellori, positions[a], grid.n, gridd);
                        let distance2 = nim[0] * nim[0] + nim[1] * nim[1] + nim[2] * nim[2];
                        if distance2 >= r2plus + TINY {
                            continue;
                        }
                        if distance2 <= TINY {
                            *fld = [0.0; 3];
                            continue;
                        }
                        let dist3i = 1.0 / (distance2 * distance2.sqrt());
                        let dfact = scaled[aa] * (dist3i - frfs);
                        let overlap = if distance2 < r2minus + TINY && distance2 > TINY {
                            1.0
                        } else {
                            let mut cin = 0i64;
                            for iii in -nc..=nc {
                                for jjj in -nc..=nc {
                                    for kkk in -nc..=nc {
                                        let dx = nim[0] + iii as f64 * ddx;
                                        let dy = nim[1] + jjj as f64 * ddy;
                                        let dz = nim[2] + kkk as f64 * ddz;
                                        if dx * dx + dy * dy + dz * dz < r2 + TINY {
                                            cin += 1;
                                        }
                                    }
                                }
                            }
                            cin as f64 * inverse_cubelets
                        };
                        let efact = dfact * overlap;
                        for c in 0..3 {
                            fld[c] += efact * nim[c];
                        }
                    }
                }
            }
            out
        })
        .collect();
    for c in 0..3 {
        for (d, s) in v[c].iter_mut().zip(&field) {
            *d = Complex64::new(s[c], 0.0);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hat_matches_the_formula() {
        let eps0 = 1.0 / (138.935484611 * 4.0 * std::f64::consts::PI);
        for &(x, y, z, alpha) in &[
            (1.0f64, 0.0f64, 0.0f64, 0.15f64),
            (3.0, -4.0, 12.0, 0.05),
            (0.2, 0.1, 0.0, 0.15),
        ] {
            let k = (x * x + y * y + z * z).sqrt();
            let ak = alpha * k;
            let by_hand =
                12.0 * (2.0 - 2.0 * ak.cos() - ak * ak.sin()) / ak.powi(4) / (eps0 * k * k);
            let got = hat(x, y, z, alpha, eps0);
            assert!(
                (got - by_hand).abs() <= 1e-12 * by_hand.abs(),
                "{got} {by_hand}"
            );
        }
        // Small αk: ψ·ε₀k² = 1 − (αk)²/15 + O((αk)⁴), the point charge in the limit. (gromos++'s
        // expression cancels catastrophically for αk ≲ 1e-3; the grids never get there.)
        let (k, alpha) = (0.1 / 0.15, 0.15);
        let got = hat(k, 0.0, 0.0, alpha, eps0) * eps0 * k * k;
        assert!((got - (1.0 - 0.01 / 15.0)).abs() < 1e-6, "{got}");
        // Below the k² guard it is zero.
        assert_eq!(hat(1e-6, 0.0, 0.0, 0.15, eps0), 0.0);
    }

    fn single_atom(pos: Vec3, r: f64) -> PbSystem {
        PbSystem {
            positions: vec![pos],
            charges: vec![1.0],
            radii: vec![r],
            atoms: vec![0],
            readout: vec![0],
            fpepsi: 138.935484611,
        }
    }

    #[test]
    fn inside_sol_on_a_single_atom() {
        // 10³ grid of 0.1 nm, one atom of radius 0.28 nm on grid point (5,5,5); dr = 0.0866.
        let g = FftGrid::new([10; 3], [1.0; 3], 1, std::f64::consts::PI);
        let sys = single_atom(Vec3::new(0.5, 0.5, 0.5), 0.28);
        let inside = inside_sol(&g, 1, &sys.positions, &sys.atoms, &sys.radii);
        // The atom's own point and its neighbours lie deep inside: full cells.
        assert_eq!(inside[g.index(5, 5, 5)], 1.0);
        assert_eq!(inside[g.index(6, 5, 5)], 1.0);
        // Far away: empty.
        assert_eq!(inside[g.index(0, 0, 0)], 0.0);
        // At 0.3 nm (between r − dr and r + dr): partially inside, a count of the 27 cubelets
        // (centred on the grid point, 0.0333 nm apart): only the nearest layer reaches in.
        let v = inside[g.index(8, 5, 5)];
        let ctoi = (1.0f64 / 3.0) * (1.0 / 3.0) * (1.0 / 3.0);
        assert!(v > 0.0 && v < 1.0, "{v}");
        assert_eq!((v / ctoi).round() * ctoi, v);
        assert!(v <= 9.0 * ctoi);
        // Beyond r + dr: empty.
        assert_eq!(inside[g.index(9, 5, 5)], 0.0);
        // Periodic: the point at 0.0 sees the image of an atom at 0.95.
        let sys = single_atom(Vec3::new(0.95, 0.5, 0.5), 0.28);
        let inside = inside_sol(&g, 1, &sys.positions, &sys.atoms, &sys.radii);
        assert_eq!(inside[g.index(0, 5, 5)], 1.0);
        // A zero radius occupies nothing.
        let sys = single_atom(Vec3::new(0.5, 0.5, 0.5), 0.0);
        let inside = inside_sol(&g, 1, &sys.positions, &sys.atoms, &sys.radii);
        assert!(inside.iter().all(|&v| v == 0.0));
    }
}
