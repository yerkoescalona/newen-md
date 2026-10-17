//! Finite-difference linearised Poisson–Boltzmann: gromos++ `pb::FDPoissonBoltzmann`
//! (`src/pb/FDPoissonBoltzmann.cc`) with its ICCG solvers.
//!
//! Grid conventions, all from gromos++: the flat index is `x + GPX·(y + GPY·z)` with only `z`
//! wrapped under PBC; a position's grid coordinate is `(r − gridstart)/h` and point `n` sits at
//! grid coordinate `n + 1`; the ε maps are face maps (the I face stored at point `n` lies between
//! points `n` and `n + 1`) holding ε/(4π). Charges are spread trilinearly as q/h, the operator is
//! `Σ ε_face (φ_i − φ_j)`, and φ is in units of e/nm; `dGelec` multiplies by 1/(4πε₀).
//!
//! The non-periodic grid carries a Coulomb (Debye–Hückel with κ = 0) Dirichlet boundary folded
//! into the right-hand side; the periodic one a neutralising background and the average potential
//! removed afterwards.

use super::iccg::{Iccg, IccgNpbc, IccgPbc};
use super::params::PbParameters;
use super::{g, PbSystem};

/// The boundary of one finite-difference solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FdBoundary {
    /// Macroscopic, non-periodic: solvent permittivity on the grid faces, Coulomb potential of the
    /// charges one step outside them.
    Npbc,
    /// Periodic with period `GP·h` on every axis (not the simulation box).
    Pbc,
}

/// The result of one finite-difference solve.
#[derive(Debug, Clone)]
pub struct FdSolution {
    /// Electrostatic potential at every read-out atom, kJ mol⁻¹ e⁻¹ (the total grid potential,
    /// grid self-potential included).
    pub potentials: Vec<f64>,
    /// ½ Σ qᵢ φᵢ over the read-out atoms, kJ/mol.
    pub dg: f64,
    /// ICCG iterations taken (gromos++ counts the initial residual as iteration 1).
    pub iterations: usize,
    /// Whether the L1 residual fell below `convergence_fd` × the initial one.
    pub converged: bool,
    /// Centre of the solute's bounding box, nm.
    pub gridcenter: [f64; 3],
    /// Position of grid coordinate 0 (point 0 sits one spacing further), nm.
    pub gridstart: [f64; 3],
}

/// One finite-difference solve, as `dGslv_pbsolv`'s `fd_ls_{npbc,pbc}_{vac,slv}` run it: a fresh
/// solver object, `setupGrid`, `solveforpotential_{npbc,pbc}` and `dGelec`, logging everything
/// gromos++ logs.
pub fn solve_fd(
    sys: &PbSystem,
    grid: [usize; 3],
    h: f64,
    boundary: FdBoundary,
    eps_solvent: f64,
    maxiter: usize,
    log: &mut String,
) -> Result<FdSolution, String> {
    if grid.contains(&0) {
        return Err(format!(
            "FD grid {grid:?}: every axis needs at least one point"
        ));
    }
    if !(h > 0.0 && h.is_finite()) {
        return Err(format!("FD grid spacing {h} must be positive"));
    }
    if maxiter == 0 {
        return Err("the iteration limit must be positive".into());
    }
    let mut fd =
        FdPoissonBoltzmann::new(sys, grid, h, boundary == FdBoundary::Pbc, eps_solvent, log);
    fd.setup_grid(log)?;
    let acceptance = fd.ppp.convergence_fd;
    let (converged, iterations) = match boundary {
        FdBoundary::Npbc => fd.solve(
            maxiter,
            acceptance,
            &IccgNpbc::new(grid[0], grid[1], grid[2]),
            log,
        ),
        FdBoundary::Pbc => fd.solve(
            maxiter,
            acceptance,
            &IccgPbc::new(grid[0], grid[1], grid[2]),
            log,
        ),
    };
    let mut potentials = Vec::with_capacity(sys.readout.len());
    let dg = fd.dg_elec(log, &mut potentials)?;
    Ok(FdSolution {
        potentials,
        dg,
        iterations,
        converged,
        gridcenter: fd.gridcenter,
        gridstart: fd.gridstart,
    })
}

/// C++ `(int)` of a grid coordinate: truncation toward zero.
#[inline]
fn trunc(x: f64) -> i64 {
    x as i32 as i64
}

struct FdPoissonBoltzmann<'a> {
    sys: &'a PbSystem,
    pbc: bool,
    epssolvent: f64,
    epssolute: f64,
    ppp: PbParameters,
    gpx: i64,
    gpy: i64,
    gpz: i64,
    n: usize,
    h: f64,
    /// Radius of every `atoms` entry, in order.
    radii: Vec<f64>,
    gridcenter: [f64; 3],
    gridstart: [f64; 3],
    phi: Vec<f64>,
    rho: Vec<f64>,
    eps_i: Vec<f64>,
    eps_j: Vec<f64>,
    eps_k: Vec<f64>,
    eps_c: Vec<f64>,
}

impl<'a> FdPoissonBoltzmann<'a> {
    // FDPoissonBoltzmann.cc:43-78
    fn new(
        sys: &'a PbSystem,
        grid: [usize; 3],
        h: f64,
        pbc: bool,
        epssolvent: f64,
        log: &mut String,
    ) -> Self {
        let ppp = PbParameters::new_logged(epssolvent, sys.fpepsi, log);
        let n = grid[0] * grid[1] * grid[2];
        log.push_str(&format!("# gridspacing {}\n", g(h)));
        log.push_str(&format!("# GPX {}\n", grid[0]));
        log.push_str(&format!("# GPY {}\n", grid[1]));
        log.push_str(&format!("# GPZ {}\n", grid[2]));
        log.push_str(&format!("# epssolvent {}\n", g(epssolvent)));
        FdPoissonBoltzmann {
            sys,
            pbc,
            epssolvent,
            epssolute: ppp.epssolute,
            ppp,
            gpx: grid[0] as i64,
            gpy: grid[1] as i64,
            gpz: grid[2] as i64,
            n,
            h,
            radii: sys.atoms.iter().map(|&a| sys.radii[a]).collect(),
            gridcenter: [0.0; 3],
            gridstart: [0.0; 3],
            phi: vec![0.0; n],
            rho: vec![0.0; n],
            eps_i: vec![0.0; n],
            eps_j: vec![0.0; n],
            eps_k: vec![0.0; n],
            eps_c: vec![0.0; n],
        }
    }

    /// `index(x, y, z)`: only `z` wraps, and only under PBC (`.cc:709-716`).
    #[inline]
    fn index(&self, x: i64, y: i64, z: i64) -> i64 {
        let z = if self.pbc { z % self.gpz } else { z };
        x + self.gpx * (y + self.gpy * z)
    }

    /// `index(x, y, z)` for a data-dependent point; gromos++ reads or writes out of bounds there
    /// (undefined behaviour), we refuse.
    fn slot(&self, x: i64, y: i64, z: i64) -> Result<usize, String> {
        let f = self.index(x, y, z);
        if f < 0 || f as usize >= self.n {
            return Err(format!(
                "grid point ({x}, {y}, {z}) is outside the {}x{}x{} grid: an atom lies too close to \
                 (or beyond) the grid edge",
                self.gpx, self.gpy, self.gpz
            ));
        }
        Ok(f as usize)
    }

    fn grid_coordinates(&self, atom: usize) -> [f64; 3] {
        let p = self.sys.positions[atom];
        [
            (p[0] - self.gridstart[0]) / self.h,
            (p[1] - self.gridstart[1]) / self.h,
            (p[2] - self.gridstart[2]) / self.h,
        ]
    }

    // setupGrid(newphi = true, gridstart = 0): .cc:81-231
    fn setup_grid(&mut self, log: &mut String) -> Result<(), String> {
        self.phi.fill(0.0);
        self.rho.fill(0.0);

        // Both branches centre on the solute's bounding box, with the same log text.
        let (lo, hi) = self.sys.extent();
        log.push_str("# Calculating gridcenter for npbc based on cog\n");
        log.push_str(&format!(
            "# mincoordinates of all atoms: {} {} {}\n",
            g(lo[0]),
            g(lo[1]),
            g(lo[2])
        ));
        log.push_str(&format!(
            "# maxcoordinates of all atoms: {} {} {}\n",
            g(hi[0]),
            g(hi[1]),
            g(hi[2])
        ));
        for d in 0..3 {
            self.gridcenter[d] = (lo[d] + hi[d]) * 0.5;
        }
        log.push_str(&format!("# GPX {}\n", self.gpx));
        log.push_str(&format!("# GPY {}\n", self.gpy));
        log.push_str(&format!("# GPZ {}\n", self.gpz));

        log.push_str("# Calculating new gridstart\n");
        let gp = [self.gpx, self.gpy, self.gpz];
        for d in 0..3 {
            self.gridstart[d] = self.gridcenter[d] - 0.5 * (self.h * gp[d] as f64) - self.h / 2.0;
            log.push_str(&format!("# gridstart[{d}] {}\n", g(self.gridstart[d])));
        }
        let (c, s) = (self.gridcenter, self.gridstart);
        log.push_str(&format!(
            "# gridcenter x,y,z = {} {} {}\n",
            g(c[0]),
            g(c[1]),
            g(c[2])
        ));
        log.push_str(&format!(
            "# gridstart x,y,z = {} {} {}\n",
            g(s[0]),
            g(s[1]),
            g(s[2])
        ));
        log.push_str(&format!("# gridspacing {}\n", g(self.h)));

        self.charge_grid_trilinear(log)?;
        if !self.pbc {
            self.debye_hueckel(log);
        }
        let mut eps = std::mem::take(&mut self.eps_i);
        self.radius_boundary(&mut eps, 0, log)?;
        self.eps_i = eps;
        let mut eps = std::mem::take(&mut self.eps_j);
        self.radius_boundary(&mut eps, 1, log)?;
        self.eps_j = eps;
        let mut eps = std::mem::take(&mut self.eps_k);
        self.radius_boundary(&mut eps, 2, log)?;
        self.eps_k = eps;

        let (gpx, gpy, gpz) = (self.gpx as usize, self.gpy as usize, self.gpz as usize);
        let gpxgpy = gpx * gpy;
        if !self.pbc {
            self.set_boundary_solvent();
            // ε_C after the boundary terms: Hünenberger & McCammon JCP 1999 eq A12.
            let mut count = 0;
            for k in 0..gpz {
                for j in 0..gpy {
                    for i in 0..gpx {
                        if i != 0 {
                            self.eps_c[count] += self.eps_i[count - 1];
                        }
                        if i != gpx - 1 {
                            self.eps_c[count] += self.eps_i[count];
                        }
                        if j != 0 {
                            self.eps_c[count] += self.eps_j[count - gpx];
                        }
                        if j != gpy - 1 {
                            self.eps_c[count] += self.eps_j[count];
                        }
                        if k != 0 {
                            self.eps_c[count] += self.eps_k[count - gpxgpy];
                        }
                        if k != gpz - 1 {
                            self.eps_c[count] += self.eps_k[count];
                        }
                        count += 1;
                    }
                }
            }
        } else {
            let at = |i: usize, j: usize, k: usize| i + gpx * (j + gpy * k);
            for k in 0..gpz {
                for j in 0..gpy {
                    for i in 0..gpx {
                        let idx = at(i, j, k);
                        let left = if i == 0 {
                            at(gpx - 1, j, k)
                        } else {
                            at(i - 1, j, k)
                        };
                        self.eps_c[idx] += self.eps_i[left];
                        self.eps_c[idx] += self.eps_i[idx];
                        let down = if j == 0 {
                            at(i, gpy - 1, k)
                        } else {
                            at(i, j - 1, k)
                        };
                        self.eps_c[idx] += self.eps_j[down];
                        self.eps_c[idx] += self.eps_j[idx];
                        let back = if k == 0 {
                            at(i, j, gpz - 1)
                        } else {
                            at(i, j, k - 1)
                        };
                        self.eps_c[idx] += self.eps_k[back];
                        self.eps_c[idx] += self.eps_k[idx];
                    }
                }
            }
        }
        Ok(())
    }

    fn add_rho(&mut self, x: i64, y: i64, z: i64, w: f64) -> Result<(), String> {
        let s = self.slot(x, y, z)?;
        self.rho[s] += w;
        Ok(())
    }

    // chargeGridtrilinear: .cc:1086-1345. Out-of-grid atoms wrap along one axis only, the first
    // of I, J, K that is out; the J < −1 and K < −1 branches put the upper weights at J−1 / K−1
    // (.cc:1238-1247, 1299-1308), as gromos++ does.
    fn charge_grid_trilinear(&mut self, log: &mut String) -> Result<(), String> {
        let sys = self.sys;
        let h = self.h;
        let (gpx, gpy, gpz) = (self.gpx, self.gpy, self.gpz);
        for (n, &a) in sys.atoms.iter().enumerate() {
            let ijk = self.grid_coordinates(a);
            let f = [
                ijk[0] - trunc(ijk[0]) as f64,
                ijk[1] - trunc(ijk[1]) as f64,
                ijk[2] - trunc(ijk[2]) as f64,
            ];
            let q = sys.charges[a];
            let mut i = trunc(ijk[0]) - 1;
            let mut j = trunc(ijk[1]) - 1;
            let mut k = trunc(ijk[2]) - 1;
            let (x1, y1, z1) = (f[0], f[1], f[2]);
            let (x0, y0, z0) = (1.0 - f[0], 1.0 - f[1], 1.0 - f[2]);
            let w = |a: f64, b: f64, c: f64| q * a * b * c / h;
            let atom = n + 1;

            if i >= gpx - 1 {
                if i >= gpx {
                    let ext = i - (gpx - 1);
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping the grid on the upper X axis by {ext} grids.\n"
                    ));
                    i = ext - 1;
                    self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
                    self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
                    self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                    self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
                } else {
                    log.push_str(&format!(
                        "# Atom {atom} is on the edge of the upper grid X and I+1 would extend by one grid.\n"
                    ));
                    self.add_rho(0, j, k + 1, w(x1, y0, z1))?;
                    self.add_rho(0, j + 1, k, w(x1, y1, z0))?;
                    self.add_rho(0, j, k, w(x1, y0, z0))?;
                    self.add_rho(0, j + 1, k + 1, w(x1, y1, z1))?;
                }
                self.add_rho(i, j, k, w(x0, y0, z0))?;
                self.add_rho(i, j + 1, k, w(x0, y1, z0))?;
                self.add_rho(i, j, k + 1, w(x0, y0, z1))?;
                self.add_rho(i, j + 1, k + 1, w(x0, y1, z1))?;
            } else if i < 0 {
                if i == -1 {
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping on the lower grid X by 1 and I+1 would be exactly 0.\n"
                    ));
                    self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                    self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
                    self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
                    self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
                    self.add_rho(gpx + i, j, k, w(x0, y0, z0))?;
                    self.add_rho(gpx + i, j + 1, k, w(x0, y1, z0))?;
                    self.add_rho(gpx + i, j, k + 1, w(x0, y0, z1))?;
                    self.add_rho(gpx + i, j + 1, k + 1, w(x0, y1, z1))?;
                } else {
                    let ext = i;
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping the grid on the lower X axis by {ext} grids.\n"
                    ));
                    i = gpx + ext;
                    self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                    self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
                    self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
                    self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
                    self.add_rho(i, j, k, w(x0, y0, z0))?;
                    self.add_rho(i, j + 1, k, w(x0, y1, z0))?;
                    self.add_rho(i, j, k + 1, w(x0, y0, z1))?;
                    self.add_rho(i, j + 1, k + 1, w(x0, y1, z1))?;
                }
            } else if j >= gpy - 1 {
                if j >= gpy {
                    let ext = j - (gpy - 1);
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping the grid on the upper Y axis by {ext} grids.\n"
                    ));
                    j = ext - 1;
                    self.add_rho(i, j + 1, k, w(x0, y1, z0))?;
                    self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
                    self.add_rho(i, j + 1, k + 1, w(x0, y1, z1))?;
                    self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
                } else {
                    log.push_str(&format!(
                        "# Atom {atom} is on the edge of the upper grid Y and J+1 would extend by one grid.\n"
                    ));
                    self.add_rho(i, 0, k, w(x0, y1, z0))?;
                    self.add_rho(i + 1, 0, k, w(x1, y1, z0))?;
                    self.add_rho(i, 0, k + 1, w(x0, y1, z1))?;
                    self.add_rho(i + 1, 0, k + 1, w(x1, y1, z1))?;
                }
                self.add_rho(i, j, k, w(x0, y0, z0))?;
                self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                self.add_rho(i, j, k + 1, w(x0, y0, z1))?;
                self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
            } else if j < 0 {
                if j == -1 {
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping on the lower grid Y by 1 and J+1 would be exactly 0.\n"
                    ));
                    self.add_rho(i, j + 1, k, w(x0, y1, z0))?;
                    self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
                    self.add_rho(i, j + 1, k + 1, w(x0, y1, z1))?;
                    self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
                    self.add_rho(i, gpy + j, k, w(x0, y0, z0))?;
                    self.add_rho(i + 1, gpy + j, k, w(x1, y0, z0))?;
                    self.add_rho(i, gpy + j, k + 1, w(x0, y0, z1))?;
                    self.add_rho(i + 1, gpy + j, k + 1, w(x1, y0, z1))?;
                } else {
                    let ext = j;
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping the grid on the lower Y axis by {ext} grids.\n"
                    ));
                    j = gpy + ext;
                    self.add_rho(i, j, k, w(x0, y0, z0))?;
                    self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                    self.add_rho(i, j - 1, k, w(x0, y1, z0))?;
                    self.add_rho(i + 1, j - 1, k, w(x1, y1, z0))?;
                    self.add_rho(i, j, k + 1, w(x0, y0, z1))?;
                    self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
                    self.add_rho(i, j - 1, k + 1, w(x0, y1, z1))?;
                    self.add_rho(i + 1, j - 1, k + 1, w(x1, y1, z1))?;
                }
            } else if k >= gpz - 1 {
                if k >= gpz {
                    let ext = k - (gpz - 1);
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping the grid on the upper Z axis by {ext} grids.\n"
                    ));
                    k = ext - 1;
                    self.add_rho(i, j, k + 1, w(x0, y0, z1))?;
                    self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
                    self.add_rho(i, j + 1, k + 1, w(x0, y1, z1))?;
                    self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
                } else {
                    log.push_str(&format!(
                        "# Atom {atom} is on the edge of the upper grid Z and K+1 would extend by one grid.\n"
                    ));
                    self.add_rho(i, j, 0, w(x0, y0, z1))?;
                    self.add_rho(i + 1, j, 0, w(x1, y0, z1))?;
                    self.add_rho(i, j + 1, 0, w(x0, y1, z1))?;
                    self.add_rho(i + 1, j + 1, 0, w(x1, y1, z1))?;
                }
                self.add_rho(i, j, k, w(x0, y0, z0))?;
                self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                self.add_rho(i, j + 1, k, w(x0, y1, z0))?;
                self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
            } else if k < 0 {
                if k == -1 {
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping on the lower grid Z by 1 and K+1 would be exactly 0.\n"
                    ));
                    self.add_rho(i, j, k + 1, w(x0, y0, z1))?;
                    self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
                    self.add_rho(i, j + 1, k + 1, w(x0, y1, z1))?;
                    self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
                    self.add_rho(i, j, gpz + k, w(x0, y0, z0))?;
                    self.add_rho(i + 1, j, gpz + k, w(x1, y0, z0))?;
                    self.add_rho(i, j + 1, gpz + k, w(x0, y1, z0))?;
                    self.add_rho(i + 1, j + 1, gpz + k, w(x1, y1, z0))?;
                } else {
                    let ext = k;
                    log.push_str(&format!(
                        "# Atom {atom} is overlapping the grid on the lower Z axis by {ext} grids.\n"
                    ));
                    k = gpz + ext;
                    self.add_rho(i, j, k, w(x0, y0, z0))?;
                    self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                    self.add_rho(i, j + 1, k, w(x0, y1, z0))?;
                    self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
                    self.add_rho(i, j, k - 1, w(x0, y0, z1))?;
                    self.add_rho(i + 1, j, k - 1, w(x1, y0, z1))?;
                    self.add_rho(i, j + 1, k - 1, w(x0, y1, z1))?;
                    self.add_rho(i + 1, j + 1, k - 1, w(x1, y1, z1))?;
                }
            } else {
                self.add_rho(i, j, k, w(x0, y0, z0))?;
                self.add_rho(i + 1, j, k, w(x1, y0, z0))?;
                self.add_rho(i, j + 1, k, w(x0, y1, z0))?;
                self.add_rho(i + 1, j + 1, k, w(x1, y1, z0))?;
                self.add_rho(i, j, k + 1, w(x0, y0, z1))?;
                self.add_rho(i + 1, j, k + 1, w(x1, y0, z1))?;
                self.add_rho(i, j + 1, k + 1, w(x0, y1, z1))?;
                self.add_rho(i + 1, j + 1, k + 1, w(x1, y1, z1))?;
            }
        }

        // Homogeneous neutralising background under PBC.
        if self.pbc {
            log.push_str("# @@@@@@@@@@@@@@@ pbc-> subtract aver \n");
            let mut sumchrg = 0.0;
            for r in &self.rho {
                sumchrg += r * h;
            }
            sumchrg /= self.n as f64 * h;
            for r in &mut self.rho {
                *r -= sumchrg;
            }
        }

        let mut znorm = 0.0;
        for r in &self.rho {
            znorm += r.abs();
        }
        log.push_str(&format!("# @chargegridtrilinear: znorm : {}\n", g(znorm)));
        Ok(())
    }

    // DebyeHueckel: .cc:1354-1525. Per atom the faces K=0, K=GPZ−1, I=0, I=GPX−1, J=0, J=GPY−1,
    // each at the virtual points one step outside the grid, with the squared distances built
    // incrementally as gromos++ builds them.
    fn debye_hueckel(&mut self, log: &mut String) {
        log.push_str("# Entering FDPoissonBoltzmann::DebyeHueckel ...\n");
        let sys = self.sys;
        let pi4 = 4.0 * self.ppp.pi;
        let kappa = self.ppp.kappa;
        let h = self.h;
        let gp = [self.gpx, self.gpy, self.gpz];
        for (n, &a) in sys.atoms.iter().enumerate() {
            let rad = self.radii[n];
            let factor = sys.charges[a] / (pi4 * (1.0 + kappa * rad));
            let pos = sys.positions[a];
            let near = [
                self.gridstart[0] - pos[0],
                self.gridstart[1] - pos[1],
                self.gridstart[2] - pos[2],
            ];
            let far = |d: usize| near[d] + (gp[d] + 1) as f64 * h;
            // (fixed axis, its point, its offset; outer axis; inner axis)
            for (fixed, point, offset, outer, inner) in [
                (2, 0, near[2], 1, 0),
                (2, gp[2] - 1, far(2), 1, 0),
                (0, 0, near[0], 2, 1),
                (0, gp[0] - 1, far(0), 2, 1),
                (1, 0, near[1], 2, 0),
                (1, gp[1] - 1, far(1), 2, 0),
            ] {
                self.debye_hueckel_face([fixed, outer, inner], point, offset, near, factor, rad);
            }
        }
    }

    fn debye_hueckel_face(
        &mut self,
        [fixed, outer, inner]: [usize; 3],
        point: i64,
        offset: f64,
        near: [f64; 3],
        factor: f64,
        rad: f64,
    ) {
        let h = self.h;
        let h_sq = h * h;
        let two_h = 2.0 * h;
        let kappa = self.ppp.kappa;
        let gp = [self.gpx, self.gpy, self.gpz];
        let mut dvec = near;
        let mut dvec_sq = [near[0] * near[0], near[1] * near[1], near[2] * near[2]];
        dvec[fixed] = offset;
        dvec_sq[fixed] = offset * offset;
        let mut coord = [0i64; 3];
        coord[fixed] = point;
        for o in 0..gp[outer] {
            dvec_sq[outer] += two_h * dvec[outer] + h_sq;
            dvec[outer] += h;
            let mut dc = dvec[inner];
            let mut dc_sq = dc * dc;
            let rest = dvec_sq[outer] + dvec_sq[fixed];
            coord[outer] = o;
            for c in 0..gp[inner] {
                dc_sq += two_h * dc + h_sq;
                dc += h;
                let dist = (dc_sq + rest).sqrt();
                coord[inner] = c;
                // NPBC only: every face point is on the grid.
                let idx = (coord[0] + self.gpx * (coord[1] + self.gpy * coord[2])) as usize;
                self.rho[idx] += factor * (kappa * (rad - dist)).exp() / dist;
            }
        }
    }

    // setboundarySolvent: .cc:719-747. ε_solvent/(4π) once per boundary face of a point.
    fn set_boundary_solvent(&mut self) {
        let epss4pi = self.epssolvent / (4.0 * self.ppp.pi);
        let (gpx, gpy, gpz) = (self.gpx as usize, self.gpy as usize, self.gpz as usize);
        let at = |i: usize, j: usize, k: usize| i + gpx * (j + gpy * k);
        for j in 0..gpy {
            for i in 0..gpx {
                self.eps_c[at(i, j, 0)] += epss4pi;
                self.eps_c[at(i, j, gpz - 1)] += epss4pi;
            }
        }
        for k in 0..gpz {
            for i in 0..gpx {
                self.eps_c[at(i, 0, k)] += epss4pi;
                self.eps_c[at(i, gpy - 1, k)] += epss4pi;
            }
        }
        for k in 0..gpz {
            for j in 0..gpy {
                self.eps_c[at(0, j, k)] += epss4pi;
                self.eps_c[at(gpx - 1, j, k)] += epss4pi;
            }
        }
    }

    // radiusboundaryEPSI/J/K: .cc:750-1084. A sharp two-valued face map: solvent everywhere, the
    // far face zeroed under NPBC (before the atoms, so an atom can overwrite it), then every face
    // point inside an atom sphere set to solute. `face` is the axis the map's faces are normal to;
    // the sqrt chain runs outer → middle → inner with the ±½ face offset on the inner (= face)
    // axis, and out-of-grid points wrap along the first out-of-range axis in that same order.
    fn radius_boundary(
        &self,
        eps: &mut [f64],
        face: usize,
        log: &mut String,
    ) -> Result<(), String> {
        let (outer, middle, inner) = match face {
            0 => (2, 1, 0),
            1 => (2, 0, 1),
            _ => (1, 0, 2),
        };
        let epsi4pi = self.epssolute / (4.0 * self.ppp.pi);
        let epss4pi = self.epssolvent / (4.0 * self.ppp.pi);
        if face == 0 {
            log.push_str(&format!("# epsi4PI: {}\n", g(epsi4pi)));
            log.push_str(&format!("# epss4PI: {}\n", g(epss4pi)));
        }

        eps.fill(epss4pi);
        let gp = [self.gpx, self.gpy, self.gpz];
        let far_value = if self.pbc { epss4pi } else { 0.0 };
        let (a, b) = match face {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };
        for u in 0..gp[a] {
            for v in 0..gp[b] {
                let mut c = [0i64; 3];
                c[face] = gp[face] - 1;
                c[a] = u;
                c[b] = v;
                let idx = self.index(c[0], c[1], c[2]) as usize;
                eps[idx] = far_value;
            }
        }

        let fuzoff = 0.5;
        for (n, &atom) in self.sys.atoms.iter().enumerate() {
            let range1 = self.radii[n] / self.h;
            let ijk = self.grid_coordinates(atom);
            let low1 = trunc(ijk[outer] - range1 + 1.0);
            let high1 = trunc(ijk[outer] + range1);
            for o in low1..=high1 {
                let d = ijk[outer] - o as f64;
                let arg = range1 * range1 - d * d;
                let mut range2 = arg.sqrt();
                if arg < 0.0 {
                    range2 = -range2;
                }
                let low2 = trunc(ijk[middle] - range2 + 1.0);
                let high2 = trunc(ijk[middle] + range2);
                for m in low2..=high2 {
                    let d = ijk[middle] - m as f64;
                    let arg = range2 * range2 - d * d;
                    let range3 = arg.sqrt();
                    // gromos++ flips range2 here (a typo for range3); harmless, range2 is only
                    // used squared from now on.
                    if arg < 0.0 {
                        range2 = -range2;
                    }
                    if range3 > 0.0 {
                        let low3 = trunc(ijk[inner] - range3 + fuzoff);
                        let high3 = trunc(ijk[inner] + range3 - fuzoff);
                        for c in low3..=high3 {
                            let mut v = [0i64; 3];
                            v[outer] = o;
                            v[middle] = m;
                            v[inner] = c;
                            let mut w = v;
                            for ax in [outer, middle, inner] {
                                if v[ax] > gp[ax] {
                                    w[ax] = v[ax] - gp[ax];
                                    break;
                                } else if v[ax] < 1 {
                                    w[ax] = gp[ax] + v[ax];
                                    break;
                                }
                            }
                            let s = self.slot(w[0] - 1, w[1] - 1, w[2] - 1)?;
                            eps[s] = epsi4pi;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    // solveforpotential_{pbc,npbc} (identical code): .cc:236-381, 387-531. Preconditioned CG
    // with an L1 residual test relative to the initial |ρ|; ρ becomes the residual.
    fn solve<S: Iccg>(
        &mut self,
        maxits: usize,
        acceptance: f64,
        iccg: &S,
        log: &mut String,
    ) -> (bool, usize) {
        let n = self.n;
        let mut z = vec![0.0; n];
        let mut p = vec![0.0; n];
        let mut l = vec![0.0; n];
        iccg.init(&mut l, &self.eps_c, &self.eps_i, &self.eps_j, &self.eps_k);

        let mut znorm = 0.0;
        for i in 0..n {
            p[i] = self.phi[i];
            znorm += self.rho[i].abs();
        }
        log.push_str(&format!("# @ ZNORM {}\n", g(znorm)));

        iccg.apply(
            &mut z,
            &p,
            &self.eps_c,
            &self.eps_i,
            &self.eps_j,
            &self.eps_k,
        );

        let bbeta0 = 0.0;
        let mut iter = 1usize;
        let mut anorm = 0.0;
        for i in 0..n {
            self.rho[i] -= z[i];
            anorm += self.rho[i].abs();
        }
        if anorm <= acceptance * znorm {
            log.push_str(&format!("# CONVERGED AFTER {iter} iterations...\n"));
            log.push_str("# Exit: solveforpotential()\n");
            return (true, iter);
        } else if iter == maxits {
            no_convergence(log, iter, anorm, acceptance * znorm);
            return (false, iter);
        }

        iccg.precondition(&mut z, &l, &self.rho, &self.eps_i, &self.eps_j, &self.eps_k);
        let mut rdotz2 = 0.0;
        for i in 0..n {
            rdotz2 += self.rho[i] * z[i];
        }
        for i in 0..n {
            p[i] = bbeta0 * p[i] + z[i];
        }
        let mut rdotz1 = rdotz2;
        iccg.apply(
            &mut z,
            &p,
            &self.eps_c,
            &self.eps_i,
            &self.eps_j,
            &self.eps_k,
        );

        let mut converged = false;
        while iter != maxits {
            iter += 1;
            let mut pdotz = 0.0;
            for i in 0..n {
                pdotz += p[i] * z[i];
            }
            let aalpha = rdotz1 / pdotz;
            anorm = 0.0;
            for i in 0..n {
                self.phi[i] += aalpha * p[i];
                self.rho[i] -= aalpha * z[i];
                anorm += self.rho[i].abs();
            }
            if iter.is_multiple_of(100) {
                log.push_str(&format!("#iteration: {iter} anorm {}\n", g(anorm)));
            }
            if anorm <= acceptance * znorm {
                log.push_str(&format!(
                    "# iter anorm acceptance*znorm {iter} {} {}\n",
                    g(anorm),
                    g(acceptance * znorm)
                ));
                log.push_str(&format!("# CONVERGED AFTER {iter} iterations...\n"));
                log.push_str("# Exit: solveforpotential()\n");
                converged = true;
                break;
            } else if iter == maxits {
                no_convergence(log, iter, anorm, acceptance * znorm);
                break;
            }

            iccg.precondition(&mut z, &l, &self.rho, &self.eps_i, &self.eps_j, &self.eps_k);
            rdotz2 = 0.0;
            for i in 0..n {
                rdotz2 += self.rho[i] * z[i];
            }
            let bbeta = rdotz2 / rdotz1;
            for i in 0..n {
                p[i] = bbeta * p[i] + z[i];
            }
            rdotz1 = rdotz2;
            iccg.apply(
                &mut z,
                &p,
                &self.eps_c,
                &self.eps_i,
                &self.eps_j,
                &self.eps_k,
            );
        }
        (converged, iter)
    }

    // dGelec: .cc:536-570. Under PBC the grid average of φ·FPEPSI is removed with a ×FPEPSI ÷FPEPSI
    // round trip, kept for bit parity.
    fn dg_elec(&mut self, log: &mut String, potentials: &mut Vec<f64>) -> Result<f64, String> {
        let n = self.n as f64;
        let fpepsi = self.ppp.fpepsi;
        let mut phitot = 0.0;
        for v in &self.phi {
            phitot += v;
        }
        log.push_str(&format!(
            "# Average potential [kJ/mol] {}\n",
            g(phitot / n * fpepsi)
        ));
        if self.pbc {
            for v in &mut self.phi {
                *v *= fpepsi;
            }
            let mut phiaver = 0.0;
            let mut phiaver2 = 0.0;
            for v in &self.phi {
                phiaver += v;
            }
            phiaver /= n;
            for v in &mut self.phi {
                *v -= phiaver;
                phiaver2 += *v;
            }
            for v in &mut self.phi {
                *v /= fpepsi;
            }
            log.push_str(&format!(
                "# AVERAGE POTENTIAL SUBTRACTED [kJ/mol]: {} -> {}\n",
                g(phiaver),
                g(phiaver2 / n)
            ));
        }
        let dg = self.dg_restricted(log, potentials)? * fpepsi;
        log.push_str(&format!("# DG RETURNED [kJ/mol]: {}\n", g(dg)));
        Ok(dg)
    }

    fn phi_at(&self, x: i64, y: i64, z: i64) -> Result<f64, String> {
        Ok(self.phi[self.slot(x, y, z)?])
    }

    // getdG_restricted: .cc:574-705. Trilinear interpolation, corners x-fastest. Under PBC
    // out-of-grid atoms are shifted by one period (the upper-K shift uses J, .cc:628, as in
    // gromos++), then the first of I, J, K on the last plane takes its +1 neighbour from plane 0.
    fn dg_restricted(&self, log: &mut String, potentials: &mut Vec<f64>) -> Result<f64, String> {
        let sys = self.sys;
        let (gpx, gpy, gpz) = (self.gpx, self.gpy, self.gpz);
        let mut potential = 0.0;
        for (n, &a) in sys.readout.iter().enumerate() {
            let ijk = self.grid_coordinates(a);
            let f = [
                ijk[0] - trunc(ijk[0]) as f64,
                ijk[1] - trunc(ijk[1]) as f64,
                ijk[2] - trunc(ijk[2]) as f64,
            ];
            let mut i = trunc(ijk[0]) - 1;
            let mut j = trunc(ijk[1]) - 1;
            let mut k = trunc(ijk[2]) - 1;
            let (x1, y1, z1) = (f[0], f[1], f[2]);
            let (x0, y0, z0) = (1.0 - f[0], 1.0 - f[1], 1.0 - f[2]);

            // Neighbour indices along each axis; the PBC branches replace one of them by 0.
            let (mut i1, mut j1, mut k1) = (i + 1, j + 1, k + 1);
            if self.pbc {
                if i >= gpx {
                    i = i - (gpx - 1) - 1;
                } else if i <= -1 {
                    i += gpx;
                }
                if j >= gpy {
                    j = j - (gpy - 1) - 1;
                } else if j <= -1 {
                    j += gpy;
                }
                if k >= gpz {
                    k = j - (gpz - 1) - 1;
                } else if k <= -1 {
                    k += gpz;
                }
                (i1, j1, k1) = (i + 1, j + 1, k + 1);
                if i == gpx - 1 {
                    i1 = 0;
                } else if j == gpy - 1 {
                    j1 = 0;
                } else if k == gpz - 1 {
                    k1 = 0;
                }
            }
            let pottmp = x0 * y0 * z0 * self.phi_at(i, j, k)?
                + x1 * y0 * z0 * self.phi_at(i1, j, k)?
                + x0 * y1 * z0 * self.phi_at(i, j1, k)?
                + x1 * y1 * z0 * self.phi_at(i1, j1, k)?
                + x0 * y0 * z1 * self.phi_at(i, j, k1)?
                + x1 * y0 * z1 * self.phi_at(i1, j, k1)?
                + x0 * y1 * z1 * self.phi_at(i, j1, k1)?
                + x1 * y1 * z1 * self.phi_at(i1, j1, k1)?;

            let q = sys.charges[a];
            let pottmpfpep = pottmp * self.ppp.fpepsi;
            potentials.push(pottmpfpep);
            log.push_str(&format!(
                "# atom  {} charge {} ele. pot. {}\n",
                n + 1,
                g(q),
                g(pottmpfpep)
            ));
            potential += 0.5 * pottmp * q;
        }
        Ok(potential)
    }
}

fn no_convergence(log: &mut String, iter: usize, anorm: f64, target: f64) {
    log.push_str(&format!("# NO CONVERGENCE AFTER {iter} iterations...\n"));
    log.push_str("# RESULTS MAY NOT BE CORRECT!\n");
    log.push_str(&format!("# Iterations {iter}\n"));
    log.push_str(&format!("# anorm {}\n", g(anorm)));
    log.push_str(&format!("# acceptance * znorm {}\n", g(target)));
    log.push_str("# Exit: solveforpotential()\n");
}
