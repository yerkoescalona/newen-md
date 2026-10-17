//! The incomplete-Cholesky preconditioner and the matrix–vector product of the finite-difference
//! Poisson operator: gromos++ `FDPoissonBoltzmann_ICCG_NPBC` and `FDPoissonBoltzmann_ICCG_PBC`.
//!
//! The operator is `A = diag(ε_C) − (ε_I, ε_J, ε_K off-diagonals)` on the flat grid index
//! `i + GPX·(j + GPY·k)`. The non-periodic variant addresses neighbours by flat offset with no row
//! checks (the zeroed far faces of the ε maps cut the coupling); the periodic one wraps every axis.
//! Both keep gromos++'s summation order term by term.

/// One ICCG flavour: the preconditioner diagonal, `z = M⁻¹ r` and `z = A p`.
pub trait Iccg {
    /// `initpciccg`: the diagonal of the incomplete-Cholesky factor.
    fn init(&self, l: &mut [f64], c: &[f64], ei: &[f64], ej: &[f64], ek: &[f64]);
    /// `pciccg`: `z = M⁻¹ r` by forward and back substitution.
    fn precondition(&self, z: &mut [f64], l: &[f64], r: &[f64], ei: &[f64], ej: &[f64], ek: &[f64]);
    /// `gqact`: `z = A p`.
    fn apply(&self, z: &mut [f64], p: &[f64], c: &[f64], ei: &[f64], ej: &[f64], ek: &[f64]);
}

/// `FDPoissonBoltzmann_ICCG_NPBC`.
#[derive(Debug, Clone, Copy)]
pub struct IccgNpbc {
    gpx: usize,
    gpxgpy: usize,
    n: usize,
}

impl IccgNpbc {
    pub fn new(gpx: usize, gpy: usize, gpz: usize) -> Self {
        IccgNpbc {
            gpx,
            gpxgpy: gpx * gpy,
            n: gpx * gpy * gpz,
        }
    }
}

impl Iccg for IccgNpbc {
    // ICCG_NPBC.cc:44-135. The `alpha = 0` "modified" terms add an exact +0 and are left out.
    fn init(&self, l: &mut [f64], c: &[f64], ei: &[f64], ej: &[f64], ek: &[f64]) {
        let (gpx, gpxgpy, n) = (self.gpx, self.gpxgpy, self.n);
        l[0] = c[0];
        for i in 1..gpx {
            l[i] = c[i] - (ei[i - 1] * ei[i - 1]) / l[i - 1];
        }
        for i in gpx..gpxgpy {
            l[i] = c[i]
                - (ei[i - 1] * ei[i - 1]) / l[i - 1]
                - (ej[i - gpx] * ej[i - gpx]) / l[i - gpx];
        }
        for i in gpxgpy..n {
            l[i] = c[i]
                - (ei[i - 1] * ei[i - 1]) / l[i - 1]
                - (ej[i - gpx] * ej[i - gpx]) / l[i - gpx]
                - (ek[i - gpxgpy] * ek[i - gpxgpy]) / l[i - gpxgpy];
        }
    }

    // ICCG_NPBC.cc:138-185.
    fn precondition(
        &self,
        z: &mut [f64],
        l: &[f64],
        r: &[f64],
        ei: &[f64],
        ej: &[f64],
        ek: &[f64],
    ) {
        let (gpx, gpxgpy, n) = (self.gpx, self.gpxgpy, self.n);
        z[0] = r[0] / l[0];
        for i in 1..gpx {
            z[i] = (r[i] + ei[i - 1] * z[i - 1]) / l[i];
        }
        for i in gpx..gpxgpy {
            z[i] = (r[i] + (ei[i - 1] * z[i - 1]) + (ej[i - gpx] * z[i - gpx])) / l[i];
        }
        for i in gpxgpy..n {
            z[i] = (r[i]
                + ei[i - 1] * z[i - 1]
                + ej[i - gpx] * z[i - gpx]
                + ek[i - gpxgpy] * z[i - gpxgpy])
                / l[i];
        }

        // The last element is not scaled (the loop stops at N−2).
        for i in 0..n - 1 {
            z[i] *= l[i];
        }

        for i in (n - gpx..n - 1).rev() {
            z[i] = (z[i] + ei[i] * z[i + 1]) / l[i];
        }
        for i in (n - gpxgpy..n - gpx).rev() {
            z[i] = (z[i] + ei[i] * z[i + 1] + ej[i] * z[i + gpx]) / l[i];
        }
        for i in (0..n - gpxgpy).rev() {
            z[i] = (z[i] + ei[i] * z[i + 1] + ej[i] * z[i + gpx] + ek[i] * z[i + gpxgpy]) / l[i];
        }
    }

    // ICCG_NPBC.cc:188-212: seven separate sweeps, so every element accumulates
    // diag, +x, +y, +z, −x, −y, −z in that order.
    fn apply(&self, z: &mut [f64], p: &[f64], c: &[f64], ei: &[f64], ej: &[f64], ek: &[f64]) {
        let (gpx, gpxgpy, n) = (self.gpx, self.gpxgpy, self.n);
        for i in 0..n {
            z[i] = c[i] * p[i];
        }
        for i in 0..n - 1 {
            z[i] -= ei[i] * p[i + 1];
        }
        for i in 0..n - gpx {
            z[i] -= ej[i] * p[i + gpx];
        }
        for i in 0..n - gpxgpy {
            z[i] -= ek[i] * p[i + gpxgpy];
        }
        for i in 1..n {
            z[i] -= ei[i - 1] * p[i - 1];
        }
        for i in gpx..n {
            z[i] -= ej[i - gpx] * p[i - gpx];
        }
        for i in gpxgpy..n {
            z[i] -= ek[i - gpxgpy] * p[i - gpxgpy];
        }
    }
}

/// `FDPoissonBoltzmann_ICCG_PBC`.
#[derive(Debug, Clone, Copy)]
pub struct IccgPbc {
    gpx: usize,
    gpy: usize,
    gpz: usize,
    gpxgpy: usize,
    n: usize,
}

impl IccgPbc {
    pub fn new(gpx: usize, gpy: usize, gpz: usize) -> Self {
        IccgPbc {
            gpx,
            gpy,
            gpz,
            gpxgpy: gpx * gpy,
            n: gpx * gpy * gpz,
        }
    }
}

impl Iccg for IccgPbc {
    // ICCG_PBC.cc:45-128.
    fn init(&self, l: &mut [f64], c: &[f64], ei: &[f64], ej: &[f64], ek: &[f64]) {
        let (gpx, gpy, gpz, gpxgpy, n) = (self.gpx, self.gpy, self.gpz, self.gpxgpy, self.n);
        for k in 0..gpz {
            for j in 0..gpy {
                for i in 0..gpx {
                    let idx = k * gpxgpy + j * gpx + i;
                    // In place, as gromos++ does: on a one-point axis the wrap partner is the
                    // point itself, already partly updated.
                    l[idx] = c[idx];
                    if i != 0 {
                        l[idx] -= (ei[idx - 1] * ei[idx - 1]) / l[idx - 1];
                    }
                    if i == gpx - 1 {
                        l[idx] -= (ei[idx] * ei[idx]) / l[idx + 1 - gpx];
                    }
                    if j != 0 {
                        l[idx] -= (ej[idx - gpx] * ej[idx - gpx]) / l[idx - gpx];
                    }
                    if j == gpy - 1 {
                        l[idx] -= (ej[idx] * ej[idx]) / l[idx + gpx - gpxgpy];
                    }
                    if k != 0 {
                        l[idx] -= (ek[idx - gpxgpy] * ek[idx - gpxgpy]) / l[idx - gpxgpy];
                    }
                    if k == gpz - 1 {
                        l[idx] -= (ek[idx] * ek[idx]) / l[idx + gpxgpy - n];
                    }
                }
            }
        }
    }

    // ICCG_PBC.cc:131-195.
    fn precondition(
        &self,
        z: &mut [f64],
        l: &[f64],
        r: &[f64],
        ei: &[f64],
        ej: &[f64],
        ek: &[f64],
    ) {
        let (gpx, gpy, gpz, gpxgpy, n) = (self.gpx, self.gpy, self.gpz, self.gpxgpy, self.n);
        for k in 0..gpz {
            for j in 0..gpy {
                for i in 0..gpx {
                    let idx = k * gpxgpy + j * gpx + i;
                    let mut t = r[idx];
                    if i != 0 {
                        t += ei[idx - 1] * z[idx - 1];
                    }
                    if i == gpx - 1 {
                        t += ei[idx] * z[idx + 1 - gpx];
                    }
                    if j != 0 {
                        t += ej[idx - gpx] * z[idx - gpx];
                    }
                    if j == gpy - 1 {
                        t += ej[idx] * z[idx + gpx - gpxgpy];
                    }
                    if k != 0 {
                        t += ek[idx - gpxgpy] * z[idx - gpxgpy];
                    }
                    if k == gpz - 1 {
                        t += ek[idx] * z[idx + gpxgpy - n];
                    }
                    z[idx] = t / l[idx];
                }
            }
        }

        for i in 0..n - 1 {
            z[i] *= l[i];
        }

        // The wrap terms (i == 0, j == 0, k == 0) multiply by ε at the point itself, where the
        // transpose of the forward sweep would use ε at the wrapped neighbour (ICCG_PBC.cc:182,
        // 184, 186), so the preconditioner is not symmetric. Kept for parity. The two ε differ
        // only where a solute sphere reaches the periodic grid faces; otherwise both are the
        // solvent value and the result is bit-identical to the symmetric form.
        for k in (0..gpz).rev() {
            for j in (0..gpy).rev() {
                for i in (0..gpx).rev() {
                    let idx = k * gpxgpy + j * gpx + i;
                    let mut t = z[idx];
                    if i != gpx - 1 {
                        t += ei[idx] * z[idx + 1];
                    }
                    if i == 0 {
                        t += ei[idx] * z[idx + gpx - 1];
                    }
                    if j != gpy - 1 {
                        t += ej[idx] * z[idx + gpx];
                    }
                    if j == 0 {
                        t += ej[idx] * z[idx + gpxgpy - gpx];
                    }
                    if k != gpz - 1 {
                        t += ek[idx] * z[idx + gpxgpy];
                    }
                    if k == 0 {
                        t += ek[idx] * z[idx + n - gpxgpy];
                    }
                    z[idx] = t / l[idx];
                }
            }
        }
    }

    // ICCG_PBC.cc:200-246: per point diag, +x, −x, +y, −y, +z, −z with wrapped neighbours.
    fn apply(&self, z: &mut [f64], p: &[f64], c: &[f64], ei: &[f64], ej: &[f64], ek: &[f64]) {
        let (gpx, gpy, gpz, gpxgpy, n) = (self.gpx, self.gpy, self.gpz, self.gpxgpy, self.n);
        for k in 0..gpz {
            for j in 0..gpy {
                for i in 0..gpx {
                    let idx = k * gpxgpy + j * gpx + i;
                    let mut t = c[idx] * p[idx];
                    if i != gpx - 1 {
                        t -= ei[idx] * p[idx + 1];
                    } else {
                        t -= ei[idx] * p[idx + 1 - gpx];
                    }
                    if i != 0 {
                        t -= ei[idx - 1] * p[idx - 1];
                    } else {
                        t -= ei[idx + gpx - 1] * p[idx + gpx - 1];
                    }
                    if j != gpy - 1 {
                        t -= ej[idx] * p[idx + gpx];
                    } else {
                        t -= ej[idx] * p[idx + gpx - gpxgpy];
                    }
                    if j != 0 {
                        t -= ej[idx - gpx] * p[idx - gpx];
                    } else {
                        t -= ej[idx + gpxgpy - gpx] * p[idx + gpxgpy - gpx];
                    }
                    if k != gpz - 1 {
                        t -= ek[idx] * p[idx + gpxgpy];
                    } else {
                        t -= ek[idx] * p[idx + gpxgpy - n];
                    }
                    if k != 0 {
                        t -= ek[idx - gpxgpy] * p[idx - gpxgpy];
                    } else {
                        t -= ek[idx + n - gpxgpy] * p[idx + n - gpxgpy];
                    }
                    z[idx] = t;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fields(n: usize) -> [Vec<f64>; 4] {
        let f = |s: f64| {
            (0..n)
                .map(|i| 1.0 + ((i as f64 * s).sin()).abs())
                .collect::<Vec<_>>()
        };
        [
            f(0.37).iter().map(|x| 10.0 * x).collect(),
            f(0.11),
            f(0.23),
            f(0.71),
        ]
    }

    /// Both operators are symmetric: pᵀ(Aq) = qᵀ(Ap).
    #[test]
    fn operators_are_symmetric() {
        let (gx, gy, gz) = (4, 3, 5);
        let n = gx * gy * gz;
        let [c, ei, ej, ek] = fields(n);
        let p: Vec<f64> = (0..n).map(|i| (i as f64 * 0.3).cos()).collect();
        let q: Vec<f64> = (0..n).map(|i| (i as f64 * 0.7).sin()).collect();
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let (mut ap, mut aq) = (vec![0.0; n], vec![0.0; n]);
        let pbc = IccgPbc::new(gx, gy, gz);
        pbc.apply(&mut ap, &p, &c, &ei, &ej, &ek);
        pbc.apply(&mut aq, &q, &c, &ei, &ej, &ek);
        assert!((dot(&p, &aq) - dot(&q, &ap)).abs() < 1e-10);
        // NPBC addresses neighbours by flat offset, so it is symmetric for any ε fields.
        let npbc = IccgNpbc::new(gx, gy, gz);
        npbc.apply(&mut ap, &p, &c, &ei, &ej, &ek);
        npbc.apply(&mut aq, &q, &c, &ei, &ej, &ek);
        assert!((dot(&p, &aq) - dot(&q, &ap)).abs() < 1e-10);
    }
}
