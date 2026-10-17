//! Rotational and translational least-squares superposition.
//!
//! Implements the **Horn (1987) quaternion method** — the standard approach used in
//! most MD analysis packages.  Unlike the newen GSL/eigendecomposition path, the
//! quaternion formulation is sign-invariant and has no degenerate-eigenvalue edge
//! cases.  The result is identical for non-degenerate cases.
//!
//! # Algorithm (Horn 1987, J. Opt. Soc. Am. A 4:629)
//! 1. Centre both sets on their weighted COG.
//! 2. Accumulate the 3×3 cross-sum `S[i][j] = Σ w_k · pos_k[i] · ref_k[j]`.
//! 3. Build the 4×4 symmetric matrix K from S.
//! 4. The eigenvector of K with the largest eigenvalue is the optimal unit quaternion.
//! 5. Convert quaternion → rotation matrix.
//! 6. Apply: pos_new = R · (pos − cog_pos) + cog_ref.

use newen_core::math::{Mat3, Vec3};

// ─── 4×4 Jacobi symmetric eigensolver ───────────────────────────────────────

type M44 = [[f64; 4]; 4];

fn eye4() -> M44 {
    let mut m = [[0.0; 4]; 4];
    for i in 0..4 {
        m[i][i] = 1.0;
    }
    m
}

/// Jacobi eigendecomposition of a 4×4 symmetric matrix.
/// Returns (eigenvalues, eigenvectors as columns).
fn jacobi4(mut a: M44) -> ([f64; 4], M44) {
    let mut v = eye4();
    for _ in 0..100 {
        let (mut p, mut q) = (0, 1);
        for i in 0..4 {
            for j in (i + 1)..4 {
                if a[i][j].abs() > a[p][q].abs() {
                    p = i;
                    q = j;
                }
            }
        }
        if a[p][q].abs() < 1e-14 {
            break;
        }

        let theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
        let t = 1.0_f64.copysign(theta) / (theta.abs() + (1.0 + theta * theta).sqrt());
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        let app = a[p][p] - t * a[p][q];
        let aqq = a[q][q] + t * a[p][q];
        a[p][p] = app;
        a[q][q] = aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for r in 0..4 {
            if r != p && r != q {
                let apr = c * a[p][r] - s * a[q][r];
                let aqr = s * a[p][r] + c * a[q][r];
                a[p][r] = apr;
                a[r][p] = apr;
                a[q][r] = aqr;
                a[r][q] = aqr;
            }
        }
        for r in 0..4 {
            let vp = c * v[r][p] - s * v[r][q];
            let vq = s * v[r][p] + c * v[r][q];
            v[r][p] = vp;
            v[r][q] = vq;
        }
    }
    ([a[0][0], a[1][1], a[2][2], a[3][3]], v)
}

// ─── Public API ──────────────────────────────────────────────────────────────

/// Weighted centre of geometry.
pub fn weighted_cog(positions: &[Vec3], indices: &[usize], weights: Option<&[f64]>) -> Vec3 {
    let mut sum = Vec3::ZERO;
    let mut w_total = 0.0f64;
    for &i in indices {
        let w = weights.map_or(1.0, |ws| ws[i]);
        sum += positions[i] * w;
        w_total += w;
    }
    if w_total > 0.0 {
        sum / w_total
    } else {
        Vec3::ZERO
    }
}

/// Build the 4×4 K matrix (Horn 1987) from the cross-sum matrix S.
///
/// S[i][j] = Σ_k w_k · pos_k[i] · ref_k[j]   (positions centred on COG)
fn build_k(s: &[[f64; 3]; 3]) -> M44 {
    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);
    [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ]
}

/// Convert unit quaternion (qw, qx, qy, qz) to 3×3 rotation matrix.
fn quat_to_mat3(q: [f64; 4]) -> Mat3 {
    let (qw, qx, qy, qz) = (q[0], q[1], q[2], q[3]);
    Mat3::from_cols(
        Vec3::new(
            qw * qw + qx * qx - qy * qy - qz * qz,
            2.0 * (qx * qy + qw * qz),
            2.0 * (qx * qz - qw * qy),
        ),
        Vec3::new(
            2.0 * (qx * qy - qw * qz),
            qw * qw - qx * qx + qy * qy - qz * qz,
            2.0 * (qy * qz + qw * qx),
        ),
        Vec3::new(
            2.0 * (qx * qz + qw * qy),
            2.0 * (qy * qz - qw * qx),
            qw * qw - qx * qx - qy * qy + qz * qz,
        ),
    )
}

/// Compute the Kabsch rotation matrix R that minimises weighted RMSD between
/// `positions[fit_indices]` and `reference[fit_indices]`.
///
/// Both sets are expected to be un-centred; this function computes their COGs
/// and centres them internally.
pub fn kabsch_rotation(
    positions: &[Vec3],
    reference: &[Vec3],
    fit_indices: &[usize],
    weights: Option<&[f64]>,
) -> (Mat3, Vec3, Vec3) {
    let pos_cog = weighted_cog(positions, fit_indices, weights);
    let ref_cog = weighted_cog(reference, fit_indices, weights);

    // S[i][j] = Σ_k w_k · p_k[i] · r_k[j]
    let mut s = [[0.0f64; 3]; 3];
    for &idx in fit_indices {
        let p = positions[idx] - pos_cog;
        let r = reference[idx] - ref_cog;
        let w = weights.map_or(1.0, |ws| ws[idx]);
        let pv = [p.x, p.y, p.z];
        let rv = [r.x, r.y, r.z];
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += w * pv[i] * rv[j];
            }
        }
    }

    let k = build_k(&s);
    let (evals, evecs) = jacobi4(k);

    // Eigenvector of largest eigenvalue = optimal quaternion
    let max_idx = evals
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    let q = [
        evecs[0][max_idx],
        evecs[1][max_idx],
        evecs[2][max_idx],
        evecs[3][max_idx],
    ];
    let rot = quat_to_mat3(q);

    (rot, pos_cog, ref_cog)
}

/// Superimpose `positions` onto `reference` using atoms in `fit_indices`.
///
/// Modifies `positions` in-place: applies the optimal rotation + translation.
/// Returns the rotation matrix applied.
pub fn superimpose(
    positions: &mut [Vec3],
    reference: &[Vec3],
    fit_indices: &[usize],
    weights: Option<&[f64]>,
) -> Mat3 {
    let (rot, pos_cog, ref_cog) = kabsch_rotation(positions, reference, fit_indices, weights);
    for p in positions.iter_mut() {
        *p = rot * (*p - pos_cog) + ref_cog;
    }
    rot
}

/// RMSD between `positions[indices]` and `reference[indices]` (no fitting).
pub fn rmsd(positions: &[Vec3], reference: &[Vec3], indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }
    let sum: f64 = indices
        .iter()
        .map(|&i| (positions[i] - reference[i]).length_squared())
        .sum();
    (sum / indices.len() as f64).sqrt()
}

/// Superimpose then compute RMSD over `rmsd_indices` (may differ from `fit_indices`).
pub fn rmsd_after_fit(
    positions: &mut [Vec3],
    reference: &[Vec3],
    fit_indices: &[usize],
    rmsd_indices: &[usize],
    weights: Option<&[f64]>,
) -> f64 {
    superimpose(positions, reference, fit_indices, weights);
    rmsd(positions, reference, rmsd_indices)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    fn triangle() -> Vec<Vec3> {
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ]
    }

    fn indices(n: usize) -> Vec<usize> {
        (0..n).collect()
    }

    #[test]
    fn identity_fit_zero_rmsd() {
        let reference = triangle();
        let mut positions = reference.clone();
        superimpose(&mut positions, &reference, &indices(3), None);
        assert!(rmsd(&positions, &reference, &indices(3)) < 1e-10);
    }

    #[test]
    fn pure_translation_zero_rmsd() {
        let reference = triangle();
        let mut positions: Vec<Vec3> = reference
            .iter()
            .map(|&p| p + Vec3::new(5., -3., 2.))
            .collect();
        superimpose(&mut positions, &reference, &indices(3), None);
        assert!(rmsd(&positions, &reference, &indices(3)) < 1e-10);
    }

    #[test]
    fn rotation_90deg_z_zero_rmsd() {
        let reference = triangle();
        // R90z: x→y, y→-x
        let mut positions: Vec<Vec3> = reference
            .iter()
            .map(|&p| Vec3::new(-p.y, p.x, p.z))
            .collect();
        superimpose(&mut positions, &reference, &indices(3), None);
        assert!(
            rmsd(&positions, &reference, &indices(3)) < 1e-8,
            "RMSD = {}",
            rmsd(&positions, &reference, &indices(3))
        );
    }

    #[test]
    fn rotation_45deg_y_zero_rmsd() {
        let reference = vec![
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
            Vec3::new(0., 0., 1.),
            Vec3::new(1., 1., 1.),
        ];
        let a = std::f64::consts::FRAC_PI_4;
        let mut positions: Vec<Vec3> = reference
            .iter()
            .map(|&p| {
                Vec3::new(
                    p.x * a.cos() + p.z * a.sin(),
                    p.y,
                    -p.x * a.sin() + p.z * a.cos(),
                )
            })
            .collect();
        superimpose(&mut positions, &reference, &indices(4), None);
        assert!(
            rmsd(&positions, &reference, &indices(4)) < 1e-8,
            "RMSD = {}",
            rmsd(&positions, &reference, &indices(4))
        );
    }

    #[test]
    fn rotation_and_translation_zero_rmsd() {
        let reference = vec![
            Vec3::new(0., 0., 0.),
            Vec3::new(1., 0., 0.),
            Vec3::new(0.5, 0.866, 0.),
            Vec3::new(0.5, 0.289, 0.816),
        ];
        let angle = std::f64::consts::PI / 6.;
        let (c, s) = (angle.cos(), angle.sin());
        let mut positions: Vec<Vec3> = reference
            .iter()
            .map(|&p| Vec3::new(c * p.x - s * p.y + 2., s * p.x + c * p.y - 1., p.z + 0.5))
            .collect();
        superimpose(&mut positions, &reference, &indices(4), None);
        assert!(
            rmsd(&positions, &reference, &indices(4)) < 1e-7,
            "RMSD = {}",
            rmsd(&positions, &reference, &indices(4))
        );
    }

    #[test]
    fn rmsd_known_shift() {
        let reference = triangle();
        let positions: Vec<Vec3> = reference
            .iter()
            .map(|&p| p + Vec3::new(0.1, 0., 0.))
            .collect();
        assert!((rmsd(&positions, &reference, &indices(3)) - 0.1).abs() < 1e-10);
    }

    #[test]
    fn fit_subset_indices() {
        // Fit on atoms 0-1 only; atom 2 moves along but isn't in fit set
        let reference = vec![
            Vec3::new(0., 0., 0.),
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
        ];
        let mut positions: Vec<Vec3> = reference
            .iter()
            .map(|&p| Vec3::new(-p.y, p.x, p.z))
            .collect();
        let fit_idx = vec![0, 1];
        superimpose(&mut positions, &reference, &fit_idx, None);
        assert!(rmsd(&positions, &reference, &fit_idx) < 1e-8);
    }
}
