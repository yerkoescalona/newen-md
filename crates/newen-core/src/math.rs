//! Mathematical primitives for molecular dynamics

use glam::{DMat3, DVec3};

/// 3D vector - double precision (f64) for accurate MD simulation
pub type Vec3 = DVec3;

/// 3x3 matrix - double precision (f64)
pub type Mat3 = DMat3;

/// Boundary condition types for periodic systems
pub trait BoundaryCondition: Send + Sync {
    /// Calculate nearest image vector from j to i (r = ri - rj)
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3;

    /// Apply periodic boundary conditions to a position
    fn put_into_box(&self, pos: Vec3) -> Vec3;
}

/// Vacuum boundary conditions (no periodicity)
#[derive(Debug, Clone, Copy)]
pub struct Vacuum;

/// Nearest image without periodicity: the plain difference `ri − rj`.
///
/// spec: periodicity.vacuum_nearest_image
#[inline(always)]
pub fn vacuum_nearest_image(ri: Vec3, rj: Vec3) -> Vec3 {
    ri - rj
}

impl BoundaryCondition for Vacuum {
    #[inline(always)]
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        vacuum_nearest_image(ri, rj)
    }

    #[inline(always)]
    fn put_into_box(&self, pos: Vec3) -> Vec3 {
        pos
    }
}

/// Rectangular periodic boundary conditions
#[derive(Debug, Clone, Copy)]
pub struct Rectangular {
    /// Box edge lengths (nm).
    pub box_size: Vec3,
    /// Precomputed half-box lengths used in the minimum-image test.
    pub half_box: Vec3,
    /// Precomputed reciprocal box lengths (1/box_size).
    pub inv_box: Vec3,
}

impl Rectangular {
    /// Create from box edge lengths; precomputes half-box and reciprocal lengths.
    ///
    /// spec: periodicity.rectangular_new
    pub fn new(box_size: Vec3) -> Self {
        Self {
            box_size,
            half_box: box_size * 0.5,
            inv_box: box_size.recip(),
        }
    }
}

/// Minimum image `ri − rj` in a rectangular box of edges `box_size` with reciprocals `inv_box`:
/// each component minus its nearest multiple of the edge, branchless.
///
/// spec: periodicity.rectangular_nearest_image
#[inline(always)]
pub fn rectangular_nearest_image(box_size: Vec3, inv_box: Vec3, ri: Vec3, rj: Vec3) -> Vec3 {
    let r = ri - rj;

    // Ties-to-even, not ties-away-from-zero (`rint`, a single `roundsd` on x86):
    // `f64::round()` has no native encoding and expands to several extra
    // instructions per component, which measured as ~47% of the whole nonbonded
    // kernel on `water_216_box`. The two differ only at exactly half a box length.
    let k = r * inv_box;
    let n = Vec3::new(
        k.x.round_ties_even(),
        k.y.round_ties_even(),
        k.z.round_ties_even(),
    );

    r - box_size * n
}

/// `pos` moved by whole edges into `[0, L)` per component (rectangular box).
///
/// spec: periodicity.rectangular_put_into_box
#[inline(always)]
pub fn rectangular_put_into_box(box_size: Vec3, inv_box: Vec3, pos: Vec3) -> Vec3 {
    pos - box_size * (pos * inv_box).floor()
}

impl BoundaryCondition for Rectangular {
    #[inline(always)]
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        rectangular_nearest_image(self.box_size, self.inv_box, ri, rj)
    }

    #[inline(always)]
    fn put_into_box(&self, pos: Vec3) -> Vec3 {
        rectangular_put_into_box(self.box_size, self.inv_box, pos)
    }
}

/// Triclinic periodic boundary conditions
#[derive(Debug, Clone, Copy)]
pub struct Triclinic {
    /// Box vectors as columns (lower-triangular GROMOS convention).
    pub box_matrix: Mat3,
    /// Inverse box matrix used for fractional-coordinate wrapping.
    pub inv_box_matrix: Mat3,
}

impl Triclinic {
    /// Create from a box matrix; precomputes the inverse.
    ///
    /// spec: periodicity.triclinic_new
    pub fn new(box_matrix: Mat3) -> Self {
        Self {
            box_matrix,
            inv_box_matrix: box_matrix.inverse(),
        }
    }
}

/// Iterative z -> y -> x lattice reduction of the fractional coordinates, as gromosXX's
/// `Boundary_Implementation<triclinic>::nearest_image` (`math/boundary_implementation.cc`) is
/// written to do. gromosXX's `Box(i, j)` is component `i` of edge `j`, so there the skew terms
/// read `b.z`, `a.y`, `a.z` (zero) and the reduction is into the brick `|x| < a.x/2`,
/// `|y| < b.y/2`, `|z| < c.z/2`; both give the same, shortest image for every vector shorter than
/// half the smallest box height (Lean `image_shortest_of_lt_half_height`). Requires the GROMOS
/// "lower triangular" box vector convention (a=(ax,0,0), b=(bx,by,0),
/// c=(cx,cy,cz)), as guaranteed by GENBOX / `truncoct_triclinic_box`.
///
/// `nim.z / c.z`, `nim_y / b.y` and `nim_x / a.x` are the fractional coordinates
/// of `nim` along c, b and a; each pair of loops brings one into `(-1/2, 1/2]`.
/// The result is the image in that parallelepiped, which in a skewed box is not
/// always the shortest one. Over ℝ, `frac - frac.round()` gives the same image
/// except at exactly half a box vector, where it takes the opposite face.
///
/// spec: periodicity.triclinic_nearest_image
#[inline(always)]
pub fn triclinic_nearest_image(box_matrix: Mat3, ri: Vec3, rj: Vec3) -> Vec3 {
    let a = box_matrix.x_axis;
    let b = box_matrix.y_axis;
    let c = box_matrix.z_axis;
    let mut nim = ri - rj;

    while nim.z >= 0.5 * c.z {
        nim -= c;
    }
    while nim.z <= -0.5 * c.z {
        nim += c;
    }

    let mut nim_y = nim.y - (c.y / c.z) * nim.z;
    while nim_y >= 0.5 * b.y {
        nim -= b;
        nim_y = nim.y - (c.y / c.z) * nim.z;
    }
    while nim_y <= -0.5 * b.y {
        nim += b;
        nim_y = nim.y - (c.y / c.z) * nim.z;
    }

    let mut nim_x = nim.x - b.x * (nim.y - c.y * nim.z / c.z) / b.y - c.x * nim.z / c.z;
    while nim_x >= 0.5 * a.x {
        nim -= a;
        nim_x = nim.x - b.x * (nim.y - c.y * nim.z / c.z) / b.y - c.x * nim.z / c.z;
    }
    while nim_x <= -0.5 * a.x {
        nim += a;
        nim_x = nim.x - b.x * (nim.y - c.y * nim.z / c.z) / b.y - c.x * nim.z / c.z;
    }

    nim
}

/// `pos` moved by whole box vectors into the cell of fractional coordinates `[0, 1)`.
///
/// spec: periodicity.triclinic_put_into_box
#[inline(always)]
pub fn triclinic_put_into_box(box_matrix: Mat3, inv_box_matrix: Mat3, pos: Vec3) -> Vec3 {
    let frac = inv_box_matrix * pos;
    let frac_wrapped = frac - frac.floor();
    box_matrix * frac_wrapped
}

impl BoundaryCondition for Triclinic {
    #[inline(always)]
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        triclinic_nearest_image(self.box_matrix, ri, rj)
    }

    #[inline(always)]
    fn put_into_box(&self, pos: Vec3) -> Vec3 {
        triclinic_put_into_box(self.box_matrix, self.inv_box_matrix, pos)
    }
}

/// Rotation matrix between the truncated-octahedron ("cubic") frame and the
/// lower-triangular triclinic frame produced by `truncoct_triclinic_box`,
/// ported from GROMOS `math::truncoct_triclinic_rotmat`
/// (`math/transformation.cc`). `forward` rotates truncoct -> triclinic
/// (applied to positions/velocities on read for NTB=-1); `!forward` is its
/// transpose (triclinic -> truncoct, applied when writing forces back out).
pub fn truncoct_triclinic_rotmat(forward: bool) -> Mat3 {
    let sq3i = 1.0 / 3.0_f64.sqrt();
    let sq2i = 1.0 / 2.0_f64.sqrt();

    // GROMOS builds `rot` from rows (sq3i,-2*sq2i*sq3i,0), (sq3i,
    // sq3i*sq2i,-sq2i), (sq3i,sq2i*sq3i,sq2i). `rot` below is that same
    // matrix `M` (its columns are GROMOS's rows, per glam's column-major
    // storage). GROMOS applies the rotation via `product(rot, v)`, which
    // (per gmath.h's `product`) computes `rot^T * v`, not `rot * v`. So for
    // `forward=true` (GROMOS's `rot` = `M`, untransposed) the applied
    // rotation is `M^T`; for `forward=false` (GROMOS transposes `rot` to
    // `M^T`, then `product` transposes again) the applied rotation is `M`.
    // See P1.4 derivation.
    let rot = Mat3::from_cols(
        Vec3::new(sq3i, sq3i, sq3i),
        Vec3::new(-2.0 * sq2i * sq3i, sq3i * sq2i, sq2i * sq3i),
        Vec3::new(0.0, -sq2i, sq2i),
    );

    if forward {
        rot.transpose()
    } else {
        rot
    }
}

/// Converts a box matrix between the legacy "cube edge length L" convention
/// (NTB=-1 BOX block, `forward=true`, where `box_matrix.x_axis.x = L`) and
/// the lower-triangular truncated-octahedron triclinic box vectors
/// (`forward=false` recovers the cubic box from the triclinic `a` vector),
/// ported from GROMOS `math::truncoct_triclinic_box`
/// (`math/transformation.cc`).
pub fn truncoct_triclinic_box(box_matrix: Mat3, forward: bool) -> Mat3 {
    if forward {
        let l = box_matrix.x_axis.x;
        let d = 0.5 * 3.0_f64.sqrt() * l;
        let a = Vec3::new(d, 0.0, 0.0);
        let b = Vec3::new(d / 3.0, 2.0 * 2.0_f64.sqrt() / 3.0 * d, 0.0);
        let c = Vec3::new(-d / 3.0, 2.0_f64.sqrt() / 3.0 * d, 6.0_f64.sqrt() / 3.0 * d);
        Mat3::from_cols(a, b, c)
    } else {
        let d = 2.0 * box_matrix.x_axis.x / 3.0_f64.sqrt();
        Mat3::from_diagonal(Vec3::splat(d))
    }
}

/// Rotates every position/velocity vector between the truncated-octahedron
/// and triclinic frames in place, ported from GROMOS
/// `math::truncoct_triclinic` (`math/transformation.cc`). `forward` rotates
/// truncoct -> triclinic (applied on read for NTB=-1); `!forward` is the
/// inverse (applied when writing forces back out).
pub fn truncoct_triclinic(positions: &mut [Vec3], forward: bool) {
    let rot = truncoct_triclinic_rotmat(forward);
    for p in positions.iter_mut() {
        *p = rot * *p;
    }
}

/// The box vectors of a GENBOX block: edge lengths `a, b, c` and angles `α` (b,c), `β` (a,c),
/// `γ` (a,b) in degrees, as gromosXX builds them (`in_configuration.cc:1977-2005`): `a` along x,
/// `b` in the xy plane, entries below 10⁻¹² set to zero (`:2038-2043`, `math::epsilon`). Right
/// angles give the diagonal box. The Euler angles of the block's third line are not applied.
pub fn genbox_box(lengths: Vec3, angles_deg: Vec3) -> Mat3 {
    let to_rad = |d: f64| std::f64::consts::PI * d / 180.0;
    let (alpha, beta, gamma) = (
        to_rad(angles_deg.x),
        to_rad(angles_deg.y),
        to_rad(angles_deg.z),
    );
    let cos_delta = (alpha.cos() - beta.cos() * gamma.cos()) / (beta.sin() * gamma.sin());
    let sin_delta = (1.0 - cos_delta * cos_delta).sqrt();
    let a = Vec3::new(lengths.x, 0.0, 0.0);
    let b = Vec3::new(lengths.y * gamma.cos(), lengths.y * gamma.sin(), 0.0);
    let c = Vec3::new(
        lengths.z * beta.cos(),
        lengths.z * cos_delta * beta.sin(),
        lengths.z * sin_delta * beta.sin(),
    );
    let zero_small = |v: Vec3| {
        let z = |x: f64| if x.abs() < 1e-12 { 0.0 } else { x };
        Vec3::new(z(v.x), z(v.y), z(v.z))
    };
    Mat3::from_cols(zero_small(a), zero_small(b), zero_small(c))
}

/// The edge lengths and angles (degrees) a GENBOX block records for `box_matrix`, as gromosXX's
/// `_print_box` computes them (`out_configuration.cc:1830-1845`): `|a|, |b|, |c|` and the
/// angles between the vectors. A zero edge gives zero angles (`:1879-1883`).
pub fn genbox_lengths_angles(box_matrix: Mat3) -> (Vec3, Vec3) {
    let (a, b, c) = (box_matrix.x_axis, box_matrix.y_axis, box_matrix.z_axis);
    let lengths = Vec3::new(a.length(), b.length(), c.length());
    if lengths.x == 0.0 || lengths.y == 0.0 || lengths.z == 0.0 {
        return (lengths, Vec3::ZERO);
    }
    let angle = |u: Vec3, v: Vec3, lu: f64, lv: f64| {
        (u.dot(v) / (lu * lv)).clamp(-1.0, 1.0).acos() * 180.0 / std::f64::consts::PI
    };
    let angles = Vec3::new(
        angle(b, c, lengths.y, lengths.z),
        angle(a, c, lengths.x, lengths.z),
        angle(a, b, lengths.x, lengths.y),
    );
    (lengths, angles)
}

/// The heights of the cell spanned by the columns `a, b, c` of `box_matrix`: the distance between
/// opposite faces, `V/|b×c|`, `V/|c×a|`, `V/|a×b|` (gromosXX `boundary_check_cutoff`,
/// `math/boundary_checks.cc:76-84`, with the same quantities written through the angles).
pub fn box_heights(box_matrix: Mat3) -> Vec3 {
    let (a, b, c) = (box_matrix.x_axis, box_matrix.y_axis, box_matrix.z_axis);
    let v = a.dot(b.cross(c)).abs();
    Vec3::new(
        v / b.cross(c).length(),
        v / c.cross(a).length(),
        v / a.cross(b).length(),
    )
}

/// Generic periodicity wrapper
#[derive(Debug, Clone)]
pub enum Periodicity {
    /// No periodic boundary — suitable for vacuum / gas-phase systems.
    Vacuum(Vacuum),
    /// Orthorhombic periodic box.
    Rectangular(Rectangular),
    /// Triclinic periodic box (lower-triangular GROMOS convention).
    Triclinic(Triclinic),
}

impl Periodicity {
    /// Apply the minimum-image convention and return the image vector `ri − rj`.
    #[inline]
    pub fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        match self {
            Periodicity::Vacuum(bc) => bc.nearest_image(ri, rj),
            Periodicity::Rectangular(bc) => bc.nearest_image(ri, rj),
            Periodicity::Triclinic(bc) => bc.nearest_image(ri, rj),
        }
    }

    /// Wrap a position back into the primary box.
    #[inline]
    pub fn put_into_box(&self, pos: Vec3) -> Vec3 {
        match self {
            Periodicity::Vacuum(bc) => bc.put_into_box(pos),
            Periodicity::Rectangular(bc) => bc.put_into_box(pos),
            Periodicity::Triclinic(bc) => bc.put_into_box(pos),
        }
    }
}

impl BoundaryCondition for Periodicity {
    #[inline]
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        match self {
            Periodicity::Vacuum(bc) => bc.nearest_image(ri, rj),
            Periodicity::Rectangular(bc) => bc.nearest_image(ri, rj),
            Periodicity::Triclinic(bc) => bc.nearest_image(ri, rj),
        }
    }

    #[inline]
    fn put_into_box(&self, pos: Vec3) -> Vec3 {
        match self {
            Periodicity::Vacuum(bc) => bc.put_into_box(pos),
            Periodicity::Rectangular(bc) => bc.put_into_box(pos),
            Periodicity::Triclinic(bc) => bc.put_into_box(pos),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vacuum_boundary() {
        let bc = Vacuum;
        let ri = Vec3::new(0.0, 0.0, 0.0);
        let rj = Vec3::new(1.0, 2.0, 3.0);

        let r = bc.nearest_image(ri, rj);
        // r = ri - rj = (0,0,0) - (1,2,3) = (-1,-2,-3)
        assert_relative_eq!(r.x, -1.0);
        assert_relative_eq!(r.y, -2.0);
        assert_relative_eq!(r.z, -3.0);
    }

    #[test]
    fn test_rectangular_boundary() {
        let bc = Rectangular::new(Vec3::splat(10.0));
        let ri = Vec3::new(9.5, 0.0, 0.0);
        let rj = Vec3::new(0.5, 0.0, 0.0);

        let r = bc.nearest_image(ri, rj);
        // ri - rj = 9.5 - 0.5 = 9.0, nearest image wraps to -1.0
        assert_relative_eq!(r.x, -1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_triclinic_truncoct_tie_boundary_diverges_from_frac_round() {
        // The gromosXX audit's finding #1 (`.claude/porting.md`): the textbook
        // `frac - frac.round()` reduction and GROMOS's while-loop z->y->x
        // reduction are NOT equivalent for strongly triclinic cells. At an
        // exact half-lattice-vector displacement they land on opposite sides
        // of the cell, a full lattice vector apart.
        let cubic = Mat3::from_diagonal(Vec3::splat(3.767055681));
        let box_matrix = truncoct_triclinic_box(cubic, true);
        let bc = Triclinic::new(box_matrix);
        let c = box_matrix.z_axis;

        let ri = 0.5 * c;
        let rj = Vec3::ZERO;

        // Old textbook reduction: frac = (0, 0, 0.5) -> round() -> (0, 0, 1)
        // -> result = -0.5 * c.
        let frac = bc.inv_box_matrix * (ri - rj);
        let old_result = box_matrix * (frac - frac.round());
        assert_relative_eq!(old_result.x, -0.5 * c.x, epsilon = 1e-9);
        assert_relative_eq!(old_result.y, -0.5 * c.y, epsilon = 1e-9);
        assert_relative_eq!(old_result.z, -0.5 * c.z, epsilon = 1e-9);

        // New while-loop reduction bounces back to +0.5 * c: z is "in range"
        // on both the `>=` and `<=` sides of the boundary.
        let new_result = bc.nearest_image(ri, rj);
        assert_relative_eq!(new_result.x, 0.5 * c.x, epsilon = 1e-9);
        assert_relative_eq!(new_result.y, 0.5 * c.y, epsilon = 1e-9);
        assert_relative_eq!(new_result.z, 0.5 * c.z, epsilon = 1e-9);

        // The two reductions differ by exactly one lattice vector c.
        assert_relative_eq!((new_result - old_result).x, c.x, epsilon = 1e-9);
        assert_relative_eq!((new_result - old_result).y, c.y, epsilon = 1e-9);
        assert_relative_eq!((new_result - old_result).z, c.z, epsilon = 1e-9);
    }

    #[test]
    fn test_triclinic_truncoct_nearest_image_reduces_large_displacement() {
        let cubic = Mat3::from_diagonal(Vec3::splat(3.767055681));
        let box_matrix = truncoct_triclinic_box(cubic, true);
        let bc = Triclinic::new(box_matrix);
        let a = box_matrix.x_axis;
        let b = box_matrix.y_axis;
        let c = box_matrix.z_axis;

        let ri = Vec3::new(10.0, 8.0, 7.0);
        let rj = Vec3::ZERO;

        let nim = bc.nearest_image(ri, rj);

        // The result is `ri - rj` shifted by an integer combination of the
        // box vectors (a valid periodic image).
        let shift = (ri - rj) - nim;
        let frac_shift = bc.inv_box_matrix * shift;
        assert_relative_eq!(frac_shift.x, frac_shift.x.round(), epsilon = 1e-9);
        assert_relative_eq!(frac_shift.y, frac_shift.y.round(), epsilon = 1e-9);
        assert_relative_eq!(frac_shift.z, frac_shift.z.round(), epsilon = 1e-9);

        // z is reduced into [-0.5*c.z, 0.5*c.z).
        assert!(nim.z < 0.5 * c.z + 1e-9);
        assert!(nim.z >= -0.5 * c.z - 1e-9);

        // y is reduced into [-0.5*b.y, 0.5*b.y) after the c-correction.
        let nim_y = nim.y - (c.y / c.z) * nim.z;
        assert!(nim_y < 0.5 * b.y + 1e-9);
        assert!(nim_y >= -0.5 * b.y - 1e-9);

        // x is reduced into [-0.5*a.x, 0.5*a.x) after the b/c-correction.
        let nim_x = nim.x - b.x * (nim.y - c.y * nim.z / c.z) / b.y - c.x * nim.z / c.z;
        assert!(nim_x < 0.5 * a.x + 1e-9);
        assert!(nim_x >= -0.5 * a.x - 1e-9);
    }

    #[test]
    fn test_truncoct_triclinic_rotmat_forward_backward_are_transposes() {
        // GROMOS: `truncoct_triclinic_rotmat(false) == transpose(truncoct_triclinic_rotmat(true))`,
        // and the matrix is orthogonal, so backward is also the inverse.
        let fwd = truncoct_triclinic_rotmat(true);
        let bwd = truncoct_triclinic_rotmat(false);
        let fwd_t = fwd.transpose();
        for axis in [
            (bwd.x_axis, fwd_t.x_axis),
            (bwd.y_axis, fwd_t.y_axis),
            (bwd.z_axis, fwd_t.z_axis),
        ] {
            assert_relative_eq!(axis.0.x, axis.1.x, epsilon = 1e-12);
            assert_relative_eq!(axis.0.y, axis.1.y, epsilon = 1e-12);
            assert_relative_eq!(axis.0.z, axis.1.z, epsilon = 1e-12);
        }

        let identity = fwd.transpose() * fwd;
        assert_relative_eq!(identity.x_axis.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(identity.y_axis.y, 1.0, epsilon = 1e-12);
        assert_relative_eq!(identity.z_axis.z, 1.0, epsilon = 1e-12);
        assert_relative_eq!(identity.x_axis.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(identity.x_axis.z, 0.0, epsilon = 1e-12);
        assert_relative_eq!(identity.y_axis.z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_truncoct_triclinic_box_round_trip() {
        // forward (cubic -> triclinic) then backward (triclinic -> cubic)
        // recovers the original cube edge length L, per GROMOS
        // `math::truncoct_triclinic_box`.
        let l = 3.767055681;
        let cubic = Mat3::from_diagonal(Vec3::splat(l));
        let triclinic = truncoct_triclinic_box(cubic, true);
        let back = truncoct_triclinic_box(triclinic, false);
        assert_relative_eq!(back.x_axis.x, l, epsilon = 1e-12);
        assert_relative_eq!(back.y_axis.y, l, epsilon = 1e-12);
        assert_relative_eq!(back.z_axis.z, l, epsilon = 1e-12);
    }

    #[test]
    fn test_truncoct_triclinic_rotates_positions() {
        let mut positions = vec![Vec3::new(1.0, 2.0, 3.0)];
        let original = positions[0];

        truncoct_triclinic(&mut positions, true);
        let rotated = positions[0];
        assert!((rotated - original).length() > 1e-9);

        truncoct_triclinic(&mut positions, false);
        assert_relative_eq!(positions[0].x, original.x, epsilon = 1e-12);
        assert_relative_eq!(positions[0].y, original.y, epsilon = 1e-12);
        assert_relative_eq!(positions[0].z, original.z, epsilon = 1e-12);
    }

    #[test]
    fn test_truncoct_triclinic_forward_matches_gromosxx_reference() {
        // GROMOS's `product(rot, v)` (gmath.h) computes `rot^T * v`, not
        // `rot * v` - the forward/backward branches in
        // `truncoct_triclinic_rotmat` must account for that transpose.
        // Reference value from a GROMOS debug build: atom 21 of
        // crates/newen-md/tests/aladip.conf, rotated via
        // `math::truncoct_triclinic(pos, true)` after
        // `math::truncoct_triclinic_box(box, true)` (PLAN 1.4 / FUTURE Dim
        // 11 #1), gives pos(21) = (1.88114910011, 0.267543648603,
        // 0.423126867767).
        let mut positions = vec![Vec3::new(0.867633465, 0.896110299, 1.494502054)];
        truncoct_triclinic(&mut positions, true);
        assert_relative_eq!(positions[0].x, 1.88114910011, epsilon = 1e-9);
        assert_relative_eq!(positions[0].y, 0.267543648603, epsilon = 1e-9);
        assert_relative_eq!(positions[0].z, 0.423126867767, epsilon = 1e-9);
    }

    #[test]
    fn genbox_right_angles_give_the_diagonal_box() {
        let m = genbox_box(Vec3::new(2.0, 2.2, 2.4), Vec3::splat(90.0));
        assert_eq!(m, Mat3::from_diagonal(Vec3::new(2.0, 2.2, 2.4)));
    }

    #[test]
    fn genbox_box_round_trips_through_lengths_and_angles() {
        let (lengths, angles) = (Vec3::new(2.0, 2.2, 2.4), Vec3::new(75.0, 80.0, 65.0));
        let m = genbox_box(lengths, angles);
        // a along x, b in the xy plane
        assert_eq!((m.x_axis.y, m.x_axis.z, m.y_axis.z), (0.0, 0.0, 0.0));
        let (l, a) = genbox_lengths_angles(m);
        assert_relative_eq!(l.x, 2.0, epsilon = 1e-14);
        assert_relative_eq!(l.y, 2.2, epsilon = 1e-14);
        assert_relative_eq!(l.z, 2.4, epsilon = 1e-14);
        assert_relative_eq!(a.x, 75.0, epsilon = 1e-12);
        assert_relative_eq!(a.y, 80.0, epsilon = 1e-12);
        assert_relative_eq!(a.z, 65.0, epsilon = 1e-12);
    }
}
