//! Generic 3D vector type for molecular dynamics
//!
//! This module provides `Vector3<F>` - a generic 3D vector that works with
//! both f32 and f64 precision. For SIMD-accelerated operations, see the
//! `simd` module which uses glam types.

use crate::float::Float;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// Generic 3D vector for molecular dynamics
///
/// This type is generic over the floating-point precision (f32 or f64).
/// For most MD simulations, f64 is recommended for energy conservation.
/// Use f32 for GPU calculations or when memory bandwidth is critical.
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Vector3<F: Float> {
    /// X component.
    pub x: F,
    /// Y component.
    pub y: F,
    /// Z component.
    pub z: F,
}

impl<F: Float> Vector3<F> {
    /// Create a new vector
    #[inline]
    pub const fn new(x: F, y: F, z: F) -> Self {
        Self { x, y, z }
    }

    /// Zero vector
    #[inline]
    pub fn zero() -> Self {
        Self::new(F::ZERO, F::ZERO, F::ZERO)
    }

    /// Unit vector along X axis
    #[inline]
    pub fn unit_x() -> Self {
        Self::new(F::ONE, F::ZERO, F::ZERO)
    }

    /// Unit vector along Y axis
    #[inline]
    pub fn unit_y() -> Self {
        Self::new(F::ZERO, F::ONE, F::ZERO)
    }

    /// Unit vector along Z axis
    #[inline]
    pub fn unit_z() -> Self {
        Self::new(F::ZERO, F::ZERO, F::ONE)
    }

    /// Create vector with all components equal
    #[inline]
    pub fn splat(v: F) -> Self {
        Self::new(v, v, v)
    }

    /// Create from f64 values (convenience for initialization)
    #[inline]
    pub fn from_f64(x: f64, y: f64, z: f64) -> Self {
        Self::new(F::from_f64(x), F::from_f64(y), F::from_f64(z))
    }

    /// Dot product
    #[inline]
    pub fn dot(self, other: Self) -> F {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product
    #[inline]
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Squared length (magnitude squared)
    #[inline]
    pub fn length_squared(self) -> F {
        self.dot(self)
    }

    /// Length (magnitude)
    #[inline]
    pub fn length(self) -> F {
        self.length_squared().sqrt()
    }

    /// Normalize to unit vector
    #[inline]
    pub fn normalize(self) -> Self {
        let len = self.length();
        if len > F::EPSILON {
            self / len
        } else {
            Self::zero()
        }
    }

    /// Normalize, returning None if length is zero
    #[inline]
    pub fn try_normalize(self) -> Option<Self> {
        let len = self.length();
        if len > F::EPSILON {
            Some(self / len)
        } else {
            None
        }
    }

    /// Distance to another vector
    #[inline]
    pub fn distance(self, other: Self) -> F {
        (self - other).length()
    }

    /// Squared distance to another vector
    #[inline]
    pub fn distance_squared(self, other: Self) -> F {
        (self - other).length_squared()
    }

    /// Linear interpolation
    #[inline]
    pub fn lerp(self, other: Self, t: F) -> Self {
        self + (other - self) * t
    }

    /// Component-wise minimum
    #[inline]
    pub fn min(self, other: Self) -> Self {
        Self::new(
            self.x.min(other.x),
            self.y.min(other.y),
            self.z.min(other.z),
        )
    }

    /// Component-wise maximum
    #[inline]
    pub fn max(self, other: Self) -> Self {
        Self::new(
            self.x.max(other.x),
            self.y.max(other.y),
            self.z.max(other.z),
        )
    }

    /// Component-wise absolute value
    #[inline]
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// Component-wise floor
    #[inline]
    pub fn floor(self) -> Self {
        Self::new(self.x.floor(), self.y.floor(), self.z.floor())
    }

    /// Check if approximately equal
    #[inline]
    pub fn approx_eq(self, other: Self, epsilon: F) -> bool {
        (self.x - other.x).abs() < epsilon
            && (self.y - other.y).abs() < epsilon
            && (self.z - other.z).abs() < epsilon
    }

    /// Convert to array
    #[inline]
    pub fn to_array(self) -> [F; 3] {
        [self.x, self.y, self.z]
    }

    /// Create from array
    #[inline]
    pub fn from_array(arr: [F; 3]) -> Self {
        Self::new(arr[0], arr[1], arr[2])
    }

    /// Convert to tuple
    #[inline]
    pub fn to_tuple(self) -> (F, F, F) {
        (self.x, self.y, self.z)
    }

    /// Sum of components
    #[inline]
    pub fn sum(self) -> F {
        self.x + self.y + self.z
    }

    /// Product of components
    #[inline]
    pub fn product(self) -> F {
        self.x * self.y * self.z
    }

    /// Minimum component
    #[inline]
    pub fn min_element(self) -> F {
        self.x.min(self.y).min(self.z)
    }

    /// Maximum component
    #[inline]
    pub fn max_element(self) -> F {
        self.x.max(self.y).max(self.z)
    }
}

impl<F: Float> Default for Vector3<F> {
    fn default() -> Self {
        Self::zero()
    }
}

// Arithmetic operations

impl<F: Float> Add for Vector3<F> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<F: Float> AddAssign for Vector3<F> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<F: Float> Sub for Vector3<F> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<F: Float> SubAssign for Vector3<F> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl<F: Float> Mul<F> for Vector3<F> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: F) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<F: Float> MulAssign<F> for Vector3<F> {
    #[inline]
    fn mul_assign(&mut self, rhs: F) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<F: Float> Div<F> for Vector3<F> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: F) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl<F: Float> DivAssign<F> for Vector3<F> {
    #[inline]
    fn div_assign(&mut self, rhs: F) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl<F: Float> Neg for Vector3<F> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

// Component-wise multiplication
impl<F: Float> Mul<Vector3<F>> for Vector3<F> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::new(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)
    }
}

// Index access
impl<F: Float> Index<usize> for Vector3<F> {
    type Output = F;

    #[inline]
    fn index(&self, index: usize) -> &F {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vector3 index out of bounds: {}", index),
        }
    }
}

impl<F: Float> IndexMut<usize> for Vector3<F> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut F {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Vector3 index out of bounds: {}", index),
        }
    }
}

// Conversions between precisions
impl Vector3<f32> {
    /// Convert to f64 precision
    #[inline]
    pub fn to_f64(self) -> Vector3<f64> {
        Vector3::new(self.x as f64, self.y as f64, self.z as f64)
    }
}

impl Vector3<f64> {
    /// Convert to f32 precision
    #[inline]
    pub fn to_f32(self) -> Vector3<f32> {
        Vector3::new(self.x as f32, self.y as f32, self.z as f32)
    }
}

/// Type alias for single precision 3D vector
pub type Vec3f = Vector3<f32>;

/// Type alias for double precision 3D vector
pub type Vec3d = Vector3<f64>;

/// ZERO constant for Vec3f
pub const VEC3F_ZERO: Vec3f = Vector3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};

/// ZERO constant for Vec3d
pub const VEC3D_ZERO: Vec3d = Vector3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_ops_f32() {
        let a: Vec3f = Vector3::new(1.0, 2.0, 3.0);
        let b: Vec3f = Vector3::new(4.0, 5.0, 6.0);

        let c = a + b;
        assert!((c.x - 5.0).abs() < 1e-6);
        assert!((c.y - 7.0).abs() < 1e-6);
        assert!((c.z - 9.0).abs() < 1e-6);
    }

    #[test]
    fn test_basic_ops_f64() {
        let a: Vec3d = Vector3::new(1.0, 2.0, 3.0);
        let b: Vec3d = Vector3::new(4.0, 5.0, 6.0);

        let c = a + b;
        assert!((c.x - 5.0).abs() < 1e-12);
        assert!((c.y - 7.0).abs() < 1e-12);
        assert!((c.z - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_dot_product() {
        let a: Vec3d = Vector3::new(1.0, 2.0, 3.0);
        let b: Vec3d = Vector3::new(4.0, 5.0, 6.0);

        // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
        assert!((a.dot(b) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_cross_product() {
        let x: Vec3d = Vector3::unit_x();
        let y: Vec3d = Vector3::unit_y();
        let z = x.cross(y);

        assert!(z.approx_eq(Vector3::unit_z(), 1e-12));
    }

    #[test]
    fn test_normalize() {
        let v: Vec3d = Vector3::new(3.0, 4.0, 0.0);
        let n = v.normalize();

        assert!((n.length() - 1.0).abs() < 1e-12);
        assert!((n.x - 0.6).abs() < 1e-12);
        assert!((n.y - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_precision_conversion() {
        let v_f64: Vec3d = Vector3::new(1.0, 2.0, 3.0);
        let v_f32 = v_f64.to_f32();
        let v_back = v_f32.to_f64();

        assert!(v_back.approx_eq(v_f64, 1e-6));
    }
}
