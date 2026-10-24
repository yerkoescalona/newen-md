//! Generic 3x3 matrix type for molecular dynamics
//!
//! This module provides `Matrix3<F>` - a generic 3x3 matrix that works with
//! both f32 and f64 precision. Used for rotation matrices, stress tensors,
//! box vectors, and other matrix operations in MD.

use crate::float::Float;
use crate::vector::Vector3;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

/// Generic 3x3 matrix for molecular dynamics
///
/// Stored in column-major order (like OpenGL/glam) for compatibility.
/// Each column represents a basis vector.
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Matrix3<F: Float> {
    /// Column 0 (x basis vector)
    pub x_axis: Vector3<F>,
    /// Column 1 (y basis vector)  
    pub y_axis: Vector3<F>,
    /// Column 2 (z basis vector)
    pub z_axis: Vector3<F>,
}

impl<F: Float> Matrix3<F> {
    /// Create matrix from column vectors
    #[inline]
    pub const fn from_cols(x_axis: Vector3<F>, y_axis: Vector3<F>, z_axis: Vector3<F>) -> Self {
        Self {
            x_axis,
            y_axis,
            z_axis,
        }
    }

    /// Create matrix from row values (row-major input, stored column-major)
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub fn from_rows(
        m00: F,
        m01: F,
        m02: F,
        m10: F,
        m11: F,
        m12: F,
        m20: F,
        m21: F,
        m22: F,
    ) -> Self {
        Self {
            x_axis: Vector3::new(m00, m10, m20),
            y_axis: Vector3::new(m01, m11, m21),
            z_axis: Vector3::new(m02, m12, m22),
        }
    }

    /// Zero matrix
    #[inline]
    pub fn zero() -> Self {
        Self::from_cols(Vector3::zero(), Vector3::zero(), Vector3::zero())
    }

    /// Identity matrix
    #[inline]
    pub fn identity() -> Self {
        Self::from_cols(Vector3::unit_x(), Vector3::unit_y(), Vector3::unit_z())
    }

    /// Create diagonal matrix
    #[inline]
    pub fn from_diagonal(diag: Vector3<F>) -> Self {
        Self::from_rows(
            diag.x,
            F::ZERO,
            F::ZERO,
            F::ZERO,
            diag.y,
            F::ZERO,
            F::ZERO,
            F::ZERO,
            diag.z,
        )
    }

    /// Create uniform scale matrix
    #[inline]
    pub fn from_scale(scale: F) -> Self {
        Self::from_diagonal(Vector3::splat(scale))
    }

    /// Create rotation matrix around X axis
    #[inline]
    pub fn from_rotation_x(angle: F) -> Self {
        let (sin, cos) = (angle.sin(), angle.cos());
        Self::from_rows(
            F::ONE,
            F::ZERO,
            F::ZERO,
            F::ZERO,
            cos,
            -sin,
            F::ZERO,
            sin,
            cos,
        )
    }

    /// Create rotation matrix around Y axis
    #[inline]
    pub fn from_rotation_y(angle: F) -> Self {
        let (sin, cos) = (angle.sin(), angle.cos());
        Self::from_rows(
            cos,
            F::ZERO,
            sin,
            F::ZERO,
            F::ONE,
            F::ZERO,
            -sin,
            F::ZERO,
            cos,
        )
    }

    /// Create rotation matrix around Z axis
    #[inline]
    pub fn from_rotation_z(angle: F) -> Self {
        let (sin, cos) = (angle.sin(), angle.cos());
        Self::from_rows(
            cos,
            -sin,
            F::ZERO,
            sin,
            cos,
            F::ZERO,
            F::ZERO,
            F::ZERO,
            F::ONE,
        )
    }

    /// Get element at (row, col)
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> F {
        match col {
            0 => self.x_axis[row],
            1 => self.y_axis[row],
            2 => self.z_axis[row],
            _ => panic!("Matrix3 column out of bounds: {}", col),
        }
    }

    /// Set element at (row, col)
    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: F) {
        match col {
            0 => self.x_axis[row] = value,
            1 => self.y_axis[row] = value,
            2 => self.z_axis[row] = value,
            _ => panic!("Matrix3 column out of bounds: {}", col),
        }
    }

    /// Transpose
    #[inline]
    pub fn transpose(&self) -> Self {
        Self::from_rows(
            self.x_axis.x,
            self.x_axis.y,
            self.x_axis.z,
            self.y_axis.x,
            self.y_axis.y,
            self.y_axis.z,
            self.z_axis.x,
            self.z_axis.y,
            self.z_axis.z,
        )
    }

    /// Determinant
    #[inline]
    pub fn determinant(&self) -> F {
        self.x_axis.x * (self.y_axis.y * self.z_axis.z - self.z_axis.y * self.y_axis.z)
            - self.y_axis.x * (self.x_axis.y * self.z_axis.z - self.z_axis.y * self.x_axis.z)
            + self.z_axis.x * (self.x_axis.y * self.y_axis.z - self.y_axis.y * self.x_axis.z)
    }

    /// Inverse (returns None if singular)
    pub fn inverse(&self) -> Option<Self> {
        let det = self.determinant();
        if det.abs() < F::EPSILON {
            return None;
        }

        let inv_det = F::ONE / det;

        Some(Self::from_rows(
            (self.y_axis.y * self.z_axis.z - self.z_axis.y * self.y_axis.z) * inv_det,
            (self.z_axis.x * self.y_axis.z - self.y_axis.x * self.z_axis.z) * inv_det,
            (self.y_axis.x * self.z_axis.y - self.z_axis.x * self.y_axis.y) * inv_det,
            (self.z_axis.y * self.x_axis.z - self.x_axis.y * self.z_axis.z) * inv_det,
            (self.x_axis.x * self.z_axis.z - self.z_axis.x * self.x_axis.z) * inv_det,
            (self.z_axis.x * self.x_axis.y - self.x_axis.x * self.z_axis.y) * inv_det,
            (self.x_axis.y * self.y_axis.z - self.y_axis.y * self.x_axis.z) * inv_det,
            (self.y_axis.x * self.x_axis.z - self.x_axis.x * self.y_axis.z) * inv_det,
            (self.x_axis.x * self.y_axis.y - self.y_axis.x * self.x_axis.y) * inv_det,
        ))
    }

    /// Trace (sum of diagonal elements)
    #[inline]
    pub fn trace(&self) -> F {
        self.x_axis.x + self.y_axis.y + self.z_axis.z
    }

    /// Get diagonal as vector
    #[inline]
    pub fn diagonal(&self) -> Vector3<F> {
        Vector3::new(self.x_axis.x, self.y_axis.y, self.z_axis.z)
    }

    /// Get column by index
    #[inline]
    pub fn col(&self, index: usize) -> Vector3<F> {
        match index {
            0 => self.x_axis,
            1 => self.y_axis,
            2 => self.z_axis,
            _ => panic!("Matrix3 column out of bounds: {}", index),
        }
    }

    /// Get row by index
    #[inline]
    pub fn row(&self, index: usize) -> Vector3<F> {
        Vector3::new(self.x_axis[index], self.y_axis[index], self.z_axis[index])
    }

    /// Transform a vector (matrix * vector)
    #[inline]
    pub fn mul_vec(&self, v: Vector3<F>) -> Vector3<F> {
        self.x_axis * v.x + self.y_axis * v.y + self.z_axis * v.z
    }

    /// Frobenius norm squared
    #[inline]
    pub fn norm_squared(&self) -> F {
        self.x_axis.length_squared() + self.y_axis.length_squared() + self.z_axis.length_squared()
    }

    /// Frobenius norm
    #[inline]
    pub fn norm(&self) -> F {
        self.norm_squared().sqrt()
    }
}

impl<F: Float> Default for Matrix3<F> {
    fn default() -> Self {
        Self::identity()
    }
}

// Matrix multiplication
impl<F: Float> Mul for Matrix3<F> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::from_cols(
            self.mul_vec(rhs.x_axis),
            self.mul_vec(rhs.y_axis),
            self.mul_vec(rhs.z_axis),
        )
    }
}

impl<F: Float> MulAssign for Matrix3<F> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

// Scalar multiplication
impl<F: Float> Mul<F> for Matrix3<F> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: F) -> Self {
        Self::from_cols(self.x_axis * rhs, self.y_axis * rhs, self.z_axis * rhs)
    }
}

// Matrix-vector multiplication
impl<F: Float> Mul<Vector3<F>> for Matrix3<F> {
    type Output = Vector3<F>;

    #[inline]
    fn mul(self, rhs: Vector3<F>) -> Vector3<F> {
        self.mul_vec(rhs)
    }
}

// Matrix addition
impl<F: Float> Add for Matrix3<F> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::from_cols(
            self.x_axis + rhs.x_axis,
            self.y_axis + rhs.y_axis,
            self.z_axis + rhs.z_axis,
        )
    }
}

impl<F: Float> AddAssign for Matrix3<F> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x_axis += rhs.x_axis;
        self.y_axis += rhs.y_axis;
        self.z_axis += rhs.z_axis;
    }
}

// Matrix subtraction
impl<F: Float> Sub for Matrix3<F> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::from_cols(
            self.x_axis - rhs.x_axis,
            self.y_axis - rhs.y_axis,
            self.z_axis - rhs.z_axis,
        )
    }
}

impl<F: Float> SubAssign for Matrix3<F> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x_axis -= rhs.x_axis;
        self.y_axis -= rhs.y_axis;
        self.z_axis -= rhs.z_axis;
    }
}

// Column indexing
impl<F: Float> Index<usize> for Matrix3<F> {
    type Output = Vector3<F>;

    #[inline]
    fn index(&self, index: usize) -> &Vector3<F> {
        match index {
            0 => &self.x_axis,
            1 => &self.y_axis,
            2 => &self.z_axis,
            _ => panic!("Matrix3 column out of bounds: {}", index),
        }
    }
}

impl<F: Float> IndexMut<usize> for Matrix3<F> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Vector3<F> {
        match index {
            0 => &mut self.x_axis,
            1 => &mut self.y_axis,
            2 => &mut self.z_axis,
            _ => panic!("Matrix3 column out of bounds: {}", index),
        }
    }
}

// Precision conversions
impl Matrix3<f32> {
    /// Convert to f64 precision
    #[inline]
    pub fn to_f64(self) -> Matrix3<f64> {
        Matrix3::from_cols(
            self.x_axis.to_f64(),
            self.y_axis.to_f64(),
            self.z_axis.to_f64(),
        )
    }
}

impl Matrix3<f64> {
    /// Convert to f32 precision
    #[inline]
    pub fn to_f32(self) -> Matrix3<f32> {
        Matrix3::from_cols(
            self.x_axis.to_f32(),
            self.y_axis.to_f32(),
            self.z_axis.to_f32(),
        )
    }
}

/// Type alias for single precision 3x3 matrix
pub type Mat3f = Matrix3<f32>;

/// Type alias for double precision 3x3 matrix
pub type Mat3d = Matrix3<f64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity() {
        let m: Mat3d = Matrix3::identity();
        let v = Vector3::new(1.0, 2.0, 3.0);
        let result = m * v;

        assert!(result.approx_eq(v, 1e-12));
    }

    #[test]
    fn test_determinant() {
        let m: Mat3d = Matrix3::identity();
        assert!((m.determinant() - 1.0).abs() < 1e-12);

        let m2 = Matrix3::from_scale(2.0);
        assert!((m2.determinant() - 8.0).abs() < 1e-12);
    }

    #[test]
    fn test_inverse() {
        let m: Mat3d = Matrix3::from_rows(1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0);

        let inv = m.inverse().unwrap();
        let product: Mat3d = m * inv;
        let identity: Mat3d = Matrix3::identity();

        for i in 0..3 {
            for j in 0..3 {
                assert!((product.get(i, j) - identity.get(i, j)).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_rotation() {
        let m: Mat3d = Matrix3::from_rotation_z(std::f64::consts::FRAC_PI_2);
        let v = Vector3::new(1.0, 0.0, 0.0);
        let result = m * v;

        // Rotating (1,0,0) by 90 degrees around Z should give (0,1,0)
        assert!(result.approx_eq(Vector3::new(0.0, 1.0, 0.0), 1e-12));
    }

    #[test]
    fn test_precision_conversion() {
        let m_f64: Mat3d = Matrix3::from_rotation_z(1.0);
        let m_f32 = m_f64.to_f32();
        let m_back = m_f32.to_f64();

        assert!((m_back.determinant() - m_f64.determinant()).abs() < 1e-5);
    }
}
