//! Floating-point precision trait for generic MD calculations
//!
//! This module provides a `Float` trait that abstracts over f32 and f64,
//! allowing users to choose their desired precision at compile time.
//!
//! # Usage
//!
//! ```rust
//! use newen_core::{Float, Vec3};
//!
//! // Create a position vector
//! let pos: Vec3 = Vec3::new(1.0, 2.0, 3.0);
//!
//! // The Float trait can be used for generic calculations
//! fn compute<T: Float>(a: T, b: T) -> T {
//!     a + b
//! }
//! let result = compute(1.0_f64, 2.0_f64);
//! ```

use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Trait for floating-point types used in MD calculations
///
/// This trait provides all the operations needed for molecular dynamics,
/// abstracting over f32 and f64 to allow compile-time precision selection.
pub trait Float:
    Copy
    + Clone
    + Debug
    + Display
    + Default
    + PartialEq
    + PartialOrd
    + Send
    + Sync
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
{
    /// The constant π
    const PI: Self;
    /// The constant 2π
    const TWO_PI: Self;
    /// The constant π/2
    const FRAC_PI_2: Self;
    /// The machine epsilon
    const EPSILON: Self;
    /// Zero
    const ZERO: Self;
    /// One
    const ONE: Self;
    /// Two
    const TWO: Self;
    /// Half
    const HALF: Self;

    /// Create from f64 (for constants and initialization)
    fn from_f64(v: f64) -> Self;

    /// Convert to f64 (for output and interoperability)
    fn to_f64(self) -> f64;

    /// Create from f32
    fn from_f32(v: f32) -> Self;

    /// Convert to f32
    fn to_f32(self) -> f32;

    /// Create from usize
    fn from_usize(v: usize) -> Self;

    /// Create from i32
    fn from_i32(v: i32) -> Self;

    /// Square root
    fn sqrt(self) -> Self;

    /// Sine
    fn sin(self) -> Self;

    /// Cosine
    fn cos(self) -> Self;

    /// Tangent
    fn tan(self) -> Self;

    /// Arc cosine
    fn acos(self) -> Self;

    /// Arc sine
    fn asin(self) -> Self;

    /// Arc tangent of y/x
    fn atan2(self, other: Self) -> Self;

    /// Exponential
    fn exp(self) -> Self;

    /// Natural logarithm
    fn ln(self) -> Self;

    /// Power
    fn powf(self, n: Self) -> Self;

    /// Integer power
    fn powi(self, n: i32) -> Self;

    /// Absolute value
    fn abs(self) -> Self;

    /// Floor
    fn floor(self) -> Self;

    /// Ceiling
    fn ceil(self) -> Self;

    /// Round
    fn round(self) -> Self;

    /// Minimum
    fn min(self, other: Self) -> Self;

    /// Maximum
    fn max(self, other: Self) -> Self;

    /// Clamp between min and max
    fn clamp(self, min: Self, max: Self) -> Self {
        self.max(min).min(max)
    }

    /// Check if approximately equal within epsilon
    fn approx_eq(self, other: Self, epsilon: Self) -> bool {
        (self - other).abs() < epsilon
    }

    /// Reciprocal (1/x)
    fn recip(self) -> Self {
        Self::ONE / self
    }

    /// Square (x^2)
    fn squared(self) -> Self {
        self * self
    }

    /// Cube (x^3)
    fn cubed(self) -> Self {
        self * self * self
    }

    /// Fused multiply-add: self * a + b
    fn mul_add(self, a: Self, b: Self) -> Self;

    /// Hyperbolic sine
    fn sinh(self) -> Self;

    /// Hyperbolic cosine
    fn cosh(self) -> Self;

    /// Hyperbolic tangent
    fn tanh(self) -> Self;

    /// Check if NaN
    fn is_nan(self) -> bool;

    /// Check if infinite
    fn is_infinite(self) -> bool;

    /// Check if finite
    fn is_finite(self) -> bool;

    /// Copysign: returns a value with magnitude of self and sign of other
    fn copysign(self, other: Self) -> Self;
}

impl Float for f32 {
    const PI: Self = std::f32::consts::PI;
    const TWO_PI: Self = 2.0 * std::f32::consts::PI;
    const FRAC_PI_2: Self = std::f32::consts::FRAC_PI_2;
    const EPSILON: Self = f32::EPSILON;
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const TWO: Self = 2.0;
    const HALF: Self = 0.5;

    #[inline]
    fn from_f64(v: f64) -> Self {
        v as f32
    }

    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }

    #[inline]
    fn from_f32(v: f32) -> Self {
        v
    }

    #[inline]
    fn to_f32(self) -> f32 {
        self
    }

    #[inline]
    fn from_usize(v: usize) -> Self {
        v as f32
    }

    #[inline]
    fn from_i32(v: i32) -> Self {
        v as f32
    }

    #[inline]
    fn sqrt(self) -> Self {
        f32::sqrt(self)
    }

    #[inline]
    fn sin(self) -> Self {
        f32::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        f32::cos(self)
    }

    #[inline]
    fn tan(self) -> Self {
        f32::tan(self)
    }

    #[inline]
    fn acos(self) -> Self {
        f32::acos(self)
    }

    #[inline]
    fn asin(self) -> Self {
        f32::asin(self)
    }

    #[inline]
    fn atan2(self, other: Self) -> Self {
        f32::atan2(self, other)
    }

    #[inline]
    fn exp(self) -> Self {
        f32::exp(self)
    }

    #[inline]
    fn ln(self) -> Self {
        f32::ln(self)
    }

    #[inline]
    fn powf(self, n: Self) -> Self {
        f32::powf(self, n)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        f32::powi(self, n)
    }

    #[inline]
    fn abs(self) -> Self {
        f32::abs(self)
    }

    #[inline]
    fn floor(self) -> Self {
        f32::floor(self)
    }

    #[inline]
    fn ceil(self) -> Self {
        f32::ceil(self)
    }

    #[inline]
    fn round(self) -> Self {
        f32::round(self)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        f32::min(self, other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        f32::max(self, other)
    }

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self {
        f32::mul_add(self, a, b)
    }

    #[inline]
    fn sinh(self) -> Self {
        f32::sinh(self)
    }

    #[inline]
    fn cosh(self) -> Self {
        f32::cosh(self)
    }

    #[inline]
    fn tanh(self) -> Self {
        f32::tanh(self)
    }

    #[inline]
    fn is_nan(self) -> bool {
        f32::is_nan(self)
    }

    #[inline]
    fn is_infinite(self) -> bool {
        f32::is_infinite(self)
    }

    #[inline]
    fn is_finite(self) -> bool {
        f32::is_finite(self)
    }

    #[inline]
    fn copysign(self, other: Self) -> Self {
        f32::copysign(self, other)
    }
}

impl Float for f64 {
    const PI: Self = std::f64::consts::PI;
    const TWO_PI: Self = 2.0 * std::f64::consts::PI;
    const FRAC_PI_2: Self = std::f64::consts::FRAC_PI_2;
    const EPSILON: Self = f64::EPSILON;
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const TWO: Self = 2.0;
    const HALF: Self = 0.5;

    #[inline]
    fn from_f64(v: f64) -> Self {
        v
    }

    #[inline]
    fn to_f64(self) -> f64 {
        self
    }

    #[inline]
    fn from_f32(v: f32) -> Self {
        v as f64
    }

    #[inline]
    fn to_f32(self) -> f32 {
        self as f32
    }

    #[inline]
    fn from_usize(v: usize) -> Self {
        v as f64
    }

    #[inline]
    fn from_i32(v: i32) -> Self {
        v as f64
    }

    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }

    #[inline]
    fn sin(self) -> Self {
        f64::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        f64::cos(self)
    }

    #[inline]
    fn tan(self) -> Self {
        f64::tan(self)
    }

    #[inline]
    fn acos(self) -> Self {
        f64::acos(self)
    }

    #[inline]
    fn asin(self) -> Self {
        f64::asin(self)
    }

    #[inline]
    fn atan2(self, other: Self) -> Self {
        f64::atan2(self, other)
    }

    #[inline]
    fn exp(self) -> Self {
        f64::exp(self)
    }

    #[inline]
    fn ln(self) -> Self {
        f64::ln(self)
    }

    #[inline]
    fn powf(self, n: Self) -> Self {
        f64::powf(self, n)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        f64::powi(self, n)
    }

    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }

    #[inline]
    fn floor(self) -> Self {
        f64::floor(self)
    }

    #[inline]
    fn ceil(self) -> Self {
        f64::ceil(self)
    }

    #[inline]
    fn round(self) -> Self {
        f64::round(self)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        f64::min(self, other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        f64::max(self, other)
    }

    #[inline]
    fn mul_add(self, a: Self, b: Self) -> Self {
        f64::mul_add(self, a, b)
    }

    #[inline]
    fn sinh(self) -> Self {
        f64::sinh(self)
    }

    #[inline]
    fn cosh(self) -> Self {
        f64::cosh(self)
    }

    #[inline]
    fn tanh(self) -> Self {
        f64::tanh(self)
    }

    #[inline]
    fn is_nan(self) -> bool {
        f64::is_nan(self)
    }

    #[inline]
    fn is_infinite(self) -> bool {
        f64::is_infinite(self)
    }

    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }

    #[inline]
    fn copysign(self, other: Self) -> Self {
        f64::copysign(self, other)
    }
}

/// Type alias for single precision (default for GPU)
pub type F32 = f32;

/// Type alias for double precision (default for CPU MD)
pub type F64 = f64;

#[cfg(test)]
mod tests {
    use super::*;

    fn test_float_ops<F: Float>() {
        let a = F::from_f64(2.0);
        let b = F::from_f64(3.0);

        assert!((a + b).to_f64() - 5.0 < F::EPSILON.to_f64());
        assert!((a * b).to_f64() - 6.0 < F::EPSILON.to_f64());
        assert!((a.sqrt().to_f64() - std::f64::consts::SQRT_2).abs() < 0.0001);
    }

    #[test]
    fn test_f32_ops() {
        test_float_ops::<f32>();
    }

    #[test]
    fn test_f64_ops() {
        test_float_ops::<f64>();
    }

    #[test]
    fn test_constants() {
        assert!((f32::PI - std::f32::consts::PI).abs() < f32::EPSILON);
        assert!((f64::PI - std::f64::consts::PI).abs() < f64::EPSILON);
    }
}
