use std::{f64::consts::{PI, SQRT_2}, ops::{Add, Mul, Sub}};

use num_complex::Complex64;

use crate::medium::ALPHA;

/// A point in 4-Dimensional Complex Space (x, y, z, w).
/// 'w' is the spatial dimension required for the Horn Torus closure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spatial4D {
    pub x: Complex64,
    pub y: Complex64,
    pub z: Complex64,
    pub w: Complex64, 
}

use num_traits::{ConstOne, ConstZero, Inv, MulAdd, Num, One, Pow, Signed, Zero};

/// Alias for a [`Complex<f64>`]
impl Zero for Spatial4D {
    
    fn zero() -> Self {
        Self::zero()
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero() && self.z.is_zero() && self.w.is_zero()
    }
}

impl Spatial4D {
    pub fn new(x: Complex64, y: Complex64, z: Complex64, w: Complex64) -> Self {
        Self { x, y, z, w }
    }

    pub fn zero() -> Self {
        Self::new(Complex64::default(), Complex64::default(), Complex64::default(), Complex64::default())
    }

    pub fn one() -> Self {
        Self {
            x: Complex64::one(),
            y: Complex64::one(),
            z: Complex64::one(),
            w: Complex64::one(),
        }
    }

    /// Geometric magnitude (Complex Euclidean Norm).
    pub fn magnitude(&self) -> Complex64 {
        let sum_sq = (self.x * self.x) + (self.y * self.y) + (self.z * self.z) + (self.w * self.w);
        sum_sq.sqrt()
    }

    /// Dot Product: Calculates scalar projection (Energy/Work).
    /// A . B = (xa*xb + ya*yb + za*zb + wa*wb)
    pub fn dot(&self, other: &Self) -> Complex64 {
        (self.x * other.x) + 
        (self.y * other.y) + 
        (self.z * other.z) + 
        (self.w * other.w)
    }

    /// 3D Cross Product (ignoring W).
    /// Essential for Orbital Mechanics (r x v).
    /// The 'w' component is usually zeroed or treated as a separate phase interaction.
    pub fn cross_3d(&self, other: &Self) -> Self {
        Self {
            x: (self.y * other.z) - (self.z * other.y),
            y: (self.z * other.x) - (self.x * other.z),
            z: (self.x * other.y) - (self.y * other.x),
            w: Complex64::default(), // Cross product doesn't produce 'W' torque
        }
    }

    /// Normalizes the vector (Direction only).
    pub fn normalize(&self) -> Self {
        let mag = self.magnitude(); // Complex Magnitude
        // Avoid division by zero
        if mag.norm() < 1e-30 {
            return Self::zero();
        }
        
        // We divide by the magnitude
        Self {
            x: self.x / mag,
            y: self.y / mag,
            z: self.z / mag,
            w: self.w / mag,
        }
    }

}

// Vector Addition: a + b
impl Add for Spatial4D {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w,
        }
    }
}

// Vector Subtraction: a - b
impl Sub for Spatial4D {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w,
        }
    }
}

// Scalar Multiplication: a * f64
impl Mul<f64> for Spatial4D {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
            w: self.w * scalar,
        }
    }
}

// Complex Scalar Multiplication: a * Complex64 (for phase/force scaling)
impl Mul<Complex64> for Spatial4D {
    type Output = Self;
    fn mul(self, scalar: Complex64) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
            w: self.w * scalar,
        }
    }
}

/// The Horn Torus Geometry ($R = r = S$).
/// Source: <https://mathworld.wolfram.com/HornTorus.html>
#[derive(Debug, Clone)]
struct HornTorus {
    /// The fundamental radius S (Major and Minor radius are equal).
    pub r: f64,
}

impl Default for HornTorus {
    fn default() -> Self {
        let s = SQRT_2 * PI.powf(0.25);
        Self { r: s }
    }
}

impl HornTorus {
    pub fn new(r: f64) -> Self {
        Self { r }
    }

    /// Exact Volume: $2 \pi^2 S^3$
    pub fn volume(&self) -> f64 {
        2.0 * PI.powi(2) * self.r.powi(3)
    }

    /// Exact Surface Area: $4 \pi^2 S^2$
    /// (Note: Wolfram MathWorld definition for Horn Torus is 4*pi^2*R*r. Since R=r, it is 4*pi^2*r^2).
    pub fn surface_area(&self) -> f64 {
        4.0 * PI.powi(2) * self.r.powi(2)
    }

    /// The "Volume Mismatch" driving Action flow.
    pub fn volume_mismatch(&self) -> f64 {
        self.volume() * ALPHA
    }

    /// Implicit Equation evaluation.
    /// Formula: $(x^2 + y^2 + z^2)^2 = 4S^2(x^2 + y^2)$
    /// Returns the residual (should be 0 on surface).
    pub fn implicit_residual(&self, x: f64, y: f64, z: f64) -> f64 {
        let lhs = (x*x + y*y + z*z).powi(2);
        let rhs = 4.0 * self.r.powi(2) * (x*x + y*y);
        lhs - rhs
    }

    /// Parametric Equations ($u, v \in [0, 2\pi)$).
    /// Wolfram: $x = (S + S \cos u) \cos v$
    pub fn parametric(&self, u: f64, v: f64) -> (f64, f64, f64) {
        let h = self.r * (1.0 + u.cos());
        (
            h * v.cos(),
            h * v.sin(),
            self.r * u.sin()
        )
    }
    
    /// Gaussian Curvature $K$.
    /// Formula: $K = \frac{\cos u}{r(R + r \cos u)}$ -> $K = \frac{\cos u}{S^2(1+\cos u)}$
    pub fn gaussian_curvature(&self, u: f64) -> f64 {
        let denom = self.r.powi(2) * (1.0 + u.cos());
        if denom.abs() < 1e-15 { return f64::INFINITY; }
        u.cos() / denom
    }
}

pub trait GemSurface {
    fn radius_a(&self) -> f64;
    fn volume(&self) -> f64;
    fn surface_area(&self) -> f64;
    fn parametric_surface(&self, u: f64, v: f64) -> [f64; 3];
    fn implicit_equation(&self, x: f64, y: f64, z: f64) -> f64;
    fn metric_tensor(&self, v: f64) -> (f64, f64) {
        let a = self.radius_a();
        let cos_half_v = (v / 2.0).cos();
        let g_uu = 4.0 * a.powi(2) * cos_half_v.powi(4);
        let g_vv = a.powi(2);
        (g_uu, g_vv)
    }
    fn gaussian_curvature(&self, v: f64) -> f64 {
        let a = self.radius_a();
        let denom = a.powi(2) * (1.0 + v.cos());
        if denom.abs() < 1e-9 { return 0.0; }
        v.cos() / denom
    }
}

impl GemSurface for HornTorus {
    fn radius_a(&self) -> f64 { self.r }
    fn volume(&self) -> f64 { 2.0 * PI.powi(2) * self.r.powi(3) }
    fn surface_area(&self) -> f64 { 4.0 * PI.powi(2) * self.r.powi(2) }
    fn parametric_surface(&self, u: f64, v: f64) -> [f64; 3] {
        let tube_factor = 1.0 + v.cos();
        [self.r * u.cos() * tube_factor, self.r * tube_factor * u.sin(), self.r * v.sin()]
    }
    fn implicit_equation(&self, x: f64, y: f64, z: f64) -> f64 {
        let sum_sq = x*x + y*y + z*z;
        sum_sq.powi(2) - (4.0 * self.r.powi(2) * (x*x + y*y))
    }
}