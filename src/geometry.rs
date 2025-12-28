use std::f64::consts::PI;
use num_complex::Complex64;
use crate::constants::{S_RADIUS, ALPHA};

/// A point in 4-Dimensional Complex Space (x, y, z, w).
/// 'w' is the spatial dimension required for the Horn Torus closure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spatial4D {
    pub x: Complex64,
    pub y: Complex64,
    pub z: Complex64,
    pub w: Complex64, 
}

impl Spatial4D {
    pub fn new(x: Complex64, y: Complex64, z: Complex64, w: Complex64) -> Self {
        Self { x, y, z, w }
    }

    pub fn zero() -> Self {
        Self::new(Complex64::default(), Complex64::default(), Complex64::default(), Complex64::default())
    }

    /// Geometric magnitude (Complex Euclidean Norm).
    pub fn magnitude(&self) -> Complex64 {
        let sum_sq = (self.x * self.x) + (self.y * self.y) + (self.z * self.z) + (self.w * self.w);
        sum_sq.sqrt()
    }
}

/// The Horn Torus Geometry ($R = r = S$).
/// Source: <https://mathworld.wolfram.com/HornTorus.html>
#[derive(Debug, Clone)]
pub struct HornTorus {
    /// The fundamental radius S (Major and Minor radius are equal).
    pub s: f64,
}

impl Default for HornTorus {
    fn default() -> Self {
        Self { s: S_RADIUS }
    }
}

impl HornTorus {
    pub fn new(s: f64) -> Self {
        Self { s }
    }

    /// Exact Volume: $2 \pi^2 S^3$
    pub fn volume(&self) -> f64 {
        2.0 * PI.powi(2) * self.s.powi(3)
    }

    /// Exact Surface Area: $4 \pi^2 S^2$
    /// (Note: Wolfram MathWorld definition for Horn Torus is 4*pi^2*R*r. Since R=r, it is 4*pi^2*r^2).
    pub fn surface_area(&self) -> f64 {
        4.0 * PI.powi(2) * self.s.powi(2)
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
        let rhs = 4.0 * self.s.powi(2) * (x*x + y*y);
        lhs - rhs
    }

    /// Parametric Equations ($u, v \in [0, 2\pi)$).
    /// Wolfram: $x = (S + S \cos u) \cos v$
    pub fn parametric(&self, u: f64, v: f64) -> (f64, f64, f64) {
        let h = self.s * (1.0 + u.cos());
        (
            h * v.cos(),
            h * v.sin(),
            self.s * u.sin()
        )
    }
    
    /// Gaussian Curvature $K$.
    /// Formula: $K = \frac{\cos u}{r(R + r \cos u)}$ -> $K = \frac{\cos u}{S^2(1+\cos u)}$
    pub fn gaussian_curvature(&self, u: f64) -> f64 {
        let denom = self.s.powi(2) * (1.0 + u.cos());
        if denom.abs() < 1e-15 { return f64::INFINITY; }
        u.cos() / denom
    }
}