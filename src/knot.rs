use std::f64::consts::PI;

use crate::{geometry::{GemSurface, Spatial4D}, medium::{ALPHA, GAMMA_P, GeometricEncodedMedium}};
// use crate::{GeometricEncodedMedium, Spatial4D};

#[derive(Debug, Clone)]
pub struct GeometricKnot {
    pub name: String,
    pub mass: f64,
    // This represents the "Net Twist" or winding number.
    // -1.0 = Electron
    // +1.0 = Proton
    //  0.0 = Neutron
    pub topology: f64, 
    
    pub geometric_radius_a: f64,
    pub base_length: f64,
}

impl Default for GeometricKnot {
    fn default() -> Self {

        let r = (4.0 * PI).powf(0.25);
        
        Self {
            name: "Vacuum".to_string(),
            mass: GAMMA_P / r,
            topology: 0.0,
            geometric_radius_a: r / ALPHA.powi(2), 
            base_length: r,
        }
    }
}

impl GeometricKnot {
    pub fn new(med: GeometricEncodedMedium, mass: f64, sub_topologies: &[f64], phys_radius: f64, name: &str) -> Self {
        if mass <= 0.0 { return Self::default(); }

        // Calculate the Base Length (Wavelength)
        // l = Gamma_p / (m * alpha)
        let base_len = GAMMA_P / (mass * ALPHA);

        let distribution = if mass < med.m_p {
            med.gamma / (mass * med.alpha.powi(2))
        } else {
            (2.0 * med.g * mass) / med.c.powi(2)
        };
        let ref_length = if phys_radius > 0.0 { phys_radius } else { distribution * med.alpha.powi(2) };
        let net_topology: f64 = sub_topologies.iter().sum();
        GeometricKnot {
            name: name.to_string(),
            mass,
            topology: net_topology,
            geometric_radius_a: distribution,
            base_length: base_len,
        }
    }
}

impl GemSurface for GeometricKnot {
    fn radius_a(&self) -> f64 {
        self.geometric_radius_a
    }

    // It uses the Horn Torus formulas defined in geometry.rs
    fn volume(&self) -> f64 {
        // V = 2 * pi^2 * r^3
        2.0 * std::f64::consts::PI.powi(2) * self.geometric_radius_a.powi(3)
    }

    fn surface_area(&self) -> f64 {
        // A = 4 * pi^2 * r^2
        // This confirms your intuition: The "Light Container" surface 
        // IS the surface of the particle.
        4.0 * std::f64::consts::PI.powi(2) * self.geometric_radius_a.powi(2)
    }

    fn parametric_surface(&self, u: f64, v: f64) -> [f64; 3] {
        // Delegate to the math model logic, but using the Knot's radius
        let r = self.geometric_radius_a;
        let tube_factor = 1.0 + v.cos();
        [r * u.cos() * tube_factor, r * tube_factor * u.sin(), r * v.sin()]
    }

    fn implicit_equation(&self, x: f64, y: f64, z: f64) -> f64 {
        let r = self.geometric_radius_a;
        let sum_sq = x*x + y*y + z*z;
        sum_sq.powi(2) - (4.0 * r.powi(2) * (x*x + y*y))
    }
}

#[derive(Debug, Clone)]
pub struct ImpedanceField {
    pub id: usize,
    pub core: GeometricKnot, 
    pub position: Spatial4D,      
    pub velocity: Spatial4D,      
    pub phase_accum: f64,        
}

// CRITICAL FIX: Manual initialization prevents Stack Overflow
impl Default for ImpedanceField {
    fn default() -> Self {
        Self {
            id: 0,
            core: GeometricKnot::default(),
            position: Spatial4D::zero(),
            velocity: Spatial4D::zero(),
            phase_accum: 0.0,
        }
    }
}

impl ImpedanceField {
    pub fn new(id: usize, core: GeometricKnot, pos: Spatial4D, vel: Spatial4D) -> Self {
        Self {
            id,
            core,
            position: pos,
            velocity: vel,
            phase_accum: 0.0,
        }
    }
}