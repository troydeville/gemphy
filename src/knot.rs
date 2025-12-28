use num_complex::Complex64;
use crate::{GeometricEncodedMedium, Spatial4D};

#[derive(Debug, Clone)]
pub struct GeometricKnot {
    pub name: String,
    pub mass: f64,
    pub charge: Complex64,
    pub geometric_radius_a: f64, 
    pub base_length: f64,
}

impl Default for GeometricKnot {
    fn default() -> Self {
        Self {
            name: "Vacuum".to_string(),
            mass: 0.0,
            charge: Complex64::default(),
            geometric_radius_a: 0.0,
            base_length: 0.0,
        }
    }
}

impl GeometricKnot {
    pub fn new(med: GeometricEncodedMedium, mass: f64, charge: Complex64, phys_radius: f64, name: &str) -> Self {
        if mass <= 0.0 { return Self::default(); }

        let distribution = if mass < med.m_p {
            med.gamma / (mass * med.alpha.powi(2))
        } else {
            (2.0 * med.g * mass) / med.c.powi(2)
        };
        let ref_length = if phys_radius > 0.0 { phys_radius } else { distribution * med.alpha.powi(2) };

        GeometricKnot {
            name: name.to_string(),
            mass,
            charge,
            geometric_radius_a: distribution,
            base_length: ref_length,
        }
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