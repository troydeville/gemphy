use std::f64::consts::PI;

use num_complex::Complex64;

use crate::{geometry::Spatial4D, medium::GeometricEncodedMedium};
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
            mass: 0.0,
            topology: 0.0,
            geometric_radius_a: r, 
            base_length: r,
        }
    }
}

impl GeometricKnot {
    pub fn new(med: GeometricEncodedMedium, mass: f64, sub_topologies: &[f64], phys_radius: f64, name: &str) -> Self {
        if mass <= 0.0 { return Self::default(); }

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