use num_complex::Complex64;

use crate::GeometricEncodedMedium;

#[derive(Debug, Clone)]
pub struct GeometricKnot {
    pub name: String,
    pub mass: f64,
    pub charge: Complex64,
    pub geometric_radius_a: f64, 
    pub physical_radius: f64,
}

impl GeometricKnot {
    pub fn new(med: GeometricEncodedMedium, mass: f64, charge: Complex64, phys_radius: f64, name: &str) -> Self {
        let distribution = if mass < med.m_p {
            med.gamma / (mass * med.alpha.powi(2))
        } else {
            (2.0 * med.g * mass) / med.c.powi(2)
        };
        let r_phys = if phys_radius > 0.0 { phys_radius } else { distribution * med.alpha.powi(2) };

        GeometricKnot {
            name: name.to_string(),
            mass,
            charge,
            geometric_radius_a: distribution,
            physical_radius: r_phys,
        }
    }
}