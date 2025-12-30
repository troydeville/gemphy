use std::f64::consts::PI;

use crate::{calculate_mass_frequency, medium, s_constant};
use crate::geometry::Spatial4D;
use crate::knot::GeometricKnot;
use crate::medium::{ALPHA, C, GemInteractionResult}; 
use num_complex::{Complex64, ComplexFloat};

#[derive(Debug, Clone)]
pub struct DynamicKnot {
    pub knot: GeometricKnot,
    pub position: Spatial4D,
    pub velocity: Spatial4D,
    pub acceleration: Spatial4D,
    pub force_accumulator: Spatial4D, 
}

impl DynamicKnot {
    pub fn new(knot: GeometricKnot, position: Spatial4D) -> Self {
        Self {
            knot,
            position,
            velocity: Spatial4D::zero(),
            acceleration: Spatial4D::zero(),
            force_accumulator: Spatial4D::zero(),
        }
    }

    pub fn distance(&self, knot: &DynamicKnot) -> Spatial4D {
        Spatial4D { 
            x: knot.position.x - self.position.x,
            y: knot.position.y - self.position.y,
            z: knot.position.z - self.position.z,
            w: knot.position.w - self.position.w
        }
    }

    pub fn clear_forces(&mut self) {
        self.force_accumulator = Spatial4D::zero();
    }

    pub fn integrate(&mut self, dt: f64) {
        let inv_mass = 1.0 / self.knot.mass;
        
        // 1. Calculate Acceleration (Complex)
        self.acceleration = self.force_accumulator * inv_mass;
        
        // 2. Symplectic Euler Integration
        self.velocity = self.velocity + (self.acceleration * dt);
        self.position = self.position + (self.velocity * dt);
    }

    pub fn apply_interaction(&mut self, result: &GemInteractionResult, ref_knot: &DynamicKnot) {
        let displacement = ref_knot.position - self.position;
        let dist_val = displacement.magnitude().norm(); 
        
        if dist_val < 1.0e-30 { return; }

        let direction = displacement.normalize();

        // FIX: Use the calculated force magnitude directly from medium.rs.
        let force_newtonian = result.force.norm();

        // Radial Component (Attraction towards ref_knot)
        let force_radial = direction * force_newtonian;

        // Tangential Component (Torque from Phase)
        let twist_ratio = if result.force.re.abs() > 1.0e-30 {
            result.force.im / result.force.re
        } else {
            0.0
        };
        
        let up_vector = if direction.x.norm() < 1e-9 && direction.y.norm() < 1e-9 {
             Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0), Complex64::new(0.0,0.0))
        } else {
             Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0))
        };
        
        let tangent = direction.cross_3d(&up_vector).normalize();
        let force_tangential = tangent * (force_newtonian * twist_ratio.abs());

        // Accumulate Total Force
        self.force_accumulator = self.force_accumulator + force_radial + force_tangential;
        
        // Removed internal integrate() call. 
        // Integration is now controlled by the main simulation loop.
    }
}