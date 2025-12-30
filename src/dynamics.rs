use crate::geometry::Spatial4D;
use crate::knot::GeometricKnot;
use crate::medium::GemInteractionResult; // Import Result Struct
use num_complex::Complex64;

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

    pub fn clear_forces(&mut self) {
        self.force_accumulator = Spatial4D::zero();
    }

    pub fn integrate(&mut self, dt: f64) {
        let inv_mass = 1.0 / self.knot.mass;
        self.acceleration = self.force_accumulator * inv_mass;
        self.velocity = self.velocity + (self.acceleration * dt);
        self.position = self.position + (self.velocity * dt);
    }

    /// UPDATED: Accepts full Interaction Result to derive Physics from Geometry
    pub fn apply_interaction(&mut self, result: &GemInteractionResult, target_pos: Spatial4D) {
        let displacement = target_pos - self.position;
        let dist_val = displacement.magnitude().norm();
        
        // Avoid singularity
        if dist_val < 1.0e-30 { return; }

        let direction = displacement.normalize();

        // --- DERIVING FORCE FROM GEOMETRY ---
        // 1. Get Binding Energy (Joules)
        // We use the Imaginary component because in GEM, the active binding force
        // (Coulomb/Strong) usually lives on the imaginary axis.
        let energy_joules = result.binding_energy.im; 

        // 2. Derive Newtonian Force Magnitude (F = E / r)
        // Energy is Potential (J). Force is Gradient (J/m).
        // This scales the geometric potential down to the physical force experienced by the mass.
        let force_newtonian = energy_joules.abs() / dist_val;

        // 3. Apply Radial Force (Gravity/Coulomb)
        // Direction: Points to target (Attraction)
        let force_radial = direction * force_newtonian;

        // 4. Apply Tangential Force (Spin/Torque)
        // We calculate the "Twist Ratio" from the raw Geometric Tension.
        // If the manifold has high Imaginary Tension vs Real Tension, we apply more torque.
        let twist_ratio = if result.force.re.abs() > 1.0e-30 {
            result.force.im / result.force.re
        } else {
            0.0
        };
        
        // Define Rotation Plane (Z-axis up)
        let up_vector = if direction.x.norm() < 1e-9 && direction.y.norm() < 1e-9 {
             Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0), Complex64::new(0.0,0.0))
        } else {
             Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0))
        };
        
        let tangent = direction.cross_3d(&up_vector).normalize();
        
        // Scale the Torque by the same Newtonian magnitude
        let force_tangential = tangent * (force_newtonian * twist_ratio.abs());

        // Accumulate Total Force
        self.force_accumulator = self.force_accumulator + force_radial + force_tangential;
    }
}