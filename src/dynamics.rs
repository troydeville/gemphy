use crate::geometry::Spatial4D;
use crate::knot::GeometricKnot;
use crate::medium::{ALPHA, C, GemInteractionResult}; // Import Result Struct
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
        // self.acceleration = self.acceleration * (1.0/ALPHA.powi(2));
        self.velocity = self.velocity + (self.acceleration * dt);
        self.position = self.position + (self.velocity * dt);
    }

    pub fn distance(&self, knot: &DynamicKnot) -> Spatial4D {
        Spatial4D { 
            x: knot.position.x-self.position.x,
            y: knot.position.y-self.position.y,
            z: knot.position.z-self.position.z,
            w: knot.position.w-self.position.w
        }
    }

    pub fn apply_interaction(&mut self, result: &GemInteractionResult, ref_knot: &DynamicKnot) {
        // let a1 = result.g1;
        let displacement = ref_knot.position - self.position;
        let direction = displacement.normalize();
        let distance = self.distance(ref_knot);
        // self.force_accumulator * result.f1.norm();
        self.acceleration = Spatial4D::new(Complex64::new(result.g1.re, 0.0),Complex64::new(0.0, result.g2.im),Complex64::ZERO,Complex64::ZERO);
        // let a = result.g_o * self.knot.mass / result.distance_natural;
        // self.velocity = self.velocity+ self.acceleration;
        // self.position = self.position + self.velocity;
        result.f1;
        let force_radial = direction * result.f1;
        let force_tangential = direction * result.f1.im;
        // Define Rotation Plane (Z-axis up)
        // let up_vector = if direction.x.norm() <  direction.y.norm()  {
        //      Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0), Complex64::new(0.0,0.0))
        // } else {
        //      Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0))
        // };

        // let tangent = direction.cross_3d(&up_vector).normalize();
        // // result.force.im / result.force.re
        // // let force_tangential = (force_radial.dot(&tangent));
        // let force_tangential = Spatial4D { 
        //     x: force_radial.x*tangent.x, 
        //     y: force_radial.y*tangent.y, 
        //     z: force_radial.z*tangent.z, 
        //     w: force_radial.w*tangent.w, 
        // };
        // self.acceleration =
        self.force_accumulator = force_tangential.normalize();

        self.integrate(1e-28);
        self.clear_forces();
    }

    // pub fn apply_interaction(&mut self, result: &GemInteractionResult, target_pos: Spatial4D) {
    //     let displacement = target_pos - self.position;
    //     let dist_val = displacement.magnitude().norm();
        
    //     // Avoid singularity
    //     if dist_val < 1.0e-30 { return; }

    //     let direction = displacement.normalize();

    //     // --- DERIVING FORCE FROM GEOMETRY ---
    //     // 1. Get Binding Energy (Joules)
    //     // We use the Imaginary component because in GEM, the active binding force
    //     // (Coulomb/Strong) usually lives on the imaginary axis.
    //     let energy_joules = result.binding_energy.im; 

    //     // 2. Derive Newtonian Force Magnitude (F = E / r)
    //     // Energy is Potential (J). Force is Gradient (J/m).
    //     // This scales the geometric potential down to the physical force experienced by the mass.
    //     let force_newtonian = energy_joules.abs() / dist_val;

    //     // 3. Apply Radial Force (Gravity/Coulomb)
    //     // Direction: Points to target (Attraction)
    //     let force_radial = direction * force_newtonian;

    //     // 4. Apply Tangential Force (Spin/Torque)
    //     // We calculate the "Twist Ratio" from the raw Geometric Tension.
    //     // If the manifold has high Imaginary Tension vs Real Tension, we apply more torque.
    //     let twist_ratio = if result.force.re.abs() > 1.0e-30 {
    //         result.force.im / result.force.re
    //     } else {
    //         0.0
    //     };
        
        // // Define Rotation Plane (Z-axis up)
        // let up_vector = if direction.x.norm() < 1e-9 && direction.y.norm() < 1e-9 {
        //      Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0), Complex64::new(0.0,0.0))
        // } else {
        //      Spatial4D::new(Complex64::new(0.0,0.0), Complex64::new(0.0,0.0), Complex64::new(1.0,0.0), Complex64::new(0.0,0.0))
        // };
        
    //     let tangent = direction.cross_3d(&up_vector).normalize();
        
    //     // Scale the Torque by the same Newtonian magnitude
    //     let force_tangential = tangent * (force_newtonian * twist_ratio.abs());

    //     // Accumulate Total Force
    //     self.force_accumulator = self.force_accumulator + force_radial + force_tangential;
    // }
}