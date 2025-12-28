use num_complex::{Complex64};
use crate::{GeometricEncodedMedium, ImpedanceField, Spatial4D};

/// The GEM System Solver.
/// Accumulates Action instead of stepping Time.
pub struct GemSystem {
    pub medium: GeometricEncodedMedium,
    pub particles: Vec<ImpedanceField>,
    pub action_accum: f64, 
}

impl Default for GemSystem {
    fn default() -> Self { Self::new() }
}

impl GemSystem {
    pub fn new() -> Self {
        Self {
            medium: GeometricEncodedMedium::new(),
            particles: Vec::new(),
            action_accum: 0.0,
        }
    }

    pub fn step(&mut self, dt: f64) {
        self.resolve_action(dt);
        self.action_accum += dt;
    }

    pub fn add_particle(&mut self, p: ImpedanceField) {
        self.particles.push(p);
    }

    /// Resolves one interval of Action Flow.
    pub fn resolve_action(&mut self, dt: f64) {
        let n = self.particles.len();
        
        // Forces are stored as 4D Complex Vectors
        let mut forces = vec![Spatial4D::zero(); n];

        for i in 0..n {
            for j in (i + 1)..n {
                let p1 = &self.particles[i];
                let p2 = &self.particles[j];

                // 1. Compute 4D Separation
                let dx = p2.position.x - p1.position.x;
                let dy = p2.position.y - p1.position.y;
                let dz = p2.position.z - p1.position.z;
                let dw = p2.position.w - p1.position.w; 
                
                // 2. Complex Distance (Preserves Phase)
                let dist_sq = dx*dx + dy*dy + dz*dz + dw*dw;
                let dist = dist_sq.sqrt();

                // 3. Horn Torus Constraint (No singularity)
                if dist.norm() < 1e-35 { continue; }

                // 4. Calculate Unified Interaction (Gravity + EM via Xi)
                // let interaction = self.medium.calculate_interaction(
                //     p1.core.mass, 
                //     p2.core.mass, 
                //     dist
                // );
                let interaction = self.medium.calculate_interaction(&p1.core, &p2.core, dist);

                // 5. Apply Force to 4D Manifold
                let f = interaction.force; // Complex Force
                let ux = dx / dist;
                let uy = dy / dist;
                let uz = dz / dist;
                let uw = dw / dist;

                let fx = f * ux;
                let fy = f * uy;
                let fz = f * uz;
                let fw = f * uw;

                forces[i].x += fx; forces[i].y += fy; forces[i].z += fz; forces[i].w += fw;
                forces[j].x -= fx; forces[j].y -= fy; forces[j].z -= fz; forces[j].w -= fw;
            }
        }

        // Apply Action to Particles
        for i in 0..n {
            let p = &mut self.particles[i];
            let m = Complex64::from(p.core.mass);
            let dt_c = Complex64::from(dt);

            // Symplectic integration on Complex Manifold
            p.velocity.x += (forces[i].x / m) * dt_c;
            p.velocity.y += (forces[i].y / m) * dt_c;
            p.velocity.z += (forces[i].z / m) * dt_c;
            p.velocity.w += (forces[i].w / m) * dt_c;

            p.position.x += p.velocity.x * dt_c;
            p.position.y += p.velocity.y * dt_c;
            p.position.z += p.velocity.z * dt_c;
            p.position.w += p.velocity.w * dt_c;
        }

        self.action_accum += dt;
    }
}
// Stub structures
#[derive(Debug, Clone, Default)] pub struct EnergyLedger { pub total_energy: f64 }
#[derive(Debug, Clone)] pub struct InteractionRecord {}