use std::f64::consts::PI;

use num_complex::ComplexFloat;

use crate::{GeometricEncodedMedium, GeometricKnot};


pub struct GemSystem {
    pub medium: GeometricEncodedMedium,
    pub particles: Vec<DynamicParticle>,
    pub time: f64,
}

impl GemSystem {
    pub fn new() -> Self {
        Self {
            medium: GeometricEncodedMedium::new(),
            particles: Vec::new(),
            time: 0.0,
        }
    }

    pub fn add_particle(&mut self, p: DynamicParticle) {
        self.particles.push(p);
    }

    pub fn step(&mut self, dt: f64) {
        let n = self.particles.len();
        let mut forces_real = vec![[0.0; 3]; n];
        let mut phases_imag = vec![0.0; n];

        for i in 0..n {
            for j in (i + 1)..n {
                let p1 = &self.particles[i];
                let p2 = &self.particles[j];

                let dx = p2.position[0] - p1.position[0];
                let dy = p2.position[1] - p1.position[1];
                let dz = p2.position[2] - p1.position[2];
                let dist_sq = dx*dx + dy*dy + dz*dz;
                let dist = dist_sq.sqrt();

                if dist < 1e-35 { continue; } 

                let result = self.medium.calculate_interaction(p1.core.mass, p2.core.mass, dist);

                let f_linear = result.force.re;
                let f_phase = result.force.im;
                let ux = dx / dist;
                let uy = dy / dist;
                let uz = dz / dist;

                forces_real[i][0] += f_linear * ux;
                forces_real[i][1] += f_linear * uy;
                forces_real[i][2] += f_linear * uz;

                forces_real[j][0] -= f_linear * ux;
                forces_real[j][1] -= f_linear * uy;
                forces_real[j][2] -= f_linear * uz;

                phases_imag[i] += f_phase;
                phases_imag[j] += f_phase;
            }
        }

        for i in 0..n {
            let p = &mut self.particles[i];
            let mass = p.core.mass;

            let ax = forces_real[i][0] / mass;
            let ay = forces_real[i][1] / mass;
            let az = forces_real[i][2] / mass;

            p.velocity[0] += ax * dt;
            p.velocity[1] += ay * dt;
            p.velocity[2] += az * dt;

            p.position[0] += p.velocity[0] * dt;
            p.position[1] += p.velocity[1] * dt;
            p.position[2] += p.velocity[2] * dt;

            p.phase_accum += (phases_imag[i] / mass) * dt;
        }

        self.time += dt;
    }

    /// Calculates the full energy state of the system.
    /// Returns a Ledger (Totals) and a detailed list of all pair interactions.
    pub fn account_for_all_energy(&self) -> (EnergyLedger, Vec<InteractionRecord>) {
        let mut ledger = EnergyLedger::default();
        let mut records = Vec::new();

        // A. Kinetic Energy (Sum of 0.5 * m * v^2)
        for p in &self.particles {
            let v_sq = p.velocity[0].powi(2) + p.velocity[1].powi(2) + p.velocity[2].powi(2);
            ledger.kinetic_energy += 0.5 * p.core.mass * v_sq;
        }

        // B. Potential Energy (Sum of all Pairs)
        let n = self.particles.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let p1 = &self.particles[i];
                let p2 = &self.particles[j];

                let dx = p2.position[0] - p1.position[0];
                let dy = p2.position[1] - p1.position[1];
                let dz = p2.position[2] - p1.position[2];
                let dist = (dx*dx + dy*dy + dz*dz).sqrt();

                if dist < 1e-35 { continue; }

                // 1. Standard Physics View (Separate Forces)
                // Grav: -G m1 m2 / r
                let u_grav = -self.medium.g * p1.core.mass * p2.core.mass / dist;
                // Elec: k q1 q2 / r
                let k_c = 1.0 / (4.0 * PI * self.medium.epsilon_o);
                // Note: We use the Real part of charge for standard Coulomb approximation
                let u_elec = k_c * p1.core.charge.re * p2.core.charge.re / dist;

                ledger.potential_gravity += u_grav;
                ledger.potential_electric += u_elec;

                // 2. GEM Unified View (Geometric Binding Energy)
                // We use your verified formula for magnitude: 
                // U_gem = (Go * m1 * m2 * ...) / (dist_factor...) 
                // For dynamic potential at distance r, GEM is effectively:
                // U_gem = - (Force_GEM * dist) integration, or roughly Unified Potential.
                // For exact accounting, we use the interaction result we verified:
                let interact = self.medium.calculate_interaction(p1.core.mass, p2.core.mass, dist);
                
                // Unified Potential U = - Integrate[F dr] ~ -F*r for simple 1/r^2 forces
                // Since GEM Force includes the Phase shift, we take the Magnitude.
                // Note: Attractive forces are negative potential.
                let u_gem = -(interact.force * dist).abs(); 

                ledger.potential_gem_unified += u_gem;

                // 3. Record this specific interaction combination
                records.push(InteractionRecord {
                    pair_name: format!("{} <-> {}", p1.core.name, p2.core.name),
                    distance: dist,
                    force_mag: interact.force.norm(),
                    energy_contribution: u_gem,
                });
            }
        }

        // Total Energy (Hamiltonian)
        ledger.total_energy = ledger.kinetic_energy + ledger.potential_gem_unified;

        (ledger, records)
    }
}


// ==============================================================================
// PART 5: Dynamic Engine
// ==============================================================================

#[derive(Debug, Clone)]
pub struct DynamicParticle {
    pub id: usize,
    pub core: GeometricKnot, 
    pub position: [f64; 3],      
    pub velocity: [f64; 3],      
    pub phase_accum: f64,        
}

impl DynamicParticle {
    pub fn new(id: usize, core: GeometricKnot, pos: [f64; 3], vel: [f64; 3]) -> Self {
        Self {
            id,
            core,
            position: pos,
            velocity: vel,
            phase_accum: 0.0,
        }
    }
}

// 1. New Structs for Energy Accounting
#[derive(Debug, Clone, Default)]
pub struct EnergyLedger {
    pub kinetic_energy: f64,
    pub potential_gravity: f64,
    pub potential_electric: f64,
    pub potential_gem_unified: f64, // The complex magnitude
    pub total_energy: f64,
}

#[derive(Debug, Clone)]
pub struct InteractionRecord {
    pub pair_name: String,
    pub distance: f64,
    pub force_mag: f64,
    pub energy_contribution: f64,
}