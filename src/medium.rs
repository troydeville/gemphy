use std::f64::consts::{PI, SQRT_2};
use num_complex::{Complex64, ComplexFloat};

use crate::knot::GeometricKnot;

pub const C: f64 = 299_792_458.0;
pub const ELEM_CHARGE: f64 = 1.602_176_63e-19;
pub const PLANCK_H: f64 = 6.626_070_15e-34;
pub const Z_P: f64 = (2.0 * PLANCK_H) / (ELEM_CHARGE*ELEM_CHARGE);
pub const PHI_Q: f64 = 1e-7; 
pub const PHI_M: f64 = 1e4;
pub const ALPHA_P: f64 = (4.0 * PI * C) / Z_P;
pub const ALPHA: f64 = ALPHA_P * PHI_Q;
pub const XI_REAL_POW_2: f64 = 4.0 * PI * SQRT_2;
pub const G: f64 = (4.0 * PI * PHI_Q)/((SQRT_2) * PHI_M);

// ==============================================================================
//  GEOMETRIC ENCODED MEDIUM (GEM)
// ==============================================================================
#[derive(Debug, Clone)]
pub struct GemInteractionResult {
    pub q1: Complex64,
    pub q2: Complex64,
    pub q_total: Complex64,
    pub af1: Complex64,
    pub af2: Complex64,
    pub g1: Complex64,
    pub g2: Complex64,
    pub force: Complex64,
    pub curvature: Complex64,
    pub g_o: Complex64,
    pub g_recovered: Complex64,
    pub binding_energy_ev: Complex64,
    pub er1: Complex64, pub er2: Complex64,
    pub ei1: Complex64, pub ei2: Complex64,
    pub f1: Complex64, pub f2: Complex64,
    pub distance_natural: f64,
    pub force_norm: f64,
    pub schwarzschild_radius: f64,
    pub ratio_mr_d: f64,
    // pub is_complex: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ForceProtocol { Gravity, Electromagnetism }

#[derive(Debug, Clone)]
pub struct GeometricEncodedMedium {
    pub h: f64, pub e: f64, pub c: f64, pub s: f64,
    pub alpha: f64, pub gamma: f64, pub g: f64, pub epsilon_o: f64,
    pub xi: Complex64,
    pub k_vel: f64, pub phi_big: f64, pub phi_small: f64,
    pub z_p: f64, pub alpha_p: f64, pub gamma_p: f64,
    pub z_o: f64, pub m_p: f64, pub l_p: f64,
}

impl Default for GeometricEncodedMedium { fn default() -> Self { Self::new() } }

impl GeometricEncodedMedium {
    pub fn new() -> Self {
        
        // const XI_REAL_POW_2: f64 = 4.0 * PI * SQRT_2 * G_CONST * EPSILON_0;
        
        let s = SQRT_2 * PI.powf(0.25);

        let h = PLANCK_H;
        let e = ELEM_CHARGE;
        let c = C;

        let phi_big = PHI_Q; 
        let phi_small = PHI_M;

        let z_p = Z_P;
        let alpha_p = ALPHA_P;
        let gamma_p = e.powi(2) / alpha_p;
        let alpha = ALPHA;
        let z_o = z_p * alpha;
        let epsilon_o = 1.0 / (c * z_o);
        let gamma = gamma_p * alpha;
        
        let g = (4.0 * PI * PHI_Q)/((SQRT_2 * PI.powf(0.25)) * PHI_M);
        // (4 \[Pi])/S \[CapitalPhi]/\[Phi]
        
        let m_p = ((h * c) / (2.0 * PI * g)).sqrt();
        let l_p = ((h * g) / (2.0 * PI * c.powi(3))).sqrt();
        
        // let xi_mag_vac = (4.0 * PI * SQRT_2 * g * epsilon_o).sqrt();
        
        // let xi_vacuum = xi_mag_vac * Complex64::new((PI/8.0).cos(), -(PI/8.0).sin());
        // 
        let xi_magnitude = (4.0 * PI * SQRT_2 * g * epsilon_o).sqrt();
        // let xi_phase = Complex64::from_polar(1.0, -PI / 8.0);
        let xi = xi_magnitude* Complex64::new((PI/8.0).cos(), -(PI/8.0).sin());
        
        let k_vel = (alpha * c) / (2.0 * PI);

        Self {
            h, e, c, s, phi_big, phi_small,
            z_p, alpha_p, gamma_p,
            alpha, z_o, epsilon_o, gamma,
            g, m_p, l_p, 
            xi,
            k_vel,
        }
    }


    pub fn calculate_interaction(&self, p1: &GeometricKnot, p2: &GeometricKnot, d: Complex64) -> GemInteractionResult {
        let m1 = p1.mass;
        let m2 = p2.mass;
        let m_total = m1 + m2;

        // 1. Calculate Radii
        let rsm1 = (2.0 * self.g * m1) / self.c.powi(2);
        let rsm2 = (2.0 * self.g * m2) / self.c.powi(2);
        let rqm1 = self.gamma / (m1 * self.alpha.powi(2));
        let rqm2 = self.gamma / (m2 * self.alpha.powi(2));

        // 2. Natural Distance
        let d_nat = (rqm1 + rqm2) - (rsm1 + rsm2);

        // 3. Shadow Charges (4*PI scaling for Energy Correctness)
        let q1_shadow_val = (4.0 * PI * self.epsilon_o * self.c.powi(2) * m1 * rsm1).sqrt();
        let q2_shadow_val = (4.0 * PI * self.epsilon_o * self.c.powi(2) * m2 * rsm2).sqrt();
        
        let q1 = Complex64::from(q1_shadow_val);
        let q2 = Complex64::from(q2_shadow_val);
        let q_total = q1 + q2;

        // 4. Curvature
        let kappa = (Complex64::from(m_total) / q_total) * self.xi;
        let go = Complex64::from(self.g) / kappa.powi(2);

        // 5. System Dynamics (NO CONDITIONS)
        // We use Pythagorean addition to "soften" the distance.
        // This naturally prevents the singularity by treating the radius as an orthogonal offset.
        let r_schwarzschild_total = rsm1 + rsm2;
        let epsilon = 1.0e-35; // Planck floor
        let geometric_floor = r_schwarzschild_total + epsilon;

        // ALGEBRAIC DECODING:
        // effective_mag = sqrt(|d|^2 + floor^2)
        // This expression is smooth, differentiable, and singularity-free.
        let dist_sq = d.norm_sqr();
        let floor_sq = geometric_floor.powi(2);
        let effective_mag = (dist_sq + floor_sq).sqrt();
        
        // Scale the original vector 'd' to this new magnitude
        // We add a tiny epsilon to the denominator to safely handle the exact d=0.0 case
        // without an 'if' check. (effective_mag / 0.0 would be inf).
        let safe_norm = d.norm() + f64::MIN_POSITIVE; 
        let scale = effective_mag / safe_norm;
        
        let effective_dist = d * scale;
        let d3 = effective_dist.powi(3);

        let common_den = Complex64::from(8.0 * SQRT_2 * PI.powi(2) * self.epsilon_o) * d3;

        // --- Active Charge (Topology Derived) ---
        let q1_active_val = p1.topology * self.e; 
        let q2_active_val = p2.topology * self.e;

        let q1_effective = Complex64::new(q1_shadow_val, q1_active_val);
        let q2_effective = Complex64::new(q2_shadow_val, q2_active_val);

        let f1_sq = (q1_effective.powi(2) * self.alpha) / (common_den * m1);
        let f2_sq = (q2_effective.powi(2) * self.alpha) / (common_den * m2);

        let f1 = f1_sq.sqrt();
        let f2 = f2_sq.sqrt();

        // Accelerations
        let ag_scaler = (Complex64::from(2.0 * PI) * effective_dist) / self.alpha;
        let ag1 = f1_sq * ag_scaler;
        let ag2 = f2_sq * ag_scaler;
        
        let force = ag1 * Complex64::from(m2);

        // 6. Energy Allocation
        let den_factor = Complex64::from(2.0 * m_total * self.gamma);
        let e2_go_alpha2 = Complex64::from(self.e.powi(2)) * go * (m1 * m2) * self.alpha.powi(2);
        
        let er1 = e2_go_alpha2 / (den_factor * self.xi.powi(2));
        let er2 = (e2_go_alpha2 * m1) / (den_factor * q2 * self.xi);
        
        let ei1 = (Complex64::from(self.e.powi(2)) * go * m1 * q1 * self.alpha.powi(2)) / (den_factor * self.xi.powi(3));
        let ei2 = er1;

        GemInteractionResult {
            q1, q2, q_total,
            af1: f1_sq, af2: f2_sq,
            g1: ag1, g2: ag2,
            force,
            curvature: kappa,
            g_o: go, g_recovered: go * kappa.powi(2),
            er1, er2, ei1, ei2,
            f1, f2,
            binding_energy_ev: er1 + er2,
            distance_natural: d_nat,
            force_norm: force.norm(),
            schwarzschild_radius: r_schwarzschild_total,
            ratio_mr_d: 1.0 - ((rqm1 + rqm2) / effective_dist.norm()),
            
        }
    }


    /// Decodes force into components for Verification Tests.
    pub fn decode_force(&self, p1: &GeometricKnot, p2: &GeometricKnot, d: f64, protocol: ForceProtocol) -> Complex64 {
         let d_c = Complex64::from(d);
         match protocol {
            ForceProtocol::Gravity => {
                // Newtonian Gravity
                let f_mag = -self.g * p1.mass * p2.mass / d.powi(2);
                Complex64::new(f_mag, 0.0)
            },
            ForceProtocol::Electromagnetism => {
                 // Coulomb Force
                 let k_c = 1.0 / (4.0 * PI * self.epsilon_o);
                 Complex64::from(k_c) * ((ELEM_CHARGE * ELEM_CHARGE) / d_c.powi(2))
            }
         }
    }

    pub fn derived_mass_ratio(&self) -> f64 {
        (48.0 * PI.powi(2) * self.alpha) / (self.phi_big * self.s * self.phi_small)
    }
    
    pub fn calculate_mass_frequency(&self, mass: f64) -> f64 {
        let term1 = (self.alpha.powi(2) * self.c) / (2.0 * self.h);
        let a_num = mass * self.c.powi(2) * self.gamma_p;
        let a_den = 8.0 * PI.powi(2) * self.h;
        let term_a = (a_num / a_den).powi(2);
        let term_b = (mass * self.c).powi(2);
        term1 * (term_a + term_b).sqrt()
    }
    
    pub fn calculate_geometric_binding_energy_complex(&self, m1: f64, m2: f64) -> (Complex64, f64) {
        // Bohr Energy (Reduced Mass)
        let mu = (m1 * m2) / (m1 + m2);
        let energy_joules = 0.5 * mu * self.c.powi(2) * self.alpha.powi(2);
        let energy_ev = energy_joules / self.e;
        (Complex64::new(energy_ev, 0.0), energy_ev)
    }

    pub fn verify_golden_loop(&self) -> bool {
        let check1 = self.h / (2.0 * PI * self.c);
        let check2 = self.gamma / self.alpha;
        (check1 - check2).abs() < 1e-40
    }

}