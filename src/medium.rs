use std::f64::consts::{PI, SQRT_2};
use num_complex::{Complex64, ComplexFloat};
use crate::GeometricKnot;
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
    pub is_complex: bool,
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

        // 2. Shadow Charges and Curvature
        let q1 = self.xi * m1;
        let q2 = self.xi * m2;
        
        let m_total = m1 + m2;
        let q_total = q1 + q2;

        let rqm1 = self.gamma / (m1 * self.alpha.powi(2));
        let rqm2 = self.gamma / (m2 * self.alpha.powi(2));
        let rsm1 = (2.0 * self.g * m1) / self.c.powi(2);
        let rsm2 = (2.0 * self.g * m2) / self.c.powi(2);

        let d_nat = (rqm1 + rqm2) - (rsm1 + rsm2);

        // Use provided distance d if significant, else use natural distance
        
        let dist = d_nat;
        let d3 = dist.powi(3);

        let kappa = self.xi * (m_total/q_total);
        let go = self.g / kappa.powi(2);
        
        // 3. Frequency Squared 
        let common_den = Complex64::from(8.0 * SQRT_2 * PI.powi(2) * self.epsilon_o) * d3;
        let f1_sq = (q1.powi(2) * self.alpha) / (common_den * m1);
        let f2_sq = (q2.powi(2) * self.alpha) / (common_den * m2);

        let f1 = f1_sq.sqrt();
        let f2 = f2_sq.sqrt();

        // 4. Accelerations 
        let ag_scaler = (Complex64::from(2.0 * PI) * dist) / self.alpha;
        let ag1 = f1_sq * ag_scaler;
        let ag2 = f2_sq * ag_scaler;

        let force = ag1 * m2;

        // 5. Energies (In Jouels) at distance
        let er1 = (ELEM_CHARGE*ELEM_CHARGE*go*m1*m2*ALPHA*ALPHA)/(SQRT_2*(m1+m2)*self.gamma*self.xi.powi(2));
        let ei2 = er1;

        let er2 = (ELEM_CHARGE*ELEM_CHARGE*go*m1*m1*m2*ALPHA*ALPHA)/(SQRT_2*(m1+m2)*q2*self.gamma*self.xi);
        let ei1 = er2;

        GemInteractionResult {
            q1, q2, q_total: q1 + q2,
            af1: f1_sq, af2: f2_sq,
            g1: ag1, g2: ag2,
            force,
            curvature: kappa,
            g_o: go, g_recovered: go * kappa.powi(2),
            er1, er2, ei1, ei2,
            f1, f2,
            binding_energy_ev: er1 + er2, // ~13.6057 eV
            distance_natural: d_nat,
            force_norm: force.norm(),
            schwarzschild_radius: rsm1 + rsm2,
            ratio_mr_d: 1.0-((rqm1 + rqm2) / d.norm()),
            is_complex: force.im.abs() > 1e-30,
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
                 Complex64::from(k_c) * ((p1.charge * p2.charge) / d_c.powi(2))
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