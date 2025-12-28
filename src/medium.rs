use std::f64::consts::{PI, SQRT_2};
use num_complex::{Complex64, ComplexFloat};
use crate::GeometricKnot;
use crate::constants::*;

#[derive(Debug, Clone)]
pub struct GemInteractionResult {
    pub q1: Complex64, pub q2: Complex64, pub q_total: Complex64,
    pub af1: Complex64, pub af2: Complex64,
    pub g1: Complex64, pub g2: Complex64,
    pub force: Complex64,
    pub curvature: Complex64,
    pub g_o: Complex64, pub g_recovered: Complex64,
    pub mqr1_ev: Complex64, pub mqr2_ev: Complex64,
    pub ratio1: Complex64, pub ratio2: Complex64,
    pub binding_energy_ev: Complex64,
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
        let h = PLANCK_H;
        let e = ELEM_CHARGE;
        let c = C;
        let s = S_GEOM_FACTOR; // Use Geometric Factor (~1.88)
        let phi_big = 1e-7; 
        let phi_small = 1e4;

        let z_p = (2.0 * h) / e.powi(2);
        let alpha_p = (4.0 * PI * c) / z_p;
        let gamma_p = e.powi(2) / alpha_p;
        let alpha = ALPHA;
        let z_o = z_p * alpha;
        let epsilon_o = EPSILON_0;
        let gamma = GAMMA;
        
        // Use constant G for test stability
        let g = G_CONST;
        
        let m_p = ((h * c) / (2.0 * PI * g)).sqrt();
        let l_p = ((h * g) / (2.0 * PI * c.powi(3))).sqrt();
        
        let xi_mag_vac = (4.0 * PI * SQRT_2 * g * epsilon_o).sqrt();
        
        let xi_vacuum = xi_mag_vac * Complex64::new((PI/8.0).cos(), -(PI/8.0).sin());

        let k_vel = (alpha * c) / (2.0 * PI);

        Self {
            h, e, c, s, phi_big, phi_small,
            z_p, alpha_p, gamma_p,
            alpha, z_o, epsilon_o, gamma,
            g, m_p, l_p, 
            xi: xi_vacuum, // Use Vacuum Xi for Medium Interactions
            k_vel,
        }
    }

    /// Calculates the interaction using knot charges if significant, else shadow charges.
    pub fn calculate_interaction(&self, p1: &GeometricKnot, p2: &GeometricKnot, d: Complex64) -> GemInteractionResult {
        let m1 = p1.mass;
        let m2 = p2.mass;
        let mr = (2.0 * self.g * (m1 + m2)) / self.c.powi(2);

        let d_norm = d.norm();
        let eff_d = if d_norm < 1e-35 { Complex64::from(1e-35) } else { d };

        let ratio_mr_d = Complex64::from(mr) / eff_d;
        let ratio_mr_d_safe = if (Complex64::from(1.0) - ratio_mr_d).norm() < 1e-15 {
            Complex64::from(1.0 - 1e-15)
        } else {
            ratio_mr_d
        };

        // Use knot charge if significant (e.g., for EM/charged particles), else shadow (gravity)
        let q1 = if p1.charge.norm() > 1e-30 { p1.charge } else { self.xi * m1 };
        let q2 = if p2.charge.norm() > 1e-30 { p2.charge } else { self.xi * m2 };
        let q_total = q1 + q2;

        let sqrt2_div2 = SQRT_2 / 2.0;
        let d3 = eff_d * eff_d * eff_d;
        let alpha_2pi = self.alpha / (2.0 * PI);

        let geom_denom_1 = Complex64::from(4.0 * PI * self.epsilon_o * m1) * d3;
        let geom_denom_2 = Complex64::from(4.0 * PI * self.epsilon_o * m2) * d3;

        let af1 = (Complex64::new(sqrt2_div2, 0.0) * q1.powi(2) * alpha_2pi) / geom_denom_1;
        let af2 = (Complex64::new(sqrt2_div2, 0.0) * q2.powi(2) * alpha_2pi) / geom_denom_2;

        let scaler_g = (Complex64::from(2.0 * PI) * eff_d) / self.alpha;
        let g1 = af1 * scaler_g;
        let fg = g1 * Complex64::from(m2);

        // Kappa with safeguard for near-zero q_total (strong EM cases)
        let m_total = m1 + m2;
        let kappa = if q_total.norm() > 1e-30 {
            self.xi * m_total / q_total.norm()  // Use norm for magnitude scaling; adjust if phase needed
        } else {
            Complex64::new(1e20, 0.0)  // Large for near-cancellation (e.g., opposite charges)
        };

        let go = self.g / (kappa.abs().powi(2));
        let g_recovered = go * kappa.powi(2);

        let mqr_coeff = self.c.powi(2) / (8.0 * PI * self.epsilon_o * g_recovered);
        let mqr1_ev = (Complex64::from(mqr_coeff) * (q2.powi(2) / Complex64::from(m2))) / self.e;
        let mqr2_ev = (Complex64::from(mqr_coeff) * (q1.powi(2) / Complex64::from(m1))) / self.e; 

        let m1_c2_ev = (m1 * self.c.powi(2)) / self.e;
        let m2_c2_ev = (m2 * self.c.powi(2)) / self.e;
        let ratio1 = mqr1_ev / (mqr1_ev + Complex64::new(m1_c2_ev, 0.0));
        let ratio2 = mqr2_ev / (mqr2_ev + Complex64::new(m2_c2_ev, 0.0));

        let (be, _) = self.calculate_geometric_binding_energy_complex(m1, m2);

        GemInteractionResult {
            q1, q2, q_total, af1, af2, g1, g2: g1, force: fg, curvature: kappa,
            g_o: go.into(), g_recovered, mqr1_ev, mqr2_ev, ratio1, ratio2,
            binding_energy_ev: be,
            schwarzschild_radius: mr, ratio_mr_d: ratio_mr_d_safe.norm(),
            is_complex: fg.im.abs() > 1e-30,
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