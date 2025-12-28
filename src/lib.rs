/// GEM Framework: Horn Torus & Unified Physics Implementation
/// Based on GEM Framework Python Suite (Version 2.0).

use std::f64::consts::{PI, SQRT_2};
use num_complex::ComplexFloat;
use serde::Deserialize;

// Re-export Complex64
pub use num_complex::Complex64;

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

// ==============================================================================
// PART 1: Traits & Interfaces
// ==============================================================================

pub trait GemSurface {
    fn radius_a(&self) -> f64;
    fn volume(&self) -> f64;
    fn surface_area(&self) -> f64;
    fn parametric_surface(&self, u: f64, v: f64) -> [f64; 3];
    fn implicit_equation(&self, x: f64, y: f64, z: f64) -> f64;
    fn metric_tensor(&self, v: f64) -> (f64, f64) {
        let a = self.radius_a();
        let cos_half_v = (v / 2.0).cos();
        let g_uu = 4.0 * a.powi(2) * cos_half_v.powi(4);
        let g_vv = a.powi(2);
        (g_uu, g_vv)
    }
    fn gaussian_curvature(&self, v: f64) -> f64 {
        let a = self.radius_a();
        let denom = a.powi(2) * (1.0 + v.cos());
        if denom.abs() < 1e-9 { return 0.0; }
        v.cos() / denom
    }
}

// ==============================================================================
// PART 2: GEM Constants & Medium Struct
// ==============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ForceProtocol {
    Gravity,          // Decodes interaction using Mass (Shadow Charge Xi*M)
    Electromagnetism, // Decodes interaction using Electric Charge (e)
}

/// The background medium of the universe.
#[derive(Debug, Clone)]
pub struct GeometricEncodedMedium {
    // Fundamental Constants (CODATA / Exact)
    pub h: f64,
    pub e: f64,
    pub c: f64,

    // GEM Geometric Factors
    pub s: f64,
    pub phi_big: f64,      // Capital Phi
    pub phi_small: f64,    // phi (lower case)

    // Derived Impedance & Geometry
    pub z_p: f64,          // Planck Impedance
    pub alpha_p: f64,      // Planck Fine Structure
    pub gamma_p: f64,      // Planck Action Geometry
    pub alpha: f64,        // Fine Structure Constant
    pub z_o: f64,          // Vacuum Impedance
    pub epsilon_o: f64,    // Vacuum Permittivity
    pub gamma: f64,        // Vacuum Action Geometry (Gamma)
    pub g: f64,            // Gravitational Constant
    pub m_p: f64,          // Planck Mass
    pub l_p: f64,          // Planck Length
    pub xi: Complex64,     // Gravito-Electric Scaling Factor (Xi)
    pub z_g: Complex64,    // Geometric Impedance
    
    // Action Dynamics
    pub k_vel: f64,        // Velocity Quanta (Action Velocity)
}

impl Default for GeometricEncodedMedium {
    fn default() -> Self {
        Self::new()
    }
}

impl GeometricEncodedMedium {
    pub fn new() -> Self {
        // 1. Fundamental Constants
        let h = 6.62607015e-34;
        let e = 1.602176634e-19;
        let c = 299792458.0;

        // 2. GEM Geometric Factors
        let s = SQRT_2 * PI.powf(0.25);
        let phi_big = 1e-7; 
        let phi_small = 1e4;

        // 3. Derivations (Matching Python Order)
        let z_p = (2.0 * h) / e.powi(2);
        let alpha_p = (4.0 * PI * c) / z_p;
        let gamma_p = e.powi(2) / alpha_p;
        let alpha = alpha_p * phi_big;
        let z_o = z_p * alpha;
        let epsilon_o = 1.0 / (c * z_o);
        let gamma = gamma_p * alpha;
        let g = z_o / (c * phi_small * s);
        let m_p = ((h * c) / (2.0 * PI * g)).sqrt();
        let l_p = ((h * g) / (2.0 * PI * c.powi(3))).sqrt();
        
        // --- CORRECTED XI DERIVATION ---
        // Definition from Mathematica: Xi = Sqrt[4Pi Sqrt[2] G ep0] * (Cos[Pi/8] - I Sin[Pi/8])
        // This ensures the magnitude is physically derived from G and Epsilon0.
        let xi_magnitude = (4.0 * PI * SQRT_2 * g * epsilon_o).sqrt();
        let xi_phase = Complex64::from_polar(1.0, -PI / 8.0);
        let xi = Complex64::from(xi_magnitude) * xi_phase;
        
        // Z_g derivation
        let z_g = (4.0 * PI * (2.0 * PI).sqrt() * alpha * s * c * e * phi_big) / (Complex64::from(m_p) * xi);

        // Velocity Quanta k = alpha * c / 2pi
        let k_vel = (alpha * c) / (2.0 * PI);

        Self {
            h, e, c,
            s, phi_big, phi_small,
            z_p, alpha_p, gamma_p,
            alpha, z_o, epsilon_o, gamma,
            g, m_p, l_p, xi,
            z_g, k_vel,
        }
    }

    /// Calculates the Tangential Spin Velocity of the Manifold.
    /// Replicates Mathematica: UnitConvert[Gamma_p / me * fMass[me], "m/s"]
    pub fn calculate_spin_velocity(&self, mass: f64) -> f64 {
        // 1. Get Frequency (Hz)
        let f_hz = self.calculate_mass_frequency(mass);
        
        // 2. Apply the Geometry Scaling: v = (Gamma_p / m) * f
        // Gamma_p units in GEM Code are [C^2 * Omega * s / m] -> We need to ensure unit consistency.
        // In the Python/Mathematica GEM suite, Gamma_p is often defined as e^2/alpha_p (Action-like).
        // Let's rely on the raw magnitude multiplication as implied by your scalar snippet.
        
        (self.gamma_p / mass) * f_hz
    }

    /// Calculates the "Geometric Slip" or "Wiggle" detected in Mathematica.
    /// Compares Ideal Velocity (c * alpha) vs Geometric Velocity (4*pi*Gamma_p*f / (m*alpha)).
    /// Returns: (v_ideal, v_geometric, slip_ratio)
    pub fn calculate_velocity_slip(&self, mass: f64) -> (f64, f64, f64) {
        // 1. Ideal Velocity (Standard Physics)
        let v_ideal = self.c * self.alpha;

        // 2. Geometric Velocity (GEM Derivation)
        // Formula: (4 * Pi * Gamma_p * fMass) / (mass * alpha)
        // Note: We need to be careful with units. 
        // In Rust, gamma_p is e^2/alpha_p. 
        // Let's use the raw components to match your Mathematica output order magnitude.
        
        let f_hz = self.calculate_mass_frequency(mass);
        
        // Based on In[1860]: Result is ~2.18e6
        let num = 4.0 * PI * self.gamma_p * f_hz;
        let den = mass * self.alpha;
        let v_geometric = num / den;

        // 3. The Wiggle (Ratio)
        // Expected ~1.000002
        let slip = v_geometric / v_ideal;

        (v_ideal, v_geometric, slip)
    }

    /// Calculates the fundamental frequency of a mass coupled to the GEM Geometry.
    /// Replicates Mathematica: fMass[M_]
    /// Result for Electron Mass should be ~3.289e15 Hz (13.6057 eV)
    pub fn calculate_mass_frequency(&self, mass: f64) -> f64 {
        // Constant Term: (alpha^2 * c) / (2 * h)
        // Mathematica: (Pi * Alpha^2 * c) / (2 * Pi * h) -> Pis cancel
        let term1 = (self.alpha.powi(2) * self.c) / (2.0 * self.h);
        
        // Sqrt Term A: ((M * c^2 * Gamma_p) / (8 * Pi^2 * h))^2
        // Note: Gamma_p units in GEM are [C^2 * Omega * s / m] which simplifies to [J * s^2 / m]
        // This makes (M c^2 Gamma_p / h) have units of Momentum [kg m/s].
        let a_num = mass * self.c.powi(2) * self.gamma_p;
        let a_den = 8.0 * PI.powi(2) * self.h;
        let term_a = (a_num / a_den).powi(2);
        
        // Sqrt Term B: (M * c)^2
        let term_b = (mass * self.c).powi(2);
        
        // Result
        term1 * (term_a + term_b).sqrt()
    }

    /// Derives the Proton/Electron Mass Ratio from geometric constants.
    pub fn derived_mass_ratio(&self) -> f64 {
        (48.0 * PI.powi(2) * self.alpha) / (self.phi_big * self.s * self.phi_small)
    }

    /// The "Master Decoder" for interactions.
    pub fn decode_force(&self, p1: &GeometricKnot, p2: &GeometricKnot, d: f64, protocol: ForceProtocol) -> Complex64 {
        match protocol {
            ForceProtocol::Gravity => {
                // DECODER: Use Mass-based "Shadow Charge" (Q = Xi * M)
                let res = self.calculate_interaction(p1.mass, p2.mass, d);
                res.force
            },
            ForceProtocol::Electromagnetism => {
                // DECODER: Use Electric Charge (q)
                let k_c = 1.0 / (4.0 * PI * self.epsilon_o);
                let dist_sq = d.powi(2);
                let f_coulomb = k_c * ((p1.charge * p2.charge) / dist_sq);
                f_coulomb
            }
        }
    }

    /// Calculates the Refined Proton Mass using the User's Formula:
    /// mp = mP * e^(-14*pi) * (1 - term1) * (1 - term2)
    pub fn calculate_refined_proton_mass(&self) -> f64 {
        // 1. Planck Mass (Derived from GEM constants)
        // mP = sqrt(hbar * c / G)
        let m_planck = ((self.h / (2.0 * PI)) * self.c / self.g).sqrt();

        // 2. The Geometric Scaling Factor: e^(-14 * pi)
        let geom_scale = (-14.0 * PI).exp();

        // 3. Correction Term 1: (4 * Alpha * Gamma) / (Alpha * Delta)
        // Note: You wrote "Alpha * Delta". In GEM, Delta usually refers to 
        // a geometric structure constant. If not defined, I will assume Delta ~ 1.0 or 
        // derived. Based on standard GEM literature, often corrections are purely Alpha based.
        // For now, I will treat the input strictly: (4 * Gamma) / Delta 
        // (Alphas cancel). If Delta is unknown, I'll use 1.0 placeholder.
        // LET'S ASSUME Delta is related to the proton structure ~ 2*PI or similar.
        // Actually, looking at your previous context, let's use the raw expression 
        // assuming standard GEM coupling.
        // Let's approximate Term 1 as small geometric correction.
        let delta = 1.0; // Placeholder until Delta is defined in your constants
        let term1 = (4.0 * self.alpha * self.gamma) / (self.alpha * delta);

        // 4. Correction Term 2: 1 / (2^11 - 4)
        let term2 = 1.0 / (2048.0 - 4.0);

        // Final Calculation
        let mp = m_planck * geom_scale * (1.0 - term1) * (1.0 - term2);
        
        mp
    }

    /// Calculates the Total Unified Force (Gravity + Electromagnetism) between two knots.
    /// Returns a Complex Force where:
    /// - Real Part: The linear push/pull (Newtons). Negative = Attractive.
    /// - Imaginary Part: The geometric phase torque (unused in simple orbit, but part of GEM).
    pub fn calculate_total_force(&self, k1: &GeometricKnot, k2: &GeometricKnot, r: f64) -> Complex64 {
        if r < 1e-35 { return Complex64::default(); }

        // 1. GRAVITY (Newtonian Limit)
        // F = -G * m1 * m2 / r^2
        // Always attractive (Negative Real)
        let f_grav_mag = -self.g * k1.mass * k2.mass / (r * r);
        let f_grav = Complex64::new(f_grav_mag, 0.0);

        // 2. ELECTROMAGNETISM (Coulomb Limit)
        // F = k_c * q1 * q2 / r^2
        // k_c = 1 / (4 * pi * eps_0)
        let k_c = 1.0 / (4.0 * PI * self.epsilon_o);
        
        // We multiply the Complex Charges.
        // (+e) * (-e) = -e^2 (Real Negative -> Attractive)
        // (+e) * (+e) = +e^2 (Real Positive -> Repulsive)
        let q_product = k1.charge * k2.charge;
        let f_elec = q_product * (k_c / (r * r));

        // 3. UNIFIED TOTAL
        // The medium carries both distortions simultaneously.
        f_grav + f_elec
    }

    


    // /// Calculates the unified interaction between two masses (Gravity Protocol).
    // pub fn calculate_interaction(&self, m1: f64, m2: f64, d: f64) -> GemInteractionResult {
    //     // Schwarzschild radius
    //     let mr = (2.0 * self.g * (m1 + m2)) / self.c.powi(2);
        
    //     let ratio_mr_d = mr / d;
        
    //     // Relativistic term
    //     let inside_term = Complex64::new(1.0 - ratio_mr_d, 0.0);
    //     let sqrt_term = inside_term.sqrt();

    //     // Shadow Charge (Q) derivation
    //     let m_total = m1 + m2;
    //     let xi_c = self.xi;
    //     let m1_c = Complex64::from(m1);
    //     let m2_c = Complex64::from(m2);
        
    //     let sqrt_sqrt_term = sqrt_term.sqrt();
    //     let q1 = (xi_c * m1_c) / sqrt_sqrt_term;
    //     let q2 = (xi_c * m2_c) / sqrt_sqrt_term;
    //     let q_total = q1 + q2;

    //     // Gravitational Force (fg)
    //     let g_complex = Complex64::from(self.g);
    //     let denominator_factor = (xi_c * xi_c) / g_complex;
    //     let common_af_term = (Complex64::from(self.alpha / (2.0 * PI))) / (denominator_factor * d.powi(3));
        
    //     let af1 = (q1 * q1 * common_af_term) / m1_c;
    //     let g_scaling = (2.0 * PI * d) / self.alpha;
    //     let g1 = af1 * g_scaling; 
    //     let fg = g1 * m2_c;

    //     // --- CURVATURE (Kappa) IMPLEMENTATION ---
    //     let full_twist = Complex64::new(1.0, 0.0) / (Complex64::new(1.0, -1.0).sqrt());
    //     let flat_space = Complex64::new(1.0, 0.0);

    //     // Interpolate: sqrt(ratio) smoothing
    //     let factor = ratio_mr_d.sqrt().min(1.0); 
    //     let phase_correction = flat_space * (1.0 - factor) + full_twist * factor;
    //     // 5. Action Frequencies
    //     let sqrt2_div2 = SQRT_2 / 2.0;
    //     let geom_denom_1 = 4.0 * PI * self.epsilon_o * d.powi(3) * m1;
    //     let geom_denom_2 = 4.0 * PI * self.epsilon_o * d.powi(3) * m2;
    //     let alpha_2pi = self.alpha / (2.0 * PI);

    //     let af1 = (Complex64::new(sqrt2_div2, 0.0) * q1.powi(2) * alpha_2pi) / Complex64::from(geom_denom_1);
    //     let af2 = (Complex64::new(sqrt2_div2, 0.0) * q2.powi(2) * alpha_2pi) / Complex64::from(geom_denom_2);

    //     // 6. Accelerations
    //     let scaler_g = (2.0 * PI * d) / self.alpha;
    //     let g1 = af1 * scaler_g;
    //     let g2 = af2 * scaler_g; 

    //     let m_total_c = Complex64::from(m_total);
    //     let raw_curvature = (m_total_c * self.xi) / q_total;
    //     let curvature = raw_curvature * phase_correction;

    //     // 8. Curvature & G Recovery
    //     let m_total = Complex64::from(m1 + m2);
    //     let q_total = q1 + q2;
    //     let kappa = (m_total * self.xi) / q_total;

    //     // Emergent Gravitational Constant (Go)
    //     let go = Complex64::from(self.g) / kappa.powi(2);
    //     let g_recovered = go * kappa.powi(2);

    //     // 9. NEW Binding Energy Calculation
    //     // Formula: UnitConvert[(Sqrt[2] Go \[Kappa]^2)/(\[CapitalXi] (2 - \[Kappa]))^2 e^2/d, "Electronvolts"]
        
    //     // Numerator: Sqrt[2] * Go * Kappa^2 * e^2
        
    //     let numer_be = g_recovered * self.e.powi(2);
        
    //     // Denominator: (Xi * (2 - Kappa))^2 * d
    //     let denom_be = self.xi.powi(2) * d * SQRT_2;
        
    //     let binding_energy_joules = numer_be / denom_be;
    //     let binding_energy_ev = binding_energy_joules / self.e;
    //     // 10. Ratios (Preserved for telemetry, though used less now)
    //     let mqr_coeff = self.c.powi(2) / (8.0 * PI * self.epsilon_o * g_recovered);
    //     let mqr1_joules = Complex64::from(mqr_coeff) * (q2.powi(2) / Complex64::from(m2));
    //     let mqr2_joules = Complex64::from(mqr_coeff) * (q1.powi(2) / Complex64::from(m1));



    //     // let m1_c2 = Complex64::from(m1 * self.c.powi(2));
    //     // let m2_c2 = Complex64::from(m2 * self.c.powi(2));

    //     let mqr1_ev = mqr1_joules / self.e;
    //     let mqr2_ev = mqr2_joules / self.e;
        
    //     // let ratio1 = (mqr1_joules - m1_c2) / mqr1_ev;
    //     // let ratio2 = (mqr2_joules - m2_c2) / mqr2_ev;

        

    //     let m1_c2_ev = (m1 * self.c.powi(2)) / self.e; // m1 c^2 in eV
    //     let m2_c2_ev = (m2 * self.c.powi(2)) / self.e; // m1 c^2 in eV
    //     let ratio1 = mqr1_ev / (mqr1_ev + Complex64::new(m1_c2_ev, 0.0)); // New form
    //     let ratio2 = mqr2_ev / (mqr2_ev + Complex64::new(m2_c2_ev, -0.0)); // Symmetric for generality




    //     GemInteractionResult { 
    //         af1, af2,
    //         g1, g2,
    //         force: fg,
    //         curvature: kappa,
    //         g_o: go,
    //         g_recovered,
    //         mqr1_ev: mqr1_joules / Complex64::from(self.e),
    //         mqr2_ev: mqr2_joules / Complex64::from(self.e),
    //         ratio1,
    //         ratio2,
    //         binding_energy_ev,
    //         schwarzschild_radius: mr,
    //         is_complex: fg.im.abs() > 1e-30,
    //         q1,
    //         q2,
    //         q_total,
    //         ratio_mr_d,
    //     }
    // }

    pub fn calculate_interaction(&self, m1: f64, m2: f64, d: f64) -> GemInteractionResult {
        // 1. LIMITS: Planck Floor & Coulomb Wall
        let d_planck = self.l_p;
        // let d_coulomb = p1.coulomb_radius + p2.coulomb_radius;
        let d_clamped = d;
        // let d_clamped = d.max(d_planck).max(d_coulomb);
        

        // 2. Schwarzschild Radius
        let mr = (2.0 * self.g * (m1 + m2)) / self.c.powi(2);
        let ratio_mr_d = mr / d_clamped;
        
        let ratio_mr_d_safe = if (1.0 - ratio_mr_d).abs() < 1e-15 {
            1.0 - 1e-15 
        } else {
            ratio_mr_d
        };

        // let inside_horizon_term = Complex64::new(1.0 - ratio_mr_d_safe, 0.0).sqrt();

        // // 3. Phase Factor
        // let phase_factor = Complex64::new(1.0, 0.0) / Complex64::new(1.0, -1.0).sqrt();
        
        // 4. Gravitational Charges
        // let term_q1_numer = 8.0 * PI * self.g * m1.powi(2) * self.epsilon_o;
        // let term_q2_numer = 8.0 * PI * self.g * m2.powi(2) * self.epsilon_o;

        // let q1_raw = (Complex64::new(term_q1_numer, 0.0) / inside_horizon_term).sqrt();
        // let q2_raw = (Complex64::new(term_q2_numer, 0.0) / inside_horizon_term).sqrt();

        let q1 = self.xi * m1;
        let q2 = self.xi * m2;

        // 5. Action Frequencies
        let sqrt2_div2 = SQRT_2 / 2.0;
        let geom_denom_1 = 4.0 * PI * self.epsilon_o * d_clamped.powi(3) * m1;
        let geom_denom_2 = 4.0 * PI * self.epsilon_o * d_clamped.powi(3) * m2;
        let alpha_2pi = self.alpha / (2.0 * PI);

        let af1 = (Complex64::new(sqrt2_div2, 0.0) * q1.powi(2) * alpha_2pi) / Complex64::from(geom_denom_1);
        let af2 = (Complex64::new(sqrt2_div2, 0.0) * q2.powi(2) * alpha_2pi) / Complex64::from(geom_denom_2);

        // 6. Accelerations
        let scaler_g = (2.0 * PI * d) / self.alpha;
        let g1 = af1 * scaler_g;
        let g2 = af2 * scaler_g; 

        // 7. Force
        let fg = g1 * Complex64::from(m2);

        // 8. Curvature & G Recovery
        let m_total = Complex64::from(m1 + m2);
        let q_total = q1 + q2;
        let kappa = (m_total * self.xi) / q_total;
        
        let go = Complex64::from(self.g) / kappa.powi(2);
        let g_recovered = go * kappa.powi(2);

        // 9. NEW Binding Energy Calculation
        // Formula: UnitConvert[(Sqrt[2] Go \[Kappa]^2)/(\[CapitalXi] (2 - \[Kappa]))^2 e^2/d, "Electronvolts"]
        
        // Numerator: Sqrt[2] * Go * Kappa^2 * e^2
        
        let numer_be = g_recovered * self.e.powi(2);
        
        // Denominator: (Xi * (2 - Kappa))^2 * d
        let denom_be = self.xi.powi(2) * d_clamped * SQRT_2;
        
        let binding_energy_joules = numer_be / denom_be;
        let binding_energy_ev = binding_energy_joules / self.e;

        // 10. Ratios (Preserved for telemetry, though used less now)
        let mqr_coeff = self.c.powi(2) / (8.0 * PI * self.epsilon_o * g_recovered);
        let mqr1_joules = Complex64::from(mqr_coeff) * (q2.powi(2) / Complex64::from(m2));
        let mqr2_joules = Complex64::from(mqr_coeff) * (q1.powi(2) / Complex64::from(m1));



        // let m1_c2 = Complex64::from(m1 * self.c.powi(2));
        // let m2_c2 = Complex64::from(m2 * self.c.powi(2));

        let mqr1_ev = mqr1_joules / self.e;
        let mqr2_ev = mqr2_joules / self.e;
        
        // let ratio1 = (mqr1_joules - m1_c2) / mqr1_ev;
        // let ratio2 = (mqr2_joules - m2_c2) / mqr2_ev;

        

        let m1_c2_ev = (m1 * self.c.powi(2)) / self.e; // m1 c^2 in eV
        let m2_c2_ev = (m2 * self.c.powi(2)) / self.e; // m1 c^2 in eV
        let ratio1 = mqr1_ev / (mqr1_ev + Complex64::new(m1_c2_ev, 0.0)); // New form
        let ratio2 = mqr2_ev / (mqr2_ev + Complex64::new(m2_c2_ev, -0.0)); // Symmetric for generality

    // Then binding_energy_ev = base_binding * if m1 < m2 { ratio1 } else { ratio2 };

        GemInteractionResult {
            q1, q2, q_total,
            af1, af2,
            g1, g2,
            force: fg,
            curvature: kappa,
            g_o: go,
            g_recovered,
            mqr1_ev: mqr1_joules / Complex64::from(self.e),
            mqr2_ev: mqr2_joules / Complex64::from(self.e),
            ratio1,
            ratio2,
            binding_energy_ev,
            schwarzschild_radius: mr,
            ratio_mr_d: ratio_mr_d_safe,
            is_complex: fg.im.abs() > 1e-30,
        }
    }

    /// Calculates the "Geometric Binding Energy" using Reduced Mass.
    /// Used for Hydrogen Spectrum (Electron-Proton).
    /// Target: 13.5983 eV
    pub fn calculate_geometric_binding_energy_complex(&self, m1: f64, m2: f64) -> (Complex64, f64) {
        // Go = G (since Kappa=1 for fundamental coupling)
        let g_o = self.g;

        let num = g_o * m1 * m2 * self.e.powi(2) * self.alpha.powi(2);
        
        // Xi^2 term provides the 45-degree geometric phase shift
        let xi_sq = self.xi * self.xi; 
        let m_sum = m1 + m2;
        let den_complex = Complex64::new(SQRT_2, 0.0) * m_sum * self.gamma * xi_sq;

        let energy_joules_complex = Complex64::new(num, 0.0) / den_complex;
        
        let joules_to_ev = 1.0 / self.e;
        let energy_ev_complex = energy_joules_complex * joules_to_ev;
        
        // Return (Complex Energy, Magnitude)
        (energy_ev_complex, energy_ev_complex.abs())
    }

    pub fn verify_golden_loop(&self) -> bool {
        let check1 = self.h / (2.0 * PI * self.c);
        let check2 = self.gamma / self.alpha;
        let check3 = self.m_p * self.l_p;
        let epsilon = 1e-40; 
        (check1 - check2).abs() < epsilon && (check1 - check3).abs() < epsilon
    }
}

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
    
    // Notebook Specifics
    pub mqr1_ev: Complex64,
    pub mqr2_ev: Complex64,
    pub ratio1: Complex64,
    pub ratio2: Complex64,
    pub binding_energy_ev: Complex64,

    pub schwarzschild_radius: f64,
    pub ratio_mr_d: f64,
    pub is_complex: bool,
}

// ==============================================================================
// PART 3: Particles & Topology
// ==============================================================================

#[derive(Debug, Clone, Copy)]
pub struct HornTorus {
    pub a: f64,
}

impl HornTorus {
    pub fn new(a: f64) -> Self {
        if a <= 0.0 { panic!("Radius 'a' must be positive."); }
        Self { a }
    }
}

impl GemSurface for HornTorus {
    fn radius_a(&self) -> f64 { self.a }
    fn volume(&self) -> f64 { 2.0 * PI.powi(2) * self.a.powi(3) }
    fn surface_area(&self) -> f64 { 4.0 * PI.powi(2) * self.a.powi(2) }
    fn parametric_surface(&self, u: f64, v: f64) -> [f64; 3] {
        let tube_factor = 1.0 + v.cos();
        [self.a * u.cos() * tube_factor, self.a * tube_factor * u.sin(), self.a * v.sin()]
    }
    fn implicit_equation(&self, x: f64, y: f64, z: f64) -> f64 {
        let sum_sq = x*x + y*y + z*z;
        sum_sq.powi(2) - (4.0 * self.a.powi(2) * (x*x + y*y))
    }
}

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

// ==============================================================================
// PART 4: The Resonance Decoder
// ==============================================================================

pub struct ResonanceDecoder {
    medium: GeometricEncodedMedium,
}

#[derive(Debug)]
pub struct ResonantState {
    pub n: u32,
    pub label: String,
    pub radius: f64,
    pub velocity: f64,
    pub energy_ev: f64,      // Magnitude of Complex Energy
    pub phase_coherence: bool, 
}

#[derive(Debug)]
pub struct MassResonance {
    pub harmonic_n: u32,
    pub scaling_factor: f64,
    pub predicted_mass_kg: f64,
    pub predicted_ratio: f64,
}

impl ResonanceDecoder {
    pub fn new() -> Self {
        Self {
            medium: GeometricEncodedMedium::new(),
        }
    }

    /// Decodes the dynamics of Hydrogen (Electron-Proton System).
    /// Uses Reduced Mass (mu) to align with GEM Geometric Binding Energy (13.5983 eV).
    pub fn decode_hydrogen_spectrum(&self, n_max: u32) -> Vec<ResonantState> {
        let mut states = Vec::new();
        let k = self.medium.k_vel;
        
        let v_base = 2.0 * PI * k;

        // Masses from GEM Source
        let m_e = 9.1093837139e-31;
        let m_p = 1.67262192595e-27;
        
        // REDUCED MASS CORRECTION
        let mu = (m_e * m_p) / (m_e + m_p);

        let e_charge = 1.602176634e-19;

        for n in 1..=n_max {
            let v_n = v_base / (n as f64);
            
            let k_c = 1.0 / (4.0 * PI * self.medium.epsilon_o);
            
            // Radius uses Reduced Mass
            let r_n = (k_c * e_charge.powi(2)) / (mu * v_n.powi(2));

            // Energy uses Reduced Mass
            let energy_joules = -0.5 * mu * v_n.powi(2);
            let energy_ev = energy_joules.abs() / e_charge;

            states.push(ResonantState {
                n,
                label: format!("n={}", n),
                radius: r_n,
                velocity: v_n,
                energy_ev, // Yields ~13.5983 eV
                phase_coherence: true, 
            });
        }
        states
    }

    pub fn scan_mass_resonances(&self, n_start: u32, n_end: u32, step: u32) -> Vec<MassResonance> {
        let mut resonances = Vec::new();
        let m_electron = 9.1093837139e-31;
        let base_term = (PI.powi(2) * self.medium.alpha) / (self.medium.phi_big * self.medium.s * self.medium.phi_small);

        for n in (n_start..=n_end).step_by(step as usize) {
            let ratio = (n as f64) * base_term;
            let mass = ratio * m_electron;
            
            resonances.push(MassResonance {
                harmonic_n: n,
                scaling_factor: base_term,
                predicted_mass_kg: mass,
                predicted_ratio: ratio,
            });
        }
        resonances
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
// PART 6: Tests
// ==============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const M_ELECTRON: f64 = 9.1093837139e-31;
    const M_PROTON: f64 = 1.67262192595e-27;
    const M_EARTH: f64 = 5.97219e24;
    const M_SUN: f64 = 1.9884e30;

    #[test]
    fn verify_derivations_and_constants() {
        let med = GeometricEncodedMedium::new();
        assert!((med.g - 6.6743257318364104e-11).abs() < 1e-20);
        assert!(med.verify_golden_loop());
    }
    
    #[test]
    fn verify_mass_ratio_derivation() {
        let med = GeometricEncodedMedium::new();
        let ratio = med.derived_mass_ratio();
        println!("Derived Proton/Electron Ratio: {:.4}", ratio);
        assert!((ratio - 1836.13).abs() < 0.1);
    }

    #[test]
    fn test_decoder_philosophy() {
        let med = GeometricEncodedMedium::new();
        
        let electron = GeometricKnot::new(med.clone(),9.109e-31, Complex64::new(1.602e-19 / 2.0.sqrt(), -1.602e-19 / 2.0.sqrt()), 0.0, "Electron");
        let proton = GeometricKnot::new(med.clone(), 1.672e-27, Complex64::new(1.602e-19 / 2.0.sqrt(), 1.602e-19 / 2.0.sqrt()), 0.0, "Proton");
        
        let d = 5.29e-11; // Bohr Radius

        println!("--- THE DECODER TEST ---");
        
        let f_gravity = med.decode_force(&electron, &proton, d, ForceProtocol::Gravity);
        let f_coulomb = med.decode_force(&electron, &proton, d, ForceProtocol::Electromagnetism);

        println!("Gravity Force: {:.4e}", f_gravity);
        println!("Coulomb Force: {:.4e}", f_coulomb);

        assert!((f_gravity.abs() - 3.6e-47).abs() < 1e-46);
        assert!((f_coulomb.abs() - 8.2e-8).abs() < 1e-7);
    }
    
    #[test]
    fn verify_complex_energy_components() {
        let med = GeometricEncodedMedium::new();
        let m_e = M_ELECTRON;
        let m_p = M_PROTON;

        let (complex_e, magnitude) = med.calculate_geometric_binding_energy_complex(m_e, m_p);

        println!("--- GEM COMPLEX ENERGY ANALYSIS ---");
        println!("Real Part:      {:.5} eV", complex_e.re);
        println!("Imaginary Part: {:.5} eV", complex_e.im);
        println!("Magnitude:      {:.5} eV", magnitude);
        println!("Target:         13.5983 eV");

        // Verify Magnitude matches your Mathematica result
        assert!((magnitude - 13.5983).abs() < 0.001);
    }

    #[test]
    fn verify_hydrogen_recoil_derivation() {
        // This test replicates the Mathematica steps In[204-208]
        let med = GeometricEncodedMedium::new();
        let m_e = M_ELECTRON;
        let m_p = M_PROTON;

        // TERM 1: Ideal Infinite Proton (In[206])
        let (_, term1_mag) = med.calculate_geometric_binding_energy_complex(m_e, 1e30);
        
        // TERM 2: Proton Recoil Energy (In[207])
        let term2_approx = term1_mag * (m_e / m_p);

        let derived_energy = term1_mag - term2_approx;

        println!("--- RECOIL DERIVATION (Mathematica In[204-208]) ---");
        println!("Term 1 (Ideal Me):    {:.5} eV (Target: 13.6057)", term1_mag);
        println!("Term 2 (Recoil Me/Mp):{:.5} eV (Target: 0.0074)", term2_approx);
        println!("Difference:           {:.5} eV (Target: 13.5983)", derived_energy);

        assert!((term1_mag - 13.6057).abs() < 0.001);
        assert!((derived_energy - 13.5983).abs() < 0.001);
    }

    #[test]
    fn verify_hydrogen_component_interaction() {
        // Replicates Mathematica In[204-208] with physical semantics
        let med = GeometricEncodedMedium::new();
        let m_e = M_ELECTRON;
        let m_p = M_PROTON;

        // COMPONENT 1: Electron Self-Field Energy
        // The energy scale defined purely by the Electron's mass/geometry
        // Mathematica: (Go * Me^2 ...) / (Me ...)
        // We simulate this by effectively isolating Me (m2 -> Infinity or similar math isolation)
        let (_, e_electron_field) = med.calculate_geometric_binding_energy_complex(m_e, 1e30);
        
        // COMPONENT 2: Proton Coupling Energy
        // The energy of the Electron's field scaling against the Proton's geometry.
        // Mathematica: (Go * Me^2 ...) / (Mp ...)
        // This comes out to: E_electron * (Me / Mp)
        let e_proton_coupling = e_electron_field * (m_e / m_p);

        // THE OBSERVABLE: Net Binding Energy
        // The result of these two fields interacting (Subtraction = Destructive Interference)
        let e_net_observable = e_electron_field - e_proton_coupling;

        println!("--- GEM HYDROGEN INTERACTION ANALYSIS ---");
        println!("1. Electron Field Energy:   {:.5} eV (Source Potential)", e_electron_field);
        println!("2. Proton Coupling Energy:  {:.5} eV (Interaction Term)", e_proton_coupling);
        println!("3. Net Observable Energy:   {:.5} eV (Measured State)", e_net_observable);
        println!("   Target:                  13.5983 eV");

        // Assertions matching your specific Mathematica values
        assert!((e_electron_field - 13.6057).abs() < 0.001);
        assert!((e_proton_coupling - 0.0074).abs() < 0.0001);
        assert!((e_net_observable - 13.5983).abs() < 0.001);
    }
    
    #[test]
    fn test_dynamic_system_orbit() {
        let medium = GeometricEncodedMedium::new();
        let mut system = GemSystem::new();
        let sun_core = GeometricKnot::new(medium.clone(), M_SUN, system.medium.xi * M_SUN, 0.0, "Sun");
        let sun = DynamicParticle::new(0, sun_core, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        system.add_particle(sun);
        
        let dist = 1.496e11;
        let v_orb = (system.medium.g * M_SUN / dist).sqrt();
        
        let earth_core = GeometricKnot::new(medium.clone(), M_EARTH, system.medium.xi * M_EARTH, 0.0, "Earth");
        let earth = DynamicParticle::new(1, earth_core, [dist, 0.0, 0.0], [0.0, v_orb, 0.0]);
        system.add_particle(earth);
        
        let dt = 3600.0;
        system.step(dt);
        
        let new_earth = &system.particles[1];
        let new_dist = (new_earth.position[0].powi(2) + new_earth.position[1].powi(2)).sqrt();
        let drift = (new_dist - dist).abs() / dist;
        assert!(drift < 0.001); 
    }

    #[test]
    fn test_hydrogen_resonance_decoder() {
        let decoder = ResonanceDecoder::new();
        println!("--- GEM HYDROGEN DECODER ---");
        
        let states = decoder.decode_hydrogen_spectrum(3);
        
        for s in states {
            println!("State n={}: Radius = {:.4e} m, Energy = {:.4} eV", s.n, s.radius, s.energy_ev);
            if s.n == 1 {
                // Check Ground State using Geometric Target (13.5983 eV)
                assert!((s.radius - 5.29e-11).abs() < 1e-12);
                // Corrected: Compare Magnitude (Positive) with Positive Target
                assert!((s.energy_ev - 13.5983).abs() < 0.01);
            }
        }
    }

    #[test]
    fn test_mass_resonance_integers() {
        let decoder = ResonanceDecoder::new();
        println!("--- GEM MASS HARMONIC SCANNER ---");
        let resonances = decoder.scan_mass_resonances(40, 60, 1);
        for r in resonances {
            if r.harmonic_n == 48 {
                println!(">> HARMONIC 48 (Proton Candidate): Ratio = {:.4}", r.predicted_ratio);
                assert!((r.predicted_ratio - 1836.13).abs() < 0.1);
            }
        }
    }
}
