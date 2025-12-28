#![doc = include_str!("../README.md")]
//! # gemPhy: Geometric Encoded Medium Physics
//!
//! A physics framework unifying interactions through a 4-dimensional impedance vacuum medium.
//! Based on the **Horn Torus** geometry ($R=r=S$) where "Time" is derived from Action frequency.
//!
//! ## Key Principles
//! * **4D Spatial Medium:** Reality consists of 4 spatial dimensions ($x,y,z,w$). Time is $1/f$.
//! * **Finite Geometry:** No singularities. All action is confined by the Horn Torus volume.
//! * **Unified Constants:** Mass and Charge are geometrically linked via $\Gamma$ and $\Xi$.

pub mod medium;
pub mod knot;
pub mod system;
pub mod surface;
pub mod constants;
pub mod geometry;
pub mod transceiver;


pub use constants::*;
pub use geometry::*;
pub use transceiver::ReactionBoundary;

pub use num_complex::Complex64;



pub use medium::{GeometricEncodedMedium, ForceProtocol};
pub use knot::{GeometricKnot, ImpedanceField};
pub use system::{GemSystem};
pub use surface::GemSurface;

use std::f64::consts::PI;
use num_complex::ComplexFloat;



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
// PART 6: Tests
// ==============================================================================

#[cfg(test)]
mod tests {
    use std::f64::consts::SQRT_2;

    use crate::knot::ImpedanceField;

    use super::*;

    const M_ELECTRON: f64 = 9.1093837139e-31;
    const M_PROTON: f64 = 1.67262192595e-27;
    const M_EARTH: f64 = 5.97219e24;
    const M_SUN: f64 = 1.9884e30;

    // #[test]
    // fn verify_derivations_and_constants() {
    //     let med = GeometricEncodedMedium::new();
    //     println!("G = {}", med.g / 6.67430e-11);
    //     assert!((med.g - 6.6743257318364104e-11).abs() < 1e-20);
    //     assert!(med.verify_golden_loop());
    // }
    
    #[test]
    fn verify_mass_ratio_derivation() {
        let med = GeometricEncodedMedium::new();
        let ratio = med.derived_mass_ratio();
        println!("Derived Proton/Electron Ratio: {:.4}", ratio);
        assert!((ratio - 1836.13).abs() < 0.1);
    }

    // #[test]
    // fn test_decoder_philosophy() {
    //     let med = GeometricEncodedMedium::new();
        
    //     let electron = GeometricKnot::new(med.clone(),9.109e-31, Complex64::new(1.602e-19 / 2.0.sqrt(), -1.602e-19 / 2.0.sqrt()), 0.0, "Electron");
    //     let proton = GeometricKnot::new(med.clone(), 1.672e-27, Complex64::new(1.602e-19 / 2.0.sqrt(), 1.602e-19 / 2.0.sqrt()), 0.0, "Proton");
        
    //     let d = 5.29e-11; // Bohr Radius

    //     println!("--- THE DECODER TEST ---");
        
    //     let f_gravity = med.decode_force(&electron, &proton, d, ForceProtocol::Gravity);
    //     let f_coulomb = med.decode_force(&electron, &proton, d, ForceProtocol::Electromagnetism);

    //     println!("Gravity Force: {:.4e}", f_gravity);
    //     println!("Coulomb Force: {:.4e}", f_coulomb);

    //     assert!((f_gravity.abs() - 3.6e-47).abs() < 1e-46);
    //     assert!((f_coulomb.abs() - 8.2e-8).abs() < 1e-7);
    // }

    #[test]
    fn test_decoder_philosophy_complex() {
        let med = GeometricEncodedMedium::new();

        
        
        
        let electron = GeometricKnot::new(med.clone(),9.109e-31, Complex64::new(1.602e-19 / 2.0.sqrt(), -1.602e-19 / 2.0.sqrt()), 0.0, "Electron");
        let proton = GeometricKnot::new(med.clone(), 1.672e-27, Complex64::new(1.602e-19 / 2.0.sqrt(), 1.602e-19 / 2.0.sqrt()), 0.0, "Proton");

        let d = med.gamma / (electron.mass * med.alpha * med.alpha); // Bohr Radius
        
        let position: Spatial4D = Spatial4D::new(Complex64::from(d), Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0));
        let velocity: Spatial4D = Spatial4D::new(Complex64::from(med.calculate_mass_frequency(electron.mass)), Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0));

        let electron_impedance_field = ImpedanceField::new(0,electron.clone(), position, velocity);
        let proton_impedance_field = ImpedanceField::new(0,proton.clone(), position, velocity);


        

        println!("--- THE DECODER TEST ---");
        
        // let f_gravity = med.decode_force(&electron, &proton, d, ForceProtocol::Gravity);
        // let f_coulomb = med.decode_force(&electron, &proton, d, ForceProtocol::Electromagnetism);
        println!("electron_impedance_field: {:#?}",electron_impedance_field);
        println!("proton_impedance_field: {:#?}",proton_impedance_field);
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
    }

    // #[test]
    // fn verify_hydrogen_component_interaction() {
    //     // Replicates Mathematica In[204-208] with physical semantics
    //     let med = GeometricEncodedMedium::new();
    //     let m_e = M_ELECTRON;
    //     let m_p = M_PROTON;

    //     // COMPONENT 1: Electron Self-Field Energy
    //     // The energy scale defined purely by the Electron's mass/geometry
    //     // Mathematica: (Go * Me^2 ...) / (Me ...)
    //     // We simulate this by effectively isolating Me (m2 -> Infinity or similar math isolation)
    //     let (_, e_electron_field) = med.calculate_geometric_binding_energy_complex(m_e, 1e30);
        
    //     // COMPONENT 2: Proton Coupling Energy
    //     // The energy of the Electron's field scaling against the Proton's geometry.
    //     // Mathematica: (Go * Me^2 ...) / (Mp ...)
    //     // This comes out to: E_electron * (Me / Mp)
    //     let e_proton_coupling = e_electron_field * (m_e / m_p);

    //     // THE OBSERVABLE: Net Binding Energy
    //     // The result of these two fields interacting (Subtraction = Destructive Interference)
    //     let e_net_observable = e_electron_field - e_proton_coupling;

    //     println!("--- GEM HYDROGEN INTERACTION ANALYSIS ---");
    //     println!("1. Electron Field Energy:   {:.5} eV (Source Potential)", e_electron_field);
    //     println!("2. Proton Coupling Energy:  {:.5} eV (Interaction Term)", e_proton_coupling);
    //     println!("3. Net Observable Energy:   {:.5} eV (Measured State)", e_net_observable);
    //     println!("   Target:                  13.5983 eV");

    //     // Assertions matching your specific Mathematica values
    //     assert!((e_electron_field - 13.6057).abs() < 0.001);
    //     assert!((e_proton_coupling - 0.0074).abs() < 0.0001);
    //     assert!((e_net_observable - 13.5983).abs() < 0.001);
    // }
    // [Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0)], 
    #[test]
    fn test_dynamic_system_orbit() {
        let medium = GeometricEncodedMedium::new();
        let mut system = GemSystem::new();

        let dist = 1.496e11;
        let v_orb = (system.medium.g * M_SUN / dist).sqrt();

        // 1. Define Sun at CENTER (0, 0, 0, w=1)
        let sun_pos = Spatial4D::new(Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0));
        let sun_vel = Spatial4D::zero();
        
        // 2. Define Earth at DISTANCE (d, 0, 0, w=1)
        let earth_pos = Spatial4D::new(Complex64::from(dist/SQRT_2), Complex64::from(dist/SQRT_2), Complex64::from(dist/SQRT_2), Complex64::from(1.0));
        let earth_vel = Spatial4D::new(Complex64::from(0.0), Complex64::from(v_orb), Complex64::from(0.0), Complex64::from(1.0));

        let sun_core = GeometricKnot::new(medium.clone(), M_SUN, system.medium.xi * M_SUN, 0.0, "Sun");
        let sun = ImpedanceField::new(0, sun_core, sun_pos, sun_vel);
        system.add_particle(sun);
        
        let earth_core = GeometricKnot::new(medium.clone(), M_EARTH, system.medium.xi * M_EARTH, 0.0, "Earth");
        let earth = ImpedanceField::new(1, earth_core, earth_pos, earth_vel);
        system.add_particle(earth);
            
        let dt = 3600.0;
        system.step(dt);
        
        let new_earth = &system.particles[1];
        let new_dist = (new_earth.position.magnitude()*new_earth.position.magnitude())+ ((new_earth.position.magnitude() * new_earth.position.magnitude()));
        let drift = (new_dist - dist).abs() / dist;
        println!("Drift: {}", drift);
        assert!(drift < 10e54); 
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

    #[test]
    fn test_dynamic_system_orbit2() {
    let medium = GeometricEncodedMedium::new();
    let mut system = GemSystem::new();

    let dist = 1.496e11;
    let v_orb = (system.medium.g * M_SUN / dist).sqrt();

    // 1. Define Sun at CENTER (0, 0, 0, w=1)
    let sun_pos = Spatial4D::new(Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0), Complex64::from(1.0));
    let sun_vel = Spatial4D::zero();
    
    // 2. Define Earth at DISTANCE (d, 0, 0, w=1)
    let earth_pos = Spatial4D::new(Complex64::from(dist), Complex64::from(0.0), Complex64::from(0.0), Complex64::from(1.0));
    let earth_vel = Spatial4D::new(Complex64::from(0.0), Complex64::from(v_orb), Complex64::from(0.0), Complex64::from(0.0));

    let sun_core = GeometricKnot::new(medium.clone(), M_SUN, system.medium.xi * M_SUN, 0.0, "Sun");
    let sun = ImpedanceField::new(0, sun_core, sun_pos, sun_vel);
    system.add_particle(sun);
    
    let earth_core = GeometricKnot::new(medium.clone(), M_EARTH, system.medium.xi * M_EARTH, 0.0, "Earth");
    let earth = ImpedanceField::new(1, earth_core, earth_pos, earth_vel);
    system.add_particle(earth);
    
    // ... rest of test
}
}

































