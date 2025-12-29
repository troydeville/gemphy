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
pub mod geometry;

// pub use geometry::*;
// pub use num_complex::Complex64;
// pub use medium::{GeometricEncodedMedium, ForceProtocol};
// pub use knot::{GeometricKnot, ImpedanceField};
// pub use system::{GemSystem};


// ==============================================================================
// Tests
// ==============================================================================

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use num_complex::Complex64;

    use crate::{knot::GeometricKnot, medium::{ELEM_CHARGE, GeometricEncodedMedium}};

    use super::*;

    #[test]
    fn verify_mass_ratio_derivation() {
        let med = GeometricEncodedMedium::new();
        let ratio = med.derived_mass_ratio();
        println!("{}", ratio);

        println!("{}", (1.0/(54.0 * PI.powi(2)) * med.alpha) / (med.phi_big * med.s * med.phi_small));
        println!("Derived Proton/Electron Ratio: {:.4}", ratio);
        assert!((ratio - 1836.13).abs() < 0.1);
    }

    #[test]
    fn view_hydrogen_energy_allocation() {

        const M_ELECTRON: f64 = 9.1093837139e-31;
        const M_PROTON: f64 = 1.67262192595e-27;

        let medium = GeometricEncodedMedium::new();
        let m1 = M_ELECTRON;
        let m2 = M_PROTON;

        let q1 = Complex64::new(-medium.e, -medium.e);
        let q2 = Complex64::new(medium.e, medium.e);

        // let electron = GeometricKnot::new(m1, Complex64::new(-medium.e / SQRT_2, -medium.e / SQRT_2), 0.0, "Electron");
        // let proton = GeometricKnot::new(m2, Complex64::new(medium.e / SQRT_2, medium.e / SQRT_2), 0.0, "Proton");
        let electron = GeometricKnot::new(medium.clone(), m1, q1, 0.0, "Electron");
        let proton = GeometricKnot::new(medium.clone(), m2, q2, 0.0, "Proton");

        
        let rg1 = (medium.gamma_p / (electron.mass * medium.alpha)).powi(2);
        let rg2 = (medium.gamma_p / (proton.mass * medium.alpha)).powi(2);
        let d = (rg1+rg2).sqrt();

        let result = medium.calculate_interaction(&electron, &proton, d.into());

        println!("**");
        println!("{:-<20} Energy Allocation {:-<20} ", "", "");
        println!("m1 allocation at d:     {:.5e} eV (Target: 13.5983)", result.er1.norm()/ELEM_CHARGE);
        println!("m2 allocation at d:     {:.5e} eV (Target: 0.0074)", result.er2.norm()/ELEM_CHARGE);
        println!("{:-<60}", "");
        println!("Total energy allocated: {:.5e} eV (Target: 13.6057)", (result.er1 + result.er2).norm()/ELEM_CHARGE);
        println!("{:-<60}", "");
        println!("{:-<60}", "");
        
    }
}
