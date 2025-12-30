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
//! 

use std::f64::consts::PI;

use physical_constants::PLANCK_CONSTANT;

use crate::medium::{ALPHA, C, GAMMA_P};

pub mod medium;
pub mod knot;
pub mod system;
pub mod geometry;
pub mod dynamics;

pub fn calculate_mass_frequency(mass: f64) -> f64 {
    let term1 = (ALPHA * C) / (2.0 * PLANCK_CONSTANT);
    let a_num = mass * C.powi(2) * GAMMA_P;
    let a_den = 8.0 * PI.powi(2) * PLANCK_CONSTANT;
    let term_a = (a_num / a_den).powi(2);
    let term_b = (mass * C).powi(2);
    term1 * (term_a + term_b).sqrt()
}

pub fn s_constant() -> f64 {
    (4.0 * PI).powf(0.25)
}

// ==============================================================================
// Tests
// ==============================================================================

#[cfg(test)] mod tests;
