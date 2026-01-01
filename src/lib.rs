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
//! 

// use gemphy::prelude::*;
use serde::Serialize;
use wasm_bindgen::prelude::*;

// Inside gemphy/src/lib.rs
#[wasm_bindgen]
pub fn wasm_calculate_interaction(m1: f64, top1: f64, m2: f64, top2: f64, dist: f64) -> JsValue {
    let medium = medium::GeometricEncodedMedium::new();
    let knot1 = knot::GeometricKnot::new(medium.clone(), m1, &[top1], 0.0, "p1");
    let knot2 = knot::GeometricKnot::new(medium.clone(), m2, &[top2], 0.0, "p2");
    
    let result = medium.calculate_interaction(&knot1, &knot2, num_complex::Complex64::from(dist));

    let freq = medium.calculate_mass_frequency(m1);

    #[derive(Serialize)]
    struct WasmResponse {
        result: GemInteractionResult,
        frequency: f64,
    }

    serde_wasm_bindgen::to_value(&WasmResponse { result, frequency: freq }).unwrap()
}

pub mod complex_serde {
    use num_complex::Complex64;
    use serde::{self, Serializer, ser::SerializeStruct};

    pub fn serialize<S>(val: &Complex64, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Complex64", 2)?;
        state.serialize_field("re", &val.re)?;
        state.serialize_field("im", &val.im)?;
        state.end()
    }
}


use std::f64::consts::PI;

use physical_constants::PLANCK_CONSTANT;

use crate::medium::{ALPHA, C, GAMMA_P, GemInteractionResult};

pub mod medium;
pub mod knot;
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
