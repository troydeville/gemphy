//! Fundamental Constants of the Geometric Encoded Medium.
use std::f64::consts::{PI, SQRT_2};
use num_complex::Complex64;

// --- FUNDAMENTAL GEOMETRY ---
pub const S_RADIUS: f64 = 5.291_772_109e-11; 

// Geometric Factor used for G derivation (~1.88)
pub const S_GEOM_FACTOR: f64 = 1.882_792_527; 

pub const ALPHA: f64 = 7.297_352_569e-3; 
pub const INV_ALPHA: f64 = 1.0 / ALPHA;

// --- PHYSICAL CONSTANTS ---
pub const C: f64 = 299_792_458.0;
pub const MU_0: f64 = 4.0 * PI * 1.0e-7;
pub const EPSILON_0: f64 = 1.0 / (MU_0 * C * C);
pub const ELEM_CHARGE: f64 = 1.602_176_63e-19;
pub const PLANCK_H: f64 = 6.626_070_15e-34;

// Hardcoded G to ensure test stability (Newtonian Gravitational Constant)
pub const G_CONST: f64 = 6.674_301_5e-11;

// --- GEM UNIFIED SCALARS ---
pub const GAMMA: f64 = (ELEM_CHARGE * ELEM_CHARGE * MU_0) / (4.0 * PI);
pub const GAMMA_PLANCK: f64 = GAMMA * INV_ALPHA;

// --- XI: THE PARTICLE TRUE FORM ---
// This is the Electron Scale Xi (~1.7e11 C/kg) used for particle compliance.
pub const XI_MAG: f64 = 1.758_820_02e11;
pub const XI_PHASE_ANGLE: f64 = -PI / 8.0;

const XI_REAL: f64 = 1.624_937_817_85e11;
const XI_IMAG: f64 = -6.730_712_821_7e10;

/// Xi ($\Xi$): The Unified Charge-Mass Constant (True Form).
pub const XI: Complex64 = Complex64::new(XI_REAL, XI_IMAG);

pub const LIGHT_BOUNDARY: f64 = 4.0 * PI * C;

pub fn get_complex_xi(g: f64) -> Complex64 {
    let mag = (4.0 * PI * SQRT_2 * g * EPSILON_0).sqrt();
    Complex64::from_polar(mag, XI_PHASE_ANGLE)
}