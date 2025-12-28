use num_complex::{Complex64};
use crate::geometry::Spatial4D;
use crate::constants::{GAMMA, XI_MAG, LIGHT_BOUNDARY, S_RADIUS};

/// The Reaction Boundary.
/// Represents a finite packet of Action within the medium.
#[derive(Debug, Clone)]
pub struct ReactionBoundary {
    pub location: Spatial4D,
    pub energy: f64,
}

impl ReactionBoundary {
    pub fn new(location: Spatial4D, energy: f64) -> Self {
        Self { location, energy }
    }

    /// Frequency derived from Energy ($E = hf$).
    pub fn frequency(&self) -> f64 {
        const H: f64 = 6.626e-34; 
        self.energy / H
    }

    /// **Derived Time**: $T = 1/f$.
    /// Time is an emergent property of Action, not a dimension.
    pub fn derived_time(&self) -> f64 {
        let f = self.frequency();
        if f.abs() < 1e-50 { 0.0 } else { 1.0 / f }
    }

    /// Calculates Geometric Mass using $\Gamma$.
    /// **Enforces Finite Geometry**: If length < S, length is clamped to S.
    pub fn geometric_mass(&self) -> Complex64 {
        let length = self.location.magnitude();
        
        // GEM Core Principle: No Singularities.
        // We cannot have a length smaller than the fundamental radius S.
        let eff_len = if length.norm() < S_RADIUS { 
            Complex64::from(S_RADIUS) 
        } else { 
            length 
        };
        
        Complex64::from(GAMMA) / eff_len
    }

    /// Validates if the Action fits the Charge-Mass Unity ($\Xi$).
    pub fn is_xi_compliant(charge: f64, mass: f64) -> bool {
        let ratio = charge / mass;
        (ratio - XI_MAG).abs() < 1.0
    }

    /// Checks if action is within the volumetric Light boundary ($4\pi c$).
    pub fn is_within_light_boundary(&self) -> bool {
        self.location.magnitude().norm() < LIGHT_BOUNDARY
    }
}