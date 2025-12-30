use crate::{knot::GeometricKnot, medium::{ELEM_CHARGE, GAMMA_P, GeometricEncodedMedium}};
use num_complex::Complex64;
// We import the external crate directly, not via 'crate::'
use physical_constants::{self, ELEMENTARY_CHARGE}; 

#[test]
fn test_muonic_hydrogen_prediction() {
    // 1. Correct Initialization
    let medium = GeometricEncodedMedium::new();

    // 2. Initialize Particles using physical constants
    let muon = GeometricKnot::new(
        medium.clone(), 
        physical_constants::MUON_MASS, 
        &[-1.0], // Topology (Negative Charge)
        0.0, 
        "Muon"
    );

    let proton = GeometricKnot::new(
        medium.clone(), 
        physical_constants::PROTON_MASS, 
        &[1.0],  // Topology (Positive Charge)
        0.0, 
        "Proton"
    );

    // 3. Calculate Quadrature Distance
    // rg = (Gamma_p / (mass * alpha))^2
    let rg1 = (GAMMA_P / (muon.mass * medium.alpha)).powi(2);
    let rg2 = (GAMMA_P / (proton.mass * medium.alpha)).powi(2);
    
    // d = sqrt(rg1 + rg2)
    let d_val = Complex64::from((rg1 + rg2).sqrt());

    // 4. Run Interaction
    let result = medium.calculate_interaction(&muon, &proton, d_val);

    // 5. Convert Joules to eV
    let binding_ev = result.binding_energy / ELEM_CHARGE;
    let er1_ev = result.er1.norm() / ELEM_CHARGE;
    let ei1_ev = result.ei1.norm() / ELEM_CHARGE;

    println!("--- Muonic Hydrogen Results ---");
    println!("Calculated Distance: {:.6e} m", d_val);
    println!("Total Binding Energy: {:.4} eV", binding_ev);
    println!("Primary Potential (Er1): {:.4} eV", er1_ev);
    println!("Impedance Potential (Ei1): {:.4} eV", ei1_ev);

    // Validation
    assert!(er1_ev > 2528.0 && er1_ev < 2529.0, "Primary Potential should match ~2528 eV");
    assert!(binding_ev.norm()  > 2750.0 && binding_ev.norm() < 2752.0, "Total Binding Energy should match ~2751 eV");
}

use crate::dynamics::DynamicKnot;
use crate::geometry::Spatial4D;


#[test]
fn test_orbit_generation() {
    let medium = GeometricEncodedMedium::new();

    let proton_knot = GeometricKnot::new(
        medium.clone(),
        physical_constants::PROTON_MASS,
        &[1.0], 
        0.0,
        "Proton"
    );

    let electron_knot = GeometricKnot::new(
        medium.clone(),
        physical_constants::ELECTRON_MASS,
        &[-1.0], 
        0.0,
        "Electron"
    );

    // FIX 1: Proton is immutable (Fixed Anchor)
    let proton = DynamicKnot::new(
        proton_knot, 
        Spatial4D::zero()
    );

    let r_bohr = physical_constants::BOHR_RADIUS;
    let mut electron = DynamicKnot::new(
        electron_knot, 
        Spatial4D::new(
            Complex64::new(r_bohr, 0.0), 
            Complex64::ZERO,           
            Complex64::ZERO,           
            Complex64::ZERO           
        )
    );

    println!("--- Orbit Generation Test ---");
    let dt = 1e-19; // Smaller timestep for stability
    
    for i in 0..10 {
        electron.clear_forces();

        let displacement = proton.position - electron.position;
        let dist_val = displacement.magnitude(); 
        
        let interaction = medium.calculate_interaction(&proton.knot, &electron.knot, dist_val);

        // Apply corrected force
        electron.apply_interaction(&interaction, &proton);
        electron.integrate(dt);

        if i % 2 == 0 {
             // Print Force to debug sign
             println!("Step {}: ForceRe={:.2e} PosX={:.4e} VelY={:.4e}", 
                i, 
                interaction.force.re,
                electron.position.x.re, 
                electron.velocity.y.re
            );
        }
    }

    // FIX 2: Check for Inward motion (X should be LESS than start)
    assert!(electron.position.x.re < r_bohr, 
        "Electron moved AWAY! Current: {:.4e}, Start: {:.4e}", 
        electron.position.x.re, r_bohr);

    // Check for Tangential motion
    let tangential_vel = electron.velocity.y.re.abs();
    assert!(tangential_vel > 0.0, "No tangential velocity generated!");
    
    println!("Orbit Start Confirmed. Tangential Vel: {:.4e}", tangential_vel);
}