use gemphy::dynamics::DynamicKnot;
use gemphy::geometry::Spatial4D;
use gemphy::knot::GeometricKnot;
use gemphy::medium::{ALPHA, C, GAMMA, GeometricEncodedMedium};
use num_complex::Complex64;
use num_traits::Zero;
use physical_constants::ELEMENTARY_CHARGE; 
use std::fs::File;
use std::io::Write;
use std::f64::consts::PI;

fn main() -> std::io::Result<()> {
    let medium = GeometricEncodedMedium::new();

    // 1. Setup Particles
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

    // 2. Geometric Time Configuration
    // Frequency of the knot
    let freq = medium.calculate_mass_frequency(electron_knot.mass);
    // let f1 = medium.calculate_mass_frequency(electron_knot.mass);
    // let f2 = medium.calculate_mass_frequency(proton_knot.mass);
    let period = 1.0 / freq;
    
    // User Formula: dt based on 720-degree (4*PI) revolution gaps
    let dt = period / (4.0 * PI.powi(2));

    // 3. Geometric Space Configuration
    // User Formula: Start at the Base Length

    let r_start = electron_knot.base_length; // Reduced radius for stable orbit

    // User Formula: Velocity coupled to Distance and Frequency
    // v = r * omega (rolling condition)
    let v_start = ALPHA * r_start * freq;

    println!("--- Geometric Configuration ---");
    println!("Radius (r):    {:.4e} m", r_start);
    println!("Frequency (f): {:.4e} Hz", freq);
    println!("Velocity (v):  {:.4e} m/s", v_start);
    println!("Time Step (dt): {:.4e} s", dt);

    // 4. Initialize Dynamics
    let proton = DynamicKnot::new(proton_knot, Spatial4D::zero());
    
    let mut electron = DynamicKnot::new(
        electron_knot, 
        Spatial4D::new(
            Complex64::from(r_start),
            Complex64::zero(),           
            Complex64::zero(),           
            Complex64::zero()           
        )
    );

    // Apply Calculated Velocity (Tangential)
    electron.velocity = Spatial4D::new(
        Complex64::zero(),
        Complex64::new(v_start, 0.0), // Y-axis
        Complex64::zero(),
        Complex64::zero()
    );

    // 5. Run Simulation
    let mut file = File::create("orbit_trajectory.csv")?;
    writeln!(file, "step,time,x,y,vx,vy,force_real,force_imag,energy_real,energy_imag")?;

    let steps = 10 * 1024; 

    for i in 0..steps {
        electron.clear_forces();

        let displacement = proton.position - electron.position;
        let dist_val = displacement.magnitude(); 
        
        let interaction = medium.calculate_interaction(&proton.knot, &electron.knot, dist_val);

        // Apply Forces
        electron.apply_interaction(&interaction, proton.position);
        electron.integrate(dt);

        

        if i % 10 == 0 {
            let current_time = dt * i as f64;
            
            // Calculate Effective Force directly from Joules
            let energy_joules = interaction.binding_energy;
            let effective_force = energy_joules/ dist_val;

            writeln!(file, "{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
                i, current_time,
                electron.position.x.re, electron.position.y.re,
                electron.velocity.x.re, electron.velocity.y.re,
                
                effective_force, // Now correctly ~10^-8 N
                interaction.force.im, // Geometric Torque
                
                interaction.binding_energy.re, 
                interaction.binding_energy.im // Joules
            )?;
        }
    }

    println!("Simulation Complete. Data saved to 'orbit_trajectory.csv'");
    Ok(())
}