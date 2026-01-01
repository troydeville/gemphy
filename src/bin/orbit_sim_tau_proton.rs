use gemphy::dynamics::DynamicKnot;
use gemphy::geometry::Spatial4D;
use gemphy::knot::GeometricKnot;
use gemphy::medium::{ALPHA, GAMMA, GeometricEncodedMedium};
use num_complex::{Complex64};
use num_traits::Zero;
use physical_constants::ELEMENTARY_CHARGE; 
use std::fs::File;
use std::io::Write;

fn main() -> std::io::Result<()> {
    let medium = GeometricEncodedMedium::new();

    let m1 = physical_constants::TAU_MASS;
    let m2 = physical_constants::PROTON_MASS;

    // 1. Setup Particles
    let proton_knot = GeometricKnot::new(
        medium.clone(),
        m2,
        &[1.0], 
        GAMMA / (m2 * ALPHA.powi(2)),
        "Proton"
    );

    let tau_knot = GeometricKnot::new(
        medium.clone(),
        m1, 
        &[-1.0], 
        GAMMA / (m1 * ALPHA.powi(2)),
        "Tau"
    );

    // 2. Geometric Configuration
    let freq = medium.calculate_mass_frequency(tau_knot.mass);
    let r_start = tau_knot.binding_radius; 

    // --- DYNAMIC STABILITY FIX ---
    let initial_dist = Spatial4D::new(
        Complex64::from(r_start), Complex64::zero(), Complex64::zero(), Complex64::zero()
    );
    
    // FIX: Swap order to (Electron, Proton) so m1=Electron. 
    // This aligns with the Mathematica derivation for Energy formulas.
    // Force magnitude remains the same (Newton's 3rd Law).
    let interaction_init = medium.calculate_interaction(&tau_knot, &proton_knot, initial_dist.magnitude());
    
    let force_magnitude = interaction_init.force.norm();
    
    // Velocity v = sqrt(F * r / m_electron)
    let v_stable = (force_magnitude * r_start / m1).sqrt();

    // Use the Period from Orbit Sim (Rydberg scale ~1e-16s)
    let period =  1.0 / freq;
    // Divide period into steps for smooth integration (e.g. 1000 steps per orbit)
    let dt = period / 1000.0; 
    
    println!("--- Initial Geometric Configuration ---");
    println!("Radius (r):    {:.4e} m", r_start);
    println!("Force (F):     {:.4e} N", force_magnitude);
    println!("Velocity (v):  {:.4e} m/s (Stable Orbit)", v_stable);
    println!("Time Step (dt): {:.4e} s", dt);

    let mut proton = DynamicKnot::new(
        proton_knot.clone(), Spatial4D::zero()
    );
    
    let mut electron = DynamicKnot::new(
        tau_knot.clone(), 
        initial_dist
    );

    electron.velocity = Spatial4D::new(
        Complex64::zero(),
        Complex64::from(v_stable), 
        Complex64::zero(),
        Complex64::zero()
    );

    let mut file = File::create("orbit_trajectory_tau_proton.csv")?;
    writeln!(file, "step,time,x,y,vx,vy,F,Er1,Ei1,Et")?;

    let steps = 50000; 
    let log_step = 100;

    for i in 0..steps {
        let dist_val = electron.distance(&proton);
        
        // FIX: Swap order here as well to maintain consistency
        let interaction = medium.calculate_interaction(&electron.knot, &proton.knot, dist_val.magnitude());
        
        if i % log_step == 0 {
            let time = i as f64 * dt;
            
            writeln!(file, "{},{:.4e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}, {:.6e}, {:.6e}",
                i,
                time,
                electron.position.x.norm(), electron.position.y.norm(),
                electron.velocity.x.norm(), electron.velocity.y.norm(),
                // (interaction.f1).norm(), 
                (interaction.force.norm()), 
                interaction.er1.norm() / ELEMENTARY_CHARGE,
                interaction.ei1.norm() / ELEMENTARY_CHARGE,
                (interaction.er1 + interaction.ei1).norm() / ELEMENTARY_CHARGE
            )?;
        }       
    
        // 1. Accumulate Forces
        electron.apply_interaction(&interaction, &proton);
        
        // 2. Integrate with Main Loop Time Step
        electron.integrate(dt);
        
        // 3. Clear Forces
        electron.clear_forces();
        proton.clear_forces();
    }

    println!("Simulation Complete. Data saved to 'orbit_trajectory.csv'");
    Ok(())
}