// bin/generate_plot.rs
use std::error::Error;
use std::fs::File;
use std::io::Write;
use gemphy::{GeometricEncodedMedium, GeometricKnot}; // Use your lib types

const M_SUN: f64 = 1.9884e30;

fn main() -> Result<(), Box<dyn Error>> {
    // 1. Initialize your exact Medium
    let med = GeometricEncodedMedium::new();
    
    // 2. Define the Sun using your Knot logic (Shadow Charge = Xi * Mass)
    // We assume the Sun is the central anchor.
    let sun = GeometricKnot::new(med.clone(), M_SUN, med.xi * M_SUN, 0.0, "Sun");

    // 3. Prepare Output File
    let mut file = File::create("gem_plot_data.csv")?;
    writeln!(file, "distance_m,newton_force,gem_force_real,gem_force_imag,kappa_real,kappa_imag")?;

    // 4. Simulation Loop (Logarithmic scale from 1km to 100 AU)
    let start_dist: f64 = 1_000.0; // 1 km
    let end_dist: f64 = 1.5e13;    // ~100 AU
    let steps = 1000;
    
    // Use a test mass (e.g., Earth Mass) to visualize the field strength at distance d
    let m_earth = 5.972e24;
    let earth = GeometricKnot::new(med.clone(), m_earth, med.xi * m_earth, 0.0, "EarthProbe");

    for i in 0..steps {
        // Logarithmic step generation
        let log_start= start_dist.ln();
        let log_end = end_dist.ln();
        let factor = (log_end - log_start) / (steps as f64);
        let d = (log_start + factor * (i as f64)).exp();

        // A. Standard Newton Calculation (for comparison)
        let f_newton = (med.g * sun.mass * earth.mass) / d.powi(2);

        // B. GEM Unified Calculation (Using your Library)
        // This automatically applies Xi, Phase Rotation, and Complex projection
        let res = med.calculate_interaction(sun.mass, earth.mass, d.into());
        
        // Write data: Distance, Newton, GEM(Real), GEM(Imag), Kappa(Real), Kappa(Imag)
        writeln!(file, "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}", 
            d, 
            f_newton, 
            res.force.re, 
            res.force.im, 
            res.curvature.re, 
            res.curvature.im
        )?;
    }
    println!(" {:#?}", earth);
    println!(" {:#?}", sun);
    println!("Data generation complete. Saved to 'gem_plot_data.csv'.");
    Ok(())
}
// use std::fs::File;
// use std::io::Write;
// use gemphy::GeometricEncodedMedium;

// fn main() -> std::io::Result<()> {
//     let medium = GeometricEncodedMedium::new();
    
//     // 1. SETUP: Unified Interaction (Planck Scale -> Macro Scale)
//     // We use Planck Mass particles to define the baseline Schwarzschild radius.
//     let mass = medium.m_p;
//     let r_s = (2.0 * medium.g * (mass + mass)) / medium.c.powi(2);

//     println!("--- GENERATING GEM UNITY DATASET ---");
//     println!("Mass Base: {:.4e} kg", mass);
//     println!("Horizon:   {:.4e} m", r_s);

//     let mut file = File::create("gem_unity_data.csv")?;
//     // CSV Header
//     writeln!(file, "distance_m,ratio_rs_d,re_go,im_go,re_force,im_force,re_kappa")?;

//     // 2. LOGARITHMIC SWEEP
//     // Start deep inside the horizon (0.1 * rs)
//     // End at Macro scale (10^20 * rs)
//     // This covers everything from Singularity -> Nucleus -> Earth -> Galaxy
    
//     let start_dist = 0.5 * r_s; 
//     let end_dist = 1.0e15 * r_s; 
    
//     let mut current_d = start_dist;
    
//     // We use a multiplier to step logarithmically (multiply by 1.1 each step)
//     // This gives us high resolution near the start and covers huge distances quickly.
//     while current_d < end_dist {
//         let result = medium.calculate_interaction(mass, mass, current_d);
        
//         // Log the raw data
//         writeln!(file, "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
//             current_d,
//             r_s / current_d,
//             result.g_o.re,
//             result.g_o.im,
//             result.force.re,
//             result.force.im,
//             result.curvature.re
//         )?;

//         // Step size: Increase distance by 5% each step for a smooth log curve
//         current_d *= 1.05;
//     }

//     println!("Success: 'gem_unity_data.csv' created.");
//     println!("Plot this file using Log-Log axes to see the Asymptotic Fade.");
//     Ok(())
// }