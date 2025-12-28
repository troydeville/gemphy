use gemphy::{GeometricKnot, GeometricEncodedMedium};
use num_complex::ComplexFloat;

fn main() -> std::io::Result<()> {

    const M_MUON: f64 = 1.883531627e-28;
    const M_PROTON: f64 = 1.67262192595e-27;

    let medium = GeometricEncodedMedium::new();
    let m1 = M_MUON;
    let m2 = M_PROTON;
    let q1 = medium.xi * m1;
    let q2 = medium.xi * m2;

    let medium = GeometricEncodedMedium::new();

    // let electron = GeometricKnot::new(m1, Complex64::new(-medium.e / SQRT_2, -medium.e / SQRT_2), 0.0, "Electron");
    // let proton = GeometricKnot::new(m2, Complex64::new(medium.e / SQRT_2, medium.e / SQRT_2), 0.0, "Proton");
    let muon = GeometricKnot::new(medium.clone(), m1, q1, 0.0, "Muon");
    let proton = GeometricKnot::new(medium.clone(), m2, q2, 0.0, "Proton");

    
    let rg1 = (medium.gamma_p / (muon.mass * medium.alpha)).powi(2);
    let rg2 = (medium.gamma_p / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let result = medium.calculate_interaction(m1, m2, d.into());

    println!("electron: {:#?}", muon);
    println!("proton: {:#?}", proton);
    println!("Result: {:#?}", result);
    let r1 = result.ratio1.abs();
    println!("Abs[binding_energy_ev] = {:#?}",  r1*(result.binding_energy_ev).abs());
    Ok(())
}
