use gemphy::{knot::GeometricKnot, medium::{ELEM_CHARGE, GeometricEncodedMedium}};
use num_complex::{Complex64, ComplexFloat};

fn main() -> std::io::Result<()> {

    const M_ELECTRON: f64 = 9.1093837139e-31;
    const M_PROTON: f64 = 1.67262192595e-27;

    let medium = GeometricEncodedMedium::new();
    let m1 = M_ELECTRON;
    let m2 = M_PROTON;

    // let q1 = Complex64::new(-medium.e, -medium.e);
    // let q2 = Complex64::new(medium.e, medium.e);

    let topology1: &[f64] = &[-1.0];
    let topology2: &[f64] = &[1.0];

    // let electron = GeometricKnot::new(m1, Complex64::new(-medium.e / SQRT_2, -medium.e / SQRT_2), 0.0, "Electron");
    // let proton = GeometricKnot::new(m2, Complex64::new(medium.e / SQRT_2, medium.e / SQRT_2), 0.0, "Proton");
    let electron = GeometricKnot::new(medium.clone(), m1, topology1, 0.0, "Electron");
    let proton = GeometricKnot::new(medium.clone(), m2, topology2, 0.0, "Proton");

    
    let rg1 = (medium.gamma_p / (electron.mass * medium.alpha)).powi(2);
    let rg2 = (medium.gamma_p / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let result = medium.calculate_interaction(&electron, &proton, d.into());

    println!("electron: {:#?}", electron);
    println!("proton: {:#?}", proton);

    println!("Result: {:#?}", result);
    println!("(result.q1 / m2).abs(): {:#?}", (result.q1 / m2).abs());
    println!("(result.q2 / m1).abs(): {:#?}", (result.q2 / m1).abs());
    
    println!(" ((result.q2 / m1).abs()/(result.q1 / m2).abs()).sqrt(): {:#?}", ((result.q2 / m2)/(result.q1 / m2)).abs());
    // println!(" ((result.q2 / m1).abs()/(result.q1 / m2).abs()).sqrt(): {:#?}", ((result.q2 / m1).abs()/(result.q1 / m2).abs()).sqrt());
    // println!(" ((result.q2 / m1).abs()/(result.q1 / m2).abs()).sqrt(): {:#?}", (m2/m1).powi(2));
    println!("Abs[binding_energy_ev] = {:#?}", result.er1.norm() / ELEM_CHARGE

);
    Ok(())
}