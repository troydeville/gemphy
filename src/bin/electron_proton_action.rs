use gemphy::{knot::GeometricKnot, medium::{GAMMA_P, GeometricEncodedMedium}};
use physical_constants::ELEMENTARY_CHARGE;

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

    
    let rg1 = (GAMMA_P / (electron.mass * medium.alpha)).powi(2);
    let rg2 = (GAMMA_P / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let result = medium.calculate_interaction(&electron, &proton, d.into());

    println!("muon:                {:#?}", electron);
    println!("proton:              {:#?}", proton);
    println!("Result:              {:#?}", result);
    println!("er1 (eV):            {:#?}", result.er1.norm()/ ELEMENTARY_CHARGE);
    println!("ei1 (eV):            {:#?}", result.ei1.norm()/ ELEMENTARY_CHARGE);
    println!("binding_energy (eV): {:#?}", result.binding_energy.norm()/ ELEMENTARY_CHARGE);
    println!("Go :                 {:#?}", result.g_o);   

    Ok(())
}