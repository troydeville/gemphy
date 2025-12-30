
use gemphy::{knot::GeometricKnot, medium::{GAMMA_P, GeometricEncodedMedium}};
use physical_constants::ELEMENTARY_CHARGE;

fn main() -> std::io::Result<()> {
    let medium = GeometricEncodedMedium::new();

    const M_MUON: f64 = 1.883531627e-28;
    const M_PROTON: f64 = 1.67262192595e-27;

    let m1 = M_MUON;
    let m2 = M_PROTON;
    // let q1 = Complex64::new(-ELEM_CHARGE, 0.0);
    // let q2 = Complex64::new(ELEM_CHARGE, 0.0);

    let muon = GeometricKnot::new(medium.clone(), m1, &[-1.0], 0.0, "Muon");
    let proton = GeometricKnot::new(medium.clone(), m2, &[1.0], 0.0, "Proton");

    let rg1 = (GAMMA_P/ (muon.mass * medium.alpha)).powi(2);
    let rg2 = (GAMMA_P / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let result = medium.calculate_interaction(&muon, &proton, d.into());

    println!("muon:                {:#?}", muon);
    println!("proton:              {:#?}", proton);
    println!("Result:              {:#?}", result);
    println!("er1 (eV):            {:#?}", result.er1.norm()/ ELEMENTARY_CHARGE);
    println!("ei1 (eV):            {:#?}", result.ei1.norm()/ ELEMENTARY_CHARGE);
    println!("binding_energy (eV): {:#?}", result.binding_energy.norm()/ ELEMENTARY_CHARGE);
    println!("Go :                 {:#?}", result.g_o);
    Ok(())
}
