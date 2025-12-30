
use gemphy::{knot::GeometricKnot, medium::{ELEM_CHARGE, GeometricEncodedMedium}};

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

    let rg1 = (medium.gamma_p / (muon.mass * medium.alpha)).powi(2);
    let rg2 = (medium.gamma_p / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let result = medium.calculate_interaction(&muon, &proton, d.into());

    println!("electron:            {:#?}", muon);
    println!("proton:              {:#?}", proton);
    println!("Result:              {:#?}", result);
    println!("er1 (eV):            {:#?}", result.er1.norm() / ELEM_CHARGE);
    println!("ei1 (eV):            {:#?}", result.ei1.norm() / ELEM_CHARGE);
    println!("binding_energy (eV): {:#?}", result.binding_energy.norm() / ELEM_CHARGE);
    println!("Go :                 {:#?}", result.g1.norm());
    Ok(())
}
