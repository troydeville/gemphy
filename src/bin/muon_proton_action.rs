
use gemphy::{knot::GeometricKnot, medium::{ELEM_CHARGE, GeometricEncodedMedium}};
use num_complex::{ComplexFloat};

fn main() -> std::io::Result<()> {
    let medium = GeometricEncodedMedium::new();

    const M_MUON: f64 = 1.883531627e-28;
    const M_PROTON: f64 = 1.67262192595e-27;

    let m1 = M_MUON;
    let m2 = M_PROTON;
    let q1 = medium.xi * m1;
    let q2 = medium.xi * m2;

    let muon = GeometricKnot::new(medium.clone(), m1, q1, 0.0, "Muon");
    let proton = GeometricKnot::new(medium.clone(), m2, q2, 0.0, "Proton");

    let rg1 = (medium.gamma_p / (muon.mass * medium.alpha)).powi(2);
    let rg2 = (medium.gamma_p / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let result = medium.calculate_interaction(&muon, &proton, d.into());

    println!("electron: {:#?}", muon);
    println!("proton: {:#?}", proton);
    println!("Result: {:#?}", result);
    println!("energy (eV): {:#?}", result.binding_energy_ev.norm() / ELEM_CHARGE);
    println!("energy (eV): {:#?}", result.er1.norm() / ELEM_CHARGE);
    println!("energy (eV): {:#?}", result.ei1.norm() / ELEM_CHARGE);

    println!("{:.12e}", medium.xi.abs());
    Ok(())
}
