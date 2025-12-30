use gemphy::{knot::GeometricKnot, medium::{ELEM_CHARGE, GeometricEncodedMedium}};
use num_complex::{ComplexFloat};

fn main() -> std::io::Result<()> {
    let medium = GeometricEncodedMedium::new();

    const M_TAU: f64 = 3.16754e-27;
    const M_PROTON: f64 = 1.67262192595e-27;

    let m1 = M_TAU;
    let m2 = M_PROTON;

    let tau = GeometricKnot::new(medium.clone(), m1, &[-1.0], 0.0, "Tau");
    let proton = GeometricKnot::new(medium.clone(), m2, &[1.0], 0.0, "Proton");

    let rg1 = (medium.gamma_p / (tau.mass * medium.alpha)).powi(2);
    let rg2 = (medium.gamma_p / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let result = medium.calculate_interaction(&tau, &proton, d.into());

    println!("electron: {:#?}", tau);
    println!("proton: {:#?}", proton);
    println!("Result: {:#?}", result);
    println!("energy (eV): {:#?}", result.binding_energy_ev.norm() / ELEM_CHARGE);
    println!("energy (eV): {:#?}", result.er1.norm() / ELEM_CHARGE);
    println!("energy (eV): {:#?}", result.ei1.norm() / ELEM_CHARGE);

    println!("{:.12e}", medium.xi.abs());
    Ok(())
}
