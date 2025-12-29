
use gemphy::{knot::GeometricKnot, medium::{GeometricEncodedMedium}};
use num_complex::ComplexFloat;

fn main() -> std::io::Result<()> {

    const M_TEST: f64 = 1.0;
    const M_NEUTRON_STAR: f64 = 2.78376e30;
    

    let medium = GeometricEncodedMedium::new();
    let m1 = M_TEST;
    let m2 = M_NEUTRON_STAR;

    let test_mass = GeometricKnot::new(medium.clone(), m1, medium.xi * m1, 0.0, "Test Mass");
    let neutron_star = GeometricKnot::new(medium.clone(), m2, medium.xi * m2, 0.0, "Neutron Star");

    
    let d: f64 = 1e4;

    let result = medium.calculate_interaction(&test_mass, &neutron_star, d.into());

    println!("{:#?}", test_mass);
    println!("{:#?}", neutron_star);

    println!("Result: {:#?}", result);
    
println!("Result: {:12e}", result.g1);
println!("Result: {:12e}", result.g2);
    Ok(())
}