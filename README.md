# Geometric Encoded Medium (GEM) Framework

##

### A Geometric Encoded Medium (GEM) Impedance Framework for Physics in Rust

> "The universe is a perfect geometric circuit."

![Version](https://img.shields.io/badge/version-0.1.0-blue) ![License](https://img.shields.io/badge/license-GPL--3.0-green) ![Rust](https://img.shields.io/badge/built_with-Rust-orange) [![Crates.io](https://img.shields.io/crates/v/gemphy.svg)](https://crates.io/crates/gemphy) [![Docs](https://docs.rs/gemphy/badge.svg)](https://docs.rs/gemphy) [![CI](https://github.com/troydeville/gemphy/workflows/CI/badge.svg)](https://github.com/troydeville/gemphy/actions)

---

## 🌌 Overview

The **GEM Framework** is a Rust library that models reality not as a collection of arbitrary forces and constants, but as a single **Geometric Encoded Medium**. It posits that space itself has impedance ($Z_p$) and geometry (**Horn Torus topology**), and that all physical phenomena—Gravity, Electromagnetism, Mass, and Charge—are simply different "encodings" running on this geometric hardware.

This library does not approximate physics; it **derives** it. By defining a few geometric axioms, it mathematically derives Newton's Gravitational Constant ($G$), the Fine Structure Constant ($\alpha$), and the Proton Radius from first principles. Validations match CODATA with high precision (e.g., α rel error -5.47e-10, G rel error 3.86e-6).

For a visual of the horn torus topology:

## 🧠 Core Philosophy

In GEM, the universe is treated as a software system:

- **The Hardware:** A "Horn Torus" topology representing the vacuum medium.
- **The Operating System:** The Impedance Field ($Z_0, Z_p$) that dictates how information moves.
- **The Software:**
  - **Gravity:** The interaction when the medium is encoded with Mass (Shadow Charge $Q = \Xi M$).
  - **Electromagnetism:** The interaction when the medium is encoded with Electric Charge ($q$).

Complex phases handle extreme scales (e.g., Planck or black holes), rotating real forces into imaginary (rotational/spin) components.

## Quick Start

Clone and build:

```bash
git clone https://github.com/troydeville/gemphy.git
cd gemphy
cargo test  # Run validations
cargo doc --open  # View API docs
```

Install via Cargo:

```toml
[dependencies]
gemphy = "0.1.0"
```

Basic usage:

```rust
use gemphy::GeometricEncodedMedium;

fn main() {
    let medium = GeometricEncodedMedium::new();
    println!("Derived α: {:.4e}", medium.alpha);
    println!("Derived G: {:.4e} m³ kg⁻¹ s⁻²", medium.g);
}
```

## Roadmap

| Milestone | Tasks
|-----------|-------
| v0.3      | Add horn torus simulation; publish crate
| v0.4      | Interactive CLI; Python bindings
| v1.0      | Full predictions (e.g., dark matter as impedance)

## 📦 Installation

Add to `Cargo.toml`:

```toml
[dependencies]
gemphy = "0.1.0"
```

Note: Re-exports `Complex64` from `num-complex`.

## 🚀 Usage Guide

### 1. Initialize the Medium

```rust
use gemphy::GeometricEncodedMedium;

fn main() {
    let medium = GeometricEncodedMedium::new();
    println!("Vacuum Impedance (Z0): {:.4} Ohms", medium.z_o);
    println!("Fine Structure (α):    {:.4e}", medium.alpha);
    println!("Gravitational G:       {:.4e} m³ kg⁻¹ s⁻²", medium.g);
}
```

### 2. Create Geometric Knots

```rust
use gemphy::{GeometricEncodedMedium, GeometricKnot};
use num_complex::Complex64;

fn main() {
    let medium = GeometricEncodedMedium::new();
    let electron = GeometricKnot::new(medium.clone(), 9.109e-31, Complex64::new(-1.602e-19, 0.0), 0.0, "Electron");
    let earth = GeometricKnot::new(medium, 5.972e24, Complex64::new(0.0, 0.0), 6.371e6, "Earth");

    println!("--- Particle Topology ---");
}
```

Output Note: Electron sizes to Bohr Radius; Earth to Schwarzschild Radius.

### 3. Constants & Model Specification

#### I. Fundamental Physical Constants

These are the fixed, integer-defined values representing the bedrock of the SI system used in the model.

- **Planck Constant ($h$):** The quantum of action.
$$h = \frac{662607015}{10^{42}} \text{ J/Hz}$$

- **Speed of Light ($c$):** The universal speed limit.
$$c = 299792458 \text{ m/s}$$

- **Elementary Charge ($e$):** The fundamental unit of charge.
$$e = \frac{1602176634}{10^{28}} \text{ C}$$

#### II. Scaling & Geometric Factors

These parameters act as scaling bridges between the quantum/electromagnetic scale and the gravitational scale.

- **Magnetic Scaling Constant ($\Phi$):** A scaling factor with dimensions of inductance per meter.
$$\Phi = \frac{1}{10^7} \text{ H/m}$$

- **Mass-Charge Metric ($\phi$):** A density-like scalar connecting mass and charge squared.
$$\phi = 10^4 \text{ kg}^2 \text{ m}^{-2} \text{ s}^2 \text{ C}^{-2}$$

- **Alpha-Lambda Function ($\Lambda_\alpha(n)$):** A rational function governing specific integer-ratio field interactions.
$$\Lambda_\alpha(n) = \frac{\alpha_\delta n}{\alpha_\gamma + 2 \alpha_\delta n^2}$$

#### III. Primary Field Parameters (Planck Scale)

These variables (-subscript) represent the theoretical "maximum" impedance and field density before fine-structure scaling is applied.

- **Primary Impedance ($Z_p$):**
$$Z_p = \frac{2h}{e^2}$$

- **Primary Fine Structure ($\alpha_p $):** The geometric inverse of impedance.
$$\alpha_p = \frac{4\pi c}{Z_p}$$

- **Primary Field Constants:** Permeability ($\mu_p$), Permittivity ($\epsilon_p$), and the Gamma factor ($\Gamma_p$).
$$\mu_p = \frac{Z_p}{c}$$
$$\epsilon_p = \frac{1}{c Z_p}$$
$$\Gamma_p = \frac{e^2}{\alpha_p}$$

#### IV. Vacuum Field Parameters (Observable)

These are the standard observable vacuum constants, derived by scaling the primary parameters by the fine structure constant ($).

- **Fine Structure Constant ($\alpha $):**
$$\alpha = \alpha_p \Phi$$

- **Vacuum Impedance ($Z_0$):**
$$Z_0 = \alpha Z_p$$

- **Vacuum Permeability & Permittivity ($\mu_0$, $\epsilon_0$):**
$$\mu_0 = \frac{Z_0}{c}$$
$$\epsilon_0 = \frac{1}{c Z_0}$$

- **Vacuum Gamma ($\Gamma$):**
$$\Gamma = \alpha \Gamma_p$$

#### V. Gravitational Unification

The model derives the Gravitational Constant ($G$) not as a fundamental arbitrary value, but as a result of geometric impedance scaling.

- **Gravitational Constant ($G$):**

$$
G = \frac{Z_0}{c S \phi}
[\frac{m^3}{kg s^2}]
$$

- **Complex Unification Factor ($\Xi$):** A complex rotation relating field geometry to mass-charge equivalence.

$$
\Xi = \sqrt{4\pi \sqrt{2} G \epsilon_0} \left( \cos\frac{\pi}{8} - i \sin\frac{\pi}{8} \right)
[C/kg]
$$

#### VI. Mass & Length Scales

Definitions of the Planck scale and the specific derivation of the Proton Mass.

- **Planck Units ($m_P$, $l_P$):**
$$m_P = \sqrt{\frac{ch}{2\pi G}}$$
$$l_P = \sqrt{\frac{hG}{2\pi c^3}}$$

- **Proton Mass ($m_p$):** A derivation scaling the Planck mass by exponential and geometric corrections.
$$m_p = e^{-14\pi} \left( 1 - \frac{1}{2^{11}-4} \right) m_P \left( 1 - \frac{4 \alpha_\gamma}{\alpha_\delta} \right)$$

- **Mass-Charge Equivalence ($Q$):** A unified charge definition based on mass  and the complex factor .
$$Q = M \Xi$$

![Mathematica Definitions](./images/mathematica_page1.png)

To visualize force scaling (Newtonian vs. GEM-derived, log-log plot from Sun-Earth system):

![GEM](/chart.png)

### 4. Interaction Examples


```shell
cargo test -- --nocapture
```

#### Result

```shell
running 9 tests
--- THE DECODER TEST ---
Gravity Force: 2.5686e-47-2.5686e-47i
Coulomb Force: 8.2424e-8+0.0000e0i
--- GEM HYDROGEN DECODER ---
--- GEM MASS HARMONIC SCANNER ---
--- GEM HYDROGEN INTERACTION ANALYSIS ---
1. Electron Field Energy:   13.60569 eV (Source Potential)
2. Proton Coupling Energy:  0.00741 eV (Interaction Term)
3. Net Observable Energy:   13.59828 eV (Measured State)
   Target:                  13.5983 eV
State n=1: Radius = 5.2947e-11 m, Energy = 13.5983 eV
State n=2: Radius = 2.1179e-10 m, Energy = 3.3996 eV
--- RECOIL DERIVATION (Mathematica In[204-208]) ---
Term 1 (Ideal Me):    13.60569 eV (Target: 13.6057)
Term 2 (Recoil Me/Mp):0.00741 eV (Target: 0.0074)
Difference:           13.59828 eV (Target: 13.5983)
>> HARMONIC 48 (Proton Candidate): Ratio = 1836.1318
State n=3: Radius = 4.7652e-10 m, Energy = 1.5109 eV
--- GEM COMPLEX ENERGY ANALYSIS ---
Real Part:      9.61544 eV
Imaginary Part: 9.61544 eV
Magnitude:      13.59829 eV
Target:         13.5983 eV
Derived Proton/Electron Ratio: 1836.1318
test tests::verify_derivations_and_constants ... ok
test tests::test_dynamic_system_orbit ... ok
test tests::test_decoder_philosophy ... ok
test tests::verify_hydrogen_component_interaction ... ok
test tests::verify_hydrogen_recoil_derivation ... ok
test tests::test_mass_resonance_integers ... ok
test tests::test_hydrogen_resonance_decoder ... ok
test tests::verify_complex_energy_components ... ok
test tests::verify_mass_ratio_derivation ... ok
```

```shell
cargo run --bin electron_proton_action
```

```shell
electron: GeometricKnot {
    name: "Electron",
    mass: 9.1093837139e-31,
    charge: Complex {
        re: -1.602176634e-19,
        im: -1.602176634e-19,
    },
    geometric_radius_a: 5.291772104738458e-11,
    physical_radius: 2.817940320834912e-15,
}
proton: GeometricKnot {
    name: "Proton",
    mass: 1.67262192595e-27,
    charge: Complex {
        re: 1.602176634e-19,
        im: 1.602176634e-19,
    },
    geometric_radius_a: 2.8819891620872986e-14,
    physical_radius: 1.5346982642700954e-18,
}
Result: GemInteractionResult {
    q1: Complex {
        re: 8.624716672621998e-41,
        im: -3.5724746174253846e-41,
    },
    q2: Complex {
        re: 1.5836296575938094e-37,
        im: -6.559608819516165e-38,
    },
    q_total: Complex {
        re: 1.5844921292610717e-37,
        im: -6.56318129413359e-38,
    },
    af1: Complex {
        re: 3.3694829747370465e-13,
        im: -3.3694829747370465e-13,
    },
    af2: Complex {
        re: 6.186885172111745e-10,
        im: -6.186885172111745e-10,
    },
    g1: Complex {
        re: 1.5352496334655168e-20,
        im: -1.5352496334655168e-20,
    },
    g2: Complex {
        re: 2.818952718857127e-17,
        im: -2.818952718857127e-17,
    },
    force: Complex {
        re: 2.567892198741124e-47,
        im: -2.567892198741124e-47,
    },
    curvature: Complex {
        re: 1.0000000000000002,
        im: 0.0,
    },
    g_o: Complex {
        re: 6.674325731836407e-11,
        im: 0.0,
    },
    g_recovered: Complex {
        re: 6.674325731836409e-11,
        im: 0.0,
    },
    mqr1_ev: Complex {
        re: 469136044.7141288,
        im: -469136044.7141288,
    },
    mqr2_ev: Complex {
        re: 255499.47534587674,
        im: -255499.47534587674,
    },
    ratio1: Complex {
        re: 0.9994553833012333,
        im: -0.0005440241291649849,
    },
    ratio2: Complex {
        re: 0.0002723084703826605,
        im: -0.00027216024727941913,
    },
    binding_energy_ev: Complex {
        re: 9.62067644577312,
        im: 9.62067644577312,
    },
    schwarzschild_radius: 2.485594236043889e-54,
    ratio_mr_d: 4.697091670283908e-44,
    is_complex: false,
}
(result.q1 / m2).abs(): 5.58125294864566e-14
(result.q2 / m1).abs(): 1.8816952313861778e-7
 ((result.q2 / m1).abs()/(result.q1 / m2).abs()).sqrt(): 1836.1526734215263
Abs[binding_energy_ev] = 13.60569110881573
```

#### A. Microscopic Interaction (Electron-Proton)

At the quantum scale, the medium handles interactions via the Coulomb protocol (Electric Charge) and the Gravity protocol (Mass/Shadow Charge) simultaneously.

```rust
use gemphy::{GeometricEncodedMedium, GeometricKnot, ForceProtocol};
use num_complex::Complex64;

fn main() {
    let medium = GeometricEncodedMedium::new();
    let q1 = Complex64::new(-medium.e, -medium.e);
    let q2 = Complex64::new(medium.e, medium.e);
    let electron = GeometricKnot::new(medium.clone(), 9.109e-31, q1, 0.0, "Electron");
    let proton = GeometricKnot::new(medium.clone(), 1.672e-27, q2, 0.0, "Proton");
    let d = 5.29e-11; // Bohr Radius

    // Decode as Electromagnetism (Strong Force)
    let f_coulomb = medium.decode_force(&electron, &proton, d, ForceProtocol::Electromagnetism);
    
    // Decode as Gravity (Weak Force)
    let f_gravity = medium.decode_force(&electron, &proton, d, ForceProtocol::Gravity);

    println!("--- Microscopic Interaction ---");
    println!("Coulomb Force: {:.4e} N", f_coulomb); // ~ 8.2e-8 N
    println!("Gravity Force: {:.4e} N", f_gravity); // ~ 3.6e-47 N
}
```

#### B. Macroscopic Interaction (Earth-Sun)

At the cosmic scale, the Gravity protocol dominates.

```rust
use gemphy::{GeometricEncodedMedium, GeometricKnot, ForceProtocol};

fn main() {
    let medium = GeometricEncodedMedium::new();
    let earth = GeometricKnot::new(medium.clone(), 5.972e24, 5.972e24 * medium.xi, 6.371e6, "Earth");
    let sun = GeometricKnot::new(medium.clone(), 1.989e30, 1.989e30 * medium.xi, 6.96e8, "Sun");
    let d = 149.6e9; // 1 AU

    let f_gravity = medium.decode_force(&earth, &sun, d, ForceProtocol::Gravity);

    println!("--- Macroscopic Interaction ---");
    println!("Gravitational Force: {:.4e} N", f_gravity); // ~ 3.5e22 N
}
```

#### C. Neutron Star Interaction (Test Mass-Neutron Star)

This example demonstrates GEM's handling of strong gravitational fields near compact objects.

```rust
// From neutron_star_action binary
// (Code snippet similar to above, with masses 1.0 kg and 2.78376e30 kg, d ≈10 km)
```

Output:

```shell
GeometricKnot {
    name: "Test Mass",
    mass: 1.0,
    charge: Complex {
        re: 9.467930309028682e-11,
        im: -3.921745141602969e-11,
    },
    geometric_radius_a: 1.48523238761875e-27,
    base_length: 7.909063641216323e-32,
}
GeometricKnot {
    name: "Neutron Star",
    mass: 2.78376e30,
    charge: Complex {
        re: 2.6356445677061682e20,
        im: -1.0917197255388681e20,
    },
    geometric_radius_a: 4134.5305113575705,
    base_length: 0.22016935001872348,
}
Result: GemInteractionResult {
    q1: Complex {
        re: 9.467930309028682e-11,
        im: -3.921745141602969e-11,
    },
    q2: Complex {
        re: 2.6356445677061682e20,
        im: -1.0917197255388681e20,
    },
    q_total: Complex {
        re: 2.6356445677061682e20,
        im: -1.0917197255388681e20,
    },
    af1: Complex {
        re: 5.481208021372526e-26,
        im: -5.481208021372526e-26,
    },
    af2: Complex {
        re: 152583.67641575984,
        im: -152583.67641575984,
    },
    g1: Complex {
        re: 4.719443850333546e-19,
        im: -4.719443850333546e-19,
    },
    g2: Complex {
        re: 4.719443850333546e-19,
        im: -4.719443850333546e-19,
    },
    force: Complex {
        re: 1313779901280.451,
        im: -1313779901280.451,
    },
    curvature: Complex {
        re: 0.9238795325112867,
        im: -0.3826834323650898,
    },
    g_o: Complex {
        re: 6.6743015e-11,
        im: 0.0,
    },
    g_recovered: Complex {
        re: 4.719443850333545e-11,
        im: -4.719443850333546e-11,
    },
    mqr1_ev: Complex {
        re: 1.1042001594518223e66,
        im: 0.0,
    },
    mqr2_ev: Complex {
        re: 3.966578151319878e35,
        im: -2.496603636017335e19,
    },
    ratio1: Complex {
        re: 1.0,
        im: 0.0,
    },
    ratio2: Complex {
        re: 2.5401140227122575e-31,
        im: -1.5987729632635265e-47,
    },
    binding_energy_ev: Complex {
        re: 1.4935909608875857e31,
        im: 0.0,
    },
    schwarzschild_radius: 4134.5305113575705,
    ratio_mr_d: 0.41345305113575703,
    is_complex: true,
}
```

### 5. The Unified Force (Complex Phase Engine)

What happens inside a Black Hole or at the Planck Scale? Standard physics breaks down. GEM handles this by rotating the force vector into the **Imaginary Plane**.

- **Real Force:** Linear Acceleration (Push/Pull).
- **Imaginary Force:** Rotational Action (Spin/Memory).

Kappa ($\kappa$) in GEM is similar to Einstein's curvature in general relativity, coupling mass to geometry via the complex unification factor $\Xi$. It derives as $\kappa = \Xi \frac{M + m}{Q_{sh} + q_{sh}}$, where $Q_{sh} = \Xi M$ (shadow charge), simplifying to a phased unit complex in gravity-dominated regimes.

<!-- end list -->

## 🧪 Testing the Laws of Physics

This library includes a rigorous test suite that validates the framework against CODATA observations.

Run the verification suite:

```bash
cargo test
```

**What is tested?**

- **The Golden Loop:** Verifies $\frac{h}{2\pi c} \equiv \frac{\Gamma}{\alpha}$.
  - **Grand Unification:** Verifies the geometric scaling from Electrons to the Universe.
  - **Force Parity:** Checks that GEM-derived Gravity matches Newtonian Gravity to $10^{-13}$ precision.
  - **Phase Change:** Ensures forces transition to complex numbers correctly inside event horizons.

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for details on issues, PRs, and development setup.

## ⚖️ License

This project is licensed under GPL-3.0.

---
