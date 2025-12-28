# Geometric Encoded Medium (GEM) Framework

##

### A Unified Impedance Framework for Physics in Rust

> "The universe is a perfect geometric circuit."

![Version](https://img.shields.io/badge/version-0.1.0-blue) ![License](https://img.shields.io/badge/license-MIT%2FApache-green) ![Rust](https://img.shields.io/badge/built_with-Rust-orange)

---

## 🌌 Overview

- The **GEM Framework** is a Rust library that models reality not as a collection of arbitrary forces and constants, but as a single **Geometric Encoded Medium**.

- It posits that space itself has impedance ($Z_p$) and geometry (**Horn Torus topology**), and that all physical phenomena—Gravity, Electromagnetism, Mass, and Charge—are simply different "encodings" running on this geometric hardware.

- This library does not approximate physics; it **derives** it. By defining a few geometric axioms, it mathematically derives Newton's Gravitational Constant ($G$), the Fine Structure Constant ($\alpha$), and the Proton Radius from first principles.

## 🧠 Core Philosophy

In GEM, the universe is treated as a software system:

- **The Hardware:** A "Horn Torus" topology representing the vacuum medium.
- **The Operating System:** The Impedance Field ($Z_0, Z_p$) that dictates how information moves.
- **The Software:**
- **Gravity:** The interaction when the medium is encoded with Mass (Shadow Charge $Q = \Xi M$).
- **Electromagnetism:** The interaction when the medium is encoded with Electric Charge ($q$).

---

## Roadmap

Milestone | v0.3 | v0.4 | v1.0
--- | --- | --- | ---
Tasks | Add horn torus simulation; publish crate | Interactive CLI; Python bindings | Full predictions (e.g., dark matter as impedance)

---

## 📦 Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
gemometry = "0.1.0"
```

> **Note:** The library re-exports `Complex64` from `num-complex`, so you do not need to install it separately unless otherwise needed.

---

## 🚀 Usage Guide

### 1\. Initialize the Medium

The `GeometricEncodedMedium` struct initializes the fundamental constants derived from the GEM geometric axioms.

```rust
use gemometry::GeometricEncodedMedium;

fn main() {
    let medium = GeometricEncodedMedium::new();

    println!("--- GEM Fundamental Derivations ---");
    println!("Vacuum Impedance (Z0): {:.4} Ohms", medium.z_o);
    println!("Fine Structure (α):    {:.4e}", medium.alpha);
    println!("Gravitational G:       {:.4e}", medium.g); // Derived, not hardcoded!
}
```

### 2\. Create Geometric Knots

The `GeometricKnot` represents any entity in the universe. It automatically detects its scale (Quantum vs. Macro) and sizes its geometric "vacuum knot" accordingly.

```rust
use gemometry::GeometricKnot;

fn main() {
    // 1. Create an Electron (Quantum Scale)
    // Mass: 9.109e-31 kg, Charge: -1.6e-19 C
    let electron = GeometricKnot::new(9.109e-31, -1.602e-19, 0.0, "Electron");
    
    // 2. Create the Earth (Macro Scale)
    // Mass: 5.972e24 kg, Charge: 0 (Neutral, calculates Shadow Charge internally)
    let earth = GeometricKnot::new(5.972e24, 0.0, 6.371e6, "Earth");

    println!("--- Particle Topology ---");
    electron.print_analysis();
    earth.print_analysis();
}
```

> **Output Note:** The Electron will automatically size to the **Bohr Radius**, while the Earth sizes to its **Schwarzschild Radius**.

Here is the text converted into a descriptive, documented format. I have grouped the equations into logical sections (Fundamental, Field Theory, Gravity, and Mass) to explain the physical significance of each variable in the context of this unified model.

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
[\frac{m^3}{kg s^2}]
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

### Next Step

Would you like me to generate a dimensional analysis table to verify that the units in the equation  correctly resolve to ?
<!-- 
### 3\. Constants

$$h = \frac{662607015}{10^{42}} \text{ J/Hz}$$

$$c = 299792458 \text{ m/s}$$

$$e = \frac{1602176634}{10^{28}} \text{ C}$$

$$\Phi = \frac{1}{10^7} \text{ H/m}$$

$$\phi = 10^4 \text{ kg}^2 \text{ m}^{-2} \text{ s}^2 \text{ C}^{-2}$$

$$\Lambda_\alpha(n) = \frac{\alpha_\delta n}{\alpha_\gamma + 2 \alpha_\delta n^2}$$

$$Z_p = \frac{2h}{e^2}$$

$$\alpha_p = \frac{4\pi c}{Z_p}$$

$$\mu_p = \frac{Z_p}{c}$$

$$\epsilon_p = \frac{1}{c Z_p}$$

$$\Gamma_p = \frac{e^2}{\alpha_p}$$

$$\alpha = \alpha_p \Phi$$

$$Z_0 = \alpha Z_p$$

$$\mu_0 = \frac{Z_0}{c}$$

$$\epsilon_0 = \frac{1}{c Z_0}$$

$$\Gamma = \alpha \Gamma_p$$

$$G = \frac{Z_0}{c S \phi}$$

$$\Xi = \sqrt{4\pi \sqrt{2} G \epsilon_0} \left( \cos\frac{\pi}{8} - i \sin\frac{\pi}{8} \right)$$

$$m_P = \sqrt{\frac{ch}{2\pi G}}$$

$$l_P = \sqrt{\frac{hG}{2\pi c^3}}$$

$$m_p = e^{-14\pi} \left( 1 - \frac{1}{2^{11}-4} \right) m_P \left( 1 - \frac{4 \alpha_\gamma}{\alpha_\delta} \right)$$

$$Q = M \Xi$$ -->

### 4\. Interaction Examples

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
use gemometry::{GeometricEncodedMedium, GeometricKnot, ForceProtocol};

fn main() {
    let medium = GeometricEncodedMedium::new();
    let electron = GeometricKnot::new(9.109e-31, -1.602e-19, 0.0, "Electron");
    let proton = GeometricKnot::new(1.672e-27, 1.602e-19, 0.0, "Proton");
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
use gemometry::{GeometricEncodedMedium, GeometricKnot, ForceProtocol};

fn main() {
    let medium = GeometricEncodedMedium::new();
    let earth = GeometricKnot::new(5.972e24, 0.0, 6.371e6, "Earth");
    let sun = GeometricKnot::new(1.989e30, 0.0, 6.96e8, "Sun");
    let d = 149.6e9; // 1 AU

    let f_gravity = medium.decode_force(&earth, &sun, d, ForceProtocol::Gravity);

    println!("--- Macroscopic Interaction ---");
    println!("Gravitational Force: {:.4e} N", f_gravity); // ~ 3.5e22 N
}
```

### 5\. The Unified Force (Complex Phase Engine)

What happens inside a Black Hole or at the Planck Scale? Standard physics breaks down. GEM handles this by rotating the force vector into the **Imaginary Plane**.

- **Real Force:** Linear Acceleration (Push/Pull).
- **Imaginary Force:** Rotational Action (Spin/Memory).

<!-- end list -->

```rust
use gemometry::{GeometricEncodedMedium, GeometricKnot};

fn main() {
    let medium = GeometricEncodedMedium::new();
    let mp = medium.m_p; // Planck Mass
    
    // Create two particles at the Planck Scale
    let p1 = GeometricKnot::new(mp, 0.0, 0.0, "Planck A");
    let p2 = GeometricKnot::new(mp, 0.0, 0.0, "Planck B");
    
    // Distance inside the Schwarzschild radius
    let d = 0.5 * medium.l_p; 

    let f_total = medium.calculate_total_force(&p1, &p2, d);

    println!("--- Unified Force (Planck Scale) ---");
    println!("Force Vector: {:.4e} + {:.4e}i N", f_total.re, f_total.im);
    
    if f_total.im.abs() > 0.0 {
        println!("Status: Phase Transition (Information Storage / Spin)");
    }
}
```

---

## 🧪 Testing the Laws of Physics

This library includes a rigorous test suite that validates the framework against CODATA observations.

Run the verification suite:

```bash
cargo test
```

**What is tested?**

- **The Golden Loop:** Verifies $h/2\pi c \equiv \Gamma/\alpha$.
  - **Grand Unification:** Verifies the geometric scaling from Electrons to the Universe.
  - **Force Parity:** Checks that GEM-derived Gravity matches Newtonian Gravity to $10^{-13}$ precision.
  - **Phase Change:** Ensures forces transition to complex numbers correctly inside event horizons.

---

## ⚖️ License

This project is an implementation of the GEM Framework.
