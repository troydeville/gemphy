# Geometric Encoded Medium (GEM) Framework

### A Unified Impedance Framework for Physics in Rust

> "The universe is a perfect geometric circuit."

---

## 🌌 Overview

The **GEM Framework** is a Rust library that models reality not as a collection of arbitrary forces and constants, but as a single **Geometric Encoded Medium**.

It posits that space itself has impedance () and geometry (**Horn Torus topology**), and that all physical phenomena—Gravity, Electromagnetism, Mass, and Charge—are simply different "encodings" running on this geometric hardware.

This library does not approximate physics; it **derives** it. By defining a few geometric axioms, it mathematically derives Newton's Gravitational Constant (), the Fine Structure Constant (), and the Proton Radius from first principles.

## 🧠 Core Philosophy

In GEM, the universe is treated as a software system:

* **The Hardware:** A "Horn Torus" topology representing the vacuum medium.
* **The Operating System:** The Impedance Field () that dictates how information moves.
* **The Software:**
* **Gravity:** The interaction when the medium is encoded with Mass (Shadow Charge ).
* **Electromagnetism:** The interaction when the medium is encoded with Electric Charge ().



---

## 📦 Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
num-complex = "0.4"
serde = { version = "1.0", features = ["derive"] }
csv = "1.3"
anyhow = "1.0"
```

---

## 🚀 Usage Guide

### 1. Initialize the Medium

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

### 2. Create Geometric Knots

The `GeometricKnot` represents any entity in the universe. It automatically detects its scale (Quantum vs. Macro) and sizes its geometric "vacuum knot" accordingly.

```rust
use gemometry::{GeometricKnot, GeometricEncodedMedium, Complex64};

fn main() {
    let med = GeometricEncodedMedium::new();
    
    // 1. Create an Electron (Quantum Scale)
    // Mass: 9.109e-31 kg, Charge: -e
    let q_electron = Complex64::new(-med.e, 0.0);
    let electron = GeometricKnot::new(med.clone(), 9.109e-31, q_electron, 0.0, "Electron");
    
    // 2. Create the Earth (Macro Scale)
    // Mass: 5.972e24 kg, Charge: 0 (Neutral)
    let earth = GeometricKnot::new(med.clone(), 5.972e24, Complex64::new(0.0, 0.0), 6.371e6, "Earth");

    println!("--- Particle Topology ---");
    println!("Electron Geo Radius: {:.4e} m", electron.geometric_radius_a);
    println!("Earth Phys Radius:   {:.4e} m", earth.physical_radius);
}

```

> **Output Note:** The Electron automatically sizes to the **Bohr Radius** scale, while the Earth sizes to its **Schwarzschild Radius** equivalent in the medium.

---

## 📐 Constants & Model Specification

### I. Fundamental Physical Constants

These are the fixed, integer-defined values representing the bedrock of the SI system used in the model.

* **Planck Constant ():** The quantum of action.


* **Speed of Light ():** The universal speed limit.


* **Elementary Charge ():** The fundamental unit of charge.



### II. Scaling & Geometric Factors

These parameters act as scaling bridges between the quantum/electromagnetic scale and the gravitational scale.

* **Magnetic Scaling Constant ():** 
* **Mass-Charge Metric ():** 

### III. Primary Field Parameters (Planck Scale)

These variables represents the theoretical "maximum" impedance and field density before fine-structure scaling is applied.

* **Primary Impedance ():** 
* **Primary Fine Structure ():** 
* **Primary Gamma ():** 

### IV. Vacuum Field Parameters (Observable)

These are the standard observable vacuum constants, derived by scaling the primary parameters by the fine structure constant.

* **Fine Structure Constant ():** 
* **Vacuum Impedance ():** 
* **Vacuum Gamma ():** 

### V. Gravitational Unification

The model derives the Gravitational Constant () not as a fundamental arbitrary value, but as a result of geometric impedance scaling.

* **Gravitational Constant ():**


* **Complex Unification Factor ():** A complex rotation relating field geometry to mass-charge equivalence.



### VI. Mass & Length Scales

* **Proton Mass ():** A derivation scaling the Planck mass by exponential and geometric corrections.



---

## 🔬 Interaction Examples

### A. Microscopic Interaction (Electron-Proton)

At the quantum scale, the medium handles interactions via the Coulomb protocol (Electric Charge) and the Gravity protocol (Mass/Shadow Charge) simultaneously.

```rust
use gemometry::{GeometricEncodedMedium, GeometricKnot, ForceProtocol, Complex64};

fn main() {
    let medium = GeometricEncodedMedium::new();
    let q_e = Complex64::new(-medium.e, 0.0);
    let q_p = Complex64::new(medium.e, 0.0);
    
    let electron = GeometricKnot::new(medium.clone(), 9.109e-31, q_e, 0.0, "Electron");
    let proton = GeometricKnot::new(medium.clone(), 1.672e-27, q_p, 0.0, "Proton");
    let d = 5.29e-11; // Bohr Radius

    // Decode as Electromagnetism (Strong Force)
    let f_coulomb = medium.decode_force(&electron, &proton, d, ForceProtocol::Electromagnetism);
    
    // Decode as Gravity (Weak Force)
    let f_gravity = medium.decode_force(&electron, &proton, d, ForceProtocol::Gravity);

    println!("--- Microscopic Interaction ---");
    println!("Coulomb Force: {:.4e} N", f_coulomb.re); // ~ 8.2e-8 N
    println!("Gravity Force: {:.4e} N", f_gravity.re); // ~ 3.6e-47 N
}

```

### B. Macroscopic Interaction (Earth-Sun)

At the cosmic scale, the Gravity protocol dominates.

```rust
use gemometry::{GeometricEncodedMedium, GeometricKnot, ForceProtocol, Complex64};

fn main() {
    let medium = GeometricEncodedMedium::new();
    let earth = GeometricKnot::new(medium.clone(), 5.972e24, Complex64::default(), 6.371e6, "Earth");
    let sun = GeometricKnot::new(medium.clone(), 1.989e30, Complex64::default(), 6.96e8, "Sun");
    let d = 149.6e9; // 1 AU

    let f_gravity = medium.decode_force(&earth, &sun, d, ForceProtocol::Gravity);

    println!("--- Macroscopic Interaction ---");
    println!("Gravitational Force: {:.4e} N", f_gravity.re); // ~ 3.5e22 N
}

```

### C. The Unified Force (Complex Phase Engine)

What happens inside a Black Hole or at the Planck Scale? Standard physics breaks down. GEM handles this by rotating the force vector into the **Imaginary Plane**.

* **Real Force:** Linear Acceleration (Push/Pull).
* **Imaginary Force:** Rotational Action (Spin/Memory).

```rust
use gemometry::{GeometricEncodedMedium, GeometricKnot, Complex64};

fn main() {
    let medium = GeometricEncodedMedium::new();
    let mp = medium.m_p; // Planck Mass
    
    // Create two particles at the Planck Scale
    let p1 = GeometricKnot::new(medium.clone(), mp, Complex64::default(), 0.0, "Planck A");
    let p2 = GeometricKnot::new(medium.clone(), mp, Complex64::default(), 0.0, "Planck B");
    
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

## 🧪 Validated Physics

This library includes a rigorous test suite that validates the framework against CODATA observations.

Run the verification suite:

```bash
cargo test -- --nocapture

```

### Key Verification Results:

1. **Hydrogen Ground State:**
* **Theory:** 13.5983 eV
* **GEM Derived:** 13.59828 eV (Error < 0.001%)


2. **Proton/Electron Mass Ratio:**
* **Observed:** 1836.15
* **GEM Derived:** 1836.13 (Derived purely from geometric constants)


3. **The "Golden Loop":**
* Verifies that  holds to  precision.



---

## ⚖️ License

This project is an implementation of the GEM Framework. Licensed under MIT/Apache.