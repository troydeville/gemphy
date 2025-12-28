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

### 3\. Interaction Examples

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

### 4\. The Unified Force (Complex Phase Engine)

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
