# Geometric Encoded Medium (GEM) Physics Framework

![Version](https://img.shields.io/badge/version-0.2.2-blue) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![Rust](https://img.shields.io/badge/built_with-Rust-orange) [![Crates.io](https://img.shields.io/crates/v/gemphy.svg)](https://crates.io/crates/gemphy) [![Docs](https://docs.rs/gemphy/badge.svg)](https://docs.rs/gemphy) [![CI](https://github.com/troydeville/gemphy/workflows/CI/badge.svg)](https://github.com/troydeville/gemphy/actions)

---

## 🌌 Overview

The **GEM Framework** is a Rust library that models reality as a single **Geometric Encoded Medium**. It posits that space-time is a medium with intrinsic impedance ($Z_p$) and a **Horn Torus topology**. All physical phenomena—Gravity, Electromagnetism, Mass, and Charge—are derived as specific "encodings" on this geometric hardware.

Unlike traditional physics engines that rely on curve-fitting or arbitrary constants, GEM **derives** fundamental values like Newton's Gravitational Constant ($G$), the Fine Structure Constant ($\alpha$), and the Proton Radius from first-principles geometric axioms.

## 🧠 Core Philosophy: The Unified Phase Engine

In GEM, the universe operates as a geometric circuit:

* **The Hardware:** A Horn Torus manifold ($R = r$) where the singularity at the center acts as a topological pump.
* **The Phase:** Interactions are not just magnitudes; they are complex vectors. Extreme scales (Planck, Black Holes, Muonic states) rotate real linear forces into **Imaginary/Rotational components**.
* **The Whirl:** Particle "Spin" is visualized as a **Traveling Wave** of energy density orbiting the singularity, rather than a point-particle rotating in a void.

## Roadmap

| Milestone | Tasks |
|-----------|-------|
| v0.3      | Add horn torus simulation; enhanced orbit sims (hydrogen, muonic, neutron star) |
| v0.4      | Interactive CLI; Python bindings; initial multi-body support |
| v0.5      | gemphy-web deployment with live visualizations; basic n-body examples (e.g., Sun-Earth-Moon) |
| v1.0      | Full predictions (e.g., dark matter as impedance, novel deviations in fine structure); arXiv preprint with derivations |

| Milestone | Tasks
|-----------|-------
| v0.3      | Add horn torus simulation
| v0.4      | Interactive CLI
| v1.0      | Full predictions (e.g., dark matter as impedance)

---

## 🚀 Quick Start

## 📦 Installation

```toml
[dependencies]
gemphy = "0.2.2"
```

Note: Re-exports `Complex64` from `num-complex`.

---

## Example Test Usage

## View tests

```bash
cargo run --bin orbit_sim
```

```bash
cargo run --bin electron_proton_action
```

```bash
cargo run --bin muon_proton_action
```

```bash
cargo run --bin neutron_star_action
```

---

### Stable Orbit Simulation

GEMPHY handles the "Dynamic Stability" of orbits by calculating the interaction between geometric knots.

```rust
use gemphy::{knot::GeometricKnot, medium::{GAMMA_P, GeometricEncodedMedium}};
use physical_constants::{ELECTRON_MASS, PROTON_MASS, ELEMENTARY_CHARGE};

fn main() {
    let m1 = ELECTRON_MASS;
    let m2 = PROTON_MASS;
    
    let medium = GeometricEncodedMedium::new();

    // Setup Particles as Geometric Knots
    let electron = GeometricKnot::new(medium.clone(), m1, &[-1.0], 0.0, "Electron");
    let proton = GeometricKnot::new(medium.clone(), m2, &[1.0], 0.0, "Proton");

    let rg1 = (GAMMA_P / (electron.mass * medium.alpha)).powi(2);
    let rg2 = (GAMMA_P / (proton.mass * medium.alpha)).powi(2);
    let d = (rg1+rg2).sqrt();

    let interaction = medium.calculate_interaction(&electron, &proton, d.into());
    
    println!("Result:             {:#?}", interaction);
    println!("Er (eV):            {:#?}", interaction.er1.norm()/ ELEMENTARY_CHARGE);
    println!("Ei (eV):            {:#?}", interaction.ei1.norm()/ ELEMENTARY_CHARGE);
    println!("E  (eV):            {:#?}", interaction.binding_energy.norm()/ ELEMENTARY_CHARGE);
}

```

---

## 📐 The Geometric Model

### I. Fundamental Scaling

The framework uses a fundamental geometric normalization constant to bridge the subatomic and cosmic scales:

* **Normalization Constant ($S$):**
$$({4 \pi})^{1/4} \approx 1.8827925275534296$$

* **Mass-Charge Metric ($\phi$):**
$$\phi = 10^4 \text{ kg}^2 \text{ m}^{-2} \text{ s}^2 \text{ C}^{-2}$$

* **Magnetic Scaling ($\Phi$):**
$$\Phi = \frac{1}{10^7} \text{ H/m}$$

* **Primary Impedance ($Z_p$):**
$$Z_p = \frac{2h}{e^2} \Omega$$

#### Fine Structure ($\alpha$) Scaling Relationships

* **Primary Fine Structure ($\alpha_p $):**
$$\alpha_p = \frac{4\pi c}{Z_p}$$

* **Fine Structure ($\alpha $):**
$$\alpha = \frac{4\pi c}{Z_p} \Phi$$

* **Impedance ($Z_p$, $Z_o$):**
$$\alpha Z_p \implies Z_0$$

<!-- * **Primary Impedance ($Z_p$):**
$$Z_p = \frac{2h}{e^2} \implies \alpha Z_p = Z_0$$

* **Primary Fine Structure ($\alpha_p $):**
$$\alpha_p = \frac{4\pi c}{Z_p}$$ -->

* **Permeability ($\mu_p $):**
$$\mu_p = \frac{Z_p}{c} \implies \mu_0 = \alpha \mu_p$$
<!-- Scaled with ($\alpha$):
$$\alpha \mu_p \Rightarrow \mu_0$$ -->
<!-- $$\frac{1}{c} (\alpha Z_p) \implies \mu_0$$ -->

* **Permittivity ($\epsilon_p$)**
$$\epsilon_p = \frac{1}{c Z_p} \implies \epsilon_0 = \frac{1}{\alpha} \epsilon_p $$

<!-- \implies \epsilon_0$$
$$\frac{1}{c} (\frac{1}{\alpha Z_p}) \implies \epsilon_0$$ -->

* **Gamma factor ($\Gamma_p$)**
$$\Gamma_p = \frac{e^2}{\alpha_p} \implies \Gamma = \frac{\Gamma_p}{\alpha}$$

### II. Gravitational Unification

GEM derives $G$ as a result of geometric impedance scaling rather than an empirical measurement:
Where $Z_0$ is the vacuum impedance and $S$ is the geometric shape factor:
$$
G = \frac{Z_0}{c S \phi}
[\frac{m^3}{kg s^2}]
$$

<!-- **Complex Unification Factor ($\Xi$)**

A complex rotation relating field geometry to mass-charge equivalence. -->

### III. Complex Geometry

**A complex rotation relating field geometry to mass-charge equivalence. ($\Xi$)**
$$\Xi = \sqrt{4\pi \sqrt{2} G \epsilon_0} \left( \cos\frac{\pi}{8} - i \sin\frac{\pi}{8} \right)[C/kg]$$

Force is calculated as a complex vector.

#### Linked Complex Potentials ($E_r$, $E_i$)

Energy interaction in GEM is calculated as the sum of two body-specific complex potentials. Each potential represents the geometric "tension" localized to that body within the medium:

* **$E_r$:** Complex potential of Body 1 (e.g., the orbiting mass).
* **$E_i$:** Complex potential of Body 2 (e.g., the central mass).
* **Total Interaction Energy ($E_r + E_i$):** .

The relationship between $E_r$ and $E_i$ is intrinsically linked by the medium's impedance. As distance or energy density changes, these values rotate in the complex plane, representing the transition from linear work to orbital/spin action.

---

## 🛠 Visualization: The Horn Torus Manifold

The provided `HornTorusManifold` component (for Three.js/React) is a **Raw Scientific Viewer** designed to mirror the Rust engine's output:

* **Traveling Waves:** Total Energy ($E$) is treated as a complex phase. This phase drives a wave that "chases its tail" around the torus, visually representing particle spin and momentum.
* **Singularity Flow:** The mesh geometry respects the Horn Torus topology ($R=r$), creating a natural "topological pump" at the center that breathes the vacuum.
* **No Fakes:** All surface deformations are driven by the `Complex64` results of the physical interaction. If the engine calculates zero energy, the manifold remains stagnant.

---

## ⚖️ License

Licensed under the **MIT License**.

> **Scientific Attribution:** If you use GEMPHY in a research paper or commercial simulation, please cite the framework to preserve the geometric integrity of the medium.

---
