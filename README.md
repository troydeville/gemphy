# Geometric Encoded Medium (GEM) Framework

##

### A Geometric Encoded Medium (GEM) Impedance Framework for Physics in Rust

> "The universe is a perfect geometric circuit."

![Version](https://img.shields.io/badge/version-0.1.1-blue) ![License](https://img.shields.io/badge/license-GPL--3.0-green) ![Rust](https://img.shields.io/badge/built_with-Rust-orange) [![Crates.io](https://img.shields.io/crates/v/gemphy.svg)](https://crates.io/crates/gemphy) [![Docs](https://docs.rs/gemphy/badge.svg)](https://docs.rs/gemphy) [![CI](https://github.com/troydeville/gemphy/workflows/CI/badge.svg)](https://github.com/troydeville/gemphy/actions)

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

## Roadmap

| Milestone | Tasks
|-----------|-------
| v0.3      | Add horn torus simulation; publish crate
| v0.4      | Interactive CLI; Python bindings
| v1.0      | Full predictions (e.g., dark matter as impedance)

---

## 📦 Installation

```bash
cargo add gemphy
```

Install via Cargo:

```toml
[dependencies]
gemphy = "0.1.1"
```

Note: Re-exports `Complex64` from `num-complex`.

---

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

### 4. Interaction Examples

### 5. The Unified Force (Complex Phase Engine)

What happens inside a Black Hole or at the Planck Scale? Standard physics breaks down. GEM handles this by rotating the force vector into the **Imaginary Plane**.

- **Real Force:** Linear Acceleration (Push/Pull).
- **Imaginary Force:** Rotational Action (Spin/Memory).

Kappa ($\kappa$) in GEM is similar to Einstein's curvature in general relativity, coupling mass to geometry via the complex unification factor $\Xi$. It derives as $\kappa = \Xi \frac{M + m}{Q_{sh} + q_{sh}}$, where $Q_{sh} = \Xi M$ (shadow charge), simplifying to a phased unit complex in gravity-dominated regimes.

<!-- end list -->

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for details on issues, PRs, and development setup.

## ⚖️ License

This project is licensed under GPL-3.0.

---
