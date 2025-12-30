# Geometric Encoded Medium (GEM) Framework

##

### A Geometric Encoded Medium (GEM) Impedance Framework for Physics in Rust

#### The universe is a perfect geometric circuit.

---

## 🌌 Overview

The **GEM Framework** is a Rust library that models reality not as a collection of arbitrary forces and constants, but as a single **Geometric Encoded Medium**. It posits that space itself has impedance () and geometry (**Horn Torus topology**), and that all physical phenomena—Gravity, Electromagnetism, Mass, and Charge—are simply different "encodings" running on this geometric hardware.

This library does not approximate physics; it **derives** it. By defining a few geometric axioms, it mathematically derives Newton's Gravitational Constant (), the Fine Structure Constant (), and the Proton Radius from first principles. Validations match CODATA with high precision (e.g.,  rel error -5.47e-10,  rel error 3.86e-6).

For a visual of the horn torus topology:

## 🧠 Core Philosophy

In GEM, the universe is treated as a software system:

* **The Hardware:** A "Horn Torus" topology representing the vacuum medium.
* **The Operating System:** The Impedance Field () that dictates how information moves.
* **The Software:**
* **Gravity:** The interaction when the medium is encoded with Mass (Shadow Charge ).
* **Electromagnetism:** The interaction when the medium is encoded with Electric Charge ().



Complex phases handle extreme scales (e.g., Planck or black holes), rotating real forces into imaginary (rotational/spin) components.

## Roadmap

| Milestone | Tasks |
| --- | --- |
| v0.3 | Add horn torus simulation; publish crate |
| v0.4 | Interactive CLI; Python bindings |
| v1.0 | Full predictions (e.g., dark matter as impedance) |

---

## 📦 Installation

```bash
cargo add gemphy
```

Install via Cargo:

```toml
[dependencies]
gemphy = "0.2.0"
```

Note: Re-exports `Complex64` from `num-complex`.

---

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

### 3. Constants & Model Specification

#### I. Fundamental Physical Constants

These are the fixed, integer-defined values representing the bedrock of the SI system used in the model.

* **Planck Constant ():** The quantum of action.


* **Speed of Light ():** The universal speed limit.


* **Elementary Charge ():** The fundamental unit of charge.



#### II. Scaling & Geometric Factors

These parameters act as scaling bridges between the quantum/electromagnetic scale and the gravitational scale.

* **Magnetic Scaling Constant ():** A scaling factor with dimensions of inductance per meter.


* **Mass-Charge Metric ():** A density-like scalar connecting mass and charge squared.


* **Alpha-Lambda Function ():** A rational function governing specific integer-ratio field interactions.



#### III. Primary Field Parameters (Planck Scale)

These variables (-subscript) represent the theoretical "maximum" impedance and field density before fine-structure scaling is applied.

* **Primary Impedance ():**


* **Primary Fine Structure ():** The geometric inverse of impedance.


* **Primary Field Constants:** Permeability (), Permittivity (), and the Gamma factor ().





#### IV. Vacuum Field Parameters (Observable)

These are the standard observable vacuum constants, derived by scaling the primary parameters by the fine structure constant ().

* **Fine Structure Constant ():**


* **Vacuum Impedance ():**


* **Vacuum Permeability & Permittivity (, ):**



* **Vacuum Gamma ():**



#### V. Gravitational Unification

The model derives the Gravitational Constant () not as a fundamental arbitrary value, but as a result of geometric impedance scaling.

* **Gravitational Constant ():**

* **Complex Unification Factor ():** A complex rotation relating field geometry to mass-charge equivalence.


#### VI. Mass & Length Scales

Definitions of the Planck scale and the specific derivation of the Proton Mass.

* **Planck Units (, ):**



* **Proton Mass ():** A derivation scaling the Planck mass by exponential and geometric corrections.


* **Mass-Charge Equivalence ():** A unified charge definition based on mass  and the complex factor .



### 4. Interaction Examples

### 5. The Unified Force (Complex Phase Engine)

What happens inside a Black Hole or at the Planck Scale? Standard physics breaks down. GEM handles this by rotating the force vector into the **Imaginary Plane**.

* **Real Force:** Linear Acceleration (Push/Pull).
* **Imaginary Force:** Rotational Action (Spin/Memory).

Kappa () in GEM is similar to Einstein's curvature in general relativity, coupling mass to geometry via the complex unification factor . It derives as , where  (shadow charge), simplifying to a phased unit complex in gravity-dominated regimes.

---

## Contributing

See [CONTRIBUTING.md](https://www.google.com/search?q=CONTRIBUTING.md) for details on issues, PRs, and development setup.

## ⚖️ License

This project is licensed under GPL-3.0.

---