# CA-PDE Tumor Growth Simulator

A hybrid **Cellular Automata (CA)–Partial Differential Equation (PDE)** framework
for simulating tumor growth driven by nutrient diffusion, consumption, and
stochastic cell behavior.

This project was developed as part of a **Data Structures & Algorithms (DSA)**
course project.

---

## 🧬 Overview

Tumor growth is modeled on a two-dimensional computational grid representing
biological tissue, where each grid cell corresponds to a biological cell.
The tumor evolves through local interaction rules while being globally regulated
by nutrient diffusion.

The combination of discrete cellular automata and continuous nutrient dynamics
allows the model to reproduce realistic tumor morphology.

---

## ⚙️ Model Components

### Cellular Automata (CA)
Each grid cell can exist in one of the following states:
- **Empty**
- **Tumor (Proliferating)**
- **Quiescent**
- **Necrotic**

State transitions depend on:
- Local nutrient concentration
- Moore (8-neighbor) neighborhood
- Probabilistic cell division

### Nutrient Dynamics (PDE)
Nutrient availability is modeled using a reaction–diffusion PDE:
- Nutrients diffuse through tissue
- Tumor and quiescent cells consume nutrients
- Boundary conditions represent blood supply

This naturally generates nutrient gradients that regulate tumor growth.

---

## 🧪 Emergent Behavior

Without explicitly defining tumor shape, the model exhibits:
- Proliferating outer rim
- Quiescent intermediate layer
- Necrotic core
- Irregular growth patterns due to stochastic effects
These structures emerge purely from local rules and diffusion dynamics.

---

## 👥 Contributors

- **Sidhhesh Phadke**
- **Sahasra Oleti**
- **Atharva Shetwe**
