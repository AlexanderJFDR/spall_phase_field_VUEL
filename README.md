# Spall Phase-Field Simulation Using Abaqus VUEL

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This project implements a **phase-field model for spall fracture** using Abaqus VUEL (User-defined Element) subroutines in FORTRAN. It simulates dynamic spall fracture under shock loading conditions, providing a robust numerical framework for predicting damage initiation and evolution in ductile materials.

The theoretical basis is the phase-field approach to spall fracture, as described in:

> Zhang, H., Peng, H., Pei, X. Y., Li, P., Tang, T. G., & Cai, L. C. (2021). A phase-field model for spall fracture. *Journal of Applied Physics*, 129(12).

---

## Features

- **Phase-field fracture modeling** for spallation under high-velocity impact / shock loading
- **Abaqus VUEL subroutine** implementation with full user-element capabilities
- **Plasticity coupling** — the phase-field evolution is coupled with plastic deformation
- **Double precision** solving for numerical stability and accuracy in explicit dynamics
- **Three test cases** are provided:
  - **1D** — 1D bar spall simulation
  - **2D_A** — 2D plate impact (configuration A)
  - **2D_B** — 2D plate impact (configuration B)

---

## Repository Structure

```
spall_phase_field_VUEL/
├── Code Files/
│   ├── 1D/
│   │   ├── pfczm_dynamic_vuel_plastic.for   # VUEL subroutine
│   │   └── 1D_plastic.inp                   # Abaqus input file
│   ├── 2D_A/
│   │   ├── pfczm_dynamic_vuel_plastic.for   # VUEL subroutine
│   │   └── 2D_plastic.inp                   # Abaqus input file
│   └── 2D_B/
│       ├── pfczm_dynamic_vuel_plastic.for   # VUEL subroutine
│       └── 2D_plastic1.inp                  # Abaqus input file
├── Free Surface Velocity Data/
│   ├── 1D/   ── 速度.csv (velocity data)
│   ├── 2D_A/ ── 速度.csv (velocity data)
│   └── 2D_B/ ── 速度.csv (velocity data)
├── Videos/
│   ├── 2D_A/
│   │   ├── 压强.avi (pressure)
│   │   ├── 损伤.avi (damage)
│   │   └── 速度.avi (velocity)
│   └── 2D_B/
│       ├── 压强.avi (pressure)
│       ├── 损伤.avi (damage)
│       └── 速度.avi (velocity)
├── LICENSE
└── README.md
```

### Code Files

Each case folder contains:

- **`pfczm_dynamic_vuel_plastic.for`** — The FORTRAN source code implementing the phase-field cohesive zone model (PFCZM) as an Abaqus VUEL subroutine, with dynamic plastic coupling.
- **`*.inp`** — Abaqus input files defining the mesh, boundary conditions, material properties, and analysis settings for each specific test configuration.

### Free Surface Velocity Data

The `速度.csv` files contain the computed free surface velocity history data for each case, which is a key output for validating spall simulations against experimental results (e.g., VISAR measurements).

### Videos

Visualization results for the 2D cases, including:

- **Pressure** (`压强.avi`) — evolution of pressure fields
- **Damage** (`损伤.avi`) — phase-field damage variable evolution
- **Velocity** (`速度.avi`) — velocity field evolution

---

## Prerequisites

- **Abaqus** (version 2016 or newer recommended) with a FORTRAN compiler configured for user subroutines
- A FORTRAN compiler compatible with your Abaqus version (e.g., Intel Fortran Compiler)

---

## Usage

### 1. Running in Abaqus/CAE

1. Open or create a job in Abaqus/CAE
2. Go to **Job → Edit → Precision** tab
3. Under **Abaqus/Explicit precision**, select **Double – XX** 
4. In the **General** tab, specify the VUEL subroutine file (`pfczm_dynamic_vuel_plastic.for`)
5. Submit the job

### 2. Running via Command Line

```bash
abaqus job=<job_name> user=pfczm_dynamic_vuel_plastic.for cpus=<N> int double=both
```

For example, to run the 2D_A case:

```bash
abaqus job=2D_plastic user=pfczm_dynamic_vuel_plastic.for cpus=4 int double=both
```

> **⚠️ Important:** The `double=both` flag is **required**. Using double precision in both the packaging and solving stages is critical for the numerical stability and accuracy of the explicit phase-field spall fracture simulation.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Author

- **HeliangYou** (https://github.com/AlexanderJFDR)

---

## Acknowledgments

- The phase-field formulation follows the framework established by Zhang et al. (2021) in the *Journal of Applied Physics*.
- This project leverages the Abaqus VUEL interface for user-defined element implementation.
