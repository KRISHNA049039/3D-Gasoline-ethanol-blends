# Ricardo Engine Combustion Analysis Suite

> **A full C++17 simulation framework for 3D meshing, physics, chemistry, combustion modelling, and emissions analysis of ethanol–gasoline blends on the Ricardo E6 research engine.**

---

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Quick Start](#quick-start)
- [Module Descriptions](#module-descriptions)
- [Physics & Models](#physics--models)
- [Results Summary](#results-summary)
- [Output Files](#output-files)
- [References](#references)

---

## Project Overview

This simulator implements a research-grade 0D/quasi-3D engine combustion analysis pipeline in modern C++17. It is built around the **Ricardo E6** variable-compression research engine (76 × 111 mm, CR 8:1, 503.5 cm³) and covers:

| Capability | Models Used |
|---|---|
| 3D Cylinder Mesh | O-grid structured hexahedral, tanh-stretched |
| Turbulence | Standard k-ε with rapid distortion |
| Fuel Chemistry | NASA 7-coefficient thermodynamics, Metghalchi-Keck / Gülder flame speed |
| Combustion Cycle | Two-zone first-law ODE, Wiebe heat release |
| Knock Prediction | Livengood-Wu integral, Chen-Flynn auto-ignition delay |
| Emissions | Westbrook-Zeldovich NOx, CO oxidation, HC from combustion efficiency |
| Wall Heat Transfer | Woschni extended model |
| Blend Sweep | E0 → E85 parametric analysis (8 blend points) |

---

## Repository Structure

```
ricardo_engine/
├── CMakeLists.txt
├── README.md
├── CONCEPTS.md               ← Physics & mathematics deep-dive
├── RESULTS.md                ← Full numerical results & analysis
├── ARCHITECTURE.md           ← Code design & module guide
│
├── include/
│   ├── utils/
│   │   └── types.hpp         ← Core types: Vec3, EngineGeometry, FuelBlend, ThermoState
│   ├── mesh/
│   │   └── cylinder_mesh.hpp ← 3D O-grid mesh generator
│   ├── physics/
│   │   └── turbulence.hpp    ← k-ε turbulence + Woschni heat transfer
│   ├── chemistry/
│   │   └── blend_chemistry.hpp ← NASA thermo, flame speed, NOx, knock
│   ├── combustion/
│   │   └── combustion_cycle.hpp ← Two-zone cycle ODE solver
│   └── analysis/
│       └── blend_analysis.hpp   ← Parametric sweep & report generation
│
├── src/
│   ├── combustion_analysis_main.cpp  ← Main entry point
│   └── main.cpp
│
└── results/
    ├── blend_summary.csv         ← All engine parameters per blend
    ├── flame_speed_map.csv       ← Su vs φ, T, blend
    ├── P_CA_E0.csv               ← Pressure trace, E0
    ├── P_CA_E5.csv
    ├── P_CA_E10.csv
    ├── P_CA_E15.csv
    ├── P_CA_E20.csv
    ├── P_CA_E30.csv
    ├── P_CA_E50.csv
    └── P_CA_E85.csv
```

---

## Quick Start

### Requirements

- GCC ≥ 9 or Clang ≥ 10 (C++17)
- CMake ≥ 3.16 (optional — direct g++ also works)
- OpenMP (optional, for parallel sweeps)

### Build with CMake

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
./combustion_analysis
```

### Build directly with g++

```bash
g++ -std=c++17 -O3 -march=native -ffast-math \
    -I include \
    src/combustion_analysis_main.cpp \
    -o combustion_analysis -lm
./combustion_analysis
```

### Expected runtime

```
Sweep (8 blends × full cycle):  ~5 ms
Mesh generation (5,000 cells):  ~1 ms per crank angle
Total wall time:                < 50 ms
```

---

## Module Descriptions

### `utils/types.hpp` — Core Data Structures (220 lines)

Defines all fundamental types used across the project:

- `Constants` — Physical/mathematical constants (R, π, σ, etc.)
- `Vec3` — 3D vector with dot, cross, norm, unit operations
- `EngineGeometry` — Ricardo E6 geometry with derived kinematics (`cylinder_volume`, `dVdtheta`, `piston_position`)
- `FuelBlend` — Blend properties via linear interpolation (LHV, AFR, RON, density, C/H/O atom counts)
- `OperatingConditions` — RPM, lambda, spark timing, intake conditions
- `ThermoState` — Full thermodynamic + species mass fraction state

### `mesh/cylinder_mesh.hpp` — 3D Mesh (341 lines)

Generates a structured hexahedral O-grid mesh of the instantaneous cylinder volume:

- Radial tanh-stretching for boundary layer resolution
- 5,720 nodes / 5,000 cells at default resolution (nr=10, nθ=20, nz=25)
- Face connectivity with boundary zone classification (fluid, liner, piston, head)
- Mesh updates dynamically as crank angle changes

### `physics/turbulence.hpp` — k-ε Model (203 lines)

Single-zone turbulence model tracking TKE evolution through the engine cycle:

- Standard k-ε with Cμ=0.09, C₁ε=1.44, C₂ε=1.92
- Rapid distortion source term during compression
- Combustion-induced TKE boost from thermal expansion
- Woschni (1967) extended heat transfer coefficient

### `chemistry/blend_chemistry.hpp` — Fuel Chemistry (289 lines)

Complete chemical kinetics for ethanol–gasoline blends:

- NASA 7-coefficient thermodynamics for O₂, N₂, CO₂, H₂O, CO, NO, iC₈H₁₈, C₂H₅OH
- Metghalchi-Keck / Gülder laminar flame speed with T/P power-law corrections
- Zimont turbulent flame speed correlation
- Westbrook-Zeldovich thermal NOx (equilibrium O-atom)
- CO oxidation (Westbrook & Dryer)
- Livengood-Wu knock integral with Chen-Flynn auto-ignition delay
- Equilibrium burned gas composition via atom balance

### `combustion/combustion_cycle.hpp` — Two-Zone ODE (228 lines)

Solves the compression-combustion-expansion cycle step by step:

- Wiebe function for heat release rate (calibrated to flame speed)
- First-law pressure update: `dP = (γ-1)/V·dQ - γ·P/V·dV`
- Zone temperatures from ideal gas law
- Woschni heat loss applied at each step
- Simultaneous Livengood-Wu knock accumulation
- Per-step NOx and CO mass accumulation

### `analysis/blend_analysis.hpp` — Parametric Sweep (207 lines)

Orchestrates the full blend sweep and report generation:

- Sweeps 8 ethanol fractions (E0 → E85)
- MBT spark adjustment per blend
- Euro 6 emissions classification
- Text report with performance / flame speed / knock / emissions tables
- CSV export for all traces

---

## Physics & Models

See **[CONCEPTS.md](CONCEPTS.md)** for full mathematical derivations of every model.

---

## Results Summary

See **[RESULTS.md](RESULTS.md)** for full tables, trends, and engineering interpretation.

---

## Output Files

| File | Content | Rows |
|---|---|---|
| `blend_summary.csv` | 22 engine parameters × 8 blends | 8 |
| `flame_speed_map.csv` | Su, δ_L × 5 blends × 7 φ × 5 T | 175 |
| `P_CA_E*.csv` | P, T_u, T_b, xb, ROHR, NOx_ppm, knock vs CA | 2,161 each |

---

## References

1. Heywood, J.B. — *Internal Combustion Engine Fundamentals*, McGraw-Hill, 1988
2. Metghalchi & Keck — *Burning velocities of mixtures of air with methanol, isooctane, and indolene*, Combustion and Flame, 1982
3. Gülder, Ö.L. — *Correlations of laminar combustion data for alternative S.I. engine fuels*, SAE 841000, 1984
4. Woschni, G. — *A universally applicable equation for the instantaneous heat transfer coefficient in the internal combustion engine*, SAE 670931, 1967
5. Livengood & Wu — *Correlation of autoignition phenomena in internal combustion engines*, 5th Int. Symp. on Combustion, 1955
6. Peters, N. — *Turbulent Combustion*, Cambridge University Press, 2000
7. Westbrook & Dryer — *Simplified reaction mechanisms for the oxidation of hydrocarbon fuels*, Combustion Science and Technology, 1981
8. Turns, S.R. — *An Introduction to Combustion*, McGraw-Hill, 3rd ed., 2012
9. Curran et al. — *A comprehensive modeling study of iso-octane oxidation*, Combustion and Flame, 2002
10. Marinov, N.M. — *A detailed chemical kinetic model for high temperature ethanol oxidation*, Int. J. Chem. Kinet., 1999
