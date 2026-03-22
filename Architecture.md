# Code Architecture Guide

> Design decisions, data flow, extension points, and API reference for the Ricardo Engine Simulation framework.

---

## Table of Contents

1. [Design Philosophy](#1-design-philosophy)
2. [Namespace Structure](#2-namespace-structure)
3. [Data Flow](#3-data-flow)
4. [Module API Reference](#4-module-api-reference)
5. [Key Algorithms — Pseudocode](#5-key-algorithms--pseudocode)
6. [Extension Guide](#6-extension-guide)
7. [Build System](#7-build-system)
8. [Numerical Stability Notes](#8-numerical-stability-notes)

---

## 1. Design Philosophy

The framework follows these principles:

**Header-only library** — all physics is in `.hpp` files under `include/`. No separate compilation of physics code required. This allows the compiler to inline small methods and enables easy modification.

**Zero external dependencies** — only the C++17 standard library and `<cmath>`. No Boost, Eigen, or third-party numerics. This makes the code portable to embedded/HPC environments.

**Value semantics** — structs are passed by value or const-reference. No heap allocation in the hot loop. The cycle solver allocates `std::vector` once at construction and appends during the sweep.

**Separation of concerns:**

```
types.hpp        ← data structures only, no computation
cylinder_mesh.hpp ← geometry only, no physics
turbulence.hpp   ← turbulence + heat transfer, no combustion
blend_chemistry.hpp ← chemistry only, no cycle integration
combustion_cycle.hpp ← ODE solver, calls chemistry
blend_analysis.hpp   ← orchestration + reporting
```

---

## 2. Namespace Structure

```cpp
namespace Ricardo {

    namespace Constants { ... }     // Physical constants

    struct Vec3 { ... }             // 3D vector
    struct EngineGeometry { ... }   // Engine dimensions + kinematics
    struct FuelBlend { ... }        // Fuel blend properties
    struct OperatingConditions { ... }
    struct ThermoState { ... }

    namespace Mesh {
        struct Node { ... }
        struct HexCell { ... }
        struct Face { ... }
        struct Zone { ... }
        class  CylinderMeshGenerator { ... }
        struct MeshStats { ... }
    }

    namespace Physics {
        struct KEpsConstants { ... }
        struct TurbState { ... }
        class  EngineTurbulence { ... }
        class  WoschniHeatTransfer { ... }
        struct SprayModel { ... }
    }

    namespace Chemistry {
        struct NASACoeffs { ... }
        class  SpeciesDB { ... }
        struct ArrheniusReaction { ... }
        class  BlendChemistry { ... }
    }

    namespace Combustion {
        class  WiebeFunction { ... }
        class  CombustionCycle { ... }
    }

    namespace Analysis {
        struct BlendSweepPoint { ... }
        class  BlendSweep { ... }
        struct EmissionsReport { ... }
        class  ReportGenerator { ... }
    }
}
```

---

## 3. Data Flow

```
EngineGeometry ──┐
FuelBlend ───────┼──→ CombustionCycle::run() ──→ CycleResult
OperatingConds ──┘         │
                           ├── WiebeFunction       (heat release shape)
                           ├── BlendChemistry      (Su, St, NOx, knock)
                           └── WoschniHeatTransfer (heat loss)

BlendSweep::run()
    └── for each ethanol fraction:
            └── CombustionCycle::run() ──→ BlendSweepPoint

BlendSweepPoint[]
    └── ReportGenerator::text_report()  ──→ stdout
    └── ReportGenerator::pressure_trace_csv() ──→ file
```

### Crank Angle Loop (inside `CombustionCycle::run`)

```
θ from -360° to +180° in steps of dθ=0.25°:
    1. V = EngineGeometry::cylinder_volume(θ)
    2. xb_new = WiebeFunction::xb(θ)
    3. dQ = dxb × Q_fuel × η_comb
    4. dP = (γ-1)/V·dQ - γ·P/V·dV   [first law]
    5. Apply Woschni heat loss
    6. Update T_u (isentropic or ideal-gas)
    7. Update T_b (adiabatic mixing)
    8. Accumulate work: W += P·dV
    9. Livengood-Wu increment
    10. NOx increment (if T_b > 1600 K)
    11. Record trace point
```

---

## 4. Module API Reference

### `EngineGeometry`

```cpp
// Instantaneous cylinder volume [m³]
double cylinder_volume(double theta_deg) const;

// Volume derivative dV/dθ [m³/rad]
double dVdtheta(double theta_deg) const;

// Piston displacement from TDC [m]
double piston_position(double theta_rad) const;

// Clearance, BDC, TDC volumes [m³]
double clearance_volume() const;
double bdc_volume() const;
double tdc_volume() const;

// Displacement [m³]
double displacement() const;

// Intake valve lift profile [m]
double intake_valve_lift_at(double cam_angle_deg) const;
```

### `FuelBlend`

```cpp
double lhv() const;              // Lower heating value [J/kg]
double stoich_afr() const;       // Stoichiometric AFR
double molar_mass() const;       // [kg/mol]
double density_liquid() const;   // [kg/m³]
double octane_number() const;    // RON
double carbon_atoms() const;     // atoms per molecule
double hydrogen_atoms() const;
double oxygen_atoms() const;
std::string name() const;        // "E0", "E10", "E85", ...
```

### `BlendChemistry`

```cpp
// Laminar flame speed [m/s] — Metghalchi-Keck/Gülder
double laminar_flame_speed(double T_u, double P, double phi) const;

// Laminar flame thickness [m]
double laminar_flame_thickness(double T_u, double P, double phi) const;

// Turbulent flame speed [m/s] — Zimont correlation
double turbulent_flame_speed(double Su, double u_prime,
                              double l_turb, double delta_L) const;

// Thermal NOx rate [mol/(m³·s)] — Westbrook-Zeldovich
double NOx_formation_rate(double T, double P, double phi) const;

// CO net rate [mol/(m³·s)]
double CO_rate(double T, double P, double phi, double Y_CO) const;

// Livengood-Wu increment [dimensionless]
double knock_integral_increment(double T, double P, double phi,
                                 double dt) const;

// Equilibrium burned gas species mass fractions
struct BurnedComposition { double Y_CO2, Y_H2O, Y_N2, Y_O2, Y_CO, Y_NO; };
BurnedComposition burned_composition(double phi) const;

// Heat release rate [W/m³]
double heat_release_rate(double T, double P, double phi, double xb) const;
```

### `CombustionCycle`

```cpp
// Constructor
CombustionCycle(const EngineGeometry& geom,
                const FuelBlend&      blend,
                const OperatingConditions& ops,
                const Config& cfg);

// Run one complete cycle, return results
CycleResult run();
```

**`CycleResult` fields:**

```cpp
struct CycleResult {
    double IMEP_bar;             // Indicated mean effective pressure
    double W_indicated_J;        // Net indicated work
    double thermal_efficiency;   // W/Q_fuel
    double comb_efficiency;      // Q_released/Q_fuel
    double BSFC_g_kWh;          // Brake-specific fuel consumption
    double peak_pressure_bar;
    double peak_pressure_CA;     // Crank angle of peak P [deg]
    double peak_temperature_K;   // Max of T_u and T_b
    double MFB10, MFB50, MFB90; // Mass fraction burned landmarks [deg]
    double CA10_90;              // Burn duration MFB90-MFB10 [deg]
    double spark_CA;
    double Su_lam_m_s;           // Laminar flame speed at spark conditions
    double St_turb_m_s;          // Turbulent flame speed
    double delta_L_um;           // Flame thickness [μm]
    double knock_index;          // Max Livengood-Wu integral
    bool   knock_detected;       // true if integral ≥ 1.0
    double NOx_g_kWh;
    double CO_g_kWh;
    double HC_g_kWh;
    // Crank angle traces:
    std::vector<double> CA_deg, P_bar, T_u_K, T_b_K;
    std::vector<double> xb_vec, ROHR_J_deg, NOx_ppm, knock_I;
};
```

### `CylinderMeshGenerator`

```cpp
// Generate mesh at given crank angle
void generate(double crank_angle_deg);

// Accessors
const std::vector<Node>&    nodes()  const;
const std::vector<HexCell>& cells()  const;
const std::vector<Face>&    faces()  const;
const std::vector<Zone>&    zones()  const;
int total_cells() const;
int total_nodes() const;
```

### `EngineTurbulence`

```cpp
// Initialise TKE from intake flow
void initialize(double k0 = -1.0);

// Advance one crank-angle step
void advance(double ca_prev, double ca_now, double density);

// TKE boost from combustion heat release
void combustion_boost(double dxb, double T_b, double T_u);

// Current state
const TurbState& state() const;

// Flame wrinkling factor Ξ
double flame_wrinkling_factor() const;
```

---

## 5. Key Algorithms — Pseudocode

### First-Law Pressure Integrator

```
for each θ in [-360°, +180°]:
    V     = cylinder_volume(θ)
    dV    = V - V_prev

    xb    = wiebe.xb(θ)
    dxb   = xb - xb_prev
    dQ    = dxb * m_fuel * LHV * eta_comb     // [J]

    gamma = mass_weighted_gamma(m_u, m_b, gamma_u, gamma_b)

    dP    = (gamma - 1) / V * dQ              // heat release
          - gamma * P / V * dV               // compression/expansion

    P     = max(P + dP, P_floor)             // floor at 500 Pa

    // Woschni correction
    h     = woschni_htc(theta, P, T_u)
    A     = wall_area(theta)
    dQ_loss = h * A * (T_u - T_wall) * dt
    P    -= (gamma - 1) / V * dQ_loss

    W    += 0.5 * (P + P_prev) * dV          // trapezoidal integration
```

### Livengood-Wu Accumulation

```
knock_integral = 0
for each θ after spark timing:
    if m_unburned > threshold:
        tau = A * P^n * exp(Ea/T_u) * f_octane / phi
        knock_integral += dtheta / (deg_per_second * tau)
        if knock_integral >= 1.0:
            KNOCK_DETECTED = true
            break
```

### Burned Zone Temperature

```
when dxb > 0:
    h_comb = LHV / (1 + AFR)               // [J/kg mixture]
    T_adiabatic = T_u + dm_b * h_comb / (m_b * Cv_b)
    T_b = (T_b_old * (m_b - dm_b) + dm_b * T_adiabatic) / m_b
    T_b = clamp(T_b, T_u, 3000 K)

    V_b = m_b * R_b * T_b / P             // burned volume
    V_u = V - V_b                          // unburned volume
    T_u = P * V_u / (m_u * R_u)           // update unburned T
```

---

## 6. Extension Guide

### Adding a new fuel species

1. Add a `static NASACoeffs species_name()` method to `SpeciesDB` in `blend_chemistry.hpp`
2. Add the species to the `FuelBlend` struct with its properties
3. Update `BlendChemistry::laminar_flame_speed()` with a new Su₀ correlation

### Adding a new emission species (e.g. N₂O)

1. Add a `double N2O_formation_rate(double T, double P, double phi)` to `BlendChemistry`
2. Add accumulator variables in `CombustionCycle::run()` following the NOx pattern
3. Add the field to `CycleResult`
4. Include in the CSV export in `blend_analysis.hpp`

### Adding a new combustion model (e.g. fractal flame)

1. Create `include/combustion/fractal_flame.hpp` with a class implementing:
   ```cpp
   double xb(double ca) const;
   double dxb_dtheta(double ca) const;
   ```
2. Replace `WiebeFunction wiebe` in `CombustionCycle::run()` with your model

### Extending the mesh (e.g. intake port)

1. Add port geometry parameters to `EngineGeometry`
2. In `CylinderMeshGenerator`, add an `generate_intake_port()` private method
3. Append port nodes/cells to the main arrays in `generate()`
4. Add a new `Zone` with type="inlet"

### Running a lambda sweep

```cpp
EngineGeometry geom; // ... configure
FuelBlend blend;     // e.g. E10
blend.ethanol_fraction = 0.10;
blend.gasoline_fraction = 0.90;

for (double lambda : {0.8, 0.9, 1.0, 1.1, 1.2, 1.3}) {
    OperatingConditions ops;
    ops.lambda = lambda;
    ops.rpm = 2000;
    Combustion::CombustionCycle::Config cfg;
    auto result = Combustion::CombustionCycle(geom, blend, ops, cfg).run();
    std::cout << lambda << "  " << result.IMEP_bar << "  " << result.NOx_g_kWh << "\n";
}
```

### Running a compression ratio sweep

```cpp
for (double CR : {8.0, 9.0, 10.0, 11.0, 12.0}) {
    EngineGeometry g;
    g.compression_ratio = CR;
    FuelBlend b; b.ethanol_fraction = 0.85; b.gasoline_fraction = 0.15;
    OperatingConditions ops; ops.rpm = 2000; ops.lambda = 1.0;
    Combustion::CombustionCycle::Config cfg;
    auto r = Combustion::CombustionCycle(g, b, ops, cfg).run();
    std::cout << "CR=" << CR
              << "  IMEP=" << r.IMEP_bar
              << "  knock=" << r.knock_index
              << (r.knock_detected ? " KNOCK" : " OK") << "\n";
}
```

---

## 7. Build System

### CMakeLists.txt targets

| Target | Description |
|---|---|
| `ricardo_engine_lib` | Static library of all physics modules |
| `ricardo_sim` | Minimal entry point (calls combustion_analysis) |
| `combustion_analysis` | Full analysis suite with all modules |

### Compiler flags

| Flag | Purpose |
|---|---|
| `-O3` | Full optimisation |
| `-march=native` | Use host CPU instruction set (AVX2 etc.) |
| `-ffast-math` | Allow floating-point reassociation (safe here) |
| `-fopenmp` | Parallelise blend sweep (optional) |

### Platform notes

The code is portable C++17. No platform-specific headers. Tested on:
- Ubuntu 22.04 / GCC 11, 13
- macOS 14 / Clang 15
- Windows 11 / MSVC 19.38 (with minor pragma adjustments)

---

## 8. Numerical Stability Notes

### Pressure floor

```cpp
P = std::max(500.0, P + dP);   // 0.005 bar minimum
```

Prevents division-by-zero in temperature calculations during the early exhaust stroke.

### Temperature clamping

```cpp
T_u = std::clamp(T_u, 200.0, 2500.0);
T_b = std::clamp(T_b, T_u, 3000.0);
```

Prevents runaway in the rare case of non-physical flame temperature (e.g. very rich misfire).

### Volume partition guard

```cpp
V_b = std::clamp(V_b, 0.0, V * 0.999);
V_u = V - V_b;
```

Ensures V_u > 0 even when xb → 1.0, preventing negative unburned volume.

### Knock integral guard

```cpp
knock_I_val += chem_.knock_integral_increment(T_u, P, phi, dt);
```

The auto-ignition delay `τ` is floored at `1e-9 s` inside `knock_integral_increment()` to prevent division by zero at very high temperatures.

### Step size

The default `dθ = 0.25°` corresponds to:
- 0.208 ms at 2000 rpm
- Captures the entire 540°-range (BDC→TDC→BDC) in 2,161 steps

Reducing to `dθ = 0.1°` gives ~2× more accurate work integration at 5× compute cost. For production sweeps `dθ = 0.5°` is adequate.
