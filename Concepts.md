# Physics & Mathematical Concepts

> Complete derivation of every physical model implemented in the Ricardo Engine Combustion Analysis Suite.

---

## Table of Contents

1. [Engine Kinematics](#1-engine-kinematics)
2. [3D Mesh Generation](#2-3d-mesh-generation)
3. [k-ε Turbulence Model](#3-k-ε-turbulence-model)
4. [Fuel Blend Thermochemistry](#4-fuel-blend-thermochemistry)
5. [Laminar Flame Speed](#5-laminar-flame-speed)
6. [Turbulent Flame Speed](#6-turbulent-flame-speed)
7. [Wiebe Heat Release Function](#7-wiebe-heat-release-function)
8. [Two-Zone First-Law Combustion Model](#8-two-zone-first-law-combustion-model)
9. [Wall Heat Transfer — Woschni Model](#9-wall-heat-transfer--woschni-model)
10. [Knock Prediction — Livengood-Wu Integral](#10-knock-prediction--livengood-wu-integral)
11. [Thermal NOx — Extended Zeldovich Mechanism](#11-thermal-nox--extended-zeldovich-mechanism)
12. [CO Formation and Oxidation](#12-co-formation-and-oxidation)
13. [Equilibrium Burned Gas Composition](#13-equilibrium-burned-gas-composition)
14. [Performance Metrics](#14-performance-metrics)
15. [Sign Conventions & Coordinate System](#15-sign-conventions--coordinate-system)

---

## 1. Engine Kinematics

### Piston Position

The instantaneous piston displacement from TDC is computed using the exact slider-crank formula (no approximation):

```
s(θ) = a(1 - cosθ) + l(1 - √(1 - λ²sin²θ))
```

where:
- `a = stroke/2` — crank radius [m]
- `l` — connecting-rod length [m]
- `λ = a/l` — crank-to-rod ratio (≈ 0.292 for Ricardo E6)
- `θ` — crank angle from TDC [rad]

### Instantaneous Cylinder Volume

```
V(θ) = V_c + A_p · s(θ)
```

where:
- `V_c = V_d / (CR - 1)` — clearance volume
- `A_p = π/4 · B²` — piston crown area
- `B` — bore diameter
- `CR` — compression ratio

### Volume Derivative

```
dV/dθ = A_p · a · sinθ · [1 + λcosθ / √(1 - λ²sin²θ)]
```

Used in the first-law pressure equation. Accuracy is critical near TDC where `sinθ → 0` and the term under the radical approaches 1.

### Ricardo E6 Geometry

| Parameter | Value |
|---|---|
| Bore B | 76.0 mm |
| Stroke S | 111.0 mm |
| Con-rod L | 190.0 mm |
| Crank-rod ratio λ | 0.2921 |
| Compression ratio CR | 8.0 |
| Displacement V_d | 503.5 cm³ |
| Clearance volume V_c | 71.9 cm³ |
| BDC volume | 575.5 cm³ |

---

## 2. 3D Mesh Generation

### O-Grid Topology

The cylinder is discretised using a structured **hexahedral O-grid**. This avoids the singularity at the cylinder axis that plagues polar grids.

```
Radial layers  : nr = 10
Circumferential: nθ = 20
Axial          : nz = 25
Total nodes    : (nr+1) × nθ × (nz+1) = 5,720
Total cells    : nr × nθ × nz         = 5,000
```

### Radial Stretching

Nodes are clustered toward the cylinder wall using hyperbolic tangent stretching:

```
r_i = R · tanh(β · i/nr) / tanh(β),   β = 2.5
```

This concentrates ~40% of nodes within the inner 20% of the radius, resolving the thermal and velocity boundary layers near the liner.

### Node Indexing

A flat 3D index maps `(i, j, k)` → node_id:

```
node_id = k · nθ · (nr+1)  +  j · (nr+1)  +  i
```

where `i` = radial, `j` = circumferential (periodic), `k` = axial.

### Boundary Zone Assignment

| Zone | Condition | BC Type |
|---|---|---|
| Cylinder fluid | interior | fluid |
| Liner wall | i = nr | wall (no-slip) |
| Piston crown | k = 0 | wall (moving) |
| Cylinder head | k = nz | wall |
| Intake port | user-specified face set | inlet |
| Exhaust port | user-specified face set | outlet |

### Cell Volume

Each hexahedral cell volume is approximated by the scalar triple product of three edge vectors from one corner:

```
V_cell = |a⃗ · (b⃗ × c⃗)|
```

where `a⃗`, `b⃗`, `c⃗` are edge vectors from node 0 to nodes 1, 3, 4.

---

## 3. k-ε Turbulence Model

### Transport Equations

The standard two-equation k-ε model solves:

**TKE equation:**
```
Dk/Dt = P_k + P_rapid - ε
```

**Dissipation equation:**
```
Dε/Dt = (C₁ε · ε/k) · (P_k + P_rapid) - C₂ε · ε²/k
```

### Model Constants

| Constant | Value | Physical meaning |
|---|---|---|
| C_μ | 0.09 | Turbulent viscosity coefficient |
| C₁ε | 1.44 | Production-to-dissipation ratio |
| C₂ε | 1.92 | Dissipation decay rate |
| σ_k | 1.0 | TKE Prandtl number |
| σ_ε | 1.3 | Dissipation Prandtl number |

### Production Terms

Mean-flow shear production:
```
P_k = 2 μ_t · S_ij · S_ij
```

Rapid distortion (compression/expansion):
```
P_rapid = -(2/3) k · (1/V) · dV/dt
```

During compression `dV/dt < 0`, so `P_rapid > 0` — TKE is generated. This is the dominant source during the compression stroke.

### Turbulent Viscosity

```
μ_t = C_μ · ρ · k² / ε
```

### Derived Quantities

| Quantity | Formula |
|---|---|
| RMS velocity fluctuation | `u' = √(2k/3)` |
| Integral length scale | `l_t = C_μ^(3/4) · k^(3/2) / ε` |
| Turbulent time scale | `τ_t = k/ε` |
| Kolmogorov length scale | `η = (ν³/ε)^(1/4)` |

### Initialisation

At intake-valve closing (BDC, θ = -360°):

```
k₀ = 1.5 · (I · ū)²      (I = turbulence intensity ≈ 5%)
ε₀ = C_μ^(3/4) · k₀^(3/2) / l₀      (l₀ ≈ 0.1 · B)
```

---

## 4. Fuel Blend Thermochemistry

### Linear Blending

All blend properties are computed by linear volume-fraction interpolation between neat ethanol and neat gasoline (iso-octane surrogate):

```
X_blend = x_eth · X_eth + x_gas · X_gas
```

### Species Properties — NASA 7-Coefficient Polynomials

Heat capacity, enthalpy, and entropy as functions of temperature:

```
Cp/R = a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴

H/RT = a₁ + a₂T/2 + a₃T²/3 + a₄T³/4 + a₅T⁴/5 + a₆/T

S/R  = a₁lnT + a₂T + a₃T²/2 + a₄T³/3 + a₅T⁴/4 + a₇
```

Two sets of coefficients (low: 200–1000 K, high: 1000–5000 K) are used with a switch at T_switch = 1000 K.

Species implemented: **O₂, N₂, CO₂, H₂O, CO, NO, iC₈H₁₈ (gasoline surrogate), C₂H₅OH (ethanol)**

### Fuel Blend Key Properties

| Property | E0 | E10 | E30 | E50 | E85 |
|---|---|---|---|---|---|
| LHV [MJ/kg] | 43.50 | 41.83 | 38.49 | 35.15 | 29.30 |
| Stoich. AFR | 14.70 | 14.13 | 12.99 | 11.85 | 9.86 |
| RON | 95.0 | 96.4 | 99.1 | 101.8 | 106.6 |
| Density [kg/m³] | 730 | 736 | 748 | 760 | 780 |
| LHV × density [MJ/L] | 31.8 | 30.8 | 28.8 | 26.7 | 22.9 |

---

## 5. Laminar Flame Speed

### Metghalchi-Keck / Gülder Correlation

The laminar burning velocity is modelled as:

```
Su(T_u, P, φ) = Su₀(φ) · (T_u/T_ref)^α · (P/P_ref)^β
```

with:
- `T_ref = 298 K`, `P_ref = 101.325 kPa`

### Reference Flame Speed Su₀(φ)

**Gasoline (iso-octane):**
```
Su₀ = 0.311 · exp(-2.29 · (φ - 1.10)²)   [m/s]
```

**Ethanol:**
```
Su₀ = 0.448 · exp(-2.29 · (φ - 1.08)²)   [m/s]
```

**Blend:**
```
Su₀_blend = x_eth · Su₀_eth + x_gas · Su₀_gas
```

### Exponents (φ-dependent)

```
α = 2.18 - 0.80·(φ - 1)
β = -0.16 + 0.22·(φ - 1)
```

### Physical Interpretation

At stoichiometric conditions and 700 K / 15 bar (representative of compression conditions at spark timing):

| Blend | Su [m/s] | Increase vs E0 |
|---|---|---|
| E0 | 1.271 | — |
| E10 | 1.328 | +4.5% |
| E30 | 1.443 | +13.5% |
| E50 | 1.558 | +22.6% |
| E85 | 1.760 | +38.5% |

Ethanol's higher Su comes from its lower C:H ratio and built-in oxygen atom, which shortens the chemical induction length.

### Laminar Flame Thickness

```
δ_L = α_th / Su
```

where `α_th = λ/(ρ·Cp)` is the thermal diffusivity of the unburned mixture.

| Blend | δ_L [μm] at 700K / 15 bar / φ=1 |
|---|---|
| E0 | 15.4 |
| E30 | 13.6 |
| E85 | 11.2 |

Thinner flames with ethanol reduce stretch sensitivity and quenching probability.

---

## 6. Turbulent Flame Speed

### Zimont Correlation

```
St = Su · [1 + 0.62 · (u'/Su)^0.41]
```

capped at `5 × Su` to avoid physically unrealistic values in the high-turbulence TDC region.

### Physical Basis

The ratio `u'/Su` is the **turbulent Damköhler number complement** — it quantifies how much the turbulence wrinkles the flame front. When `u' >> Su`, the flame is in the **corrugated flamelets regime** and wrinkling dominates.

### Typical Values at Spark Timing

| Quantity | Value |
|---|---|
| Mean piston speed | 3.7 m/s |
| u' (estimated) | ~0.5 × ū_piston = 1.85 m/s |
| l_t | 6–8 mm |
| Su (E0) | 1.11 m/s |
| St (E0) | 2.23 m/s |
| Su (E85) | 1.53 m/s |
| St (E85) | 2.90 m/s |

---

## 7. Wiebe Heat Release Function

### Mass Fraction Burned

```
x_b(θ) = 1 - exp[-a · ξ^(m+1)]
```

where `ξ = (θ - θ₀) / Δθ` is the normalised crank angle, with:
- `θ₀` — start of combustion (spark timing)
- `Δθ` — combustion duration
- `a = 5.0` — efficiency parameter (≈ 5 for SI engines, ensures 99.3% completion)
- `m = 2.0` — form factor (controls shape: m=0 is exponential, m=2 gives bell-shaped ROHR)

### Rate of Heat Release

```
dxb/dθ = a·(m+1)/Δθ · ξ^m · (1 - xb)
```

### Combustion Duration Calibration

Duration `Δθ` is calibrated from the mean flame speed:

```
Δθ [deg] = (flame_travel_distance / S_avg) · (RPM × 6)
```

where:
- `flame_travel_distance ≈ 0.75 · B` (from spark to farthest wall)
- `S_avg = (Su + St) / 2`

This couples the Wiebe parameterisation to the actual chemistry and turbulence.

### MFB Landmarks

| Landmark | x_b |
|---|---|
| MFB10 (end of flame kernel) | 0.10 |
| MFB50 (centre of combustion) | 0.50 |
| MFB90 (end of main burn) | 0.90 |
| CA10-90 (burn duration) | MFB90 − MFB10 |

---

## 8. Two-Zone First-Law Combustion Model

### Zone Definition

At any crank angle after spark timing, the cylinder contains two zones:

```
Unburned zone (u): T_u, m_u, V_u — compressed unburned mixture
Burned zone   (b): T_b, m_b, V_b — combustion products
```

with the constraints:
```
m_u + m_b = m_total = const.
V_u + V_b = V(θ)
P_u = P_b = P  (pressure equilibrium)
```

### First-Law Pressure Update

Applying the first law to the entire cylinder:

```
dP/dθ = [(γ_eff - 1)/V] · dQ/dθ  -  [γ_eff · P / V] · dV/dθ
```

where `γ_eff` is the mass-weighted effective ratio of specific heats:

```
γ_eff = (m_u·Cv_u·γ_u + m_b·Cv_b·γ_b) / (m_u·Cv_u + m_b·Cv_b)
```

| Zone | Cv [J/kg·K] | γ |
|---|---|---|
| Unburned | 740 | 1.35 |
| Burned | 850 | 1.28 |

### Zone Temperature Updates

**Unburned zone** — isentropic compression/expansion:
```
T_u,new = T_u,old · (V_old/V_new)^(γ_u - 1)
```

**Burned zone** — from energy balance when mass `dm_b` burns:
```
T_ad = T_u + dm_b · h_comb / (m_b · Cv_b)
T_b,new = [T_b,old · (m_b - dm_b) + dm_b · T_ad] / m_b
```

where `h_comb = LHV / (1 + AFR)` is the heat release per kg of mixture.

### Volume Partition

After temperature update:
```
V_b = m_b · R_b · T_b / P
V_u = V - V_b
T_u = P · V_u / (m_u · R_u)
```

### Indicated Work

```
W = ∫ P dV  ≈  Σ 0.5·(P_i + P_{i-1})·(V_i - V_{i-1})
```

Positive W = net expansion work (energy delivered to piston).

---

## 9. Wall Heat Transfer — Woschni Model

### Heat Transfer Coefficient

```
h = 3.26 · B^{-0.2} · P^{0.8} · T^{-0.55} · w^{0.8}   [W/m²·K]
```

where:
- `B` — bore [m]
- `P` — pressure [kPa]
- `T` — gas temperature [K]
- `w` — characteristic gas velocity [m/s]:

```
w = C₁ · ū_piston + C₂ · (V_d · T_ref / P_ref) · |P - P_motored|
```

with `C₁ = 2.28` (gas exchange), `C₂ = 0.00324` (combustion enhancement).

### Heat Loss Rate

```
Q̇_loss = h · A_wall · (T_gas - T_wall)
```

**Wall temperatures used:**

| Surface | Temperature |
|---|---|
| Liner | T_coolant = 353 K |
| Piston crown | T_oil + 30 K = 393 K |
| Cylinder head | T_coolant + 20 K = 373 K |

### Effect on Cycle

The Woschni correction typically reduces peak pressure by 3–8% and thermal efficiency by 2–5 percentage points compared to an adiabatic simulation.

---

## 10. Knock Prediction — Livengood-Wu Integral

### Knock Criterion

Engine knock occurs when the end-gas auto-ignites before the flame front reaches it. The **Livengood-Wu integral** tracks the fractional progress toward auto-ignition:

```
I = ∫ dt / τ(T_u, P, φ)
```

Knock occurs when `I ≥ 1.0`.

### Auto-Ignition Delay Correlation (Chen-Flynn)

```
τ = A · P^n · exp(E_a / T_u) · f(ON) / φ
```

Blend-weighted parameters:

| Parameter | Gasoline | Ethanol |
|---|---|---|
| A | 11.2 | 17.15 |
| n (pressure exp.) | −4.1 | −3.3 |
| E_a/R [K] | 13,400 | 12,700 |

Octane number correction factor:
```
f(ON) = (100 - ON)/100 · 0.4 + 0.6
```

Higher-octane fuels (ethanol blends) have a longer `τ` → lower `I` → better knock resistance.

### Knock Results at 2000 RPM, 20° BTDC, λ=1

| Blend | RON | Knock Index | Status |
|---|---|---|---|
| E0 | 95.0 | 1.003 | **KNOCK** |
| E10 | 96.4 | 1.011 | **KNOCK** |
| E20 | 97.7 | 1.001 | **KNOCK** |
| E30 | 99.1 | 0.790 | Borderline |
| E50 | 101.8 | 0.399 | Safe |
| E85 | 106.6 | 0.120 | Safe (large margin) |

At CR=8 with 20° BTDC spark, E0–E20 are knock-limited at 2000 rpm. Increasing ethanol content to E30+ provides knock safety margin, allowing higher compression ratios or more advanced spark timing for better efficiency.

---

## 11. Thermal NOx — Extended Zeldovich Mechanism

### Reaction Scheme

The dominant thermal NOx pathway at high temperatures (T > 1800 K):

```
R1:  O  + N₂  →  NO + N      (rate-limiting)
R2:  N  + O₂  →  NO + O
R3:  N  + OH  →  NO + H
```

### Rate Expression

Using the Westbrook/Turns formulation with equilibrium O-atom concentration:

```
d[NO]/dt = 6×10¹⁰ / √T · exp(-69090/T) · [N₂] · √[O₂]   [mol/(m³·s)]
```

where concentrations are in mol/m³.

### Equilibrium O-atom

```
[O] = K_eq(T) · √[O₂]
```

with `K_eq` from thermodynamic equilibrium of `½O₂ ⇌ O`.

### Key Feature

The exponential temperature dependence (`exp(-69090/T)`) means NOx formation is extremely sensitive to peak temperature:

- At 2400 K: rate is ~100× higher than at 2000 K
- Hence, retarding spark timing (lowering T_peak) is the primary NOx reduction strategy
- Ethanol's lower adiabatic flame temperature (due to lower LHV) reduces NOx vs gasoline

---

## 12. CO Formation and Oxidation

### Reaction

```
CO + OH  →  CO₂ + H
```

### Rate (Westbrook & Dryer)

```
k_CO = 6.324×10⁷ · T^1.5 · exp(-15098/T)   [m³/(mol·s)]

r_CO = k_CO · [CO] · [OH]   [mol/(m³·s)]
```

### Physical Behaviour

- CO forms during rich combustion (φ > 1) from partial oxidation of the fuel
- CO is subsequently oxidised in the expansion stroke as temperature decreases
- Lean operation (λ = 1 in this simulation) minimises CO
- Ethanol's built-in oxygen atom (−OH group) helps oxidise CO during flame passage

---

## 13. Equilibrium Burned Gas Composition

### Atom Balance for Lean/Stoichiometric (φ ≤ 1)

For a fuel `C_a H_b O_γ` burned with λ × stoichiometric air:

```
CₐHᵦOᵧ + (a + b/4 − γ/2)/λ · (O₂ + 3.76N₂)  →
    a·CO₂ + b/2·H₂O + (a + b/4 − γ/2)·(1−λ)/λ · O₂ + 3.76·(a + b/4 − γ/2)/λ · N₂
```

### For Rich Combustion (φ > 1)

CO and H₂ appear in the products:

```
x·CO₂ + (1−x)·CO + y·H₂O + (1−y)/2·H₂ + N₂
```

where `x` and `y` are solved from O and H atom balances.

---

## 14. Performance Metrics

### Indicated Mean Effective Pressure (IMEP)

```
IMEP = W_net / V_d   [Pa]
```

### Thermal Efficiency

```
η_th = W_net / Q_fuel
```

where `Q_fuel = m_fuel · LHV`.

### Brake-Specific Fuel Consumption

```
BSFC = m_fuel / W_net   [g/kWh]
```

Note: BSFC increases with ethanol fraction despite similar efficiency, because ethanol's lower LHV requires more fuel mass to deliver the same energy.

### Combustion Efficiency

```
η_comb = Q_released / Q_fuel
```

For this simulation, set at 99.5% (Wiebe xb → 1 with 0.5% unburned HC).

### Emissions — g/kWh Conversion

```
NOx [g/kWh] = m_NO_per_cycle [kg] × 1000 / (W_net [J] / 3,600,000)
```

---

## 15. Sign Conventions & Coordinate System

### Crank Angle

```
θ = 0°     : TDC (top dead centre, firing)
θ = -360°  : BDC (start of compression, IVC)
θ = +180°  : BDC (end of expansion, EVO)
Negative θ : before TDC (BTDC)
Positive θ : after TDC (ATDC)
```

### 3D Coordinate System

```
Z-axis : cylinder axis (positive upward toward head)
X, Y   : radial plane
Origin : cylinder axis at piston crown at TDC
```

### Thermodynamic Sign Convention

```
dQ > 0 : heat added to working fluid (combustion)
dW > 0 : work done by working fluid (expansion)
First law: dU = dQ - dW
```
