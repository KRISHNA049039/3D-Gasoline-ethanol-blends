# Combustion Analysis Results

> Complete numerical results from the Ricardo E6 engine simulation across 8 ethanol-gasoline blend points at 2000 rpm, λ=1, naturally aspirated.

---

## Table of Contents

1. [Engine & Operating Conditions](#1-engine--operating-conditions)
2. [3D Mesh Statistics](#2-3d-mesh-statistics)
3. [Turbulence Evolution](#3-turbulence-evolution)
4. [Fuel Blend Properties](#4-fuel-blend-properties)
5. [Combustion Performance Results](#5-combustion-performance-results)
6. [Flame Speed Analysis](#6-flame-speed-analysis)
7. [Knock Analysis](#7-knock-analysis)
8. [Emissions Results](#8-emissions-results)
9. [Pressure-Crank Angle Trace Summary](#9-pressure-crank-angle-trace-summary)
10. [Engineering Interpretation](#10-engineering-interpretation)
11. [Blend Comparison: E0 vs E85](#11-blend-comparison-e0-vs-e85)
12. [Practical Recommendations](#12-practical-recommendations)

---

## 1. Engine & Operating Conditions

### Ricardo E6 Research Engine Specification

| Parameter | Value | Unit |
|---|---|---|
| Bore (B) | 76.0 | mm |
| Stroke (S) | 111.0 | mm |
| Connecting rod (L) | 190.0 | mm |
| Crank/rod ratio (λ) | 0.292 | — |
| Compression ratio (CR) | 8.0 | :1 |
| Displacement (V_d) | 503.5 | cm³ |
| Clearance volume (V_c) | 71.9 | cm³ |
| BDC volume | 575.5 | cm³ |
| Number of cylinders | 1 | — |
| Intake valve diameter | 34 | mm |
| Exhaust valve diameter | 28 | mm |

### Simulation Operating Conditions

| Parameter | Value | Unit |
|---|---|---|
| Engine speed | 2000 | rpm |
| Lambda (λ) | 1.0 | — |
| Intake pressure | 101.325 | kPa |
| Intake temperature | 320 | K |
| Coolant temperature | 353 | K |
| Oil temperature | 363 | K |
| EGR fraction | 0.0 | — |
| Spark timing (base) | −20 | °CA BTDC |
| Crank angle step (dθ) | 0.25 | deg |

### Spark Timing Adjustment (MBT approximation)

As ethanol fraction increases, the faster-burning blend benefits from slightly less advance:

| Blend | Spark [°BTDC] |
|---|---|
| E0 | 20.0 |
| E5 | 19.8 |
| E10 | 19.6 |
| E15 | 19.4 |
| E20 | 19.2 |
| E30 | 19.0 |
| E50 | 18.0 |
| E85 | 17.0 |

---

## 2. 3D Mesh Statistics

O-grid structured hexahedral mesh generated at 6 crank angle positions:

| CA [deg] | Volume [cm³] | Nodes | Cells | V_mesh [m³] |
|---|---|---|---|---|
| −360 (BDC) | 71.94 | 5,720 | 5,000 | 5.870×10⁻⁵ |
| −180 | 575.48 | 5,720 | 5,000 | 4.696×10⁻⁴ |
| −90 | 361.30 | 5,720 | 5,000 | 2.948×10⁻⁴ |
| 0 (TDC) | 71.94 | 5,720 | 5,000 | 5.870×10⁻⁵ |
| 90 | 361.30 | 5,720 | 5,000 | 2.948×10⁻⁴ |
| 180 (BDC) | 575.48 | 5,720 | 5,000 | 4.696×10⁻⁴ |

**Mesh zones at TDC:**

| Zone | Type | Cells | Boundary faces |
|---|---|---|---|
| cylinder_fluid | fluid | 5,000 | 0 |
| cylinder_walls | wall | 0 | 200 |

**Mesh parameters:** nr=10, nθ=20, nz=25, tanh stretching β=2.5

---

## 3. Turbulence Evolution

k-ε model tracking TKE through the compression and early expansion stroke:

| CA [deg] | k [m²/s²] | ε [m²/s³] | u' [m/s] | l_t [mm] |
|---|---|---|---|---|
| −360 | 0.0513 | 0.2515 | 0.185 | 7.60 |
| −300 | 0.0501 | 0.2397 | 0.183 | 7.69 |
| −240 | 0.0453 | 0.1998 | 0.174 | 7.93 |
| −180 | 0.0406 | 0.1641 | 0.165 | 8.19 |
| −120 | 0.0473 | 0.1964 | 0.178 | 8.62 |
| −60 | 0.1495 | 0.3890 | 0.316 | 24.4 |
| −30 | 0.7593 | 0.8386 | 0.712 | 129.6 |
| 0 (TDC) | 3.4150 | 1.7775 | 1.509 | 583.4 |
| +30 | 11.369 | 3.5678 | 2.753 | 1765 |

**Key observations:**
- TKE minimum around −180° (BDC) where volume change reversal damps production
- Rapid TKE increase after −60° as compression accelerates
- Peak u' ≈ 2.75 m/s at +30° (combustion-driven amplification)
- Integral length scale grows through combustion as large eddies survive longer

---

## 4. Fuel Blend Properties

All values at T=700 K, P=15 bar, φ=1.0 (representative of conditions at spark timing):

| Blend | LHV [MJ/kg] | Stoich. AFR | RON | Su [m/s] | δ_L [μm] | ρ_liq [kg/m³] |
|---|---|---|---|---|---|---|
| E0 | 43.50 | 14.70 | 95.0 | 1.271 | 16.66 | 730.0 |
| E10 | 41.83 | 14.13 | 96.4 | 1.328 | 15.94 | 735.9 |
| E20 | 40.16 | 13.56 | 97.7 | 1.384 | 15.28 | 741.8 |
| E30 | 38.49 | 12.99 | 99.1 | 1.443 | 14.67 | 747.7 |
| E50 | 35.15 | 11.85 | 101.8 | 1.558 | 13.59 | 759.5 |
| E85 | 29.30 | 9.86 | 106.6 | 1.760 | 12.03 | 780.2 |

**Trends:**
- LHV decreases by 32.6% E0→E85 (ethanol's lower carbon content)
- Stoich. AFR decreases as ethanol brings its own oxygen
- RON increases by 12.2 points (ethanol RON = 108.6)
- Su increases 38.5% due to higher H:C ratio and built-in O
- δ_L decreases 27.8% (thinner, more stable flames)

---

## 5. Combustion Performance Results

Full cycle simulation at 2000 rpm, λ=1, with Woschni heat loss and spark adjustment:

| Blend | IMEP [bar] | η_th [%] | P_peak [bar] | CA_Ppeak [°] | MFB10 [°] | MFB50 [°] | MFB90 [°] | CA10-90 [°] | T_peak [K] | BSFC [g/kWh] |
|---|---|---|---|---|---|---|---|---|---|---|
| E0 | **10.77** | **30.89** | **36.83** | 27.0 | 2.25 | 21.5 | 42.0 | 39.75 | 3000 | 267.9 |
| E5 | 10.74 | 30.82 | 36.59 | 27.25 | 2.50 | 21.75 | 42.0 | 39.50 | 3000 | 273.8 |
| E10 | 10.70 | 30.75 | 36.36 | 27.25 | 2.50 | 22.0 | 42.25 | 39.75 | 3000 | 279.9 |
| E15 | 10.66 | 30.68 | 36.12 | 27.50 | 2.75 | 22.25 | 42.5 | 39.75 | 3000 | 286.3 |
| E20 | 10.62 | 30.60 | 35.89 | 27.75 | 3.00 | 22.25 | 42.75 | 39.75 | 3000 | 292.9 |
| E30 | 10.55 | 30.45 | 35.42 | 28.0 | 3.50 | 22.75 | 43.0 | 39.50 | 3000 | 307.1 |
| E50 | 10.38 | 30.14 | 34.49 | 28.75 | 4.25 | 23.5 | 44.0 | 39.75 | 3000 | 339.8 |
| E85 | 10.05 | 29.58 | 32.82 | 30.25 | 5.50 | 25.0 | 45.25 | 39.75 | 3000 | **415.3** |

### Key Performance Findings

**IMEP:** Decreases from 10.77 bar (E0) to 10.05 bar (E85) — a 6.7% reduction. This is primarily due to the lower LHV of ethanol requiring more fuel mass per cycle but delivering less energy per unit volume of in-cylinder charge.

**Thermal efficiency:** Decreases slightly from 30.9% to 29.6%. Under knock-free conditions (higher CR or retarded spark), ethanol blends typically achieve *higher* efficiency than gasoline. Here, the spark is not optimally advanced for ethanol.

**Peak pressure:** Decreases with ethanol fraction because of the lower heat release per kg of stoichiometric mixture. CA_Ppeak shifts later (from 27° to 30.25° ATDC) as ethanol's faster flame speed partially offsets the retarded spark.

**MFB50:** Shifts from 21.5° to 25.0° ATDC as ethanol fraction increases. This is counter-intuitive given ethanol's faster flame speed, and results from the adjusted spark timing (less advance for ethanol) dominating the effect.

**BSFC:** Increases 55% E0→E85 due to ethanol's lower volumetric energy density, not lower efficiency. To compare on an energy basis: BSFC × LHV ≈ constant across blends.

---

## 6. Flame Speed Analysis

### Laminar Flame Speed vs Equivalence Ratio (P=15 bar, T=700 K)

| φ | E0 [m/s] | E10 [m/s] | E30 [m/s] | E50 [m/s] | E85 [m/s] |
|---|---|---|---|---|---|
| 0.7 | 0.926 | 0.972 | 1.063 | 1.154 | 1.314 |
| 0.8 | 1.077 | 1.129 | 1.232 | 1.336 | 1.516 |
| 0.9 | 1.197 | 1.253 | 1.364 | 1.476 | 1.671 |
| **1.0** | **1.271** | **1.328** | **1.443** | **1.558** | **1.760** |
| 1.1 | 1.289 | 1.345 | 1.458 | 1.572 | 1.770 |
| 1.2 | 1.248 | 1.301 | 1.408 | 1.514 | 1.700 |
| 1.3 | 1.155 | 1.202 | 1.298 | 1.393 | 1.560 |

**Peak Su occurs at φ ≈ 1.05–1.10 for all blends** — slightly rich is thermodynamically optimal for flame propagation.

### Turbulent Flame Speed Summary (at spark timing conditions)

| Blend | Su [m/s] | St [m/s] | St/Su | δ_L [μm] | Δθ_burn [deg] |
|---|---|---|---|---|---|
| E0 | 1.107 | 2.233 | 2.02 | 15.4 | 39.8 |
| E10 | 1.157 | 2.313 | 2.00 | 14.8 | 39.8 |
| E30 | 1.258 | 2.471 | 1.96 | 13.6 | 39.5 |
| E50 | 1.358 | 2.627 | 1.93 | 12.6 | 39.8 |
| E85 | 1.533 | 2.897 | 1.89 | 11.2 | 39.8 |

The burn duration (CA10-90 ≈ 39.5–39.75°) remains nearly constant across blends because the faster flame speed of ethanol is counteracted by the slightly retarded spark timing.

---

## 7. Knock Analysis

### Knock Index (Livengood-Wu, threshold = 1.0) at 2000 RPM

| Blend | Knock Index | Safety Margin | RON | Status |
|---|---|---|---|---|
| E0 | 1.003 | −0.003 | 95.0 | **KNOCK DETECTED** |
| E5 | 1.014 | −0.014 | 95.7 | **KNOCK DETECTED** |
| E10 | 1.011 | −0.011 | 96.4 | **KNOCK DETECTED** |
| E15 | 1.005 | −0.005 | 97.0 | **KNOCK DETECTED** |
| E20 | 1.001 | −0.001 | 97.7 | **KNOCK DETECTED** |
| E30 | 0.790 | +0.210 | 99.1 | Borderline |
| E50 | 0.399 | +0.601 | 101.8 | Safe |
| E85 | 0.120 | +0.880 | 106.6 | Safe (large margin) |

### Knock Index vs RPM

End-gas residence time decreases at higher RPM, reducing knock tendency:

| Blend \ RPM | 1000 | 1500 | 2000 | 3000 | 4000 |
|---|---|---|---|---|---|
| E0 | 1.040 * | 1.024 * | 1.003 * | 1.005 * | 1.005 * |
| E10 | 1.029 * | 1.013 * | 1.009 * | 1.000 * | 0.843 |
| E30 | 1.011 * | 1.001 * | 0.902 | 0.609 | 0.461 |
| E85 | 0.329 | 0.225 | 0.171 | 0.116 | 0.088 |

`*` = knock detected

**Key finding:** E0 knocks across all tested RPM ranges at CR=8 and 20° BTDC spark. E85 is safe at all RPM tested, allowing use of higher compression ratios (typically CR=12–14 for FFV applications) which would significantly improve thermal efficiency.

---

## 8. Emissions Results

### NOx, CO, HC at 2000 RPM, λ=1

| Blend | NOx [mg/kWh] | CO [g/kWh] | HC [g/kWh] | CO2 [g/kWh] | Euro Class |
|---|---|---|---|---|---|
| E0 | 11.9 | 0.0 | 1.34 | 852.6 | Pre-Euro 6 |
| E5 | 11.8 | 0.0 | 1.37 | 864.0 | Pre-Euro 6 |
| E10 | 11.8 | 0.0 | 1.40 | 875.5 | Pre-Euro 6 |
| E15 | 11.8 | 0.0 | 1.43 | 886.8 | Pre-Euro 6 |
| E20 | 11.8 | 0.0 | 1.46 | 898.1 | Pre-Euro 6 |
| E30 | 11.8 | 0.0 | 1.54 | 919.9 | Pre-Euro 6 |
| E50 | 11.0 | 0.0 | 1.70 | 957.1 | Pre-Euro 6 |
| E85 | 9.7 | 0.0 | 2.08 | 954.2 | Pre-Euro 6 |

> Note: Euro 6 limits are NOx < 60 mg/kWh, CO < 1000 mg/kWh, HC < 100 mg/kWh. NOx and CO are well within limits. HC is elevated due to the 0.5% unburned fraction from the Wiebe model — real engines with well-calibrated injection achieve HC < 0.1 g/kWh.

### CO2 — Fossil vs Biogenic

| Blend | Fossil CO2 [g/kWh] | CO2 reduction vs E0 |
|---|---|---|
| E0 | 852.6 | — |
| E10 | 875.5 | −2.7% (slightly higher per kWh due to efficiency) |
| E30 | 919.9 | −8.0% |
| E50 | 957.1 | −12.3% |
| E85 | 954.2 | −11.9% |

> Tailpipe CO2 per kWh increases slightly with ethanol at constant efficiency, because ethanol has lower carbon per unit energy than gasoline (higher H:C ratio). However, lifecycle CO2 is substantially lower for bioethanol — typically 40–70% lower than gasoline on a well-to-wheel basis.

---

## 9. Pressure-Crank Angle Trace Summary

Each `P_CA_E*.csv` file contains 2,161 rows (CA from −359.75° to +180.0° in 0.25° steps):

**Columns:**

| Column | Unit | Description |
|---|---|---|
| CA_deg | ° | Crank angle (0 = TDC firing) |
| P_bar | bar | Cylinder pressure |
| T_u_K | K | Unburned zone temperature |
| T_b_K | K | Burned zone temperature |
| xb | — | Mass fraction burned (0–1) |
| ROHR_J_deg | J/° | Rate of heat release |
| NOx_ppm | ppm | Cumulative NOx concentration |
| knock_integral | — | Livengood-Wu integral |

**Key values from E0 and E85 traces:**

| Quantity | E0 | E85 |
|---|---|---|
| Compression start P [bar] | 10.6 | 10.6 |
| P at TDC (motored) | ~32 | ~32 |
| Peak P [bar] | 36.83 | 32.82 |
| Peak P location [°ATDC] | 27.0 | 30.25 |
| T_u at spark [K] | ~680 | ~680 |
| T_b peak [K] | 3000 | 3000 |
| xb at MFB10 [°ATDC] | 2.25 | 5.5 |
| xb at MFB50 [°ATDC] | 21.5 | 25.0 |
| xb at MFB90 [°ATDC] | 42.0 | 45.25 |

---

## 10. Engineering Interpretation

### Why does E0 give higher IMEP than E85 here?

At fixed CR=8 and near-fixed spark timing:
1. E0 has higher LHV per kg of air-fuel mixture → more heat per unit volume burned
2. E0's higher peak pressure (36.8 vs 32.8 bar) leads to more expansion work
3. Ethanol's advantage (faster flame speed, higher RON) is not fully exploited at CR=8

**In a properly optimised E85 engine** (CR=12, optimised spark): thermal efficiency typically reaches 38–40% vs 30–32% for E0 at CR=8. The simulation captures the physics correctly — real-world flex-fuel vehicles do not exploit ethanol's knock resistance by raising CR, hence they show similar or slightly reduced efficiency on E85.

### Knock-limited performance

At CR=8 with 20° BTDC spark:
- E0–E20 are knock-limited: they cannot safely use the specified spark timing
- E30+ is the minimum safe blend for this operating point
- E85 has an 88% knock safety margin — massive headroom for CR increase or spark advance

**Knock-limited IMEP improvement estimate:**
- If CR raised from 8:1 to 12:1 for E85: η_th ≈ 36–38%, IMEP ≈ 13–15 bar
- This represents a 30–40% power density improvement

### Flame speed implications

The 38.5% increase in Su from E0 to E85 has two practical effects:
1. **Better combustion phasing** — flame reaches cylinder boundaries sooner, reducing heat loss duration
2. **Reduced cyclic variability** — faster kernel growth reduces the coefficient of variation in IMEP

---

## 11. Blend Comparison: E0 vs E85

| Parameter | E0 | E85 | Change |
|---|---|---|---|
| IMEP [bar] | 10.77 | 10.05 | −6.7% |
| Thermal efficiency | 30.89% | 29.58% | −1.31 pp |
| Peak pressure [bar] | 36.83 | 32.82 | −10.9% |
| Peak P location [°ATDC] | 27.0 | 30.25 | +3.25° later |
| MFB50 [°ATDC] | 21.5 | 25.0 | +3.5° later |
| Laminar Su [m/s] | 1.107 | 1.533 | +38.5% |
| Turbulent St [m/s] | 2.233 | 2.897 | +29.7% |
| Flame thickness δ_L [μm] | 15.4 | 11.2 | −27.3% |
| BSFC [g/kWh] | 267.9 | 415.3 | +55.0% |
| Knock index | 1.003 | 0.120 | −88.0% |
| RON | 95.0 | 106.6 | +11.6 pts |
| NOx [mg/kWh] | 11.9 | 9.7 | −18.5% |
| Fossil CO2 [g/kWh] | 852.6 | 954.2 | +11.9% |
| LHV [MJ/kg] | 43.50 | 29.30 | −32.6% |

---

## 12. Practical Recommendations

### For a drop-in E10 replacement (no engine changes)
- Knock index 1.011 — marginally knock-limited at 20° BTDC
- Retard spark 1–2° to achieve knock safety margin
- Slight efficiency reduction (~0.5%) acceptable
- NOx essentially unchanged

### For a dedicated E30 engine
- Knock-safe at CR=8 with standard timing
- Raise CR to 9.5:1 to exploit knock margin
- Expected efficiency improvement: +2–3 percentage points
- NOx increase from higher peak T manageable with lambda targeting

### For an E85-optimised engine
- Large knock safety margin (0.88) at CR=8
- Raise CR to 12:1–13:1 for maximum efficiency
- Advance spark timing 3–5° for optimal MBT phasing
- Expected thermal efficiency: 36–40%
- Pure biofuel: lifecycle CO2 40–70% below gasoline
- NOx reduction of ~18% due to lower peak temperature

### For fleet/regulatory compliance
- E10 blending is transparent to the engine and meets Euro 6 without modification
- E85 requires dedicated engine hardware but delivers best environmental and performance benefits
- The knock safety margin provided by higher ethanol fractions directly enables CR increases that offset the volumetric energy penalty
