#pragma once
#include "utils/types.hpp"
#include <cmath>
#include <array>
#include <string>
#include <vector>

namespace Ricardo {
namespace Chemistry {

// ─────────────────────────────────────────────
//  Species thermodynamic data (NASA 7-coefficient)
//  Cp/R = a1 + a2*T + a3*T² + a4*T³ + a5*T⁴
//  H/RT = a1 + a2/2*T + ... + a6/T
//  S/R  = a1*ln(T) + a2*T + ...
// ─────────────────────────────────────────────
struct NASACoeffs {
    std::string name;
    double molar_mass;             // kg/mol
    std::array<double,7> lo;       // 200–1000 K
    std::array<double,7> hi;       // 1000–5000 K
    double T_switch = 1000.0;

    double Cp_R(double T) const {
        const auto& a = (T < T_switch) ? lo : hi;
        return a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*a[4])));
    }
    double H_RT(double T) const {
        const auto& a = (T < T_switch) ? lo : hi;
        return a[0] + T*(a[1]/2 + T*(a[2]/3 + T*(a[3]/4 + T*a[4]/5))) + a[5]/T;
    }
    double S_R(double T) const {
        const auto& a = (T < T_switch) ? lo : hi;
        return a[0]*std::log(T) + T*(a[1] + T*(a[2]/2 + T*(a[3]/3 + T*a[4]/4))) + a[6];
    }
    double Cp(double T) const { return Cp_R(T) * Constants::R_GAS / molar_mass; }
    double H(double T)  const { return H_RT(T) * Constants::R_GAS * T / molar_mass; }
    double G(double T)  const { return (H_RT(T) - S_R(T)) * Constants::R_GAS * T; }
};

// ─────────────────────────────────────────────
//  Species database
// ─────────────────────────────────────────────
class SpeciesDB {
public:
    static NASACoeffs O2() {
        NASACoeffs s; s.name = "O2"; s.molar_mass = 0.031999;
        s.lo = {3.78245636,-2.99673416e-3,9.84730201e-6,-9.68129509e-9,3.24372837e-12,-1063.94356,3.65767573};
        s.hi = {3.61219356,7.48531660e-4,-1.82255654e-7,2.27977781e-11,-1.13556012e-15,-1200.97885,3.67768058};
        return s;
    }
    static NASACoeffs N2() {
        NASACoeffs s; s.name = "N2"; s.molar_mass = 0.028014;
        s.lo = {3.53100528,-1.23660988e-4,-5.02999433e-7,2.43530612e-9,-1.40881235e-12,-1046.97628,2.96747038};
        s.hi = {2.95257637,1.39690040e-3,-4.92631603e-7,7.86010195e-11,-4.60755204e-15,-923.948688,5.87188762};
        return s;
    }
    static NASACoeffs CO2() {
        NASACoeffs s; s.name = "CO2"; s.molar_mass = 0.044010;
        s.lo = {2.35681300,8.98412990e-3,-7.12206320e-6,2.45730080e-9,-1.42885480e-13,-48371.9697,9.90090350};
        s.hi = {3.85746029,4.41437026e-3,-2.21481404e-6,5.23490188e-10,-4.72084164e-14,-48759.166,2.27163806};
        return s;
    }
    static NASACoeffs H2O() {
        NASACoeffs s; s.name = "H2O"; s.molar_mass = 0.018015;
        s.lo = {4.19864056,-2.03643410e-3,6.52040211e-6,-5.48797062e-9,1.77197250e-12,-30293.7267,-0.849032208};
        s.hi = {3.03399249,2.17691804e-3,-1.64072518e-7,-9.70419870e-11,1.68200992e-14,-30004.2971,4.96675409};
        return s;
    }
    static NASACoeffs CO() {
        NASACoeffs s; s.name = "CO"; s.molar_mass = 0.028010;
        s.lo = {3.57953347,-6.10353680e-4,1.01681433e-6,9.07005884e-10,-9.04424499e-13,-14344.086,3.50840928};
        s.hi = {2.71518561,2.06252743e-3,-9.98825771e-7,2.30053008e-10,-2.03647716e-14,-14151.8724,7.81868772};
        return s;
    }
    static NASACoeffs NO() {
        NASACoeffs s; s.name = "NO"; s.molar_mass = 0.030006;
        s.lo = {4.21859896,-4.63988124e-3,1.10443049e-5,-9.34055507e-9,2.80554874e-12,9845.09964,2.28061001};
        s.hi = {3.26071234,1.19101135e-3,-4.29122646e-7,6.94481463e-11,-4.03295681e-15,9921.43132,6.36900469};
        return s;
    }
    // Surrogate for gasoline (iso-octane C8H18)
    static NASACoeffs isoOctane() {
        NASACoeffs s; s.name = "IC8H18"; s.molar_mass = 0.11423;
        s.lo = {-4.20868,0.0780376,-4.87548e-5,1.71876e-8,-2.65836e-12,-2.75730e4,3.49674e1};
        s.hi = {1.75479e1,0.0473780,-1.63727e-5,2.54337e-9,-1.49300e-13,-3.20805e4,-7.04530e1};
        return s;
    }
    // Ethanol C2H5OH
    static NASACoeffs ethanol() {
        NASACoeffs s; s.name = "C2H5OH"; s.molar_mass = 0.046069;
        s.lo = {4.85887168,-3.74017174e-3,6.95553926e-5,-8.86547318e-8,3.51688350e-11,-2.99966764e4,4.80186978};
        s.hi = {6.56243650,1.52406230e-2,-5.38385050e-6,8.62536440e-10,-5.12120880e-14,-3.15256880e4,-9.46871510};
        return s;
    }
};

// ─────────────────────────────────────────────
//  Reaction rate (Arrhenius)
// ─────────────────────────────────────────────
struct ArrheniusReaction {
    std::string name;
    double A;       // pre-exponential factor
    double n;       // temperature exponent
    double Ea_R;    // activation energy / R_gas  [K]
    double rate(double T) const {
        return A * std::pow(T, n) * std::exp(-Ea_R / T);
    }
};

// ─────────────────────────────────────────────
//  Reduced chemistry: 12-step ethanol-gasoline blend
//  Based on Curran et al. + Marinov et al. mechanisms
// ─────────────────────────────────────────────
class BlendChemistry {
public:
    explicit BlendChemistry(const FuelBlend& blend) : blend_(blend) {}

    // ── Global heat release rate  [W/m³] ─────────
    // Wiebe-based approach combined with Arrhenius for knocking
    double heat_release_rate(double T, double P, double phi, double xb) const {
        // Unburned zone temperature dependent rate
        double rho = P / (287.0 * T);   // approximate
        double lhv = blend_.lhv();
        // Modified Arrhenius for blend
        double A_glob = blend_.ethanol_fraction * 6.5e9 + blend_.gasoline_fraction * 5.0e9;
        double Ea_R   = blend_.ethanol_fraction * 14500.0 + blend_.gasoline_fraction * 15100.0;
        double k      = A_glob * std::exp(-Ea_R / T);
        double Yf     = std::max(0.0, (1.0 - xb) * 0.06); // approx fuel mass fraction
        double YO2    = 0.233 * std::max(0.0, 1.0 - phi * 0.06 / blend_.stoich_afr());
        return rho * k * Yf * std::pow(YO2, 1.5) * lhv;
    }

    // ── Laminar flame speed [m/s] ─────────────────
    //  Metghalchi-Keck + ethanol correction (Gülder model)
    double laminar_flame_speed(double T_u, double P, double phi) const {
        // Reference conditions
        double T_ref = 298.0, P_ref = 101325.0;
        double T_rat = T_u / T_ref;
        double P_rat = P  / P_ref;

        // Base flame speed at stoichiometry for each component
        double Su0_eth = 0.448 * std::exp(-2.29 * std::pow(phi - 1.08, 2)); // m/s
        double Su0_gas = 0.311 * std::exp(-2.29 * std::pow(phi - 1.10, 2)); // m/s

        // Blend
        double Su0 = blend_.ethanol_fraction * Su0_eth + blend_.gasoline_fraction * Su0_gas;

        // Temperature & pressure corrections
        double alpha = 2.18 - 0.80 * (phi - 1.0);
        double beta  = -0.16 + 0.22 * (phi - 1.0);
        double Su    = Su0 * std::pow(T_rat, alpha) * std::pow(P_rat, beta);

        // EGR dilution correction (Baratta et al.)
        // already included via phi for simplicity
        return std::max(0.01, Su);
    }

    // ── Turbulent flame speed [m/s] ──────────────
    //  Peters correlation
    double turbulent_flame_speed(double Su, double u_prime, double l_turb,
                                  double delta_L) const {
        // Karlowitz number based
        double Da  = Su * l_turb / (u_prime * delta_L + 1e-9); // Damköhler
        double Ka  = std::pow(u_prime / Su, 1.5) *
                     std::pow(delta_L / l_turb, 0.5);
        (void)Da; (void)Ka;
        // Peters (1999) correlation
        double ratio = u_prime / Su;
        // Zimont correlation (more stable for SI engines)
        double St = Su * (1.0 + 0.62 * std::pow(ratio, 0.41));
        return std::max(Su, std::min(St, 5.0 * Su));  // cap at 5x laminar
    }

    // ── Flame thickness (laminar) [m] ────────────
    double laminar_flame_thickness(double T_u, double P, double phi) const {
        double Su   = laminar_flame_speed(T_u, P, phi);
        // Thermal diffusivity of mixture
        double lambda_mix = 0.026 + 2e-4 * T_u;  // W/(m·K) approx
        double rho   = P / (287.0 * T_u);
        double Cp    = 1050.0;
        double alpha_th = lambda_mix / (rho * Cp);
        return alpha_th / std::max(Su, 1e-6);
    }

    // ── NOx formation (Zeldovich mechanism) ──────
    // Returns d[NO]/dt in mol/(m³·s)
    double NOx_formation_rate(double T, double P, double phi) const {
        if (T < 1800.0) return 0.0;  // Thermal NOx threshold
        double rho  = P / (287.0 * T);
        // O₂ and N₂ concentrations [mol/m³]
        double XO2  = (phi < 1.0) ? 0.21 * (1.0 - phi) / (1.0 + 0.21*(phi-1.0)+0.001) : 0.001;
        double CO2m = rho * 1000.0 / 29.0;  // approximate total molar conc [mol/m³]
        double CO2  = CO2m * XO2;
        double CN2  = CO2m * 0.79;

        // Westbrook/Turns Zeldovich thermal NOx (SI units)
        // rNO [mol/(m3*s)] = 6e10/sqrt(T)*exp(-69090/T)*[N2]*sqrt([O2])
        //  where concentrations in mol/m3
        // This gives physically correct ppm-scale NOx for SI engines
        double rNO = 6.0e10 / std::sqrt(T) * std::exp(-69090.0 / T)
                     * CN2 * std::sqrt(std::max(CO2, 1e-4));
        return std::max(0.0, rNO);
    }

    // ── CO formation/oxidation rate [mol/(m³·s)] ─
    double CO_rate(double T, double P, double phi, double Y_CO) const {
        double rho = P / (287.0 * T);
        double CCO = rho * 1000.0 * Y_CO / 0.028;   // mol/m³
        double CO2_conc = rho * 1000.0 * 0.15 / 0.044;

        // CO oxidation: CO + OH -> CO2 + H  (Westbrook & Dryer)
        double k_CO = 6.324e7 * std::pow(T, 1.5) * std::exp(-15098.0 / T);
        double COH  = rho * 1000.0 / 29.0 * 0.001;  // rough OH estimate
        double rate_ox  = k_CO * CCO * COH;

        // Rich-condition CO formation from partial oxidation
        double rate_form = (phi > 1.0) ? 1.5e6 * std::exp(-8000.0 / T) *
                            rho * 0.01 : 0.0;
        (void)CO2_conc;
        return rate_form - rate_ox;
    }

    // ── Knock integral (Livengood-Wu) ──────────
    // Returns knock integral; engine knocks when integral >= 1
    double knock_integral_increment(double T, double P, double phi, double dt) const {
        // Auto-ignition delay time (Chen correlation for blends)
        double A_knock  = blend_.ethanol_fraction * 17.15 + blend_.gasoline_fraction * 11.2;
        double n_knock  = blend_.ethanol_fraction * (-3.3) + blend_.gasoline_fraction * (-4.1);
        double Ea_knock = blend_.ethanol_fraction * 12700.0 + blend_.gasoline_fraction * 13400.0;
        double ON_factor = (100.0 - blend_.octane_number()) / 100.0 * 0.4 + 0.6;

        double tau = A_knock * std::pow(P / Constants::ATM_PA, n_knock) *
                     std::exp(Ea_knock / T) * ON_factor / phi;
        return dt / std::max(tau, 1e-9);
    }

    // ── Equilibrium burned gas composition ───────
    struct BurnedComposition {
        double Y_CO2, Y_H2O, Y_N2, Y_O2, Y_CO, Y_NO;
    };
    BurnedComposition burned_composition(double phi) const {
        BurnedComposition bc{};
        // Simple atom-balance for CₐHᵦOᵧ fuel + air
        double a = blend_.carbon_atoms();
        double b = blend_.hydrogen_atoms();
        double g = blend_.oxygen_atoms();
        double stoich_O = a + b/4.0 - g/2.0; // O atoms needed stoich

        if (phi <= 1.0) {
            // Lean/stoich – complete combustion
            double f = phi;
            // Per mole of fuel
            double n_CO2 = a;
            double n_H2O = b / 2.0;
            double n_O2_excess = stoich_O * (1.0 - f) / f;
            double n_N2  = stoich_O / f * (79.0 / 21.0);
            double total  = n_CO2 + n_H2O + n_O2_excess + n_N2;
            bc.Y_CO2 = n_CO2 * 0.044 / total;
            bc.Y_H2O = n_H2O * 0.018 / total;
            bc.Y_O2  = n_O2_excess * 0.032 / total;
            bc.Y_N2  = n_N2 * 0.028 / total;
            bc.Y_CO  = 0.0;
            bc.Y_NO  = 0.001 * phi; // rough
        } else {
            // Rich – partial combustion
            double n_CO2 = a * 1.0/phi;
            double n_H2O = b / 2.0 / phi;
            double n_CO  = a * (1.0 - 1.0/phi);
            double n_H2  = (b / 2.0) * (1.0 - 1.0/phi);
            double n_N2  = stoich_O * (79.0 / 21.0);
            double total  = n_CO2 + n_H2O + n_CO + n_H2 + n_N2;
            bc.Y_CO2 = n_CO2 * 0.044 / total;
            bc.Y_H2O = n_H2O * 0.018 / total;
            bc.Y_CO  = n_CO  * 0.028 / total;
            bc.Y_N2  = n_N2  * 0.028 / total;
            bc.Y_O2  = 0.0;
            bc.Y_NO  = 0.0;
        }
        return bc;
    }

private:
    FuelBlend blend_;
    double Schmidt_ = 0.7; // turbulent Schmidt number
};

} // namespace Chemistry
} // namespace Ricardo
