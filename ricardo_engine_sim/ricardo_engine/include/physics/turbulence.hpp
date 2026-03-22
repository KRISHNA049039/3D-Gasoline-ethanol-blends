#pragma once
#include "utils/types.hpp"
#include <cmath>
#include <vector>

namespace Ricardo {
namespace Physics {

// ─────────────────────────────────────────────
//  Standard k-ε turbulence model constants
// ─────────────────────────────────────────────
struct KEpsConstants {
    double Cmu   = 0.09;
    double C1eps = 1.44;
    double C2eps = 1.92;
    double sigma_k   = 1.0;
    double sigma_eps = 1.3;
};

// ─────────────────────────────────────────────
//  Turbulent kinetic energy state
// ─────────────────────────────────────────────
struct TurbState {
    double k   = 0.0;  // TKE [m²/s²]
    double eps = 0.0;  // dissipation [m²/s³]
    double mu_t = 0.0; // turbulent viscosity [Pa·s]
    double l_int = 0.0;// integral length scale [m]
    double u_prime = 0.0; // rms velocity fluctuation [m/s]
};

// ─────────────────────────────────────────────
//  Single-zone k-ε model for engine cylinder
//  Tracks TKE generation during compression/combustion
// ─────────────────────────────────────────────
class EngineTurbulence {
public:
    explicit EngineTurbulence(const EngineGeometry& geom,
                               const OperatingConditions& ops,
                               const KEpsConstants& consts = {})
        : geom_(geom), ops_(ops), C_(consts) {}

    void initialize(double k0 = -1.0) {
        // Initial TKE from intake flow
        double u_bar = ops_.rpm / 60.0 * geom_.stroke;  // mean piston speed [m/s]
        state_.k     = (k0 > 0) ? k0 : 1.5 * std::pow(0.05 * u_bar, 2.0);
        // Initial length scale = fraction of bore
        state_.l_int = 0.1 * geom_.bore;
        state_.eps   = std::pow(C_.Cmu, 0.75) *
                       std::pow(state_.k, 1.5) / state_.l_int;
        update_derived();
    }

    // Advance model one crank-angle step
    void advance(double ca_prev, double ca_now, double density) {
        double dCA = ca_now - ca_prev;
        if (std::abs(dCA) < 1e-9) return;
        double dt = dCA / ops_.crank_deg_per_second();

        double V     = geom_.cylinder_volume(ca_now);
        double V_prev = geom_.cylinder_volume(ca_prev);
        double dV_dt = (V - V_prev) / (dt + 1e-12);

        // Strain rate from mean piston motion
        double S_ij  = -dV_dt / (3.0 * V + 1e-12); // volumetric rate

        // Production term: P_k = 2 * mu_t * Sij * Sij
        double P_k   = 2.0 * state_.mu_t * std::pow(S_ij, 2.0);

        // Compression work on TKE (rapid distortion)
        double P_rapid = -2.0 / 3.0 * state_.k * dV_dt / (V + 1e-12);

        // k-equation
        double dk_dt = P_k + P_rapid - state_.eps;
        state_.k     = std::max(1e-6, state_.k + dk_dt * dt);

        // ε-equation
        double deps_dt = (C_.C1eps * state_.eps / (state_.k + 1e-9)) *
                         (P_k + P_rapid) -
                         C_.C2eps * state_.eps * state_.eps / (state_.k + 1e-9);
        state_.eps   = std::max(1e-8, state_.eps + deps_dt * dt);

        update_derived();
        (void)density;
    }

    // Combustion-induced TKE boost (from thermal expansion)
    void combustion_boost(double dxb, double T_b, double T_u) {
        double gamma_boost = (T_b / (T_u + 1e-9) - 1.0) * dxb;
        state_.k   *= (1.0 + 0.5 * gamma_boost);
        state_.eps *= (1.0 + 0.3 * gamma_boost);
        update_derived();
    }

    const TurbState& state() const { return state_; }

    // Turbulent Prandtl flame thickness enhancement factor
    double flame_wrinkling_factor() const {
        // Damköhler: ratio of integral to flame scales
        return 1.0 + std::sqrt(state_.u_prime / (state_.u_prime + 0.1));
    }

private:
    EngineGeometry      geom_;
    OperatingConditions ops_;
    KEpsConstants       C_;
    TurbState           state_;

    void update_derived() {
        state_.mu_t   = C_.Cmu * state_.k * state_.k / (state_.eps + 1e-9) * 1.2;
        state_.u_prime = std::sqrt(2.0 / 3.0 * state_.k);
        state_.l_int  = std::pow(C_.Cmu, 0.75) *
                         std::pow(state_.k, 1.5) / (state_.eps + 1e-9);
    }
};

// ─────────────────────────────────────────────
//  Wall heat transfer (Woschni extended model)
// ─────────────────────────────────────────────
class WoschniHeatTransfer {
public:
    struct Config {
        double C1 = 2.28;   // Woschni constants
        double C2 = 0.00324;
        Config() = default;
    };

    WoschniHeatTransfer(const EngineGeometry& g, const OperatingConditions& o) : geom_(g), ops_(o), cfg_(Config()) {}
    WoschniHeatTransfer(const EngineGeometry& g, const OperatingConditions& o,
                         const Config& cfg)
        : geom_(g), ops_(o), cfg_(cfg) {}

    // Heat transfer coefficient [W/(m²·K)]
    double htc(double ca, double P_Pa, double T_K,
               double P_motored_Pa = 0.0) const {
        double bore = geom_.bore;
        double w_mean = 2.0 * geom_.stroke * ops_.rpm / 60.0; // mean piston speed

        double V   = geom_.cylinder_volume(ca);
        double V_d = geom_.displacement();
        double T_ref = ops_.intake_temp;
        double P_ref = ops_.intake_pressure;

        double c1 = cfg_.C1;
        double c2 = (P_motored_Pa > 0)
            ? cfg_.C2 * V_d * T_ref / P_ref * std::abs(P_Pa - P_motored_Pa)
            : 0.0;

        double w = c1 * w_mean + c2;
        // Woschni formula: h = 3.26 * bore^-0.2 * P^0.8 * T^-0.55 * w^0.8
        double P_kPa = P_Pa / 1000.0;
        return 3.26 * std::pow(bore, -0.2) *
               std::pow(P_kPa, 0.8) *
               std::pow(T_K, -0.55) *
               std::pow(w, 0.8);
    }

    // Total heat loss [W]
    double Q_loss(double ca, double P_Pa, double T_K,
                  double T_wall_K = 453.0) const {
        double h   = htc(ca, P_Pa, T_K);
        double V   = geom_.cylinder_volume(ca);
        double bore = geom_.bore;
        double h_piston = V / (Constants::PI / 4.0 * bore * bore);
        // Surface area (simplified cylinder)
        double A = 2.0 * Constants::PI / 4.0 * bore * bore +
                   Constants::PI * bore * h_piston;
        return h * A * (T_K - T_wall_K);
    }

private:
    EngineGeometry      geom_;
    OperatingConditions ops_;
    Config              cfg_;
};

// ─────────────────────────────────────────────
//  Spray / mixture formation model (simplified)
// ─────────────────────────────────────────────
struct SprayModel {
    // Sauter mean diameter [m] from empirical correlation
    static double SMD(double injection_pressure_MPa,
                      double density_liquid,
                      double viscosity_liquid,
                      double surface_tension) {
        // Modified Hiroyasu-Arai correlation
        double dp  = injection_pressure_MPa * 1e6;
        double rho = density_liquid;
        double We  = rho * dp / surface_tension;
        double Re  = std::sqrt(2.0 * dp / rho) * 1e-3 / viscosity_liquid;
        return 23.9 * std::pow(We, -0.261) * std::pow(Re, -0.33) * 1e-3;
    }

    // Evaporation time [s] for a droplet (D² law)
    static double evap_time(double D0, double T_gas, double T_boil,
                             double lambda_g = 0.026) {
        double delta_T = std::max(1.0, T_gas - T_boil);
        double Kv = 8.0 * lambda_g / (800.0 * 2.0 * 250000.0 / delta_T);
        return D0 * D0 / (Kv + 1e-15);
    }
};

} // namespace Physics
} // namespace Ricardo
