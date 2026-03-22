#pragma once
#include <cmath>
#include <array>
#include <vector>
#include <string>

namespace Ricardo {

// ─────────────────────────────────────────────
//  Physical & mathematical constants
// ─────────────────────────────────────────────
namespace Constants {
    constexpr double PI          = 3.14159265358979323846;
    constexpr double R_UNIVERSAL = 8314.46261815324;   // J/(kmol·K)
    constexpr double R_GAS       = 8.314462618;        // J/(mol·K)
    constexpr double AVOGADRO    = 6.02214076e23;
    constexpr double BOLTZMANN   = 1.380649e-23;       // J/K
    constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;
    constexpr double ATM_PA      = 101325.0;           // Pa
    constexpr double T_STP       = 298.15;             // K
    constexpr double P_STP       = 101325.0;           // Pa

    // Engine-specific
    constexpr double DEG_TO_RAD  = PI / 180.0;
    constexpr double RAD_TO_DEG  = 180.0 / PI;
    constexpr double RPM_TO_RPS  = 1.0 / 60.0;
}

// ─────────────────────────────────────────────
//  3D Vector
// ─────────────────────────────────────────────
struct Vec3 {
    double x = 0, y = 0, z = 0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s)       const { return {x*s,   y*s,   z*s};   }
    Vec3 operator/(double s)       const { return {x/s,   y/s,   z/s};   }
    Vec3& operator+=(const Vec3& o) { x+=o.x; y+=o.y; z+=o.z; return *this; }
    Vec3& operator*=(double s)      { x*=s;   y*=s;   z*=s;   return *this; }

    double dot(const Vec3& o)  const { return x*o.x + y*o.y + z*o.z; }
    Vec3   cross(const Vec3& o) const {
        return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
    }
    double norm()  const { return std::sqrt(dot(*this)); }
    double norm2() const { return dot(*this); }
    Vec3   unit()  const { double n = norm(); return (n > 1e-15) ? *this/n : Vec3{}; }
};

inline Vec3 operator*(double s, const Vec3& v) { return v * s; }

// ─────────────────────────────────────────────
//  Engine geometry specification (Ricardo E6 style)
// ─────────────────────────────────────────────
struct EngineGeometry {
    // Cylinder dimensions
    double bore            = 0.0760;   // m  (76 mm)
    double stroke          = 0.1110;   // m (111 mm)
    double con_rod_length  = 0.1900;   // m (190 mm)
    double compression_ratio = 8.0;
    int    num_cylinders   = 1;

    // Valve geometry
    double intake_valve_dia   = 0.034; // m
    double exhaust_valve_dia  = 0.028; // m
    double intake_valve_lift  = 0.010; // m (max)
    double exhaust_valve_lift = 0.009; // m (max)

    // Port geometry
    double intake_port_dia    = 0.032; // m
    double exhaust_port_dia   = 0.026; // m

    // Combustion chamber (hemispherical)
    double squish_height      = 0.0010; // m
    double bowl_depth         = 0.008;  // m

    // Derived quantities
    double displacement() const {
        return Constants::PI / 4.0 * bore * bore * stroke;
    }
    double clearance_volume() const {
        return displacement() / (compression_ratio - 1.0);
    }
    double tdc_volume() const { return clearance_volume(); }
    double bdc_volume() const { return clearance_volume() + displacement(); }

    // Piston position from crank angle (radians)
    double piston_position(double theta) const {
        double a = stroke / 2.0;
        double l = con_rod_length;
        double lambda = a / l;
        return a * (1.0 - std::cos(theta)) +
               l * (1.0 - std::sqrt(1.0 - lambda*lambda * std::sin(theta)*std::sin(theta)));
    }

    // Instantaneous cylinder volume (m³)
    double cylinder_volume(double theta_deg) const {
        double theta = theta_deg * Constants::DEG_TO_RAD;
        double s = piston_position(theta);
        double piston_area = Constants::PI / 4.0 * bore * bore;
        return clearance_volume() + piston_area * s;
    }

    // dV/dtheta (m³/rad)
    double dVdtheta(double theta_deg) const {
        double theta = theta_deg * Constants::DEG_TO_RAD;
        double a  = stroke / 2.0;
        double l  = con_rod_length;
        double lm = a / l;
        double sinT = std::sin(theta);
        double cosT = std::cos(theta);
        double piston_area = Constants::PI / 4.0 * bore * bore;
        double ds_dtheta = a * sinT * (1.0 + lm * cosT /
            std::sqrt(1.0 - lm*lm * sinT*sinT));
        return piston_area * ds_dtheta;
    }

    // Intake valve lift profile (simplified cosine)
    double intake_valve_lift_at(double cam_angle_deg) const {
        double phi = cam_angle_deg * Constants::DEG_TO_RAD;
        double open  = 0.0 * Constants::DEG_TO_RAD;
        double close = 240.0 * Constants::DEG_TO_RAD;
        if (phi < open || phi > close) return 0.0;
        double x = (phi - open) / (close - open);
        return intake_valve_lift * std::sin(Constants::PI * x);
    }
};

// ─────────────────────────────────────────────
//  Fuel blend specification
// ─────────────────────────────────────────────
struct FuelBlend {
    double ethanol_fraction = 0.0;    // 0–1 volume fraction
    double gasoline_fraction = 1.0;

    // Fuel properties (mass basis)
    double lhv() const {               // Lower heating value J/kg
        return ethanol_fraction * 26.8e6 + gasoline_fraction * 43.5e6;
    }
    double stoich_afr() const {        // Stoichiometric air–fuel ratio
        return ethanol_fraction * 9.0 + gasoline_fraction * 14.7;
    }
    double molar_mass() const {        // kg/mol
        return ethanol_fraction * 0.04607 + gasoline_fraction * 0.1142;
    }
    double density_liquid() const {    // kg/m³
        return ethanol_fraction * 789.0 + gasoline_fraction * 730.0;
    }
    double octane_number() const {     // RON
        return ethanol_fraction * 108.6 + gasoline_fraction * 95.0;
    }
    double reid_vapor_pressure() const { // kPa
        return ethanol_fraction * 16.0 + gasoline_fraction * 62.0;
    }

    // Carbon/hydrogen/oxygen atoms (for a representative molecule)
    double carbon_atoms()   const { return ethanol_fraction * 2.0 + gasoline_fraction * 8.26; }
    double hydrogen_atoms() const { return ethanol_fraction * 6.0 + gasoline_fraction * 15.5; }
    double oxygen_atoms()   const { return ethanol_fraction * 1.0 + gasoline_fraction * 0.0;  }

    std::string name() const {
        int pct = static_cast<int>(ethanol_fraction * 100.0 + 0.5);
        return "E" + std::to_string(pct);
    }
};

// ─────────────────────────────────────────────
//  Operating conditions
// ─────────────────────────────────────────────
struct OperatingConditions {
    double rpm              = 2000.0;  // rev/min
    double lambda           = 1.0;     // relative air-fuel ratio (1 = stoich)
    double spark_timing_deg = -20.0;   // crank deg BTDC (negative = before TDC)
    double intake_pressure  = 101325.0; // Pa
    double intake_temp      = 300.0;   // K
    double coolant_temp     = 353.0;   // K
    double oil_temp         = 363.0;   // K
    double egr_fraction     = 0.0;     // 0–1

    double omega() const {             // rad/s
        return rpm * 2.0 * Constants::PI * Constants::RPM_TO_RPS;
    }
    double cycle_time() const {        // s for 4-stroke
        return 2.0 * 60.0 / rpm;
    }
    double crank_deg_per_second() const {
        return rpm * 360.0 / 60.0;
    }
};

// ─────────────────────────────────────────────
//  Thermodynamic state
// ─────────────────────────────────────────────
struct ThermoState {
    double pressure    = 0.0;  // Pa
    double temperature = 0.0;  // K
    double density     = 0.0;  // kg/m³
    double mass        = 0.0;  // kg
    double volume      = 0.0;  // m³
    double internal_energy = 0.0; // J
    double enthalpy    = 0.0;  // J

    // Composition (mass fractions)
    double Y_fuel   = 0.0;
    double Y_O2     = 0.0;
    double Y_N2     = 0.0;
    double Y_CO2    = 0.0;
    double Y_H2O    = 0.0;
    double Y_CO     = 0.0;
    double Y_NOx    = 0.0;
    double Y_unburned = 0.0; // UHC

    double burned_fraction = 0.0; // 0–1
};

} // namespace Ricardo
