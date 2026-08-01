#include "cprocess/oxidation.hpp"
#include "cprocess/param_db.hpp"
#include <cmath>

namespace cp {

// Arrhenius: k = A_pre * exp(-Ea / (kB * T))
static double arrhenius(double A_pre, double Ea_eV, double T_celsius) {
  const double kB = 8.617333e-5;
  const double T_K = T_celsius + 273.15;
  return A_pre * std::exp(-Ea_eV / (kB * T_K));
}

// Applies the pressure/HCl/orientation scaling common to dry and wet paths.
// B scales linearly with O2 partial pressure; B/A scales as P^0.75
// (Deal & Grove 1965 high-pressure extension). HCl grows both B and B/A by
// (1 + gain*hcl_frac). orient_factor multiplies B/A only (interfacial
// reaction rate is orientation-dependent; the parabolic/diffusion-limited
// B term is not).
static void apply_ambient_scaling(double& B, double& B_A, double pressure_atm,
                                  double hcl_frac, double orient_factor) {
  const double hcl_gain = ParamDB::instance().get("ox.hcl.gain", 6.0);
  const double hcl_mult = 1.0 + hcl_gain * hcl_frac;
  B *= pressure_atm * hcl_mult;
  B_A *= std::pow(pressure_atm, 0.75) * hcl_mult * orient_factor;
}

// [C-3] Deal-Grove B/A Arrhenius constants, ParamDB-ified. Defaults match
// the original hardcoded Deal & Grove (1965) Table I values exactly, so
// leaving the keys unset reproduces the original numbers bit-identically.
//
// Dry O2, <100> Si:
//   B (parabolic) = 7.72e2 * exp(-1.23/kT)  um²/min
//   B/A (linear)  = 6.23e6 * exp(-2.00/kT)  um/min
DealGroveParams deal_grove_dry(double T_celsius, double pressure_atm,
                               double hcl_frac, double orient_factor) {
  const ParamDB& db = ParamDB::instance();
  const double b0 = db.get("ox.dry.b0", 7.72e2);
  const double be = db.get("ox.dry.be", 1.23);
  const double a0 = db.get("ox.dry.a0", 6.23e6);
  const double ae = db.get("ox.dry.ae", 2.00);
  double B    = arrhenius(b0, be, T_celsius);
  double B_A  = arrhenius(a0, ae, T_celsius);
  apply_ambient_scaling(B, B_A, pressure_atm, hcl_frac, orient_factor);
  const double A    = B / B_A;
  // tau: time to grow initial 25 Å native oxide
  const double x_i  = 0.0025;  // um (25 Å, native oxide)
  const double tau  = (x_i * x_i + A * x_i) / B;  // hours
  return {A, B, tau};
}

// Wet O2 (H2O), <100> Si:
//   B  = 3.86e2 * exp(-0.78/kT)  um²/min
//   B/A= 1.63e8 * exp(-2.05/kT)  um/min
DealGroveParams deal_grove_wet(double T_celsius, double pressure_atm,
                               double hcl_frac, double orient_factor) {
  const ParamDB& db = ParamDB::instance();
  const double b0 = db.get("ox.wet.b0", 3.86e2);
  const double be = db.get("ox.wet.be", 0.78);
  const double a0 = db.get("ox.wet.a0", 1.63e8);
  const double ae = db.get("ox.wet.ae", 2.05);
  double B   = arrhenius(b0, be, T_celsius);
  double B_A = arrhenius(a0, ae, T_celsius);
  apply_ambient_scaling(B, B_A, pressure_atm, hcl_frac, orient_factor);
  const double A   = B / B_A;
  const double x_i = 0.0;  // wet: no significant native oxide offset
  const double tau = (x_i * x_i + A * x_i) / B;
  return {A, B, tau};
}

double deal_grove_step(double x0_um, double dt_min, double T_celsius,
                       bool wet, double pressure_atm, double hcl_frac,
                       double orient_factor) {
  const DealGroveParams p = wet
      ? deal_grove_wet(T_celsius, pressure_atm, hcl_frac, orient_factor)
      : deal_grove_dry(T_celsius, pressure_atm, hcl_frac, orient_factor);
  // Arrhenius constants are per hour; convert dt from minutes to hours.
  const double dt_hr = dt_min / 60.0;
  const double t_eq = (x0_um * x0_um + p.A * x0_um) / p.B;
  const double t_new = t_eq + dt_hr;
  // x = A/2 * (sqrt(1 + 4*B/A^2 * t) - 1)
  const double A2 = p.A / 2.0;
  const double arg = 1.0 + t_new * p.B / (A2 * A2);
  if (arg <= 0) return x0_um;
  return A2 * (std::sqrt(arg) - 1.0);
}

double deal_grove_step_massoud(double x0_um, double dt_min, double T_celsius,
                               bool wet, double massoud_c, double massoud_l,
                               double pressure_atm, double hcl_frac,
                               double orient_factor) {
  if (massoud_c == 0.0) {
    return deal_grove_step(x0_um, dt_min, T_celsius, wet, pressure_atm,
                           hcl_frac, orient_factor);
  }
  const DealGroveParams p = wet
      ? deal_grove_wet(T_celsius, pressure_atm, hcl_frac, orient_factor)
      : deal_grove_dry(T_celsius, pressure_atm, hcl_frac, orient_factor);
  // Numerically integrate dx/dt = B / (2x + A) * (1 + C*exp(-x/L)) with
  // fixed-step RK4 (N substeps -- dt_min is typically a small oxidize()
  // sub-step already, so N=20 gives ample accuracy without being slow).
  const int N = 20;
  const double h = (dt_min / 60.0) / N;  // hours (B/A are per-hour Arrhenius rates)
  auto rate = [&](double x) {
    const double base = p.B / (2.0 * x + p.A);
    const double enh = 1.0 + massoud_c * std::exp(-x / massoud_l);
    return base * enh;
  };
  double x = x0_um;
  for (int i = 0; i < N; ++i) {
    const double k1 = rate(x);
    const double k2 = rate(x + 0.5 * h * k1);
    const double k3 = rate(x + 0.5 * h * k2);
    const double k4 = rate(x + h * k3);
    x += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    if (x < 0) x = 0;
  }
  return x;
}

}  // namespace cp
