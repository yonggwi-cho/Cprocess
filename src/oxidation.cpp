#include "cprocess/oxidation.hpp"
#include <cmath>

namespace cp {

// Arrhenius: k = A_pre * exp(-Ea / (kB * T))
// kB = 8.617e-5 eV/K
static double arrhenius(double A_pre, double Ea_eV, double T_celsius) {
  const double kB = 8.617333e-5;
  const double T_K = T_celsius + 273.15;
  return A_pre * std::exp(-Ea_eV / (kB * T_K));
}

// Standard Deal-Grove parameters from Deal & Grove (1965), Table I.
// Dry O2, <100> Si:
//   B (parabolic) = 7.72e2 * exp(-1.23/kT)  um²/min
//   B/A (linear)  = 6.23e6 * exp(-2.00/kT)  um/min
DealGroveParams deal_grove_dry(double T_celsius) {
  const double B    = arrhenius(7.72e2, 1.23, T_celsius);
  const double B_A  = arrhenius(6.23e6, 2.00, T_celsius);
  const double A    = B / B_A;
  // tau: time to grow initial 25 Å native oxide
  const double x_i  = 0.0025;  // um (25 Å, native oxide)
  const double tau  = (x_i * x_i + A * x_i) / B;  // hours
  return {A, B, tau};
}

// Wet O2 (H2O), <100> Si:
//   B  = 3.86e2 * exp(-0.78/kT)  um²/min
//   B/A= 1.63e8 * exp(-2.05/kT)  um/min
DealGroveParams deal_grove_wet(double T_celsius) {
  const double B   = arrhenius(3.86e2, 0.78, T_celsius);
  const double B_A = arrhenius(1.63e8, 2.05, T_celsius);
  const double A   = B / B_A;
  const double x_i = 0.0;  // wet: no significant native oxide offset
  const double tau = (x_i * x_i + A * x_i) / B;
  return {A, B, tau};
}

double deal_grove_step(double x0_um, double dt_min, double T_celsius, bool wet) {
  const DealGroveParams p = wet ? deal_grove_wet(T_celsius)
                                : deal_grove_dry(T_celsius);
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

}  // namespace cp
