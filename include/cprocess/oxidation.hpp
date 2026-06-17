#pragma once

namespace cp {

struct DealGroveParams {
  double A;    // linear rate constant (um)
  double B;    // parabolic rate constant (um²/min)
  double tau;  // time offset to account for initial oxide (min)
};

// Standard Deal-Grove rate constants for <100> Si.
// T_celsius: oxidation temperature.
DealGroveParams deal_grove_dry(double T_celsius);
DealGroveParams deal_grove_wet(double T_celsius);

// Integrate Deal-Grove for one step.
// x0_um: initial oxide thickness (um).
// dt_min: time step (min).
// Returns new oxide thickness (um).
double deal_grove_step(double x0_um, double dt_min, double T_celsius,
                       bool wet = false);

}  // namespace cp
