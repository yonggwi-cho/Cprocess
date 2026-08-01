#pragma once

namespace cp {

struct DealGroveParams {
  double A;    // linear rate constant (um)
  double B;    // parabolic rate constant (um²/min)
  double tau;  // time offset to account for initial oxide (min)
};

// [C-3] Standard Deal-Grove rate constants for <100> Si. The base Arrhenius
// pre-exponentials/activation energies are ParamDB-backed
// (ox.dry.b0/be/a0/ae, ox.wet.b0/be/a0/ae) with defaults matching the
// original hardcoded Deal & Grove (1965) Table I values exactly, so leaving
// them unset reproduces the original numbers bit-identically.
//
// pressure_atm: B scales linearly with O2 partial pressure (B ∝ P), B/A
// scales as P^0.75 (Deal & Grove 1965, high-pressure oxidation extension).
// Defaults to 1.0 atm (no scaling, bit-identical).
//
// hcl_frac: fractional HCl in the oxidizing ambient (e.g. 0.03 for 3%
// HCl); grows both B and B/A by (1 + ox.hcl.gain * hcl_frac). Defaults to
// 0.0 (no scaling, bit-identical). ox.hcl.gain ParamDB key, default 6.0
// (empirical; ~3% HCl gives ~15-20% growth-rate enhancement).
//
// orient_factor: multiplies B/A (the linear, reaction-rate-limited term)
// to account for crystal-orientation dependence of the Si-SiO2 interfacial
// reaction rate; <111> Si oxidizes ~1.68x faster (linearly) than <100> Si
// (the parabolic/diffusion-limited B term is orientation-independent).
// Defaults to 1.0 (i.e. <100>, bit-identical).
DealGroveParams deal_grove_dry(double T_celsius, double pressure_atm = 1.0,
                               double hcl_frac = 0.0,
                               double orient_factor = 1.0);
DealGroveParams deal_grove_wet(double T_celsius, double pressure_atm = 1.0,
                               double hcl_frac = 0.0,
                               double orient_factor = 1.0);

// Integrate Deal-Grove for one step (pure Deal-Grove -- no Massoud
// thin-film correction; see deal_grove_step_massoud below for that).
// x0_um: initial oxide thickness (um).
// dt_min: time step (min).
// Returns new oxide thickness (um).
double deal_grove_step(double x0_um, double dt_min, double T_celsius,
                       bool wet = false, double pressure_atm = 1.0,
                       double hcl_frac = 0.0, double orient_factor = 1.0);

// [C-3] Deal-Grove step with an optional Massoud (1985) thin-oxide growth
// enhancement: the local growth rate dx/dt is multiplied by
// (1 + massoud_c * exp(-x / massoud_l)) before being integrated, which
// boosts the initial (thin-film, <~30 nm) dry-oxidation rate above the
// classic Deal-Grove parabolic/linear law without perturbing thick-film
// growth (the exponential decays to ~0 well before typical Tier-A
// benchmark thicknesses). massoud_c == 0 (default) skips the numerical
// integration entirely and delegates to deal_grove_step() for an
// identical result to the pre-C-3 code path.
double deal_grove_step_massoud(double x0_um, double dt_min,
                               double T_celsius, bool wet,
                               double massoud_c, double massoud_l,
                               double pressure_atm = 1.0,
                               double hcl_frac = 0.0,
                               double orient_factor = 1.0);

}  // namespace cp
