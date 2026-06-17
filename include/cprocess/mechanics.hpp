#pragma once
#include <vector>

namespace cp {

// Per-cell Maxwell viscoelastic state: stress tensor (Voigt: xx,yy,zz,xy,yz,xz).
struct StressState {
  double s[6] = {};  // MPa
};

struct ViscoElasticParams {
  double E_GPa = 130.0;   // Young's modulus (Si: ~130 GPa)
  double nu    = 0.28;    // Poisson's ratio
  double tau_relax_min = 1e6;  // stress relaxation time (min); large = elastic
};

// Update stress state for one time step using Maxwell model:
//   σ_new = σ_old * exp(-dt/τ) + C : dε
// dstrain[6]: incremental strain tensor (Voigt), same units as σ/E.
// Returns updated stress state.
StressState maxwell_update(const StressState& s_old,
                           const double dstrain[6],
                           double dt_min,
                           const ViscoElasticParams& p);

// Update all cell stress states in place.
// dstrain_cells: per-cell incremental strain (nc * 6 values, row-major).
void mechanics_step(std::vector<StressState>& states,
                    const std::vector<double>& dstrain_cells,
                    double dt_min,
                    const ViscoElasticParams& p);

}  // namespace cp
