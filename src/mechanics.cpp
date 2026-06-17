#include "cprocess/mechanics.hpp"
#include <cmath>

namespace cp {

StressState maxwell_update(const StressState& s_old,
                           const double dstrain[6],
                           double dt_min,
                           const ViscoElasticParams& p) {
  const double decay = std::exp(-dt_min / p.tau_relax_min);
  const double E = p.E_GPa * 1e3;  // MPa
  const double nu = p.nu;

  // Isotropic elastic stiffness (Voigt): C * dε components.
  // C_ijkl for isotropic: λ*δ_ij*δ_kl + μ*(δ_ik*δ_jl + δ_il*δ_jk)
  const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const double mu  = E / (2.0 * (1.0 + nu));
  const double tr  = dstrain[0] + dstrain[1] + dstrain[2];

  StressState s_new;
  // Normal components.
  for (int i = 0; i < 3; ++i)
    s_new.s[i] = s_old.s[i] * decay + lam * tr + 2.0 * mu * dstrain[i];
  // Shear components (factor 2 in Voigt convention for shear strains).
  for (int i = 3; i < 6; ++i)
    s_new.s[i] = s_old.s[i] * decay + mu * dstrain[i];

  return s_new;
}

void mechanics_step(std::vector<StressState>& states,
                    const std::vector<double>& dstrain_cells,
                    double dt_min,
                    const ViscoElasticParams& p) {
  const int nc = static_cast<int>(states.size());
  for (int i = 0; i < nc; ++i)
    states[i] = maxwell_update(states[i], &dstrain_cells[i * 6], dt_min, p);
}

}  // namespace cp
