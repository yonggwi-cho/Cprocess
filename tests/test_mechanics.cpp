#include <cstdio>
#include <cmath>

#include "cprocess/mechanics.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // ---- Stress relaxes to zero with no new strain --------------------------
  {
    ViscoElasticParams p;
    p.tau_relax_min = 10.0;

    StressState s{};
    s.s[0] = 1000.0;  // 1 GPa initial stress

    const double dstrain[6] = {};
    for (int i = 0; i < 50; ++i)
      s = maxwell_update(s, dstrain, 10.0, p);

    // After 50 * 10 min = 500 min >> tau=10 min: stress ~ 0
    std::printf("relaxed stress = %.2e MPa\n", s.s[0]);
    CHECK(std::fabs(s.s[0]) < 1e-5);
  }

  // ---- Elastic response: no relaxation (tau very large) -------------------
  {
    ViscoElasticParams p;
    p.tau_relax_min = 1e12;  // essentially infinite
    p.E_GPa = 130.0;
    p.nu = 0.0;  // simplified: no Poisson coupling

    // Uniaxial strain in x, small step.
    const double eps = 1e-4;
    const double dstrain[6] = {eps, 0, 0, 0, 0, 0};
    StressState s{};
    s = maxwell_update(s, dstrain, 1.0, p);

    // σ_xx = E * ε (nu=0, lam=0 → only 2*mu = E contribution for normal).
    // With nu=0: lam=0, mu=E/2, σ_xx = 2*mu*eps = E*eps
    const double expected = 130e3 * eps;  // MPa
    std::printf("elastic σ_xx = %.2f MPa (expected %.2f)\n", s.s[0], expected);
    CHECK_NEAR(s.s[0], expected, 0.01 * expected);
  }

  // ---- mechanics_step: operates on vector of states -----------------------
  {
    ViscoElasticParams p;
    p.tau_relax_min = 1e12;
    p.E_GPa = 130.0;
    p.nu = 0.28;

    const int nc = 4;
    std::vector<StressState> states(nc);
    std::vector<double> dstrain(nc * 6, 0.0);
    // Apply uniaxial strain to first cell only.
    dstrain[0] = 1e-4;

    mechanics_step(states, dstrain, 1.0, p);
    CHECK(states[0].s[0] > 0.0);    // first cell has stress
    CHECK(states[1].s[0] == 0.0);   // others still zero
    std::printf("mechanics_step: s[0].xx=%.2f  s[1].xx=%.2f\n",
                states[0].s[0], states[1].s[0]);
  }

  std::printf("mechanics tests passed\n");
  return 0;
}
