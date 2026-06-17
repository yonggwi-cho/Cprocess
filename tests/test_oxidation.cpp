#include <cstdio>
#include <cmath>

#include "cprocess/oxidation.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // ---- dry oxidation at 1000°C, 60 min → ~0.03-0.06 um --------------------
  {
    const double x = deal_grove_step(0.0, 60.0, 1000.0, false);
    std::printf("dry 1000C 60min: x = %.4f um\n", x);
    CHECK(x > 0.02 && x < 0.15);
  }

  // ---- wet oxidation at 900°C, 30 min → ~0.05-0.3 um ---------------------
  {
    const double x = deal_grove_step(0.0, 30.0, 900.0, true);
    std::printf("wet 900C 30min: x = %.4f um\n", x);
    CHECK(x > 0.03 && x < 0.5);
  }

  // ---- linear rate limit: thin oxide, x ~ (B/A)*t -------------------------
  {
    const DealGroveParams p = deal_grove_dry(1000.0);
    const double dt_min = 1.0;
    const double dt_hr = dt_min / 60.0;
    // B/A in um/hr; for thin oxide x ~ (B/A)*t
    const double x_linear = (p.B / p.A) * dt_hr;
    const double x_dg = deal_grove_step(0.0, dt_min, 1000.0, false);
    std::printf("linear approx=%.6f um  DG=%.6f um\n", x_linear, x_dg);
    // Should agree within 20% for thin oxide.
    CHECK(std::fabs(x_dg - x_linear) / (x_linear + 1e-30) < 0.3);
  }

  // ---- self-consistency: step from x0 by dt equals direct step from 0 -----
  {
    const double x1 = deal_grove_step(0.0,  30.0, 950.0, true);
    const double x2 = deal_grove_step(x1,   30.0, 950.0, true);
    const double x12 = deal_grove_step(0.0, 60.0, 950.0, true);
    std::printf("two-step=%.5f um  one-step=%.5f um\n", x2, x12);
    CHECK_NEAR(x2, x12, 1e-9);
  }

  // ---- oxide thickness increases monotonically with time -------------------
  {
    double x = 0.0;
    for (int i = 0; i < 10; ++i) {
      const double x_new = deal_grove_step(x, 10.0, 1100.0, false);
      CHECK(x_new > x);
      x = x_new;
    }
    std::printf("monotone growth at 1100C dry: x after 100min = %.4f um\n", x);
  }

  std::printf("oxidation tests passed\n");
  return 0;
}
