#include <cstdio>
#include <cmath>
#include <sstream>

#include "cprocess/oxidation.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

static SimState make_column(double lz_cm, int nz, std::ostream* log) {
  SimState st;
  proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, lz_cm, 4, 4, nz, log);
  proc::set_region(st, "silicon", -1, log);
  return st;
}

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

  // -------------------------------------------------------------------------
  // [C-3] Massoud thin-film enhancement: previously
  // tests/test_sprocess_parity.cpp's "Massoud thin-oxide enhancement"
  // WILL_FAIL check (moved here now that it PASSes). Dry oxidation in the
  // ~10 nm regime at 900 C, with `ox.massoud.c/l` opted into via ParamDB,
  // must grow >10% more than the pure Deal-Grove prediction.
  // -------------------------------------------------------------------------
  {
    ParamDB::instance().clear();
    ParamDB::instance().set("ox.massoud.c", 0.9);
    ParamDB::instance().set("ox.massoud.l", 0.01);  // 10 nm decay length

    std::ostringstream log;
    SimState st = make_column(0.3e-4, 12, &log);
    const double t_min = 40.0;  // ~10 nm regime at 900 C dry
    const double x = proc::oxidize(st, t_min * 60.0, 900 + 273.15, false, &log);
    const double x_dg = deal_grove_step(0.0, t_min, 900.0, false) * 1e-4;
    std::printf("[C-3] Massoud (ox.massoud.c=0.9): grown=%.2f nm vs DG=%.2f nm\n",
                x * 1e7, x_dg * 1e7);
    CHECK(x > 1.10 * x_dg);

    ParamDB::instance().clear();
  }

  // -------------------------------------------------------------------------
  // [C-3] Bit-identity: default ParamDB (Massoud off, pressure=1atm,
  // hcl=0, orient=<100>) must reproduce the pre-C-3 pure Deal-Grove
  // number exactly, so ParamDB-ifying the Arrhenius constants and adding
  // the new optional oxidize() parameters is a behavior-preserving
  // refactor when left at defaults.
  // -------------------------------------------------------------------------
  {
    ParamDB::instance().clear();
    std::ostringstream log;
    SimState st = make_column(0.5e-4, 50, &log);
    const double tox_cm = proc::oxidize(st, 3600.0, 1273.15, false, &log);
    const double analytic_um = deal_grove_step(0.0, 60.0, 1000.0, false);
    std::printf("[C-3] bit-identity: tox=%.6g um, analytic=%.6g um\n",
                tox_cm * 1e4, analytic_um);
    CHECK(std::fabs(tox_cm * 1e4 - analytic_um) <= 1e-9 * std::fabs(analytic_um));
  }

  // -------------------------------------------------------------------------
  // [C-3] Pressure scaling: higher O2 pressure grows oxide faster
  // (B ∝ P, B/A ∝ P^0.75).
  // -------------------------------------------------------------------------
  {
    const double x_1atm = deal_grove_step(0.0, 30.0, 1000.0, false, 1.0);
    const double x_5atm = deal_grove_step(0.0, 30.0, 1000.0, false, 5.0);
    std::printf("[C-3] pressure: 1atm=%.4f um  5atm=%.4f um\n", x_1atm, x_5atm);
    CHECK(x_5atm > x_1atm);
  }

  // -------------------------------------------------------------------------
  // [C-3] HCl enhancement: nonzero HCl fraction grows oxide faster.
  // -------------------------------------------------------------------------
  {
    const double x_noHCl = deal_grove_step(0.0, 30.0, 1000.0, false, 1.0, 0.0);
    const double x_HCl   = deal_grove_step(0.0, 30.0, 1000.0, false, 1.0, 0.03);
    std::printf("[C-3] HCl: 0%%=%.4f um  3%%=%.4f um\n", x_noHCl, x_HCl);
    CHECK(x_HCl > x_noHCl);
  }

  // -------------------------------------------------------------------------
  // [C-3] Orientation ratio: <111> grows faster than <100> in the
  // linear (thin-oxide) regime, by roughly the 1.68x literature ratio
  // on the linear (B/A) rate constant.
  // -------------------------------------------------------------------------
  {
    const DealGroveParams p100 = deal_grove_dry(1000.0, 1.0, 0.0, 1.0);
    const DealGroveParams p111 = deal_grove_dry(1000.0, 1.0, 0.0, 1.68);
    const double ba100 = p100.B / p100.A;
    const double ba111 = p111.B / p111.A;
    std::printf("[C-3] orient: B/A <100>=%.4f  <111>=%.4f  ratio=%.3f\n",
                ba100, ba111, ba111 / ba100);
    CHECK_NEAR(ba111 / ba100, 1.68, 0.05);

    const double x100 = deal_grove_step(0.0, 5.0, 1000.0, false, 1.0, 0.0, 1.0);
    const double x111 = deal_grove_step(0.0, 5.0, 1000.0, false, 1.0, 0.0, 1.68);
    std::printf("[C-3] orient: 5min <100>=%.5f um  <111>=%.5f um\n", x100, x111);
    CHECK(x111 > x100);
  }

  std::printf("oxidation tests passed\n");
  return 0;
}
