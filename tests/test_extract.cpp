// C-2: proc::sheet_resistance / proc::junction_depth tests.
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "cprocess/materials.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

static SimState make_column(double zmax_cm, int nz) {
  SimState st;
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, zmax_cm, 3, 3, nz);
  proc::set_region(st, "silicon", -1);
  return st;
}

int main() {
  // 1. sheet_resistance is positive, finite, and decreases with dose (more
  //    active carriers -> lower Rs) for a fixed shallow implant + anneal.
  double rs_lo, rs_hi;
  {
    SimState st = make_column(1.0e-4, 80);
    proc::implant_gauss(st, "B", 1e14, 0, 0.03e-4, 0.015e-4, 0, false, 0, 0, 0, 0,
                        false, "gauss");
    DiffuseOpts d; d.temp = 1273.15; d.time = 600; d.verbosity = 0;
    proc::diffuse(st, d);
    rs_lo = proc::sheet_resistance(st, "B");
    CHECK(std::isfinite(rs_lo) && rs_lo > 0);
    std::printf("Rs(B, 1e14, 1000C/10min) = %.4g Ohm/sq\n", rs_lo);
  }
  {
    SimState st = make_column(1.0e-4, 80);
    proc::implant_gauss(st, "B", 1e16, 0, 0.03e-4, 0.015e-4, 0, false, 0, 0, 0, 0,
                        false, "gauss");
    DiffuseOpts d; d.temp = 1273.15; d.time = 600; d.verbosity = 0;
    proc::diffuse(st, d);
    rs_hi = proc::sheet_resistance(st, "B");
    std::printf("Rs(B, 1e16, 1000C/10min) = %.4g Ohm/sq\n", rs_hi);
  }
  CHECK(rs_hi < rs_lo);  // 100x the dose -> substantially lower sheet Rs

  // 2. sheet_resistance over a depth window is >= the whole-column value
  //    (a narrower window captures less integrated conductance -> higher Rs).
  {
    SimState st = make_column(1.0e-4, 80);
    proc::implant_gauss(st, "B", 1e15, 0, 0.03e-4, 0.015e-4, 0, false, 0, 0, 0, 0,
                        false, "gauss");
    DiffuseOpts d; d.temp = 1273.15; d.time = 600; d.verbosity = 0;
    proc::diffuse(st, d);
    const double rs_full = proc::sheet_resistance(st, "B");
    const double rs_window = proc::sheet_resistance(st, "B", 0.0, 0.05e-4);
    CHECK(rs_window >= rs_full);
    std::printf("Rs full=%.4g, Rs[0,0.05um]=%.4g\n", rs_full, rs_window);
  }

  // 3. junction_depth: single-species fallback vs. background level.
  {
    SimState st = make_column(1.0e-4, 100);
    proc::implant_gauss(st, "B", 1e15, 0, 0.03e-4, 0.015e-4, 0, false, 0, 0, 0, 0,
                        false, "gauss");
    DiffuseOpts d; d.temp = 1273.15; d.time = 900; d.verbosity = 0;
    proc::diffuse(st, d);
    const double xj = proc::junction_depth(st, "B", 1e15);
    CHECK(xj > 0 && xj < 1.0e-4);
    std::printf("Xj(B, bg=1e15) = %.4g um\n", xj * 1e4);
    // A higher background level is reached at a shallower depth.
    const double xj_hi_bg = proc::junction_depth(st, "B", 1e17);
    CHECK(xj_hi_bg <= xj);
  }

  // 4. junction_depth with an opposing dopant: crossing lands between the two
  //    peaks for a simple P-well/B-implant style pair.
  {
    SimState st = make_column(1.0e-4, 100);
    proc::init(st, "P", 1e15);  // n-type background
    proc::implant_gauss(st, "B", 5e15, 0, 0.03e-4, 0.015e-4, 0, false, 0, 0, 0, 0,
                        false, "gauss");
    DiffuseOpts d; d.temp = 1273.15; d.time = 300; d.verbosity = 0;
    proc::diffuse(st, d);
    const double xj = proc::junction_depth(st, "B");
    CHECK(xj > 0 && xj < 1.0e-4);
    std::printf("Xj(B/P well) = %.4g um\n", xj * 1e4);
  }

  // 5. Errors on unknown species.
  {
    SimState st = make_column(1.0e-4, 10);
    bool threw = false;
    try { proc::sheet_resistance(st, "As"); }
    catch (const std::exception&) { threw = true; }
    CHECK(threw);
  }

  std::printf("extract tests passed\n");
  return 0;
}
