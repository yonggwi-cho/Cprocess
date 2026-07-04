// Tests for P2-3: oxidation-enhanced diffusion (OED) via interstitial
// injection at the growing Si/SiO2 interface, wired into proc::oxidize().

#include <cmath>
#include <cstdio>
#include <sstream>

#include "cprocess/materials.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// Mass-weighted std-dev of a 1-D profile along z (cf. tests/test_ted.cpp),
// restricted to silicon cells only (so mesh regridding/oxide conversion
// during oxidize() doesn't skew the comparison).
static double profile_spread_si(const SimState& st, const std::vector<double>& c) {
  const std::vector<char> mask = proc::silicon_mask(st);
  double m = 0, mz = 0, mz2 = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    if (!mask[i]) continue;
    const double w = c[i] * st.mesh.cell_vol[i];
    const double z = st.mesh.cell_cent[i].z;
    m += w; mz += w * z; mz2 += w * z * z;
  }
  if (m <= 0) return 0;
  const double mean = mz / m;
  return std::sqrt(std::max(0.0, mz2 / m - mean * mean));
}

static double total_mass_all(const SimState& st, const std::string& sym) {
  const auto& f = st.fields.at(sym);
  double m = 0;
  for (std::size_t i = 0; i < f.size(); ++i) m += f[i] * st.mesh.cell_vol[i];
  return m;
}

// Common setup: 0.2x0.2x0.8 um box (cell height 0.01 um, same as
// tests/test_oxidize_flow.cpp), B 1e18 cm^-3 Gaussian at Rp=0.3 um.
// NOTE: the spec's own recipe (0.5 um box, Rp=0.1 um) is too shallow for a
// wet-1000C-30min anneal here: Deal-Grove growth for those conditions is
// tox~0.29 um (matches literature for wet O2 at 1000 C), so 0.44*0.29 =
// ~0.13 um of Si is consumed -- deeper than Rp=0.1 um, meaning the *entire*
// implant peak would be swallowed into the frozen oxide-converted band
// regardless of theta. That makes the "theta=0 barely changes the profile"
// comparison (test 2) fail for a purely geometric reason (mass deletion),
// not a TED bug. Rp is deepened to 0.3 um (box to 0.8 um) so the Si
// consumption stays well short of the implant peak.
static void setup(SimState& st) {
  proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.8e-4, 4, 4, 80);
  proc::set_region(st, "silicon", -1);
  proc::implant_gauss(st, "B", 1e18 * 0.02e-4, 0.0, 0.3e-4, 0.02e-4, 0.0,
                      false, 0, 0, 0, 0, false);
}

int main() {
  ParamDB::instance().clear();

  // --- 1. OED enhancement: oxidizing (theta default 0.01) spreads B more
  //        than an inert (theta=0) anneal of the same thermal budget. ---
  double spread_A, spread_B;
  {
    SimState stA;
    setup(stA);
    ParamDB::instance().set("oed.theta", 0.01);
    proc::oxidize(stA, 1800.0, 1273.15, true);  // 30 min, 1000 C, wet
    spread_A = profile_spread_si(stA, stA.fields.at("B"));
  }
  {
    SimState stB;
    setup(stB);
    ParamDB::instance().set("oed.theta", 0.0);  // inert: geometry-only
    proc::oxidize(stB, 1800.0, 1273.15, true);
    spread_B = profile_spread_si(stB, stB.fields.at("B"));
  }
  std::printf("OED spread: A(theta=0.01)=%.4g um, B(theta=0, inert)=%.4g um\n",
              spread_A * 1e4, spread_B * 1e4);
  CHECK(spread_A > 1.2 * spread_B);
  std::printf("OED enhancement factor: %.3fx\n", spread_A / spread_B);

  // --- 2. theta=0 disables OED entirely and reproduces the P1-6 legacy
  //        oxidize() behavior: geometry (Deal-Grove + retag + nearest-
  //        centroid regridding) ONLY, no thermal diffusion at all -- so the
  //        Si-side profile spread is essentially unchanged by the call
  //        (DoD: "theta=0 matches P1-6's traditional (geometry-only)
  //        behavior"). Compare against the spread immediately after implant,
  //        before any oxidize() call.
  double spread0;
  {
    SimState st0;
    setup(st0);
    spread0 = profile_spread_si(st0, st0.fields.at("B"));
  }
  std::printf("theta=0 (inert) vs pre-oxidize spread: B=%.4g um, pre=%.4g um\n",
              spread_B * 1e4, spread0 * 1e4);
  CHECK(std::fabs(spread_B - spread0) <= 0.02 * spread0);

  // --- 3. Dose conservation: state A's total B (Si + oxide-converted cells)
  //        is conserved within 5% across the oxidize() call. ---
  {
    SimState st;
    setup(st);
    const double mass0 = total_mass_all(st, "B");
    ParamDB::instance().set("oed.theta", 0.01);
    proc::oxidize(st, 1800.0, 1273.15, true);
    const double mass1 = total_mass_all(st, "B");
    std::printf("dose conservation: before=%.6e after=%.6e (%.3g%% change)\n",
                mass0, mass1, 100.0 * std::fabs(mass1 - mass0) / mass0);
    CHECK(std::fabs(mass1 - mass0) <= 0.05 * mass0);
  }

  // --- 4. The "I" field is positive after oxidize(), peaking near the new
  //        Si/SiO2 interface. ---
  {
    SimState st;
    setup(st);
    ParamDB::instance().set("oed.theta", 0.01);
    const double tox_cm = proc::oxidize(st, 1800.0, 1273.15, true);
    CHECK(st.fields.count("I") == 1);
    const auto& I = st.fields.at("I");
    double maxI = 0;
    int argmax = -1;
    for (std::size_t i = 0; i < I.size(); ++i) {
      if (I[i] > maxI) { maxI = I[i]; argmax = static_cast<int>(i); }
    }
    std::printf("I field: max=%.4g cm^-3 at cell %d\n", maxI, argmax);
    CHECK(maxI > 0);
    CHECK(argmax >= 0);

    // Interface height: bbox top - grown oxide thickness.
    const BBox bb = st.mesh.bbox();
    const double z_if = bb.hi.z - tox_cm;
    const double h = (bb.hi.z - bb.lo.z + tox_cm) / 80.0;  // approx cell height
    std::printf("  z_if=%.4g um, argmax centroid z=%.4g um, 2h=%.4g um\n",
                z_if * 1e4, st.mesh.cell_cent[argmax].z * 1e4, 2 * h * 1e4);
    CHECK(st.mesh.cell_cent[argmax].z >= z_if - 2 * h);
  }

  ParamDB::instance().clear();
  std::printf("\nall oed tests passed\n");
  return 0;
}
