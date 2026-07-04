#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/materials.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// Standard deviation (spread) of a 1-D profile along z, mass-weighted.
static double profile_spread(const SimState& st, const std::vector<double>& c) {
  double m = 0, mz = 0, mz2 = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    const double w = c[i] * st.mesh.cell_vol[i];
    const double z = st.mesh.cell_cent[i].z;
    m += w; mz += w * z; mz2 += w * z * z;
  }
  if (m <= 0) return 0;
  const double mean = mz / m;
  return std::sqrt(std::max(0.0, mz2 / m - mean * mean));
}

static double total_dose(const SimState& st, const std::string& sym) {
  const auto& c = st.fields.at(sym);
  double tot = 0;
  for (std::size_t i = 0; i < c.size(); ++i) tot += c[i] * st.mesh.cell_vol[i];
  return tot;
}

static void setup_mesh(SimState& st) {
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40);
  proc::set_region(st, "silicon", -1);
}

static DiffuseOpts opts(double temp_k, double time_s) {
  DiffuseOpts o;
  o.temp = temp_k;
  o.time = time_s;
  o.verbosity = 0;
  o.nonortho = false;  // see test_dopant_clusters.cpp: orthogonal box mesh,
                       // avoids a deferred-correction mass-conservation
                       // artifact unrelated to the physics under test
  return o;
}

int main() {
  // --- 1. Basic operation: implant + diffuse for each of In/C/F/Ge. ---
  for (const char* sym : {"In", "C", "F", "Ge"}) {
    SimState st;
    setup_mesh(st);
    proc::implant_gauss(st, sym, 1e14, 50.0, 0, 0, 0, false, 0, 0, 0, 0,
                        /*damage=*/false);
    const double dose0 = total_dose(st, sym);
    CHECK(dose0 > 0);
    proc::diffuse(st, opts(1273.15, 600));  // 1000 C, 10 min
    const double dose1 = total_dose(st, sym);
    const double rel = std::fabs(dose1 - dose0) / dose0;
    std::printf("new_dopants: %s basic op dose change=%.4g%%\n", sym,
               rel * 100.0);
    CHECK(rel < 5e-3);  // < 0.5%
    for (double v : st.fields.at(sym)) CHECK(std::isfinite(v));
  }

  // --- 2. MC implant works for In and Ge. ---
  for (const char* sym : {"In", "Ge"}) {
    SimState st;
    setup_mesh(st);
    McImplantStats s = proc::implant_mc(st, sym, 1e14, 30.0, 20000, 0.0, 0.0,
                                        12345ULL, 1, /*channeling=*/false,
                                        false, 0, 0, 0, 0, /*damage=*/false);
    std::printf("new_dopants: %s MC deposited=%lld\n", sym, s.deposited);
    CHECK(s.deposited > 0);
  }

  // --- 3. Neutral species doesn't shift charge neutrality (nni). ---
  {
    SimState st1, st2;
    setup_mesh(st1);
    setup_mesh(st2);
    proc::init(st1, "P", 1e18);
    proc::init(st2, "P", 1e18);
    proc::init(st2, "Ge", 1e20);
    const DiffuseOpts o = opts(1273.15, 300);
    proc::diffuse(st1, o);
    proc::diffuse(st2, o);
    const auto& p1 = st1.fields.at("P");
    const auto& p2 = st2.fields.at("P");
    CHECK(p1.size() == p2.size());
    double maxrel = 0;
    for (std::size_t i = 0; i < p1.size(); ++i) {
      const double denom = std::fabs(p1[i]) + 1.0;
      maxrel = std::max(maxrel, std::fabs(p1[i] - p2[i]) / denom);
    }
    std::printf("new_dopants: P with/without Ge background maxrel=%.4g\n",
               maxrel);
    CHECK(maxrel < 1e-9);
  }

  // --- 4. Ge is immobile. ---
  {
    SimState st;
    setup_mesh(st);
    proc::implant_gauss(st, "Ge", 1e15, 50.0, 0, 0, 0, false, 0, 0, 0, 0, false);
    std::vector<double> before = st.fields.at("Ge");
    proc::diffuse(st, opts(1273.15, 600));
    const auto& after = st.fields.at("Ge");
    double maxrel = 0;
    for (std::size_t i = 0; i < before.size(); ++i) {
      const double denom = std::fabs(before[i]) + 1.0;
      maxrel = std::max(maxrel, std::fabs(after[i] - before[i]) / denom);
    }
    std::printf("new_dopants: Ge immobility maxrel=%.4g\n", maxrel);
    CHECK(maxrel < 1e-12);
  }

  // --- 5. C suppresses TED (B+C co-implant vs. B alone). ---
  {
    double spread_B_alone, spread_B_with_C;
    const DiffuseOpts o = opts(1173.15, 60);  // 900 C, 1 min
    {
      SimState st;
      setup_mesh(st);
      proc::implant_gauss(st, "B", 1e14, 0.0, 0.05e-4, 0.02e-4, 0.0, false,
                          0, 0, 0, 0, /*damage=*/true);
      const double s0 = profile_spread(st, st.fields.at("B"));
      proc::diffuse_ted(st, o);
      spread_B_alone = profile_spread(st, st.fields.at("B")) - s0;
    }
    {
      SimState st;
      setup_mesh(st);
      proc::implant_gauss(st, "B", 1e14, 0.0, 0.05e-4, 0.02e-4, 0.0, false,
                          0, 0, 0, 0, /*damage=*/true);
      proc::init(st, "C", 1e19);
      const double s0 = profile_spread(st, st.fields.at("B"));
      proc::diffuse_ted(st, o);
      spread_B_with_C = profile_spread(st, st.fields.at("B")) - s0;
      const auto& ccl = st.fields.at("C_cl");
      double ccl_sum = 0;
      for (std::size_t i = 0; i < ccl.size(); ++i)
        ccl_sum += ccl[i] * st.mesh.cell_vol[i];
      std::printf("new_dopants: C_cl accumulated dose=%.4g\n", ccl_sum);
      CHECK(ccl_sum > 0);
    }
    std::printf("new_dopants: TED spread increment B alone=%.4g um, "
               "B+C=%.4g um\n",
               spread_B_alone * 1e4, spread_B_with_C * 1e4);
    CHECK(spread_B_with_C <= 0.8 * spread_B_alone);  // >= 20% reduction
  }

  // --- 6. In diffuses more slowly than B. ---
  {
    const DiffuseOpts o = opts(1273.15, 1800);  // 1000 C, 30 min
    double spread_B, spread_In;
    {
      SimState st;
      setup_mesh(st);
      proc::implant_gauss(st, "B", 1e14, 0.0, 0.05e-4, 0.02e-4, 0.0, false,
                          0, 0, 0, 0, false);
      const double s0 = profile_spread(st, st.fields.at("B"));
      proc::diffuse(st, o);
      spread_B = profile_spread(st, st.fields.at("B")) - s0;
    }
    {
      SimState st;
      setup_mesh(st);
      proc::implant_gauss(st, "In", 1e14, 0.0, 0.05e-4, 0.02e-4, 0.0, false,
                          0, 0, 0, 0, false);
      const double s0 = profile_spread(st, st.fields.at("In"));
      proc::diffuse(st, o);
      spread_In = profile_spread(st, st.fields.at("In")) - s0;
    }
    std::printf("new_dopants: 1000C/30min spread increment B=%.4g um, "
               "In=%.4g um\n",
               spread_B * 1e4, spread_In * 1e4);
    CHECK(spread_In < spread_B);
  }

  std::printf("new dopants tests passed\n");
  return 0;
}
