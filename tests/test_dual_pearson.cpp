// [C-1] Tests for the extended (1 keV .. 3 MeV) implant moment table and the
// dual-Pearson (primary Pearson-IV + channeling-tail) analytic implant
// profile. See docs/tasks/C1_implant_moments.md for the full spec and data
// provenance.

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "cprocess/implant.hpp"
#include "cprocess/materials.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

static std::string mat_of(const SimState& st, int ci) {
  auto it = st.region_material.find(st.mesh.cell_region[ci]);
  return it == st.region_material.end() ? std::string("silicon") : it->second;
}

static double si_peak_depth(const SimState& st, const std::string& sp,
                            double z_ref) {
  const auto& f = st.fields.at(sp);
  double best = -1, zbest = z_ref;
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (mat_of(st, (int)i) != "silicon") continue;
    if (f[i] > best) { best = f[i]; zbest = st.mesh.cell_cent[i].z; }
  }
  return z_ref - zbest;
}

static SimState make_column(double zmax, int nz, std::ostream* log) {
  SimState st;
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, zmax, 3, 3, nz, log);
  proc::set_region(st, "silicon", -1, log);
  return st;
}

static double conc_at(const SimState& s, const std::string& sp, double z0,
                      double depth) {
  const auto& f = s.fields.at(sp);
  double best = 1e300, c = 0;
  for (std::size_t i = 0; i < f.size(); ++i) {
    const double d = std::fabs((z0 - s.mesh.cell_cent[i].z) - depth);
    if (d < best) { best = d; c = f[i]; }
  }
  return c;
}

}  // namespace

int main() {
  // --- 1: [C-1] the parity-flip check itself -- B 40 keV dual profile must
  // not underestimate the MC channeling tail concentration @2*Rp by more
  // than 10x. (Moved here from tests/test_sprocess_parity.cpp once it began
  // passing; see docs/tasks/C1_implant_moments.md for the calibration.)
  {
    std::ostringstream log;
    SimState an = make_column(0.8e-4, 80, &log);
    SimState mc = make_column(0.8e-4, 80, &log);
    proc::implant_gauss(an, "B", 1e14, 40.0, 0, 0, 0, false, 0, 0, 0, 0, false,
                        "dual", &log);
    proc::implant_mc(mc, "B", 1e14, 40.0, 200000, 0, 0, 3, 0,
                     /*channeling=*/true, false, 0, 0, 0, 0, false, false,
                     &log);
    const double z0 = an.mesh.bbox().hi.z;
    const double rp = si_peak_depth(an, "B", z0);
    const double ca = conc_at(an, "B", z0, 2 * rp);
    const double cm = conc_at(mc, "B", z0, 2 * rp);
    CHECK(cm > 0);
    CHECK(ca >= 0.1 * cm);
  }

  // --- 2: profile="dual" conserves total dose (primary + tail sum to the
  // requested dose, same tolerance as the existing plain-pearson test). ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    proc::set_region(st, "silicon", -1);
    const double dose = 1e14;
    const double area = 0.2e-4 * 0.2e-4;
    const double atoms = proc::implant_gauss(st, "B", dose, 80.0, 0, 0, 0,
                                             false, 0, 0, 0, 0, false, "dual");
    CHECK(std::fabs(atoms / (dose * area) - 1.0) < 0.01);
  }

  // --- 3: dual falls back to plain pearson4 (zero-tail) behavior when
  // dp_frac/dp_l are zero (ImplantParams defaults) -- bit-identical to
  // pearson4 with the same moments. ---
  {
    SimState st1, st2;
    proc::mesh_box(st1, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    proc::mesh_box(st2, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    const Dopant* b = find_dopant("B");
    CHECK(b != nullptr);
    ImplantParams p1;
    p1.dopant = b; p1.dose = 1e14; p1.rp = 531e-7; p1.drp = 94e-7;
    p1.profile = ImplantParams::Profile::pearson4;
    p1.gamma = -1.5; p1.beta = 8.0;
    ImplantParams p2 = p1;
    p2.profile = ImplantParams::Profile::dual;  // dp_frac=dp_l=0 -> no tail
    std::vector<double> c1, c2;
    std::vector<char> mask;
    const double a1 = apply_implant(st1.mesh, mask, p1, c1);
    const double a2 = apply_implant(st2.mesh, mask, p2, c2);
    CHECK(a1 == a2);
    CHECK(c1.size() == c2.size());
    for (std::size_t i = 0; i < c1.size(); ++i) CHECK(c1[i] == c2[i]);
  }

  // --- 4: ParamDB overridability of the tail: "B.dp.frac"=0 must reproduce
  // the plain-pearson4 (no-tail) dose profile through the process.cpp path,
  // and a larger explicit frac must raise the @1.5*Rp tail concentration. ---
  {
    ParamDB::instance().set("B.dp.frac", 0.0);
    SimState st0;
    proc::mesh_box(st0, 0, 0.3e-4, 0, 0.3e-4, 0, 0.8e-4, 3, 3, 160);
    proc::set_region(st0, "silicon", -1);
    proc::implant_gauss(st0, "B", 1e14, 40.0, 0, 0, 0, false, 0, 0, 0, 0,
                        false, "dual");
    const double z0 = st0.mesh.bbox().hi.z;
    const double rp0 = si_peak_depth(st0, "B", z0);
    const double c_notail = conc_at(st0, "B", z0, 1.5 * rp0);

    ParamDB::instance().set("B.dp.frac", 0.20);
    ParamDB::instance().set("B.dp.decay_mult", 2.0);
    SimState st1;
    proc::mesh_box(st1, 0, 0.3e-4, 0, 0.3e-4, 0, 0.8e-4, 3, 3, 160);
    proc::set_region(st1, "silicon", -1);
    proc::implant_gauss(st1, "B", 1e14, 40.0, 0, 0, 0, false, 0, 0, 0, 0,
                        false, "dual");
    const double c_tail = conc_at(st1, "B", z0, 1.5 * rp0);
    CHECK(c_tail > c_notail);

    ParamDB::instance().clear();
  }

  // --- 5: wider (1 keV .. 3 MeV) range table spot-checks, LSS-derivation
  // sanity per docs/tasks/C1_implant_moments.md. All new points preserve the
  // existing 10-200 keV anchor rows bit-identically (tested via the
  // permanent Tier-A benchmarks in test_benchmarks.cpp); here we only check
  // qualitative LSS trends the extension must satisfy at the new points. ---
  {
    const Dopant* b = find_dopant("B");
    CHECK(b != nullptr);
    double rp, drp, gamma, beta;

    // (a) monotonic increase in Rp with energy across the full 1 keV..3 MeV
    // span (both the retained core and the new low/high extensions).
    double prev_rp = -1;
    const double kEnergies[] = {1, 2, 5, 10, 30, 80, 200, 400, 1000, 3000};
    for (double e : kEnergies) {
      CHECK(implant_moments(*b, e, rp, drp, gamma, beta));
      CHECK(rp > prev_rp);
      prev_rp = rp;
    }

    // (b) high-energy asymptote: log-log slope of Rp(E) near 3 MeV should be
    // close to the LSS electronic-stopping value 0.5 (sqrt(E) scaling), i.e.
    // sub-linear and well below the near-linear low-energy slope.
    double rp_hi1, drp_hi1, g_, be_;
    double rp_hi2, drp_hi2;
    CHECK(implant_moments(*b, 2000, rp_hi1, drp_hi1, g_, be_));
    CHECK(implant_moments(*b, 3000, rp_hi2, drp_hi2, g_, be_));
    const double n_hi = std::log(rp_hi2 / rp_hi1) / std::log(3000.0 / 2000.0);
    CHECK(n_hi > 0.3 && n_hi < 0.7);

    // (c) gamma/beta relax toward the Gaussian point (0, 3) at the 3 MeV
    // extension endpoint (documented approximation, see materials.cpp).
    CHECK(implant_moments(*b, 3000, rp, drp, gamma, beta));
    CHECK(std::fabs(gamma) < 0.05);
    CHECK(std::fabs(beta - 3.0) < 0.05);

    // (d) below-table extension (1 keV) stays positive/finite and below the
    // 10 keV anchor.
    double rp_lo, drp_lo;
    CHECK(implant_moments(*b, 1, rp_lo, drp_lo, gamma, beta));
    CHECK(rp_lo > 0 && rp_lo < 33e-7 /* 33 nm, the 10 keV anchor */);
  }

  std::printf("test_dual_pearson: OK\n");
  return 0;
}
