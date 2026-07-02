// Tests for the P1-1 Pearson-IV analytic implant depth profile.

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "cprocess/implant.hpp"
#include "cprocess/materials.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

// Mass-weighted moments of a depth field: mean, sigma, and 3rd central
// moment mu3, all with d = ztop - z.
struct Moments { double mean, sigma, mu3; };

Moments depth_moments(const SimState& st, const std::vector<double>& c) {
  const BBox bb = st.mesh.bbox();
  double m = 0, md = 0, md2 = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    const double w = c[i] * st.mesh.cell_vol[i];
    const double d = bb.hi.z - st.mesh.cell_cent[i].z;
    m += w; md += w * d; md2 += w * d * d;
  }
  const double mean = (m > 0) ? md / m : 0.0;
  const double sigma2 = (m > 0) ? std::max(0.0, md2 / m - mean * mean) : 0.0;
  const double sigma = std::sqrt(sigma2);
  double m3 = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    const double w = c[i] * st.mesh.cell_vol[i];
    const double d = bb.hi.z - st.mesh.cell_cent[i].z;
    m3 += w * (d - mean) * (d - mean) * (d - mean);
  }
  const double mu3 = (m > 0) ? m3 / m : 0.0;
  return {mean, sigma, mu3};
}

double field_sum(const SimState& st, const std::vector<double>& c) {
  double s = 0;
  for (std::size_t i = 0; i < c.size(); ++i) s += c[i] * st.mesh.cell_vol[i];
  return s;
}

template <typename F>
void expect_throw(const char* what, F&& f) {
  bool threw = false;
  try {
    f();
  } catch (const std::exception&) {
    threw = true;
  }
  if (!threw) {
    std::fprintf(stderr, "FAILED: expected throw for %s\n", what);
    std::exit(1);
  }
}

}  // namespace

int main() {
  const double area = 0.2e-4 * 0.2e-4;

  // --- 1 & 2: dose conservation + moment match (B, 200 keV). ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    proc::set_region(st, "silicon", -1);
    const double dose = 1e14;
    const double atoms =
        proc::implant_gauss(st, "B", dose, 200.0, 0, 0, 0, false, 0, 0, 0, 0,
                            false, "pearson");
    CHECK(std::fabs(atoms / (dose * area) - 1.0) < 0.005);

    const auto& f = st.fields.at("B");
    const double fsum = field_sum(st, f);
    CHECK(std::fabs(fsum / (dose * area) - 1.0) < 0.005);

    const Moments mom = depth_moments(st, f);
    const double rp = 531e-7, drp = 94e-7;
    CHECK(std::fabs(mom.mean - rp) < 0.02 * rp);
    CHECK(std::fabs(mom.sigma - drp) < 0.05 * drp);

    // --- 3: negative skewness for pearson; near-zero for gauss. ---
    CHECK(mom.mu3 < 0.0);

    SimState stg;
    proc::mesh_box(stg, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    proc::set_region(stg, "silicon", -1);
    proc::implant_gauss(stg, "B", dose, 200.0, 0, 0, 0, false, 0, 0, 0, 0,
                        false, "gauss");
    const auto& fg = stg.fields.at("B");
    const Moments momg = depth_moments(stg, fg);
    const double skew_g = std::fabs(momg.mu3) / (momg.sigma * momg.sigma * momg.sigma);
    CHECK(skew_g < 0.1);
  }

  // --- 4: invalid Pearson-IV moments fall back to Gaussian, bit-identical. ---
  {
    SimState st1, st2;
    proc::mesh_box(st1, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    proc::mesh_box(st2, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);

    const Dopant* b = find_dopant("B");
    CHECK(b != nullptr);

    ImplantParams p1;
    p1.dopant = b;
    p1.dose = 1e14;
    p1.rp = 531e-7; p1.drp = 94e-7;
    p1.profile = ImplantParams::Profile::pearson4;
    p1.gamma = -2.0;
    p1.beta = 4.0;  // beta=4.0 <= 1+gamma^2=5.0 -> invalid Type-IV

    ImplantParams p2 = p1;
    p2.profile = ImplantParams::Profile::gauss;

    std::vector<double> c1, c2;
    std::vector<char> mask;  // empty = all cells
    const double atoms1 = apply_implant(st1.mesh, mask, p1, c1);
    const double atoms2 = apply_implant(st2.mesh, mask, p2, c2);

    CHECK(atoms1 == atoms2);
    CHECK(c1.size() == c2.size());
    for (std::size_t i = 0; i < c1.size(); ++i) CHECK(c1[i] == c2[i]);
  }

  // --- 5: throw on pearson w/o energy, and on unknown profile string. ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    proc::set_region(st, "silicon", -1);
    expect_throw("pearson without energy", [&] {
      proc::implant_gauss(st, "B", 1e13, 0.0, 0.1e-4, 0.03e-4, 0, false, 0, 0,
                          0, 0, false, "pearson");
    });
    expect_throw("unknown profile", [&] {
      proc::implant_gauss(st, "B", 1e13, 50.0, 0, 0, 0, false, 0, 0, 0, 0,
                          false, "foo");
    });
  }

  // --- 6: default profile arg -> non-regression Gaussian. ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240);
    proc::set_region(st, "silicon", -1);
    const double atoms =
        proc::implant_gauss(st, "B", 1e14, 50.0, 0, 0, 0, false, 0, 0, 0, 0);
    CHECK(atoms > 0);
  }

  std::printf("test_pearson: OK\n");
  return 0;
}
