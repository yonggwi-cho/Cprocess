// W-8: material-aware implant transport regression tests.
//
// Moved from tests/test_sprocess_parity.cpp per the WILL_FAIL workflow
// (docs/tasks/README.md "SProcess パリティテスト"): all three [W-8] checks
// now PASS, so they live here as permanent regression tests instead of the
// executable-spec suite.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

static std::string mat_of(const SimState& st, int ci) {
  auto it = st.region_material.find(st.mesh.cell_region[ci]);
  return it == st.region_material.end() ? std::string("silicon") : it->second;
}

// Peak-concentration depth of `sp` below z_ref, silicon cells only (cm).
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

static SimState make_column(double zmax, int nz) {
  SimState st;
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, zmax, 3, 3, nz);
  proc::set_region(st, "silicon", -1);
  return st;
}

int main() {
  // ---- MC screen-oxide attenuation ----------------------------------------
  {
    SimState bare = make_column(0.5e-4, 50);
    const double z0 = bare.mesh.bbox().hi.z;
    proc::implant_mc(bare, "B", 1e13, 30.0, 40000, 0, 0, 11, 1, true, false, 0,
                     0, 0, 0, false, false);
    const double d_bare = si_peak_depth(bare, "B", z0);

    SimState ox = make_column(0.5e-4, 50);
    const double z_si = ox.mesh.bbox().hi.z;
    proc::deposit(ox, "oxide", 0.05e-4, 5, {});  // 50nm
    proc::implant_mc(ox, "B", 1e13, 30.0, 40000, 0, 0, 11, 1, true, false, 0, 0,
                     0, 0, false, false);
    const double d_ox = si_peak_depth(ox, "B", z_si);

    std::printf("MC screen-oxide: bare=%.1f nm, ox50=%.1f nm\n", d_bare * 1e7,
                d_ox * 1e7);
    CHECK(d_ox < 0.95 * d_bare);
  }

  // ---- All-Si MC bit-identity (legacy fast path unaffected) --------------
  {
    SimState a = make_column(0.5e-4, 50);
    SimState b = make_column(0.5e-4, 50);
    proc::implant_mc(a, "B", 1e13, 30.0, 40000, 0, 0, 11, 1, true, false, 0, 0,
                     0, 0, false, false);
    proc::implant_mc(b, "B", 1e13, 30.0, 40000, 0, 0, 11, 1, true, false, 0, 0,
                     0, 0, false, false);
    const auto& fa = a.fields.at("B");
    const auto& fb = b.fields.at("B");
    CHECK(fa.size() == fb.size());
    bool bitexact = true;
    for (std::size_t i = 0; i < fa.size(); ++i)
      if (fa[i] != fb[i]) { bitexact = false; break; }
    std::printf("all-Si MC bit-identity (same seed): %s\n",
                bitexact ? "yes" : "no");
    CHECK(bitexact);
  }

  // ---- STI shielding: oxide-filled trench blocks dopant reaching Si ------
  {
    SimState st;
    proc::mesh_box(st, 0, 0.6e-4, 0, 0.3e-4, 0, 0.6e-4, 6, 3, 30);
    proc::set_region(st, "silicon", -1);
    std::vector<std::pair<double, double>> tr = {
        {0, 0}, {0.3e-4, 0}, {0.3e-4, 0.3e-4}, {0, 0.3e-4}};
    proc::etch(st, 0.2e-4, tr, "");
    for (auto& [t, m] : st.region_material)
      if (m == "gas") m = "oxide";
    proc::implant_mc(st, "P", 1e13, 100.0, 60000, 0, 0, 5, 0, true, false, 0, 0,
                     0, 0, false, false);

    const auto& f = st.fields.at("P");
    double du = 0, db = 0;
    for (std::size_t i = 0; i < f.size(); ++i) {
      if (mat_of(st, (int)i) != "silicon") continue;
      const double q = f[i] * st.mesh.cell_vol[i];
      if (st.mesh.cell_cent[i].x < 0.3e-4) du += q;
      else                                 db += q;
    }
    std::printf("STI shielding: total P dose under-oxide=%.4e bare=%.4e atoms "
                "(ratio=%.3f)\n", du, db, du / db);
    CHECK(du < 0.5 * db);
  }

  // ---- Analytic screening offset ------------------------------------------
  {
    SimState bare = make_column(0.5e-4, 50);
    const double z0 = bare.mesh.bbox().hi.z;
    proc::implant_gauss(bare, "B", 1e13, 30.0, 0, 0, 0, false, 0, 0, 0, 0);
    const double d_bare = si_peak_depth(bare, "B", z0);

    SimState ox = make_column(0.5e-4, 50);
    const double z_si = ox.mesh.bbox().hi.z;
    proc::deposit(ox, "oxide", 0.05e-4, 5, {});
    proc::implant_gauss(ox, "B", 1e13, 30.0, 0, 0, 0, false, 0, 0, 0, 0);
    const double d_ox = si_peak_depth(ox, "B", z_si);

    std::printf("analytic screen-oxide: bare=%.1f nm, ox50=%.1f nm\n",
                d_bare * 1e7, d_ox * 1e7);
    CHECK(d_ox < 0.8 * d_bare);
  }

  // ---- All-Si analytic bit-identity ---------------------------------------
  {
    SimState a = make_column(0.5e-4, 50);
    SimState b = make_column(0.5e-4, 50);
    proc::implant_gauss(a, "B", 1e13, 30.0, 0, 0, 0, false, 0, 0, 0, 0);
    proc::implant_gauss(b, "B", 1e13, 30.0, 0, 0, 0, false, 0, 0, 0, 0);
    const auto& fa = a.fields.at("B");
    const auto& fb = b.fields.at("B");
    bool bitexact = true;
    for (std::size_t i = 0; i < fa.size(); ++i)
      if (fa[i] != fb[i]) { bitexact = false; break; }
    std::printf("all-Si analytic bit-identity: %s\n", bitexact ? "yes" : "no");
    CHECK(bitexact);
  }

  std::printf("test_implant_materials: all checks passed\n");
  return 0;
}
