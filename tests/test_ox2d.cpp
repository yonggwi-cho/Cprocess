// Tests for P2-4 Stage B: 2D LOCOS oxidation (proc::oxidize_2d), bird's-beak
// lateral encroachment under a nitride mask.

#include <cmath>
#include <cstdio>
#include <sstream>

#include "cprocess/deck.hpp"
#include "cprocess/oxidation.hpp"
#include "cprocess/process.hpp"
#include "cprocess/remesh.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

// Oxide thickness at a given x (averaged over y, over cells whose xy
// centroid lies within `tol_x` of `x_um`): (max z of an oxide cell there)
// minus (max z of a silicon cell there).
double oxide_thickness_at_x(const SimState& st, double x_um, double tol_x_um) {
  const double x_cm = x_um * 1e-4, tol_cm = tol_x_um * 1e-4;
  double z_ox = -1e300, z_si = -1e300;
  bool any_ox = false, any_si = false;
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    const Vec3& c = st.mesh.cell_cent[i];
    if (std::fabs(c.x - x_cm) > tol_cm) continue;
    auto it = st.region_material.find(st.mesh.cell_region[i]);
    const std::string mat = (it == st.region_material.end()) ? "silicon" : it->second;
    const auto& cell = st.mesh.cells[i];
    if (mat == "oxide" || mat == "sio2") {
      any_ox = true;
      for (int k = 0; k < 4; ++k) z_ox = std::max(z_ox, st.mesh.nodes[cell[k]].z);
    } else if (mat == "silicon" || mat == "si") {
      any_si = true;
      for (int k = 0; k < 4; ++k) z_si = std::max(z_si, st.mesh.nodes[cell[k]].z);
    }
  }
  if (!any_ox || !any_si) return 0.0;
  return (z_ox - z_si) * 1e4;  // um
}

double total_mass(const SimState& st, const std::string& sym) {
  const auto& f = st.fields.at(sym);
  double m = 0;
  for (std::size_t i = 0; i < f.size(); ++i) m += f[i] * st.mesh.cell_vol[i];
  return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Test 1: planar (blanket) sanity -- oxidize_2d requires a nitride mask, so
// force a 2D path with a dummy nitride stripe far from the measurement point
// and check the open-field center thickness against deal_grove_step.
// ---------------------------------------------------------------------------
static void test_open_field_deal_grove() {
  std::printf("test_open_field_deal_grove\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.1e-4, 0, 0.5e-4, 8, 2, 40, &log);
  proc::set_region(st, "silicon", -1, &log);
  // Nitride stripe over x in [0.3, 0.4] um; open field elsewhere.
  proc::deposit(st, "nitride", 0.05e-4, 2, {}, &log);
  proc::etch(st, 0.05e-4, {{0.3e-4, 0}, {0.4e-4, 0}, {0.4e-4, 0.1e-4}, {0.3e-4, 0.1e-4}},
            "nitride", &log);

  const double tox_cm = proc::oxidize_2d(st, 30.0 * 60.0, 1273.15, true, &log);
  const double tox_um = tox_cm * 1e4;
  const double analytic_um = deal_grove_step(0, 30.0, 1000.0, /*wet=*/true);

  std::printf("  open-field tox=%.4g um analytic=%.4g um\n", tox_um, analytic_um);
  const double rel = std::fabs(tox_um - analytic_um) / analytic_um;
  CHECK(rel < 0.15);
}

// ---------------------------------------------------------------------------
// Test 2/3/4: bird's beak taper, mesh sanity, conservation.
// ---------------------------------------------------------------------------
static void test_birds_beak() {
  std::printf("test_birds_beak\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.2e-4, 0, 0.5e-4, 8, 4, 50, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e18, -1, &log);
  const double mass_before = total_mass(st, "B");

  proc::deposit(st, "nitride", 0.05e-4, 2, {}, &log);
  // Open the field for x < 0.2 um (nitride remains for x in [0.2, 0.4]).
  proc::etch(st, 0.05e-4,
            {{0, 0}, {0.2e-4, 0}, {0.2e-4, 0.2e-4}, {0, 0.2e-4}},
            "nitride", &log);

  const double vol_before = st.mesh.total_volume();

  proc::oxidize_2d(st, 60.0 * 60.0, 1273.15, true, &log);

  QualityStats qs = mesh_quality(st.mesh);
  std::printf("  mesh min_q after oxidize_2d = %.4g\n", qs.min_q);
  CHECK(qs.min_q > 0.02);

  const double t_open = oxide_thickness_at_x(st, 0.1, 0.05);
  const double t_edge = oxide_thickness_at_x(st, 0.3, 0.05);
  std::printf("  t(open, x=0.1um)=%.4g um  t(under mask, x=0.3um)=%.4g um\n",
             t_open, t_edge);
  CHECK(t_open > 0);
  const double ratio = t_edge / t_open;
  std::printf("  ratio = %.4g\n", ratio);
  CHECK(ratio > 0.20 && ratio < 0.80);

  double prev = 1e300;
  bool monotone = true;
  for (double x : {0.22, 0.26, 0.30, 0.34}) {
    const double t = oxide_thickness_at_x(st, x, 0.03);
    std::printf("    t(x=%.2f um) = %.4g um\n", x, t);
    if (t > prev * 1.05) monotone = false;
    prev = t;
  }
  CHECK(monotone);

  const double mass_after = total_mass(st, "B");
  const double rel_mass = std::fabs(mass_after - mass_before) / mass_before;
  std::printf("  B mass rel change = %.4g\n", rel_mass);
  CHECK(rel_mass < 0.05);

  const double vol_after = st.mesh.total_volume();
  const double rel_vol = std::fabs(vol_after - vol_before) / vol_before;
  std::printf("  volume rel change = %.4g (expected some growth from oxide swell)\n", rel_vol);
  CHECK(vol_after > vol_before);
}

int main() {
  test_open_field_deal_grove();
  test_birds_beak();
  std::printf("test_ox2d: all tests passed\n");
  return 0;
}
