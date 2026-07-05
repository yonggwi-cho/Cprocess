// Test for P2-4 Stage A: steady-state oxidant diffusion-reaction solver
// (cp::solve_oxidant) against the planar Deal-Grove limit.

#include <cmath>
#include <cstdio>
#include <sstream>

#include "cprocess/mesh.hpp"
#include "cprocess/oxidant_solver.hpp"
#include "cprocess/oxidation.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

constexpr double kNOx = 2.25e22;   // SiO2 molecular density, cm^-3
constexpr double kCgas = 5.2e16;   // dissolved oxidant conc, dry 1000C, cm^-3

// Deal-Grove B / (B/A), converted to cm^2/s and cm/s (Arrhenius constants in
// oxidation.cpp are per-hour, um: B [um^2/hr], B/A [um/hr]).
void deal_grove_rates(double temp_c, bool wet, double& d_ox, double& ks) {
  const DealGroveParams p = wet ? deal_grove_wet(temp_c) : deal_grove_dry(temp_c);
  const double B_cm2_s = p.B * 1e-8 / 3600.0;      // um^2/hr -> cm^2/s
  const double BA_cm_s = (p.B / p.A) * 1e-4 / 3600.0;  // um/hr -> cm/s
  d_ox = B_cm2_s * kNOx / (2.0 * kCgas);
  ks = BA_cm_s * kNOx / kCgas;
}

// Build a planar oxide-over-silicon mesh: total column height = n_ox+n_si
// layers of thickness h = x_o_cm / n_ox, oxide occupying the top x_o_cm.
Mesh make_planar(double x_o_cm, int n_ox, int n_si) {
  const double h = x_o_cm / n_ox;
  const double z_total = h * (n_ox + n_si);
  Mesh m = make_box_mesh(0, 1e-4, 0, 1e-4, 0, z_total, 1, 1, n_ox + n_si);
  const double z_if = z_total - x_o_cm;
  for (std::size_t i = 0; i < m.cells.size(); ++i)
    m.cell_region[i] = (m.cell_cent[i].z > z_if) ? 1 : 0;
  return m;
}

// Area-weighted average interface velocity (individual Kuhn-tet interface
// faces vary in area/orientation even on a nominally flat planar interface,
// so a plain arithmetic mean over faces is not representative).
double avg_iface_vn(const Mesh& m, const OxidantResult& r) {
  double num = 0, den = 0;
  for (const auto& [f, v] : r.iface_vn) {
    const double area = norm(m.faces[f].S);
    num += v * area;
    den += area;
  }
  return den > 0 ? num / den : 0.0;
}

}  // namespace

// ---------------------------------------------------------------------------
// Planar Deal-Grove agreement at 3 oxide thicknesses (dry, 1000 C).
// ---------------------------------------------------------------------------
static void test_planar_deal_grove() {
  std::printf("test_planar_deal_grove\n");
  const double temp_c = 1000.0;
  double d_ox, ks;
  deal_grove_rates(temp_c, /*wet=*/false, d_ox, ks);

  for (double x_o_um : {0.01, 0.05, 0.2}) {
    const double x_o_cm = x_o_um * 1e-4;
    Mesh m = make_planar(x_o_cm, 10, 10);

    std::vector<char> oxide_mask(m.cells.size()), si_mask(m.cells.size());
    for (std::size_t i = 0; i < m.cells.size(); ++i) {
      oxide_mask[i] = (m.cell_region[i] == 1);
      si_mask[i] = (m.cell_region[i] == 0);
    }

    std::ostringstream log;
    const OxidantResult res =
        solve_oxidant(m, oxide_mask, si_mask, d_ox, ks, kCgas, &log);
    CHECK(!res.iface_vn.empty());
    const double vn_num = avg_iface_vn(m, res);

    // Analytic reference: dx/dt from a small finite difference of
    // deal_grove_step around x_o_um.
    const double dt_min = 1e-4;  // small step, min
    const double x1 = deal_grove_step(x_o_um, dt_min, temp_c, /*wet=*/false);
    const double vn_ref_um_min = (x1 - x_o_um) / dt_min;
    const double vn_ref_cm_s = vn_ref_um_min * 1e-4 / 60.0;

    const double rel = std::fabs(vn_num - vn_ref_cm_s) / vn_ref_cm_s;
    std::printf("  x_o=%.3g um: v_n(solver)=%.4g cm/s v_n(ref)=%.4g cm/s rel=%.3g\n",
               x_o_um, vn_num, vn_ref_cm_s, rel);
    CHECK(rel < 0.10);
  }
}

int main() {
  test_planar_deal_grove();
  std::printf("test_oxidant: all tests passed\n");
  return 0;
}
