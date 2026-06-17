#include <cstdio>
#include <cmath>

#include "cprocess/mesh.hpp"
#include "cprocess/field_transfer.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // ---- Transfer a uniform field: should be exact --------------------------
  {
    Mesh m_old = make_box_mesh(0, 1, 0, 1, 0, 1, 3, 3, 3);
    Mesh m_new = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);

    const int nc_old = static_cast<int>(m_old.cells.size());
    std::vector<double> field_old(nc_old, 1e18);

    const auto field_new = transfer_field_nearest(m_old, field_old, m_new);

    // Uniform field: every transferred value should be 1e18.
    for (double v : field_new) CHECK_NEAR(v, 1e18, 1e10);

    const double err = check_mass_conservation(m_old, field_old, m_new, field_new);
    std::printf("uniform field mass error: %.2e\n", err);
    CHECK(err < 0.01);
  }

  // ---- Transfer a Gaussian profile: mass conservation < 5% ----------------
  {
    Mesh m_old = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    Mesh m_new = make_box_mesh(0, 1, 0, 1, 0, 1, 5, 5, 5);

    const int nc_old = static_cast<int>(m_old.cells.size());
    std::vector<double> field_old(nc_old);
    for (int i = 0; i < nc_old; ++i) {
      const Vec3& c = m_old.cell_cent[i];
      const double r2 = (c.x - 0.5) * (c.x - 0.5) +
                        (c.y - 0.5) * (c.y - 0.5) +
                        (c.z - 0.5) * (c.z - 0.5);
      field_old[i] = 1e18 * std::exp(-r2 / (2.0 * 0.1 * 0.1));
    }

    const auto field_new = transfer_field_nearest(m_old, field_old, m_new);
    const double err = check_mass_conservation(m_old, field_old, m_new, field_new);
    std::printf("Gaussian field mass error: %.2f%%\n", err * 100.0);
    CHECK(err < 0.10);  // nearest-centroid: within 10%
  }

  std::printf("field_transfer tests passed\n");
  return 0;
}
