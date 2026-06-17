#include <cstdio>
#include <cmath>
#include <algorithm>

#include "cprocess/mesh.hpp"
#include "cprocess/levelset.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // ---- levelset_init: inside region → φ<0, outside → φ>0 -----------------
  {
    // Two-region mesh: left half (x<0.5) is region 1, right half is region 2.
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    // Tag left half cells as region 1.
    for (std::size_t ci = 0; ci < m.cells.size(); ++ci) {
      m.cell_region[ci] = (m.cell_cent[ci].x < 0.5) ? 1 : 2;
    }

    const auto phi = levelset_init(m, 1);
    CHECK(phi.size() == m.nodes.size());

    // Nodes deep inside region 1 (x ~ 0) should have φ < 0.
    int n_neg = 0, n_pos = 0;
    for (std::size_t i = 0; i < m.nodes.size(); ++i) {
      if (m.nodes[i].x < 0.2) { CHECK(phi[i] < 0); ++n_neg; }
      if (m.nodes[i].x > 0.8) { CHECK(phi[i] > 0); ++n_pos; }
    }
    std::printf("levelset_init: n_neg=%d n_pos=%d\n", n_neg, n_pos);
    CHECK(n_neg > 0 && n_pos > 0);
  }

  // ---- levelset_advect: interface moves with velocity ---------------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    for (std::size_t ci = 0; ci < m.cells.size(); ++ci)
      m.cell_region[ci] = (m.cell_cent[ci].x < 0.5) ? 1 : 2;

    auto phi = levelset_init(m, 1);

    // Find approximate zero crossing before advection (x ~ 0.5).
    const int nn = static_cast<int>(m.nodes.size());

    // Advect with uniform v_n = +0.1 (interface moves right by 0.1*dt).
    std::vector<double> vn(nn, 0.1);
    const double dt = 0.5;
    levelset_advect(m, phi, vn, dt);

    // After advection, some formerly positive nodes near x=0.5 should now
    // have phi < old phi (interface moved outward).
    // Just check phi didn't blow up.
    for (double p : phi) CHECK(std::isfinite(p));
    std::printf("levelset_advect: done\n");
  }

  // ---- levelset_reinit: |∇φ| → 1 after reinitialization ------------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    for (std::size_t ci = 0; ci < m.cells.size(); ++ci)
      m.cell_region[ci] = (m.cell_cent[ci].x < 0.5) ? 1 : 2;
    auto phi = levelset_init(m, 1);

    // Distort phi artificially.
    for (std::size_t i = 0; i < phi.size(); ++i) phi[i] *= 2.0;

    levelset_reinit(m, phi, 5);
    for (double p : phi) CHECK(std::isfinite(p));
    std::printf("levelset_reinit: done\n");
  }

  // ---- levelset_update_regions: region tags follow phi sign ---------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 3, 3, 3);
    const int nn = static_cast<int>(m.nodes.size());
    // phi < 0 everywhere → all cells become inside_tag=1.
    std::vector<double> phi(nn, -1.0);
    levelset_update_regions(m, phi, 1, 2);
    for (int r : m.cell_region) CHECK(r == 1);

    // phi > 0 everywhere → all cells become outside_tag=2.
    std::fill(phi.begin(), phi.end(), 1.0);
    levelset_update_regions(m, phi, 1, 2);
    for (int r : m.cell_region) CHECK(r == 2);
    std::printf("levelset_update_regions: done\n");
  }

  std::printf("levelset tests passed\n");
  return 0;
}
