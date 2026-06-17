#include <cstdio>
#include <cmath>

#include "cprocess/mesh.hpp"
#include "cprocess/topology.hpp"
#include "cprocess/ale_mover.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // ---- compute_node_normals: boundary nodes get unit normals ---------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 2, 2, 2);
    MeshTopology topo;
    topo.build(m);
    const auto normals = compute_node_normals(m, topo);

    // Every boundary node should have a non-zero unit normal.
    int bnd_count = 0;
    for (std::size_t i = 0; i < m.nodes.size(); ++i) {
      if (topo.node_boundary[i]) {
        const double len = norm(normals[i]);
        CHECK(len > 0.5);  // unit normal
        ++bnd_count;
      }
    }
    CHECK(bnd_count > 0);
    std::printf("boundary nodes with normals: %d\n", bnd_count);
  }

  // ---- ale_move: moving zmax surface outward increases bbox ----------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 3, 3, 3);
    MeshTopology topo;
    topo.build(m);
    const auto normals = compute_node_normals(m, topo);

    const double zmax_before = m.bbox().hi.z;

    // Only move nodes on the zmax patch (normal pointing in +z direction).
    // v_n = +1 (outward), dt = 0.05.
    auto mask_fn = [&](int i) {
      // Block all boundary nodes except those on zmax (z ~ 1.0).
      if (!topo.node_boundary[i]) return true;
      return m.nodes[i].z < 0.99;
    };
    const auto res = ale_move(m, topo, normals, 1.0, 0.05, mask_fn, 2);
    m.finalize();

    std::printf("ale_move: moved=%d skipped=%d\n", res.n_moved, res.n_skipped);
    CHECK(res.n_moved > 0);

    const double zmax_after = m.bbox().hi.z;
    std::printf("zmax before=%.4f after=%.4f\n", zmax_before, zmax_after);
    CHECK(zmax_after > zmax_before + 0.01);

    // All cell volumes remain positive.
    for (double v : m.cell_vol) CHECK(v > 0);
  }

  // ---- rescale_fields_for_volume_change: conserves total mass --------------
  {
    // Two cells, field = {1e18, 2e18}, volumes double.
    std::vector<double> old_vol = {1.0e-12, 2.0e-12};
    std::vector<double> new_vol = {2.0e-12, 4.0e-12};
    std::vector<std::vector<double>> fields = {{1e18, 2e18}};

    const double mass_before = fields[0][0] * old_vol[0] + fields[0][1] * old_vol[1];
    rescale_fields_for_volume_change(fields, old_vol, new_vol);
    const double mass_after = fields[0][0] * new_vol[0] + fields[0][1] * new_vol[1];

    std::printf("mass before=%.6e after=%.6e\n", mass_before, mass_after);
    CHECK_NEAR(mass_after, mass_before, 1e-6 * mass_before);
  }

  std::printf("ale tests passed\n");
  return 0;
}
