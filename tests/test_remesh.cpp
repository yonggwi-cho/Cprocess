#include <cstdio>

#include "cprocess/mesh.hpp"
#include "cprocess/remesh.hpp"
#include "cprocess/topology.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // ---- tet_quality: regular tet ~ 1, sliver ~ 0 ----------------------------
  {
    // A near-regular tetrahedron.
    const Vec3 a{0, 0, 0}, b{1, 0, 0}, c{0.5, 0.8660254, 0},
        d{0.5, 0.2886751, 0.8164966};
    const double q = tet_quality(a, b, c, d);
    std::printf("regular tet quality = %.4f\n", q);
    CHECK(q > 0.95 && q <= 1.0 + 1e-9);

    // A sliver: nearly coplanar.
    const double qs = tet_quality({0, 0, 0}, {1, 0, 0}, {0, 1, 0},
                                  {0.3, 0.3, 1e-4});
    std::printf("sliver quality = %.6f\n", qs);
    CHECK(qs < 0.05);
  }

  // ---- MeshTopology on a box mesh -----------------------------------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 2, 2, 2);
    MeshTopology topo;
    topo.build(m);

    // Every cell appears in exactly its 6 edges.
    std::size_t edge_cell_refs = 0;
    for (const auto& [k, cells] : topo.edge_cells) {
      (void)k;
      edge_cell_refs += cells.size();
    }
    CHECK(edge_cell_refs == m.cells.size() * 6);

    // Corner node (0,0,0) is a boundary node.
    int corner = -1;
    for (std::size_t i = 0; i < m.nodes.size(); ++i)
      if (norm(m.nodes[i]) < 1e-12) { corner = static_cast<int>(i); break; }
    CHECK(corner >= 0);
    CHECK(topo.node_boundary[corner] == 1);

    // The dead-center node of a 2x2x2 box is interior.
    int center = -1;
    for (std::size_t i = 0; i < m.nodes.size(); ++i)
      if (norm(m.nodes[i] - Vec3{0.5, 0.5, 0.5}) < 1e-12)
        { center = static_cast<int>(i); break; }
    CHECK(center >= 0);
    CHECK(topo.node_boundary[center] == 0);
  }

  // ---- edge_split preserves total volume and grows the cell count ----------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 2, 2, 2);
    const double vol0 = m.total_volume();
    const std::size_t nc0 = m.cells.size();

    // Pick an existing edge: scan the first cell's first edge.
    const int a = m.cells[0][0], b = m.cells[0][1];
    MeshTopology topo;
    topo.build(m);
    const std::size_t ninc = topo.edge_incident(a, b).size();
    CHECK(ninc > 0);

    const int mid = edge_split(m, a, b);
    CHECK(mid >= 0);
    m.finalize();

    CHECK(m.cells.size() == nc0 + ninc);
    CHECK_NEAR(m.total_volume(), vol0, 1e-12 * vol0);
    // New node sits at the edge midpoint.
    // (geometry only; volumes already validated by finalize()).
    for (double v : m.cell_vol) CHECK(v > 0);
  }

  // ---- laplacian_smooth: volume preserved, boundary fixed, no inversion ----
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    const double vol0 = m.total_volume();

    // Record boundary node positions.
    MeshTopology topo;
    topo.build(m);
    std::vector<Vec3> bnd_before;
    std::vector<int> bnd_ids;
    for (std::size_t i = 0; i < m.nodes.size(); ++i)
      if (topo.node_boundary[i]) { bnd_ids.push_back(i); bnd_before.push_back(m.nodes[i]); }

    const int moved = laplacian_smooth(m, 3, 0.5);
    m.finalize();
    std::printf("laplacian moved %d nodes\n", moved);

    // Total volume of a closed box is conserved exactly (boundary pinned).
    CHECK_NEAR(m.total_volume(), vol0, 1e-10 * vol0);
    for (double v : m.cell_vol) CHECK(v > 0);
    // Boundary nodes did not move.
    for (std::size_t k = 0; k < bnd_ids.size(); ++k)
      CHECK(norm(m.nodes[bnd_ids[k]] - bnd_before[k]) < 1e-15);
  }

  std::printf("remesh tests passed\n");
  return 0;
}
