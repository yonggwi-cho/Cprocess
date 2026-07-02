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

  // ---- split_edges: finalize is included, single edge --------------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    const int a = m.cells[0][0], b = m.cells[0][1];
    SplitResult r = split_edges(m, {{a, b}});
    CHECK(r.n_split == 1);
    CHECK(r.n_skipped == 0);
    CHECK(m.faces.size() > 0);
    CHECK(m.cell_vol.size() == m.cells.size());
    CHECK(mesh_quality(m).min_q > 0);
  }

  // ---- split_edges: volume preservation -----------------------------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    const double vol0 = m.total_volume();
    const int a = m.cells[0][0], b = m.cells[0][1];
    split_edges(m, {{a, b}});
    CHECK_NEAR(m.total_volume(), vol0, 1e-12 * vol0);
  }

  // ---- split_edges: batch split with conflict skip ------------------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    // Two edges of the same cell: they share incident cells.
    const int v0 = m.cells[0][0], v1 = m.cells[0][1], v2 = m.cells[0][2];
    SplitResult r = split_edges(m, {{v0, v1}, {v0, v2}});
    CHECK(r.n_split == 1);
    CHECK(r.n_skipped == 1);
  }

  // ---- redistribute_field: value and mass preservation --------------------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    const double vol0 = m.total_volume();
    double mass0 = 0;
    for (double v : m.cell_vol) mass0 += 1e15 * v;

    const int v0 = m.cells[0][0], v1 = m.cells[0][1];
    const int v2 = m.cells[3][0], v3 = m.cells[3][1];
    std::vector<double> conc0(m.cells.size(), 1e15);
    SplitResult r = split_edges(m, {{v0, v1}, {v2, v3}});
    CHECK(r.n_split >= 1);

    std::vector<double> conc1 = redistribute_field(conc0, r.cell_parent);
    CHECK(conc1.size() == m.cells.size());
    for (double c : conc1) CHECK(std::fabs(c - 1e15) < 1e-6);

    double mass1 = 0;
    for (std::size_t i = 0; i < conc1.size(); ++i) mass1 += conc1[i] * m.cell_vol[i];
    CHECK_NEAR(m.total_volume(), vol0, 1e-12 * vol0);
    CHECK_NEAR(mass1, mass0, 1e-12 * mass0);
  }

  // ---- split_edges: many random edges, quality and volume preserved -------
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 4, 4, 4);
    CHECK(mesh_quality(m).min_q > 0.1);
    const double vol0 = m.total_volume();

    MeshTopology topo;
    topo.build(m);
    std::vector<std::pair<int, int>> all_edges;
    for (const auto& [key, cells] : topo.edge_cells) {
      (void)cells;
      const int a = static_cast<int>(key >> 32);
      const int b = static_cast<int>(key & 0xffffffffu);
      all_edges.push_back({a, b});
    }

    std::srand(42);
    std::vector<std::pair<int, int>> chosen;
    for (int i = 0; i < 50 && !all_edges.empty(); ++i) {
      const std::size_t idx = static_cast<std::size_t>(std::rand()) % all_edges.size();
      chosen.push_back(all_edges[idx]);
    }

    SplitResult r = split_edges(m, chosen);
    std::printf("random split: n_split=%d n_skipped=%d\n", r.n_split, r.n_skipped);
    CHECK(r.n_split > 0);
    CHECK(mesh_quality(m).min_q > 0);
    CHECK_NEAR(m.total_volume(), vol0, 1e-12 * vol0);
  }

  std::printf("remesh tests passed\n");
  return 0;
}
