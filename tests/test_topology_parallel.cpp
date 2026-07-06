// Tests for PA-1: OpenMP parallelization of MeshTopology::build. The
// requirement is bit-identical results between 1 and 4 threads for every
// derived structure (edge_cells, node_cells, node_adj, node_boundary,
// node_interface), since the parallel merge discipline (thread-local
// collection + fixed thread-order sequential merge, plus an explicit
// node_adj sort) is designed to be thread-count independent.

#include <chrono>
#include <cstdio>

#include "cprocess/mesh.hpp"
#include "cprocess/topology.hpp"
#include "test_util.hpp"

#ifdef _OPENMP
#  include <omp.h>
#else
static inline void omp_set_num_threads(int) {}
#endif

using namespace cp;

int main() {
  Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 24, 24, 24);

  // ---- 1: determinism of MeshTopology::build across thread counts ---------
  {
    omp_set_num_threads(1);
    auto t0 = std::chrono::steady_clock::now();
    MeshTopology t1;
    t1.build(m);
    auto t1_end = std::chrono::steady_clock::now();

    omp_set_num_threads(4);
    auto t2 = std::chrono::steady_clock::now();
    MeshTopology t4;
    t4.build(m);
    auto t4_end = std::chrono::steady_clock::now();
    omp_set_num_threads(1);

    const double ms1 = std::chrono::duration<double, std::milli>(t1_end - t0).count();
    const double ms4 = std::chrono::duration<double, std::milli>(t4_end - t2).count();
    std::printf("MeshTopology::build 24^3: 1thr=%.3fms 4thr=%.3fms (info only)\n",
                ms1, ms4);

    // Edge count and per-edge incident-cell lists match exactly.
    CHECK(t1.edge_cells.size() == t4.edge_cells.size());
    for (const auto& [key, cells1] : t1.edge_cells) {
      auto it = t4.edge_cells.find(key);
      CHECK(it != t4.edge_cells.end());
      CHECK(cells1 == it->second);
    }

    CHECK(t1.node_cells.size() == t4.node_cells.size());
    CHECK(t1.node_adj.size() == t4.node_adj.size());
    CHECK(t1.node_boundary.size() == t4.node_boundary.size());
    CHECK(t1.node_interface.size() == t4.node_interface.size());

    for (std::size_t i = 0; i < t1.node_cells.size(); ++i) {
      CHECK(t1.node_cells[i] == t4.node_cells[i]);
      CHECK(t1.node_adj[i] == t4.node_adj[i]);
      CHECK(t1.node_boundary[i] == t4.node_boundary[i]);
      CHECK(t1.node_interface[i] == t4.node_interface[i]);
    }
  }

  // ---- 2: node_adj is sorted (order-normalized regardless of thread count) -
  {
    MeshTopology topo;
    topo.build(m);
    for (const auto& adj : topo.node_adj) {
      for (std::size_t i = 1; i < adj.size(); ++i) CHECK(adj[i - 1] <= adj[i]);
    }
  }

  // ---- 3: sanity -- edge_incident() lookups still work after parallel build
  {
    MeshTopology topo;
    topo.build(m);
    CHECK(!topo.edge_cells.empty());
    // Every cell contributes 6 (edge, cell) incidences.
    std::size_t total = 0;
    for (const auto& [key, cells] : topo.edge_cells) { (void)key; total += cells.size(); }
    CHECK(total == m.cells.size() * 6);
  }

  std::printf("topology_parallel tests passed\n");
  return 0;
}
