#include "cprocess/topology.hpp"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#else
static inline int omp_get_thread_num()  { return 0; }
static inline int omp_get_max_threads() { return 1; }
#endif

namespace cp {

namespace {
// The six edges of a tet as local-vertex index pairs.
constexpr int kEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
}  // namespace

// Parallelization (PA-1): the per-cell work below (computing 4 node-incidence
// pairs and 6 edge-incidence pairs) is embarrassingly parallel and touches no
// shared state, so it is done in a first pass with #pragma omp parallel for
// schedule(static) into thread-local buffers. A second, strictly sequential
// pass then merges those buffers in fixed thread order (0, 1, 2, ...); since
// schedule(static) hands out contiguous cell-index chunks to threads in
// increasing order, this fixed-order merge reproduces the same cell-ascending
// order per edge/node as the old serial scan, for any thread count -- so
// edge_cells[key] and node_cells[node] are bit-identical between 1 and 4
// threads. node_adj's insertion order depends on unordered_map iteration
// order, which is not itself guaranteed parallel-invariant, so node_adj is
// explicitly sorted at the end to normalize it (this also matches the
// contract that edge_incident() callers only key-lookup edge_cells and never
// rely on edge_cells' own iteration order).
void MeshTopology::build(const Mesh& m) {
  const int nn = static_cast<int>(m.nodes.size());
  const int nc = static_cast<int>(m.cells.size());

  edge_cells.clear();
  edge_cells.reserve(static_cast<std::size_t>(nc) * 6);
  node_cells.assign(nn, {});
  node_adj.assign(nn, {});
  node_boundary.assign(nn, 0);
  node_interface.assign(nn, 0);

  const int nth = std::max(1, omp_get_max_threads());

  // Pass 1 (parallel): gather (edge_key, cell) and (node, cell) pairs into
  // thread-local buffers. schedule(static) gives each thread a contiguous
  // range of cell indices, so within a buffer cell ids are ascending, and
  // buffer 0's cells all precede buffer 1's, etc.
  std::vector<std::vector<std::pair<std::uint64_t, int>>> tl_edges(nth);
  std::vector<std::vector<std::pair<int, int>>> tl_nodes(nth);
#pragma omp parallel for schedule(static)
  for (int ci = 0; ci < nc; ++ci) {
    const auto& c = m.cells[ci];
    const int tid = omp_get_thread_num();
    for (int v = 0; v < 4; ++v) tl_nodes[tid].push_back({c[v], ci});
    for (const auto& e : kEdges)
      tl_edges[tid].push_back({edge_key(c[e[0]], c[e[1]]), ci});
  }

  // Pass 2 (sequential merge, fixed thread order): reproduces the serial
  // cell-ascending insertion order regardless of thread count.
  for (int t = 0; t < nth; ++t)
    for (const auto& [key, ci] : tl_edges[t]) edge_cells[key].push_back(ci);
  for (int t = 0; t < nth; ++t)
    for (const auto& [node, ci] : tl_nodes[t]) node_cells[node].push_back(ci);

  // Node adjacency from unique edges; sorted below to normalize any
  // dependence on unordered_map iteration order.
  for (const auto& [key, cells] : edge_cells) {
    (void)cells;
    const int a = static_cast<int>(key >> 32);
    const int b = static_cast<int>(key & 0xffffffffu);
    node_adj[a].push_back(b);
    node_adj[b].push_back(a);
  }
#pragma omp parallel for schedule(static)
  for (int i = 0; i < nn; ++i)
    std::sort(node_adj[i].begin(), node_adj[i].end());

  // Boundary nodes: any node on a face with no neighbour. Each face writes
  // disjoint-or-idempotent flag bits (node_boundary[v] = 1), so this is safe
  // to parallelize over faces even though multiple faces may touch the same
  // node.
  const int nf = static_cast<int>(m.faces.size());
#pragma omp parallel for schedule(static)
  for (int fi = 0; fi < nf; ++fi) {
    const auto& f = m.faces[fi];
    if (f.neigh < 0)
      for (int v : f.n) node_boundary[v] = 1;
  }

  // Interface nodes: incident cells span more than one region tag. Each
  // iteration writes only node_interface[i], so this is safe per-node.
#pragma omp parallel for schedule(static)
  for (int i = 0; i < nn; ++i) {
    if (node_cells[i].empty()) continue;
    const int r0 = m.cell_region[node_cells[i].front()];
    for (int ci : node_cells[i])
      if (m.cell_region[ci] != r0) { node_interface[i] = 1; break; }
  }
}

const std::vector<int>& MeshTopology::edge_incident(int a, int b) const {
  static const std::vector<int> kEmpty;
  auto it = edge_cells.find(edge_key(a, b));
  return it == edge_cells.end() ? kEmpty : it->second;
}

}  // namespace cp
